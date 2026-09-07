#!/usr/bin/env python3
"""
PSE_NVAUTO_SEGMENTATION_FINAL.py

Publication-facing implementation of the FINAL acute infarct segmentation
pipeline used in the PSE thalamic diaschisis pilot.

Algorithm
---------
The segmentation model is NVAUTO as distributed through the
``stroke_segmentor`` Python package. This is NOT the complete DeepISLES
ensemble; historical development folders/scripts in the local project used
"deepisles" in some names, but the final model used here is NVAUTO.

Final preprocessing frozen for the pilot
-----------------------------------------
1. Load high-b DWI and ADC after rigid diffusion->T1 coregistration.
2. Use the DWI grid as the reference grid.
3. Resample ADC linearly to the exact DWI grid and enforce identical
   qform/sform geometry for both channels.
4. Project the FreeSurfer native-T1 aseg to the DWI grid with nearest-neighbour
   interpolation and derive a conservative brain mask.
5. Set extracerebral DWI/ADC signal to zero; clip negative intensities to zero.
6. Run NVAUTO/stroke_segmentor on CPU.
7. Save the binary lesion mask, lesion volume, preprocessing audit values and
   a DWI/ADC visual-QC montage.

The original DWI, ADC and FreeSurfer files are never modified.

Usage
-----
    python PSE_NVAUTO_SEGMENTATION_FINAL.py all
    python PSE_NVAUTO_SEGMENTATION_FINAL.py P003
    python PSE_NVAUTO_SEGMENTATION_FINAL.py sub-P003 --overwrite

Environment
-----------
PSE_ROOT may be set to the project root. If unset, ~/PSE is used.

Expected inputs
---------------
<PSE_ROOT>/derivatives/coreg_diffusion_v02/sub-P00X/acute/
    DWI_b1000_coreg_T1.nii
    ADC_coreg_T1.nii

<PSE_ROOT>/derivatives/roi_masks_native/sub-P00X/acute/
    aseg_in_native_T1.nii.gz

Outputs
-------
<PSE_ROOT>/derivatives/lesion_nvauto_inputs/sub-P00X/acute/
    DWI_for_nvauto.nii.gz
    ADC_for_nvauto.nii.gz

<PSE_ROOT>/derivatives/lesion_nvauto_brainmasked/sub-P00X/acute/
    brain_mask_in_DWI.nii.gz
    DWI_brainmasked.nii.gz
    ADC_brainmasked.nii.gz
    lesion_nvauto_brainmasked.nii.gz
    sub-P00X_nvauto_brainmasked_QC.png

Global status table:
<PSE_ROOT>/derivatives/lesion_nvauto_brainmasked/
    PSE_NVAUTO_FINAL_STATUS.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import traceback
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy import ndimage
from stroke_segmentor.inferer import Inferer

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit("matplotlib is required") from exc


ROOT = Path(os.environ.get("PSE_ROOT", str(Path.home() / "PSE"))).expanduser().resolve()
INPUT_ROOT = ROOT / "derivatives" / "coreg_diffusion_v02"
ASEG_ROOT = ROOT / "derivatives" / "roi_masks_native"
HARM_ROOT = ROOT / "derivatives" / "lesion_nvauto_inputs"
OUT_ROOT = ROOT / "derivatives" / "lesion_nvauto_brainmasked"
STATUS_CSV = OUT_ROOT / "PSE_NVAUTO_FINAL_STATUS.csv"


def normalize_subject(arg: str) -> str:
    arg = arg.strip()
    if arg.startswith("sub-P"):
        return arg
    if arg.startswith("P"):
        return f"sub-{arg}"
    raise ValueError(f"Unrecognized subject format: {arg}")


def discover_subjects() -> list[str]:
    subjects = []
    if INPUT_ROOT.is_dir():
        for p in sorted(INPUT_ROOT.glob("sub-P*")):
            if (p / "acute").is_dir():
                subjects.append(p.name)
    return subjects


def save_like_ref(data: np.ndarray, ref_img, path: Path, dtype=np.float32) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    hdr = ref_img.header.copy()
    hdr.set_data_dtype(dtype)
    out = nib.Nifti1Image(np.asarray(data, dtype=dtype), ref_img.affine.copy(), header=hdr)
    out.set_qform(ref_img.affine.copy(), code=1)
    out.set_sform(ref_img.affine.copy(), code=1)
    nib.save(out, str(path))


def harmonize_adc_to_dwi(dwi_img, adc_img) -> tuple[np.ndarray, np.ndarray]:
    dwi = dwi_img.get_fdata(dtype=np.float32)
    same_shape = tuple(adc_img.shape) == tuple(dwi_img.shape)
    same_affine = np.allclose(adc_img.affine, dwi_img.affine, atol=1e-5, rtol=0)

    if same_shape and same_affine:
        adc = adc_img.get_fdata(dtype=np.float32)
    else:
        adc_rs = resample_from_to(
            adc_img,
            (dwi_img.shape, dwi_img.affine),
            order=1,
            mode="constant",
            cval=0.0,
        )
        adc = adc_rs.get_fdata(dtype=np.float32)

    return dwi, adc


def build_brain_mask(aseg_img, ref_img) -> np.ndarray:
    aseg_res = resample_from_to(
        aseg_img,
        (ref_img.shape, ref_img.affine),
        order=0,
        mode="constant",
        cval=0.0,
    )
    aseg = np.rint(aseg_res.get_fdata(dtype=np.float32)).astype(np.int32)
    mask = aseg > 0

    # In-plane only: diffusion slices are relatively thick.
    structure_2d = np.array(
        [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
        dtype=bool,
    )
    cleaned = np.zeros_like(mask, dtype=bool)

    for z in range(mask.shape[2]):
        sl = mask[:, :, z]
        if not np.any(sl):
            continue
        sl = ndimage.binary_fill_holes(sl)
        sl = ndimage.binary_closing(sl, structure=structure_2d, iterations=1)
        sl = ndimage.binary_dilation(sl, structure=structure_2d, iterations=1)
        cleaned[:, :, z] = sl

    labels, n_labels = ndimage.label(cleaned)
    if n_labels > 0:
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        keep = np.where(sizes >= 100)[0]
        cleaned = np.isin(labels, keep)

    return cleaned


def robust_window(x: np.ndarray) -> tuple[float, float]:
    vals = x[np.isfinite(x) & (x != 0)]
    if vals.size == 0:
        return 0.0, 1.0
    lo, hi = np.percentile(vals, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(vals.min()), float(vals.max())
        if hi <= lo:
            hi = lo + 1.0
    return float(lo), float(hi)


def lesion_volume_ml(mask: np.ndarray, img) -> float:
    voxel_mm3 = float(abs(np.linalg.det(img.affine[:3, :3])))
    return float(np.count_nonzero(mask) * voxel_mm3 / 1000.0)


def pick_slices(mask: np.ndarray, n: int = 6) -> list[int]:
    zz = np.where(np.any(mask, axis=(0, 1)))[0]
    if zz.size == 0:
        z0 = max(0, int(round(mask.shape[2] * 0.20)))
        z1 = min(mask.shape[2] - 1, int(round(mask.shape[2] * 0.80)))
    else:
        z0 = max(0, int(zz.min()) - 1)
        z1 = min(mask.shape[2] - 1, int(zz.max()) + 1)
    return sorted(set(np.linspace(z0, z1, n).round().astype(int).tolist()))


def make_qc(sub: str, dwi: np.ndarray, adc: np.ndarray, lesion: np.ndarray, out_png: Path) -> None:
    zlist = pick_slices(lesion, n=6)
    dwi_w = robust_window(dwi)
    adc_w = robust_window(adc)

    fig, axes = plt.subplots(2, len(zlist), figsize=(3.0 * len(zlist), 6.0))
    if len(zlist) == 1:
        axes = axes.reshape(2, 1)

    for j, z in enumerate(zlist):
        D = np.rot90(dwi[:, :, z])
        A = np.rot90(adc[:, :, z])
        L = np.rot90(lesion[:, :, z])

        axes[0, j].imshow(D, cmap="gray", vmin=dwi_w[0], vmax=dwi_w[1])
        if np.any(L):
            axes[0, j].contour(L.astype(float), levels=[0.5], linewidths=1.4)
        axes[0, j].set_title(f"DWI z={z}")
        axes[0, j].axis("off")

        axes[1, j].imshow(A, cmap="gray", vmin=adc_w[0], vmax=adc_w[1])
        if np.any(L):
            axes[1, j].contour(L.astype(float), levels=[0.5], linewidths=1.4)
        axes[1, j].set_title("ADC + NVAUTO")
        axes[1, j].axis("off")

    fig.suptitle(f"{sub}: final brain-masked NVAUTO segmentation")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def resample_binary_to_ref(img, ref_img) -> np.ndarray:
    same_shape = tuple(img.shape) == tuple(ref_img.shape)
    same_affine = np.allclose(img.affine, ref_img.affine, atol=1e-5, rtol=0)
    if same_shape and same_affine:
        data = img.get_fdata(dtype=np.float32)
    else:
        res = resample_from_to(
            img,
            (ref_img.shape, ref_img.affine),
            order=0,
            mode="constant",
            cval=0.0,
        )
        data = res.get_fdata(dtype=np.float32)
    return data > 0.5


def process_subject(sub: str, overwrite: bool) -> dict[str, str]:
    input_dir = INPUT_ROOT / sub / "acute"
    aseg_dir = ASEG_ROOT / sub / "acute"
    harm_dir = HARM_ROOT / sub / "acute"
    out_dir = OUT_ROOT / sub / "acute"
    harm_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    dwi_file = input_dir / "DWI_b1000_coreg_T1.nii"
    adc_file = input_dir / "ADC_coreg_T1.nii"
    aseg_file = aseg_dir / "aseg_in_native_T1.nii.gz"
    if not aseg_file.is_file():
        aseg_file = aseg_dir / "aseg_in_native_T1.nii"

    dwi_harm_file = harm_dir / "DWI_for_nvauto.nii.gz"
    adc_harm_file = harm_dir / "ADC_for_nvauto.nii.gz"
    brain_file = out_dir / "brain_mask_in_DWI.nii.gz"
    dwi_masked_file = out_dir / "DWI_brainmasked.nii.gz"
    adc_masked_file = out_dir / "ADC_brainmasked.nii.gz"
    seg_file = out_dir / "lesion_nvauto_brainmasked.nii.gz"
    qc_file = out_dir / f"{sub}_nvauto_brainmasked_QC.png"

    row = {
        "Subject": sub,
        "Status": "",
        "LesionVolume_mL": "",
        "BrainMaskNvox": "",
        "BrainMaskFractionOfGrid": "",
        "DWIOutsideNonzeroBefore": "",
        "ADCOutsideNonzeroBefore": "",
        "QCPNG": str(qc_file),
        "Message": "",
    }

    try:
        for f in (dwi_file, adc_file, aseg_file):
            if not f.is_file():
                raise FileNotFoundError(str(f))

        if seg_file.is_file() and not overwrite:
            seg_img = nib.load(str(seg_file))
            lesion = seg_img.get_fdata(dtype=np.float32) > 0.5
            row["LesionVolume_mL"] = f"{lesion_volume_ml(lesion, seg_img):.6f}"
            row["Status"] = "SKIPPED_EXISTING"
            row["Message"] = "Existing final NVAUTO mask reused."
            return row

        dwi_img = nib.load(str(dwi_file))
        adc_img = nib.load(str(adc_file))
        aseg_img = nib.load(str(aseg_file))

        dwi, adc = harmonize_adc_to_dwi(dwi_img, adc_img)
        save_like_ref(dwi, dwi_img, dwi_harm_file, np.float32)
        save_like_ref(adc, dwi_img, adc_harm_file, np.float32)

        brain = build_brain_mask(aseg_img, dwi_img)
        row["BrainMaskNvox"] = str(int(np.count_nonzero(brain)))
        row["BrainMaskFractionOfGrid"] = f"{np.count_nonzero(brain) / brain.size:.8f}"
        row["DWIOutsideNonzeroBefore"] = str(
            int(np.count_nonzero((~brain) & np.isfinite(dwi) & (dwi != 0)))
        )
        row["ADCOutsideNonzeroBefore"] = str(
            int(np.count_nonzero((~brain) & np.isfinite(adc) & (adc != 0)))
        )

        dwi_clean = np.where(brain & np.isfinite(dwi), dwi, 0.0).astype(np.float32)
        adc_clean = np.where(brain & np.isfinite(adc), adc, 0.0).astype(np.float32)
        dwi_clean[dwi_clean < 0] = 0
        adc_clean[adc_clean < 0] = 0

        save_like_ref(brain.astype(np.uint8), dwi_img, brain_file, np.uint8)
        save_like_ref(dwi_clean, dwi_img, dwi_masked_file, np.float32)
        save_like_ref(adc_clean, dwi_img, adc_masked_file, np.float32)

        if overwrite and seg_file.exists():
            seg_file.unlink()

        inferer = Inferer(force_cpu=True)
        inferer.infer(
            adc_path=str(adc_masked_file),
            dwi_path=str(dwi_masked_file),
            segmentation_path=str(seg_file),
        )
        if not seg_file.is_file():
            raise RuntimeError("NVAUTO returned without creating a segmentation.")

        lesion = resample_binary_to_ref(nib.load(str(seg_file)), dwi_img)
        volume = lesion_volume_ml(lesion, dwi_img)
        row["LesionVolume_mL"] = f"{volume:.6f}"
        make_qc(sub, dwi_clean, adc_clean, lesion, qc_file)

        row["Status"] = "OK"
        row["Message"] = "Final brain-masked NVAUTO inference completed."

    except Exception as exc:
        row["Status"] = "ERROR"
        row["Message"] = f"{type(exc).__name__}: {exc}"
        print(traceback.format_exc(limit=5))

    return row


def write_status(rows: list[dict[str, str]]) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    fields = [
        "Subject",
        "Status",
        "LesionVolume_mL",
        "BrainMaskNvox",
        "BrainMaskFractionOfGrid",
        "DWIOutsideNonzeroBefore",
        "ADCOutsideNonzeroBefore",
        "QCPNG",
        "Message",
    ]
    with STATUS_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("subject", help="P003, sub-P003, or all")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.subject.lower() == "all":
        subjects = discover_subjects()
        if not subjects:
            raise SystemExit(f"No subjects found under {INPUT_ROOT}")
    else:
        subjects = [normalize_subject(args.subject)]

    rows = []
    print("=" * 76)
    print("PSE FINAL NVAUTO LESION SEGMENTATION")
    print(f"Project root: {ROOT}")
    print("=" * 76)

    for sub in subjects:
        print(f"\n[{sub}]")
        row = process_subject(sub, args.overwrite)
        rows.append(row)
        print(
            f"  {row['Status']} | lesion volume={row['LesionVolume_mL']} mL | "
            f"{row['Message']}"
        )

    write_status(rows)
    print(f"\nStatus CSV: {STATUS_CSV}")
    return 1 if any(r["Status"] == "ERROR" for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
