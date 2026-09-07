#!/usr/bin/env python3
"""
PSE_FINALIZE_LESIONS_AND_ROI_OVERLAP_01.py

Finalize the 7 brain-masked NVAUTO lesion masks and quantify direct lesion
overlap with left/right thalamus and hippocampus.

This script is intended for the current PSE acute pilot.

Pipeline frozen here
--------------------
DWI b1000 + ADC
    -> ADC harmonized to exact DWI grid
    -> FreeSurfer-derived brain mask
    -> extracerebral signal set to zero
    -> NVAUTO / stroke_segmentor
    -> visual QC
    -> final binary lesion mask
    -> lesion-to-deep-ROI overlap

Important methodological choice
-------------------------------
The FINAL lesion volume is measured in the DWI grid.

For thalamus/hippocampus overlap, the binary lesion mask is resampled with
nearest-neighbour interpolation INTO native T1/aseg space. Overlap is then
computed against the original FreeSurfer aseg labels:
    Left thalamus      10
    Left hippocampus   17
    Right thalamus     49
    Right hippocampus  53

This avoids downsampling the anatomical ROIs to the relatively coarse DWI
grid. Only a BINARY lesion mask is resampled; no quantitative PET/ASL/CTP
image is modified.

Inputs
------
Brain-masked NVAUTO:
  <PSE_ROOT>/derivatives/lesion_nvauto_brainmasked/
    sub-P00X/acute/lesion_nvauto_brainmasked.nii.gz
    sub-P00X/acute/DWI_brainmasked.nii.gz

Native FreeSurfer aseg:
  <PSE_ROOT>/derivatives/roi_masks_native/
    sub-P00X/acute/aseg_in_native_T1.nii.gz

Stroke side:
  <PSE_ROOT>/derivatives/roi_values/final/
    PSE_STROKE_SIDE_TEMPLATE.csv

Outputs
-------
Final lesion masks:
  <PSE_ROOT>/derivatives/lesion_masks_final/
    sub-P00X/acute/lesion_final_in_DWI.nii
    sub-P00X/acute/sub-P00X_final_lesion_deepROI_QC.png

Tables:
  <PSE_ROOT>/derivatives/lesion_masks_final/
    PSE_FINAL_LESION_STATUS.csv
    PSE_LESION_ROI_OVERLAP_LONG_v01.csv
    PSE_LESION_ROI_OVERLAP_IPSICONTRA_v01.csv
    PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv

Usage
-----
    conda activate pse-stroke
    cd <PSE_ROOT>/code
    python PSE_FINALIZE_LESIONS_AND_ROI_OVERLAP_01.py

Dependencies
------------
numpy, scipy, nibabel, matplotlib
(already present in the current pse-stroke environment)

Safety
------
Existing final masks are NOT overwritten by default.
To intentionally regenerate them, set OVERWRITE_FINAL = True below.
"""

from __future__ import annotations

import os

import csv
import shutil
import sys
import traceback
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required. Install it with:\n"
        "  conda install matplotlib -y"
    ) from exc


# ============================================================================
# Configuration

ROOT = Path(os.environ.get("PSE_ROOT", str(Path.home() / "PSE"))).expanduser().resolve()

SOURCE_ROOT = ROOT / "derivatives" / "lesion_nvauto_brainmasked"
ASEG_ROOT = ROOT / "derivatives" / "roi_masks_native"
SIDE_FILE = (
    ROOT
    / "derivatives"
    / "roi_values"
    / "final"
    / "PSE_STROKE_SIDE_TEMPLATE.csv"
)

FINAL_ROOT = ROOT / "derivatives" / "lesion_masks_final"

OVERWRITE_FINAL = False

SUBJECTS = [f"sub-P{i:03d}" for i in range(1, 8)]

# Visual QC decisions from the completed pilot review.
# These labels document the decision; the script does not infer QC quality.
QC_DECISIONS = {
    "sub-P001": "ACCEPTED_PILOT_WITH_RESERVE",
    "sub-P002": "ACCEPTED_PILOT",
    "sub-P003": "ACCEPTED_PILOT",
    "sub-P004": "ACCEPTED_PILOT",
    "sub-P005": "ACCEPTED_PILOT",
    "sub-P006": "ACCEPTED_PILOT",
    "sub-P007": "ACCEPTED_PILOT",
}

ASEG_LABELS = {
    ("Thalamus", "L"): 10,
    ("Hippocampus", "L"): 17,
    ("Thalamus", "R"): 49,
    ("Hippocampus", "R"): 53,
}

# We DO NOT force a categorical biological interpretation from one arbitrary
# threshold. Instead, the output contains:
#   - any overlap
#   - >=1% ROI involvement
#   - >=5% ROI involvement
# These are sensitivity descriptors, not definitions of diaschisis.
THRESHOLDS_PCT = (1.0, 5.0)


# ============================================================================
# Helpers

def read_stroke_sides(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Stroke-side CSV not found: {path}")

    out: dict[str, str] = {}

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames:
            raise RuntimeError("Stroke-side CSV has no header.")

        # Tolerate minor capitalization/naming variation.
        field_map = {str(x).strip().lower(): x for x in reader.fieldnames}

        subject_field = field_map.get("subject")
        side_field = field_map.get("strokeside") or field_map.get("stroke_side")

        if subject_field is None or side_field is None:
            raise RuntimeError(
                "Stroke-side CSV must contain Subject and StrokeSide columns. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            sub = str(row.get(subject_field, "")).strip()
            side = str(row.get(side_field, "")).strip().upper()

            if not sub:
                continue

            if not sub.startswith("sub-") and sub.startswith("P"):
                sub = f"sub-{sub}"

            if side in {"L", "R"}:
                out[sub] = side

    return out


def find_aseg(sub: str) -> Path:
    base = ASEG_ROOT / sub / "acute"
    gz = base / "aseg_in_native_T1.nii.gz"
    nii = base / "aseg_in_native_T1.nii"

    if gz.is_file():
        return gz
    if nii.is_file():
        return nii

    raise FileNotFoundError(f"Native aseg not found for {sub}: {base}")


def robust_window(x: np.ndarray) -> tuple[float, float]:
    vals = x[np.isfinite(x) & (x != 0)]

    if vals.size == 0:
        return 0.0, 1.0

    lo, hi = np.percentile(vals, [2, 98])

    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.min(vals))
        hi = float(np.max(vals))

        if hi <= lo:
            hi = lo + 1.0

    return float(lo), float(hi)


def voxel_volume_mm3(img) -> float:
    return float(abs(np.linalg.det(img.affine[:3, :3])))


def save_binary_like(
    mask: np.ndarray,
    ref_img,
    out_path: Path,
) -> None:
    header = ref_img.header.copy()
    header.set_data_dtype(np.uint8)

    out = nib.Nifti1Image(
        np.asarray(mask, dtype=np.uint8),
        ref_img.affine.copy(),
        header=header,
    )

    # Preserve exact reference geometry in both forms.
    out.set_qform(ref_img.affine.copy(), code=1)
    out.set_sform(ref_img.affine.copy(), code=1)

    nib.save(out, str(out_path))


def load_binary(path: Path) -> tuple[object, np.ndarray]:
    img = nib.load(str(path))
    data = img.get_fdata(dtype=np.float32)
    mask = np.isfinite(data) & (data > 0.5)
    return img, mask


def resample_binary_to_ref(mask_img, ref_img) -> np.ndarray:
    same_shape = tuple(mask_img.shape) == tuple(ref_img.shape)
    same_affine = np.allclose(
        mask_img.affine,
        ref_img.affine,
        atol=1e-5,
        rtol=0,
    )

    if same_shape and same_affine:
        data = mask_img.get_fdata(dtype=np.float32)
    else:
        res = resample_from_to(
            mask_img,
            (ref_img.shape, ref_img.affine),
            order=0,
            mode="constant",
            cval=0.0,
        )
        data = res.get_fdata(dtype=np.float32)

    return np.isfinite(data) & (data > 0.5)


def resample_label_mask_to_ref(
    aseg_img,
    label: int,
    ref_img,
) -> np.ndarray:
    """
    For DISPLAY ONLY: project a native aseg label to DWI grid.
    Quantitative overlap is computed in native aseg space instead.
    """
    label_native = (
        np.rint(aseg_img.get_fdata(dtype=np.float32)).astype(np.int32) == label
    )

    label_img = nib.Nifti1Image(
        label_native.astype(np.uint8),
        aseg_img.affine.copy(),
    )
    label_img.set_qform(aseg_img.affine.copy(), code=1)
    label_img.set_sform(aseg_img.affine.copy(), code=1)

    return resample_binary_to_ref(label_img, ref_img)


def pick_qc_slices(
    lesion: np.ndarray,
    thal_l: np.ndarray,
    thal_r: np.ndarray,
    hip_l: np.ndarray,
    hip_r: np.ndarray,
    n: int = 6,
) -> list[int]:
    support = lesion | thal_l | thal_r | hip_l | hip_r
    zz = np.where(np.any(support, axis=(0, 1)))[0]

    if zz.size == 0:
        z0 = max(0, int(round(lesion.shape[2] * 0.20)))
        z1 = min(lesion.shape[2] - 1, int(round(lesion.shape[2] * 0.80)))
    else:
        z0 = max(0, int(zz.min()) - 1)
        z1 = min(lesion.shape[2] - 1, int(zz.max()) + 1)

    return sorted(
        set(np.linspace(z0, z1, n).round().astype(int).tolist())
    )


def make_final_qc(
    sub: str,
    dwi: np.ndarray,
    lesion: np.ndarray,
    thal_l: np.ndarray,
    thal_r: np.ndarray,
    hip_l: np.ndarray,
    hip_r: np.ndarray,
    out_png: Path,
    stroke_side: str,
    lesion_volume_ml: float,
) -> None:

    zlist = pick_qc_slices(
        lesion,
        thal_l,
        thal_r,
        hip_l,
        hip_r,
        n=6,
    )

    win = robust_window(dwi)

    ncols = len(zlist)
    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=(3.0 * ncols, 6.1),
    )

    if ncols == 1:
        axes = axes.reshape(2, 1)

    for j, z in enumerate(zlist):
        D = np.rot90(dwi[:, :, z])

        L = np.rot90(lesion[:, :, z])
        TL = np.rot90(thal_l[:, :, z])
        TR = np.rot90(thal_r[:, :, z])
        HL = np.rot90(hip_l[:, :, z])
        HR = np.rot90(hip_r[:, :, z])

        # Row 1: final lesion.
        ax = axes[0, j]
        ax.imshow(D, cmap="gray", vmin=win[0], vmax=win[1])

        if np.any(L):
            ax.contour(
                L.astype(float),
                levels=[0.5],
                linewidths=1.6,
                linestyles="-",
            )

        ax.set_title(f"DWI + final lesion z={z}")
        ax.axis("off")

        # Row 2: final lesion + deep anatomical ROIs.
        ax = axes[1, j]
        ax.imshow(D, cmap="gray", vmin=win[0], vmax=win[1])

        if np.any(L):
            ax.contour(
                L.astype(float),
                levels=[0.5],
                linewidths=1.5,
                linestyles="-",
            )

        if np.any(TL):
            ax.contour(
                TL.astype(float),
                levels=[0.5],
                linewidths=1.0,
                linestyles="--",
            )

        if np.any(TR):
            ax.contour(
                TR.astype(float),
                levels=[0.5],
                linewidths=1.0,
                linestyles="--",
            )

        if np.any(HL):
            ax.contour(
                HL.astype(float),
                levels=[0.5],
                linewidths=1.0,
                linestyles=":",
            )

        if np.any(HR):
            ax.contour(
                HR.astype(float),
                levels=[0.5],
                linewidths=1.0,
                linestyles=":",
            )

        ax.set_title("Lesion solid | thalamus dashed | hippocampus dotted")
        ax.axis("off")

    fig.suptitle(
        f"{sub} | stroke side {stroke_side} | "
        f"final lesion {lesion_volume_ml:.2f} mL",
        fontsize=13,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Main

def main() -> int:
    FINAL_ROOT.mkdir(parents=True, exist_ok=True)

    stroke_sides = read_stroke_sides(SIDE_FILE)

    status_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []

    print("=" * 76)
    print("PSE FINAL LESION MASKS + THALAMUS/HIPPOCAMPUS OVERLAP")
    print("=" * 76)
    print(f"Source root : {SOURCE_ROOT}")
    print(f"Final root  : {FINAL_ROOT}")
    print(f"Stroke side : {SIDE_FILE}")
    print("=" * 76)

    for i, sub in enumerate(SUBJECTS, start=1):

        print()
        print("=" * 76)
        print(f"[{i}/{len(SUBJECTS)}] {sub}")
        print("=" * 76)

        side = stroke_sides.get(sub)

        source_dir = SOURCE_ROOT / sub / "acute"
        source_seg = source_dir / "lesion_nvauto_brainmasked.nii.gz"
        source_dwi = source_dir / "DWI_brainmasked.nii.gz"

        final_dir = FINAL_ROOT / sub / "acute"
        final_dir.mkdir(parents=True, exist_ok=True)

        final_seg = final_dir / "lesion_final_in_DWI.nii"
        qc_png = final_dir / f"{sub}_final_lesion_deepROI_QC.png"

        status = {
            "Subject": sub,
            "StrokeSide": side or "",
            "QDecision": QC_DECISIONS.get(sub, "UNSPECIFIED"),
            "SourceSegmentation": str(source_seg),
            "FinalSegmentation": str(final_seg),
            "LesionVolume_mL_DWIgrid": "",
            "LesionNvox_DWIgrid": "",
            "FinalMaskCreated": False,
            "QCPNG": str(qc_png),
            "Status": "",
            "Message": "",
        }

        try:
            if side not in {"L", "R"}:
                raise RuntimeError(
                    f"Missing/invalid stroke side for {sub}: {side}"
                )

            if not source_seg.is_file():
                raise FileNotFoundError(
                    f"Brain-masked NVAUTO segmentation not found: {source_seg}"
                )

            if not source_dwi.is_file():
                raise FileNotFoundError(
                    f"Brain-masked DWI not found: {source_dwi}"
                )

            aseg_file = find_aseg(sub)

            # ----------------------------------------------------------------
            # Load and validate source final candidate.

            lesion_img, lesion = load_binary(source_seg)

            if lesion.ndim != 3:
                raise RuntimeError(
                    f"Lesion mask must be 3D, got {lesion.shape}"
                )

            if np.count_nonzero(lesion) == 0:
                raise RuntimeError(
                    "Brain-masked NVAUTO lesion is empty; refusing promotion."
                )

            lesion_nvox = int(np.count_nonzero(lesion))
            lesion_vol_ml = (
                lesion_nvox * voxel_volume_mm3(lesion_img) / 1000.0
            )

            status["LesionNvox_DWIgrid"] = lesion_nvox
            status["LesionVolume_mL_DWIgrid"] = f"{lesion_vol_ml:.6f}"

            # ----------------------------------------------------------------
            # Promote to final mask, safely.

            if final_seg.exists() and not OVERWRITE_FINAL:
                # Validate existing final mask instead of silently overwriting.
                existing_img, existing = load_binary(final_seg)

                same_shape = tuple(existing.shape) == tuple(lesion.shape)
                same_affine = np.allclose(
                    existing_img.affine,
                    lesion_img.affine,
                    atol=1e-5,
                    rtol=0,
                )
                same_data = (
                    same_shape
                    and same_affine
                    and np.array_equal(existing, lesion)
                )

                if not same_data:
                    raise RuntimeError(
                        "A different final mask already exists. "
                        "Set OVERWRITE_FINAL=True only if intentional."
                    )

                print("  Existing final mask matches source; kept unchanged.")

            else:
                save_binary_like(lesion, lesion_img, final_seg)
                print(f"  Final mask written: {final_seg}")

            status["FinalMaskCreated"] = True

            # ----------------------------------------------------------------
            # Quantitative overlap in native T1 / aseg space.

            aseg_img = nib.load(str(aseg_file))
            aseg = (
                np.rint(
                    aseg_img.get_fdata(dtype=np.float32)
                ).astype(np.int32)
            )

            # Resample BINARY lesion to native T1 / aseg grid.
            lesion_native = resample_binary_to_ref(
                nib.load(str(final_seg)),
                aseg_img,
            )

            native_voxel_mm3 = voxel_volume_mm3(aseg_img)
            lesion_native_nvox = int(np.count_nonzero(lesion_native))

            subject_structure_rows: dict[
                tuple[str, str], dict[str, object]
            ] = {}

            for (structure, hemi), label in ASEG_LABELS.items():

                roi = aseg == label

                roi_nvox = int(np.count_nonzero(roi))
                overlap = lesion_native & roi
                overlap_nvox = int(np.count_nonzero(overlap))

                roi_vol_ml = roi_nvox * native_voxel_mm3 / 1000.0
                overlap_vol_ml = (
                    overlap_nvox * native_voxel_mm3 / 1000.0
                )

                roi_overlap_frac = (
                    overlap_nvox / roi_nvox if roi_nvox > 0 else np.nan
                )

                lesion_overlap_frac = (
                    overlap_nvox / lesion_native_nvox
                    if lesion_native_nvox > 0
                    else np.nan
                )

                relation = "Ipsi" if hemi == side else "Contra"

                row = {
                    "Subject": sub,
                    "StrokeSide": side,
                    "Structure": structure,
                    "Hemisphere": hemi,
                    "RelationToStroke": relation,
                    "AsegLabel": label,
                    "ROINvox_nativeT1": roi_nvox,
                    "ROIVolume_mL_nativeT1": f"{roi_vol_ml:.6f}",
                    "OverlapNvox_nativeT1": overlap_nvox,
                    "OverlapVolume_mL_nativeT1": f"{overlap_vol_ml:.6f}",
                    "ROIOverlapFraction": (
                        "" if np.isnan(roi_overlap_frac)
                        else f"{roi_overlap_frac:.8f}"
                    ),
                    "ROIOverlapPercent": (
                        "" if np.isnan(roi_overlap_frac)
                        else f"{100.0 * roi_overlap_frac:.6f}"
                    ),
                    "LesionOverlapFraction_nativeT1": (
                        "" if np.isnan(lesion_overlap_frac)
                        else f"{lesion_overlap_frac:.8f}"
                    ),
                    "AnyOverlap": overlap_nvox > 0,
                    "ROIOverlap_GE1pct": (
                        False
                        if np.isnan(roi_overlap_frac)
                        else 100.0 * roi_overlap_frac >= 1.0
                    ),
                    "ROIOverlap_GE5pct": (
                        False
                        if np.isnan(roi_overlap_frac)
                        else 100.0 * roi_overlap_frac >= 5.0
                    ),
                }

                overlap_rows.append(row)
                subject_structure_rows[(structure, relation)] = row

            # ----------------------------------------------------------------
            # Ipsi/contra paired summary.

            for structure in ("Thalamus", "Hippocampus"):

                ipsi = subject_structure_rows[(structure, "Ipsi")]
                contra = subject_structure_rows[(structure, "Contra")]

                pair_rows.append({
                    "Subject": sub,
                    "StrokeSide": side,
                    "Structure": structure,

                    "IpsiHemisphere": ipsi["Hemisphere"],
                    "IpsiOverlapVolume_mL": ipsi[
                        "OverlapVolume_mL_nativeT1"
                    ],
                    "IpsiROIOverlapPercent": ipsi["ROIOverlapPercent"],
                    "IpsiAnyOverlap": ipsi["AnyOverlap"],
                    "IpsiROIOverlap_GE1pct": ipsi[
                        "ROIOverlap_GE1pct"
                    ],
                    "IpsiROIOverlap_GE5pct": ipsi[
                        "ROIOverlap_GE5pct"
                    ],

                    "ContraHemisphere": contra["Hemisphere"],
                    "ContraOverlapVolume_mL": contra[
                        "OverlapVolume_mL_nativeT1"
                    ],
                    "ContraROIOverlapPercent": contra[
                        "ROIOverlapPercent"
                    ],
                    "ContraAnyOverlap": contra["AnyOverlap"],
                    "ContraROIOverlap_GE1pct": contra[
                        "ROIOverlap_GE1pct"
                    ],
                    "ContraROIOverlap_GE5pct": contra[
                        "ROIOverlap_GE5pct"
                    ],
                })

            # ----------------------------------------------------------------
            # QC figure in DWI grid.

            dwi_img = nib.load(str(source_dwi))

            if (
                tuple(dwi_img.shape) != tuple(lesion_img.shape)
                or not np.allclose(
                    dwi_img.affine,
                    lesion_img.affine,
                    atol=1e-5,
                    rtol=0,
                )
            ):
                # Resample DWI only for display if geometry somehow differs.
                dwi_display_img = resample_from_to(
                    dwi_img,
                    (lesion_img.shape, lesion_img.affine),
                    order=1,
                    mode="constant",
                    cval=0.0,
                )
            else:
                dwi_display_img = dwi_img

            dwi = dwi_display_img.get_fdata(dtype=np.float32)

            thal_l = resample_label_mask_to_ref(
                aseg_img,
                10,
                lesion_img,
            )
            hip_l = resample_label_mask_to_ref(
                aseg_img,
                17,
                lesion_img,
            )
            thal_r = resample_label_mask_to_ref(
                aseg_img,
                49,
                lesion_img,
            )
            hip_r = resample_label_mask_to_ref(
                aseg_img,
                53,
                lesion_img,
            )

            make_final_qc(
                sub=sub,
                dwi=dwi,
                lesion=lesion,
                thal_l=thal_l,
                thal_r=thal_r,
                hip_l=hip_l,
                hip_r=hip_r,
                out_png=qc_png,
                stroke_side=side,
                lesion_volume_ml=lesion_vol_ml,
            )

            status["Status"] = "OK"
            status["Message"] = (
                "Final lesion promoted and deep-ROI overlap computed."
            )

            print(f"  Lesion volume: {lesion_vol_ml:.3f} mL")
            print(f"  QC: {qc_png}")

        except Exception as exc:
            status["Status"] = "ERROR"
            status["Message"] = f"{type(exc).__name__}: {exc}"

            print(f"  ERROR: {status['Message']}")
            print(traceback.format_exc(limit=4))

        status_rows.append(status)

    # =========================================================================
    # Write tables.

    status_fields = [
        "Subject",
        "StrokeSide",
        "QDecision",
        "SourceSegmentation",
        "FinalSegmentation",
        "LesionVolume_mL_DWIgrid",
        "LesionNvox_DWIgrid",
        "FinalMaskCreated",
        "QCPNG",
        "Status",
        "Message",
    ]

    overlap_fields = [
        "Subject",
        "StrokeSide",
        "Structure",
        "Hemisphere",
        "RelationToStroke",
        "AsegLabel",
        "ROINvox_nativeT1",
        "ROIVolume_mL_nativeT1",
        "OverlapNvox_nativeT1",
        "OverlapVolume_mL_nativeT1",
        "ROIOverlapFraction",
        "ROIOverlapPercent",
        "LesionOverlapFraction_nativeT1",
        "AnyOverlap",
        "ROIOverlap_GE1pct",
        "ROIOverlap_GE5pct",
    ]

    pair_fields = [
        "Subject",
        "StrokeSide",
        "Structure",

        "IpsiHemisphere",
        "IpsiOverlapVolume_mL",
        "IpsiROIOverlapPercent",
        "IpsiAnyOverlap",
        "IpsiROIOverlap_GE1pct",
        "IpsiROIOverlap_GE5pct",

        "ContraHemisphere",
        "ContraOverlapVolume_mL",
        "ContraROIOverlapPercent",
        "ContraAnyOverlap",
        "ContraROIOverlap_GE1pct",
        "ContraROIOverlap_GE5pct",
    ]

    status_csv = FINAL_ROOT / "PSE_FINAL_LESION_STATUS.csv"
    overlap_csv = FINAL_ROOT / "PSE_LESION_ROI_OVERLAP_LONG_v01.csv"
    pair_csv = (
        FINAL_ROOT
        / "PSE_LESION_ROI_OVERLAP_IPSICONTRA_v01.csv"
    )

    with status_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=status_fields)
        writer.writeheader()
        writer.writerows(status_rows)

    with overlap_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=overlap_fields)
        writer.writeheader()
        writer.writerows(overlap_rows)

    with pair_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=pair_fields)
        writer.writeheader()
        writer.writerows(pair_rows)

    # =========================================================================
    # Wide subject-level direct-involvement summary.

    direct_rows: list[dict[str, object]] = []

    for sub in SUBJECTS:
        side = stroke_sides.get(sub, "")

        thal = next(
            (
                x
                for x in pair_rows
                if x["Subject"] == sub
                and x["Structure"] == "Thalamus"
            ),
            None,
        )

        hip = next(
            (
                x
                for x in pair_rows
                if x["Subject"] == sub
                and x["Structure"] == "Hippocampus"
            ),
            None,
        )

        status = next(
            (x for x in status_rows if x["Subject"] == sub),
            {},
        )

        if thal is None or hip is None:
            direct_rows.append({
                "Subject": sub,
                "StrokeSide": side,
                "QDecision": QC_DECISIONS.get(sub, ""),
                "LesionVolume_mL_DWIgrid": status.get(
                    "LesionVolume_mL_DWIgrid",
                    "",
                ),
                "Thalamus_IpsiOverlapPercent": "",
                "Thalamus_IpsiAnyOverlap": "",
                "Thalamus_Ipsi_GE1pct": "",
                "Thalamus_Ipsi_GE5pct": "",
                "Hippocampus_IpsiOverlapPercent": "",
                "Hippocampus_IpsiAnyOverlap": "",
                "Hippocampus_Ipsi_GE1pct": "",
                "Hippocampus_Ipsi_GE5pct": "",
                "EitherDeepROI_IpsiAnyOverlap": "",
                "EitherDeepROI_Ipsi_GE1pct": "",
                "EitherDeepROI_Ipsi_GE5pct": "",
                "Status": status.get("Status", "ERROR"),
            })
            continue

        th_any = bool(thal["IpsiAnyOverlap"])
        th_1 = bool(thal["IpsiROIOverlap_GE1pct"])
        th_5 = bool(thal["IpsiROIOverlap_GE5pct"])

        hi_any = bool(hip["IpsiAnyOverlap"])
        hi_1 = bool(hip["IpsiROIOverlap_GE1pct"])
        hi_5 = bool(hip["IpsiROIOverlap_GE5pct"])

        direct_rows.append({
            "Subject": sub,
            "StrokeSide": side,
            "QDecision": QC_DECISIONS.get(sub, ""),
            "LesionVolume_mL_DWIgrid": status.get(
                "LesionVolume_mL_DWIgrid",
                "",
            ),

            "Thalamus_IpsiOverlapPercent": thal[
                "IpsiROIOverlapPercent"
            ],
            "Thalamus_IpsiAnyOverlap": th_any,
            "Thalamus_Ipsi_GE1pct": th_1,
            "Thalamus_Ipsi_GE5pct": th_5,

            "Hippocampus_IpsiOverlapPercent": hip[
                "IpsiROIOverlapPercent"
            ],
            "Hippocampus_IpsiAnyOverlap": hi_any,
            "Hippocampus_Ipsi_GE1pct": hi_1,
            "Hippocampus_Ipsi_GE5pct": hi_5,

            "EitherDeepROI_IpsiAnyOverlap": th_any or hi_any,
            "EitherDeepROI_Ipsi_GE1pct": th_1 or hi_1,
            "EitherDeepROI_Ipsi_GE5pct": th_5 or hi_5,

            "Status": status.get("Status", ""),
        })

    direct_fields = [
        "Subject",
        "StrokeSide",
        "QDecision",
        "LesionVolume_mL_DWIgrid",

        "Thalamus_IpsiOverlapPercent",
        "Thalamus_IpsiAnyOverlap",
        "Thalamus_Ipsi_GE1pct",
        "Thalamus_Ipsi_GE5pct",

        "Hippocampus_IpsiOverlapPercent",
        "Hippocampus_IpsiAnyOverlap",
        "Hippocampus_Ipsi_GE1pct",
        "Hippocampus_Ipsi_GE5pct",

        "EitherDeepROI_IpsiAnyOverlap",
        "EitherDeepROI_Ipsi_GE1pct",
        "EitherDeepROI_Ipsi_GE5pct",

        "Status",
    ]

    direct_csv = (
        FINAL_ROOT
        / "PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv"
    )

    with direct_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=direct_fields)
        writer.writeheader()
        writer.writerows(direct_rows)

    # =========================================================================
    # Console summary.

    n_ok = sum(x["Status"] == "OK" for x in status_rows)

    print()
    print("=" * 76)
    print("FINALIZATION COMPLETE")
    print("=" * 76)
    print(f"Successful subjects: {n_ok}/{len(SUBJECTS)}")
    print(f"Status:      {status_csv}")
    print(f"Long overlap:{overlap_csv}")
    print(f"Ipsi/contra: {pair_csv}")
    print(f"Direct ROI:  {direct_csv}")
    print()
    print(
        "Interpretation reminder: AnyOverlap / >=1% / >=5% are "
        "descriptive sensitivity flags. The script does NOT label "
        "a subject as 'diaschisis'."
    )
    print("=" * 76)

    return 0 if n_ok == len(SUBJECTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
