#!/usr/bin/env python3
"""
PSE_JHU_WM_TOPOGRAPHY_01.py

Exploratory white-matter lesion-topography analysis using the
JHU ICBM-DTI-81 white-matter labels atlas.

Pipeline
--------
For each acute subject:
1) Load FreeSurfer brain.mgz.
2) Resample the FINAL binary ischemic lesion mask to the FreeSurfer brain grid
   with nearest-neighbour interpolation.
3) Register native FreeSurfer brain.mgz -> MNI152 1-mm brain using
   FreeSurfer mri_synthmorph (default joint affine+deformable registration).
4) Apply the SAME transform to the binary lesion mask with nearest-neighbour
   interpolation.
5) Visually QC the native->MNI registration and lesion placement.
6) Quantify ipsilesional lesion load in selected JHU white-matter regions:
      ALIC  = anterior limb of internal capsule
      PLIC  = posterior limb of internal capsule
      RLIC  = retrolenticular part of internal capsule
      ACR   = anterior corona radiata
      SCR   = superior corona radiata
      PCR   = posterior corona radiata
      PTR   = posterior thalamic radiation
      EC    = external capsule
7) Merge tract lesion load with thalamic FDG-PET / ASL / CTP asymmetry.
8) Produce descriptive Spearman rho values in:
      ALL
      THALAMUS_SPARED_LT1PCT
      STRICT_DEEP_ROI_SPARED
   No p-values / confirmatory inference.

IMPORTANT
---------
Large acute lesions may affect nonlinear normalization. Therefore the QC PNG
for EVERY subject must be reviewed before interpreting tract lesion loads.

Run from the existing pse-stroke environment:
    conda activate pse-stroke
    cd <PSE_ROOT>/code
    python PSE_JHU_WM_TOPOGRAPHY_01.py
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy.stats import spearmanr

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required. Install it with:\n"
        "  conda install matplotlib -y"
    ) from exc


# ============================================================================
# Paths

ROOT = Path(os.environ.get("PSE_ROOT", str(Path.home() / "PSE"))).expanduser().resolve()

FS_HOME = Path(os.environ.get("FREESURFER_HOME", "/Applications/freesurfer/8.2.0")).expanduser().resolve()
SYNTHMORPH = FS_HOME / "bin" / "mri_synthmorph"
FS_DIR = ROOT / "derivatives" / "freesurfer"

LESION_DIR = ROOT / "derivatives" / "lesion_masks_final"
LESION_SUMMARY = (
    LESION_DIR / "PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv"
)

THALAMUS_FILE = (
    ROOT / "derivatives" / "acute_pilot_thalamus_lesion_aware"
    / "PSE_ACUTE_THALAMUS_PATIENT_LEVEL_v01.csv"
)

CTP_FILE = (
    ROOT / "derivatives" / "acute_pilot_thalamus_ctp_extended"
    / "PSE_ACUTE_THALAMUS_CTP_PATIENT_LEVEL_v01.csv"
)

FSL_DATA_ENV = Path(os.environ.get("PSE_FSL_DATA_ENV", str(Path.home() / "miniforge3/envs/pse-fsl-data"))).expanduser().resolve()

MNI_TEMPLATE = (
    FSL_DATA_ENV / "data" / "standard" / "MNI152_T1_1mm_brain.nii.gz"
)

JHU_ATLAS = (
    FSL_DATA_ENV / "data" / "atlases" / "JHU"
    / "JHU-ICBM-labels-1mm.nii.gz"
)

JHU_XML = (
    FSL_DATA_ENV / "data" / "atlases" / "JHU-labels.xml"
)

OUT_DIR = (
    ROOT / "derivatives" / "acute_pilot_jhu_wm_topography"
)

REG_DIR = OUT_DIR / "registration"
QC_DIR = OUT_DIR / "qc"

SUBJECTS = [f"sub-P{i:03d}" for i in range(1, 8)]


# ============================================================================
# Regions and outcomes

TRACT_PATTERNS = {
    "ALIC": "anterior limb of internal capsule",
    "PLIC": "posterior limb of internal capsule",
    "RLIC": "retrolenticular part of internal capsule",
    "ACR": "anterior corona radiata",
    "SCR": "superior corona radiata",
    "PCR": "posterior corona radiata",
    "PTR": "posterior thalamic radiation",
    "EC": "external capsule",
}

TRACT_ORDER = list(TRACT_PATTERNS)

OUTCOME_FIELDS = [
    "PET_FDG_AI",
    "ASL_CBF_AI",
    "CTP_CBF_AI",
    "CTP_CBV_AI",
    "CTP_MTT_AI",
    "CTP_Tmax_AI",
    "CTP_TTP_AI",
]


# ============================================================================
# Helpers

def to_float(x):
    if x is None:
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def to_bool(x):
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def fmt(x, nd=6):
    try:
        x = float(x)
    except Exception:
        return ""
    return "" if not np.isfinite(x) else f"{x:.{nd}f}"


def read_csv_by_subject(path):
    data = {}
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            data[row["Subject"].strip()] = row
    return data


def find_fs_subject_dir(subject):
    direct = FS_DIR / subject
    if direct.is_dir():
        return direct

    code = subject.replace("sub-", "")
    candidates = sorted(
        p for p in FS_DIR.glob(f"*{code}*")
        if p.is_dir()
    )

    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(
        f"Could not uniquely locate FreeSurfer directory for {subject}. "
        f"Candidates: {candidates}"
    )


def find_lesion_mask(subject):
    acute = LESION_DIR / subject / "acute"

    preferred = [
        acute / "lesion_final_in_DWI.nii",
        acute / "lesion_final_in_DWI.nii.gz",
    ]

    for p in preferred:
        if p.is_file():
            return p

    candidates = sorted(
        set(
            list(acute.glob("*final*lesion*.nii"))
            + list(acute.glob("*final*lesion*.nii.gz"))
            + list(acute.glob("*lesion*final*.nii"))
            + list(acute.glob("*lesion*final*.nii.gz"))
        )
    )

    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(
        f"Could not uniquely locate final lesion mask for {subject}. "
        f"Candidates: {candidates}"
    )


def run_command(cmd):
    print("  RUN:", " ".join(str(x) for x in cmd))

    env = os.environ.copy()
    env["FREESURFER_HOME"] = str(FS_HOME)
    env["FS_ALLOW_DEEP"] = "1"

    proc = subprocess.run(
        [str(x) for x in cmd],
        env=env,
        text=True,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed with return code "
            f"{proc.returncode}: {' '.join(str(x) for x in cmd)}"
        )


def normalize_name(name):
    s = str(name).lower()
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def detect_side(name):
    """
    Returns L/R/None from common JHU/FSL label naming conventions.
    """
    s = normalize_name(name)

    left_patterns = [
        r"\bleft\b",
        r"\bl$",
        r"\(l\)$",
    ]
    right_patterns = [
        r"\bright\b",
        r"\br$",
        r"\(r\)$",
    ]

    for pat in left_patterns:
        if re.search(pat, s):
            return "L"

    for pat in right_patterns:
        if re.search(pat, s):
            return "R"

    return None


def parse_jhu_xml():
    """
    Read label index/name pairs from FSL's JHU-labels.xml.
    The ICBM-DTI-81 atlas is a discrete Label atlas.
    """
    tree = ET.parse(JHU_XML)
    root = tree.getroot()

    labels = []

    for elem in root.findall(".//label"):
        idx = int(elem.attrib["index"])
        name = (elem.text or "").strip()

        labels.append(
            {
                "value": idx,
                "name": name,
                "normalized": normalize_name(name),
                "side": detect_side(name),
            }
        )

    return labels


def find_jhu_label(labels, phrase, side):
    phrase_n = normalize_name(phrase)

    candidates = [
        item for item in labels
        if phrase_n in item["normalized"]
        and item["side"] == side
    ]

    if len(candidates) == 1:
        return candidates[0]

    nearby = [
        item for item in labels
        if phrase_n.split()[0] in item["normalized"]
    ]

    raise RuntimeError(
        f"Could not uniquely identify JHU label for '{phrase}', side {side}.\n"
        f"Candidates matching full phrase: {candidates}\n"
        f"Related labels: {nearby}"
    )


def resample_binary_to_reference(src_path, ref_path, out_path):
    src = nib.load(str(src_path))
    ref = nib.load(str(ref_path))

    rs = resample_from_to(
        src,
        ref,
        order=0,
    )

    data = (np.asarray(rs.dataobj) > 0.5).astype(np.uint8)

    out = nib.Nifti1Image(
        data,
        ref.affine,
        ref.header,
    )
    out.set_data_dtype(np.uint8)

    nib.save(out, str(out_path))


def create_registration_qc(
    subject,
    mni_path,
    moved_t1_path,
    lesion_mni_path,
    out_path,
):
    ref_img = nib.load(str(mni_path))
    mov_img = nib.load(str(moved_t1_path))
    lesion_img = nib.load(str(lesion_mni_path))

    ref = np.asarray(ref_img.dataobj, dtype=float)
    mov = np.asarray(mov_img.dataobj, dtype=float)
    lesion = np.asarray(lesion_img.dataobj) > 0.5

    # Use lesion centre if available; otherwise MNI volume centre.
    coords = np.argwhere(lesion)

    if coords.size:
        center = np.round(coords.mean(axis=0)).astype(int)
    else:
        center = np.array(ref.shape) // 2

    center = np.clip(center, 0, np.array(ref.shape) - 1)

    slices = [
        (
            "Sagittal",
            ref[center[0], :, :],
            mov[center[0], :, :],
            lesion[center[0], :, :],
        ),
        (
            "Coronal",
            ref[:, center[1], :],
            mov[:, center[1], :],
            lesion[:, center[1], :],
        ),
        (
            "Axial",
            ref[:, :, center[2]],
            mov[:, :, center[2]],
            lesion[:, :, center[2]],
        ),
    ]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(13, 8),
    )

    for col, (title, r, m, l) in enumerate(slices):
        axes[0, col].imshow(
            np.rot90(r),
            cmap="gray",
        )
        axes[0, col].set_title(
            f"MNI template — {title}"
        )
        axes[0, col].axis("off")

        axes[1, col].imshow(
            np.rot90(r),
            cmap="gray",
        )

        mrot = np.rot90(m)

        if np.isfinite(mrot).any():
            finite = mrot[np.isfinite(mrot)]
            threshold = np.percentile(finite, 20)

            brain_like = mrot > threshold

            if np.any(brain_like):
                axes[1, col].contour(
                    brain_like.astype(float),
                    levels=[0.5],
                    linewidths=1.0,
                )

        if np.any(l):
            axes[1, col].contour(
                np.rot90(l.astype(float)),
                levels=[0.5],
                linewidths=1.5,
            )

        axes[1, col].set_title(
            f"Registered T1 boundary + lesion — {title}"
        )
        axes[1, col].axis("off")

    fig.suptitle(
        f"{subject}: native T1 → MNI152 registration QC"
    )

    fig.tight_layout()
    fig.savefig(
        out_path,
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)


def merge_outcomes(subject, thalamus_data, ctp_data):
    out = {k: np.nan for k in OUTCOME_FIELDS}

    t = thalamus_data.get(subject, {})

    for field in [
        "PET_FDG_AI",
        "ASL_CBF_AI",
        "CTP_CBF_AI",
    ]:
        out[field] = to_float(
            t.get(field)
        )

    c = ctp_data.get(subject, {})

    for field in [
        "CTP_CBF_AI",
        "CTP_CBV_AI",
        "CTP_MTT_AI",
        "CTP_Tmax_AI",
        "CTP_TTP_AI",
    ]:
        value = to_float(
            c.get(field)
        )

        if np.isfinite(value):
            out[field] = value

    return out


def tract_load(atlas, label_value, lesion, voxel_ml):
    tract = atlas == int(label_value)

    n_tract = int(
        np.count_nonzero(tract)
    )

    if n_tract == 0:
        return np.nan, np.nan

    overlap = tract & lesion

    n_overlap = int(
        np.count_nonzero(overlap)
    )

    load_pct = (
        100.0
        * n_overlap
        / n_tract
    )

    overlap_ml = (
        n_overlap
        * voxel_ml
    )

    return load_pct, overlap_ml


# ============================================================================
# Main processing

def process_subjects(label_lookup):
    lesion_summary = read_csv_by_subject(
        LESION_SUMMARY
    )
    thalamus_data = read_csv_by_subject(
        THALAMUS_FILE
    )
    ctp_data = read_csv_by_subject(
        CTP_FILE
    )

    atlas_img = nib.load(
        str(JHU_ATLAS)
    )
    atlas = np.asarray(
        atlas_img.dataobj
    ).astype(np.int32)

    mni_img = nib.load(
        str(MNI_TEMPLATE)
    )

    # The JHU 1-mm label atlas should match the MNI152 1-mm template grid.
    if (
        atlas_img.shape != mni_img.shape
        or not np.allclose(
            atlas_img.affine,
            mni_img.affine,
            atol=1e-4,
        )
    ):
        raise RuntimeError(
            "JHU atlas and MNI152 template do not have identical grids. "
            "Stop and inspect before proceeding."
        )

    voxel_ml = (
        abs(
            np.linalg.det(
                atlas_img.affine[:3, :3]
            )
        )
        / 1000.0
    )

    rows = []

    for subject in SUBJECTS:
        print()
        print("=" * 72)
        print(subject)
        print("=" * 72)

        meta = lesion_summary[subject]
        side = meta["StrokeSide"].strip().upper()

        if side not in {"L", "R"}:
            raise RuntimeError(
                f"Unexpected stroke side for {subject}: {side}"
            )

        subject_reg = (
            REG_DIR / subject
        )
        subject_reg.mkdir(
            parents=True,
            exist_ok=True,
        )

        fs_subject = find_fs_subject_dir(
            subject
        )

        brain_native = (
            fs_subject
            / "mri"
            / "brain.mgz"
        )

        lesion_native = find_lesion_mask(
            subject
        )

        lesion_fsgrid = (
            subject_reg
            / f"{subject}_lesion_in_FSbrain_grid.nii.gz"
        )

        moved_t1 = (
            subject_reg
            / f"{subject}_brain_in_MNI152_1mm.nii.gz"
        )

        transform = (
            subject_reg
            / f"{subject}_native_to_MNI152_1mm_transform.nii.gz"
        )

        lesion_mni = (
            subject_reg
            / f"{subject}_lesion_in_MNI152_1mm.nii.gz"
        )

        print(f"Stroke side: {side}")
        print(f"Native brain: {brain_native}")
        print(f"Final lesion: {lesion_native}")

        # First place lesion on exact native FreeSurfer brain grid.
        if not lesion_fsgrid.is_file():
            print(
                "Resampling binary lesion to FreeSurfer native brain grid..."
            )
            resample_binary_to_reference(
                lesion_native,
                brain_native,
                lesion_fsgrid,
            )

        # Estimate native T1 -> MNI transform.
        if (
            not moved_t1.is_file()
            or not transform.is_file()
        ):
            print(
                "Running SynthMorph native T1 -> MNI registration..."
            )

            run_command(
                [
                    SYNTHMORPH,
                    "register",
                    "-o",
                    moved_t1,
                    "-t",
                    transform,
                    brain_native,
                    MNI_TEMPLATE,
                ]
            )
        else:
            print(
                "Existing SynthMorph registration found; reusing it."
            )

        # Apply transform to lesion with NN interpolation.
        if not lesion_mni.is_file():
            print(
                "Applying transform to lesion mask with nearest-neighbour interpolation..."
            )

            run_command(
                [
                    SYNTHMORPH,
                    "apply",
                    "-m",
                    "nearest",
                    transform,
                    lesion_fsgrid,
                    lesion_mni,
                ]
            )
        else:
            print(
                "Existing MNI lesion found; reusing it."
            )

        # Check geometry.
        lesion_mni_img = nib.load(
            str(lesion_mni)
        )

        if (
            lesion_mni_img.shape != atlas_img.shape
            or not np.allclose(
                lesion_mni_img.affine,
                atlas_img.affine,
                atol=1e-4,
            )
        ):
            print(
                "MNI lesion grid differs from JHU atlas grid; "
                "resampling lesion once to exact JHU grid with nearest neighbour."
            )

            rs = resample_from_to(
                lesion_mni_img,
                atlas_img,
                order=0,
            )

            data = (
                np.asarray(rs.dataobj)
                > 0.5
            ).astype(np.uint8)

            corrected = nib.Nifti1Image(
                data,
                atlas_img.affine,
                atlas_img.header,
            )
            corrected.set_data_dtype(
                np.uint8
            )

            nib.save(
                corrected,
                str(lesion_mni),
            )

            lesion_mni_img = nib.load(
                str(lesion_mni)
            )

        lesion = (
            np.asarray(
                lesion_mni_img.dataobj
            )
            > 0.5
        )

        QC_DIR.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_path = (
            QC_DIR
            / f"{subject}_native_to_MNI_JHU_QC.png"
        )

        create_registration_qc(
            subject,
            MNI_TEMPLATE,
            moved_t1,
            lesion_mni,
            qc_path,
        )

        row = {
            "Subject": subject,
            "StrokeSide": side,
            "LesionVolume_mL_DWIgrid": to_float(
                meta.get(
                    "LesionVolume_mL_DWIgrid"
                )
            ),
            "Thalamus_IpsiOverlapPercent": to_float(
                meta.get(
                    "Thalamus_IpsiOverlapPercent"
                )
            ),
            "Hippocampus_IpsiOverlapPercent": to_float(
                meta.get(
                    "Hippocampus_IpsiOverlapPercent"
                )
            ),
            "ThalamusDirect_GE1pct": (
                to_float(
                    meta.get(
                        "Thalamus_IpsiOverlapPercent"
                    )
                )
                >= 1.0
            ),
            "DeepROI_AnyOverlap": to_bool(
                meta.get(
                    "EitherDeepROI_IpsiAnyOverlap"
                )
            ),
        }

        print("JHU ipsilesional tract lesion load:")

        for code in TRACT_ORDER:
            info = label_lookup[
                (code, side)
            ]

            load_pct, overlap_ml = tract_load(
                atlas,
                info["value"],
                lesion,
                voxel_ml,
            )

            row[
                f"{code}_LesionLoadPct"
            ] = load_pct

            row[
                f"{code}_OverlapVolume_mL"
            ] = overlap_ml

            print(
                f"  {code:4s} | "
                f"{info['name']:<55s} | "
                f"{fmt(load_pct, 3)}% "
                f"({fmt(overlap_ml, 3)} mL)"
            )

        row.update(
            merge_outcomes(
                subject,
                thalamus_data,
                ctp_data,
            )
        )

        rows.append(
            row
        )

    return rows


# ============================================================================
# Outputs

def write_patient_csv(rows):
    path = (
        OUT_DIR
        / "PSE_JHU_WM_TOPOGRAPHY_PATIENT_LEVEL_v01.csv"
    )

    fields = [
        "Subject",
        "StrokeSide",
        "LesionVolume_mL_DWIgrid",
        "Thalamus_IpsiOverlapPercent",
        "Hippocampus_IpsiOverlapPercent",
        "ThalamusDirect_GE1pct",
        "DeepROI_AnyOverlap",
    ]

    for code in TRACT_ORDER:
        fields += [
            f"{code}_LesionLoadPct",
            f"{code}_OverlapVolume_mL",
        ]

    fields += OUTCOME_FIELDS

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
        )
        writer.writeheader()

        for row in rows:
            clean = {}

            for field in fields:
                value = row.get(
                    field,
                    "",
                )

                if isinstance(
                    value,
                    (float, np.floating),
                ):
                    clean[field] = fmt(
                        value
                    )
                else:
                    clean[field] = value

            writer.writerow(
                clean
            )

    return path


def group_rows(rows, group):
    if group == "ALL":
        return list(rows)

    if group == "THALAMUS_SPARED_LT1PCT":
        return [
            r for r in rows
            if not r["ThalamusDirect_GE1pct"]
        ]

    if group == "STRICT_DEEP_ROI_SPARED":
        return [
            r for r in rows
            if not r["DeepROI_AnyOverlap"]
        ]

    raise ValueError(group)


def build_spearman_table(rows):
    groups = [
        "ALL",
        "THALAMUS_SPARED_LT1PCT",
        "STRICT_DEEP_ROI_SPARED",
    ]

    out = []

    for group in groups:
        subset = group_rows(
            rows,
            group,
        )

        for code in TRACT_ORDER:
            xfield = (
                f"{code}_LesionLoadPct"
            )

            for outcome in OUTCOME_FIELDS:
                valid = [
                    r for r in subset
                    if np.isfinite(
                        r.get(
                            xfield,
                            np.nan,
                        )
                    )
                    and np.isfinite(
                        r.get(
                            outcome,
                            np.nan,
                        )
                    )
                ]

                n_nonzero = sum(
                    r[xfield] > 0
                    for r in valid
                )

                rho = np.nan

                # Avoid manufacturing rho from almost all-zero predictors.
                if (
                    len(valid) >= 3
                    and n_nonzero >= 2
                    and len(
                        set(
                            round(
                                float(r[xfield]),
                                10,
                            )
                            for r in valid
                        )
                    ) >= 2
                ):
                    rho = float(
                        spearmanr(
                            [
                                r[xfield]
                                for r in valid
                            ],
                            [
                                r[outcome]
                                for r in valid
                            ],
                        ).statistic
                    )

                out.append(
                    {
                        "Group": group,
                        "Tract": code,
                        "Outcome": outcome,
                        "N": len(valid),
                        "N_TractNonZero": n_nonzero,
                        "SpearmanRho": rho,
                        "Subjects": ";".join(
                            r["Subject"]
                            for r in valid
                        ),
                    }
                )

    return out


def write_spearman_csv(rows):
    path = (
        OUT_DIR
        / "PSE_JHU_WM_TOPOGRAPHY_SPEARMAN_v01.csv"
    )

    fields = [
        "Group",
        "Tract",
        "Outcome",
        "N",
        "N_TractNonZero",
        "SpearmanRho",
        "Subjects",
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
        )
        writer.writeheader()

        for row in rows:
            clean = dict(row)
            clean["SpearmanRho"] = fmt(
                row["SpearmanRho"]
            )
            writer.writerow(
                clean
            )

    return path


def make_heatmap(rows):
    path = (
        OUT_DIR
        / "PSE_JHU_WM_TOPOGRAPHY_HEATMAP_v01.png"
    )

    matrix = np.asarray(
        [
            [
                r[
                    f"{code}_LesionLoadPct"
                ]
                for code in TRACT_ORDER
            ]
            for r in rows
        ],
        dtype=float,
    )

    fig, ax = plt.subplots(
        figsize=(10.5, 6)
    )

    im = ax.imshow(
        matrix,
        aspect="auto",
    )

    ax.set_xticks(
        np.arange(
            len(TRACT_ORDER)
        )
    )
    ax.set_xticklabels(
        TRACT_ORDER
    )

    ax.set_yticks(
        np.arange(
            len(rows)
        )
    )
    ax.set_yticklabels(
        [
            r["Subject"].replace(
                "sub-",
                "",
            )
            for r in rows
        ]
    )

    ax.set_xlabel(
        "JHU white-matter region"
    )
    ax.set_ylabel(
        "Subject"
    )
    ax.set_title(
        "Ipsilesional JHU white-matter lesion load"
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
    )
    cbar.set_label(
        "Region occupied by infarct (%)"
    )

    for i in range(
        matrix.shape[0]
    ):
        for j in range(
            matrix.shape[1]
        ):
            value = matrix[i, j]

            if np.isfinite(value):
                text = (
                    "0"
                    if abs(value) < 0.05
                    else f"{value:.1f}"
                )

                ax.text(
                    j,
                    i,
                    text,
                    ha="center",
                    va="center",
                    fontsize=8,
                )

    fig.tight_layout()
    fig.savefig(
        path,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    return path


# ============================================================================
# Entrypoint

def main():
    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    REG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    QC_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 80)
    print("PSE — JHU WHITE-MATTER LESION TOPOGRAPHY")
    print("=" * 80)

    required = [
        SYNTHMORPH,
        MNI_TEMPLATE,
        JHU_ATLAS,
        JHU_XML,
        LESION_SUMMARY,
        THALAMUS_FILE,
        CTP_FILE,
    ]

    for path in required:
        if not path.is_file():
            raise FileNotFoundError(
                path
            )

    print(f"SynthMorph : {SYNTHMORPH}")
    print(f"MNI 1 mm   : {MNI_TEMPLATE}")
    print(f"JHU atlas  : {JHU_ATLAS}")
    print(f"JHU XML    : {JHU_XML}")

    labels = parse_jhu_xml()

    print()
    print("Selected JHU labels:")

    label_lookup = {}

    for code, phrase in TRACT_PATTERNS.items():
        for side in ["L", "R"]:
            info = find_jhu_label(
                labels,
                phrase,
                side,
            )

            label_lookup[
                (code, side)
            ] = info

            print(
                f"  {code}-{side}: "
                f"value={info['value']} | "
                f"{info['name']}"
            )

    rows = process_subjects(
        label_lookup
    )

    patient_csv = write_patient_csv(
        rows
    )

    stats = build_spearman_table(
        rows
    )

    spearman_csv = write_spearman_csv(
        stats
    )

    heatmap = make_heatmap(
        rows
    )

    print()
    print("=" * 80)
    print("OUTPUTS")
    print("=" * 80)
    print(
        f"Patient table  : {patient_csv}"
    )
    print(
        f"Spearman table : {spearman_csv}"
    )
    print(
        f"Heatmap        : {heatmap}"
    )
    print(
        f"Registration   : {REG_DIR}"
    )
    print(
        f"QC folder      : {QC_DIR}"
    )

    print()
    print(
        "Largest absolute descriptive rho values "
        "in THALAMUS-SPARED subgroup:"
    )

    finite = [
        r for r in stats
        if r["Group"]
        == "THALAMUS_SPARED_LT1PCT"
        and np.isfinite(
            r["SpearmanRho"]
        )
    ]

    finite.sort(
        key=lambda r: abs(
            r["SpearmanRho"]
        ),
        reverse=True,
    )

    for r in finite[:20]:
        print(
            f"  {r['Tract']:4s} | "
            f"{r['Outcome']:12s} | "
            f"N={r['N']} | "
            f"nonzero={r['N_TractNonZero']} | "
            f"rho={r['SpearmanRho']:.3f}"
        )

    print()
    print("IMPORTANT:")
    print(
        "  1) Review every native->MNI QC image before interpreting tract loads."
    )
    print(
        "  2) N<=7: rho values are descriptive only; no p-values."
    )
    print(
        "  3) The thalamus-spared subgroup is the main topographic analysis."
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
