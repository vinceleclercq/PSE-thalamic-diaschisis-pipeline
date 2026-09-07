#!/usr/bin/env python3
from __future__ import annotations

import os

import csv
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from scipy.stats import spearmanr

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit("matplotlib is required. Install it with: conda install matplotlib -y") from exc

ROOT = Path(os.environ.get("PSE_ROOT", str(Path.home() / "PSE"))).expanduser().resolve()
FS_DIR = ROOT / "derivatives" / "freesurfer"
LESION_DIR = ROOT / "derivatives" / "lesion_masks_final"
LESION_SUMMARY = LESION_DIR / "PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv"
THALAMUS_FILE = ROOT / "derivatives" / "acute_pilot_thalamus_lesion_aware" / "PSE_ACUTE_THALAMUS_PATIENT_LEVEL_v01.csv"
CTP_FILE = ROOT / "derivatives" / "acute_pilot_thalamus_ctp_extended" / "PSE_ACUTE_THALAMUS_CTP_PATIENT_LEVEL_v01.csv"
OUT_DIR = ROOT / "derivatives" / "acute_pilot_lesion_topography"
QC_DIR = OUT_DIR / "qc"

SUBJECTS = [f"sub-P{i:03d}" for i in range(1, 8)]

SUBCORTICAL_LABELS = {
    "L": {"Thalamus":[10], "Caudate":[11], "Putamen":[12], "Pallidum":[13], "Hippocampus":[17], "Accumbens":[26], "VentralDC":[28]},
    "R": {"Thalamus":[49], "Caudate":[50], "Putamen":[51], "Pallidum":[52], "Hippocampus":[53], "Accumbens":[58], "VentralDC":[60]},
}

CORTICAL_GROUPS_BASE = {
    "FrontalCortex":[3,12,14,17,18,19,20,24,27,28,32],
    "ParietalCortex":[8,22,25,29,31],
    "TemporalCortex":[1,6,7,9,15,16,30,33,34],
    "OccipitalCortex":[5,11,13,21],
    "CingulateCortex":[2,10,23,26],
    "Insula":[35],
}

REGION_ORDER = [
    "Caudate","Putamen","Pallidum","Lentiform","Thalamus","Hippocampus",
    "Accumbens","VentralDC","Insula","FrontalCortex","ParietalCortex",
    "TemporalCortex","OccipitalCortex","CingulateCortex"
]

OUTCOME_FIELDS = [
    "PET_FDG_AI","ASL_CBF_AI","CTP_CBF_AI","CTP_CBV_AI",
    "CTP_MTT_AI","CTP_Tmax_AI","CTP_TTP_AI"
]

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
    candidates = sorted(p for p in FS_DIR.glob(f"*{subject.replace('sub-','')}*") if p.is_dir())
    if len(candidates) == 1:
        return candidates[0]
    raise FileNotFoundError(f"Could not uniquely locate FreeSurfer directory for {subject}. Candidates: {candidates}")

def find_lesion_mask(subject):
    acute = LESION_DIR / subject / "acute"
    for p in [acute/"lesion_final_in_DWI.nii", acute/"lesion_final_in_DWI.nii.gz"]:
        if p.is_file():
            return p
    globs = sorted(set(
        list(acute.glob("*final*lesion*.nii")) +
        list(acute.glob("*final*lesion*.nii.gz")) +
        list(acute.glob("*lesion*final*.nii")) +
        list(acute.glob("*lesion*final*.nii.gz"))
    ))
    if len(globs) == 1:
        return globs[0]
    raise FileNotFoundError(f"Could not uniquely locate final lesion mask for {subject}. Candidates: {globs}")

def cortical_codes(side, base_ids):
    prefix = 1000 if side == "L" else 2000
    return [prefix+x for x in base_ids]

def region_mask(aparc, side, region):
    if region == "Lentiform":
        return np.isin(aparc, SUBCORTICAL_LABELS[side]["Putamen"] + SUBCORTICAL_LABELS[side]["Pallidum"])
    if region in SUBCORTICAL_LABELS[side]:
        return np.isin(aparc, SUBCORTICAL_LABELS[side][region])
    if region in CORTICAL_GROUPS_BASE:
        return np.isin(aparc, cortical_codes(side, CORTICAL_GROUPS_BASE[region]))
    raise KeyError(region)

def lesion_load(mask, lesion, voxel_ml):
    n_region = int(mask.sum())
    if n_region == 0:
        return np.nan, np.nan
    n_overlap = int((mask & lesion).sum())
    return 100.0*n_overlap/n_region, n_overlap*voxel_ml

def merge_outcomes(subject, thal, ctp):
    out = {k:np.nan for k in OUTCOME_FIELDS}
    t = thal.get(subject,{})
    out["PET_FDG_AI"] = to_float(t.get("PET_FDG_AI"))
    out["ASL_CBF_AI"] = to_float(t.get("ASL_CBF_AI"))
    out["CTP_CBF_AI"] = to_float(t.get("CTP_CBF_AI"))
    c = ctp.get(subject,{})
    for f in ["CTP_CBF_AI","CTP_CBV_AI","CTP_MTT_AI","CTP_Tmax_AI","CTP_TTP_AI"]:
        v = to_float(c.get(f))
        if np.isfinite(v):
            out[f] = v
    return out

def save_qc(subject, brain_img, lesion, out_path):
    brain = np.asarray(brain_img.dataobj, dtype=float)
    coords = np.argwhere(lesion)
    center = np.array(brain.shape)//2 if coords.size == 0 else np.round(coords.mean(axis=0)).astype(int)
    center = np.clip(center, 0, np.array(brain.shape)-1)
    slices = [
        ("Sagittal", brain[center[0],:,:], lesion[center[0],:,:]),
        ("Coronal", brain[:,center[1],:], lesion[:,center[1],:]),
        ("Axial", brain[:,:,center[2]], lesion[:,:,center[2]]),
    ]
    fig, axes = plt.subplots(1,3,figsize=(13,4.5))
    for ax,(title,bg,msk) in zip(axes,slices):
        ax.imshow(np.rot90(bg), cmap="gray")
        if np.any(msk):
            ax.contour(np.rot90(msk.astype(float)), levels=[0.5], linewidths=1.5)
        ax.set_title(title)
        ax.axis("off")
    fig.suptitle(f"{subject}: final lesion in FreeSurfer aparc+aseg space")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

def build_patient_table():
    lesion_summary = read_csv_by_subject(LESION_SUMMARY)
    thalamus_data = read_csv_by_subject(THALAMUS_FILE)
    ctp_data = read_csv_by_subject(CTP_FILE)
    rows = []
    for subject in SUBJECTS:
        print(f"\n[{subject}]")
        sr = lesion_summary[subject]
        side = sr["StrokeSide"].strip().upper()
        fs_subject = find_fs_subject_dir(subject)
        aparc_path = fs_subject/"mri"/"aparc+aseg.mgz"
        brain_path = fs_subject/"mri"/"brain.mgz"
        lesion_path = find_lesion_mask(subject)
        if not aparc_path.is_file(): raise FileNotFoundError(aparc_path)
        if not brain_path.is_file(): raise FileNotFoundError(brain_path)

        lesion_img = nib.load(str(lesion_path))
        aparc_img = nib.load(str(aparc_path))
        brain_img = nib.load(str(brain_path))
        lesion_rs_img = resample_from_to(lesion_img, aparc_img, order=0)
        lesion_rs = np.asarray(lesion_rs_img.dataobj) > 0.5
        aparc = np.asarray(aparc_img.dataobj).astype(np.int32)
        voxel_ml = abs(np.linalg.det(aparc_img.affine[:3,:3]))/1000.0

        row = {
            "Subject":subject,
            "StrokeSide":side,
            "LesionVolume_mL_DWIgrid":to_float(sr.get("LesionVolume_mL_DWIgrid")),
            "Thalamus_IpsiOverlapPercent_previous":to_float(sr.get("Thalamus_IpsiOverlapPercent")),
            "Hippocampus_IpsiOverlapPercent_previous":to_float(sr.get("Hippocampus_IpsiOverlapPercent")),
        }

        for region in REGION_ORDER:
            load, vol = lesion_load(region_mask(aparc, side, region), lesion_rs, voxel_ml)
            row[f"{region}_LesionLoadPct"] = load
            row[f"{region}_OverlapVolume_mL"] = vol
            print(f"  {region:16s}: {fmt(load,3)}% ({fmt(vol,3)} mL)")

        row.update(merge_outcomes(subject, thalamus_data, ctp_data))
        QC_DIR.mkdir(parents=True, exist_ok=True)
        save_qc(subject, brain_img, lesion_rs, QC_DIR/f"{subject}_lesion_topography_QC.png")
        rows.append(row)
    return rows

def write_patient_csv(rows):
    path = OUT_DIR/"PSE_LESION_TOPOGRAPHY_PATIENT_LEVEL_v01.csv"
    fields = ["Subject","StrokeSide","LesionVolume_mL_DWIgrid",
              "Thalamus_IpsiOverlapPercent_previous","Hippocampus_IpsiOverlapPercent_previous"]
    for r in REGION_ORDER:
        fields += [f"{r}_LesionLoadPct", f"{r}_OverlapVolume_mL"]
    fields += OUTCOME_FIELDS
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            clean = {k:(fmt(row.get(k)) if isinstance(row.get(k),(float,np.floating)) else row.get(k,"")) for k in fields}
            w.writerow(clean)
    return path

def build_spearman(rows):
    out = []
    for region in REGION_ORDER:
        rf = f"{region}_LesionLoadPct"
        for outcome in OUTCOME_FIELDS:
            valid = [r for r in rows if np.isfinite(r.get(rf,np.nan)) and np.isfinite(r.get(outcome,np.nan))]
            nonzero = sum(r[rf] > 0 for r in valid)
            rho = np.nan
            if len(valid) >= 3 and nonzero >= 2:
                rho = float(spearmanr([r[rf] for r in valid], [r[outcome] for r in valid]).statistic)
            out.append({"Region":region,"Outcome":outcome,"N":len(valid),
                        "N_RegionNonZero":nonzero,"SpearmanRho":rho,
                        "Subjects":";".join(r["Subject"] for r in valid)})
    return out

def write_spearman(rows):
    path = OUT_DIR/"PSE_LESION_TOPOGRAPHY_SPEARMAN_v01.csv"
    fields = ["Region","Outcome","N","N_RegionNonZero","SpearmanRho","Subjects"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            x = dict(row); x["SpearmanRho"] = fmt(row["SpearmanRho"]); w.writerow(x)
    return path

def make_heatmap(rows):
    path = OUT_DIR/"PSE_LESION_TOPOGRAPHY_HEATMAP_v01.png"
    matrix = np.asarray([[r[f"{reg}_LesionLoadPct"] for reg in REGION_ORDER] for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(14,6))
    im = ax.imshow(matrix, aspect="auto")
    ax.set_xticks(np.arange(len(REGION_ORDER)))
    ax.set_xticklabels(REGION_ORDER, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([r["Subject"].replace("sub-","") for r in rows])
    ax.set_title("Ipsilesional regional lesion load")
    ax.set_xlabel("FreeSurfer anatomical region")
    ax.set_ylabel("Subject")
    cbar = fig.colorbar(im, ax=ax); cbar.set_label("Region occupied by infarct (%)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix[i,j]
            if np.isfinite(v):
                ax.text(j,i,"0" if abs(v)<0.05 else f"{v:.1f}",ha="center",va="center",fontsize=7)
    fig.tight_layout()
    fig.savefig(path,dpi=200,bbox_inches="tight")
    plt.close(fig)
    return path

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    print("="*80)
    print("PSE ACUTE PILOT — EXPLORATORY LESION TOPOGRAPHY")
    print("="*80)
    print("FreeSurfer aparc+aseg-based regional lesion load.")
    print("Internal capsule is NOT included in this version.")
    print("Nearest-neighbour resampling is used for the binary lesion mask.")
    print("="*80)

    for req in [LESION_SUMMARY, THALAMUS_FILE, CTP_FILE]:
        if not req.is_file():
            raise FileNotFoundError(req)

    rows = build_patient_table()
    patient_csv = write_patient_csv(rows)
    spearman_rows = build_spearman(rows)
    spearman_csv = write_spearman(spearman_rows)
    heatmap = make_heatmap(rows)

    print("\nOUTPUTS")
    print(f"  Patient-level table : {patient_csv}")
    print(f"  Spearman table      : {spearman_csv}")
    print(f"  Heatmap             : {heatmap}")
    print(f"  QC folder           : {QC_DIR}")

    finite = [r for r in spearman_rows if np.isfinite(r["SpearmanRho"])]
    finite.sort(key=lambda r: abs(r["SpearmanRho"]), reverse=True)
    print("\nLargest absolute descriptive rho values:")
    for r in finite[:20]:
        print(f"  {r['Region']:16s} | {r['Outcome']:12s} | N={r['N']} | "
              f"nonzero={r['N_RegionNonZero']} | rho={r['SpearmanRho']:.3f}")
    print("\nIMPORTANT: descriptive only, N<=7, no p-values.")
    print("First inspect the heatmap and QC images before interpreting associations.")
    print("="*80)

if __name__ == "__main__":
    main()
