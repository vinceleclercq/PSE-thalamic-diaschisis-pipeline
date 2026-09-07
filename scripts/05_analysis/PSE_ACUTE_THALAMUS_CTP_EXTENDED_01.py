#!/usr/bin/env python3
"""
PSE_ACUTE_THALAMUS_CTP_EXTENDED_01.py

Extended lesion-aware acute thalamic CT-perfusion analysis.

Purpose
-------
Extend the previous CTP-CBF analysis to all available CTP maps in the
analysis-ready ROI dataset, especially:
    CBF, CBV, MTT, Tmax, TTP

For every map, the script extracts ipsilateral and contralateral thalamic
median values and computes:
    Delta = ipsilateral - contralateral
    AI    = (ipsilateral - contralateral) / (ipsilateral + contralateral)

Interpretation
--------------
For CBF / CBV:
    negative Delta / AI -> lower ipsilateral value

For MTT / Tmax / TTP:
    positive Delta / AI -> longer ipsilateral transit/delay

The script reports BOTH AI and absolute ipsi-contra difference.
No inferential tests are performed.

Groups
------
1) ALL
2) THALAMUS_SPARED_LT1PCT
3) STRICT_DEEP_ROI_SPARED
   = no ipsilateral thalamic OR hippocampal lesion overlap at all
"""

from __future__ import annotations

import os

import csv
from collections import defaultdict
from pathlib import Path
import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required. Install it with:\n"
        "  conda install matplotlib -y"
    ) from exc

ROOT = Path(os.environ.get("PSE_ROOT", str(Path.home() / "PSE"))).expanduser().resolve()

ROI_FILE = (
    ROOT / "derivatives" / "roi_values" / "final"
    / "PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv"
)

LESION_FILE = (
    ROOT / "derivatives" / "lesion_masks_final"
    / "PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv"
)

OUT_DIR = (
    ROOT / "derivatives" / "acute_pilot_thalamus_ctp_extended"
)

SUBJECTS = [f"sub-P{i:03d}" for i in range(1, 8)]
DIRECT_THAL_THRESHOLD_PCT = 1.0
CANONICAL_ORDER = ["CTP_CBF", "CTP_CBV", "CTP_MTT", "CTP_Tmax", "CTP_TTP"]
TIME_MAPS = {"CTP_MTT", "CTP_Tmax", "CTP_TTP"}


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


def canonical_modality(name):
    s = str(name).strip()
    low = s.lower().replace("-", "_")
    mapping = {
        "ctp_cbf": "CTP_CBF",
        "ctp_cbv": "CTP_CBV",
        "ctp_mtt": "CTP_MTT",
        "ctp_tmax": "CTP_Tmax",
        "ctp_ttp": "CTP_TTP",
    }
    return mapping.get(low, s)


def read_lesion_summary():
    out = {}
    with LESION_FILE.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sub = row["Subject"].strip()
            thal_pct = to_float(row["Thalamus_IpsiOverlapPercent"])
            out[sub] = {
                "StrokeSide": row["StrokeSide"].strip(),
                "LesionVolume_mL": to_float(row["LesionVolume_mL_DWIgrid"]),
                "Thalamus_IpsiOverlapPercent": thal_pct,
                "Hippocampus_IpsiOverlapPercent": to_float(
                    row["Hippocampus_IpsiOverlapPercent"]
                ),
                "ThalamusDirect_GE1pct": (
                    np.isfinite(thal_pct)
                    and thal_pct >= DIRECT_THAL_THRESHOLD_PCT
                ),
                "DeepROI_AnyOverlap": to_bool(
                    row["EitherDeepROI_IpsiAnyOverlap"]
                ),
            }
    return out


def read_ctp_thalamus():
    data = defaultdict(dict)
    available = set()

    with ROI_FILE.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        for row in reader:
            if row["Structure"].strip() != "Thalamus":
                continue

            mod = canonical_modality(row["Modality"])

            if not mod.lower().startswith("ctp_"):
                continue
            if mod not in CANONICAL_ORDER:
                continue

            usable = str(row.get("PrimaryPairUsable", "")).strip()
            flag = str(row.get("AnalysisFlag", "")).strip()

            if usable not in {"1", "True", "true"}:
                continue
            if flag != "OK":
                continue

            sub = row["Subject"].strip()
            ipsi = to_float(row["IpsiMedian"])
            contra = to_float(row["ContraMedian"])
            ai = to_float(row["Median_AI_IpsiContra"])
            ratio = to_float(row["Median_Ratio_IpsiContra"])

            delta = np.nan
            if np.isfinite(ipsi) and np.isfinite(contra):
                delta = ipsi - contra

            data[sub][mod] = {
                "IpsiMedian": ipsi,
                "ContraMedian": contra,
                "Delta_IpsiMinusContra": delta,
                "AI": ai,
                "Ratio": ratio,
            }
            available.add(mod)

    ordered = [m for m in CANONICAL_ORDER if m in available]
    return data, ordered


def build_patient_rows(lesion, ctp, modalities):
    rows = []

    for sub in SUBJECTS:
        l = lesion.get(sub, {})
        row = {
            "Subject": sub,
            "StrokeSide": l.get("StrokeSide", ""),
            "LesionVolume_mL": l.get("LesionVolume_mL", np.nan),
            "Thalamus_IpsiOverlapPercent": l.get(
                "Thalamus_IpsiOverlapPercent", np.nan
            ),
            "Hippocampus_IpsiOverlapPercent": l.get(
                "Hippocampus_IpsiOverlapPercent", np.nan
            ),
            "ThalamusDirect_GE1pct": l.get(
                "ThalamusDirect_GE1pct", False
            ),
            "DeepROI_AnyOverlap": l.get(
                "DeepROI_AnyOverlap", False
            ),
        }

        for mod in modalities:
            d = ctp.get(sub, {}).get(mod, {})
            row[f"{mod}_IpsiMedian"] = d.get("IpsiMedian", np.nan)
            row[f"{mod}_ContraMedian"] = d.get("ContraMedian", np.nan)
            row[f"{mod}_Delta"] = d.get("Delta_IpsiMinusContra", np.nan)
            row[f"{mod}_AI"] = d.get("AI", np.nan)
            row[f"{mod}_Ratio"] = d.get("Ratio", np.nan)

        rows.append(row)

    return rows


def members(rows, group):
    if group == "ALL":
        return list(rows)
    if group == "THALAMUS_SPARED_LT1PCT":
        return [r for r in rows if not r["ThalamusDirect_GE1pct"]]
    if group == "STRICT_DEEP_ROI_SPARED":
        return [r for r in rows if not r["DeepROI_AnyOverlap"]]
    raise ValueError(group)


def safe_median(vals):
    a = np.asarray([v for v in vals if np.isfinite(v)], dtype=float)
    return np.nan if a.size == 0 else float(np.median(a))


def build_summary(rows, modalities):
    groups = ["ALL", "THALAMUS_SPARED_LT1PCT", "STRICT_DEEP_ROI_SPARED"]
    out = []

    for group in groups:
        subset = members(rows, group)

        for mod in modalities:
            ais = [r[f"{mod}_AI"] for r in subset if np.isfinite(r[f"{mod}_AI"])]
            deltas = [
                r[f"{mod}_Delta"] for r in subset
                if np.isfinite(r[f"{mod}_Delta"])
            ]
            ratios = [
                r[f"{mod}_Ratio"] for r in subset
                if np.isfinite(r[f"{mod}_Ratio"])
            ]

            out.append(
                {
                    "Group": group,
                    "Modality": mod,
                    "N": len(ais),
                    "Subjects": ";".join(
                        r["Subject"]
                        for r in subset
                        if np.isfinite(r[f"{mod}_AI"])
                    ),
                    "Median_AI": safe_median(ais),
                    "N_AI_negative": sum(v < 0 for v in ais),
                    "N_AI_positive": sum(v > 0 for v in ais),
                    "Median_Delta_IpsiMinusContra": safe_median(deltas),
                    "N_Delta_negative": sum(v < 0 for v in deltas),
                    "N_Delta_positive": sum(v > 0 for v in deltas),
                    "Median_Ratio": safe_median(ratios),
                }
            )

    return out


def write_patient_csv(rows, modalities):
    out = OUT_DIR / "PSE_ACUTE_THALAMUS_CTP_PATIENT_LEVEL_v01.csv"

    fields = [
        "Subject",
        "StrokeSide",
        "LesionVolume_mL",
        "Thalamus_IpsiOverlapPercent",
        "Hippocampus_IpsiOverlapPercent",
        "ThalamusDirect_GE1pct",
        "DeepROI_AnyOverlap",
    ]

    for mod in modalities:
        fields += [
            f"{mod}_IpsiMedian",
            f"{mod}_ContraMedian",
            f"{mod}_Delta",
            f"{mod}_AI",
            f"{mod}_Ratio",
        ]

    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            clean = {}
            for k in fields:
                v = row.get(k, "")
                clean[k] = fmt(v) if isinstance(v, float) else v
            writer.writerow(clean)

    return out


def write_summary_csv(rows):
    out = OUT_DIR / "PSE_ACUTE_THALAMUS_CTP_SUBGROUP_SUMMARY_v01.csv"

    fields = [
        "Group",
        "Modality",
        "N",
        "Subjects",
        "Median_AI",
        "N_AI_negative",
        "N_AI_positive",
        "Median_Delta_IpsiMinusContra",
        "N_Delta_negative",
        "N_Delta_positive",
        "Median_Ratio",
    ]

    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            clean = dict(row)
            for k in ["Median_AI", "Median_Delta_IpsiMinusContra", "Median_Ratio"]:
                clean[k] = fmt(row[k])
            writer.writerow(clean)

    return out


def subject_label(row):
    label = row["Subject"].replace("sub-", "")
    if row["ThalamusDirect_GE1pct"]:
        label += "†"
    return label


def plot_patient_ai(rows, modalities):
    out = OUT_DIR / "PSE_ACUTE_THALAMUS_CTP_PATIENT_LEVEL_AI_v01.png"

    fig, ax = plt.subplots(figsize=(12, 6.5))
    x = np.arange(len(rows), dtype=float)
    offsets = np.linspace(-0.28, 0.28, max(1, len(modalities)))
    markers = ["o", "s", "^", "D", "v"]

    for j, mod in enumerate(modalities):
        xs, ys = [], []
        for i, row in enumerate(rows):
            y = row[f"{mod}_AI"]
            if np.isfinite(y):
                xs.append(x[i] + offsets[j])
                ys.append(y)

        ax.scatter(
            xs, ys,
            s=75,
            marker=markers[j % len(markers)],
            label=mod.replace("CTP_", ""),
        )

    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels([subject_label(r) for r in rows])
    ax.set_xlabel("Subject")
    ax.set_ylabel("Thalamic asymmetry index")
    ax.set_title(
        "Acute thalamic asymmetry across CT-perfusion maps\n"
        "† = direct ipsilateral thalamic lesion overlap ≥1%"
    )
    ax.legend(ncol=min(5, len(modalities)))
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_time_delta(rows, modalities):
    time_mods = [m for m in modalities if m in TIME_MAPS]
    if not time_mods:
        return None

    out = OUT_DIR / "PSE_ACUTE_THALAMUS_CTP_TIME_MAP_DELTA_v01.png"

    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(rows), dtype=float)
    offsets = np.linspace(-0.22, 0.22, len(time_mods))
    markers = ["o", "s", "^"]

    for j, mod in enumerate(time_mods):
        xs, ys = [], []
        for i, row in enumerate(rows):
            y = row[f"{mod}_Delta"]
            if np.isfinite(y):
                xs.append(x[i] + offsets[j])
                ys.append(y)

        ax.scatter(
            xs, ys,
            s=80,
            marker=markers[j % len(markers)],
            label=mod.replace("CTP_", ""),
        )

    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels([subject_label(r) for r in rows])
    ax.set_xlabel("Subject")
    ax.set_ylabel("Ipsilateral − contralateral thalamic value")
    ax.set_title(
        "Acute thalamic side-to-side differences in CTP time maps\n"
        "Positive values = longer ipsilateral transit/delay"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_group_medians(summary, modalities):
    out = OUT_DIR / "PSE_ACUTE_THALAMUS_CTP_SUBGROUP_MEDIANS_AI_v01.png"

    groups = ["ALL", "THALAMUS_SPARED_LT1PCT", "STRICT_DEEP_ROI_SPARED"]
    labels = ["All", "Thalamus spared <1%", "No thal./hipp. overlap"]
    lookup = {(r["Group"], r["Modality"]): r for r in summary}

    fig, ax = plt.subplots(figsize=(11.5, 6.5))
    x = np.arange(len(groups), dtype=float)
    offsets = np.linspace(-0.30, 0.30, max(1, len(modalities)))
    markers = ["o", "s", "^", "D", "v"]

    for j, mod in enumerate(modalities):
        xs, ys = [], []
        for i, group in enumerate(groups):
            row = lookup[(group, mod)]
            y = row["Median_AI"]
            if np.isfinite(y):
                xs.append(x[i] + offsets[j])
                ys.append(y)

        ax.scatter(
            xs, ys,
            s=90,
            marker=markers[j % len(markers)],
            label=mod.replace("CTP_", ""),
        )

    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Median thalamic asymmetry index")
    ax.set_title("Lesion-aware sensitivity analysis across CTP maps")
    ax.legend(ncol=min(5, len(modalities)))
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("PSE ACUTE THALAMUS — EXTENDED CT PERFUSION ANALYSIS")
    print("=" * 78)
    print(f"ROI file    : {ROI_FILE}")
    print(f"Lesion file : {LESION_FILE}")
    print(f"Output dir  : {OUT_DIR}")
    print("=" * 78)

    if not ROI_FILE.is_file():
        raise FileNotFoundError(ROI_FILE)
    if not LESION_FILE.is_file():
        raise FileNotFoundError(LESION_FILE)

    lesion = read_lesion_summary()
    ctp, modalities = read_ctp_thalamus()

    if not modalities:
        raise RuntimeError(
            "No usable CTP thalamic modalities were found in the analysis-ready CSV."
        )

    print("Usable CTP modalities found:")
    for mod in modalities:
        print(f"  - {mod}")

    patient_rows = build_patient_rows(lesion, ctp, modalities)
    summary_rows = build_summary(patient_rows, modalities)

    patient_csv = write_patient_csv(patient_rows, modalities)
    summary_csv = write_summary_csv(summary_rows)
    fig_ai = plot_patient_ai(patient_rows, modalities)
    fig_time = plot_time_delta(patient_rows, modalities)
    fig_groups = plot_group_medians(summary_rows, modalities)

    print()
    print("Outputs:")
    print(f"  {patient_csv}")
    print(f"  {summary_csv}")
    print(f"  {fig_ai}")
    if fig_time is not None:
        print(f"  {fig_time}")
    print(f"  {fig_groups}")

    print()
    print("Subgroup summary:")
    for row in summary_rows:
        print(
            f"  {row['Group']} | {row['Modality']} | "
            f"N={row['N']} | "
            f"median AI={fmt(row['Median_AI'], 4)} | "
            f"AI<0={row['N_AI_negative']}/{row['N']} | "
            f"AI>0={row['N_AI_positive']}/{row['N']} | "
            f"median delta={fmt(row['Median_Delta_IpsiMinusContra'], 4)}"
        )

    print()
    print("Interpretation reminder:")
    print("  CBF/CBV: negative AI or delta = lower ipsilateral value.")
    print("  MTT/Tmax/TTP: positive AI or delta = longer ipsilateral delay.")
    print("  All results are descriptive; no inferential tests are performed.")
    print("=" * 78)


if __name__ == "__main__":
    main()
