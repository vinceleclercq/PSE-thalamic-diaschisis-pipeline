#!/usr/bin/env python3
"""
PSE_ACUTE_THALAMUS_LESION_AWARE_ANALYSIS_01.py

Lesion-aware revision of the acute thalamic pilot analysis.

Purpose
-------
Combine:
1) acute thalamic PET/ASL/CTP asymmetry measures, and
2) direct ischemic overlap with the ipsilateral thalamus/hippocampus.

Primary thalamic modalities:
    PET_FDG
    ASL_CBF
    CTP_CBF

Primary quantitative measure:
    Median_AI_IpsiContra

Direct thalamic involvement:
    ipsilateral thalamic ROI overlap >= 1%

Strict deep-ROI-spared subgroup:
    no ipsilateral thalamic OR hippocampal overlap at all

Outputs
-------
<PSE_ROOT>/derivatives/acute_pilot_thalamus_lesion_aware/

    PSE_ACUTE_THALAMUS_PATIENT_LEVEL_v01.csv
    PSE_ACUTE_THALAMUS_SUBGROUP_SUMMARY_v01.csv
    PSE_ACUTE_THALAMUS_PATIENT_LEVEL_AI_v01.png
    PSE_ACUTE_THALAMUS_SUBGROUP_MEDIANS_v01.png

The patient-level plot marks subjects with direct thalamic involvement >=1%
by adding a dagger (†) to the subject label.

No inferential hypothesis tests are performed.
"""

from __future__ import annotations

import os

import csv
import math
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
    ROOT
    / "derivatives"
    / "roi_values"
    / "final"
    / "PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv"
)

LESION_FILE = (
    ROOT
    / "derivatives"
    / "lesion_masks_final"
    / "PSE_LESION_DIRECT_INVOLVEMENT_SUMMARY_v01.csv"
)

OUT_DIR = (
    ROOT
    / "derivatives"
    / "acute_pilot_thalamus_lesion_aware"
)

PRIMARY_MODALITIES = ["PET_FDG", "ASL_CBF", "CTP_CBF"]
SUBJECTS = [f"sub-P{i:03d}" for i in range(1, 8)]

DIRECT_THALAMUS_THRESHOLD_PCT = 1.0


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
    s = str(x).strip().lower()
    return s in {"true", "1", "yes", "y"}


def safe_median(values):
    arr = np.asarray(
        [v for v in values if np.isfinite(v)],
        dtype=float,
    )

    if arr.size == 0:
        return np.nan

    return float(np.median(arr))


def fmt(x, nd=6):
    return "" if not np.isfinite(x) else f"{x:.{nd}f}"


def read_lesion_summary():
    lesion = {}

    with LESION_FILE.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:
        reader = csv.DictReader(f)

        for row in reader:
            sub = row["Subject"].strip()

            thal_pct = to_float(
                row["Thalamus_IpsiOverlapPercent"]
            )
            hip_pct = to_float(
                row["Hippocampus_IpsiOverlapPercent"]
            )

            lesion[sub] = {
                "StrokeSide": row["StrokeSide"].strip(),
                "QDecision": row["QDecision"].strip(),
                "LesionVolume_mL": to_float(
                    row["LesionVolume_mL_DWIgrid"]
                ),
                "Thalamus_IpsiOverlapPercent": thal_pct,
                "Hippocampus_IpsiOverlapPercent": hip_pct,
                "ThalamusDirect_GE1pct": (
                    np.isfinite(thal_pct)
                    and thal_pct
                    >= DIRECT_THALAMUS_THRESHOLD_PCT
                ),
                "DeepROI_AnyOverlap": to_bool(
                    row["EitherDeepROI_IpsiAnyOverlap"]
                ),
            }

    return lesion


def read_primary_thalamic_data():
    data = defaultdict(dict)

    with ROI_FILE.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as f:
        reader = csv.DictReader(f)

        for row in reader:
            sub = row["Subject"].strip()
            modality = row["Modality"].strip()
            structure = row["Structure"].strip()

            if structure != "Thalamus":
                continue

            if modality not in PRIMARY_MODALITIES:
                continue

            usable = str(
                row.get("PrimaryPairUsable", "")
            ).strip()

            analysis_flag = str(
                row.get("AnalysisFlag", "")
            ).strip()

            if usable not in {"1", "True", "true"}:
                continue

            if analysis_flag != "OK":
                continue

            data[sub][modality] = {
                "Median_AI_IpsiContra": to_float(
                    row["Median_AI_IpsiContra"]
                ),
                "Median_Ratio_IpsiContra": to_float(
                    row["Median_Ratio_IpsiContra"]
                ),
                "IpsiMedian": to_float(
                    row["IpsiMedian"]
                ),
                "ContraMedian": to_float(
                    row["ContraMedian"]
                ),
                "IpsiCoverage": to_float(
                    row["Ipsi_StatsCoverageFraction"]
                ),
                "ContraCoverage": to_float(
                    row["Contra_StatsCoverageFraction"]
                ),
            }

    return data


def build_patient_table(lesion, roi):
    rows = []

    for sub in SUBJECTS:
        l = lesion.get(sub, {})

        row = {
            "Subject": sub,
            "StrokeSide": l.get("StrokeSide", ""),
            "QDecision": l.get("QDecision", ""),
            "LesionVolume_mL": l.get(
                "LesionVolume_mL",
                np.nan,
            ),
            "Thalamus_IpsiOverlapPercent": l.get(
                "Thalamus_IpsiOverlapPercent",
                np.nan,
            ),
            "Hippocampus_IpsiOverlapPercent": l.get(
                "Hippocampus_IpsiOverlapPercent",
                np.nan,
            ),
            "ThalamusDirect_GE1pct": l.get(
                "ThalamusDirect_GE1pct",
                False,
            ),
            "DeepROI_AnyOverlap": l.get(
                "DeepROI_AnyOverlap",
                False,
            ),
        }

        for mod in PRIMARY_MODALITIES:
            d = roi.get(sub, {}).get(mod, {})

            row[f"{mod}_AI"] = d.get(
                "Median_AI_IpsiContra",
                np.nan,
            )
            row[f"{mod}_Ratio"] = d.get(
                "Median_Ratio_IpsiContra",
                np.nan,
            )
            row[f"{mod}_IpsiMedian"] = d.get(
                "IpsiMedian",
                np.nan,
            )
            row[f"{mod}_ContraMedian"] = d.get(
                "ContraMedian",
                np.nan,
            )

        rows.append(row)

    return rows


def write_patient_table(rows):
    out = OUT_DIR / "PSE_ACUTE_THALAMUS_PATIENT_LEVEL_v01.csv"

    fields = [
        "Subject",
        "StrokeSide",
        "QDecision",
        "LesionVolume_mL",
        "Thalamus_IpsiOverlapPercent",
        "Hippocampus_IpsiOverlapPercent",
        "ThalamusDirect_GE1pct",
        "DeepROI_AnyOverlap",
    ]

    for mod in PRIMARY_MODALITIES:
        fields.extend(
            [
                f"{mod}_AI",
                f"{mod}_Ratio",
                f"{mod}_IpsiMedian",
                f"{mod}_ContraMedian",
            ]
        )

    with out.open(
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

            for k in fields:
                v = row.get(k, "")

                if isinstance(v, float):
                    clean[k] = fmt(v)
                else:
                    clean[k] = v

            writer.writerow(clean)

    return out


def group_members(rows, group_name):
    if group_name == "ALL":
        return list(rows)

    if group_name == "THALAMUS_SPARED_LT1PCT":
        return [
            r
            for r in rows
            if not r["ThalamusDirect_GE1pct"]
        ]

    if group_name == "STRICT_DEEP_ROI_SPARED":
        return [
            r
            for r in rows
            if not r["DeepROI_AnyOverlap"]
        ]

    raise ValueError(group_name)


def build_group_summary(rows):
    groups = [
        "ALL",
        "THALAMUS_SPARED_LT1PCT",
        "STRICT_DEEP_ROI_SPARED",
    ]

    out_rows = []

    for group in groups:
        members = group_members(rows, group)

        for mod in PRIMARY_MODALITIES:
            ais = [
                r[f"{mod}_AI"]
                for r in members
                if np.isfinite(r[f"{mod}_AI"])
            ]

            ratios = [
                r[f"{mod}_Ratio"]
                for r in members
                if np.isfinite(r[f"{mod}_Ratio"])
            ]

            neg = sum(v < 0 for v in ais)

            out_rows.append(
                {
                    "Group": group,
                    "Modality": mod,
                    "SubjectsInGroup": ";".join(
                        r["Subject"]
                        for r in members
                    ),
                    "N": len(ais),
                    "Median_AI_IpsiContra": safe_median(
                        ais
                    ),
                    "N_AI_negative": neg,
                    "Fraction_AI_negative": (
                        np.nan
                        if len(ais) == 0
                        else neg / len(ais)
                    ),
                    "Median_Ratio_IpsiContra": safe_median(
                        ratios
                    ),
                }
            )

    return out_rows


def write_group_summary(rows):
    out = (
        OUT_DIR
        / "PSE_ACUTE_THALAMUS_SUBGROUP_SUMMARY_v01.csv"
    )

    fields = [
        "Group",
        "Modality",
        "SubjectsInGroup",
        "N",
        "Median_AI_IpsiContra",
        "N_AI_negative",
        "Fraction_AI_negative",
        "Median_Ratio_IpsiContra",
    ]

    with out.open(
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

            for k in [
                "Median_AI_IpsiContra",
                "Fraction_AI_negative",
                "Median_Ratio_IpsiContra",
            ]:
                clean[k] = fmt(
                    to_float(clean[k])
                )

            writer.writerow(clean)

    return out


def plot_patient_level(rows):
    out = (
        OUT_DIR
        / "PSE_ACUTE_THALAMUS_PATIENT_LEVEL_AI_v01.png"
    )

    x = np.arange(len(rows), dtype=float)

    offsets = {
        "PET_FDG": -0.18,
        "ASL_CBF": 0.0,
        "CTP_CBF": 0.18,
    }

    markers = {
        "PET_FDG": "o",
        "ASL_CBF": "s",
        "CTP_CBF": "^",
    }

    fig, ax = plt.subplots(
        figsize=(11.5, 6.5),
    )

    for mod in PRIMARY_MODALITIES:
        xs = []
        ys = []

        for i, row in enumerate(rows):
            ai = row[f"{mod}_AI"]

            if np.isfinite(ai):
                xs.append(x[i] + offsets[mod])
                ys.append(ai)

        ax.scatter(
            xs,
            ys,
            marker=markers[mod],
            s=75,
            label=mod,
        )

    ax.axhline(0, linewidth=1)

    labels = []

    for row in rows:
        short = row["Subject"].replace(
            "sub-",
            "",
        )

        if row["ThalamusDirect_GE1pct"]:
            short += "†"

        labels.append(short)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(
        "Thalamic median AI (ipsi − contra) / (ipsi + contra)"
    )
    ax.set_xlabel("Subject")
    ax.set_title(
        "Acute thalamic asymmetry by subject\n"
        "† = direct ipsilateral thalamic lesion overlap ≥1%"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        out,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    return out


def plot_subgroup_medians(summary_rows):
    out = (
        OUT_DIR
        / "PSE_ACUTE_THALAMUS_SUBGROUP_MEDIANS_v01.png"
    )

    group_order = [
        "ALL",
        "THALAMUS_SPARED_LT1PCT",
        "STRICT_DEEP_ROI_SPARED",
    ]

    group_labels = {
        "ALL": "All",
        "THALAMUS_SPARED_LT1PCT": "Thalamus spared <1%",
        "STRICT_DEEP_ROI_SPARED": "No deep ROI overlap",
    }

    x = np.arange(len(group_order), dtype=float)

    offsets = {
        "PET_FDG": -0.18,
        "ASL_CBF": 0.0,
        "CTP_CBF": 0.18,
    }

    markers = {
        "PET_FDG": "o",
        "ASL_CBF": "s",
        "CTP_CBF": "^",
    }

    lookup = {
        (r["Group"], r["Modality"]): r
        for r in summary_rows
    }

    fig, ax = plt.subplots(
        figsize=(10.5, 6.5),
    )

    for mod in PRIMARY_MODALITIES:
        xs = []
        ys = []

        for i, group in enumerate(group_order):
            row = lookup[(group, mod)]
            y = row["Median_AI_IpsiContra"]

            if np.isfinite(y):
                xs.append(x[i] + offsets[mod])
                ys.append(y)

        ax.scatter(
            xs,
            ys,
            marker=markers[mod],
            s=90,
            label=mod,
        )

    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [group_labels[g] for g in group_order]
    )
    ax.set_ylabel(
        "Median thalamic AI"
    )
    ax.set_title(
        "Lesion-aware sensitivity analysis of acute thalamic asymmetry"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        out,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    return out


def main():
    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 72)
    print("PSE ACUTE THALAMUS LESION-AWARE ANALYSIS")
    print("=" * 72)
    print(f"ROI file    : {ROI_FILE}")
    print(f"Lesion file : {LESION_FILE}")
    print(f"Output dir  : {OUT_DIR}")
    print("=" * 72)

    if not ROI_FILE.is_file():
        raise FileNotFoundError(ROI_FILE)

    if not LESION_FILE.is_file():
        raise FileNotFoundError(LESION_FILE)

    lesion = read_lesion_summary()
    roi = read_primary_thalamic_data()

    patient_rows = build_patient_table(
        lesion,
        roi,
    )

    patient_csv = write_patient_table(
        patient_rows
    )

    summary_rows = build_group_summary(
        patient_rows
    )

    summary_csv = write_group_summary(
        summary_rows
    )

    patient_png = plot_patient_level(
        patient_rows
    )

    subgroup_png = plot_subgroup_medians(
        summary_rows
    )

    print()
    print("Patient-level table:")
    print(f"  {patient_csv}")

    print("Subgroup summary:")
    print(f"  {summary_csv}")

    print("Patient-level figure:")
    print(f"  {patient_png}")

    print("Subgroup-median figure:")
    print(f"  {subgroup_png}")

    print()
    print("Key patient-level values:")

    for row in patient_rows:
        direct = (
            "DIRECT>=1%"
            if row["ThalamusDirect_GE1pct"]
            else "SPARED<1%"
        )

        print(
            f"  {row['Subject']} | {direct} | "
            f"PET={fmt(row['PET_FDG_AI'], 4)} | "
            f"ASL={fmt(row['ASL_CBF_AI'], 4)} | "
            f"CTP={fmt(row['CTP_CBF_AI'], 4)}"
        )

    print()
    print("Subgroup summary:")

    for row in summary_rows:
        print(
            f"  {row['Group']} | "
            f"{row['Modality']} | "
            f"N={row['N']} | "
            f"median AI={fmt(row['Median_AI_IpsiContra'], 4)} | "
            f"negative={row['N_AI_negative']}/{row['N']}"
        )

    print("=" * 72)
    print(
        "No inferential tests performed; "
        "this remains an exploratory pilot analysis."
    )
    print("=" * 72)


if __name__ == "__main__":
    main()
