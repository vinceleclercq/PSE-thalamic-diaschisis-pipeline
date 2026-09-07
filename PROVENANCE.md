# Provenance and curation notes

This repository is a cleaned, publication-facing snapshot of the code used to develop the acute PSE thalamic diaschisis pilot analysis.

## What was changed for public release

The public copy was curated to remove workstation-specific paths and development-only clutter while preserving the analytical logic of the final pipeline.

Changes include:

- replaced local absolute paths with environment variables (`PSE_ROOT`, `SPM_DIR`, `FREESURFER_HOME`, `PSE_FSL_DATA_ENV`);
- retained only the final/relevant versions of duplicate MATLAB/QC scripts;
- consolidated the two final NVAUTO preprocessing stages into `PSE_NVAUTO_SEGMENTATION_FINAL.py`;
- renamed public NVAUTO output folders to avoid the historical but inaccurate implication that the full DeepISLES ensemble was used;
- added a reproducible ADC reconstruction script based on the b=0/b=1000 reconstruction used when vendor ADC maps were absent;
- excluded patient data, images, lesion masks and generated analysis tables.

## NVAUTO development provenance

During development, two scripts were run sequentially:

1. `PSE_DEEPISLES_INFER_03.py`
   - harmonized ADC to the exact DWI grid;
   - forced identical DWI/ADC qform/sform geometry;
   - ran an initial `stroke_segmentor` inference.

2. `PSE_NVAUTO_BRAINMASK_TEST_01.py`
   - projected FreeSurfer aseg to DWI space;
   - created an in-plane cleaned/dilated brain mask;
   - zeroed extracerebral DWI/ADC signal;
   - reran NVAUTO on CPU;
   - generated the visually accepted final lesion masks.

The public script `PSE_NVAUTO_SEGMENTATION_FINAL.py` consolidates those final operations into a single reproducible implementation. The model used is NVAUTO from `stroke_segmentor`, not the complete DeepISLES ensemble.

## Development scripts intentionally not included in the main pipeline

The following local scripts represented intermediate or rejected methods and should not be cited as the final analysis method:

- `PSE_LESION_CANDIDATE_DWI_ADC_01.m`
- `PSE_LESION_CANDIDATE_DWI_ADC_02.m`
- `PSE_LESION_CANDIDATE_DWI_ADC_03.m`
- `PSE_LESION_CANDIDATE_DWI_ADC_04.m`
- early `PSE_DEEPISLES_INFER_0x.py` variants
- manual lesion-edit preparation helpers
- patient-specific artifact investigation scripts
- superseded v01 ROI/FreeSurfer scripts when a corrected v02/v03 version exists

These development steps were useful for QC and method selection but did not define the final lesion masks used in the lesion-aware analysis.

## Pilot-specific analytical decisions

Some scripts intentionally preserve pilot-specific analytical decisions, including:

- subject identifiers are discovered from the local project structure or local side/QC tables rather than hard-coded in the public scripts;
- ADC reconstruction is triggered for subjects without a vendor ADC map, using the b=0/b=1000 file mapping only after source-DICOM verification of b-values;
- direct ipsilateral thalamic overlap >=1% as the main descriptive sensitivity threshold;
- no inferential p-values in the 7-subject exploratory analyses;
- ASL interpreted through relative/asymmetry measures because absolute units were not formally verified.

These decisions should be revisited before applying the repository unchanged to a larger/final cohort.
