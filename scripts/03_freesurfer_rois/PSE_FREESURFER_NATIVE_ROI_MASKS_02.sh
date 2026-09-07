#!/bin/zsh
# Map FreeSurfer aseg.mgz back to the native T1 grid and create bilateral
# thalamic/hippocampal binary masks for every completed acute FreeSurfer subject.

export FREESURFER_HOME="${FREESURFER_HOME:-/Applications/freesurfer/8.2.0}"
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"

ROOT="${PSE_ROOT:-$HOME/PSE}"
export SUBJECTS_DIR="$ROOT/derivatives/freesurfer"
OUTROOT="$ROOT/derivatives/roi_masks_native"
LOGROOT="$ROOT/derivatives/roi_masks_native_logs"
FSBIN="$FREESURFER_HOME/bin"
mkdir -p "$OUTROOT" "$LOGROOT"

subjects=(${(f)"$(find "$SUBJECTS_DIR" -maxdepth 1 -type d -name 'sub-P*_acute' -exec basename {} \; 2>/dev/null | sed 's/_acute$//' | sort)"})
if (( ${#subjects[@]} == 0 )); then
    echo "ERROR: no completed/partial sub-P*_acute directories found under $SUBJECTS_DIR"
    exit 1
fi

echo "============================================================"
echo "PSE FREESURFER -> NATIVE ROI MASKS"
echo "Subjects: ${subjects[*]}"
echo "============================================================"

for SUB in "${subjects[@]}"; do
    FSID="${SUB}_acute"
    FSDIR="$SUBJECTS_DIR/$FSID"
    ASEG="$FSDIR/mri/aseg.mgz"
    RAWAVG="$FSDIR/mri/rawavg.mgz"
    OUTDIR="$OUTROOT/$SUB/acute"
    LOGFILE="$LOGROOT/${SUB}_acute.log"
    mkdir -p "$OUTDIR"

    if [[ ! -f "$ASEG" || ! -f "$RAWAVG" ]]; then
        echo "$SUB: missing aseg.mgz or rawavg.mgz -> skip" | tee "$LOGFILE"
        continue
    fi

    ASEG_NATIVE="$OUTDIR/aseg_in_native_T1.nii.gz"
    "$FSBIN/mri_label2vol" --seg "$ASEG" --temp "$RAWAVG" --o "$ASEG_NATIVE" --regheader "$ASEG" 2>&1 | tee "$LOGFILE"

    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 10 --o "$OUTDIR/Left_Thalamus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"
    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 17 --o "$OUTDIR/Left_Hippocampus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"
    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 49 --o "$OUTDIR/Right_Thalamus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"
    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 53 --o "$OUTDIR/Right_Hippocampus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"
    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 10 49 --o "$OUTDIR/Bilateral_Thalamus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"
    "$FSBIN/mri_binarize" --i "$ASEG_NATIVE" --match 17 53 --o "$OUTDIR/Bilateral_Hippocampus_mask.nii.gz" --uchar 2>&1 | tee -a "$LOGFILE"

    echo "$SUB: DONE" | tee -a "$LOGFILE"
done
