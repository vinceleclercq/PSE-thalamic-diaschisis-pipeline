#!/bin/zsh
# Run cross-sectional FreeSurfer recon-all for every acute T1 found under
# $PSE_ROOT/derivatives/nifti/sub-P*/acute/T1/T1.nii.
# Existing completed subjects are skipped. Partial subjects stop the batch for
# review rather than being overwritten.

set -u

export FREESURFER_HOME="${FREESURFER_HOME:-/Applications/freesurfer/8.2.0}"
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"

ROOT="${PSE_ROOT:-$HOME/PSE}"
export SUBJECTS_DIR="$ROOT/derivatives/freesurfer"
LOGDIR="$ROOT/derivatives/freesurfer_batch_logs"
mkdir -p "$SUBJECTS_DIR" "$LOGDIR"

if [[ -z "${FS_LICENSE:-}" || ! -f "$FS_LICENSE" ]]; then
    echo "ERROR: FS_LICENSE is not set or invalid."
    exit 1
fi

subjects=(${(f)"$(find "$ROOT/derivatives/nifti" -maxdepth 1 -type d -name 'sub-P*' -exec basename {} \; 2>/dev/null | sort)"})
if (( ${#subjects[@]} == 0 )); then
    echo "ERROR: no sub-P* folders found under $ROOT/derivatives/nifti"
    exit 1
fi

MEM_BYTES=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
MEM_GB=$(( MEM_BYTES / 1024 / 1024 / 1024 ))
if (( MEM_GB > 0 && MEM_GB < 32 )); then
    export FS_V8_XOPTS=0
else
    unset FS_V8_XOPTS 2>/dev/null || true
fi

echo "============================================================"
echo "PSE FREESURFER ACUTE BATCH"
echo "FreeSurfer: $(recon-all -version 2>&1 | tail -1)"
echo "SUBJECTS_DIR: $SUBJECTS_DIR"
echo "Subjects: ${subjects[*]}"
echo "============================================================"

for SUB in "${subjects[@]}"; do
    FSID="${SUB}_acute"
    T1="$ROOT/derivatives/nifti/$SUB/acute/T1/T1.nii"
    MASTERLOG="$LOGDIR/${FSID}_batch.log"
    DONEFILE="$SUBJECTS_DIR/$FSID/scripts/recon-all.done"

    echo "\n------------------------------------------------------------"
    echo "$FSID"
    echo "------------------------------------------------------------"

    if [[ -f "$DONEFILE" ]]; then
        echo "Already complete -> skip"
        continue
    fi
    if [[ ! -f "$T1" ]]; then
        echo "Missing T1 -> skip: $T1"
        continue
    fi
    if [[ -d "$SUBJECTS_DIR/$FSID" ]]; then
        echo "ERROR: partial subject directory already exists: $SUBJECTS_DIR/$FSID"
        echo "Review it before rerunning this subject."
        exit 1
    fi

    FOV=$("$FREESURFER_HOME/bin/mri_info" "$T1" 2>/dev/null | awk '/^[[:space:]]*fov:/ {print $2; exit}')
    if [[ -z "$FOV" ]]; then
        echo "ERROR: could not determine FOV for $T1"
        exit 1
    fi

    EXTRA=()
    if awk -v f="$FOV" 'BEGIN {exit !(f > 256.0)}'; then
        EXTRA=(-cw256)
        echo "FOV > 256 mm -> adding -cw256"
    fi

    echo "START: $(date '+%Y-%m-%d %H:%M:%S')" | tee "$MASTERLOG"
    recon-all -s "$FSID" -i "$T1" "${EXTRA[@]}" -all 2>&1 | tee -a "$MASTERLOG"
    RC=${pipestatus[1]}
    echo "END: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$MASTERLOG"
    echo "EXIT CODE: $RC" | tee -a "$MASTERLOG"

    if [[ $RC -ne 0 ]]; then
        echo "ERROR: recon-all failed for $FSID; stopping for review."
        exit $RC
    fi
done

echo "\nBATCH FINISHED"
