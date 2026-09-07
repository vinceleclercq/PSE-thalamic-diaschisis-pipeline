#!/bin/sh
# Copy to a private local file (e.g. config.local.sh), edit the paths, then:
#   source config.local.sh

export PSE_ROOT="$HOME/PSE"
export SPM_DIR="$HOME/spm25"
export FREESURFER_HOME="/Applications/freesurfer/8.2.0"
export PSE_FSL_DATA_ENV="$HOME/miniforge3/envs/pse-fsl-data"

# FreeSurfer also requires a valid FS_LICENSE path in your environment.
# export FS_LICENSE="$HOME/.freesurfer_license.txt"
