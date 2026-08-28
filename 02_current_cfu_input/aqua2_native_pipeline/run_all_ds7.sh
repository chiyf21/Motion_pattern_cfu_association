#!/usr/bin/env bash
set -euo pipefail

# Portable placeholder (edit for another machine):
BASE_DIR="/path/to/your/aqua_validation"
# Current server example:
BASE_DIR="/home/cyf/wbi/wbi_code/experiments/_trash/roi_ca_dependency/aqua_validation"

# Portable placeholder (edit for another machine):
AQUA2_DIR="/path/to/your/AQuA2"
# Current server example:
AQUA2_DIR="/home/cyf/wbi/wbi_code/experiments/_trash/roi_ca_dependency/AQuA2"

MATLAB="matlab"
WORKERS=8
FORCE="${FORCE:-0}"
SLICES=(1 2 3 4 5 6 7 8 9 10 11 12)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for z in "${SLICES[@]}"; do
  export AQUA_Z="$z"
  export AQUA_BASE_DIR="$BASE_DIR"
  export AQUA2_DIR
  export AQUA_POOL_WORKERS="$WORKERS"
  "$MATLAB" -nodisplay -nosplash -nodesktop -r "run('$SCRIPT_DIR/02_run_event_detection_ds7.m'); exit;"
  if [[ "$FORCE" == "1" || ! -f "$BASE_DIR/output/slices/cfu_ds7_all_ot030_min5_group5/slice_Z$(printf '%02d' "$z")_ds7_native_CFU_ot030_min5_group5.mat" ]]; then
    "$MATLAB" -nodisplay -nosplash -nodesktop -r "run('$SCRIPT_DIR/03_run_cfu_aggregation_ds7.m'); exit;"
  else
    echo "Skipping existing CFU output for Z$(printf '%02d' "$z")"
  fi
done
