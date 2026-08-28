#!/usr/bin/env bash
set -euo pipefail

BASE=/home/cyf/wbi/wbi_code
ROOT=$BASE/experiments/motion_pattern_cfu_association
cd "$BASE"

python "$ROOT/03_cfu_pattern_spatial_overlap/render_spatial_overlap_stats.py"
python "$ROOT/03_cfu_pattern_spatial_overlap/render_spatial_overlap_couple_overlays.py"
