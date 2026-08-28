#!/usr/bin/env bash
set -euo pipefail

BASE=/home/cyf/wbi/wbi_code
ROOT=$BASE/experiments/motion_pattern_cfu_association
LOG=$ROOT/logs
mkdir -p "$LOG"
cd "$BASE"

echo "[$(date -Is)] Stage 01: motion patterns omega=0.5 mu=0.5" | tee "$LOG/pipeline.log"
pattern_pids=()
python "$ROOT/01_motion_pattern_extraction_omega05_mu05/rerun_patterns_omega05_mu05.py" 1 5 9 > "$LOG/pattern_group_1.log" 2>&1 & pattern_pids+=("$!")
python "$ROOT/01_motion_pattern_extraction_omega05_mu05/rerun_patterns_omega05_mu05.py" 2 6 10 > "$LOG/pattern_group_2.log" 2>&1 & pattern_pids+=("$!")
python "$ROOT/01_motion_pattern_extraction_omega05_mu05/rerun_patterns_omega05_mu05.py" 3 7 11 > "$LOG/pattern_group_3.log" 2>&1 & pattern_pids+=("$!")
python "$ROOT/01_motion_pattern_extraction_omega05_mu05/rerun_patterns_omega05_mu05.py" 4 8 12 > "$LOG/pattern_group_4.log" 2>&1 & pattern_pids+=("$!")
for pid in "${pattern_pids[@]}"; do wait "$pid"; done

echo "[$(date -Is)] Temporal co-occurrence" | tee -a "$LOG/pipeline.log"
python "$ROOT/04_all_pattern_cfu_lag_cooccurrence/run_all_cfu_pattern_lag8_w3_pairwise.py" > "$LOG/temporal.log" 2>&1

echo "[$(date -Is)] Spatial modules, network and gallery" | tee -a "$LOG/pipeline.log"
python "$ROOT/05_local_mechanical_modules_distributed_ca_network/run_module_network.py" > "$LOG/spatial_and_network.log" 2>&1
python "$ROOT/05_local_mechanical_modules_distributed_ca_network/render_module_cfu_spatial_gallery.py" > "$LOG/gallery.log" 2>&1
python "$ROOT/03_cfu_pattern_spatial_overlap/render_spatial_overlap_stats.py" > "$LOG/fig4_style_stats.log" 2>&1
python "$ROOT/03_cfu_pattern_spatial_overlap/render_spatial_overlap_couple_overlays.py" > "$LOG/fig4_style_overlays.log" 2>&1

echo "[$(date -Is)] COMPLETE" | tee -a "$LOG/pipeline.log"
