# Slice 7: spatial + DTW versus spatial + lag correlation

## Purpose

This controlled comparison changes only the temporal mode-distance component. It does not re-extract motion modes, alter production patterns, or use CFUs.

## Inputs and method

- Stage-06 Slice 7 cache: `objects.pkl` and `distance_matrix.npz`.
- 1,034 cached modes; the runner verifies object/matrix ordering and reproduces the source partition at cut 0.45.
- Shared compatibility gate: different episodes and mask IoU ≥ 0.08. The saved run retained 28,321 of 534,061 unordered pairs.
- Canonical distance: `0.5 × D_h + 0.5 × D_b`, where `D_h` is sign-aware DTW and `D_b` is the unchanged cached response-field spatial distance.
- Alternative: `0.5 × D_lag + 0.5 × D_b`, where `D_lag = 1 − max |Pearson r|` over integer lags −8…+8. Minimum overlap is `max(3, ceil(0.8 × shorter trace))`.
- Both use complete linkage and cut 0.30. A 0.20–0.45 threshold sensitivity is also saved. Equal numeric cuts need not imply equal effective stringency.

## Reproduce on the current server

```bash
ROOT=/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055/06_methodological_review_q2_q4
python "$ROOT/Q2/code/run_slice07_dual_clustering.py" \
  --cache /home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice07_velocity_decomp/06_patterns \
  --membrane-dir /mnt/data21T_2/cyf/f338/f338_registrated_0530/membrane \
  --out "$ROOT/Q2/results_lag_plus_spatial_reproduced"
```

The output directory must be empty. Add `--overwrite-results` only when intentionally regenerating an existing output. The Q2 source module `src_registration/motion_correlation_pattern_v2.py` is currently at `/home/cyf/wbi/wbi_code/src_registration/`.

## Outputs

`results_lag_plus_spatial/` contains the pair-level distances and lag information, mode/cluster assignments, threshold sensitivity, summary validation, and pair-review overlays. `pair_review/` includes all 252 pairs co-clustered by both methods and the 100 largest distance disagreements in each one-sided category. The older pure lag-only diagnostic is not part of this corrected comparison folder.
