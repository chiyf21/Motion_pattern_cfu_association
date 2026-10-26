# Slice 7 pair review folders

Each category folder contains a representative native-membrane overlay and a pair table. Batch-rendered pair figures are stored in the nested folders below:

- `01_both_methods/all_pair_overlays/`: individual PNGs for all 252 pairs co-clustered by both methods.
- `02_spatial_dtw_only/top_distance_disagreements/`: the 100 spatial/DTW-only pairs with the largest absolute difference between the two composite distances.
- `03_lag_plus_spatial_only/top_distance_disagreements/`: the 100 lag-time + spatial-only pairs ranked the same way.

The two disagreement sets are ranked by `|D_DTW+space - D_lag+space|`; each image uses the same format as the original representative examples. `render_slice07_pair_review.py` regenerates this collection on the analysis server; pass `--category both_methods`, `--category spatial_dtw_only`, or `--category lag_plus_spatial_only` to render one group at a time. The complete pair-level tables remain in each category folder as `all_pairs_in_category.csv` and in the parent as `../all_gated_mode_pairs.csv`.
