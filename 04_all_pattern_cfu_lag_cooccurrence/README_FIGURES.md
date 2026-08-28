# motion_pattern_cfu_association critical visual outputs

All scripts in this folder use the omega=0.5, mu=0.5 v0827 patterns and the
linked current CFUs. They never read the legacy `Fig5/23` pattern objects.

## All significant temporal pairs

Run `render_significant_pairs_overlay_fast.py`. It renders every pair in
`empirical_FDR_significant_pairs.csv` to
`04_all_pattern_cfu_lag_cooccurrence/figures/significant_pairs_spatiotemporal_original_overlay_fast/`.
Each figure contains both original-reference spatial panels plus the complete
0-based timeline.

## One article-style pair

`render_article_pair.py` is generic. For example:

```bash
PATTERN_SLICE=8 PATTERN_ID=200 CFU_SLICE=8 CFU_ID=2 \
python render_article_pair.py
```

It saves separate pattern/CFU spatial-strength images and aligned timelines.
The requested IDs must exist in the v0827 pattern/CFU tables; the script is a
renderer and does not alter statistical results.

## One native-resolution overlap image per slice

Run `render_all_slice_original_overlap.py`. It creates 12 original-pixel-size
reference overlays in the v0827 spatial result folder. Magenta is the union of
selected pattern masks, cyan the union of selected paired CFU masks, and yellow
their intersection.
