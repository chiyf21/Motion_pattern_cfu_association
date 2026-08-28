# 03 CFU–pattern spatial overlap

This experiment is a spatial analysis independent of temporal significance. For each slice, it compares each eligible pattern mask with each CFU mask and records spatial metrics in `spatial_all_overlap_slice*.csv`. The selected module pairs are in `spatial_final_slice*.csv`, with a cross-slice summary in `summary_all_slices.csv`.

Selection requires pattern members ≥5, `ratio≤3`, and `coverage≥0.5`. `run_spatial_overlap.py` recomputes the tables. `render_all_slice_original_overlap.py` renders original-resolution spatial overviews; `render_spatial_overlap_stats.py` renders the all-slice Fig4-style summary; `render_spatial_overlap_couple_overlays.py` renders module motion-arrow overlays. Current output: 52 spatial modules. These modules are not pre-filtered by experiment 04 p-values.

The per-slice module overview figures are in `figures/original_resolution_pattern_cfu_overlap/` (`slice01...slice12...png`). Each figure uses the original field-of-view geometry and displays the union of the selected pattern masks, the union of their matched CFU masks, and their spatial intersections. The accompanying `manifest.csv` records the selected module count for each slice. The detailed arrow overlays under `figures/couple_overlays/` are retained on the server for inspection but are intentionally excluded from Git.
