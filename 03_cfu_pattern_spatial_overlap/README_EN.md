# 03: CFU–pattern spatial overlap

This experiment compares motion-pattern masks with CFU masks within each slice. Spatially corresponding pairs are defined as candidate local mechanical modules. It is parallel to the lag-co-occurrence experiment and does not use temporal significance for selection.

Patterns are restricted to those with at least five members, then filtered by `ratio≤3` and `coverage≥0.5`. `spatial_all_overlap_slice*.csv` contains all candidates and metrics; `spatial_final_slice*.csv` contains final modules; `summary_all_slices.csv` is the summary. Use `run_spatial_overlap.py` to recompute and the three `render_*` scripts to draw figures. The current result contains 52 modules under `figures/`.
