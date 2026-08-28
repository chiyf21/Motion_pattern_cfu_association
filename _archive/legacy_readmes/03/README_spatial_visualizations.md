# motion_pattern_cfu_association spatial module and Fig4-style annotation outputs

## Stage 03: spatial module extraction

`03_cfu_pattern_spatial_overlap/run_spatial_overlap.py` recomputes spatial
pattern-CFU modules only. It uses the v0827 omega=0.5/mu=0.5 pattern files,
linked current CFUs, pattern members >=5, ratio <=3 and coverage >=0.5.

## Stage 03.1: Fig4-style spatial visualizations

`03_cfu_pattern_spatial_overlap/render_spatial_overlap_stats.py` produces the
all-slice relative-area versus coverage scatter plot in PNG/PDF plus its CSV.
`render_spatial_overlap_couple_overlays.py` produces the per-module motion-arrow
overlays and a manifest. Run both with `run_spatial_overlap_visualizations.sh`.

These figures use only v0827 spatial modules and v0827 pattern objects; they
do not use archived Fig5 spatial or pattern data.
