# Fig5_v0827 experimental overview

This is an independent comparison using `omega=0.5` and `mu=0.5`. It does not overwrite the legacy Fig5 patterns, CFU files, statistical results, or figures. Displayed slice numbers in filenames are 1-based; the raw motion H5 `slice` argument is 0-based.

## Experimental chain

1. `01_motion_pattern_extraction_omega05_mu05`: extract motion patterns from registered motion data.
2. `02_current_cfu_input`: freeze the current AQuA2 CFU outputs as read-only inputs.
3. `03_cfu_pattern_spatial_overlap`: identify spatial pattern–CFU modules.
4. `04_all_pattern_cfu_lag_cooccurrence`: test lagged co-occurrence for all eligible pattern×CFU pairs.
5. `05_local_mechanical_modules_distributed_ca_network`: connect spatial modules to temporally significant CFU associations.

CFU extraction is not rerun in this experiment directory. The complete motion pipeline is documented under `01_motion_pattern_extraction_omega05_mu05/pipeline_from_raw_omega05_mu05/`.

## Current results

- 2,915 patterns across 12 slices; 288 patterns have at least 5 members.
- 724 CFUs across 12 slices.
- 52 spatial modules using `pattern members≥5`, `ratio≤3`, and `coverage≥0.5`.
- Lag analysis uses frames `-8…8`, window=3, and 500 global circular-shift nulls; 168 pairs have global q<0.05.
- The module network contains 34 q<0.05 module–CFU edges involving 8 modules.
