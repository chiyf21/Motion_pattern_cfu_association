# Results index

This repository separates reproducible code from large experimental inputs and caches.

## Inputs to provide separately

- Patch-motion arrays for each 0-based slice, containing `motion_delta`, `motion_abs`, and `mask_patched` in `arrays.npz` format.
- The 12 current AQuA2 CFU MAT files referenced by `02_current_cfu_input/cfu/`. The symbolic links are server-specific; provide the resolved MAT files separately.
- The reference TIFF only when rendering original-resolution figures.
- Pattern object caches under `01_motion_pattern_extraction_omega05_mu05/patterns/SliceXX_velocity_decomp/06_patterns/objects.pkl` if the collaborator should reuse the exact current run without recomputing.

## Pipeline and outputs

1. `01_motion_pattern_extraction_omega05_mu05/` creates motion patterns. The `members>=5` overview figures are tracked under `patterns/_overview_members_ge5/`.
2. `03_cfu_pattern_spatial_overlap/` computes spatial pattern–CFU modules. The 12 tracked original-resolution module overviews are under `figures/original_resolution_pattern_cfu_overlap/`.
3. `04_all_pattern_cfu_lag_cooccurrence/` tests all eligible pattern×CFU pairs at lags -8…+8 with a three-frame window. Large result tables and overlays are intentionally external; an example can be regenerated with the supplied render scripts.
4. `05_local_mechanical_modules_distributed_ca_network/` joins spatial modules and temporal associations.

Run the complete analysis with `bash run_all_pipeline.sh` after editing `config.py` and supplying the external inputs. CFU extraction itself is not rerun by this repository.

All paths and parameters that commonly need adjustment are documented in `config.py`; the scripts also expose stage-specific command-line options where appropriate.
