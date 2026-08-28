# Motion–calcium spatiotemporal association analysis

This repository contains the independent `omega=0.5, mu=0.5` analysis branch for testing associations between registered motion patterns and AQuA2 calcium functional units (CFUs). It is organized as a reproducible analysis package, not as a single figure-generation folder.

## What is included

`01_motion_pattern_extraction_omega05_mu05/` contains the motion-pattern pipeline and derived pattern objects. The complete workflow starts from patch-level motion outputs and proceeds through motion units, episodes, sparse-compact modes, and mode-based patterns.

`02_current_cfu_input/` contains the frozen CFU inputs used by every downstream analysis. The links currently point to the laboratory server; the target MAT files are not stored in Git.

`03_cfu_pattern_spatial_overlap/` tests spatial correspondence between pattern masks and CFU masks. Its outputs define candidate local mechanical modules using pattern member count, ratio, and coverage thresholds.

`04_all_pattern_cfu_lag_cooccurrence/` tests temporal association independently of spatial modules. It evaluates all eligible pattern×CFU pairs across lags -8 to +8 with a three-frame window and empirical global FDR.

`05_local_mechanical_modules_distributed_ca_network/` joins the spatial modules from 03 with significant temporal CFU associations from 04 to examine distributed Ca signals associated with each local mechanical module.

`wholistic_registration/` is the pinned Git submodule containing the registration and motion-analysis implementation used by the complete motion pipeline.

## Current result snapshot

The current derived results contain 2,915 patterns across 12 slices, 288 patterns with at least five members, and 724 CFUs. Spatial filtering gives 52 candidate modules. The lag analysis gives 168 globally significant pairs at q<0.05. The module-centric network contains 34 q<0.05 module–CFU edges involving 8 modules.

## Reproducibility and data boundary

The repository deliberately excludes large generated caches, raw motion H5 files, pickle objects, NPZ arrays, MAT inputs, and lag result tables. Raw registered motion is approximately 801 GB; it must be provided through a mounted data location or a separate data archive. The current pattern objects are approximately 25 GB and the lag results approximately 343 MB.

## Files to provide to collaborators

The intended handoff starts from patch-level motion, not from raw motion H5 files. For each displayed slice, provide `01_motion_pattern_extraction_omega05_mu05/patch_motion/SliceXX_velocity_decomp/01_patch_motion/arrays.npz` and its `metadata.json`. The NPZ must contain `motion_delta`, `motion_abs`, and `mask_patched`; these are patch-grid arrays, not full-resolution movies. The legacy convention is that `mask_patched=True` denotes background/invalid patches, so downstream code uses its complement as the valid analysis mask.

For users who only need to analyze the existing patterns, also provide `patterns/SliceXX_velocity_decomp/06_patterns/objects.pkl`, `distance_matrix.npz` when present, and the matching `metadata.json`. For users who need to regenerate patterns from patch motion, provide the stage caches `02_motion_units/`, `03_episodes/`, and `04_modes/`, or run the stage scripts in `pipeline_from_raw_omega05_mu05/` after supplying a patch-motion-compatible stage-1 cache. The complete pattern objects are large and are intentionally outside Git.

The CFU handoff consists of the 12 actual MAT files represented by the links in `02_current_cfu_input/cfu/`. The links themselves are server-specific and are not sufficient for another machine. These CFUs are the ds7 AQuA2 native outputs with a minimum of five events per CFU.

The spatial analysis reads pattern objects plus CFU MAT files and writes `03_cfu_pattern_spatial_overlap/spatial_all_overlap_slice*.csv`, `spatial_final_slice*.csv`, and `summary_all_slices.csv`. The temporal analysis reads the same pattern objects and CFU MAT files and writes the tables under `04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset/`. The module-network analysis reads the spatial final tables and temporal significant-pair table, then writes outputs under `05_local_mechanical_modules_distributed_ca_network/`.

In short: patch motion is the motion input; CFU MAT files are the calcium input; 03 and 04 are parallel analyses; 05 consumes both. No raw motion H5, raw calcium movie, or AQuA2 rerun is required to reproduce the downstream analyses once these derived inputs are supplied.

The code should be run in the order 01 → 03 and 04 → 05. Experiment 03 and experiment 04 are parallel analyses; 05 depends on both. See the English and Chinese README in each directory for exact inputs, outputs, parameters, and command entry points.

For a single place to check the handoff and current tracked figures, see `RESULTS_INDEX.md` (Chinese: `RESULTS_INDEX_CN.md`). Before running on another machine, edit the active paths and parameters in `config.py`; the file contains both portable placeholder lines and the current server examples. `run_all_pipeline.sh` is repository-relative and now includes the spatial-module computation explicitly.
