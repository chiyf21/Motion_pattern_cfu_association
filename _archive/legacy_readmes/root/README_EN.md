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

The code should be run in the order 01 → 03 and 04 → 05. Experiment 03 and experiment 04 are parallel analyses; 05 depends on both. See the English and Chinese README in each directory for exact inputs, outputs, parameters, and command entry points.
