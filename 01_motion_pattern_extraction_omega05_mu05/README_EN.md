# 01: Motion pattern extraction

## Purpose

Starting from registered raw motion displacement fields, this experiment generates patch motion, motion units, motion episodes, sparse-compact motion modes, and finally motion patterns by clustering modes. With `omega=0.5` and `mu=0.5`, the clustering distance incorporates both spatial and motion features.

## Results and code

The directory contains `06_patterns/objects.pkl` for 12 slices: 2,915 patterns in total, including 288 patterns with at least five members. `rerun_patterns_omega05_mu05.py` starts from existing `04_modes`; the complete raw-motion workflow is under `pipeline_from_raw_omega05_mu05/`, with resumable stage-wise caches.

Key parameters: patch size 7; velocity is frame-to-frame motion difference; modes use `Kmax=8`, `svd_target_r2=0.90`, `lambda_sc=0.05`, `rho=1`, and `kappa=4`; pattern clustering uses complete linkage, `cluster_dist_thresh=0.45`, `min_iou=0.08`, and best-connected-component unified masks.
