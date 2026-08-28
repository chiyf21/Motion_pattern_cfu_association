# Canonical motion-pattern extraction

This replaces the previous Stage-6-only rerun. It starts with registered
`motion_*.h5` displacement fields and an invalid/background mask in the same raw H×W
geometry, then stores every reusable intermediate result for a **0-based**
slice.

The code imports the pinned clone:
`Fig5/wholistic_registration` commit `2b3c4e611ca194d391a31b56f6bc93a28ca90b13`.

## Stages

1. `01_patch_motion`: patch-averaged displacement, velocity, and patch mask.
2. `02_motion_units`: local-MAD velocity events and filtered motion units.
3. `03_episodes`: spatial-temporal episode grouping and artifact filtering.
4. `04_modes`: sparse-compact velocity modes (`Kmax=8`, `lambda_sc=0.05`,
   `rho=1`, `kappa=4`).
5. `06_patterns`: direct mode clustering with `omega=mu=0.5`, complete linkage,
   `cluster_dist_thresh=0.45`, and best-connected-component unified masks.

Stage 05 is deliberately absent: the canonical analysis clusters modes directly,
not regions.

## Run

```bash
python run_motion_to_patterns.py \
  --raw-motion-dir /path/to/registered/motion \
  --invalid-mask-npy /path/to/invalid_background_mask.npy \
  --slice 6 \
  --output-root /path/to/output \
  --through-stage 06
```

The raw H5 `slice` is 0-based. Output is resumable: re-running skips a completed
stage. Use `--overwrite` only when intentionally regenerating the requested and
downstream stages. `metadata.json` records paths, clone commit, and parameters.

The mask is **True for invalid/background pixels**. The current analysis uses
`valid_mask = ~mask_patched` because legacy cached `mask_patched` records
invalid/background patches; this compatibility convention is preserved here.
