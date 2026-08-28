# 01 Motion pattern extraction

This is the only stage that creates motion patterns. `rerun_patterns_omega05_mu05.py` is a fast comparison script that starts from existing mode caches. For a fresh run from motion data, use `pipeline_from_raw_omega05_mu05/run_motion_to_patterns.py`; it writes resumable stage directories `01_patch_motion`, `02_motion_units`, `03_episodes`, `04_modes`, and `06_patterns` for one 0-based slice.

The canonical branch uses patch size 7, frame-to-frame velocity for event detection and mode decomposition, local-MAD rest motion (`window_t=21`, `window_xy=3`), episode artifact filtering, sparse-compact modes (`Kmax=8`, SVD target R²=.90, λ=.05, ρ=1, κ=4), and direct mode clustering with complete linkage, `min_iou=.08`, `cluster_dist_thresh=.45`, `omega=mu=.5`, and best-connected-component unified masks. The pinned implementation is the submodule at `../wholistic_registration/`.

The current final objects are under `patterns/SliceXX_velocity_decomp/06_patterns/objects.pkl`. There are 2,915 patterns across 12 slices and 288 patterns with at least five members. Large objects are intentionally ignored by Git; provide them separately if downstream users should analyze the existing run.

For a compact visual check of the eligible patterns, run `render_pattern_overview_members_ge5.py`. It reads the current pattern objects, keeps only patterns with `n_members >= 5`, and writes one original-resolution overview per slice plus a contact sheet and `manifest.csv` to `patterns/_overview_members_ge5/`. These images are the canonical `members>=5` pattern overview outputs; the underlying pattern caches remain excluded from Git.
