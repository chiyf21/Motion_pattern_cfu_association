# Complete motion-to-pattern pipeline

`pipeline_from_raw_omega05_mu05/run_motion_to_patterns.py` is the complete resumable pipeline from patch-level motion to patterns. It is not a self-contained 801-GB raw-motion distribution: the raw motion fields must be supplied through an external data location, or converted to patch motion using the script's raw-H5 input interface.

Slice arguments are 0-based. A run requires a motion-file directory, an invalid/background mask, a slice index, and an output root. It writes `01_patch_motion`, `02_motion_units`, `03_episodes`, `04_modes`, and `06_patterns`. Default parameters match the canonical pattern branch and are recorded in each cache's metadata.json.
