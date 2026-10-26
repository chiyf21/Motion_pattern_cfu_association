# Q2–Q4 analyses for the cluster055 motion–Ca project

This folder contains three distinct analyses and their current code/results. The detailed, version-specific methods and input checklist are in [`Q2_Q4_METHODS_EN.md`](Q2_Q4_METHODS_EN.md). Start there before rerunning any analysis.

## Current entry points

- **Q2 — mode metric comparison:** `Q2/results_lag_plus_spatial/`. It compares canonical spatial+DTW clustering with spatial+lag-correlation clustering for the already-retained Slice 7 modes.
- **Q3 — local membrane appearance mismatch:** `Q3/results/` and `Q3/figures/`. It evaluates CFU-region membrane MSE around significant pair-hit onsets, other Ca onsets, and event-free controls. The MSE is an appearance-mismatch proxy, not ground-truth registration error.
- **Q4 — continuous response (primary M2):** `Q4/continuous_response_20261026/M2_main_results_png/`. This is the curated PNG-only figure set: one M2-vs-M0 significance overview, one M2-kernel overview, and one observed-response/kernel PNG for each q-significant pair.

Q4's saved tables, run configuration, and runner remain under `Q4/continuous_response_20261026/`. The older Q4 pair-figure gallery and stability-diagnostic folder were removed from this curated package; they are not needed to regenerate the current PNG set from the saved tables.

## Inputs needed to reproduce

The scripts and figures alone are not a complete rerun package. A collaborator needs the Slice 7 Stage-06 mode cache, `/home/cyf/wbi/wbi_code/src_registration/motion_correlation_pattern_v2.py`, and registered membrane images for Q2; the relevant cluster055 motion-pattern caches and compatible `wholistic_registration` source; the current AQuA2 CFU MAT files; the cluster055 temporal significant-pair table/code; the 1,598 registered membrane-channel volumes and saved reference volumes for Q3; and the frozen spatial module–CFU association table plus Q3 metrics table for Q4. The CFU entries in the project are symlinks, so transfer their target MAT files rather than the links alone. Full server paths and each file's role are itemized in `Q2_Q4_METHODS_EN.md` and the relevant `run_config.json` files.

Server paths currently embedded in historical configs are machine-specific. For a collaborator, provide the data at an agreed location and pass its paths explicitly to the runners. Do not infer that all Q2–Q4 inputs are included in this code/results folder.

## Reproducibility status

The runners pass syntax and CLI checks on the analysis server. The Q2 missing helper import and Q3 full-frame runner's implicit CFU-root selection were corrected in this folder. Q4's historical CLI defaults do not all match the paths recorded for the saved run; reruns must use the explicit paths in its `run_config.json`. A full rerun has not been performed as part of this packaging/cleanup pass.
