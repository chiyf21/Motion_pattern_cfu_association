# 02 CFU extraction pipeline

This directory documents and stores the code for the current frozen CFU input. The canonical result is the DS7 patch, native AQuA2 2-D pipeline, with at least five events per CFU. The existing `.mat` files in `../cfu/` are outputs, not scripts.

## Exact processing chain

```text
patch calcium input
  -> 01_prepare_aqua_input.py
  -> AQuA2 event detection
       burst.prep1
       pre.baselineRemoveAndNoiseEstimation
       act.acDetect
       se.seDetection
       evt.se2evtTop
       fea.getFeaturesTop
  -> 03_run_cfu_aggregation_ds7.m
       cfu.CFUdetectScript
       cfu.calAllDependencyScript
       cfu.groupCFUscript
  -> CFU MAT output
```

The AQuA2 event-detection code is `02_run_event_detection_ds7.m`. It uses a bounded MATLAB pool and reads the 0-based slice from `AQUA_Z`. The CFU code is `03_run_cfu_aggregation_ds7.m`; it reuses the event result and does not rerun detection.

## Current parameters

Event detection: `regMaskGap=5`, `frameRate=2.0`, `thrARScl=3`, `minDur=5`, `minSize=10`, `maxDelay=0.6`, `needRefine=0`, `needGrow=0`, `sourceSensitivity=8`, `whetherExtend=1`, `cut=800`, `movAvgWin=25`, `correctTrend=1`, `gapExt=5`.

CFU aggregation: `overlapThr1=0.30`, `overlapThr2=0.50`, `minNumEvt1=5`, `minNumEvt2=3`, `maxDist=10`, `shift=0`, `pValueThr=1e-5`, `cfuNumThr=5`.

## Running one slice

1. Edit the active `baseDir` and AQuA2 startup path in `02_run_event_detection_ds7.m` and `03_run_cfu_aggregation_ds7.m`. Each file contains a `/path/to/your/...` placeholder followed by the current server example.
2. Set the displayed slice index. `AQUA_Z=1` means the first displayed slice in the current file convention; the underlying AQuA2 slice index is passed explicitly to the script.
3. Run event detection, then CFU aggregation:

```bash
export AQUA_Z=6
matlab -nodisplay -nosplash -nodesktop -r "run('02_run_event_detection_ds7.m'); exit;"
matlab -nodisplay -nosplash -nodesktop -r "run('03_run_cfu_aggregation_ds7.m'); exit;"
```

The scripts print input, output, event count, CFU count, worker count and elapsed seconds. Change only the parameter block when testing a new AQuA2 variant; write its output to a new directory instead of overwriting the frozen result.

## Running all slices

`run_all_ds7.sh` runs the two MATLAB stages for slices 1–12. It is intentionally explicit and easy to edit: change `MATLAB`, `BASE_DIR`, `AQUA2_DIR`, `WORKERS`, or the slice list at the top. It does not rerun CFU extraction when a matching output already exists unless `FORCE=1` is set.

## Data boundary

The raw calcium movie, intermediate event `.mat` files and final CFU `.mat` files are not committed to GitHub. The 12 frozen CFU files in `../cfu/` must be supplied separately for downstream analyses. The AQuA2 project itself is also an external dependency and must be supplied or cloned separately.
