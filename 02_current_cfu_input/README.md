# 02 Current CFU input

This directory freezes the AQuA2 CFU files used by experiments 03–05. It contains one link per displayed slice (slice01–slice12). The links are not portable because they currently point to MAT files outside this repository on the laboratory server.

The inputs come from the ds7 AQuA2 native event-detection and CFU aggregation run, with a minimum of five events per CFU. No AQuA2 code is run here. To reproduce this input on another machine, copy the target MAT files into a local data directory and update the links or a local path configuration; do not commit private data to GitHub.

The reproducible source pipeline is documented in `aqua2_native_pipeline/`. It contains the calcium preparation script, the native AQuA2 DS7 event-detection script, the CFU aggregation script, and an explicit all-slice runner. The runner reuses existing event outputs and skips existing CFU files unless `FORCE=1` is set.
