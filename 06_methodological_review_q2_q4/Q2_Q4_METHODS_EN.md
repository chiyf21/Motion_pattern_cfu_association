# Q2–Q4 methodological record (cluster055)

This note records the current implementations in `06_methodological_review_q2_q4/`. The three questions are separate diagnostics; none is a causal analysis. Q2 and Q3 are the current controlled comparison and registration-consistency analyses. Q4 is the primary continuous-response analysis and is presented through `Q4/continuous_response_20261026/M2_main_results_png/`.

## Q2 — DTW versus lagged-correlation mode clustering

### Question

Do two temporal similarity definitions group the same already-gated Slice 7 motion modes? This is a sensitivity comparison of the mode-clustering metric. It does not re-extract modes, alter the production pattern set, or use downstream CFUs.

### Inputs and population

- The Stage-06 Slice 7 cache: `objects.pkl` (retained mode objects and cached pair metadata/distance matrix) and `distance_matrix.npz`.
- The cached object order and matrix order are checked. The runner also verifies that the cached canonical labels are reproduced at the original composite-distance cut of 0.45 before comparing alternatives.
- The cache contains 1,034 retained modes, hence 534,061 unordered mode pairs. A pair is eligible only when its episodes differ and its spatial-mask IoU is at least 0.08; 28,321 pairs pass this shared gate in the saved run.
- The membrane images are used only to render representative spatial examples; the clustering distances are computed from cached mode data.

### Distances and clustering

The canonical arm uses the cached spatial/temporal composite

`D_DTW+space = 0.5 × D_h + 0.5 × D_b`,

where `D_h` is sign-aware DTW distance between mode activation traces and `D_b` is the cached response-field spatial correlation distance. The cached `D_b`, including its canonical sign alignment, is kept unchanged.

The comparison arm replaces only `D_h`:

`D_lag = 1 − max |Pearson r(lag)|`, for integer lags −8…+8,

`D_lag+space = 0.5 × D_lag + 0.5 × D_b`.

At each lag, the two traces are aligned over their overlapping samples. At least `max(3, ceil(0.8 × shorter trace length))` samples must overlap. The maximum absolute correlation defines the distance; signed correlation and best lag are retained for interpretation. Both arms use complete-linkage hierarchical clustering and a strict cut of 0.30 in the saved comparison. A 0.20–0.45 cut sensitivity table is also produced. The same numeric cut is not guaranteed to have equal effective stringency for the two different distance distributions.

### Outputs

`Q2/results_lag_plus_spatial/` contains mode assignments, cluster memberships/intersections, the eligible-pair table, threshold sensitivity, summary diagnostics, and representative spatial overlays. The three example categories are “both methods,” “spatial-DTW only,” and “lag-plus-spatial only.” The review images overlay mode masks on the registered Slice 7 membrane image.

### Reproduction inputs

Provide the Stage-06 Slice 7 cache, `/home/cyf/wbi/wbi_code/src_registration/motion_correlation_pattern_v2.py` (the Q2 runner imports this exact module), and registered Slice 7 membrane frames for the representative overlays. The runner accepts `--cache`, `--membrane-dir`, and `--out`; the output directory must be empty unless `--overwrite-results` is explicitly supplied.

## Q3 — local membrane appearance mismatch around CFUs

### Question and interpretation

Does membrane-image mismatch inside the CFU region tend to be higher at Ca onsets, especially onsets temporally hit by a significant motion-pattern/CFU pair, than at event-free control frames? The measured quantity is an image-appearance residual, **not** a ground-truth registration error: no independent landmarks or true displacement field are used in this diagnostic.

### Pair set and event/control sampling

- The pair list is every cluster055 temporal-co-occurrence pair with global empirical-FDR `q < 0.05` (260 pairs; 55 unique CFUs in the saved run). It is not the smaller 70-edge Q4 spatial-module set.
- In the saved pair table, the used slice IDs are 1–12; the current pattern/CFU input folders are named `Slice01`…`Slice12` and `slice_Z01`…`slice_Z12`. No `Slice00`/`slice_Z00` result is present in this project input set. The scripts use the stored numeric IDs directly; do not add or subtract one when mapping those files.
- For each pair, the temporal analysis supplies the pattern's motion-event peak and best lag. An onset is classified as a pair hit when a CFU onset falls in the three-frame window beginning at `motion peak + best lag` (the existing ±8-frame search and 3-frame onset window are inherited from that temporal test).
- For each CFU, at most 40 Ca-onset frames are sampled with a fixed seed. Pair-hit onsets are included first; remaining slots are filled by non-hit onsets. Event-free controls are sampled from the same temporal quartile, outside ±8 frames of every onset, without reusing a control frame.
- The saved run uses the 1,598-frame movie, DS7 membrane images, CFU masks, and the time-block reference actually assigned to each frame. When reference intervals share a boundary, the interval with the later start is selected.

### Spatial mapping and metrics

Registered membrane TIFF planes are averaged over non-overlapping 7×7 native pixels to the DS7 grid (the trailing native column that does not fit the grid is dropped). The CFU mask is mapped into the registration patch grid using the recorded 5-pixel gap. For each fixed CFU mask and sampled frame:

- `masked MSE = mean((frame − assigned reference)^2)` over CFU-mask pixels.
- `normalized MSE = masked MSE / spatial variance of the same reference pixels`.

Per-CFU onset-minus-control differences are summarized with bootstrap 95% confidence intervals (5,000 resamples) and paired Wilcoxon tests. Pair-hit onset values are also compared descriptively with other onset values. A separate full-frame proxy runner writes the same local metrics for every frame and unique CFU. The Q3 trace figures combine motion RMS activation, CFU dF/F, and the local membrane proxy on a shared time axis; the “error” label should always be read as an appearance-mismatch proxy.

### Reproduction inputs

Provide the cluster055 significant-pair table and its temporal-analysis code, the 12 DS7 AQuA2 CFU MAT files, the registered membrane channel (`membrane/vol_ch1_*.tif`, 1,598 3-D frames), and the saved reference TIFFs (`reference/vol_ref_<start>_<end>.tif`). The sampled runner accepts `--repo-root`, `--cfu-root`, `--registration-root`, `--reference-dir`, `--out`, `--max-onset-frames-per-cfu` (default 40), and `--seed` (default 20261015). The full-frame runner has explicit `--cfu-root` support in the audited version.

## Q4 — pair-specific continuous motion–Ca response (primary model M2)

### Candidate set and occurrence definition

Q4 uses the 70 frozen spatial-module edges in the cluster055 module-to-CFU association table, not every temporal pair. Candidate patterns have at least five members. For each pattern, one occurrence is formed per episode by merging all its modes in that episode. Its interval is the union of member-mode intervals; its peak is the frame with the largest spatial RMS of the summed reconstructed motion contribution. The M2 amplitude predictor is the peak RMS magnitude, divided by the median positive occurrence amplitude. Motion amplitude remains in input units and is not converted to a calibrated physical unit.

The response is continuous AQuA2 CFU dF/F (`cfuInfo1` column 6 in MATLAB's 1-based description; the raw weighted fluorescence field is not used). The model is fitted to the continuous dF/F trace; it is **not** pre-subtracted by an event-specific baseline. For observed event-triggered response curves and summaries, each occurrence's baseline is its mean dF/F over lags −16…−9. The response window and modelled lag range are −8…+8 frames around the motion peak. Negative lag means Ca precedes the motion peak; positive lag means Ca follows it.

### Models

All models include an intercept, linear time drift, and quadratic time drift (`M0`). Motion-related predictors are represented by three L2-normalized Gaussian lag bases centred at −8, 0, and +8 frames, with σ = 4.5 frames. Ridge-penalized linear regression uses λ = 1 on motion basis coefficients; baseline drift terms are not penalized.

- `M0`: baseline drift only.
- `M1`: binary motion-event impulses convolved with the three lag bases.
- `M2` (primary): impulses weighted by each motion occurrence's normalized peak amplitude, convolved with the lag bases.
- `M3`: continuous reconstructed-motion envelope convolved with the same bases; exploratory comparator.

The displayed M2 kernel is the full-record descriptive fitted lag-response coefficient curve. It is not a causal impulse response and should be interpreted alongside the individual observed, motion-peak-aligned Ca traces.

### Validation and inference

- Predictive validation uses five contiguous time-block folds. A 16-frame guard is excluded around held-out blocks to limit leakage across the ±8-frame response window.
- An additional event-group validation merges occurrence windows that overlap within ±8 frames, holds one merged event group out, and excludes its ±8-frame halo from training. This is a stability check, not a second hypothesis test.
- The primary statistic is `CV MSE(M0) − CV MSE(M2)`; positive values favor M2.
- For each pair, the Ca series is circularly shifted 500 times by offsets of at least 17 frames, preserving its internal temporal structure while disrupting alignment to the fixed motion design. The empirical one-sided p-value is `(1 + number of null statistics ≥ observed statistic)/(1 + 500)`.
- Benjamini–Hochberg correction is applied across the 70 tested candidate pairs. The saved 2026-10-26 result has 23/70 pairs at `q ≤ 0.05`; 35/70 have lower mean time-block CV MSE under M2. The mean `MSE(M0) − MSE(M2)` across pairs is −0.0002697, so M2 does not show a general average improvement across all candidates. Thus statistical evidence and the sign of the average CV difference are reported separately. Pair discovery and testing use the same recording, so these are exploratory, selection-conditioned associations—not independent confirmation or causal evidence.

### Current figures and outputs

The curated PNG-only entry point is `Q4/continuous_response_20261026/M2_main_results_png/`:

- `M2_vs_M0_significance_all_pairs.png`: all 70 pairs, their M0−M2 CV difference, and BH-adjusted shift-test q-values.
- `M2_kernel_overview_qsig_pairs.png`: the 23 q-significant M2 kernels on a common color scale.
- `significant_pair_png/<pair>_M2.png`: one image per q-significant pair, with observed motion-aligned occurrences/mean response next to its fitted M2 lag curve.

The tables and code under `tables/` and `code/` are retained for audit and figure regeneration. The main PNG renderer does not refit the models; it reads the saved tables. Existing exploratory PDFs are not part of the curated figure set.

### Reproduction inputs

Provide (i) the relevant pattern `objects.pkl` caches and the pinned `wholistic_registration` source needed to interpret/reconstruct cached motion objects; (ii) the current DS7 AQuA2 CFU MAT files for candidate CFU slices (the saved candidate table uses zero-based CFU slice IDs 1–12); (iii) the frozen module-CFU association CSV defining the 70 candidates; and (iv) the sampled Q3 membrane-metric table used only to attach Q3 context to Q4 pair summaries. The current candidate table uses zero-based motion slice IDs 1–6. Q4's primary M2 fit and circular-shift test do not require raw registered TIFFs once these derived inputs are available. The Q4 runner accepts explicit `--pattern-root`, `--cfu-root`, `--module-table`, `--q3-table`, `--motion-repo`, `--out`, `--shifts`, and `--seed` paths/options. Use the exact paths recorded in `run_config.json`; the script's historical defaults do not all match the recorded run configuration.

## Data checklist for a collaborator

The analysis code/results folder alone is insufficient to rerun all analyses. On the current server, the source paths recorded/used are:

1. **Q2 mode cache:** `/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice07_velocity_decomp/06_patterns/` (`objects.pkl`, `distance_matrix.npz`).
2. **Registered membrane for Q2 overlays and Q3:** `/mnt/data21T_2/cyf/f338/f338_registrated_0530/membrane/` (`vol_ch1_*.tif`, 1,598 3-D frames).
3. **Time-block reference volumes for Q3:** `/mnt/data21T_2/cyf/f338/f338_registrated_0530/reference/` (`vol_ref_<start>_<end>.tif`).
4. **Cluster055 motion project, temporal pair results, and motion caches:** `/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055/`, especially `01_motion_pattern_extraction_omega05_mu05/patterns/` and `04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset/empirical_FDR_significant_pairs.csv`; the temporal runner, `/home/cyf/wbi/wbi_code/src_registration/motion_correlation_pattern_v2.py` for Q2, and the project's `wholistic_registration/` source for Q4 are also required.
5. **Current AQuA2 CFU MAT files:** `/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055/02_current_cfu_input/cfu/` (the saved Q3 significant-pair set uses slice IDs 1–12).

6. **Frozen Q4 module–CFU table:** `/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055/05_local_mechanical_modules_distributed_ca_network/02_module_cfu_network/module_cfu_associations_fdr_q005.csv`.

7. **Q4's auxiliary Q3 metrics table (context only; not a predictor or test input):** `/home/cyf/wbi/wbi_code/experiments/methodological_review_q2_q4/03_q3_local_membrane_mse/results/q3_membrane_qc_reference_matched_mse/q3_membrane_frame_metrics.csv`.

All 12 CFU MAT files in the original and cluster055 project roots were checksum-identical. The Q3 full-frame runner previously relied on the original-root default; it now has an explicit `--cfu-root` option and defaults relative to the selected project root.

The files under the cluster055 CFU input directory are symlinks into `/home/cyf/wbi/wbi_code/experiments/Fig5/23_all_cfu_pattern_lag_cooccurrence_current/input_links/cfu/`. For transfer to a collaborator, copy the target `.mat` files themselves (or recreate the links against the collaborator's data root); copying only the symlinks will not make a portable dataset.

## Command templates (current server)

Run into new, empty output directories; do not point `--out` at the curated source results unless intentionally replacing them. Replace server-specific paths with the collaborator's locations.

```bash
PROJECT=/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055
REVIEW=$PROJECT/06_methodological_review_q2_q4
REG=/mnt/data21T_2/cyf/f338/f338_registrated_0530

# Q2: rerun the Slice 7 clustering comparison
python "$REVIEW/Q2/code/run_slice07_dual_clustering.py" \
  --cache /home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice07_velocity_decomp/06_patterns \
  --membrane-dir "$REG/membrane" \
  --out "$REVIEW/Q2/results_lag_plus_spatial_reproduced"

# Q3 sampled event/control membrane assessment
python "$REVIEW/Q3/code/03_q3_local_membrane_mse/run_q3_membrane_qc.py" \
  --repo-root "$PROJECT" \
  --cfu-root "$PROJECT/02_current_cfu_input/cfu" \
  --registration-root "$REG" --reference-dir "$REG/reference" \
  --out "$REVIEW/Q3/results/q3_membrane_qc_reproduced"

# Q3 all-frame local proxy (can be storage/I/O intensive)
python "$REVIEW/Q3/code/05_q3_significant_pair_traces/run_q3_full_membrane_proxy.py" \
  --repo-root "$PROJECT" \
  --cfu-root "$PROJECT/02_current_cfu_input/cfu" \
  --registration-root "$REG" --reference-dir "$REG/reference" \
  --out "$REVIEW/Q3/results/q3_framewise_local_mse_reproduced"

# Q4 model/table rerun: use all input paths from the saved run_config.json.
# Set --out to a new directory. Then render the curated PNG-only figures:
python "$REVIEW/Q4/continuous_response_20261026/code/run_q4_continuous_response.py" \
  --analysis-root "$PROJECT" \
  --pattern-root "$PROJECT/01_motion_pattern_extraction_omega05_mu05/patterns" \
  --cfu-root "$PROJECT/02_current_cfu_input/cfu" \
  --module-table "$PROJECT/05_local_mechanical_modules_distributed_ca_network/02_module_cfu_network/module_cfu_associations_fdr_q005.csv" \
  --q3-table /home/cyf/wbi/wbi_code/experiments/methodological_review_q2_q4/03_q3_local_membrane_mse/results/q3_membrane_qc_reference_matched_mse/q3_membrane_frame_metrics.csv \
  --motion-repo "$PROJECT/wholistic_registration" \
  --out "$REVIEW/Q4/continuous_response_20261026/reproduced_tables"
python "$REVIEW/Q4/continuous_response_20261026/code/render_m2_main_figures.py" \
  --q4-dir "$REVIEW/Q4/continuous_response_20261026/reproduced_tables" \
  --out "$REVIEW/Q4/continuous_response_20261026/reproduced_M2_main_results_png"
```

The Q4 Q3-metrics path above is auxiliary metadata for pair summaries, not an input to M2 predictors or its circular-shift test. The exact ridge, fold, and shift parameters are recorded in `Q4/continuous_response_20261026/run_config.json`.

Before rerunning, replace server-specific paths with the collaborator's mounted data locations and pass them explicitly to the CLI. See each subfolder README and the saved `run_config.json` for the exact run record.
