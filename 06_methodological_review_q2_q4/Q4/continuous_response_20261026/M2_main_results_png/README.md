# Focused Q4 M2 figures

This folder contains PNGs only. It answers two questions for the 70 frozen Q4 candidate pairs: (1) whether M2 improves over M0 under the existing circular-shift test with BH correction; (2) what lag-response relationship M2 fits for the q≤0.05 pairs.

- `M2_vs_M0_significance_all_pairs.png`: all candidate pairs; vertical 0 means no mean time-block CV improvement, and horizontal line is q=0.05. Green points are significant.
- `M2_kernel_overview_qsig_pairs.png`: M2 kernels for the 23 q-significant pairs on a shared color scale.
- `significant_pair_png/`: one PNG per q-significant pair, showing observed motion-peak-aligned Ca occurrences and the M2 fitted kernel.

The event-response curves are baseline-corrected observed dF/F. The kernel is a full-record descriptive fit, not a fold-specific estimate; it describes predictive association, not causation. Existing q values come from the 500-shift test and BH correction over the candidate set. No model was refit here.
