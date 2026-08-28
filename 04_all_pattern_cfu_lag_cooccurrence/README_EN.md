# 04: All pattern–CFU lag co-occurrence

This experiment does not use spatial-module filtering. It directly tests temporal co-occurrence between every eligible motion pattern and every CFU. Lags are -8…8 frames with a three-frame window; p-values are computed using 500 global circular-shift nulls. The 17 lag tests are corrected within each pair first, followed by global correction across all pairs using the best-lag p-values.

Run `run_all_cfu_pattern_lag8_w3_pairwise.py`. The `results/` directory stores lag-level and pair-level outputs. `render_significant_pairs_overlay_fast.py` renders significant pairs in batch, while `render_article_pair.py` renders one article-style pair. There are currently 168 pairs with global q<0.05; figures are under `figures/`.
