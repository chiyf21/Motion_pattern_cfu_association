# 04 All pattern–CFU lag co-occurrence

This is the parallel temporal analysis. It tests all eligible motion-pattern × CFU pairs without using spatial-module selection. For every pair it evaluates 17 lags (-8…+8 frames) with a three-frame window. The null distribution uses 500 global circular shifts applied consistently to CFU timelines. The implementation first handles the lag search within a pair, then applies the global pair-level empirical FDR procedure to the best-lag statistics.

Run `run_all_cfu_pattern_lag8_w3_pairwise.py`. Pair-level and lag-level tables are written under `results/`; significant-pair overlays are rendered by `render_significant_pairs_overlay_fast.py`; `render_article_pair.py` creates a publication-style figure for one requested pair. The current run has 168 pairs with global q<0.05. It does not use ratio, coverage, or module membership from experiment 03.
