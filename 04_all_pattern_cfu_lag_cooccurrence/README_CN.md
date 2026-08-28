# 04 全部 pattern–CFU lag co-occurrence

这是与空间 module 并列的时间分析。它不使用实验 03 的空间筛选，而是对所有 eligible motion pattern×CFU pair 进行检验。每个 pair 检验 -8 到 +8 帧共 17 个 lag，window=3 帧；null 使用 500 次对 CFU timeline 一致施加的 global circular shift。先处理 pair 内的 lag 搜索，再对 best-lag 统计量进行全局 empirical FDR 校正。

运行入口是 `run_all_cfu_pattern_lag8_w3_pairwise.py`。逐 lag 和逐 pair 表位于 `results/`；`render_significant_pairs_overlay_fast.py` 批量绘制显著 pair；`render_article_pair.py` 绘制单个文章级 pair。当前全局 q<0.05 的 pair 为 168 个。本实验不使用实验 03 的 ratio、coverage 或 module membership。
