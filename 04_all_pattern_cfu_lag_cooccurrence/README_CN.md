# 04：All pattern–CFU lag co-occurrence（中文）

本实验不使用空间 module 筛选，而是对所有 eligible motion pattern 与所有 CFU 直接检验时间 co-occurrence。lag 为 -8…8 帧，window=3；每个 pair 使用 500 次全局 circular-shift null 计算 p-value。先对同一 pair 的 17 个 lag 校正，再对全部 pair 的 best-lag p-value 做全局校正。

代码入口是 `run_all_cfu_pattern_lag8_w3_pairwise.py`；`results/` 保存逐 lag 和逐 pair 结果；`render_significant_pairs_overlay_fast.py` 批量生成显著 pair 图；`render_article_pair.py` 生成单个文章级 pair 图。当前全局 q<0.05 有 168 个 pair，图在 `figures/`。
