# 04：All pattern–CFU lag co-occurrence

## 目的

不使用空间 module 筛选，而是对全部 eligible motion pattern 与全部 CFU 直接检验时间上的 co-occurrence，回答哪些 Ca CFU 与哪些 mechanical pattern 在时间上相互关联。

## 方法

使用 lag=-8…8 帧、window=3 帧。每个 pattern×CFU pair 在各 lag 上计算经验统计量，并使用 500 次全局 circular-shift null 生成 p-value；随后先对每个 pair 的 17 个 lag 做一次校正，再对全部 pair 的 best-lag p-value 做全局校正。显著性阈值为 q<0.05。

## 结果与代码

- `run_all_cfu_pattern_lag8_w3_pairwise.py`：完整计算入口。
- `results/`：逐 pair、逐 lag、校正后的结果。
- `figures/significant_pairs_spatiotemporal_original_overlay_fast/`：显著 pair 的原始尺寸空间/时间 overlay。
- `render_significant_pairs_overlay_fast.py`：批量绘图。
- `render_article_pair.py`：单个 pair 的文章级空间图与 timeline。

当前全局 q<0.05 有 168 个 pair。该实验不使用前一个目录的 ratio、coverage 或 module 结果。
