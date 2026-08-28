# 03：CFU–pattern spatial overlap（中文）

本实验在每个 slice 内比较 motion pattern mask 与 CFU mask。空间上高度对应的 pair 定义为候选 local mechanical module；本实验与 lag co-occurrence 并列，不使用时间显著性筛选。

仅保留 pattern members≥5，并使用 `ratio≤3` 与 `coverage≥0.5`。`spatial_all_overlap_slice*.csv` 保存全部候选及指标，`spatial_final_slice*.csv` 保存最终 module，`summary_all_slices.csv` 为汇总。`run_spatial_overlap.py` 负责计算；`render_all_slice_original_overlap.py`、`render_spatial_overlap_stats.py` 和 `render_spatial_overlap_couple_overlays.py` 负责绘图。当前共 52 个 module，图位于 `figures/`。
