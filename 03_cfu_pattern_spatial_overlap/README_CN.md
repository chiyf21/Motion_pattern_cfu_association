# 03 CFU–pattern 空间 overlap

这是独立于时间显著性的空间分析。每个 slice 中，将 eligible pattern mask 与 CFU mask 两两比较，全部指标写入 `spatial_all_overlap_slice*.csv`；最终 module pair 写入 `spatial_final_slice*.csv`；`summary_all_slices.csv` 为跨 slice 汇总。

筛选条件是 pattern members≥5、`ratio≤3`、`coverage≥0.5`。`run_spatial_overlap.py` 重新计算表格；`render_all_slice_original_overlap.py` 绘制原始尺寸空间 overview；`render_spatial_overlap_stats.py` 绘制 Fig4 风格汇总；`render_spatial_overlap_couple_overlays.py` 绘制 module 的 motion-arrow overlay。当前共 52 个空间 module，不使用实验 04 的 p-value 预筛选。
