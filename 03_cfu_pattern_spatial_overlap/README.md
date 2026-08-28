# 03：CFU–pattern spatial overlap

## 目的

在每个 slice 内比较 motion pattern 的空间 mask 与 CFU 的空间 mask。空间上高度对应的 pair 被定义为候选 local mechanical module；这个实验与 lag co-occurrence 并列，不使用时间显著性筛选。

## 判据

只保留 pattern members≥5 的 pattern，并使用 `ratio≤3` 与 `coverage≥0.5`。其中 ratio 表示相对面积关系，coverage 表示 pattern/CFU 的空间覆盖程度，具体字段保存在逐 slice CSV 中。

## 结果与代码

- `spatial_all_overlap_slice*.csv`：全部候选 pair 及其空间指标。
- `spatial_final_slice*.csv`：最终选中的 module pair。
- `summary_all_slices.csv`：12 个 slice 汇总。
- `run_spatial_overlap.py`：重新计算空间 overlap。
- `render_all_slice_original_overlap.py`：原始尺寸空间 overview。
- `render_spatial_overlap_stats.py`：Fig4-style 全 slice 统计图。
- `render_spatial_overlap_couple_overlays.py`：module–CFU motion arrow overlay。

当前共得到 52 个 spatial module。图和 overlay 位于 `figures/`。
