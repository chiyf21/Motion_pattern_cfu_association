# 05：Local mechanical module embedded in distributed Ca network

## 目的

以 03 中的 spatial module 为起点，连接 04 中时间上显著的 CFU，观察一个局部 mechanical module 是否对应多个分布式 Ca 上游信号。这是 module-centric 的后续网络分析。

## 输入

- 03 的 `spatial_final_slice*.csv`：module 定义。
- 04 的全局显著 pair 表：CFU–pattern 时间关系。
- 01 的 pattern masks 与 activation：module 空间和 pattern 身份。
- 02 的当前 CFU 文件：CFU 空间位置。

## 结果与代码

- `run_module_network.py`：构建 module–CFU 网络并生成汇总表。
- `render_module_cfu_spatial_gallery.py`：绘制 module 及其 q<0.05 CFU 的空间 gallery。
- `00_input_version_audit.md/csv`：记录各输入文件的版本来源。
- `02_module_cfu_network/`：网络边表、q<0.05/q<0.001 筛选结果。
- `03_module_cfu_spatial_gallery_q005/`：module-centric 空间图。

当前 q<0.05 module–CFU 边为 34 条，涉及 8 个 module；当前没有额外把正 lag/负 lag解释为不同方向，方向解释留给后续分析。
