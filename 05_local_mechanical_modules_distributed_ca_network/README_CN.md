# 05 局部 mechanical module 嵌入分布式 Ca network

本实验从 03 的空间 module 出发，使用 04 的时间关联结果，分析每个 module 中的 pattern 与哪些 CFU 在时间上相关。它把空间身份和时间关联连接起来，但不会重新运行 pattern 提取、CFU 提取或 lag 检验。

`run_module_network.py` 连接 module 表、pattern objects、CFU 输入和显著 pair 表；`render_module_cfu_spatial_gallery.py` 绘制每个 module 及其 q<0.05 CFU 的空间位置；`00_input_version_audit.md/csv` 记录精确输入来源。网络表位于 `02_module_cfu_network/`，图位于 `03_module_cfu_spatial_gallery_q005/`。当前有 34 条 q<0.05 module–CFU 边，涉及 8 个 module。
