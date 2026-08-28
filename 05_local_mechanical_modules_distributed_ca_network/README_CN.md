# 05：Local mechanical module embedded in distributed Ca network（中文）

本实验以 03 中的空间 module 为起点，连接 04 中时间上显著的 CFU，观察一个局部 mechanical module 是否对应多个分布式 Ca 上游信号。

输入包括 03 的 spatial module 表、04 的显著 pair 表、01 的 pattern mask/activation，以及 02 的 CFU 文件。`run_module_network.py` 构建网络，`render_module_cfu_spatial_gallery.py` 绘制空间 gallery，`00_input_version_audit.md/csv` 记录输入版本。当前有 34 条 q<0.05 module–CFU 边，涉及 8 个 module。
