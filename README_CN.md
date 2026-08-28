# Motion–calcium 时空关联分析

本仓库是 `omega=0.5, mu=0.5` 的独立分析版本，用于研究注册后的 motion pattern 与 AQuA2 calcium functional unit（CFU）之间的空间和时间关联。它是一个完整的分析包，不只是单张图的绘图目录。

`01_motion_pattern_extraction_omega05_mu05/` 保存 motion pattern 流程和结果，完整流程从 patch-level motion 开始，依次经过 motion units、episodes、sparse-compact modes，最后得到基于 mode 的 patterns。

`02_current_cfu_input/` 保存所有下游分析共用的固定 CFU 输入。当前链接指向实验室服务器，实际 MAT 文件没有提交到 Git。

`03_cfu_pattern_spatial_overlap/` 检验 pattern mask 与 CFU mask 的空间对应关系，并根据 pattern 成员数、ratio 和 coverage 定义候选 local mechanical module。

`04_all_pattern_cfu_lag_cooccurrence/` 独立检验全部 pattern×CFU 的时间关联，不使用空间 module 筛选。使用 lag=-8 到 +8 帧、window=3 和经验全局 FDR。

`05_local_mechanical_modules_distributed_ca_network/` 将 03 的空间 module 与 04 的时间显著 CFU 关联起来，分析每个局部 mechanical module 可能对应的分布式 Ca 信号。

当前结果包括 12 个 slice 的 2915 个 pattern，其中成员数至少为 5 的有 288 个；CFU 共 724 个；空间 module 共 52 个；lag 分析全局 q<0.05 的 pair 为 168 个；module 网络有 34 条 q<0.05 的 module–CFU 边，涉及 8 个 module。

仓库不包含大型缓存、原始 motion H5、pickle、NPZ、MAT 和大型 lag 结果。完整 motion 原始数据约 801 GB，应通过服务器挂载或独立数据包提供。推荐顺序为 01 → 03/04 → 05，具体输入、输出、参数和命令见各目录的中英文 README。
