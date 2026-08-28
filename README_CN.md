# Motion–calcium 时空关联分析

本仓库是 `omega=0.5, mu=0.5` 的独立分析版本，用于研究注册后的 motion pattern 与 AQuA2 calcium functional unit（CFU）之间的空间和时间关联。它是一个完整的分析包，不只是单张图的绘图目录。

`01_motion_pattern_extraction_omega05_mu05/` 保存 motion pattern 流程和结果，完整流程从 patch-level motion 开始，依次经过 motion units、episodes、sparse-compact modes，最后得到基于 mode 的 patterns。

`02_current_cfu_input/` 保存所有下游分析共用的固定 CFU 输入。当前链接指向实验室服务器，实际 MAT 文件没有提交到 Git。

`03_cfu_pattern_spatial_overlap/` 检验 pattern mask 与 CFU mask 的空间对应关系，并根据 pattern 成员数、ratio 和 coverage 定义候选 local mechanical module。

`04_all_pattern_cfu_lag_cooccurrence/` 独立检验全部 pattern×CFU 的时间关联，不使用空间 module 筛选。使用 lag=-8 到 +8 帧、window=3 和经验全局 FDR。

`05_local_mechanical_modules_distributed_ca_network/` 将 03 的空间 module 与 04 的时间显著 CFU 关联起来，分析每个局部 mechanical module 可能对应的分布式 Ca 信号。

当前结果包括 12 个 slice 的 2915 个 pattern，其中成员数至少为 5 的有 288 个；CFU 共 724 个；空间 module 共 52 个；lag 分析全局 q<0.05 的 pair 为 168 个；module 网络有 34 条 q<0.05 的 module–CFU 边，涉及 8 个 module。

仓库不包含大型缓存、原始 motion H5、pickle、NPZ、MAT 和大型 lag 结果。完整 motion 原始数据约 801 GB，应通过服务器挂载或独立数据包提供。推荐顺序为 01 → 03/04 → 05，具体输入、输出、参数和命令见各目录的中英文 README。

## 提供给合作者的文件

本项目的交付起点是 patch-level motion，而不是原始 motion H5。每个展示编号的 slice 应提供 `01_motion_pattern_extraction_omega05_mu05/patch_motion/SliceXX_velocity_decomp/01_patch_motion/arrays.npz` 及对应的 `metadata.json`。NPZ 必须包含 `motion_delta`、`motion_abs` 和 `mask_patched`；这些是 patch 网格数据，不是原始分辨率视频。历史约定中 `mask_patched=True` 表示背景/无效 patch，后续分析使用其补集作为有效区域。

如果合作者只分析已有 pattern，还需要提供 `patterns/SliceXX_velocity_decomp/06_patterns/objects.pkl`、存在时的 `distance_matrix.npz` 和对应 `metadata.json`。如果要从 patch motion 重新生成 pattern，则还应提供 `02_motion_units/`、`03_episodes/`、`04_modes/`，或者在提供兼容 Stage 1 缓存后运行 `pipeline_from_raw_omega05_mu05/` 中的脚本。pattern objects 体积较大，故被 Git 忽略。

CFU 交付内容是 `02_current_cfu_input/cfu/` 中链接所代表的 12 个实际 MAT 文件。只提供 symbolic link 不够，因为它们目前指向服务器外部路径。这些 CFU 来自 ds7 的 AQuA2 原生结果，每个 CFU 至少包含 5 个 event。

空间分析读取 pattern objects 和 CFU MAT，输出 `03_cfu_pattern_spatial_overlap/spatial_all_overlap_slice*.csv`、`spatial_final_slice*.csv` 与 `summary_all_slices.csv`。时间分析读取同一批 pattern 和 CFU，输出位于 `04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset/`。module network 分析读取空间 final 表和时间显著 pair 表，输出位于 `05_local_mechanical_modules_distributed_ca_network/`。

简而言之：patch motion 是 motion 输入，CFU MAT 是 calcium 输入；03 和 04 是并列分析；05 同时依赖二者。只要这些派生输入齐全，就不需要提供原始 motion H5、原始 calcium movie，也不需要重新运行 AQuA2。
