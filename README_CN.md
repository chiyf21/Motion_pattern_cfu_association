# motion_pattern_cfu_association 实验总览（中文）

这是 `omega=0.5, mu=0.5` 的独立比较版本，不覆盖旧 Fig5 的 pattern、CFU、统计结果或图。文件名中的展示编号为 1-based；原始 motion H5 的 `slice` 参数为 0-based。

## 实验链

1. `01_motion_pattern_extraction_omega05_mu05`：从注册后的 motion 数据提取 motion pattern。
2. `02_current_cfu_input`：固定当前 AQuA2 CFU 结果，作为只读输入。
3. `03_cfu_pattern_spatial_overlap`：根据空间位置寻找 pattern–CFU module。
4. `04_all_pattern_cfu_lag_cooccurrence`：对全部 eligible pattern×CFU 进行 lag co-occurrence 检验。
5. `05_local_mechanical_modules_distributed_ca_network`：把空间 module 与时间显著的 CFU 关联连接起来。

CFU 提取不在此目录重新执行。完整 motion pipeline 见实验 01 下的 `pipeline_from_raw_omega05_mu05/`。

## 当前结果

- 12 个 slice，共 2915 个 pattern；成员数≥5 的 pattern 为 288 个。
- 12 个 slice，共 724 个 CFU。
- 空间 module 52 个，准则为 `pattern members≥5, ratio≤3, coverage≥0.5`。
- lag 分析使用 `-8…8` 帧、window=3、500 次 global circular-shift null；全局 q<0.05 有 168 个 pair。
- module 网络有 34 条 q<0.05 的 module–CFU 边，涉及 8 个 module。
