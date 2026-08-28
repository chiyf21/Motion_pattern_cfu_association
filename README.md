# motion_pattern_cfu_association 实验总览

这是 omega=0.5、mu=0.5 的独立比较版本。它不覆盖旧 Fig5 的 pattern、CFU、统计结果或图。所有 slice 编号在文件名中使用 1-based 展示编号；原始 motion H5 的 slice 参数仍明确是 0-based。

## 实验链

1. `01_motion_pattern_extraction_omega05_mu05`：从 motion 数据得到 motion pattern。
2. `02_current_cfu_input`：固定当前 AQuA2 CFU 结果，作为只读输入。
3. `03_cfu_pattern_spatial_overlap`：寻找空间上对应的 pattern–CFU module。
4. `04_all_pattern_cfu_lag_cooccurrence`：对全部 eligible pattern×CFU 做 lag co-occurrence。
5. `05_local_mechanical_modules_distributed_ca_network`：将空间 module 与时间显著 CFU 关联连接起来。

CFU 提取不在本目录重新运行；本版本使用 `02_current_cfu_input/cfu/` 中的链接文件。完整 motion pipeline 的运行入口见实验 01 的 `pipeline_from_raw_omega05_mu05/README.md`。

## 当前结果摘要

- pattern：12 个 slice，共 2915 个 pattern；成员数≥5 的 pattern 为 288 个。
- CFU：12 个 slice，共 724 个 CFU 输入。
- 空间 module：52 个通过 `ratio≤3, coverage≥0.5, pattern members≥5` 的 pattern–CFU 对。
- 时间分析：lag=-8…8、window=3、500 次全局循环移位 null；全局校正后 q<0.05 有 168 个 pair。
- module 网络：q<0.05 的 module–CFU 边 34 条，涉及 8 个 module。

推荐运行顺序：先确认 01 的所有 pattern 文件，再运行 03 和 04，最后运行 05。`run_all_pipeline.sh` 是批量入口。
