# 完整 motion→pattern pipeline

`pipeline_from_raw_omega05_mu05/run_motion_to_patterns.py` 是从 patch-level motion 开始的完整可恢复流程。它不是从原始 801 GB motion H5 开始；原始 motion 需要先由外部数据流程转换为 patch motion，或由脚本的 raw-H5 输入接口生成。

每个 slice 使用 0-based 索引。运行时需要提供 motion 文件目录、invalid/background mask、slice 编号和输出根目录。输出阶段为 `01_patch_motion`、`02_motion_units`、`03_episodes`、`04_modes` 和 `06_patterns`。默认参数与当前 canonical pattern 版本一致，具体参数写入每个 cache 的 metadata.json。
