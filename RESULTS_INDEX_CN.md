# 结果索引

本仓库将可复现代码与大型输入数据、缓存结果分开管理。

## 需要单独提供的输入

- 每个 0-based slice 的 patch-motion 数组，`arrays.npz` 中应包含 `motion_delta`、`motion_abs` 和 `mask_patched`。
- `02_current_cfu_input/cfu/` 中链接所指向的 12 个当前 AQuA2 CFU MAT 文件。符号链接依赖服务器，不能单独交付，需要提供实际 MAT 文件。
- 只有在绘制原始尺寸图时才需要 reference TIFF。
- 如果希望直接复用当前 pattern 结果而不重新计算，还需要单独提供 `06_patterns/objects.pkl`。

## 流程和结果

1. `01_motion_pattern_extraction_omega05_mu05/` 提取 motion pattern；`members>=5` 的 overview 位于 `patterns/_overview_members_ge5/`，这些图已纳入 Git。
2. `03_cfu_pattern_spatial_overlap/` 计算 pattern–CFU 空间 module；12 个原始尺寸 module overview 位于 `figures/original_resolution_pattern_cfu_overlap/`，这些图已纳入 Git。
3. `04_all_pattern_cfu_lag_cooccurrence/` 对所有符合条件的 pattern×CFU 在 -8…+8 帧、window=3 下进行检验。大型结果表和 overlay 不纳入 Git，可用已有绘图脚本重新生成。
4. `05_local_mechanical_modules_distributed_ca_network/` 合并空间 module 与时间关联结果。

修改 `config.py` 并准备外部输入后，可运行 `bash run_all_pipeline.sh` 完整执行分析。本仓库不重新运行 CFU 提取。

常用路径和参数集中写在 `config.py` 中；需要局部调试时，各阶段脚本也提供相应命令行参数。
