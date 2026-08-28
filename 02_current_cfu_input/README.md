# 02：Current CFU input

## 目的

集中保存本版本所使用的 AQuA2 CFU 输入，确保空间和时间分析使用同一套 CFU，而不因复制文件产生版本漂移。

## 内容

`cfu/` 下有 12 个 symbolic link，分别指向 Fig5/23 中的 `slice_Z01` 到 `slice_Z12` 的 CFU MAT 文件。CFU 使用 ds7、AQuA2 原生 event detection/aggregation 结果，event 数量阈值为 5。

本目录不包含 AQuA2 重跑代码，也不应在这里修改 CFU。若要改变 CFU 参数，应在 `Fig5/22_aqua2_cfu_pipeline_ds7` 建立新版本并更新输入审计。
