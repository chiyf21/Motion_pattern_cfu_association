# 02：Current CFU input（中文）

此目录集中保存本版本使用的 AQuA2 CFU 输入，避免不同分析读取不同版本。`cfu/` 下有 12 个 symbolic link，分别指向 Fig5/23 中 slice01–slice12 的 CFU MAT 文件。

这些 CFU 来自 ds7 的 AQuA2 原生 event detection/aggregation，event 数量阈值为 5。本目录不重新运行 AQuA2，也不应直接修改链接目标；如需改变 CFU 参数，应在 Fig5/22 中建立新版本并更新输入审计。
