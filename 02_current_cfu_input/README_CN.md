# 02 当前 CFU 输入

本目录固定保存实验 03–05 共用的 AQuA2 CFU 文件，按展示编号 slice01–slice12 各有一个链接。当前链接指向仓库外部的实验室服务器，因此不是可移植路径。

这些输入来自 ds7 的 AQuA2 原生 event detection 和 CFU 聚合结果，每个 CFU 至少包含 5 个 event。本目录不运行 AQuA2。若在其他机器复现，应将实际 MAT 文件复制到本地数据目录，再修改链接或本地路径配置，不要把私人数据提交到 GitHub。

可复现的源代码整理在 `aqua2_native_pipeline/` 中，包含 calcium 输入准备、原生 AQuA2 DS7 event detection、CFU 聚合以及明确的全 slice 批处理脚本。批处理脚本会复用已有 event 结果；如果对应 CFU 文件已存在则跳过，设置 `FORCE=1` 才会强制重新提取。
