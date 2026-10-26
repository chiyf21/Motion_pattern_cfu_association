# cluster055 补充分析：Q2–Q4

本目录包含三个彼此独立的分析。英文详细方法、参数、复现命令及数据清单以 [`Q2_Q4_METHODS_EN.md`](Q2_Q4_METHODS_EN.md) 为准；最外层 `README.md` 是面向合作者的英文入口。

- **Q2**：比较 Slice 7 已保留 motion modes 的“空间距离 + DTW”和“相同空间距离 + lag correlation”两种聚类。当前结果在 `Q2/results_lag_plus_spatial/`。
- **Q3**：使用全部 q<0.05 的时间显著 pattern–CFU pairs（当前 260 对、55 个 CFU），检查 pair-hit onset、其他 Ca onset 和 event-free control 附近的局部 membrane mismatch。这里的 MSE/RMSE 是图像外观差异 proxy，不是真值配准误差。结果和逐 pair 三曲线图位于 `Q3/`。
- **Q4**：对冻结的 70 条空间 module–CFU 候选关系拟合连续 Motion–Ca 响应，主要关注 M2 相对 M0 的检验和 M2 lag-response。当前主图只从 `Q4/continuous_response_20261026/M2_main_results_png/` 查看：两张总览 PNG 和每个 q 显著 pair 一张 PNG。

Q4 中已按本次约定删除旧 `pair_figures/` 和 `M2_stability_20261026/`；表格、runner、参数记录和当前 PNG 保留。Q4 `response_metrics/` 是早先单独计算的响应幅度/延迟/一致性指标，并非这版 M2 主图的输入；`synthetic_tests/` 当前为空。

注意：项目中的 CFU MAT 路径是符号链接，移交数据时应复制链接目标文件，而不是只复制链接。Q2/Q3/Q4 所需数据位置、slice ID 和脚本参数见英文详细说明。完整的模型重跑没有在本次整理中执行；本次做了脚本检查与现存结果验证。
