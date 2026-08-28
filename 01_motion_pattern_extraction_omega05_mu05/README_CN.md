# 01：Motion pattern extraction（中文）

## 目的

从注册后的原始 motion 位移场开始，依次生成 patch motion、motion units、motion episodes、sparse-compact motion modes，最后将 mode 聚类为 motion patterns。`omega=0.5, mu=0.5` 使聚类距离同时考虑空间和运动特征。

## 结果与代码

当前保存了 12 个 slice 的 `06_patterns/objects.pkl`，共 2915 个 pattern，其中成员数≥5 的有 288 个。`rerun_patterns_omega05_mu05.py` 只从已有 `04_modes` 开始；`pipeline_from_raw_omega05_mu05/` 中的 `run_motion_to_patterns.py` 才是从 raw motion 到 pattern 的完整流程，并按阶段保存缓存。

关键参数：patch size=7；velocity 为相邻帧 motion difference；mode 为 `Kmax=8, svd_target_r2=0.90, lambda_sc=0.05, rho=1, kappa=4`；pattern 使用 complete linkage、`cluster_dist_thresh=0.45`、`min_iou=0.08` 和 best-connected-component unified mask。
