# 01：Motion pattern extraction

## 目的

从注册后的原始 motion 位移场开始，提取 patch-level motion、motion units、motion episodes、sparse-compact motion modes，最后将 mode 聚类为 motion patterns。这里使用 `omega=0.5, mu=0.5`，因此 pattern 距离同时考虑空间与运动特征。

## 当前结果

当前 v0827 目录保存了 12 个 slice 的最终 `06_patterns/objects.pkl`。总 pattern 数为 2915，成员数≥5 的 pattern 为 288 个。

## 代码

- `rerun_patterns_omega05_mu05.py`：从已有 `04_modes` 缓存开始，只重新执行最终 mode→pattern 聚类。
- `pipeline_from_raw_omega05_mu05/run_motion_to_patterns.py`：完整 raw motion→pattern 流程。
- `01_patch_motion.sh`、`02_motion_units.sh`、`03_episodes.sh`、`04_modes_to_patterns.sh`：可分阶段运行的入口。

完整流程使用 Fig5 下 clone 的 `wholistic_registration`，commit 为 `2b3c4e611ca194d391a31b56f6bc93a28ca90b13`。中间结果按 `01_patch_motion`、`02_motion_units`、`03_episodes`、`04_modes`、`06_patterns` 保存，并在 metadata 中记录参数。

## 关键参数

patch size=7；velocity 使用 frame-to-frame motion difference；mode 使用 `Kmax=8, svd_target_r2=0.90, lambda_sc=0.05, rho=1, kappa=4`；pattern 使用 complete linkage、`cluster_dist_thresh=0.45`、`min_iou=0.08`、best-connected-component unified mask。
