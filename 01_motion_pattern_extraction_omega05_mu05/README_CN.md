# 01 Motion pattern 提取

这是唯一产生 motion pattern 的阶段。`rerun_patterns_omega05_mu05.py` 是从已有 mode 缓存开始的快速重聚类脚本；如果需要从 motion 数据重新开始，应使用 `pipeline_from_raw_omega05_mu05/run_motion_to_patterns.py`。它会为一个 0-based slice 保存可恢复的 `01_patch_motion`、`02_motion_units`、`03_episodes`、`04_modes` 和 `06_patterns`。

当前版本使用 patch size=7，相邻帧差分作为 velocity，用 local-MAD（时间窗 21、空间窗 3）检测 motion unit，并进行 episode artifact filtering。mode 参数为 `Kmax=8`、SVD target R²=0.90、`lambda_sc=.05, rho=1, kappa=4`；最终直接对 mode 聚类，使用 complete linkage、`min_iou=.08`、`cluster_dist_thresh=.45`、`omega=mu=.5` 和 best-connected-component unified mask。实现来自 `../wholistic_registration/` submodule。

结果位于 `patterns/SliceXX_velocity_decomp/06_patterns/objects.pkl`。12 个 slice 共 2915 个 pattern，其中成员数至少为 5 的有 288 个。大型对象被 Git 忽略，如需复用当前结果需要单独提供。
