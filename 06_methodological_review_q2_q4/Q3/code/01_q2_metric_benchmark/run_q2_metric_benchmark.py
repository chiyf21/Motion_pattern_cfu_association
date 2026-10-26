#!/usr/bin/env python3
"""Compare the current sign-aware DTW with lagged correlation on synthetic and real modes.

The real-data arm is a sensitivity audit, not an independent validation: the
existing pattern labels were themselves produced using DTW.
"""
from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

EXP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXP_ROOT))

from review_common import PATTERN_ROOT, REPO_ROOT, load_patterns, write_csv, write_json


def dtw(x, y, radius=None):
    x = np.asarray(x, dtype=np.float32).reshape(-1)
    y = np.asarray(y, dtype=np.float32).reshape(-1)
    n, m = len(x), len(y)
    if not n or not m:
        return np.inf
    constrained = radius is not None
    radius = max(n, m) if radius is None else max(1, int(radius))
    dp = np.full((n + 1, m + 1), np.inf, dtype=np.float32)
    dp[0, 0] = 0.0
    for i in range(1, n + 1):
        if constrained:
            center = i * m / n
            lo, hi = max(1, int(np.ceil(center - radius))), min(m, int(np.floor(center + radius)))
        else:
            lo, hi = 1, m
        for j in range(lo, hi + 1):
            dp[i, j] = abs(x[i - 1] - y[j - 1]) + min(
                dp[i - 1, j], dp[i, j - 1], dp[i - 1, j - 1]
            )
    return float(dp[n, m] / (n + m + 1e-8))


def sign_dtw(x, y, radius=None):
    return min(dtw(x, y, radius), dtw(x, -np.asarray(y), radius))


def lagged_corr(x, y, max_lag=8, min_overlap_fraction=0.8):
    """Maximum absolute Pearson correlation over bounded integer shifts."""
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    best_r, best_lag = 0.0, 0
    min_overlap = max(3, int(np.ceil(min(len(x), len(y)) * min_overlap_fraction)))
    for lag in range(-max_lag, max_lag + 1):
        xs, ys = (lag, 0) if lag >= 0 else (0, -lag)
        n = min(len(x) - xs, len(y) - ys)
        if n < min_overlap:
            continue
        a, b = x[xs:xs + n], y[ys:ys + n]
        a, b = a - a.mean(), b - b.mean()
        den = np.linalg.norm(a) * np.linalg.norm(b)
        r = float(np.dot(a, b) / den) if den > 1e-12 else 0.0
        if abs(r) > abs(best_r):
            best_r, best_lag = r, lag
    return float(np.clip(abs(best_r), 0.0, 1.0)), int(best_lag)


def gaussian_trace(n, centers, widths, amplitudes):
    t = np.arange(n, dtype=float)
    y = np.zeros(n, dtype=float)
    for c, w, a in zip(centers, widths, amplitudes):
        y += a * np.exp(-0.5 * ((t - c) / w) ** 2)
    return y


def resample_trace(x, n):
    return np.interp(np.linspace(0, len(x) - 1, n), np.arange(len(x)), x)


def synthetic_cases():
    base = gaussian_trace(120, [22, 54, 91], [4, 7, 5], [1.0, 0.65, 0.9])
    shifted = np.zeros_like(base)
    shifted[7:] = base[:-7]
    slow = resample_trace(base, 165)
    # Nonuniform time warp: early part is compressed and late part expanded.
    source_t = np.arange(len(base), dtype=float)
    warped_t = np.interp(np.linspace(0, len(base) - 1, 120), [0, 30, 60, 119], [0, 55, 70, 119])
    piecewise = np.interp(warped_t, source_t, base)
    different = gaussian_trace(120, [22, 55, 114], [4, 7, 5], [1.0, 0.65, 0.9])
    return [
        ("translation_only", base, shifted),
        ("uniform_speed_change", base, slow),
        ("piecewise_speed_change", base, piecewise),
        ("similar_peaks_different_global_timing", base, different),
    ]


def metric_row(case, x, y):
    corr, lag = lagged_corr(x, y)
    radius = max(1, int(round(0.10 * max(len(x), len(y)))))
    return {
        "case": case,
        "length_x": len(x),
        "length_y": len(y),
        "unconstrained_sign_dtw": sign_dtw(x, y),
        "constrained_sign_dtw_10pct_band": sign_dtw(x, y, radius),
        "max_abs_lagged_pearson_r_lag_pm8": corr,
        "best_lag_frames": lag,
        "lagged_correlation_distance": 1.0 - corr,
    }


def region_mask(region):
    m = getattr(region, "region_mask", None)
    if m is None:
        m = getattr(region, "support_mask", None)
    if m is None:
        a = getattr(region, "response_strength", None)
        if a is None:
            return None
        m = np.asarray(a) > 0
    m = np.asarray(m, dtype=bool)
    return m if m.ndim == 2 and np.any(m) else None


def sample_pattern_members(patterns, max_regions, seed):
    rng = np.random.default_rng(seed)
    groups = []
    for p in patterns:
        if int(getattr(p, "n_members", 0)) < 5:
            continue
        members = getattr(p, "components", None) or getattr(p, "regions", [])
        members = [m for m in members if region_mask(m) is not None and getattr(m, "activation", None) is not None]
        if len(members) >= 2:
            groups.append((int(p.pattern_id), members))
    rng.shuffle(groups)
    chosen = []
    for pid, members in groups:
        if len(chosen) >= max_regions:
            break
        k = min(2, len(members), max_regions - len(chosen))
        for idx in rng.choice(len(members), size=k, replace=False):
            chosen.append((pid, members[int(idx)]))
    if len(chosen) < max_regions:
        have = {(pid, id(r)) for pid, r in chosen}
        rest = [(pid, r) for pid, members in groups for r in members if (pid, id(r)) not in have]
        rng.shuffle(rest)
        chosen.extend(rest[: max_regions - len(chosen)])
    records = []
    for pid, r in chosen:
        h = np.asarray(r.activation, dtype=np.float32).reshape(-1)
        if h.size < 3:
            continue
        mask = region_mask(r)
        indices = np.flatnonzero(mask.reshape(-1)).astype(np.int32)
        if not indices.size:
            continue
        rf = getattr(r, "response_field", None)
        rs = getattr(r, "response_strength", None)
        records.append({
            "pattern_id": pid,
            "episode_id": int(getattr(r, "episode_id", -1)),
            "activation": h,
            "mask_shape": mask.shape,
            "mask_indices": indices,
            "response_field": None if rf is None else np.asarray(rf, dtype=np.float32),
            "response_strength": None if rs is None else np.asarray(rs, dtype=np.float32),
        })
    return records


def sparse_iou(a, b):
    if a["mask_shape"] != b["mask_shape"]:
        return 0.0
    inter = np.intersect1d(a["mask_indices"], b["mask_indices"], assume_unique=True).size
    return float(inter / (len(a["mask_indices"]) + len(b["mask_indices"]) - inter + 1e-12))


def response_distance(a, b, iou):
    # Match the active pipeline's sign-insensitive correlation distance on overlap.
    if iou <= 0 or a["response_field"] is None or b["response_field"] is None:
        return np.nan
    m = np.zeros(a["mask_shape"][0] * a["mask_shape"][1], dtype=bool)
    common = np.intersect1d(a["mask_indices"], b["mask_indices"], assume_unique=True)
    if common.size < 3:
        return np.nan
    m[common] = True
    m = m.reshape(a["mask_shape"])
    x = a["response_field"][m].reshape(-1).astype(float)
    y = b["response_field"][m].reshape(-1).astype(float)
    if x.size != y.size or x.size < 6:
        return np.nan
    x -= x.mean(); y -= y.mean()
    den = np.linalg.norm(x) * np.linalg.norm(y)
    r = float(np.dot(x, y) / den) if den > 1e-12 else 0.0
    return 1.0 - abs(float(np.clip(r, -1, 1)))


def lower_distance_auc(values, same):
    values = np.asarray(values, dtype=float)
    same = np.asarray(same, dtype=bool)
    p, n = values[same], values[~same]
    if not len(p) or not len(n):
        return np.nan
    # Mann-Whitney U for the probability that a same-pattern pair has smaller distance.
    u = mannwhitneyu(p, n, alternative="two-sided", method="asymptotic").statistic
    return float(1.0 - u / (len(p) * len(n)))


def run_real(out_dir, repo_root, pattern_root, slices, max_regions, seed):
    pair_rows, summary_rows = [], []
    for z in slices:
        payload = load_patterns(z, pattern_root, repo_root)
        records = sample_pattern_members(payload["patterns"], max_regions, seed + z)
        del payload
        gc.collect()
        same_flags, d_dtw, d_cdtw, d_xcorr, d_resp = [], [], [], [], []
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                a, b = records[i], records[j]
                if a["episode_id"] == b["episode_id"]:
                    continue
                iou = sparse_iou(a, b)
                if iou < 0.08:
                    continue
                xcorr, xlag = lagged_corr(a["activation"], b["activation"], max_lag=8)
                radius = max(1, int(round(0.10 * max(len(a["activation"]), len(b["activation"])))))
                dd = sign_dtw(a["activation"], b["activation"])
                dc = sign_dtw(a["activation"], b["activation"], radius)
                dr = response_distance(a, b, iou)
                same = a["pattern_id"] == b["pattern_id"]
                same_flags.append(same); d_dtw.append(dd); d_cdtw.append(dc); d_xcorr.append(1.0 - xcorr); d_resp.append(dr)
                pair_rows.append({
                    "slice_id": z, "pattern_i": a["pattern_id"], "pattern_j": b["pattern_id"],
                    "episode_i": a["episode_id"], "episode_j": b["episode_id"],
                    "same_current_pattern": int(same), "spatial_iou": iou,
                    "current_sign_aware_dtw": dd, "constrained_sign_dtw_10pct": dc,
                    "lagged_xcorr_distance_pm8": 1.0 - xcorr, "xcorr_best_lag_frames": xlag,
                    "response_field_distance": dr,
                })
        same = np.asarray(same_flags, dtype=bool)
        for method, distances in (("current_dtw", d_dtw), ("constrained_dtw", d_cdtw), ("lagged_xcorr", d_xcorr)):
            arr = np.asarray(distances, dtype=float)
            auc = lower_distance_auc(arr, same) if len(arr) else np.nan
            summary_rows.append({
                "slice_id": z, "method": method, "sampled_modes": len(records),
                "spatially_eligible_cross_episode_pairs": int(len(arr)),
                "same_current_pattern_pairs": int(np.sum(same)),
                "within_vs_between_auc_lower_is_same": auc,
                "median_within_distance": float(np.nanmedian(arr[same])) if np.any(same) else np.nan,
                "median_between_distance": float(np.nanmedian(arr[~same])) if np.any(~same) else np.nan,
                "spearman_vs_current_dtw": float(spearmanr(arr, d_dtw).statistic) if len(arr) > 2 else np.nan,
            })
        print(f"Q2 real slice {z:02d}: sampled {len(records)} modes; eligible pairs {len(same_flags)}", flush=True)
        del records
        gc.collect()
    write_csv(out_dir / "q2_real_pair_distances.csv", pair_rows)
    write_csv(out_dir / "q2_real_slice_summary.csv", summary_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results/q2_metric_benchmark")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    ap.add_argument("--pattern-root", type=Path, default=PATTERN_ROOT)
    ap.add_argument("--slices", type=int, nargs="+", default=list(range(1, 13)))
    ap.add_argument("--max-regions-per-slice", type=int, default=180)
    ap.add_argument("--seed", type=int, default=20261015)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    write_json(args.out / "run_config.json", {
        "experiment": "Q2 DTW versus bounded lagged cross-correlation",
        "repo_root": args.repo_root, "pattern_root": args.pattern_root,
        "slices": args.slices, "max_regions_per_slice": args.max_regions_per_slice,
        "seed": args.seed, "spatial_iou_gate": 0.08, "max_correlation_lag_frames": 8,
        "constrained_dtw_band_fraction": 0.10,
        "real_data_note": "Descriptive sensitivity audit; labels were originally generated using DTW.",
    })
    synth = [metric_row(*case) for case in synthetic_cases()]
    write_csv(args.out / "q2_synthetic_metrics.csv", synth)
    run_real(args.out, args.repo_root, args.pattern_root, args.slices, args.max_regions_per_slice, args.seed)
    print(f"Q2 outputs: {args.out}")


if __name__ == "__main__":
    main()
