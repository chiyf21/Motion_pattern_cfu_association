#!/usr/bin/env python3
"""Test whether Ca onsets / significant pair hits coincide with poorer local membrane consistency."""
from __future__ import annotations

import argparse
import gc
import re
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy.stats import spearmanr, wilcoxon

EXP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXP_ROOT))

from review_common import (
    CFU_ROOT, REPO_ROOT, TEMPORAL_ROOT, T, WINDOW, load_temporal_module,
    read_cfu, read_csv, write_csv, write_json,
)


PS = 7
GAP = 5
H_PATCH, W_PATCH = 244, 329
DEFAULT_REGISTRATION = Path("/mnt/data21T_2/cyf/f338/f338_registrated_0530")
DEFAULT_REFERENCE_DIR = DEFAULT_REGISTRATION / "reference"


def patch_grid(mask):
    mask = np.asarray(mask, dtype=bool)
    yy = np.arange(H_PATCH) - GAP
    xx = np.arange(W_PATCH) - GAP
    valid_y = (yy >= 0) & (yy < mask.shape[0])
    valid_x = (xx >= 0) & (xx < mask.shape[1])
    out = mask[np.ix_(np.clip(yy, 0, mask.shape[0] - 1), np.clip(xx, 0, mask.shape[1] - 1))]
    out[~valid_y, :] = False
    out[:, ~valid_x] = False
    return out


def downsample7(image):
    image = np.asarray(image, dtype=np.float32)
    if image.shape[0] < H_PATCH * PS or image.shape[1] < W_PATCH * PS:
        raise ValueError(f"Registered plane too small for DS7 grid: {image.shape}")
    # Match the current patch pipeline: use floor(H/7) x floor(W/7), dropping
    # the one trailing native-resolution column (2304 is not divisible by 7).
    image = image[:H_PATCH * PS, :W_PATCH * PS]
    return image.reshape(H_PATCH, PS, W_PATCH, PS).mean(axis=(1, 3), dtype=np.float32)


def sample_ca_and_control(onsets, pair_hits, max_events, seed):
    rng = np.random.default_rng(seed)
    all_onsets = np.asarray(sorted(set(map(int, onsets))), dtype=int)
    hit = np.intersect1d(all_onsets, np.asarray(sorted(pair_hits), dtype=int))
    n_pair_hit_total = int(len(hit))
    other = np.setdiff1d(all_onsets, hit)
    if len(hit) > max_events:
        hit = np.sort(rng.choice(hit, size=max_events, replace=False))
    remaining = max(0, max_events - len(hit))
    if len(other) > remaining:
        other = np.sort(rng.choice(other, size=remaining, replace=False))
    event_rows = [(int(t), "ca_onset_pair_hit") for t in hit]
    event_rows += [(int(t), "ca_onset_not_pair_hit") for t in other]

    excluded = np.zeros(T, dtype=bool)
    for t in all_onsets:
        excluded[max(0, t - 8):min(T, t + 9)] = True
    candidates_by_quartile = []
    for q in range(4):
        lo, hi = q * T // 4, (q + 1) * T // 4
        candidates_by_quartile.append(np.flatnonzero(~excluded[lo:hi]) + lo)
    controls, used = [], set()
    for t, _ in event_rows:
        q = min(3, int(4 * t / T))
        choices = np.asarray([x for x in candidates_by_quartile[q] if int(x) not in used], dtype=int)
        if not len(choices):
            choices = np.asarray([x for x in np.flatnonzero(~excluded) if int(x) not in used], dtype=int)
        if len(choices):
            c = int(rng.choice(choices))
            used.add(c)
            controls.append((c, "event_free_baseline"))
    return event_rows, controls, int(len(all_onsets)), n_pair_hit_total


def extract_pair_hit_frames(rows, fdr, repo_root):
    by_slice = {}
    for row in rows:
        z = int(row["motion_slice_0based"])
        by_slice.setdefault(z, {})
    for z in by_slice:
        peaks = {pid: ts for pid, _, ts in fdr.motion_event_peaks(z)}
        by_slice[z] = peaks
    hits = {}
    pair_meta = []
    for row in rows:
        mz, pid = int(row["motion_slice_0based"]), int(row["pattern_id"])
        cz, cid = int(row["cfu_slice_0based"]), int(row["cfu_id"])
        lag = int(row["best_lag"])
        motion_peaks = by_slice.get(mz, {}).get(pid, np.empty(0, dtype=int))
        hit_frames = set()
        for start in np.asarray(motion_peaks, dtype=int) + lag:
            if 0 <= start < T - WINDOW + 1:
                hit_frames.update(range(int(start), int(start + WINDOW)))
        key = (cz, cid)
        hits.setdefault(key, set()).update(hit_frames)
        pair_meta.append((row, key, hit_frames))
    return hits, pair_meta


def mem_patch(frame, z, registration_root):
    path = registration_root / f"membrane/vol_ch1_{frame:06d}.tif"
    mapped = tifffile.memmap(path)
    try:
        return downsample7(np.asarray(mapped[z], dtype=np.float32))
    finally:
        del mapped


def reference_by_frame(reference_dir, n_frames):
    """Map each movie frame to its matching saved registration reference.

    Reference intervals are encoded as vol_ref_<start>_<end>.tif. If two
    intervals touch/overlap at a boundary, use the one with the later start,
    matching the forward-moving registration block at that boundary.
    """
    intervals = []
    for path in sorted(Path(reference_dir).glob("vol_ref_*_*.tif")):
        match = re.fullmatch(r"vol_ref_(\d{6})_(\d{6})", path.stem)
        if match:
            intervals.append((int(match.group(1)), int(match.group(2)), path))
    if not intervals:
        raise FileNotFoundError(f"No vol_ref_<start>_<end>.tif files in {reference_dir}")
    out = {}
    for frame in range(n_frames):
        candidates = [(a, b, p) for a, b, p in intervals if a <= frame <= b]
        if candidates:
            out[frame] = max(candidates, key=lambda item: item[0])[2]
    missing = sorted(set(range(n_frames)) - set(out))
    if missing:
        raise ValueError(f"No saved registration reference for frames: {missing[:10]}")
    return out


def metric(ref, image, mask):
    n = int(mask.sum())
    if n < 8:
        return {"n_mask_pixels": n, "masked_mse": np.nan, "reference_variance": np.nan, "normalized_mse": np.nan}
    ref_values = ref[mask].astype(np.float64)
    image_values = image[mask].astype(np.float64)
    mse = float(np.mean(np.square(image_values - ref_values)))
    ref_variance = float(np.var(ref_values))
    normalized_mse = mse / ref_variance if ref_variance > 1e-12 else np.nan
    return {"n_mask_pixels": n, "masked_mse": mse, "reference_variance": ref_variance, "normalized_mse": normalized_mse}


def bootstrap_mean_ci(values, seed=1, n_boot=5000):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = np.mean(rng.choice(x, size=(n_boot, len(x)), replace=True), axis=1)
    return float(np.mean(x)), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results/q3_membrane_qc_reference_matched_mse")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    ap.add_argument("--cfu-root", type=Path, default=CFU_ROOT)
    ap.add_argument("--registration-root", type=Path, default=DEFAULT_REGISTRATION)
    ap.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR,
                    help="Directory of time-block references; used frame-by-frame by default.")
    ap.add_argument("--reference", type=Path, default=None,
                    help="Optional single-reference override for legacy comparisons only.")
    ap.add_argument("--slices", type=int, nargs="+", default=list(range(1, 13)),
                    help="Slice labels to include; a pair is kept only if both endpoints are selected.")
    ap.add_argument("--max-onset-frames-per-cfu", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261015)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    write_json(args.out / "run_config.json", {
        "experiment": "Membrane-channel local registration-consistency proxy",
        "repo_root": args.repo_root, "cfu_root": args.cfu_root,
        "registration_root": args.registration_root,
        "reference_dir": args.reference_dir,
        "reference_override": args.reference,
        "slices": args.slices,
        "registered_channel": "membrane/vol_ch1_{frame:06d}.tif", "reference_version": "f338_registrated_0530",
        "frame_count": T, "patch_size": PS, "reg_mask_gap": GAP,
        "max_onset_frames_per_cfu": args.max_onset_frames_per_cfu,
        "global_q_threshold": 0.05, "lag_window_from_existing_test": [-8, 8, WINDOW],
        "metrics": ["masked_mse", "normalized_mse"],
        "normalized_mse_definition": "masked_mse / spatial variance of the matching reference pixels inside the CFU mask; NaN if variance <= 1e-12.",
        "reference_policy": "Match each sampled frame to its saved time-block reference; later-start interval wins at shared boundaries.",
        "interpretation": "Appearance-consistency proxy, not a direct deformation-field registration error.",
    })

    fdr = load_temporal_module(args.repo_root)
    result_dir = args.repo_root / "04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset"
    sig_rows = [
        r for r in read_csv(result_dir / "empirical_FDR_significant_pairs.csv")
        if int(r["motion_slice_0based"]) in args.slices and int(r["cfu_slice_0based"]) in args.slices
    ]
    if not sig_rows:
        raise RuntimeError("No q<0.05 pairs remain for the requested slices")
    cfus_by_slice = {}
    for z in sorted({int(r["cfu_slice_0based"]) for r in sig_rows}):
        ids = [int(r["cfu_id"]) for r in sig_rows if int(r["cfu_slice_0based"]) == z]
        cfus_by_slice[z] = read_cfu(z, ids, args.cfu_root)
    pair_hits, pair_meta = extract_pair_hit_frames(sig_rows, fdr, args.repo_root)

    sampled = {}
    sample_plans = {}
    for z, data in cfus_by_slice.items():
        for cid, c in data.items():
            key = (z, cid)
            events, controls, n_onsets, n_pair_hits = sample_ca_and_control(
                np.flatnonzero(c["onset"]), pair_hits.get(key, set()),
                args.max_onset_frames_per_cfu, args.seed + z * 10000 + cid,
            )
            for frame, label in events:
                sampled[(z, cid, frame)] = label
            sample_plans[key] = {
                "events": events, "controls": controls,
                "n_total_onsets": n_onsets, "n_pair_hit_onsets": n_pair_hits,
            }
            for frame, _ in controls:
                sampled[(z, cid, frame)] = "event_free_baseline"

    # Load each unique full 3-D TIFF frame only once, then retain only its DS7 grid.
    frames_by_slice = {}
    for z, cid, frame in sampled:
        frames_by_slice.setdefault(z, set()).add(frame)
    if args.reference is None:
        ref_for_frame = reference_by_frame(args.reference_dir, T)
    else:
        ref_for_frame = {frame: args.reference for frame in range(T)}
    needed_refs = {ref_for_frame[frame] for frames in frames_by_slice.values() for frame in frames}
    needed_slices = set(frames_by_slice)
    ref_by_path_slice = {}
    for ref_path in sorted(needed_refs):
        ref_stack = tifffile.memmap(ref_path)
        for z in needed_slices:
            ref_by_path_slice[(ref_path, z)] = downsample7(np.asarray(ref_stack[z], dtype=np.float32))
        del ref_stack

    mask_by_cfu = {
        (z, cid): patch_grid(c["mask"])
        for z, data in cfus_by_slice.items() for cid, c in data.items()
    }
    rows = []
    frame_cache = {}
    total = sum(len(frames) for frames in frames_by_slice.values())
    done = 0
    for z, frames in frames_by_slice.items():
        for frame in sorted(frames):
            frame_cache[(z, frame)] = mem_patch(frame, z, args.registration_root)
            done += 1
            if done == 1 or done % 100 == 0 or done == total:
                print(f"Q3 membrane frame reads {done}/{total}", flush=True)
    for (z, cid, frame), label in sorted(sampled.items()):
        ref = ref_by_path_slice[(ref_for_frame[frame], z)]
        image, mask = frame_cache[(z, frame)], mask_by_cfu[(z, cid)]
        m = metric(ref, image, mask)
        ca = cfus_by_slice[z][cid]
        rows.append({
            "slice_id": z, "cfu_id": cid, "frame_0based": frame, "sample_class": label,
            "reference_file": ref_for_frame[frame].name,
            "dff_at_frame": float(ca["dff"][frame]), "n_member_events": ca["n_member_events"],
            **m,
        })
    write_csv(args.out / "q3_membrane_frame_metrics.csv", rows)

    per_cfu = []
    for (z, cid), plan in sample_plans.items():
        own = [r for r in rows if int(r["slice_id"]) == z and int(r["cfu_id"]) == cid]
        hit = [r for r in own if r["sample_class"] == "ca_onset_pair_hit"]
        nonhit = [r for r in own if r["sample_class"] == "ca_onset_not_pair_hit"]
        baseline = [r for r in own if r["sample_class"] == "event_free_baseline"]
        all_event = hit + nonhit
        raw_mse = lambda xs: np.asarray([r["masked_mse"] for r in xs], dtype=float)
        norm_mse = lambda xs: np.asarray([r["normalized_mse"] for r in xs], dtype=float)
        dffs = lambda xs: np.asarray([r["dff_at_frame"] for r in xs], dtype=float)
        h_mse, b_mse = raw_mse(all_event), raw_mse(baseline)
        h_nmse, b_nmse = norm_mse(all_event), norm_mse(baseline)
        per_cfu.append({
            "slice_id": z, "cfu_id": cid,
            "n_total_onsets": plan["n_total_onsets"], "n_pair_hit_onsets": plan["n_pair_hit_onsets"],
            "n_sampled_ca_onsets": len(all_event), "n_sampled_pair_hit_onsets": len(hit),
            "n_sampled_baselines": len(baseline),
            "mean_event_masked_mse": float(np.nanmean(h_mse)) if len(h_mse) else np.nan,
            "mean_baseline_masked_mse": float(np.nanmean(b_mse)) if len(b_mse) else np.nan,
            "event_minus_baseline_masked_mse": float(np.nanmean(h_mse) - np.nanmean(b_mse)) if len(h_mse) and len(b_mse) else np.nan,
            "mean_event_normalized_mse": float(np.nanmean(h_nmse)) if len(h_nmse) else np.nan,
            "mean_baseline_normalized_mse": float(np.nanmean(b_nmse)) if len(b_nmse) else np.nan,
            "event_minus_baseline_normalized_mse": float(np.nanmean(h_nmse) - np.nanmean(b_nmse)) if len(h_nmse) and len(b_nmse) else np.nan,
            "within_cfu_spearman_dff_vs_normalized_mse": float(spearmanr(dffs(all_event), h_nmse).statistic) if len(h_nmse) >= 3 and np.nanstd(dffs(all_event)) > 0 else np.nan,
        })
    write_csv(args.out / "q3_cfu_summary.csv", per_cfu)

    lookup = {(int(r["slice_id"]), int(r["cfu_id"]), int(r["frame_0based"])): r for r in rows}
    pair_rows = []
    for row, key, hitframes in pair_meta:
        z, cid = key
        plan = sample_plans[key]
        actual_onsets = set(map(int, np.flatnonzero(cfus_by_slice[z][cid]["onset"])))
        hit_onsets = sorted(actual_onsets.intersection(hitframes))
        sampled_hit = [lookup[(z, cid, t)] for t in hit_onsets if (z, cid, t) in lookup and lookup[(z, cid, t)]["sample_class"] == "ca_onset_pair_hit"]
        sampled_nonhit = [r for r in rows if int(r["slice_id"]) == z and int(r["cfu_id"]) == cid and r["sample_class"] == "ca_onset_not_pair_hit"]
        pair_rows.append({
            "motion_slice": int(row["motion_slice_0based"]), "pattern_id": int(row["pattern_id"]),
            "cfu_slice": z, "cfu_id": cid, "best_lag": int(row["best_lag"]),
            "q_empirical_global": float(row["q_empirical_global"]),
            "n_pair_hit_onsets_total": len(hit_onsets), "n_pair_hit_onsets_sampled": len(sampled_hit),
            "mean_pair_hit_normalized_mse": float(np.nanmean([r["normalized_mse"] for r in sampled_hit])) if sampled_hit else np.nan,
            "mean_other_onset_normalized_mse": float(np.nanmean([r["normalized_mse"] for r in sampled_nonhit])) if sampled_nonhit else np.nan,
            "pair_hit_minus_other_onset_normalized_mse": float(np.nanmean([r["normalized_mse"] for r in sampled_hit]) - np.nanmean([r["normalized_mse"] for r in sampled_nonhit])) if sampled_hit and sampled_nonhit else np.nan,
        })
    write_csv(args.out / "q3_significant_pair_membrane_metrics.csv", pair_rows)

    raw_d = np.asarray([r["event_minus_baseline_masked_mse"] for r in per_cfu], dtype=float)
    norm_d = np.asarray([r["event_minus_baseline_normalized_mse"] for r in per_cfu], dtype=float)
    mean_raw, lo_raw, hi_raw = bootstrap_mean_ci(raw_d, args.seed)
    mean_norm, lo_norm, hi_norm = bootstrap_mean_ci(norm_d, args.seed + 1)
    raw_valid, norm_valid = raw_d[np.isfinite(raw_d)], norm_d[np.isfinite(norm_d)]
    summary = [
        {"metric": "event_minus_baseline_masked_mse", "n_cfUs": len(raw_valid), "mean_difference": mean_raw, "bootstrap_ci95_low": lo_raw, "bootstrap_ci95_high": hi_raw, "paired_wilcoxon_p": float(wilcoxon(raw_valid).pvalue) if len(raw_valid) >= 2 and np.any(raw_valid != 0) else np.nan},
        {"metric": "event_minus_baseline_normalized_mse", "n_cfUs": len(norm_valid), "mean_difference": mean_norm, "bootstrap_ci95_low": lo_norm, "bootstrap_ci95_high": hi_norm, "paired_wilcoxon_p": float(wilcoxon(norm_valid).pvalue) if len(norm_valid) >= 2 and np.any(norm_valid != 0) else np.nan},
    ]
    write_csv(args.out / "q3_population_summary.csv", summary)

    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
        axes[0].boxplot([
            [r["normalized_mse"] for r in rows if r["sample_class"].startswith("ca_onset")],
            [r["normalized_mse"] for r in rows if r["sample_class"] == "event_free_baseline"],
        ], tick_labels=["Ca onset", "Event-free\nbaseline"], showfliers=False)
        axes[0].set_ylabel("Normalized masked MSE\n(lower = closer to reference)")
        axes[1].hist(norm_valid, bins=16, color="#557A95", edgecolor="white")
        axes[1].axvline(0, color="black", lw=1)
        axes[1].set_xlabel("Per-CFU event − baseline normalized MSE")
        axes[1].set_ylabel("CFUs")
        fig.tight_layout()
        fig.savefig(args.out / "q3_membrane_qc_summary.png", dpi=220)
        fig.savefig(args.out / "q3_membrane_qc_summary.pdf")
        plt.close(fig)
    except Exception as exc:
        print(f"Q3 plotting skipped: {exc}", flush=True)
    del frame_cache
    gc.collect()
    print(f"Q3 outputs: {args.out}")


if __name__ == "__main__":
    main()
