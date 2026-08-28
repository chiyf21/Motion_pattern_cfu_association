#!/usr/bin/env python3
"""Compare pattern counts across cluster_dist_thresh values.

This is a read-only comparison: it loads existing 04_modes caches and does
not overwrite the current pattern objects.  Results are written to a small
CSV summary for threshold selection.
"""
from __future__ import annotations

import argparse
import csv
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, "/home/cyf/wbi/wbi_code")
# Use the same implementation as rerun_patterns_omega05_mu05.py so that the
# comparison is compatible with the existing pickled mode caches.
import src_registration.motion_correlation_pattern_v2 as mcp
# Portable placeholder: sys.path.insert(0, "/path/to/your/wbi_code")

# Current mode-cache location used by the frozen pattern run.
MODE_ROOT = Path("/home/cyf/wbi/wbi_code/data/f338_velocity_decomp")
# Portable placeholder: MODE_ROOT = Path("/path/to/your/mode_cache")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--thresholds", nargs="+", type=float,
                   default=[0.25, 0.35, 0.45, 0.55, 0.65])
    p.add_argument("--slices", nargs="+", type=int, default=list(range(1, 13)))
    p.add_argument("--output", type=Path,
                   default=Path(__file__).resolve().parent / "cluster_threshold_sweep")
    return p.parse_args()


def load_modes(z):
    path = MODE_ROOT / f"Slice{z:02d}_velocity_decomp" / "04_modes" / "objects.pkl"
    with path.open("rb") as f:
        data = pickle.load(f)
    episodes = data.get("episodes", [])
    modes = data.get("modes", [])
    by_id = {}
    for ep in episodes:
        by_id[getattr(ep, "episode_id", -1)] = ep
        ep.modes = []
    for mode in modes:
        ep = by_id.get(getattr(mode, "episode_id", -1))
        if ep is not None:
            ep.modes.append(mode)
    return episodes, modes


def run_one(episodes, threshold):
    params = dict(
        min_strength=0.0, min_area=5, min_duration=1,
        min_iou=0.08, omega=0.5, mu=0.5,
        b_distance="correlation", spatial_rule="iou",
        cluster_dist_thresh=threshold, linkage_method="complete",
        incompatible_dist=1e6, compute_unified=True,
        unified_mask_mode="best_cc", unified_sign_method="correlation",
        min_pattern_members=2, min_unified_area=50,
    )
    patterns, *_ = mcp.getMotionPattern(
        episodes, unit_type="mode", verbose=False, **params
    )
    sizes = np.asarray([int(getattr(p, "n_members", 0)) for p in patterns], dtype=int)
    return {
        "pattern_count": int(len(sizes)),
        "members_ge5": int(np.sum(sizes >= 5)),
        "members_ge8": int(np.sum(sizes >= 8)),
        "members_ge10": int(np.sum(sizes >= 10)),
        "max_members": int(sizes.max()) if sizes.size else 0,
        "median_members": float(np.median(sizes)) if sizes.size else 0.0,
        "mean_members": float(sizes.mean()) if sizes.size else 0.0,
        "singleton_or_pair": int(np.sum(sizes <= 2)),
    }


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    started = time.time()
    for z in args.slices:
        episodes, modes = load_modes(z)
        for threshold in args.thresholds:
            t0 = time.time()
            stats = run_one(episodes, threshold)
            row = {"slice": z, "cluster_dist_thresh": threshold,
                   "n_episodes": len(episodes), "n_modes": len(modes),
                   "elapsed_seconds": round(time.time() - t0, 2), **stats}
            rows.append(row)
            print(row, flush=True)
    fields = list(rows[0]) if rows else []
    with (args.output / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "README.txt").open("w") as f:
        f.write("Read-only cluster_dist_thresh sweep from existing 04_modes caches.\n")
        f.write("Current pattern objects were not overwritten.\n")
        f.write(f"thresholds={args.thresholds}\nslices={args.slices}\n")
    print(f"Wrote {args.output / 'summary.csv'} in {time.time()-started:.1f}s")


if __name__ == "__main__":
    main()
