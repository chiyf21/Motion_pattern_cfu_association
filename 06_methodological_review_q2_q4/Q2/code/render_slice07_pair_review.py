#!/usr/bin/env python3
"""Render pair-level overlays for agreement and strongest method disagreements."""
from __future__ import annotations

import argparse
import csv
import pickle
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import run_slice07_dual_clustering as analysis

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results_lag_plus_spatial"
CACHE = analysis.DEFAULT_CACHE
FLOAT_FIELDS = (
    "spatial_iou", "D_h", "D_b", "composite_distance",
    "lagged_correlation_distance", "lag_plus_spatial_distance",
    "signed_r_at_best_lag",
)
INT_FIELDS = (
    "mode_index_i_0based", "mode_index_j_0based", "mode_id_i", "mode_id_j",
    "episode_id_i", "episode_id_j", "best_lag_frames",
    "overlap_frames_at_best_lag", "spatial_dtw_cluster_i",
    "lag_plus_spatial_cluster_i",
)
CATEGORIES = (
    ("both_methods", "01_both_methods", None),
    ("spatial_dtw_only", "02_spatial_dtw_only", "spatial_dtw_only"),
    ("lag_plus_spatial_only", "03_lag_plus_spatial_only", "lag_plus_spatial_only"),
)


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for field in FLOAT_FIELDS:
            row[field] = float(row[field])
        for field in INT_FIELDS:
            row[field] = int(row[field])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit-per-disagreement", type=int, default=100)
    parser.add_argument("--category", choices=("all", *(item[0] for item in CATEGORIES)), default="all")
    parser.add_argument("--out", type=Path, default=RESULTS / "pair_review")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite-existing", action="store_true")
    args = parser.parse_args()
    if args.limit_per_disagreement < 1:
        parser.error("--limit-per-disagreement must be >= 1")

    rows = read_rows(RESULTS / "all_gated_mode_pairs.csv")
    selected = {}
    categories = [item for item in CATEGORIES if args.category == "all" or item[0] == args.category]
    for key, folder, exclusive_key in categories:
        group = [row for row in rows if row["co_cluster_category"] == key]
        if exclusive_key is None:
            group.sort(key=lambda r: (r["mode_index_i_0based"], r["mode_index_j_0based"]))
        else:
            group.sort(key=lambda r: (
                -abs(r["composite_distance"] - r["lag_plus_spatial_distance"]),
                r["mode_index_i_0based"], r["mode_index_j_0based"],
            ))
            group = group[:args.limit_per_disagreement]
        selected[key] = (folder, group)

    print("Selected pair figures:", flush=True)
    for key, (folder, group) in selected.items():
        print(f"  {folder}: {len(group)}", flush=True)
    if args.dry_run:
        return

    with (RESULTS / "mode_cluster_assignments.csv").open(newline="", encoding="utf-8-sig") as stream:
        assignments = list(csv.DictReader(stream))
    labels_comp = np.array([int(r["spatial_dtw_cluster_at_cut"]) - 1 for r in assignments], dtype=int)
    labels_lag = np.array([int(r["lag_plus_spatial_cluster_at_cut"]) - 1 for r in assignments], dtype=int)
    with (CACHE / "objects.pkl").open("rb") as stream:
        modes = pickle.load(stream)["kept_regions"]
    if len(modes) != len(labels_comp) or len(modes) != len(labels_lag):
        raise RuntimeError("Mode cache and cluster-assignment table are not aligned.")

    for category, (folder, group) in selected.items():
        target = args.out / folder / ("all_pair_overlays" if category == "both_methods"
                                      else "top_distance_disagreements")
        target.mkdir(parents=True, exist_ok=True)
        for rank, pair in enumerate(group, start=1):
            i, j = pair["mode_index_i_0based"], pair["mode_index_j_0based"]
            stem = (
                f"rank{rank:03d}_ep{pair['episode_id_i']}_m{pair['mode_id_i']}_idx{i}"
                f"__ep{pair['episode_id_j']}_m{pair['mode_id_j']}_idx{j}.png"
            )
            path = target / stem
            if path.exists() and not args.overwrite_existing:
                continue
            fig = analysis.make_example_figure(pair, modes, labels_comp, labels_lag,
                                               analysis.DEFAULT_MEMBRANE)
            try:
                fig.savefig(path, dpi=220, bbox_inches="tight")
            finally:
                plt.close(fig)
            if rank % 25 == 0 or rank == len(group):
                print(f"{folder}: rendered {rank}/{len(group)}", flush=True)
        gap = "all co-clustered pairs" if category == "both_methods" else (
            f"top {len(group)} by |D_DTW+space - D_lag+space|"
        )
        (target / "selection.txt").write_text(
            f"Category: {category}\nSelected: {len(group)} pairs\nRanking: {gap}\n",
            encoding="utf-8",
        )
    print(f"Pair figures saved under {args.out}", flush=True)


if __name__ == "__main__":
    main()
