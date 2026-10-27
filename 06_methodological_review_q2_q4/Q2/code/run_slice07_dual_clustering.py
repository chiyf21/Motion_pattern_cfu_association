#!/usr/bin/env python3
"""Compare spatial/DTW with lag-time + unchanged spatial-loss clusters for S07."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgb
import numpy as np
import tifffile
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

WBI_CODE = Path("/home/cyf/wbi/wbi_code")
if WBI_CODE.exists():
    sys.path.insert(0, str(WBI_CODE))
EXP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXP_ROOT))
sys.path.insert(0, str(EXP_ROOT / "01_q2_metric_benchmark"))

import src_registration.motion_correlation_pattern_v2 as mcp  # noqa: E402


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


DEFAULT_CACHE = Path(
    "/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association/"
    "01_motion_pattern_extraction_omega05_mu05/patterns/"
    "Slice07_velocity_decomp/06_patterns"
)
DEFAULT_MEMBRANE = Path(
    "/mnt/data21T_2/cyf/f338/f338_registrated_0530/membrane"
)
T = 1598
PATCH_SIZE = 7
INCOMPATIBLE = 1e6


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def mode_frames(mode, n: int) -> np.ndarray:
    interval = getattr(mode, "time_range", None)
    if interval is None:
        return np.arange(n, dtype=int)
    start, stop = map(int, interval)
    if stop - start == n:
        return np.arange(start, stop, dtype=int)
    if stop - start + 1 == n:
        return np.arange(start, stop + 1, dtype=int)
    return np.arange(start, start + n, dtype=int)


def signed_corr_at_lag(x: np.ndarray, y: np.ndarray, lag: int) -> tuple[float, int]:
    x, y = np.asarray(x, float).reshape(-1), np.asarray(y, float).reshape(-1)
    sx, sy = (lag, 0) if lag >= 0 else (0, -lag)
    n = min(len(x) - sx, len(y) - sy)
    if n < 3:
        return float("nan"), max(n, 0)
    a, b = x[sx:sx + n], y[sy:sy + n]
    if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
        return float("nan"), n
    return float(np.corrcoef(a, b)[0, 1]), n


def same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.array_equal(a[:, None] == a[None, :], b[:, None] == b[None, :]))


def cluster_counts(labels: np.ndarray, pairs: list[tuple[int, int]]) -> tuple[int, int, int]:
    _, counts = np.unique(labels, return_counts=True)
    non_singleton = int(np.sum(counts > 1))
    coassigned = sum(int(labels[i] == labels[j]) for i, j in pairs)
    return len(counts), non_singleton, coassigned


def mode_identifier(mode, index: int):
    value = getattr(mode, "mode_id", None)
    return int(value) if value is not None else index


def mask_on_native_grid(mask: np.ndarray, image_shape: tuple[int, int]) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    out = np.repeat(np.repeat(mask, PATCH_SIZE, axis=0), PATCH_SIZE, axis=1)
    h, w = image_shape
    if out.shape[0] < h or out.shape[1] < w:
        out = np.pad(out, ((0, max(0, h - out.shape[0])),
                           (0, max(0, w - out.shape[1]))))
    return out[:h, :w]


def draw_mode_overlay(ax, image: np.ndarray, mask: np.ndarray, color: str,
                      title: str, box: tuple[int, int, int, int],
                      mask_lowres: np.ndarray) -> None:
    x0, x1, y0, y1 = box
    raw = image[y0:y1, x0:x1]
    lo, hi = np.nanpercentile(raw, [1, 99.5])
    if not np.isfinite(lo + hi) or hi <= lo:
        lo, hi = float(np.nanmin(raw)), float(np.nanmax(raw) + 1)
    ax.imshow(raw, cmap="gray", origin="upper", vmin=lo, vmax=hi,
              interpolation="nearest", extent=(x0, x1, y1, y0))
    local = mask[y0:y1, x0:x1]
    rgba = np.zeros((*local.shape, 4), dtype=np.uint8)
    rgba[..., :3] = np.round(np.asarray(to_rgb(color)) * 255).astype(np.uint8)
    rgba[..., 3] = local.astype(np.uint8) * 102
    ax.imshow(rgba, origin="upper", interpolation="nearest",
              extent=(x0, x1, y1, y0))
    gx0, gy0 = x0 // PATCH_SIZE, y0 // PATCH_SIZE
    gx1 = min(mask_lowres.shape[1], (x1 + PATCH_SIZE - 1) // PATCH_SIZE)
    gy1 = min(mask_lowres.shape[0], (y1 + PATCH_SIZE - 1) // PATCH_SIZE)
    contour_mask = mask_lowres[gy0:gy1, gx0:gx1]
    if np.any(contour_mask):
        ax.contour(contour_mask.astype(np.uint8), levels=[0.5], colors=[color],
                   linewidths=0.7, origin="upper",
                   extent=(gx0 * PATCH_SIZE, min(gx1 * PATCH_SIZE, image.shape[1]),
                           min(gy1 * PATCH_SIZE, image.shape[0]), gy0 * PATCH_SIZE))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_xlabel("Membrane x (native pixels)")
    ax.set_ylabel("Membrane y (native pixels; top = 0)")
    ax.set_title(title, fontsize=9)


def make_example_figure(pair: dict, modes: list, labels_comp: np.ndarray,
                        labels_lag_space: np.ndarray, membrane_dir: Path):
    i = int(pair["mode_index_i_0based"])
    j = int(pair["mode_index_j_0based"])
    mi, mj = modes[i], modes[j]
    hi = np.asarray(mcp._get_region_activation(mi), dtype=float).reshape(-1)
    hj = np.asarray(mcp._get_region_activation(mj), dtype=float).reshape(-1)
    mask_i_low = np.asarray(mcp._get_region_mask_simple(mi), dtype=bool)
    mask_j_low = np.asarray(mcp._get_region_mask_simple(mj), dtype=bool)
    mask_i = mask_on_native_grid(mask_i_low, (1708, 2304))
    mask_j = mask_on_native_grid(mask_j_low, (1708, 2304))

    def peak_frame(mode, trace):
        frames = mode_frames(mode, len(trace))
        keep = (frames >= 0) & (frames < T)
        frames, trace = frames[keep], trace[keep]
        if not len(frames):
            return 0
        return int(frames[np.argmax(np.abs(trace))])

    f_i, f_j = peak_frame(mi, hi), peak_frame(mj, hj)
    img_i = np.asarray(tifffile.memmap(membrane_dir / f"vol_ch1_{f_i:06d}.tif")[7],
                       dtype=np.float32).copy()
    img_j = np.asarray(tifffile.memmap(membrane_dir / f"vol_ch1_{f_j:06d}.tif")[7],
                       dtype=np.float32).copy()
    yy, xx = np.nonzero(mask_i | mask_j)
    pad = 70
    box = (max(0, int(xx.min()) - pad), min(2304, int(xx.max()) + pad + 1),
           max(0, int(yy.min()) - pad), min(1708, int(yy.max()) + pad + 1))

    fig = plt.figure(figsize=(14, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.25, 0.8), hspace=0.38, wspace=0.24)
    ax1, ax2 = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    draw_mode_overlay(ax1, img_i, mask_i, "#2878B5",
                      f"Mode {mode_identifier(mi, i)} (episode {getattr(mi, 'episode_id', -1)}), peak frame {f_i}", box,
                      mask_i_low)
    draw_mode_overlay(ax2, img_j, mask_j, "#E07A1F",
                      f"Mode {mode_identifier(mj, j)} (episode {getattr(mj, 'episode_id', -1)}), peak frame {f_j}", box,
                      mask_j_low)
    ax2.set_ylabel("")
    ax2.tick_params(axis="y", labelleft=False)
    ax = fig.add_subplot(grid[1, :])
    lag = int(pair["best_lag_frames"])
    ax.plot(np.arange(len(hi)), hi, color="#2878B5", lw=1.5,
            label=f"Mode {mode_identifier(mi, i)}")
    ax.plot(np.arange(len(hj)) + lag, hj, color="#E07A1F", lw=1.5,
            label=f"Mode {mode_identifier(mj, j)}, shifted by best lag {lag:+d}")
    ax.axhline(0, color="0.75", lw=0.7)
    ax.set_xlabel("Within-mode frame index (mode 2 shifted for lag alignment)")
    ax.set_ylabel("Signed activation")
    ax.legend(fontsize=9, ncol=2)
    ax.grid(alpha=0.18)

    category = pair["co_cluster_category"]
    cat_title = {
        "both_methods": "Co-clustered by both methods",
        "spatial_dtw_only": "Co-clustered by spatial/DTW only",
        "lag_plus_spatial_only": "Co-clustered by lag-time + spatial only",
    }[category]
    fig.suptitle(f"Slice 7 | {cat_title}", y=0.985, fontsize=12)
    fig.text(
        0.5, 0.945,
        f"IoU={pair['spatial_iou']:.3f}; D_h={pair['D_h']:.3f}; D_b={pair['D_b']:.3f}; "
        f"D_DTW+space={pair['composite_distance']:.3f}; D_lag={pair['lagged_correlation_distance']:.3f}; "
        f"D_lag+space={pair['lag_plus_spatial_distance']:.3f}\n"
        f"Signed r={pair['signed_r_at_best_lag']:.3f} at lag {lag:+d}; "
        f"clusters: spatial/DTW={int(labels_comp[i])+1}, lag+space={int(labels_lag_space[i])+1}",
        ha="center", va="top", fontsize=9)
    fig.subplots_adjust(top=0.87, bottom=0.09, left=0.08, right=0.98,
                        hspace=0.36, wspace=0.18)
    return fig


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--membrane-dir", type=Path, default=DEFAULT_MEMBRANE)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results_lag_plus_spatial")
    ap.add_argument("--composite-cut", type=float, default=0.30,
                    help="Strict complete-linkage cut for 0.5*D_h + 0.5*D_b; canonical cut is 0.45.")
    ap.add_argument("--lag-composite-cut", type=float, default=0.30,
                    help="Strict complete-linkage cut for 0.5*D_lag + 0.5*D_b.")
    ap.add_argument("--overwrite-results", action="store_true",
                    help="Regenerate this runner's outputs after an interrupted/failed run.")
    ap.add_argument("--min-iou", type=float, default=0.08)
    args = ap.parse_args()
    if args.out.exists() and any(args.out.iterdir()) and not args.overwrite_results:
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {args.out}")
    args.out.mkdir(parents=True, exist_ok=True)

    import pickle
    object_path = args.cache / "objects.pkl"
    matrix_path = args.cache / "distance_matrix.npz"
    with object_path.open("rb") as stream:
        data = pickle.load(stream)
    modes = data["kept_regions"]
    n_modes = len(modes)
    info = data["info"]
    composite = np.asarray(info["distance_matrix"], dtype=np.float32)
    saved_matrix = np.load(matrix_path)["distance_matrix"]
    if composite.shape != (n_modes, n_modes) or not np.allclose(composite, saved_matrix):
        raise RuntimeError("Cached mode ordering/matrix mismatch; refusing to cluster misaligned modes.")
    if len(info["pair_info"]) != n_modes * (n_modes - 1) // 2:
        raise RuntimeError("Pair-info table is incomplete for the cached mode list.")

    # Reproduce the saved partition at the canonical threshold as an ordering check.
    original = np.asarray(data["labels"], dtype=int)
    _, original_check = mcp.cluster_regions_hierarchical(
        composite, cluster_dist_thresh=0.45, linkage_method="complete",
        incompatible_dist=INCOMPATIBLE, verbose=False)
    if not same_partition(original, original_check):
        raise RuntimeError("The saved matrix does not reproduce the original 0.45 partition.")

    labels_comp = np.asarray(mcp.cluster_regions_hierarchical(
        composite, cluster_dist_thresh=args.composite_cut,
        linkage_method="complete", incompatible_dist=INCOMPATIBLE,
        verbose=False)[1], dtype=int)

    lag_space_matrix = np.full((n_modes, n_modes), INCOMPATIBLE, dtype=np.float32)
    np.fill_diagonal(lag_space_matrix, 0.0)
    eligible: list[tuple[int, int]] = []
    pair_rows: list[dict] = []
    lag_cache: dict[int, np.ndarray] = {}
    for (i, j), gate in info["pair_info"].items():
        if not gate.get("compatible", False):
            continue
        if i not in lag_cache:
            lag_cache[i] = np.asarray(mcp._get_region_activation(modes[i]), dtype=float).reshape(-1)
        if j not in lag_cache:
            lag_cache[j] = np.asarray(mcp._get_region_activation(modes[j]), dtype=float).reshape(-1)
        hi, hj = lag_cache[i], lag_cache[j]
        abs_r, best_lag = lagged_corr(hi, hj, max_lag=8, min_overlap_fraction=0.8)
        signed_r, n_overlap = signed_corr_at_lag(hi, hj, best_lag)
        lag_dist = 1.0 - float(abs_r)
        lag_plus_spatial_dist = 0.5 * lag_dist + 0.5 * float(gate["D_b"])
        lag_space_matrix[i, j] = lag_space_matrix[j, i] = lag_plus_spatial_dist
        eligible.append((i, j))
        if len(eligible) % 5000 == 0:
            print(f"Computed lag distances for {len(eligible)} gated pairs", flush=True)
        pair_rows.append({
            "mode_index_i_0based": i, "mode_index_j_0based": j,
            "mode_id_i": mode_identifier(modes[i], i),
            "mode_id_j": mode_identifier(modes[j], j),
            "episode_id_i": int(getattr(modes[i], "episode_id", -1)),
            "episode_id_j": int(getattr(modes[j], "episode_id", -1)),
            "spatial_iou": float(gate["iou"]),
            "D_h": float(gate["D_h"]), "D_b": float(gate["D_b"]),
            "composite_distance": float(gate["distance"]),
            "lagged_correlation_distance": lag_dist,
            "lag_plus_spatial_distance": lag_plus_spatial_dist,
            "best_lag_frames": int(best_lag), "signed_r_at_best_lag": signed_r,
            "overlap_frames_at_best_lag": n_overlap,
        })

    if not pair_rows:
        raise RuntimeError("No spatially gated pairs were found.")
    if any(not np.isclose(
            row["lag_plus_spatial_distance"],
            0.5 * row["lagged_correlation_distance"] + 0.5 * row["D_b"],
            rtol=0, atol=1e-7) for row in pair_rows):
        raise RuntimeError("Lag replacement changed the intended 0.5*D_lag + 0.5*cached D_b formula.")
    labels_lag_space = np.asarray(mcp.cluster_regions_hierarchical(
        lag_space_matrix, cluster_dist_thresh=args.lag_composite_cut,
        linkage_method="complete", incompatible_dist=INCOMPATIBLE,
        verbose=False)[1], dtype=int)

    category_counts = Counter()
    for row in pair_rows:
        i, j = int(row["mode_index_i_0based"]), int(row["mode_index_j_0based"])
        a_same, b_same = labels_comp[i] == labels_comp[j], labels_lag_space[i] == labels_lag_space[j]
        category = ("both_methods" if a_same and b_same else
                    "spatial_dtw_only" if a_same else
                    "lag_plus_spatial_only" if b_same else "neither")
        row.update({
            "same_spatial_dtw_cluster": int(a_same),
            "same_lag_plus_spatial_cluster": int(b_same),
            "spatial_dtw_cluster_i": int(labels_comp[i]) + 1,
            "spatial_dtw_cluster_j": int(labels_comp[j]) + 1,
            "lag_plus_spatial_cluster_i": int(labels_lag_space[i]) + 1,
            "lag_plus_spatial_cluster_j": int(labels_lag_space[j]) + 1,
            "co_cluster_category": category,
        })
        category_counts[category] += 1

    mode_rows = []
    for i, mode in enumerate(modes):
        mode_rows.append({
            "mode_index_0based": i,
            "mode_id": mode_identifier(mode, i),
            "episode_id": int(getattr(mode, "episode_id", -1)),
            "cached_cluster_at_canonical_cut_0.45": int(original[i]) + 1,
            "spatial_dtw_cluster_at_cut": int(labels_comp[i]) + 1,
            "lag_plus_spatial_cluster_at_cut": int(labels_lag_space[i]) + 1,
            "activation_frames": len(mcp._get_region_activation(mode)),
            "support_area_patches": int(np.asarray(mcp._get_region_mask_simple(mode), bool).sum()),
        })

    cluster_sizes = {}
    cluster_rows = []
    for name, labels in (("spatial_dtw", labels_comp), ("lag_plus_spatial", labels_lag_space)):
        unique, counts = np.unique(labels, return_counts=True)
        cluster_sizes[name] = counts.tolist()
        for label, size in zip(unique, counts):
            indices = np.flatnonzero(labels == label)
            cluster_rows.append({
                "method": name, "cluster_id": int(label) + 1,
                "n_modes": int(size),
                "mode_indices_0based": ";".join(map(str, indices.tolist())),
                "mode_ids": ";".join(str(mode_identifier(modes[k], k)) for k in indices),
                "episode_ids": ";".join(str(getattr(modes[k], "episode_id", -1)) for k in indices),
            })

    intersections = defaultdict(list)
    for i in range(n_modes):
        intersections[(int(labels_comp[i]), int(labels_lag_space[i]))].append(i)
    size_comp = Counter(map(int, labels_comp))
    size_lag = Counter(map(int, labels_lag_space))
    intersection_rows = []
    for (ca, cb), indices in sorted(intersections.items()):
        overlap = len(indices)
        intersection_rows.append({
            "spatial_dtw_cluster_id": ca + 1, "lag_plus_spatial_cluster_id": cb + 1,
            "n_shared_modes": overlap,
            "spatial_dtw_cluster_size": size_comp[ca],
            "lag_plus_spatial_cluster_size": size_lag[cb],
            "jaccard": overlap / (size_comp[ca] + size_lag[cb] - overlap),
            "mode_indices_0based": ";".join(map(str, indices)),
            "mode_ids": ";".join(str(mode_identifier(modes[k], k)) for k in indices),
        })

    sensitivity_rows = []
    eligible_index_pairs = [(int(r["mode_index_i_0based"]), int(r["mode_index_j_0based"])) for r in pair_rows]
    for method, matrix, thresholds in (
        ("spatial_dtw", composite, (0.20, 0.25, 0.30, 0.35, 0.40, 0.45)),
        ("lag_plus_spatial", lag_space_matrix, (0.20, 0.25, 0.30, 0.35, 0.40, 0.45)),
    ):
        tree = linkage(squareform(matrix, checks=False), method="complete")
        for cut in thresholds:
            current = fcluster(tree, t=cut, criterion="distance")
            groups, non_singletons, coassigned = cluster_counts(current, eligible_index_pairs)
            sensitivity_rows.append({
                "method": method, "distance_cut": cut,
                "equivalent_abs_correlation_minimum": "",
                "n_clusters": groups, "n_non_singleton_clusters": non_singletons,
                "gated_pairs_coassigned": coassigned,
            })

    by_category = defaultdict(list)
    for row in pair_rows:
        by_category[row["co_cluster_category"]].append(row)
    examples = []
    if by_category["both_methods"]:
        examples.append(min(by_category["both_methods"], key=lambda r:
            r["composite_distance"] / args.composite_cut +
            r["lag_plus_spatial_distance"] / args.lag_composite_cut))
    if by_category["spatial_dtw_only"]:
        examples.append(min(by_category["spatial_dtw_only"], key=lambda r:
            r["composite_distance"] / args.composite_cut -
            r["lag_plus_spatial_distance"] / args.lag_composite_cut))
    if by_category["lag_plus_spatial_only"]:
        examples.append(min(by_category["lag_plus_spatial_only"], key=lambda r:
            r["lag_plus_spatial_distance"] / args.lag_composite_cut -
            r["composite_distance"] / args.composite_cut))

    write_csv(args.out / "all_gated_mode_pairs.csv", pair_rows)
    write_csv(args.out / "mode_cluster_assignments.csv", mode_rows)
    write_csv(args.out / "cluster_memberships.csv", cluster_rows)
    write_csv(args.out / "cluster_intersections.csv", intersection_rows)
    write_csv(args.out / "threshold_sensitivity.csv", sensitivity_rows)

    example_rows = []
    if examples:
        with PdfPages(args.out / "representative_cluster_pairs.pdf") as pdf:
            for number, pair in enumerate(examples, start=1):
                fig = make_example_figure(pair, modes, labels_comp, labels_lag_space, args.membrane_dir)
                stem = f"example_{number:02d}_{pair['co_cluster_category']}"
                fig.savefig(args.out / f"{stem}.png", dpi=220, bbox_inches="tight")
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                example_rows.append({"figure": f"{stem}.png", **pair})
    write_csv(args.out / "representative_pairs.csv", example_rows)

    colors = {
        "both_methods": "#2A8C68", "spatial_dtw_only": "#2878B5",
        "lag_plus_spatial_only": "#E07A1F", "neither": "#A0A0A0",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax = axes[0, 0]
    for name, labels, color in (("Spatial/DTW", labels_comp, "#2878B5"),
                               ("Lag time + spatial", labels_lag_space, "#E07A1F")):
        sizes = np.sort(np.unique(labels, return_counts=True)[1])[::-1]
        ax.plot(np.arange(1, len(sizes) + 1), sizes, lw=1.6, color=color,
                label=f"{name} (clusters={len(sizes)})")
    ax.set_yscale("log")
    ax.set_xlabel("Cluster rank (largest first)")
    ax.set_ylabel("Modes per cluster (log scale)")
    ax.set_title("Cluster-size profiles")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.18)

    ax = axes[0, 1]
    category_order = ("both_methods", "spatial_dtw_only", "lag_plus_spatial_only", "neither")
    vals = [category_counts[k] for k in category_order]
    ax.bar(["Both", "Spatial/DTW\nonly", "Lag time + spatial\nonly", "Neither"], vals,
           color=[colors[k] for k in category_order])
    for x, value in enumerate(vals):
        ax.text(x, value, f"{value:,}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Gated mode pairs")
    ax.set_title("Pairwise co-clustering agreement")
    ax.grid(axis="y", alpha=0.18)

    ax = axes[1, 0]
    for category in category_order:
        subset = [r for r in pair_rows if r["co_cluster_category"] == category]
        if not subset:
            continue
        ax.scatter([r["composite_distance"] for r in subset],
                   [r["lag_plus_spatial_distance"] for r in subset],
                   s=8 if category == "neither" else 13,
                   alpha=0.16 if category == "neither" else 0.65,
                   color=colors[category], label=f"{category.replace('_', ' ')} (n={len(subset):,})",
                   rasterized=True)
    ax.axvline(args.composite_cut, color="#2878B5", ls="--", lw=1)
    ax.axhline(args.lag_composite_cut, color="#E07A1F", ls="--", lw=1)
    ax.set_xlabel("Original distance (0.5 D_h + 0.5 D_b)")
    ax.set_ylabel("Replacement distance (0.5 D_lag + 0.5 D_b)")
    ax.set_title("Distances for all spatially gated pairs")
    ax.legend(fontsize=7, loc="best")
    ax.grid(alpha=0.15)

    ax = axes[1, 1]
    for method, color in (("spatial_dtw", "#2878B5"), ("lag_plus_spatial", "#E07A1F")):
        rows = [r for r in sensitivity_rows if r["method"] == method]
        ax.plot([r["distance_cut"] for r in rows], [r["n_clusters"] for r in rows],
                marker="o", color=color,
                label="Spatial/DTW" if method == "spatial_dtw" else "Lag time + spatial")
    ax.axvline(args.composite_cut, color="#2878B5", ls=":", lw=1)
    ax.axvline(args.lag_composite_cut, color="#E07A1F", ls=":", lw=1)
    ax.set_xlabel("Method-specific complete-linkage distance cut")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Threshold sensitivity (cuts are not scale-equivalent)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.18)
    fig.suptitle(
        f"Slice 7 strict clustering comparison | {n_modes:,} modes; "
        f"{len(pair_rows):,} gated pairs; DTW+space cut={args.composite_cut:.2f}; "
        f"lag-time+space cut={args.lag_composite_cut:.2f}",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out / "slice07_dtw_vs_lag_plus_spatial_summary.png", dpi=240, bbox_inches="tight")
    fig.savefig(args.out / "slice07_dtw_vs_lag_plus_spatial_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    reasons = Counter(r.get("reason", "missing") for r in info["pair_info"].values())
    ncomp, ncomp_multi, comp_edges = cluster_counts(labels_comp, eligible_index_pairs)
    nlag, nlag_multi, lag_edges = cluster_counts(labels_lag_space, eligible_index_pairs)
    coassigned_union = (category_counts["both_methods"] +
                        category_counts["spatial_dtw_only"] + category_counts["lag_plus_spatial_only"])
    intersection_multi = [r for r in intersection_rows if r["n_shared_modes"] >= 2]
    exact_multi = [r for r in intersection_multi
                   if r["n_shared_modes"] == r["spatial_dtw_cluster_size"] == r["lag_plus_spatial_cluster_size"]]
    overlap_diagnostics = {}
    for category in ("both_methods", "spatial_dtw_only", "lag_plus_spatial_only", "neither"):
        vals = [int(r["overlap_frames_at_best_lag"]) for r in by_category[category]]
        overlap_diagnostics[category] = {
            "median_overlap_frames": float(np.median(vals)) if vals else None,
            "pairs_with_overlap_le_5": sum(v <= 5 for v in vals),
        }
    summary = {
        "slice_label_1based": 7, "n_gated_modes": n_modes,
        "all_unordered_mode_pairs": len(info["pair_info"]),
        "pair_gate_reasons": dict(reasons),
        "n_spatially_gated_cross_episode_pairs": len(pair_rows),
        "spatial_gate": f"IoU >= {args.min_iou}; same-episode pairs excluded",
        "spatial_dtw_method": {
            "distance": "cached source distance = 0.5 * sign-aware activation DTW + 0.5 * response-field correlation distance",
            "linkage": "complete", "distance_cut": args.composite_cut,
            "cluster_count": ncomp, "non_singleton_cluster_count": ncomp_multi,
            "gated_pairs_coassigned": comp_edges,
        },
        "lag_plus_spatial_method": {
            "distance": "0.5 * (1 - max absolute Pearson r over integer lags -8..+8) + 0.5 * unchanged cached response-field spatial distance D_b; minimum temporal overlap is max(3, ceil(0.8 * shorter trace))",
            "linkage": "complete", "distance_cut": args.lag_composite_cut,
            "spatial_component": "unchanged cached D_b from the canonical pipeline",
            "cluster_count": nlag, "non_singleton_cluster_count": nlag_multi,
            "gated_pairs_coassigned": lag_edges,
        },
        "pairwise_co_cluster_categories": dict(category_counts),
        "pairwise_co_cluster_jaccard": (category_counts["both_methods"] / coassigned_union
                                         if coassigned_union else None),
        "cluster_intersection_diagnostics": {
            "nonempty_cluster_intersections": len(intersection_rows),
            "intersections_sharing_2plus_modes": len(intersection_multi),
            "exact_matching_non_singleton_cluster_pairs": len(exact_multi),
            "largest_shared_intersection": max(r["n_shared_modes"] for r in intersection_rows),
        },
        "lag_overlap_frame_diagnostics": overlap_diagnostics,
        "lag_plus_spatial_only_negative_signed_r_pairs": sum(
            float(r["signed_r_at_best_lag"]) < 0
            for r in by_category["lag_plus_spatial_only"]
        ),
        "original_cache_partition_reproduced_at_cut_0.45": True,
        "interpretation_note": "Only the temporal component was replaced: D_new=0.5*D_lag+0.5*the unchanged cached D_b. Lag correlation uses absolute r, so strong anticorrelation can still count as temporally similar. The hard IoU/episode gate is unchanged.",
    }
    with (args.out / "summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    print(f"Results written to {args.out}", flush=True)


if __name__ == "__main__":
    main()
