#!/usr/bin/env python3
"""Make inspectable Q2 mode examples and Q3 significant-pair trace pages.

Run on the analysis server from the report directory. Q3 uses the corrected
time-block-matched membrane reference; its third panel shows a framewise
normalized-MSE proxy with sampled onset/control points, not ground-truth error.
"""
from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(EXP / "01_q2_metric_benchmark"))
sys.path.insert(0, str(EXP / "03_q3_local_membrane_mse"))

from review_common import (  # noqa: E402
    CFU_ROOT, PATTERN_ROOT, REPO_ROOT, T, WINDOW, configure_imports,
    load_patterns, read_cfu, read_csv, write_csv, write_json,
)
from run_q2_metric_benchmark import lagged_corr  # noqa: E402


Q2_CASES = [
    # Same current cluster, close response fields, but activation traces have
    # very weak bounded-lag absolute correlation.
    {"slice": 5, "pattern_i": 166, "episode_i": 375, "pattern_j": 166, "episode_j": 332},
    # Different current clusters; nearly perfect bounded-lag correlation, but
    # the spatial response field dissimilarity pushes the combined distance up.
    {"slice": 8, "pattern_i": 142, "episode_i": 355, "pattern_j": 153, "episode_j": 299},
]

PS = 7
GAP = 5


def _members(pattern):
    return getattr(pattern, "components", None) or getattr(pattern, "regions", [])


def find_region(z, pid, eid, cache):
    if z not in cache:
        cache[z] = load_patterns(z, PATTERN_ROOT, REPO_ROOT)["patterns"]
    for p in cache[z]:
        if int(p.pattern_id) == pid:
            matches = [m for m in _members(p) if int(getattr(m, "episode_id", -1)) == eid]
            if len(matches) != 1:
                raise RuntimeError(f"Expected one mode for S{z}/P{pid}/episode{eid}; found {len(matches)}")
            return p, matches[0]
    raise KeyError(f"No current pattern S{z}/P{pid}")


def signed_corr_at_lag(x, y, lag):
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    xs, ys = (lag, 0) if lag >= 0 else (0, -lag)
    n = min(len(x) - xs, len(y) - ys)
    a, b = x[xs:xs + n], y[ys:ys + n]
    if n < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def q2_figures(out_dir):
    configure_imports(REPO_ROOT)
    from wholistic_registration.utils.motion_correlation_pattern import (  # noqa: PLC0415
        _get_region_activation,
        _get_region_mask_simple,
        compute_region_distance_matrix_simple,
    )

    source_csv = EXP / "01_q2_metric_benchmark/results/q2_metric_benchmark/q2_real_pair_distances.csv"
    if not source_csv.exists():
        source_csv = HERE / "data/q2_real_pair_distances.csv"
    rows = read_csv(source_csv)
    by_key = {
        (int(r["slice_id"]), int(r["pattern_i"]), int(r["episode_i"]),
         int(r["pattern_j"]), int(r["episode_j"])): r
        for r in rows
    }
    cache = {}
    records = []
    with PdfPages(out_dir / "q2_mode_pair_examples.pdf") as pdf:
        for idx, case in enumerate(Q2_CASES, start=1):
            z = case["slice"]
            p1, a = find_region(z, case["pattern_i"], case["episode_i"], cache)
            p2, b = find_region(z, case["pattern_j"], case["episode_j"], cache)
            row = by_key[(z, case["pattern_i"], case["episode_i"], case["pattern_j"], case["episode_j"])]
            x, y = _get_region_activation(a), _get_region_activation(b)
            ma, mb = _get_region_mask_simple(a), _get_region_mask_simple(b)
            if x is None or y is None or ma is None or mb is None:
                raise RuntimeError(f"Incomplete mode data in Q2 case {idx}")

            matrix, pair_info = compute_region_distance_matrix_simple(
                [a, b], min_iou=0.08, omega=0.5, mu=0.5,
                b_distance="correlation", spatial_rule="iou", verbose=False,
            )
            exact = pair_info[(0, 1)]
            corr_abs, lag = lagged_corr(x, y, max_lag=8)
            signed_r = signed_corr_at_lag(x, y, lag)
            same_label = int(p1.pattern_id) == int(p2.pattern_id)

            fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1, 1.6]})
            labels = ma.astype(np.uint8) + 2 * mb.astype(np.uint8)
            axes[0].imshow(labels, origin="upper", interpolation="nearest",
                           cmap=ListedColormap(["#f7f7f7", "#2878B5", "#E07A1F", "#8F65A8"]),
                           vmin=0, vmax=3)
            axes[0].set_title("Spatial support zoom (absolute DS7 coordinates)")
            axes[0].set_xlabel("patch x")
            axes[0].set_ylabel("patch y; image row order")
            yy, xx = np.nonzero(ma | mb)
            pad = max(8, int(0.15 * max(xx.max() - xx.min() + 1, yy.max() - yy.min() + 1)))
            axes[0].set_xlim(max(-0.5, xx.min() - pad), min(ma.shape[1] - 0.5, xx.max() + pad))
            axes[0].set_ylim(min(ma.shape[0] - 0.5, yy.max() + pad), max(-0.5, yy.min() - pad))
            axes[0].set_aspect("equal")
            axes[0].legend(handles=[
                plt.Line2D([0], [0], color="#2878B5", lw=7, label=f"P{case['pattern_i']} / E{case['episode_i']}"),
                plt.Line2D([0], [0], color="#E07A1F", lw=7, label=f"P{case['pattern_j']} / E{case['episode_j']}"),
            ], loc="upper right", fontsize=8, framealpha=0.9)

            axes[1].plot(np.arange(len(x)), x, color="#2878B5", lw=1.8,
                         label=f"P{case['pattern_i']} E{case['episode_i']} (n={len(x)})")
            axes[1].plot(np.arange(len(y)) + lag, y, color="#E07A1F", lw=1.8,
                         label=f"P{case['pattern_j']} E{case['episode_j']} (lag-aligned, n={len(y)})")
            axes[1].axvline(0, color="0.8", lw=0.7)
            axes[1].set_title("Raw signed activation traces; shifted by best ±8-frame lag")
            axes[1].set_xlabel("within-mode frame index")
            axes[1].set_ylabel("activation")
            axes[1].legend(fontsize=8, loc="best")
            axes[1].grid(alpha=0.2)

            title = (
                f"Q2 example {idx} - slice {z:02d}: "
                f"{'same current pattern' if same_label else 'different current patterns'}\n"
                f"IoU={exact['iou']:.3f}; D_h={exact['D_h']:.3f}; D_b={exact['D_b']:.3f}; "
                f"D=.5D_h+.5D_b={matrix[0, 1]:.3f} (cut=.45); "
                f"bounded-lag r={signed_r:.3f} at lag {lag:+d} (distance={1-corr_abs:.3f})"
            )
            fig.suptitle(title, fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 0.84))
            pdf.savefig(fig, bbox_inches="tight")
            fig.savefig(out_dir / f"q2_example_{idx}.png", dpi=220, bbox_inches="tight")
            plt.close(fig)

            records.append({
                "example": idx, "slice": z,
                "pattern_i": case["pattern_i"], "episode_i": case["episode_i"],
                "pattern_j": case["pattern_j"], "episode_j": case["episode_j"],
                "same_current_pattern": same_label, "source_audit_dtw": row["current_sign_aware_dtw"],
                "source_audit_lagged_xcorr_distance_pm8": row["lagged_xcorr_distance_pm8"],
                "source_audit_best_lag": row["xcorr_best_lag_frames"],
                "exact_pipeline_iou": exact["iou"], "exact_pipeline_D_h": exact["D_h"],
                "exact_pipeline_D_b": exact["D_b"], "exact_pipeline_combined_D": matrix[0, 1],
                "combined_cut": 0.45, "signed_lagged_r": signed_r,
                "activation_len_i": len(x), "activation_len_j": len(y),
                "warning": "Different/same cluster assignment is the observed current label; threshold applies to the weighted combined distance plus complete-linkage clustering.",
            })
    write_csv(out_dir / "q2_mode_pair_examples.csv", records)


def load_current_significant_pairs():
    path = REPO_ROOT / "04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset/empirical_FDR_significant_pairs.csv"
    return [r for r in read_csv(path) if float(r["q_empirical_global"]) < 0.05]


def rms_pattern_trace(pattern, temporal):
    sq = np.zeros(T, dtype=float)
    count = np.zeros(T, dtype=np.int32)
    for component in _members(pattern):
        frames, h = temporal.component_frames(component)
        if len(frames):
            sq[frames] += h * h
            count[frames] += 1
    trace = np.zeros(T, dtype=float)
    active = count > 0
    trace[active] = np.sqrt(sq[active] / count[active])
    return trace


def q3_figures(out_dir, q3_metrics_dir, q3_full_proxy,
               selected_pair_keys=None, shared_mse_scale=False):
    configure_imports(REPO_ROOT)
    from run_q3_membrane_qc import WINDOW as HIT_WINDOW  # noqa: PLC0415
    from review_common import load_temporal_module  # noqa: PLC0415

    metrics_path = q3_metrics_dir / "q3_membrane_frame_metrics.csv"
    metrics = read_csv(metrics_path)
    by_cfu = {}
    for r in metrics:
        by_cfu.setdefault((int(r["slice_id"]), int(r["cfu_id"])), []).append(r)
    dense_proxy = {}
    if q3_full_proxy.exists():
        for r in read_csv(q3_full_proxy):
            dense_proxy.setdefault((int(r["slice_id"]), int(r["cfu_id"])), []).append(r)

    pairs = load_current_significant_pairs()
    pairs.sort(key=lambda r: float(r["q_empirical_global"]))
    if selected_pair_keys:
        pairs = [r for r in pairs if (
            int(r["motion_slice_0based"]), int(r["pattern_id"]),
            int(r["cfu_slice_0based"]), int(r["cfu_id"])
        ) in selected_pair_keys]
        found = {
            (int(r["motion_slice_0based"]), int(r["pattern_id"]),
             int(r["cfu_slice_0based"]), int(r["cfu_id"])) for r in pairs
        }
        missing = selected_pair_keys - found
        if missing:
            raise ValueError(f"Requested Q3 pairs are missing from current q<0.05 results: {sorted(missing)}")
    by_slice = {}
    for r in pairs:
        z = int(r["motion_slice_0based"])
        by_slice.setdefault(z, set()).add(int(r["pattern_id"]))
    temporal = load_temporal_module(REPO_ROOT)
    patterns_by_key = {}
    for z, ids in by_slice.items():
        # Preserve the exact slice label convention used by the existing
        # temporal runner's motion_event_peaks(z) call.
        payload = load_patterns(z, PATTERN_ROOT, REPO_ROOT)
        patterns_by_key[z] = {int(p.pattern_id): p for p in payload["patterns"] if int(p.pattern_id) in ids}
        del payload
        gc.collect()

    cfus_by_slice = {}
    for z in sorted({int(r["cfu_slice_0based"]) for r in pairs}):
        ids = [int(r["cfu_id"]) for r in pairs if int(r["cfu_slice_0based"]) == z]
        cfus_by_slice[z] = read_cfu(z, ids, CFU_ROOT)

    page_rows = []
    common_rmse_ylim = None
    shared_rmse_info = None
    if shared_mse_scale:
        unique_cfus = {(int(r["cfu_slice_0based"]), int(r["cfu_id"])) for r in pairs}
        all_mse = np.asarray([
            float(x["masked_mse"])
            for key in sorted(unique_cfus)
            for x in dense_proxy.get(key, [])
            if np.isfinite(float(x["masked_mse"])) and float(x["masked_mse"]) >= 0
        ], dtype=float)
        if not len(all_mse):
            raise RuntimeError("Cannot determine a shared RMSE scale: no dense local-MSE values found")
        all_rmse = np.sqrt(all_mse)
        upper = float(np.quantile(all_rmse, 0.995))
        clipped = int(np.count_nonzero(all_rmse > upper))
        common_rmse_ylim = (0.0, upper)
        shared_rmse_info = {
            "metric": "sqrt(framewise raw MSE) = RMSE",
            "roi": "existing CFU mask; no ROI or registration-reference recomputation",
            "reference_policy": "frame-matched saved block reference, as in source metrics",
            "axis_rule": "shared y-axis from 0 to the pooled 99.5th percentile of all finite framewise RMSE values across unique CFUs represented by significant pairs",
            "ymin": 0.0,
            "ymax": upper,
            "percentile": 99.5,
            "n_values": int(len(all_rmse)),
            "global_max_rmse": float(np.max(all_rmse)),
            "n_above_axis": clipped,
            "fraction_above_axis": float(clipped / len(all_rmse)),
            "exceedances_are_marked_at_axis_top": True,
        }
        write_json(out_dir / "q3_shared_rmse_axis.json", shared_rmse_info)
        print(
            f"Shared Q3 RMSE axis 0–{upper:.4g} (P99.5); "
            f"{clipped}/{len(all_rmse)} values ({clipped / len(all_rmse):.3%}) exceed it; "
            f"global max={np.max(all_rmse):.4g}", flush=True
        )
    pdf_path = out_dir / "q3_significant_pairs_3traces_reference_matched_rmse_shared_axis.pdf"
    with PdfPages(pdf_path) as pdf:
        for page, r in enumerate(pairs, start=1):
            mz, pid = int(r["motion_slice_0based"]), int(r["pattern_id"])
            cz, cid = int(r["cfu_slice_0based"]), int(r["cfu_id"])
            lag = int(r["best_lag"])
            pattern = patterns_by_key[mz][pid]
            trace = rms_pattern_trace(pattern, temporal)
            prominence = max(float(temporal.MOTION_PROMINENCE_FACTOR) * float(np.std(trace)), 1e-8)
            peaks, _ = temporal.find_peaks(trace, prominence=prominence,
                                           distance=temporal.MOTION_MIN_DISTANCE)
            ca = cfus_by_slice[cz][cid]
            onset = np.flatnonzero(ca["onset"])
            hit_frames = set()
            for start in peaks + lag:
                if 0 <= start < T - HIT_WINDOW + 1:
                    hit_frames.update(range(int(start), int(start + HIT_WINDOW)))
            hit_onsets = set(map(int, onset)).intersection(hit_frames)

            points = by_cfu.get((cz, cid), [])
            xs, ys, classes = [], [], []
            for m in points:
                frame = int(m["frame_0based"])
                if m["sample_class"].startswith("ca_onset"):
                    label = "pair-hit onset" if frame in hit_onsets else "other sampled onset"
                else:
                    label = "event-free control"
                value = float(m["masked_mse"])
                if np.isfinite(value) and value >= 0:
                    xs.append(frame); ys.append(float(np.sqrt(value))); classes.append(label)

            fig, axes = plt.subplots(3, 1, figsize=(11.0, 7.0), sharex=True,
                                     gridspec_kw={"height_ratios": [1.0, 1.25, 1.0]})
            t = np.arange(T)
            axes[0].plot(t, trace, color="#3D6F91", lw=1.0)
            axes[0].scatter(peaks, trace[peaks], color="#C44E52", s=12, zorder=3, label="pattern peak")
            for p in peaks + lag:
                if 0 <= p < T - HIT_WINDOW + 1:
                    axes[0].axvspan(p, p + HIT_WINDOW - 1, color="#C44E52", alpha=0.08, lw=0)
            axes[0].set_ylabel("RMS activation")
            axes[0].legend(loc="upper right", fontsize=8)
            # These CSV columns are mislabeled "0based": the temporal runner
            # loops over z=1..12 and reads Slice{z}/slice_Z{z} directly.
            # Keep figure labels consistent with those actual source IDs.
            axes[0].set_title(f"Q3 significant pair {page}/{len(pairs)}: S{mz:02d} P{pid} → S{cz:02d} CFU{cid:03d}; lag={lag:+d} frames; q={float(r['q_empirical_global']):.3g}", fontsize=10)

            axes[1].plot(t, ca["dff"], color="#2A8C68", lw=0.85)
            axes[1].scatter(onset, ca["dff"][onset], color="#C44E52", s=10,
                            label="CFU onset", zorder=3)
            axes[1].set_ylabel("CFU dF/F")
            axes[1].legend(loc="upper right", fontsize=8)

            colors = {"pair-hit onset": "#C44E52", "other sampled onset": "#777777",
                      "event-free control": "#4C72B0"}
            for label in colors:
                idx = [i for i, x in enumerate(classes) if x == label]
                if idx:
                    axes[2].scatter(np.asarray(xs)[idx], np.asarray(ys)[idx], s=20,
                                    color=colors[label], alpha=0.8, label=label)
            axes[2].set_ylabel("CFU-mask membrane RMSE")
            dense_rows = dense_proxy.get((cz, cid), [])
            if dense_rows:
                dense_rows.sort(key=lambda x: int(x["frame_0based"]))
                dense_x = np.asarray([int(x["frame_0based"]) for x in dense_rows], dtype=int)
                dense_mse = np.asarray([float(x["masked_mse"]) for x in dense_rows], dtype=float)
                dense_y = np.sqrt(np.maximum(dense_mse, 0))
                axes[2].plot(dense_x, dense_y, color="#222222", lw=0.65, alpha=0.75,
                             label="framewise raw masked RMSE", zorder=1)
                if common_rmse_ylim is not None:
                    clipped_here = np.isfinite(dense_y) & (dense_y > common_rmse_ylim[1])
                    if clipped_here.any():
                        axes[2].scatter(dense_x[clipped_here],
                                        np.full(int(clipped_here.sum()), common_rmse_ylim[1]),
                                        marker="^", s=12, color="#B22222", zorder=4,
                                        label=f"above shared P99.5 (n={int(clipped_here.sum())})")
            if common_rmse_ylim is not None:
                axes[2].set_ylim(*common_rmse_ylim)
            axes[2].set_xlabel("frame (0-based)")
            axes[2].legend(loc="upper right", fontsize=8, ncol=3)
            axes[2].text(0.01, 0.02,
                         "Local image MSE proxy; not ground-truth registration error.",
                         transform=axes[2].transAxes, fontsize=8, color="#7A1F1F")
            for ax in axes:
                ax.grid(alpha=0.18)
                ax.set_xlim(0, T - 1)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            page_rows.append({
                "page": page, "motion_slice_0based": mz, "motion_slice_file_id": mz,
                "pattern_id": pid, "cfu_slice_0based": cz,
                "cfu_slice_file_id": cz, "cfu_id": cid, "best_lag": lag,
                "q_empirical_global": r["q_empirical_global"],
                "cfu_rmse_axis_ymin": common_rmse_ylim[0] if common_rmse_ylim else "",
                "cfu_rmse_axis_ymax": common_rmse_ylim[1] if common_rmse_ylim else "",
                "n_pattern_peaks": len(peaks), "n_pair_hit_ca_onsets": len(hit_onsets),
                "n_sampled_proxy_markers": len(xs),
                "n_framewise_proxy_values": len(dense_proxy.get((cz, cid), [])),
                "membrane_metric": "sqrt(raw masked MSE) within existing CFU mask versus frame-matched block reference",
                "n_framewise_values_above_shared_axis": int(np.count_nonzero(
                    np.isfinite(dense_y) & (dense_y > common_rmse_ylim[1])
                )) if dense_rows and common_rmse_ylim is not None else 0,
                "proxy_caveat": "Local image mismatch proxy, not ground-truth registration error.",
            })
    write_csv(out_dir / "q3_significant_pair_trace_index_rmse_shared_axis.csv", page_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=None,
                    help="Optional shared output directory. By default Q2 and Q3 figures go to their own experiment folders.")
    ap.add_argument("--q3-metrics", type=Path,
                    default=EXP / "03_q3_local_membrane_mse/results/q3_membrane_qc_reference_matched_mse")
    ap.add_argument("--q3-full-proxy", type=Path,
                    default=EXP / "05_q3_significant_pair_traces/results/q3_membrane_qc_reference_matched_full_proxy_mse/q3_local_proxy_timeseries.csv")
    ap.add_argument("--q3-only", action="store_true",
                    help="Regenerate only the Q3 figure; leave existing Q2 figures untouched.")
    ap.add_argument("--q3-pair", action="append", default=[], metavar="MOTION_Z,PATTERN,CFU_Z,CFU_ID",
                    help="Limit Q3 output to a pair key; repeat this option for multiple pairs.")
    ap.add_argument("--shared-mse-y", action="store_true",
            help="Use one shared CFU-mask RMSE axis (pooled P99.5) across all selected Q3 pages; exceedances are marked.")
    args = ap.parse_args()
    selected_pairs = set()
    for value in args.q3_pair:
        parts = tuple(int(x) for x in value.split(","))
        if len(parts) != 4:
            ap.error(f"Invalid --q3-pair {value!r}; expected MOTION_Z,PATTERN,CFU_Z,CFU_ID")
        selected_pairs.add(parts)
    if args.out is not None:
        args.out.mkdir(parents=True, exist_ok=True)
        if not args.q3_only:
            q2_figures(args.out)
        q3_figures(args.out, args.q3_metrics, args.q3_full_proxy,
                   selected_pairs or None, args.shared_mse_y)
        print(f"Follow-up figures written to {args.out}")
        return

    if not args.q3_only:
        q2_out = EXP / "01_q2_metric_benchmark/figures/mode_pair_examples"
        q2_out.mkdir(parents=True, exist_ok=True)
        q2_figures(q2_out)
        print(f"Q2 figures written to {q2_out}")
    q3_out = EXP / "05_q3_significant_pair_traces/figures"
    q3_out.mkdir(parents=True, exist_ok=True)
    q3_figures(q3_out, args.q3_metrics, args.q3_full_proxy,
               selected_pairs or None, args.shared_mse_y)
    print(f"Q3 figures written to {q3_out}")


if __name__ == "__main__":
    main()
