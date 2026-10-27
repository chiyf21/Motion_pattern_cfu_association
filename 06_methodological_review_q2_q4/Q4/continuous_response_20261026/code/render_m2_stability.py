#!/usr/bin/env python3
"""Summarize existing Q4 M2 validation and render pair-wise diagnostics."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q4-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tab = args.q4_dir / "tables"
    args.out.mkdir(parents=True, exist_ok=True)

    pairs = [r for r in read_csv(tab / "pair_summary.csv") if r.get("status") == "complete"]
    pair_by_key = {r["pair_key"]: r for r in pairs}
    time_rows = read_csv(tab / "time_block_cv_scores.csv")
    event_rows = [r for r in read_csv(tab / "event_group_cv_scores.csv")
                  if r.get("model") in {"M0_baseline", "M1_binary_events", "M2_amplitude_weighted", "M3_continuous_motion"}]
    curves = read_csv(tab / "event_triggered_response_curves.csv")
    kernels = read_csv(tab / "kernel_coefficients_descriptive.csv")
    response_path = args.q4_dir / "response_metrics" / "pair_response_metrics.csv"
    response = {r["pair_key"]: r for r in read_csv(response_path)} if response_path.exists() else {}

    time_score = {(r["pair_key"], r["fold"], r["model"]): float(r["test_mse"]) for r in time_rows}
    event_score = {(r["pair_key"], r["fold"], r["model"]): float(r["test_mse"]) for r in event_rows}
    time_folds = sorted({r["fold"] for r in time_rows})
    curves_by_pair = defaultdict(dict)
    for r in curves:
        k = (r["episode_id"], int(r["motion_peak_frame_0based"]))
        curves_by_pair[r["pair_key"]].setdefault(k, {})[int(r["lag_frames"])] = float(r["baseline_corrected_dff"])
    kernels_by_pair = defaultdict(lambda: defaultdict(dict))
    for r in kernels:
        kernels_by_pair[r["pair_key"]][r["model"]][int(r["lag_frames"])] = float(r["kernel_coefficient_dff_per_standardized_predictor"])

    summary = []
    for key, pair in pair_by_key.items():
        fold_d = [time_score[(key, f, "M0_baseline")] - time_score[(key, f, "M2_amplitude_weighted")]
                  for f in time_folds if (key, f, "M0_baseline") in time_score and (key, f, "M2_amplitude_weighted") in time_score]
        ev_folds = sorted({r["fold"] for r in event_rows if r["pair_key"] == key and
                           (key, r["fold"], "M0_baseline") in event_score and
                           (key, r["fold"], "M2_amplitude_weighted") in event_score})
        ev_d = [event_score[(key, f, "M0_baseline")] - event_score[(key, f, "M2_amplitude_weighted")]
                for f in ev_folds]
        rr = response.get(key, {})
        summary.append({
            "pair_key": key,
            "q_bh_circular_shift": pair.get("q_bh_circular_shift", ""),
            "q_le_0p05": int(float(pair.get("q_bh_circular_shift", 1)) <= 0.05),
            "n_occurrences": pair.get("n_motion_occurrences", ""),
            "n_event_groups": pair.get("n_event_groups", ""),
            "M0_mean_time_block_mse": float(np.mean([time_score[(key, f, "M0_baseline")] for f in time_folds if (key, f, "M0_baseline") in time_score])),
            "M2_mean_time_block_mse": float(np.mean([time_score[(key, f, "M2_amplitude_weighted")] for f in time_folds if (key, f, "M2_amplitude_weighted") in time_score])),
            "M2_winning_time_blocks": sum(d > 0 for d in fold_d),
            "n_time_blocks": len(fold_d),
            "mean_delta_mse_time_blocks_M0_minus_M2": float(np.mean(fold_d)) if fold_d else np.nan,
            "M2_winning_occurrence_groups": sum(d > 0 for d in ev_d),
            "n_occurrence_group_folds": len(ev_d),
            "fraction_occurrence_group_folds_M2_wins": float(np.mean(np.asarray(ev_d) > 0)) if ev_d else np.nan,
            "median_delta_mse_occurrence_groups_M0_minus_M2": float(np.median(ev_d)) if ev_d else np.nan,
            "response_magnitude_median_delta_dff": rr.get("response_magnitude_median_peak_delta_dff", ""),
            "response_latency_median_frames": rr.get("response_latency_median_frames_ca_minus_motion", ""),
            "occurrence_curve_consistency_median_r": rr.get("curve_consistency_median_pairwise_pearson_r", ""),
        })
    summary.sort(key=lambda r: (int(r["q_le_0p05"]) == 0, float(r["q_bh_circular_shift"] or 1), r["pair_key"]))
    write_csv(args.out / "m2_pair_stability_summary.csv", summary)

    pdf_path = args.out / "m2_pair_diagnostics_all_70.pdf"
    with PdfPages(pdf_path) as pdf:
        for page, row in enumerate(summary, 1):
            key = row["pair_key"]
            p = pair_by_key[key]
            fig, ax = plt.subplots(4, 1, figsize=(10.5, 12.5), constrained_layout=True)
            occ = curves_by_pair.get(key, {})
            lags = np.arange(-8, 9)
            vals = []
            for occ_key in sorted(occ):
                trace = np.asarray([occ[occ_key].get(int(lag), np.nan) for lag in lags], dtype=float)
                vals.append(trace)
                ax[0].plot(lags, trace, color="0.68", alpha=.38, lw=.8)
            if vals:
                mat = np.asarray(vals)
                mean = np.nanmean(mat, axis=0)
                sem = np.nanstd(mat, axis=0) / np.maximum(np.sqrt(np.isfinite(mat).sum(axis=0)), 1)
                ax[0].plot(lags, mean, color="#126e82", lw=2, label="mean ± SEM")
                ax[0].fill_between(lags, mean - sem, mean + sem, color="#126e82", alpha=.18)
            ax[0].axvline(0, color="black", lw=.8, ls="--")
            ax[0].set(xlabel="Lag from motion peak (frames)", ylabel="Baseline-corrected dF/F",
                      title=f"Observed occurrences ({len(occ)}); negative lag = Ca before motion")
            ax[0].legend(frameon=False, loc="best")

            lag_k = np.arange(-8, 9)
            for model, color, label in (("M1_binary_events", "#f2a541", "M1 binary"),
                                        ("M2_amplitude_weighted", "#126e82", "M2 amplitude-weighted"),
                                        ("M3_continuous_motion", "#7a5195", "M3 continuous")):
                values = kernels_by_pair.get(key, {}).get(model, {})
                if values:
                    ax[1].plot(lag_k, [values.get(int(l), np.nan) for l in lag_k], color=color, lw=1.8, label=label)
            ax[1].axhline(0, color="0.5", lw=.7)
            ax[1].axvline(0, color="black", lw=.8, ls="--")
            ax[1].set(xlabel="Lag from motion peak (frames)", ylabel="Estimated kernel coefficient",
                      title="Full-record descriptive kernels (not fold-specific stability estimates)")
            ax[1].legend(frameon=False, ncol=3, loc="best")

            folds = [f for f in time_folds if (key, f, "M0_baseline") in time_score and (key, f, "M2_amplitude_weighted") in time_score]
            delta = [time_score[(key, f, "M0_baseline")] - time_score[(key, f, "M2_amplitude_weighted")] for f in folds]
            ax[2].bar(np.arange(len(folds)), delta, color=["#4c9f70" if d > 0 else "#bb6b6b" for d in delta])
            ax[2].axhline(0, color="black", lw=.8)
            ax[2].set_xticks(np.arange(len(folds)), [f.replace("block_", "B") for f in folds])
            ax[2].set(ylabel="M0 MSE − M2 MSE", title="Contiguous time-block CV; positive favors M2")

            eg = sorted({r["fold"] for r in event_rows if r["pair_key"] == key and
                         (key, r["fold"], "M0_baseline") in event_score and
                         (key, r["fold"], "M2_amplitude_weighted") in event_score})
            egd = [event_score[(key, f, "M0_baseline")] - event_score[(key, f, "M2_amplitude_weighted")] for f in eg]
            ax[3].axhline(0, color="black", lw=.8)
            if egd:
                ax[3].scatter(np.arange(len(egd)), egd, color=["#4c9f70" if d > 0 else "#bb6b6b" for d in egd], s=22)
                ax[3].plot(np.arange(len(egd)), egd, color="0.65", lw=.7, zorder=0)
            ax[3].set_xticks(np.arange(len(eg)), [f"G{i+1}" for i in range(len(eg))] if eg else [])
            ax[3].set(xlabel="Held-out occurrence group", ylabel="M0 MSE − M2 MSE",
                      title="Occurrence-group CV; positive favors M2")

            q = float(p.get("q_bh_circular_shift", 1))
            fig.suptitle(f"M2 diagnostic {page}/{len(summary)}: {key} | q={q:.3g} | "
                         f"time blocks won {row['M2_winning_time_blocks']}/{row['n_time_blocks']} | "
                         f"event groups won {row['M2_winning_occurrence_groups']}/{row['n_occurrence_group_folds']}", fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)

    n_sig = sum(int(r["q_le_0p05"]) for r in summary)
    n_time_stable = sum(int(r["M2_winning_time_blocks"]) >= 3 for r in summary)
    sig_time_stable = sum(int(r["q_le_0p05"]) and int(r["M2_winning_time_blocks"]) >= 3 for r in summary)
    n_event_stable = sum(float(r["fraction_occurrence_group_folds_M2_wins"]) > .5 for r in summary)
    sig_event_stable = sum(int(r["q_le_0p05"]) and float(r["fraction_occurrence_group_folds_M2_wins"]) > .5 for r in summary)
    (args.out / "README.md").write_text(
        "# Q4 M2 stability diagnostics\n\n"
        "Derived only from the saved cluster055 Q4 outputs; no model was refit. The PDF gives one page per pair: individual motion-peak-aligned Ca responses, full-record M1/M2/M3 kernels, time-block fold-wise M0−M2 MSE, and occurrence-group fold-wise M0−M2 MSE. Positive M0−M2 favors M2.\n\n"
        f"- Complete pairs: {len(summary)}; circular-shift BH q≤0.05: {n_sig}.\n"
        f"- M2 wins at least 3/5 contiguous time blocks: {n_time_stable}/{len(summary)}; among q-significant pairs: {sig_time_stable}/{n_sig}.\n"
        f"- M2 wins a majority of held-out occurrence groups: {n_event_stable}/{len(summary)}; among q-significant pairs: {sig_event_stable}/{n_sig}.\n\n"
        "Important: prediction stability is not kernel-shape stability. Saved kernel coefficients were fit on the full recording and are descriptive; fold-specific kernels or occurrence bootstrap intervals would be needed to claim the kernel shape itself is stable. Q-values are inherited from the original 500-shift test and BH correction, not recalculated here. Candidate selection was within the same recording, so this remains exploratory, not causal or independent confirmation.\n",
        encoding="utf-8",
    )
    print(f"pairs={len(summary)}; q<=.05={n_sig}; time-block stable={n_time_stable}; "
          f"q-significant time-block stable={sig_time_stable}; event-group majority wins={n_event_stable}; "
          f"q-significant event-group majority wins={sig_event_stable}")
    print(pdf_path)


if __name__ == "__main__":
    main()
