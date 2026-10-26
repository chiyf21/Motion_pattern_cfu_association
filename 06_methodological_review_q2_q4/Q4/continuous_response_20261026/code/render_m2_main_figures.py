#!/usr/bin/env python3
"""Render the focused, PNG-only Q4 M2 results from saved tables."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q4-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tables = args.q4_dir / "tables"
    args.out.mkdir(parents=True, exist_ok=True)
    pair_dir = args.out / "significant_pair_png"
    pair_dir.mkdir(exist_ok=True)

    pairs = [r for r in read_csv(tables / "pair_summary.csv") if r.get("status") == "complete"]
    pairs.sort(key=lambda r: float(r.get("q_bh_circular_shift", 1)))
    sig = [r for r in pairs if float(r.get("q_bh_circular_shift", 1)) <= 0.05]
    curves = read_csv(tables / "event_triggered_response_curves.csv")
    curve_by_pair = defaultdict(dict)
    for r in curves:
        key = (r["episode_id"], int(r["motion_peak_frame_0based"]))
        curve_by_pair[r["pair_key"]].setdefault(key, {})[int(r["lag_frames"])] = float(r["baseline_corrected_dff"])
    kernel_rows = read_csv(tables / "kernel_coefficients_descriptive.csv")
    kernel_by_pair = defaultdict(lambda: defaultdict(dict))
    for r in kernel_rows:
        kernel_by_pair[r["pair_key"]][r["model"]][int(r["lag_frames"])] = float(r["kernel_coefficient_dff_per_standardized_predictor"])

    # Figure 1: one dot per candidate; q-significant pairs are labeled.
    fig, ax = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
    for r in pairs:
        q = float(r["q_bh_circular_shift"])
        delta = float(r["delta_mse_M0_minus_M2"])
        passed = q <= 0.05
        y = -np.log10(max(q, 1e-300))
        ax.scatter(delta, y, s=35 if passed else 22,
                   color="#167a53" if passed else "#a9adb2",
                   edgecolor="black" if passed else "none", linewidth=.35,
                   alpha=.9 if passed else .65, zorder=3 if passed else 2)
        if passed:
            ax.annotate(r["pair_key"], (delta, y), xytext=(4, 3), textcoords="offset points", fontsize=6.5)
    ax.axvline(0, color="black", lw=.8)
    ax.axhline(-np.log10(.05), color="#b22222", lw=1, ls="--", label="BH q = 0.05")
    ax.set(xlabel="Mean time-block CV MSE improvement (M0 − M2; positive favors M2)",
           ylabel="−log10(BH-adjusted circular-shift q)",
           title=f"M2 versus M0: {len(sig)}/{len(pairs)} candidate pairs pass q ≤ 0.05")
    ax.grid(alpha=.18)
    ax.legend(frameon=False)
    fig.savefig(args.out / "M2_vs_M0_significance_all_pairs.png", dpi=240)
    plt.close(fig)

    # Figure 2: directly compare the learned M2 lag-response kernels of q-significant pairs.
    lags = np.arange(-8, 9)
    kernels = np.asarray([[kernel_by_pair[r["pair_key"]]["M2_amplitude_weighted"].get(int(l), np.nan)
                           for l in lags] for r in sig], dtype=float)
    lim = float(np.nanmax(np.abs(kernels))) if kernels.size else 1.0
    lim = max(lim, 1e-12)
    fig_h = max(7.0, 0.26 * len(sig) + 2.2)
    fig, ax = plt.subplots(figsize=(11.0, fig_h), constrained_layout=True)
    im = ax.imshow(kernels, aspect="auto", interpolation="nearest", cmap="RdBu_r",
                   vmin=-lim, vmax=lim, extent=[-8.5, 8.5, len(sig)-0.5, -0.5])
    ax.set_xticks(lags)
    ax.set_yticks(np.arange(len(sig)), [r["pair_key"] for r in sig], fontsize=7)
    ax.axvline(0, color="black", lw=.7, alpha=.7)
    ax.set(xlabel="Lag from motion peak (frames; positive = after motion)",
           ylabel="Pair (sorted by q)",
           title="M2 learned lag-response kernels for q-significant pairs\n(full-record descriptive fits; shared color scale)")
    cb = fig.colorbar(im, ax=ax, shrink=.8)
    cb.set_label("Kernel coefficient (ΔdF/F per standardized predictor)")
    fig.savefig(args.out / "M2_kernel_overview_qsig_pairs.png", dpi=240)
    plt.close(fig)

    # Individual pair panels: measured occurrence responses beside the fitted M2 relation.
    for r in sig:
        key = r["pair_key"]
        occ = curve_by_pair.get(key, {})
        traces = []
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
        for occ_key in sorted(occ):
            y = np.asarray([occ[occ_key].get(int(l), np.nan) for l in lags], dtype=float)
            traces.append(y)
            ax[0].plot(lags, y, color="0.65", alpha=.45, lw=.8)
        if traces:
            mat = np.asarray(traces)
            mean = np.nanmean(mat, axis=0)
            sem = np.nanstd(mat, axis=0) / np.maximum(np.sqrt(np.isfinite(mat).sum(axis=0)), 1)
            ax[0].plot(lags, mean, color="#176b87", lw=2, label="mean observed response")
            ax[0].fill_between(lags, mean - sem, mean + sem, color="#176b87", alpha=.18, label="± SEM")
            ax[0].legend(frameon=False, fontsize=8)
        ax[0].axhline(0, color="0.5", lw=.7)
        ax[0].axvline(0, color="black", lw=.8, ls="--")
        ax[0].set(xlabel="Lag from motion peak (frames)", ylabel="Baseline-corrected dF/F",
                  title=f"Observed occurrences (n={len(traces)})")

        m2 = kernel_by_pair.get(key, {}).get("M2_amplitude_weighted", {})
        kval = np.asarray([m2.get(int(l), np.nan) for l in lags])
        ax[1].plot(lags, kval, color="#167a53", lw=2.2, marker="o", ms=3)
        ax[1].fill_between(lags, 0, kval, color="#167a53", alpha=.12)
        ax[1].axhline(0, color="0.5", lw=.7)
        ax[1].axvline(0, color="black", lw=.8, ls="--")
        ax[1].set(xlabel="Lag from motion peak (frames)",
                  ylabel="M2 kernel coefficient\n(ΔdF/F per standardized predictor)",
                  title="M2 fitted lag relationship")
        q = float(r["q_bh_circular_shift"])
        delta = float(r["delta_mse_M0_minus_M2"])
        fig.suptitle(f"{key} | M2 vs M0: BH q={q:.3g}; mean CV ΔMSE={delta:.4g} (positive favors M2)", fontsize=10)
        safe = key.replace("/", "_")
        fig.savefig(pair_dir / f"{safe}_M2.png", dpi=240)
        plt.close(fig)

    (args.out / "README.md").write_text(
        "# Focused Q4 M2 figures\n\n"
        f"This folder contains PNGs only. It answers two questions for the {len(pairs)} frozen Q4 candidate pairs: (1) whether M2 improves over M0 under the existing circular-shift test with BH correction; (2) what lag-response relationship M2 fits for the q≤0.05 pairs.\n\n"
        f"- `M2_vs_M0_significance_all_pairs.png`: all candidate pairs; vertical 0 means no mean time-block CV improvement, and horizontal line is q=0.05. Green points are significant.\n"
        f"- `M2_kernel_overview_qsig_pairs.png`: M2 kernels for the {len(sig)} q-significant pairs on a shared color scale.\n"
        f"- `significant_pair_png/`: one PNG per q-significant pair, showing observed motion-peak-aligned Ca occurrences and the M2 fitted kernel.\n\n"
        "The event-response curves are baseline-corrected observed dF/F. The kernel is a full-record descriptive fit, not a fold-specific estimate; it describes predictive association, not causation. Existing q values come from the 500-shift test and BH correction over the candidate set. No model was refit here.\n",
        encoding="utf-8",
    )
    print(f"candidate_pairs={len(pairs)} significant={len(sig)} pngs={2+len(sig)} output={args.out}")


if __name__ == "__main__":
    main()
