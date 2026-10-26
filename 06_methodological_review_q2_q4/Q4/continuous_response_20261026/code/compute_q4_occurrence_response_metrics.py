#!/usr/bin/env python3
"""Summarize observed response magnitude, latency, and curve consistency per Q4 pair."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

LAGS = np.arange(-8, 9)


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def qstats(values):
    a = np.asarray(values, dtype=float)
    return (float(np.median(a)), float(np.quantile(a, .25)), float(np.quantile(a, .75)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q4-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tables = args.q4_dir / "tables"
    curve_rows = read_csv(tables / "event_triggered_response_curves.csv")
    source_rows = read_csv(tables / "occurrence_level_response.csv")

    curves = defaultdict(dict)
    for r in curve_rows:
        key = (r["pair_key"], r["episode_id"], int(r["motion_peak_frame_0based"]))
        lag = int(r["lag_frames"])
        if lag in curves[key]:
            raise ValueError(f"Duplicate lag {lag} for occurrence {key}")
        curves[key][lag] = float(r["baseline_corrected_dff"])
    source = {
        (r["pair_key"], r["episode_id"]): r
        for r in source_rows
    }

    occurrence_rows = []
    grouped = defaultdict(list)
    for (pair, episode, motion_peak), trace_by_lag in sorted(curves.items()):
        if sorted(trace_by_lag) != LAGS.tolist():
            raise ValueError(f"Expected all 17 lags for {pair}, episode {episode}")
        y = np.asarray([trace_by_lag[int(lag)] for lag in LAGS], dtype=float)
        peak_i = int(np.argmax(y))  # earliest lag wins ties
        magnitude = float(y[peak_i])
        src = source[(pair, episode)]
        if abs(magnitude - float(src["event_triggered_peak_delta_dff_pm8"])) > 1e-6:
            raise ValueError(f"Saved occurrence peak does not match lag curve: {pair}/{episode}")
        row = {
            "pair_key": pair,
            "episode_id": episode,
            "motion_peak_frame_0based": motion_peak,
            "response_magnitude_peak_delta_dff": magnitude,
            "response_peak_lag_frames_ca_minus_motion": int(LAGS[peak_i]),
            "response_peak_is_positive": int(magnitude > 0),
            "response_sum_delta_dff_pm8": float(np.sum(y)),
            "cfu_onset_hit_at_frozen_lag_window": src["cfu_onset_hit_at_temporal_best_lag_window"],
        }
        occurrence_rows.append(row)
        grouped[pair].append(y)

    pair_rows = []
    for pair in sorted(grouped):
        occurrences = [r for r in occurrence_rows if r["pair_key"] == pair]
        mags = [float(r["response_magnitude_peak_delta_dff"]) for r in occurrences]
        lats = [int(r["response_peak_lag_frames_ca_minus_motion"]) for r in occurrences]
        traces = grouped[pair]
        correlations = []
        for i in range(len(traces)):
            for j in range(i + 1, len(traces)):
                if np.std(traces[i]) > 1e-12 and np.std(traces[j]) > 1e-12:
                    correlations.append(float(np.corrcoef(traces[i], traces[j])[0, 1]))
        mag_med, mag_q25, mag_q75 = qstats(mags)
        lat_med, lat_q25, lat_q75 = qstats(lats)
        if correlations:
            corr_med, corr_q25, corr_q75 = qstats(correlations)
        else:
            corr_med = corr_q25 = corr_q75 = float("nan")
        pair_rows.append({
            "pair_key": pair,
            "n_occurrences": len(traces),
            "response_magnitude_median_peak_delta_dff": mag_med,
            "response_magnitude_q25_peak_delta_dff": mag_q25,
            "response_magnitude_q75_peak_delta_dff": mag_q75,
            "response_latency_median_frames_ca_minus_motion": lat_med,
            "response_latency_q25_frames": lat_q25,
            "response_latency_q75_frames": lat_q75,
            "curve_consistency_median_pairwise_pearson_r": corr_med,
            "curve_consistency_q25_pairwise_pearson_r": corr_q25,
            "curve_consistency_q75_pairwise_pearson_r": corr_q75,
            "n_valid_pairwise_curve_correlations": len(correlations),
            "fraction_occurrences_with_positive_peak": float(np.mean(np.asarray(mags) > 0)),
        })

    if not pair_rows or not occurrence_rows:
        raise ValueError(f"Unexpected result size: {len(pair_rows)} pairs, {len(occurrence_rows)} occurrences")
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "occurrence_response_metrics.csv", occurrence_rows)
    write_csv(args.out / "pair_response_metrics.csv", pair_rows)
    (args.out / "metric_definitions.md").write_text(
        "# Q4 occurrence response metrics\n\n"
        "Computed from the existing baseline-corrected AQuA2 dF/F curves at lags -8..+8; no model was refit.\n\n"
        "- Response magnitude: maximum baseline-corrected dF/F in the 17-frame window.\n"
        "- Response latency: lag of that maximum, in frames, defined as Ca frame minus motion peak; negative means Ca peak precedes motion. Earliest lag wins ties.\n"
        "- Response consistency: median and IQR of Pearson correlations for every pair of occurrence curves from the same pair, calculated over the 17 aligned lag samples. This measures shape similarity while being insensitive to amplitude scaling.\n"
        "- Pair-level magnitude and latency are reported as median and 25th/75th percentiles across occurrences. A non-positive peak is retained as a signed maximum and flagged; its timing should be interpreted cautiously.\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(occurrence_rows)} occurrence rows and {len(pair_rows)} pair rows to {args.out}")


if __name__ == "__main__":
    main()
