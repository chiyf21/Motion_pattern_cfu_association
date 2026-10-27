#!/usr/bin/env python3
"""Pair-specific continuous Motion–Ca response analysis (independent Q4).

This script reads frozen, user-generated motion-pattern caches and filtered
AQuA2 CFUs. It never reruns registration, motion extraction, or AQuA2.
The pickle inputs are Python object caches produced by the checked-in pipeline;
only use them from a trusted source.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANALYSIS = Path("/mnt/data21T_3/cyf/f338_ds7_features_integrated_3prs_20261020/analysis_membrane_qc_z3p5_20261021")
PATTERN_ROOT = ANALYSIS / "01_motion_pattern_extraction_omega05_mu05/patterns"
CFU_ROOT = ANALYSIS / "02_current_cfu_input/cfu"
MODULE_TABLE = ANALYSIS / "05_local_mechanical_modules_distributed_ca_network/02_module_cfu_network/module_cfu_associations_fdr_q005.csv"
Q3_TABLE = ANALYSIS / "06_supplemental/q3_all_significant_pairs/q3_membrane_frame_metrics.csv"
MOTION_REPO = Path("/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association_cluster055/wholistic_registration")
OUT = ANALYSIS / "07_Q2_Q3_Q4_current_CFU/Q4_continuous_response"
N_FRAMES = 1598
MAX_LAG = 8
N_BASIS = 3
RIDGE = 1.0
N_SHIFTS = 500
SEED = 20261022


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def frame_range(component) -> tuple[int, int] | None:
    tr = getattr(component, "time_range", None)
    if tr is None or len(tr) < 2:
        return None
    a, b = int(tr[0]), int(tr[1])
    h = np.asarray(getattr(component, "activation", [])).reshape(-1)
    if b - a + 1 == h.size:
        end = b
    elif b - a == h.size:
        end = b - 1
    else:
        end = a + h.size - 1
    a, end = max(0, a), min(N_FRAMES - 1, end)
    return (a, end) if end >= a else None


def load_patterns(z: int, pattern_root: Path, motion_repo: Path) -> list:
    path = pattern_root / f"Slice{z:02d}_velocity_decomp/06_patterns/objects.pkl"
    if not path.is_file():
        raise FileNotFoundError(path)
    if motion_repo.exists() and str(motion_repo / "src") not in sys.path:
        sys.path.insert(0, str(motion_repo / "src"))
    # Trusted user-produced cache; its generation metadata is preserved next
    # to the object file and the producer is the pinned local motion pipeline.
    with path.open("rb") as f:
        payload = pickle.load(f)
    return payload["patterns"]


def pattern_occurrences(z: int, pattern_root: Path, motion_repo: Path) -> list[dict]:
    """One occurrence per pattern×episode, merging all modes in that episode.

    Amplitude is the peak spatial-RMS magnitude of the sum of the modes'
    reconstructed vector contributions, not the unit-normalized latent h(t).
    Values stay in the input motion units (not converted to physical units).
    """
    out: list[dict] = []
    for pat in load_patterns(z, pattern_root, motion_repo):
        if int(getattr(pat, "n_members", 0)) < 5:
            continue
        comps = list(getattr(pat, "components", None) or getattr(pat, "regions", []) or [])
        by_episode: dict[int, list] = defaultdict(list)
        for c in comps:
            er = frame_range(c)
            if er is not None and getattr(c, "episode_id", None) is not None:
                by_episode[int(c.episode_id)].append(c)
        for ep, modes in by_episode.items():
            ranges = [frame_range(c) for c in modes]
            lo = min(r[0] for r in ranges if r is not None)
            hi = max(r[1] for r in ranges if r is not None)
            vec = None
            for c in modes:
                rm = getattr(c, "reconstructed_motion", None)
                coords = getattr(c, "valid_coords", None)
                er = frame_range(c)
                if rm is None or coords is None or er is None:
                    continue
                rm = np.asarray(rm, dtype=np.float32)
                coords = np.asarray(coords)
                if rm.ndim != 3 or rm.shape[-1] != 2 or coords.ndim != 2 or coords.shape[0] != rm.shape[1]:
                    continue
                ca, cb = er
                expected = cb - ca + 1
                if rm.shape[0] != expected:
                    continue
                if vec is None:
                    vec = np.zeros((N_FRAMES, rm.shape[1], 2), dtype=np.float32)
                # Modes from one episode share valid_coords in the canonical
                # decomposition; reject inconsistent coordinates explicitly.
                if vec.shape[1] != rm.shape[1]:
                    continue
                vec[ca:cb + 1] += rm
            motion = np.zeros(N_FRAMES, dtype=np.float32)
            if vec is not None:
                motion[lo:hi + 1] = np.sqrt(np.mean(np.sum(vec[lo:hi + 1] ** 2, axis=-1), axis=1))
            # Fallback is transparent and flagged: use the mode activation
            # envelope only if a reconstructed vector contribution is absent.
            amp_source = "reconstructed_motion_rms"
            if not np.any(motion[lo:hi + 1] > 0):
                sq = np.zeros(N_FRAMES, dtype=np.float64)
                n = np.zeros(N_FRAMES, dtype=np.int32)
                for c in modes:
                    er = frame_range(c)
                    h = np.asarray(getattr(c, "activation", []), dtype=float).reshape(-1)
                    if er is None or not h.size:
                        continue
                    a, b = er
                    hh = h[:b - a + 1]
                    sq[a:a + len(hh)] += hh * hh
                    n[a:a + len(hh)] += 1
                ok = n > 0
                motion[ok] = np.sqrt(sq[ok] / n[ok]).astype(np.float32)
                amp_source = "latent_activation_rms_fallback"
            peak = lo + int(np.argmax(motion[lo:hi + 1])) if hi >= lo else lo
            amp = float(motion[peak])
            if not np.isfinite(amp) or amp <= 0:
                continue
            out.append({"slice_id": z, "pattern_id": int(pat.pattern_id),
                        "n_members": int(pat.n_members), "episode_id": ep,
                        "start_frame": lo, "end_frame": hi, "peak_frame": peak,
                        "duration_frames": hi - lo + 1, "amplitude": amp,
                        "amplitude_source": amp_source, "motion_envelope": motion})
    return out


def read_cfus(z: int, cfu_root: Path) -> dict[int, dict]:
    path = cfu_root / f"slice_Z{z:02d}_ds7_native_CFU_membrane_qc_z3p5.mat"
    if not path.exists():
        # Input folder may link to the filtered source with its original name.
        hits = sorted(cfu_root.glob(f"slice_Z{z:02d}*.mat"))
        if not hits:
            raise FileNotFoundError(path)
        path = hits[0]
    result = {}
    with h5py.File(path, "r") as f:
        cells = f["cfuInfo1"]
        fields, nc = cells.shape
        layout = "field_first" if fields in (9, 10) else "cfu_first"
        for cid in range(1, nc + 1):
            def get(field0: int):
                ref = cells[field0, cid - 1] if layout == "field_first" else cells[cid - 1, field0]
                return np.asarray(f[ref])
            member = get(1).reshape(-1)
            onset = get(3).reshape(-1).astype(bool)  # MATLAB cfuInfo1 column 4
            dff = get(5).reshape(-1).astype(np.float64)  # MATLAB column 6; column 5 is raw weighted Ca
            if onset.size != N_FRAMES or dff.size != N_FRAMES:
                raise ValueError(f"Slice {z} CFU{cid}: got {onset.size}/{dff.size}, expected {N_FRAMES}")
            result[cid] = {"onset": onset, "dff": dff, "n_member_events": int(member.size)}
    return result


def rbf_basis() -> np.ndarray:
    lags = np.arange(-MAX_LAG, MAX_LAG + 1, dtype=float)
    centers = np.linspace(-MAX_LAG, MAX_LAG, N_BASIS)
    sigma = 4.5
    b = np.column_stack([np.exp(-0.5 * ((lags - c) / sigma) ** 2) for c in centers])
    b /= np.maximum(np.sqrt(np.sum(b * b, axis=0, keepdims=True)), 1e-12)
    return b


BASIS = rbf_basis()
LAGS = np.arange(-MAX_LAG, MAX_LAG + 1)


def convolve_events(events: list[dict], amplitude_weighted: bool) -> np.ndarray:
    x = np.zeros((N_FRAMES, N_BASIS), dtype=np.float64)
    if not events:
        return x
    amps = np.asarray([e["amplitude"] for e in events], float)
    scale = np.median(amps[amps > 0]) if amplitude_weighted and np.any(amps > 0) else 1.0
    for e in events:
        a = (e["amplitude"] / scale) if amplitude_weighted else 1.0
        for lag_i, lag in enumerate(LAGS):
            frame = int(e["peak_frame"]) + int(lag)
            if 0 <= frame < N_FRAMES:
                x[frame] += a * BASIS[lag_i]
    return x


def convolve_continuous(signal: np.ndarray) -> np.ndarray:
    x = np.zeros((N_FRAMES, N_BASIS), dtype=np.float64)
    for lag_i, lag in enumerate(LAGS):
        if lag > 0:
            x[lag:] += signal[:-lag, None] * BASIS[lag_i]
        elif lag < 0:
            x[:lag] += signal[-lag:, None] * BASIS[lag_i]
        else:
            x += signal[:, None] * BASIS[lag_i]
    return x


def baseline_design() -> np.ndarray:
    t = np.linspace(-1, 1, N_FRAMES)
    return np.column_stack([np.ones(N_FRAMES), t, t * t])


def ridge_fit_predict(xtr: np.ndarray, ytr: np.ndarray, xte: np.ndarray,
                      penalize_from: int, ridge: float = RIDGE) -> tuple[np.ndarray, np.ndarray]:
    mu = np.mean(xtr[:, penalize_from:], axis=0) if xtr.shape[1] > penalize_from else np.array([])
    sd = np.std(xtr[:, penalize_from:], axis=0) if xtr.shape[1] > penalize_from else np.array([])
    sd[sd < 1e-10] = 1.0
    a, b = xtr.copy(), xte.copy()
    if len(mu):
        a[:, penalize_from:] = (a[:, penalize_from:] - mu) / sd
        b[:, penalize_from:] = (b[:, penalize_from:] - mu) / sd
    p = xtr.shape[1]
    penalty = np.zeros((p, p), dtype=float)
    if p > penalize_from:
        penalty[penalize_from:, penalize_from:] = np.eye(p - penalize_from) * ridge
    beta = np.linalg.solve(a.T @ a + penalty + np.eye(p) * 1e-10, a.T @ ytr)
    return b @ beta, beta


def block_folds(nfolds: int = 5):
    for i, block in enumerate(np.array_split(np.arange(N_FRAMES), nfolds)):
        lo, hi = int(block[0]), int(block[-1]) + 1
        test = np.zeros(N_FRAMES, bool)
        test[min(hi, lo + MAX_LAG * 2):max(lo, hi - MAX_LAG * 2)] = True
        train = np.ones(N_FRAMES, bool)
        train[max(0, lo - MAX_LAG * 2):min(N_FRAMES, hi + MAX_LAG * 2)] = False
        yield f"block_{i+1}", train, test


def grouped_event_folds(events: list[dict]):
    spans = sorted((max(0, e["peak_frame"] - MAX_LAG), min(N_FRAMES - 1, e["peak_frame"] + MAX_LAG), j)
                    for j, e in enumerate(events))
    groups: list[list[int]] = []
    end = -1
    for lo, hi, j in spans:
        if not groups or lo > end:
            groups.append([j]); end = hi
        else:
            groups[-1].append(j); end = max(end, hi)
    for k, group in enumerate(groups):
        test = np.zeros(N_FRAMES, bool)
        for j in group:
            frame = events[j]["peak_frame"]
            test[max(0, frame - MAX_LAG):min(N_FRAMES, frame + MAX_LAG + 1)] = True
        halo = np.convolve(test.astype(int), np.ones(2 * MAX_LAG + 1, dtype=int), mode="same") > 0
        train = ~halo
        yield f"event_group_{k+1}", train, test, group


def model_matrices(events: list[dict], continuous: np.ndarray) -> dict[str, np.ndarray]:
    base = baseline_design()
    binary = convolve_events(events, False)
    weighted = convolve_events(events, True)
    cont = convolve_continuous(continuous)
    return {"M0_baseline": base,
            "M1_binary_events": np.column_stack([base, binary]),
            "M2_amplitude_weighted": np.column_stack([base, weighted]),
            "M3_continuous_motion": np.column_stack([base, cont])}


def score_models(y: np.ndarray, matrices: dict[str, np.ndarray], fold_iter) -> tuple[list[dict], dict[str, list]]:
    rows: list[dict] = []
    pred_store: dict[str, list] = defaultdict(list)
    for fold_name, train, test, *rest in fold_iter:
        if train.sum() < 100 or test.sum() < 10:
            continue
        for model, x in matrices.items():
            pred, _ = ridge_fit_predict(x[train], y[train], x[test], penalize_from=3)
            err = y[test] - pred
            rows.append({"fold": fold_name, "model": model, "n_train_frames": int(train.sum()),
                         "n_test_frames": int(test.sum()), "test_mse": float(np.mean(err ** 2)),
                         "test_rmse": float(np.sqrt(np.mean(err ** 2))),
                         "test_r2": float(1 - np.sum(err ** 2) / max(np.sum((y[test] - np.mean(y[train])) ** 2), 1e-12))})
            pred_store[model].append((test.copy(), pred.copy(), y[test].copy()))
    return rows, pred_store


def aggregate_prediction_score(rows: list[dict], model: str) -> float:
    r = [x["test_mse"] for x in rows if x["model"] == model]
    return float(np.mean(r)) if r else float("nan")


def q3_summary(q3_path: Path, z: int, cid: int) -> dict:
    if not q3_path.exists():
        return {"q3_n_frames": 0, "q3_median_normalized_mse": np.nan, "q3_p95_normalized_mse": np.nan,
                "q3_pair_hit_median_normalized_mse": np.nan, "q3_other_median_normalized_mse": np.nan}
    vals, hit, other = [], [], []
    for r in read_csv(q3_path):
        if int(r.get("slice_id", -1)) != z or int(r.get("cfu_id", -1)) != cid:
            continue
        v = float(r.get("normalized_mse", "nan"))
        if not np.isfinite(v):
            continue
        vals.append(v)
        cls = r.get("sample_class", "").lower()
        (hit if "hit" in cls else other).append(v)
    med = lambda a: float(np.median(a)) if a else float("nan")
    return {"q3_n_frames": len(vals), "q3_median_normalized_mse": med(vals),
            "q3_p95_normalized_mse": float(np.quantile(vals, .95)) if vals else np.nan,
            "q3_pair_hit_median_normalized_mse": med(hit), "q3_other_median_normalized_mse": med(other)}


def bh_qvalues(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def event_hit(onset: np.ndarray, peak: int, best_lag: int, window: int) -> bool:
    lo = peak + best_lag
    hi = lo + max(1, window)
    lo, hi = max(0, lo), min(N_FRAMES, hi)
    return bool(np.any(onset[lo:hi]))


def event_response_rows(events: list[dict], dff: np.ndarray, onset: np.ndarray, best_lag: int, window: int,
                        pair_key: str) -> tuple[list[dict], list[np.ndarray], list[dict]]:
    rows, curves, curve_rows = [], [], []
    for e in events:
        pk = int(e["peak_frame"])
        base_lo, base_hi = pk - 16, pk - 8
        lag_frames = np.arange(-MAX_LAG, MAX_LAG + 1)
        idx = pk + lag_frames
        curve = np.full(len(lag_frames), np.nan)
        ok = (idx >= 0) & (idx < N_FRAMES)
        curve[ok] = dff[idx[ok]]
        if base_lo >= 0:
            baseline = float(np.mean(dff[base_lo:base_hi]))
            delta = curve - baseline
            curves.append(delta)
            peak_resp = float(np.nanmax(delta))
            area = float(np.nansum(delta))
        else:
            baseline, peak_resp, area = np.nan, np.nan, np.nan
            curves.append(curve)
        for lag, value in zip(lag_frames, curves[-1]):
            curve_rows.append({"pair_key": pair_key, "motion_peak_frame_0based": pk,
                               "episode_id": e["episode_id"], "lag_frames": int(lag),
                               "ca_frame_0based": int(pk + lag),
                               "baseline_corrected_dff": float(value) if np.isfinite(value) else np.nan})
        rows.append({"pair_key": pair_key, "motion_slice_0based": e["slice_id"],
                     "pattern_id": e["pattern_id"], "episode_id": e["episode_id"],
                     "cfu_peak_frame_0based": pk, "motion_start_frame_0based": e["start_frame"],
                     "motion_end_frame_0based": e["end_frame"], "motion_duration_frames": e["duration_frames"],
                     "motion_amplitude": e["amplitude"], "amplitude_source": e["amplitude_source"],
                     "cfu_prebaseline_dff_mean_frames_minus16_to_minus9": baseline,
                     "event_triggered_peak_delta_dff_pm8": peak_resp,
                     "event_triggered_sum_delta_dff_pm8": area,
                     "cfu_onset_hit_at_temporal_best_lag_window": int(event_hit(onset, pk, best_lag, window)),
                     "temporal_best_lag_frames": best_lag, "temporal_window_frames": window})
    return rows, curves, curve_rows


def circular_shift_test(y: np.ndarray, matrices: dict[str, np.ndarray], n_shifts: int,
                        rng: np.random.Generator) -> tuple[float, float, np.ndarray]:
    # Primary statistic is block-CV MSE improvement M0 minus M2.
    obs_rows, _ = score_models(y, matrices, list(block_folds()))
    obs = aggregate_prediction_score(obs_rows, "M0_baseline") - aggregate_prediction_score(obs_rows, "M2_amplitude_weighted")
    offsets = rng.integers(MAX_LAG * 2 + 1, N_FRAMES - MAX_LAG * 2, size=n_shifts)
    # All folds use fixed motion designs. Precompute the ridge least-squares
    # maps once, then evaluate all shifted Ca traces as a batch. This is
    # numerically equivalent to refitting every fold for every circular shift.
    shifted = np.stack([np.roll(y, int(s)) for s in offsets], axis=0)
    null_m0 = np.zeros(n_shifts, dtype=float)
    null_m2 = np.zeros(n_shifts, dtype=float)
    fold_count = 0
    for _, train, test in block_folds():
        if train.sum() < 100 or test.sum() < 10:
            continue
        fold_count += 1
        for model, target_acc in (("M0_baseline", null_m0), ("M2_amplitude_weighted", null_m2)):
            x = matrices[model]
            xtr, xte = x[train].copy(), x[test].copy()
            if xtr.shape[1] > 3:
                mu, sd = xtr[:, 3:].mean(axis=0), xtr[:, 3:].std(axis=0)
                sd[sd < 1e-10] = 1.0
                xtr[:, 3:] = (xtr[:, 3:] - mu) / sd
                xte[:, 3:] = (xte[:, 3:] - mu) / sd
            penalty = np.zeros((xtr.shape[1], xtr.shape[1]), dtype=float)
            if xtr.shape[1] > 3:
                penalty[3:, 3:] = np.eye(xtr.shape[1] - 3) * RIDGE
            inv_map = np.linalg.solve(xtr.T @ xtr + penalty + np.eye(xtr.shape[1]) * 1e-10, xtr.T)
            beta_batch = inv_map @ shifted[:, train].T
            pred_batch = xte @ beta_batch
            err = shifted[:, test].T - pred_batch
            target_acc += np.mean(err * err, axis=0) / 5.0
    null = null_m0 - null_m2
    p = float((1 + np.sum(null >= obs)) / (1 + n_shifts))
    return obs, p, null


def plot_pair(pair_dir: Path, key: str, events: list[dict], dff: np.ndarray, onset: np.ndarray,
              event_curves: list[np.ndarray], cv_rows: list[dict], null: np.ndarray,
              delta_mse: float, p_shift: float, beta_kernels: dict[str, np.ndarray], title: str) -> None:
    pair_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(4, 1, figsize=(12, 13), constrained_layout=True)
    xlag = np.arange(-MAX_LAG, MAX_LAG + 1)
    if event_curves:
        mat = np.asarray(event_curves)
        for curve in mat:
            ax[0].plot(xlag, curve, color="#9aa0a6", alpha=.22, lw=.8)
        mean = np.nanmean(mat, axis=0)
        se = np.nanstd(mat, axis=0) / max(np.sqrt(np.sum(np.isfinite(mat[:, 0]))), 1)
        ax[0].plot(xlag, mean, color="#126e82", lw=2.5, label="Mean baseline-corrected dF/F")
        ax[0].fill_between(xlag, mean - se, mean + se, color="#126e82", alpha=.18)
    ax[0].axvline(0, color="black", lw=1, ls="--")
    ax[0].set(xlabel="Lag from motion peak (frames; positive = Ca after motion)", ylabel="Δ dF/F", title="Occurrence-centered continuous Ca response")
    ax[0].legend(loc="best", frameon=False)
    t = np.arange(N_FRAMES)
    ax[1].plot(t, dff, color="#167a53", lw=.8, label="AQuA2 CFU dF/F")
    ax[1].scatter(t[onset], dff[onset], color="#d24b4b", s=12, label="AQuA2 CFU onset")
    for e in events:
        ax[1].axvspan(e["start_frame"], e["end_frame"], color="#4c78a8", alpha=.08)
        ax[1].axvline(e["peak_frame"], color="#4c78a8", alpha=.22, lw=.6)
    ax[1].set(xlabel="Frame (0-based)", ylabel="dF/F", title="Full-recording CFU trace and pattern occurrence intervals")
    ax[1].legend(loc="upper right", ncol=2, frameon=False)
    models = ["M0_baseline", "M1_binary_events", "M2_amplitude_weighted", "M3_continuous_motion"]
    vals = [aggregate_prediction_score(cv_rows, m) for m in models]
    ax[2].bar(np.arange(4), vals, color=["#888888", "#f2a541", "#126e82", "#7a5195"])
    ax[2].set_xticks(np.arange(4), ["M0\nbaseline", "M1\nbinary", "M2\namplitude-weighted", "M3\ncontinuous"])
    ax[2].set_ylabel("Blocked-CV MSE (lower is better)")
    ax[2].set_title(f"Held-out full-recording MSE; M0−M2={delta_mse:.4g}; circular-shift p={p_shift:.4g}")
    if len(null):
        ax[3].hist(null, bins=30, color="#b7c9d6", edgecolor="white")
        ax[3].axvline(delta_mse, color="#c33", lw=2, label="Observed M0−M2")
        ax[3].axvline(0, color="black", lw=.8)
        ax[3].legend(frameon=False)
    ax[3].set(xlabel="Blocked-CV MSE improvement under circular shifts", ylabel="Shifts", title="Autocorrelation-preserving circular-shift null")
    fig.suptitle(title, fontsize=13)
    fig.savefig(pair_dir / f"{key}_continuous_response.png", dpi=200)
    fig.savefig(pair_dir / f"{key}_continuous_response.pdf")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis-root", type=Path, default=ANALYSIS)
    ap.add_argument("--pattern-root", type=Path, default=PATTERN_ROOT)
    ap.add_argument("--cfu-root", type=Path, default=CFU_ROOT)
    ap.add_argument("--module-table", type=Path, default=MODULE_TABLE)
    ap.add_argument("--q3-table", type=Path, default=Q3_TABLE)
    ap.add_argument("--motion-repo", type=Path, default=MOTION_REPO)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--shifts", type=int, default=N_SHIFTS)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    started = time.time()
    args.out.mkdir(parents=True, exist_ok=True)
    for d in ("pair_figures", "tables", "diagnostics", "synthetic_tests"):
        (args.out / d).mkdir(exist_ok=True)
    config = {"experiment": "pair-specific continuous Motion-Ca response; independent from binary-onset Q4",
              "analysis_root": str(args.analysis_root), "pattern_root": str(args.pattern_root),
              "cfu_root": str(args.cfu_root), "module_table": str(args.module_table),
              "q3_table": str(args.q3_table), "motion_repo": str(args.motion_repo), "out": str(args.out), "n_frames": N_FRAMES,
              "candidate_scope": "frozen spatial-module edges that are globally significant in onset co-occurrence (module table)",
              "pattern_filter": "n_members >= 5", "motion_occurrence": "one unique pattern x episode_id; all modes in an episode are combined",
              "occurrence_start_end": "union of member mode time_range within that episode; peak is max spatial RMS of summed reconstructed_motion",
              "amplitude": "peak spatial-RMS magnitude of summed reconstructed vector contribution; input motion units, not physical-unit calibrated",
              "response": "AQuA2 cfuInfo1 column 6 (1-based), continuous dF/F; column 5 is raw weighted fluorescence and is not used",
              "baseline": "per-occurrence mean dF/F at lags -16 through -9; response display window -8 through +8",
              "lag_sign": "negative: Ca precedes motion; positive: Ca follows motion",
              "basis": {"type": "three L2-normalized Gaussian radial basis functions", "centers_frames": [-8, 0, 8], "sigma_frames": 4.5},
              "models": ["M0 baseline intercept+linear+quadratic drift", "M1 binary event impulse convolution", "M2 amplitude-weighted event impulse convolution (primary)", "M3 continuous reconstructed-motion envelope convolution (exploratory)"],
              "validation": "5 contiguous time-block folds with 16-frame guards; plus occurrence-window group folds merged when ±8-frame windows overlap",
              "circular_shift": {"n_shifts": args.shifts, "random_seed": args.seed, "min_offset_frames": 17, "statistic": "block-CV MSE(M0)-MSE(M2); positive favors M2"},
              "ridge_lambda": RIDGE, "q_adjustment": "BH over tested candidate pairs, separately for circular-shift p values",
              "selection_caveat": "Candidate pairs were discovered in this recording, so this is exploratory conditional analysis, not independent confirmation or a causal test."}
    (args.out / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    module_rows = read_csv(args.module_table)
    q3_cache: dict[tuple[int, int], dict] = {}
    pattern_cache: dict[int, list[dict]] = {}
    cfu_cache: dict[int, dict[int, dict]] = {}
    for z in sorted({int(r["module_slice_0based"]) for r in module_rows} | {int(r["source_cfu_slice_0based"]) for r in module_rows}):
        pattern_cache[z] = pattern_occurrences(z, args.pattern_root, args.motion_repo)
        cfu_cache[z] = read_cfus(z, args.cfu_root)
        print(f"Loaded slice {z:02d}: {len(pattern_cache[z])} pattern-episode occurrences; {len(cfu_cache[z])} CFUs", flush=True)
    by_pat = {(z, p): [e for e in evs if e["pattern_id"] == p]
              for z, evs in pattern_cache.items() for p in {e["pattern_id"] for e in evs}}

    pair_rows, occurrence_rows, response_curve_rows, cv_rows_all, group_rows_all, kernel_rows, null_rows = [], [], [], [], [], [], []
    rng = np.random.default_rng(args.seed)
    seen = set()
    for ix, edge in enumerate(module_rows, 1):
        pz, pid = int(edge["module_slice_0based"]), int(edge["module_pattern_id"])
        cz, cid = int(edge["source_cfu_slice_0based"]), int(edge["source_cfu_id"])
        key_tuple = (pz, pid, cz, cid)
        if key_tuple in seen:
            continue
        seen.add(key_tuple)
        key = f"S{pz:02d}_P{pid}__S{cz:02d}_CFU{cid:03d}"
        events = sorted(by_pat.get((pz, pid), []), key=lambda e: e["peak_frame"])
        if not events:
            pair_rows.append({"pair_key": key, "status": "no_frozen_occurrences_for_pattern", **edge})
            continue
        if cid not in cfu_cache[cz]:
            pair_rows.append({"pair_key": key, "status": "CFU_missing", **edge})
            continue
        cf = cfu_cache[cz][cid]
        dff, onset = cf["dff"], cf["onset"]
        best_lag, window = int(float(edge.get("best_lag", 0))), int(float(edge.get("window_frames", 3)))
        evrows, evcurves, evcurve_rows = event_response_rows(events, dff, onset, best_lag, window, key)
        occurrence_rows.extend(evrows)
        response_curve_rows.extend(evcurve_rows)
        motion = np.zeros(N_FRAMES, dtype=np.float64)
        for e in events:
            motion += e["motion_envelope"]
        matrices = model_matrices(events, motion)
        folds = list(block_folds())
        cv, _ = score_models(dff, matrices, folds)
        cv_rows_all.extend([{"pair_key": key, **r} for r in cv])
        group_folds = list(grouped_event_folds(events))
        if len(group_folds) >= 3:
            gr, _ = score_models(dff, matrices, [(n, tr, te) for n, tr, te, _ in group_folds])
            group_rows_all.extend([{"pair_key": key, "n_occurrence_groups": len(group_folds), **r} for r in gr])
            event_cv_status = "evaluated"
        else:
            gr, event_cv_status = [], "insufficient_independent_groups_lt3"
        for gname, tr, te, group in group_folds:
            group_rows_all.append({"pair_key": key, "fold": gname, "model": "group_membership",
                                   "n_train_frames": int(tr.sum()), "n_test_frames": int(te.sum()),
                                   "test_mse": np.nan, "test_rmse": np.nan, "test_r2": np.nan,
                                   "occurrence_indices": ";".join(map(str, group))})
        obs, p_shift, null = circular_shift_test(dff, matrices, args.shifts, rng)
        null_rows.extend([{"pair_key": key, "shift_index": j + 1, "delta_mse_m0_minus_m2": float(v)} for j, v in enumerate(null)])
        m0, m1, m2, m3 = [aggregate_prediction_score(cv, m) for m in ("M0_baseline", "M1_binary_events", "M2_amplitude_weighted", "M3_continuous_motion")]
        # Descriptive full-data basis kernels, never used as CV scores.
        for model in ("M1_binary_events", "M2_amplitude_weighted", "M3_continuous_motion"):
            x = matrices[model]
            _, beta = ridge_fit_predict(x, dff, x, penalize_from=3)
            b = beta[3:]
            kernel = BASIS @ b.reshape(N_BASIS, 1)
            for lag, kval in zip(LAGS, kernel[:, 0]):
                kernel_rows.append({"pair_key": key, "model": model, "lag_frames": int(lag), "kernel_coefficient_dff_per_standardized_predictor": float(kval)})
        if (cz, cid) not in q3_cache:
            q3_cache[(cz, cid)] = q3_summary(args.q3_table, cz, cid)
        q3 = q3_cache[(cz, cid)]
        pair_rows.append({"pair_key": key, "status": "complete", "motion_slice_0based": pz,
                          "pattern_id": pid, "pattern_n_members": int(events[0]["n_members"]),
                          "cfu_slice_0based": cz, "cfu_id": cid, "cfu_n_member_events": cf["n_member_events"],
                          "n_motion_occurrences": len(events), "n_cfu_onsets": int(onset.sum()),
                          "n_onset_hits_at_temporal_best_lag": int(sum(r["cfu_onset_hit_at_temporal_best_lag_window"] for r in evrows)),
                          "median_motion_amplitude": float(np.median([e["amplitude"] for e in events])),
                          "amplitude_sources": ";".join(sorted(set(e["amplitude_source"] for e in events))),
                          "best_lag_onset_analysis_frames": best_lag, "onset_window_frames": window,
                          "M0_block_cv_mse": m0, "M1_block_cv_mse": m1, "M2_block_cv_mse": m2, "M3_block_cv_mse": m3,
                          "delta_mse_M0_minus_M2": m0 - m2, "delta_mse_M1_minus_M2": m1 - m2,
                          "circular_shift_p_m0_vs_m2": p_shift, "event_group_cv_status": event_cv_status,
                          "n_event_groups": len(group_folds), "q_empirical_global_onset": edge.get("q_empirical_global", ""),
                          "module_ratio": edge.get("module_ratio", ""), "module_coverage": edge.get("module_coverage", ""), **q3})
        plot_pair(args.out / "pair_figures", key, events, dff, onset, evcurves, cv, null, obs, p_shift,
                  {}, f"{key} | continuous Motion–Ca response")
        if ix % 5 == 0 or ix == len(module_rows):
            print(f"Q4 continuous pair {ix}/{len(module_rows)}: {key}; occurrences={len(events)}; elapsed={time.time()-started:.1f}s", flush=True)

    # Pair-level empirical p values are corrected only over the tested frozen candidates.
    completed = [r for r in pair_rows if r.get("status") == "complete"]
    if completed:
        pvals = np.asarray([float(r["circular_shift_p_m0_vs_m2"]) for r in completed])
        qvals = bh_qvalues(pvals)
        for r, q in zip(completed, qvals):
            r["q_bh_circular_shift"] = float(q)
    write_csv(args.out / "tables/pair_summary.csv", pair_rows)
    write_csv(args.out / "tables/occurrence_level_response.csv", occurrence_rows)
    write_csv(args.out / "tables/event_triggered_response_curves.csv", response_curve_rows)
    write_csv(args.out / "tables/time_block_cv_scores.csv", cv_rows_all)
    write_csv(args.out / "tables/event_group_cv_scores.csv", group_rows_all)
    write_csv(args.out / "tables/kernel_coefficients_descriptive.csv", kernel_rows)
    write_csv(args.out / "tables/circular_shift_null.csv", null_rows)
    # Q3-linked subset: keep metric type explicitly descriptive/proxy.
    write_csv(args.out / "tables/q3_registration_error_linkage.csv",
              [{"pair_key": r.get("pair_key"), "cfu_slice_0based": r.get("cfu_slice_0based"),
                "cfu_id": r.get("cfu_id"), "q3_n_frames": r.get("q3_n_frames"),
                "q3_median_normalized_mse": r.get("q3_median_normalized_mse"),
                "q3_p95_normalized_mse": r.get("q3_p95_normalized_mse"),
                "q3_pair_hit_median_normalized_mse": r.get("q3_pair_hit_median_normalized_mse"),
                "q3_other_median_normalized_mse": r.get("q3_other_median_normalized_mse"),
                "interpretation": "local/full-frame MSE ratio proxy from Q3; not ground-truth registration error"}
               for r in completed])
    summary = {"n_candidate_pairs": len(seen), "n_completed_pairs": len(completed),
               "n_pairs_with_event_groups_ge3": sum(r.get("event_group_cv_status") == "evaluated" for r in completed),
               "n_pairs_delta_mse_positive": sum(float(r["delta_mse_M0_minus_M2"]) > 0 for r in completed),
               "n_pairs_circular_shift_q_bh_le_0p05": sum(float(r.get("q_bh_circular_shift", 1)) <= .05 for r in completed),
               "mean_delta_mse_M0_minus_M2": float(np.mean([float(r["delta_mse_M0_minus_M2"]) for r in completed])) if completed else None,
               "elapsed_seconds": time.time() - started,
               "caveat": "Within-recording candidate discovery makes all Q4 inference exploratory; no causal claim."}
    (args.out / "diagnostics/summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
