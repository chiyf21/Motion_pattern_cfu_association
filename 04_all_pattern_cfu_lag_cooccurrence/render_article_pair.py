#!/usr/bin/env python3
"""Article-style rendering for slice08 P200 x slice08 CFU002.

This is a presentation-only renderer. It reuses the current co-occurrence
inputs and AQuA2 CFU fields; it does not recompute event detection or change
the statistical result.
"""
from __future__ import annotations

import csv
import importlib.util
import os
import pickle
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.ticker import MaxNLocator


BASE = Path(__file__).resolve().parents[1]
EXP = BASE / "04_all_pattern_cfu_lag_cooccurrence"
RESULT = EXP / "results/global_shift_empirical_fdr_onset"
CFU_BASE = BASE / "02_current_cfu_input/cfu"
PATTERN_SLICE = int(os.environ.get("PATTERN_SLICE", "8"))
PATTERN_ID = int(os.environ.get("PATTERN_ID", "200"))
CFU_SLICE = int(os.environ.get("CFU_SLICE", "8"))
CFU_ID = int(os.environ.get("CFU_ID", "2"))
OUT = EXP / (
    f"figures/article_pair_slice{PATTERN_SLICE:02d}_P{PATTERN_ID:03d}__"
    f"slice{CFU_SLICE:02d}_CFU{CFU_ID:03d}"
)
T = 1598
PS = 7
REG_MASK_GAP = 5
FRAME_SECONDS = 3.80843294434
TIMELINE_FIGSIZE = (8.4, 3.2)
TIMELINE_XTICKS = np.arange(0, 101, 20)
sys.path.insert(0, str(BASE))
from config import REFERENCE_TIF
REF_TIF = REFERENCE_TIF

spec = importlib.util.spec_from_file_location(
    "fdr", EXP / "run_all_cfu_pattern_lag8_w3_pairwise.py"
)
fdr = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(fdr)


def cmap_from_hue(light: str, dark: str, name: str):
    return mcolors.LinearSegmentedColormap.from_list(name, [light, dark])


PATTERN_CMAP = cmap_from_hue("#F5DCEB", "#8D145D", "pattern_strength")
CFU_CMAP = cmap_from_hue("#D9EEF2", "#087F8C", "cfu_strength")


def load_pattern():
    path = (
        BASE
        / f"experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{PATTERN_SLICE:02d}_velocity_decomp"
        / "06_patterns/objects.pkl"
    )
    with path.open("rb") as f:
        patterns = pickle.load(f)["patterns"]
    return next(p for p in patterns if int(p.pattern_id) == PATTERN_ID)


def load_cfu_fields():
    path = (
        CFU_BASE
        / f"slice_Z{CFU_SLICE:02d}_ds7_native_CFU_ot030_min5_group5.mat"
    )
    with h5py.File(path, "r") as f:
        cells = f["cfuInfo1"]
        layout, n = fdr._detect_cfu_cell_layout(cells)
        if CFU_ID > n:
            raise ValueError(f"CFU{CFU_ID:03d} is not present; n={n}")

        def read_field(field_zero_based: int):
            ref = fdr._get_cfu_field_ref(cells, layout, CFU_ID - 1, field_zero_based)
            return np.asarray(f[ref])

        # AQuA2 fields: 3=weighted spatial region, 4=onset, 6=dF/F.
        return (
            read_field(2).T.astype(float),
            read_field(3).reshape(-1).astype(bool),
            read_field(5).reshape(-1).astype(float),
        )


def load_background(z):
    stack = tifffile.imread(REF_TIF).astype(np.float32)
    bg = stack[z]
    finite = bg[np.isfinite(bg)]
    lo, hi = np.percentile(finite, [1, 99])
    return np.clip((bg - lo) / max(hi - lo, 1e-6), 0, 1) ** 0.65


def cfu_to_reference(mask: np.ndarray, ref_shape: tuple[int, int]) -> np.ndarray:
    """Map the cropped DS7 weighted CFU map into the full reference image."""
    h, w = ref_shape
    yp = (np.arange(h) // PS).astype(int) - REG_MASK_GAP
    xp = (np.arange(w) // PS).astype(int) - REG_MASK_GAP
    valid_y = (yp >= 0) & (yp < mask.shape[0])
    valid_x = (xp >= 0) & (xp < mask.shape[1])
    ysafe = np.clip(yp, 0, mask.shape[0] - 1)
    xsafe = np.clip(xp, 0, mask.shape[1] - 1)
    out = mask[np.ix_(ysafe, xsafe)].copy()
    out[~valid_y, :] = 0
    out[:, ~valid_x] = 0
    return out


def robust_limits(values: np.ndarray, mask: np.ndarray | None = None):
    x = np.asarray(values, dtype=float)
    if mask is not None:
        x = x[np.asarray(mask, dtype=bool)]
    x = x[np.isfinite(x)]
    x = x[x > 0]
    if not len(x):
        return 0.0, 1.0
    lo, hi = np.percentile(x, [1, 99])
    if hi <= lo:
        lo, hi = float(np.min(x)), float(np.max(x))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def save_both(fig, stem: str, tight: bool = True):
    bbox = "tight" if tight else None
    fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white", bbox_inches=bbox)
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white", bbox_inches=bbox)
    plt.close(fig)


def style_spatial_axis(ax, bg):
    ax.imshow(bg, origin="lower", cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xlim(0, bg.shape[1])
    # Keep the established AQuA2 display convention: y increases downward.
    ax.set_ylim(bg.shape[0], 0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_spatial(pattern, cfu_weight, bg_pattern, bg_cfu):
    pmask = np.asarray(pattern.unified_mask, dtype=bool)
    response = np.asarray(pattern.unified_response_field, dtype=float)
    strength = np.linalg.norm(response, axis=-1)
    pattern_strength = np.where(pmask, strength, 0.0)
    pattern_lo, pattern_hi = robust_limits(pattern_strength, pmask)

    cfu_full = cfu_to_reference(cfu_weight, bg_cfu.shape)
    cfu_mask = cfu_full > 0
    cfu_lo, cfu_hi = robust_limits(cfu_full, cfu_mask)

    # Pattern is defined on the DS7 grid. Expand its per-patch strength to
    # native reference pixels, matching the existing PS=7 coordinate map.
    pattern_full = np.repeat(np.repeat(pattern_strength, PS, axis=0), PS, axis=1)
    pattern_full = pattern_full[: bg_pattern.shape[0], : bg_pattern.shape[1]]
    pattern_mask_full = np.repeat(np.repeat(pmask, PS, axis=0), PS, axis=1)
    pattern_mask_full = pattern_mask_full[: bg_pattern.shape[0], : bg_pattern.shape[1]]

    pnorm = mcolors.Normalize(pattern_lo, pattern_hi)
    pimg = np.ma.masked_where(~pattern_mask_full, pattern_full)
    cnorm = mcolors.Normalize(cfu_lo, cfu_hi)
    cimg = np.ma.masked_where(~cfu_mask, cfu_full)

    # Save the two spatial panels independently so they can be composed in
    # the manuscript without a built-in legend/colorbar.
    figp, axp = plt.subplots(figsize=(6.1, 5.1), constrained_layout=True)
    style_spatial_axis(axp, bg_pattern)
    axp.imshow(pimg, origin="lower", cmap=PATTERN_CMAP, norm=pnorm,
               alpha=0.90, interpolation="nearest")
    axp.contour(pattern_mask_full, levels=[0.5], colors=["#6D1049"], linewidths=0.8)
    save_both(figp, f"pattern_P{PATTERN_ID:03d}_spatial_strength")

    figc, axc = plt.subplots(figsize=(6.1, 5.1), constrained_layout=True)
    style_spatial_axis(axc, bg_cfu)
    axc.imshow(cimg, origin="lower", cmap=CFU_CMAP, norm=cnorm,
               alpha=0.90, interpolation="nearest")
    axc.contour(cfu_mask, levels=[0.5], colors=["#05616B"], linewidths=0.8)
    save_both(figc, f"CFU{CFU_ID:03d}_spatial_strength")

    return {
        "pattern_strength_limits": (pattern_lo, pattern_hi),
        "cfu_weight_limits": (cfu_lo, cfu_hi),
        "pattern_mask_shape": pmask.shape,
        "cfu_mask_shape": cfu_weight.shape,
    }


def pattern_activation_sequence(pattern):
    """Return the binary union of all member activation durations."""
    sequence = np.zeros(T, dtype=np.uint8)
    components = getattr(pattern, "components", None) or getattr(pattern, "regions", [])
    for component in components:
        frames, _ = fdr.component_frames(component)
        if len(frames):
            sequence[frames] = 1
    return sequence


def render_pattern_timeline(pattern, peaks):
    peaks = np.asarray(peaks, dtype=int)
    peaks = peaks[(peaks >= 0) & (peaks < T)]
    sequence = pattern_activation_sequence(pattern)
    frames = np.arange(T)
    x = frames * FRAME_SECONDS / 60.0
    y = sequence.astype(float)
    xmax = (T - 1) * FRAME_SECONDS / 60.0
    # Fixed canvas and axes geometry are shared with the CFU timeline below,
    # so the two PDFs can be stacked without horizontal misalignment.
    fig = plt.figure(figsize=TIMELINE_FIGSIZE, facecolor="white")
    ax = fig.add_axes([0.24, 0.30, 0.70, 0.66])
    ax.step(x, y, where="post", color="#8D145D", linewidth=2.5)
    ax.fill_between(x, 0, y, step="post", color="#C979A7", alpha=0.42, linewidth=0)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.08, 1.55)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"])
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("activation")
    ax.set_xticks(TIMELINE_XTICKS)
    ax.tick_params(axis="both", labelsize=24, width=1.5, length=6)
    ax.xaxis.label.set_size(28)
    ax.yaxis.label.set_size(28)
    ax.grid(axis="x", color="0.84", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_both(fig, f"pattern_P{PATTERN_ID:03d}_activation_timeline", tight=False)
    return peaks, sequence


def render_cfu_dff(dff, onsets):
    dff = np.asarray(dff, dtype=float)
    dff = dff[:T]
    frames = np.arange(len(dff))
    minutes = frames * FRAME_SECONDS / 60.0
    xmax = (T - 1) * FRAME_SECONDS / 60.0
    onsets = np.asarray(onsets, dtype=int)
    onsets = onsets[(onsets >= 0) & (onsets < len(dff))]
    fig = plt.figure(figsize=TIMELINE_FIGSIZE, facecolor="white")
    ax = fig.add_axes([0.24, 0.30, 0.70, 0.66])
    ax.plot(minutes, dff, color="#087F8C", linewidth=2.3, rasterized=True)
    ax.axhline(0, color="0.45", linewidth=1.2)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("")
    ax.set_ylabel("dF/F")
    ax.set_xticks(TIMELINE_XTICKS)
    ax.tick_params(axis="both", labelsize=24, width=1.5, length=6)
    ax.tick_params(axis="x", labelbottom=False, bottom=False)
    ax.xaxis.label.set_size(28)
    ax.yaxis.label.set_size(28)
    ax.grid(axis="x", color="0.88", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_both(fig, f"CFU{CFU_ID:03d}_dff_timeline", tight=False)
    return onsets


def read_pair_row():
    path = RESULT / "empirical_FDR_pairs_q001.csv"
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if (
                int(row["motion_slice_0based"]) == PATTERN_SLICE
                and int(row["pattern_id"]) == PATTERN_ID
                and int(row["cfu_slice_0based"]) == CFU_SLICE
                and int(row["cfu_id"]) == CFU_ID
            ):
                return row
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pattern = load_pattern()
    cfu_weight, onsets, dff = load_cfu_fields()
    bg_pattern = load_background(PATTERN_SLICE)
    bg_cfu = load_background(CFU_SLICE)
    spatial_meta = render_spatial(pattern, cfu_weight, bg_pattern, bg_cfu)
    peaks = fdr.motion_event_peaks(PATTERN_SLICE)
    peak_map = {pid: p for pid, _, p in peaks}
    # A user may request a valid spatial pattern that was excluded from the
    # lag test (e.g. too few detected motion peaks).  Render it faithfully as
    # an all-zero activation sequence instead of failing or implying a test.
    pattern_peaks, pattern_sequence = render_pattern_timeline(
        pattern, peak_map.get(PATTERN_ID, np.array([], dtype=int))
    )
    cfu_onsets = render_cfu_dff(dff, np.flatnonzero(onsets))

    row = read_pair_row()
    metadata = {
        "pattern_slice_0based": PATTERN_SLICE,
        "pattern_id": PATTERN_ID,
        "cfu_slice_0based": CFU_SLICE,
        "cfu_id": CFU_ID,
        "pattern_activation_peaks_n": len(pattern_peaks),
        "pattern_active_frames_n": int(pattern_sequence.sum()),
        "pattern_entered_lag_test": PATTERN_ID in peak_map,
        "cfu_onsets_n": len(cfu_onsets),
        "frame_seconds": FRAME_SECONDS,
        "pattern_min_max_display": spatial_meta["pattern_strength_limits"],
        "cfu_weight_min_max_display": spatial_meta["cfu_weight_limits"],
    }
    if row:
        metadata.update({f"pair_{k}": v for k, v in row.items()})
    with (OUT / "metadata.txt").open("w") as f:
        for k, v in metadata.items():
            f.write(f"{k}: {v}\n")
    (OUT / "README.md").write_text(
        f"Article-style rendering for slice{PATTERN_SLICE:02d} P{PATTERN_ID:03d} x "
        f"slice{CFU_SLICE:02d} CFU{CFU_ID:03d}.\n\n"
        "The pattern and CFU spatial maps are saved as separate figures without colorbars. "
        "Within each same-hue map, brightness encodes the pattern unified_response_field vector norm "
        "or the AQuA2 cfuInfo1 field 3 weighted spatial map. "
        "The pattern timeline is a binary 0/1 union of all member activation durations; dots mark "
        "the same motion peaks used by the current lag co-occurrence run. "
        "The CFU curve is AQuA2 cfuInfo1 field 6 (native dF/F); onset field 4 is retained in metadata but not drawn. "
        "Coordinates use the current PS=7 and regMaskGap=5 mapping, with the established inverted-y display convention.\n"
    )
    print(f"output: {OUT}")
    print(f"pattern peaks: {len(pattern_peaks)}; CFU onsets: {len(cfu_onsets)}")


if __name__ == "__main__":
    main()
