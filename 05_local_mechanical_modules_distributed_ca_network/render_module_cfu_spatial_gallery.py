#!/usr/bin/env python3
"""Spatial gallery: Fig5/24 modules and all q<0.05 associated CFUs.

Each page is module-centric.  The first panel shows the target local module
(pattern response magnitude + its local CFU spatial weight) on the target
slice. Remaining panels show every non-local significant CFU on its native
slice reference. Direction/lag is displayed only as metadata, not a filter.

This is an exploratory gallery, so it writes PNG only for fast, restartable
batch rendering. Article-quality PDFs can be produced later for selected
modules.
"""
from __future__ import annotations

import csv
import pickle
import sys
from collections import defaultdict
from functools import lru_cache
from math import ceil
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import tifffile


BASE = Path(__file__).resolve().parents[1]
EXP23 = BASE / "04_all_pattern_cfu_lag_cooccurrence"
EXP24 = BASE / "05_local_mechanical_modules_distributed_ca_network"
MODULES = EXP24 / "02_module_cfu_network/current_module_table.csv"
EDGES = EXP24 / "02_module_cfu_network/module_cfu_associations_fdr_q005.csv"
CFU_LINKS = BASE / "02_current_cfu_input/cfu"
OUT = EXP24 / "03_module_cfu_spatial_gallery_q005"
sys.path.insert(0, str(BASE))
from config import REFERENCE_TIF
REF_TIF = REFERENCE_TIF

T = 1598
PS = 7
REG_MASK_GAP = 5
PATTERN_SHAPE = (244, 329)
MASK_THRESHOLD = 0.1

sys.path.insert(0, str(BASE))

PATTERN_CMAP = mcolors.LinearSegmentedColormap.from_list("pattern", ["#F5DCEB", "#8D145D"])
CFU_CMAP = mcolors.LinearSegmentedColormap.from_list("cfu", ["#D9EEF2", "#087F8C"])


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def field_ref(cells, field0: int, cfu0: int):
    if cells.shape[0] == 9:
        return cells[field0, cfu0]
    if cells.shape[1] == 9:
        return cells[cfu0, field0]
    raise ValueError(f"Unknown cfuInfo1 layout {cells.shape}")


@lru_cache(maxsize=None)
def cfu_weight(z: int, cid: int):
    path = CFU_LINKS / f"slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat"
    with h5py.File(path, "r") as f:
        raw = np.asarray(f[field_ref(f["cfuInfo1"], 2, cid - 1)])
    return raw.T.astype(float)


@lru_cache(maxsize=None)
def pattern_object(z: int, pid: int):
    path = (
        BASE
        / f"experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{z:02d}_velocity_decomp"
        / "06_patterns/objects.pkl"
    )
    with path.open("rb") as f:
        patterns = {int(p.pattern_id): p for p in pickle.load(f)["patterns"]}
    return patterns[pid]


def load_backgrounds():
    stack = tifffile.imread(REF_TIF).astype(np.float32)
    out = {}
    for z in range(stack.shape[0]):
        bg = stack[z]
        lo, hi = np.percentile(bg[np.isfinite(bg)], [1, 99])
        out[z] = np.clip((bg - lo) / max(hi - lo, 1e-6), 0, 1) ** 0.65
    return out


def cfu_to_reference(weight: np.ndarray, shape: tuple[int, int]):
    h, w = shape
    yp = (np.arange(h) // PS).astype(int) - REG_MASK_GAP
    xp = (np.arange(w) // PS).astype(int) - REG_MASK_GAP
    valid_y = (yp >= 0) & (yp < weight.shape[0])
    valid_x = (xp >= 0) & (xp < weight.shape[1])
    ysafe = np.clip(yp, 0, weight.shape[0] - 1)
    xsafe = np.clip(xp, 0, weight.shape[1] - 1)
    out = weight[np.ix_(ysafe, xsafe)].copy()
    out[~valid_y, :] = 0
    out[:, ~valid_x] = 0
    return out


def robust_norm(values: np.ndarray, valid: np.ndarray):
    x = np.asarray(values, float)[valid]
    x = x[np.isfinite(x) & (x > 0)]
    if not len(x):
        return mcolors.Normalize(0, 1)
    lo, hi = np.percentile(x, [1, 99])
    if hi <= lo:
        hi = lo + 1e-6
    return mcolors.Normalize(lo, hi)


def setup(ax, bg):
    ax.imshow(bg, origin="lower", cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xlim(0, bg.shape[1])
    ax.set_ylim(bg.shape[0], 0)  # same established AQuA2 display convention
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def add_cfu(ax, bg, z, cid):
    full = cfu_to_reference(cfu_weight(z, cid), bg.shape)
    valid = full > MASK_THRESHOLD
    ax.imshow(
        np.ma.masked_where(~valid, full), origin="lower", cmap=CFU_CMAP,
        norm=robust_norm(full, valid), alpha=.88, interpolation="nearest",
    )
    ax.contour(valid, levels=[.5], colors=["#05616B"], linewidths=.65)


def add_module(ax, bg, z, pid, local_cid):
    p = pattern_object(z, pid)
    pmask = np.asarray(p.unified_mask, bool)
    response = np.asarray(p.unified_response_field, float)
    strength = np.linalg.norm(response, axis=-1)
    full_strength = np.repeat(np.repeat(np.where(pmask, strength, 0), PS, axis=0), PS, axis=1)
    full_mask = np.repeat(np.repeat(pmask, PS, axis=0), PS, axis=1)
    full_strength = full_strength[:bg.shape[0], :bg.shape[1]]
    full_mask = full_mask[:bg.shape[0], :bg.shape[1]]
    ax.imshow(
        np.ma.masked_where(~full_mask, full_strength), origin="lower", cmap=PATTERN_CMAP,
        norm=robust_norm(full_strength, full_mask), alpha=.74, interpolation="nearest",
    )
    ax.contour(full_mask, levels=[.5], colors=["#6D1049"], linewidths=.65)
    add_cfu(ax, bg, z, local_cid)


def fmt_q(q: float):
    return f"q={q:.2g}" if q >= 1e-3 else f"q={q:.1e}"


def render_one(module, edges, backgrounds):
    module_id = module["module_id"]
    z = int(module["module_slice_0based"])
    pid = int(module["module_pattern_id"])
    local_cid = int(module["local_cfu_id"])
    # The local CFU is shown in the target module panel, not duplicated.
    sources = [
        e for e in sorted(edges, key=lambda x: (
            float(x["q_empirical_global"]),
            int(x["source_cfu_slice_0based"]),
            int(x["source_cfu_id"]),
        ))
        if e["source_scope"] != "intramodule_local_cfu"
    ]
    panels = 1 + len(sources)
    ncol = min(4, panels)
    nrow = ceil(panels / ncol)
    stem = f"{module_id}_q005_CFU_locations"
    png_path = OUT / "figures" / f"{stem}.png"
    if png_path.exists():
        return len(sources), stem
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.85 * ncol, 3.25 * nrow), facecolor="white")
    axs = np.atleast_1d(axs).ravel()

    setup(axs[0], backgrounds[z])
    add_module(axs[0], backgrounds[z], z, pid, local_cid)
    axs[0].set_title(
        f"Target module\nS{z:02d} P{pid:03d} + local C{local_cid:03d}",
        fontsize=10, pad=5,
    )
    for ax, edge in zip(axs[1:], sources):
        source_z = int(edge["source_cfu_slice_0based"])
        source_cid = int(edge["source_cfu_id"])
        setup(ax, backgrounds[source_z])
        add_cfu(ax, backgrounds[source_z], source_z, source_cid)
        ax.set_title(
            f"Associated CFU\nS{source_z:02d} C{source_cid:03d} | {fmt_q(float(edge['q_empirical_global']))}, lag={int(edge['best_lag'])}",
            fontsize=9, pad=5,
        )
    for ax in axs[panels:]:
        ax.axis("off")
    fig.suptitle(f"{module_id} | q<0.05 associated CFUs: {len(edges)}", fontsize=13, y=.995)
    fig.tight_layout(rect=[0, 0, 1, .97])
    fig.savefig(png_path, dpi=180, facecolor="white")
    plt.close(fig)
    return len(sources), stem


def main():
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    modules = read_csv(MODULES)
    by_module = defaultdict(list)
    for edge in read_csv(EDGES):
        by_module[edge["module_id"]].append(edge)
    modules = [m for m in modules if m["module_id"] in by_module]
    modules.sort(key=lambda m: (int(m["module_slice_0based"]), int(m["module_pattern_id"]), int(m["local_cfu_id"])))
    backgrounds = load_backgrounds()
    manifest = []
    for i, module in enumerate(modules, 1):
        n_nonlocal, stem = render_one(module, by_module[module["module_id"]], backgrounds)
        manifest.append({
            "module_id": module["module_id"],
            "module_slice_0based": module["module_slice_0based"],
            "module_pattern_id": module["module_pattern_id"],
            "local_cfu_id": module["local_cfu_id"],
            "n_q005_associated_cfus": len(by_module[module["module_id"]]),
            "n_q005_nonlocal_cfus_shown": n_nonlocal,
            "figure_png": f"figures/{stem}.png",
            "figure_pdf": "",
        })
        if i == 1 or i % 10 == 0 or i == len(modules):
            print(f"rendered {i}/{len(modules)}", flush=True)
    fields = list(manifest[0]) if manifest else ["module_id"]
    with (OUT / "module_gallery_manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest)
    (OUT / "README.md").write_text(
        "Each module-centric PNG shows the target current spatial module (pattern response magnitude in magenta and local CFU spatial weight in teal) plus every empirical-FDR q<0.05 CFU associated with the module pattern. Associated CFUs are plotted on their own native slice reference. No lag-direction filter is applied; lag is displayed only for interpretation. The source tables are Fig5/24/02_module_cfu_network/module_cfu_associations_fdr_q005.csv.\n"
    )
    print(f"output: {OUT}; module figures: {len(manifest)}")


if __name__ == "__main__":
    main()
