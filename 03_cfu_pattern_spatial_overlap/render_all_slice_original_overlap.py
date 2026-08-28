#!/usr/bin/env python3
"""Render one native-resolution spatial-overlap image for every v0827 slice.

Each output keeps the registered reference image's original pixel dimensions.
Magenta marks the union of selected pattern masks, cyan marks the union of
their paired CFU masks, and yellow marks their spatial intersection.  These
are display overlays only: membership and selection are read exclusively from
the v0827 spatial-final CSV tables.
"""
from __future__ import annotations

import csv
import pickle
import sys
from pathlib import Path

import h5py
import numpy as np
import tifffile
from PIL import Image

BASE = Path("/home/cyf/wbi/wbi_code")
sys.path.insert(0, str(BASE))  # required to unpickle MotionPattern objects
SPATIAL = BASE / "experiments/Fig5_v0827/03_cfu_pattern_spatial_overlap"
PATTERN_ROOT = BASE / "experiments/Fig5_v0827/01_motion_pattern_extraction_omega05_mu05/patterns"
CFU_BASE = BASE / "experiments/Fig5_v0827/02_current_cfu_input/cfu"
OUT = SPATIAL / "figures/original_resolution_pattern_cfu_overlap"
REFERENCE = Path("/mnt/data21T_2/cyf/f338/f338_registrated_0530/reference/vol_ref_000599_000999.tif")

PATCH = 7
GAP = 5
PATCH_SHAPE = (244, 329)
CFU_WEIGHT_THRESHOLD = 0.1


def cfu_ref(cells, field0: int, cfu0: int):
    if cells.shape[0] == 9:
        return cells[field0, cfu0]
    if cells.shape[1] == 9:
        return cells[cfu0, field0]
    raise ValueError(f"Unexpected cfuInfo1 layout: {cells.shape}")


def full_cfu_mask(z: int, cid: int) -> np.ndarray:
    path = CFU_BASE / f"slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat"
    with h5py.File(path, "r") as f:
        raw = np.asarray(f[cfu_ref(f["cfuInfo1"], 2, cid - 1)]).T
    if raw.shape != (234, 319):
        raise ValueError(f"Z{z:02d} CFU{cid:03d}: unexpected cropped shape {raw.shape}")
    out = np.zeros(PATCH_SHAPE, dtype=bool)
    out[GAP:-GAP, GAP:-GAP] = raw > CFU_WEIGHT_THRESHOLD
    return out


def full_pattern_mask(z: int, pid: int) -> np.ndarray:
    path = PATTERN_ROOT / f"Slice{z:02d}_velocity_decomp/06_patterns/objects.pkl"
    with path.open("rb") as f:
        patterns = pickle.load(f)["patterns"]
    pattern = next(p for p in patterns if int(p.pattern_id) == pid)
    mask = np.asarray(pattern.unified_mask, dtype=bool)
    if mask.shape != PATCH_SHAPE:
        raise ValueError(f"Z{z:02d} P{pid:03d}: unexpected pattern shape {mask.shape}")
    return mask


def normalized_background(reference: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(reference, [1, 99])
    x = np.clip((reference.astype(float) - lo) / (hi - lo + 1e-12), 0, 1) ** 0.65
    return np.round(x * 255).astype(np.uint8)


def blend(base: np.ndarray, mask: np.ndarray, color: tuple[int, int, int], alpha: float) -> np.ndarray:
    out = base.astype(float).copy()
    color_arr = np.asarray(color, float)
    out[mask] = (1 - alpha) * out[mask] + alpha * color_arr
    return np.round(np.clip(out, 0, 255)).astype(np.uint8)


def fit_to_reference(mask: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Zero-pad the DS7 grid mapping at the bottom/right reference border."""
    out = np.zeros(shape, dtype=bool)
    h = min(shape[0], mask.shape[0])
    w = min(shape[1], mask.shape[1])
    out[:h, :w] = mask[:h, :w]
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stack = tifffile.imread(REFERENCE)
    manifest = []
    for z in range(1, 13):
        rows_path = SPATIAL / f"spatial_final_slice{z}.csv"
        with rows_path.open(newline="") as f:
            rows = list(csv.DictReader(f))
        pattern_union = np.zeros(PATCH_SHAPE, dtype=bool)
        cfu_union = np.zeros(PATCH_SHAPE, dtype=bool)
        for row in rows:
            pattern_union |= full_pattern_mask(z, int(row["pattern_id"]))
            cfu_union |= full_cfu_mask(z, int(row["cfu_id"]))
        overlap = pattern_union & cfu_union
        background = normalized_background(stack[z])
        rgb = np.repeat(background[..., None], 3, axis=2)
        pfull = np.repeat(np.repeat(pattern_union, PATCH, axis=0), PATCH, axis=1)
        cfull = np.repeat(np.repeat(cfu_union, PATCH, axis=0), PATCH, axis=1)
        ofull = np.repeat(np.repeat(overlap, PATCH, axis=0), PATCH, axis=1)
        h, w = rgb.shape[:2]
        pfull = fit_to_reference(pfull, (h, w))
        cfull = fit_to_reference(cfull, (h, w))
        ofull = fit_to_reference(ofull, (h, w))
        rgb = blend(rgb, pfull & ~ofull, (218, 57, 128), 0.55)
        rgb = blend(rgb, cfull & ~ofull, (15, 164, 180), 0.55)
        rgb = blend(rgb, ofull, (249, 202, 47), 0.78)
        path = OUT / f"slice{z:02d}_pattern_cfu_overlap_original.png"
        Image.fromarray(rgb).save(path, compress_level=1)
        manifest.append({
            "slice_0based": z,
            "n_modules": len(rows),
            "pattern_union_patch_pixels": int(pattern_union.sum()),
            "cfu_union_patch_pixels": int(cfu_union.sum()),
            "overlap_union_patch_pixels": int(overlap.sum()),
            "image": str(path),
        })
        print(f"slice{z:02d}: modules={len(rows)} -> {path}", flush=True)
    with (OUT / "manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    (OUT / "README.md").write_text(
        "# v0827 native-resolution pattern-CFU overlap\n\n"
        "One PNG per 0-based slice on the original registered reference grid. "
        "Magenta=union of selected pattern masks; cyan=union of selected CFU masks; "
        "yellow=their union-level intersection. Pattern and CFU masks are mapped from "
        "the shared DS7 244x329 grid using patch size 7 and CFU regMaskGap=5.\n"
    )


if __name__ == "__main__":
    main()
