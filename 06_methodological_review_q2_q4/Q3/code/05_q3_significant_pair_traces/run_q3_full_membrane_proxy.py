#!/usr/bin/env python3
"""Compute a framewise local membrane-reference appearance proxy for Q3.

This is intentionally not named registration error: without independent
landmarks/ground truth, image mismatch is only an appearance residual proxy.
The script reads each registered membrane TIFF once and evaluates all CFU masks
on the slices represented in the current q<0.05 pair table.
"""
from __future__ import annotations

import argparse
import gc
import os
import sys
from pathlib import Path

# This point-wise image metric does not benefit from BLAS fan-out.
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"

import numpy as np
import tifffile

EXP_ROOT = Path(__file__).resolve().parents[1]
Q3_DIR = EXP_ROOT / "03_q3_local_membrane_mse"
sys.path.insert(0, str(EXP_ROOT))
sys.path.insert(0, str(Q3_DIR))

from review_common import (
    REPO_ROOT, T, read_csv, read_cfu, write_csv, write_json,
)
from run_q3_membrane_qc import (
    DEFAULT_REGISTRATION, GAP, T as Q3_T, downsample7,
    patch_grid, reference_by_frame,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results/q3_membrane_qc_reference_matched_full_proxy_mse")
    ap.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    ap.add_argument("--cfu-root", type=Path, default=None,
                    help="CFU MAT directory; defaults to <repo-root>/02_current_cfu_input/cfu.")
    ap.add_argument("--registration-root", type=Path, default=DEFAULT_REGISTRATION)
    ap.add_argument("--reference-dir", type=Path, default=DEFAULT_REGISTRATION / "reference")
    args = ap.parse_args()
    cfu_root = args.cfu_root or (args.repo_root / "02_current_cfu_input" / "cfu")
    if Q3_T != T:
        raise RuntimeError(f"Inconsistent frame counts: {T} vs {Q3_T}")
    args.out.mkdir(parents=True, exist_ok=True)

    sig_path = args.repo_root / "04_all_pattern_cfu_lag_cooccurrence/results/global_shift_empirical_fdr_onset/empirical_FDR_significant_pairs.csv"
    sig = [r for r in read_csv(sig_path) if float(r["q_empirical_global"]) < 0.05]
    cfus_by_slice = {}
    for z in sorted({int(r["cfu_slice_0based"]) for r in sig}):
        ids = [int(r["cfu_id"]) for r in sig if int(r["cfu_slice_0based"]) == z]
        cfus_by_slice[z] = read_cfu(z, ids, cfu_root)
    masks = {
        (z, cid): patch_grid(c["mask"])
        for z, cfus in cfus_by_slice.items() for cid, c in cfus.items()
    }
    refs = reference_by_frame(args.reference_dir, T)
    slices = sorted(cfus_by_slice)
    ref_cache = {}
    for path in sorted(set(refs.values())):
        stack = tifffile.memmap(path)
        for z in slices:
            ref = downsample7(np.asarray(stack[z], dtype=np.float32))
            ref_cache[(path, z)] = ref
        del stack

    rows = []
    for frame in range(T):
        path = args.registration_root / f"membrane/vol_ch1_{frame:06d}.tif"
        mapped = tifffile.memmap(path)
        ref_path = refs[frame]
        try:
            for z in slices:
                image = downsample7(np.asarray(mapped[z], dtype=np.float32))
                ref = ref_cache[(ref_path, z)]
                for cid, c in cfus_by_slice[z].items():
                    mask = masks[(z, cid)]
                    n = int(mask.sum())
                    rv = ref[mask].astype(np.float64) if n else np.empty(0)
                    iv = image[mask].astype(np.float64) if n else np.empty(0)
                    mse = float(np.mean(np.square(iv - rv))) if n else np.nan
                    ref_variance = float(np.var(rv)) if n else np.nan
                    normalized_mse = mse / ref_variance if n and ref_variance > 1e-12 else np.nan
                    rows.append({
                        "slice_id": z, "cfu_id": cid, "frame_0based": frame,
                        "reference_file": ref_path.name, "n_mask_pixels": n,
                        "masked_mse": mse, "reference_variance": ref_variance,
                        "normalized_mse": normalized_mse,
                    })
        finally:
            del mapped
        if frame == 0 or (frame + 1) % 100 == 0 or frame + 1 == T:
            print(f"Framewise local membrane MSE: {frame + 1}/{T}", flush=True)

    out_csv = args.out / "q3_local_proxy_timeseries.csv"
    write_csv(out_csv, rows)
    write_json(args.out / "run_config.json", {
        "experiment": "Framewise local membrane-reference MSE for current significant-pair CFUs",
        "repo_root": args.repo_root, "cfu_root": cfu_root,
        "registration_root": args.registration_root,
        "reference_dir": args.reference_dir, "frame_count": T,
        "slices": slices, "n_unique_cfus": sum(len(x) for x in cfus_by_slice.values()),
        "n_rows": len(rows), "patch_size": 7, "reg_mask_gap": GAP,
        "reference_policy": "Saved interval containing frame; later-start interval wins on a shared boundary.",
        "metrics": ["masked_mse", "normalized_mse"],
        "normalized_mse_definition": "masked_mse / spatial variance of matching reference pixels inside the CFU mask; NaN if variance <= 1e-12.",
        "interpretation": "Dense framewise local image mismatch proxy; not a ground-truth registration error.",
    })
    print(f"Saved {len(rows)} local time-series rows to {out_csv}", flush=True)
    gc.collect()


if __name__ == "__main__":
    main()
