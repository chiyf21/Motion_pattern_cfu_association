#!/usr/bin/env python3
"""Canonical Fig5 motion pipeline: registered motion H5 -> MotionPatterns.

Each invocation writes a resumable cache for one 0-based slice:
01_patch_motion, 02_motion_units, 03_episodes, 04_modes and 06_patterns.
The implementation imports only the pinned Fig5/wholistic_registration clone.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy import ndimage

CLONE = Path('/home/cyf/wbi/wbi_code/experiments/Fig5/wholistic_registration')
sys.path.insert(0, str(CLONE / 'src'))
from wholistic_registration.utils import motion_correlation_pattern as mcp
from wholistic_registration.utils import motion_stage_cache as cache

CANONICAL = {
    'patch_size': 7, 'median_filter_size': (1, 3, 3, 1),
    'rest_window_t': 21, 'rest_window_xy': 3,
    'episode': dict(tolerant_time=1, min_total_area=30, expand_frames=1,
                    min_cc_area=8, global_motion_mode='median'),
    'artifact': dict(max_fov_fraction=.5, min_duration=3, max_duration=None,
                     max_global_corr=.90, max_edge_fraction=.80, edge_width=3),
    # This is the sparse-compact mode version used by the current v0827 results.
    'mode': dict(Kmax=8, K_selection_method='svd', svd_target_r2=.90,
                 lambda_sc=.05, rho=1., kappa=4., support_rel_thresh=.08,
                 max_iter=200, use_velocity=True),
    'pattern': dict(min_strength=0., min_area=5, min_duration=1, min_iou=.08,
                    omega=.5, mu=.5, b_distance='correlation', spatial_rule='iou',
                    cluster_dist_thresh=.45, linkage_method='complete',
                    incompatible_dist=1e6, compute_unified=True,
                    unified_mask_mode='best_cc', unified_sign_method='correlation',
                    min_pattern_members=2, min_unified_area=50),
}

def motion_files(raw_dir: Path):
    fs = list(raw_dir.glob('motion_*.h5'))
    return sorted(fs, key=lambda p: int(re.search(r'motion_(\d+)\.h5$', p.name).group(1)))

def stage1(raw_dir: Path, mask: np.ndarray, slice_id: int, out: Path, overwrite: bool):
    target = out / '01_patch_motion' / 'arrays.npz'
    if target.exists() and not overwrite: return
    fs = motion_files(raw_dir)
    if not fs: raise FileNotFoundError(f'No motion_*.h5 under {raw_dir}')
    H, W = mask.shape; ps = CANONICAL['patch_size']; nr, nc = H // ps, W // ps
    abs_motion = np.empty((len(fs), nr, nc, 2), np.float32)
    with h5py.File(fs[0], 'r') as f: key = next(iter(f.keys()))
    for t, p in enumerate(fs):
        with h5py.File(p, 'r') as f:
            frame = np.asarray(f[key][slice_id, :nr*ps, :nc*ps, :2], np.float32)
        abs_motion[t] = frame.reshape(nr, ps, nc, ps, 2).mean(axis=(1, 3))
    delta, aligned = abs_motion[1:] - abs_motion[:-1], abs_motion[1:]
    # The original current cache used a 1x3x3x1 median filter before detection.
    delta = ndimage.median_filter(delta, size=CANONICAL['median_filter_size'])
    aligned = ndimage.median_filter(aligned, size=CANONICAL['median_filter_size'])
    mp = mask[:nr*ps, :nc*ps].reshape(nr, ps, nc, ps).mean(axis=(1,3)) > .3
    cache.save_patch_motion_stage(out, delta, aligned, mp, params={
        'slice_index_0based': slice_id, 'raw_motion_dir': str(raw_dir),
        'n_input_frames': len(fs), 'patch_size': ps,
        'median_filter_size': CANONICAL['median_filter_size']})

def stage2(out: Path, overwrite: bool):
    if (out/'02_motion_units/objects.pkl').exists() and not overwrite: return
    a = cache.load_patch_motion_stage(out); d = a['motion_delta']; valid = ~a['mask_patched']
    mag = np.linalg.norm(d, axis=-1)
    rest = mcp.estimate_rest_state_motion(mag, window_size_t=21, window_size_xy=3, use_gpu='auto')
    units, active = mcp.getMotionUnit(d, mag * valid[None], rest, save_motion=False, use_gpu='auto')
    filtered = mcp.filterMotionUnits(units, active, use_gpu='auto')
    cache.save_motion_units_stage(out, rest, units, active, filtered, params={
        'event_detection_field': 'motion_delta (velocity)', 'rest_window_t':21,
        'rest_window_xy':3, 'use_abs_dev': True})

def stage3(out: Path, overwrite: bool):
    if (out/'03_episodes/objects.pkl').exists() and not overwrite: return
    a = cache.load_patch_motion_stage(out); u = cache.load_motion_units_stage(out)
    d, valid = a['motion_delta'], ~a['mask_patched']
    raw = mcp.getMotionEpisode(u['motion_units_filtered'], d, motion_full_abs=d,
        global_valid_mask=valid, **CANONICAL['episode'])
    eps = mcp.filter_episodes_artifacts(raw, valid_mask=valid, verbose=False, **CANONICAL['artifact'])
    cache.save_episodes_stage(out, eps, params={'episode':CANONICAL['episode'],
        'artifact_filter':CANONICAL['artifact'], 'n_episodes_raw':len(raw),
        'n_episodes_filtered':len(eps)})

def stage4(out: Path, overwrite: bool):
    if (out/'04_modes/objects.pkl').exists() and not overwrite: return
    eps = cache.load_episodes_stage(out)
    modes = mcp.getMotionModes(eps, verbose=False, **CANONICAL['mode'])
    cache.save_modes_stage(out, eps, modes, params=CANONICAL['mode'])

def stage5(out: Path, overwrite: bool):
    if (out/'06_patterns/objects.pkl').exists() and not overwrite: return
    data = cache.load_modes_stage(out); eps, modes = data['episodes'], data['modes']
    by_id = {getattr(e, 'episode_id', -1): e for e in eps}
    for e in eps: e.modes = []
    for m in modes:
        if getattr(m, 'episode_id', -1) in by_id: by_id[m.episode_id].modes.append(m)
    patterns, kept, groups, labels, info = mcp.getMotionPattern(eps, unit_type='mode', verbose=False,
                                                                  **CANONICAL['pattern'])
    cache.save_patterns_stage(out, patterns, kept, groups, labels, info,
                              params=CANONICAL['pattern'])

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--raw-motion-dir', type=Path)
    p.add_argument('--patch-motion-cache', type=Path,
                   help='Existing 01_patch_motion directory containing arrays.npz and metadata.json; skips raw H5 input.')
    p.add_argument('--invalid-mask-npy', type=Path,
                   help='Boolean invalid/background mask in the raw HxW geometry (legacy-compatible convention).')
    p.add_argument('--slice', type=int, required=True, help='0-based H5 slice index.')
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--through-stage', choices=['01','02','03','04','06'], default='06')
    p.add_argument('--overwrite', action='store_true')
    args = p.parse_args()
    out = args.output_root / f'Slice{args.slice:02d}_velocity_decomp'; out.mkdir(parents=True, exist_ok=True)
    cache.write_cache_metadata(out, project_name='Fig5_canonical_motion_pipeline', extra={
        'slice_index_0based':args.slice, 'clone_commit':'2b3c4e611ca194d391a31b56f6bc93a28ca90b13',
        'canonical_parameters':CANONICAL})
    if args.patch_motion_cache is not None:
        source = args.patch_motion_cache
        if not (source / 'arrays.npz').exists():
            raise FileNotFoundError(f'Missing arrays.npz in {source}')
        (out / '01_patch_motion').mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / 'arrays.npz', out / '01_patch_motion' / 'arrays.npz')
        if (source / 'metadata.json').exists():
            shutil.copy2(source / 'metadata.json', out / '01_patch_motion' / 'metadata.json')
    else:
        if args.raw_motion_dir is None:
            p.error('provide either --raw-motion-dir or --patch-motion-cache')
        if args.invalid_mask_npy is None:
            p.error('--invalid-mask-npy is required with --raw-motion-dir')
        mask = np.load(args.invalid_mask_npy).astype(bool)
        stage1(args.raw_motion_dir, mask, args.slice, out, args.overwrite)
    if args.through_stage == '01': return
    stage2(out, args.overwrite)
    if args.through_stage == '02': return
    stage3(out, args.overwrite)
    if args.through_stage == '03': return
    stage4(out, args.overwrite)
    if args.through_stage == '04': return
    stage5(out, args.overwrite)

if __name__ == '__main__': main()
