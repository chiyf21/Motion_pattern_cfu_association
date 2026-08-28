#!/usr/bin/env python3
"""Render Fig4-style motion-arrow overlays for Fig5_v0827 spatial modules."""

from pathlib import Path
import csv
import pickle
import sys
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/home/cyf/wbi/wbi_code')
sys.path.insert(0, str(BASE))
CFU_SPATIAL = BASE / 'experiments/Fig5_v0827/03_cfu_pattern_spatial_overlap'
OUT = BASE / 'experiments/Fig5_v0827/03_cfu_pattern_spatial_overlap/figures/couple_overlays'
REF = Path('/mnt/data21T_2/cyf/f338/f338_registrated_0530/reference/vol_ref_000599_000999.tif')
PS = 7
GRID_STEP = 3

OUT.mkdir(parents=True, exist_ok=True)
manifest = []

for csv_path in sorted(CFU_SPATIAL.glob('spatial_final_slice*.csv')):
    slice_id = int(csv_path.stem.split('slice')[-1])
    pattern_path = BASE / f'experiments/Fig5_v0827/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{slice_id:02d}_velocity_decomp/06_patterns/objects.pkl'
    patterns = {p.pattern_id: p for p in pickle.load(pattern_path.open('rb'))['patterns']}
    reference = tifffile.imread(REF)[slice_id].astype(float)
    rows = list(csv.DictReader(csv_path.open()))
    output_slice = OUT / f'slice{slice_id:02d}'
    output_slice.mkdir(parents=True, exist_ok=True)

    for row in rows:
        pattern_id = int(row['pattern_id'])
        cfu_id = int(row['cfu_id'])
        p = patterns[pattern_id]
        mask = np.asarray(p.unified_mask, bool)
        response = np.asarray(p.unified_response_field, float)
        magnitude = np.linalg.norm(response, axis=2)
        ys, xs = np.where(mask)
        height, width = mask.shape
        ylo, yhi = max(0, ys.min() - 8), min(height, ys.max() + 9)
        xlo, xhi = max(0, xs.min() - 8), min(width, xs.max() + 9)

        yy, xx = np.mgrid[ylo:yhi:GRID_STEP, xlo:xhi:GRID_STEP]
        u = response[ylo:yhi:GRID_STEP, xlo:xhi:GRID_STEP, 0]
        v = response[ylo:yhi:GRID_STEP, xlo:xhi:GRID_STEP, 1]
        mag = magnitude[ylo:yhi:GRID_STEP, xlo:xhi:GRID_STEP]
        valid = mask[ylo:yhi:GRID_STEP, xlo:xhi:GRID_STEP] & np.isfinite(mag)
        norm = np.sqrt(u * u + v * v)
        scale_ref = np.nanpercentile(mag[valid], 98) if np.any(valid) else 1
        arrow_len = np.minimum(1.15, 1.8 * np.sqrt(np.clip(mag, 0, scale_ref) / (scale_ref + 1e-12)))
        u_norm = u / (norm + 1e-12) * arrow_len
        v_norm = v / (norm + 1e-12) * arrow_len
        u_norm[~valid] = np.nan
        v_norm[~valid] = np.nan

        background = reference[ylo * PS:yhi * PS, xlo * PS:xhi * PS]
        lo, hi = np.percentile(background, [1, 99])
        background = np.clip((background - lo) / (hi - lo + 1e-12), 0, 1) ** 0.65

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(background, cmap='gray', origin='lower',
                  extent=[xlo * PS, xhi * PS, ylo * PS, yhi * PS])
        ax.quiver(xx * PS, yy * PS, u_norm * PS, v_norm * PS,
                  color='#FF1493', angles='xy', scale=80, width=.013,
                  headwidth=3.5, headlength=5, headaxislength=4.5, alpha=.97)
        ax.invert_yaxis()
        ax.set_axis_off()
        fig.subplots_adjust(0, 0, 1, 1)
        output = output_slice / f'P{pattern_id:03d}_CFU{cfu_id:03d}.png'
        fig.savefig(output, dpi=220, bbox_inches='tight', pad_inches=0, facecolor='white')
        plt.close(fig)
        manifest.append([slice_id, pattern_id, cfu_id, row.get('ratio', ''),
                         row.get('coverage', ''), str(output)])
    print(f'slice{slice_id:02d}: {len(rows)} figures -> {output_slice}')

with (OUT / 'current_cfu_couple_overlay_manifest.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['slice', 'pattern_id', 'cfu_id', 'ratio', 'coverage', 'output'])
    writer.writerows(manifest)
print(f'total overlays: {len(manifest)}')
