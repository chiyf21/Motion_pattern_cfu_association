#!/usr/bin/env python3
"""Render canonical motion-pattern overviews for patterns with >=5 members."""
from pathlib import Path
import pickle
import sys
import numpy as np
import tifffile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import REFERENCE_TIF
PATTERN_ROOT = ROOT / '01_motion_pattern_extraction_omega05_mu05/patterns'
OUT = PATTERN_ROOT / '_overview_members_ge5'
REFERENCE = REFERENCE_TIF
PATCH = 7

def norm(x):
    lo, hi = np.percentile(x, [1, 99])
    return np.clip((x.astype(float)-lo)/(hi-lo+1e-12), 0, 1)

def draw(z, ref, save_path):
    with (PATTERN_ROOT / f'Slice{z:02d}_velocity_decomp/06_patterns/objects.pkl').open('rb') as f:
        patterns = pickle.load(f)['patterns']
    pats = [p for p in patterns if int(getattr(p, 'n_members', len(getattr(p, 'regions', [])))) >= 5]
    h, w = ref.shape
    ph, pw = h // PATCH, w // PATCH
    cmap = plt.cm.hsv(np.linspace(0, .96, max(1, len(pats))))
    rgb = np.repeat(norm(ref)[..., None], 3, axis=2)
    claimed = np.zeros((ph, pw), bool)
    for p, color in zip(sorted(pats, key=lambda x: x.pattern_id), cmap):
        m = np.asarray(p.unified_mask, bool)
        mh, mw = min(ph, m.shape[0]), min(pw, m.shape[1])
        new = m[:mh, :mw] & ~claimed[:mh, :mw]
        px = np.repeat(np.repeat(new, PATCH, 0), PATCH, 1)
        hh, ww = min(h, px.shape[0]), min(w, px.shape[1])
        rgb[:hh, :ww][px[:hh, :ww]] = .35 * rgb[:hh, :ww][px[:hh, :ww]] + .65 * np.asarray(color[:3])
        claimed[:mh, :mw] |= m[:mh, :mw]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.imshow(rgb, origin='upper'); ax.axis('off')
    ax.set_title(f'Slice {z:02d}: {len(pats)} motion patterns (>=5 members)')
    fig.savefig(save_path, dpi=220, bbox_inches='tight')
    plt.close(fig)
    return len(patterns), len(pats)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref = tifffile.imread(REFERENCE)
    rows = []
    for z in range(1, 13):
        n_all, n_ge5 = draw(z, ref[z-1], OUT / f'Slice{z:02d}_overview_members_ge5.png')
        rows.append((z, n_all, n_ge5))
    with (OUT / 'manifest.csv').open('w') as f:
        f.write('slice_display,slice_0based,n_patterns,n_members_ge5,output\n')
        for z, n, k in rows:
            f.write(f'{z},{z-1},{n},{k},Slice{z:02d}_overview_members_ge5.png\n')
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    for ax, (z, _, k) in zip(axes.ravel(), rows):
        img = plt.imread(OUT / f'Slice{z:02d}_overview_members_ge5.png')
        ax.imshow(img); ax.set_title(f'Slice {z:02d} | >=5: {k}'); ax.axis('off')
    fig.savefig(OUT / 'all_slices_overview_members_ge5.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'created {len(rows)} slice overviews in {OUT}')

if __name__ == '__main__': main()
