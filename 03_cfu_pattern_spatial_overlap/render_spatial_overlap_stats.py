#!/usr/bin/env python3
"""Create the Fig4-style spatial-module summary for Fig5_v0827."""

from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/home/cyf/wbi/wbi_code')
RES = BASE / 'experiments/Fig5_v0827/03_cfu_pattern_spatial_overlap'
OUT = BASE / 'experiments/Fig5_v0827/03_cfu_pattern_spatial_overlap/figures/fig4_style_spatial_overlap'
OUT.mkdir(parents=True, exist_ok=True)

all_rows = []
all_x = []
all_y = []
selected = []
for z in range(1, 13):
    with (RES / f'spatial_all_overlap_slice{z}.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    with (RES / f'spatial_final_slice{z}.csv').open(newline='') as f:
        final = {(int(r['pattern_id']), int(r['cfu_id'])) for r in csv.DictReader(f)}
    for r in rows:
        if int(r['pattern_n_members']) < 5:
            continue
        ratio = float(r['ratio'])
        coverage = float(r['coverage'])
        if not (np.isfinite(ratio) and np.isfinite(coverage) and ratio > 0):
            continue
        is_selected = (int(r['pattern_id']), int(r['cfu_id'])) in final
        all_rows.append([z, int(r['pattern_id']), int(r['cfu_id']), ratio, coverage, is_selected])
        all_x.append(ratio)
        all_y.append(coverage)
        selected.append(is_selected)

x = np.asarray(all_x, float)
y = np.asarray(all_y, float)
selected = np.asarray(selected, bool)
fig = plt.figure(figsize=(10.5, 8.0))
ax = fig.add_axes([0.16, 0.17, 0.79, 0.75])
ax.scatter(x, y, s=28, c='#5A5A5A', alpha=0.48, edgecolors='none')
ax.scatter(x[selected], y[selected], s=90, c='#8B1E3F', alpha=0.92,
           edgecolor='white', linewidth=0.3)
ax.axvline(3, color='black', linestyle='--', linewidth=1.3)
ax.axhline(0.5, color='black', linestyle='--', linewidth=1.3)
ax.set_xscale('log')
ax.set_xlim(0.7, max(x) * 1.1)
ax.set_ylim(0, 1)
ax.set_xticks([1, 3, 10, 100, 1000])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xlabel('Relative area (log)', fontsize=29)
ax.set_ylabel('IoM', fontsize=29)
ax.tick_params(labelsize=27)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.savefig(OUT / 'fig5_all_slices_fig3_stats.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT / 'fig5_all_slices_fig3_stats.pdf', bbox_inches='tight', facecolor='white')
plt.close(fig)

with (OUT / 'fig5_all_slices_fig3_stats_points.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['slice', 'pattern_id', 'cfu_id', 'ratio', 'coverage', 'selected'])
    writer.writerows(all_rows)
print(f'candidates={len(x)} selected={int(selected.sum())} slices=12')
