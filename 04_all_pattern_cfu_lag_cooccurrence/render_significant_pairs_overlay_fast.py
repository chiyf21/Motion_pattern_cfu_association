#!/usr/bin/env python3
"""Fast raster overlay renderer for significant pattern-CFU pairs."""
import csv
import importlib.util
import pickle
from pathlib import Path

import h5py
import numpy as np
import tifffile
from PIL import Image, ImageDraw

BASE = Path('/home/cyf/wbi/wbi_code')
EXP = BASE / 'experiments/motion_pattern_cfu_association/04_all_pattern_cfu_lag_cooccurrence'
RESULT = EXP / 'results/global_shift_empirical_fdr_onset'
OUT = EXP / 'figures/significant_pairs_spatiotemporal_original_overlay_fast'
CFU_BASE = BASE / 'experiments/motion_pattern_cfu_association/02_current_cfu_input/cfu'
OUT.mkdir(parents=True, exist_ok=True)
T, WINDOW, PS, GAP = 1598, 3, 7, 5
REF_TIF = Path('/mnt/data21T_2/cyf/f338/f338_registrated_0530/reference/vol_ref_000599_000999.tif')
PANEL_W, PANEL_H = 680, 500

spec = importlib.util.spec_from_file_location('fdr', EXP / 'run_all_cfu_pattern_lag8_w3_pairwise.py')
fdr = importlib.util.module_from_spec(spec); spec.loader.exec_module(fdr)

def load_patterns(z):
    p = BASE / f'experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{z:02d}_velocity_decomp/06_patterns/objects.pkl'
    with p.open('rb') as f: return {int(x.pattern_id): x for x in pickle.load(f)['patterns']}

def load_cfu_masks(z):
    p = CFU_BASE / f'slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat'
    out = {}
    with h5py.File(p, 'r') as f:
        cells=f['cfuInfo1']; layout,n=fdr._detect_cfu_cell_layout(cells)
        for i in range(n):
            ref=fdr._get_cfu_field_ref(cells,layout,i,2); out[i+1]=np.asarray(f[ref]).T > 0
    return out

def load_timing(z): return {cid:on for cid,_,on in fdr.cfu_event_onsets(z)}
def load_peaks(z): return {pid:pk for pid,_,pk in fdr.motion_event_peaks(z)}

def norm_bg(x):
    finite=x[np.isfinite(x)]; lo,hi=np.percentile(finite,[1,99]); return np.clip((x-lo)/max(hi-lo,1e-6),0,1)**.65

def cfu_to_ref(mask, shape):
    h,w=shape; yp=np.arange(h)//PS-GAP; xp=np.arange(w)//PS-GAP
    vy=(yp>=0)&(yp<mask.shape[0]); vx=(xp>=0)&(xp<mask.shape[1])
    out=mask[np.ix_(np.clip(yp,0,mask.shape[0]-1),np.clip(xp,0,mask.shape[1]-1))]
    out[~vy,:]=False; out[:,~vx]=False; return out

def pattern_to_ref(mask, shape):
    x=np.kron(mask.astype(bool),np.ones((PS,PS),bool)); out=np.zeros(shape,bool)
    out[:min(shape[0],x.shape[0]),:min(shape[1],x.shape[1])]=x[:shape[0],:shape[1]]
    return out

def make_bg_panel(bg):
    return Image.fromarray(np.uint8(np.clip(bg*255,0,255)), 'L').resize((PANEL_W,PANEL_H), Image.Resampling.BILINEAR).convert('RGBA')

def panel(bg_panel, mask, color):
    im=bg_panel.copy()
    m=Image.fromarray(np.uint8(mask*255), 'L').resize((PANEL_W,PANEL_H), Image.Resampling.NEAREST)
    overlay=Image.new('RGBA',(PANEL_W,PANEL_H),color+(0,)); overlay.putalpha(m.point(lambda v:int(v*.62)))
    return Image.alpha_composite(im,overlay)

def draw_pair(row, patterns, masks, timing, peaks, bgs, bg_panels):
    ms,pid,cs,cid=(int(row['motion_slice_0based']),int(row['pattern_id']),int(row['cfu_slice_0based']),int(row['cfu_id']))
    lag=int(row['best_lag']); q=float(row['q_empirical_global']); refp=bgs[ms]; refc=bgs[cs]
    pm=pattern_to_ref(np.asarray(patterns[ms][pid].unified_mask),refp.shape); cm=cfu_to_ref(masks[cs][cid],refc.shape)
    canvas=Image.new('RGB',(1400,1030),'white'); canvas.paste(panel(bg_panels[ms],pm,(235,45,85)),(20,55)); canvas.paste(panel(bg_panels[cs],cm,(20,150,235)),(700,55))
    d=ImageDraw.Draw(canvas); d.text((20,15),f'Pattern slice{ms:02d} P{pid:03d}  ×  CFU slice{cs:02d} CFU{cid:03d} | q={q:.4g} | best lag={lag} frames',fill='black')
    d.text((25,565),f'Pattern mask on reference: slice{ms:02d}',fill=(180,20,50)); d.text((705,565),f'CFU mask on reference: slice{cs:02d}',fill=(10,90,180))
    x0,y0,x1,y1=60,650,1340,980; d.rectangle((x0,y0,x1,y1),outline=(100,100,100),width=1); d.line((x0,780,x1,780),fill=(210,210,210),width=1)
    pon=np.asarray(peaks[ms][pid],int); con=np.asarray(timing[cs][cid],int); sx=(x1-x0)/T
    for x in pon:
        xx=int(x0+x*sx); d.line((xx,705,xx,735),fill=(194,55,55),width=2)
        st=int(x+lag); en=st+WINDOW
        if len(con[(con>=st)&(con<en)]):
            xa=int(x0+max(0,st)*sx); xb=int(x0+min(T,en)*sx); d.rectangle((xa,790,xb,815),fill=(242,177,52))
    for x in con:
        xx=int(x0+x*sx); d.line((xx,820,xx,850),fill=(20,110,180),width=2)
    d.text((70,685),'pattern motion peaks',fill=(194,55,55)); d.text((70,855),f'CFU onsets | {len(pon)} peaks, {len(con)} onsets | orange=best-lag hit window',fill=(20,90,150))
    for tick in range(0,T+1,200):
        xx=int(x0+tick*sx); d.line((xx,960,xx,968),fill=(80,80,80)); d.text((xx-10,970),str(tick),fill=(60,60,60))
    name=f'slice{ms:02d}_P{pid:03d}__slice{cs:02d}_CFU{cid:03d}.png'; canvas.save(OUT/name,compress_level=1); return name

def main():
    with (RESULT/'empirical_FDR_significant_pairs.csv').open() as f: rows=list(csv.DictReader(f))
    slices=sorted({int(r['motion_slice_0based']) for r in rows}|{int(r['cfu_slice_0based']) for r in rows})
    stack=tifffile.imread(REF_TIF).astype(np.float32); bgs={z:norm_bg(stack[z]) for z in slices}; bg_panels={z:make_bg_panel(bgs[z]) for z in slices}
    patterns={z:load_patterns(z) for z in slices}; masks={z:load_cfu_masks(z) for z in slices}; timing={z:load_timing(z) for z in slices}; peaks={z:load_peaks(z) for z in slices}
    manifest=[]
    for i,r in enumerate(rows,1):
        name=draw_pair(r,patterns,masks,timing,peaks,bgs,bg_panels); x=dict(r); x['figure']=name; manifest.append(x)
        if i==1 or i%50==0 or i==len(rows): print(f'rendered {i}/{len(rows)}',flush=True)
    with (OUT/'manifest.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
    top=sorted(rows,key=lambda r:float(r['q_empirical_global']))[:20]; sheet=Image.new('RGB',(1400,1030),'white')
    for i,r in enumerate(top):
        im=Image.open(OUT/draw_pair(r,patterns,masks,timing,peaks,bgs,bg_panels)).convert('RGB').resize((350,257),Image.Resampling.BILINEAR); sheet.paste(im,((i%4)*350,(i//4)*257))
    sheet.save(OUT/'top20_q_montage.png',compress_level=1); (OUT/'README.md').write_text('motion_pattern_cfu_association (omega=mu=0.5) fast raster overlay version. Original registered reference is shown in both spatial panels. Pattern masks are mapped by PS=7; CFU masks are transposed from HDF5 and mapped with regMaskGap=5. Timeline uses 0-based frames; red=pattern motion peaks, blue=CFU onsets, orange=best-lag hit windows.\n')
    print('rendered',len(rows),'pairs'); print('output',OUT)

if __name__=='__main__': main()
