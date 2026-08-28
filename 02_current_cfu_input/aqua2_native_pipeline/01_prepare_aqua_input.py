#!/usr/bin/env python3
"""
Stage 1: Prepare Ca patch data for AQuA2 analysis.

Loads the Ca patch stack from .npz, transposes to AQuA2's expected
(H, W, Z, T) layout, optionally computes dF/F, and saves as .mat files.
Also writes the AQuA2 batch parameter CSV.

Output:
    output/ca_patch_raw.mat       — raw fluorescence (H,W,Z,T)
    output/ca_patch_dff.mat       — dF/F (H,W,Z,T) [optional]
    output/parameters_for_batch.csv — AQuA2 config
"""

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
# Current server input used for the frozen DS7 CFU run (directly runnable here):
CA_CACHE = Path("/home/cyf/wbi/wbi_code/data/ca_cache/ca_patch_f338_full_ps7.npz")
# Portable placeholder (uncomment and edit for another machine):
# CA_CACHE = Path("/path/to/yourdata/ca_patch_f338_full_ps7.npz")


def load_ca_patch_mmap(path):
    """Load Ca patch stack in mmap mode (no memory copy)."""
    data = np.load(path, mmap_mode='r')
    return data['ca_patch_stack']  # (T, Z, H, W) float32


def save_as_hdf5(data, path, var_name='datOrg1'):
    """
    Save a 4D (H,W,Z,T) array as HDF5 for later MATLAB import.

    MATLAB can't directly load() arbitrary HDF5, so Stage 2's
    run_aqua_batch.m will convert this to .mat via h5read+save.
    """
    import h5py
    print(f"  Saving HDF5 to {path} ...")
    size_gb = data.nbytes / (1024**3)
    print(f"    Array shape: {data.shape}, size: {size_gb:.2f} GB")

    with h5py.File(str(path), 'w') as f:
        ds = f.create_dataset(
            var_name,
            data=data,
            compression='gzip',
            compression_opts=4,
            shuffle=True,
            chunks=True,
        )
        ds.attrs['dtype'] = str(data.dtype)
        ds.attrs['shape'] = data.shape

    print(f"  Done: {size_gb:.2f} GB written")


def transpose_for_aqua(ca_tzhw):
    """
    Convert from Python (T, Z, H, W) to AQuA2 (H, W, Z, T).

    AQuA2 expects datOrg1 with size [H, W, L, T] where L is layers (Z).
    burst.prep1.m does: if ndims==3 → permute([1,2,4,3]), i.e. H×W×T → H×W×1×T.
    For 4D data we simply permute to H×W×Z×T directly.
    """
    T, Z, H, W = ca_tzhw.shape
    print(f"  Input shape (T,Z,H,W): ({T},{Z},{H},{W})")
    # Transpose: (T,Z,H,W) -> (H,W,Z,T)
    ca_hwzt = np.transpose(ca_tzhw, (2, 3, 1, 0))
    print(f"  Output shape (H,W,Z,T): {ca_hwzt.shape}")
    return np.asarray(ca_hwzt, dtype=np.float32)


def compute_dff_global(ca_tzhw, baseline_percentile=20, spatial_chunk_size=50000):
    """
    Compute dF/F using a static global per-pixel baseline.

    F₀ = baseline_percentile-th percentile of the full 1598-frame trace for each pixel.
    dF/F = (F - F₀) / F₀

    This is fast (percentile along axis=0) and AQuA2's internal baseline
    removal handles residual drift. Processes spatial chunks for memory.

    Parameters
    ----------
    ca_tzhw : ndarray (T, Z, H, W)
    baseline_percentile : float
        Percentile for baseline F₀.

    Returns
    -------
    dff : ndarray (T, Z, H, W), float32
    """
    T, Z, H, W = ca_tzhw.shape
    N_spatial = Z * H * W

    # Step 1: Compute F₀ per pixel (process spatially in chunks)
    print(f"    Computing F₀ ({baseline_percentile}th percentile) per pixel...")
    ca_flat = np.asarray(ca_tzhw).reshape(T, N_spatial)
    f0 = np.zeros(N_spatial, dtype=np.float32)

    for s_start in range(0, N_spatial, spatial_chunk_size):
        s_end = min(s_start + spatial_chunk_size, N_spatial)
        chunk = ca_flat[:, s_start:s_end]
        f0[s_start:s_end] = np.percentile(chunk, baseline_percentile, axis=0)
        if (s_start // spatial_chunk_size) % 20 == 0:
            print(f"      F₀: {s_start}/{N_spatial} pixels ({100*s_start/N_spatial:.0f}%)")
    print(f"      F₀: {N_spatial}/{N_spatial} pixels (100%)")

    # Step 2: Compute dF/F in temporal chunks
    print(f"    Computing dF/F (T={T} frames)...")
    dff = np.zeros_like(ca_tzhw, dtype=np.float32)
    f0_3d = f0.reshape(Z, H, W)
    f0_3d = np.maximum(f0_3d, 1e-6)

    chunk_t = 200
    for t_start in range(0, T, chunk_t):
        t_end = min(t_start + chunk_t, T)
        slab = np.asarray(ca_tzhw[t_start:t_end], dtype=np.float32)
        dff[t_start:t_end] = (slab - f0_3d[np.newaxis, :, :, :]) / f0_3d[np.newaxis, :, :, :]
        if (t_start // chunk_t) % 10 == 0:
            print(f"      dF/F: {t_start}/{T} frames ({100*t_start/T:.0f}%)")
    print(f"      dF/F: {T}/{T} frames (100%)")

    return dff


def write_parameters_csv(csv_path):
    """Write AQuA2 batch parameters CSV (one row for our dataset)."""
    # Column order matches AQuA2's parameters.csv header
    # Values chosen for GCaMP patch data (500ms/frame, ~13 min recording)
    # Must match AQuA2's expected columns: Name, Variable, Type, Default, ...
    header = "Name,Variable,Type,Default,Long-duration signal,Notes"
    rows = [
        header,
        "Registration mode,registrateCorrect,preprocessing,1,1,",
        "Bleach correction mode,bleachCorrect,preprocessing,1,1,",
        "Median filter radius (For salt and pepper noise),medSmo,preprocessing,0,0,",
        "Spatial smoothing level,smoXY,preprocessing,0.5,0.5,Spatial smoothing filter size",
        ",,,,,",
        "Active voxels threshold scale,thrARScl,activeregion,3,3,",
        "Minimum duration of signal,minDur,activeregion,5,5,",
        "Minimum size,minSize,activeregion,20,20,Minimum event size in voxels",
        "Maximum size,maxSize,activeregion,inf,inf,",
        "Circularity threshold,circularityThr,activeregion,0,0,0=no limit",
        "Allowed distance in the same signal,spaMergeDist,activeregion,0,0,",
        ",,,,,",
        "Whether need temporal segmentation or not,needTemp,tempSeg,1,1,",
        "Seed size ratio compared to region ,seedSzRatio,tempSeg,0.01,0.01,",
        "Significance Threshold,sigThr,tempSeg,3.5,3.5,",
        "Delay score for merging,maxDelay,tempSeg,0.6,0.6,",
        "Whether signals are too close so that further refine is needed,needRefine,tempSeg,0,0,",
        "Whether signals are not large enough so that further growing is needed,needGrow,tempSeg,0,0,",
        ",,,,,",
        "Whether need spatial segmentation or not,needSpa,spaSeg,1,1,",
        "Source size relative to super event,sourceSzRatio,spaSeg,0.01,0.01,",
        "Sensitivity to detect source (Level 1 to 10),sourceSensitivity,spaSeg,8,8,",
        "Whether to extend event temporally,whetherExtend,spaSeg,1,1,",
        ",,,,,",
        "Whether detect global signals,detectGlo,glo,0,0,",
        "Minimum duration of global signal,gloDur,glo,20,20,",
        ",,,,,",
        "Ignore decay tau calculation,ignoreTau,post,1,1,",
        "Propagation metric in different directions,propMetric,post,0,0,Disabled for 3D",
        "Network features,networkFeatures,post,0,0,",
        ",,,,,",
        "Propagation smoothness,gtwSmo,latent,0.2,0.2,",
        "Intensity step for checking signal source,ratio,latent,0.5,0.5,",
        "Remove pixels close to image boundary,regMaskGap,latent,5,5,",
        "Frames per segment,cut,latent,200,800,Larger for long recordings",
        "Baseline window,movAvgWin,latent,25,25,",
        "Event show threshold on raw data,minShow1,latent,0.2,0.2,",
        "Correct baseline trend,correctTrend,latent,1,1,",
        "Propagation threshold minimum,propthrmin,latent,0.5,0.5,",
        "Propagation threshold step,propthrstep,latent,0.1,0.1,",
        "Propagation threshold maximum,propthrmax,latent,0.5,0.5,",
        "How long duration ratio should be seen as footprint,compress,latent,0,0,",
        "Check more time,gapExt,latent,5,5,",
        "Downsample duration,TPatch,latent,20,20,",
        "Max Downsample Spatial,maxSpaScale,latent,7,7,",
        "Min Downsample Spatial,minSpaScale,latent,3,3,",
        ",,,,,",
        "Frame rate,frameRate,from data,2.0,2.0,500ms per frame",
        "Spatial resolution,spatialRes,from data,1.0,1.0,patch resolution",
        "Estimated noise variance,varEst,from data,0.02,0.02,",
        "Foreground threshold,fgFluo,from data,0,0,auto-detect",
        "Background threshold,bgFluo,from data,0,0,auto-detect",
        "X cooridante for north vector,northx,from data,0,0,",
        "Y cooridante for north vector,northy,from data,1,1,",
    ]
    with open(csv_path, 'w') as f:
        f.write('\n'.join(rows))
    print(f"  Parameters CSV written to {csv_path}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Stage 1: Preparing Ca data for AQuA2")
    print("=" * 60)

    # 1. Load Ca patch data
    print("\n[1/4] Loading Ca patch stack (mmap)...")
    ca_tzhw = load_ca_patch_mmap(CA_CACHE)
    print(f"  Shape: {ca_tzhw.shape}, dtype: {ca_tzhw.dtype}")
    print(f"  Size: {ca_tzhw.nbytes / (1024**3):.2f} GB")

    # 2. Save raw Ca as .mat
    print("\n[2/4] Converting raw Ca to AQuA2 format...")
    ca_raw_hwzt = transpose_for_aqua(ca_tzhw)
    save_as_hdf5(ca_raw_hwzt, OUTPUT_DIR / "ca_patch_raw.h5")
    del ca_raw_hwzt

    # 3. Compute dF/F and save
    print("\n[3/4] Computing dF/F ...")
    dff_tzhw = compute_dff_global(ca_tzhw, baseline_percentile=20)

    # dF/F can have extreme outliers; mild clip for stability
    dff_low, dff_high = np.percentile(dff_tzhw, [0.01, 99.99])
    print(f"  dF/F range (p0.01-p99.99): [{dff_low:.4f}, {dff_high:.4f}]")
    # Only clip extreme outliers (beyond ±3σ of typical Ca dF/F)
    clip_low = max(dff_low, -2.0)
    clip_high = min(dff_high, 5.0)
    dff_tzhw = np.clip(dff_tzhw, clip_low, clip_high)
    print(f"  Clipped to [{clip_low}, {clip_high}]")
    # Do NOT rescale — AQuA2's burst.prep1 normalizes internally

    print("\n[3b/4] Converting dF/F to AQuA2 format...")
    ca_dff_hwzt = transpose_for_aqua(dff_tzhw)
    save_as_hdf5(ca_dff_hwzt, OUTPUT_DIR / "ca_patch_dff.h5")
    del dff_tzhw, ca_dff_hwzt

    # 4. Write config
    print("\n[4/4] Writing AQuA2 parameter CSV...")
    write_parameters_csv(OUTPUT_DIR / "parameters_for_batch.csv")

    print("\n" + "=" * 60)
    print("Stage 1 complete. Output files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        if f.is_file():
            size_mb = f.stat().st_size / (1024**2)
            print(f"  {f.name}  ({size_mb:.1f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main()
