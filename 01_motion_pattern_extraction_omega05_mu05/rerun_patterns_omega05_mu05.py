#!/usr/bin/env python3
"""Independent Stage 6 rerun with omega=0.5, mu=0.5 + best_cc unified.
Usage: python rerun_patterns_spatial.py [slice1 slice2 ...]  (default: 1..12)
Saves under this experiment's own patterns/SliceXX_velocity_decomp/06_patterns/.
"""
import sys, pickle, time, warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings('ignore')

sys.path.insert(0, "/home/cyf/wbi/wbi_code")
import src_registration.motion_correlation_pattern_v2 as mcp_v2
import src_registration.motion_stage_cache as msc

DATA_BASE = Path("/home/cyf/wbi/wbi_code/data/f338_velocity_decomp")
OUT_BASE = Path("/home/cyf/wbi/wbi_code/experiments/Fig5_v0827/01_motion_pattern_extraction_omega05_mu05/patterns")
STAGE = "06_patterns"


def run_slice(slice_idx):
    cache = DATA_BASE / f"Slice{slice_idx:02d}_velocity_decomp"
    t0 = time.time()
    with open(cache / "04_modes" / "objects.pkl", "rb") as f:
        data = pickle.load(f)
    modes = data.get("modes", [])
    episodes = data.get("episodes", [])
    ep_by_id = {}
    for ep in episodes:
        eid = getattr(ep, "episode_id", -1)
        ep_by_id[eid] = ep
        ep.modes = []
    for m in modes:
        eid = getattr(m, "episode_id", -1)
        if eid in ep_by_id:
            ep_by_id[eid].modes.append(m)

    pp = dict(
        min_strength=0.0, min_area=5, min_duration=1,
        min_iou=0.08, omega=0.5, mu=0.5,
        b_distance="correlation",
        spatial_rule="iou",
        cluster_dist_thresh=0.45,
        linkage_method="complete", incompatible_dist=1e6,
        compute_unified=True,
        unified_mask_mode="best_cc",
        unified_sign_method="correlation",
        min_pattern_members=2,
        min_unified_area=50,
    )
    patterns, kept, groups, labels, info = mcp_v2.getMotionPattern(
        episodes, unit_type="mode", verbose=False, **pp)
    out_cache = OUT_BASE / f"Slice{slice_idx:02d}_velocity_decomp"
    msc.save_patterns_stage(str(out_cache), patterns=patterns,
                            kept_regions=kept, groups=groups,
                            labels=labels, info=info, params=pp,
                            stage_name=STAGE)
    sizes = [p.n_members for p in patterns]
    print(f"  Slice {slice_idx:02d}: {len(episodes)} ep, {len(modes)} modes -> "
          f"{len(patterns)} patterns (>=5: {sum(1 for s in sizes if s >= 5)}, "
          f"max={max(sizes) if sizes else 0}) in {time.time()-t0:.0f}s", flush=True)
    return len(patterns)


if __name__ == "__main__":
    slices = [int(x) for x in sys.argv[1:]] or list(range(1, 13))
    for s in slices:
        run_slice(s)
