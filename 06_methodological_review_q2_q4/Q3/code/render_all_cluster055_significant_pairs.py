#!/usr/bin/env python3
"""Render Q3 traces for every q<0.05 pair in the cluster055 temporal results."""
import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, required=True)
    ap.add_argument("--review-root", type=Path, required=True)
    ap.add_argument("--q3-metrics", type=Path, required=True)
    ap.add_argument("--q3-full-proxy", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    sys.path[:0] = [str(args.review_root / "report"), str(args.review_root)]
    import make_followup_figures as figures
    import review_common

    figures.REPO_ROOT = args.project_root
    figures.PATTERN_ROOT = args.project_root / "01_motion_pattern_extraction_omega05_mu05" / "patterns"
    figures.CFU_ROOT = args.project_root / "02_current_cfu_input" / "cfu"
    review_common.REPO_ROOT = figures.REPO_ROOT
    review_common.PATTERN_ROOT = figures.PATTERN_ROOT
    review_common.CFU_ROOT = figures.CFU_ROOT

    pair_file = (args.project_root / "04_all_pattern_cfu_lag_cooccurrence" / "results" /
                 "global_shift_empirical_fdr_onset" / "empirical_FDR_significant_pairs.csv")

    def cluster055_pairs():
        return [r for r in review_common.read_csv(pair_file)
                if float(r["q_empirical_global"]) < 0.05]

    figures.load_current_significant_pairs = cluster055_pairs
    args.out.mkdir(parents=True, exist_ok=True)
    figures.q3_figures(args.out, args.q3_metrics, args.q3_full_proxy,
                       selected_pair_keys=None, shared_mse_scale=True)
    print(f"Rendered all {len(cluster055_pairs())} cluster055 q<0.05 pairs to {args.out}")


if __name__ == "__main__":
    main()
