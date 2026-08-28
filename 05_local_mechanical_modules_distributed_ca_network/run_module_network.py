#!/usr/bin/env python3
"""Build current spatial CFU-pattern modules and their distributed Ca network.

This analysis intentionally recomputes spatial modules from the same current
CFU HDF5 files used by Fig5/23 global co-occurrence.  It does not rerun
AQuA2 event detection or the all-pattern x all-CFU statistical test.
"""
from __future__ import annotations

import csv
import os
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np


BASE = Path("/home/cyf/wbi/wbi_code")
EXP23 = BASE / "experiments/motion_pattern_cfu_association/04_all_pattern_cfu_lag_cooccurrence"
OUT = BASE / "experiments/motion_pattern_cfu_association/05_local_mechanical_modules_distributed_ca_network"
# The current spatial result is a standalone, reusable analysis stage.  Fig5/24
# consumes it together with Fig5/23's temporal associations.
SPATIAL_OUT = BASE / "experiments/motion_pattern_cfu_association/03_cfu_pattern_spatial_overlap"
NETWORK_OUT = OUT / "02_module_cfu_network"
CFU_LINKS = BASE / "experiments/motion_pattern_cfu_association/02_current_cfu_input/cfu"
FDR_RESULT = (
    EXP23
    / "results/global_shift_empirical_fdr_onset/empirical_FDR_significant_pairs.csv"
)
Q001_RESULT = (
    EXP23
    / "results/global_shift_empirical_fdr_onset/empirical_FDR_pairs_q001.csv"
)

SLICES = range(1, 13)
PATTERN_SHAPE = (244, 329)
REG_MASK_GAP = 5
MASK_THRESHOLD = 0.1
RATIO_THRESHOLD = 3.0
COVERAGE_THRESHOLD = 0.5
MIN_PATTERN_MEMBERS = 5

sys.path.insert(0, str(BASE))


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def cfu_field_ref(cells, field_zero_based: int, cfu_zero_based: int):
    """Read a MATLAB cell field reference independent of HDF5 cell layout."""
    if cells.shape[0] == 9:
        return cells[field_zero_based, cfu_zero_based]
    if cells.shape[1] == 9:
        return cells[cfu_zero_based, field_zero_based]
    raise ValueError(f"Unrecognized cfuInfo1 shape: {cells.shape}")


def load_current_cfus(z: int):
    """Return current CFU masks on the shared 244 x 329 pattern grid."""
    path = CFU_LINKS / f"slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat"
    out = []
    with h5py.File(path, "r") as f:
        cells = f["cfuInfo1"]
        n_cfu = cells.shape[1] if cells.shape[0] == 9 else cells.shape[0]
        for cid0 in range(n_cfu):
            # HDF5 reverses MATLAB spatial dimensions.  Transpose back to
            # AQuA2's y,x convention before adding regMaskGap=5 pixels.
            weight = np.asarray(f[cfu_field_ref(cells, 2, cid0)]).T
            if weight.shape != (234, 319):
                raise ValueError(f"Z{z:02d} CFU{cid0 + 1}: unexpected mask {weight.shape}")
            mask = np.zeros(PATTERN_SHAPE, dtype=bool)
            mask[REG_MASK_GAP:-REG_MASK_GAP, REG_MASK_GAP:-REG_MASK_GAP] = (
                weight > MASK_THRESHOLD
            )
            event_ids = np.asarray(f[cfu_field_ref(cells, 1, cid0)]).reshape(-1)
            out.append({"cfu_id": cid0 + 1, "mask": mask, "n_events": len(event_ids)})
    return out


def load_patterns(z: int):
    path = (
        BASE
        / f"experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{z:02d}_velocity_decomp"
        / "06_patterns/objects.pkl"
    )
    with path.open("rb") as f:
        patterns = pickle.load(f)["patterns"]
    return sorted(patterns, key=lambda p: int(p.pattern_id))


def write_input_version_audit():
    """Record the exact current CFU and pattern inputs used by Fig5/24."""
    rows = []
    for z in SLICES:
        cfu_link = CFU_LINKS / f"slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat"
        cfu_real = Path(os.path.realpath(cfu_link))
        pattern_path = (
            BASE
            / f"experiments/motion_pattern_cfu_association/01_motion_pattern_extraction_omega05_mu05/patterns/Slice{z:02d}_velocity_decomp"
            / "06_patterns/objects.pkl"
        )
        with h5py.File(cfu_link, "r") as f:
            cells = f["cfuInfo1"]
            n_cfu = cells.shape[1] if cells.shape[0] == 9 else cells.shape[0]
            n_events = int(np.asarray(f["nEvents"]).squeeze())
        rows.append({
            "slice_0based": z,
            "cfu_link": str(cfu_link),
            "cfu_resolved_path": str(cfu_real),
            "cfu_file_mtime_utc": cfu_real.stat().st_mtime_ns,
            "aqua2_n_events": n_events,
            "aqua2_n_cfus": n_cfu,
            "pattern_path": str(pattern_path),
            "pattern_file_mtime_utc": pattern_path.stat().st_mtime_ns,
        })
    fields = list(rows[0])
    write_csv(OUT / "00_input_version_audit.csv", rows, fields)
    (OUT / "00_input_version_audit.md").write_text(
        "# Fig5/24 input-version audit\n\n"
        "Spatial modules are computed from the Fig5/23 `input_links/cfu` files listed in the CSV. "
        "The temporal network is joined from Fig5/23 `empirical_FDR_significant_pairs.csv`, "
        "whose runner reads the same `input_links/cfu` directory. Pattern masks and activations "
        "come from the listed motion_pattern_cfu_association pattern `06_patterns/objects.pkl` files.\n"
    )


SPATIAL_FIELDS = [
    "slice_0based", "pattern_id", "cfu_id", "pattern_n_members", "cfu_n_events",
    "Ap_pattern_pixels", "Acfu_pixels", "overlap_pixels", "ratio", "coverage",
    "passes_overlap", "passes_ratio", "passes_coverage", "passes_members",
]


def compute_current_spatial_modules():
    all_final = []
    summaries = []
    for z in SLICES:
        cfus = load_current_cfus(z)
        patterns = load_patterns(z)
        rows = []
        for pattern in patterns:
            pid = int(pattern.pattern_id)
            pmask = np.asarray(pattern.unified_mask, dtype=bool)
            if pmask.shape != PATTERN_SHAPE:
                raise ValueError(f"Z{z:02d} P{pid}: unexpected mask {pmask.shape}")
            ap = int(pmask.sum())
            members = int(pattern.n_members)
            for cfu in cfus:
                cmask = cfu["mask"]
                ac = int(cmask.sum())
                overlap = int(np.logical_and(pmask, cmask).sum())
                if not overlap or not ap or not ac:
                    continue
                ratio = max(ap, ac) / min(ap, ac)
                coverage = overlap / min(ap, ac)
                rows.append({
                    "slice_0based": z,
                    "pattern_id": pid,
                    "cfu_id": cfu["cfu_id"],
                    "pattern_n_members": members,
                    "cfu_n_events": cfu["n_events"],
                    "Ap_pattern_pixels": ap,
                    "Acfu_pixels": ac,
                    "overlap_pixels": overlap,
                    "ratio": ratio,
                    "coverage": coverage,
                    "passes_overlap": 1,
                    "passes_ratio": int(ratio <= RATIO_THRESHOLD),
                    "passes_coverage": int(coverage >= COVERAGE_THRESHOLD),
                    "passes_members": int(members >= MIN_PATTERN_MEMBERS),
                })
        rows.sort(key=lambda r: (-r["coverage"], r["pattern_id"], r["cfu_id"]))
        final = [
            r for r in rows
            if r["passes_ratio"] and r["passes_coverage"] and r["passes_members"]
        ]
        write_csv(SPATIAL_OUT / f"spatial_all_overlap_slice{z}.csv", rows, SPATIAL_FIELDS)
        write_csv(SPATIAL_OUT / f"spatial_final_slice{z}.csv", final, SPATIAL_FIELDS)
        all_final.extend(final)
        summaries.append({
            "slice_0based": z,
            "n_current_cfus": len(cfus),
            "n_patterns": len(patterns),
            "n_overlap_pairs": len(rows),
            "n_modules": len(final),
            "n_module_patterns": len({r["pattern_id"] for r in final}),
            "n_local_cfus": len({r["cfu_id"] for r in final}),
        })
        print(f"Z{z:02d}: current CFUs={len(cfus)}, modules={len(final)}", flush=True)
    write_csv(
        SPATIAL_OUT / "summary_all_slices.csv", summaries,
        list(summaries[0]) if summaries else ["slice_0based"],
    )
    return all_final, summaries


def read_rows(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


NETWORK_FIELDS = [
    "module_id", "module_slice_0based", "module_pattern_id", "local_cfu_id",
    "module_ratio", "module_coverage", "source_cfu_slice_0based", "source_cfu_id",
    "source_scope", "temporal_direction", "candidate_upstream", "best_lag",
    "window_frames", "q_empirical_global", "best_lag_binomial_p", "pair_score",
    "n_motion_events", "k_hit_motion_events_at_best_lag", "n_valid_motion_events_at_best_lag",
]


def module_rows(spatial_final: list[dict]):
    modules = []
    for r in spatial_final:
        z, pid, cid = int(r["slice_0based"]), int(r["pattern_id"]), int(r["cfu_id"])
        modules.append({
            "module_id": f"M_S{z:02d}_P{pid:03d}_C{cid:03d}",
            "module_slice_0based": z,
            "module_pattern_id": pid,
            "local_cfu_id": cid,
            "module_ratio": float(r["ratio"]),
            "module_coverage": float(r["coverage"]),
            "pattern_n_members": int(r["pattern_n_members"]),
            "local_cfu_n_events": int(r["cfu_n_events"]),
        })
    return modules


def join_network(modules: list[dict], association_rows: list[dict]):
    by_pattern = defaultdict(list)
    for row in association_rows:
        by_pattern[(int(row["motion_slice_0based"]), int(row["pattern_id"]))].append(row)
    joined = []
    for module in modules:
        key = (module["module_slice_0based"], module["module_pattern_id"])
        for row in by_pattern.get(key, []):
            source_slice = int(row["cfu_slice_0based"])
            source_cfu = int(row["cfu_id"])
            lag = int(row["best_lag"])
            local = source_slice == module["module_slice_0based"] and source_cfu == module["local_cfu_id"]
            if local:
                scope = "intramodule_local_cfu"
            elif source_slice == module["module_slice_0based"]:
                scope = "same_slice_nonlocal_cfu"
            else:
                scope = "cross_slice_cfu"
            direction = "candidate_upstream" if lag < 0 else ("synchronous" if lag == 0 else "candidate_downstream")
            joined.append({
                **{k: module[k] for k in [
                    "module_id", "module_slice_0based", "module_pattern_id", "local_cfu_id",
                    "module_ratio", "module_coverage",
                ]},
                "source_cfu_slice_0based": source_slice,
                "source_cfu_id": source_cfu,
                "source_scope": scope,
                "temporal_direction": direction,
                "candidate_upstream": int((not local) and lag < 0),
                "best_lag": lag,
                "window_frames": int(row["window_frames"]),
                "q_empirical_global": float(row["q_empirical_global"]),
                "best_lag_binomial_p": float(row["best_lag_binomial_p"]),
                "pair_score": float(row["pair_score"]),
                "n_motion_events": int(row["n_motion_events"]),
                "k_hit_motion_events_at_best_lag": int(row["k_hit_motion_events_at_best_lag"]),
                "n_valid_motion_events_at_best_lag": int(row["n_valid_motion_events_at_best_lag"]),
            })
    return joined


def summarize(modules, links, label: str):
    covered = {r["module_id"] for r in links}
    upstream = [r for r in links if r["candidate_upstream"]]
    counts = Counter((r["source_scope"], r["temporal_direction"]) for r in links)
    return {
        "threshold_set": label,
        "n_modules": len(modules),
        "n_modules_with_any_cfu_association": len(covered),
        "n_module_cfu_edges": len(links),
        "n_intramodule_edges": sum(r["source_scope"] == "intramodule_local_cfu" for r in links),
        "n_nonlocal_edges": sum(r["source_scope"] != "intramodule_local_cfu" for r in links),
        "n_candidate_upstream_edges": len(upstream),
        "n_modules_with_candidate_upstream": len({r["module_id"] for r in upstream}),
        "n_unique_upstream_cfus": len({(r["source_cfu_slice_0based"], r["source_cfu_id"]) for r in upstream}),
        "scope_direction_counts": "; ".join(
            f"{scope}|{direction}={n}" for (scope, direction), n in sorted(counts.items())
        ),
    }


def main(spatial_only: bool = False):
    SPATIAL_OUT.mkdir(parents=True, exist_ok=True)
    NETWORK_OUT.mkdir(parents=True, exist_ok=True)
    write_input_version_audit()
    spatial_final, spatial_summary = compute_current_spatial_modules()
    modules = module_rows(spatial_final)
    write_csv(
        NETWORK_OUT / "current_module_table.csv", modules,
        list(modules[0]) if modules else ["module_id"],
    )

    if spatial_only:
        (SPATIAL_OUT / "README.md").write_text(
            "# motion_pattern_cfu_association current spatial CFU-pattern modules\n\n"
            "Computed from the motion_pattern_cfu_association omega=0.5, mu=0.5 pattern files and "
            "the linked current CFUs. Criteria: pattern members>=5, ratio<=3, coverage>=0.5. "
            "Run `run_spatial_overlap.py` in this directory to recompute this stage only.\n"
        )
        print(f"current spatial modules: {len(modules)}")
        return

    fdr_rows = read_rows(FDR_RESULT)
    # The pairwise runner writes the q<=0.05 table directly.  q<=0.001 is a
    # strict subset, so derive and persist its compatibility table here when
    # the standalone filtering helper has not been run.
    if Q001_RESULT.exists():
        q001_rows = read_rows(Q001_RESULT)
    else:
        q001_rows = [r for r in fdr_rows if float(r["q_empirical_global"]) <= 0.001]
        write_csv(
            Q001_RESULT,
            q001_rows,
            list(fdr_rows[0]) if fdr_rows else ["q_empirical_global"],
        )
    fdr_links = join_network(modules, fdr_rows)
    q001_links = join_network(modules, q001_rows)
    write_csv(NETWORK_OUT / "module_cfu_associations_fdr_q005.csv", fdr_links, NETWORK_FIELDS)
    write_csv(NETWORK_OUT / "module_cfu_associations_q001.csv", q001_links, NETWORK_FIELDS)
    upstream = [r for r in fdr_links if r["candidate_upstream"]]
    write_csv(NETWORK_OUT / "candidate_upstream_cfu_to_module_edges_q005.csv", upstream, NETWORK_FIELDS)

    summaries = [summarize(modules, fdr_links, "empirical-FDR q<=0.05"), summarize(modules, q001_links, "q<=0.001")]
    write_csv(NETWORK_OUT / "network_summary.csv", summaries, list(summaries[0]))
    (NETWORK_OUT / "README.md").write_text(
        "# Local mechanical modules embedded in the distributed Ca network\n\n"
        "Spatial modules were recomputed from the same current CFU HDF5 inputs used by Fig5/23 global co-occurrence, with ratio<=3, coverage>=0.5, and pattern members>=5.\n\n"
        "A module is a local pattern-CFU spatial pair.  For every module pattern, all significant co-occurrence CFUs are joined as candidate network inputs.  A non-local CFU is labelled candidate_upstream only when best_lag<0, because the Fig5/23 convention is CFU onset = motion activation peak + best_lag.  Lag 0 is synchronous and positive lag is candidate_downstream; these are association labels, not causal proof.\n"
    )
    (OUT / "README.md").write_text(
        "This experiment connects current spatial CFU-pattern modules to the current all-pattern x all-CFU empirical-FDR co-occurrence results.  It recomputes only spatial overlap; AQuA2 event detection and global co-occurrence are reused unchanged.\n"
    )
    print(f"current modules: {len(modules)}")
    print(f"FDR q<=0.05 module-CFU edges: {len(fdr_links)}; candidate upstream: {len(upstream)}")
    print(f"output: {OUT}")


if __name__ == "__main__":
    main(spatial_only=("--spatial-only" in sys.argv))
