"""Small adapters for the existing motion-pattern/CFU result files."""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import h5py
import numpy as np


SERVER_ROOT = Path("/home/cyf/wbi/wbi_code")
REPO_ROOT = SERVER_ROOT / "experiments/motion_pattern_cfu_association"
PATTERN_ROOT = REPO_ROOT / "01_motion_pattern_extraction_omega05_mu05/patterns"
CFU_ROOT = REPO_ROOT / "02_current_cfu_input/cfu"
TEMPORAL_ROOT = REPO_ROOT / "04_all_pattern_cfu_lag_cooccurrence"
T = 1598
LAGS = np.arange(-8, 9, dtype=int)
WINDOW = 3


def configure_imports(repo_root: Path = REPO_ROOT) -> None:
    for root in (repo_root, repo_root.parents[1] if len(repo_root.parents) > 1 else repo_root):
        if root.exists() and str(root) not in sys.path:
            sys.path.insert(0, str(root))
    source = repo_root / "wholistic_registration/src"
    if source.exists() and str(source) not in sys.path:
        sys.path.insert(0, str(source))


def load_temporal_module(repo_root: Path = REPO_ROOT):
    script = repo_root / "04_all_pattern_cfu_lag_cooccurrence/run_all_cfu_pattern_lag8_w3_pairwise.py"
    configure_imports(repo_root)
    spec = importlib.util.spec_from_file_location("current_lag_analysis", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import current analysis: {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_patterns(z: int, pattern_root: Path = PATTERN_ROOT, repo_root: Path = REPO_ROOT):
    configure_imports(repo_root)
    path = pattern_root / f"Slice{z:02d}_velocity_decomp/06_patterns/objects.pkl"
    import pickle

    with path.open("rb") as f:
        return pickle.load(f)


def cfu_path(z: int, cfu_root: Path = CFU_ROOT) -> Path:
    return cfu_root / f"slice_Z{z:02d}_ds7_native_CFU_ot030_min5_group5.mat"


def _cell_layout(cells):
    if cells.ndim != 2:
        raise ValueError(f"Unexpected cfuInfo1 shape {cells.shape}")
    if cells.shape[0] in (9, 10):
        return "field_first", cells.shape[1]
    if cells.shape[1] in (9, 10):
        return "cfu_first", cells.shape[0]
    raise ValueError(f"Cannot orient cfuInfo1 shape {cells.shape}")


def read_cfu(z: int, ids=None, cfu_root: Path = CFU_ROOT):
    """Read current AQuA2 CFU mask, onset, dF/F, and member-event count."""
    path = cfu_path(z, cfu_root)
    out = {}
    with h5py.File(path, "r") as f:
        cells = f["cfuInfo1"]
        layout, n = _cell_layout(cells)
        wanted = range(1, n + 1) if ids is None else sorted(set(map(int, ids)))
        for cid in wanted:
            if cid < 1 or cid > n:
                continue

            def get(field0):
                ref = cells[field0, cid - 1] if layout == "field_first" else cells[cid - 1, field0]
                return np.asarray(f[ref])

            member_ids = get(1).reshape(-1)
            weight = get(2).T.astype(np.float32, copy=False)
            onset = get(3).reshape(-1).astype(bool, copy=False)
            dff = get(5).reshape(-1).astype(np.float32, copy=False)
            if onset.size != T or dff.size != T:
                raise ValueError(f"Slice {z} CFU {cid}: expected {T} frames, got {onset.size}/{dff.size}")
            out[cid] = {
                "weight": weight,
                "mask": weight > 0,
                "onset": onset,
                "dff": dff,
                "n_member_events": int(member_ids.size),
            }
    return out


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=str)
