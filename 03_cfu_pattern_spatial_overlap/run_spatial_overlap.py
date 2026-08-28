#!/usr/bin/env python3
"""Recompute only the motion_pattern_cfu_association spatial pattern-CFU module stage.

The implementation lives in the shared network runner; --spatial-only keeps
this stage independent of temporal FDR tables and does not write network data.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path("/home/cyf/wbi/wbi_code/experiments/motion_pattern_cfu_association")
RUNNER = ROOT / "05_local_mechanical_modules_distributed_ca_network/run_module_network.py"
subprocess.run([sys.executable, str(RUNNER), "--spatial-only"], check=True)
