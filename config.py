"""Local configuration for the motion-pattern/CFU association pipeline.

Edit the active paths below when moving the project to another machine.
All analysis scripts should derive repository-relative paths from REPO_ROOT.
"""
from pathlib import Path
import os

REPO_ROOT = Path(__file__).resolve().parent

# Current server data root.  For a collaborator, replace this with the local
# data directory, e.g. Path('/path/to/yourdata').
DATA_ROOT = Path('/home/cyf/wbi/wbi_code/data')
# DATA_ROOT = Path('/path/to/yourdata')

# Current server reference image used only for publication figures.
REFERENCE_TIF = Path('/mnt/data21T_2/cyf/f338/f338_registrated_0530/reference/vol_ref_000599_000999.tif')
# REFERENCE_TIF = Path('/path/to/your/reference/vol_ref_000599_000999.tif')

PATTERN_DIR = REPO_ROOT / '01_motion_pattern_extraction_omega05_mu05' / 'patterns'
CFU_DIR = REPO_ROOT / '02_current_cfu_input' / 'cfu'
SPATIAL_DIR = REPO_ROOT / '03_cfu_pattern_spatial_overlap'
TEMPORAL_DIR = REPO_ROOT / '04_all_pattern_cfu_lag_cooccurrence'
NETWORK_DIR = REPO_ROOT / '05_local_mechanical_modules_distributed_ca_network'

# Analysis parameters are intentionally centralized and easy to edit.
MIN_PATTERN_MEMBERS = 5
RATIO_THRESHOLD = 3.0
COVERAGE_THRESHOLD = 0.5
LAG_MIN, LAG_MAX = -8, 8
WINDOW_W = 3
