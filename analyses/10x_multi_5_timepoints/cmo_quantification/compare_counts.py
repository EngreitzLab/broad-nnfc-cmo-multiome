#!/usr/bin/env python3
"""Compare new kb count output against the reference .bak files."""

import sys
import anndata
import numpy as np
from pathlib import Path

CMO_COUNTS = Path("/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts")
REFERENCE = Path(__file__).resolve().parents[3] / "tmp"
CHANNELS = ["channel1", "channel2"]

all_ok = True
for ch in CHANNELS:
    base = CMO_COUNTS / ch / "counts_unfiltered"
    old_path = REFERENCE / ch / "adata.h5ad"
    new_path = base / "adata.h5ad"

    print(f"\n=== {ch} ===")

    if not old_path.exists():
        print(f"  SKIP: no .bak file at {old_path}")
        continue
    if not new_path.exists():
        print(f"  SKIP: new file not yet written at {new_path}")
        continue

    old = anndata.read_h5ad(old_path)
    new = anndata.read_h5ad(new_path)

    shape_ok = old.shape == new.shape
    print(f"  shape  old={old.shape}  new={new.shape}  {'OK' if shape_ok else 'MISMATCH'}")

    obs_ok = list(old.obs_names) == list(new.obs_names)
    print(f"  obs    {'OK' if obs_ok else 'MISMATCH'}")

    var_ok = list(old.var_names) == list(new.var_names)
    print(f"  var    {'OK' if var_ok else 'MISMATCH'}")

    if shape_ok:
        diff = (old.X - new.X)
        counts_ok = diff.nnz == 0
        print(f"  counts {'IDENTICAL' if counts_ok else f'DIFFER: {diff.nnz} non-zero differences, max abs={np.abs(diff).max():.0f}'}")
    else:
        counts_ok = False
        print("  counts SKIP (shape mismatch)")

    if all([shape_ok, obs_ok, var_ok, counts_ok]):
        print(f"  PASS")
    else:
        print(f"  FAIL")
        all_ok = False

print(f"\n{'ALL CHANNELS PASS' if all_ok else 'SOME CHANNELS FAILED'}")
sys.exit(0 if all_ok else 1)
