"""The QM-only arm: 25 cluster-disjoint folds, ridge, no free descriptors.

Notebook 46's runner is IMPORTED and reused, not copied. Only four things are overridden: the
feature directory, the three output directories, and the NaN-fill policy lookup (neither new
block has a NaN, and notebook 46's one NaN-carrying block is not run here). The fitting path --
`run_tabular_one`'s per-fold median fill, `make_model`'s 1st/99th clip plus StandardScaler plus
RidgeCV, `Ctx.population`'s `assign_screen_split`, the content-based cache, the bootstrap seeding
and the anti-clobber timing rule -- is notebook 46's own code, executing unmodified.

THREE ARMS, all ridge on the frozen Butina 25-fold partition:

    free__ridge        the reproduction check. Its feature block is byte-identical to notebook
                       46's and its seeds come from the same manifest, so its out-of-fold column
                       must reproduce notebook 46's stored one. Nothing else is trusted until
                       it does.
    qm13__ridge        13 molecule-level QM columns, alone
    qm13_atom4__ridge  + notebook 45's four atom descriptors, 17 columns, alone

No blind prediction is produced here. `outputs/49_qm_only_arm/gate.json` has to pass first, and
the gate was frozen before this script was ever run.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "outputs" / "49_qm_only_arm"
NB46 = REPO_ROOT / "outputs" / "46_qm_arm"

# Notebook 34 established OMP_NUM_THREADS changes a training trajectory. Nothing here trains
# chemprop, but the environment is held identical to the one every stored score was measured
# under. KMP_DUPLICATE_LIB_OK is set by `import sklearn` in this environment, so it is reported.
assert os.environ.get("OMP_NUM_THREADS") is None, "OMP_NUM_THREADS must be unset"
_KMP = os.environ.get("KMP_DUPLICATE_LIB_OK")

BLOCKS = ["FREE", "QM13", "QM13_ATOM4"]
MSTAR = "ridge"


def load_nb46_runner():
    spec = importlib.util.spec_from_file_location(
        "nb46_run_arms", REPO_ROOT / "scripts" / "46_run_arms.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    # repoint at this script's own tree; Ctx reads these module globals when constructed
    m.FEAT_DIR = OUT / "features"
    m.OUT = OUT
    m.PRED_DIR = OUT / "predictions"
    m.SCORE_DIR = OUT / "scores"
    m.OOF_DIR = OUT / "oof_long"
    m.T1 = OUT / "tier1_unused"
    # notebook 46's prereg is consulted only for the NaN-fill block list; neither new block has
    # a NaN and FREE_ATOM4_HAL is not run here, so the list is empty by construction.
    m.prereg = lambda: {"nan_fill_policy": {"blocks": []}}
    return m


def main() -> int:
    gate = OUT / "gate.json"
    assert gate.exists(), "the gate must be frozen before this arm is trained"
    g = json.loads(gate.read_text())
    print(f"gate frozen {g['written_at_utc']}, thresholds "
          + ", ".join(f"{k} {v:.4f}" for k, v in g["rule"]["threshold_per_isoform"].items()))

    m = load_nb46_runner()
    ctx = m.Ctx("confirm")
    print(f"\nnotebook 46's runner imported and repointed; M* = {MSTAR}, "
          f"{len(BLOCKS)} blocks x 25 folds")

    t0, n = time.time(), 0
    for block in BLOCKS:
        for repeat in range(5):
            for fold in range(5):
                seed = ctx.fold_seed[(repeat, fold)]
                w = m.run_unit(ctx, block, MSTAR, "butina", repeat, fold, seed)
                n += 1
                print(f"  [{n:3d}/{len(BLOCKS) * 25}] {m.arm_name(block, MSTAR):20s} "
                      f"r{repeat}f{fold}  {w:5.2f} min")
    wall = (time.time() - t0) / 60
    print(f"\n{n} units, {wall:.2f} min total")

    # ---- reproduction check: free__ridge must reproduce notebook 46's stored out-of-fold column
    print("\n--- reproduction check: free__ridge against notebook 46's stored predictions ---")
    worst, checked = 0.0, 0
    for repeat in range(5):
        for fold in range(5):
            seed = ctx.fold_seed[(repeat, fold)]
            tag = f"free__ridge__r{repeat}f{fold}"
            a = pd.read_csv(OUT / "oof_long" / f"{tag}.csv")
            b = pd.read_csv(NB46 / "oof_long" / f"{tag}.csv")
            k = ["arm", "repeat", "fold", "seed", "inchikey", "isoform"]
            j = a.merge(b, on=k, suffixes=("_new", "_nb46"), validate="one_to_one")
            assert len(j) == len(a) == len(b), tag
            worst = max(worst, float(np.abs(j.y_pred_new - j.y_pred_nb46).max()))
            checked += len(j)
    print(f"  {checked} (fold, compound, isoform) predictions compared, "
          f"max abs diff {worst:.3e}")
    assert worst == 0.0, "free__ridge does not reproduce notebook 46 bit-identically"
    print("  BIT-IDENTICAL -- the imported runner is notebook 46's, executing unmodified.")

    json.dump({"ran_at_utc": datetime.now(timezone.utc).isoformat(),
               "blocks": BLOCKS, "mstar": MSTAR, "n_units": n, "wall_min": wall,
               "free_ridge_reproduces_nb46_max_abs_diff": worst,
               "n_predictions_compared": checked,
               "KMP_DUPLICATE_LIB_OK_at_start": _KMP,
               "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
              open(OUT / "run_meta.json", "w"), indent=1)
    print(f"\nwrote {OUT / 'run_meta.json'}")
    print(f"KMP_DUPLICATE_LIB_OK at start: {_KMP!r} (reported, not asserted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
