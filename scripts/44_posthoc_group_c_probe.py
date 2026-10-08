#!/usr/bin/env python
"""POST-HOC probe: does the group (c) exclusion rule track predictive value?

**THIS IS POST-HOC AND CANNOT CARRY A CLAIM.** It was written after the pre-registered results
were read, because one pre-registered sensitivity arm (`qm_full_hmin`, which adds `hirshfeld_min`)
came back improving `qm_full` on 25/25 folds — larger than any primary effect, from a column the
frozen rule had EXCLUDED as "reproducible from the free 2D block at out-of-sample R^2 >= 0.90".

The rule that excluded `hirshfeld_min` (at 0.885, the one judgement call) excluded six others on
the same logic, two of them at far higher reproducibility (`hirshfeld_max` 0.993,
`hirshfeld_abs_mean` 0.951). If those also improve `qm_full`, the group (c) threshold does not
track predictive value — which is exactly the caveat the plan wrote down in Part 2.3 and Part 7
entry 13, and this is its first measurement.

Everything here is reported SEPARATELY from the pre-registered results and is barred from
changing any verdict. `prereg.json` is not touched. Writes to a separate directory so no
pre-registered output can be overwritten.

Usage:  python scripts/44_posthoc_group_c_probe.py
"""
from __future__ import annotations
import os
assert "OMP_NUM_THREADS" not in os.environ, "OMP_NUM_THREADS is set -- unset it."
import json, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
import importlib.util
_s = importlib.util.spec_from_file_location("drv", REPO_ROOT / "scripts" / "44_run_qm_arms.py")
drv = importlib.util.module_from_spec(_s); _s.loader.exec_module(drv)

OUT = REPO_ROOT / "outputs" / "44_qm_molecular_block"
PH = OUT / "posthoc_group_c"; PH.mkdir(parents=True, exist_ok=True)
FEAT = OUT / "features"

# the group (c) columns, from prereg.json's own forbidden list rather than retyped
PRE = json.loads((OUT / "prereg.json").read_text())
GROUP_C = ["hirshfeld_min", "hirshfeld_max", "hirshfeld_abs_mean", "hirshfeld_range",
           "formal_charge", "dft_energy_eh", "energy_per_atom"]
QM_FULL = PRE["blocks"]["QM_FULL"]["columns"]


def main() -> int:
    qd = pd.read_csv(REPO_ROOT / "outputs" / "qm_descriptors" / "descriptors.csv")
    qf = pd.read_csv(FEAT / "QM_FULL.csv")
    qd = qd.assign(hirshfeld_range=qd.hirshfeld_max - qd.hirshfeld_min,
                   energy_per_atom=qd.dft_energy_eh / qd.total_atoms_with_h)
    base = qf.merge(qd[["Molecule_Name"] + GROUP_C], on="Molecule_Name", validate="1:1")
    assert len(base) == 5655 and not base[GROUP_C].isna().any().any()

    blocks = {}
    for c in GROUP_C:
        blocks[f"PH_QMFULL_{c.upper()}"] = QM_FULL + [c]
    blocks["PH_QMFULL_ALLC"] = QM_FULL + GROUP_C
    written = []
    for nm, cols in blocks.items():
        df = base[["Molecule_Name", "inchikey", "split"] + cols]
        p = FEAT / f"{nm}.csv"
        assert not p.exists(), f"{p} exists -- refusing to overwrite a pre-registered block"
        df.to_csv(p, index=False)
        written.append(nm)
        print(f"  built {nm:34s} {len(cols):3d} cols")

    ctx = drv.Ctx("confirm")
    ctx.scores_path = PH / "run_scores.csv"
    ctx.timings_path = PH / "run_timings.csv"
    ctx.pred_dir = PH / "predictions"; ctx.score_dir = PH / "scores"; ctx.oof_dir = PH / "oof_long"
    for d in [ctx.pred_dir, ctx.score_dir, ctx.oof_dir]:
        d.mkdir(parents=True, exist_ok=True)
    m_star = json.loads((OUT / "model_selection.json").read_text())["m_star"]
    print(f"\n{len(written)} post-hoc arms x 25 folds at M* = {m_star}")
    t0 = time.time()
    for r in range(5):
        for f in range(drv.N_FOLDS):
            seed = ctx.fold_seed[(r, f)]
            for nm in written:
                drv.run_tabular_unit(ctx, nm, m_star, "butina", r, f, seed)
    print(f"\ndone in {(time.time()-t0)/60:.1f} min")
    (PH / "POSTHOC_README.json").write_text(json.dumps({
        "status": "POST-HOC -- written after the pre-registered results were read",
        "cannot_carry_a_claim": True,
        "why": ("the pre-registered qm_full_hmin sensitivity arm improved qm_full on 25/25 folds, "
                "from a column the frozen group (c) rule had excluded; this probes whether the "
                "rule tracks predictive value at all"),
        "arms": written, "group_c_columns": GROUP_C, "m_star": m_star,
        "prereg_untouched": True,
        "at_utc": datetime.now(timezone.utc).isoformat()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
