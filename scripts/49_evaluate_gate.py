"""Evaluate the frozen decorrelation gate against the QM-only arm's out-of-fold predictions.

Reads `outputs/49_qm_only_arm/gate.json` for the thresholds rather than recomputing them, so the
rule cannot drift between freezing and evaluation. Writes `gate_result.json` with the per-isoform
verdict and the evidence behind it.

Accuracy is reported alongside but is NOT part of the rule. Notebook 46 already measured accuracy
and the point of this arm is decorrelation: a block at half its materiality bar that disagrees
with the pool is worth more to an ensemble than one at 95% of bar that tracks `plain`.
"""

from __future__ import annotations

import glob
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "outputs" / "49_qm_only_arm"
NB37 = REPO_ROOT / "outputs" / "37_butina_5x5_complete"
NB40 = REPO_ROOT / "outputs" / "40_single_object_ensemble"
NB46 = REPO_ROOT / "outputs" / "46_qm_arm"
ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
EP = {i: f"{i}_pIC50_direct_inhibition" for i in ISOFORMS}
NEW_ARMS = ["qm13__ridge", "qm13_atom4__ridge"]

import sys  # noqa: E402
sys.path.insert(0, str(REPO_ROOT))
from src.vendor.openadmet_eval.custom_scoring_functions import (  # noqa: E402
    rae_soft_threshold_absolute_error as _st)


def strae(y, p, lo, hi):
    return float(_st(np.asarray(y, float), np.asarray(p, float),
                     y_true_upper=np.asarray(hi, float), y_true_lower=np.asarray(lo, float)))


def long_oof(oof_dir: Path, pattern: str) -> pd.DataFrame:
    files = sorted(glob.glob(str(oof_dir / pattern)))
    assert files, f"no files matching {pattern} under {oof_dir}"
    return pd.concat([pd.read_csv(p) for p in files], ignore_index=True)


def main() -> int:
    g = json.loads((OUT / "gate.json").read_text())
    thr = g["rule"]["threshold_per_isoform"]
    pool_of = {i: g["context_not_part_of_the_rule"][i]["pool_members"] for i in ISOFORMS}

    stab = pd.read_csv(NB40 / "part8_stabilised_member" / "xgboost_oof_10seeds.csv")
    sa = stab.groupby(["repeat", "fold", "inchikey"], as_index=False)[
        [EP[i] for i in ISOFORMS]].mean()
    qm_new = long_oof(OUT / "oof_long", "qm13*.csv")
    qm_nb46 = long_oof(NB46 / "oof_long", "*.csv")
    assert sorted(qm_new.arm.unique()) == sorted(NEW_ARMS), sorted(qm_new.arm.unique())

    rows, acc = [], []
    for iso in ISOFORMS:
        d = pd.read_csv(NB37 / "oof_long" / f"oof_long_{iso}_25fold.csv")
        n0 = len(d)
        d = d.merge(sa[["repeat", "fold", "inchikey", EP[iso]]].rename(
            columns={EP[iso]: "ecfp4_narrow__xgboost_stab"}),
            on=["repeat", "fold", "inchikey"], validate="one_to_one")
        for src in (qm_nb46, qm_new):
            w = src[src.isoform == iso].pivot_table(
                index=["repeat", "fold", "inchikey"], columns="arm", values="y_pred")
            d = d.merge(w.reset_index(), on=["repeat", "fold", "inchikey"], validate="one_to_one")
        assert len(d) == n0 and d.notna().all().all(), iso

        pool = pool_of[iso]
        C = pd.DataFrame({a: d[a] - d["y_true"] for a in pool + NEW_ARMS
                          + ["free__ridge", "free_qm13__ridge"]}).corr()
        for arm in NEW_ARMS + ["free__ridge", "free_qm13__ridge"]:
            v = float(C.loc[arm, pool].mean())
            rows.append(dict(isoform=iso, arm=arm, mean_vs_pool=v,
                             min_vs_pool=float(C.loc[arm, pool].min()),
                             threshold_aid=thr[iso], below_aid=bool(v < thr[iso]),
                             vs_free_ridge=float(C.loc[arm, "free__ridge"]),
                             vs_free_qm13=float(C.loc[arm, "free_qm13__ridge"]),
                             under_test=arm in NEW_ARMS))
        y, lo, hi = (d.y_true.to_numpy(float), d.conf_low.to_numpy(float),
                     d.conf_high.to_numpy(float))
        for arm in NEW_ARMS + ["free__ridge", "free_qm13__ridge", "plain", "aid", "deadzone"]:
            acc.append(dict(isoform=iso, arm=arm, st_rae=strae(y, d[arm].to_numpy(float), lo, hi),
                            spearman=float(spearmanr(y, d[arm]).statistic),
                            pred_sd=float(d[arm].std(ddof=1))))
    CORR, ACC = pd.DataFrame(rows), pd.DataFrame(acc)
    CORR.to_csv(OUT / "gate_correlations.csv", index=False)
    ACC.to_csv(OUT / "arm_accuracy.csv", index=False)

    print("### THE GATE -- mean residual correlation against notebook 48's deployable pool")
    print("### PASS requires BELOW `aid` on at least one isoform.\n")
    print(CORR[CORR.under_test][["isoform", "arm", "mean_vs_pool", "threshold_aid",
                                 "below_aid", "vs_free_ridge"]].round(4).to_string(index=False))
    passes = {a: sorted(CORR[(CORR.arm == a) & CORR.below_aid].isoform) for a in NEW_ARMS}
    overall = any(passes.values())
    print(f"\nisoforms passed, per arm: " + ", ".join(f"{a}: {v or 'none'}"
                                                      for a, v in passes.items()))
    print(f"\nGATE VERDICT: {'PASS' if overall else 'FAIL'}")

    print("\n### accuracy, reported alongside and NOT part of the rule")
    print(ACC.pivot(index="arm", columns="isoform", values="st_rae").round(4).to_string())

    json.dump({"evaluated_at_utc": datetime.now(timezone.utc).isoformat(),
               "gate_frozen_at_utc": g["written_at_utc"],
               "thresholds": thr,
               "arms_under_test": NEW_ARMS,
               "isoforms_passed_per_arm": passes,
               "verdict": "PASS" if overall else "FAIL",
               "consequence": ("blind columns may be built for the passing arm"
                               if overall else
                               "no blind column is built; the question is closed"),
               "mean_vs_pool": {a: {i: float(CORR[(CORR.arm == a) & (CORR.isoform == i)]
                                             .mean_vs_pool.iloc[0]) for i in ISOFORMS}
                                for a in NEW_ARMS},
               }, open(OUT / "gate_result.json", "w"), indent=1)
    print(f"\nwrote {OUT / 'gate_result.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
