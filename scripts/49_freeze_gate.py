"""Freeze the decorrelation gate BEFORE the QM-only arm is trained.

Every quantity here is a function of predictions that already existed: notebook 37's 25-fold
out-of-fold table, notebook 40's seed-averaged xgboost out-of-fold table, and notebook 46's six
arms. None of them involves the QM-only arm, which has not been trained at the time this runs,
so the threshold cannot be tuned to its result.

THE RULE. Ensemble value is decorrelation, not accuracy. `aid` is the member whose diversity the
pool currently depends on -- notebook 48 found that dropping it takes the greedy blend's gain to
exactly 0.0000 on all four isoforms -- and it is also the member the real board contradicts on
4 of 4 isoforms. So the bar a replacement has to clear is `aid`'s own decorrelation, not some
absolute number:

    PASS  if the arm's mean residual correlation against notebook 48's pool is BELOW `aid`'s
          own mean residual correlation against the rest of that pool, on at least one isoform.

    FAIL  otherwise. On a FAIL, no blind column is built and the question is closed.

A second, weaker reference is recorded for context only and is not part of the rule: the pool's
own internal mean, and the most decorrelated non-`aid` member (`mordred_pca__xgboost`), which is
already in the pool and already available for weighting.
"""

from __future__ import annotations

import glob
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "outputs" / "49_qm_only_arm"
NB37 = REPO_ROOT / "outputs" / "37_butina_5x5_complete"
NB40 = REPO_ROOT / "outputs" / "40_single_object_ensemble"
NB46 = REPO_ROOT / "outputs" / "46_qm_arm"

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
EP = {i: f"{i}_pIC50_direct_inhibition" for i in ISOFORMS}

# notebook 48's deployable pool per isoform: an arm needs an out-of-fold column AND a full-data
# blind column. Reproduced from outputs/48_board_ensemble/all_fitted_weights.csv's own pools.
POOL = {
    "CYP1A2": ["plain", "ecfp4_narrow__rf", "ecfp4_narrow__xgboost", "ecfp4_narrow__lightgbm",
               "mordred_pca__xgboost", "mordred_pca__lightgbm", "aid", "deadzone", "mae_only",
               "ecfp4_narrow__xgboost_stab"],
    "CYP2C9": ["plain", "ecfp4_narrow__rf", "ecfp4_narrow__xgboost", "ecfp4_narrow__lightgbm",
               "mordred_pca__lightgbm", "aid", "deadzone", "mae_only",
               "ecfp4_narrow__xgboost_stab"],
    "CYP2D6": ["plain", "ecfp4_narrow__rf", "ecfp4_narrow__xgboost", "ecfp4_narrow__lightgbm",
               "mordred_pca__rf", "mordred_pca__xgboost", "mordred_pca__lightgbm", "aid",
               "deadzone", "mae_only", "ecfp4_narrow__xgboost_stab"],
    "CYP3A4": ["plain", "ecfp4_narrow__rf", "ecfp4_narrow__xgboost", "ecfp4_narrow__lightgbm",
               "mordred_pca__lightgbm", "aid", "deadzone", "mae_only",
               "ecfp4_narrow__xgboost_stab"],
}
# The full 12-arm pool-A table, which is the reference the approved plan's `aid` figures were
# computed against. It is NOT the gate reference: notebook 48's ensemble is built from the
# deployable pool above, so that is what a replacement member has to beat. Both are recorded,
# the plan's values are asserted against the reference they were actually computed on, and the
# difference is reported rather than hidden -- it does not change any verdict, because every QM
# arm sits far above `aid` under either reference.
FULL_POOL_A = ["plain", "ecfp4_narrow__rf", "ecfp4_narrow__xgboost", "ecfp4_narrow__lightgbm",
               "mordred_pca__rf", "mordred_pca__xgboost", "mordred_pca__lightgbm", "single",
               "aid", "deadzone", "mae_only", "ecfp4_narrow__xgboost_stab"]
AID_EXPECTED_FULL_A = {"CYP1A2": 0.8321, "CYP2C9": 0.7866, "CYP2D6": 0.8838, "CYP3A4": 0.7396}


def merged_oof():
    """Notebook 37's 25-fold table + notebook 40's stabilised xgboost + notebook 46's six arms."""
    stab = pd.read_csv(NB40 / "part8_stabilised_member" / "xgboost_oof_10seeds.csv")
    sa = stab.groupby(["repeat", "fold", "inchikey"], as_index=False)[
        [EP[i] for i in ISOFORMS]].mean()
    q = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(str(NB46 / "oof_long" / "*.csv")))],
                  ignore_index=True)
    out = {}
    for iso in ISOFORMS:
        d = pd.read_csv(NB37 / "oof_long" / f"oof_long_{iso}_25fold.csv")
        n0 = len(d)
        d = d.merge(sa[["repeat", "fold", "inchikey", EP[iso]]].rename(
            columns={EP[iso]: "ecfp4_narrow__xgboost_stab"}),
            on=["repeat", "fold", "inchikey"], validate="one_to_one")
        w = q[q.isoform == iso].pivot_table(index=["repeat", "fold", "inchikey"],
                                            columns="arm", values="y_pred")
        d = d.merge(w.reset_index(), on=["repeat", "fold", "inchikey"], validate="one_to_one")
        assert len(d) == n0 and d.notna().all().all(), iso
        out[iso] = d
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    gate_path = OUT / "gate.json"
    assert not gate_path.exists(), (
        f"{gate_path} already exists -- refusing to rewrite a frozen gate. Notebook 44's audit "
        f"found an unrecoverable prereg provenance exactly this way.")
    # the arm must not exist yet, or the gate is not frozen before the result
    trained = sorted(glob.glob(str(OUT / "oof_long" / "qm13*")))
    assert not trained, f"QM-only predictions already exist ({len(trained)} files); gate is late"

    M = merged_oof()
    ctx, aid_bar = {}, {}
    for iso in ISOFORMS:
        d = M[iso]
        pool = POOL[iso]
        C = pd.DataFrame({a: d[a] - d["y_true"] for a in pool}).corr()
        aid_bar[iso] = float(C.loc["aid", [a for a in pool if a != "aid"]].mean())
        off = C.to_numpy()[np.triu_indices(len(pool), 1)]
        per_member = {a: float(C.loc[a, [x for x in pool if x != a]].mean()) for a in pool}
        best_non_aid = min((v, k) for k, v in per_member.items() if k != "aid")

        Cf = pd.DataFrame({a: d[a] - d["y_true"] for a in FULL_POOL_A}).corr()
        aid_full = float(Cf.loc["aid", [a for a in FULL_POOL_A if a != "aid"]].mean())
        d_exp = abs(aid_full - AID_EXPECTED_FULL_A[iso])

        ctx[iso] = {"pool_members": pool, "pool_internal_mean": float(off.mean()),
                    "aid_mean_vs_rest_DEPLOYABLE_this_is_the_threshold": aid_bar[iso],
                    "aid_mean_vs_rest_full_pool_A_the_plan_quoted_this": aid_full,
                    "members_excluded_from_deployable": sorted(set(FULL_POOL_A) - set(pool)),
                    "best_non_aid_member": best_non_aid[1],
                    "best_non_aid_value": best_non_aid[0],
                    "every_member": per_member}
        print(f"  {iso}: aid vs deployable {aid_bar[iso]:.6f} (THRESHOLD)   "
              f"vs full pool A {aid_full:.6f} (plan recorded "
              f"{AID_EXPECTED_FULL_A[iso]:.4f}, diff {d_exp:.1e})")
        print(f"          pool internal {off.mean():.4f}   best non-aid "
              f"{best_non_aid[1]} {best_non_aid[0]:.4f}   "
              f"excluded from deployable: {sorted(set(FULL_POOL_A) - set(pool))}")
        assert d_exp < 5e-5, (
            f"{iso}: the plan's `aid` figure does not reproduce on the reference it was "
            f"computed against")

    json.dump({
        "written_at_utc": datetime.now(timezone.utc).isoformat(),
        "written_by": "scripts/49_freeze_gate.py",
        "status": ("FROZEN BEFORE THE QM-ONLY ARM WAS TRAINED. Every input is a prediction that "
                   "already existed; none involves the arm under test."),
        "question": ("Can a QM-only arm carry the diversity `aid` is currently carrying in "
                     "notebook 48's ensemble pool?"),
        "metric": "mean Pearson correlation of residuals (prediction - truth), pooled over the "
                  "25 cluster-disjoint folds, against notebook 48's deployable pool",
        "rule": {
            "PASS": "mean residual correlation against the pool BELOW `aid`'s own mean against "
                    "the rest of the pool, on at least one isoform",
            "FAIL": "otherwise -- no blind column is built and the question is closed",
            "threshold_per_isoform": aid_bar,
            "threshold_reference": "`aid` against notebook 48's DEPLOYABLE pool for that isoform, which is the pool the ensemble is actually built from. The approved plan quoted `aid` against the full 12-arm pool-A table instead; both are recorded per isoform in `context_not_part_of_the_rule`, and the choice does not change any verdict.",
        },
        "context_not_part_of_the_rule": ctx,
        "arms_under_test": ["qm13__ridge", "qm13_atom4__ridge"],
        "reproduction_check_arm": "free__ridge, which must reproduce notebook 46's stored "
                                  "out-of-fold column before anything else is trusted",
    }, open(gate_path, "w"), indent=1)
    print(f"\nfrozen: {gate_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
