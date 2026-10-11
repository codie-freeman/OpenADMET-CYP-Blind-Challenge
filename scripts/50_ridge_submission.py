"""The shared implementation behind notebooks 50 and 51. ONE copy, deliberately.

THE PAIRING IS THE POINT. Notebooks 50 and 51 are a matched pair sent about 12 hours apart, and
the only thing that may differ between them is the feature block. If anything else differed, the
board difference would not be attributable to the quantum data.

The guarantee is structural rather than asserted after the fact: there is exactly one
implementation of the fit, the correction, the validation and the output format, it lives here,
and both notebooks call `build(block=...)` with nothing else changed. `build` writes a
`pipeline_fingerprint.json` recording every parameter it used, and notebook 51 asserts its
fingerprint equals notebook 50's on every field except the block name and the output paths.

THE RIDGE PIPELINE IS NOTEBOOK 46'S OWN FUNCTION, imported rather than retyped, so the alphas,
the 1st/99th clip and the scaler cannot drift from the arms whose out-of-fold numbers both
notebooks quote as reference.

WHAT IS FITTED. One ridge per endpoint, on the training compounds that carry a label for that
endpoint. That is notebook 46's `run_tabular_one` structure and it is forced by the labelling:
the four endpoints are 26-48% sparse and no single-output ridge can take a NaN target. It is NOT
per-isoform model selection -- the block, the estimator, the alpha grid and the clip are fixed
and identical across endpoints and across the two notebooks. The alpha itself is chosen inside
each fit by RidgeCV's own leave-one-out criterion on training rows only, the same way in both.

FULL-DATA SCOPE. There is no held-out portion in a full-data fit, so the clip percentiles are
computed on each endpoint's own labelled training rows. That follows notebooks 30/32/40's
full-data precedent (fit-portion statistics computed on everything available) and is the only
thing that differs from notebook 46's per-fold version, identically in both notebooks.

CORRECTION. The two layers this project has already measured on the board, applied as separate
reported layers and never silently:
  layer 1, CYP2D6 only -- placement onto notebook 16's deployed solved blind population, using
           THIS model's own pooled out-of-fold Pearson rho. The chain NB16/NB19/NB27-widened/
           NB32-A/27b/47 all use.
  layer 2, all four    -- notebook 36's widening to the spread ratio of the column holding that
           isoform's best board ST-RAE, applied only where the relative gap reaches 10%.
The procedure is identical in both notebooks. The parameters are fitted per model and so differ,
which is the point of a correction.

Nothing here sends anything. The notebooks hold the SEND flag.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
FEATURES = REPO_ROOT / "outputs" / "49_qm_only_arm" / "features_blind_extended"
NB46 = REPO_ROOT / "outputs" / "46_qm_arm"
CURATED = REPO_ROOT / "data" / "processed" / "train_inhibition_curated.csv"
BLIND = REPO_ROOT / "data" / "processed" / "test_blinded_curated.csv"
LEADERBOARD = REPO_ROOT / "docs" / "leaderboard_submissions.md"

ISOFORMS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
EP = {i: f"{i}_pIC50_direct_inhibition" for i in ISOFORMS}
# ridge is deterministic at these settings; the seed is recorded for the fingerprint and the
# determinism of the fit is checked rather than assumed.
SEED = 2684470948
WIDEN_BAR = 0.10
NB16_MU, NB16_SD = 3.155762171801848, 1.511412192589019

SENT = {
    "04b": "outputs/submissions/activity_submission_v1.csv",
    "NB10": "outputs/10_final_retrain_predict/submission_candidate.csv",
    "10c": "outputs/10c_control_submission/submission_candidate.csv",
    "NB12": "outputs/12_caruana_retrain_predict/submission_candidate.csv",
    "NB13-calib": "outputs/13_aid1851_blind_population_calibration/submission_candidate.csv",
    "NB16": "outputs/board_solved_calibration/submission_candidate.csv",
    "NB19": "outputs/19_cyp2c9_revert/submission_candidate.csv",
    "NB27-widened": "outputs/27_calibration/cyp2d6_widened/submission_candidate.csv",
    "NB30": "outputs/30_aid_full_retrain/submission_candidate.csv",
    "27b": "outputs/27b_aid_cyp2d6_corrected/submission_candidate.csv",
    "NB32-A": "outputs/32_deadzone_aid_retrain/submission_candidate_A_allfour.csv",
    "NB36-widened": "outputs/36_spread_corrected_nb32a/submission_candidate_spread_corrected.csv",
}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def nb46_module():
    """Notebook 46's runner, imported so `make_model` and its constants cannot drift."""
    spec = importlib.util.spec_from_file_location(
        "nb46_run_arms_for_submission", REPO_ROOT / "scripts" / "46_run_arms.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def board_strae() -> dict:
    md = LEADERBOARD.read_text()
    tbl = md[md.index("## All twelve submissions to date"):md.index("## UNCONFIRMED SUBMISSION")]
    lines = [l for l in tbl.splitlines() if l.strip().startswith("|")]
    subs = [c.strip().strip("*") for c in lines[0].strip("|").split("|")][2:]
    out, cur = {}, None
    for l in lines[2:]:
        cells = [c.strip().strip("*") for c in l.strip("|").split("|")]
        if cells[0]:
            cur = cells[0].replace(" (Overall)", "")
        if cells[1] != "ST-RAE":
            continue
        for s, v in zip(subs, cells[2:]):
            out[(cur, s)] = float(v.replace("−", "-"))
    for s in subs:
        assert abs(np.mean([out[(i, s)] for i in ISOFORMS]) - out[("MA", s)]) < 5e-4, s
    return out


def widening_targets(trsd: dict) -> dict:
    """Notebook 36's rule, re-derived from the sent files plus the parsed board table."""
    key = pd.read_csv(BLIND)[["Molecule_Name"]]
    board = board_strae()
    ratios = {}
    for tag, rel in SENT.items():
        d = key.merge(pd.read_csv(REPO_ROOT / rel), on="Molecule_Name", validate="one_to_one")
        for iso in ISOFORMS:
            ratios[(tag, iso)] = float(d[EP[iso]].std(ddof=1) / trsd[iso])
    target = {}
    for iso in ISOFORMS:
        best = min(board[(iso, s)] for s in SENT)
        ties = sorted(s for s in SENT if abs(board[(iso, s)] - best) < 1e-12)
        vals = np.array([ratios[(t, iso)] for t in ties])
        assert vals.max() - vals.min() < 1e-9, (iso, ties)   # rule must not depend on the tie
        target[iso] = float(vals[0])
    nb47 = pd.read_csv(REPO_ROOT / "outputs/47_mae_only_full_retrain/widening_applied.csv"
                       ).set_index("isoform")["target_ratio"]
    for iso in ISOFORMS:
        want = 0.85 if iso == "CYP2D6" else float(nb47.loc[iso])
        assert abs(target[iso] - want) < 1e-9, (iso, target[iso], want)
    return target


def oof_pearson_rho(arm: str) -> dict:
    """This model's own pooled out-of-fold Pearson rho, from notebook 46's saved columns."""
    import glob
    from scipy.stats import pearsonr
    q = pd.concat([pd.read_csv(p) for p in
                   sorted(glob.glob(str(NB46 / "oof_long" / f"{arm}__r*.csv")))],
                  ignore_index=True)
    assert set(q.arm.unique()) == {arm} and q.groupby(["repeat", "fold"]).ngroups == 25
    cur = pd.read_csv(CURATED)
    out = {}
    for iso in ISOFORMS:
        s = q[q.isoform == iso].merge(cur[["inchikey", EP[iso]]], on="inchikey", validate="many_to_one")
        s = s[s[EP[iso]].notna()]
        out[iso] = float(pearsonr(s[EP[iso]], s.y_pred).statistic)
    return out


def build(block: str, arm: str, out_dir: Path, model_tag: str) -> dict:
    """Fit, predict, correct, validate, write. The ONLY argument that differs between the two
    notebooks is `block` (and the names that follow from it)."""
    import sys
    sys.path.insert(0, str(REPO_ROOT))
    from src.vendor.validation.activity_validation import validate_activity_submission

    out_dir = Path(out_dir)
    if not out_dir.is_absolute():
        out_dir = REPO_ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    m46 = nb46_module()
    assert os.environ.get("OMP_NUM_THREADS") is None, "OMP_NUM_THREADS must be unset"

    curated = pd.read_csv(CURATED)
    blind = pd.read_csv(BLIND)
    assert len(curated) == 4905 and len(blind) == 750
    feat_path = FEATURES / f"{block}.csv"
    F = pd.read_csv(feat_path)
    cols = [c for c in F.columns if c not in ("Molecule_Name", "inchikey")]
    assert len(F) == 5655 and F[cols].notna().all().all(), f"{block}: unexpected NaN"
    tr = F[F.Molecule_Name.isin(set(curated.Molecule_Name))].reset_index(drop=True)
    bl = blind[["Molecule_Name"]].merge(F, on="Molecule_Name", validate="one_to_one")
    assert len(tr) == 4905 and len(bl) == 750
    y = curated.set_index("Molecule_Name")
    Xtr_all = tr[cols].to_numpy(float)
    Xbl = bl[cols].to_numpy(float)

    # ---- fit: one ridge per endpoint, notebook 46's own make_model, on that endpoint's labels
    raw, alphas, nfit = {}, {}, {}
    for iso in ISOFORMS:
        lab = y.loc[tr.Molecule_Name, EP[iso]].to_numpy(float)
        sel = ~np.isnan(lab)
        model = m46.make_model("ridge", SEED, Xtr_all[sel])
        model.fit(Xtr_all[sel], lab[sel])
        raw[iso] = np.asarray(model.predict(Xbl)).reshape(-1)
        alphas[iso] = float(model.named_steps["ridge"].alpha_)
        nfit[iso] = int(sel.sum())
        # determinism, checked rather than assumed: ridge takes a seed it does not use
        again = m46.make_model("ridge", SEED + 1, Xtr_all[sel])
        again.fit(Xtr_all[sel], lab[sel])
        d = float(np.abs(np.asarray(again.predict(Xbl)).reshape(-1) - raw[iso]).max())
        assert d == 0.0, f"{iso}: ridge is not deterministic ({d:.3e})"

    trsd = {i: float(curated[EP[i]].dropna().std(ddof=1)) for i in ISOFORMS}
    trmin = {i: float(curated[EP[i]].min()) for i in ISOFORMS}
    target = widening_targets(trsd)
    rho = oof_pearson_rho(arm)

    # ---- correction, two reported layers
    corr, chain, widen = {}, [], []
    for iso in ISOFORMS:
        v = raw[iso].astype(float).copy()
        m0, s0 = float(v.mean()), float(v.std(ddof=1))
        if iso == "CYP2D6":
            slope = rho[iso] * NB16_SD / s0
            v = NB16_MU + slope * (v - m0)
            m1, s1 = float(v.mean()), float(v.std(ddof=1))
            assert abs(m1 - NB16_MU) < 1e-9 and abs(s1 - rho[iso] * NB16_SD) < 1e-9
            chain.append(dict(isoform=iso, rho=rho[iso], target_mean=NB16_MU,
                              target_sd=NB16_SD, slope=slope, intercept=NB16_MU - slope * m0,
                              mean_before=m0, sd_before=s0, ratio_before=s0 / trsd[iso],
                              mean_after=m1, sd_after=s1, ratio_after=s1 / trsd[iso]))
            m0, s0 = m1, s1
        gap = (target[iso] - s0 / trsd[iso]) / target[iso]
        do = bool(gap >= WIDEN_BAR)
        mult = target[iso] * trsd[iso] / s0 if do else 1.0
        v = m0 + mult * (v - m0)
        corr[iso] = v
        widen.append(dict(isoform=iso, ratio_before=s0 / trsd[iso], target_ratio=target[iso],
                          relative_gap=gap, bar=WIDEN_BAR, applied=do, multiplier=mult,
                          mean_before=m0, mean_after=float(v.mean()), sd_before=s0,
                          sd_after=float(v.std(ddof=1)),
                          ratio_after=float(v.std(ddof=1)) / trsd[iso],
                          total_sd_multiplier=float(v.std(ddof=1)) / float(raw[iso].std(ddof=1)),
                          spearman_vs_raw=float(pd.Series(raw[iso]).corr(
                              pd.Series(v), method="spearman")),
                          min_after=float(v.min()), max_after=float(v.max()),
                          n_below_zero=int((v < 0).sum()),
                          n_below_train_min=int((v < trmin[iso]).sum())))
    W = pd.DataFrame(widen)
    for _, r in W.iterrows():
        assert abs(r.spearman_vs_raw - 1.0) < 1e-12, f"{r.isoform}: rank order moved"
        assert abs(r.ratio_after - (r.target_ratio if r.applied else r.ratio_before)) < 1e-9
    W.to_csv(out_dir / "widening_applied.csv", index=False)
    pd.DataFrame(chain).to_csv(out_dir / "cyp2d6_correction_chain.csv", index=False)

    # ---- write and validate
    paths = {}
    for tag, colset in (("raw", raw), ("corrected", corr)):
        df = pd.DataFrame({"SMILES": blind["SMILES"].to_numpy(),
                           "Molecule_Name": blind["Molecule_Name"].to_numpy()})
        for iso in ISOFORMS:
            df[EP[iso]] = colset[iso]
        p = out_dir / (f"submission_candidate_{tag}.csv" if tag == "raw"
                       else "submission_candidate.csv")
        df.to_csv(p, index=False)
        back = pd.read_csv(p)
        assert list(back.columns) == list(df.columns) and len(back) == 750
        assert np.abs(back[[EP[i] for i in ISOFORMS]].to_numpy()
                      - df[[EP[i] for i in ISOFORMS]].to_numpy()).max() < 1e-9
        ok, errs = validate_activity_submission(p, expected_ids=set(blind.Molecule_Name))
        assert ok, errs
        paths[tag] = dict(path=str(p.relative_to(REPO_ROOT)), sha256=sha256(p))

    spread = pd.DataFrame([
        dict(isoform=i, raw_mean=float(raw[i].mean()), raw_sd=float(raw[i].std(ddof=1)),
             raw_ratio=float(raw[i].std(ddof=1) / trsd[i]),
             corrected_sd=float(corr[i].std(ddof=1)),
             corrected_ratio=float(corr[i].std(ddof=1) / trsd[i]),
             target_ratio=target[i], train_sd=trsd[i]) for i in ISOFORMS])
    spread.to_csv(out_dir / "spread_ratios.csv", index=False)

    # ---- the fingerprint notebook 51 checks against notebook 50's
    fp = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "implementation": "scripts/50_ridge_submission.py",
        "implementation_sha256": sha256(Path(__file__)),
        "nb46_runner_sha256": sha256(REPO_ROOT / "scripts" / "46_run_arms.py"),
        "VARIES_BETWEEN_NOTEBOOKS": {"block": block, "arm": arm,
                                     "out_dir": str(out_dir.relative_to(REPO_ROOT)),
                                     "model_tag": model_tag, "n_features": len(cols)},
        "MUST_BE_IDENTICAL": {
            "estimator": "notebook 46 make_model('ridge'): Pipeline[clip p1-p99 on the fit rows,"
                         " StandardScaler, RidgeCV]",
            "ridge_alphas": [float(a) for a in m46.RIDGE_ALPHAS],
            "ridge_clip_pct": list(m46.RIDGE_CLIP_PCT),
            "alpha_selection": "RidgeCV internal leave-one-out, per endpoint, training rows only",
            "seed": SEED,
            "fit_rows_per_endpoint": nfit,
            "n_train": 4905, "n_blind": 750,
            "nan_fill": "none required; both blocks are complete (asserted)",
            "correction_layer1": "CYP2D6 placement onto NB16 solved population "
                                 f"mean={NB16_MU} sd={NB16_SD} using this model's own OOF rho",
            "correction_layer2": "notebook 36 widening to the best-board-column ratio, "
                                 f"applied where relative gap >= {WIDEN_BAR}",
            "widening_targets": target,
            "output_columns": ["SMILES", "Molecule_Name"] + [EP[i] for i in ISOFORMS],
            "validator": "src/vendor/validation/activity_validation.validate_activity_submission",
        },
        "FITTED_PER_MODEL": {"alpha_per_endpoint": alphas, "oof_pearson_rho": rho,
                             "widening_multiplier": {r.isoform: float(r.multiplier)
                                                     for _, r in W.iterrows()}},
        "files": paths,
    }
    json.dump(fp, open(out_dir / "pipeline_fingerprint.json", "w"), indent=1)
    return {"spread": spread, "widening": W, "chain": pd.DataFrame(chain),
            "fingerprint": fp, "raw": raw, "corrected": corr, "alphas": alphas, "rho": rho}
