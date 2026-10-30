# Pre-registration: predictions to be checked against the released full-set labels

**Written and committed 2026-10-30, before the OpenADMET CYP450 blind labels are released.**
Submissions close **2026-11-03**. Nothing in this file may be edited after the labels are
available; outcomes are to be recorded in a separate, clearly-dated section appended below, with
the original predictions left exactly as written. A prediction written after the reveal is worth
nothing, which is the entire reason this file exists and is committed now.

**Every board figure quoted as a starting point is a half-set measurement.** The live leaderboard
scores only half of the blind test set as submissions arrive. One full-set reading exists — the
interim board of 2026-09-21, which scored `27b` — and it disagrees with the live half materially
and in opposite directions by isoform (`docs/leaderboard_submissions.md`, "The interim (full-set)
leaderboard"). That disagreement is the main reason any of these predictions could fail.

## The object under test

Two submissions, **NB56b-s050** and **NB56b-s100**, that differ in **exactly one scalar**: the
shrinkage fraction `s` applied to the affine placement-and-scale correction on CYP1A2 and CYP2C9.
Same base file (NB36-widened), same solved blind populations, same recovered ρ, same code path.
CYP2D6 and CYP3A4 are carried byte-identical into both.

| | file | sha256 |
|---|---|---|
| s = 0.5 | `outputs/56b_cyp1a2_corrected/submission_candidate_s050.csv` | `a1da00e209975f4b82fdaeae1bc0210eb83953ef13872d22aa04d354c438aa8e` |
| s = 1.0 | `outputs/56b_cyp1a2_corrected/submission_candidate_s100.csv` | `b2c53f136f1d2684643126492a98d2b718b696e55b216a17d60bf7ff6f79a228` |

Both hashes are verified against the files on disk and against each send's own server receipt.
**Checking any prediction below requires these exact files**, not rebuilt ones: the project has a
documented float round-trip hazard and 13 rows whose CXSMILES field contains a comma, so a
regenerated file may not be byte-identical.

**Half-set readings, as the starting point** (`docs/leaderboard_submissions.md`):

| | CYP1A2 | CYP2C9 | CYP2D6 | CYP3A4 | macro |
|---|---:|---:|---:|---:|---:|
| ST-RAE, s = 0.5 | 0.5631 | 0.5408 | 0.7502 | 0.4904 | 0.5861 |
| ST-RAE, s = 1.0 | 0.5767 | 0.5562 | 0.7502 | 0.4904 | 0.5933 |
| R², s = 0.5 | 0.5070 | 0.5595 | 0.3768 | 0.6563 | 0.5249 |
| R², s = 1.0 | 0.5467 | 0.5667 | 0.3768 | 0.6563 | 0.5366 |
| MAE, s = 0.5 | 0.7988 | 0.5212 | 0.9319 | 0.5373 | 0.6973 |
| MAE, s = 1.0 | 0.7575 | 0.5190 | 0.9319 | 0.5373 | 0.6864 |
| Spearman, both | 0.7565 | 0.7554 | 0.4438 | 0.8020 | 0.6894 |
| Kendall, both | 0.5611 | 0.5652 | 0.3076 | 0.6132 | 0.5118 |

## How to check these

With labels in hand, score **both** files on the full 750 compounds using the unmodified vendored
evaluator (`src/vendor/openadmet_eval/`), the same code every board comparison in this project has
used. Both files can be scored locally once labels exist, so these predictions do **not** depend on
both having been sent — only one will be standing at the close.

---

## Prediction B — rank order does not move

**Claim.** Scored on the full 750, **NB56b-s050 and NB56b-s100 return identical Spearman and
identical Kendall on all four isoforms and at macro level.**

**Pass.** All eight isoform cells (4 isoforms × 2 rank metrics) plus both macro cells agree to
`|Δ| ≤ 1e-12` at full float precision — equivalently, exact equality at the four decimal places
the board displays.

**Fail.** Any one of those ten cells differing by more than `1e-12`.

**This is a construction guarantee, not an estimate, and the construction was verified on the
actual files before the labels existed.** The correction is a strictly increasing affine map, so
the two files induce the same ordering on every subset of compounds, the scored half included.
Checked directly on the sent files:

- Spearman between the s = 0.5 and s = 1.0 columns is **1.000000000000000** on all four isoforms.
- The **rank vectors are exactly equal** on all four isoforms (max absolute rank difference 0.0).
- **Zero tied values** in any column of either file.
- The map is strictly increasing with a minimum step of **5.4e-06** (CYP1A2) and **3.5e-06**
  (CYP2C9), against a written precision of 1e-07 — so rounding cannot create a tie or an inversion.

**So this prediction should be regarded as near-certain, and that is stated rather than dressed up
as a bold forecast.** It is still worth registering, because it is the project's central
methodological claim — an affine correction cannot change which compound ranks above which — and
because it is falsifiable in three specific ways, each of which would be a real finding:

1. The files scored are not the files named above (wrong archive, or a rebuild).
2. The evaluator's Spearman is not a pure rank statistic, or handles near-ties by value rather
   than by rank.
3. The full-set scoring applies some per-submission transform before scoring.

**What a pass licenses.** Any claim in this project that a placement or spread correction improved
*ranking* is false by construction, confirmed on real labels. It follows that the ordering quality
of every affine descendant in this record is a property of the model underneath it — which is what
makes finding 3's board-versus-cross-validation transfer disagreement a statement about models and
not about calibration.

**What it does not license.** Nothing about whether either correction improved the submission
overall (that is prediction C), and nothing about rank order under a non-affine change. Both
ensembles in this project's record (NB40-ensemble, NB48-ensemble) did move Spearman, as they must.

### B2 — the byte-identical control

**Claim.** On the full set, **all ten CYP2D6 and CYP3A4 cells** (2 isoforms × 5 metrics) are
identical between the two submissions.

**Pass.** `|Δ| ≤ 1e-12` on all ten. **Fail.** Any cell differing by more than that.

The two columns are byte-identical between the files (max absolute difference `0.00e+00`, verified
directly), so a failure here means the scoring is not a deterministic function of the submitted
column — which would invalidate every resend comparison in this project's record, of which there
are seven.

---

## Prediction C — ST-RAE and R² disagree about which submission is better

**Claim.** On the two isoforms the lever moved, going from s = 0.5 to s = 1.0 on the full set,
**ST-RAE gets worse while R² gets better** — on both isoforms, as it did on the half.

**C1, primary — the four signs.** On CYP1A2 and CYP2C9:

- `ST-RAE(s=1.0) > ST-RAE(s=0.5)` (worse), and
- `R²(s=1.0) > R²(s=0.5)` (better).

**Pass.** All four cells hold, each with `|Δ| ≥ 0.001` so that a null is not scored as a pass.
**Partial.** Two or three of four. **Fail.** One or none.

**This one is a genuine risk and is not a safe bet.** The half-set ST-RAE deltas are **+0.0136**
(CYP1A2) and **+0.0154** (CYP2C9). The only half-to-full comparison this project has — `27b` on the
interim board — shows per-isoform ST-RAE *level* shifts of −0.0510, −0.0258, **+0.0797** and
−0.0084, so the largest observed half-to-full movement is five times the gap being predicted. A
sign flip on ST-RAE is entirely possible. R² has more margin on CYP1A2 (+0.0397 against observed
level shifts of at most 0.0587) and very little on CYP2C9 (+0.0072).

**C2, secondary — MAE moves with R², not with ST-RAE.** `MAE(s=1.0) < MAE(s=0.5)` on both isoforms.
**Pass** on both; this is the claim that makes ST-RAE the odd one out rather than merely one of two
metrics disagreeing, and it is the form of the finding worth carrying.

**C3, secondary — macro-level divergence.** Macro ST-RAE worsens and macro R² improves from
s = 0.5 to s = 1.0.

**C4, secondary and scale-free — the divergence is asymmetric on CYP1A2.** `|ΔR²| > |ΔST-RAE|` on
CYP1A2 (half-set ratio 0.0397 / 0.0136 ≈ 2.9). Registered because it survives a uniform
rescaling of either metric, where C1's magnitudes do not.

**C5, weak — magnitudes.** Each of the four C1 deltas stays within **±0.0797** (ST-RAE) and
**±0.0587** (R²) of its half-set value. These bounds are the largest half-to-full *level* shifts
actually observed on this board (`27b`, CYP2D6 on each metric), borrowed for want of anything
tighter. They are deliberately generous and **wide enough to contain zero**, so C5 can pass while
C1 fails; it is recorded as the weakest item here and should not be quoted as support for C1.

**C6 — the mechanism, which is the part worth being wrong about.** The proposed explanation for C1
is that ST-RAE scores a prediction as zero error anywhere inside its credible interval and as
distance to the nearest edge outside it; those intervals are wide at low pIC50 and narrow at high,
so over-prediction is cheap and under-prediction is punished, and a correction that shifts
predictions *down* onto a low-centred blind population moves mass from the cheap side to the
expensive side. R² and MAE carry no such asymmetry. If that is right, the ST-RAE-optimal uniform
offset should sit **above** the R²-optimal one.

**Claim.** Computing both optima per isoform on the full set, by sweeping a uniform additive
offset applied to the s = 0.5 column:

- the ST-RAE-optimal offset is **greater** than the R²-optimal offset on **all four** isoforms, and
- the gap lies in **[0.02, 0.25] pIC50** on at least three of four.

**Pass.** Sign holds on 4 of 4 **and** the gap is in range on ≥ 3 of 4. **Partial.** Sign on 4 of 4,
gap in range on fewer. **Fail.** Sign fails on any isoform.

The predicted range is anchored on notebook 21, which measured this gap at **+0.06 to +0.12 pIC50**
on all four isoforms using this project's own out-of-fold predictions
(`outputs/21_strae_offset_curve/offset_curve.csv`); the registered window is widened around that
because notebook 21's measurement is out-of-fold and notebook 21 itself established that out-of-fold
optima need not transfer to the blind population. **C6 is the only prediction here that would, if it
passes, explain C1 rather than merely restate it** — and the only one whose failure would leave C1
standing as an unexplained observation.

---

## What none of this tests

- **The base model.** Three of NB36-widened's four columns are unchanged in s = 0.5 and two in
  s = 1.0, so neither submission says anything about whether `deadzone + aid` was the right model.
- **The optimum in `s`.** It was sampled at two points. The ST-RAE optimum lies at or below 0.5 and
  the R² optimum at or above 1.0; neither is bracketed, and nothing here locates either.
- **Whether choosing s = 0.5 as the final submission was right.** It wins on ST-RAE and loses on R²
  and MAE on the half set, so the choice is conditional on ST-RAE being the metric that decides.
  **That condition is assumed, not confirmed**: ST-RAE is first in the vendored evaluator's
  `ACTIVITY_METRICS` and the exemplar `primary_metric` in its docstring, but no default is set
  there, the live backend's sort key is not in this repository, and the recorded board ranks are
  not monotone in macro ST-RAE across dates. A full-set reveal in which C1 flips, **or** a
  confirmation that the board ranks on R² or MAE, would each mean the choice was wrong — for
  different reasons — and either should be recorded plainly.

## Provenance

Board readings from `docs/leaderboard_submissions.md`. Half-to-full shift figures from its
"interim (full-set) leaderboard" section, which is a **user-supplied board capture not derivable
from anything in this repository**. The offset-gap figure is this project's own (notebook 21). The
pre-release file checks reported under prediction B were run on 2026-10-30 against the two sha256
files named above. Full write-up of the half-set comparison is `docs/findings.md` finding 5.
