# Provenance

Source repo: https://github.com/OpenADMET/CYP-Challenge-Tutorial
Source directory: `validation/`
Commit hash: `9d4925eb4a0fb914256da1b27d110593bcbe3cf0`
  (`main` branch HEAD of the local clone at
  `/Users/codiefreeman/CYP-tutorial/CYP-Challenge-Tutorial` at time of pull -- a single
  "Initial commit"; the source repo's history does not currently contain the
  `858ae63c...` commit referenced by `src/vendor/openadmet_eval/PROVENANCE.md`, so that
  commit is not usable as a pin here)
Date pulled: 2026-08-26
License: Apache License 2.0 (as declared by the source repository)

## Files vendored (unmodified, byte-for-byte)

| File | Size (bytes) | Blob SHA (`git hash-object`) |
|---|---|---|
| `__init__.py` | 196 | `ab5934ebd950cc6774fbe519c0e1b10125632db5` |
| `activity_validation.py` | 3051 | `27bbcb2639f2e79abf984f58aaef810f77010dd6` |
| `tdi_validation.py` | 4046 | `a0eaa08e105b1a4d2239b564bb65022e2abf2699` |

Copied directly from the local clone above and verified byte-identical via
`git hash-object` against the clone's own tracked blobs before being placed here.
`tdi_validation.py` is vendored alongside `activity_validation.py` only because
`__init__.py` imports both (`from .tdi_validation import validate_tdi_submission`) --
this project's activity-track submission code only ever calls
`validate_activity_submission`.

## Why vendored here instead of referenced via sys.path to the external clone

`scripts/train_final_submission_multitask.py` needs
`validate_activity_submission` per the tutorial's own import pattern
(`from validation.activity_validation import validate_activity_submission`) rather than
a reimplementation. Vendoring here (matching the existing `src/vendor/openadmet_eval/`
convention) keeps that import reproducible from a fresh clone of *this* repo, rather
than depending on an external clone's absolute path that exists only on this one
machine. `sys.path` is pointed at `src/vendor/` (the parent of this `validation/`
package) so the import line itself is byte-identical to the tutorial notebook's own.

## Pre-vendoring review notes

`validate_activity_submission` (the only function this project calls) requires
`SMILES` and `Molecule_Name` id columns plus one `{CYP}_pIC50_direct_inhibition`
column per isoform in `ACTIVITY_VALUE_COLUMNS`; when `expected_ids` is passed it checks
the submitted `Molecule_Name` set against it exactly (missing/extra), otherwise it
falls back to a fixed `ACTIVITY_DATASET_SIZE = 750` row-count check. No test-set
labels are read or required -- it validates submission *format*, not prediction
accuracy. No modification needed for this project's use.
