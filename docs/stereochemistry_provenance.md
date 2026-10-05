# Stereochemistry provenance — what an unspecified centre means in this dataset

This records the primary source behind one modelling-relevant fact: that a compound with an
unspecified stereocentre in this challenge's data was **assayed as a mixture**, rather than having a
known configuration omitted from the record. Notebook `42` relies on it, and the dissertation will.

## What can and cannot be verified here

This repository's convention for verifiable quotation is the `*_quote_verification.json` pattern —
`outputs/37_butina_5x5_complete/ash_quote_verification.json` and
`outputs/39_ambiguity_and_shared_blend/kv_quote_verification.json`. Both pin a PDF by `sha256`,
re-extract the text at runtime, and emit `"status": "VERIFIED"`.

**That pattern does not apply here and is deliberately not imitated.** A Discord message has no DOI,
no stable citable URL, and no artefact that can be hashed and re-read. Emitting a `VERIFIED` status
against it would be a fabricated guarantee. The status this project records for it is
`CONSISTENT_WITH_LOCAL_DATA`, which is a weaker and accurate claim.

What *is* checkable is whether the statement agrees with this repository's own independently
computed structural data. It does, on every row — see the cross-check below, which notebook `42`
re-runs programmatically on every execution and which raises if it ever stops agreeing.

## Capture

| field | value |
|---|---|
| source | OpenADMET community Discord |
| poster | Scott Simpkins (OpenADMET) |
| message timestamp, as displayed | `01/10/2026, 19:57` |
| capture route | pasted verbatim into a Claude Code session by the repository owner, 2026-10-03 |
| cross-check re-run | `notebooks/42_orca_timing_run.ipynb`, Part 6.1 |

Two honesty notes on the timestamp. Discord renders it in the *reader's* local timezone, so it is
not an absolute UTC instant. And `01/10/2026` is ambiguous between 1 October and 10 January under
different locale conventions; it is recorded exactly as displayed. Given this project's
September–October 2026 activity, 1 October 2026 is the likely reading — stated here as an inference,
not as a fact.

## The message, verbatim

```
Scott Simpkins — 01/10/2026, 19:57
Hi @Kyrylo ! After looking back at the documentation from Enamine, I would say your assumption is
correct (assuming that Enamine's information is correct, which we more or less have to):

| OCNT Number  | Enamine Z Number | Stereochem                          |
|--------------|------------------|-------------------------------------|
| OCNT-1965934 |        Z57479431 | Racemic or presumed racemic or meso |
| OCNT-2328651 |      Z2242128693 | Single known enantiomer             |
| OCNT-2328781 |      Z1501485356 | Racemic or presumed racemic or meso |
| OCNT-2328911 |      Z2756516862 | Single known enantiomer             |
| OCNT-2395305 |      Z2972407088 | Single known enantiomer             |
| OCNT-2395406 |      Z2972818747 | Racemic or presumed racemic or meso |
```

## Cross-check against this repository's own data

Source: `outputs/41_dataset_audit_3d/per_compound_audit.csv`, computed by notebook `41` from the
curated SMILES with RDKit's `FindPotentialStereo`, entirely independently of the statement above.

| OCNT | Enamine Z | Claim | assigned centres | unassigned centres | matched partner | agrees |
|---|---|---|---:|---:|---|---|
| `OCNT-1965934` | `Z57479431` | Racemic or presumed racemic or meso | 0 | 1 | `OCNT-2328651` | yes |
| `OCNT-2328651` | `Z2242128693` | Single known enantiomer | 1 | 0 | `OCNT-1965934` | yes |
| `OCNT-2328781` | `Z1501485356` | Racemic or presumed racemic or meso | 0 | 1 | `OCNT-2328911` | yes |
| `OCNT-2328911` | `Z2756516862` | Single known enantiomer | 1 | 0 | `OCNT-2328781` | yes |
| `OCNT-2395305` | `Z2972407088` | Single known enantiomer | 1 | 0 | `OCNT-2395406` | yes |
| `OCNT-2395406` | `Z2972818747` | Racemic or presumed racemic or meso | 0 | 1 | `OCNT-2395305` | yes |

**6/6 agree exactly.** All six are in the curated **training** set. "Racemic or
presumed racemic or meso" corresponds in every case to one *unassigned* centre and zero assigned;
"single known enantiomer" to one *assigned* centre and zero unassigned.

**They form 3 matched constitutional pairs** — confirmed two independent ways, by
stereo-stripped canonical SMILES and by the InChIKey constitution block, which agree. Enamine lists
the racemate and the single enantiomer of the *same structure* as separate OCNT entries, and this
project's SMILES encode exactly that distinction.

That is the substantive point: **the unspecified/specified split in this dataset is deliberate and
load-bearing, not sloppy record-keeping.**

## What this licenses, and what it does not

**Licenses**: treating an unassigned centre in this dataset as genuinely unspecified in the source
material — the compound was assayed as a racemate or meso form — rather than as a curation defect.
Consequently, for such compounds the chemically correct representation is the *mixture*, so
conformer generation for a 3D descriptor should enumerate the unspecified centres and average rather
than silently embedding one arbitrary stereoisomer.

**Does not license**: imputing a configuration for any individual compound; nor extending the claim
to the blind test set, where no analogous statement exists. The statement covers six named compounds
and is treated as evidence about the dataset's *convention*, not about any unnamed compound.

## Prior Discord-sourced facts on this project

Recorded here for discoverability. Neither is retrospectively transcribed — this register begins
with the message above and does not reconstruct verbatim text it does not have.

| poster | date | fact | currently recorded in |
|---|---|---|---|
| Sean Colby | 2026-09-03 | CYP2D6 uses a different assay (Echo-MS vs fluorescence) and shares its test compounds with the other three isoforms rather than being separately hit-expansion-constructed | `CLAUDE.md`, notebook 04c status entry |
| Pat Walters | 2026-09-10 | CYP3A4 test-set construction: 2048-bit Morgan radius=3, Tanimoto *distance* cutoff 0.65, pairing each potent hit with a cluster-mate | `CLAUDE.md`, notebook 33 status entry; `notebooks/33_butina_cluster_split.ipynb` |

## Reproduce the cross-check

```python
import pandas as pd
a = pd.read_csv("outputs/41_dataset_audit_3d/per_compound_audit.csv")
ids = ["OCNT-1965934", "OCNT-2328651", "OCNT-2328781",
       "OCNT-2328911", "OCNT-2395305", "OCNT-2395406"]
print(a[a.Molecule_Name.isin(ids)][
    ["set", "Molecule_Name", "n_stereocentres_assigned",
     "n_stereocentres_unassigned", "SMILES"]].to_string(index=False))
```
