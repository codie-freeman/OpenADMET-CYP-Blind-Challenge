# Data Provenance

Source: https://huggingface.co/datasets/openadmet/cyp-challenge-train-test

- Downloaded: 2026-08-19
- Dataset commit SHA (revision pinned at download time): `85f8b358d0a2056a98b990dd75d3b3ec9247862b`
- Dataset `lastModified` (per HF API at download time): 2026-08-06 20:25:31 UTC
- Downloaded via: `huggingface_hub.hf_hub_download`, `repo_type="dataset"`, `revision="85f8b358d0a2056a98b990dd75d3b3ec9247862b"`

## Files

| File | Rows | Columns |
|---|---|---|
| cyp-challenge-TRAIN_inhibition.csv | 4905 | 18 |
| cyp-challenge-TEST-BLINDED.csv | 750 | 2 |
| cyp-challenge-TRAIN_TDI.csv | 6145 | 36 |
| cyp-challenge-single-concentration-TRAIN.csv | 17504 | 12 |
| cyp-challenge-TRAIN_Emax.csv | 6146 | 30 |

No cleaning, filtering, deduplication, or modification was applied to any of the above files — they are byte-identical to what was retrieved from the Hugging Face Hub at the pinned revision.

To reproduce this exact snapshot later (e.g. if OpenADMET revises the dataset):

```python
from huggingface_hub import hf_hub_download
hf_hub_download(
    repo_id="openadmet/cyp-challenge-train-test",
    filename="<filename>",
    repo_type="dataset",
    revision="85f8b358d0a2056a98b990dd75d3b3ec9247862b",
)
```
