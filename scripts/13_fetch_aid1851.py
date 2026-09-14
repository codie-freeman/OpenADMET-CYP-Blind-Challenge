"""Fetch PubChem AID 1851 (NCATS qHTS cytochrome panel) full per-SID assay data.

Batched under the PUG REST 10,000-SID cap, resumable via on-disk per-batch
CSV caching (a partial run skips already-fetched batches). Uses the full
(non-"concise") CSV export so the continuous Fit_LogAC50 dose-response fit
value is present, not just the binary active/inactive outcome.

This is a calibration-target fetch only (see CLAUDE.md / notebook 13) — the
data is never used as training/fine-tuning signal.
"""

import json
import sys
import time
from pathlib import Path

import requests

AID = 1851
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "raw" / "aid1851"
BATCH_SIZE = 2000
REQUEST_PAUSE_S = 1.0
MAX_RETRIES = 3
TIMEOUT_S = 120

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def fetch_sid_list() -> list[int]:
    sids_path = OUT_DIR / "sids.json"
    if sids_path.exists():
        sids = json.loads(sids_path.read_text())
        print(f"[sids] using cached list: {len(sids)} SIDs ({sids_path})")
        return sids

    url = f"{PUG_BASE}/assay/aid/{AID}/sids/JSON?sid_type=all"
    resp = requests.get(url, timeout=TIMEOUT_S)
    resp.raise_for_status()
    sids = resp.json()["InformationList"]["Information"][0]["SID"]
    sids_path.write_text(json.dumps(sids))
    print(f"[sids] fetched fresh: {len(sids)} SIDs -> {sids_path}")
    return sids


def fetch_batch_with_retry(sid_batch: list[int]) -> str:
    url = f"{PUG_BASE}/assay/aid/{AID}/CSV"
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.post(
                url,
                data={"sid": ",".join(str(s) for s in sid_batch)},
                timeout=TIMEOUT_S,
            )
            resp.raise_for_status()
            return resp.text
        except requests.RequestException as e:
            last_err = e
            wait = 5 * attempt
            print(f"    attempt {attempt}/{MAX_RETRIES} failed ({e}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"batch fetch failed after {MAX_RETRIES} attempts: {last_err}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sids = fetch_sid_list()

    batches = [sids[i : i + BATCH_SIZE] for i in range(0, len(sids), BATCH_SIZE)]
    print(f"[batches] {len(batches)} batches of up to {BATCH_SIZE} SIDs each")

    for i, batch in enumerate(batches):
        batch_path = OUT_DIR / f"batch_{i:04d}.csv"
        if batch_path.exists():
            print(f"[batch {i:04d}] already fetched, skipping ({batch_path})")
            continue
        print(f"[batch {i:04d}] fetching {len(batch)} SIDs...")
        csv_text = fetch_batch_with_retry(batch)
        batch_path.write_text(csv_text)
        n_lines = csv_text.count("\n")
        print(f"[batch {i:04d}] saved -> {batch_path} ({n_lines} lines)")
        time.sleep(REQUEST_PAUSE_S)

    done_marker = OUT_DIR / "_FETCH_COMPLETE"
    done_marker.write_text(f"{len(batches)} batches, {len(sids)} SIDs, AID {AID}\n")
    print(f"[done] all {len(batches)} batches present -> {done_marker}")


if __name__ == "__main__":
    sys.exit(main())
