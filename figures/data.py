"""Read result tables out of ``outputs/``.

Figures read what the notebooks already wrote. This module is the only door to
disk in the figure layer, it never writes, and it never imports from a notebook.

Directory layout, surveyed across all 10,875 CSVs under ``outputs/``:

* 527 summary tables sit at ``outputs/<nn>_<slug>/<name>.csv`` -- the convention
  documented in ``outputs/README.md``.
* Nine top-level directories are not ``<nn>_``-prefixed (``population_moments/``,
  ``qm_parsed/``, ``qm_descriptors/``, ``board_solved_calibration/``,
  ``calibration_blind_population/``, ``final_submission/``, ``mordred_pca/``,
  ``submissions/``) and one is prefixed but is not a notebook
  (``55_superseded_preserved/``). ``notebook_dir`` is taken verbatim, so all of
  them are reachable.
* The remainder are per-run artefacts one level deeper (``predictions/``,
  ``scores/``, ``chemprop_runs/``, ``oof_long/``, ``tier1/`` ...). Pass a
  subdirectory in ``notebook_dir`` to reach those.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from figures.style.base import repo_root


def path_to(notebook_dir: str, filename: str) -> Path:
    """Resolve ``outputs/<notebook_dir>/<filename>`` without reading it."""
    return repo_root() / "outputs" / notebook_dir / filename


def load(notebook_dir: str, filename: str, **read_csv_kwargs) -> pd.DataFrame:
    """Load a result CSV from ``outputs/<notebook_dir>/<filename>``.

    ``float_precision="round_trip"`` is not decoration. Notebook 56b established
    that pandas' default C float parser is not correctly rounded: a plain
    read/write round-trip of a submission column perturbed 25 of 750 values by up
    to 8.9e-16, and two checks disagreeing is what exposed it. Every table here is
    numeric result data, so it is read exactly.

    No ``index_col`` is assumed -- ``submitted_means.csv``'s index header is blank,
    and guessing would be wrong elsewhere. Callers that want an index pass
    ``index_col=0``.
    """
    path = path_to(notebook_dir, filename)
    if not path.is_file():
        raise FileNotFoundError(
            f"no such result table: {path}\n"
            f"  (notebook_dir={notebook_dir!r}, filename={filename!r})"
        )
    read_csv_kwargs.setdefault("float_precision", "round_trip")
    return pd.read_csv(path, **read_csv_kwargs)
