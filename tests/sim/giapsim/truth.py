"""Readers of the M7-12 truth tables shared by the calibration tests (giapsim.scenario_truth).

Each returns a pandas DataFrame with the file's columns in order (strings), integer coordinates
and counts, and ``rbs_planted`` as bool; a header-only file gives an empty frame with the columns.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

__all__ = ["read_alt_sites", "read_indel_sites", "read_novel_orfs"]


def _read(path: Path, ints: Sequence[str], bools: Sequence[str] = ()) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ints:
        if c in df.columns:
            df[c] = df[c].astype("int64")
    for c in bools:
        if c in df.columns:
            df[c] = df[c].str.strip().str.lower().isin(("true", "1")).astype(bool)
    return df


def read_alt_sites(outdir: Path | str) -> pd.DataFrame:
    """``truth/alt_sites.tsv``: one row per non-synonymous SNV codon on the reference."""
    return _read(Path(outdir) / "truth" / "alt_sites.tsv", ("codon_start", "n_true_peptides"))


def read_indel_sites(outdir: Path | str) -> pd.DataFrame:
    """``truth/indel_sites.tsv``: in-frame indels inside genes (lead decision 2)."""
    return _read(Path(outdir) / "truth" / "indel_sites.tsv", ("pos_sample", "pos_reference"))


def read_novel_orfs(outdir: Path | str) -> pd.DataFrame:
    """``truth/novel_orfs.tsv`` with the M7-12 columns (``stop``, ``true_start``, ...)."""
    ints = ("start", "end", "stop", "true_start", "n_true_peptides")
    return _read(Path(outdir) / "truth" / "novel_orfs.tsv", ints, ("rbs_planted",))
