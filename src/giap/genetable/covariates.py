"""Protein covariates per gene (docs/specs/M4-10-protein-covariates.md; spec §7.3.3, §7.3.4, G6).

Six covariates per gene: ``length_aa``, ``n_observable_peptides`` (distinct in-silico peptides of
the configured enzyme, 0 missed cleavages, length in ``abundance.observable_len``), ``gravy``,
``pi`` and ``mw_da`` (Biopython ``ProteinAnalysis`` on the cleaned sequence) and ``tm_segments``
(Kyte-Doolittle sliding windows, G6 proxy). Pseudogenes without an annotated amino-acid sequence
are translated from ``nt`` (table 11, first codon read as M, stop before the first stop codon).

Lead decisions: a protein with nothing left after cleaning has ``length_aa = 0``, NaN ``gravy``,
``pi`` and ``mw_da`` and 0 peptides / TM segments; letters outside the 20 standard amino acids
score 0 in TM windows; leading and trailing ``*`` are stripped; ``covariate_table`` keeps the
input row order. Everything here is pure.
"""

from __future__ import annotations

import math
import re
from dataclasses import astuple, dataclass, fields
from typing import TYPE_CHECKING, Literal, get_args

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqUtils.ProtParam import ProteinAnalysis

if TYPE_CHECKING:
    from giap.config import GiapConfig

__all__ = [
    "Covariates",
    "Enzyme",
    "covariate_table",
    "covariates",
    "digest",
    "n_observable_peptides",
    "protein_sequence",
    "tm_segments",
]

Enzyme = Literal["trypsin_p", "trypsin"]
_ENZYMES: frozenset[str] = frozenset(get_args(Enzyme))
_CLEAVE: dict[str, re.Pattern[str]] = {
    "trypsin_p": re.compile(r"(?<=[KR])"),  # after K/R, also before P (DIA-NN default)
    "trypsin": re.compile(r"(?<=[KR])(?!P)"),  # after K/R, not before P
}

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
_NON_STANDARD = re.compile(f"[^{STANDARD_AA}]")
_VALID = re.compile(r"[A-Z*]*")

# Kyte & Doolittle (1982) hydropathy, scaled by 10 so window sums are exact integers.
_KD10: dict[str, int] = {
    "A": 18, "R": -45, "N": -35, "D": -35, "C": 25, "Q": -35, "E": -35, "G": -4, "H": -32,
    "I": 45, "L": 38, "K": -39, "M": 19, "F": 28, "P": -16, "S": -8, "T": -7, "W": -9,
    "Y": -13, "V": 42,
}  # fmt: skip
_KD10_LUT = np.zeros(256, dtype=np.int64)
for _aa, _v in _KD10.items():
    _KD10_LUT[ord(_aa)] = _v

COVARIATE_COLUMNS: tuple[str, ...] = (
    "length_aa",
    "n_observable_peptides",
    "gravy",
    "pi",
    "mw_da",
    "tm_segments",
)


@dataclass(frozen=True, slots=True)
class Covariates:
    """Protein covariates of one gene (spec §7.3.3, §9.3)."""

    length_aa: int
    n_observable_peptides: int
    gravy: float
    pi: float
    mw_da: float
    tm_segments: int


def _check_enzyme(enzyme: str) -> None:
    if enzyme not in _ENZYMES:
        raise ValueError(f"unknown enzyme {enzyme!r}; expected one of {sorted(_ENZYMES)}")


def digest(aa: str, enzyme: Enzyme, missed: int = 0) -> list[str]:
    """In-silico digest: fragments in sequence order, then concatenations of 2..missed+1."""
    _check_enzyme(enzyme)
    if missed < 0:
        raise ValueError(f"missed must be >= 0, got {missed}")
    frags = [f for f in _CLEAVE[enzyme].split(aa) if f]
    out = list(frags)
    for span in range(2, missed + 2):
        out.extend("".join(frags[i : i + span]) for i in range(len(frags) - span + 1))
    return out


def n_observable_peptides(aa: str, enzyme: Enzyme, length: tuple[int, int] = (7, 30)) -> int:
    """Distinct peptides (as written; I and L distinct) with 0 missed cleavages and length in the
    inclusive range ``length``."""
    lo, hi = length
    return len({p for p in digest(aa, enzyme) if lo <= len(p) <= hi})


def tm_segments(aa: str, window: int = 19, threshold: float = 1.6, cap: int = 10) -> int:
    """Kyte-Doolittle TM proxy (G6): maximal runs of consecutive windows with mean >= threshold.

    Letters outside the 20 standard amino acids score 0; sequences shorter than the window give 0.
    """
    if window < 1:
        raise ValueError(f"window must be >= 1, got {window}")
    n = len(aa)
    if n < window:
        return 0
    codes = np.frombuffer(aa.encode("latin-1", errors="replace"), dtype=np.uint8)
    csum = np.concatenate(([0], np.cumsum(_KD10_LUT[codes])))
    sums = csum[window:] - csum[:-window]
    # mean >= threshold  <=>  sum(10 * KD) >= 10 * threshold * window (exact integer sums).
    hit = sums >= round(10.0 * threshold * window, 9)
    starts = int(hit[0]) + int(np.count_nonzero(hit[1:] & ~hit[:-1]))
    return min(cap, starts)


def _translation_maps(table: int) -> tuple[dict[str, str], frozenset[str]]:
    ct = CodonTable.unambiguous_dna_by_id[table]
    return dict(ct.forward_table), frozenset(ct.stop_codons)


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def protein_sequence(aa: str | None, nt: str, feature_type: str, table: int = 11) -> str:
    """Sequence the covariates are computed on (spec §7.3.3 pseudogene rule).

    ``aa`` when present (even empty); otherwise the translation of ``nt`` with ``table`` from the
    first codon, whole codons only, first codon read as M, stopping before the first stop codon.
    ``feature_type`` is accepted for the contract; the rule applies whenever ``aa`` is missing.
    """
    del feature_type
    if not _is_missing(aa):
        return str(aa)
    forward, stops = _translation_maps(table)
    seq = str(nt).upper()
    out: list[str] = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        if i == 0:
            out.append("M")
            continue
        if codon in stops:
            break
        out.append(forward.get(codon, "X"))
    return "".join(out)


def _protparam(clean: str) -> tuple[float, float, float]:
    """(GRAVY, pI, average MW in Da) from Biopython ProtParam on a non-empty standard sequence."""
    pa = ProteinAnalysis(clean)  # type: ignore[no-untyped-call]
    return (
        float(pa.gravy()),  # type: ignore[no-untyped-call]
        float(pa.isoelectric_point()),  # type: ignore[no-untyped-call]
        float(pa.molecular_weight()),  # type: ignore[no-untyped-call]
    )


def covariates(aa: str, enzyme: Enzyme, length: tuple[int, int]) -> Covariates:
    """Covariates of one protein; ``aa`` upper case, leading/trailing ``*`` stripped."""
    _check_enzyme(enzyme)
    if _VALID.fullmatch(aa) is None:
        raise ValueError("protein sequence must be upper-case letters (and '*')")
    seq = aa.strip("*")
    clean = _NON_STANDARD.sub("", seq)
    if not clean:
        nan = math.nan
        return Covariates(0, 0, nan, nan, nan, 0)
    gravy, pi, mw = _protparam(clean)
    return Covariates(
        length_aa=len(seq),
        n_observable_peptides=n_observable_peptides(seq, enzyme, length),
        gravy=gravy,
        pi=pi,
        mw_da=mw,
        tm_segments=tm_segments(seq),
    )


def covariate_table(genes: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """``locus_tag`` + the six covariate columns, one row per C2 gene, input row order kept."""
    enzyme = cfg.abundance.enzyme
    lo, hi = cfg.abundance.observable_len
    table = cfg.annotation.translation_table
    rows = [
        astuple(covariates(protein_sequence(aa, nt, ft, table), enzyme, (lo, hi)))
        for aa, nt, ft in zip(
            genes["aa"].tolist(),
            genes["nt"].tolist(),
            genes["feature_type"].tolist(),
            strict=True,
        )
    ]
    names = [f.name for f in fields(Covariates)]
    out = pd.DataFrame.from_records(rows, columns=names) if rows else pd.DataFrame(columns=names)
    out = out.astype(
        {
            "length_aa": "int64",
            "n_observable_peptides": "int64",
            "gravy": "float64",
            "pi": "float64",
            "mw_da": "float64",
            "tm_segments": "int64",
        }
    )
    out.insert(0, "locus_tag", genes["locus_tag"].to_numpy())
    return out
