"""Independent oracle and builders for M4-14 gene-level abundance (docs/specs/M4-14-abundance.md).

Written from the spec text only (GIPPA_SPEC.md §7.3.4, §7.2.6, §9.3), with plain Python loops so it
shares no code with ``gippa.genetable.abundance``:

- peptide route: accepted rows with ``abundance_eligible`` and a finite intensity > 0;
  iBAQ = sum / max(1, n_observable); Top3 = arithmetic mean of the (up to) 3 largest linear values;
- protein route: C10 joined through the protein -> gene map, values on the log2 scale (linear values
  -> log2, values <= 0 -> missing), the larger value per gene;
- imputation: fill statistic over observed log2 values of detected eligible (non-partial) genes
  (median, mean, min, zero = min - 1, exclude = none), filled for every detected gene with a
  missing value, partial included (PR #82 lead decision 5).

Also a small peptide-level simulator on top of gippasim layouts/expression (M4-07/M4-08), since the
M4-09 evidence simulator is not merged: each detected gene's observable tryptic peptides (0 missed
cleavages, length 7-30) are quantified independently with probability ``p_obs`` and intensity
``ln I = abundance_log2_true * ln 2 + N(0, sd^2)`` (the M4-09 intensity model).
"""

from __future__ import annotations

import math
import re
import statistics
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

ABUNDANCE_COLUMNS = (
    "abundance_raw",
    "abundance_log2",
    "abundance_imputed",
    "imputation_policy",
    "abundance_type",
    "abundance_scale",
    "top3_log2",
)
POLICIES = ("median", "mean", "min", "zero", "exclude")
ACCEPTED_COLUMNS = (
    "gene", "peptide", "q_value", "pep", "score", "intensity", "uniqueness", "route",
    "locus_class", "abundance_eligible", "start", "end", "contig", "crosses_origin",
    "is_contaminant",
)  # fmt: skip
COVARIATE_COLUMNS = ("length_aa", "n_observable_peptides", "gravy", "pi", "mw_da", "tm_segments")


def is_na(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, (list, tuple, np.ndarray)):
        return False
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


# ------------------------------------------------------------------------------------ builders


def acc(
    gene: str,
    peptide: str,
    intensity: float | None,
    *,
    uniqueness: str = "unique",
    route: str = "exact",
    locus_class: str = "cds_inframe",
    eligible: bool | None = None,
    start: int = 0,
) -> dict[str, Any]:
    """One accepted (gene, peptide) row in the M4-11 ``accepted_peptides`` shape."""
    if eligible is None:
        eligible = uniqueness == "unique" and locus_class not in (
            "nterm_extension",
            "cterm_extension",
        )
    return {
        "gene": gene,
        "peptide": peptide,
        "q_value": 0.001,
        "pep": 0.01,
        "score": 50.0,
        "intensity": math.nan if intensity is None else float(intensity),
        "uniqueness": uniqueness,
        "route": route,
        "locus_class": locus_class,
        "abundance_eligible": bool(eligible),
        "start": int(start),
        "end": int(start) + 3 * len(peptide),
        "contig": "chr1",
        "crosses_origin": False,
        "is_contaminant": False,
        "placements": (
            ("chr1", int(start), int(start) + 3 * len(peptide), route, uniqueness, locus_class),
        ),
    }


def accepted_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    rows = list(rows)
    cols = (*ACCEPTED_COLUMNS, "placements")
    if not rows:
        df = pd.DataFrame({c: pd.Series([], dtype=object) for c in cols})
    else:
        df = pd.DataFrame(rows, columns=list(cols))
    return df.astype(
        {
            "gene": object,
            "peptide": object,
            "q_value": "float64",
            "pep": "float64",
            "score": "float64",
            "intensity": "float64",
            "abundance_eligible": bool,
            "start": "int64",
            "end": "int64",
            "crosses_origin": bool,
            "is_contaminant": bool,
        }
    )


def covariates_frame(n_obs: Mapping[str, int], order: Sequence[str] | None = None) -> pd.DataFrame:
    """M4-10 ``covariate_table`` shape; only ``n_observable_peptides`` matters here."""
    tags = list(order) if order is not None else list(n_obs)
    return pd.DataFrame(
        {
            "locus_tag": tags,
            "length_aa": np.full(len(tags), 300, dtype=np.int64),
            "n_observable_peptides": np.array([int(n_obs.get(t, 0)) for t in tags], np.int64),
            "gravy": np.zeros(len(tags)),
            "pi": np.full(len(tags), 7.0),
            "mw_da": np.full(len(tags), 33000.0),
            "tm_segments": np.zeros(len(tags), dtype=np.int64),
        }
    )


@dataclass(frozen=True)
class DetectionStub:
    """Stand-in for M4-12 ``Detection`` (only ``table`` is used by M4-14)."""

    tier: str
    table: pd.DataFrame
    summary_rows: pd.DataFrame
    warnings: tuple[str, ...] = ()


def detection(
    tags: Sequence[str], detected: Iterable[str], partial: Iterable[str] = ()
) -> DetectionStub:
    det, part = set(detected), set(partial)
    status = [
        "partial_excluded" if t in part else ("detected" if t in det else "not_detected")
        for t in tags
    ]
    table = pd.DataFrame(
        {
            "locus_tag": list(tags),
            "detected": [t in det for t in tags],
            "detection_status": status,
            "gene_q_value": [0.001 if t in det else math.nan for t in tags],
            "detection_tier": ["decoy_picked"] * len(tags),
        }
    )
    return DetectionStub("decoy_picked", table, pd.DataFrame())


def c10(values: Mapping[str, float | None]) -> pd.DataFrame:
    """C10 ``abundance.norm.parquet`` rows: ``protein_id``, ``value`` and NO ``scale`` column (M4-06
    decision on PR #51; the measurement scale is ``abundance.scale`` or preflight, passed to
    ``gene_abundance`` as ``c10_scale``; PR #82 lead decision 1)."""
    return pd.DataFrame(
        {
            "protein_id": list(values),
            "value": [math.nan if v is None else float(v) for v in values.values()],
        }
    )


def pmap_frame(rows: Iterable[tuple[str, str | None, str]]) -> pd.DataFrame:
    """M4-12 protein -> gene map: only ``protein_id``, ``gene``, ``method`` (PR #82 decision 2)."""
    rows = list(rows)
    return pd.DataFrame(
        {
            "protein_id": pd.Series([r[0] for r in rows], dtype=object),
            "gene": pd.Series([r[1] for r in rows], dtype=object),
            "method": pd.Series([r[2] for r in rows], dtype=object),
        }
    )


PROTEIN_HITS_COLUMNS = (
    "protein_id", "protein_len_aa", "route", "contig", "start", "end", "strand", "rank",
    "identity", "query_coverage", "reciprocal_overlap", "class", "gene", "flags",
)  # fmt: skip


def protein_hits_frame(rows: Iterable[tuple[Any, ...]]) -> pd.DataFrame:
    """C12 ``protein_hits`` (M3-16) in the M4Inputs shape: one (protein_id, gene, flags[, rank])
    per hit, ``flags`` a list (e.g. ``cross_strain_protein``); other columns are placeholders."""
    recs = []
    for i, r in enumerate(rows):
        rank = int(r[3]) if len(r) > 3 else 1
        recs.append(
            {
                "protein_id": r[0], "protein_len_aa": 300, "route": "miniprot", "contig": "chr1",
                "start": 1000 * i, "end": 1000 * i + 900, "strand": "+", "rank": rank,
                "identity": 0.95, "query_coverage": 0.95, "reciprocal_overlap": 0.95,
                "class": "near", "gene": r[1], "flags": list(r[2]),
            }
        )  # fmt: skip
    if not recs:
        return pd.DataFrame({c: pd.Series([], dtype=object) for c in PROTEIN_HITS_COLUMNS})
    return pd.DataFrame(recs, columns=list(PROTEIN_HITS_COLUMNS))


# -------------------------------------------------------------------------------------- oracle


def peptide_values(accepted: pd.DataFrame) -> dict[str, list[float]]:
    """Quantified linear intensities per gene (eligible, finite, > 0)."""
    out: dict[str, list[float]] = {}
    for row in accepted.to_dict("records"):
        if not bool(row["abundance_eligible"]):
            continue
        x = row["intensity"]
        if is_na(x) or not math.isfinite(float(x)) or float(x) <= 0:
            continue
        out.setdefault(str(row["gene"]), []).append(float(x))
    return out


def ibaq(values: Sequence[float], n_observable: int) -> float:
    return math.fsum(values) / max(1, int(n_observable))


def top3(values: Sequence[float]) -> float:
    top = sorted(values, reverse=True)[:3]
    return math.fsum(top) / len(top)


def fill_value(observed: Sequence[float], policy: str) -> float | None:
    if policy == "exclude" or not observed:
        return None
    if policy == "median":
        return float(statistics.median(observed))
    if policy == "mean":
        return math.fsum(observed) / len(observed)
    if policy == "min":
        return float(min(observed))
    if policy == "zero":
        return float(min(observed)) - 1.0
    raise ValueError(policy)


def impute(
    log2: Mapping[str, float],
    detected: Mapping[str, bool],
    policy: str,
    eligible: Mapping[str, bool] | None = None,
) -> tuple[dict[str, float], dict[str, bool], float | None]:
    """Fill statistic over observed values of detected AND eligible (non-partial; default all)
    genes; filled: every detected gene with a missing value, partial or not (PR #82 decision 5)."""
    observed = [
        v
        for t, v in log2.items()
        if detected.get(t, False) and (eligible is None or eligible.get(t, True)) and not is_na(v)
    ]
    fill = fill_value(observed, policy)
    filled: dict[str, float] = {}
    mask: dict[str, bool] = {}
    for t, v in log2.items():
        do = fill is not None and detected.get(t, False) and is_na(v)
        filled[t] = fill if do else v  # type: ignore[assignment]
        mask[t] = bool(do)
    return filled, mask, fill


def to_log2(value: float | None, scale: str) -> float:
    if is_na(value):
        return math.nan
    v = float(value)  # type: ignore[arg-type]
    if not math.isfinite(v):
        return math.nan  # +-inf never survives (PR #82 decision 3)
    if scale == "log2":
        return v
    return math.log2(v) if v > 0 else math.nan


# ----------------------------------------------------------------------- peptide-level simulator

_CLEAVE = re.compile(r"(?<=[KR])(?!P)")


def observable_peptides(aa: str, lo: int = 7, hi: int = 30) -> list[str]:
    """Trypsin (no cleavage before P), 0 missed cleavages, length lo..hi (oracle copy)."""
    seq = aa[:-1] if aa.endswith("*") else aa
    return [p for p in _CLEAVE.split(seq) if lo <= len(p) <= hi]


@dataclass(frozen=True)
class SimPeptides:
    accepted: pd.DataFrame
    n_observable: dict[str, int]
    truth: dict[str, float]  # abundance_log2_true per gene
    detected: list[str]


def simulate_peptides(
    layout: Any,
    expression: Any,
    rng: np.random.Generator,
    *,
    p_obs: float = 0.6,
    sd_ln: float = 0.5,
) -> SimPeptides:
    """Quantified unique exact peptides of every ``detected_true`` gene (see module docstring)."""
    genes = layout.genes
    ex = expression.genes.set_index("locus_tag")
    rows: list[dict[str, Any]] = []
    n_obs: dict[str, int] = {}
    truth: dict[str, float] = {}
    detected: list[str] = []
    for tag, aa in zip(genes["locus_tag"].astype(str), genes["aa"], strict=True):
        peps = observable_peptides(aa) if isinstance(aa, str) and aa else []
        n_obs[tag] = len(peps)
        a = float(ex.at[tag, "abundance_log2_true"])
        truth[tag] = a
        if not bool(ex.at[tag, "detected_true"]):
            continue
        detected.append(tag)
        seen = rng.random(len(peps)) < p_obs
        noise = rng.normal(0.0, sd_ln, len(peps))
        used: set[str] = set()
        for j, pep in enumerate(peps):
            if not seen[j] or pep in used:
                continue
            used.add(pep)
            rows.append(acc(tag, pep, math.exp(a * math.log(2.0) + noise[j]), start=3 * j))
    return SimPeptides(accepted_frame(rows), n_obs, truth, detected)
