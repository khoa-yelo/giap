"""Detection summary: basic O - E rows, QC rows and the ``detection_summary`` assembler.

GIAP_SPEC.md §7.3.5(a), §8.2 ("resample TUs, not genes"), §9.3 (``detection_summary.tsv``);
docs/specs/M4-16-detection-summary.md; docs/methods/detection-summary.md. Pure functions: the
caller passes the generator from ``giap.rng.for_step("context_stats", seed)``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from giap.stats.bootstrap import tu_block_bootstrap

__all__ = [
    "FEW_UNITS_FLAG",
    "SECTIONS",
    "SUMMARY_COLUMNS",
    "basic_rows",
    "detection_summary",
    "qc_rows",
]

SUMMARY_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n",
                   "flags")  # fmt: skip
CORE_COLUMNS = SUMMARY_COLUMNS[:7]
SECTIONS = ("tier", "basic", "detectability", "qc")
FEW_UNITS_FLAG = "ci_few_units"  # docs/methods/bootstrap.md: CI from < 10 resampling units
FEW_UNITS_MIN = 10
FEATURE_CLASSES = ("cds", "sorf", "pseudogene")
MISMATCH_FLAG = "engine_protein_mismatch"

Row = tuple[str, str, str, str, float, float, float, Any, str]


def _frame(rows: Sequence[Row]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(SUMMARY_COLUMNS))
    for c in ("section", "group_type", "group", "metric", "flags"):
        df[c] = df[c].astype(object)
    for c in ("value", "ci_low", "ci_high"):
        df[c] = df[c].astype("float64")
    df["n"] = df["n"].astype("Int64")
    return df


def _cog_letters(v: object) -> list[str]:
    if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA or str(v) == "":
        return ["none"]
    return sorted(set(str(v)))


def _groups(genes: pd.DataFrame) -> list[tuple[str, str, np.ndarray]]:
    """Sorted (group_type, group, mask) over eligible genes; empty groups are dropped."""
    contig = genes["contig"].astype(str).to_numpy()
    ftype = genes["feature_type"].astype(str).to_numpy()
    letters = [_cog_letters(v) for v in genes["cog_category"]]
    out: list[tuple[str, str, np.ndarray]] = []
    out += [("replicon", c, contig == c) for c in sorted(set(contig))]
    out += [("feature_class", f, ftype == f) for f in FEATURE_CLASSES]
    out.append(("feature_class", "hypothetical", genes["hypothetical"].astype(bool).to_numpy()))
    for letter in sorted({x for ls in letters for x in ls}):
        out.append(("cog_category", letter, np.array([letter in ls for ls in letters], dtype=bool)))
    return sorted((g for g in out if g[2].any()), key=lambda g: (g[0], g[1]))


def _sum_v(d: pd.DataFrame) -> float:
    return float(d["v"].sum())


def basic_rows(genes: pd.DataFrame, B: int, rng: np.random.Generator) -> pd.DataFrame:
    """Basic O - E rows (section ``basic``) per replicon, feature class and COG category.

    Eligible genes only (``detection_status != partial_excluded``). Metrics per group:
    ``n_genes``, ``n_detected``, ``pct_detected``, ``expected`` (sum of ``p_detect``),
    ``o_minus_e`` (95 % TU-block bootstrap CI over the group's TUs, flag ``ci_few_units`` with
    fewer than 10 TUs), ``o_over_e`` (NA when E = 0) and ``n_tu``. ``n`` = n_genes on every row.
    Groups are processed in sorted order so the draws are reproducible.
    """
    elig = genes[genes["detection_status"].astype(str) != "partial_excluded"]
    elig = elig.reset_index(drop=True)
    det = elig["detected"].astype(bool).to_numpy()
    p = elig["p_detect"].astype(float).to_numpy()
    tu = elig["tu_id"].astype(str).to_numpy()
    rows: list[Row] = []
    for gtype, group, mask in _groups(elig):
        n = int(mask.sum())
        o, e = float(det[mask].sum()), float(p[mask].sum())
        resid = pd.Series(det[mask].astype(float) - p[mask])
        sums = resid.groupby(tu[mask], sort=True).sum()
        tus = pd.DataFrame({"tu": sums.index.astype(str), "v": sums.to_numpy()})
        ci = tu_block_bootstrap(tus, "tu", _sum_v, B, rng)
        n_tu = len(tus)
        flag = FEW_UNITS_FLAG if n_tu < FEW_UNITS_MIN else ""
        nan = math.nan
        vals: list[tuple[str, float, float, float, str]] = [
            ("n_genes", float(n), nan, nan, ""),
            ("n_detected", o, nan, nan, ""),
            ("pct_detected", 100.0 * o / n, nan, nan, ""),
            ("expected", e, nan, nan, ""),
            ("o_minus_e", o - e, ci.ci_low, ci.ci_high, flag),
            ("o_over_e", o / e if e > 0 else nan, nan, nan, ""),
            ("n_tu", float(n_tu), nan, nan, ""),
        ]
        rows += [("basic", gtype, group, m, v, lo, hi, n, f) for m, v, lo, hi, f in vals]
    return _frame(rows)


def _has_flag(v: object, flag: str) -> bool:
    if isinstance(v, list | tuple | set | np.ndarray):
        return flag in {str(x) for x in v}
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return False
    return flag in str(v).split(";")


def qc_rows(
    peptides: pd.DataFrame | None, peptide_hits: pd.DataFrame, unmapped: pd.DataFrame | None
) -> pd.DataFrame:
    """QC rows (section ``qc``, group ``run``/``all``) from C9, C11 and C13.

    ``id_mixups``: distinct C11 peptides flagged ``engine_protein_mismatch``.
    ``genome_inconsistent_rate``: share of C9 target peptides with ``search_db = same`` without an
    exact C11 placement (NA without such peptides). C13 counts by ``label`` and ``kind`` (NA
    without C13).
    """
    hits = peptide_hits
    pep_col = hits["peptide"].astype(str)
    mixed = [_has_flag(f, MISMATCH_FLAG) for f in hits["flags"]] if "flags" in hits else []
    id_mixups = float(pep_col[np.asarray(mixed, dtype=bool)].nunique()) if len(hits) else 0.0

    rate, n_same = math.nan, pd.NA
    if peptides is not None and "search_db" in peptides:
        target = ~peptides["is_decoy"].astype(bool)
        same = peptides[target & (peptides["search_db"].astype(str) == "same")]
        if len(same):
            exact = set(pep_col[hits["route"].astype(str) == "exact"]) if len(hits) else set()
            uniq = same["peptide"].astype(str).unique()
            rate = float(np.mean([pe not in exact for pe in uniq]))
            n_same = len(uniq)

    def count(label: str, kind: str) -> float:
        if unmapped is None:
            return math.nan
        sel = (unmapped["label"].astype(str) == label) & (unmapped["kind"].astype(str) == kind)
        return float(sel.sum())

    nan = math.nan
    vals: list[tuple[str, float, Any]] = [
        ("id_mixups", id_mixups, pd.NA),
        ("genome_inconsistent_rate", rate, n_same),
        ("contaminant_peptides", count("contaminant", "peptide"), pd.NA),
        ("contaminant_proteins", count("contaminant", "protein"), pd.NA),
        ("unmapped_peptides", count("unmapped", "peptide"), pd.NA),
        ("unmapped_proteins", count("unmapped", "protein"), pd.NA),
    ]
    return _frame([("qc", "run", "all", m, v, nan, nan, n, "") for m, v, n in vals])


def _flags_str(v: object) -> str:
    if isinstance(v, list | tuple | set | np.ndarray):
        return ";".join(str(x) for x in v if str(x))
    if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)):
        return ""
    return str(v)


def detection_summary(parts: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate summary parts (tier, basic, detectability, qc), validate and sort.

    Parts need the seven core columns; a missing ``n`` becomes NA and missing ``flags`` empty.
    Rows sort by section order, then ``group_type``, ``group``, ``metric`` as strings.

    Raises:
        ValueError: on a part without a core column or a row with an unknown section.
    """
    frames = []
    for i, part in enumerate(parts):
        missing = [c for c in CORE_COLUMNS if c not in part.columns]
        if missing:
            raise ValueError(f"summary part {i} lacks columns {missing}")
        df = part.copy()
        if "n" not in df:
            df["n"] = pd.NA
        if "flags" not in df:
            df["flags"] = ""
        frames.append(df[list(SUMMARY_COLUMNS)])
    rows: list[Row] = []
    for df in frames:
        for r in df.itertuples(index=False):
            n = None if pd.isna(r.n) else int(r.n)
            rows.append((str(r.section), str(r.group_type), str(r.group), str(r.metric),
                         float(r.value), float(r.ci_low), float(r.ci_high), n,
                         _flags_str(r.flags)))  # fmt: skip
    bad = sorted({r[0] for r in rows} - set(SECTIONS))
    if bad:
        raise ValueError(f"unknown detection_summary sections {bad}")
    order = {s: i for i, s in enumerate(SECTIONS)}
    rows.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))
    return _frame(rows)
