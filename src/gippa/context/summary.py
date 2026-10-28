"""Detection summary: basic O - E rows, QC rows and the ``detection_summary`` assembler.

GIPPA_SPEC.md §7.3.5(a), §8.2 ("resample TUs, not genes"), §9.3 (``detection_summary.tsv``);
docs/specs/M4-16-detection-summary.md; docs/methods/detection-summary.md. Pure functions: the
caller passes the generator from ``gippa.rng.for_step("context_stats", seed)``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from gippa.stats.bootstrap import tu_block_bootstrap

__all__ = [
    "FEW_UNITS_FLAG",
    "REQUIRES_VARIANTS_FLAG",
    "SECTIONS",
    "STRATUM0_METRIC",
    "SUMMARY_COLUMNS",
    "basic_rows",
    "detection_summary",
    "placement_kind_rows",
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
REQUIRES_VARIANTS_FLAG = "requires_variants"  # id_mixups (§7.5.1) needs the variant module (M7)
STRATUM0_METRIC = "stratum0_unmapped_share"  # ADR-0019 decision 4 (§7.5.1 genome-inconsistent)
NO_COG = frozenset({"", "-"})

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
    if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA:
        return ["none"]
    s = str(v).strip()
    if s in NO_COG:
        return ["none"]
    return sorted(set(s))


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


def _is_missing(v: object) -> bool:
    if v is None or v is pd.NA:
        return True
    if isinstance(v, float) and math.isnan(v):
        return True
    return isinstance(v, str) and v.strip() == ""


def _detected(col: pd.Series) -> np.ndarray:
    """``detected`` as a bool array; strings or other non-boolean values raise ValueError."""
    if pd.api.types.is_bool_dtype(col.dtype):
        if col.isna().any():
            raise ValueError("detected has missing values on eligible genes")
        return np.asarray(col.to_numpy(dtype=bool), dtype=bool)
    if col.dtype == object and all(isinstance(v, bool | np.bool_) for v in col):
        return np.asarray(col.to_numpy(dtype=bool), dtype=bool)
    raise ValueError(f"detected must be boolean, got dtype {col.dtype}")


def _check_eligible(elig: pd.DataFrame) -> None:
    """Fail closed: every eligible gene needs a finite ``p_detect`` and a ``tu_id``."""
    tags = elig["locus_tag"].astype(str).to_numpy()
    p = pd.to_numeric(elig["p_detect"], errors="coerce").astype(float).to_numpy()
    bad_p = sorted(tags[~np.isfinite(p)].tolist())
    no_tu = np.array([_is_missing(v) for v in elig["tu_id"]], dtype=bool)
    bad_tu = sorted(tags[no_tu].tolist()) if len(tags) else []
    msgs = []
    if bad_p:
        msgs.append(f"eligible genes with NaN p_detect: {bad_p}")
    if bad_tu:
        msgs.append(f"eligible genes without tu_id: {bad_tu}")
    if msgs:
        raise ValueError("; ".join(msgs))


def _few_units(n_tu: int, expected: float, n_genes: int) -> bool:
    """CI flag: < 10 TUs, or min(E, n - E) < 10 (few expected detections or non-detections)."""
    return n_tu < FEW_UNITS_MIN or min(expected, n_genes - expected) < FEW_UNITS_MIN


def basic_rows(genes: pd.DataFrame, B: int, rng: np.random.Generator) -> pd.DataFrame:
    """Basic O - E rows (section ``basic``) per replicon, feature class and COG category.

    Eligible genes only (``detection_status != partial_excluded``). Metrics per group:
    ``n_genes``, ``n_detected``, ``pct_detected``, ``expected`` (sum of ``p_detect``),
    ``o_minus_e`` (95 % TU-block bootstrap CI over the group's TUs, flag ``ci_few_units`` with
    fewer than 10 TUs or min(E, n_genes - E) < 10), ``o_over_e`` (NA when E = 0) and ``n_tu``.
    ``n`` = n_genes on every row. Groups are processed in sorted order so the draws are
    reproducible. COG ``-`` and blank count in ``none``.

    Raises:
        ValueError: ``detected`` is not boolean, or an eligible gene has NaN ``p_detect`` or no
            ``tu_id`` (the genes are named).
    """
    elig = genes[genes["detection_status"].astype(str) != "partial_excluded"]
    elig = elig.reset_index(drop=True)
    det = _detected(elig["detected"])
    _check_eligible(elig)
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
        flag = FEW_UNITS_FLAG if _few_units(n_tu, e, n) else ""
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


def _stratum0_rows(stratum0: pd.DataFrame) -> list[Row]:
    """ADR-0019 decision 4 QC rows: the stratum-0 unmapped share k/n, overall and per class.

    ``stratum0`` is ``gippa.mapping.placement_fdr.stratum0_pi_classes``'s table (one row per
    length class with ``n``, ``k`` and ``upper``). The overall row is group ``run``/``all``, the
    per-class rows group_type ``length_class`` with the class label; ``value`` = k/n (NA for an
    empty class, never 0) and ``n`` = n. A per-class row reports the table's ``upper`` as
    ``ci_high``: the one-sided 95 % Clopper-Pearson bound of that same k/n (``ci_low`` stays NA,
    the interval is one-sided). The overall row has no such bound in the table and leaves both NA.
    """
    nan = math.nan
    n_all, k_all = int(stratum0["n"].sum()), int(stratum0["k"].sum())
    share = k_all / n_all if n_all else nan
    rows: list[Row] = [("qc", "run", "all", STRATUM0_METRIC, share, nan, nan, n_all, "")]
    for r in stratum0.itertuples(index=False):
        n, k = int(r.n), int(r.k)
        upper = float(r.upper) if n else nan
        rows.append(("qc", "length_class", str(r.length_class), STRATUM0_METRIC,
                     k / n if n else nan, nan, upper, n, ""))  # fmt: skip
    return rows


def qc_rows(
    peptides: pd.DataFrame | None,
    peptide_hits: pd.DataFrame,
    unmapped: pd.DataFrame | None,
    stratum0: pd.DataFrame | None = None,
    id_mixups: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """QC rows (section ``qc``, group ``run``/``all``) from C9, C11 and C13 (§7.5.1, §9.3).

    ``id_mixups``: the variant module's ``(count, n)`` (M7-14 rule 3: same-db ALT sites flagged
    neither ``start_codon`` nor ``stop_codon``, of ``n`` same-db ALT sites); without it (no
    ``variants_qc``, e.g. no peptides) NA with flag ``requires_variants`` (M4-16).
    ``engine_protein_mismatch_rate``: distinct C11 peptides flagged ``engine_protein_mismatch``
    over distinct placed (C11) peptides; NA without C11 rows.
    ``genome_inconsistent_rate``: share of distinct C9 target, non-contaminant (``is_contaminant``)
    peptides with ``search_db = other`` without an exact C11 placement; ``n`` = that denominator
    (NA without such peptides, e.g. all ``same``). C13 counts by ``label`` and ``kind`` (NA
    without C13).

    With the ADR-0019 class table ``stratum0`` (the ``placement`` step's ``stratum0_classes``)
    the rows of ``_stratum0_rows`` follow. They measure another population than
    ``genome_inconsistent_rate`` and replace no row of it (M3-25 decision 2).
    """
    hits = peptide_hits
    pep_col = hits["peptide"].astype(str)
    mm_rate: float = math.nan
    n_placed: Any = pd.NA
    if len(hits):
        if "flags" in hits:
            mixed = np.array([_has_flag(f, MISMATCH_FLAG) for f in hits["flags"]], dtype=bool)
        else:
            mixed = np.zeros(len(hits), dtype=bool)
        n_placed = int(pep_col.nunique())
        mm_rate = float(pep_col[mixed].nunique()) / n_placed

    rate: float = math.nan
    n_other: Any = pd.NA
    if peptides is not None and "search_db" in peptides:
        keep = ~peptides["is_decoy"].astype(bool)
        if "is_contaminant" in peptides:
            keep &= ~peptides["is_contaminant"].fillna(False).astype(bool)
        other = peptides[keep & (peptides["search_db"].astype(str) == "other")]
        if len(other):
            exact = set(pep_col[hits["route"].astype(str) == "exact"]) if len(hits) else set()
            uniq = other["peptide"].astype(str).unique()
            rate = float(np.mean([pe not in exact for pe in uniq]))
            n_other = len(uniq)

    def count(label: str, kind: str) -> float:
        if unmapped is None:
            return math.nan
        sel = (unmapped["label"].astype(str) == label) & (unmapped["kind"].astype(str) == kind)
        return float(sel.sum())

    nan = math.nan
    vals: list[tuple[str, float, Any, str]] = [
        ("id_mixups", nan, pd.NA, REQUIRES_VARIANTS_FLAG)
        if id_mixups is None
        else ("id_mixups", float(id_mixups[0]), int(id_mixups[1]), ""),
        ("engine_protein_mismatch_rate", mm_rate, n_placed, ""),
        ("genome_inconsistent_rate", rate, n_other, ""),
        ("contaminant_peptides", count("contaminant", "peptide"), pd.NA, ""),
        ("contaminant_proteins", count("contaminant", "protein"), pd.NA, ""),
        ("unmapped_peptides", count("unmapped", "peptide"), pd.NA, ""),
        ("unmapped_proteins", count("unmapped", "protein"), pd.NA, ""),
    ]
    rows: list[Row] = [("qc", "run", "all", m, v, nan, nan, n, f) for m, v, n, f in vals]
    if stratum0 is not None:
        rows += _stratum0_rows(stratum0)
    return _frame(rows)


KIND_GROUPS = ("exact", "tolerant", "frameshift")
DIFF_GROUPS = ("substitution", "swap", "block", "near_block")
FLAG_GROUPS = ("readthrough", "pyl_site", "sec_site")  # rule 17; 16c adds pyl_site, sec_site
REASON_GROUPS = ("no_hit", "rejected_placement_fdr", "too_short", "contaminant_id")


def placement_kind_rows(peptide_hits: pd.DataFrame, unmapped: pd.DataFrame | None) -> pd.DataFrame:
    """M7-16 rule 17: distinct peptides by placement kind, tolerant diff type, placement flag and
    unmapped reason (section ``qc``, metric ``n_peptides``). Counts with their population ``n``,
    never rates; every group is written. Without C13 the ``unmapped_reason`` rows are NA."""
    kind: dict[str, str] = {}
    diffs: dict[str, set[str]] = {}
    flagged: dict[str, set[str]] = {}
    h = peptide_hits
    for pep, route, dt, fl in zip(
        h["peptide"].astype(str), h["route"].astype(str), h["diff_types"], h["flags"], strict=True
    ):
        fs = _has_flag(dt, "frameshift")
        k = "exact" if route == "exact" else "frameshift" if fs else "tolerant"
        if pep not in kind or KIND_GROUPS.index(k) < KIND_GROUPS.index(kind[pep]):
            kind[pep] = k
        if k == "tolerant":
            diffs.setdefault(pep, set()).update(g for g in DIFF_GROUPS if _has_flag(dt, g))
        flagged.setdefault(pep, set()).update(g for g in FLAG_GROUPS if _has_flag(fl, g))
    tol = [p for p, k in kind.items() if k == "tolerant"]
    nan = math.nan
    rows: list[Row] = []

    def add(gt: str, g: str, v: float, n: Any) -> None:
        rows.append(("qc", gt, g, "n_peptides", v, nan, nan, n, ""))

    for g in KIND_GROUPS:
        add("placement_kind", g, float(sum(k == g for k in kind.values())), len(kind))
    for g in DIFF_GROUPS:
        add("tolerant_diff_type", g, float(sum(g in diffs.get(p, ()) for p in tol)), len(tol))
    for g in FLAG_GROUPS:
        add("placement_flag", g, float(sum(g in f for f in flagged.values())), len(kind))
    pep_rows = None if unmapped is None else unmapped[unmapped["kind"].astype(str) == "peptide"]
    for g in REASON_GROUPS:
        if pep_rows is None:
            add("unmapped_reason", g, nan, pd.NA)
        else:
            add("unmapped_reason", g, float((pep_rows["reason"].astype(str) == g).sum()),
                len(pep_rows))  # fmt: skip
    return _frame(rows)


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
