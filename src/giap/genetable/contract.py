"""M4 input contract: the tables M4 reads from M1, M2 and M3 (docs/specs/M4-06-input-contract.md).

Contract tables (key: file, producer; invariants checked here):

- C1 ``replicons``: replicons.parquet (M2-04); ``contig`` unique.
- C2 ``genes``: genes.parquet (M2-04); ``locus_tag`` unique; ``0 <= start < end``;
  ``end <= length`` unless ``crosses_origin`` (then circular and ``length < end < 2 length``);
  ``aa`` NA only for pseudogenes and without a trailing ``*``.
- C3 ``xrefs``: xrefs.parquet (M2-04); long format, columns only.
- C5 ``tus``: tus.parquet (M2-12); every C2 gene in exactly one TU.
- C6 ``layers``: layers.parquet (M2-20); one row per C2 gene; ``tu_id`` agrees with C5.
- C7 ``regions``: regions.parquet (M2-19); AMR rows have NA coordinates.
- C8 ``region_genes``: region_genes.parquet (M2-19); ``subset`` non-NA iff prophage member.
- C9 ``peptides``: peptides.norm.parquet (M1-12); one row per peptide; ``protein_ids`` a native
  list of non-empty strings (lead decision 1; ``;``-joined strings are dtype errors).
- C10 ``abundance``: abundance.norm.parquet (M1-12, optional); ``raw_id`` optional, no scale
  column (the loader reads ``routing.abundance_scale`` from preflight.json).
- C11 ``peptide_hits``: mapping/peptide_hits.tsv (M3-16); ``gene`` set for assigning classes,
  NA for non-assigning ones (pseudogene: either), always a C2 tag; no C9 decoy placed.
- C12 ``protein_hits``: mapping/protein_hits.tsv (M3-16); ``gene`` a C2 tag when set.
- C13 ``unmapped``: mapping/unmapped.tsv (M3-16); C14 ``id_map``: annotation/id_map.tsv (opt).

C4 (``annotation_meta.json``) is read by the loader (``giap.io.work_tables``). Frames hold the
internal representation: 0-based half-open coordinates and list cells (``list``, ``tuple`` or numpy
arrays) for list columns. Every check is pure and returns M1-04 ``Issue`` objects with codes
``contract.<table>.<rule>``; nothing here raises on content.
"""

from __future__ import annotations

import numbers
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np
import pandas as pd

from giap.io.issues import MAX_EXAMPLES, Issue, Severity

__all__ = [
    "ASSIGNING_CLASSES",
    "CONTRACT_COLUMNS",
    "TABLES",
    "Column",
    "InputTables",
    "validate_inputs",
    "validate_table",
]

Dtype = Literal["str", "int", "float", "bool", "list"]


@dataclass(frozen=True, slots=True)
class Column:
    """One contract column: name, logical dtype, NA allowed, allowed values, may be absent."""

    name: str
    dtype: Dtype
    nullable: bool = False
    enum: tuple[str, ...] | None = None
    optional: bool = False


def _s(name: str, nullable: bool = False, enum: tuple[str, ...] | None = None) -> Column:
    return Column(name, "str", nullable, enum)


def _i(name: str, nullable: bool = False) -> Column:
    return Column(name, "int", nullable)


def _f(name: str, nullable: bool = False) -> Column:
    return Column(name, "float", nullable)


def _b(name: str, nullable: bool = False) -> Column:
    return Column(name, "bool", nullable)


def _l(name: str) -> Column:
    return Column(name, "list")


STRAND = ("+", "-")
LOCUS_CLASSES = (
    "rna_overlap", "pseudogene", "sorf", "cds_inframe", "nterm_extension", "cterm_extension",
    "cds_outframe", "antisense", "intergenic",
)  # fmt: skip
# Classes whose placements must carry a gene (``cterm_extension``: ADR-0013); ``pseudogene`` may
# (lead decision 2) but need not.
ASSIGNING_CLASSES: tuple[str, ...] = (
    "cds_inframe", "sorf", "nterm_extension", "cterm_extension",
)  # fmt: skip
_MAY_ASSIGN = (*ASSIGNING_CLASSES, "pseudogene")

_LAYER_FLOATS = (
    "promoter_score", "promoter_score_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
    "terminator_conf_downstream", "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg",
)  # fmt: skip
_CODON_FLOATS = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10")

# fmt: off
TABLES: Mapping[str, tuple[Column, ...]] = {
    "replicons": (
        _s("contig"), _i("length"), _s("topology", enum=("circular", "linear")),
        _s("kind", enum=("chromosome", "plasmid", "unknown")),
    ),
    "genes": (
        _s("locus_tag"), _s("gene", True), _s("product", True),
        _s("feature_type", enum=("cds", "sorf", "pseudogene")), _s("contig"), _i("start"),
        _i("end"), _s("strand", enum=STRAND), _b("partial"), _b("crosses_origin"),
        _b("hypothetical"), _s("aa", True), _s("nt"), _s("start_codon", True),
        Column("ribosomal_protein", "bool", nullable=True, optional=True),
        Column("amr", "bool", nullable=True, optional=True),
    ),
    "xrefs": (_s("locus_tag"), _s("db"), _s("id")),
    "tus": (
        _s("tu_id"), _s("contig"), _i("start"), _i("end"), _s("strand", enum=STRAND),
        _l("genes"),
    ),
    "layers": (
        _s("locus_tag"), _s("tu_id"), _i("tu_position"), _i("tu_size"),
        _f("uniop_p_prev", True), _f("uniop_p_next", True), _l("region_ids"),
        *(_f(c, True) for c in _LAYER_FLOATS[:2]), _i("n_promoter_windows_tu", True),
        *(_f(c, True) for c in _LAYER_FLOATS[2:4]), _i("tss_distance_bp", True),
        *(_f(c, True) for c in _LAYER_FLOATS[4:]), _i("spacer_bp", True),
        _s("start_codon", True), *(_f(c, True) for c in _CODON_FLOATS),
        _i("rare_codons_2_10", True), _l("flags"),
    ),
    "regions": (
        _s("region_id"), _s("region_type", enum=("prophage", "plasmid", "bgc", "crispr", "amr")),
        _s("contig"), _i("start", True), _i("end", True), _s("source_tool"), _f("score", True),
        Column("genes", "list", optional=True), Column("notes", "str", True, optional=True),
    ),
    "region_genes": (
        _s("region_id"), _s("locus_tag"), _f("overlap_frac"),
        _s("subset", True, ("integrase", "repressor", "structural", "lysis", "cargo")),
        _b("core_biosynthetic"),
    ),
    "peptides": (
        _s("peptide"), _f("q_value", True), _f("pep", True), _f("score", True), _b("is_decoy"),
        _b("is_contaminant"), _l("protein_ids"), _l("modified_peptides"), _i("n_rows"),
        _i("n_runs"), _f("intensity", True), _s("source", True), _s("search_db", True),
    ),
    "abundance": (_s("protein_id"), Column("raw_id", "str", optional=True), _f("value", True)),
    "peptide_hits": (
        _s("peptide"), _i("peptide_len"), _s("source", True), _f("best_q_input", True),
        _f("pep_input", True), _f("score_input", True), _s("contig"), _i("start"), _i("end"),
        _s("strand", enum=STRAND), _i("frame", True), _b("crosses_origin"),
        _s("route", enum=("exact", "tolerant")), _i("n_diff", True), _l("diff_types"),
        _f("placement_q", True), _s("length_class", True), _s("locus_class", enum=LOCUS_CLASSES),
        _s("uniqueness", enum=("unique", "multi_identical", "shared")), _s("gene", True),
        _s("novel_candidate_id", True), _l("flags"),
    ),
    "protein_hits": (
        _s("protein_id"), _i("protein_len_aa", True), _s("route"), _s("contig"), _i("start"),
        _i("end"), _s("strand", enum=STRAND), _i("rank"), _f("identity", True),
        _f("query_coverage", True), _f("reciprocal_overlap", True),
        _s("class", enum=("identical", "near", "divergent")), _s("gene", True), _l("flags"),
    ),
    "unmapped": (
        _s("sequence_id"), _s("sequence"), _s("kind"), _s("label", True), _s("reason"),
        _s("contaminant_match", True),
    ),
    "id_map": (_s("source_id"), _s("locus_tag")),
}
# fmt: on

#: Column names per table (what giapsim's literal copy must equal).
CONTRACT_COLUMNS: Mapping[str, tuple[str, ...]] = {
    k: tuple(c.name for c in cols) for k, cols in TABLES.items()
}


class InputTables(Protocol):
    """The fields of ``giap.io.work_tables.M4Inputs`` that the cross-table checks read."""

    @property
    def replicons(self) -> pd.DataFrame: ...
    @property
    def genes(self) -> pd.DataFrame: ...
    @property
    def tus(self) -> pd.DataFrame: ...
    @property
    def layers(self) -> pd.DataFrame: ...
    @property
    def regions(self) -> pd.DataFrame: ...
    @property
    def region_genes(self) -> pd.DataFrame: ...
    @property
    def peptides(self) -> pd.DataFrame | None: ...
    @property
    def peptide_hits(self) -> pd.DataFrame: ...
    @property
    def protein_hits(self) -> pd.DataFrame: ...


# --------------------------------------------------------------------------- helpers


def _issue(
    table: str,
    rule: str,
    message: str,
    rows: Sequence[int] = (),
    values: Iterable[object] = (),
    severity: Severity = "error",
) -> Issue:
    """Issue ``contract.<table>.<rule>``; ``rows`` are 0-based positions, ``values`` are named."""
    vals = list(dict.fromkeys(str(v) for v in values))
    if vals:
        shown = ", ".join(vals[:MAX_EXAMPLES]) + (" …" if len(vals) > MAX_EXAMPLES else "")
        message = f"{message}: {shown}"
    if rows:
        message = f"{message} ({len(rows)} rows)"
    examples = tuple(str(r + 1) for r in rows[:MAX_EXAMPLES])
    return Issue(f"contract.{table}.{rule}", severity, message, len(rows), examples)


def _is_listlike(v: object) -> bool:
    return isinstance(v, list | tuple | np.ndarray)


def _is_na(v: object) -> bool:
    if _is_listlike(v):
        return False
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def _na_mask(s: pd.Series) -> np.ndarray:
    return np.asarray(s.isna(), dtype=bool)


def _positions(mask: Iterable[bool]) -> list[int]:
    return [i for i, m in enumerate(mask) if m]


def _is_int_value(v: object) -> bool:
    if isinstance(v, bool | np.bool_):
        return False
    if isinstance(v, numbers.Integral):
        return True
    return isinstance(v, numbers.Real) and np.isfinite(float(v)) and float(v).is_integer()


def _is_float_value(v: object) -> bool:
    return isinstance(v, numbers.Real) and not isinstance(v, bool | np.bool_)


_CELL_OK: dict[str, Callable[[object], bool]] = {
    "str": lambda v: isinstance(v, str),
    "int": _is_int_value,
    "float": _is_float_value,
    "bool": lambda v: isinstance(v, bool | np.bool_),
    "list": _is_listlike,
}


def _bad_dtype_rows(s: pd.Series, dtype: str, na: np.ndarray) -> list[int]:
    """Positions of non-NA cells whose type does not fit ``dtype``."""
    kind = s.dtype
    if dtype != "list" and not pd.api.types.is_object_dtype(kind):
        is_bool = pd.api.types.is_bool_dtype(kind)
        if dtype == "str":
            ok = pd.api.types.is_string_dtype(kind)
        elif dtype == "bool":
            ok = is_bool
        elif dtype == "float":
            ok = pd.api.types.is_numeric_dtype(kind) and not is_bool
        elif pd.api.types.is_integer_dtype(kind):
            ok = True
        elif pd.api.types.is_float_dtype(kind):
            vals = s.to_numpy(dtype="float64", na_value=np.nan)
            return _positions(~na & ~(np.isfinite(vals) & (np.floor(vals) == vals)))
        else:
            ok = False
        return [] if ok else _positions(~na)
    check = _CELL_OK[dtype]
    return [i for i, v in enumerate(s.tolist()) if not na[i] and not check(v)]


def _bad_list_cells(s: pd.Series, na: np.ndarray) -> list[int]:
    """Positions of list cells holding an NA, empty or non-string element."""
    bad = [_is_listlike(v) and any(not isinstance(x, str) or x == "" for x in list(v)) for v in s]
    return _positions(~na & np.asarray(bad, dtype=bool))


def _check_column(table: str, col: Column, s: pd.Series) -> list[Issue]:
    na = _na_mask(s)
    out: list[Issue] = []
    bad = _bad_dtype_rows(s, col.dtype, na)
    if bad:
        got = sorted({type(v).__name__ for v in s.iloc[bad[:50]].tolist()})
        msg = f"column {col.name!r} must be {col.dtype} (dtype {s.dtype}, cell types {got})"
        out.append(_issue(table, "dtype", msg, bad))
    if not col.nullable and na.any():
        out.append(_issue(table, "na", f"column {col.name!r} has missing values", _positions(na)))
    if col.enum is not None:
        enum = set(col.enum)
        rows = [i for i, v in enumerate(s.tolist()) if isinstance(v, str) and v not in enum]
        if rows:
            values = [s.iloc[i] for i in rows]
            msg = f"column {col.name!r} values not in {list(col.enum)}"
            out.append(_issue(table, "enum", msg, rows, values))
    if col.dtype == "list":
        rows = _bad_list_cells(s, na)
        if rows:
            msg = f"column {col.name!r} list elements must be non-empty, non-NA strings"
            out.append(_issue(table, "list_element", msg, rows))
    return out


def _duplicates(table: str, rule: str, s: pd.Series, what: str) -> list[Issue]:
    keys = s.astype(object).where(~_na_mask(s), None)
    dup = np.asarray(keys.duplicated(keep=False) & keys.notna(), dtype=bool)
    if not dup.any():
        return []
    return [_issue(table, rule, f"duplicate {what}", _positions(dup), keys[dup].tolist())]


def _num(s: pd.Series) -> np.ndarray:
    """Numeric view of a column (non-numeric cells become NaN)."""
    return np.asarray(pd.to_numeric(s, errors="coerce"), dtype="float64")


def _truthy(s: pd.Series) -> np.ndarray:
    return np.array([bool(v) if isinstance(v, bool | np.bool_) else False for v in s.tolist()])


def _has(df: pd.DataFrame | None, *cols: str) -> bool:
    return df is not None and all(c in df.columns for c in cols)


# --------------------------------------------------------------------------- per table


def _rules_genes(df: pd.DataFrame) -> list[Issue]:
    out: list[Issue] = []
    if _has(df, "start", "end"):
        start, end = _num(df["start"]), _num(df["end"])
        rows = _positions(~((start >= 0) & (start < end)) & ~np.isnan(start) & ~np.isnan(end))
        if rows:
            out.append(_issue("genes", "interval", "need 0 <= start < end", rows))
    if _has(df, "aa"):
        aa = df["aa"].tolist()
        stop = [i for i, v in enumerate(aa) if isinstance(v, str) and v.endswith("*")]
        if stop:
            out.append(_issue("genes", "aa_stop", "aa must not end with '*'", stop))
        if _has(df, "feature_type"):
            ft = df["feature_type"].tolist()
            miss = [i for i, v in enumerate(aa) if _is_na(v) and ft[i] != "pseudogene"]
            if miss:
                out.append(_issue("genes", "aa_missing", "aa may be NA only for pseudogenes", miss))
    return out


def _rules_regions(df: pd.DataFrame) -> list[Issue]:
    if not _has(df, "region_type", "start", "end"):
        return []
    amr = np.asarray(df["region_type"].astype(object) == "amr", dtype=bool)
    rows = _positions(amr & ~(_na_mask(df["start"]) & _na_mask(df["end"])))
    msg = "AMR regions have NA coordinates"
    return [_issue("regions", "amr_coordinates", msg, rows)] if rows else []


# Columns whose values are unique within the table (rule ``<column>_unique``).
_UNIQUE = dict(replicons="contig", genes="locus_tag", layers="locus_tag", peptides="peptide")
_TABLE_RULES: dict[str, Callable[[pd.DataFrame], list[Issue]]] = {
    "genes": _rules_genes,
    "regions": _rules_regions,
}


def validate_table(name: str, df: pd.DataFrame) -> list[Issue]:
    """Columns, dtypes, enums, NA and single-table invariants of contract table ``name``.

    Missing required columns are errors, extra columns warnings; ``df`` is not modified.
    """
    if name not in TABLES:
        raise ValueError(f"unknown contract table {name!r}; expected one of {sorted(TABLES)}")
    cols = TABLES[name]
    known = {c.name for c in cols}
    out: list[Issue] = []
    for c in cols:
        if c.name not in df.columns:
            if not c.optional:
                out.append(_issue(name, "missing_column", f"missing required column {c.name!r}"))
            continue
        out += _check_column(name, c, df[c.name])
    extra = [str(c) for c in df.columns if c not in known]
    if extra:
        out.append(
            _issue(name, "extra_column", "extra columns ignored", values=extra, severity="warn")
        )
    key = _UNIQUE.get(name)
    if key is not None and key in df.columns:
        out += _duplicates(name, f"{key}_unique", df[key], key)
    rules = _TABLE_RULES.get(name)
    if rules is not None:
        out += rules(df)
    return out


# --------------------------------------------------------------------------- cross-table


def _unknown_refs(
    table: str, rule: str, df: pd.DataFrame, col: str, known: set[str], what: str
) -> list[Issue]:
    if not _has(df, col):
        return []
    vals = df[col].tolist()
    rows = [i for i, v in enumerate(vals) if isinstance(v, str) and v not in known]
    if not rows:
        return []
    return [_issue(table, rule, f"{col} not a known {what}", rows, [vals[i] for i in rows])]


def _check_bounds(genes: pd.DataFrame, reps: dict[str, tuple[float, str]]) -> list[Issue]:
    """``end <= length`` unless ``crosses_origin`` (then circular, ``length < end < 2·length``)."""
    if not _has(genes, "contig", "end", "crosses_origin"):
        return []
    end, cross = _num(genes["end"]), _truthy(genes["crosses_origin"])
    over, origin = [], []
    for i, c in enumerate(genes["contig"].tolist()):
        if not isinstance(c, str) or c not in reps or np.isnan(end[i]):
            continue
        length, topology = reps[c]
        if cross[i] and (topology != "circular" or not length < end[i] < 2 * length):
            origin.append(i)
        elif not cross[i] and end[i] > length:
            over.append(i)
    out = []
    if over:
        out.append(_issue("genes", "bounds", "end beyond replicon length", over))
    if origin:
        msg = "crosses_origin needs a circular replicon and length < end < 2*length"
        out.append(_issue("genes", "origin", msg, origin))
    return out


def _check_tus(t: InputTables, genes: set[str]) -> list[Issue]:
    if not _has(t.tus, "genes"):
        return []
    count: dict[str, int] = {}
    unknown: list[str] = []
    for cell in t.tus["genes"].tolist():
        if not _is_listlike(cell):
            continue
        for g in list(cell):
            if not isinstance(g, str):
                continue
            count[g] = count.get(g, 0) + 1
            if g not in genes:
                unknown.append(g)
    out = []
    missing = sorted(g for g in genes if g not in count)
    multi = sorted(g for g, n in count.items() if n > 1 and g in genes)
    if missing:
        out.append(_issue("tus", "gene_coverage", "genes in no TU", values=missing))
    if multi:
        out.append(_issue("tus", "gene_coverage", "genes in more than one TU", values=multi))
    if unknown:
        out.append(_issue("tus", "unknown_gene", "TU lists unknown genes", values=sorted(unknown)))
    return out


def _check_layers(t: InputTables, genes: set[str]) -> list[Issue]:
    lay = t.layers
    if not _has(lay, "locus_tag"):
        return []
    out = []
    tags = [v for v in lay["locus_tag"].tolist() if isinstance(v, str)]
    missing = sorted(genes - set(tags))
    extra = sorted(set(tags) - genes)
    if missing:
        out.append(_issue("layers", "gene_coverage", "genes without a layers row", values=missing))
    if extra:
        out.append(_issue("layers", "gene_coverage", "layers rows for unknown genes", values=extra))
    if _has(lay, "tu_id") and _has(t.tus, "tu_id", "genes"):
        tu_of = {
            g: tu
            for tu, cell in zip(t.tus["tu_id"].tolist(), t.tus["genes"].tolist(), strict=True)
            if _is_listlike(cell)
            for g in list(cell)
            if isinstance(g, str)
        }
        tus = lay["tu_id"].tolist()
        rows = [
            i
            for i, g in enumerate(lay["locus_tag"].tolist())
            if isinstance(g, str) and g in tu_of and tus[i] != tu_of[g]
        ]
        if rows:
            tags_bad = lay["locus_tag"].iloc[rows].tolist()
            out.append(
                _issue("layers", "tu_id", "tu_id disagrees with tus.parquet", rows, tags_bad)
            )
    return out


def _check_subset(t: InputTables) -> list[Issue]:
    reg, rg = t.regions, t.region_genes
    if not (_has(reg, "region_id", "region_type") and _has(rg, "region_id", "subset")):
        return []
    region_type = dict(zip(reg["region_id"].tolist(), reg["region_type"].tolist(), strict=True))
    subset_na = _na_mask(rg["subset"])
    rows = [
        i
        for i, rid in enumerate(rg["region_id"].tolist())
        if rid in region_type and (region_type[rid] == "prophage") == bool(subset_na[i])
    ]
    msg = "subset is set exactly for prophage members"
    return [_issue("region_genes", "subset", msg, rows)] if rows else []


def _check_peptide_hits(t: InputTables, genes: set[str]) -> list[Issue]:
    h = t.peptide_hits
    out = _unknown_refs("peptide_hits", "unknown_gene", h, "gene", genes, "C2 locus tag")
    if _has(h, "locus_class", "gene"):
        gene_na = _na_mask(h["gene"])
        rows = []
        for i, lc in enumerate(h["locus_class"].tolist()):
            needs = lc in ASSIGNING_CLASSES
            forbids = isinstance(lc, str) and lc not in _MAY_ASSIGN
            if (needs and gene_na[i]) or (forbids and not gene_na[i]):
                rows.append(i)
        if rows:
            msg = f"gene must be set exactly for assigning classes {list(ASSIGNING_CLASSES)}"
            out.append(_issue("peptide_hits", "assignment", msg, rows))
    pep = t.peptides
    if _has(h, "peptide") and _has(pep, "peptide", "is_decoy"):
        assert pep is not None
        decoys = {
            p for p, d in zip(pep["peptide"].tolist(), _truthy(pep["is_decoy"]), strict=True) if d
        }
        vals = h["peptide"].tolist()
        rows = [i for i, p in enumerate(vals) if p in decoys]
        if rows:
            msg = "decoy peptides placed"
            out.append(_issue("peptide_hits", "decoy", msg, rows, [vals[i] for i in rows]))
    return out


def validate_inputs(t: InputTables) -> list[Issue]:
    """Cross-table invariants (C2 vs C1, TU partition, layers, region subsets, mapping genes)."""
    rep = t.replicons
    reps: dict[str, tuple[float, str]] = {}
    if _has(rep, "contig", "length", "topology"):
        for c, n, topo in zip(
            rep["contig"].tolist(), _num(rep["length"]), rep["topology"].tolist(), strict=True
        ):
            if isinstance(c, str) and c not in reps and not np.isnan(n):
                reps[c] = (float(n), str(topo))
    tags = t.genes["locus_tag"].tolist() if _has(t.genes, "locus_tag") else []
    genes = {g for g in tags if isinstance(g, str)}
    out: list[Issue] = []
    out += _unknown_refs("genes", "unknown_contig", t.genes, "contig", set(reps), "replicon")
    out += _check_bounds(t.genes, reps)
    out += _check_tus(t, genes)
    out += _check_layers(t, genes)
    out += _check_subset(t)
    out += _check_peptide_hits(t, genes)
    out += _unknown_refs(
        "protein_hits", "unknown_gene", t.protein_hits, "gene", genes, "C2 locus tag"
    )
    return out
