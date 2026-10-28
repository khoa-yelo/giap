"""``genes.tsv`` assembly: one row per C2 gene (docs/specs/M4-15-gene-table-assembly.md).

GIPPA_SPEC.md §7.3 (one gene table), §9.3 column groups in order, §9.4 Gene flags (closed
vocabulary), §6.5 sort order. Inputs are the M4-06 contract tables (``M4Inputs``) and the outputs of
M4-10 (covariates), M4-11 (evidence), M4-12 (detection), M4-13 (detectability, optional: ``None``
leaves ``p_detect`` NA) and M4-14 (abundance). The frame holds the internal representation:
0-based half-open coordinates, Python lists in list columns, ``None``/NaN/``<NA>`` for missing;
``gippa.io.genes_tsv`` converts it for the TSV. Lead decisions 4, 5, 9: ``gc_local_dev`` copied from
C6 when present; ``cog_category`` from the COG 2020 table (``cog_def``, else NA); several COG IDs ->
the smallest; no C14 id map -> ``source_id = locus_tag``; M5 columns all NA. Everything is pure.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from gippa.genetable.abundance import AbundanceResult
    from gippa.genetable.detectability import DetectabilityResult
    from gippa.genetable.detection import Detection
    from gippa.io.work_tables import M4Inputs

__all__ = ["GENES_COLUMNS", "GENES_TYPES", "GENE_FLAGS", "assemble_genes", "function_columns",
           "gene_flags"]  # fmt: skip

ColumnType = Literal["str", "int", "float", "bool", "list"]

_IDENTITY = (
    "locus_tag", "source_id", "gene", "product", "feature_type", "hypothetical", "contig", "start",
    "end", "strand", "length_aa", "partial", "crosses_origin", "replicon_type",
)  # fmt: skip
_CONTEXT = (
    "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids", "island_id",
    "replichore", "replichore_position", "gc_local_dev",
)  # fmt: skip
_LAYER_VALUES = (
    "promoter_score", "promoter_score_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
    "terminator_conf_downstream", "tir_ostir", "sd_dg", "mrna_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10",
)  # fmt: skip
_PROTEIN = ("mw_da", "pi", "gravy", "tm_segments", "n_observable_peptides")
_DETECTION = ("detected", "detection_tier", "detection_status", "gene_q_value")
_EVIDENCE = (
    "detected_2pep", "n_peptides_unique", "n_peptides_multi_identical", "n_peptides_shared",
    "n_peptides_tolerant", "coverage_frac", "best_peptide_q",
)  # fmt: skip
_ABUNDANCE = (
    "abundance_raw", "abundance_log2", "abundance_imputed", "imputation_policy", "abundance_type",
    "abundance_scale", "top3_log2",
)  # fmt: skip
_FUNCTION = ("cog_id", "cog_category", "go_terms", "kegg_ko", "ec")

#: §9.3 ``tables/genes.tsv`` columns in order (identity, context, regulatory and sequence features,
#: evidence, abundance, function, flags).
GENES_COLUMNS: tuple[str, ...] = (
    *_IDENTITY, *_CONTEXT, *_LAYER_VALUES, *_PROTEIN, *_DETECTION, *_EVIDENCE, "p_detect",
    *_ABUNDANCE, *_FUNCTION, "annotation_source", "flags",
)  # fmt: skip

#: §9.4 flags of scope Gene, in vocabulary order.
GENE_FLAGS: tuple[str, ...] = (
    "partial", "crosses_origin", "pseudogene", "multi_identical", "shared_evidence",
    "two_proteins_one_gene", "cross_strain_protein", "imputed_abundance", "tolerant_evidence",
    "strain_difference", "nterm_extension_evidence", "cai_reference", "top3_lt3",
    "detectability_constant", "abundance_not_explained_by_sequence",
)  # fmt: skip
_FLAG_RANK = {f: i for i, f in enumerate(GENE_FLAGS)}

_LISTS = ("region_ids", "go_terms", "kegg_ko", "ec", "flags")
_INTS = (
    "start", "end", "length_aa", "tu_position", "tu_size", "spacer_bp", "rare_codons_2_10",
    "tm_segments", "n_observable_peptides", "n_peptides_unique", "n_peptides_multi_identical",
    "n_peptides_shared", "n_peptides_tolerant",
)  # fmt: skip
_BOOLS = ("hypothetical", "partial", "crosses_origin", "detected", "detected_2pep",
          "abundance_imputed")  # fmt: skip
_STRS = (
    "locus_tag", "source_id", "gene", "product", "feature_type", "contig", "strand",
    "replicon_type", "tu_id", "island_id", "replichore", "start_codon", "detection_tier",
    "detection_status", "imputation_policy", "abundance_type", "abundance_scale", "cog_id",
    "cog_category", "annotation_source",
)  # fmt: skip


def _type_of(c: str) -> ColumnType:
    if c in _LISTS:
        return "list"
    if c in _INTS:
        return "int"
    if c in _BOOLS:
        return "bool"
    return "str" if c in _STRS else "float"


#: Logical type of every ``GENES_COLUMNS`` column (TSV writer, reader and Parquet schema).
GENES_TYPES: Mapping[str, ColumnType] = {c: _type_of(c) for c in GENES_COLUMNS}

_COG_ID = re.compile(r"^COG\d+$")
_LIST_DBS = {"go_terms": "GO", "kegg_ko": "KEGG", "ec": "EC"}


def gene_flags(row_parts: Mapping[str, Sequence[str]]) -> list[str]:
    """Union of the flag lists in ``row_parts``, unique, in ``GENE_FLAGS`` order.

    A word outside the §9.4 Gene vocabulary raises ``ValueError`` (the vocabulary is closed).
    """
    seen: set[str] = set()
    for part, flags in row_parts.items():
        for f in flags:
            if f not in _FLAG_RANK:
                raise ValueError(f"unknown gene flag {f!r} (from {part}); see GIPPA_SPEC.md §9.4")
            seen.add(f)
    return sorted(seen, key=_FLAG_RANK.__getitem__)


def _is_listlike(v: object) -> bool:
    return isinstance(v, list | tuple | np.ndarray)


def _missing(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v))


def _as_list(v: object) -> list[Any]:
    return list(v) if _is_listlike(v) else []  # type: ignore[call-overload]


def function_columns(xrefs: pd.DataFrame, cog_def: pd.DataFrame | None) -> pd.DataFrame:
    """Pivot of C3 to one row per gene: ``cog_id`` (the smallest ``COG<digits>`` ID; bare
    category xrefs such as ``COG:L`` are ignored), ``cog_category`` (from ``cog_def`` columns
    ``cog_id``/``cog_category``; NA without the table or for an unknown ID) and the sorted distinct
    ``go_terms`` / ``kegg_ko`` / ``ec``. Rows are sorted by ``locus_tag``; other databases ignored.
    """
    cols = ["locus_tag", "cog_id", "cog_category", *_LIST_DBS]
    tags = sorted({str(t) for t in xrefs["locus_tag"]}) if len(xrefs) else []
    db = xrefs["db"].astype(str) if len(xrefs) else pd.Series([], dtype=object)
    ids = xrefs["id"].astype(str) if len(xrefs) else pd.Series([], dtype=object)
    lt = xrefs["locus_tag"].astype(str) if len(xrefs) else pd.Series([], dtype=object)
    cog_mask = (db == "COG") & ids.str.match(_COG_ID.pattern)
    cog_id = ids[cog_mask].groupby(lt[cog_mask]).min()
    cat: dict[str, str] = {}
    if cog_def is not None:
        for k, v in zip(cog_def["cog_id"], cog_def["cog_category"], strict=True):
            if isinstance(k, str) and isinstance(v, str) and v:
                cat.setdefault(k, v)
    data: dict[str, list[Any]] = {"locus_tag": tags}
    data["cog_id"] = [cog_id.get(t) for t in tags]
    data["cog_category"] = [cat.get(c) if c is not None else None for c in data["cog_id"]]
    for col, name in _LIST_DBS.items():
        m = db == name
        per = ids[m].groupby(lt[m]).agg(lambda s: sorted(set(s)))
        data[col] = [list(per.get(t, [])) for t in tags]
    return pd.DataFrame(data, columns=cols)


def _by_tag(df: pd.DataFrame, tags: pd.Index, cols: Sequence[str]) -> pd.DataFrame:
    """``cols`` of ``df`` aligned to ``tags`` by ``locus_tag`` (absent columns become NA)."""
    keyed = df.set_index(df["locus_tag"].astype(str))
    out = {c: keyed[c].reindex(tags) if c in keyed.columns else pd.Series(None, index=tags,
           dtype=object) for c in cols}  # fmt: skip
    return pd.DataFrame(out, index=tags)


def _typed(values: pd.Series, kind: ColumnType) -> pd.Series:
    vals = values.reset_index(drop=True)
    if kind == "list":
        return pd.Series([_as_list(v) for v in vals], dtype=object)
    if kind == "int":
        return pd.to_numeric(vals, errors="coerce").astype("Int64")
    if kind == "float":
        return pd.to_numeric(vals, errors="coerce").astype("float64")
    if kind == "bool":
        return vals.astype("boolean").fillna(False).astype(bool)
    return pd.Series([None if _missing(v) else str(v) for v in vals], dtype=object)


def _sort_order(genes: pd.DataFrame, contig_order: Sequence[str]) -> np.ndarray:
    """§6.5 order: contig (C1 order), start, end, strand, locus_tag."""
    rank = {c: i for i, c in enumerate(contig_order)}
    keys = (
        genes["locus_tag"].astype(str).to_numpy(),
        genes["strand"].astype(str).to_numpy(),
        genes["end"].astype("int64").to_numpy(),
        genes["start"].astype("int64").to_numpy(),
        np.array([rank.get(str(c), len(rank)) for c in genes["contig"]], dtype=np.int64),
    )
    return np.lexsort(keys)


def assemble_genes(
    inputs: M4Inputs,
    cov: pd.DataFrame,
    evidence: pd.DataFrame,
    detection: Detection,
    detectability: DetectabilityResult | None,
    abundance: AbundanceResult,
    *,
    cog_def: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """The §9.3 gene table (columns ``GENES_COLUMNS``), one row per C2 gene, §6.5-sorted.

    Inputs are matched by ``locus_tag`` (their row order does not matter). ``detectability`` None
    leaves ``p_detect`` NA and sets no ``detectability_constant`` flag; with ``constant`` set, the
    flag goes on eligible (non-partial) genes.
    """
    genes = inputs.genes.iloc[_sort_order(inputs.genes, inputs.replicons["contig"].astype(str))]
    tags = pd.Index(genes["locus_tag"].astype(str).to_numpy())
    g = genes.set_index(tags)
    cols: dict[str, pd.Series] = {c: g[c] for c in _IDENTITY if c in g.columns}
    ids = inputs.id_map
    src = pd.Series(tags, index=tags, dtype=object)
    if ids is not None and len(ids):
        m = ids.astype({"locus_tag": str}).drop_duplicates("locus_tag").set_index("locus_tag")
        mapped = m["source_id"].reindex(tags).astype(object)
        src = mapped.where(mapped.notna(), src)
    cols["source_id"] = src
    kind = inputs.replicons.set_index(inputs.replicons["contig"].astype(str))["kind"]
    cols["replicon_type"] = g["contig"].astype(str).map(kind)
    cols["length_aa"] = _by_tag(cov, tags, ["length_aa"])["length_aa"]
    lay = _by_tag(inputs.layers, tags, [*_CONTEXT, *_LAYER_VALUES, "flags"])
    for c in (*_CONTEXT, *_LAYER_VALUES):
        cols[c] = lay[c]
    for c in ("island_id", "replichore", "replichore_position"):  # M5
        cols[c] = pd.Series(None, index=tags, dtype=object)
    for src_df, names in ((cov, _PROTEIN), (detection.table, _DETECTION), (evidence, _EVIDENCE),
                          (abundance.table, _ABUNDANCE)):  # fmt: skip
        part = _by_tag(src_df, tags, names)
        for c in names:
            cols[c] = part[c]
    p = pd.Series(np.nan, index=tags) if detectability is None else detectability.p_detect
    cols["p_detect"] = p.set_axis(p.index.astype(str)).reindex(tags)
    fn = _by_tag(function_columns(inputs.xrefs, cog_def), tags, _FUNCTION)
    for c in _FUNCTION:
        cols[c] = fn[c]
    cols["annotation_source"] = pd.Series(inputs.annotation_source or None, index=tags,
                                          dtype=object)  # fmt: skip
    ev_flags = _by_tag(evidence, tags, ["evidence_flags"])["evidence_flags"]
    ab_flags = _by_tag(abundance.table, tags, ["abundance_flags"])["abundance_flags"]
    constant = detectability is not None and detectability.constant
    partial = g["partial"].astype(bool)
    cols["flags"] = pd.Series(
        [
            gene_flags(
                {
                    "layers": _as_list(a),
                    "evidence": _as_list(b),
                    "abundance": _as_list(c),
                    "detectability": ["detectability_constant"] if constant and not pa else [],
                }
            )
            for a, b, c, pa in zip(lay["flags"], ev_flags, ab_flags, partial, strict=True)
        ],
        index=tags,
        dtype=object,
    )
    return pd.DataFrame({c: _typed(cols[c], GENES_TYPES[c]) for c in GENES_COLUMNS})
