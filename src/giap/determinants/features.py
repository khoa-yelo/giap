"""The determinants feature matrix: one row per eligible gene, five feature blocks
(GIAP_SPEC.md §7.8.1, §7.8.2 step 1, §7.8.3; docs/specs/M8-07-determinants-features.md).

``feature_frame`` builds ``tables/determinants_features.tsv`` (raw values, no transforms);
``model_matrix`` encodes it for M8-08 (natural logs, one-hot levels ``<feature>_<level>``, all-NA
feature columns dropped). Coordinates are 0-based half-open; pure functions, no file I/O.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd

from giap.annotate.bakta_layers import is_ribosomal_protein
from giap.annotate.operons import Tu
from giap.model.coords import Replicon
from giap.model.genes import Gene

__all__ = ["BLOCKS", "FEATURE_COLUMNS", "dist_to_tu_start", "feature_frame", "model_matrix",
           "nterm_residue"]  # fmt: skip

Kind = Literal["str", "int", "float", "bool"]

IDENTITY_COLUMNS: tuple[str, ...] = (
    "locus_tag", "tu_id", "detected", "abundance_log2", "abundance_imputed", "ribosomal_protein",
)  # fmt: skip
#: block -> features of the TSV, §7.8.1 order
FEATURE_BLOCKS: Mapping[str, tuple[str, ...]] = {
    "transcription": (
        "promotech_max", "promotech_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
        "n_promoters_tu_to_gene", "tu_position", "dist_to_tu_start_bp",
    ),
    "translation": ("tir_ostir", "sd_dg", "spacer_bp", "start_codon", "mrna_dg"),
    "codon": ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10"),
    "protein": ("length_aa", "mw", "pi", "gravy", "nterm_residue", "tm_segments"),
    "context": ("replicon_type", "replichore_position", "in_prophage", "in_bgc", "gc_local_dev"),
}  # fmt: skip
#: the TSV columns: identity columns, then the 31 features in §7.8.1 order
FEATURE_COLUMNS: tuple[str, ...] = (
    *IDENTITY_COLUMNS, *(f for fs in FEATURE_BLOCKS.values() for f in fs),
)  # fmt: skip

#: categorical feature -> one-hot levels (every level kept)
LEVELS: Mapping[str, tuple[str, ...]] = {
    "start_codon": ("ATG", "GTG", "TTG", "other"),
    "nterm_residue": tuple("ACDEFGHIKLMNPQRSTVWY"),
    "replicon_type": ("chromosome", "plasmid", "unknown"),
}
LOG_FEATURES: tuple[str, ...] = ("tir_ostir", "length_aa")


def _encoded(feature: str) -> tuple[str, ...]:
    levels = LEVELS.get(feature)
    return (feature,) if levels is None else tuple(f"{feature}_{lv}" for lv in levels)


#: block -> model columns after encoding (the attribution blocks ``shap_<block>``)
BLOCKS: Mapping[str, tuple[str, ...]] = {
    b: tuple(c for f in fs for c in _encoded(f)) for b, fs in FEATURE_BLOCKS.items()
}

#: feature (or identity column) -> gene-table column it is copied from (spec step 2)
GENE_TABLE_SOURCES: Mapping[str, str] = {
    "tu_id": "tu_id", "detected": "detected", "abundance_log2": "abundance_log2",
    "abundance_imputed": "abundance_imputed", "promotech_max": "promoter_score",
    "promotech_tu": "promoter_score_tu", "tx_rate_max_tu": "tx_rate_max_tu",
    "tx_rate_sum_tu": "tx_rate_sum_tu", "tu_position": "tu_position", "tir_ostir": "tir_ostir",
    "sd_dg": "sd_dg", "spacer_bp": "spacer_bp", "start_codon": "start_codon",
    "mrna_dg": "mrna_dg", **{c: c for c in FEATURE_BLOCKS["codon"]}, "length_aa": "length_aa",
    "mw": "mw_da", "pi": "pi", "gravy": "gravy", "tm_segments": "tm_segments",
    "replicon_type": "replicon_type", "gc_local_dev": "gc_local_dev",
}  # fmt: skip

_INTS = ("n_promoters_tu_to_gene", "tu_position", "dist_to_tu_start_bp", "spacer_bp",
         "rare_codons_2_10", "length_aa", "tm_segments")  # fmt: skip
_BOOLS = ("detected", "abundance_imputed", "ribosomal_protein", "in_prophage", "in_bgc")
_STRS = ("locus_tag", "tu_id", "start_codon", "nterm_residue", "replicon_type")
#: TSV column -> logical type
KINDS: Mapping[str, Kind] = {
    c: "int" if c in _INTS else "bool" if c in _BOOLS else "str" if c in _STRS else "float"
    for c in FEATURE_COLUMNS
}

_MET_REMOVED = frozenset("ACGPSTV")
_REGION_FLAGS = {"in_prophage": "prophage", "in_bgc": "bgc"}


def nterm_residue(aa: str) -> str | None:
    """The mature N-terminal residue (§7.8.3, decision 19): residue 2 when residue 1 is M and
    residue 2 ∈ {A, C, G, P, S, T, V}, else residue 1; None below 2 aa."""
    if len(aa) < 2:
        return None
    return aa[1] if aa[0] == "M" and aa[1] in _MET_REMOVED else aa[0]


def _start_codon_pos(g: Gene) -> int:
    iv = g.interval
    return int(iv.start) if iv.strand == "+" else int(iv.end) - 1


def dist_to_tu_start(
    genes: Sequence[Gene], tus: Sequence[Tu], replicons: Mapping[str, Replicon]
) -> dict[str, int]:
    """bp from each TU's first start codon to each member's start codon, in the direction of
    transcription (0 for the first gene; wraps on circular replicons). Genes in no TU (or whose
    TU's first gene is not in ``genes``) have no entry."""
    by_tag = {g.locus_tag: g for g in genes}
    out: dict[str, int] = {}
    for tu in tus:
        first = by_tag.get(tu.genes[0]) if tu.genes else None
        if first is None:
            continue
        rep = replicons[tu.contig]
        p0 = _start_codon_pos(first)
        for tag in tu.genes:
            g = by_tag.get(tag)
            if g is None:
                continue
            d = _start_codon_pos(g) - p0
            if g.interval.strand == "-":
                d = -d
            out[tag] = d % rep.length if rep.topology == "circular" else d
    return out


def _typed(values: Any, kind: Kind, n: int) -> pd.Series:
    vals = pd.Series(values).reset_index(drop=True) if n else pd.Series([], dtype=object)
    if kind == "int":
        return pd.to_numeric(vals, errors="coerce").astype("Int64")
    if kind == "float":
        return pd.to_numeric(vals, errors="coerce").astype("float64")
    if kind == "bool":
        return vals.astype("boolean")
    return pd.Series([v if isinstance(v, str) else None for v in vals], dtype=object)


def _flags(v: Any) -> list[str]:
    if isinstance(v, str):
        return [x for x in v.split(";") if x]
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return []
    return [str(x) for x in v]


def _eligible(tbl: pd.DataFrame) -> np.ndarray:
    """§7.3 regions-table rule: not partial, finite ``p_detect`` (dedupe with M8-03's
    ``giap.genetable.eligible`` once it merges; docs/backlog.md)."""
    p = pd.to_numeric(tbl["p_detect"], errors="coerce").to_numpy(float)
    part = tbl["partial"].astype("boolean").fillna(False).to_numpy(bool)
    return np.asarray(~part & np.isfinite(p), dtype=bool)


def _sorted(tbl: pd.DataFrame, contig_order: Sequence[str]) -> pd.DataFrame:
    """§6.5 order: contig (C1 order), start, end, strand, locus_tag."""
    rank = {c: i for i, c in enumerate(contig_order)}
    keys = (
        tbl["locus_tag"].astype(str).to_numpy(), tbl["strand"].astype(str).to_numpy(),
        tbl["end"].astype("int64").to_numpy(), tbl["start"].astype("int64").to_numpy(),
        np.array([rank.get(str(c), len(rank)) for c in tbl["contig"]], dtype=np.int64),
    )  # fmt: skip
    return tbl.iloc[np.lexsort(keys)].reset_index(drop=True)


def _lookup(frame: pd.DataFrame | None, col: str, tags: pd.Series) -> pd.Series:
    if frame is None or col not in frame.columns:
        return pd.Series([None] * len(tags), dtype=object)
    m = frame.drop_duplicates("locus_tag").set_index(frame["locus_tag"].astype(str))[col]
    return m.reindex(tags.to_numpy()).reset_index(drop=True)


def _ribosomal(tbl: pd.DataFrame, by_tag: Mapping[str, Gene]) -> list[bool]:
    """C2 ``ribosomal_protein`` when set, else product / gene name or the ``cai_reference``
    flag (lead decision 2)."""
    given = tbl["ribosomal_protein"] if "ribosomal_protein" in tbl.columns else None
    out = []
    for i, (tag, flags) in enumerate(zip(tbl["locus_tag"], tbl["flags"], strict=True)):
        v = None if given is None else given.iloc[i]
        if v is not None and not pd.isna(v):
            out.append(bool(v))
            continue
        g = by_tag.get(str(tag))
        out.append((g is not None and is_ribosomal_protein(g)) or "cai_reference" in _flags(flags))
    return out


def feature_frame(
    genes_tbl: pd.DataFrame,
    genes: Sequence[Gene],
    tus: Sequence[Tu],
    layers: pd.DataFrame,
    regions: pd.DataFrame | None,
    gene_context: pd.DataFrame | None,
    replicons: Mapping[str, Replicon],
) -> pd.DataFrame:
    """``determinants_features.tsv`` (columns ``FEATURE_COLUMNS``): one row per eligible gene
    (§7.3) in §6.5 order, raw values. ``regions`` is C8 with the C7 ``region_type`` joined on
    (None: no C7/C8, membership NA); ``gene_context`` is M5-09's (None: ``replichore_position``
    NA); ``replicons`` in C1 order."""
    tbl = _sorted(genes_tbl.loc[_eligible(genes_tbl)], list(replicons))
    tags = tbl["locus_tag"].astype(str)
    by_tag = {g.locus_tag: g for g in genes}
    dist = dist_to_tu_start(genes, tus, replicons)
    cols: dict[str, Any] = {"locus_tag": tags}
    for feat, src in GENE_TABLE_SOURCES.items():
        cols[feat] = tbl[src] if src in tbl.columns else [None] * len(tbl)
    cols["ribosomal_protein"] = _ribosomal(tbl, by_tag)
    cols["n_promoters_tu_to_gene"] = _lookup(layers, "n_promoter_windows_tu", tags)
    cols["dist_to_tu_start_bp"] = [dist.get(t) for t in tags]
    cols["nterm_residue"] = [nterm_residue(g.aa) if (g := by_tag.get(t)) else None for t in tags]
    for feat, rtype in _REGION_FLAGS.items():
        if regions is None:
            cols[feat] = [None] * len(tbl)
        else:
            members = set(regions.loc[regions["region_type"] == rtype, "locus_tag"].astype(str))
            cols[feat] = [t in members for t in tags]
    cols["replichore_position"] = _lookup(gene_context, "replichore_position", tags)
    n = len(tbl)
    return pd.DataFrame({c: _typed(cols[c], KINDS[c], n) for c in FEATURE_COLUMNS})


def _numeric(features: pd.DataFrame, feature: str) -> np.ndarray:
    s = features[feature]
    if KINDS[feature] == "bool":
        s = s.astype("boolean").astype("Float64")
    x = pd.to_numeric(s, errors="coerce").to_numpy(float, na_value=np.nan)
    if feature in LOG_FEATURES:
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(x > 0, np.log(np.where(x > 0, x, 1.0)), np.nan)
    return np.asarray(x, dtype=float)


def _one_hot(features: pd.DataFrame, feature: str) -> dict[str, np.ndarray]:
    levels = LEVELS[feature]
    vals = [v if isinstance(v, str) else None for v in features[feature]]
    if feature == "start_codon":
        vals = [None if v is None else v if v in levels else "other" for v in vals]
    out: dict[str, np.ndarray] = {}
    for lv in levels:
        out[f"{feature}_{lv}"] = np.array(
            [np.nan if v is None else float(v == lv) for v in vals], dtype=float
        )
    return out


def model_matrix(features: pd.DataFrame) -> tuple[pd.DataFrame, Mapping[str, str]]:
    """Encoded model columns (float64, ``BLOCKS`` order) and column -> block.

    ``tir_ostir`` and ``length_aa`` -> ln (≤ 0 -> NaN); ``start_codon``, ``nterm_residue`` and
    ``replicon_type`` one-hot (NA -> NaN in every level); booleans 0/1. A feature whose TSV
    column is NA for every row is dropped with all its encoded columns (spec step 4)."""
    data: dict[str, np.ndarray] = {}
    block_of: dict[str, str] = {}
    for block, feats in FEATURE_BLOCKS.items():
        for feat in feats:
            if bool(features[feat].isna().all()):
                continue
            enc = _one_hot(features, feat) if feat in LEVELS else {feat: _numeric(features, feat)}
            data |= enc
            block_of |= dict.fromkeys(enc, block)
    return pd.DataFrame(data, index=features.index, columns=list(data)), block_of
