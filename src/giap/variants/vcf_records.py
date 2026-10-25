"""VCF records and the ``variants.tsv`` table from M7-02 sites (GIAP_SPEC.md §7.5.5, G14, B.7;
docs/specs/M7-03-variant-tables-and-vcf.md rules 1-8).

Substitutions and equal-length two-codon events (swap, 2<->2 block) become records; 1<->2 and
2<->1 blocks and frameshifts stay in the table only (``alt_codon_vcf`` NA). A REF allele crossing
the origin of a circular replicon is flagged ``vcf_unrepresentable`` and gets no record.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from functools import cache
from itertools import product
from typing import Any

import pandas as pd

from giap.io.vcf import VcfHeaderField, VcfRecord
from giap.mapping.sixframe import translate11_text
from giap.model.coords import Replicon

__all__ = [
    "FORMAT_FIELDS",
    "INFO_FIELDS",
    "VARIANTS_COLUMNS",
    "VARIANT_FLAGS",
    "minimal_alt_codon",
    "vcf_records",
]

VARIANTS_COLUMNS: tuple[str, ...] = (
    "variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon", "ref_aa",
    "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "nt_changes", "alt_codon_vcf",
    "class_error_rate", "class_error_n", "source", "flags",
)  # fmt: skip
_FLAG_ORDER: tuple[str, ...] = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon", "vcf_unrepresentable",
)  # fmt: skip
VARIANT_FLAGS: frozenset[str] = frozenset(_FLAG_ORDER)
UNESTIMATED = "unestimated"

_F = VcfHeaderField
INFO_FIELDS: tuple[VcfHeaderField, ...] = (
    _F("GENE", "1", "String", "Locus tag of the gene framing the codon"),
    _F("AA_POS", "1", "Integer", "1-based residue position of the codon in the gene"),
    _F("AA_REF", "1", "String", "Reference residue(s) on the gene strand"),
    _F("AA_ALT", "A", "String", "Alternative residue(s) observed in peptides (J = I or L)"),
    _F("CALL", "1", "String", "Codon call: ALT or mixed"),
    _F("CONF", "A", "String", "Confidence of each ALT allele: high or candidate"),
    _F("NT_CHANGES", "A", "Integer", "Nucleotide changes of the minimal-change ALT codon(s)"),
    _F("ALT_TYPE", "A", "String", "ALT type: substitution, swap or block"),
    _F("FLAGS", ".", "String", "Variant flags of the alleles, deduplicated, in GIAP §9.4 order"),
    _F("CLASS_ERR", "A", "String", "Class-specific error rate of each allele, or unestimated"),
)
FORMAT_FIELDS: tuple[VcfHeaderField, ...] = (
    _F("GT", "1", "String", "Haploid genotype: the best-supported ALT allele; . for mixed or tie"),
    _F("AD", "R", "Integer", "Distinct peptides supporting REF, then each ALT allele"),
    _F("DP", "1", "Integer", "Distinct eligible peptides covering the codon"),
)

_COMP = str.maketrans("ACGTN", "TGCAN")
_RECODED = {"U": "TGA", "O": "TAG"}  # M7-16 rule 11: only for the zero-change check


def _rc(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


@cache
def _codons_of(aa: str) -> tuple[str, ...]:
    """Table-11 codons encoding ``aa`` (``J``: every I or L codon), sorted."""
    want = "IL" if aa == "J" else aa
    return tuple(c for c in map("".join, product("ACGT", repeat=3)) if translate11_text(c) in want)


def minimal_alt_codon(ref_codon: str, alt_aa: str) -> tuple[str, int]:
    """G14: the codon for ``alt_aa`` with the fewest differences from ``ref_codon`` (gene strand),
    ties to the lexicographically smallest; returns ``(codon, n_changes)``."""
    codons = (_RECODED[alt_aa],) if alt_aa in _RECODED else _codons_of(alt_aa)
    if not codons:
        raise ValueError(f"no table-11 codon encodes {alt_aa!r}")
    n, codon = min((sum(a != b for a, b in zip(c, ref_codon, strict=True)), c) for c in codons)
    return codon, n


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _flags(v: object) -> tuple[str, ...]:
    items: tuple[str, ...] = () if _na(v) else tuple(v)  # type: ignore[arg-type]
    bad = sorted(set(items) - VARIANT_FLAGS)
    if bad:
        raise ValueError(f"flags outside the Variants vocabulary (§9.4): {bad}")
    return _ordered(items)


def _ordered(items: Iterable[str]) -> tuple[str, ...]:
    """Deduplicated flags in the §9.4 order (M7-05 rule 7)."""
    seen = set(items)
    return tuple(f for f in _FLAG_ORDER if f in seen)


def _class_err(v: object) -> str:
    if _na(v) or v == UNESTIMATED:
        return UNESTIMATED
    return f"{float(v):.6g}"  # type: ignore[arg-type]


def _dp(pile: pd.DataFrame) -> dict[tuple[str, int, str], int]:
    """Distinct eligible peptides per pileup codon: ``n_ref + n_il + Σ n_alt`` (decision 4)."""
    out = {}
    for r in pile[["contig", "codon_start", "strand", "n_ref", "n_il", "alt_counts"]].itertuples():
        alts = "" if _na(r.alt_counts) else str(r.alt_counts)
        n = sum(int(x.rpartition(":")[2]) for x in alts.split(";") if x)
        out[(str(r.contig), int(r.codon_start), str(r.strand))] = int(r.n_ref) + int(r.n_il) + n
    return out


def _record(rows: list[dict[str, Any]], ref: str, pos0: int, dp: int) -> VcfRecord:
    rows = sorted(rows, key=lambda r: (-int(r["n_alt"]), str(r["alt_aa"])))
    first = rows[0]
    n_alt = [int(r["n_alt"]) for r in rows]
    best = max(n_alt)
    gt = "." if first["call"] == "mixed" or n_alt.count(best) > 1 else str(n_alt.index(best) + 1)
    gene_ref = _rc(ref) if first["strand"] == "-" else ref
    info: dict[str, Any] = {}
    if not _na(first["gene"]):
        info["GENE"] = str(first["gene"])
    if not _na(first["aa_pos"]):
        info["AA_POS"] = int(first["aa_pos"])
    info |= {
        "AA_REF": translate11_text(gene_ref),
        "AA_ALT": [str(r["alt_aa"]) for r in rows],
        "CALL": str(first["call"]),
        "CONF": [str(r["confidence"]) for r in rows],
        "NT_CHANGES": [int(r["nt_changes"]) for r in rows],
        "ALT_TYPE": [str(r["alt_type"]) for r in rows],
    }
    flags = list(_ordered(f for r in rows for f in r["flags"]))
    if flags:
        info["FLAGS"] = flags
    info["CLASS_ERR"] = [str(r["class_error_rate"]) for r in rows]
    alts = tuple(str(r["alt_codon_vcf"]) for r in rows)
    ad = (int(first["n_ref"]), *n_alt)
    return VcfRecord(str(first["contig"]), pos0, ref, alts, info, gt, ad, dp)


def vcf_records(
    sites: pd.DataFrame,
    pile: pd.DataFrame,
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
) -> tuple[list[VcfRecord], pd.DataFrame]:
    """VCF records (rules 2-8) and the ``variants.tsv`` table, one row per ``sites`` row in order.

    Optional ``sites`` columns from M7-05/M7-10: ``flags`` (tuples), ``class_error_rate`` (float
    or ``unestimated``) and ``class_error_n``.
    """
    dp = _dp(pile)
    groups: dict[tuple[str, int, str, int, int], list[dict[str, Any]]] = {}
    rows: list[dict[str, Any]] = []
    for s in sites.to_dict("records"):
        row = {c: s.get(c, pd.NA) for c in VARIANTS_COLUMNS}
        row["flags"] = _flags(s.get("flags"))
        row["class_error_rate"] = _class_err(s.get("class_error_rate"))
        ref_nt, alt = str(s["ref_codon"]), str(s["alt_aa"])
        width = len(ref_nt)
        row["nt_changes"] = row["alt_codon_vcf"] = pd.NA
        rows.append(row)
        if s["alt_type"] == "frameshift" or width != 3 * len(alt):
            continue
        cods = [minimal_alt_codon(ref_nt[3 * i : 3 * i + 3], a) for i, a in enumerate(alt)]
        row["nt_changes"] = sum(n for _, n in cods)
        contig, codon, strand = str(s["contig"]), int(s["codon_start"]), str(s["strand"])
        rep = replicons[contig]
        low = codon if strand == "+" or width == 3 else codon - 3
        low = low % rep.length if rep.topology == "circular" else low
        if low + width > rep.length or row["nt_changes"] == 0:  # origin, or U/O over its stop
            row["flags"] = _ordered((*row["flags"], "vcf_unrepresentable"))
            continue
        gene_alt = "".join(c for c, _ in cods)
        row["alt_codon_vcf"] = _rc(gene_alt) if strand == "-" else gene_alt
        groups.setdefault((contig, codon, strand, width, low), []).append(row)
    recs = []
    for (contig, codon, strand, width, low), grp in groups.items():
        ref = genome[contig][low : low + width].upper()
        recs.append(_record(grp, ref, low, dp[(contig, codon, strand)]))
    table = pd.DataFrame(rows, columns=list(VARIANTS_COLUMNS))
    for col in ("aa_pos", "nt_changes", "class_error_n"):
        table[col] = table[col].astype("Int64")
    return recs, table
