"""Independent oracle and frame builders for M4-11 accepted evidence (docs/specs/M4-11-...).

Written from the spec text (GIPPA_SPEC.md §7.3.1, §7.2.6, §7.2.7, §9.3, §9.4 and the M4-11 public
interface), not from ``gippa.genetable.evidence``. Coverage is the union of ALL accepted placement
intervals of a gene, repeats of one peptide inside the gene included (M4-11 lead decision 5 after
the PR #74 review); it uses explicit position sets taken
modulo the replicon length, so origin-crossing genes and placements (``end > length``, §6.5) need
no interval splitting here. Nothing in this module imports ``gippa``.

Builders return frames in the M4-06 internal representation (0-based half-open, list cells).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

# M4-06 C11 (M3-16 peptide_hits.tsv) column order.
HIT_COLUMNS: tuple[str, ...] = (
    "peptide",
    "peptide_len",
    "source",
    "best_q_input",
    "pep_input",
    "score_input",
    "contig",
    "start",
    "end",
    "strand",
    "frame",
    "crosses_origin",
    "route",
    "n_diff",
    "diff_types",
    "placement_q",
    "length_class",
    "locus_class",
    "uniqueness",
    "gene",
    "novel_candidate_id",
    "flags",
)
# M1-06 AGG_COLUMNS (C9).
PEPTIDE_COLUMNS: tuple[str, ...] = (
    "peptide",
    "q_value",
    "pep",
    "score",
    "is_decoy",
    "is_contaminant",
    "protein_ids",
    "modified_peptides",
    "n_rows",
    "n_runs",
    "intensity",
    "source",
    "search_db",
)
GENE_COLUMNS: tuple[str, ...] = (
    "locus_tag",
    "gene",
    "product",
    "feature_type",
    "contig",
    "start",
    "end",
    "strand",
    "partial",
    "crosses_origin",
    "hypothetical",
    "aa",
    "nt",
    "start_codon",
)
EVIDENCE_COLUMNS: tuple[str, ...] = (
    "n_peptides_unique",
    "n_peptides_multi_identical",
    "n_peptides_shared",
    "n_peptides_tolerant",
    "coverage_frac",
    "best_peptide_q",
    "detected_2pep",
    "n_accepted",
    "evidence_flags",
)
ACCEPTED_COLUMNS: tuple[str, ...] = (
    "gene",
    "peptide",
    "q_value",
    "pep",
    "score",
    "intensity",
    "uniqueness",
    "route",
    "locus_class",
    "abundance_eligible",
    "start",
    "end",
    "contig",
    "crosses_origin",
)
# §7.2.6 assigning classes; ``pseudogene`` (M4-index lead decision 1) and ``cterm_extension``
# (ADR-0013, M4-06 lead decision 2 after review) included.
ASSIGNING = frozenset({"cds_inframe", "sorf", "nterm_extension", "cterm_extension", "pseudogene"})


# ----------------------------------------------------------------------------- builders


def hit(
    peptide: str,
    gene: str | None,
    start: int,
    end: int | None = None,
    *,
    contig: str = "chr1",
    strand: str = "+",
    locus_class: str = "cds_inframe",
    uniqueness: str = "unique",
    route: str = "exact",
    n_diff: int | None = None,
    placement_q: float | None = 0.0,
    flags: Sequence[str] = (),
    best_q: float | None = 0.001,
    crosses_origin: bool = False,
) -> dict[str, Any]:
    """One C11 row; ``end`` defaults to ``start + 3 * len(peptide)``."""
    e = start + 3 * len(peptide) if end is None else end
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "best_q_input": best_q,
        "pep_input": 0.01,
        "score_input": 50.0,
        "contig": contig,
        "start": start,
        "end": e,
        "strand": strand,
        "frame": start % 3,
        "crosses_origin": crosses_origin,
        "route": route,
        "n_diff": (0 if route == "exact" else 1) if n_diff is None else n_diff,
        "diff_types": [] if route == "exact" else ["substitution"],
        "placement_q": placement_q,
        "length_class": "short",
        "locus_class": locus_class,
        "uniqueness": uniqueness,
        "gene": gene,
        "novel_candidate_id": None,
        "flags": list(flags),
    }


def hits_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    rows = [dict(r) for r in rows]
    if not rows:
        df = pd.DataFrame({c: pd.Series([], dtype=object) for c in HIT_COLUMNS})
    else:
        df = pd.DataFrame(rows, columns=list(HIT_COLUMNS))
    return df.astype(
        {
            "peptide": object,
            "peptide_len": "int64",
            "best_q_input": "float64",
            "pep_input": "float64",
            "score_input": "float64",
            "start": "int64",
            "end": "int64",
            "frame": "Int64",
            "crosses_origin": bool,
            "n_diff": "Int64",
            "placement_q": "float64",
        }
    )


def pep(
    peptide: str,
    q: float | None = 0.001,
    *,
    pep_value: float | None = 0.01,
    score: float | None = 50.0,
    intensity: float | None = 1.0e6,
    is_decoy: bool = False,
    is_contaminant: bool = False,
    protein_ids: Sequence[str] = ("P1",),
    n_rows: int = 1,
    n_runs: int = 1,
    modified: Sequence[str] | None = None,
) -> dict[str, Any]:
    """One C9 row (one stripped peptide after aggregation)."""
    return {
        "peptide": peptide,
        "q_value": q,
        "pep": pep_value,
        "score": score,
        "is_decoy": is_decoy,
        "is_contaminant": is_contaminant,
        "protein_ids": list(protein_ids),
        "modified_peptides": list(modified) if modified is not None else [peptide],
        "n_rows": n_rows,
        "n_runs": n_runs,
        "intensity": intensity,
        "source": "database",
        "search_db": "related",
    }


def peptides_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    rows = [dict(r) for r in rows]
    if not rows:
        df = pd.DataFrame({c: pd.Series([], dtype=object) for c in PEPTIDE_COLUMNS})
    else:
        df = pd.DataFrame(rows, columns=list(PEPTIDE_COLUMNS))
    return df.astype(
        {
            "q_value": "float64",
            "pep": "float64",
            "score": "float64",
            "is_decoy": bool,
            "is_contaminant": bool,
            "n_rows": "int64",
            "n_runs": "int64",
            "intensity": "float64",
        }
    )


def gene(
    locus_tag: str,
    start: int,
    end: int,
    *,
    contig: str = "chr1",
    strand: str = "+",
    feature_type: str = "cds",
    partial: bool = False,
    crosses_origin: bool = False,
) -> dict[str, Any]:
    n_codons = max((end - start) // 3, 2)
    return {
        "locus_tag": locus_tag,
        "gene": None,
        "product": None,
        "feature_type": feature_type,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "partial": partial,
        "crosses_origin": crosses_origin,
        "hypothetical": True,
        "aa": None if feature_type == "pseudogene" else "M" + "A" * (n_codons - 2),
        "nt": "ATG" + "GCT" * (n_codons - 2) + "TAA",
        "start_codon": "ATG",
    }


def genes_frame(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    rows = [dict(r) for r in rows]
    if not rows:
        df = pd.DataFrame({c: pd.Series([], dtype=object) for c in GENE_COLUMNS})
    else:
        df = pd.DataFrame(rows, columns=list(GENE_COLUMNS))
    return df.astype(
        {
            "start": "int64",
            "end": "int64",
            "partial": bool,
            "crosses_origin": bool,
            "hypothetical": bool,
        }
    )


def replicons_frame(
    rows: Sequence[tuple[str, int, str]] = (("chr1", 100_000, "circular"),),
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "contig": [r[0] for r in rows],
            "length": pd.Series([r[1] for r in rows], dtype="int64"),
            "topology": [r[2] for r in rows],
            "kind": ["chromosome"] + ["plasmid"] * (len(rows) - 1),
        }
    )


# ----------------------------------------------------------------------------- oracle


def _na(x: Any) -> bool:
    if x is None:
        return True
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


@dataclass(frozen=True)
class Rules:
    q_value_max: float = 0.01
    alpha_placement: float = 0.01
    genome_relation: str = "related"


def _flags(cell: Any) -> list[str]:
    if isinstance(cell, (list, tuple, np.ndarray)):
        return [str(f) for f in cell]
    return [] if _na(cell) else str(cell).split(";")


def accepted_rows(
    hits: pd.DataFrame, peptides: pd.DataFrame | None, rules: Rules
) -> list[dict[str, Any]]:
    """Every C11 row that is accepted evidence (§7.3.1), with the C9 q joined (C9 wins)."""
    q_of: dict[str, float | None] = {}
    if peptides is not None:
        for p, q in zip(peptides["peptide"], peptides["q_value"], strict=True):
            q_of[str(p)] = None if _na(q) else float(q)
    out = []
    for r in hits.to_dict("records"):
        p = str(r["peptide"])
        if peptides is not None and p not in q_of:
            raise KeyError(p)
        q = (
            q_of.get(p)
            if peptides is not None
            else (None if _na(r["best_q_input"]) else float(r["best_q_input"]))
        )
        if q is not None and not q <= rules.q_value_max:
            continue
        pq = r["placement_q"]
        if not (
            (not _na(pq) and float(pq) <= rules.alpha_placement) or "same_db" in _flags(r["flags"])
        ):
            continue
        if _na(r["gene"]) or r["locus_class"] not in ASSIGNING:
            continue
        if r["uniqueness"] not in ("unique", "multi_identical"):
            continue
        if rules.genome_relation == "same" and r["route"] != "exact":
            continue
        out.append({**r, "q": q})
    return out


def accepted_pairs(
    hits: pd.DataFrame, peptides: pd.DataFrame | None, rules: Rules
) -> set[tuple[str, str]]:
    return {(str(r["gene"]), str(r["peptide"])) for r in accepted_rows(hits, peptides, rules)}


def abundance_eligible(row: Mapping[str, Any], rules: Rules) -> bool:
    return (
        row["uniqueness"] == "unique"
        and row["locus_class"] not in ("nterm_extension", "cterm_extension")  # M4-11 decision 2
        and (
            row["route"] == "exact"
            or (row["route"] == "tolerant" and rules.genome_relation == "related")
        )
    )


def _positions(start: int, end: int, length: int) -> set[int]:
    return {p % length for p in range(int(start), int(end))}


def gene_evidence(
    genes: pd.DataFrame,
    hits: pd.DataFrame,
    peptides: pd.DataFrame | None,
    replicons: pd.DataFrame,
    rules: Rules,
    min_peptides: int = 2,
) -> dict[str, dict[str, Any]]:
    """Per-gene evidence columns (§9.3) and flags (§9.4), keyed by locus tag."""
    length = dict(
        zip(replicons["contig"].astype(str), replicons["length"].astype(int), strict=True)
    )
    acc = accepted_rows(hits, peptides, rules)
    by_gene: dict[str, list[dict[str, Any]]] = {}
    for r in acc:
        by_gene.setdefault(str(r["gene"]), []).append(r)
    shared: dict[str, set[str]] = {}
    for r in hits.to_dict("records"):
        if not _na(r["gene"]) and r["uniqueness"] == "shared":
            shared.setdefault(str(r["gene"]), set()).add(str(r["peptide"]))
    out: dict[str, dict[str, Any]] = {}
    for g in genes.to_dict("records"):
        tag = str(g["locus_tag"])
        rows = by_gene.get(tag, [])
        L = length[str(g["contig"])]
        gpos = _positions(g["start"], g["end"], L)
        cov: set[int] = set()
        for r in rows:  # every accepted placement, not one row per peptide (decision 5)
            if str(r["contig"]) == str(g["contig"]):
                cov |= _positions(r["start"], r["end"], L) & gpos
        n_acc = len({r["peptide"] for r in rows})
        qs = [r["q"] for r in rows if r["q"] is not None]
        n_uni = len({r["peptide"] for r in rows if r["uniqueness"] == "unique"})
        n_multi = len({r["peptide"] for r in rows if r["uniqueness"] == "multi_identical"})
        n_tol = len({r["peptide"] for r in rows if r["route"] == "tolerant"})
        n_sh = len(shared.get(tag, set()))
        flags = set()
        if n_multi:
            flags.add("multi_identical")
        if n_sh:
            flags.add("shared_evidence")
        if n_tol and rules.genome_relation == "related":
            flags |= {"tolerant_evidence", "strain_difference"}
        if any(r["locus_class"] == "nterm_extension" for r in rows):
            flags.add("nterm_extension_evidence")
        out[tag] = {
            "n_peptides_unique": n_uni,
            "n_peptides_multi_identical": n_multi,
            "n_peptides_shared": n_sh,
            "n_peptides_tolerant": n_tol,
            "coverage_frac": len(cov) / (int(g["end"]) - int(g["start"])),
            "best_peptide_q": min(qs) if qs else math.nan,
            "detected_2pep": (not bool(g["partial"])) and n_acc >= min_peptides,
            "n_accepted": n_acc,
            "evidence_flags": flags,
        }
    return out


# ----------------------------------------------------------------------------- realistic fixture


def letters(n: int) -> str:
    """Base-17 code of ``n`` in residue letters (no K/R, so no internal cleavage site)."""
    alphabet = "ACDEFGHILMNPQSTVY"  # no W: W separates codes
    out = alphabet[n % 17]
    n //= 17
    while n:
        out += alphabet[n % 17]
        n //= 17
    return out


def peptide_span(g: Mapping[str, Any], aa_from: int, n_aa: int) -> tuple[int, int]:
    """Genome interval (unwrapped like the gene) of residues ``[aa_from, aa_from + n_aa)``."""
    if g["strand"] == "+":
        s = int(g["start"]) + 3 * aa_from
    else:
        s = int(g["end"]) - 3 * (aa_from + n_aa)
    return s, s + 3 * n_aa


def perfect_evidence(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    detected: Mapping[str, bool],
    rng: np.random.Generator,
    n_false: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """C9 and C11 frames as if search and placement were perfect (M4-09 rules, simplified).

    Every ``detected_true`` gene gets 1-4 distinct peptides (substrings of its protein, or a
    synthetic sequence for pseudogenes without ``aa``) placed at their true locus with
    ``placement_q = 0`` and q <= 0.01; ``n_false`` high-q false targets land on random genes and
    as many decoys (C9 only) are added.
    """
    length = dict(
        zip(replicons["contig"].astype(str), replicons["length"].astype(int), strict=True)
    )
    cls = {"cds": "cds_inframe", "sorf": "sorf", "pseudogene": "pseudogene"}
    hit_rows: list[dict[str, Any]] = []
    pep_rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    recs = genes.to_dict("records")
    for k, g in enumerate(recs):
        aa = g["aa"] if isinstance(g["aa"], str) else None
        n_codons = (int(g["end"]) - int(g["start"])) // 3
        many = int(rng.integers(1, 5)) if detected.get(str(g["locus_tag"]), False) else 0
        for j in range(many):
            n_aa = int(rng.integers(7, 16))
            top = max(min(len(aa) if aa else n_codons, n_codons) - n_aa, 1)
            a0 = int(rng.integers(0, top))
            seq = (aa[a0 : a0 + n_aa] if aa else "") + "W" + letters(k) + "W" + letters(j) + "K"
            if seq in seen:
                continue
            seen.add(seq)
            s, e = peptide_span(g, a0, n_aa)
            L = length[str(g["contig"])]
            if s >= L:
                s, e = s - L, e - L
            hit_rows.append(
                hit(
                    seq,
                    str(g["locus_tag"]),
                    s,
                    e,
                    contig=str(g["contig"]),
                    strand=str(g["strand"]),
                    locus_class=cls[str(g["feature_type"])],
                    crosses_origin=e > L,
                )
            )
            pep_rows.append(
                pep(seq, float(rng.uniform(0.0, 0.009)), protein_ids=[str(g["locus_tag"])])
            )
    for i in range(n_false):
        g = recs[int(rng.integers(0, len(recs)))]
        seq = "FALSEW" + letters(i) + "PEPK"
        s, e = peptide_span(g, 0, len(seq))
        L = length[str(g["contig"])]
        if s >= L:
            s, e = s - L, e - L
        hit_rows.append(
            hit(
                seq,
                str(g["locus_tag"]),
                s,
                e,
                contig=str(g["contig"]),
                strand=str(g["strand"]),
                crosses_origin=e > L,
            )
        )
        pep_rows.append(pep(seq, float(rng.uniform(0.05, 1.0))))
        pep_rows.append(
            pep(
                "DECQYW" + letters(i) + "PEPK",
                float(rng.uniform(0.0, 1.0)),
                is_decoy=True,
                protein_ids=["rev_" + str(g["locus_tag"])],
            )
        )
    return peptides_frame(pep_rows), hits_frame(hit_rows)
