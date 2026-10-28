"""Scenario builder for M4-12 detection tiers (docs/specs/M4-12-detection-tiers.md).

Inputs are built in the M4-06 internal representation (C2 genes, C9 peptides, C11 peptide hits,
C12 protein hits) with the M4-11 builders of ``evidence_oracle``; ``accepted`` and ``evidence``
come from the merged producers (``gippa.genetable.evidence.accepted_peptides`` /
``gene_evidence``, which also carries ``contaminant_evidence``) and the covariates from M4-10
``covariate_table``, so ``detect`` sees exactly the columns its producers write. Nothing here
re-implements detection.

Gene ``i`` of a scenario lies at ``[OFFSET + STEP i, OFFSET + STEP i + GENE_NT)`` on circular
``chr1``. Target peptides get a C9 row and one exact, unique, ``cds_inframe`` C11 placement in their
gene (``placement_q`` 0, so they pass the five conditions of §7.3.1 whenever q <= q_value_max or
q is NA); the C11 ``best_q_input`` / ``pep_input`` / ``score_input`` mirror the C9 values.
Decoys are C9 rows only (no decoy is ever placed, M4-06 C11 invariant).
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abundance_oracle as ao
import evidence_oracle as eo

GENE_NT = 600
STEP = 700
OFFSET = 100
NAN = math.nan
TIERS = ("decoy_picked", "mayu", "two_peptide", "protein_list")
STATUSES = ("detected", "not_detected", "single_peptide_unestimated", "partial_excluded")
TABLE_COLUMNS = ("locus_tag", "detected", "detection_status", "gene_q_value", "detection_tier")
SUMMARY_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n",
                   "flags")  # fmt: skip
CURVE_THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
_ALPHABET = "ACDEFGHILMNPQSTVWY"  # no K/R inside the body: every generated peptide is tryptic


def is_na(x: Any) -> bool:
    if x is None or x is pd.NA:
        return True
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


def cfg(**kw: Any) -> Any:
    """``GippaConfig`` with nested overrides, e.g. ``cfg(stats={"alpha": {...}})``."""
    from gippa.config import GippaConfig

    return GippaConfig(**kw)


def tag(i: int) -> str:
    return f"GSIM_{5 * (i + 1):05d}"


def peptide_name(k: int, last: str = "K") -> str:
    """Distinct tryptic-looking peptide number ``k`` (10 residues + ``last``)."""
    body = []
    for _ in range(10):
        body.append(_ALPHABET[k % len(_ALPHABET)])
        k //= len(_ALPHABET)
    return "".join(body) + last


@dataclass
class Run:
    """One ``detect`` call with everything it was given."""

    det: Any
    genes: pd.DataFrame
    peptides: pd.DataFrame | None
    peptide_hits: pd.DataFrame
    protein_hits: pd.DataFrame
    accepted: pd.DataFrame
    evidence: pd.DataFrame
    covariates: pd.DataFrame

    def row(self, locus_tag: str) -> dict[str, Any]:
        t = self.det.table
        sel = t[t["locus_tag"].astype(str) == locus_tag]
        assert len(sel) == 1, f"{locus_tag}: {len(sel)} rows in Detection.table"
        return sel.iloc[0].to_dict()

    def q(self, locus_tag: str) -> float:
        v = self.row(locus_tag)["gene_q_value"]
        return NAN if is_na(v) else float(v)

    def status(self, locus_tag: str) -> str:
        return str(self.row(locus_tag)["detection_status"])

    def detected(self, locus_tag: str) -> bool:
        return bool(self.row(locus_tag)["detected"])


@dataclass
class Builder:
    """Hand-built M4 inputs; ``n`` genes ``tag(0..n-1)`` unless ``tags`` is given."""

    n: int = 0
    tags: list[str] | None = None
    partial: set[str] = field(default_factory=set)
    peptides: list[dict[str, Any]] = field(default_factory=list)
    hits: list[dict[str, Any]] = field(default_factory=list)
    protein_hit_rows: list[tuple[Any, ...]] = field(default_factory=list)
    _k: int = 0
    _slot: dict[str, int] = field(default_factory=dict)
    _pos: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.tags is None:
            self.tags = [tag(i) for i in range(self.n)]
        self.n = len(self.tags)

    def _next(self, last: str = "K") -> str:
        self._k += 1
        return peptide_name(self._k, last)

    def gene_start(self, locus_tag: str) -> int:
        assert self.tags is not None
        if len(self._pos) != len(self.tags):
            self._pos = {t: i for i, t in enumerate(self.tags)}
        return OFFSET + STEP * self._pos[locus_tag]

    def target(
        self,
        gene: str | None,
        q: float | None = 0.001,
        pep: float | None = 0.01,
        score: float | None = 50.0,
        *,
        peptide: str | None = None,
        protein_ids: Sequence[str] | None = None,
        contaminant: bool = False,
        uniqueness: str = "unique",
        placed: bool = True,
    ) -> str:
        """A target peptide (C9 row, plus a C11 placement in ``gene`` when ``placed``)."""
        p = peptide if peptide is not None else self._next()
        ids = list(protein_ids) if protein_ids is not None else ([gene] if gene else ["OTHER"])
        self.peptides.append(
            eo.pep(p, q, pep_value=pep, score=score, is_contaminant=contaminant, protein_ids=ids)
        )
        if placed and gene is not None:
            j = self._slot.get(gene, 0)
            self._slot[gene] = j + 1
            h = eo.hit(p, gene, self.gene_start(gene) + 3 * (j % 150), best_q=q,
                       uniqueness=uniqueness)  # fmt: skip
            h["pep_input"], h["score_input"] = pep, score
            self.hits.append(h)
        return p

    def decoy(
        self,
        protein_ids: Sequence[str],
        q: float | None = 0.001,
        pep: float | None = 0.01,
        score: float | None = 50.0,
        *,
        peptide: str | None = None,
    ) -> str:
        p = peptide if peptide is not None else self._next()
        self.peptides.append(
            eo.pep(p, q, pep_value=pep, score=score, is_decoy=True, protein_ids=list(protein_ids))
        )
        return p

    def protein_hit(self, protein_id: str, gene: str | None, rank: int = 1) -> None:
        self.protein_hit_rows.append((protein_id, gene, (), rank))

    # ------------------------------------------------------------------------------- frames

    def genes_frame(self) -> pd.DataFrame:
        assert self.tags is not None
        rows = [
            eo.gene(t, self.gene_start(t), self.gene_start(t) + GENE_NT, partial=t in self.partial)
            for t in self.tags
        ]
        return eo.genes_frame(rows)

    def replicons(self) -> pd.DataFrame:
        return eo.replicons_frame((("chr1", OFFSET + STEP * (self.n + 1), "circular"),))

    def peptides_frame(self) -> pd.DataFrame:
        return eo.peptides_frame(self.peptides)

    def run(
        self,
        config: Any = None,
        *,
        sizes: Mapping[str, float] | float | None = None,
        no_peptides: bool = False,
        shuffle: np.random.Generator | None = None,
    ) -> Run:
        """``detect`` on the scenario; ``no_peptides`` passes C9 = None (proteins-only input);
        ``sizes`` overrides ``n_observable_peptides``; ``shuffle`` permutes C9/C11/C12 rows."""
        from gippa.genetable.covariates import covariate_table
        from gippa.genetable.detection import detect
        from gippa.genetable.evidence import accepted_peptides, gene_evidence

        c = config if config is not None else cfg()
        genes = self.genes_frame()
        peps = None if no_peptides else self.peptides_frame()
        hits = eo.hits_frame(self.hits)
        prot = ao.protein_hits_frame(self.protein_hit_rows)
        if shuffle is not None:
            if peps is not None and len(peps):
                peps = peps.iloc[shuffle.permutation(len(peps))].reset_index(drop=True)
            if len(hits):
                hits = hits.iloc[shuffle.permutation(len(hits))].reset_index(drop=True)
            if len(prot):
                prot = prot.iloc[shuffle.permutation(len(prot))].reset_index(drop=True)
        acc = accepted_peptides(hits, peps, c)
        ev = gene_evidence(genes, hits, acc, self.replicons(),
                           min_peptides=c.detection.min_peptides_sequences_only)  # fmt: skip
        cov = covariate_table(genes, c)
        if sizes is not None:
            if isinstance(sizes, Mapping):
                cov["n_observable_peptides"] = [sizes.get(t, 1) for t in cov["locus_tag"]]
            else:
                cov["n_observable_peptides"] = sizes
        det = detect(genes, ev, acc, peps, hits, prot, cov, c)
        return Run(det, genes, peps, hits, prot, acc, ev, cov)


# ---------------------------------------------------------------------------------- summary


def summary(det: Any) -> pd.DataFrame:
    s = det.summary_rows
    missing = set(SUMMARY_COLUMNS) - set(s.columns)
    assert not missing, f"summary_rows lacks §9.3 columns {sorted(missing)}"
    assert set(s["section"].astype(str)) <= {"tier"}, "M4-12 writes section 'tier' only"
    return s


def tier_metric(det: Any, metric: str) -> float:
    s = summary(det)
    sel = s[(s["group_type"] == "tier") & (s["metric"] == metric)]
    assert len(sel) == 1, f"tier metric {metric}: {len(sel)} rows"
    assert str(sel["group"].iloc[0]) == det.tier
    return float(sel["value"].iloc[0])


def curve_rows(det: Any) -> pd.DataFrame:
    s = summary(det)
    return s[s["group_type"] == "q_threshold"]


def curve_metric(det: Any, threshold: float, metric: str) -> float:
    c = curve_rows(det)
    sel = c[(np.isclose(c["group"].astype(float), threshold)) & (c["metric"] == metric)]
    assert len(sel) == 1, f"curve {threshold}/{metric}: {len(sel)} rows"
    return float(sel["value"].iloc[0])


def mapped(pmap: pd.DataFrame, protein_id: str) -> tuple[str | None, str | None]:
    """(gene, method) of ``protein_id`` in the protein -> gene map; (None, None) if unmapped
    (absent, or present with gene NA)."""
    assert {"protein_id", "gene", "method"} <= set(pmap.columns), list(pmap.columns)
    sel = pmap[pmap["protein_id"].astype(str) == protein_id]
    assert len(sel) <= 1, f"{protein_id}: {len(sel)} rows in the protein -> gene map"
    if sel.empty or is_na(sel["gene"].iloc[0]):
        return None, None
    return str(sel["gene"].iloc[0]), str(sel["method"].iloc[0])


def flat_ids(rows: Iterable[Any]) -> list[str]:
    return [str(x) for r in rows for x in r]
