"""Hand-built inputs for the M7-09 evidence and ranking tests
(docs/specs/M7-09-novel-evidence-and-ranking.md).

The ORFs come from the merged M7-06 ``reconstruct`` and the M7-07 ``orf_table`` on the genomes of
tests/helpers/novel_outputs_fixture.py, so ``apply_evidence`` is tested on real table rows:

* ``ex1``: intergenic ORF ``NOVEL_c1_43_p``, starts ATG [42, 45) and GTG [60, 63), stop TGA
  [120, 123), two peptides; ``start_best_ostir = 60`` gives 20 codons (spec example 1);
* ``nterm2``: M7-06 example 6 (gene G + [300, 900), TTG [249, 252) upstream of the annotated ATG,
  peptide [255, 285)) with a second distinct peptide [264, 285) (spec example 7);
* ``pseudo``: a Bakta pseudogene P + [540, 900) with an internal TGA [600, 603) and a peptide
  [570, 597) classed ``pseudogene`` (spec example 8 with a homopolymer placed by the test).

Interface pinned by the test author (``giap.novel.evidence``):

* ``apply_evidence(table, orfs, *, gene_table, terminators, gm, genome, cfg,
  peptide_intervals=None)`` -> the table (same rows, same order, ``NOVEL_ORFS_COLUMNS``) with
  ``rbs_supported``, ``promoter_supported``, the terminator columns, ``called``, ``rank`` and
  ``flags`` filled. ``table`` is the M7-07 ``orf_table`` frame; ``gene_table`` needs only
  ``partial``, ``detected`` and ``tir_ostir``; ``terminators`` is the run's terminator layer as
  ``Terminator`` records, ``None`` when the layer is missing; ``peptide_intervals`` maps an
  ``orf_id`` to the genome intervals of its peptides (the pseudogene / cterm labels read only
  these coordinates).
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from giap.annotate.terminators import Terminator
from giap.config.model import GiapConfig
from giap.model.coords import Interval, Replicon
from helpers.novel_oracle import Case, Seq, cds, hit_row, translate
from helpers.novel_outputs_fixture import (
    C1,
    L,
    ex1,
    ex1_hits,
    ex1_seq,
    ex6_put,
    g_gene,
    orfs_of,
    scored,
    table,
)

ID1, ID6 = "NOVEL_c1_43_p", "NOVEL_c1_250_p"
# §9.4 "Novel ORFs" first, then the milestone additions in the order of the novel_orfs schema
# pattern (M7-06/07: no_start_codon, partial, ostir_missing; M7-09: rbs_threshold_unestimated,
# junction_unknown; ADR-0012: promoter_supported, terminator_supported; ADR-0013:
# readthrough_candidate; M7-17 / ADR-0050: denovo_only_unsupported).
FLAG_ORDER: tuple[str, ...] = (
    "junction_spanned", "homopolymer_junction", "class_fdr_unestimated", "rbs_supported",
    "no_start_codon", "partial", "ostir_missing", "rbs_threshold_unestimated", "junction_unknown",
    "promoter_supported", "terminator_supported", "readthrough_candidate",
    "denovo_only_unsupported",
)  # fmt: skip
TIRS_30 = [100.0 * k for k in range(1, 31)]  # spec examples: 10 % quantile = 390
THRESHOLD_30 = 390.0


def ev() -> Any:
    import giap.novel.evidence as m

    return m


# the threshold the M7-09 examples were written for, not the calibrated default (ADR-0022)
_TERMINATORS = {"terminators": {"terminatornet_yes": 0.5}}


def cfg(**novel: Any) -> GiapConfig:
    if novel:
        return GiapConfig(seed=1, novel_orfs=novel, annotation=_TERMINATORS)
    return GiapConfig(seed=1, annotation=_TERMINATORS)


def gene_table(
    tirs: Sequence[float | None],
    *,
    detected: bool | None = True,
    partial: bool = False,
    extra: Sequence[tuple[float | None, bool | None, bool]] = (),
) -> pd.DataFrame:
    """Gene-table rows ``(tir_ostir, detected, partial)``: ``tirs`` with the shared flags, then
    ``extra`` rows as given."""
    rows = [(t, detected, partial) for t in tirs] + list(extra)
    return pd.DataFrame({
        "locus_tag": [f"G{i:04d}" for i in range(len(rows))],
        "partial": pd.array([r[2] for r in rows], dtype="boolean"),
        "detected": pd.array([r[1] for r in rows], dtype="boolean"),
        "tir_ostir": pd.array([float("nan") if r[0] is None else r[0] for r in rows],
                              dtype="float64"),
    })  # fmt: skip


GT30 = gene_table(TIRS_30)
GT15 = gene_table(TIRS_30[:15])


def term(contig: str, start: int, end: int, strand: str, conf: float) -> Terminator:
    return Terminator(Interval(contig, start, end, strand), conf, "TerminatorNet")  # type: ignore[arg-type]


def flags_of(row: Mapping[str, Any]) -> list[str]:
    v = row["flags"]
    return [] if v is None or v is pd.NA or (isinstance(v, float) and v != v) else str(v).split(";")


def is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and v != v)


# ---------------------------------------------------------------- class_evidence rows


def row(**kw: Any) -> dict[str, Any]:
    """A ``novel_orfs`` row (M7-07 table cells, 0-based positions) of spec example 1, plus the two
    keys ``class_evidence`` reads beyond the table: ``start_upstream_nt`` (nt from
    ``start_longest`` to the annotated start of ``overlapping_gene`` in transcription direction,
    positive = upstream; NA without a start or a gene) and ``spans_annotated_stop`` (a peptide
    interval contains the annotated stop codon of ``overlapping_gene``)."""
    base: dict[str, Any] = {
        "orf_id": ID1, "class": "intergenic", "contig": "c1", "strand": "+", "frame": 1,
        "start_longest": 42, "start_best_ostir": 60, "end": 122, "length_aa_longest": 26,
        "length_aa_best": 20, "start_codon_best": "GTG", "n_peptides": 2,
        "peptides": "ALEKRDSPHW;SPHWFNQ", "best_peptide_q": 0.001, "ostir_tir": 500.0,
        "promoter_score": None, "promoter_supported": None, "terminator_downstream": None,
        "terminator_conf": None, "terminator_distance_nt": None, "rbs_supported": None,
        "class_fdr": "unestimated", "called": False, "rank": None, "overlapping_gene": None,
        "flags": "class_fdr_unestimated", "start_upstream_nt": None,
        "spans_annotated_stop": False,
    }  # fmt: skip
    base.update(kw)
    return base


# ---------------------------------------------------------------- apply_evidence cases


def ex1_case(cls: str = "intergenic", *, n_hits: int = 2, gene: str | None = None) -> Case:
    seq = ex1_seq().text()
    hits = ex1_hits(seq, cls=cls)[:n_hits]
    if gene is not None:
        hits = [{**h, "gene": gene} for h in hits]
    genes = [cds(gene, "c1", 42, 75, "+")] if gene is not None else []
    return Case({"c1": seq}, dict(C1), hits, genes=genes)


def nterm2() -> Case:
    """Spec example 7: M7-06 example 6 with two distinct peptides."""
    seq = ex6_put(Seq(L)).text()
    hits = [
        hit_row(translate(seq[255:285]), "c1", 255, 285, "+", "nterm_extension", L, gene="G"),
        hit_row(translate(seq[264:285]), "c1", 264, 285, "+", "nterm_extension", L, gene="G"),
    ]
    return Case({"c1": seq}, dict(C1), hits, genes=[g_gene()])


PSEUDO_JUNCTION = 600


def pseudo_case(run: tuple[int, str] | None = (595, "AAAAA"), *, premature: bool = True) -> Case:
    """Pseudogene P + [540, 900): ATG [540, 543), filler, an internal TGA [600, 603) (or none),
    TAA [897, 900); peptide [570, 597) classed ``pseudogene`` (gene P); ``run`` = a homopolymer
    ``(start, bases)`` written over the genome."""
    s = Seq(L).fill(543, 897).put(540, "ATG").put(897, "TAA")
    if premature:
        s.put(PSEUDO_JUNCTION, "TGA")
    if run is not None:
        s.put(run[0], run[1])
    seq = s.text()
    hits = [
        hit_row(translate(seq[570:597]), "c1", 570, 597, "+", "pseudogene", L, gene="P"),
        hit_row(translate(seq[546:570]), "c1", 546, 570, "+", "pseudogene", L, gene="P"),
    ]
    return Case({"c1": seq}, dict(C1), hits, genes=[cds("P", "c1", 540, 900, "+", "pseudogene")])


def evidence_table(
    case: Case,
    *,
    best: Mapping[str, int | None] | None = None,
    tir: Mapping[str, float] | None = None,
    promoter: Mapping[str, float] | None = None,
    gt: pd.DataFrame = GT30,
    terminators: Sequence[Terminator] | None = (),
    config: GiapConfig | None = None,
    peptide_intervals: Mapping[str, Sequence[Interval]] | None = None,
    orfs: list[Any] | None = None,
) -> pd.DataFrame:
    """``orf_table`` of the case (with OSTIR / promoter scores) through ``apply_evidence``."""
    orfs = orfs if orfs is not None else orfs_of(case)
    cols: dict[str, dict[str, Any]] = {}
    if tir is not None:
        cols["ostir_tir"] = dict(tir)
    if promoter is not None:
        cols["promoter_score"] = dict(promoter)
    t = table(case, orfs, scored(orfs, best=dict(best) if best is not None else None, **cols))
    return ev().apply_evidence(
        t, orfs, gene_table=gt, terminators=None if terminators is None else list(terminators),
        gm=case.gm, genome=case.genome, cfg=config or cfg(),
        peptide_intervals=peptide_intervals,
    )  # fmt: skip


def ex1_table(**kw: Any) -> dict[str, Any]:
    """The single row of example 1 (``start_best_ostir`` GTG 60, TIR 500 unless overridden)."""
    case = kw.pop("case", None) or ex1()
    kw.setdefault("best", {ID1: 60})
    kw.setdefault("tir", {ID1: 500.0})
    (r,) = evidence_table(case, **kw).to_dict("records")
    return r


def replace_orf(orf: Any, **kw: Any) -> Any:
    return dataclasses.replace(orf, **kw)


def circ(n: int = L) -> Replicon:
    return Replicon("c1", n, "circular", "chromosome")


def lin(n: int = L) -> Replicon:
    return Replicon("c1", n, "linear", "plasmid")
