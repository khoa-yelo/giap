"""Builders for the M7-11b class-call tests (docs/specs/M7-11b-class-fdr-wiring.md).

The M7-11 frame layouts and score encodings are reused from ``helpers.novel_fdr_fixture``
(``PickedWorld``: rank r -> ``q_value = r * 1e-5``, ``score = 1000 - r``). On top of them:

* a world genome (seeded random ACGT, four circular 200-kb contigs ``c1..c4``) and its
  ``GeneModel``, so ``class_fdr`` can read replicon lengths and ``nterm_entries``;
* ``gene_table()``: 25 detected, non-partial genes with ``tir_ostir`` 1..25, so the M7-09 RBS
  threshold exists; ORFs get ``ostir_tir = RBS_TIR`` (far above it) and so meet the RBS part of
  their class rule ("every ORF meets its class rule unless stated", spec examples);
* ``assemble``: the step body of spec "novel_orfs step" items 3-5 on in-memory frames:
  ``class_fdr`` -> ``minimal_call`` with ``class_fdr = class_fdr_cell(...)`` -> ``orf_table`` ->
  ``class_rules`` -> ``apply_evidence``; it returns the result, the final table and the rules;
* ``costs_frame``: a ``placement_costs`` frame (``best_costs`` columns, ``PeptideCosts`` order).

Test-author readings (spec gaps, see the PR):

* ``class_fdr(..., denovo=None)`` gives every ``denovo_subgroup`` ORF ``reason =
  denovo_inputs_missing`` (the step passes ``None`` when an input is missing, step item 2);
* the side-table / ``OrfClassFdr.fdr_class`` of an ORF that never entered a picked list is its
  own class.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from giap.config.model import GiapConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel
from giap.novel.candidates import Orf

LENGTH = 200_000
CONTIGS = ("c1", "c2", "c3", "c4")
RBS_TIR = 1.0e6
ALPHA = 0.01
COST_COLUMNS = ("peptide", "length", "source", "search_db", "m_target", "m_decoy", "e_target",
                "e_decoy", "pep", "q_value", "pi_denovo", "is_contaminant")  # fmt: skip
SUMMARY_FIELDS = ("regime", "cls", "estimated", "reason", "n_targets", "n_target_wins",
                  "n_decoys", "paired_share", "c", "q_floor", "floor_above_alpha")  # fmt: skip
NTERM_GENE_AT = 50_000  # b0002 on c1, + strand, [50000, 50300)


def cc() -> Any:
    import giap.novel.class_call as m  # fails (ModuleNotFoundError) until M7-11b lands

    return m


def world_genome(
    contigs: Sequence[str] = CONTIGS, length: int = LENGTH, seed: int = 0
) -> dict[str, str]:
    rng = np.random.default_rng(seed)
    return {c: "".join(rng.choice(list("ACGT"), length)) for c in contigs}


def gene_model(
    genome: Mapping[str, str], genes: Sequence[Gene] = (), topology: str = "circular"
) -> GeneModel:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return GeneModel(replicons=reps, genes=tuple(genes), features=())


def nterm_world() -> tuple[dict[str, str], Gene]:
    """The world genome with gene ``b0002`` on c1 (+) whose in-frame upstream walk finds a stop
    19 codons up: ``TAA (GCT)x19 | ATG (GCT)x98 TAA``."""
    g = world_genome()
    p = NTERM_GENE_AT
    up = "TAA" + "GCT" * 19
    body = "ATG" + "GCT" * 98 + "TAA"
    s = g["c1"]
    g["c1"] = s[: p - len(up)] + up + body + s[p + len(body) :]
    gene = Gene(locus_tag="b0002", feature_type="cds", interval=Interval("c1", p, p + 300, "+"),
                aa="M" + "A" * 98, nt=body, start_codon="ATG", product="test protein")  # fmt: skip
    return g, gene


def config(alpha: float | None = None) -> GiapConfig:
    if alpha is None:
        return GiapConfig()
    return GiapConfig.model_validate({"stats": {"alpha": {"novel_orfs": alpha}}})


def gene_table(n: int = 25) -> pd.DataFrame:
    return pd.DataFrame({
        "locus_tag": [f"g{i:04d}" for i in range(n)],
        "partial": pd.Series([False] * n, dtype="boolean"),
        "detected": pd.Series([True] * n, dtype="boolean"),
        "tir_ostir": np.arange(1.0, n + 1.0),
    })  # fmt: skip


def orf_at(
    orf_id: str,
    cls: str,
    contig: str,
    strand: str,
    stop: int,
    peptides: Sequence[str],
    *,
    start: int | None = None,
    overlapping_gene: str | None = None,
    end_3p_open: bool = False,
) -> Orf:
    """An ORF with one start codon (300 nt upstream of the stop unless ``start`` is given)."""
    if start is None:
        start = stop - 300 if strand == "+" else stop + 300
    return Orf(
        orf_id=orf_id, cls=cls, contig=contig, strand=strand, frame=0,  # type: ignore[arg-type]
        stop=stop, starts=(start,), start_longest=start, end_5p_open=False,
        end_3p_open=end_3p_open,
        peptides=tuple(sorted({p.upper().replace("I", "L") for p in peptides})), best_q=None,
        overlapping_gene=overlapping_gene, first_peptide=start,
    )  # fmt: skip


def costs_frame(
    peptides: Sequence[str],
    *,
    e_target: Sequence[float] | float = 1e-3,
    e_decoy: Sequence[float] | float = math.nan,
    pi: Sequence[float] | float = math.nan,
    source: Sequence[str] | str = "denovo",
) -> pd.DataFrame:
    """A ``placement_costs`` frame (``best_costs`` columns): ``m_* = 0`` where ``e_*`` is finite."""
    n = len(peptides)

    def col(v: Any) -> list[Any]:
        return list(v) if isinstance(v, list | tuple | np.ndarray) else [v] * n

    et, ed = np.asarray(col(e_target), float), np.asarray(col(e_decoy), float)
    df = pd.DataFrame({
        "peptide": list(peptides), "length": [len(p) for p in peptides], "source": col(source),
        "search_db": ["other"] * n, "m_target": np.where(np.isnan(et), np.nan, 0.0),
        "m_decoy": np.where(np.isnan(ed), np.nan, 0.0), "e_target": et, "e_decoy": ed,
        "pep": [math.nan] * n, "q_value": [math.nan] * n, "pi_denovo": col(pi),
        "is_contaminant": [False] * n,
    }, columns=list(COST_COLUMNS))  # fmt: skip
    df["length"] = df["length"].astype("int64")
    for c in ("m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo"):
        df[c] = df[c].astype("float64")
    return df.sort_values("peptide", kind="stable").reset_index(drop=True)


def decoy_hits(peptides: Sequence[str], intervals: Sequence[Interval]) -> pd.DataFrame:
    """``HIT_COLUMNS`` decoy rows (no E column: the mapping hit tables carry none)."""
    from helpers.novel_fdr_fixture import hit_rows

    df = hit_rows(peptides, intervals, [0.0] * len(peptides), search="decoy")
    return df.drop(columns=["e_decoy"])


def target_hits(peptides: Sequence[str], intervals: Sequence[Interval]) -> pd.DataFrame:
    from helpers.novel_fdr_fixture import hit_rows

    df = hit_rows(peptides, intervals, [0.0] * len(peptides), search="target")
    return df.drop(columns=["e_target"])


@dataclass
class Assembled:
    res: Any  # ClassFdrResult
    table: pd.DataFrame  # after apply_evidence
    rules: dict[str, bool]
    pre: pd.DataFrame  # orf_table before apply_evidence

    def row(self, orf_id: str) -> dict[str, Any]:
        rows = self.table[self.table["orf_id"] == orf_id].to_dict("records")
        assert len(rows) == 1, orf_id
        return rows[0]

    def flags(self, orf_id: str) -> list[str]:
        f = self.row(orf_id)["flags"]
        return [] if f is None or f is pd.NA or (isinstance(f, float) and math.isnan(f)) \
            else str(f).split(";")  # fmt: skip

    def cell(self, orf_id: str) -> str:
        return str(self.row(orf_id)["class_fdr"])

    def called(self, orf_id: str) -> bool:
        return bool(self.row(orf_id)["called"])

    def summary(self, regime: str, cls: str) -> Any:
        got = [s for s in self.res.classes if s.regime == regime and s.cls == cls]
        assert len(got) == 1, (regime, cls, self.res.classes)
        return got[0]


def assemble(
    orfs: Sequence[Orf],
    peptides: pd.DataFrame,
    *,
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GiapConfig | None = None,
    denovo: Any = None,
    tir: float = RBS_TIR,
    gt: pd.DataFrame | None = None,
) -> Assembled:
    """Step items 3-5 in memory (no terminators, no peptide intervals)."""
    from giap.novel.candidates import minimal_call
    from giap.novel.evidence import apply_evidence, class_rules
    from giap.novel.outputs import orf_table

    m = cc()
    cfg = cfg or config()
    gt = gene_table() if gt is None else gt
    res = m.class_fdr(orfs, peptides, gm=gm, genome=genome, cfg=cfg, denovo=denovo)
    scored = minimal_call(orfs, cfg)
    scored["class_fdr"] = pd.Series([m.class_fdr_cell(res.orfs[i]) for i in scored["orf_id"]],
                                     dtype="object")  # fmt: skip
    by_id = {o.orf_id: o for o in orfs}
    scored["start_best_ostir"] = pd.Series(
        [by_id[i].starts[0] if by_id[i].starts else None for i in scored["orf_id"]], dtype="object"
    )
    scored["ostir_tir"] = tir
    pre = orf_table(orfs, scored, genome=genome, replicons=gm.replicons)
    rules = class_rules(pre, orfs, gene_table=gt, gm=gm, genome=genome, cfg=cfg,
                        peptide_intervals={})  # fmt: skip
    table = apply_evidence(pre, orfs, gene_table=gt, terminators=None, gm=gm, genome=genome,
                           cfg=cfg, peptide_intervals={})  # fmt: skip
    return Assembled(res, table, dict(rules), pre)
