"""The golden GSEA fixture of M8-13 §1.2 and its test-side oracles (ADR-0038).

- :func:`golden_world`: ``gsea_table``'s inputs from ``tests/data/golden/r/inputs``. One gene
  ``GENE_<tu>`` per TU, detected and eligible, ``abundance_log2 = score`` (so the TU mean is the
  score exactly); one collection whose sets are the golden sets mapped to genes; the
  ``TU_ABSENT_*`` members as genes *absent* from the gene table (they drop as R's name matching
  drops them); one linear replicon.
- :class:`CapturingEngine` records every ``EngineCall`` and delegates (default: the replay).
- :class:`GoldenReplayEngine` returns the golden ``outputs/fgsea.tsv.gz`` rows for one job.
- :class:`NormalisingEngine` (G3a, lead decision 3): sets the seed to GOLDEN_SEED and reorders
  the TUs inside each tie block to golden file order, then delegates. Nothing else changes.
- :func:`oracle_es`: ES and leading edge from the GSEA definition (gseaParam 1, ``std``): the
  ranking is the scores sorted decreasing, ties in input order; a hit steps up by |r| / Σ|r_hits|,
  a miss steps down by 1 / (N - m); ES is the extreme of the running sum (0 if max = -min). No
  fgsea source.
- :func:`log10p_check`: the §3 p rule with the statistician's NA classes (question 2).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from helpers.golden_r_gate import GOLDEN, num, read_gz

GENE = "GENE_"
CONTIG = "golden_chr"


@dataclass(frozen=True)
class GoldenWorld:
    genes: pd.DataFrame
    colls: list[Any]
    replicons: pd.DataFrame
    stats: dict[str, float]  # tu_id -> score, file order
    sets: dict[str, list[str]]  # set_id -> member TUs, file order, absent members kept


def golden_stats() -> dict[str, float]:
    return {r["tu_id"]: float(r["score"]) for r in read_gz(GOLDEN / "inputs" / "gsea_stats.tsv.gz")}


def golden_sets() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for r in read_gz(GOLDEN / "inputs" / "gsea_sets.tsv.gz"):
        out.setdefault(r["set_id"], []).append(r["member"])
    return out


def golden_fgsea() -> dict[str, dict[str, str]]:
    return {r["pathway"]: r for r in read_gz(GOLDEN / "outputs" / "fgsea.tsv.gz")}


def golden_file_order() -> dict[str, int]:
    """TU -> its row index in ``inputs/gsea_stats.tsv.gz`` (the order R saw)."""
    return {t: i for i, t in enumerate(golden_stats())}


def golden_world(collection: str = "golden", shift: float = 0.0) -> GoldenWorld:
    """``shift`` is added to every ``abundance_log2`` (the mutation of example 9); ``stats``
    stays the golden reference."""
    from giap.function.collections import GeneSets

    stats, sets = golden_stats(), golden_sets()
    tus = list(stats)
    genes = pd.DataFrame(
        {
            "locus_tag": [GENE + t for t in tus],
            "tu_id": tus,
            "tu_position": 0,
            "contig": CONTIG,
            "start": np.arange(len(tus), dtype=np.int64) * 1000,
            "end": np.arange(len(tus), dtype=np.int64) * 1000 + 900,
            "strand": "+",
            "partial": False,
            "p_detect": 0.5,
            "detected": True,
            "abundance_log2": [stats[t] + shift for t in tus],
            "abundance_scale": "log2",
        }
    )
    gene_sets = {s: frozenset(GENE + m for m in ms) for s, ms in sets.items()}
    coll = GeneSets(collection, gene_sets, {}, {}, 1.0, 0)  # type: ignore[arg-type]
    replicons = pd.DataFrame({"contig": [CONTIG], "topology": ["linear"]})
    return GoldenWorld(genes, [coll], replicons, stats, sets)


def run_gsea(world: GoldenWorld, engine: Any, seed: int) -> pd.DataFrame:
    """``gsea_table`` on the fixture as the function step calls it (sizes [10, 500], defaults)."""
    from giap.config.model import GseaConfig
    from giap.function.gsea import gsea_table
    from giap.rng import for_step

    return gsea_table(world.genes, world.colls, (10, 500), GseaConfig(), 0.05,
                      for_step("function.gsea", seed), world.replicons, engine=engine)  # fmt: skip


@dataclass
class CapturingEngine:
    """Records every call; ``inner`` serves it (None: the golden replay)."""

    inner: Callable[[Any], Any] | None = None
    method: str = "fgsea_multilevel"
    calls: list[Any] = field(default_factory=list)

    def __call__(self, call: Any) -> Any:
        self.calls.append(call)
        return (self.inner or GoldenReplayEngine())(call)


@dataclass
class GoldenReplayEngine:
    """The golden fgsea rows, in the sent set order, for every job (one job expected)."""

    override: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    method: str = "fgsea_multilevel"
    calls: list[Any] = field(default_factory=list)

    def __call__(self, call: Any) -> Any:
        from dataclasses import replace

        from giap.function.gsea import EngineRow

        self.calls.append(call)
        gold = golden_fgsea()
        out = {}
        for job in call.jobs:
            rows = []
            for term in job.sets:
                g = gold[term]
                edge = tuple(g["leadingEdge"].split(";")) if g["leadingEdge"] else ()
                row = EngineRow(term, int(g["size"]), num(g["ES"]), num(g["NES"]),
                                num(g["pval"]), num(g["log2err"]), edge)  # fmt: skip
                rows.append(replace(row, **dict(self.override.get(term, {}))))
            out[job.job_id] = tuple(rows)
        return out


@dataclass
class NormalisingEngine:
    """G3a's wrapper: ``seed`` replaces the call's seed and every tie block (TUs with equal
    scores) is put in ``order`` (golden file order) on the positions it occupies; scores, jobs
    and parameters are passed through unchanged."""

    inner: Callable[[Any], Any]
    seed: int
    order: Mapping[str, int]
    method: str = "fgsea_multilevel"
    calls: list[Any] = field(default_factory=list)

    def __call__(self, call: Any) -> Any:
        from dataclasses import replace

        tus = list(call.tu_ids)
        blocks: dict[float, list[int]] = {}
        for i, v in enumerate(np.asarray(call.scores, dtype=np.float64).tolist()):
            blocks.setdefault(v, []).append(i)
        for idx in blocks.values():
            if len(idx) > 1:
                block = sorted((tus[i] for i in idx), key=self.order.__getitem__)
                for i, t in zip(idx, block, strict=True):
                    tus[i] = t
        norm = replace(call, tu_ids=tuple(tus), seed=self.seed)
        self.calls.append(norm)
        return self.inner(norm)


def oracle_es(
    tu_ids: Sequence[str], scores: Sequence[float], members: Sequence[str]
) -> tuple[float, frozenset[str]]:
    """(ES, leading-edge TUs) of the set ``members`` on the ranking ``tu_ids``/``scores``."""
    s = np.asarray(scores, dtype=np.float64)
    order = np.argsort(-s, kind="stable")
    ranked = [tu_ids[i] for i in order]
    r = np.abs(s[order])
    hit = np.isin(ranked, list(set(members)))
    m, n = int(hit.sum()), len(ranked)
    nr = float(r[hit].sum())
    step = np.where(hit, r / nr if nr > 0 else 1.0 / m, -1.0 / (n - m))
    run = np.cumsum(step)
    top, bottom = float(run.max()), float(min(run.min(), 0.0))
    if top == -bottom:
        return 0.0, frozenset()
    pos = np.flatnonzero(hit)
    if top > -bottom:
        edge = pos[pos <= int(np.argmax(run))]
        return top, frozenset(ranked[i] for i in edge)
    edge = pos[pos > int(np.argmin(run))]
    return bottom, frozenset(ranked[i] for i in edge)


def log10p_check(p: float, e: float, p_r: float, e_r: float, k: float) -> tuple[str, float]:
    """(class, |Δlog10 p| / bound) of one row under §3 and question 2:

    - ``na_fail``: p NA on exactly one side (a gate failure, no cap);
    - ``na_pval_both``: p NA on both sides, not checked;
    - ``na_log2err_both``: log2err NA on both sides, not checked (ratio NaN);
    - ``na_log2err_one``: one log2err NA, bound k·e_finite / log2(10);
    - ``ok`` / ``exceed``: bound k·sqrt(e² + e_R²) / log2(10)."""
    if math.isnan(p) != math.isnan(p_r):
        return "na_fail", math.nan
    if math.isnan(p):
        return "na_pval_both", math.nan
    if math.isnan(e) and math.isnan(e_r):
        return "na_log2err_both", math.nan
    one = math.isnan(e) or math.isnan(e_r)
    combined = (e_r if math.isnan(e) else e) if one else math.hypot(e, e_r)
    bound = k * combined / math.log2(10.0)
    ratio = abs(math.log10(p) - math.log10(p_r)) / bound if bound > 0 else math.inf
    if one:
        return "na_log2err_one", ratio
    return ("ok" if ratio <= 1.0 else "exceed"), ratio


def na_log2err_over_cap(classes: Sequence[str], compared: int) -> bool:
    """The question-2 alert: both log2err NA classes together exceed ⌊0.05 · compared⌋."""
    n_na = sum(c in ("na_log2err_one", "na_log2err_both") for c in classes)
    return n_na > math.floor(0.05 * compared)
