"""Attack nulls for the abundance GSEA on the R fgsea engine (PR #214 statistician review;
GIAP_SPEC.md §8.1 "GSEA: BH per collection", §8.5 attack 3; ADR-0038; M8-17 §11 T11).
Statistician-owned.

``test_gsea_r_engine_calibration.py`` re-runs the three single-set level claims on
``RFgseaEngine``. This file adds what the switch leaves unchecked on the shipped engine:

1. **Teeth on the R path** (confirmed fails through ``calibrate``, private log):
   a. attack 3 ranked at gene level: the operon-packed set of ``_packed`` with every gene its own
      TU, so the TU-label permutation treats the 5 correlated genes of an operon as exchangeable
      units (port pilot: 0.375 at 0.05);
   b. ties handed to R in TU-ID order: the ``_tied`` world sent as one ``EngineCall`` in natural
      TU order instead of the seeded order; R's stable sort keeps that order, so a cluster of
      consecutive tied TUs sits together (port pilot: 0.165-0.22 at 0.05).
   Without these, a level pass of the R arms could not be told from an arm that cannot fail.
2. **BH per collection without the 1/(1 + nPermSimple) floor.** The family evidence in
   docs/methods/function-gsea.md (FWER 0.0015 under the complete null; FDR 0.017 / 0.053 with
   shifted children) was measured on ``fgsea_simple``, whose floor keeps BH from rejecting; the
   note defers the re-check to M8-11, which ADR-0038 replaces. GO-like nested collection
   (1,000 TUs, E. coli operon lengths, rho 0.8; 10 parents = union of 3 children of 10-25 TUs,
   20 standalone sets of 10-60 TUs):
   a. complete null: P(any q <= alpha) at level (FWER = FDR under the complete null);
   b. the first child of 3 parents shifted by +1.5: mean FDP <= q at each level. These parents and
      children are significant together, so the second (GO re-test) R call with the shared seed
      runs live; the test asserts it did.
3. **The gated structure arms** of ``test_attack_gsea_structure.py`` (median imputation, size
   biased membership, scattered between-TU AR(1)) on the R engine, same generators and criterion.

Slurm only (``slow``, ``requires_r``); the ``giap-fgsea`` env must be built.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.gsea_data import Gene, frame, gene_sets
from test_attack_gsea_structure import _first_genes, _frame, _world
from test_gsea_tu_calibration import _packed, _tied

from giap.config.model import GseaConfig
from giap.rng import for_step

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calsuite("gsea"),
    pytest.mark.slow,
    pytest.mark.requires_r,
]

SUITE = "gsea"
CFG = GseaConfig()
SIZE = (10, 500)
ALPHA = 0.05
FAST_FULL = [
    pytest.param(False, id="fast"),
    pytest.param(True, id="full", marks=pytest.mark.calibration_full),
]


@pytest.fixture(scope="module")
def engine(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    return r_fgsea_engine(
        resolve_launcher(os.environ), tmp_root=Path(tmp_path_factory.mktemp("fgsea_attack"))
    )


class _Counting:
    """Wraps an engine and counts its calls that carry GO re-test jobs."""

    def __init__(self, inner: Any) -> None:
        self.inner = inner
        self.method = inner.method
        self.retest_calls = 0

    def __call__(self, call: Any) -> Any:
        if any(j.job_id.startswith("retest:") for j in call.jobs):
            self.retest_calls += 1
        return self.inner(call)


def _table(engine: Any, genes: pd.DataFrame, colls: list[Any], rng: np.random.Generator) -> Any:
    from giap.function.gsea import gsea_table

    seed = int(rng.integers(0, 2**31 - 1))
    out = gsea_table(genes, colls, SIZE, CFG, ALPHA, for_step("function.gsea", seed),
                     engine=engine)  # fmt: skip
    assert len(out) > 0 and set(out["method"]) == {engine.method}, out
    return out


def _p(engine: Any, genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> float:
    out = _table(engine, genes, [gene_sets("kegg", {"S": members})], rng)
    assert len(out) == 1, out
    p = float(out["p_value"].iloc[0])
    assert np.isfinite(p) and 0 < p <= 1, p
    return p


# --------------------------------------------------------------------------- 1. teeth


def _gene_level(engine: Any) -> Any:
    def rep(rng: np.random.Generator) -> float:
        genes, members = _packed(rng)
        pairs = zip(genes["locus_tag"], genes["abundance_log2"], strict=True)
        solo = {f"TU_{k + 1}": [Gene(str(g), float(v))] for k, (g, v) in enumerate(pairs)}
        return _p(engine, frame(solo), members, rng)

    return rep


def _tied_id_order(engine: Any) -> Any:
    def rep(rng: np.random.Generator) -> float:
        from giap.function.gsea import EngineCall, GseaJob, tu_ranking

        genes, members = _tied(rng)
        rank = tu_ranking(genes)
        pos = {t: i for i, t in enumerate(rank.tu_ids)}
        tu_of = dict(zip(genes["locus_tag"], genes["tu_id"], strict=True))
        s = tuple(sorted({tu_of[m] for m in members}, key=pos.__getitem__))
        seed = int(rng.integers(1, 2**31 - 1))
        job = GseaJob("main:kegg", {"S": s}, SIZE[0], SIZE[1])
        call = EngineCall(tuple(rank.tu_ids), np.asarray(rank.scores, dtype=np.float64), (job,),
                          CFG, seed)  # fmt: skip
        rows = engine(call)["main:kegg"]
        assert len(rows) == 1, rows
        return float(rows[0].p)

    return rep


def test_teeth_r_attack3_gene_level_ranking_fails(engine: Any, teeth_log: Path) -> None:
    rec = C.calibrate(SUITE, "level", ALPHA, _gene_level(engine), C.replicates(SUITE, full=False),
                      stream=f"{SUITE}.r.teeth_attack3_gene_level", log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", f"gene-level ranking on R is not a confirmed fail: {rec}"


@pytest.mark.calibration_full
def test_teeth_r_ties_in_tu_id_order_fail(engine: Any, teeth_log: Path) -> None:
    rec = C.calibrate(SUITE, "level", ALPHA, _tied_id_order(engine), C.replicates(SUITE, full=True),
                      stream=f"{SUITE}.r.teeth_ties_id_order", log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", f"TU-ID tie order on R is not a confirmed fail: {rec}"


# --------------------------------------------------------------------------- 2. BH per collection


def _go_world(rng: np.random.Generator, shift: float) -> tuple[pd.DataFrame, Any, set[str]]:
    n = 1000
    tu, x = _world(rng, n, "ecoli", 0.8)
    x = x.copy()
    order = rng.permutation(n)
    taken = 0
    sets: dict[str, list[str]] = {}
    parents: dict[str, list[str]] = {}
    shifted: set[str] = set()
    for p in range(10):
        kids = []
        for c in range(3):
            k = int(rng.integers(10, 26))
            pick, taken = order[taken : taken + k], taken + k
            name = f"GO:{p:02d}{c}"
            sets[name] = _first_genes(tu, pick)
            parents[name] = [f"GO:{p:02d}P"]
            kids.append(name)
            if shift and c == 0 and p < 3:
                shifted |= set(sets[name])
        sets[f"GO:{p:02d}P"] = [g for kid in kids for g in sets[kid]]
        parents[f"GO:{p:02d}P"] = []
    for s in range(20):
        pick = rng.choice(n, int(rng.integers(10, 61)), replace=False)
        sets[f"GO:S{s:02d}"] = _first_genes(tu, pick)
        parents[f"GO:S{s:02d}"] = []
    if shift:
        idx = np.array(sorted(int(g[1:]) - 1 for g in shifted))
        x[idx] += shift
    null = {t for t, gs in sets.items() if not set(gs) & shifted}
    return _frame(tu, x), gene_sets("go_bp", sets, parents=parents), null


def _fdp(engine: Any, shift: float, levels: list[float]) -> Any:
    def rep(rng: np.random.Generator) -> list[float]:
        genes, coll, null = _go_world(rng, shift)
        out = _table(engine, genes, [coll], rng)
        q = out["q_value"].to_numpy(dtype=float)
        terms = list(out["term"])
        v = []
        for a in levels:
            rej = q <= a  # NaN q never rejects
            fd = sum(1 for t, r in zip(terms, rej, strict=True) if r and t in null)
            v.append(fd / max(int(rej.sum()), 1))
        return v

    return rep


@pytest.mark.parametrize("full", FAST_FULL)
def test_r_bh_per_collection_complete_null_fwer(engine: Any, full: bool) -> None:
    levels = list(C.load_tolerances().levels)
    reps = C.replicates(SUITE, full=full)
    calibrate_levels(SUITE, "fwer", levels, _fdp(engine, 0.0, levels), reps,
                     stream=f"{SUITE}.r.attack.bh_go_null").assert_ok()  # fmt: skip


@pytest.mark.parametrize("full", FAST_FULL)
def test_r_bh_per_collection_fdr_with_shifted_children(engine: Any, full: bool) -> None:
    levels = list(C.load_tolerances().levels)
    reps = C.replicates(SUITE, full=full)
    counted = _Counting(engine)
    gate = calibrate_levels(SUITE, "fdr", levels, _fdp(counted, 1.5, levels), reps,
                            stream=f"{SUITE}.r.attack.bh_go_shift1.5")  # fmt: skip
    gate.assert_ok()
    if gate.values:  # replicates ran here (not a skipped confirm-mode check)
        assert counted.retest_calls > 0, "no replicate reached the GO re-test call"


# --------------------------------------------------------------------------- 3. structure arms


def _median_imputed(engine: Any) -> Any:
    def rep(rng: np.random.Generator) -> float:
        n = 600
        tu, x = _world(rng, n, "1-4", 0.8)
        imp = rng.random(tu.size) < 0.2
        x = x.copy()
        x[imp] = np.median(x[~imp])
        a = int(rng.integers(0, n - 15))
        return _p(engine, _frame(tu, x, imp=imp), _first_genes(tu, np.arange(a, a + 15)), rng)

    return rep


def _size_biased(engine: Any) -> Any:
    def rep(rng: np.random.Generator) -> float:
        tu, x = _world(rng, 1500, "ecoli", 0.2)
        g = rng.choice(tu.size, 40, replace=False)
        return _p(engine, _frame(tu, x), [f"g{int(i) + 1}" for i in g], rng)

    return rep


def _ar1_scattered(engine: Any) -> Any:
    def rep(rng: np.random.Generator) -> float:
        n = 1500
        tu, x = _world(rng, n, "ecoli", 0.8, phi=0.5)
        pick = rng.choice(n, 15, replace=False)
        return _p(engine, _frame(tu, x), _first_genes(tu, pick), rng)

    return rep


ARMS = {"median_imputed": _median_imputed, "size_biased": _size_biased,
        "ar1_scattered": _ar1_scattered}  # fmt: skip


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("arm", sorted(ARMS))
def test_r_attack_gsea_structure_at_level(engine: Any, arm: str, full: bool) -> None:
    levels = C.load_tolerances().levels
    reps = C.replicates(SUITE, full=full)
    calibrate_levels(SUITE, "level", levels, ARMS[arm](engine), reps,
                     stream=f"{SUITE}.r.attack.{arm}").assert_ok()  # fmt: skip
