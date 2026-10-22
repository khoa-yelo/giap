"""Level of the TU-collapsed GSEA (issue M8-05; GIAP_SPEC.md §7.7.4, §7.7.5, §8.2, §8.4, §8.5
attack 3; M8-02 lead decision 7).

Claims checked through ``giap.function.gsea.gsea_table`` (fgsea simple stage, pass 1), with the
§8.4 level criterion at every ``levels`` alpha, R = ``replicates_fast.gsea`` (full arms:
``replicates_full.gsea``):

1. ``gsea_null`` - abundance independent of the sets, operons of 1-4 genes sharing a TU effect,
   one set of 12 uniformly drawn TUs (one detected gene each): the set p is at level.
2. attack 3 ``gsea_attack_operon_packed`` - every gene of 12 random 5-gene operons in one set
   (60 genes, 12 TUs), genes of one operon strongly correlated: at level after the TU collapse.
   Teeth: the same set ranked at gene level (genes permuted as if independent) is a confirmed
   fail, i.e. permuting gene labels is anti-conservative (§8.2 "resample TUs, not genes").
3. ties - 40 % of TUs imputed to one value (``abundance.missing = min``) and a set of 15
   consecutive TU IDs (a genomic cluster): at level because ties enter fgsea in a seeded random
   order (M8-02 lead decision 7). Teeth: the same data ranked in TU-ID order is a confirmed fail
   (the statistician measured 0.165).

The fgsea stream of each replicate is ``giap.rng.for_step("function.gsea", seed)`` with ``seed``
drawn from the replicate generator.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.gsea_data import Gene, frame, gene_sets

from giap.config.model import GseaConfig
from giap.rng import for_step
from giap.stats.fgsea import fgsea_simple, prepare

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("gsea")]

SUITE = "gsea"
CFG = GseaConfig()  # n_perm_simple 1000, sample_size 101, eps 0
SIZE = (10, 500)
ALPHA = 0.05


def _g() -> Any:
    from giap.function import gsea

    return gsea


def _p(genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> float:
    seed = int(rng.integers(0, 2**31 - 1))
    out = _g().gsea_table(
        genes,
        [gene_sets("kegg", {"S": members})],
        SIZE,
        CFG,
        ALPHA,
        for_step("function.gsea", seed),
    )
    assert len(out) == 1, out
    p = float(out["p_value"].iloc[0])
    assert np.isfinite(p) and 1 / (CFG.n_perm_simple + 1) <= p <= 1, p
    return p


# --------------------------------------------------------------------------- 1. gsea_null


def _operon_world(rng: np.random.Generator, n_tu: int = 300) -> pd.DataFrame:
    tus = {}
    k = 0
    effect = rng.standard_normal(n_tu) + 20.0
    for i in range(n_tu):
        genes = []
        for _ in range(int(rng.integers(1, 5))):
            k += 1
            genes.append(Gene(f"g{k}", float(effect[i] + 0.3 * rng.standard_normal())))
        tus[f"TU_{i + 1}"] = genes
    return frame(tus)


def _null_rep(rng: np.random.Generator) -> float:
    genes = _operon_world(rng)
    tus = sorted(set(genes["tu_id"]))
    pick = set(rng.choice(tus, 12, replace=False))
    first = genes[genes["tu_id"].isin(pick)].groupby("tu_id")["locus_tag"].first()
    return _p(genes, list(first), rng)


def _levels(name: str, fn: Any, full: bool) -> None:
    R = C.replicates(SUITE, full=full)
    levels = C.load_tolerances().levels
    calibrate_levels(SUITE, "level", levels, fn, R, stream=f"{SUITE}.tu.{name}").assert_ok()


def test_gsea_null_at_level() -> None:
    _levels("null", _null_rep, full=False)


@pytest.mark.calibration_full
def test_gsea_null_at_level_full() -> None:
    _levels("null", _null_rep, full=True)


# --------------------------------------------------------------------------- 2. attack 3

N_OP, OP_LEN, N_SET = 200, 5, 12


def _packed(rng: np.random.Generator) -> tuple[pd.DataFrame, list[str]]:
    effect = rng.standard_normal(N_OP) + 20.0
    tus = {
        f"TU_{i + 1}": [
            Gene(f"g{i * OP_LEN + j + 1}", float(effect[i] + 0.3 * rng.standard_normal()))
            for j in range(OP_LEN)
        ]
        for i in range(N_OP)
    }
    sel = rng.choice(N_OP, N_SET, replace=False)
    members = [f"g{int(i) * OP_LEN + j + 1}" for i in sel for j in range(OP_LEN)]
    return frame(tus), members


def _packed_tu_rep(rng: np.random.Generator) -> float:
    genes, members = _packed(rng)
    return _p(genes, members, rng)


def _packed_gene_level_rep(rng: np.random.Generator) -> float:
    """The arm GIAP rejects: genes ranked individually and permuted as if independent."""
    genes, members = _packed(rng)
    x = genes["abundance_log2"].to_numpy(float)
    prep = prepare(x - np.median(x), list(genes["locus_tag"]), {"S": members}, 10, 500)
    seed = int(rng.integers(0, 2**31 - 1))
    res = fgsea_simple(prep, CFG.n_perm_simple, CFG.sample_size, for_step("function.gsea", seed))
    return float(res.p[0])


def test_attack3_operon_packed_sets_at_level_after_tu_collapse() -> None:
    _levels("attack3_packed", _packed_tu_rep, full=False)


@pytest.mark.calibration_full
def test_attack3_operon_packed_sets_at_level_after_tu_collapse_full() -> None:
    _levels("attack3_packed", _packed_tu_rep, full=True)


def test_teeth_attack3_gene_level_ranking_fails(teeth_log: Path) -> None:
    """Gene-label permutation on operon-packed sets is anti-conservative (pilot: 0.375 at 0.05)."""
    rec = C.calibrate(
        SUITE,
        "level",
        ALPHA,
        _packed_gene_level_rep,
        C.replicates(SUITE, full=True),
        stream=f"{SUITE}.tu.teeth_attack3_gene_level",
        log_dir=teeth_log,
    )
    assert rec["status"] == "fail", f"the gene-level arm is not a confirmed fail: {rec}"


# --------------------------------------------------------------------------- 3. ties

N_TIE = 300


def _tied(rng: np.random.Generator) -> tuple[pd.DataFrame, list[str]]:
    x = rng.standard_normal(N_TIE) + 20.0
    imp = rng.random(N_TIE) < 0.4
    x[imp] = x[~imp].min()
    genes = frame(
        {
            f"TU_{i + 1}": [Gene(f"g{i + 1}", float(x[i]), imputed=bool(imp[i]))]
            for i in range(N_TIE)
        }
    )
    start = int(rng.integers(0, N_TIE - 15))
    return genes, [f"g{i + 1}" for i in range(start, start + 15)]


def _tied_rep(rng: np.random.Generator) -> float:
    genes, members = _tied(rng)
    return _p(genes, members, rng)


def _tied_id_order_rep(rng: np.random.Generator) -> float:
    """The order the lead rejected: the TU ranking fed to fgsea in TU-ID order."""
    genes, members = _tied(rng)
    rank = _g().tu_ranking(genes)
    labels = list(rank.tu_ids)
    tu_of = dict(zip(genes["locus_tag"], genes["tu_id"], strict=True))
    prep = prepare(np.asarray(rank.scores), labels, {"S": [tu_of[m] for m in members]}, 10, 500)
    seed = int(rng.integers(0, 2**31 - 1))
    res = fgsea_simple(prep, CFG.n_perm_simple, CFG.sample_size, for_step("function.gsea", seed))
    return float(res.p[0])


def test_ties_at_level_with_seeded_random_order() -> None:
    _levels("ties_min", _tied_rep, full=False)


@pytest.mark.calibration_full
def test_ties_at_level_with_seeded_random_order_full() -> None:
    _levels("ties_min", _tied_rep, full=True)


def test_teeth_ties_in_tu_id_order_fail(teeth_log: Path) -> None:
    """TU-ID tie order with a cluster of consecutive TUs is anti-conservative (pilot: 0.22)."""
    rec = C.calibrate(
        SUITE,
        "level",
        ALPHA,
        _tied_id_order_rep,
        C.replicates(SUITE, full=True),
        stream=f"{SUITE}.tu.teeth_ties_id_order",
        log_dir=teeth_log,
    )
    assert rec["status"] == "fail", f"TU-ID tie order is not a confirmed fail: {rec}"
