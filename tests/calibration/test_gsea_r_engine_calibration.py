"""M8-17b T11 and A10: the M8-15 ``gsea`` level claims on the R fgsea engine
(docs/specs/M8-17-r-fgsea-engine.md §11 T11; ADR-0038 Consequences: replaces M8-11's multilevel
calibration).

The three level claims of ``test_gsea_tu_calibration.py`` (``gsea_null``, attack 3 operon-packed
sets after the TU collapse, ties at level with the seeded random order) are re-run through
``gsea_table(..., engine=RFgseaEngine)``: the same replicate generators, the same §8.4 level
criterion at every ``levels`` alpha, ``replicates_fast.gsea`` / ``replicates_full.gsea`` from
``tolerances.yaml`` (protected, unchanged). Each replicate makes one R call (one main job, no GO
collection, so no re-test call).

Slurm only: ``slow`` and ``requires_r`` keep these out of GitHub CI. Run with
``pytest tests/calibration -m "calibration and requires_r" -n 8`` inside
``scripts/slurm/calibration_full.sbatch gsea_r`` (8 CPUs, < 2 h for the full suite); the
``giap-fgsea`` env must be built (``sbatch scripts/slurm/build_envs.sbatch fgsea``).
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
from helpers.gsea_data import gene_sets
from test_gsea_tu_calibration import _operon_world, _packed, _tied

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


@pytest.fixture(scope="module")
def engine(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    return r_fgsea_engine(
        resolve_launcher(os.environ), tmp_root=Path(tmp_path_factory.mktemp("fgsea_cal"))
    )


def _p(engine: Any, genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> float:
    from giap.function.gsea import gsea_table

    seed = int(rng.integers(0, 2**31 - 1))
    out = gsea_table(
        genes,
        [gene_sets("kegg", {"S": members})],
        SIZE,
        CFG,
        ALPHA,
        for_step("function.gsea", seed),
        engine=engine,
    )
    assert len(out) == 1, out
    assert set(out["method"]) == {"fgsea_multilevel"}
    p = float(out["p_value"].iloc[0])
    assert np.isfinite(p) and 0 < p <= 1, p
    return p


def _null(engine: Any):
    def rep(rng: np.random.Generator) -> float:
        genes = _operon_world(rng)
        tus = sorted(set(genes["tu_id"]))
        pick = set(rng.choice(tus, 12, replace=False))
        first = genes[genes["tu_id"].isin(pick)].groupby("tu_id")["locus_tag"].first()
        return _p(engine, genes, list(first), rng)

    return rep


def _packed_rep(engine: Any):
    def rep(rng: np.random.Generator) -> float:
        genes, members = _packed(rng)
        return _p(engine, genes, members, rng)

    return rep


def _tied_rep(engine: Any):
    def rep(rng: np.random.Generator) -> float:
        genes, members = _tied(rng)
        return _p(engine, genes, members, rng)

    return rep


def _levels(name: str, fn: Any, full: bool) -> None:
    R = C.replicates(SUITE, full=full)
    levels = C.load_tolerances().levels
    calibrate_levels(SUITE, "level", levels, fn, R, stream=f"{SUITE}.r.{name}").assert_ok()


def test_r_engine_gsea_null_at_level(engine: Any) -> None:
    _levels("null", _null(engine), full=False)


@pytest.mark.calibration_full
def test_r_engine_gsea_null_at_level_full(engine: Any) -> None:
    _levels("null", _null(engine), full=True)


def test_r_engine_attack3_operon_packed_at_level(engine: Any) -> None:
    _levels("attack3_packed", _packed_rep(engine), full=False)


@pytest.mark.calibration_full
def test_r_engine_attack3_operon_packed_at_level_full(engine: Any) -> None:
    _levels("attack3_packed", _packed_rep(engine), full=True)


def test_r_engine_ties_at_level(engine: Any) -> None:
    _levels("ties_min", _tied_rep(engine), full=False)


@pytest.mark.calibration_full
def test_r_engine_ties_at_level_full(engine: Any) -> None:
    _levels("ties_min", _tied_rep(engine), full=True)
