"""M8-16 review round, code review B2c (lead decision 5): the ``determinants`` family takes the
``detected`` calls from the run's gene table, not from the M8-07 features table
(docs/specs/M8-16-sensitivity-function-determinants.md item 2 "with ``detected``,
``abundance_log2`` and ``abundance_imputed`` taken from the run's gene table").

By gap 11, ``single_peptide_detected`` is one of the two paths on which the determinants family
can change: newly detected genes must enter the model rows. Here 60 genes have an abundance and
are detected in the gene table but not in the features table (the primary's calls); the rerun
must fit on all of them.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_function_data as sf

N_FLIP = 60


def _rng() -> np.random.Generator:
    from gippa.rng import for_step

    return for_step(f"sensitivity/{sf.DET}", sf.SEED)


@pytest.fixture(scope="module")
def world() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes, features): every gene detected in ``genes``; the first ``N_FLIP`` not detected
    in ``features`` (their abundance kept)."""
    genes, feats = sf.planted_world()
    feats = feats.copy()
    det = feats["detected"].astype("boolean").copy()
    det.iloc[:N_FLIP] = False
    feats["detected"] = det
    return genes, feats


def test_with_gene_calls_takes_detected_from_genes(world: tuple[pd.DataFrame, Any]) -> None:
    from gippa.determinants.sensitivity import with_gene_calls

    genes, feats = world
    out = with_gene_calls(feats, genes)
    assert out["detected"].astype(bool).all()
    off = genes.copy()
    off["detected"] = [i % 3 != 0 for i in range(len(off))]
    got = with_gene_calls(feats, off)["detected"].astype(bool).tolist()
    assert got == off["detected"].tolist()


def test_newly_detected_genes_enter_the_rerun(world: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """The result equals the M8-10 chain on the features with the gene table's calls (all 180
    genes) and differs from the chain on the features' own calls (120 genes)."""
    genes, feats = world
    cfg = sf.base_config()
    got = sf.det_result(genes, feats, cfg, _rng())
    want = sf.det_oracle(genes, feats, cfg, _rng())
    assert got is not None
    assert got.flagged == want.flagged
    assert set(got.dominant) == set(want.explained)
    assert sf.non_none(got.dominant) == want.dominant
    assert got.cv_r2 == pytest.approx(want.cv_r2, rel=1e-9, abs=1e-12)
    stale = genes.copy()
    stale["detected"] = feats["detected"].astype(bool).to_numpy()
    old = sf.det_oracle(stale, feats, cfg, _rng())
    assert set(got.dominant) != set(old.explained) or got.cv_r2 != pytest.approx(old.cv_r2)
    flipped = set(feats["locus_tag"].iloc[:N_FLIP].astype(str))
    assert set(got.dominant) & flipped, "no newly detected gene among the rerun's explained genes"


def test_undetected_in_genes_leaves_the_model_rows() -> None:
    """The other direction: a gene detected in the features but not in the gene table is not a
    model row (here: too few rows left, so no result)."""
    genes, feats = sf.planted_world()
    off = genes.copy()
    off["detected"] = [i < 90 for i in range(len(off))]
    assert sf.det_result(off, feats, sf.base_config(), _rng()) is None
