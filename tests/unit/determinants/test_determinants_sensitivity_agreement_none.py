"""M8-16 review round, lead decision 5 (statistician B1, code review N2): ``agreement`` leaves out
the genes whose dominant block is ``None`` in both runs, is NA when no gene is left, and counts
``None`` in only one run as a disagreement; a run compared with itself gives 1
(docs/specs/M8-16-sensitivity-function-determinants.md "Lead decision 5"; GIPPA_SPEC.md §8.6
"plus ``agreement`` for the dominant determinant block").

This reverses decision 4's "a block of None in both runs does not count as agreeing". The
identity case on a real rerun follows the statistician's attack: a planted world with imputed
detected genes, whose ``missing=min`` rerun equals the primary (imputed values never enter the
determinants models), must read ``jaccard 1, agreement 1``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import determinants_models_data as dm
from helpers import sensitivity_function_data as sf
from hypothesis import given, settings
from hypothesis import strategies as st


def _S() -> Any:
    from gippa.context import sensitivity as S

    return S


def _agreement(a: dict[str, str | None], b: dict[str, str | None]) -> Any:
    return _S().agreement(sf.det_sens([], a, 0.1), sf.det_sens([], b, 0.1))


def _row(a: dict[str, str | None], b: dict[str, str | None]) -> Any:
    S = _S()
    rows = S.sensitivity_rows({sf.DET: sf.det_sens(["x"], a, 0.1)},
                              {(sf.DET, "v"): sf.det_sens(["x"], b, 0.1)})  # fmt: skip
    return rows.iloc[0]


# --------------------------------------------------------------------------- examples


def test_run_against_itself_with_a_none_block_is_one() -> None:
    """B1: a ``None`` block in both runs is left out, so r vs r is 1 (not 1/2)."""
    r = {"g1": "codon", "g2": None}
    assert _agreement(r, r) == pytest.approx(1.0)
    assert float(_row(r, r)["agreement"]) == pytest.approx(1.0)


def test_only_none_blocks_shared_is_na() -> None:
    """No gene left after dropping the both-``None`` genes: agreement is NA (not 0)."""
    r = {"g1": None, "g2": None}
    assert _agreement(r, r) is None
    assert pd.isna(_row(r, r)["agreement"])
    assert pd.isna(_row({"g1": None, "g5": "codon"}, {"g1": None, "g6": "codon"})["agreement"])


@pytest.mark.parametrize(("a", "b"), [
    ({"g1": "codon", "g2": None}, {"g1": "codon", "g2": "protein"}),
    ({"g1": "codon", "g2": "protein"}, {"g1": "codon", "g2": None}),
])  # fmt: skip
def test_none_in_one_run_is_a_disagreement(a: dict[str, Any], b: dict[str, Any]) -> None:
    """``None`` in only one run stays in the denominator and does not agree: 1 of 2."""
    assert _agreement(a, b) == pytest.approx(0.5)
    assert float(_row(a, b)["agreement"]) == pytest.approx(0.5)


def test_mixed_none_cases() -> None:
    """g1 agrees, g2 is ``None`` in both (left out), g3 is ``None`` in one run (disagrees), g4
    disagrees: 1 of 3."""
    a = {"g1": "codon", "g2": None, "g3": None, "g4": "codon"}
    b = {"g1": "codon", "g2": None, "g3": "protein", "g4": "context"}
    assert _agreement(a, b) == pytest.approx(1 / 3)
    assert _agreement(b, a) == pytest.approx(1 / 3)


# --------------------------------------------------------------------------- properties

_BLOCKS = st.sampled_from(["transcription", "translation", "codon", "protein", "context"])
_DOM = st.dictionaries(st.sampled_from([f"g{i}" for i in range(12)]),
                       st.one_of(st.none(), _BLOCKS), max_size=10)  # fmt: skip


@settings(max_examples=100, deadline=None)
@given(r=_DOM)
def test_identity_property(r: dict[str, str | None]) -> None:
    """agreement(r, r) is 1 when any gene of r has a block, else NA."""
    got = _agreement(r, r)
    if any(v is not None for v in r.values()):
        assert got == pytest.approx(1.0)
    else:
        assert got is None


@settings(max_examples=100, deadline=None)
@given(p=_DOM, v=_DOM)
def test_agreement_property_with_none(p: dict[str, str | None], v: dict[str, str | None]) -> None:
    """agreement = |same non-None block| / |shared genes not None in both| (NA when that set is
    empty); symmetric; in [0, 1]; the row carries the same value."""
    kept = [g for g in set(p) & set(v) if not (p[g] is None and v[g] is None)]
    ab, ba = _agreement(p, v), _agreement(v, p)
    if not kept:
        assert ab is None and ba is None
        assert pd.isna(_row(p, v)["agreement"])
        return
    want = sum(p[g] is not None and p[g] == v[g] for g in kept) / len(kept)
    assert ab == pytest.approx(want)
    assert ba == pytest.approx(want)
    assert 0.0 <= ab <= 1.0
    assert float(_row(p, v)["agreement"]) == pytest.approx(want)


# --------------------------------------------------------------------------- the identity attack


def _imputed_world(n: int, seed: int, signal: float, n_imputed: int) -> tuple[Any, Any]:
    """``sf.planted_world`` with ``n_imputed`` detected genes imputed at the observed median (in
    the gene table and the features table)."""
    r = np.random.default_rng(seed)
    tu = dm.tu_ids(dm.mixed_sizes(r, n, 1, 3))
    tg = [f"g{i + 1}" for i in range(n)]
    cols = {c: r.normal(0.0, 1.0, n).tolist() for c in dm.PLANTED_BLOCKS}
    y = signal * np.asarray(cols["cai"]) + r.normal(0.0, 1.0, n)
    imp = np.zeros(n, dtype=bool)
    imp[r.choice(n, n_imputed, replace=False)] = True
    ab = np.where(imp, float(np.median(y[~imp])), y)
    feats = dm.features(cols, tu=tu, tags=tg, abundance=ab.tolist(), imputed=imp.tolist())
    genes = sf.gene_table(ab.tolist(), [True] * n, [0.5] * n, imputed=imp.tolist())
    genes["tu_id"] = tu
    return genes, feats


def test_missing_min_on_imputed_genes_reads_jaccard_1_agreement_1() -> None:
    """Statistician B1 attack (n = 400, signal 0.6, 40 detected genes imputed; seed 3 has a
    top-decile gene without a positive block): the ``missing=min`` rerun equals the primary,
    so its row reads jaccard 1 and agreement 1. Since M8-16c (|SHAP| rule) a None block needs an
    all-zero row, so the fixture only requires explained genes with a block."""
    from gippa.rng import for_step

    S = _S()
    cfg = sf.base_config()
    genes, feats = _imputed_world(400, 3, 0.6, 40)
    stream = f"sensitivity/{sf.DET}"
    p = sf.det_result(genes, feats, cfg, for_step(stream, sf.SEED))
    var = S.Variant("missing=min", {"abundance.missing": "min"})
    g_min = S.variant_genes(genes, var, cfg)
    v = sf.det_result(g_min, feats, cfg, for_step(stream, sf.SEED))
    assert p is not None and v is not None
    assert any(b is not None for b in p.dominant.values()), "fixture: no explained gene"
    assert p == v, "fixture: missing=min should reproduce the primary"
    r = S.sensitivity_rows({sf.DET: p}, {(sf.DET, "missing=min"): v}).iloc[0]
    assert float(r["jaccard"]) == 1.0
    assert float(r["agreement"]) == 1.0
