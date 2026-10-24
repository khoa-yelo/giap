"""M8-16c: ``agreement`` over the genes the features explain (Khoa's decision of 2026-10-24,
follow-up to M8-16 #265; docs/specs/M8-16-sensitivity-function-determinants.md lead decisions
4-5; GIAP_SPEC.md §8.6 "plus ``agreement`` for the dominant determinant block").

New definition (lead):

* an **explained gene** is a model row not flagged ``abundance_not_explained_by_sequence`` in
  that run;
* the sensitivity rerun's ``DeterminantsSensitivity.dominant`` gives a block for **every
  explained model row** (not only the top decile) from that row's LightGBM SHAP block sums;
* that block is the **agreement block** (Khoa, option (a) of the #299 statistician review,
  2026-10-24): the block with the largest **|SHAP|** sum, not M8-08 rule 3's largest positive
  sum, so a gene predicted below the mean keeps its informative (negative) block instead of a
  noise block just above 0. Interface: ``giap.determinants.sensitivity.agreement_block(sums)``.
  Ties and degenerate rows mirror ``_dominant``: NaN ignored, ties go to the first block in
  ``BLOCKS`` order (whatever the signs), ``None`` when every |sum| is 0 or NaN;
* ``agreement`` is the share of the genes explained in both runs whose dominant block is the
  same: both-``None`` genes left out, ``None`` in one run a disagreement, NA when no gene is
  left, a run against itself 1;
* ``determinants_attribution.tsv`` is unchanged (``dominant_block`` still rule 3, top-decile
  only); the sensitivity.tsv columns are unchanged.

The oracles here (``sf.abs_block`` / ``sf.positive_block``) are written out, not ``_dominant``
(statistician finding 6). Calibration of the fix: tests/calibration/
test_determinants_agreement_calibration.py.

Interpretation (test author, reported to the lead): ``agreement`` excludes a gene that is in a
run's ``flagged`` set even if a hand-built ``dominant`` lists it, so "explained in both runs" is
read from ``flagged`` and not only from the ``dominant`` keys.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import sensitivity_function_data as sf
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.determinants.attribution import SHAP_COLUMNS


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def _rng() -> np.random.Generator:
    from giap.rng import for_step

    return for_step(f"sensitivity/{sf.DET}", sf.SEED)


def _agree(a: Any, b: Any) -> Any:
    return _S().agreement(a, b)


def _row_agreement(a: Any, b: Any) -> Any:
    rows = _S().sensitivity_rows({sf.DET: a}, {(sf.DET, "v"): b})
    return rows.iloc[0]["agreement"]


def _want(a: Any, b: Any) -> float | None:
    """The new definition, written out."""
    ea = {g: v for g, v in a.dominant.items() if g not in a.flagged}
    eb = {g: v for g, v in b.dominant.items() if g not in b.flagged}
    kept = [g for g in ea.keys() & eb.keys() if not (ea[g] is None and eb[g] is None)]
    if not kept:
        return None
    return sum(ea[g] is not None and ea[g] == eb[g] for g in kept) / len(kept)


# --------------------------------------------------------------------------- hand-built


def test_unexplained_gene_is_excluded_from_agreement() -> None:
    """g2 is flagged in the variant: it is not explained in both runs, so it leaves the
    denominator even though both ``dominant`` maps list it (with different blocks). g1 agrees,
    g3 disagrees: 1 of 2 (the old top-decile rule read 1 of 3)."""
    p = sf.det_sens([], {"g1": "codon", "g2": "codon", "g3": "protein"}, 0.4)
    v = sf.det_sens(["g2"], {"g1": "codon", "g2": "transcription", "g3": "context"}, 0.4)
    assert _agree(p, v) == pytest.approx(0.5)
    assert _agree(v, p) == pytest.approx(0.5)
    assert float(_row_agreement(p, v)) == pytest.approx(0.5)


def test_gene_unexplained_in_both_runs_is_excluded() -> None:
    """g9 is flagged in both runs: left out. Only g1 is explained in both (and agrees): 1."""
    p = sf.det_sens(["g9"], {"g1": "codon", "g9": "codon"}, 0.4)
    v = sf.det_sens(["g9"], {"g1": "codon", "g9": "protein"}, 0.4)
    assert _agree(p, v) == pytest.approx(1.0)


def test_only_unexplained_genes_shared_is_na() -> None:
    """Every shared gene is flagged in one run or the other: no gene left, NA."""
    p = sf.det_sens(["g1"], {"g1": "codon", "g2": "codon"}, 0.4)
    v = sf.det_sens(["g2"], {"g1": "codon", "g2": "codon"}, 0.4)
    assert _agree(p, v) is None
    assert pd.isna(_row_agreement(p, v))


def test_explained_non_top_decile_genes_count() -> None:
    """Spec-style example: three top-decile genes (g1-g3) plus five explained genes outside the
    top decile (g10-g14), all in ``dominant`` under the new rule; g4 is flagged in the primary.
    Same block on g1, g2, g10, g11, g12; different on g3, g13; g14 is ``None`` in the variant
    only (a disagreement): 5 of 8."""
    p = sf.det_sens(["g4"], {"g1": "codon", "g2": "codon", "g3": "transcription",
                             "g10": "codon", "g11": "protein", "g12": "context",
                             "g13": "codon", "g14": "codon"}, 0.41)  # fmt: skip
    v = sf.det_sens([], {"g1": "codon", "g2": "codon", "g3": "protein", "g4": "codon",
                         "g10": "codon", "g11": "protein", "g12": "context",
                         "g13": "translation", "g14": None}, 0.38)  # fmt: skip
    assert _agree(p, v) == pytest.approx(5 / 8)
    assert float(_row_agreement(p, v)) == pytest.approx(5 / 8)


def test_none_handling_on_explained_genes() -> None:
    """``None`` in both runs left out (g2); ``None`` in one run a disagreement (g3); a flagged
    gene with ``None`` in the other run is not a disagreement but excluded (g4): 1 of 2."""
    p = sf.det_sens(["g4"], {"g1": "codon", "g2": None, "g3": None, "g4": "codon"}, 0.3)
    v = sf.det_sens([], {"g1": "codon", "g2": None, "g3": "protein", "g4": None}, 0.3)
    assert _agree(p, v) == pytest.approx(0.5)
    assert _agree(v, p) == pytest.approx(0.5)


def test_run_against_itself_with_flagged_genes_is_one() -> None:
    """A run with flagged genes (listed in its hand-built ``dominant``) against itself: 1."""
    r = sf.det_sens(["g3"], {"g1": "codon", "g2": None, "g3": "protein"}, 0.3)
    assert _agree(r, r) == pytest.approx(1.0)


def test_self_agreement_na_when_only_flagged_genes_have_a_block() -> None:
    """The explained genes all have ``None``, the only block sits on a flagged gene: no gene is
    left, so a run against itself is NA (the old rule read 1)."""
    r = sf.det_sens(["g3"], {"g1": None, "g2": None, "g3": "protein"}, 0.3)
    assert _agree(r, r) is None


# --------------------------------------------------------------------------- the block rule


def _block(sums: list[float]) -> Any:
    from giap.determinants.sensitivity import agreement_block

    return agreement_block(np.asarray(sums, dtype=float))


def _B() -> tuple[str, ...]:
    from giap.determinants.attribution import BLOCKS

    return BLOCKS


def test_negative_dominant_row_picks_the_negative_block() -> None:
    """A gene predicted below the mean: block 1 carries -2.0, the others noise around 0. The
    agreement block is block 1; rule 3 (largest positive) would pick block 2 (0.3)."""
    from giap.determinants.attribution import _dominant

    B = _B()
    row = [0.1, -2.0, 0.3, -0.05, 0.0]
    assert _block(row) == B[1]
    assert _dominant(np.asarray(row)) == B[2], "rule 3 itself must stay as it is"


def test_all_negative_row_has_a_block() -> None:
    """Every sum negative: rule 3 gives ``None``, the agreement block the largest |sum|."""
    B = _B()
    assert _block([-0.2, -0.1, -0.9, -0.3, -0.4]) == B[2]


def test_ties_mirror_dominant() -> None:
    """Ties at the largest |sum| go to the first block in ``BLOCKS`` order, as ``_dominant``
    breaks a positive tie; a +x / -x tie the same way."""
    from giap.determinants.attribution import _dominant

    B = _B()
    assert _dominant(np.asarray([0.0, 1.0, 1.0, 0.0, 0.0])) == B[1]
    assert _block([0.0, 1.0, 1.0, 0.0, 0.0]) == B[1]
    assert _block([0.0, -1.0, 1.0, 0.0, 0.0]) == B[1]
    assert _block([0.0, 1.0, -1.0, 0.0, 0.0]) == B[1]
    assert _block([0.5, 0.0, 0.0, -1.0, -1.0]) == B[3]


@pytest.mark.parametrize(
    "row",
    [[0.0] * 5, [-0.0] * 5, [math.nan] * 5, [math.nan, 0.0, 0.0, math.nan, 0.0]],
    ids=["zeros", "neg-zeros", "all-nan", "nan-and-zeros"],
)
def test_all_zero_or_nan_row_is_none(row: list[float]) -> None:
    """No block carries any SHAP mass: ``None`` (as ``_dominant`` on an all-zero row)."""
    from giap.determinants.attribution import _dominant

    assert _block(row) is None
    assert _dominant(np.asarray(row)) is None


def test_nan_is_ignored() -> None:
    """A NaN block (absent from the model matrix) never wins."""
    B = _B()
    assert _block([math.nan, -0.5, 0.2, math.nan, 0.1]) == B[1]


_SUM = st.one_of(
    st.floats(min_value=-5, max_value=5, allow_nan=False),
    st.just(0.0),
    st.just(math.nan),
    st.sampled_from([-1.0, 1.0]),
)


@settings(max_examples=200, deadline=None)
@given(row=st.lists(_SUM, min_size=5, max_size=5))
def test_agreement_block_property(row: list[float]) -> None:
    """Property: the agreement block equals the written-out |SHAP| rule; on a row without a
    negative sum it equals rule 3; it is ``None`` exactly when no |sum| is > 0."""
    from giap.determinants.attribution import _dominant

    got = _block(row)
    assert got == sf.abs_block(row)
    finite = [x for x in row if not math.isnan(x)]
    assert (got is None) == (not any(abs(x) > 0 for x in finite))
    if all(x >= 0 for x in finite):
        assert got == _dominant(np.asarray(row, dtype=float))


# --------------------------------------------------------------------------- properties

_BLOCKS = st.sampled_from(["transcription", "translation", "codon", "protein", "context"])
_GENES = [f"g{i}" for i in range(12)]
_DOM = st.dictionaries(st.sampled_from(_GENES), st.one_of(st.none(), _BLOCKS), max_size=10)
_FLAG = st.frozensets(st.sampled_from(_GENES), max_size=5)


@settings(max_examples=150, deadline=None)
@given(d=_DOM, fl=_FLAG)
def test_self_agreement_property(d: dict[str, str | None], fl: frozenset[str]) -> None:
    """Property: a run against itself is 1 when any explained gene has a block, else NA."""
    r = sf.det_sens(sorted(fl), d, 0.2)
    got = _agree(r, r)
    if any(v is not None for g, v in d.items() if g not in fl):
        assert got == pytest.approx(1.0)
    else:
        assert got is None


@settings(max_examples=150, deadline=None)
@given(p=_DOM, fp=_FLAG, v=_DOM, fv=_FLAG)
def test_agreement_matches_the_explained_definition(
    p: dict[str, str | None], fp: frozenset[str], v: dict[str, str | None], fv: frozenset[str]
) -> None:
    """Property: agreement equals the written-out definition, symmetric, in [0, 1], and the
    sensitivity row carries the same value."""
    a, b = sf.det_sens(sorted(fp), p, 0.2), sf.det_sens(sorted(fv), v, 0.2)
    want = _want(a, b)
    ab, ba = _agree(a, b), _agree(b, a)
    if want is None:
        assert ab is None and ba is None
        assert pd.isna(_row_agreement(a, b))
        return
    assert ab == pytest.approx(want)
    assert ba == pytest.approx(want)
    assert 0.0 <= ab <= 1.0
    assert float(_row_agreement(a, b)) == pytest.approx(want)


# --------------------------------------------------------------------------- the rerun


def _wide_oracle(
    genes: pd.DataFrame, features: pd.DataFrame, cfg: Any, rng: np.random.Generator
) -> tuple[dict[str, str | None], pd.DataFrame]:
    """The M8-10 chain without its null (as ``sf.det_oracle``) and the written-out |SHAP| rule
    (``sf.abs_block``) on the SHAP block sums of every explained model row; also the attribution
    table."""
    from giap.determinants.attribution import (
        SHAP_COLUMNS,
        UNEXPLAINED,
        attribution_table,
        shap_blocks,
    )
    from giap.determinants.features import model_matrix
    from giap.determinants.models import cross_validate, model_rows, tu_folds

    f = sf.merged_features(genes, features)
    rows = model_rows(f)
    X_all, blocks = model_matrix(f)
    X = X_all[rows]
    y = pd.to_numeric(f.loc[rows, "abundance_log2"]).to_numpy(float)
    tu = f.loc[rows, "tu_id"].astype(str).tolist()
    d = cfg.determinants
    folds = tu_folds(tu, d.outer_folds, d.repeats, rng)
    cv = cross_validate(X, y, folds, blocks, d, cfg.seed, 1, tu=tu)
    shap = shap_blocks(cv["lightgbm"], X, blocks)
    shap.index = pd.Index(f.loc[shap.index, "locus_tag"].astype(str))
    att = attribution_table(f, rows, cv["lightgbm"], shap)
    model_tags = set(f.loc[rows, "locus_tag"].astype(str))
    dom: dict[str, str | None] = {}
    for t, fl in zip(att["locus_tag"].astype(str), att["flags"], strict=True):
        if t in model_tags and UNEXPLAINED not in str(fl).split(";"):
            dom[t] = sf.abs_block(np.asarray(shap.loc[t, list(SHAP_COLUMNS)], dtype=float))
    return dom, att


_SHAP = SHAP_COLUMNS


@pytest.fixture(scope="module")
def world() -> tuple[pd.DataFrame, pd.DataFrame]:
    return sf.planted_world()


def test_dominant_covers_every_explained_model_row(
    world: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """``dominant`` has one entry per explained model row (top decile or not), none for a
    flagged gene, each the |SHAP| block of that row's SHAP sums (not rule 3's: the two differ on
    many genes predicted below the mean)."""
    genes, feats = world
    cfg = sf.base_config()
    res = sf.det_result(genes, feats, cfg, _rng())
    want, att = _wide_oracle(genes, feats, cfg, _rng())
    assert res is not None
    top = set(att.loc[att["top_decile"].astype(bool), "locus_tag"].astype(str))
    assert res.flagged, "fixture: no flagged gene"
    assert len(set(want) - top) > 50, "fixture: few explained genes outside the top decile"
    assert set(res.dominant) == set(want)
    assert not set(res.dominant) & set(res.flagged)
    assert set(res.dominant) - top, "dominant still covers the top decile only"
    rule3 = {t: sf.positive_block(att.set_index("locus_tag").loc[t, list(_SHAP)].to_numpy(float))
             for t in want}  # fmt: skip
    assert sum(want[t] != rule3[t] for t in want) > 20, "fixture: rules agree on most genes"
    assert dict(res.dominant) == want


def test_top_decile_blocks_follow_the_agreement_rule(
    world: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """On explained top-decile rows the wider map gives the |SHAP| block of the row's SHAP
    columns in the attribution table (the same sums the table's rule-3 ``dominant_block`` is read
    from); the table still leaves ``dominant_block`` NA off the top decile."""
    genes, feats = world
    cfg = sf.base_config()
    res = sf.det_result(genes, feats, cfg, _rng())
    _, att = _wide_oracle(genes, feats, cfg, _rng())
    assert res is not None
    is_top = att["top_decile"].astype(bool)
    assert att.loc[~is_top, "dominant_block"].isna().all()
    top = att[is_top]
    for t, row in zip(top["locus_tag"].astype(str), top[list(_SHAP)].to_numpy(float),
                      strict=True):  # fmt: skip
        if t in res.flagged:
            continue
        assert res.dominant[t] == sf.abs_block(row)


def test_attribution_dominant_block_keeps_rule_3() -> None:
    """``determinants_attribution.tsv``'s ``dominant_block`` is unchanged: rule 3 (largest
    positive SHAP sum) on top-decile rows, NA elsewhere, even where the |SHAP| block differs
    (n = 400, signal 0.35: about 19 top-decile rows). The rerun's ``dominant`` takes the
    |SHAP| block on those same rows."""
    genes, feats = sf.planted_world(n=400, signal=0.35)
    cfg = sf.base_config()
    res = sf.det_result(genes, feats, cfg, _rng())
    _, att = _wide_oracle(genes, feats, cfg, _rng())
    assert res is not None
    is_top = att["top_decile"].astype(bool)
    assert att.loc[~is_top, "dominant_block"].isna().all()
    top = att[is_top]
    differ = 0
    for t, b, row in zip(top["locus_tag"].astype(str), top["dominant_block"],
                         top[list(_SHAP)].to_numpy(float), strict=True):  # fmt: skip
        assert (None if pd.isna(b) else str(b)) == sf.positive_block(row)
        if sf.abs_block(row) != sf.positive_block(row):
            differ += 1
            if t not in res.flagged:
                assert res.dominant[t] == sf.abs_block(row)
    assert differ >= 5, "fixture: the rules agree on the top decile"


def test_variant_agreement_is_over_genes_explained_in_both(
    world: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    """A noisy variant (top3-like change, 15 genes losing their abundance): the row's agreement
    equals the definition on the two wider maps, and its denominator is far larger than the
    two top deciles could give."""
    genes, feats = world
    g = genes.copy()
    r = np.random.default_rng(11)
    g["abundance_log2"] = g["abundance_log2"].to_numpy(float) + r.normal(0, 0.7, len(g))
    g.loc[g.index[:15], "abundance_log2"] = math.nan
    cfg = sf.base_config()
    p = sf.det_result(genes, feats, cfg, _rng())
    v = sf.det_result(g, feats, cfg, _rng())
    assert p is not None and v is not None
    wp, _ = _wide_oracle(genes, feats, cfg, _rng())
    wv, _ = _wide_oracle(g, feats, cfg, _rng())
    both = [t for t in wp.keys() & wv.keys() if not (wp[t] is None and wv[t] is None)]
    assert len(both) > 50, "fixture: few genes explained in both runs"
    want = sum(wp[t] is not None and wp[t] == wv[t] for t in both) / len(both)
    assert _agree(p, v) == pytest.approx(want)
    assert float(_row_agreement(p, v)) == pytest.approx(want)


def test_rerun_against_itself_agrees_fully(world: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """The same rerun twice (same seed): agreement 1 over the wider map."""
    genes, feats = world
    cfg = sf.base_config()
    a = sf.det_result(genes, feats, cfg, _rng())
    b = sf.det_result(genes, feats, cfg, _rng())
    assert a is not None and b is not None
    assert len(a.dominant) > 50
    assert _agree(a, b) == pytest.approx(1.0)


# --------------------------------------------------------------------------- n_agreement


def _n_want(a: Any, b: Any) -> int:
    ea = {g: v for g, v in a.dominant.items() if g not in a.flagged}
    eb = {g: v for g, v in b.dominant.items() if g not in b.flagged}
    return sum(not (ea[g] is None and eb[g] is None) for g in ea.keys() & eb.keys())


_BLK = st.sampled_from(["codon", "protein", "transcription", None])


@settings(max_examples=80, deadline=None)
@given(
    p=st.dictionaries(st.sampled_from([f"g{i}" for i in range(10)]), _BLK, max_size=10),
    v=st.dictionaries(st.sampled_from([f"g{i}" for i in range(10)]), _BLK, max_size=10),
    fp=st.sets(st.sampled_from([f"g{i}" for i in range(10)]), max_size=3),
    fv=st.sets(st.sampled_from([f"g{i}" for i in range(10)]), max_size=3),
)
def test_n_agreement_counts_the_kept_genes(
    p: dict[str, Any], v: dict[str, Any], fp: set[str], fv: set[str]
) -> None:
    """#299 code review 1: ``n_agreement`` is the denominator of ``agreement``; 0 exactly when
    ``agreement`` is NA."""
    a, b = sf.det_sens(sorted(fp), p, 0.3), sf.det_sens(sorted(fv), v, 0.3)
    n = _S().n_agreement(a, b)
    assert n == _n_want(a, b) == _S().n_agreement(b, a)
    assert (n == 0) == (_want(a, b) is None)
