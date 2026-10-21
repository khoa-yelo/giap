"""M6-05 pass 1: discordant OP pairs — the t test and the BH family
(docs/specs/M6-05-discordance.md; GIAP_SPEC.md §7.4.2, §7.4.3, §8.1).

Integration pass (MVP-first): the both-detected branch, eligibility, the degenerate table and the
pass-1 one-detected branch (``p_miss = None`` -> ``no_detection_limit_model``). Per M6-index lead
decision 1, P(miss) never enters the BH family, so the family holds only t p-values. Example 6 and
the ``p_miss``-returns-NaN row are pass 2 (M6-05H) and are not here.

Inputs are M6-04 pair frames built with ``helpers.pairs_data`` on a linear chain of + strand
genes, all adjacent pairs OP, with the gene abundances chosen so the pair ``delta_log2`` values are
exactly the requested (dyadic-rounded) Δ.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.pairs_data import frame, gene, reps
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import giap.stats.tdist as tdist
from giap.config import GiapConfig
from giap.operons import discordance as disc
from giap.stats.multitest import bh
from giap.stats.tdist import StudentT, student_t_two_sided

ROW_COLS = ("pair_id", "test", "stat", "effect", "p_value", "q_value", "method", "flags")
_DYADIC = 2.0**20


def _cfg(**abundance: Any) -> GiapConfig:
    doc: dict[str, Any] = {"seed": 1}
    if abundance:
        doc["abundance"] = abundance
    return GiapConfig.model_validate(doc)


def _t_draws(n: int, seed: int = 7, loc: float = 0.05, scale: float = 0.4, df: float = 6.0):
    return loc + scale * np.random.default_rng(seed).standard_t(df, size=n)


def chain(
    deltas: Sequence[float],
    *,
    contig: str = "c1",
    prefix: str = "G",
    undetected: Sequence[int] = (),
    no_value: Sequence[int] = (),
    imputed: Sequence[int] = (),
    partial: bool = False,
    uniop_p: float = 0.9,
    abundance_type: str = "ibaq",
    abundance_scale: str = "molar",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pairs, genes) for a chain of len(deltas) + 1 genes; pair k (0-based) is
    (gene k, gene k + 1) with delta_log2 == deltas[k] rounded to 2**-20. ``undetected`` genes have
    no value; ``no_value`` genes are detected with NaN abundance; indices are gene indices."""
    d = np.round(np.asarray(deltas, dtype=float) * _DYADIC) / _DYADIC
    ab = np.concatenate([[20.0], 20.0 + np.cumsum(d)])
    tags = [f"{prefix}{i}" for i in range(len(ab))]
    genes = [gene(t, 100 + 400 * i, 400 + 400 * i, "+", contig) for i, t in enumerate(tags)]
    det = np.array([i not in undetected for i in range(len(ab))])
    val = np.where(det & ~np.isin(np.arange(len(ab)), list(no_value)), ab, np.nan)
    tbl = pd.DataFrame(
        {
            "locus_tag": tags,
            "detected": det,
            "p_detect": 0.9,
            "abundance_log2": val,
            "partial": partial,
            "abundance_imputed": np.isin(np.arange(len(ab)), list(imputed)),
            "abundance_type": abundance_type,
            "abundance_scale": abundance_scale,
        }
    )
    probs = {(tags[i], tags[i + 1]): uniop_p for i in range(len(tags) - 1)}
    pairs = frame(genes, reps((contig, 400 * len(ab) + 1000, "linear")), probs, tbl=tbl)
    if prefix != "G":
        pairs["pair_id"] = [f"{prefix}_{p}" for p in pairs["pair_id"]]
    return pairs, tbl


def run(pairs: pd.DataFrame, genes: pd.DataFrame, cfg: GiapConfig | None = None) -> Any:
    return disc.discordance(pairs, genes, cfg or _cfg(), None)


def rows_by_id(res: Any) -> pd.DataFrame:
    return res.rows.set_index("pair_id")


def flags_of(v: Any) -> list[str]:
    if isinstance(v, str):
        return [f for f in v.replace(",", ";").split(";") if f]
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return []
    return [str(f) for f in v]


@pytest.fixture
def fixed_fit(monkeypatch: pytest.MonkeyPatch) -> Iterator[Any]:
    """Replace the t fit by a fixed StudentT so the M6-01 constants are exact."""
    holder: dict[str, StudentT] = {}

    def fake(x: Any, *a: Any, **k: Any) -> StudentT:
        f = holder["fit"]
        n = int(np.isfinite(np.asarray(x, dtype=float)).sum())
        return StudentT(f.loc, f.scale, f.df, n, f.converged)

    monkeypatch.setattr(tdist, "fit_student_t", fake)
    if hasattr(disc, "fit_student_t"):
        monkeypatch.setattr(disc, "fit_student_t", fake)

    def set_fit(loc: float, scale: float, df: float, converged: bool = True) -> None:
        holder["fit"] = StudentT(loc, scale, df, 0, converged)

    yield set_fit


# ------------------------------------------------------------------ interface


def test_min_op_pairs_is_twenty() -> None:
    assert disc.MIN_OP_PAIRS == 20


def test_result_shape_and_row_columns() -> None:
    pairs, genes = chain(_t_draws(60))
    res = run(pairs, genes)
    assert isinstance(res, disc.Discordance)
    assert isinstance(res.fit, StudentT)
    assert set(ROW_COLS) <= set(res.rows.columns)
    assert set(pairs.loc[pairs["class"] == "OP", "pair_id"]) <= set(res.rows["pair_id"])
    assert res.rows["pair_id"].is_unique


def test_fit_rows_selects_usable_op_pairs() -> None:
    pairs, _ = chain([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], undetected=[6], no_value=[4])
    pairs = pairs.copy()
    pairs["abundance_imputed_a"] = [False, True, False, False, False, False]
    pairs["abundance_imputed_b"] = False
    pairs.loc[0, "class"] = "NO"
    got = np.sort(np.asarray(disc.fit_rows(pairs), dtype=float))
    # pair 0 not OP, pair 1 imputed, pairs 3/4 touch the NaN gene 4, pair 5 the undetected gene 6
    np.testing.assert_allclose(got, [0.3])


# ------------------------------------------------------------------ examples 1-3, 5


def test_example1_positive_delta_constant(fixed_fit: Any) -> None:
    fixed_fit(0.05, 0.4, 6.0)
    pairs, genes = chain([*_t_draws(200), 1.2])
    r = rows_by_id(run(pairs, genes)).loc["PAIR_201"]
    assert r["test"] == "both_detected_t"
    assert r["method"] == "student_t"
    assert float(r["p_value"]) == pytest.approx(0.0282, abs=5e-5)
    assert float(r["stat"]) == pytest.approx(1.2, abs=1e-6)
    assert float(r["effect"]) == pytest.approx(1.2 - 0.05, abs=1e-6)


def test_example2_negative_delta_constant(fixed_fit: Any) -> None:
    fixed_fit(0.05, 0.4, 6.0)
    pairs, genes = chain([*_t_draws(200), -1.2])
    r = rows_by_id(run(pairs, genes)).loc["PAIR_201"]
    assert float(r["p_value"]) == pytest.approx(0.0205, abs=5e-5)


def test_example3_delta_at_loc_gives_p_one_and_zero_effect(fixed_fit: Any) -> None:
    fixed_fit(0.25, 0.4, 6.0)
    pairs, genes = chain([*_t_draws(200), 0.25])
    r = rows_by_id(run(pairs, genes)).loc["PAIR_201"]
    assert float(r["p_value"]) == 1.0
    assert float(r["effect"]) == 0.0


def test_both_detected_rows_use_the_real_fit() -> None:
    pairs, genes = chain([*_t_draws(200), 1.2])
    res = run(pairs, genes)
    assert res.fit.converged and res.fit.n == 201
    ref = tdist.fit_student_t(pairs["delta_log2"].to_numpy(float))
    assert res.fit.loc == pytest.approx(ref.loc, abs=1e-3)
    assert res.fit.scale == pytest.approx(ref.scale, abs=1e-3)
    rows = rows_by_id(res)
    d = pairs.set_index("pair_id")["delta_log2"]
    np.testing.assert_allclose(
        rows.loc[d.index, "p_value"].to_numpy(float),
        student_t_two_sided(d.to_numpy(float), res.fit),
        rtol=1e-12,
    )
    np.testing.assert_allclose(rows.loc[d.index, "stat"].to_numpy(float), d.to_numpy(float))
    np.testing.assert_allclose(
        rows.loc[d.index, "effect"].to_numpy(float), d.to_numpy(float) - res.fit.loc
    )
    assert (rows["test"] == "both_detected_t").all()
    assert (rows["method"] == "student_t").all()
    if "null_model" in rows.columns:
        assert (rows["null_model"] == "within_operon_t").all()


def test_ci_is_na_and_flagged_unestimated() -> None:
    """Lead decision 2: no CI exists for a single pair's effect."""
    res = run(*chain(_t_draws(60)))
    tested = res.rows[res.rows["test"] == "both_detected_t"]
    assert len(tested) == 60
    for col in ("ci_low", "ci_high"):
        if col in tested.columns:
            assert tested[col].isna().all()
    assert all("unestimated" in flags_of(f) for f in tested["flags"])


def test_example5_imputed_pair_not_in_fit_but_tested() -> None:
    pairs, genes = chain([*_t_draws(200), 0.7][::-1], imputed=[0])
    res = run(pairs, genes)
    assert res.fit.n == 200
    r = rows_by_id(res).loc["PAIR_1"]
    assert r["test"] == "both_detected_t"
    assert np.isfinite(float(r["p_value"]))
    assert float(r["p_value"]) == pytest.approx(
        float(student_t_two_sided([0.7], res.fit)[0]), abs=1e-6
    )


def test_rna_between_pair_is_tested_normally() -> None:
    pairs, genes = chain(_t_draws(60))
    pairs = pairs.copy()
    pairs.at[3, "flags"] = ["rna_between"]
    r = rows_by_id(run(pairs, genes)).loc[pairs.at[3, "pair_id"]]
    assert r["test"] == "both_detected_t"
    assert np.isfinite(float(r["p_value"]))


# ------------------------------------------------------------------ examples 4, 7-10 and degenerate


def _all_untested(res: Any, flag: str) -> None:
    assert (res.rows["test"] == "none").all()
    assert res.rows["p_value"].isna().all()
    assert res.rows["q_value"].isna().all()
    assert all(flag in flags_of(f) for f in res.rows["flags"])


def test_example4_nineteen_usable_pairs_is_too_few() -> None:
    res = run(*chain(_t_draws(19)))
    _all_untested(res, "too_few_op_pairs")


def test_twenty_usable_pairs_is_enough() -> None:
    pairs, genes = chain(_t_draws(20, seed=3))
    assert tdist.fit_student_t(pairs["delta_log2"].to_numpy(float)).converged
    res = run(pairs, genes)
    assert (res.rows["test"] == "both_detected_t").all()
    assert res.rows["q_value"].notna().all()


def test_too_few_counts_usable_pairs_not_op_rows() -> None:
    """25 OP pairs, but the imputed genes leave only 19 usable for the fit."""
    res = run(*chain(_t_draws(25), imputed=[0, 2, 4]))
    _all_untested(res, "too_few_op_pairs")


def test_all_identical_deltas_degenerate_fit() -> None:
    res = run(*chain([0.5] * 60))
    assert not res.fit.converged
    _all_untested(res, "too_few_op_pairs")


def test_non_converged_fit_untested_and_warns(
    fixed_fit: Any, caplog: pytest.LogCaptureFixture
) -> None:
    fixed_fit(0.05, 0.4, 6.0, converged=False)
    with caplog.at_level(logging.WARNING):
        res = run(*chain(_t_draws(60)))
    _all_untested(res, "too_few_op_pairs")
    warn = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert warn and any("60" in r.getMessage() for r in warn)


def test_no_op_pair_at_all_is_an_empty_family() -> None:
    res = run(*chain(_t_draws(30), uniop_p=0.1))
    assert (res.rows["test"] == "none").all()
    assert res.rows["q_value"].isna().all()


def test_example7_one_detected_pass1_has_no_detection_limit_model() -> None:
    pairs, genes = chain(_t_draws(60), undetected=[60])
    r = rows_by_id(run(pairs, genes)).loc["PAIR_60"]
    assert r["test"] == "none"
    assert pd.isna(r["p_value"]) and pd.isna(r["q_value"])
    assert "no_detection_limit_model" in flags_of(r["flags"])


def test_neither_detected_and_excluded_missing_are_untested() -> None:
    # genes 60, 61 undetected -> PAIR_61 neither detected; gene 30 detected without a value
    pairs, genes = chain(_t_draws(61), undetected=[60, 61], no_value=[30])
    rows = rows_by_id(run(pairs, genes, _cfg(missing="exclude")))
    for pid in ("PAIR_61", "PAIR_30", "PAIR_31"):
        assert rows.loc[pid, "test"] == "none"
        assert pd.isna(rows.loc[pid, "p_value"])
    assert rows.loc["PAIR_10", "test"] == "both_detected_t"


def test_example8_one_family_of_t_p_values_only() -> None:
    """100 both-detected rows and 50 one-detected rows: P(miss) stays out (M6-index decision 1)."""
    p1, g1 = chain(_t_draws(100, seed=11))
    p2, g2 = chain(_t_draws(50, seed=12), contig="c2", prefix="H", undetected=range(1, 51, 2))
    rows = run(pd.concat([p1, p2], ignore_index=True), pd.concat([g1, g2], ignore_index=True)).rows
    one = rows["pair_id"].isin(p2["pair_id"])
    assert one.sum() == 50
    assert (rows.loc[one, "test"] == "none").all()
    assert rows.loc[one, "q_value"].isna().all()
    fam = rows["p_value"].notna()
    assert fam.sum() == 100 and (rows.loc[fam, "test"] == "both_detected_t").all()
    np.testing.assert_allclose(
        rows.loc[fam, "q_value"].to_numpy(float), bh(rows.loc[fam, "p_value"].to_numpy(float))
    )
    s = rows[fam].sort_values("p_value")
    assert (np.diff(s["q_value"].to_numpy(float)) >= -1e-15).all()


def test_example9_relative_scale_skips_every_test() -> None:
    pairs, genes = chain(_t_draws(60), abundance_type="maxlfq", abundance_scale="relative")
    before = pairs.copy(deep=True)
    res = run(pairs, genes)
    _all_untested(res, "relative_scale")
    pd.testing.assert_frame_equal(pairs, before)  # classes, distances, evidence untouched
    assert (pairs["class"] == "OP").all()


def test_example10_deterministic() -> None:
    pairs, genes = chain(_t_draws(80))
    a, b = run(pairs, genes), run(pairs.copy(), genes.copy())
    pd.testing.assert_frame_equal(a.rows, b.rows)
    assert a.fit == b.fit


def test_input_frames_not_mutated() -> None:
    pairs, genes = chain(_t_draws(40))
    p0, g0 = pairs.copy(deep=True), genes.copy(deep=True)
    run(pairs, genes)
    pd.testing.assert_frame_equal(pairs, p0)
    pd.testing.assert_frame_equal(genes, g0)


# ------------------------------------------------------------------ properties

_PROP = settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])


@_PROP
@given(
    n=st.integers(5, 45),
    seed=st.integers(0, 2**16),
    und=st.sets(st.integers(0, 45), max_size=8),
)
def test_bh_family_is_exactly_the_tested_rows(n: int, seed: int, und: set[int]) -> None:
    rows = run(*chain(_t_draws(n, seed=seed), undetected=sorted(i for i in und if i <= n))).rows
    tested = rows["test"] != "none"
    assert (rows["p_value"].notna() == tested).all()
    assert (rows["q_value"].notna() == tested).all()
    if tested.any():
        np.testing.assert_allclose(
            rows.loc[tested, "q_value"].to_numpy(float),
            bh(rows.loc[tested, "p_value"].to_numpy(float)),
        )


@_PROP
@given(seed=st.integers(0, 2**16), perm_seed=st.integers(0, 2**16))
def test_fit_invariant_to_row_order_and_pair_id(seed: int, perm_seed: int) -> None:
    pairs, genes = chain(_t_draws(40, seed=seed))
    base = run(pairs, genes)
    order = np.random.default_rng(perm_seed).permutation(len(pairs))
    shuf = pairs.iloc[order].reset_index(drop=True)
    relabel = {p: f"X{k}" for k, p in enumerate(shuf["pair_id"])}
    shuf["pair_id"] = shuf["pair_id"].map(relabel)
    got = run(shuf, genes)
    assert got.fit.converged == base.fit.converged and got.fit.n == base.fit.n
    for attr in ("loc", "scale", "df"):
        assert getattr(got.fit, attr) == pytest.approx(
            getattr(base.fit, attr), rel=1e-9, nan_ok=True
        )
    a = base.rows.set_index("pair_id")[["test", "p_value"]]
    b = got.rows.set_index("pair_id")[["test", "p_value"]].rename(
        index={v: k for k, v in relabel.items()}
    )
    pd.testing.assert_frame_equal(b.loc[a.index], a, rtol=1e-9)


@_PROP
@given(seed=st.integers(0, 2**16), n_partial=st.integers(1, 30))
def test_partial_gene_pairs_change_nothing(seed: int, n_partial: int) -> None:
    pairs, genes = chain(_t_draws(40, seed=seed))
    base = run(pairs, genes)
    p2, g2 = chain(_t_draws(n_partial, seed=seed + 1) * 5, contig="c2", prefix="Q", partial=True)
    got = run(pd.concat([pairs, p2], ignore_index=True), pd.concat([genes, g2], ignore_index=True))
    assert got.fit == base.fit
    rows = got.rows.set_index("pair_id")
    pd.testing.assert_frame_equal(
        rows.loc[base.rows["pair_id"], [c for c in ROW_COLS if c != "pair_id"]],
        base.rows.set_index("pair_id")[[c for c in ROW_COLS if c != "pair_id"]],
    )
    extra = rows.loc[p2["pair_id"]]
    assert (extra["test"] == "none").all() and extra["q_value"].isna().all()
