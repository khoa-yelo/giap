"""M6-13: the abundance-dependent scale in the discordance step (ADR-0055 option (a);
docs/specs/M6-13-discordance-abundance-sigma.md; GIAP_SPEC.md §7.4.2, §7.4.3, §8.1).

Interface pinned here:

- ``discordance(..., *, sigma_model="loglinear")``; ``"constant"`` reproduces current main bit for
  bit, with ``scaled.reason = "disabled"``;
- ``Discordance.scaled: ScaledStudentT`` is the null of every both-detected p; ``Discordance.fit``
  stays the constant ``StudentT``; ``Discordance.ratios`` (= ``fit_rows``) is unchanged;
- ``fit_pairs(pairs) -> (Delta, m)`` of the usable OP pairs, sorted lexicographically;
- an imputed tested pair (never a fit pair) takes z from its own m, clamped to the fit range
  (S-2, S-3).

Examples E5 (through the step), E6, E7, E12, E13 and E14, the ``scaled`` field, the fallback
WARNING and the bit-identity of ``sigma_model = "constant"`` on the M6-05 fixtures. The frames are
hand-built pair tables (the columns ``discordance`` reads), so Delta and the pair mean m can be
set independently; M6-04 frames from ``test_discordance.chain`` cover the bit-identity check.
"""

from __future__ import annotations

import inspect
import logging
import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.optimize
from scipy import stats as sps
from test_discordance import _t_draws
from test_discordance import chain as m605_chain

import giap.stats.tdist as tdist
from giap.config import GiapConfig
from giap.operons import discordance as disc
from giap.stats.multitest import bh
from giap.stats.tdist import ScaledStudentT, StudentT, fit_student_t, student_t_two_sided

_REAL_MINIMIZE = scipy.optimize.minimize


def _cfg() -> GiapConfig:
    return GiapConfig(seed=1)


def frame(
    delta: Sequence[float],
    m: Sequence[float],
    *,
    imputed: Sequence[int] = (),
    one_detected: Sequence[tuple[float, float]] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pairs, genes): OP pair i is (A{i}, B{i}) with a = m - Delta / 2 and b = m + Delta / 2.
    ``imputed`` pair indices get an imputed gene_b; each ``one_detected`` (a_v, _) adds an OP
    pair with gene_a detected at a_v and gene_b undetected."""
    d = np.asarray(delta, dtype=float)
    mm = np.asarray(m, dtype=float)
    a, b = mm - d / 2, mm + d / 2
    n = d.size
    rows = {
        "pair_id": [f"P{i}" for i in range(n)],
        "class": ["OP"] * n,
        "gene_a": [f"A{i}" for i in range(n)],
        "gene_b": [f"B{i}" for i in range(n)],
        "strand": ["+"] * n,
        "detected_a": [True] * n,
        "detected_b": [True] * n,
        "abundance_log2_a": a,
        "abundance_log2_b": b,
    }
    pairs = pd.DataFrame(rows)
    if one_detected:
        k = len(one_detected)
        extra = pd.DataFrame({
            "pair_id": [f"U{j}" for j in range(k)], "class": ["OP"] * k,
            "gene_a": [f"UA{j}" for j in range(k)], "gene_b": [f"UB{j}" for j in range(k)],
            "strand": ["+"] * k, "detected_a": [True] * k, "detected_b": [False] * k,
            "abundance_log2_a": [float(v) for v, _ in one_detected],
            "abundance_log2_b": [math.nan] * k,
        })  # fmt: skip
        pairs = pd.concat([pairs, extra], ignore_index=True)
    pairs["delta_log2"] = pairs["abundance_log2_b"] - pairs["abundance_log2_a"]
    pairs["flags"] = ""
    tags = [*pairs["gene_a"], *pairs["gene_b"]]
    imp = {f"B{i}" for i in imputed}
    genes = pd.DataFrame({
        "locus_tag": tags, "partial": False, "abundance_imputed": [t in imp for t in tags],
        "abundance_type": "ibaq", "abundance_scale": "absolute",
    })  # fmt: skip
    return pairs, genes


def hetero(n: int, seed: int, gamma: float = -0.5) -> tuple[np.ndarray, np.ndarray]:
    """E8's generator at size n: z ~ N(0, 1), Delta = -0.3 + 0.5 exp(gamma z) t_5, m = 20 + 2 z."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    t5 = rng.standard_t(5, n)
    return -0.3 + 0.5 * np.exp(gamma * z) * t5, 20.0 + 2.0 * z


def e8_frame() -> tuple[pd.DataFrame, pd.DataFrame]:
    """E8 (n = 2,000, default_rng(11)) as a pair table (m = z)."""
    rng = np.random.default_rng(11)
    z = rng.standard_normal(2000)
    t5 = rng.standard_t(5, 2000)
    return frame(-0.3 + 0.5 * np.exp(-0.5 * z) * t5, z)


def run(pairs: pd.DataFrame, genes: pd.DataFrame, p_miss: Any = None, **kw: Any) -> Any:
    return disc.discordance(pairs, genes, _cfg(), p_miss, **kw)


def t_rows(res: Any) -> pd.DataFrame:
    return res.rows[res.rows["test"] == "both_detected_t"]


def _flags(v: Any) -> list[str]:
    if isinstance(v, str):
        return [f for f in v.replace(",", ";").split(";") if f]
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return []
    return [str(f) for f in v]


def _usable_frame(pairs: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """The pair frame with the per-gene imputed flags joined, as ``discordance`` builds it."""
    imp = genes.drop_duplicates("locus_tag").set_index("locus_tag")["abundance_imputed"]
    out = pairs.copy()
    for side in ("a", "b"):
        out[f"abundance_imputed_{side}"] = out[f"gene_{side}"].map(imp).fillna(False).astype(bool)
    if "partial" in genes.columns:
        part = genes.drop_duplicates("locus_tag").set_index("locus_tag")["partial"]
        bad = out["gene_a"].map(part).fillna(False).astype(bool) | out["gene_b"].map(
            part).fillna(False).astype(bool)  # fmt: skip
        out = out.loc[~bad.to_numpy()]
    return out


def assert_main_constant(res: Any, pairs: pd.DataFrame, genes: pd.DataFrame) -> None:
    """Current main on the same frame: fit = fit_student_t(fit_rows), p = student_t_two_sided,
    q = BH over the t rows, effect = Delta - loc; compared bit for bit."""
    delta = disc.fit_rows(_usable_frame(pairs, genes))
    ref = fit_student_t(delta) if delta.size >= disc.MIN_OP_PAIRS else None
    if ref is not None:
        assert res.fit == ref
    np.testing.assert_array_equal(res.ratios, delta)
    t = t_rows(res)
    if t.empty:
        return
    assert ref is not None
    d = pairs.set_index("pair_id").loc[t["pair_id"], "delta_log2"].to_numpy(float)
    p = student_t_two_sided(d, ref)
    np.testing.assert_array_equal(t["p_value"].to_numpy(float), p)
    np.testing.assert_array_equal(t["q_value"].to_numpy(float), bh(p))
    np.testing.assert_array_equal(t["stat"].to_numpy(float), d)
    np.testing.assert_array_equal(t["effect"].to_numpy(float), d - ref.loc)


def _expected_p(delta: np.ndarray, m: np.ndarray, s: ScaledStudentT) -> np.ndarray:
    z = np.clip((m - s.center) / s.sd, s.z_min, s.z_max)
    sigma = s.scale0 * np.exp(s.gamma * z)
    return np.asarray(2.0 * sps.t.sf(np.abs(delta - s.loc) / sigma, s.df), dtype=float)


# --------------------------------------------------------------------------- interface


def test_sigma_model_parameter_defaults_to_loglinear() -> None:
    sig = inspect.signature(disc.discordance)
    assert sig.parameters["sigma_model"].default == "loglinear"
    assert sig.parameters["sigma_model"].kind is inspect.Parameter.KEYWORD_ONLY


def test_module_constant_sigma_model_is_loglinear() -> None:
    """Q2 default: no config key yet; a module constant ``SIGMA_MODEL = "loglinear"``."""
    import giap.steps.operons as step

    values = {getattr(mod, "SIGMA_MODEL", None) for mod in (disc, step)} - {None}
    assert values == {"loglinear"}


def test_result_has_scaled_and_keeps_fit_and_ratios() -> None:
    pairs, genes = frame(*hetero(120, 1))
    res = run(pairs, genes)
    assert isinstance(res.scaled, ScaledStudentT)
    assert isinstance(res.fit, StudentT)
    assert res.fit == fit_student_t(disc.fit_rows(pairs))
    np.testing.assert_array_equal(res.ratios, disc.fit_rows(pairs))
    assert res.scaled.n == 120
    assert res.scaled.model == "loglinear"


def test_fit_pairs_is_delta_and_mean_sorted_lexicographically() -> None:
    d = [0.5, -0.25, 0.5, 0.125, -0.25]  # dyadic, so b - a is exact and the ties are real
    m = [21.0, 19.0, 18.0, 25.0, 23.0]
    pairs, _ = frame(d, m, one_detected=[(17.0, 0.0)])
    x, mm = disc.fit_pairs(pairs)
    np.testing.assert_array_equal(x, [-0.25, -0.25, 0.125, 0.5, 0.5])
    np.testing.assert_array_equal(mm, [19.0, 23.0, 25.0, 18.0, 21.0])
    np.testing.assert_array_equal(x, disc.fit_rows(pairs))


# --------------------------------------------------------------------------- examples


def test_e6_thirty_pairs_use_the_loglinear_null() -> None:
    """S-5: no pair-count gate between MIN_OP_PAIRS = 20 and 50."""
    pairs, genes = frame(*hetero(30, 6))
    res = run(pairs, genes)
    assert res.scaled.model == "loglinear" and res.scaled.reason is None
    t = t_rows(res)
    assert len(t) == 30
    assert np.isfinite(t["p_value"].to_numpy(float)).all()
    assert not any(disc.TOO_FEW in _flags(f) for f in res.rows["flags"])
    d = pairs.set_index("pair_id").loc[t["pair_id"], "delta_log2"].to_numpy(float)
    mm = pairs.set_index("pair_id").loc[t["pair_id"], ["abundance_log2_a", "abundance_log2_b"]]
    np.testing.assert_allclose(t["p_value"].to_numpy(float),
                               _expected_p(d, mm.mean(axis=1).to_numpy(float), res.scaled),
                               rtol=1e-9)  # fmt: skip


def test_e6_constant_is_bit_identical_to_student_t() -> None:
    pairs, genes = frame(*hetero(30, 6))
    res = run(pairs, genes, sigma_model="constant")
    assert res.scaled.model == "constant" and res.scaled.reason == "disabled"
    assert_main_constant(res, pairs, genes)


def test_e7_flat_abundance_is_the_constant_null() -> None:
    delta = _t_draws(200)
    pairs, genes = frame(delta, np.full(200, 20.0))
    res = run(pairs, genes)
    assert res.scaled.model == "constant" and res.scaled.reason == "flat_abundance"
    assert res.scaled.gamma == 0.0
    t = t_rows(res)
    d = t["stat"].to_numpy(float)
    np.testing.assert_allclose(t["p_value"].to_numpy(float),
                               student_t_two_sided(d, fit_student_t(disc.fit_rows(pairs))),
                               rtol=1e-12)  # fmt: skip


def test_e12_constant_on_e8_is_main_bit_for_bit() -> None:
    pairs, genes = e8_frame()
    res = run(pairs, genes, sigma_model="constant")
    assert res.scaled.reason == "disabled" and res.scaled.model == "constant"
    assert res.scaled.gamma == 0.0
    assert_main_constant(res, pairs, genes)


def test_e8_through_the_step_recovers_gamma() -> None:
    pairs, genes = e8_frame()
    res = run(pairs, genes)
    s = res.scaled
    assert s.model == "loglinear" and s.n == 2000
    assert -0.62 <= s.gamma <= -0.38
    t = t_rows(res)
    m = t.merge(pairs, on="pair_id")[["abundance_log2_a", "abundance_log2_b"]].mean(axis=1)
    np.testing.assert_allclose(t["p_value"].to_numpy(float),
                               _expected_p(t["stat"].to_numpy(float), m.to_numpy(float), s),
                               rtol=1e-9)  # fmt: skip
    np.testing.assert_allclose(t["effect"].to_numpy(float), t["stat"].to_numpy(float) - s.loc,
                               rtol=0, atol=1e-12)  # fmt: skip
    assert (t["method"] == "student_t").all()
    assert (t["null_model"] == "within_operon_t").all()
    assert all(f == ["ci_unestimated"] for f in t["flags"].map(_flags))
    np.testing.assert_allclose(t["q_value"].to_numpy(float), bh(t["p_value"].to_numpy(float)))


def test_e8_rows_do_not_depend_on_row_order() -> None:
    pairs, genes = e8_frame()
    a = run(pairs, genes).rows.set_index("pair_id").sort_index()
    b = run(pairs.iloc[::-1].reset_index(drop=True), genes).rows.set_index("pair_id").sort_index()
    np.testing.assert_array_equal(a["p_value"].to_numpy(float), b["p_value"].to_numpy(float))


def test_e13_nineteen_pairs_is_too_few() -> None:
    pairs, genes = frame(*hetero(19, 13))
    res = run(pairs, genes)
    assert res.scaled.model == "none"
    assert (res.rows["test"] == "none").all()
    for f in res.rows["flags"]:
        assert {disc.TOO_FEW, disc.UNESTIMATED} <= set(_flags(f))


class Recorder:
    def __init__(self) -> None:
        self.ratios: list[np.ndarray] = []

    def __call__(self, u: str, a_v: float, ratios: Any, v: str, /) -> float:
        self.ratios.append(np.asarray(ratios, dtype=float).copy())
        return 0.25


def test_e14_pmiss_gets_the_raw_fit_rows_and_one_detected_rows_are_unchanged() -> None:
    """Q1: the one-detected P(miss) branch is unchanged; ``ratios`` is the raw sorted Delta."""
    rng = np.random.default_rng(11)
    z = rng.standard_normal(2000)
    t5 = rng.standard_t(5, 2000)
    one = [(16.0, 0.0), (-1.5, 0.0), (2.5, 0.0)]
    pairs, genes = frame(-0.3 + 0.5 * np.exp(-0.5 * z) * t5, z, one_detected=one)
    rec_l, rec_c = Recorder(), Recorder()
    log_rows = run(pairs, genes, rec_l).rows
    const_rows = run(pairs, genes, rec_c, sigma_model="constant").rows
    assert len(rec_l.ratios) == 3
    for r in rec_l.ratios:
        np.testing.assert_array_equal(r, disc.fit_rows(pairs))
    for a, b in zip(rec_l.ratios, rec_c.ratios, strict=True):
        np.testing.assert_array_equal(a, b)
    sel = log_rows["test"] == disc.PMISS_TEST
    assert int(sel.sum()) == 3
    cols = ["pair_id", "test", "stat", "p_value", "method", "null_model"]
    pd.testing.assert_frame_equal(log_rows.loc[sel, cols], const_rows.loc[sel, cols])
    assert log_rows.loc[sel, "q_value"].isna().all()


def test_e5_imputed_tested_pair_is_clamped_to_the_fit_range() -> None:
    """S-2 / S-3: the imputed pair is tested, never fitted; its z comes from its own m and is
    clamped to [z_min, z_max] of the fit pairs."""
    d, m = hetero(300, 5)
    d = np.append(d, 1.0)
    m = np.append(m, float(m.max()) + 5.0)
    pairs, genes = frame(d, m, imputed=[300])
    res = run(pairs, genes)
    s = res.scaled
    assert s.n == 300 and s.model == "loglinear"
    assert s.gamma < -0.2
    row = res.rows.set_index("pair_id").loc["P300"]
    assert row["test"] == "both_detected_t"
    sig = s.scale0 * math.exp(s.gamma * s.z_max)
    want = float(2.0 * sps.t.sf(abs(1.0 - s.loc) / sig, s.df))
    assert float(row["p_value"]) == pytest.approx(want, rel=1e-9)
    z_raw = (float(m[-1]) - s.center) / s.sd
    assert z_raw > s.z_max + 1.0
    raw = float(2.0 * sps.t.sf(abs(1.0 - s.loc) / (s.scale0 * math.exp(s.gamma * z_raw)), s.df))
    assert float(row["p_value"]) > raw * 1.5


# --------------------------------------------------------------------------- fallbacks


@pytest.fixture
def fail_scaled(monkeypatch: pytest.MonkeyPatch) -> None:
    """The four-parameter optimiser reports failure; the constant fit runs unchanged."""

    def fake(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        res = _REAL_MINIMIZE(fun, x0, *a, **k)
        if np.asarray(x0).size == 4:
            res = scipy.optimize.OptimizeResult(dict(res))
            res.update(success=False, status=2)
        return res

    monkeypatch.setattr(scipy.optimize, "minimize", fake)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", fake)


def test_non_converged_scaled_fit_is_the_constant_null_not_too_few(
    fail_scaled: None, caplog: pytest.LogCaptureFixture
) -> None:
    """Q3: a failed scaled fit falls back to gamma = 0 with a WARNING; a valid null exists, so it
    is not ``too_few_op_pairs``."""
    pairs, genes = frame(*hetero(200, 3))
    with caplog.at_level(logging.WARNING, logger="giap"):
        res = run(pairs, genes)
    assert res.scaled.model == "constant" and res.scaled.reason == "not_converged"
    t = t_rows(res)
    assert len(t) == 200
    assert not any(disc.TOO_FEW in _flags(f) for f in res.rows["flags"])
    ref = fit_student_t(disc.fit_rows(pairs))
    np.testing.assert_allclose(
        t["p_value"].to_numpy(float),
        student_t_two_sided(t["stat"].to_numpy(float), ref),
        rtol=1e-12,
    )
    warned = [r for r in caplog.records if r.levelno >= logging.WARNING
              and "not_converged" in r.getMessage()]  # fmt: skip
    assert warned, [r.getMessage() for r in caplog.records]


def test_constant_model_does_not_try_the_scaled_fit(
    fail_scaled: None, caplog: pytest.LogCaptureFixture
) -> None:
    pairs, genes = frame(*hetero(200, 3))
    with caplog.at_level(logging.WARNING, logger="giap"):
        res = run(pairs, genes, sigma_model="constant")
    assert res.scaled.reason == "disabled"
    assert not [r for r in caplog.records if "not_converged" in r.getMessage()]


# --------------------------------------------------------------------------- bit identity


def _m605_fixtures() -> list[tuple[str, tuple[pd.DataFrame, pd.DataFrame]]]:
    return [
        ("draws60", m605_chain(_t_draws(60))),
        ("draws200_plus_1.2", m605_chain([*_t_draws(200), 1.2])),
        (
            "mixed_eligibility",
            m605_chain(_t_draws(80), undetected=[3, 40], no_value=[10], imputed=[20, 21]),
        ),
        ("twenty", m605_chain(_t_draws(20, seed=3))),
        ("nineteen", m605_chain(_t_draws(19, seed=3))),
        ("partial", m605_chain(_t_draws(60), partial=True)),
        ("relative", m605_chain(_t_draws(60), abundance_scale="relative")),
    ]


@pytest.mark.parametrize("case", range(7),
                         ids=[n for n, _ in _m605_fixtures()])  # fmt: skip
def test_constant_is_main_bit_for_bit_on_the_m605_fixtures(case: int) -> None:
    _, (pairs, genes) = _m605_fixtures()[case]
    res = run(pairs, genes, sigma_model="constant")
    assert_main_constant(res, pairs, genes)
    other = run(pairs, genes)
    rest = res.rows["test"] != "both_detected_t"
    pd.testing.assert_series_equal(res.rows.loc[rest, "test"], other.rows.loc[rest, "test"])
    assert list(res.rows.loc[rest, "flags"].map(_flags)) == list(
        other.rows.loc[rest, "flags"].map(_flags))  # fmt: skip


def test_loglinear_changes_p_on_varied_abundance() -> None:
    """The M6-05 chain fixtures are random-walk abundances (m varies), so the default moves their
    p: the reason the lead pins ``sigma_model = "constant"`` on M6-05's constant-fit assertions."""
    pairs, genes = m605_chain([*_t_draws(200), 1.2])
    a = t_rows(run(pairs, genes))["p_value"].to_numpy(float)
    b = t_rows(run(pairs, genes, sigma_model="constant"))["p_value"].to_numpy(float)
    assert a.shape == b.shape
    assert not np.array_equal(a, b)
