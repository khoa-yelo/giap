"""M6-05H (pass 2): the one-detected P(miss) branch of the discordance test
(docs/specs/M6-05-discordance.md item 4, example 6, example 8 and the degenerate rows;
GIAP_SPEC.md §7.4.3 "One detected, one not"; M6-02 for the toy model).

Interface pinned here (test-author reading of the spec and of the merged pass-1 code):

* ``discordance(pairs, genes, cfg, p_miss, *, n_seen=None, k_min=1)`` takes ``p_miss`` as the
  callable of the Protocol ``giap.operons.discordance.PMissFn`` (renamed from ``PMissModel``,
  review N4): ``p_miss(u, a_v, ratios, v) -> float``, all positional, where ``u`` is the
  undetected gene's locus tag, ``a_v`` the detected partner's ``abundance_log2``, ``ratios`` the
  oriented within-operon log-ratios r (§7.4.3): Δ of the fit pairs when ``u`` is downstream
  (``gene_b``), -Δ when ``u`` is upstream (``gene_a``), and ``v`` the detected partner's locus
  tag. ``v`` is there so the adaptor can use the out-of-fold model of **v's** TU (M6-02 lead
  decision 9, "call ``model_for_tu`` with the TU of the detected partner"; review B1).
  ``giap.operons.detection_limit.p_miss_for_genes`` is that adaptor:
  ``p_miss(model_for_tu(model, tu(v)), a_v, ratios, observable_peptides(aa(u)))``.
* **The one-detected event** (M6-02 lead decision 9, review S1): with ``n_seen`` (locus tag ->
  number of distinct observable peptides of the gene seen in C11; absent = 0) a gene is *seen*
  when ``n_seen >= k_min`` (``k_min`` of the model in force). A one-detected row is an eligible
  OP pair with exactly one gene seen: u = the gene with fewer than ``k_min`` observable peptides
  seen, whatever its ``detected`` call or imputed value says, v = the seen partner. The
  both-detected t test is unchanged (``detected`` on both sides with finite values) and takes
  precedence: a row with a t test is never scored by P(miss). With ``n_seen=None`` (no peptide
  evidence; unit callers, M6-09 attack C) the event falls back to ``detected_a xor detected_b``.
* A one-detected row whose partner v has no ``abundance_log2`` is ``none`` with no flag, with or
  without a model (spec item 5; review N1).
* A one-detected OP row: ``test = "one_detected_pmiss"``, ``stat`` = ``a_v``, ``effect`` NA,
  ``p_value`` = the callable's value, ``method = "detection_limit"``,
  ``null_model = "p_detect_model"`` (spec item 4).
* **BH family.** The spec body (item 7, example 8) pools both kinds in one family of 150. The
  M6-index lead decision 1 (2026-10-21) and M6-02 lead decision 4 reverse that: P(miss) is
  reported **with no q-value** and stays out of the BH family unless M6-09 attack C shows the
  mixture is safe, and then an ADR records why. ``docs/methods/operons.md`` says the same. No such
  ADR exists, so these tests pin the lead decision. They do it through the single switch
  ``PMISS_IN_BH_FAMILY``. If an ADR later admits the mixture, the lead flips it to ``True`` and
  the same tests then check the pooled family of 150.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import detection_limit_data as D
from helpers.pairs_data import frame, gene, reps
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import giap.stats.tdist as tdist
from giap.config import GiapConfig
from giap.operons import detection_limit as dl
from giap.operons import discordance as disc
from giap.stats.multitest import bh
from giap.stats.tdist import StudentT

#: M6-index lead decision 1 / M6-02 lead decision 4: P(miss) gets no q and is not in the family.
PMISS_IN_BH_FAMILY = False

PMISS = "one_detected_pmiss"
TTEST = "both_detected_t"
FAMILY_TESTS = {TTEST, PMISS} if PMISS_IN_BH_FAMILY else {TTEST}
_DYADIC = 2.0**20
EX6_RATIOS = (-1.0, -0.2, 0.3, 1.5)


def _cfg(**abundance: Any) -> GiapConfig:
    doc: dict[str, Any] = {"seed": 1}
    if abundance:
        doc["abundance"] = abundance
    return GiapConfig.model_validate(doc)


def _t_draws(
    n: int, seed: int = 7, loc: float = 0.05, scale: float = 0.4, df: float = 6.0
) -> np.ndarray:
    return loc + scale * np.random.default_rng(seed).standard_t(df, size=n)


def chain(
    deltas: Sequence[float],
    *,
    base: float = 20.0,
    contig: str = "c1",
    prefix: str = "G",
    undetected: Sequence[int] = (),
    no_value: Sequence[int] = (),
    imputed: Sequence[int] = (),
    values: Mapping[int, float] | None = None,
    partial: bool = False,
    uniop_p: float = 0.9,
    abundance_type: str = "ibaq",
    abundance_scale: str = "molar",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pairs, genes) for a + strand chain of len(deltas) + 1 genes starting at abundance
    ``base``; pair k is (gene k, gene k + 1), so gene k + 1 is downstream, with
    delta_log2 == deltas[k] (rounded to 2**-20). Indices in ``undetected`` / ``no_value`` /
    ``imputed`` are gene indices; ``values`` overrides the ``abundance_log2`` of a gene index
    (e.g. a finite imputed value for an undetected gene). Pair ids are ``PAIR_<k+1>``,
    prefixed ``<prefix>_`` unless the prefix is ``G``."""
    d = np.round(np.asarray(deltas, dtype=float) * _DYADIC) / _DYADIC
    ab = np.concatenate([[base], base + np.cumsum(d)])
    tags = [f"{prefix}{i}" for i in range(len(ab))]
    genes = [gene(t, 100 + 400 * i, 400 + 400 * i, "+", contig) for i, t in enumerate(tags)]
    det = np.array([i not in undetected for i in range(len(ab))])
    val = np.where(det & ~np.isin(np.arange(len(ab)), list(no_value)), ab, np.nan)
    for i, x in (values or {}).items():
        val[i] = x
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


def cat(*parts: tuple[pd.DataFrame, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (
        pd.concat([p for p, _ in parts], ignore_index=True),
        pd.concat([g for _, g in parts], ignore_index=True),
    )


def _toy(k_min: int = 1) -> Any:
    return dl.PMissModel(
        kind="peptide", coef=dict(D.TOY_COEF), k_min=k_min, n_fit=1000, auc=0.8, usable=True
    )


def toy_pmiss(u: str, a_v: float, ratios: Any, v: str) -> float:
    """M6-02's toy model for every gene u (peptides of length 10, 20, 30)."""
    return dl.p_miss(_toy(), a_v, ratios, list(D.TOY_PEPS))


class Spy:
    """Records every call of the ``p_miss`` callable."""

    def __init__(self, fn: Callable[[str, float, Any, str], float] = toy_pmiss) -> None:
        self.fn = fn
        self.calls: list[tuple[str, float, np.ndarray, str]] = []

    def __call__(self, u: str, a_v: float, ratios: Any, v: str, /) -> float:
        r = np.asarray(ratios, dtype=float).copy()
        self.calls.append((str(u), float(a_v), r, str(v)))
        return self.fn(u, a_v, r, v)

    def by_gene(self) -> dict[str, tuple[float, np.ndarray]]:
        return {u: (a, r) for u, a, r, _ in self.calls}

    def partner(self) -> dict[str, str]:
        return {u: v for u, _, _, v in self.calls}


def run(
    pairs: pd.DataFrame,
    genes: pd.DataFrame,
    p_miss: Any = toy_pmiss,
    cfg: GiapConfig | None = None,
    **kw: Any,
) -> Any:
    return disc.discordance(pairs, genes, cfg or _cfg(), p_miss, **kw)


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
    """Replace the t fit by a fixed, converged StudentT (the fit itself is M6-01's)."""
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

    set_fit(0.15, 0.5, 6.0)
    yield set_fit


def ex6_frames(*, upstream_missing: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Example 6: 20 usable OP pairs whose Δ are the r sample {-1.0, -0.2, 0.3, 1.5} five times
    (same mean as the four values), plus one OP pair on contig c2 with the partner v at
    a_v = 17.0. By default only the upstream gene H0 is detected (u = H1 is downstream)."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    if upstream_missing:
        one = chain([0.5], base=16.5, contig="c2", prefix="H", undetected=[0])  # H1 = 17.0
    else:
        one = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])  # H0 = 17.0
    return cat(fit_part, one)


# ------------------------------------------------------------------ example 6 and item 4


def test_example6_upstream_detected_p_value_constant(fixed_fit: Any) -> None:
    pairs, genes = ex6_frames()
    r = rows_by_id(run(pairs, genes)).loc["H_PAIR_1"]
    assert r["test"] == PMISS
    assert float(r["p_value"]) == pytest.approx(0.6889, abs=5e-5)


def test_example6_row_columns_follow_item4(fixed_fit: Any) -> None:
    pairs, genes = ex6_frames()
    r = rows_by_id(run(pairs, genes)).loc["H_PAIR_1"]
    assert float(r["stat"]) == pytest.approx(17.0, abs=1e-12)  # the detected partner's a_v
    assert pd.isna(r["effect"])
    assert r["method"] == "detection_limit"
    assert r["null_model"] == "p_detect_model"
    assert "no_detection_limit_model" not in flags_of(r["flags"])
    assert "too_few_op_pairs" not in flags_of(r["flags"])


def test_example6_callable_gets_u_av_and_the_fit_ratios(fixed_fit: Any) -> None:
    """u is the undetected gene, a_v the partner's abundance, and r = Δ of the fit pairs because
    u is downstream (§7.4.3)."""
    pairs, genes = ex6_frames()
    spy = Spy()
    res = run(pairs, genes, spy)
    calls = spy.by_gene()
    assert set(calls) == {"H1"}, "only the one-detected pair calls p_miss"
    a_v, ratios = calls["H1"]
    assert a_v == pytest.approx(17.0, abs=1e-12)
    np.testing.assert_allclose(np.sort(ratios), np.sort(list(EX6_RATIOS) * 5), atol=1e-6)
    np.testing.assert_allclose(np.sort(ratios), np.sort(res.ratios), atol=0)


def test_upstream_missing_partner_gets_negated_ratios(fixed_fit: Any) -> None:
    """u upstream (gene_a undetected): r = -Δ (§7.4.3), and the p follows the toy model."""
    pairs, genes = ex6_frames(upstream_missing=True)
    spy = Spy()
    r = rows_by_id(run(pairs, genes, spy)).loc["H_PAIR_1"]
    a_v, ratios = spy.by_gene()["H0"]
    assert a_v == pytest.approx(17.0, abs=1e-12)
    want_r = -np.asarray(list(EX6_RATIOS) * 5)
    np.testing.assert_allclose(np.sort(ratios), np.sort(want_r), atol=1e-6)
    assert r["test"] == PMISS
    assert float(r["stat"]) == pytest.approx(17.0, abs=1e-12)
    want = dl.p_miss(_toy(), 17.0, want_r, list(D.TOY_PEPS))
    assert float(r["p_value"]) == pytest.approx(want, abs=1e-6)
    assert want != pytest.approx(0.6889, abs=1e-3), "orientation must matter in this example"


def test_ratios_come_from_the_fit_pairs_only(fixed_fit: Any) -> None:
    """An imputed pair is tested but is not in the fit, so its Δ is not in the r sample."""
    fit_part = chain([*_t_draws(40), 9.0], imputed=[41])  # pair 41 (Δ = 9) touches gene 41
    one = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    spy = Spy()
    res = run(*cat(fit_part, one), spy)
    _, ratios = spy.by_gene()["H1"]
    assert ratios.size == res.fit.n == 40
    assert not np.any(np.isclose(ratios, 9.0))
    np.testing.assert_allclose(np.sort(ratios), np.sort(res.ratios))


def test_p_value_is_the_callable_value_verbatim(fixed_fit: Any) -> None:
    pairs, genes = ex6_frames()
    r = rows_by_id(run(pairs, genes, lambda u, a_v, ratios, v: 0.123456789)).loc["H_PAIR_1"]
    assert r["test"] == PMISS
    assert float(r["p_value"]) == 0.123456789


# ------------------------------------------------------------------ the BH family (item 7, ex. 8)


def _example8(p_miss: Any = toy_pmiss) -> tuple[pd.DataFrame, pd.Series]:
    p1 = chain(_t_draws(100, seed=11))
    p2 = chain(_t_draws(50, seed=12), contig="c2", prefix="H", undetected=range(1, 51, 2))
    rows = run(*cat(p1, p2), p_miss).rows
    return rows, rows["pair_id"].isin(p2[0]["pair_id"])


def test_example8_fifty_one_detected_rows_carry_pmiss() -> None:
    rows, one = _example8()
    assert one.sum() == 50
    assert (rows.loc[one, "test"] == PMISS).all()
    p = rows.loc[one, "p_value"].to_numpy(float)
    assert np.isfinite(p).all() and (p >= 0).all() and (p <= 1).all()
    assert (rows.loc[~one, "test"] == TTEST).all()


def test_example8_bh_family_composition() -> None:
    """Lead decision (``PMISS_IN_BH_FAMILY``): the family is the 100 t p-values and P(miss)
    has no q; with the switch on, it is one family of 150 (the spec body's example 8)."""
    rows, one = _example8()
    fam = rows["test"].isin(FAMILY_TESTS).to_numpy()
    assert fam.sum() == (150 if PMISS_IN_BH_FAMILY else 100)
    np.testing.assert_allclose(
        rows.loc[fam, "q_value"].to_numpy(float), bh(rows.loc[fam, "p_value"].to_numpy(float))
    )
    if not PMISS_IN_BH_FAMILY:
        assert rows.loc[one, "q_value"].isna().all()
    s = rows[fam].sort_values("p_value", kind="stable")
    assert (np.diff(s["q_value"].to_numpy(float)) >= -1e-15).all()  # monotone in the pooled p


def test_t_q_values_do_not_depend_on_pmiss_values() -> None:
    """With P(miss) outside the family, the t rows' q cannot move when P(miss) changes."""
    if PMISS_IN_BH_FAMILY:
        pytest.skip("P(miss) shares the family (ADR); this invariant does not hold")
    lo, _ = _example8(lambda u, a_v, ratios, v: 1e-6)
    hi, _ = _example8(lambda u, a_v, ratios, v: 0.999)
    t_lo = lo[lo["test"] == TTEST].set_index("pair_id")["q_value"]
    t_hi = hi[hi["test"] == TTEST].set_index("pair_id")["q_value"]
    pd.testing.assert_series_equal(t_lo, t_hi.loc[t_lo.index])


# ------------------------------------------------------------------ degenerate rows with a model


def test_pmiss_nan_for_a_gene_is_untested_with_no_detection_limit_model(fixed_fit: Any) -> None:
    """Degenerate row: ``p_miss`` NaN (no observable peptide) → ``none`` +
    ``no_detection_limit_model``; the other one-detected pair is still tested."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    k = chain([0.5], base=17.0, contig="c3", prefix="K", undetected=[1])

    def fn(u: str, a_v: float, ratios: Any, v: str) -> float:
        return math.nan if u == "K1" else toy_pmiss(u, a_v, ratios, v)

    rows = rows_by_id(run(*cat(fit_part, h, k), fn))
    bad = rows.loc["K_PAIR_1"]
    assert bad["test"] == "none"
    assert pd.isna(bad["p_value"]) and pd.isna(bad["q_value"])
    assert "no_detection_limit_model" in flags_of(bad["flags"])
    assert rows.loc["H_PAIR_1", "test"] == PMISS


def test_too_few_op_pairs_leaves_one_detected_untested() -> None:
    """Spec item 2: below 20 usable pairs every OP test is ``none`` (no r sample to use)."""
    pairs, genes = cat(
        chain(_t_draws(19)), chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    )
    r = rows_by_id(run(pairs, genes)).loc["H_PAIR_1"]
    assert r["test"] == "none"
    assert pd.isna(r["p_value"]) and pd.isna(r["q_value"])
    assert "too_few_op_pairs" in flags_of(r["flags"])


def test_non_converged_fit_leaves_one_detected_untested(fixed_fit: Any) -> None:
    fixed_fit(0.05, 0.4, 6.0, converged=False)
    pairs, genes = ex6_frames()
    res = run(pairs, genes)
    assert (res.rows["test"] == "none").all()
    assert res.rows["p_value"].isna().all()
    assert "too_few_op_pairs" in flags_of(rows_by_id(res).loc["H_PAIR_1", "flags"])


def test_relative_scale_skips_one_detected(fixed_fit: Any) -> None:
    fit_part = chain(list(EX6_RATIOS) * 5, abundance_type="maxlfq", abundance_scale="relative")
    one = chain(
        [0.5], base=17.0, contig="c2", prefix="H", undetected=[1],
        abundance_type="maxlfq", abundance_scale="relative",
    )  # fmt: skip
    spy = Spy()
    res = run(*cat(fit_part, one), spy)
    assert (res.rows["test"] == "none").all()
    assert res.rows["p_value"].isna().all()
    assert "relative_scale" in flags_of(rows_by_id(res).loc["H_PAIR_1", "flags"])
    assert not spy.calls


def test_neither_detected_is_untested_and_never_scored(fixed_fit: Any) -> None:
    fit_part = chain(list(EX6_RATIOS) * 5)
    nn = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[0, 1])
    spy = Spy()
    r = rows_by_id(run(*cat(fit_part, nn), spy)).loc["H_PAIR_1"]
    assert r["test"] == "none" and pd.isna(r["p_value"])
    assert not spy.calls


def test_detected_partner_without_value_is_untested_under_exclude(fixed_fit: Any) -> None:
    """Item 5: the partner is detected but its abundance is missing → no a_v, no test."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1], no_value=[0])
    spy = Spy()
    r = rows_by_id(run(*cat(fit_part, h), spy, _cfg(missing="exclude"))).loc["H_PAIR_1"]
    assert r["test"] == "none" and pd.isna(r["p_value"])
    assert "H1" not in spy.by_gene()


def test_partial_gene_one_detected_pair_is_not_tested(fixed_fit: Any) -> None:
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1], partial=True)
    spy = Spy()
    res = run(*cat(fit_part, h), spy)
    r = rows_by_id(res).loc["H_PAIR_1"]
    assert r["test"] == "none" and pd.isna(r["p_value"])
    assert not spy.calls


def test_no_class_one_detected_pair_is_not_scored(fixed_fit: Any) -> None:
    fit_part = chain(list(EX6_RATIOS) * 5)
    no = chain([0.5], base=17.0, contig="c2", prefix="N", undetected=[1], uniop_p=0.05)
    assert (no[0]["class"] == "NO").all()
    spy = Spy()
    res = run(*cat(fit_part, no), spy)
    assert not res.rows["pair_id"].isin(no[0]["pair_id"]).any()
    assert not spy.calls


def test_p_miss_none_keeps_the_pass1_unestimated_path(fixed_fit: Any) -> None:
    """No model at all (one run, no peptides) → ``none`` + ``no_detection_limit_model``, while the
    same frames with a model give ``one_detected_pmiss`` (the branch is live, not removed)."""
    pairs, genes = ex6_frames()
    without = rows_by_id(run(pairs, genes, None)).loc["H_PAIR_1"]
    assert without["test"] == "none" and pd.isna(without["p_value"])
    assert "no_detection_limit_model" in flags_of(without["flags"])
    assert rows_by_id(run(pairs, genes)).loc["H_PAIR_1", "test"] == PMISS


def test_deterministic_with_a_model() -> None:
    pairs, genes = cat(
        chain(_t_draws(60)), chain(_t_draws(10, seed=3), contig="c2", prefix="H",
                                   undetected=range(1, 11, 2))
    )  # fmt: skip
    a, b = run(pairs, genes), run(pairs.copy(), genes.copy())
    pd.testing.assert_frame_equal(a.rows, b.rows)
    assert (a.rows["test"] == PMISS).sum() == 10


def test_input_frames_not_mutated_with_a_model(fixed_fit: Any) -> None:
    pairs, genes = ex6_frames()
    p0, g0 = pairs.copy(deep=True), genes.copy(deep=True)
    run(pairs, genes)
    pd.testing.assert_frame_equal(pairs, p0)
    pd.testing.assert_frame_equal(genes, g0)


# ------------------------------------------------------------------ review of #311 (M6-05H)


def test_pmissfn_protocol_is_exported() -> None:
    """N4: the callable's Protocol is ``PMissFn``; ``PMissModel`` is M6-02's dataclass only."""
    assert "PMissFn" in disc.__all__
    assert isinstance(disc.PMissFn, type)
    assert "PMissModel" not in disc.__all__


@pytest.mark.parametrize("upstream_missing", [False, True])
def test_callable_gets_the_detected_partner_v(fixed_fit: Any, upstream_missing: bool) -> None:
    """B1: the fourth positional argument is the detected partner v (either orientation)."""
    pairs, genes = ex6_frames(upstream_missing=upstream_missing)
    spy = Spy()
    run(pairs, genes, spy)
    u, v = ("H0", "H1") if upstream_missing else ("H1", "H0")
    assert spy.partner() == {u: v}


def test_partner_without_value_has_no_flag_with_or_without_a_model(fixed_fit: Any) -> None:
    """N1 / item 5: a one-detected row whose detected partner has no abundance is ``none`` with
    no reason flag, and the ``p_miss = None`` path agrees with the model path."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1], no_value=[0])
    pairs, genes = cat(fit_part, h)
    spy = Spy()
    with_model = rows_by_id(run(pairs, genes, spy)).loc["H_PAIR_1"]
    without = rows_by_id(run(pairs, genes, None)).loc["H_PAIR_1"]
    for r in (with_model, without):
        assert r["test"] == "none" and pd.isna(r["p_value"]) and pd.isna(r["q_value"])
        assert "no_detection_limit_model" not in flags_of(r["flags"])
        assert "unestimated" not in flags_of(r["flags"])
    assert flags_of(with_model["flags"]) == flags_of(without["flags"])
    assert not spy.calls


def test_undetected_gene_with_finite_imputed_value_is_one_detected(fixed_fit: Any) -> None:
    """S4 (#147 review): u is not detected but has a finite imputed abundance. The row is still
    one-detected (scored by P(miss) on the partner's a_v, not on u's value) and the pair's Δ
    is not part of the r sample."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1], imputed=[1],
              values={1: 13.0})  # fmt: skip
    pairs, genes = cat(fit_part, h)
    assert float(pairs.set_index("pair_id").loc["H_PAIR_1", "abundance_log2_b"]) == 13.0
    for kw in ({}, {"n_seen": _seen_all(fit_part, h, k_min=1, unseen=("H1",)), "k_min": 1}):
        spy = Spy()
        res = run(pairs, genes, spy, **kw)
        r = rows_by_id(res).loc["H_PAIR_1"]
        assert r["test"] == PMISS, kw
        assert float(r["stat"]) == pytest.approx(17.0, abs=1e-12)
        a_v, ratios = spy.by_gene()["H1"]
        assert a_v == pytest.approx(17.0, abs=1e-12)
        assert spy.partner() == {"H1": "H0"}
        assert ratios.size == 20 and not np.any(np.isclose(ratios, -4.0))
        assert float(r["p_value"]) == pytest.approx(0.6889, abs=5e-5)


def test_imputed_partner_is_scored_and_carries_abundance_imputed(fixed_fit: Any) -> None:
    """S4 / F6: the detected partner's a_v is imputed. The row is scored on it (as example 5
    tests an imputed both-detected pair) and the merged pair row carries the pair-scope
    ``abundance_imputed`` flag (M6-08 lead decision 6)."""
    from giap.operons.flags import merge_tests

    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1], imputed=[0])
    pairs, genes = cat(fit_part, h)
    d = run(pairs, genes)
    r = rows_by_id(d).loc["H_PAIR_1"]
    assert r["test"] == PMISS
    assert float(r["stat"]) == pytest.approx(17.0, abs=1e-12)
    merged = merge_tests(pairs, d.rows, d.rows.iloc[:0], genes).set_index("pair_id")
    m = merged.loc["H_PAIR_1"]
    assert m["test"] == PMISS
    assert float(m["p_value"]) == pytest.approx(float(r["p_value"]), abs=0)
    assert "abundance_imputed" in flags_of(m["flags"])
    assert "abundance_imputed" not in flags_of(merged.loc["PAIR_1", "flags"])


# --------------------------------------------- S1: the one-detected event is "< k_min seen"


def _seen_all(
    *parts: tuple[pd.DataFrame, pd.DataFrame],
    k_min: int,
    unseen: Sequence[str] = (),
    counts: Mapping[str, int] | None = None,
) -> dict[str, int]:
    """``n_seen``: every detected gene of ``parts`` has k_min observable peptides seen, every
    gene in ``unseen`` 0; ``counts`` overrides."""
    out: dict[str, int] = {}
    for _, g in parts:
        for tag, det in zip(g["locus_tag"], g["detected"], strict=True):
            out[str(tag)] = k_min if bool(det) else 0
    out.update({u: 0 for u in unseen})
    out.update(counts or {})
    return out


@pytest.mark.parametrize("k_min", [1, 2])
@pytest.mark.parametrize("u_seen", ["below", "at"])
def test_one_detected_event_boundary_at_k_min_for_u(
    fixed_fit: Any, k_min: int, u_seen: str
) -> None:
    """M6-02 lead decision 9 (S1): u is missed when it has fewer than k_min observable peptides
    seen. u (not ``detected``, no value) with k_min - 1 seen is scored; with k_min seen it was
    not missed, so the row is not ``one_detected_pmiss`` and p_miss is never asked about u."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    n_u = k_min - 1 if u_seen == "below" else k_min
    seen = _seen_all(fit_part, h, k_min=k_min, counts={"H1": n_u})
    spy = Spy()
    r = rows_by_id(run(*cat(fit_part, h), spy, n_seen=seen, k_min=k_min)).loc["H_PAIR_1"]
    if u_seen == "below":
        assert r["test"] == PMISS
        assert spy.partner() == {"H1": "H0"}
        assert float(r["p_value"]) == pytest.approx(0.6889, abs=5e-5)
    else:
        assert r["test"] != PMISS
        assert pd.isna(r["p_value"])
        assert "H1" not in spy.by_gene()


@pytest.mark.parametrize("k_min", [1, 2])
@pytest.mark.parametrize("v_seen", ["below", "at"])
def test_one_detected_event_boundary_at_k_min_for_v(
    fixed_fit: Any, k_min: int, v_seen: str
) -> None:
    """S1, the partner's side: v must itself be seen (>= k_min). v ``detected`` with a value but
    k_min - 1 peptides seen and u with none: neither gene is seen, the row is ``none``."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    n_v = k_min - 1 if v_seen == "below" else k_min
    seen = _seen_all(fit_part, h, k_min=k_min, counts={"H0": n_v, "H1": 0})
    spy = Spy()
    r = rows_by_id(run(*cat(fit_part, h), spy, n_seen=seen, k_min=k_min)).loc["H_PAIR_1"]
    if v_seen == "at":
        assert r["test"] == PMISS
        assert spy.partner() == {"H1": "H0"}
    else:
        assert r["test"] == "none" and pd.isna(r["p_value"])
        assert not spy.calls


def test_gene_absent_from_n_seen_counts_as_unseen(fixed_fit: Any) -> None:
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    seen = _seen_all(fit_part, h, k_min=1)
    del seen["H1"]
    spy = Spy()
    r = rows_by_id(run(*cat(fit_part, h), spy, n_seen=seen, k_min=1)).loc["H_PAIR_1"]
    assert r["test"] == PMISS
    assert spy.partner() == {"H1": "H0"}


def test_n_seen_leaves_the_t_rows_alone_and_t_takes_precedence(fixed_fit: Any) -> None:
    """S1: the event only re-defines the one-detected branch. A both-detected pair with finite
    values keeps its t test even when one gene has fewer than k_min peptides seen (k_min = 2,
    a gene called on one peptide), and the t rows equal the ``n_seen=None`` run."""
    fit_part = chain(list(EX6_RATIOS) * 5)
    h = chain([0.5], base=17.0, contig="c2", prefix="H", undetected=[1])
    pairs, genes = cat(fit_part, h)
    seen = _seen_all(fit_part, h, k_min=2, counts={"G5": 1})
    spy = Spy()
    a = run(pairs, genes, spy, n_seen=seen, k_min=2).rows
    b = run(pairs, genes).rows
    ta, tb = a[a["test"] == TTEST], b[b["test"] == TTEST]
    pd.testing.assert_frame_equal(ta.reset_index(drop=True), tb.reset_index(drop=True))
    assert set(ta["pair_id"]) >= {"PAIR_5", "PAIR_6"}  # the pairs of G5
    assert set(spy.by_gene()) == {"H1"}


# --------------------------------------------- B1: the adaptor uses v's out-of-fold model

_SEQ_U = "MSTEKLAAGVNPDLKAEGFTWQYPSRHHEAVDNLGTRWASDEFGHIKLMNPQR"
_SEQ_V = "MKVLAAGIDEFGNPTWKQEELAVSTGHPYR"


def _fold_model() -> Any:
    """A peptide model whose two TU folds differ: TU_U is in fold 0, TU_V in fold 1."""
    lo = dict(D.TOY_COEF) | {"intercept": D.TOY_COEF["intercept"] - 2.0}
    hi = dict(D.TOY_COEF) | {"intercept": D.TOY_COEF["intercept"] + 2.0}
    return dl.PMissModel(
        kind="peptide", coef=dict(D.TOY_COEF), k_min=1, n_fit=1000, auc=0.8, usable=True,
        fold_of_tu={"TU_U": 0, "TU_V": 1}, fold_coef=(lo, hi),
    )  # fmt: skip


def _adaptor_genes() -> pd.DataFrame:
    return pd.DataFrame(
        {"locus_tag": ["U", "V", "W"], "tu_id": ["TU_U", "TU_V", None],
         "aa": [_SEQ_U, _SEQ_V, _SEQ_V]}
    )  # fmt: skip


def test_p_miss_for_genes_uses_the_tu_of_the_detected_partner() -> None:
    """M6-02 lead decision 9 / review B1: u and v in different TUs of different folds. The value
    is ``p_miss(model_for_tu(model, tu(v)), a_v, r, observable_peptides(aa(u)))``, and it differs
    from the u's-TU reading (so this test tells the two apart)."""
    cfg, model = _cfg(), _fold_model()
    peps_u = dl.observable_peptides(_SEQ_U, cfg)
    assert peps_u, "the test protein needs observable peptides"
    fn = dl.p_miss_for_genes(model, _adaptor_genes(), cfg)
    r = np.asarray(EX6_RATIOS)
    got = fn("U", 17.0, r, "V")
    want = dl.p_miss(dl.model_for_tu(model, "TU_V"), 17.0, r, peps_u)
    wrong = dl.p_miss(dl.model_for_tu(model, "TU_U"), 17.0, r, peps_u)
    assert abs(want - wrong) > 0.05
    assert got == pytest.approx(want, abs=1e-12)


def test_p_miss_for_genes_partner_without_tu_uses_the_full_model() -> None:
    cfg, model = _cfg(), _fold_model()
    peps_u = dl.observable_peptides(_SEQ_U, cfg)
    fn = dl.p_miss_for_genes(model, _adaptor_genes(), cfg)
    r = np.asarray(EX6_RATIOS)
    assert fn("U", 17.0, r, "W") == pytest.approx(dl.p_miss(model, 17.0, r, peps_u), abs=1e-12)
    assert fn("U", 17.0, r, "unknown") == pytest.approx(
        dl.p_miss(model, 17.0, r, peps_u), abs=1e-12
    )


def test_p_miss_for_genes_u_without_sequence_is_nan() -> None:
    fn = dl.p_miss_for_genes(_fold_model(), _adaptor_genes(), _cfg())
    assert math.isnan(fn("unknown", 17.0, np.asarray(EX6_RATIOS), "V"))


# ------------------------------------- S1: the per-gene count of observable peptides seen


def test_observable_peptides_seen_counts_distinct_observable_assignments() -> None:
    """``dl.observable_peptides_seen(genes, peptide_hits, cfg)``: per gene, the number of its
    distinct observable peptides (``observable_peptides``' set: 0 missed cleavages, length in
    ``abundance.observable_len``) with a C11 assignment to it (any route, I and L the same;
    M6-02 lead decision 2 and decision 8) — the count the peptide model's training set uses."""
    cfg = _cfg()
    seqs = dl._observable_seqs(_SEQ_U, cfg)
    assert len(seqs) >= 2
    il = seqs[0].replace("L", "I") if "L" in seqs[0] else seqs[0]
    hits = pd.DataFrame(
        {"gene": ["U", "U", "U", "U", "V", None],
         "peptide": [il, seqs[0], seqs[1], "AAK", seqs[1], seqs[1]]}
    )  # fmt: skip
    genes = _adaptor_genes()
    got = dl.observable_peptides_seen(genes, hits, cfg)
    assert got.get("U") == 2  # seqs[0] (twice, once I/L-swapped) and seqs[1]; AAK not observable
    assert got.get("V", 0) == 0  # seqs[1] is not an observable peptide of V's protein
    assert got.get("W", 0) == 0


# ------------------------------------------------------------------ properties

_PROP = settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])


@settings(max_examples=200, deadline=None)
@given(
    a_v=st.floats(-10.0, 60.0),
    ratios=st.lists(st.floats(-12.0, 12.0), min_size=1, max_size=30),
    lengths=st.lists(st.integers(7, 30), min_size=1, max_size=12),
    k_min=st.sampled_from([1, 2]),
)
def test_p_miss_is_a_probability(
    a_v: float, ratios: list[float], lengths: list[int], k_min: int
) -> None:
    """M6-02's P(miss) lies in [0, 1] for any a_v, r sample and peptide set."""
    p = dl.p_miss(_toy(k_min), a_v, ratios, [(n, 0.0) for n in lengths])
    assert 0.0 <= p <= 1.0


@_PROP
@given(
    n=st.integers(25, 45),
    seed=st.integers(0, 2**16),
    n_one=st.integers(1, 12),
    shift=st.floats(-6.0, 6.0),
)
def test_family_and_pmiss_rows_property(n: int, seed: int, n_one: int, shift: float) -> None:
    """Both arms: every ``one_detected_pmiss`` p is the model's P(miss) in [0, 1]; the BH family
    is exactly the rows of ``FAMILY_TESTS``, q = BH(p) over it and q is monotone in p; every
    untested row has neither p nor q (the M4-17 rule)."""
    spy = Spy(lambda u, a_v, r, v: dl.p_miss(_toy(), a_v + shift, r, list(D.TOY_PEPS)))
    one = chain(_t_draws(2 * n_one, seed=seed + 1), contig="c2", prefix="H",
                undetected=range(1, 2 * n_one + 1, 2))  # fmt: skip
    rows = run(*cat(chain(_t_draws(n, seed=seed)), one), spy).rows
    calls = spy.by_gene()
    pm = rows[rows["test"] == PMISS]
    if (rows["test"] == TTEST).any():  # the fit converged, so every one-detected pair is scored
        assert len(pm) == 2 * n_one
    if len(pm):
        p = pm["p_value"].to_numpy(float)
        assert ((p >= 0) & (p <= 1)).all()
        und = {f"H{i}" for i in range(1, 2 * n_one + 1, 2)}
        assert set(calls) == und
        assert len(pm) == 2 * n_one
    untested = rows["test"] == "none"
    assert rows.loc[untested, "p_value"].isna().all()
    assert rows.loc[untested, "q_value"].isna().all()
    fam = rows["test"].isin(FAMILY_TESTS)
    assert rows.loc[fam, "q_value"].notna().all()
    assert rows.loc[~fam, "q_value"].isna().all()
    if fam.any():
        p_f = rows.loc[fam, "p_value"].to_numpy(float)
        q_f = rows.loc[fam, "q_value"].to_numpy(float)
        np.testing.assert_allclose(q_f, bh(p_f))
        o = np.argsort(p_f, kind="stable")
        assert (np.diff(q_f[o]) >= -1e-15).all()
