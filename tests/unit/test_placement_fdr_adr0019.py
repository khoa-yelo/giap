"""ADR-0019: stratum-0 π̂_wrong for database input from exact-match shares per length class
(docs/adr/ADR-0019-stratum0-pi-hat-exact-share.md decisions 1-4; extends ADR-0010; supersedes
GIAP_SPEC.md §7.2.5 step 4 for database input).

Interface assumed here (test-author; the ADR names no functions):

* ``giap.mapping.placement_fdr.stratum0_pi_classes(costs, classes) -> DataFrame``: pure. ``costs``
  holds the rows of the population P (the caller drops non-database rows and ``same_db`` exact
  peptides) with at least ``length``, ``m_target``, ``m_decoy`` (the ``best_costs`` schema; NaN =
  no hit). It returns one row per class of ``classes`` (= ``mapping.length_classes``), in config
  order, with columns ``length_class`` (the ``length_class_labels`` label), ``n``, ``k`` (no exact
  target hit: ``m_target`` != 0 or NaN), ``r`` (share of the n with ``m_decoy`` == 0), ``upper``
  (U, one-sided 95 % Clopper-Pearson upper bound of k/n; 1 when k = n) and ``pi_wrong`` (1 when
  r = 1, else min(1, U/(1 - r))). An empty class has n = k = 0; its other values are not pinned.
* ``giap.mapping.placement_fdr.input_prefiltered(n, k, mean_pep, r) -> bool``: decision 4's run
  warning, True iff k < the 0.001 quantile of Binomial(n, mean_pep * (1 - r)) (equivalently
  P(X <= k) < 0.001); False when ``mean_pep`` is NaN or n = 0.
* ``placement(costs, cfg)`` keeps its signature and output columns (decision 4: ``mapping_fdr``
  columns unchanged); its database stratum-0 lists use the class π̂ as decoy weights and
  c = max π̂_l over classes with n_l > 0; a peptide outside every class weighs 1 and does not
  enter c. ``pep`` and ``q_value`` are not read for it (decisions 2 and 3).
* ``giap.context.summary.qc_rows(peptides, peptide_hits, unmapped, stratum0=table)``: with the
  ``stratum0_pi_classes`` table it adds metric ``stratum0_unmapped_share`` (value k/n, ``n`` = n)
  as group ``run``/``all`` and per class as group_type ``length_class``, group = label. The M4-16
  QC rows, ``genome_inconsistent_rate`` included, stay exactly as they are (M3-25 decision 2).
* Contaminant rows (``is_contaminant``) leave the π̂ population P like ``same_db`` exact peptides
  (M3-25 decision 7), so the costs frame carries an optional ``is_contaminant`` column (absent or
  NaN = False), filled by ``best_costs`` from the peptide table. Only the exclusion from P is
  pinned here: whether a contaminant peptide's own hits stay in the stratum lists is untested.

Hand values: with k = 0, U = 1 - 0.05^(1/n) (Clopper-Pearson); for k > 0 U is pinned and checked
against its defining equation P(X <= k; n, U) = 0.05.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig
from giap.mapping.placement_fdr import placement

ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
PREFIX = "MDEFGHKLNPQRSTVWYACDEFGHKLMNPQRSTVWY"
CLASSES = [[7, 7], [8, 8], [9, 9], [10, 11], [12, 14], [15, 19], [20, 40]]
LABELS = ["7", "8", "9", "10-11", "12-14", "15-19", "20-40"]


def u0(n: int) -> float:
    """Clopper-Pearson one-sided 95 % upper bound for k = 0 of n."""
    return 1.0 - 0.05 ** (1.0 / n)


def binom_cdf(k: int, n: int, p: float) -> float:
    return math.fsum(math.comb(n, j) * p**j * (1 - p) ** (n - j) for j in range(k + 1))


# U(k, n) for k > 0 (scipy.stats.beta.ppf(0.95, k + 1, n - k)); checked by the defining equation.
U_1_100 = 0.046559811453538935
U_1_3 = 0.8646496378284161
U_10_40 = 0.3870602400920182


@pytest.mark.parametrize(("k", "n", "u"), [(1, 100, U_1_100), (1, 3, U_1_3), (10, 40, U_10_40)])
def test_pinned_upper_bounds_solve_the_clopper_pearson_equation(k: int, n: int, u: float) -> None:
    assert binom_cdf(k, n, u) == pytest.approx(0.05, abs=1e-12)


# ------------------------------------------------------------------------------------ helpers


def pep_of(i: int, length: int = 9) -> str:
    assert 7 <= length <= len(PREFIX) + 3 and i < 19**3
    digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (2, 1, 0))
    return PREFIX[: length - 3] + digits


LONG45 = ("MDEFGHKLNPQRSTVWYAC" * 3)[:45]  # outside every default class


def row(
    peptide: str, *, mt: float | None = None, et: float | None = None, md: float | None = None,
    ed: float | None = None, **kw: Any,
) -> dict[str, Any]:  # fmt: skip
    nan = math.nan
    r: dict[str, Any] = {
        "peptide": peptide, "length": len(peptide), "source": "database", "search_db": "other",
        "is_contaminant": False,
        "m_target": nan if mt is None else float(mt), "m_decoy": nan if md is None else float(md),
        "e_target": nan if et is None else float(et), "e_decoy": nan if ed is None else float(ed),
        "pep": nan, "q_value": nan, "pi_denovo": nan,
    }  # fmt: skip
    r.update(kw)
    return r


def frame(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


def targets0(n: int, *, start: int = 0, length: int = 9, e0: float = 0.001, **kw: Any) -> list:
    """n exact (m = 0) targets, E = e0 * (i + 1)."""
    return [row(pep_of(start + i, length), mt=0, et=e0 * (i + 1), **kw) for i in range(n)]


def run(rows: Sequence[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    return placement(frame(rows), GiapConfig())


def q_of(per: pd.DataFrame, peptides: Sequence[str]) -> np.ndarray:
    q = per.set_index("peptide")["placement_q"].astype(float)
    return np.asarray(q.loc[list(peptides)].to_numpy(), dtype=float)


def fdr_row(fdr: pd.DataFrame, m: int, cls: str, dt: str = "reverse") -> dict[str, Any]:
    sel = fdr[
        (fdr["section"] == "placement")
        & (fdr["decoy_type"] == dt)
        & (fdr["n_diff"].astype(int) == m)
        & (fdr["length_class"].astype(str) == cls)
    ]
    assert len(sel) == 1
    out: dict[str, Any] = sel.iloc[0].to_dict()
    return out


def classes_table(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    from giap.mapping.placement_fdr import stratum0_pi_classes

    return stratum0_pi_classes(frame(rows), CLASSES)


def class_row(tab: pd.DataFrame, label: str) -> dict[str, Any]:
    sel = tab[tab["length_class"].astype(str) == label]
    assert len(sel) == 1, tab
    out: dict[str, Any] = sel.iloc[0].to_dict()
    return out


# ------------------------------------------------------------- decision 1: the class estimator


def test_table_has_one_row_per_class_in_config_order() -> None:
    tab = classes_table(targets0(5))
    assert tab["length_class"].astype(str).tolist() == LABELS
    assert {"length_class", "n", "k", "r", "upper", "pi_wrong"} <= set(tab.columns)


def test_k0_r0_pi_is_the_upper_bound() -> None:
    """n = 100, k = 0, r = 0: π̂ = U = 1 - 0.05^(1/100) (edge case r = 0 gives π̂ = U)."""
    r9 = class_row(classes_table(targets0(100)), "9")
    assert (int(r9["n"]), int(r9["k"])) == (100, 0)
    assert float(r9["r"]) == 0.0
    assert float(r9["upper"]) == pytest.approx(u0(100), rel=1e-9)
    assert float(r9["pi_wrong"]) == pytest.approx(u0(100), rel=1e-9)


def test_decoy_share_divides_the_upper_bound() -> None:
    """n = 100, k = 0, 10 with an exact decoy hit: r = 0.1, π̂ = U/(1 - r)."""
    rows = targets0(100)
    for r in rows[:10]:
        r.update(m_decoy=0.0, e_decoy=1e-9)
    r9 = class_row(classes_table(rows), "9")
    assert float(r9["r"]) == pytest.approx(0.1, rel=1e-12)
    assert float(r9["pi_wrong"]) == pytest.approx(u0(100) / 0.9, rel=1e-9)


def test_unmapped_peptide_counts_in_k() -> None:
    """99 exact targets + 1 decoy-only peptide: n = 100, k = 1, r = 0.01."""
    rows = [*targets0(99), row(pep_of(500), md=0, ed=10.0)]
    r9 = class_row(classes_table(rows), "9")
    assert (int(r9["n"]), int(r9["k"])) == (100, 1)
    assert float(r9["r"]) == pytest.approx(0.01, rel=1e-12)
    assert float(r9["upper"]) == pytest.approx(U_1_100, rel=1e-9)
    assert float(r9["pi_wrong"]) == pytest.approx(U_1_100 / 0.99, rel=1e-9)


def test_tolerant_or_hitless_peptides_count_in_k() -> None:
    """k counts m_T != 0 (here m_T = 1 and m_T = 2) and peptides with no hit at all; r counts
    exact decoy hits only (a cost-1 decoy hit is not one)."""
    rows = [
        *targets0(30, length=12),
        *[row(pep_of(100 + i, 12), mt=1, et=0.1, md=1, ed=0.2) for i in range(6)],
        *[row(pep_of(200 + i, 12), mt=2, et=0.1) for i in range(2)],
        *[row(pep_of(300 + i, 12)) for i in range(2)],
    ]
    r12 = class_row(classes_table(rows), "12-14")
    assert (int(r12["n"]), int(r12["k"])) == (40, 10)
    assert float(r12["r"]) == 0.0
    assert float(r12["upper"]) == pytest.approx(U_10_40, rel=1e-9)
    assert float(r12["pi_wrong"]) == pytest.approx(U_10_40, rel=1e-9)


def test_r_one_gives_pi_one() -> None:
    rows = targets0(20)
    for r in rows:
        r.update(m_decoy=0.0, e_decoy=1.0)
    r9 = class_row(classes_table(rows), "9")
    assert float(r9["r"]) == 1.0
    assert float(r9["pi_wrong"]) == 1.0


def test_k_equals_n_gives_upper_one_and_pi_one() -> None:
    rows = [row(pep_of(i), md=0, ed=1.0) for i in range(2)]
    rows += [row(pep_of(10 + i)) for i in range(3)]
    r9 = class_row(classes_table(rows), "9")
    assert (int(r9["n"]), int(r9["k"])) == (5, 5)
    assert float(r9["upper"]) == 1.0
    assert float(r9["pi_wrong"]) == 1.0


def test_exact_share_below_decoy_share_is_clipped_to_one() -> None:
    """m = 1 - k/n = 0.2 < r = 0.5 ("m < r is clipped to 1")."""
    rows = [*targets0(2), *[row(pep_of(50 + i)) for i in range(8)]]
    for r in rows[:2] + rows[2:5]:
        r.update(m_decoy=0.0, e_decoy=1.0)
    r9 = class_row(classes_table(rows), "9")
    assert (int(r9["n"]), int(r9["k"]), float(r9["r"])) == (10, 8, 0.5)
    assert float(r9["pi_wrong"]) == 1.0


def test_empty_class_has_no_peptides() -> None:
    tab = classes_table(targets0(10))
    for label in LABELS:
        if label != "9":
            r = class_row(tab, label)
            assert (int(r["n"]), int(r["k"])) == (0, 0)


def test_classes_are_counted_separately() -> None:
    rows = [*targets0(100), *targets0(3, start=200, length=12), row(pep_of(300, 12))]
    tab = classes_table(rows)
    assert (int(class_row(tab, "9")["n"]), int(class_row(tab, "9")["k"])) == (100, 0)
    r12 = class_row(tab, "12-14")
    assert (int(r12["n"]), int(r12["k"])) == (4, 1)


def test_pep_and_q_value_do_not_enter_the_table() -> None:
    a = [*targets0(50), row(pep_of(400), md=0, ed=1.0)]
    b = [dict(r, pep=0.5, q_value=0.001) for r in a]
    pd.testing.assert_frame_equal(classes_table(a), classes_table(b))


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.tuples(
            st.integers(7, 39),
            st.sampled_from([0.0, 1.0, 2.0, math.nan]),
            st.sampled_from([0.0, 1.0, math.nan]),
        ),
        min_size=1,
        max_size=60,
    )
)
def test_property_table_invariants(items: list[tuple[int, float, float]]) -> None:
    """n sums to the in-class rows; 0 <= k <= n; 0 <= r <= 1; U >= k/n; π̂ in (0, 1] and
    π̂ >= min(1, (k/n)/(1 - r)) (the upper bound never lies below the point estimate)."""
    rows = [
        row(pep_of(i, L), mt=None if math.isnan(mt) else mt, et=0.1,
            md=None if math.isnan(md) else md, ed=0.2)
        for i, (L, mt, md) in enumerate(items)
    ]  # fmt: skip
    tab = classes_table(rows)
    assert int(tab["n"].sum()) == len(items)
    for r in tab.to_dict("records"):
        n, k = int(r["n"]), int(r["k"])
        assert 0 <= k <= n
        if n == 0:
            continue
        rr, u, pi = float(r["r"]), float(r["upper"]), float(r["pi_wrong"])
        assert 0.0 <= rr <= 1.0
        assert k / n <= u <= 1.0
        assert 0.0 < pi <= 1.0
        point = 1.0 if rr >= 1.0 else min(1.0, (k / n) / (1.0 - rr))
        assert pi >= point - 1e-12


# ------------------------------------------- decision 2: entry into the stratum-0 FDR (placement)


def test_placement_k0_r0_q_is_u_over_t() -> None:
    """100 exact targets, no decoys: c = π̂_9 = U, FDR(e_j) = U/j -> q = U/100 for all."""
    rows = targets0(100)
    per, fdr = run(rows)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), u0(100) / 100, rtol=1e-9)
    for cls in ("9", "all"):
        r = fdr_row(fdr, 0, cls)
        assert float(r["pi_wrong"]) == pytest.approx(u0(100), rel=1e-9)  # (W + c)/(D + 1), W = 0
        assert float(r["fdr_stratum"]) == pytest.approx(u0(100) / 100, rel=1e-9)
        assert int(r["n_accepted"]) == 100


def test_placement_decoy_hits_weigh_the_class_pi() -> None:
    """10 of the 100 also hit the decoy exactly, before every target: W = 10 π̂, c = π̂,
    π̂ = U/0.9; FDR(e_j) = 11 π̂/j -> q = 11 π̂/100."""
    rows = targets0(100)
    for r in rows[:10]:
        r.update(m_decoy=0.0, e_decoy=1e-9)
    per, fdr = run(rows)
    pi = u0(100) / 0.9
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), 11 * pi / 100, rtol=1e-9)
    r9 = fdr_row(fdr, 0, "9")
    assert int(r9["decoys"]) == 10
    assert float(r9["pi_wrong"]) == pytest.approx(pi, rel=1e-9)
    assert float(r9["fdr_stratum"]) == pytest.approx(11 * pi / 100, rel=1e-9)


def test_placement_unmapped_peptide_raises_pi() -> None:
    """99 targets + 1 decoy-only peptide after every target: π̂ = U(1, 100)/0.99, q = π̂/99."""
    rows = [*targets0(99), row(pep_of(500), md=0, ed=10.0)]
    per, fdr = run(rows)
    pi = U_1_100 / 0.99
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows[:99]]), pi / 99, rtol=1e-9)
    tot = fdr_row(fdr, 0, "all")
    assert (int(tot["targets"]), int(tot["decoys"])) == (99, 1)
    assert float(tot["pi_wrong"]) == pytest.approx(pi, rel=1e-9)
    assert float(tot["fdr_stratum"]) == pytest.approx(2 * pi / 99, rel=1e-9)


def test_placement_c_is_the_max_class_pi_and_w_the_own_class_pi() -> None:
    """Class 9: 100 targets, one with an exact decoy hit -> π̂_9 = U(0, 100)/0.99. Class 12-14: 3
    targets -> π̂_12 = U(0, 3) = c. One list: FDR(e_j) = (π̂_9 + c)/j -> q = (π̂_9 + c)/103.
    Empty classes do not enter c."""
    t9 = targets0(100)
    t9[0].update(m_decoy=0.0, e_decoy=1e-9)
    t12 = targets0(3, start=200, length=12, e0=1.0)
    per, fdr = run(t9 + t12)
    pi9, c = u0(100) / 0.99, u0(3)
    names = [r["peptide"] for r in t9 + t12]
    np.testing.assert_allclose(q_of(per, names), (pi9 + c) / 103, rtol=1e-9)
    r9 = fdr_row(fdr, 0, "9")
    assert float(r9["pi_wrong"]) == pytest.approx((pi9 + c) / 2, rel=1e-9)
    assert float(r9["fdr_stratum"]) == pytest.approx((pi9 + c) / 100, rel=1e-9)
    r12 = fdr_row(fdr, 0, "12-14")
    assert float(r12["pi_wrong"]) == pytest.approx(c, rel=1e-9)
    assert float(r12["fdr_stratum"]) == pytest.approx(c / 3, rel=1e-9)


def test_placement_peptide_outside_every_class_weighs_one_and_not_in_c() -> None:
    """A 45-mer decoy-only peptide (outside 7-40) before every target: W = 1, c = π̂_9 = U(0, 100)
    (it does not enter c, nor class 9's n, k, r): q = (1 + U)/100."""
    rows = [*targets0(100), row(LONG45, md=0, ed=1e-9)]
    per, fdr = run(rows)
    u = u0(100)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows[:100]]), (1 + u) / 100,
                               rtol=1e-9)  # fmt: skip
    tot = fdr_row(fdr, 0, "all")
    assert int(tot["decoys"]) == 1
    assert float(tot["pi_wrong"]) == pytest.approx((1 + u) / 2, rel=1e-9)


def test_placement_same_db_exact_peptides_leave_the_population() -> None:
    """50 ``same`` exact peptides with exact decoy hits are not in P: class 9 keeps n = 100, r = 0,
    so q = U(0, 100)/100 for the ``other`` peptides; the ``same`` ones get q = 0."""
    other = targets0(100)
    same = [row(pep_of(1000 + i), mt=0, et=0.5, md=0, ed=1e-9, search_db="same") for i in range(50)]
    per, _ = run(other + same)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in other]), u0(100) / 100, rtol=1e-9)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in same]), 0.0)


def test_placement_contaminant_rows_leave_the_population() -> None:
    """M3-25 decision 7: 20 hit-less contaminant rows of class 9 must not raise k (they have no
    hits, so they are in no list either): π̂ stays U(0, 100), not U(20, 120)/1."""
    rows = [*targets0(100), *[row(pep_of(700 + i), is_contaminant=True) for i in range(20)]]
    per, fdr = run(rows)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows[:100]]), u0(100) / 100,
                               rtol=1e-9)  # fmt: skip
    assert float(fdr_row(fdr, 0, "9")["pi_wrong"]) == pytest.approx(u0(100), rel=1e-9)


def test_placement_contaminant_exact_hits_do_not_lower_k() -> None:
    """The other direction: 100 class-9 contaminants with exact target hits must not pull the
    class's k/n down (with them in P: k = 1 of 200 instead of 1 of 100)."""
    base = [*targets0(99), row(pep_of(500), md=0, ed=10.0)]
    cont = [row(pep_of(700 + i), mt=0, et=0.002, is_contaminant=True) for i in range(100)]
    pi = U_1_100 / 0.99
    assert float(fdr_row(run(base + cont)[1], 0, "all")["pi_wrong"]) == pytest.approx(pi, rel=1e-9)


# ------------------------------------------------- decisions 2 and 3: PEP and q are not read


@pytest.mark.parametrize(
    "extra",
    [
        {"pep": 0.0},
        {"pep": 0.02},
        {"pep": 0.9},
        {"q_value": 0.001},
        {"q_value": 0.0098},
        {"pep": 0.3, "q_value": 0.005},
    ],
)
def test_placement_ignores_pep_and_q_value_for_database_stratum0(extra: dict[str, float]) -> None:
    base = [*targets0(99), row(pep_of(500), md=0, ed=10.0)]
    base[3].update(m_decoy=0.0, e_decoy=0.0035)
    with_extra = [dict(r, **extra) for r in base]
    per_a, fdr_a = run(base)
    per_b, fdr_b = run(with_extra)
    pd.testing.assert_frame_equal(per_a, per_b)
    pd.testing.assert_frame_equal(fdr_a, fdr_b)


def test_placement_min_collapsed_pep_no_longer_deflates_q() -> None:
    """Run collapse: PEP 10x too small (§5.3 min over runs) must not change the q-values."""
    rng = np.random.default_rng(7)
    rows = [*targets0(300, e0=1e-4)]
    for i in range(30):
        rows.append(row(pep_of(600 + i), mt=0 if i < 8 else None, et=0.03 + 1e-4 * i,
                        md=0 if i % 2 else None, ed=0.02 + 1e-4 * i))  # fmt: skip
    pep = rng.uniform(0.0, 0.2, len(rows))
    a = [dict(r, pep=float(p)) for r, p in zip(rows, pep, strict=True)]
    b = [dict(r, pep=float(p) / 10) for r, p in zip(rows, pep, strict=True)]
    pd.testing.assert_frame_equal(run(a)[0], run(b)[0])


def test_placement_without_pep_or_q_uses_the_class_estimate_not_one() -> None:
    """Previously π̂ = 1 without PEP and q-values (§7.2.5 step 4 fallback); now U(0, 2)."""
    rows = targets0(2)
    per, fdr = run(rows)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), u0(2) / 2, rtol=1e-9)
    assert float(fdr_row(fdr, 0, "all")["pi_wrong"]) == pytest.approx(u0(2), rel=1e-9)


# ----------------------------------------------------------------------- m >= 1 is unchanged


def test_stratum1_keeps_pi_one_while_its_peptides_count_in_stratum0_k() -> None:
    """100 exact targets + 50 cost-1 targets + 1 cost-1 decoy-only, all class 9. Stratum 1:
    π̂ = 1, c = 1, the decoy before every target -> q = 2/50. Stratum 0: k = 51 of n = 151, r = 0
    -> π̂_9 = U(51, 151), q = π̂_9/100."""
    from scipy.stats import beta

    t0 = targets0(100)
    t1 = [row(pep_of(300 + i), mt=1, et=0.001 * (i + 1)) for i in range(50)]
    d1 = row(pep_of(900), md=1, ed=1e-9)
    per, fdr = run([*t0, *t1, d1])
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t1]), 2 / 50, rtol=1e-12)
    assert float(fdr_row(fdr, 1, "all")["pi_wrong"]) == 1.0
    u = float(beta.ppf(0.95, 52, 100))
    assert binom_cdf(51, 151, u) == pytest.approx(0.05, abs=1e-9)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t0]), u / 100, rtol=1e-9)


def test_denovo_stratum0_still_uses_pi_denovo() -> None:
    """The class estimator is for database input only: de novo rows keep ``pi_denovo``."""
    t = [row(pep_of(i, 7), mt=0, et=float(i + 1), source="denovo", pi_denovo=0.2)
         for i in range(4)]  # fmt: skip
    per, fdr = run(t)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t]), 0.05, rtol=1e-12)
    assert float(fdr_row(fdr, 0, "all", "shuffle")["pi_wrong"]) == pytest.approx(0.2)


# -------------------------------------------------- decision 4: input_prefiltered run warning


@pytest.mark.parametrize(
    ("n", "k", "mean_pep", "r", "fires"),
    [
        (1000, 0, 0.1, 0.0, True),
        (1000, 71, 0.1, 0.0, True),  # P(X <= 71; 1000, 0.1) = 0.00086 < 0.001
        (1000, 72, 0.1, 0.0, False),  # P(X <= 72) = 0.00127
        (1000, 100, 0.1, 0.0, False),
        (1000, 54, 0.1, 0.2, True),  # Binomial(1000, 0.08): 0.001 quantile 55
        (1000, 55, 0.1, 0.2, False),
        (0, 0, 0.1, 0.0, False),
        (1000, 0, math.nan, 0.0, False),
    ],
)
def test_input_prefiltered_threshold(
    n: int, k: int, mean_pep: float, r: float, fires: bool
) -> None:
    from giap.mapping.placement_fdr import input_prefiltered

    assert input_prefiltered(n, k, mean_pep, r) is fires


def test_input_prefiltered_boundaries_are_the_binomial_quantile() -> None:
    """Hand check of the pinned boundaries above."""
    assert binom_cdf(71, 1000, 0.1) < 0.001 <= binom_cdf(72, 1000, 0.1)
    assert binom_cdf(54, 1000, 0.08) < 0.001 <= binom_cdf(55, 1000, 0.08)


# --------------------------------------------------- decision 4: detection_summary QC rows


def test_qc_rows_report_the_stratum0_unmapped_share_overall_and_per_class() -> None:
    from giap.context.summary import qc_rows

    rows = [*targets0(99), row(pep_of(500), md=0, ed=10.0), *targets0(30, start=200, length=12)]
    rows += [row(pep_of(300 + i, 12)) for i in range(10)]
    tab = classes_table(rows)
    peptides = pd.DataFrame(
        {"peptide": ["X"], "is_decoy": [False], "search_db": ["other"], "source": ["database"]}
    )
    hits = pd.DataFrame({"peptide": pd.Series([], dtype=str), "route": pd.Series([], dtype=str)})
    qc = qc_rows(peptides, hits, None, stratum0=tab)
    # decision 2: the M4-16 QC rows are untouched, each over its own population
    base = qc_rows(peptides, hits, None)
    keep = qc[qc["metric"] != "stratum0_unmapped_share"].reset_index(drop=True)
    pd.testing.assert_frame_equal(keep, base.reset_index(drop=True))
    assert "genome_inconsistent_rate" in set(qc["metric"])
    sel = qc[qc["metric"] == "stratum0_unmapped_share"]
    assert set(sel["section"]) == {"qc"}
    overall = sel[(sel["group_type"] == "run") & (sel["group"] == "all")]
    assert len(overall) == 1
    assert float(overall["value"].iloc[0]) == pytest.approx(11 / 140, rel=1e-12)
    assert int(overall["n"].iloc[0]) == 140
    per = sel[sel["group_type"] == "length_class"].set_index("group")
    assert float(per.loc["9", "value"]) == pytest.approx(1 / 100, rel=1e-12)
    assert int(per.loc["9", "n"]) == 100
    assert float(per.loc["12-14", "value"]) == pytest.approx(10 / 40, rel=1e-12)
    assert int(per.loc["12-14", "n"]) == 40
    # classes with n = 0 report no share (NA), never 0
    if "7" in per.index:
        assert pd.isna(per.loc["7", "value"])
