"""M8-04b: TU-counted detection ORA (docs/specs/M8-04-category-summary-ora.md §2, §4-§7, examples
1-10 and the degenerate-case table; ADR-0029, ADR-0035).

Interface pinned here (module ``giap.function.ora``):

- ``tu_universe(genes, folds=None) -> TuUniverse`` and ``tu_members(u, genes) -> frozenset[int]``;
- ``detection_enrichment(genes, colls, size, alpha, *, folds=None) -> DataFrame`` with the columns
  of ``ORA_COLUMNS`` in the order of the module constant ``COLUMNS``; ``folds`` is M4-13's
  ``DetectabilityResult.fold_of_tu``. Without folds the TU recalibration (decision 12) is not run
  and the raw q is used (run flag ``tu_recal_skipped``), which is how examples 1-10 keep the
  spec's raw-q numbers; ``test_ora_recalibration.py`` pins the fitted path;
- ``attrs["run_flags"]`` is a list of strings; ``flags`` cells are sequences of strings (the writer
  joins them, as M8-05); every row carries ``ci_unestimated`` (lead decision 6).

Resolutions (the spec leaves them open; reported to the lead):

- **Eligibility** = non-partial rows (``partial`` NA counts as not partial; §7.3, M4-13 decision:
  every non-partial gene gets a p_detect). A non-partial gene with NaN ``p_detect`` or no
  ``tu_id`` fails closed (§2), naming the locus tags. Set members that are not eligible genes of
  the frame are ignored by ``tu_members``.
- **Size filter and elim.** A term outside ``size`` (counted in TUs on its full set) is not
  tested, writes no row and does not enter the elim walk, so it removes nothing from its
  ancestors (topGO prunes by node size before the walk).
- ``stat`` = k of the full TU set (``k_detected_in_set``); ``effect`` = k - E[K | N] of the tested
  set (full, or elim-reduced for GO ``over``). ``effect`` is NA where ``p_value`` is NA.
- ``p_detect_clipped`` goes on the rows whose full TU set holds a TU with a clipped p̂ (0 or 1),
  and in ``run_flags``.
- Rows are ordered by collection (contiguous blocks), then term, then direction (over < under);
  the collection order itself is not pinned.
- ``null_degenerate`` (a q = 1 TU undetected) is no longer reachable from a gene frame: p̂ is
  clipped and the fitted q is a logistic value with its own complement (decision 12). It is
  tested on the primitive in M8-04a only.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.ora_data import (
    ORA_COLUMNS,
    T6,
    Gene,
    as_list,
    cpb_raw,
    frame,
    gene_sets,
    kegg,
    ora,
    row,
    singletons,
    t6_frame,
)
from scipy.stats import hypergeom

from giap.stats.multitest import bh

ALPHA = 0.05
ALL = (1, 500)


def _run(genes: pd.DataFrame, colls: list[Any], size: tuple[int, int] = ALL, **kw: Any) -> Any:
    return ora().detection_enrichment(genes, colls, size, ALPHA, **kw)


def _na(v: Any) -> bool:
    return v is None or bool(pd.isna(v))


def _flags(v: Any) -> list[str]:
    return as_list(v)


# ----------------------------------------------------------------------------- tu_universe


def test_tu_universe_t6_raw_fields() -> None:
    u = ora().tu_universe(t6_frame())
    assert u.tu_ids == ("TU_1", "TU_2", "TU_3", "TU_4", "TU_5", "TU_6")
    q = [0.95, 0.5, 0.36, 0.5, 1 - 0.7**3, 0.1]
    np.testing.assert_allclose(u.q_raw, q, rtol=1e-14)
    lc = [math.log(0.1 * 0.5), math.log(0.5), 2 * math.log(0.8), math.log(0.5),
          3 * math.log(0.7), math.log(0.9)]  # fmt: skip
    np.testing.assert_allclose(u.log_comp, lc, rtol=1e-14)
    np.testing.assert_allclose(u.q_raw, -np.expm1(u.log_comp), rtol=1e-15)
    np.testing.assert_array_equal(u.n_genes, [2, 1, 2, 1, 3, 1])
    np.testing.assert_array_equal(u.detected, [True, False, True, True, True, False])
    np.testing.assert_array_equal(u.clipped, [False] * 6)
    assert u.gene_tu["g1"] == 0 and u.gene_tu["g2"] == 0 and u.gene_tu["g10"] == 5
    assert set(u.gene_tu) == {g for gs in T6.values() for g, _ in gs}


def test_tu_universe_without_folds_uses_the_raw_q() -> None:
    u = ora().tu_universe(t6_frame())
    np.testing.assert_array_equal(u.q, u.q_raw)
    np.testing.assert_allclose(u.comp, np.exp(u.log_comp), rtol=1e-15)
    assert "tu_recal_skipped" in u.flags


def test_tu_ids_are_natural_sorted_and_row_order_free() -> None:
    g = singletons([0.3] * 12, [0, 4, 10])
    u = ora().tu_universe(g)
    assert u.tu_ids == tuple(f"TU_{i}" for i in range(1, 13))
    shuffled = g.sample(frac=1.0, random_state=7).reset_index(drop=True)
    v = ora().tu_universe(shuffled)
    assert v.tu_ids == u.tu_ids
    assert dict(v.gene_tu) == dict(u.gene_tu)
    np.testing.assert_array_equal(v.q, u.q)
    np.testing.assert_array_equal(v.detected, u.detected)


def test_partial_genes_are_not_in_the_universe() -> None:
    g = frame(
        {
            "TU_1": [Gene("g1", 10, True, p_detect=0.5), Gene("g2", 10, True, partial=True,
                                                               p_detect=math.nan)],
            "TU_2": [Gene("g3", None, False, partial=True, p_detect=math.nan)],
            "TU_3": [Gene("g4", None, False, p_detect=0.2)],
        }
    )  # fmt: skip
    u = ora().tu_universe(g)
    assert u.tu_ids == ("TU_1", "TU_3")  # TU_2 has no eligible gene
    assert "g2" not in u.gene_tu and "g3" not in u.gene_tu
    np.testing.assert_allclose(u.q_raw, [0.5, 0.2], rtol=1e-15)


def test_nan_p_detect_on_an_eligible_gene_fails_closed() -> None:
    g = t6_frame()
    g.loc[g["locus_tag"] == "g5", "p_detect"] = math.nan
    with pytest.raises(ValueError, match="g5"):
        ora().tu_universe(g)
    with pytest.raises(ValueError, match="g5"):
        _run(g, [kegg({"x": {"g1", "g4", "g7"}})])


def test_missing_tu_id_fails_closed() -> None:
    g = t6_frame()
    g["tu_id"] = g["tu_id"].astype(object)
    g.loc[g["locus_tag"] == "g8", "tu_id"] = None
    with pytest.raises(ValueError, match="g8"):
        ora().tu_universe(g)


# ------------------------------------------------------------------------------ tu_members


def test_tu_members_any_eligible_gene() -> None:
    u = ora().tu_universe(t6_frame())
    assert ora().tu_members(u, frozenset({"g1", "g4", "g7"})) == frozenset({0, 2, 4})
    assert ora().tu_members(u, frozenset({"g8", "g9"})) == frozenset({4})
    assert ora().tu_members(u, frozenset()) == frozenset()


def test_tu_members_ignores_genes_outside_the_universe() -> None:
    u = ora().tu_universe(t6_frame())
    assert ora().tu_members(u, frozenset({"g3", "nope"})) == frozenset({1})


def test_tu_members_does_not_depend_on_detection() -> None:
    """Q1 / lead decision 1: membership over eligible genes, never over detected ones."""
    u0 = ora().tu_universe(t6_frame(detected=[]))
    u1 = ora().tu_universe(t6_frame())
    s = frozenset({"g1", "g3", "g5", "g8"})
    assert ora().tu_members(u0, s) == ora().tu_members(u1, s) == frozenset({0, 1, 2, 4})


# --------------------------------------------------------------------------- column contract


def test_columns_and_constants() -> None:
    out = _run(t6_frame(), [kegg({"x": {"g1", "g4", "g7"}})])
    assert set(out.columns) == ORA_COLUMNS
    assert list(out.columns) == list(ora().COLUMNS)
    assert set(out["method"]) == {"cpb_tu_ora"}
    assert set(out["null_model"]) == {"conditional_poisson_binomial_tu"}
    assert (out["n_perm"] == 0).all()
    assert out["ci_low"].isna().all() and out["ci_high"].isna().all()
    assert all("ci_unestimated" in _flags(f) for f in out["flags"])
    for c in ("n_universe", "n_in_set", "n_detected", "k_detected_in_set", "n_in_set_genes",
              "k_detected_in_set_genes", "n_perm"):  # fmt: skip
        assert pd.api.types.is_integer_dtype(out[c]), c
    for c in ("odds_weight", "p_classic", "p_noelim", "effect", "p_value", "q_value", "coverage"):
        assert pd.api.types.is_float_dtype(out[c]), c
    assert isinstance(out.attrs.get("run_flags"), list)


def test_empty_collections_give_an_empty_frame_with_every_column() -> None:
    out = _run(t6_frame(), [kegg({})])
    assert len(out) == 0
    assert set(out.columns) == ORA_COLUMNS
    out = _run(t6_frame(), [])
    assert len(out) == 0 and set(out.columns) == ORA_COLUMNS


def test_no_calibrated_claim_is_written() -> None:
    """ADR-0029: no column or cell says "calibrated" (only the flag ``ora_uncalibrated`` may)."""
    out = _run(t6_frame(), [kegg({"x": {"g1", "g4", "g7"}})])
    texts = list(out.columns) + [str(v) for c in out.columns for v in out[c]]
    texts += [str(f) for v in out["flags"] for f in _flags(v)] + list(out.attrs["run_flags"])
    for s in texts:
        assert "calibrated" not in s.lower().replace("uncalibrated", ""), s


# -------------------------------------------------------------------------------- examples


def test_example1_t6_partial_membership() -> None:
    out = _run(t6_frame(), [kegg({"x": {"g1", "g4", "g7"}}, coverage=0.7)], size=(3, 500))
    over, under = row(out, "x", "over"), row(out, "x", "under")
    for r in (over, under):
        assert (r["n_universe"], r["n_in_set"], r["n_detected"], r["k_detected_in_set"]) == (
            6, 3, 4, 3)  # fmt: skip
        assert (r["n_in_set_genes"], r["k_detected_in_set_genes"]) == (3, 1)
        assert r["stat"] == 3
        assert r["effect"] == pytest.approx(0.60913953, abs=1e-8)
        assert r["odds_weight"] == pytest.approx(0.46666667 / 0.34285714, rel=1e-7)
        assert r["coverage"] == 0.7
        assert _na(r["parent_group"])
    assert over["p_value"] == pytest.approx(2134593 / 5159698, rel=1e-11)
    assert over["p_noelim"] == over["p_value"]
    assert under["p_value"] == pytest.approx(1.0, rel=1e-12)
    assert over["p_classic"] == pytest.approx(0.89763087, rel=1e-6)
    assert under["p_classic"] == pytest.approx(0.54832563, rel=1e-6)
    assert over["q_value"] == pytest.approx(over["p_value"], rel=1e-12)  # a family of one


def test_example2_tu_packing() -> None:
    tus = {"TU_1": [Gene(f"c{i}", 10, True, p_detect=0.5) for i in range(1, 6)]}
    for i in range(20):
        tus[f"TU_{i + 2}"] = [Gene(f"s{i + 1}", 10 if i < 10 else None, i < 10, p_detect=0.5)]
    out = _run(frame(tus), [kegg({"C": {f"c{i}" for i in range(1, 6)}})])
    r = row(out, "C", "over")
    assert (r["n_universe"], r["n_in_set"], r["n_detected"], r["k_detected_in_set"]) == (
        21, 1, 11, 1)  # fmt: skip
    assert (r["n_in_set_genes"], r["k_detected_in_set_genes"]) == (5, 5)
    assert r["p_value"] == pytest.approx(341 / 351, rel=1e-11)
    assert r["odds_weight"] == 1.0
    assert r["p_classic"] == pytest.approx(13 / 230, rel=1e-9)


def test_example3_size_filter_counts_tus() -> None:
    out = _run(t6_frame(), [kegg({"x": {"g1", "g4", "g7"}})], size=(4, 500))
    assert len(out) == 0
    # 5 genes but 3 TUs: still filtered at 4 (sizes are TUs, Q6)
    out = _run(t6_frame(), [kegg({"y": {"g1", "g2", "g7", "g8", "g9"}})], size=(3, 500))
    assert len(out) == 0
    out = _run(t6_frame(), [kegg({"x": {"g1", "g4", "g7"}})], size=(1, 2))
    assert len(out) == 0


def test_example4_equal_q_is_the_tu_hypergeometric() -> None:
    det = [0, 1, 2, 3, 5, 6, 7, 8]  # k = 4 in C = the first five, N = 8
    out = _run(singletons([0.4] * 20, det), [kegg({"C": {f"g{i}" for i in range(1, 6)}})])
    assert row(out, "C", "over")["p_value"] == pytest.approx(hypergeom(20, 5, 8).sf(3), rel=1e-11)
    assert row(out, "C", "under")["p_value"] == pytest.approx(hypergeom(20, 5, 8).cdf(4), rel=1e-11)
    assert row(out, "C", "over")["p_value"] == pytest.approx(0.05779154, abs=1e-8)
    assert row(out, "C", "under")["p_value"] == pytest.approx(0.99638803, abs=1e-8)


def _bh_families_frame() -> tuple[pd.DataFrame, list[Any]]:
    colls = [
        kegg({"x": {"g1", "g4", "g7"}, "y": {"g3", "g6"}, "z": {"g10", "g2"}}),
        gene_sets("cog", {"K": {"g1", "g3", "g5"}, "L": {"g6", "g10"}}),
    ]
    return t6_frame(), colls


def test_example6_under_rows_and_bh_per_collection_and_direction() -> None:
    g, colls = _bh_families_frame()
    out = _run(g, colls)
    for (coll, direction), sub in out.groupby(["collection", "direction"]):
        np.testing.assert_allclose(
            sub["q_value"].to_numpy(), bh(sub["p_value"].to_numpy()), rtol=1e-15, err_msg=coll
        )
        if direction == "under":
            np.testing.assert_array_equal(sub["p_value"], sub["p_noelim"])
    under = row(out, "x", "under")
    expect = cpb_raw(g, ["TU_1", "TU_3", "TU_5"]).p_under
    assert under["p_value"] == pytest.approx(expect, rel=1e-11)


def test_non_go_collections_have_no_elim_and_no_parent_group() -> None:
    g, colls = _bh_families_frame()
    out = _run(g, colls)
    np.testing.assert_array_equal(out["p_value"], out["p_noelim"])
    assert out["parent_group"].isna().all()


# ------------------------------------------------------------------- elim on TUs (examples 7, 8)

GO_SETS = {
    "A1": {"g1", "g2", "g3"},
    "A2": {"g4"},
    "A": {"g1", "g2", "g3", "g4", "g5"},
    "B": {"g6"},
    "R": {"g1", "g2", "g3", "g4", "g5", "g6"},
}
GO_PARENTS = {"A1": {"A"}, "A2": {"A"}, "A": {"R"}, "B": {"R"}, "R": set()}


def _elim_frame(p_c: float, extra: bool = False) -> pd.DataFrame:
    """Example 7's TUs a = TU_1 {g1, g2}, b = TU_2 {g3, g4}, c = TU_3 {g5}, d = TU_4 {g6}, plus
    200 unannotated singleton TUs (p 0.05, two detected) so that A1 = {a, b}, both detected, is
    elim-significant (p ~ 5e-4 < 0.01). c is detected; ``p_c`` sets whether A on {c} is too.
    ``extra`` adds TU_205 = {h1}, detected with p 0.001 and in no example-7 term."""
    tus = {
        "TU_1": [Gene("g1", 10, True, p_detect=0.03), Gene("g2", None, False, p_detect=0.02)],
        "TU_2": [Gene("g3", None, False, p_detect=0.03), Gene("g4", 10, True, p_detect=0.02)],
        "TU_3": [Gene("g5", 10, True, p_detect=p_c)],
        "TU_4": [Gene("g6", None, False, p_detect=0.05)],
    }
    for i in range(200):
        tus[f"TU_{i + 5}"] = [Gene(f"f{i + 1}", 10 if i < 2 else None, i < 2, p_detect=0.05)]
    if extra:
        tus["TU_205"] = [Gene("h1", 10, True, p_detect=0.001)]
    return frame(tus)


def _go(
    extra: dict[str, set[str]] | None = None, parents: dict[str, set[str]] | None = None
) -> Any:
    par = {**GO_PARENTS, **(parents or {})}
    return gene_sets("go_bp", {**GO_SETS, **(extra or {})}, parents=par)


def test_example7_elim_removes_tus() -> None:
    g = _elim_frame(0.05)
    out = _run(g, [_go()])
    over = {t: row(out, t, "over") for t in GO_SETS}
    a1 = cpb_raw(g, ["TU_1", "TU_2"])
    assert a1.p_over < 0.01  # A1 is elim-significant, visited first
    assert over["A1"]["p_value"] == pytest.approx(a1.p_over, rel=1e-10)
    # A2 is not an ancestor of A1: tested on its full set {b}
    assert over["A2"]["p_value"] == pytest.approx(cpb_raw(g, ["TU_2"]).p_over, rel=1e-10)
    # A loses TUs {a, b}: tested on {c}; its p_noelim keeps {a, b, c}
    a_elim = cpb_raw(g, ["TU_3"])
    assert a_elim.p_over >= 0.01
    assert over["A"]["p_value"] == pytest.approx(a_elim.p_over, rel=1e-10)
    assert over["A"]["p_noelim"] == pytest.approx(
        cpb_raw(g, ["TU_1", "TU_2", "TU_3"]).p_over, rel=1e-10
    )
    assert over["A"]["effect"] == pytest.approx(1 - a_elim.mean, rel=1e-10)
    assert over["A"]["n_in_set"] == 3 and over["A"]["k_detected_in_set"] == 3
    assert over["A"]["stat"] == 3
    # the rejected gene-removal design would test A on {g4, g5} -> TUs {b, c}
    gene_removal = cpb_raw(g, ["TU_2", "TU_3"]).p_over
    assert not math.isclose(over["A"]["p_value"], gene_removal, rel_tol=1e-3)
    # R loses {a, b} (A was not significant on {c}, so c stays): tested on {c, d}
    assert over["R"]["p_value"] == pytest.approx(cpb_raw(g, ["TU_3", "TU_4"]).p_over, rel=1e-10)
    assert over["R"]["p_noelim"] == pytest.approx(
        cpb_raw(g, ["TU_1", "TU_2", "TU_3", "TU_4"]).p_over, rel=1e-10
    )


def test_elim_is_over_only_under_rows_use_the_full_set() -> None:
    g = _elim_frame(0.05)
    out = _run(g, [_go()])
    for t in GO_SETS:
        r = row(out, t, "under")
        assert r["p_value"] == r["p_noelim"]
    assert row(out, "A", "under")["p_value"] == pytest.approx(
        cpb_raw(g, ["TU_1", "TU_2", "TU_3"]).p_under, rel=1e-10
    )


def test_p_classic_is_gene_level_wallenius_on_the_full_set() -> None:
    """``p_classic`` is never elim'd: A's over row uses all five genes (k = 3 of them: g1, g4,
    g5) with ω = mean p in A / mean p outside, as M8-01's ``wallenius_ora``."""
    from giap.stats.wallenius import wallenius_ora

    g = _elim_frame(0.05)
    out = _run(g, [_go()])
    r = row(out, "A", "over")
    tags = g["locus_tag"].to_numpy()
    in_a = np.isin(tags, sorted(GO_SETS["A"]))
    w = wallenius_ora(g["detected"].to_numpy(bool), g["p_detect"].to_numpy(float), in_a)
    assert (r["n_in_set_genes"], r["k_detected_in_set_genes"]) == (5, 3)
    assert r["odds_weight"] == pytest.approx(w.odds, rel=1e-12)
    assert r["p_classic"] == pytest.approx(w.p_over, rel=1e-9)
    assert row(out, "A", "under")["p_classic"] == pytest.approx(w.p_under, rel=1e-9)


def test_elim_size_filtered_children_remove_nothing() -> None:
    """With size [3, 500] A1 (2 TUs), A2 and B (1 TU) write no row and do not enter the walk,
    so A is tested on its full set (resolution: filtered terms are pruned before elim)."""
    g = _elim_frame(0.05)
    out = _run(g, [_go()], size=(3, 500))
    assert set(out["term"]) == {"A", "R"}
    a = row(out, "A", "over")
    assert a["p_value"] == a["p_noelim"]
    assert a["p_value"] == pytest.approx(cpb_raw(g, ["TU_1", "TU_2", "TU_3"]).p_over, rel=1e-10)


def test_elim_threshold_is_strict_and_on_the_elim_p() -> None:
    """A child whose p is 0.01 or more removes nothing: in example 7, A (p on {c} >= 0.01) does
    not remove c from R."""
    g = _elim_frame(0.05)
    out = _run(g, [_go()])
    assert row(out, "A", "over")["p_value"] >= 0.01
    assert row(out, "R", "over")["p_value"] == pytest.approx(
        cpb_raw(g, ["TU_3", "TU_4"]).p_over, rel=1e-10
    )


def test_example8_parent_group() -> None:
    """q_c = 0.001 makes A significant on {c} as well; R (tested on {d}) is not."""
    g = _elim_frame(0.001)
    out = _run(g, [_go()])
    over = {t: row(out, t, "over") for t in GO_SETS}
    assert over["A1"]["q_value"] <= ALPHA and over["A"]["q_value"] <= ALPHA
    assert over["R"]["q_value"] > ALPHA and over["B"]["q_value"] > ALPHA
    assert over["R"]["p_value"] == pytest.approx(cpb_raw(g, ["TU_4"]).p_over, rel=1e-10)
    assert over["A1"]["parent_group"] == "A"
    assert over["A"]["parent_group"] == "A"
    assert over["A2"]["q_value"] <= ALPHA and over["A2"]["parent_group"] == "A"
    assert _na(over["R"]["parent_group"]) and _na(over["B"]["parent_group"])
    for t in GO_SETS:  # no under row is significant here
        r = row(out, t, "under")
        assert r["q_value"] > ALPHA and _na(r["parent_group"])


def test_parent_group_is_the_most_general_significant_ancestor() -> None:
    """A gets a second parent S (no parents, depth 0) whose own TU h1 (q 0.001, detected) keeps
    S significant after elim removes A's TUs. A1's significant ancestors are A (depth 1) and S
    (depth 0): its group is S, as is A's and S's own; R stays not significant."""
    g = _elim_frame(0.001, extra=True)
    sets = {"S": GO_SETS["A"] | {"h1"}}
    out = _run(g, [_go(sets, {"A": {"R", "S"}, "S": set()})])
    s = row(out, "S", "over")
    assert s["q_value"] <= ALPHA
    assert s["p_value"] == pytest.approx(cpb_raw(g, ["TU_205"]).p_over, rel=1e-10)
    assert row(out, "A", "over")["q_value"] <= ALPHA
    assert row(out, "R", "over")["q_value"] > ALPHA
    assert row(out, "A1", "over")["parent_group"] == "S"
    assert row(out, "A", "over")["parent_group"] == "S"
    assert s["parent_group"] == "S"


def test_elim_empty() -> None:
    """P = {g1, g4} covers TUs {a, b} only and is a parent of A1: elim empties it."""
    g = _elim_frame(0.05)
    out = _run(g, [_go({"P": {"g1", "g4"}}, {"A1": {"A", "P"}, "P": set()})])
    p = row(out, "P", "over")
    assert _na(p["p_value"]) and _na(p["q_value"]) and _na(p["effect"])
    assert {"unestimated", "elim_empty"} <= set(_flags(p["flags"]))
    assert p["p_noelim"] == pytest.approx(cpb_raw(g, ["TU_1", "TU_2"]).p_over, rel=1e-10)
    assert math.isfinite(p["p_classic"])
    assert p["n_in_set"] == 2
    fam = out[(out["collection"] == "go_bp") & (out["direction"] == "over")]
    ok = fam["p_value"].notna()
    np.testing.assert_allclose(fam.loc[ok, "q_value"], bh(fam.loc[ok, "p_value"].to_numpy()))


def test_go_term_without_significant_descendants_keeps_p_noelim() -> None:
    g = _elim_frame(0.05)
    out = _run(g, [_go()])
    for t in ("A1", "A2", "B"):  # leaves
        r = row(out, t, "over")
        assert r["p_value"] == r["p_noelim"]


# ------------------------------------------------------------------------- degenerate rows


def test_example9_c_covers_every_tu() -> None:
    g = t6_frame()
    every = {x for gs in T6.values() for x, _ in gs}
    out = _run(g, [kegg({"all": every, "x": {"g1", "g4", "g7"}, "y": {"g3", "g6"}})])
    for direction in ("over", "under"):
        r = row(out, "all", direction)
        assert _na(r["p_value"]) and _na(r["p_noelim"]) and _na(r["q_value"])
        assert _na(r["odds_weight"]) and _na(r["effect"])
        assert "unestimated" in _flags(r["flags"])
        assert r["n_in_set"] == 6
        rest = out[(out["direction"] == direction) & (out["term"] != "all")]
        np.testing.assert_allclose(rest["q_value"], bh(rest["p_value"].to_numpy()), rtol=1e-15)


def test_no_detected_tu() -> None:
    out = _run(t6_frame(detected=[]), [kegg({"x": {"g1", "g4", "g7"}, "y": {"g3", "g6"}})])
    assert len(out) == 4
    np.testing.assert_allclose(out["p_value"], 1.0, rtol=1e-12)
    np.testing.assert_allclose(out["q_value"], 1.0, rtol=1e-12)
    assert (out["n_detected"] == 0).all() and (out["k_detected_in_set"] == 0).all()
    np.testing.assert_allclose(out["effect"], 0.0, atol=1e-12)


def test_every_tu_detected() -> None:
    every = [x for gs in T6.values() for x, _ in gs]
    out = _run(t6_frame(detected=every), [kegg({"x": {"g1", "g4", "g7"}})])
    np.testing.assert_allclose(out["p_value"], 1.0, rtol=1e-12)
    assert (out["n_detected"] == 6).all() and (out["k_detected_in_set"] == 3).all()
    np.testing.assert_allclose(out["effect"], 0.0, atol=1e-12)


@pytest.mark.parametrize(("gene", "p", "detected"), [("g10", 1.0, False), ("g3", 0.0, True)])
def test_p_detect_of_exactly_0_or_1_is_clipped_and_flagged(gene: str, p: float, detected: bool
                                                           ) -> None:  # fmt: skip
    """Decision 12: a certain-looking gene is clipped to [1e-12, 1 - 1e-12]; the run goes on, the
    rows whose TU set holds it carry ``p_detect_clipped`` and the others do not."""
    g = t6_frame()
    g.loc[g["locus_tag"] == gene, "p_detect"] = p
    g.loc[g["locus_tag"] == gene, "detected"] = detected
    sets = {"has": {gene, "g1", "g4"}, "not": {"g6", "g7"}}
    out = _run(g, [kegg(sets)])
    assert "p_detect_clipped" in out.attrs["run_flags"]
    for direction in ("over", "under"):
        has, other = row(out, "has", direction), row(out, "not", direction)
        assert "p_detect_clipped" in _flags(has["flags"])
        assert "p_detect_clipped" not in _flags(other["flags"])
        assert math.isfinite(has["p_value"]) and math.isfinite(other["p_value"])
        assert "null_degenerate" not in _flags(has["flags"])
    u = ora().tu_universe(g)
    i = u.gene_tu[gene]
    assert u.clipped[i] and int(u.clipped.sum()) == 1
    pc = min(max(p, 1e-12), 1 - 1e-12)
    others = [pp for tu in T6.values() for x, pp in tu if u.gene_tu[x] == i and x != gene]
    expect = math.log1p(-pc) + sum(math.log1p(-x) for x in others)
    assert u.log_comp[i] == pytest.approx(expect, rel=1e-13)
    assert u.q[i] > 0.0 and u.comp[i] > 0.0


def test_detectability_constant_is_not_hypergeometric_on_tus() -> None:
    """Constant p̂ = 0.3: ω = 1 so ``p_classic`` is the gene hypergeometric, but q grows with TU
    size and the TU null is the CPB, not the TU hypergeometric."""
    tus: dict[str, list[Gene]] = {}
    n = 0
    for i in range(30):
        size = 4 if i < 6 else 1
        tus[f"TU_{i + 1}"] = [
            Gene(f"g{n + j + 1}", 10 if (j == 0 and i % 3 == 0) else None, j == 0 and i % 3 == 0,
                 p_detect=0.3)
            for j in range(size)
        ]  # fmt: skip
        n += size
    g = frame(tus)
    big = {x.tag for i in range(6) for x in tus[f"TU_{i + 1}"]}
    out = _run(g, [kegg({"big": big})])
    r = row(out, "big", "over")
    assert r["odds_weight"] == 1.0
    in_big = g["locus_tag"].isin(big)
    M, nn, N, k = len(g), len(big), int(g["detected"].sum()), int(g.loc[in_big, "detected"].sum())
    assert r["p_classic"] == pytest.approx(hypergeom(M, nn, N).sf(k - 1), rel=1e-9)
    tu_big = [f"TU_{i + 1}" for i in range(6)]
    assert r["p_value"] == pytest.approx(cpb_raw(g, tu_big).p_over, rel=1e-11)
    hyp_tu = hypergeom(30, 6, int(r["n_detected"])).sf(int(r["k_detected_in_set"]) - 1)
    assert abs(r["p_value"] - hyp_tu) > 1e-3


def test_odds_out_of_range_drops_p_classic_only() -> None:
    p = [0.5] * 5 + [1e-7] * 15  # ω = 5e6 > 1e6 (ODDS_RANGE)
    out = _run(singletons(p, [0, 1, 7]), [kegg({"C": {f"g{i}" for i in range(1, 6)}})])
    for direction in ("over", "under"):
        r = row(out, "C", direction)
        assert _na(r["p_classic"]) and _na(r["odds_weight"])
        assert "odds_out_of_range" in _flags(r["flags"])
        assert math.isfinite(r["p_value"]) and math.isfinite(r["q_value"])
        assert "unestimated" not in _flags(r["flags"])


def test_odds_undefined_keeps_the_cpb_columns() -> None:
    """Every gene outside C has p̂ 0 (clipped for q): ω is undefined, the CPB p stands."""
    p = [0.5] * 5 + [0.0] * 15
    out = _run(singletons(p, [0, 1]), [kegg({"C": {f"g{i}" for i in range(1, 6)}})])
    r = row(out, "C", "over")
    assert _na(r["p_classic"]) and _na(r["odds_weight"])
    assert math.isfinite(r["p_value"])


def test_example10_deterministic_and_row_order_free() -> None:
    g, colls = _bh_families_frame()
    a = _run(g, colls)
    b = _run(g, colls)
    pd.testing.assert_frame_equal(a, b)
    c = _run(g.sample(frac=1.0, random_state=3).reset_index(drop=True), colls)
    pd.testing.assert_frame_equal(a, c)


def test_rows_are_sorted_by_collection_term_direction() -> None:
    g = _elim_frame(0.05)
    colls = [_go(), kegg({"k2": {"g1", "g6"}, "k1": {"g3", "g5"}})]
    out = _run(g, colls)
    seen: list[str] = []
    for c in out["collection"]:
        if not seen or seen[-1] != c:
            assert c not in seen, "collections are not contiguous"
            seen.append(c)
    for _, sub in out.groupby("collection", sort=False):
        keys = list(zip(sub["term"], sub["direction"], strict=True))
        assert keys == sorted(keys)
    assert set(out["direction"]) == {"over", "under"}


def test_effect_is_k_minus_conditional_mean_in_tus() -> None:
    g, colls = _bh_families_frame()
    out = _run(g, colls)
    r = row(out, "y", "over")
    law = cpb_raw(g, ["TU_2", "TU_4"])
    assert r["effect"] == pytest.approx(r["k_detected_in_set"] - law.mean, rel=1e-11)
    assert row(out, "y", "under")["effect"] == r["effect"]
