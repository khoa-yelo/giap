"""M8-04b: ``category_summary`` (docs/specs/M8-04 §3; GIAP_SPEC.md §7.7.2, §9.3; unchanged by
ADR-0029 and ADR-0035).

Interface pinned here (module ``giap.function.summary``):
``category_summary(genes, colls, B, rng) -> DataFrame`` with the §9.3 columns plus ``flags``
(``SUMMARY_COLUMNS``); one row per category of the ``cog``, ``kegg`` and ``go_slim`` collections
(others are ignored), every category whatever its size, sorted by collection in that order, then
category. ``expected`` = Σ p_detect (raw, unclipped), ``o_minus_e`` with a 95 % TU-block bootstrap
CI of the per-TU sums of d - p (M4-16 ``basic_rows``), NA with the flag ``ci_few_units`` below 10
TUs.

Resolutions (reported to the lead):

- the degenerate row "an empty collection: category_summary still reports coverage" is met by
  ``attrs["coverage"]``: {collection: coverage} for every cog / kegg / go_slim collection passed,
  empty or not (a table row needs a category);
- ``flags`` cells are sequences of strings, as in ``detection_enrichment``;
- eligibility and fail-closed rules are ``detection_enrichment``'s (non-partial genes; NaN
  p_detect on one raises);
- on singleton TUs named so that lexicographic and natural order agree, the CI is exactly
  ``tu_block_bootstrap`` of the per-TU d - p sums with the stat "sum" and the caller's generator
  (one category, so no draw-order question arises).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.ora_data import SUMMARY_COLUMNS, Gene, as_list, frame, gene_sets, summary

from giap.rng import for_step
from giap.stats.bootstrap import tu_block_bootstrap

STREAM = "function.bootstrap"
B = 200


def _run(genes: pd.DataFrame, colls: list[Any], seed: int = 1, b: int = B) -> pd.DataFrame:
    return summary().category_summary(genes, colls, b, for_step(STREAM, seed))


def _rev1_frame() -> pd.DataFrame:
    """Revision-1 universe: 10 singleton TUs, p 0.5 except g1, g2 at 0.9; detected g1, g2, g3,
    g5."""
    det = {"g1", "g2", "g3", "g5"}
    p = {"g1": 0.9, "g2": 0.9}
    return frame(
        {
            f"TU_{i}": [Gene(f"g{i}", 10.0 if f"g{i}" in det else None, f"g{i}" in det,
                             p_detect=p.get(f"g{i}", 0.5))]
            for i in range(1, 11)
        }
    )  # fmt: skip


def _row(out: pd.DataFrame, coll: str, cat: str) -> pd.Series:
    r = out[(out["collection"] == coll) & (out["category"] == cat)]
    assert len(r) == 1, (coll, cat, len(r))
    return r.iloc[0]


def test_revision1_cog_k() -> None:
    cog = gene_sets("cog", {"K": {"g1", "g2", "g3", "g4", "g6"}}, coverage=0.5)
    out = _run(_rev1_frame(), [cog])
    assert set(out.columns) == SUMMARY_COLUMNS
    r = _row(out, "cog", "K")
    assert (r["n_genes"], r["n_detected"]) == (5, 3)
    assert r["expected"] == pytest.approx(3.3, abs=1e-12)
    assert r["o_minus_e"] == pytest.approx(-0.3, abs=1e-12)
    assert pd.isna(r["ci_low"]) and pd.isna(r["ci_high"])
    assert "ci_few_units" in as_list(r["flags"])
    assert r["coverage"] == 0.5
    assert r["name"] == "name of K"


def _big_frame(n: int = 40) -> pd.DataFrame:
    rng = np.random.default_rng(5)
    p = rng.uniform(0.3, 0.7, n)
    det = rng.random(n) < p
    return frame(
        {
            f"TU_{i + 1:02d}": [Gene(f"g{i + 1:02d}", 10.0 + i if det[i] else None, bool(det[i]),
                                     p_detect=float(p[i]))]
            for i in range(n)
        }
    )  # fmt: skip


def test_large_category_has_a_tu_bootstrap_ci() -> None:
    g = _big_frame()
    tags = list(g["locus_tag"])
    out = _run(g, [gene_sets("kegg", {"map1": set(tags)})], seed=3)
    r = _row(out, "kegg", "map1")
    d = g["detected"].astype(float).to_numpy()
    p = g["p_detect"].to_numpy()
    assert r["o_minus_e"] == pytest.approx(float((d - p).sum()), abs=1e-12)
    assert "ci_few_units" not in as_list(r["flags"])
    assert r["ci_low"] <= r["o_minus_e"] <= r["ci_high"]
    tus = pd.DataFrame({"tu": g["tu_id"].astype(str), "v": d - p})
    ci = tu_block_bootstrap(tus, "tu", lambda x: float(x["v"].sum()), B, for_step(STREAM, 3))
    assert r["ci_low"] == pytest.approx(ci.ci_low, abs=1e-12)
    assert r["ci_high"] == pytest.approx(ci.ci_high, abs=1e-12)


def test_ci_resamples_tus_not_genes() -> None:
    """Twelve 2-gene TUs whose genes are perfectly anti-correlated (one detected, one not, p 0.5
    each): every TU sum of d - p is 0, so a TU bootstrap CI is exactly [0, 0], while resampling
    genes would spread it."""
    tus = {
        f"TU_{i + 1:02d}": [Gene(f"a{i + 1:02d}", 10.0, True, p_detect=0.5),
                            Gene(f"b{i + 1:02d}", None, False, p_detect=0.5)]
        for i in range(12)
    }  # fmt: skip
    g = frame(tus)
    out = _run(g, [gene_sets("cog", {"J": set(g["locus_tag"])})])
    r = _row(out, "cog", "J")
    assert r["o_minus_e"] == pytest.approx(0.0, abs=1e-12)
    assert r["ci_low"] == pytest.approx(0.0, abs=1e-12)
    assert r["ci_high"] == pytest.approx(0.0, abs=1e-12)


def test_median_abundance_over_detected_members_with_a_value() -> None:
    g = frame(
        {
            "TU_1": [Gene("g1", 10.0, True, p_detect=0.5), Gene("g2", 12.0, True, imputed=True,
                                                                 p_detect=0.5)],
            "TU_2": [Gene("g3", None, True, p_detect=0.5)],
            "TU_3": [Gene("g4", 99.0, False, p_detect=0.5)],
            "TU_4": [Gene("g5", None, False, p_detect=0.5)],
        }
    )  # fmt: skip
    out = _run(g, [gene_sets("cog", {"K": {"g1", "g2", "g3", "g4"}, "L": {"g5"}})])
    assert _row(out, "cog", "K")["median_abundance_log2"] == pytest.approx(11.0)
    assert pd.isna(_row(out, "cog", "L")["median_abundance_log2"])


def test_only_cog_kegg_go_slim_and_sorted() -> None:
    g = _rev1_frame()
    colls = [
        gene_sets("go_slim", {"GO:2": {"g1"}, "GO:1": {"g2", "g3"}}, coverage=0.3),
        gene_sets("go_bp", {"GO:9": {"g1"}}, parents={"GO:9": set()}),
        gene_sets("kegg", {"map2": {"g4"}, "map1": {"g5"}}, coverage=0.2),
        gene_sets("cog", {"K": {"g6"}, "C": {"g7"}}, coverage=0.4),
    ]
    out = _run(g, colls)
    keys = list(zip(out["collection"], out["category"], strict=True))
    assert keys == [("cog", "C"), ("cog", "K"), ("kegg", "map1"), ("kegg", "map2"),
                    ("go_slim", "GO:1"), ("go_slim", "GO:2")]  # fmt: skip
    assert list(out["coverage"]) == [0.4, 0.4, 0.2, 0.2, 0.3, 0.3]
    assert (out["n_genes"] >= 1).all()  # one-gene categories are reported


def test_empty_collection_still_reports_coverage() -> None:
    colls = [gene_sets("cog", {}, coverage=0.0), gene_sets("kegg", {"map1": {"g1"}}, coverage=0.1)]
    out = _run(_rev1_frame(), colls)
    assert set(out["collection"]) == {"kegg"}
    assert out.attrs["coverage"] == {"cog": 0.0, "kegg": 0.1}
    empty = _run(_rev1_frame(), [gene_sets("go_slim", {}, coverage=0.0)])
    assert len(empty) == 0 and set(empty.columns) == SUMMARY_COLUMNS
    assert empty.attrs["coverage"] == {"go_slim": 0.0}


def test_partial_genes_are_not_counted() -> None:
    g = _rev1_frame()
    g.loc[g["locus_tag"] == "g4", "partial"] = True
    g.loc[g["locus_tag"] == "g4", "p_detect"] = math.nan
    out = _run(g, [gene_sets("cog", {"K": {"g1", "g2", "g3", "g4", "g6"}})])
    r = _row(out, "cog", "K")
    assert (r["n_genes"], r["n_detected"]) == (4, 3)
    assert r["expected"] == pytest.approx(2.8, abs=1e-12)


def test_nan_p_detect_on_an_eligible_gene_fails_closed() -> None:
    g = _rev1_frame()
    g.loc[g["locus_tag"] == "g6", "p_detect"] = math.nan
    with pytest.raises(ValueError, match="g6"):
        _run(g, [gene_sets("cog", {"K": {"g1", "g6"}})])


def test_expected_uses_the_unclipped_p_detect() -> None:
    g = _rev1_frame()
    g.loc[g["locus_tag"] == "g6", "p_detect"] = 1.0
    out = _run(g, [gene_sets("cog", {"K": {"g6", "g7"}})])
    assert _row(out, "cog", "K")["expected"] == 1.5


def test_example10_same_seed_same_table() -> None:
    g = _big_frame()
    colls = [gene_sets("kegg", {"map1": set(g["locus_tag"]), "map2": set(g["locus_tag"][:5])})]
    a, b = _run(g, colls, seed=8), _run(g, colls, seed=8)
    pd.testing.assert_frame_equal(a, b)
    c = _run(g, colls, seed=9)
    assert not math.isclose(_row(a, "kegg", "map1")["ci_low"], _row(c, "kegg", "map1")["ci_low"])
