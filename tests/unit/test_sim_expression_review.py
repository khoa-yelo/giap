"""giapsim expression: lead decisions 5 and 6 after the PR #71 code review (M4-08).

docs/specs/M4-08-sim-expression.md, "Lead decisions after code review":

5. Tests kill the surviving mutants: islands at rho = 0.5 shift by exactly ``effect`` times the
   stationary SD ``tu_sd / sqrt(1 - rho^2)`` (decision 4); an island ending on the last TU of its
   replicon is accepted, one TU later raises; TUs straddling a region boundary keep ``tu_log2``
   unchanged (``dt == 0``) while each gene gets its own region effect; ``expression_for_simulation``
   equals ``abundance_log2_true * ln 2`` of ``simulate_expression`` with the same rng.
6. The gene table has ``island_id``: NA outside islands, the island's id (as in the islands truth
   frame) inside; several islands are distinguishable.

Expected values come from tests/helpers/expression_oracle.py (written from the spec).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo


def _m() -> Any:
    import giapsim.expression

    return giapsim.expression


def layout(seed: int = 1, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def express(lay: Any, seed: int = 1, **kw: Any) -> Any:
    m = _m()
    return m.simulate_expression(lay, m.ExpressionConfig(**kw), np.random.default_rng(seed))


def island(*args: Any) -> Any:
    return _m().IslandSpec(*args)


def by_tag(expr: Any) -> pd.DataFrame:
    return expr.genes.set_index("locus_tag")


def island_genes(lay: Any, contig: str, first: int, n: int) -> set[str]:
    mem = eo.tu_members(lay)
    return {g for t in eo.tu_order(lay)[contig][first : first + n] for g in mem[t]}


@pytest.fixture(scope="module")
def lay1() -> Any:
    return layout(1)


# ----------------------------------------------------------------------------- decision 5


@pytest.mark.parametrize(("rho", "effect"), [(0.5, 2.0), (0.5, -1.5), (0.9, 1.0)])
def test_island_shift_is_effect_times_stationary_sd(lay1: Any, rho: float, effect: float) -> None:
    """Paired runs: island genes move by exactly effect * tu_sd / sqrt(1 - rho^2), others by 0."""
    tu_sd = 2.0
    want = effect * tu_sd / math.sqrt(1.0 - rho * rho)
    base = by_tag(express(lay1, 5, ar1_rho=rho, tu_sd=tu_sd))
    e = express(lay1, 5, ar1_rho=rho, tu_sd=tu_sd, islands=(island("chr1", 30, 10, effect),))
    eff = by_tag(e)
    inside = island_genes(lay1, "chr1", 30, 10)
    for g in eff.index:
        w = want if g in inside else 0.0
        assert eff.at[g, "tu_log2"] - base.at[g, "tu_log2"] == pytest.approx(w, abs=1e-9), g
        assert eff.at[g, "abundance_log2_true"] - base.at[g, "abundance_log2_true"] == (
            pytest.approx(w, abs=1e-9)
        ), g
    assert float(e.islands["effect_log2"].iloc[0]) == pytest.approx(want, abs=1e-12)


def test_island_ending_on_last_tu_accepted(lay1: Any) -> None:
    n_tu = len(eo.tu_order(lay1)["chr1"])
    e = express(lay1, 1, islands=(island("chr1", n_tu - 5, 5, 1.0),))
    g = by_tag(e)
    inside = island_genes(lay1, "chr1", n_tu - 5, 5)
    assert inside
    assert set(g.index[g["in_island"].astype(bool)]) == inside


@pytest.mark.parametrize(("first_back", "n"), [(4, 5), (5, 6), (0, 1)])
def test_island_one_tu_past_the_end_raises(lay1: Any, first_back: int, n: int) -> None:
    """first_tu + n_tus == number of TUs + 1: the island runs past the replicon's last TU."""
    n_tu = len(eo.tu_order(lay1)["chr1"])
    first = n_tu - first_back
    assert first + n == n_tu + 1
    with pytest.raises(ValueError):
        express(lay1, 1, islands=(island("chr1", first, n, 1.0),))


REGION_TYPES = ("prophage", "bgc", "amr", "crispr")


@pytest.mark.parametrize("rtype", REGION_TYPES)
def test_straddling_tus_keep_tu_level_genes_get_own_effect(lay1: Any, rtype: str) -> None:
    """TUs with some (not all) members in the region: dt == 0; each gene's da = its own effect."""
    status = eo.tu_region_status(lay1, rtype)
    if not any(s == "some" for s in status.values()):
        pytest.skip(f"no TU straddles a {rtype} boundary in the default layout")
    rt = eo.region_types_of_genes(lay1)
    mem = eo.tu_members(lay1)
    base = by_tag(express(lay1, 1))
    eff = by_tag(express(lay1, 1, region_effect_log2={rtype: -3.0}))
    n_some = 0
    for t, s in status.items():
        if s != "some":
            continue
        for g in mem[t]:
            n_some += 1
            own = -3.0 if rtype in rt.get(g, set()) else 0.0
            dt = eff.at[g, "tu_log2"] - base.at[g, "tu_log2"]
            da = eff.at[g, "abundance_log2_true"] - base.at[g, "abundance_log2_true"]
            assert dt == pytest.approx(0.0, abs=1e-9), (t, g)
            assert da == pytest.approx(own, abs=1e-9), (t, g)
            assert eff.at[g, "region_effect_log2"] == pytest.approx(own), (t, g)
    assert n_some >= 2


def test_prophage_example3_has_a_straddling_tu(lay1: Any) -> None:
    """Example 3's layout really exercises the 'some' branch (else the test above is vacuous)."""
    assert "some" in set(eo.tu_region_status(lay1, "prophage").values())


@pytest.mark.parametrize("seed", [1, 7])
def test_expression_for_simulation_is_abundance_times_ln2(seed: int) -> None:
    """Same rng: the hook returns ln-scale abundance = abundance_log2_true * ln 2 (no regions)."""
    m = _m()
    lay = layout(2, n_genes_chromosome=300, regions=(), plasmids=(20,))
    c = m.ExpressionConfig(plasmid_copy_log2={}, ar1_rho=0.7)
    ref = m.simulate_expression(lay, c, np.random.default_rng(seed)).genes
    rng = np.random.default_rng(seed)
    got = np.asarray(m.expression_for_simulation(lay.genes, lay.tus, c, rng), dtype=float)
    tags = lay.genes["locus_tag"].astype(str)
    want = ref.set_index("locus_tag").loc[tags, "abundance_log2_true"].to_numpy() * math.log(2.0)
    np.testing.assert_allclose(got, want, rtol=0, atol=1e-12)


# ----------------------------------------------------------------------------- decision 6


def test_island_id_column_na_outside_id_inside(lay1: Any) -> None:
    e = express(lay1, 1, islands=(island("chr1", 12, 7, 1.0),))
    g = by_tag(e)
    assert "island_id" in g.columns
    iid = str(e.islands["island_id"].iloc[0])
    inside = island_genes(lay1, "chr1", 12, 7)
    for tag, v in g["island_id"].items():
        if tag in inside:
            assert not eo.is_na(v) and str(v) == iid, tag
        else:
            assert eo.is_na(v), tag
    assert (g["in_island"].astype(bool) == g["island_id"].notna()).all()


def test_island_id_all_na_without_islands(lay1: Any) -> None:
    g = express(lay1, 1).genes
    assert "island_id" in g.columns
    assert g["island_id"].isna().all()


def test_several_islands_distinguishable(lay1: Any) -> None:
    specs = (island("chr1", 0, 5, 1.0), island("chr1", 40, 20, -2.0), island("chr1", 100, 6, 0.5))
    e = express(lay1, 1, islands=specs)
    g = by_tag(e)
    ids = [str(v) for v in e.islands["island_id"]]
    assert len(set(ids)) == len(specs)
    for iid, s in zip(ids, specs, strict=True):
        inside = island_genes(lay1, s.replicon, s.first_tu, s.n_tus)
        got = {str(t) for t, v in g["island_id"].items() if not eo.is_na(v) and str(v) == iid}
        assert got == inside, iid
    assert set(g["island_id"].dropna().astype(str)) == set(ids)


def test_island_id_on_plasmid_island() -> None:
    """Islands on different replicons keep their own ids."""
    lay = layout(4, n_genes_chromosome=300, plasmids=(40,))
    order = eo.tu_order(lay)
    plasmid = next(c for c in order if c != "chr1")
    specs = (island("chr1", 2, 4, 1.0), island(plasmid, 1, 3, -1.0))
    e = express(lay, 1, islands=specs)
    g = by_tag(e)
    ids = [str(v) for v in e.islands["island_id"]]
    assert ids[0] != ids[1]
    for iid, s in zip(ids, specs, strict=True):
        inside = island_genes(lay, s.replicon, s.first_tu, s.n_tus)
        assert {str(v) for v in g.loc[sorted(inside), "island_id"]} == {iid}
