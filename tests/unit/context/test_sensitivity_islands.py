"""M5-10 ``islands`` sensitivity family: examples 1-5, edge cases, the flagged-set rule and the
signature of lead decision 38 (docs/specs/M5-10-sensitivity-islands.md; GIPPA_SPEC.md §8.6, G9,
§7.3.5(d); M4-19 lead decisions 2 and 5).

Example 4 (``output.sensitivity = off``) and the step-level B = 199 / stream checks live in
tests/integration/test_sensitivity_islands_step.py. The reference flagged set comes from the
merged M5-06 producers (``tests/helpers/sensitivity_islands_data.oracle_flagged``); examples 1
and 3 stub ``islands_frame`` so the ID rule is pinned on chosen island rows.
"""

from __future__ import annotations

import inspect
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import sensitivity_islands_data as sid


def _S() -> Any:
    from gippa.context import sensitivity as S

    return S


def _rng() -> np.random.Generator:
    from gippa.rng import for_step

    return for_step(sid.STREAM, 20261022)


def _is_na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


def _row(primary: frozenset[str], variant: frozenset[str], name: str = "missing=min") -> Any:
    df = _S().sensitivity_rows({sid.FAMILY: primary}, {(sid.FAMILY, name): variant})
    assert len(df) == 1
    return df.iloc[0]


# --------------------------------------------------------------------------- interface


def test_islands_family_is_registered() -> None:
    S = _S()
    assert "regions" in S.FAMILIES, "M4-19's family must stay registered"
    assert S.FAMILIES["islands"] is S.islands_family
    assert "islands_family" in S.__all__


def test_islands_family_signature_decision_38() -> None:
    """Q1 / decision 38: ``(genes, cfg, rng)`` positional, ``tus`` and ``replicons``
    keyword-only and required, so the step binds them by name (M4-19 decision 5)."""
    params = inspect.signature(_S().islands_family).parameters
    assert list(params)[:3] == ["genes", "cfg", "rng"]
    for name in ("tus", "replicons"):
        assert params[name].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert params[name].default is inspect.Parameter.empty, name
    assert {"regions", "region_genes"}.isdisjoint(params), "islands do not read C7 / C8"


@pytest.mark.parametrize("given", [(), ("tus",), ("replicons",)])
def test_islands_family_requires_both_tables(given: tuple[str, ...]) -> None:
    t = sid.unit_tables()
    tables = {"tus": t.tus, "replicons": t.replicons}
    with pytest.raises(TypeError):
        _S().islands_family(t.genes, sid.unit_config(), _rng(), **{k: tables[k] for k in given})


# --------------------------------------------------------------------------- flagged-set rule


def test_flagged_set_is_every_gene_of_every_island() -> None:
    """The family's set is ``<direction>:<locus_tag>`` over the best-window gene list of every
    M5-06 island, on the generator it is given (same state -> the M5-06 pipeline's islands)."""
    t = sid.unit_tables()
    cfg = sid.unit_config()
    want = sid.oracle_flagged(t.genes, t.tus, t.replicons, cfg, _rng())
    assert want, "fixture: the planted island is not found"
    got = sid.call_islands_family(t, cfg, _rng())
    assert isinstance(got, frozenset)
    assert all(isinstance(x, str) for x in got)
    assert got == want
    assert all(x.split(":", 1)[0] in ("high", "low") for x in got)


def test_best_window_genes_not_the_cluster() -> None:
    """Edge case 2 / M5-06 Q3: the set holds exactly the reported ``genes`` of each island row
    (its best window), so its size is the sum of the rows' gene counts (no duplicates here)."""
    t = sid.unit_tables()
    cfg = sid.unit_config()
    frame = sid.oracle_frame(t.genes, t.tus, t.replicons, cfg, _rng())
    got = sid.call_islands_family(t, cfg, _rng())
    assert len(got) == int(frame["n_genes"].sum())
    tags = {x.split(":", 1)[1] for x in got}
    assert tags == {str(g) for gs in frame["genes"] for g in gs}


def test_ids_keep_direction(monkeypatch: pytest.MonkeyPatch) -> None:
    """A gene in a high and in a low island gives two IDs (the gene -> island map of step 9
    would keep only one; the family uses every island, not the map)."""
    frame = sid.fake_islands([("high", ["G1", "G5"]), ("low", ["G5", "G9"])])
    sid.patch_islands_frame(monkeypatch, [frame])
    got = sid.call_islands_family(sid.tiny_tables(), sid.unit_config(), _rng())
    assert got == frozenset({"high:G1", "high:G5", "low:G5", "low:G9"})


def test_no_island_gives_the_empty_set(monkeypatch: pytest.MonkeyPatch) -> None:
    sid.patch_islands_frame(monkeypatch, [sid.fake_islands([])])
    assert sid.call_islands_family(sid.tiny_tables(), sid.unit_config(), _rng()) == frozenset()


def test_unscanned_replicon_gives_the_empty_set() -> None:
    """M5-06 decision 40: a replicon with < 2 ``min_tus`` TUs is not scanned -> no island."""
    t = sid.unscanned_tables()
    assert sid.call_islands_family(t, sid.unit_config(), _rng()) == frozenset()


def test_family_passes_genes_cfg_and_rng_to_islands_frame(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The stubbed ``islands_frame`` is called once per family call with the gene table, a
    config carrying the caller's island knobs, and the caller's generator."""
    calls = sid.patch_islands_frame(monkeypatch, [sid.fake_islands([("high", ["G1"])])])
    t = sid.tiny_tables()
    cfg = sid.unit_config(alpha=0.01)
    rng = _rng()
    sid.call_islands_family(t, cfg, rng)
    assert len(calls) == 1
    args = list(calls[0]["args"]) + list(calls[0]["kwargs"].values())
    assert any(a is rng for a in args), "islands_frame must draw from the family's generator"
    seen = [a for a in args if hasattr(a, "context") and hasattr(a, "stats")]
    assert len(seen) == 1
    assert seen[0].stats.alpha.context == 0.01
    assert seen[0].context.islands == cfg.context.islands
    frames = [a for a in args if isinstance(a, pd.DataFrame)]
    assert any(set(f["locus_tag"]) == set(t.genes["locus_tag"]) for f in frames)


# --------------------------------------------------------------------------- examples 1-3


def test_example_1(monkeypatch: pytest.MonkeyPatch) -> None:
    """primary {high:G1,G2,G3}, ``missing=min`` {high:G2,G3,G4} -> 3, 3, 2, 0.5."""
    frames = [sid.fake_islands([("high", ["G1", "G2", "G3"])]),
              sid.fake_islands([("high", ["G2", "G3", "G4"])])]  # fmt: skip
    sid.patch_islands_frame(monkeypatch, frames)
    t = sid.tiny_tables()
    cfg = sid.unit_config()
    p = sid.call_islands_family(t, cfg, _rng())
    v = sid.call_islands_family(t, cfg, _rng())
    assert p == frozenset({"high:G1", "high:G2", "high:G3"})
    r = _row(p, v)
    assert (r["family"], r["variant"]) == ("islands", "missing=min")
    assert (int(r["n_primary"]), int(r["n_variant"]), int(r["n_shared"])) == (3, 3, 2)
    assert float(r["jaccard"]) == pytest.approx(0.5)
    assert _is_na(r["agreement"])


def test_example_2_no_island_in_either(monkeypatch: pytest.MonkeyPatch) -> None:
    sid.patch_islands_frame(monkeypatch, [sid.fake_islands([])])
    t = sid.tiny_tables()
    p = sid.call_islands_family(t, sid.unit_config(), _rng())
    v = sid.call_islands_family(t, sid.unit_config(), _rng())
    r = _row(p, v)
    assert (int(r["n_primary"]), int(r["n_variant"]), int(r["n_shared"])) == (0, 0, 0)
    assert _is_na(r["jaccard"])


def test_example_3_direction_change_is_a_change(monkeypatch: pytest.MonkeyPatch) -> None:
    """primary {low:G7}, variant {high:G7} -> n_shared 0, jaccard 0.0."""
    sid.patch_islands_frame(monkeypatch, [sid.fake_islands([("low", ["G7"])]),
                                          sid.fake_islands([("high", ["G7"])])])  # fmt: skip
    t = sid.tiny_tables()
    p = sid.call_islands_family(t, sid.unit_config(), _rng())
    v = sid.call_islands_family(t, sid.unit_config(), _rng())
    assert (p, v) == (frozenset({"low:G7"}), frozenset({"high:G7"}))
    r = _row(p, v)
    assert (int(r["n_primary"]), int(r["n_variant"]), int(r["n_shared"])) == (1, 1, 0)
    assert float(r["jaccard"]) == 0.0


# --------------------------------------------------------------------------- example 5, edges


def _without_abundance(genes: pd.DataFrame, k: int) -> tuple[pd.DataFrame, list[str]]:
    """``genes`` with ``abundance_log2`` NaN on ``k`` detected chr1 genes (what ``exclude``
    leaves for genes the primary imputed); detection unchanged."""
    det = genes[(genes["contig"] == "chr1") & genes["detected"].astype(bool)]
    tags = det["locus_tag"].astype(str).tolist()[::7][:k]
    out = genes.copy()
    out.loc[out["locus_tag"].isin(tags), "abundance_log2"] = np.nan
    return out, tags


def test_example_5_exclude_drops_unscored_genes_and_still_yields_a_set() -> None:
    """Detected genes without an abundance are not scored (M5-06 step 1); the family runs on
    such a table and returns the M5-06 islands of it."""
    t = sid.unit_tables()
    genes, tags = _without_abundance(t.genes, 20)
    assert len(tags) == 20
    cfg = sid.unit_config(abundance_missing="exclude")
    got = sid.call_islands_family(sid.IslandTables(genes, t.tus, t.replicons), cfg, _rng())
    assert got == sid.oracle_flagged(genes, t.tus, t.replicons, cfg, _rng())


def test_exclude_with_nothing_missing_equals_primary() -> None:
    """Edge case 1: with every detected gene observed, ``exclude`` = primary (jaccard 1 or NA)."""
    t = sid.unit_tables()
    p = sid.call_islands_family(t, sid.unit_config(), _rng())
    v = sid.call_islands_family(t, sid.unit_config(abundance_missing="exclude"), _rng())
    assert p == v
    r = _row(p, v, "missing=exclude")
    assert _is_na(r["jaccard"]) or float(r["jaccard"]) == 1.0


# --------------------------------------------------------------------------- determinism, purity


def test_same_generator_state_same_set() -> None:
    t = sid.unit_tables()
    cfg = sid.unit_config()
    a = sid.call_islands_family(t, cfg, _rng())
    np.random.seed(1)  # noqa: NPY002 - the legacy global state must not matter
    b = sid.call_islands_family(t, cfg, _rng())
    np.random.seed(2)  # noqa: NPY002
    c = sid.call_islands_family(t, cfg, _rng())
    assert a == b == c


def test_family_draws_only_from_the_given_generator() -> None:
    """The generator advances: the scan's bootstrap draws come from it."""
    t = sid.unit_tables()
    rng = _rng()
    before = rng.bit_generator.state
    sid.call_islands_family(t, sid.unit_config(), rng)
    assert rng.bit_generator.state != before


def test_inputs_are_not_mutated() -> None:
    t = sid.unit_tables()
    snap = (t.genes.copy(deep=True), t.tus.copy(deep=True), t.replicons.copy(deep=True))
    sid.call_islands_family(t, sid.unit_config(), _rng())
    pd.testing.assert_frame_equal(t.genes, snap[0])
    pd.testing.assert_frame_equal(t.tus, snap[1])
    pd.testing.assert_frame_equal(t.replicons, snap[2])
