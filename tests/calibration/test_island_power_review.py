"""PR #188 review (code reviewer B1, B2; statistician nit 1): the power recorder's planted truth,
its one stream per replicate, and the region-free power layout."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from test_island_power import power_module

import gippasim.layout as la


def _no_scan(m: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Skip the M5-06 scan: the truth and the streams are read before it."""
    monkeypatch.setattr(m.idat, "build_series", lambda *a, **k: None)
    empty = pd.DataFrame({"direction": [], "genes": []})
    monkeypatch.setattr(m.idat, "step", lambda *a, **k: (empty, None))


@pytest.mark.parametrize("seed", [3, 11])
def test_planted_set_is_the_island_tus_of_that_replicates_layout(
    seed: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """B1: planted = chr1 TUs first_tu .. first_tu + n_tus of the replicate's own layout."""
    m = power_module()
    _no_scan(m, monkeypatch)
    pair = m.scenario(20, 2.0, "high")
    (isl,) = pair[1].islands
    _, planted = m._replicate(np.random.default_rng(seed), pair)
    lay = la.simulate_layout(pair[0], np.random.default_rng(seed))
    chr_tus = lay.tus[lay.tus["contig"] == isl.replicon].sort_values("start", kind="mergesort")
    want = frozenset(chr_tus["tu_id"].iloc[isl.first_tu : isl.first_tu + isl.n_tus].astype(str))
    assert planted == want
    assert len(planted) == 20


def test_every_replicate_draws_its_own_stream(monkeypatch: pytest.MonkeyPatch) -> None:
    """B2: replicate r uses stream ``<STREAM>.r<r>``; R distinct names, R distinct layouts."""
    m = power_module()
    names: list[str] = []
    real = m.for_step

    def spy(name: str, seed: int) -> Any:
        names.append(name)
        return real(name, seed)

    tus_seen: list[frozenset[str]] = []

    def rep(rng: Any, pair: Any) -> tuple[list[Any], frozenset[str]]:
        lay = la.simulate_layout(pair[0], rng)
        tus_seen.append(frozenset(map(str, lay.tus["start"].astype(str) + lay.tus["contig"])))
        return [], frozenset()

    monkeypatch.setattr(m, "for_step", spy)
    monkeypatch.setattr(m, "_replicate", rep)
    m.run_cell(10, 1.5, "high", 4, m.cell_seed(0))
    assert names == [f"{m.STREAM}.r{r}" for r in range(4)]
    assert len(set(tus_seen)) == 4


def test_power_layout_has_no_regions() -> None:
    """Statistician nit 1: the power scenario drops regions (it carries no region effect), so
    the planted island cannot overlap one."""
    m = power_module()
    lcfg, ecfg = m.scenario(20, 2.0, "high")
    assert tuple(lcfg.regions) == ()
    assert not getattr(ecfg, "region_effect_log2", None)
    assert not getattr(ecfg, "region_detect_shift", None)
