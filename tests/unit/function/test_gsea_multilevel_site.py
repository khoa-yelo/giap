"""M8-11 call site: ``gsea_table`` runs ``fgsea_multilevel`` after ``fgsea_simple`` (spec item 6;
docs/specs/M8-11-stats-fgsea-multilevel.md; M8-index decisions 12-13).

Pinned here:

- ``giap.function.gsea`` imports ``fgsea_multilevel`` by name (tests replace that module attribute)
  and calls it as ``fgsea_multilevel(prep, simple, cfg.sample_size, cfg.eps, rng)`` right after
  each main-run ``fgsea_simple`` of the same ``prep``;
- ``method`` is ``fgsea_multilevel`` on every row of the run, also when no set needs the stage
  (one method per run, decision Q2);
- ``p_value`` / ``log2err`` come from the multilevel result, ``q_value`` is BH of those p-values
  within the collection, and the p floor 1/(1 + modeFraction) is gone;
- a row whose multilevel ``log2err`` is NaN with a finite p (``isCpGeHalf`` false) keeps its p and
  carries the flag ``log2err_unestimated`` (not ``unestimated``).
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import as_list, gene_sets, singleton_frame
from scipy.special import polygamma

from giap.config.model import GseaConfig
from giap.rng import for_step
from giap.stats.fgsea import SimpleResult
from giap.stats.multitest import bh

STREAM = "function.gsea"


def _g() -> Any:
    from giap.function import gsea

    return gsea


def _ml_error(p: float, s: int = 101) -> float:
    single = float(polygamma(1, (s + 1) / 2) - polygamma(1, s + 1))
    return math.sqrt(math.floor(-math.log2(p) + 1) * single) / math.log(2)


def _genes(n: int = 300, seed: int = 3) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return singleton_frame(20.0 + 0.01 * rng.permutation(n))


def _ranked_tags(genes: pd.DataFrame) -> list[str]:
    return list(genes.sort_values("abundance_log2", ascending=False)["locus_tag"])


def _run(genes: pd.DataFrame, colls: Sequence[Any], cfg: GseaConfig | None = None) -> pd.DataFrame:
    cfg = cfg or GseaConfig()
    return _g().gsea_table(genes, colls, (10, 500), cfg, 0.05, for_step(STREAM, 1))


def _row(out: pd.DataFrame, term: str) -> pd.Series:
    rows = out[out["term"] == term]
    assert len(rows) == 1
    return rows.iloc[0]


def test_gsea_module_uses_fgsea_multilevel() -> None:
    from giap.stats.fgsea import fgsea_multilevel

    assert _g().fgsea_multilevel is fgsea_multilevel


def test_method_is_fgsea_multilevel_and_the_floor_is_gone() -> None:
    genes = _genes()
    tags = _ranked_tags(genes)
    coll = gene_sets("kegg", {"top": tags[:20], "mid": tags[140:160], "spread": tags[::15]})
    out = _run(genes, [coll])
    assert len(out) == 3
    assert set(out["method"]) == {"fgsea_multilevel"}
    r = _row(out, "top")
    assert 0 < r["p_value"] < 1e-6  # far below 1/(1 + modeFraction)
    assert r["log2err"] == pytest.approx(_ml_error(float(r["p_value"])), rel=1e-9)
    assert "log2err_unestimated" not in as_list(r["flags"])
    np.testing.assert_allclose(out["q_value"], bh(out["p_value"].to_numpy(float)), atol=1e-15)


def test_method_is_fgsea_multilevel_when_no_set_needs_the_stage() -> None:
    genes = _genes()
    tags = _ranked_tags(genes)
    coll = gene_sets("cog", {f"s{j}": tags[j::25] for j in range(3)})
    out = _run(genes, [coll])
    assert len(out) == 3
    assert set(out["method"]) == {"fgsea_multilevel"}
    assert (out["p_value"] > 0.01).all()


class _Spy:
    def __init__(self, real: Any, nan_log2err_for: frozenset[str] = frozenset()) -> None:
        self.real = real
        self.nan_for = nan_log2err_for
        self.calls: list[tuple[tuple[str, ...], int, float, bool]] = []

    def __call__(
        self, prep: Any, simple: SimpleResult, sample_size: int, eps: float, rng: Any
    ) -> SimpleResult:
        self.calls.append((tuple(prep.names), sample_size, eps, isinstance(simple, SimpleResult)))
        out = self.real(prep, simple, sample_size, eps, rng)
        if not self.nan_for:
            return out
        err = np.array(out.log2err, dtype=float)
        for i, name in enumerate(prep.names):
            if name in self.nan_for:
                err[i] = np.nan
        return dataclasses.replace(out, log2err=err)


def _install(monkeypatch: pytest.MonkeyPatch, **kw: Any) -> _Spy:
    from giap.stats.fgsea import fgsea_multilevel

    spy = _Spy(fgsea_multilevel, **kw)
    monkeypatch.setattr(_g(), "fgsea_multilevel", spy)
    return spy


def test_called_with_the_config_sample_size_and_eps(monkeypatch: pytest.MonkeyPatch) -> None:
    spy = _install(monkeypatch)
    genes = _genes()
    tags = _ranked_tags(genes)
    colls = [
        gene_sets("kegg", {"top": tags[:20], "spread": tags[::15]}),
        gene_sets("cog", {"bottom": tags[-25:]}),
    ]
    cfg = GseaConfig(sample_size=51, eps=1e-12)
    _run(genes, colls, cfg)
    wanted = ({"spread", "top"}, {"bottom"})
    main_calls = [c for c in spy.calls if set(c[0]) in wanted]
    assert len(main_calls) == 2, spy.calls
    for _, ss, eps, is_simple in main_calls:
        assert ss == 51 and eps == 1e-12 and is_simple


def test_eps_from_the_config_floors_the_row_p() -> None:
    genes = _genes()
    tags = _ranked_tags(genes)
    out = _run(genes, [gene_sets("kegg", {"top": tags[:20]})], GseaConfig(eps=1e-4))
    r = _row(out, "top")
    assert r["p_value"] == 1e-4
    assert math.isnan(r["log2err"])


def test_cp_ge_half_false_flags_log2err_unestimated(monkeypatch: pytest.MonkeyPatch) -> None:
    _install(monkeypatch, nan_log2err_for=frozenset({"top"}))
    genes = _genes()
    tags = _ranked_tags(genes)
    out = _run(genes, [gene_sets("kegg", {"top": tags[:20], "spread": tags[::15]})])
    r = _row(out, "top")
    assert math.isnan(r["log2err"])
    assert math.isfinite(r["p_value"]) and r["p_value"] < 1e-6  # p kept
    assert math.isfinite(r["q_value"])
    flags = as_list(r["flags"])
    assert "log2err_unestimated" in flags
    assert "unestimated" not in flags
    assert "log2err_unestimated" not in as_list(_row(out, "spread")["flags"])
