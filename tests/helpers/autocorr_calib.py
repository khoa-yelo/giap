"""Replicate machinery for the M5-05 autocorrelation calibration suite (GIPPA_SPEC.md §8.4).

One simulated genome per replicate feeds every (statistic, null_model, criterion) check of a
scenario, as in ``region_tests_calib``: each check has its own ``_criteria`` stream, and
replicate ``i`` of any stream maps back to the scenario's replicate ``i`` (seed
``suite_seed(<scenario stream>, i)``), simulated and tabulated once per process.

**Scenarios** are the M5-04 generators (``gippasim.scenarios``, merged in PR #119; lead
decision 33); this module only wraps them in a ``Scenario`` with its ``_criteria`` stream:

* ``label_null``: ``tu_sd = 0``, ``ar1_rho = 0``, no position decay, iid gene and stoichiometry
  noise, covariate-only detection -> no structure at all, so both nulls hold (§8.4 "gene-level
  noise (label null)");
* ``tu_null``: ``ar1_rho = 0``, ``tu_sd = 2``, decay 0.2, ``detection = "abundance"`` -> within-TU
  structure only, so the TU-block null holds and the gene-label null does not (§8.4
  "TU-structured noise without inter-TU structure"), and detection inherits the TU structure.

Both take ``oriter_log2`` (the ADR-0020 ori-ter cosine of the attack arms) and ``tu_null`` takes
``short_genes_in``; the scenarios carry the default regions (no planted effect) and, with a
nonzero gradient, an oriC at ``scenarios.ORIC_FRAC``.

``p_detect`` = ``p_detect_true`` and the abundance of detected genes = ``abundance_log2_true``
(spec acceptance "Calibration"); B = ``PERMUTATIONS`` permutations per null.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import autocorr_data as ad
import numpy as np
import pandas as pd
import region_tests_data as rd

#: permutations per null and replicon in calibration runs (the smallest p is 1 / (B + 1) = 0.005)
PERMUTATIONS = 199
#: the lags of ``context.moran_lags``
LAGS = tuple(range(1, 21))
SUITE = "autocorrelation"
CONTIG = "chr1"


@dataclass(frozen=True)
class Scenario:
    """One M5-04 scenario: its name (the ``_criteria`` stream), size and config pair."""

    name: str
    n_genes: int
    #: ``(LayoutConfig, ExpressionConfig)`` from ``gippasim.scenarios``
    pair: tuple[Any, Any]
    #: null models whose level this scenario can gate
    nulls: tuple[str, ...] = ("gene_label", "tu_block")

    @property
    def stream(self) -> str:
        return f"{SUITE}.{self.name}"

    @property
    def layout_cfg(self) -> Any:
        return self.pair[0]

    @property
    def expr_cfg(self) -> Any:
        return self.pair[1]


def label_null(n_genes: int, *, oriter_log2: float = 0.0, name: str | None = None) -> Scenario:
    """``scenarios.label_null``: no structure at all, so both nulls hold."""
    import gippasim.scenarios as scen

    pair = scen.label_null(n_genes, oriter_log2=oriter_log2)
    return Scenario(name or f"label_null_{n_genes}", n_genes, pair)


def tu_null(
    n_genes: int,
    *,
    short_genes_in: str | None = None,
    oriter_log2: float = 0.0,
    name: str | None = None,
) -> Scenario:
    """``scenarios.tu_null``: within-TU structure only, so the TU-block null holds."""
    import gippasim.scenarios as scen

    return Scenario(
        name or f"tu_null_{n_genes}",
        n_genes,
        scen.tu_null(n_genes, short_genes_in=short_genes_in, oriter_log2=oriter_log2),
        nulls=("tu_block",),
    )


def simulate(sc: Scenario, rng: np.random.Generator) -> Any:
    """The M4-06 tables of one draw of ``sc`` (``p_detect`` = ``p_detect_true``)."""
    import gippasim.expression as ex
    import gippasim.layout as la

    lay = la.simulate_layout(sc.layout_cfg, rng)
    expr = ex.simulate_expression(lay, sc.expr_cfg, rng)
    return rd.tables_from(lay, expr, "true", rng)


_FRAMES: dict[tuple[str, int], pd.DataFrame] = {}


def _C() -> Any:
    import _criteria as C

    return C


def replicate_rows(sc: Scenario, i: int) -> pd.DataFrame:
    """``autocorrelation_rows`` of replicate ``i`` of scenario ``sc`` (cached per process)."""
    key = (sc.stream, i)
    if key not in _FRAMES:
        rng = _C().replicate_rng(sc.stream, i)
        tables = simulate(sc, rng)
        seed = int(rng.integers(0, 2**31 - 1))
        cfg = ad.config(PERMUTATIONS, lags=(LAGS[0], LAGS[-1]))
        _FRAMES[key] = ad.rows(tables.genes, tables.replicons, cfg, seed)
    return _FRAMES[key]


def shared_fn(
    sc: Scenario, stream: str, R: int, value: Callable[[pd.DataFrame], Any]
) -> Callable[[np.random.Generator], Any]:
    """Replicate function for ``calibrate(..., stream=stream)``: maps the generator handed in by
    ``_criteria`` (seed ``suite_seed(stream, i)``) to replicate ``i`` of ``sc``."""
    C = _C()
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def fn(rng: np.random.Generator) -> Any:
        key = int(rng.bit_generator.seed_seq.entropy)  # type: ignore[attr-defined]
        return value(replicate_rows(sc, index[key]))

    return fn


def p_or_one(p: Any) -> float:
    """An NA p (degenerate null, invalid lag) does not reject."""
    x = float(p)
    return 1.0 if np.isnan(x) else x


def lag_p(
    res: pd.DataFrame, statistic: str, null_model: str, lag: int, contig: str = CONTIG
) -> float:
    """The p of one (contig, statistic, null_model, lag) row; NA -> 1."""
    return p_or_one(ad.pick(res, contig, statistic, null_model, lag)["p_value"])


def family_rejection(
    res: pd.DataFrame,
    statistic: str,
    null_model: str,
    levels: Sequence[float],
    contig: str = CONTIG,
) -> np.ndarray:
    """Per level: 1 if the BH family (contig, statistic, null_model) rejects any lag (ADR-0020
    decision 5, FWER-style bound under the complete null)."""
    q = ad.block(res, contig, statistic, null_model)["q_value"].to_numpy(float)
    q = np.where(np.isnan(q), 1.0, q)
    return np.array([float((q <= a).any()) for a in levels])


def pooled_share(
    res: pd.DataFrame,
    statistic: str,
    null_model: str,
    levels: Sequence[float],
    contig: str = CONTIG,
) -> np.ndarray:
    """Per level: the share of the 20 lags with p <= level (recorded, not gated: decision 18)."""
    p = ad.block(res, contig, statistic, null_model)["p_value"].to_numpy(float)
    p = np.where(np.isnan(p), 1.0, p)
    return np.array([float(np.count_nonzero(p <= a)) / p.size for a in levels])
