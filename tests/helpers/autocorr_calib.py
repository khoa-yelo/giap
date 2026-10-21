"""Replicate machinery for the M5-05 autocorrelation calibration suite (GIAP_SPEC.md §8.4).

One simulated genome per replicate feeds every (statistic, null_model, criterion) check of a
scenario, as in ``region_tests_calib``: each check has its own ``_criteria`` stream, and
replicate ``i`` of any stream maps back to the scenario's replicate ``i`` (seed
``suite_seed(<scenario stream>, i)``), simulated and tabulated once per process.

**Scenarios.** docs/specs/M5-04 defines ``giapsim.scenarios.label_null`` / ``tu_null``; that
module is not merged yet (PR #119), so the two null scenarios are spelled out here with the
merged M4-07 layout and M4-08 expression knobs and swapped for the simulator's versions when
#119 lands:

* ``label_null``: ``tu_sd = 0``, ``ar1_rho = 0``, ``position_decay_log2 = 0``, iid gene and
  stoichiometry noise, covariate-only detection -> no structure at all, so both nulls hold
  (§8.4 "gene-level noise (label null)");
* ``tu_null``: ``ar1_rho = 0``, ``tu_sd = 2``, decay 0.2, ``detection = "abundance"`` (beta
  ``log2_abundance`` 0.5 with the M4-08 intercept offset) -> within-TU structure only, so the
  TU-block null holds and the gene-label null is false (§8.4 "TU-structured noise without
  inter-TU structure (TU null)"), and detection inherits the TU structure.

A chromosome alone (no plasmid, no regions) keeps the nulls exact: regions and plasmids would
add no effect here but change gene sizes and the TU grid for nothing.

``p_detect`` = ``p_detect_true`` and the abundance of detected genes = ``abundance_log2_true``
(spec acceptance "Calibration"); B = ``PERMUTATIONS`` permutations per null.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
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
#: detection logit slope on log2 abundance for the TU null (M4-08 "detection = abundance")
TU_NULL_BETA_ABUNDANCE = 0.5


def _beta(**shift: float) -> dict[str, float]:
    import giapsim.expression as ex

    b = dict(ex.DEFAULT_BETA)
    for k, v in shift.items():
        b[k] = b[k] + v
    return b


@dataclass(frozen=True)
class Scenario:
    """A null generator: chromosome size and M4-08 expression knobs."""

    name: str
    n_genes: int
    expr_kw: Mapping[str, Any] = field(default_factory=dict)
    layout_kw: Mapping[str, Any] = field(default_factory=dict)
    #: (statistic, null_model) pairs whose level this scenario can gate
    nulls: tuple[str, ...] = ("gene_label", "tu_block")

    @property
    def stream(self) -> str:
        return f"{SUITE}.{self.name}"


def label_null(n_genes: int) -> Scenario:
    """No structure at all: both nulls hold (M5-04 ``scenarios.label_null``)."""
    return Scenario(
        f"label_null_{n_genes}",
        n_genes,
        {
            "tu_sd": 0.0,
            "ar1_rho": 0.0,
            "position_decay_log2": 0.0,
            "gene_sd": 0.5,
            "stoich_sd": 0.7,
            "detection": "covariates_only",
        },
    )


def tu_null(n_genes: int) -> Scenario:
    """Within-TU structure only: the TU-block null holds (M5-04 ``scenarios.tu_null``)."""
    b4 = TU_NULL_BETA_ABUNDANCE
    return Scenario(
        f"tu_null_{n_genes}",
        n_genes,
        {
            "tu_sd": 2.0,
            "ar1_rho": 0.0,
            "position_decay_log2": 0.2,
            "detection": "abundance",
            "beta": _beta(log2_abundance=b4, intercept=-b4 * 20.0),
        },
        nulls=("tu_block",),
    )


def layout_config(n_genes: int, **kw: Any) -> Any:
    """One circular chromosome, no plasmid and no region (origin crossing on, M4-07)."""
    import giapsim.layout as la

    base: dict[str, Any] = {
        "n_genes_chromosome": int(n_genes),
        "plasmids": (),
        "regions": (),
        "origin_crossing": True,
    }
    base.update(kw)
    return la.LayoutConfig(**base)


_FRAMES: dict[tuple[str, int], pd.DataFrame] = {}


def _C() -> Any:
    import _criteria as C

    return C


def replicate_rows(sc: Scenario, i: int) -> pd.DataFrame:
    """``autocorrelation_rows`` of replicate ``i`` of scenario ``sc`` (cached per process)."""
    key = (sc.stream, i)
    if key not in _FRAMES:
        rng = _C().replicate_rng(sc.stream, i)
        tables = rd.sim_tables(
            rng, layout_config(sc.n_genes, **dict(sc.layout_kw)), sc.expr_kw, p_source="true"
        )
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
