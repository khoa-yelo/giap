"""Null and power scenarios for the M5 calibration and power runs (issue M5-04; GIAP_SPEC.md §8.4,
§8.5 attack 1, §11.3).

Each scenario is a pure function returning ``(LayoutConfig, ExpressionConfig)``; the caller draws
the layout and the expression with its own generators. All use the circular ``chr1`` of ``n_genes``
genes and no other replicon (except ``draft_ar1``), the default regions with no planted region
effect, and ``detection = "covariates_only"`` unless stated. A nonzero ori-ter gradient always comes
with ``oric_frac = 0.8`` and the cosine shape (lead decision 16).

- ``label_null``: no TU, position or AR(1) structure (``tu_sd = 0``, no decay); gene abundance is
  iid normal (stoichiometry + noise), so both the gene-label and the TU-block nulls hold;
- ``tu_null``: TU-structured noise without inter-TU structure (``ar1_rho = 0``, ``tu_sd = 2``,
  decay 0.2); ``detection = "abundance"`` by default so detection inherits the TU structure;
- ``ar1_background``: the default expression model at AR(1) ``rho``, no islands;
- ``island_power``: ``ar1_background`` plus one chr1 island at a fixed TU start (lead decision 17:
  no seed; whether it overlaps a region is left to the caller);
- ``oriter_gradient``: ``ar1_background`` with an oriC and the cosine gradient;
- ``draft_ar1``: a draft assembly (500-gene chromosome + 133 linear 30-gene contigs) at ``rho``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Literal

from giapsim.expression import DEFAULT_BETA, ExpressionConfig, FrozenMap, IslandSpec
from giapsim.layout import LayoutConfig

__all__ = [
    "ABUNDANCE_BETA",
    "ORIC_FRAC",
    "ar1_background",
    "draft_ar1",
    "island_first_tu",
    "island_power",
    "label_null",
    "oriter_gradient",
    "tu_null",
]

Detection = Literal["covariates_only", "abundance"]
Pair = tuple[LayoutConfig, ExpressionConfig]

#: oriC position of every scenario with a nonzero gradient (lead decision 16).
ORIC_FRAC = 0.8
_B4 = 0.5  # detection logit per log2 abundance under ``detection = "abundance"``
_MU = ExpressionConfig().mu_log2
#: Detection coefficients for ``detection = "abundance"``: slope ``_B4`` on log2 abundance with
#: the intercept offset ``-_B4 * mu`` so the detected fraction stays near the default.
ABUNDANCE_BETA: Mapping[str, float] = FrozenMap(
    {**DEFAULT_BETA, "log2_abundance": _B4, "intercept": DEFAULT_BETA["intercept"] - _B4 * _MU}
)
_DRAFT_CONTIGS = (30,) * 133


def _beta(detection: Detection) -> Mapping[str, float]:
    if detection not in ("covariates_only", "abundance"):
        raise ValueError(f"unknown detection mode {detection!r}")
    return ABUNDANCE_BETA if detection == "abundance" else DEFAULT_BETA


def _layout(n_genes: int, oriter_log2: float, short_genes_in: str | None = None) -> LayoutConfig:
    frac = ORIC_FRAC if oriter_log2 != 0.0 else None
    return LayoutConfig(n_genes_chromosome=n_genes, oric_frac=frac, short_genes_in=short_genes_in)


def label_null(n_genes: int, *, oriter_log2: float = 0.0) -> Pair:
    """Gene-level iid noise only (both nulls hold); ``oriter_log2 != 0`` adds the cosine gradient
    (then the gradient is the only structure)."""
    ecfg = ExpressionConfig(tu_sd=0.0, ar1_rho=0.0, position_decay_log2=0.0,
                            oriter_log2=oriter_log2)  # fmt: skip
    return _layout(n_genes, oriter_log2), ecfg


def tu_null(
    n_genes: int,
    detection: Detection = "abundance",
    *,
    short_genes_in: str | None = None,
    oriter_log2: float = 0.0,
) -> Pair:
    """Within-TU correlation only: the TU-block null holds, the gene-label null does not."""
    ecfg = ExpressionConfig(tu_sd=2.0, ar1_rho=0.0, position_decay_log2=0.2,
                            detection=detection, beta=_beta(detection),
                            oriter_log2=oriter_log2)  # fmt: skip
    return _layout(n_genes, oriter_log2, short_genes_in=short_genes_in), ecfg


def ar1_background(n_genes: int, rho: float) -> Pair:
    """The default expression model at AR(1) ``rho``: no islands, no region effects."""
    return _layout(n_genes, 0.0), ExpressionConfig(ar1_rho=rho)


def island_first_tu(n_genes: int, mean_tu_size: float = LayoutConfig().mean_tu_size) -> int:
    """0-based chr1 TU index where ``island_power`` starts its island: a quarter of the expected
    TU count, so islands of up to half the TUs fit every layout drawn from the scenario."""
    return int(n_genes / (4.0 * mean_tu_size))


def island_power(
    n_genes: int,
    rho: float,
    n_tus: int,
    effect_sd: float,
    *,
    detection: Detection = "covariates_only",
) -> Pair:
    """``ar1_background`` plus one chr1 island of ``n_tus`` TUs shifted by ``effect_sd`` stationary
    SDs (the sign gives the direction), starting at ``island_first_tu(n_genes)``."""
    island = IslandSpec("chr1", island_first_tu(n_genes), n_tus, effect_sd)
    ecfg = ExpressionConfig(ar1_rho=rho, islands=(island,), detection=detection,
                            beta=_beta(detection))  # fmt: skip
    return _layout(n_genes, 0.0), ecfg


def oriter_gradient(n_genes: int, rho: float, oriter_log2: float) -> Pair:
    """``ar1_background`` with an oriC at ``ORIC_FRAC`` and the cosine ori-ter gradient."""
    lcfg = LayoutConfig(n_genes_chromosome=n_genes, oric_frac=ORIC_FRAC)
    return lcfg, ExpressionConfig(ar1_rho=rho, oriter_log2=oriter_log2, oriter_shape="cosine")


def draft_ar1(rho: float) -> Pair:
    """A draft assembly at AR(1) ``rho``: 500-gene ``chr1`` plus 133 linear 30-gene contigs."""
    lcfg = LayoutConfig(n_genes_chromosome=500, linear_contigs=_DRAFT_CONTIGS)
    return lcfg, ExpressionConfig(ar1_rho=rho)
