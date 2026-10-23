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

M8-14 function scenarios (drawn with ``giapsim.function_sim.draw``; the ``ExpressionConfig`` is a
``FunctionExpressionConfig`` carrying the scenario's ``FunctionSpec``):

- ``ora_null``: covariate-only detection, categories independent of detection;
- ``ora_attack_confounded``: one GO category's members get ``member_shift`` on the detection logit
  (members drawn at random, lead decision 3 gap 7); detection still depends on p_true only;
- ``gsea_null``: the M4-08 model without position decay, sets independent of abundance;
- ``gsea_attack_operon_packed``: ``n_sets`` KEGG sets of 1-3 packed chr1 TUs (10-30 genes) each
  with a large random TU effect; no set effect;
- ``determinants_null`` / ``determinants_planted``: TU effect (sd 1) + noise (``noise_sd``),
  plus ``sum beta_j z(feature_j)`` for ``promoter``, ``rbs``, ``cai`` in the planted one.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Literal

from giapsim.expression import DEFAULT_BETA, ExpressionConfig, FrozenMap, IslandSpec
from giapsim.function_sim import FEATURES, PACKED_TU_SIZE, FunctionExpressionConfig, FunctionSpec
from giapsim.layout import LayoutConfig

__all__ = [
    "ABUNDANCE_BETA",
    "ORIC_FRAC",
    "ar1_background",
    "determinants_null",
    "determinants_planted",
    "draft_ar1",
    "gsea_attack_operon_packed",
    "gsea_null",
    "island_first_tu",
    "island_power",
    "label_null",
    "ora_attack_confounded",
    "ora_null",
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


# --------------------------------------------------------------------------- M8-14 function

_MIN_FUNCTION_GENES = 100  # the default regions alone take 81 chr1 genes
_MAX_SETS = 20  # the 2-20 synthetic KEGG maps


def _function_n(n_genes: int) -> int:
    if isinstance(n_genes, bool) or not isinstance(n_genes, int) or n_genes < _MIN_FUNCTION_GENES:
        raise ValueError(f"n_genes must be an int >= {_MIN_FUNCTION_GENES}, got {n_genes!r}")
    return n_genes


def ora_null(n_genes: int) -> Pair:
    """Detection ~ Bernoulli(p_true) on the covariates only; categories independent of it."""
    spec = FunctionSpec("ora")
    return LayoutConfig(n_genes_chromosome=_function_n(n_genes)), FunctionExpressionConfig(
        function=spec
    )


def ora_attack_confounded(n_genes: int, member_shift: float) -> Pair:
    """``ora_null`` plus one confounded GO category whose members' detection logit moves by
    ``member_shift`` (0 reproduces ``ora_null``'s detection truth)."""
    if not math.isfinite(member_shift):
        raise ValueError(f"member_shift must be finite, got {member_shift}")
    spec = FunctionSpec("ora", attack=True, member_shift=float(member_shift))
    return LayoutConfig(n_genes_chromosome=_function_n(n_genes)), FunctionExpressionConfig(
        function=spec
    )


def gsea_null(n_genes: int) -> Pair:
    """M4-08 abundance (TU level + within-TU noise, no position decay); sets independent."""
    ecfg = FunctionExpressionConfig(position_decay_log2=0.0, function=FunctionSpec("gsea"))
    return LayoutConfig(n_genes_chromosome=_function_n(n_genes)), ecfg


def gsea_attack_operon_packed(n_genes: int, n_sets: int) -> Pair:
    """``n_sets`` sets of 1, 2, 3, 1, ... packed chr1 TUs of 10-30 genes; no set effect."""
    n = _function_n(n_genes)
    if isinstance(n_sets, bool) or not isinstance(n_sets, int) or not 1 <= n_sets <= _MAX_SETS:
        raise ValueError(f"n_sets must be an int in [1, {_MAX_SETS}], got {n_sets!r}")
    lo, hi = PACKED_TU_SIZE
    n_tus = sum(1 + k % 3 for k in range(n_sets))
    sizes = tuple(lo + (7 * j) % (hi - lo + 1) for j in range(n_tus))
    if sum(sizes) > n // 2:
        raise ValueError(f"{n_sets} packed sets need {sum(sizes)} genes, more than half of {n}")
    ecfg = FunctionExpressionConfig(position_decay_log2=0.0,
                                    function=FunctionSpec("packed", n_sets=n_sets))  # fmt: skip
    return LayoutConfig(n_genes_chromosome=n, packed_tus=sizes), ecfg


def _determinants(n_genes: int, coefs: Mapping[str, float], noise_sd: float,
                  features: str) -> Pair:  # fmt: skip
    n = _function_n(n_genes)
    bad = set(coefs) - set(FEATURES)
    if bad:
        raise ValueError(f"unknown determinants features {sorted(bad)}; use {list(FEATURES)}")
    if not all(math.isfinite(float(v)) for v in coefs.values()):
        raise ValueError(f"coefs must be finite, got {dict(coefs)}")
    if not (math.isfinite(noise_sd) and noise_sd >= 0):
        raise ValueError(f"noise_sd must be finite and >= 0, got {noise_sd}")
    if features not in ("planted", "gaussian"):
        raise ValueError(f"features must be 'planted' or 'gaussian', got {features!r}")
    spec = FunctionSpec("determinants", coefs={f: float(coefs.get(f, 0.0)) for f in FEATURES},
                        features=features)  # type: ignore[arg-type]  # fmt: skip
    ecfg = FunctionExpressionConfig(tu_sd=1.0, ar1_rho=0.0, position_decay_log2=0.0,
                                    stoich_sd=0.0, gene_sd=float(noise_sd),
                                    function=spec)  # fmt: skip
    return LayoutConfig(n_genes_chromosome=n), ecfg


def determinants_null(n_genes: int, *, noise_sd: float = 1.0, features: str = "planted") -> Pair:
    """Response = TU effect (sd 1) + noise (``noise_sd``), independent of every feature."""
    return _determinants(n_genes, {}, noise_sd, features)


def determinants_planted(
    n_genes: int, coefs: Mapping[str, float], noise_sd: float, *, features: str = "planted"
) -> Pair:
    """``determinants_null`` plus ``sum_j coefs[j] * z(feature_j)`` (features ``promoter``,
    ``rbs``, ``cai``; missing ones 0); ``features="gaussian"`` draws iid N(0, 1) features."""
    return _determinants(n_genes, coefs, noise_sd, features)
