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

Operon scenarios (M6-03; drawn with ``simulate_operons``, see ``giapsim.operons``): ``operon_null``,
``missed_operon_power``, ``discordance_power`` and the attacks ``operon_attack_ar1`` and
``operon_attack_detection``. All use the generating t ``OPERON_T`` for the within-TU increment, no
position decay, and the detection intercept ``OPERON_BETA`` (~80 % of genes detected, so 1,800 genes
give >= 300 OP, >= 200 NO, >= 100 CV and >= 30 SW pairs with both genes detected).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from typing import Any, Literal

from giapsim.expression import DEFAULT_BETA, ExpressionConfig, FrozenMap, IslandSpec
from giapsim.layout import DEFAULT_REGIONS, LayoutConfig, RegionSpec
from giapsim.operons import OperonConfig, OperonSim, simulate_operons

__all__ = [
    "ABUNDANCE_BETA",
    "ATTACK_BETA",
    "ATTACK_BLOCK_FRAC",
    "ATTACK_DETECT_B4",
    "ATTACK_DETECT_OFFSET",
    "BLOCK_DETECT_SHIFT",
    "OPERON_BETA",
    "OPERON_P_SAME_STRAND",
    "OPERON_T",
    "ORIC_FRAC",
    "OperonSim",
    "ar1_background",
    "discordance_power",
    "draft_ar1",
    "island_first_tu",
    "island_power",
    "label_null",
    "missed_operon_power",
    "operon_attack_ar1",
    "operon_attack_detection",
    "operon_null",
    "oriter_gradient",
    "simulate_operons",
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


#: Generating Student-t (loc, scale, df) of the within-TU increment (lead decision 8: the
#: methods-note simulation, Slurm 1864589).
OPERON_T: tuple[float, float, float] = (-0.3, 0.5, 5.0)
_OPERON_INTERCEPT = -7.3  # ~80 % detected on the default layout (DEFAULT_BETA gives ~55 %)
#: Detection coefficients of the operon scenarios (``covariates_only``).
OPERON_BETA: Mapping[str, float] = FrozenMap({**DEFAULT_BETA, "intercept": _OPERON_INTERCEPT})


#: Same-strand probability of adjacent TUs in the operon scenarios: half the TU boundaries flip,
#: so CV, DV and SW decoys are frequent (§11.3 "set convergent and divergent frequencies").
OPERON_P_SAME_STRAND = 0.5


#: ``operon_attack_detection`` detection (lead decision 10): logit = the merged covariate logit at
#: the operon intercept + ``ATTACK_DETECT_B4`` x (final abundance - 20) + ``ATTACK_DETECT_OFFSET``,
#: plus ``BLOCK_DETECT_SHIFT`` on the block. Detection is then mostly a soft threshold on the
#: abundance: ~94 % detected outside the block, ~49 % inside (>= 25 % of block pairs both detected).
ATTACK_DETECT_B4 = 2.0
ATTACK_DETECT_OFFSET = 7.0
#: Logit shift on the short, hydrophobic block (= -3 log2 of abundance at slope 2).
BLOCK_DETECT_SHIFT = -6.0
ATTACK_BETA: Mapping[str, float] = FrozenMap(
    {**DEFAULT_BETA, "log2_abundance": ATTACK_DETECT_B4,
     "intercept": _OPERON_INTERCEPT + ATTACK_DETECT_OFFSET - 20.0 * ATTACK_DETECT_B4}
)  # fmt: skip
#: Share of the chromosome's genes in the ``operon_attack_detection`` block: the within-bin excess
#: lives only in the block's pairs, so a 40-gene block cannot move a 1,800-gene CMH z (R = 20).
ATTACK_BLOCK_FRAC = 1.0 / 3.0


def _regions(n_genes: int) -> tuple[RegionSpec, ...]:
    """The default regions, shrunk proportionally (>= 1 gene each) when they would take more than
    a quarter of a small genome."""
    total = sum(r.n_genes for r in DEFAULT_REGIONS)
    if 4 * total <= n_genes:
        return DEFAULT_REGIONS
    f = n_genes / (4.0 * total)
    return tuple(replace(r, n_genes=max(1, int(r.n_genes * f))) for r in DEFAULT_REGIONS)


def _lcfg(n_genes: int, **kw: Any) -> LayoutConfig:
    return LayoutConfig(n_genes_chromosome=n_genes, p_same_strand=OPERON_P_SAME_STRAND,
                        regions=_regions(n_genes), **kw)  # fmt: skip


def _operon_cfg(rho: float = 0.0, **plan: Any) -> OperonConfig:
    return OperonConfig(ar1_rho=rho, position_decay_log2=0.0, beta=OPERON_BETA,
                        operon_t=OPERON_T, **plan)  # fmt: skip


def operon_null(n_genes: int, *, rho: float = 0.0) -> Pair:
    """Missed-operon and discordance null: OP Delta ~ ``OPERON_T``, NO pairs independent (TU-level
    AR(1) at ``rho``), pure UniOP classes (lead decision 4) and no internal promoter or terminator
    (lead decision 1). Small genomes are allowed; the M6-09 suites need >= 1,000 genes."""
    return _lcfg(n_genes, uniop_pure=True), _operon_cfg(rho)


def _plant(cfgs: Pair, n_planted: int, effect_sd: float) -> Pair:
    """``cfgs`` with ``n_planted`` co-abundant NO pairs (the ``missed_operon_power`` mechanism)."""
    if n_planted == 0 and effect_sd == 0.0:
        return cfgs
    lcfg, ecfg = cfgs
    assert isinstance(ecfg, OperonConfig)
    return lcfg, replace(ecfg, n_missed=n_planted, missed_effect_sd=effect_sd)


def missed_operon_power(n_genes: int, n_planted: int, effect_sd: float) -> Pair:
    """``n_planted`` true NO pairs (UniOP false negatives) with ``|Delta| ~ effect_sd`` x the
    background NO SD; UniOP keeps its own errors (default Beta classes) and convergent TU
    boundaries get -20...300 bp gaps (negative: 3' overlaps), so CV decoys fill every bin up to
    300 bp (lead decisions 11, 23). A ``ValueError`` from ``simulate_operons`` when fewer NO pairs
    exist. The M6-09 suite needs >= 1,000 genes."""
    lcfg = _lcfg(n_genes, short_convergent_gaps=True)
    return lcfg, _operon_cfg(n_missed=n_planted, missed_effect_sd=effect_sd)


def discordance_power(n_genes: int, n_planted: int, effect_log2: float) -> Pair:
    """``operon_null`` plus ``n_planted`` true OP pairs shifted by ``+effect_log2`` over a planted
    internal promoter (ceil(n/2)) or ``-effect_log2`` over an internal terminator (floor(n/2)). A
    ``ValueError`` from ``simulate_operons`` when fewer OP pairs exist; the suite needs >= 1,000
    genes."""
    lcfg = _lcfg(n_genes, uniop_pure=True)
    return lcfg, _operon_cfg(n_discordant=n_planted, discordant_effect_log2=effect_log2)


def operon_attack_ar1(
    n_genes: int, rho: float, *, n_planted: int = 0, effect_sd: float = 0.0
) -> Pair:
    """Attack 1 (§8.5): ``operon_null(n_genes, rho=rho)`` (lead decision 7), TU abundance AR(1)
    along the chromosome: neighbours look alike by position, which the distance-reweighted decoys
    must absorb. ``rho = 0`` is ``operon_null``. ``n_planted`` / ``effect_sd`` (default nothing)
    plant co-abundant NO pairs as ``missed_operon_power`` does (lead decision 13), so the FDR gate
    has calls. The suite needs >= 1,000 genes."""
    return _plant(operon_null(n_genes, rho=rho), n_planted, effect_sd)


def operon_attack_detection(
    n_genes: int, short_genes_in: str, *, n_planted: int = 0, effect_sd: float = 0.0
) -> Pair:
    """Detection-distance confounding (ADR-0018; lead decisions 10, 11): the ``short_genes_in``
    block (a third of the chromosome, ``ATTACK_BLOCK_FRAC``) has short genes and intra-TU gaps
    between all of them; detection depends on the final (post-walk, post-planting) abundance
    (``ATTACK_BETA``) with the block shift ``BLOCK_DETECT_SHIFT`` standing in for hydrophobicity;
    convergent TU boundaries get -20...300 bp gaps. The both-detected filter keeps only the block's
    high-abundance genes, so its NO pairs (the targets of the < 50 bp bins) have smaller |Delta|
    than the CV decoys from outside the block sharing those bins: the attack on the decoy
    distance reweighting. UniOP keeps its own errors. ``n_planted`` / ``effect_sd`` plant as
    ``missed_operon_power`` does (lead decision 13). The suite needs >= 1,000 genes."""
    k = int(n_genes * ATTACK_BLOCK_FRAC)
    regions = tuple(replace(r, n_genes=max(r.n_genes, k)) if r.region_type == short_genes_in
                    else r for r in _regions(n_genes))  # fmt: skip
    lcfg = LayoutConfig(n_genes_chromosome=n_genes, p_same_strand=OPERON_P_SAME_STRAND,
                        regions=regions, short_genes_in=short_genes_in,
                        short_gaps_in=short_genes_in, short_convergent_gaps=True)  # fmt: skip
    shift = FrozenMap({short_genes_in: BLOCK_DETECT_SHIFT})
    ecfg = OperonConfig(ar1_rho=0.0, position_decay_log2=0.0, beta=ATTACK_BETA, operon_t=OPERON_T,
                        detection="abundance", region_detect_shift=shift)  # fmt: skip
    return _plant((lcfg, ecfg), n_planted, effect_sd)
