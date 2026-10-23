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

M7-12 (docs/specs/M7-12-sim-variant-novel-scenarios.md; §7.5.6, §7.6.6, §8.4, §8.5 attack 6,
§11.6 TE2) adds four nucleotide scenarios returning a ``giapsim.SimConfig`` for
``giapsim.simulate(cfg, seed, outdir)``. The config does not depend on ``seed`` (the simulation
seed draws everything); ``size`` "100k" / "300k" / "1M" sets ``genome_length`` (100 kb is for
fast unit tests only; novel-ORF power and FDR gates use >= 300 kb, >= 200 genes, so GIAP's
detection q floor 1 / n is below 0.01: lead decision 13). All are simulator-v2
genomes (signals at strength 0.8) with ``m7_truth`` (``truth/alt_sites.tsv`` always, header only
without a strain; see ``giapsim.scenario_truth``):

- ``variant_f1``: sample vs reference by SNVs (1 / kb) and indels (0.2 / kb); sample-derived
  database peptides (labelled ``search_db = same``: the reference identification carries the
  ALT, which is what the pairing needs), de novo reads with every error type and paired
  ``reference_peptides.tsv``;
- ``novel_orf_power``: intergenic ORFs of 30-80 aa, out-of-frame and antisense ORFs (>= 30 aa)
  found in CDSs and N-terminal extensions, ``NOVEL_PER_CLASS[size]`` of each (fewer intergenic
  ones, with a WARNING, when the genome has too few gaps: ``truth/novel_orfs.tsv`` holds the
  realised rows and ``truth/novel_orfs.json`` the requested / realised / shortfall counts;
  ``n_per_class`` overrides the count, ``NOVEL_PER_CLASS_GATE`` for M7-13 G3/G4, decision 26);
  ``AGGAGG`` at spacer 5-10 before every planted intergenic start and N-terminal true start,
  whose ``rbs_novel`` TIRs are lognormal(``SimConfig.planted_tir``) = (8, 0.8), the
  first-round TU RBS distribution (decision 14: the matched arm; the decision 8 sanity arm sets
  (9, 0.5)); peptides with 2 % wrong identifications drawn from non-coding six-frame
  translations (``wrong_source = "sixframe_noncoding"``, decision 12: the false novel candidates
  the class-FDR gate controls) with matched decoys, de novo errors and the precomputed layers.
  Database peptides are sample-derived and labelled ``search_db = same``;
- ``same_strain_insilico``: every Trypsin/P peptide (7-30 aa, <= 1 missed cleavage) of every
  annotated protein, detected with probability 1, no wrong peptide, ``search_db = same``;
- ``isobaric_attack``: de novo errors only of the block, near-block, modification and swap
  types, no strain differences.

MVP-04 (docs/specs/MVP-04-ecoli-evidence-run.md) adds ``k12_evidence``, a
``giapsim.real_genome.RealGenomeConfig`` for ``simulate_on_genome`` on the committed K-12 genome at
``K12_EVIDENCE_SEED`` (its docstring records the tuned detection intercept b0).

Operon scenarios (M6-03; drawn with ``simulate_operons``, see ``giapsim.operons``): ``operon_null``,
``missed_operon_power``, ``discordance_power`` and the attacks ``operon_attack_ar1`` and
``operon_attack_detection``. All use the generating t ``OPERON_T`` for the within-TU increment, no
position decay, and the detection intercept ``OPERON_BETA`` (~80 % of genes detected, so 1,800 genes
give >= 300 OP, >= 200 NO, >= 100 CV and >= 30 SW pairs with both genes detected).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Literal

from giapsim.expression import DEFAULT_BETA, ExpressionConfig, FrozenMap, IslandSpec
from giapsim.layout import DEFAULT_REGIONS, LayoutConfig, RegionSpec
from giapsim.operons import OperonConfig, OperonSim, simulate_operons

if TYPE_CHECKING:
    from giapsim.real_genome import RealGenomeConfig
    from giapsim.simulate import SimConfig

__all__ = [
    "ABUNDANCE_BETA",
    "ATTACK_BETA",
    "ATTACK_BLOCK_FRAC",
    "ATTACK_DETECT_B4",
    "ATTACK_DETECT_OFFSET",
    "BLOCK_DETECT_SHIFT",
    "K12_B0",
    "K12_EVIDENCE_SEED",
    "K12_ISLANDS",
    "MATCHED_TIR",
    "NOVEL_PER_CLASS",
    "NOVEL_PER_CLASS_GATE",
    "OPERON_BETA",
    "OPERON_P_SAME_STRAND",
    "OPERON_T",
    "ORIC_FRAC",
    "SIZES",
    "OperonSim",
    "ar1_background",
    "discordance_power",
    "draft_ar1",
    "island_first_tu",
    "island_power",
    "isobaric_attack",
    "k12_evidence",
    "label_null",
    "missed_operon_power",
    "novel_orf_power",
    "operon_attack_ar1",
    "operon_attack_detection",
    "operon_null",
    "oriter_gradient",
    "same_strain_insilico",
    "simulate_operons",
    "tu_null",
    "variant_f1",
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


# ---------------------------------------------------------------- M7-12 nucleotide scenarios

Size = Literal["100k", "300k", "1M"]
#: Genome length of each scenario size.
SIZES: Mapping[str, int] = MappingProxyType({"100k": 100_000, "300k": 300_000, "1M": 1_000_000})
#: Planted / found novel ORFs and N-terminal extensions per class of ``novel_orf_power``.
NOVEL_PER_CLASS: Mapping[str, int] = MappingProxyType({"100k": 8, "300k": 16, "1M": 40})
# Lead decision 26 (class q floor): the picked class q (D + 1) / T is >= 1 / T, so a class needs
# >= 1 / alpha = 100 targets to call anything at alpha 0.01; M7-13 G3/G4 plant this many per
# class at 1M (intergenic saturates at the usable-gap supply, about 100-145; reported in
# truth/novel_orfs.json, not forced).
NOVEL_PER_CLASS_GATE = 300
#: (mu, sigma) of the planted starts' log TIR in ``novel_orf_power`` (decision 14).
MATCHED_TIR = (8.0, 0.8)
_ISOBARIC_WEIGHTS = (("block", 0.4), ("near_block", 0.2), ("modification", 0.2), ("swap", 0.2))


def _length(size: str) -> int:
    if size not in SIZES:
        raise ValueError(f"unknown size {size!r}; expected one of {sorted(SIZES)}")
    return SIZES[size]


def _v2(genome_length: int, **kw: object) -> SimConfig:
    from giapsim.simulate import SimConfig

    return SimConfig(genome_length=genome_length, rbs_strength=0.8, promoter_strength=0.8,
                     terminator_strength=0.8, m7_truth=True, **kw)  # type: ignore[arg-type]  # fmt: skip


def variant_f1(seed: int, size: Size = "100k") -> SimConfig:
    """Strain differences for the ALT-classification F1 (§7.5.6); ``seed`` is not used."""
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.strain import StrainConfig

    strain = StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=0.2)
    return _v2(_length(size), strain=strain, proteomics=ProteomicsConfig(), denovo=DenovoConfig())


def novel_orf_power(seed: int, size: Size = "100k", n_per_class: int | None = None) -> SimConfig:
    """Planted and found novel ORFs for the class-FDR / power gate (§7.6.6, §8.4); ``seed`` is
    not used. Planted RBS TIRs follow the first-round distribution (module doc). ``n_per_class``
    (lead decision 26) replaces ``NOVEL_PER_CLASS[size]`` for all four classes; M7-13 G3/G4 use
    ``NOVEL_PER_CLASS_GATE`` at 1M. ``None`` keeps the default (byte-identical outputs)."""
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig

    length = _length(size)
    n = NOVEL_PER_CLASS[size] if n_per_class is None else n_per_class
    return _v2(
        length,
        n_novel_orfs=(("intergenic", n), ("cds_outframe", n), ("antisense", n)),
        n_nterm_extensions=n,
        novel_orf_len_aa=(30, 80),
        plant_novel_rbs=True,
        planted_tir=MATCHED_TIR,
        write_precomputed=True,
        proteomics=ProteomicsConfig(pi_wrong=0.02, wrong_source="sixframe_noncoding"),
        denovo=DenovoConfig(),
    )


def same_strain_insilico(seed: int) -> SimConfig:
    """TE2's in-silico same-strain peptide set (§11.6) on a 100-kb genome; ``seed`` is not used."""
    from giapsim.proteomics import ProteomicsConfig

    prot = ProteomicsConfig(pi_wrong=0.0, detect_logit=(50.0, 0.0, 0.0), search_db="same")
    return _v2(SIZES["100k"], n_sorfs=3, proteomics=prot)


def isobaric_attack(seed: int) -> SimConfig:
    """§8.5 attack 6: isobaric de novo errors only, no strain (100 kb); ``seed`` is not used."""
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig

    dn = DenovoConfig(error_weights=dict(_ISOBARIC_WEIGHTS))
    return _v2(SIZES["100k"], proteomics=ProteomicsConfig(), denovo=dn)


# ---------------------------------------------------------------- MVP-04: the K-12 evidence run

#: The one fixed seed of the MVP-04 acceptance run and the weekly slot (lead decision 1, Q11).
K12_EVIDENCE_SEED = 20261023
#: Detection intercept b0 of the K-12 evidence run (rule 7; see ``k12_evidence``).
K12_B0 = -5.25
#: Six planted islands (step 3): (first TU, effect in stationary SDs), 8 TUs each, alternately up
#: and down, spread over the chromosome; valid for the RegulonDB (2,601 TUs) and the gap-rule
#: (2,862 TUs) partitions of the recorded annotation.
K12_ISLANDS = ((150, 3.0), (560, -3.0), (970, 3.0), (1380, -3.0), (1790, 3.0), (2200, -3.0))


def k12_evidence() -> RealGenomeConfig:
    """The MVP-04 K-12 evidence inputs (spec "giapsim on a real genome").

    SNVs only at 0.5 / kb (Q4); RegulonDB 14.5.0 TUs with the gap-rule fallback (Q3); TU-level
    AR(1) at rho = 0.5 with the six ``K12_ISLANDS`` (+-3 stationary SDs, 8 TUs each);
    ``NOVEL_PER_CLASS_GATE`` = 300 novel ORFs of 30-150 aa per class (Q5); 3 runs of the related
    sample (``search_db = "other"``) at pi_wrong = 0.01, with scored null decoys
    (``decoys = "scored_null"``, MVP-04 step 6).

    Rule 7: only the intercept b0 of the M3-19 ``detect_logit`` is tuned; b1 = 0.08, b2 = 0.9
    stay. **b0 = -5.25**: truth-detected share 0.631 (2,757 of 4,367 eligible protein genes with
    >= 1 true peptide; RegulonDB TUs, ``K12_EVIDENCE_SEED``), inside the [0.55, 0.70] band. A
    giapsim-only grid on the full genome (Slurm 1870316, 1870329) measured
    b0 -2.0: 0.944 (the default saturates, example 8), -3.0: 0.894, -4.0: 0.807, -4.5: 0.741,
    -5.0: 0.669, -5.2: 0.638, -5.25: 0.631, -5.5: 0.592, -6.0: 0.503. At -5.25 the novel classes
    have 203 intergenic, 217 cds_outframe and 213 antisense truth ORFs with a true peptide.
    """
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.real_genome import RealGenomeConfig
    from giapsim.strain import StrainConfig

    _, b1, b2 = ProteomicsConfig().detect_logit
    prot = ProteomicsConfig(detect_logit=(K12_B0, b1, b2), pi_wrong=0.01, search_db="other",
                            n_runs=3, decoys="scored_null")  # fmt: skip
    expr = ExpressionConfig(
        ar1_rho=0.5,
        islands=tuple(IslandSpec("U00096.3", first, 8, eff) for first, eff in K12_ISLANDS),
    )
    return RealGenomeConfig(
        proteomics=prot,
        strain=StrainConfig(snv_rate_per_kb=0.5, indel_rate_per_kb=0.0),
        expression=expr,
        tu_source="regulondb",
        n_novel_per_class=NOVEL_PER_CLASS_GATE,
        novel_len_aa=(30, 150),
    )


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
