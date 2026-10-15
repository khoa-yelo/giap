"""Pydantic v2 configuration model (GIAP_SPEC.md §5.2).

Defaults here MUST equal ``defaults.yaml`` (enforced by tests/unit/test_config.py).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PlainSerializer,
    SerializationInfo,
    WithJsonSchema,
    field_validator,
)

Module = Literal["core", "mobile", "bgc", "function"]
GenomeRelation = Literal["related", "same"]
PeptideSource = Literal["auto", "database", "denovo"]
SearchDb = Literal["auto", "same", "other"]
PromoterMode = Literal["targeted", "genome"]
TerminatorTool = Literal["terminatornet", "transtermhp"]
AbundanceMethod = Literal["ibaq", "top3"]
AbundanceType = Literal["ibaq", "top3", "maxlfq", "intensity", "unknown"]
AbundanceScale = Literal["auto", "log2", "linear"]
AbundanceMissing = Literal["median", "exclude", "min", "zero", "mean"]
Enzyme = Literal["trypsin_p", "trypsin"]
SensitivityLevel = Literal["basic", "full", "off"]
RegionType = Literal["prophage", "plasmid", "bgc", "crispr", "amr"]
OperonDecoy = Literal["convergent", "sandwich"]
PlotFormat = Literal["svg", "png", "pdf"]


def _encode_inf(value: float, info: SerializationInfo) -> float | str:
    """In JSON mode write ±inf as the strings ``"inf"``/``"-inf"`` (strict JSON has no Infinity).

    Pydantic's float validation parses those strings back, so the dump round-trips.
    """
    if info.mode_is_json() and math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return value


ExtendedFloat = Annotated[
    float,
    PlainSerializer(_encode_inf, return_type=Any),
    WithJsonSchema({"anyOf": [{"type": "number"}, {"enum": ["inf", "-inf"]}]}),
]
"""A float that may be ±inf; encoded as ``"inf"``/``"-inf"`` in JSON."""


class _Strict(BaseModel):
    """Base for every config section: unknown keys are errors."""

    # coerce_numbers_to_str: ``--set sample_name=2024`` parses as YAML int but means "2024".
    model_config = ConfigDict(extra="forbid", coerce_numbers_to_str=True)


# ---------------------------------------------------------------- input / preflight


class InputConfig(_Strict):
    peptide_source: PeptideSource = "auto"
    search_db: SearchDb = "auto"
    decoy_prefixes: list[str] = ["rev_", "REV_", "DECOY_", "decoy_", "XXX_", "##"]
    contaminant_prefixes: list[str] = ["Cont_", "CON__", "CONT_", "cRAP"]
    q_value_max: float = 0.01
    peptide_len: list[int] = [7, 40]
    i_equals_l: bool = True


class PreflightConfig(_Strict):
    max_contigs_warn: int = 10
    min_total_bp_warn: int = 100000
    min_n50_warn: int = 100000
    min_largest_contig_warn: int = 1000000
    small_replicon_bp: int = 500000
    gc_range_warn: list[float] = [0.25, 0.75]
    n_run_min_bp: int = Field(default=10, ge=1)
    abundance_missing_warn_frac: float = 0.20


# ---------------------------------------------------------------- annotation


class BaktaConfig(_Strict):
    db: str = "auto"
    complete: bool | Literal["auto"] = "auto"
    extra_args: list[str] = []


class UniopConfig(_Strict):
    pair_threshold: float = 0.5
    fallback: Literal["distance"] = "distance"


class PromotersConfig(_Strict):
    mode: PromoterMode = "targeted"
    window_nt: int = 40
    window_3p_range: list[int] = [-200, 20]
    promotech_yes: float = 0.5


class PromoterCalculatorConfig(_Strict):
    enabled: bool = True


class TerminatorsConfig(_Strict):
    tool: TerminatorTool = "terminatornet"
    cross_check_transtermhp: bool = True
    bactermfinder: bool = False
    terminatornet_yes: float = 0.5


class OstirConfig(_Strict):
    upstream_nt: int = 50
    downstream_nt: int = 50
    anti_sd: str = "auto"


class CodonConfig(_Strict):
    cai_reference: str = "ribosomal_proteins"


class PrecomputedConfig(_Strict):
    bakta_json: Path | None = None
    operon_pairs: Path | None = None
    promoters: Path | None = None
    tx_rates: Path | None = None
    terminators: Path | None = None
    rbs: Path | None = None
    regions: Path | None = None


class AnnotationConfig(_Strict):
    bakta: BaktaConfig = Field(default_factory=BaktaConfig)
    translation_table: int = 11
    start_codons: list[str] = ["ATG", "GTG", "TTG"]
    region_min_overlap_frac: float = 0.5
    uniop: UniopConfig = Field(default_factory=UniopConfig)
    tu_upstream_bp: list[int] = [-300, 50]
    promoters: PromotersConfig = Field(default_factory=PromotersConfig)
    promoter_calculator: PromoterCalculatorConfig = Field(default_factory=PromoterCalculatorConfig)
    terminators: TerminatorsConfig = Field(default_factory=TerminatorsConfig)
    ostir: OstirConfig = Field(default_factory=OstirConfig)
    codon: CodonConfig = Field(default_factory=CodonConfig)
    precomputed: PrecomputedConfig = Field(default_factory=PrecomputedConfig)


# ---------------------------------------------------------------- mapping


class ExactConfig(_Strict):
    min_len: int = 7
    max_len: int = 40


class ProteinsConfig(_Strict):
    min_identity: float = 0.80
    min_coverage: float = 0.80
    near_identity: float = 0.95
    miniprot_args: str = "-S -N 5 --outs=0.95"
    tblastn_args: str = "-db_gencode 11 -seg no -evalue 1e-5"


class FrameshiftConfig(_Strict):
    enabled: bool = True
    min_side_aa: int = Field(default=4, ge=1)


class TolerantConfig(_Strict):
    enabled: bool = True
    min_len: int = 10
    max_diff: dict[str, int] = {"10-11": 1, "12-40": 2}
    swaps: bool = True
    mass_equivalents: bool = True
    frameshift: FrameshiftConfig = Field(default_factory=FrameshiftConfig)


class MappingDecoysConfig(_Strict):
    database_input: Literal["reverse", "shuffle"] = "reverse"
    denovo_input: Literal["reverse", "shuffle"] = "shuffle"
    shuffle_window_aa: int = 50


class MappingConfig(_Strict):
    exact: ExactConfig = Field(default_factory=ExactConfig)
    proteins: ProteinsConfig = Field(default_factory=ProteinsConfig)
    tolerant: TolerantConfig = Field(default_factory=TolerantConfig)
    decoys: MappingDecoysConfig = Field(default_factory=MappingDecoysConfig)
    length_classes: list[list[int]] = [
        [7, 7],
        [8, 8],
        [9, 9],
        [10, 11],
        [12, 14],
        [15, 19],
        [20, 40],
    ]


# ---------------------------------------------------------------- gene table


class DetectionConfig(_Strict):
    min_peptides_sequences_only: int = 2
    mayu_size_bins: int = 10


class AbundanceConfig(_Strict):
    method: AbundanceMethod = "ibaq"
    type: AbundanceType = "unknown"
    scale: AbundanceScale = "auto"
    missing: AbundanceMissing = "median"
    sensitivity: list[AbundanceMissing] = ["min", "exclude"]
    enzyme: Enzyme = "trypsin_p"
    observable_len: list[int] = [7, 30]
    run_aggregation: Literal["median"] = "median"


class DetectabilityConfig(_Strict):
    covariates: list[str] = ["n_observable_peptides", "length_aa", "gravy", "pi", "tm_segments"]
    spline_knots: int = 5
    folds: int = 5


# ---------------------------------------------------------------- statistics and modules


class AlphaConfig(_Strict):
    identification: float = 0.01
    placement: float = 0.01
    novel_orfs: float = 0.01
    variants: float = 0.01
    context: float = 0.05
    enrichment: float = 0.05
    operon_flags: float = 0.10
    determinants: float = 0.05


class StatsConfig(_Strict):
    alpha: AlphaConfig = Field(default_factory=AlphaConfig)
    permutations: int = 1000
    bootstrap: int = 1000
    besag_clifford_h: int = 10


class IslandsConfig(_Strict):
    min_tus: int = 2
    max_tus: int = 20
    detrend_span_frac: float = 0.10
    block_len: int | Literal["auto"] = "auto"


class ContextConfig(_Strict):
    region_types: list[RegionType] = ["prophage", "plasmid", "bgc", "crispr", "amr"]
    exhaustive_shift_max_genes: int = 20000
    moran_lags: list[int] = [1, 20]
    islands: IslandsConfig = Field(default_factory=IslandsConfig)


class OperonsConfig(_Strict):
    decoys: list[OperonDecoy] = ["convergent", "sandwich"]
    distance_bins_bp: list[ExtendedFloat] = Field(
        default=[-math.inf, 0, 50, 150, 300, math.inf],
        # Pydantic encodes schema defaults without field serializers (±inf would become null).
        json_schema_extra={"default": ["-inf", 0, 50, 150, 300, "inf"]},
    )
    pi0: float = 1.0


class VariantsConfig(_Strict):
    high_min_peptides: int = 2
    mass_equivalent_tol_da: float = 0.05


class NovelOrfsConfig(_Strict):
    min_codons: int = 10
    min_peptides_without_decoys: int = 2
    rbs_min_quantile: float = 0.10


class GseaConfig(_Strict):
    sample_size: int = 101
    eps: float = 0.0
    n_perm_simple: int = 1000


class FunctionConfig(_Strict):
    go_slim: str = "generic"
    ora_size: list[int] = [5, 500]
    gsea_size: list[int] = [10, 500]
    gsea: GseaConfig = Field(default_factory=GseaConfig)


class LightgbmConfig(_Strict):
    n_estimators: int = 500
    learning_rate: float = 0.03
    num_leaves: int = 15
    min_child_samples: int = 20
    subsample: float = 0.8
    subsample_freq: int = 1
    colsample_bytree: float = 0.8
    deterministic: bool = True


class ElasticNetConfig(_Strict):
    l1_ratio: list[float] = [0.1, 0.5, 0.9]
    inner_folds: int = 3


class DeterminantsConfig(_Strict):
    outer_folds: int = 5
    repeats: int = 5
    permutations: int = 1000
    lightgbm: LightgbmConfig = Field(default_factory=LightgbmConfig)
    elastic_net: ElasticNetConfig = Field(default_factory=ElasticNetConfig)


class ReportConfig(_Strict):
    interactive: bool = True
    max_locus_plots: int = 50
    formats: list[PlotFormat] = ["svg", "png", "pdf"]
    dpi: int = 300


class OutputConfig(_Strict):
    parquet: bool = False
    sensitivity: SensitivityLevel = "basic"


# ---------------------------------------------------------------- root


class GiapConfig(_Strict):
    """The resolved GIAP configuration (root of ``defaults.yaml``)."""

    schema_version: str = "1.0"
    seed: int = 20261014
    sample_name: str = "auto"
    threads: int = 16
    modules: list[Module] = ["core", "mobile", "bgc", "function"]
    genome_relation: GenomeRelation = "related"
    input: InputConfig = Field(default_factory=InputConfig)
    preflight: PreflightConfig = Field(default_factory=PreflightConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    mapping: MappingConfig = Field(default_factory=MappingConfig)
    detection: DetectionConfig = Field(default_factory=DetectionConfig)
    abundance: AbundanceConfig = Field(default_factory=AbundanceConfig)
    detectability: DetectabilityConfig = Field(default_factory=DetectabilityConfig)
    stats: StatsConfig = Field(default_factory=StatsConfig)
    context: ContextConfig = Field(default_factory=ContextConfig)
    operons: OperonsConfig = Field(default_factory=OperonsConfig)
    variants: VariantsConfig = Field(default_factory=VariantsConfig)
    novel_orfs: NovelOrfsConfig = Field(default_factory=NovelOrfsConfig)
    function: FunctionConfig = Field(default_factory=FunctionConfig)
    determinants: DeterminantsConfig = Field(default_factory=DeterminantsConfig)
    report: ReportConfig = Field(default_factory=ReportConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    @field_validator("modules")
    @classmethod
    def _core_always_on(cls, value: list[Module]) -> list[Module]:
        return value if "core" in value else ["core", *value]
