"""Inputs for the M4-15 gene-table tests (docs/specs/M4-15-gene-table-assembly.md).

``world`` builds the hand-made M4-06 contract frames of ``tests/unit/genetable/m4_contract_fixture``
(20 genes on circular ``chr1``; gene 3 a sORF, gene 7 a pseudogene, prophage on genes 10-19, AMR on
genes 2 and 19) and runs the merged M4-10 ... M4-14 producers on them, so ``assemble_genes`` sees
exactly what its producers return. Nothing here re-implements assembly.

``COG_DEF`` is a tiny stand-in for the COG 2020 ``cog-20.def.tab`` (tab-separated, no header:
COG ID, functional category letters, name, gene, pathway, PubMed, PDB); the tests never read the
real ``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "unit" / "genetable"))
import m4_contract_fixture as mf  # noqa: E402

NAN = math.nan

# §9.3 tables/genes.tsv, in order (identity, context, regulatory/sequence, evidence, abundance,
# function, flags).
IDENTITY = (
    "locus_tag", "source_id", "gene", "product", "feature_type", "hypothetical", "contig", "start",
    "end", "strand", "length_aa", "partial", "crosses_origin", "replicon_type",
)  # fmt: skip
CONTEXT = (
    "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids", "island_id",
    "replichore", "replichore_position", "gc_local_dev",
)  # fmt: skip
REGULATORY = (
    "promoter_score", "promoter_score_tu", "tx_rate_max_tu", "tx_rate_sum_tu",
    "terminator_conf_downstream", "tir_ostir", "sd_dg", "mrna_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10",
    "mw_da", "pi", "gravy", "tm_segments", "n_observable_peptides",
)  # fmt: skip
EVIDENCE = (
    "detected", "detection_tier", "detection_status", "gene_q_value", "detected_2pep",
    "n_peptides_unique", "n_peptides_multi_identical", "n_peptides_shared", "n_peptides_tolerant",
    "coverage_frac", "best_peptide_q", "p_detect",
)  # fmt: skip
ABUNDANCE = (
    "abundance_raw", "abundance_log2", "abundance_imputed", "imputation_policy", "abundance_type",
    "abundance_scale", "top3_log2",
)  # fmt: skip
FUNCTION = ("cog_id", "cog_category", "go_terms", "kegg_ko", "ec", "annotation_source")
GENES_COLUMNS = (*IDENTITY, *CONTEXT, *REGULATORY, *EVIDENCE, *ABUNDANCE, *FUNCTION, "flags")

# §9.4 scope Gene, in vocabulary order.
GENE_FLAGS = (
    "partial", "crosses_origin", "pseudogene", "multi_identical", "shared_evidence",
    "two_proteins_one_gene", "cross_strain_protein", "imputed_abundance", "tolerant_evidence",
    "strain_difference", "nterm_extension_evidence", "cai_reference", "top3_lt3",
    "detectability_constant", "abundance_not_explained_by_sequence",
)  # fmt: skip

LIST_COLUMNS = ("region_ids", "go_terms", "kegg_ko", "ec", "flags")
INT_COLUMNS = (
    "start", "end", "length_aa", "tu_position", "tu_size", "spacer_bp", "rare_codons_2_10",
    "tm_segments", "n_observable_peptides", "n_peptides_unique", "n_peptides_multi_identical",
    "n_peptides_shared", "n_peptides_tolerant",
)  # fmt: skip
BOOL_COLUMNS = (
    "hypothetical", "partial", "crosses_origin", "detected", "detected_2pep", "abundance_imputed",
)  # fmt: skip
FLOAT_COLUMNS = (
    "uniop_p_prev", "uniop_p_next", "gc_local_dev", "promoter_score", "promoter_score_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "terminator_conf_downstream", "tir_ostir", "sd_dg",
    "mrna_dg", "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "mw_da", "pi", "gravy",
    "gene_q_value", "coverage_frac", "best_peptide_q", "p_detect", "abundance_raw",
    "abundance_log2", "top3_log2",
)  # fmt: skip
# Columns whose §9.3 type is not stated (M5 fills them); never pinned.
UNTYPED = ("island_id", "replichore", "replichore_position")
STRING_COLUMNS = tuple(
    c
    for c in GENES_COLUMNS
    if c not in (*LIST_COLUMNS, *INT_COLUMNS, *BOOL_COLUMNS, *FLOAT_COLUMNS, *UNTYPED)
)
ENUMS: dict[str, tuple[str, ...]] = {
    "feature_type": ("cds", "sorf", "pseudogene"),
    "detection_tier": ("decoy_picked", "mayu", "two_peptide", "protein_list"),
    "detection_status": ("detected", "not_detected", "single_peptide_unestimated",
                         "partial_excluded"),
    "abundance_type": ("ibaq", "top3", "maxlfq", "intensity", "unknown", "none"),
    "abundance_scale": ("molar", "relative"),
    "annotation_source": ("bakta_full", "bakta_light", "user"),
}  # fmt: skip

COG_DEF_ROWS = (
    ("COG0593", "L", "Chromosomal replication initiation ATPase DnaA", "DnaA", "", "", ""),
    ("COG0000", "J", "Translation factor", "TfX", "", "", ""),
    ("COG0002", "KT", "Signal transduction regulator", "RegX", "", "", ""),
)


def tag(i: int) -> str:
    return mf.tag(i)


def cfg(**kw: Any) -> Any:
    from giap.config import GiapConfig

    return GiapConfig(**kw)


def cog_def_frame() -> pd.DataFrame:
    """The in-memory COG table passed to ``function_columns`` (columns ``cog_id``,
    ``cog_category``; extra columns allowed)."""
    return pd.DataFrame(
        {"cog_id": [r[0] for r in COG_DEF_ROWS], "cog_category": [r[1] for r in COG_DEF_ROWS],
         "name": [r[2] for r in COG_DEF_ROWS]}
    )  # fmt: skip


def write_cog_def(path: Path) -> Path:
    """``cog-20.def.tab`` in the NCBI layout (no header, 7 tab-separated columns)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join("\t".join(r) + "\n" for r in COG_DEF_ROWS), encoding="utf-8")
    return path


@dataclass
class World:
    """Contract frames, ``M4Inputs`` and the producer outputs ``assemble_genes`` takes."""

    frames: dict[str, Any]
    inputs: Any
    cov: pd.DataFrame
    evidence: pd.DataFrame
    detection: Any
    detectability: Any
    abundance: Any
    cfg: Any

    def assemble(self, **override: Any) -> pd.DataFrame:
        from giap.genetable.assemble import assemble_genes

        parts = {
            "inputs": self.inputs, "cov": self.cov, "evidence": self.evidence,
            "detection": self.detection, "detectability": self.detectability,
            "abundance": self.abundance,
        } | override  # fmt: skip
        return assemble_genes(
            parts["inputs"], parts["cov"], parts["evidence"], parts["detection"],
            parts["detectability"], parts["abundance"],
        )  # fmt: skip

    def with_inputs(self, **frames: Any) -> Any:
        return dataclasses.replace(self.inputs, **frames)


def detectability_result(
    p: dict[str, float], constant: bool = False, tags: list[str] | None = None
) -> Any:
    """An M4-13 ``DetectabilityResult`` with the given ``p_detect`` (NaN for unlisted tags)."""
    from giap.genetable.detectability import DetectabilityResult

    idx = list(tags if tags is not None else p)
    s = pd.Series([p.get(t, NAN) for t in idx], index=idx, dtype=float, name="p_detect")
    diag = pd.DataFrame(columns=["section", "group_type", "group", "metric", "value", "ci_low",
                                 "ci_high"])  # fmt: skip
    return DetectabilityResult(s, constant, {}, (), diag)


def build_world(
    n_genes: int = 20,
    frames: dict[str, Any] | None = None,
    constant: bool = False,
    config: Any = None,
) -> World:
    """Run M4-10 ... M4-14 on the contract frames (hand fixture unless ``frames`` given).

    ``p_detect`` is 0.5 + 0.01 i for eligible (non-partial) gene i and NaN otherwise
    (``constant`` sets ``DetectabilityResult.constant``).
    """
    from giap.genetable.abundance import gene_abundance
    from giap.genetable.covariates import covariate_table
    from giap.genetable.detection import detect
    from giap.genetable.evidence import accepted_peptides, gene_evidence

    fr = frames if frames is not None else mf.build_frames(n_genes)
    inputs = mf.make_inputs(fr)
    c = config if config is not None else cfg()
    cov = covariate_table(inputs.genes, c)
    acc = accepted_peptides(inputs.peptide_hits, inputs.peptides, c)
    ev = gene_evidence(inputs.genes, inputs.peptide_hits, acc, inputs.replicons)
    det = detect(inputs.genes, ev, acc, inputs.peptides, inputs.peptide_hits,
                 inputs.protein_hits, cov, c)  # fmt: skip
    ab = gene_abundance(inputs.genes, acc, cov, det, inputs.abundance, det.protein_gene_map, c,
                        c10_scale=inputs.abundance_scale,
                        protein_hits=inputs.protein_hits)  # fmt: skip
    tags = inputs.genes["locus_tag"].astype(str).tolist()
    partial = inputs.genes["partial"].astype(bool).tolist()
    p = {t: 0.5 + 0.01 * i for i, (t, pa) in enumerate(zip(tags, partial, strict=True)) if not pa}
    dres = detectability_result(p, constant, tags)
    return World(fr, inputs, cov, ev, det, dres, ab, c)


def frames_with(n_genes: int = 20) -> dict[str, Any]:
    """A fresh copy of the hand fixture frames (callers edit them before ``build_world``)."""
    return mf.build_frames(n_genes)


def row(df: pd.DataFrame, locus_tag: str) -> pd.Series:
    sel = df[df["locus_tag"] == locus_tag]
    assert len(sel) == 1, f"{locus_tag}: {len(sel)} rows"
    return sel.iloc[0]


def is_na(x: Any) -> bool:
    if x is None or x is pd.NA:
        return True
    if isinstance(x, list | tuple | np.ndarray):
        return False
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


def as_list(x: Any) -> list[Any]:
    """A list cell as a Python list (NA -> [])."""
    if isinstance(x, list | tuple | np.ndarray):
        return list(x)
    if is_na(x):
        return []
    raise AssertionError(f"not a list cell: {x!r}")


def sort_key(df: pd.DataFrame, contig_order: list[str]) -> list[tuple[Any, ...]]:
    """§6.5 sort key per row: contig (C1 order), start, end, strand, locus_tag."""
    rank = {c: i for i, c in enumerate(contig_order)}
    return [
        (rank[c], int(s), int(e), str(st), str(t))
        for c, s, e, st, t in zip(df["contig"], df["start"], df["end"], df["strand"],
                                  df["locus_tag"], strict=True)
    ]  # fmt: skip
