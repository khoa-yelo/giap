"""Step ``gene_table`` (GIPPA_SPEC.md §6.2 "genes, layers, loci, abundance → genes_pre.parquet";
docs/specs/M4-15-gene-table-assembly.md): the M4-06 contract tables through M4-10 ... M4-14 into
the §9.3 gene table without ``p_detect``.

Inputs: M2's Parquet intermediates (C1-C6; not C7/C8) and ``annotation_meta``, ``mapping`` (the
``loci`` directory with ``peptide_hits.tsv``, ``protein_hits.tsv``, ``unmapped.tsv``), optional
C9/C10, the C14 id map, preflight.json and ``cog2020_def`` (raw NCBI ``cog-20.def.tab``; absent
→ ``cog_category`` NA and run flag ``layer_missing:cog``).
A contract ``error`` issue fails the step (``StepError``).
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from gippa.config import GippaConfig
from gippa.genetable.abundance import gene_abundance
from gippa.genetable.assemble import assemble_genes
from gippa.genetable.covariates import covariate_table
from gippa.genetable.detection import detect
from gippa.genetable.evidence import accepted_peptides, gene_evidence
from gippa.io.genes_tsv import read_cog_def, write_genes_parquet
from gippa.io.issues import Issue
from gippa.io.summary_tsv import write_summary_part
from gippa.io.work_tables import FILE_NAMES, OPTIONAL_KEYS, default_paths, load_m4_inputs
from gippa.setupdb.installed import cog_table
from gippa.steps.registry import InputPath, StepContext, run_step, step

__all__ = ["GenesPre", "build_genes_pre", "contract_paths", "m4_paths", "run_gene_table",
           "work_inputs"]  # fmt: skip

_M3_TSVS = ("peptide_hits", "protein_hits", "unmapped")
# C7/C8 are not read: region membership arrives through C6 ``region_ids``.
_REQUIRED = ("genes", "replicons", "xrefs", "annotation_meta", "tus", "layers", "mapping")
_UNREAD = frozenset({"contract.regions.missing_file", "contract.region_genes.missing_file"})
_OPTIONAL = ("peptides", "abundance", "id_map", "preflight", "cog2020_def")


def _one(p: InputPath) -> Path:
    return p[0] if isinstance(p, list) else p


def m4_paths(inputs: Mapping[str, InputPath]) -> dict[str, Path]:
    """Step input keys → ``load_m4_inputs`` paths (the M3 TSVs live in the ``mapping`` dir)."""
    paths = {k: _one(v) for k, v in inputs.items() if k in FILE_NAMES or k == "preflight"}
    if "preflight" in paths and paths["preflight"].is_dir():  # local executor: work/preflight/
        paths["preflight"] = paths["preflight"] / "preflight.json"
    if "mapping" in inputs:
        mapping = _one(inputs["mapping"])
        paths |= {k: mapping / FILE_NAMES[k] for k in _M3_TSVS}
    return paths


def check_issues(issues: list[Issue], log: logging.Logger, what: str) -> None:
    """Log warnings; raise ``ValueError`` listing every contract ``error`` issue."""
    for i in issues:
        if i.severity != "error":
            log.warning("%s: %s: %s", what, i.code, i.message)
    errors = [f"{i.code}: {i.message}" for i in issues if i.severity == "error"]
    if errors:
        raise ValueError(f"{what}: M4 input contract errors: " + "; ".join(errors))


@dataclass(frozen=True, slots=True)
class GenesPre:
    """What the step writes: the gene table without ``p_detect`` and M4-12's ``tier`` rows."""

    table: pd.DataFrame
    tier_summary: pd.DataFrame


def build_genes_pre(
    inputs: Mapping[str, InputPath], cfg: GippaConfig, log: logging.Logger
) -> GenesPre:
    """Load and check the contract tables, run M4-10 ... M4-14 and assemble (no ``p_detect``)."""
    m4, issues = load_m4_inputs(m4_paths(inputs))
    check_issues([i for i in issues if i.code not in _UNREAD], log, "gene_table")
    cog = inputs.get("cog2020_def")
    if cog is None:
        log.warning(
            "layer_missing:cog: COG 2020 table not installed (gippa setup --cog); cog_category NA",
            extra={"event": "gene_table.layer_missing", "data": {"layer": "cog",
                                                                 "flag": "layer_missing:cog"}},
        )  # fmt: skip
    cog_def = None if cog is None else read_cog_def(_one(cog))
    cov = covariate_table(m4.genes, cfg)
    acc = accepted_peptides(m4.peptide_hits, m4.peptides, cfg)
    ev = gene_evidence(m4.genes, m4.peptide_hits, acc, m4.replicons)
    det = detect(m4.genes, ev, acc, m4.peptides, m4.peptide_hits, m4.protein_hits, cov, cfg)
    for w in det.warnings:
        log.warning("gene_table: detection: %s", w,
                    extra={"event": f"gene_table.{w}", "data": {"warning": w}})  # fmt: skip
    ab = gene_abundance(m4.genes, acc, cov, det, m4.abundance, det.protein_gene_map, cfg,
                        c10_scale=m4.abundance_scale, protein_hits=m4.protein_hits)  # fmt: skip
    table = assemble_genes(m4, cov, ev, det, None, ab, cog_def=cog_def)
    return GenesPre(table.drop(columns=["p_detect"]), det.summary_rows)


@step(
    "gene_table",
    inputs=_REQUIRED,
    optional_inputs=_OPTIONAL,
    outputs=["genes_pre", "tier_summary"],
)
def gene_table(ctx: StepContext) -> None:
    """Contract tables → ``genetable/genes_pre.parquet`` (everything but ``p_detect``) and the
    ``tier`` rows of ``detection_summary``."""
    pre = build_genes_pre(ctx.inputs, ctx.cfg, ctx.log)
    write_genes_parquet(ctx.outputs["genes_pre"], pre.table)
    write_summary_part(ctx.outputs["tier_summary"], pre.tier_summary)
    n_det = int(pre.table["detected"].sum())
    ctx.log.info("gene_table: %d genes, %d detected, %d tier rows", len(pre.table), n_det,
                 len(pre.tier_summary))  # fmt: skip


def _first(*paths: Path) -> Path | None:
    return next((p for p in paths if p.exists()), None)


def work_inputs(work_dir: Path) -> dict[str, Path]:
    """Step inputs found in ``work_dir``: the ``gippa run`` layout (M3 TSVs under ``mapping/``, C14
    under ``annotation/``) or the flat ``default_paths`` layout (lead decision 8)."""
    found: dict[str, Path | None] = {
        k: _first(work_dir / FILE_NAMES[k], work_dir / k)
        for k in _REQUIRED + _OPTIONAL
        if k in FILE_NAMES
    }
    found["id_map"] = _first(work_dir / "annotation" / "id_map.tsv", work_dir / "id_map.tsv")
    found["preflight"] = _first(work_dir / "preflight.json")
    flat = all((work_dir / FILE_NAMES[k]).is_file() for k in _M3_TSVS)
    found["mapping"] = work_dir if flat else work_dir / "mapping"
    db = os.environ.get("GIPPA_DB_DIR")
    found["cog2020_def"] = cog_table(Path(db) if db else None)
    return {k: p for k, p in found.items() if p is not None and p.exists()}


def contract_paths(work_dir: Path) -> dict[str, Path]:
    """M4-06 key -> file in ``work_dir`` for every contract table, in either layout of
    ``work_inputs`` (M4-21 lead decision 3; C7/C8 as ``regions.parquet`` or the key itself). A
    required table that is not there keeps its ``default_paths`` name, so ``load_m4_inputs``
    reports it as missing; an optional one is left out."""
    found = m4_paths(work_inputs(work_dir))
    for key in ("regions", "region_genes"):
        p = _first(work_dir / FILE_NAMES[key], work_dir / key)
        if p is not None:
            found[key] = p
    paths = {k: p for k, p in default_paths(work_dir).items() if k not in OPTIONAL_KEYS}
    return paths | {k: p for k, p in found.items() if k in FILE_NAMES}


def run_gene_table(work_dir: Path, cfg: GippaConfig) -> Path:
    """Run ``gene_table`` on ``work_dir``; returns ``work_dir/genetable/genes_pre.parquet``
    (the ``tier`` rows go to ``genetable/detection_summary_parts/tier.parquet``)."""
    gt = work_dir / "genetable"
    out = gt / "genes_pre.parquet"
    outputs = {"genes_pre": out, "tier_summary": gt / "detection_summary_parts" / "tier.parquet"}
    run_step("gene_table", cfg, work_inputs(work_dir), outputs)
    return out
