"""The ``variants`` step (GIAP_SPEC.md §6.2, §7.5, §9.1; docs/specs/M7-04-variants-step-and-run.md).

Module 5 (core): ``peptide_hits.tsv`` of the LOCI directory ``mapping`` (read with the M3-16
contract reader, so 0-based half-open) → M7-01 ``align_hits`` → M7-02 ``eligible``/``pileup``/
``sites`` → M7-03 ``vcf_records`` → ``pileup`` (always gzip), ``variants`` and ``variants_vcf``.
Without ``hits_tolerant`` frameshift rows are left out with the WARNING run flag
``layer_missing:hits_tolerant``. ``reference_peptides`` is declared for M7-10 and not read yet.
One pass, no randomness. The VCF sample is ``cfg.sample_name``, which ``giap run`` resolves; an
unresolved ``auto`` fails the step (``write_vcf``).
"""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Mapping
from pathlib import Path
from typing import Literal

import pandas as pd

import giap
from giap.config import GiapConfig
from giap.io.fasta import read_genome_fasta
from giap.io.output_checks import check_tsv, check_vcf
from giap.io.variant_tables import write_pileup, write_variants
from giap.io.vcf import VcfRecord, write_vcf
from giap.io.work_tables import load_tables
from giap.model.genes import GeneModel
from giap.steps.mapping import _path, read_gene_model, resolve_search_db
from giap.steps.registry import StepContext, step
from giap.variants.align import align_hits
from giap.variants.flags import flags_for, start_codons
from giap.variants.pileup import Support, eligible, pileup, sites
from giap.variants.vcf_records import FORMAT_FIELDS, INFO_FIELDS, vcf_records

LAYER_MISSING_HITS_TOLERANT = "layer_missing:hits_tolerant"
SUMMARY_EVENT = "variants.summary"


def _has_frameshift(v: object) -> bool:
    return isinstance(v, str) and "frameshift" in v.split(";")


def without_frameshifts(peptide_hits: pd.DataFrame) -> pd.DataFrame:
    """``peptide_hits`` minus its frameshift rows (rule 1, no ``hits_tolerant``). Pure."""
    keep = ~peptide_hits["diff_types"].map(_has_frameshift).astype(bool)
    return peptide_hits[keep].reset_index(drop=True)


def support_of(
    peptide_hits: pd.DataFrame, peptides: pd.DataFrame, database_ids: Collection[str]
) -> dict[str, Support]:
    """Per peptide: C11 ``source``/``placement_q``/``uniqueness`` and C9 ``search_db``. A value
    other than ``same``/``other`` (``auto``), or none, is resolved table-wide as placement does
    (``resolve_search_db`` against ``database_ids``, the gene model's locus tags). Pure."""
    sdb: dict[str, str] = {}
    if "search_db" in peptides.columns:
        sdb = {str(p): str(s) for p, s in zip(peptides["peptide"], peptides["search_db"],
                                               strict=True)}  # fmt: skip
    resolved: str | None = None
    out: dict[str, Support] = {}
    for r in peptide_hits.to_dict("records"):
        s = sdb.get(r["peptide"])
        if s not in ("same", "other"):
            if resolved is None:
                resolved = resolve_search_db(peptides, database_ids)
            s = resolved
        db: Literal["same", "other"] = "other" if s == "other" else "same"
        out.setdefault(r["peptide"], Support(r["peptide"], r["source"], db,
                       float(r["placement_q"]), r["uniqueness"]))  # fmt: skip
    return out


def run_variants(
    peptide_hits: pd.DataFrame,
    hits_tolerant: pd.DataFrame | None,
    peptides: pd.DataFrame,
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GiapConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, list[VcfRecord]]:
    """Rule 2's chain: the pileup (with its ``attrs``), the ``variants.tsv`` frame and the VCF
    records. ``peptide_hits`` is 0-based half-open with ``;``-joined ``diff_types``. Each site
    carries its M7-05 flags; they never change a call, confidence or count. Pure."""
    support = support_of(peptide_hits, peptides, {g.locus_tag for g in gm.genes})
    aligns = align_hits(peptide_hits, hits_tolerant, genome, gm.replicons)
    pile = pileup(eligible(aligns, support, cfg), gm, genome, cfg, support=support)
    st = sites(pile, cfg)
    starts = start_codons(gm)
    st["flags"] = [
        tuple(flags_for(r, starts, cfg, {(support[p].source, support[p].search_db) for p in peps}))
        for r, peps in zip(st.to_dict("records"), st.attrs["alt_peptides"], strict=True)
    ]
    records, table = vcf_records(st, pile, genome, gm.replicons)
    return pile, table, records


def _peptide_hits(mapping: Path) -> pd.DataFrame:
    frames, issues = load_tables({"peptide_hits": mapping / "peptide_hits.tsv"}, ["peptide_hits"])
    errors = [i for i in issues if i.severity == "error"]
    if errors:
        raise ValueError(f"{mapping / 'peptide_hits.tsv'}: {errors[0]}")
    ph = frames["peptide_hits"]
    ph["diff_types"] = [";".join(v) if isinstance(v, list) else v for v in ph["diff_types"]]
    return ph


def _counts(col: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in col.value_counts().sort_index().items()}


def _write_gzip(path: Path, pile: pd.DataFrame) -> None:
    if path.name.endswith(".gz"):
        write_pileup(path, pile)
        return
    gz = path.with_name(path.name + ".gz")  # the writer gzips by suffix
    write_pileup(gz, pile)
    gz.replace(path)


@step(
    "variants",
    inputs=["mapping", "genome", "genes", "replicons", "peptides"],
    optional_inputs=["hits_tolerant", "reference_peptides"],
    outputs=["pileup", "variants", "variants_vcf"],
)
def variants_step(ctx: StepContext) -> None:
    """Codon pileup, ``variants.tsv`` and ``variants.vcf`` (Module 5)."""
    ph = _peptide_hits(_path(ctx.inputs["mapping"]))
    ht_in = ctx.inputs.get("hits_tolerant")
    ht = None if ht_in is None else pd.read_parquet(_path(ht_in))
    if ht is None:
        n = len(ph)
        ph = without_frameshifts(ph)
        ctx.log.warning("%s: hits_tolerant not available; %d frameshift row(s) skipped",
                        LAYER_MISSING_HITS_TOLERANT, n - len(ph))  # fmt: skip
    gpath = _path(ctx.inputs["genome"])
    genome = {r.id: r.sequence for r in read_genome_fasta(gpath)}
    gm = read_gene_model(_path(ctx.inputs["genes"]), _path(ctx.inputs["replicons"]))
    peptides = pd.read_parquet(_path(ctx.inputs["peptides"]))
    pile, table, records = run_variants(ph, ht, peptides, gm, genome, ctx.cfg)
    data = {
        "aligned_hits": len(ph),
        "covered_codons": len(pile),
        "sites_by_call": _counts(table["call"]),
        "sites_by_confidence": _counts(table["confidence"]),
    }
    ctx.log.info("variants: %d aligned hits, %d covered codons, %d sites (by call %s)",
                 len(ph), len(pile), len(table), data["sites_by_call"],
                 extra={"event": SUMMARY_EVENT, "data": data})  # fmt: skip
    sha = hashlib.sha256(gpath.read_bytes()).hexdigest()
    write_vcf(ctx.outputs["variants_vcf"], records, list(gm.replicons.values()),
              giap_version=giap.__version__, genome_sha256=sha, sample=ctx.cfg.sample_name,
              info_fields=INFO_FIELDS, format_fields=FORMAT_FIELDS, sequences=genome)  # fmt: skip
    _write_gzip(ctx.outputs["pileup"], pile)
    write_variants(ctx.outputs["variants"], table)
    check_vcf("variants_vcf", ctx.outputs["variants_vcf"])  # M7-15 item 2
    check_tsv("pileup", ctx.outputs["pileup"], "pileup")
    check_tsv("variants", ctx.outputs["variants"], "variants")
