"""The ``variants`` step (GIPPA_SPEC.md §6.2, §7.5, §9.1;
docs/specs/M7-04-variants-step-and-run.md).

Module 5 (core): ``peptide_hits.tsv`` of the LOCI directory ``mapping`` (read with the M3-16
contract reader, so 0-based half-open) → M7-01 ``align_hits`` → M7-02 ``eligible``/``pileup``/
``sites`` → M7-03 ``vcf_records`` → ``pileup`` (always gzip), ``variants`` and ``variants_vcf``.
Without ``hits_tolerant`` frameshift rows are left out with the WARNING run flag
``layer_missing:hits_tolerant``. With the optional C17 ``peptide_spectra`` (M7-10) every row
carries its class error rate (``gippa.variants.error_rates``), else ``unestimated``; the rates never
change a call. ``reference_peptides`` is declared only (preflight pairs it into C17). One pass,
no randomness. The VCF sample is ``cfg.sample_name``, which ``gippa run`` resolves; an unresolved
``auto`` fails the step (``write_vcf``). M7-14: the work Parquet ``variants_qc`` holds the §7.5.1
``id_mixups`` count (one row: ``metric``, ``value``, ``n``), which ``context_stats`` writes into
``detection_summary.tsv``. M7-18: when the pileup skipped any op on a B, J, O, U, X or Z
peptide residue, one INFO event ``variants.nonstandard_residues`` carries
``attrs["nonstandard"]``.
"""

from __future__ import annotations

import hashlib
from collections.abc import Collection, Mapping
from pathlib import Path
from typing import Literal

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

import gippa
from gippa.config import GippaConfig
from gippa.io.fasta import read_genome_fasta
from gippa.io.output_checks import check_tsv, check_vcf
from gippa.io.parquet import write_parquet_deterministic
from gippa.io.variant_tables import write_pileup, write_variants
from gippa.io.vcf import VcfRecord, write_vcf
from gippa.io.work_tables import load_tables
from gippa.mapping.readthrough import annotated_stops, site_keys
from gippa.model.genes import GeneModel
from gippa.steps.mapping import _path, read_gene_model, resolve_search_db
from gippa.steps.registry import InputPath, StepContext, step
from gippa.variants.align import HitAlignment, align_hits
from gippa.variants.error_rates import alt_support, annotate_sites, class_rates, paired_observations
from gippa.variants.flags import flags_for, left_block_site, start_codons
from gippa.variants.pileup import Support, eligible, pileup, sites
from gippa.variants.qc import id_mixups, site_support
from gippa.variants.vcf_records import FORMAT_FIELDS, INFO_FIELDS, vcf_records

LAYER_MISSING_HITS_TOLERANT = "layer_missing:hits_tolerant"
SUMMARY_EVENT = "variants.summary"
NONSTANDARD_EVENT = "variants.nonstandard_residues"


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
    cfg: GippaConfig,
    spectra: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[VcfRecord]]:
    """Rule 2's chain: the pileup (with its ``attrs``), the ``variants.tsv`` frame and the VCF
    records. ``peptide_hits`` is 0-based half-open with ``;``-joined ``diff_types``. Each site
    carries its M7-05 flags and its M7-10 class error rate from the C17 ``spectra`` (none:
    ``unestimated``); neither changes a call, confidence or count. Pure."""
    pile, table, records, _ = run_variants_qc(peptide_hits, hits_tolerant, peptides, gm, genome,
                                              cfg, spectra)  # fmt: skip
    return pile, table, records


def run_variants_qc(
    peptide_hits: pd.DataFrame,
    hits_tolerant: pd.DataFrame | None,
    peptides: pd.DataFrame,
    gm: GeneModel,
    genome: Mapping[str, str],
    cfg: GippaConfig,
    spectra: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[VcfRecord], pd.DataFrame]:
    """:func:`run_variants` plus the long ``(variant_id, peptide, source, search_db)`` support
    frame of each ``variants.tsv`` row (M7-14, :func:`gippa.variants.qc.id_mixups`). Pure."""
    support = support_of(peptide_hits, peptides, {g.locus_tag for g in gm.genes})
    stops = site_keys(annotated_stops(gm, genome, recoded=True)[0])  # M7-16 rules 10, 22
    aligns = align_hits(peptide_hits, hits_tolerant, genome, gm.replicons, stop_sites=stops)
    pile = pileup(eligible(aligns, support, cfg), gm, genome, cfg, support=support)
    st = sites(pile, cfg)
    starts = start_codons(gm)
    alt_peps = st.attrs["alt_peptides"]
    by_pep: dict[str, list[HitAlignment]] = {}
    for a in aligns:
        by_pep.setdefault(a.peptide, []).append(a)
    tol = cfg.variants.mass_equivalent_tol_da
    flags = []
    for r, peps in zip(st.to_dict("records"), alt_peps, strict=True):
        alt_aligns = [a for p in sorted(peps) for a in by_pep.get(p, ())]
        left = left_block_site(r, alt_aligns, genome, gm.replicons, tol)
        alt_src = {(support[p].source, support[p].search_db) for p in peps}
        flags.append(tuple(flags_for(r, starts, cfg, alt_src, left_block=left)))
    st["flags"] = flags
    rates = []
    if spectra is not None:
        sup = alt_support(st, alt_peps)
        is_dn = spectra["source"].astype(str).eq("denovo").to_numpy(dtype=bool)
        rates = class_rates(paired_observations(spectra[is_dn], spectra[~is_dn], sup))
    st = annotate_sites(st, rates)
    records, table = vcf_records(st, pile, genome, gm.replicons)
    sources = {p: (str(s.source), str(s.search_db)) for p, s in support.items()}
    return pile, table, records, site_support(table["variant_id"], alt_peps, sources)


def qc_frame(table: pd.DataFrame, support: pd.DataFrame) -> pd.DataFrame:
    """The ``variants_qc`` rows (``metric``, ``value``, ``n``): ``id_mixups`` (rule 3). Pure."""
    count, n = id_mixups(table, support)
    return pd.DataFrame({"metric": ["id_mixups"], "value": [count], "n": [n]})


def write_variants_qc(path: Path, qc: pd.DataFrame) -> None:
    """``work/variants_qc`` (byte-reproducible Parquet, TE7)."""
    cols = [("metric", pa.string()), ("value", pa.int64()), ("n", pa.int64())]
    rows = [(str(m), int(v), int(n)) for m, v, n in zip(qc["metric"], qc["value"], qc["n"],
                                                        strict=True)]  # fmt: skip
    write_parquet_deterministic(path, cols, rows)


def _peptide_hits(mapping: Path) -> pd.DataFrame:
    frames, issues = load_tables({"peptide_hits": mapping / "peptide_hits.tsv"}, ["peptide_hits"])
    errors = [i for i in issues if i.severity == "error"]
    if errors:
        raise ValueError(f"{mapping / 'peptide_hits.tsv'}: {errors[0]}")
    ph = frames["peptide_hits"]
    ph["diff_types"] = [";".join(v) if isinstance(v, list) else v for v in ph["diff_types"]]
    return ph


def _spectra(path: InputPath | None) -> pd.DataFrame | None:
    """C17 ``peptide_spectra`` through the contract reader (M7-10); None when not staged."""
    if path is None:
        return None
    frames, issues = load_tables({"peptide_spectra": _path(path)}, ["peptide_spectra"])
    errors = [i for i in issues if i.severity == "error"]
    if errors:
        raise ValueError(f"{path}: {errors[0]}")
    return frames["peptide_spectra"]


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
    optional_inputs=["hits_tolerant", "reference_peptides", "peptide_spectra"],
    outputs=["pileup", "variants", "variants_vcf", "variants_qc"],
)
def variants_step(ctx: StepContext) -> None:
    """Codon pileup, ``variants.tsv``, ``variants.vcf`` and ``variants_qc`` (Module 5)."""
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
    spectra = _spectra(ctx.inputs.get("peptide_spectra"))
    pile, table, records, support = run_variants_qc(ph, ht, peptides, gm, genome, ctx.cfg, spectra)
    data = {
        "aligned_hits": len(ph),
        "covered_codons": len(pile),
        "sites_by_call": _counts(table["call"]),
        "sites_by_confidence": _counts(table["confidence"]),
    }
    ctx.log.info("variants: %d aligned hits, %d covered codons, %d sites (by call %s)",
                 len(ph), len(pile), len(table), data["sites_by_call"],
                 extra={"event": SUMMARY_EVENT, "data": data})  # fmt: skip
    ns = pile.attrs["nonstandard"]
    if ns["skipped_sites"]:
        ctx.log.info("variants: %d site(s) with an ambiguous or rare peptide residue skipped (%s)",
                     ns["skipped_sites"], ns["by_residue"],
                     extra={"event": NONSTANDARD_EVENT, "data": ns})  # fmt: skip
    sha = hashlib.sha256(gpath.read_bytes()).hexdigest()
    write_vcf(ctx.outputs["variants_vcf"], records, list(gm.replicons.values()),
              gippa_version=gippa.__version__, genome_sha256=sha, sample=ctx.cfg.sample_name,
              info_fields=INFO_FIELDS, format_fields=FORMAT_FIELDS, sequences=genome)  # fmt: skip
    _write_gzip(ctx.outputs["pileup"], pile)
    write_variants(ctx.outputs["variants"], table)
    write_variants_qc(ctx.outputs["variants_qc"], qc_frame(table, support))
    check_vcf("variants_vcf", ctx.outputs["variants_vcf"])  # M7-15 item 2
    check_tsv("pileup", ctx.outputs["pileup"], "pileup")
    check_tsv("variants", ctx.outputs["variants"], "variants")
