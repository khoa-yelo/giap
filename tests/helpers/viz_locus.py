"""Small synthetic ``OUT/`` trees for the V3 ``locus_<id>`` tests (M9-10).

Every builder takes 1-based closed coordinates, as the published TSVs carry them (``RunView``
converts the ``start`` columns to 0-based on read); BED rows are 0-based half-open. Tables go
through ``helpers.viz_out_tree.write_schema_table`` (schema header; absent cells ``NA``);
``annotation/tus.tsv`` and ``annotation/regions.tsv`` have no schema and use the writer column
orders of ``giap.annotate.operons.TUS_COLUMNS`` and ``giap.annotate.regions.TSV_COLUMNS``.
Nothing here touches ``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from giap.io.tsv import write_tsv
from helpers.viz_out_tree import gene_row, make_out, write_schema_table

#: chr1 100 kb circular chromosome, lin1 30 kb linear contig
REPLICONS: tuple[dict[str, Any], ...] = (
    {"id": "chr1", "kind": "chromosome", "length": 100_000, "topology": "circular"},
    {"id": "lin1", "kind": "unknown", "length": 30_000, "topology": "linear"},
)
TUS_COLUMNS = (
    "tu_id", "contig", "start", "end", "strand", "genes", "upstream_start", "upstream_end",
)  # fmt: skip
REGIONS_COLUMNS = (
    "region_id", "region_type", "contig", "start", "end", "source_tool", "score", "genes", "notes",
)  # fmt: skip


def preflight(replicons: Sequence[Mapping[str, Any]] = REPLICONS) -> dict[str, Any]:
    """A minimal ``preflight.json`` (``preflight.schema.json`` fields) carrying ``replicons``."""
    return {
        "schema_version": "1.0",
        "giap_version": "0.0.1",
        "inputs": [],
        "checks": [],
        "summary": {"n_error": 0, "n_warn": 0, "n_info": 0, "n_ok": 0},
        "evidence_tier": "decoy_picked",
        "routing": {
            "n_proteins_long": 0,
            "n_proteins_peptide_path": 0,
            "n_proteins_dropped": 0,
            "n_peptides": 0,
            "abundance_used": True,
            "abundance_scale": "log2",
        },
        "replicons": [dict(r) for r in replicons],
        "contig_map": [{"id": r["id"], "original_id": r["id"]} for r in replicons],
    }


def region_stat(
    region_id: str,
    region_type: str,
    contig: str,
    start: int,
    end: int,
    q: float | None,
    *,
    stat: str = "o_minus_e",
    subset: str = "all",
    n_eligible: int = 100,
    observed: float = 12.0,
    expected: float = 31.0,
    **kw: object,
) -> dict[str, object]:
    """One ``regions_stats.tsv`` row (1-based closed ``start``/``end``)."""
    row: dict[str, object] = {
        "region_id": region_id,
        "region_type": region_type,
        "subset": subset,
        "contig": contig,
        "start": start,
        "end": end,
        "n_genes": n_eligible,
        "n_eligible": n_eligible,
        "observed": observed,
        "expected": expected,
        "o_minus_e": observed - expected,
        "o_over_e": observed / expected if expected else None,
        "stat": stat,
        "effect": observed - expected,
        "p_value": q,
        "q_value": q,
        "method": "poisson_binomial",
        "null_model": "none",
        "n_perm": 0,
        "flags": "small_region" if q is None else None,
    }
    row.update(kw)
    return row


def island(
    island_id: str,
    contig: str,
    start: int,
    end: int,
    p: float | None,
    *,
    direction: str = "high",
    mean: float = 21.5,
    **kw: object,
) -> dict[str, object]:
    """One ``islands.tsv`` row (1-based closed; ``end`` may exceed the replicon length)."""
    row: dict[str, object] = {
        "island_id": island_id,
        "direction": direction,
        "contig": contig,
        "start": start,
        "end": end,
        "tu_first": "TU0001",
        "tu_last": "TU0002",
        "n_tus": 2,
        "n_genes": 3,
        "genes": "g1;g2;g3",
        "mean_abundance_log2": mean,
        "detected_frac": 0.5,
        "scan_s": 4.0,
        "block_len": 2,
        "stat": "scan_max",
        "effect": 1.0,
        "p_value": p,
        "q_value": None,
        "method": "scan_max_stat",
        "null_model": "circular_block_bootstrap",
        "n_perm": 999,
    }
    row.update(kw)
    return row


def flag(
    pair_id: str,
    gene_a: str,
    gene_b: str,
    q: float | None,
    *,
    flag: str = "internal_promoter_supported",
) -> dict[str, object]:
    """One ``operon_flags.tsv`` row."""
    return {
        "pair_id": pair_id,
        "gene_a": gene_a,
        "gene_b": gene_b,
        "uniop_class": "OP",
        "promoter_between_yes": True,
        "terminator_between_yes": False,
        "abundance_pattern": "discordant",
        "flag": flag,
        "q_value": q,
        "lfdr": q,
        "evidence_note": "test",
    }


def gene(
    tag: str,
    contig: str,
    start: int,
    end: int,
    *,
    strand: str = "+",
    log2: float | None = None,
    status: str = "detected",
    imputed: bool = False,
    tu_id: str | None = None,
    name: str | None = None,
    **kw: object,
) -> dict[str, object]:
    """One ``genes.tsv`` row (1-based closed); ``detected`` follows ``status``."""
    row = gene_row(
        locus_tag=tag,
        gene=name,
        contig=contig,
        start=start,
        end=end,
        strand=strand,
        feature_type="cds",
        detection_status=status,
        detected=status == "detected",
        abundance_log2=log2,
        abundance_raw=None if log2 is None or imputed else 2.0**log2,
        abundance_imputed=imputed,
        abundance_type="none" if log2 is None else "ibaq",
        tu_id=tu_id,
    )
    row.update(kw)
    return row


def tu(tu_id: str, contig: str, start: int, end: int, strand: str, genes: Sequence[str]) -> tuple:
    """One ``annotation/tus.tsv`` row (1-based closed, upstream region = 150 bp before)."""
    return (tu_id, contig, start, end, strand, ";".join(genes), max(1, start - 150), start - 1)


def region(region_id: str, region_type: str, contig: str, start: int, end: int) -> tuple:
    """One ``annotation/regions.tsv`` row (1-based closed)."""
    return (region_id, region_type, contig, start, end, "test", None, "", None)


def hit(peptide: str, contig: str, start: int, end: int, strand: str = "+") -> dict[str, object]:
    """One ``mapping/peptide_hits.tsv`` placement (1-based closed)."""
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": 1,
        "crosses_origin": False,
        "route": "exact",
        "n_diff": 0,
        "uniqueness": "unique",
    }


BedRow = tuple[str, int, int, str, float, str]


def _write_bed(path: Path, rows: Sequence[BedRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join("\t".join(str(c) for c in r) + "\n" for r in rows), encoding="utf-8")


def build_out(
    tmp: Path,
    *,
    replicons: Sequence[Mapping[str, Any]] = REPLICONS,
    genes: Sequence[Mapping[str, object]] = (),
    regions_stats: Sequence[Mapping[str, object]] | None = None,
    islands: Sequence[Mapping[str, object]] | None = None,
    flags: Sequence[Mapping[str, object]] | None = None,
    tus: Sequence[tuple] | None = None,
    regions: Sequence[tuple] | None = None,
    promoters: Sequence[BedRow] | None = None,
    terminators: Sequence[BedRow] | None = None,
    hits: Sequence[Mapping[str, object]] | None = None,
    novel_orfs: Sequence[Mapping[str, object]] | None = None,
) -> Path:
    """An ``OUT/`` with ``preflight.json`` and ``genes.tsv``, plus every table given."""
    out = make_out(tmp)
    pf = out / "preflight" / "preflight.json"
    pf.parent.mkdir(parents=True, exist_ok=True)
    pf.write_text(json.dumps(preflight(replicons), sort_keys=True), encoding="utf-8")
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if regions_stats is not None:
        write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats", regions_stats)
    if islands is not None:
        write_schema_table(out / "tables" / "islands.tsv", "islands", islands)
    if flags is not None:
        write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", flags)
    if novel_orfs is not None:
        write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", novel_orfs)
    if hits is not None:
        write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits", hits)
    if tus is not None:
        (out / "annotation").mkdir(parents=True, exist_ok=True)
        write_tsv(out / "annotation" / "tus.tsv", TUS_COLUMNS, tus)
    if regions is not None:
        (out / "annotation").mkdir(parents=True, exist_ok=True)
        write_tsv(out / "annotation" / "regions.tsv", REGIONS_COLUMNS, regions)
    if promoters is not None:
        _write_bed(out / "annotation" / "promoters.bed", promoters)
    if terminators is not None:
        _write_bed(out / "annotation" / "terminators.bed", terminators)
    return out
