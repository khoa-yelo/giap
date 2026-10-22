"""Hand-built, valid M4 input tables (M4-06 contract C1-C14) for the contract tests.

Independent of ``giapsim`` (M4-07 is not merged yet) and of ``giap.genetable.contract``: the frames
below are written from the M4-06 spec table, so they check the validator rather than mirror it.

Layout: one circular contig ``chr1``; gene ``i`` (0-based) is ``GSIM_<5(i+1):05d>`` at
``[100 + 350 i, 100 + 350 i + 300)``; TUs are consecutive runs of ``tu_sizes`` (cycled), strands
alternate per TU; gene 3 is a sORF, gene 7 a pseudogene with ``aa = NA``. In-memory frames use the
internal representation (0-based half-open, Python lists for list columns); ``write_work`` writes
Parquet for C1-C10/C14-free files and M3's final 1-based TSVs for C11-C13 (spec "The contract").
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

GENE_LEN = 300
GENE_STEP = 350
OFFSET = 100

# Spec C6 (M2-20 LAYER_COLUMNS), in order.
LAYER_COLUMNS: tuple[str, ...] = (
    "locus_tag", "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids",
    "promoter_score", "promoter_score_tu", "n_promoter_windows_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream",
    "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10", "flags",
)  # fmt: skip
PROMOTER_COLUMNS: tuple[str, ...] = (
    "promoter_score", "promoter_score_tu", "n_promoter_windows_tu", "tx_rate_max_tu",
    "tx_rate_sum_tu", "tss_distance_bp",
)  # fmt: skip
_LAYER_VALUE_COLUMNS = tuple(c for c in LAYER_COLUMNS[7:] if c not in ("start_codon", "flags"))

PEPTIDE_HITS_COLUMNS: tuple[str, ...] = (
    "peptide", "peptide_len", "source", "best_q_input", "pep_input", "score_input", "contig",
    "start", "end", "strand", "frame", "crosses_origin", "route", "n_diff", "diff_types",
    "placement_q", "length_class", "locus_class", "uniqueness", "gene", "novel_candidate_id",
    "flags",
)  # fmt: skip
PROTEIN_HITS_COLUMNS: tuple[str, ...] = (
    "protein_id", "protein_len_aa", "route", "contig", "start", "end", "strand", "rank",
    "identity", "query_coverage", "reciprocal_overlap", "class", "gene", "flags",
)  # fmt: skip
#: C15 ``features`` (M5-04; the columns M2-04 ``parse_bakta`` writes), in order.
FEATURE_COLUMNS: tuple[str, ...] = (
    "feature_id", "type", "contig", "start", "end", "strand", "gene", "product", "amino_acid",
    "anticodon",
)  # fmt: skip
#: C16 ``peptide_runs`` (M6-02; the columns M1-06 ``aggregate_peptides`` writes), in order.
PEPTIDE_RUN_COLUMNS: tuple[str, ...] = ("peptide", "run", "intensity", "n_rows")
UNMAPPED_COLUMNS: tuple[str, ...] = (
    "sequence_id", "sequence", "kind", "label", "reason", "contaminant_match",
)  # fmt: skip
AGG_COLUMNS: tuple[str, ...] = (
    "peptide", "q_value", "pep", "score", "is_decoy", "is_contaminant", "protein_ids",
    "modified_peptides", "n_rows", "n_runs", "intensity", "source", "search_db",
)  # fmt: skip

# File names of the default layout (spec "The contract" column "File").
FILE_NAMES: dict[str, str] = {
    "replicons": "replicons.parquet",
    "genes": "genes.parquet",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
    "tus": "tus.parquet",
    "layers": "layers.parquet",
    "regions": "regions.parquet",
    "region_genes": "region_genes.parquet",
    "peptides": "peptides.norm.parquet",
    "abundance": "abundance.norm.parquet",
    "peptide_hits": "peptide_hits.tsv",
    "protein_hits": "protein_hits.tsv",
    "unmapped": "unmapped.tsv",
    "id_map": "id_map.tsv",
}
FRAME_KEYS: tuple[str, ...] = (
    "replicons", "genes", "xrefs", "tus", "layers", "regions", "region_genes", "peptides",
    "abundance", "peptide_hits", "protein_hits", "unmapped", "id_map",
)  # fmt: skip
TSV_KEYS = ("peptide_hits", "protein_hits", "unmapped", "id_map")


def tag(i: int) -> str:
    """Locus tag of 0-based gene ``i`` (M0 convention ``GSIM_<5·(i+1):05d>``)."""
    return f"GSIM_{5 * (i + 1):05d}"


def _modified_peptides_value(peptide: str) -> Any:
    """``modified_peptides`` encoding is not fixed by M4-06; follow ``TABLES`` when it exists."""
    try:
        from giap.genetable.contract import TABLES

        col = next(c for c in TABLES["peptides"] if c.name == "modified_peptides")
        return [peptide] if col.dtype == "list" else peptide
    except Exception:
        return peptide


def build_frames(n_genes: int = 20, tu_sizes: Sequence[int] = (1, 2, 3)) -> dict[str, Any]:
    """Valid C1-C15 frames (internal representation) plus ``annotation_source``."""
    if n_genes < 0:
        raise ValueError("n_genes must be >= 0")
    length = OFFSET + GENE_STEP * max(n_genes, 1) + 1000

    # TUs: consecutive runs, strands alternate per TU.
    tu_of: list[int] = []
    strand_of: list[str] = []
    t = 0
    while len(tu_of) < n_genes:
        size = int(tu_sizes[t % len(tu_sizes)])
        for _ in range(min(size, n_genes - len(tu_of))):
            tu_of.append(t)
            strand_of.append("+" if t % 2 == 0 else "-")
        t += 1
    n_tus = t

    starts = [OFFSET + GENE_STEP * i for i in range(n_genes)]
    ends = [s + GENE_LEN for s in starts]
    ftype = ["cds"] * n_genes
    if n_genes > 3:
        ftype[3] = "sorf"
    if n_genes > 7:
        ftype[7] = "pseudogene"
    aa: list[str | None] = ["M" + "A" * 98 for _ in range(n_genes)]
    if n_genes > 7:
        aa[7] = None
    genes = pd.DataFrame(
        {
            "locus_tag": [tag(i) for i in range(n_genes)],
            "gene": [f"gen{i}" if i % 3 == 0 else None for i in range(n_genes)],
            "product": [None if i % 4 == 1 else f"protein {i}" for i in range(n_genes)],
            "feature_type": ftype,
            "contig": ["chr1"] * n_genes,
            "start": pd.Series(starts, dtype="int64"),
            "end": pd.Series(ends, dtype="int64"),
            "strand": strand_of,
            "partial": pd.Series([False] * n_genes, dtype=bool),
            "crosses_origin": pd.Series([False] * n_genes, dtype=bool),
            "hypothetical": pd.Series([i % 5 == 2 for i in range(n_genes)], dtype=bool),
            "aa": aa,
            "nt": ["ATG" + "GCT" * 98 + "TAA" for _ in range(n_genes)],
            "start_codon": ["ATG"] * n_genes,
        }
    )
    genes = genes.astype({"start": "int64", "end": "int64"})

    replicons = pd.DataFrame(
        {
            "contig": ["chr1"],
            "length": pd.Series([length], dtype="int64"),
            "topology": ["circular"],
            "kind": ["chromosome"],
        }
    )

    xref_rows = []
    for i in range(0, n_genes, 2):
        xref_rows.append((tag(i), "COG", f"COG{i:04d}"))
        xref_rows.append((tag(i), "GO", f"GO:{i:07d}"))
    xrefs = pd.DataFrame(xref_rows, columns=["locus_tag", "db", "id"])

    # TU table and per-gene transcription order.
    tu_rows = []
    tu_position: dict[int, int] = {}
    tu_size: dict[int, int] = {}
    tu_order: dict[int, list[int]] = {}
    for u in range(n_tus):
        members = [i for i in range(n_genes) if tu_of[i] == u]
        order = members if strand_of[members[0]] == "+" else members[::-1]
        tu_order[u] = order
        for pos, g in enumerate(order, start=1):
            tu_position[g] = pos
            tu_size[g] = len(order)
        tu_rows.append(
            (
                f"TU_{u + 1}",
                "chr1",
                starts[members[0]],
                ends[members[-1]],
                strand_of[members[0]],
                [tag(g) for g in order],
            )
        )
    tus = pd.DataFrame(tu_rows, columns=["tu_id", "contig", "start", "end", "strand", "genes"])
    tus = tus.astype({"start": "int64", "end": "int64"})

    # Regions: prophage on genes 10..19 (when present), AMR scattered (genes 2 and last).
    region_rows: list[tuple[Any, ...]] = []
    rg_rows: list[tuple[Any, ...]] = []
    region_ids: list[list[str]] = [[] for _ in range(n_genes)]
    if n_genes >= 20:
        region_rows.append(("prophage_1", "prophage", "chr1", starts[10], ends[19], "genomad", 0.9))
        subsets = ("integrase", "repressor", "structural", "lysis", "cargo")
        for j, g in enumerate(range(10, 20)):
            rg_rows.append(("prophage_1", tag(g), 1.0, subsets[j % 5], False))
            region_ids[g].append("prophage_1")
    if n_genes >= 3:
        region_rows.append(("amr_1", "amr", "chr1", None, None, "amrfinderplus", None))
        for g in (2, n_genes - 1):
            rg_rows.append(("amr_1", tag(g), 1.0, None, False))
            region_ids[g].append("amr_1")
    regions = pd.DataFrame(
        region_rows,
        columns=["region_id", "region_type", "contig", "start", "end", "source_tool", "score"],
    )
    regions = regions.astype({"start": "Int64", "end": "Int64", "score": "float64"})
    region_genes = pd.DataFrame(
        rg_rows, columns=["region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic"]
    )
    region_genes = region_genes.astype({"overlap_frac": "float64", "core_biosynthetic": bool})

    # Layers (C6): structural columns filled, every regulatory/codon layer all-NA.
    prev_p: list[float] = []
    next_p: list[float] = []
    for g in range(n_genes):
        order = tu_order[tu_of[g]]
        pos = tu_position[g] - 1
        prev_p.append(0.9 if pos > 0 else math.nan)
        next_p.append(0.9 if pos < len(order) - 1 else math.nan)
    layers: dict[str, Any] = {
        "locus_tag": [tag(i) for i in range(n_genes)],
        "tu_id": [f"TU_{tu_of[i] + 1}" for i in range(n_genes)],
        "tu_position": pd.Series([tu_position[i] for i in range(n_genes)], dtype="int64"),
        "tu_size": pd.Series([tu_size[i] for i in range(n_genes)], dtype="int64"),
        "uniop_p_prev": pd.Series(prev_p, dtype="float64"),
        "uniop_p_next": pd.Series(next_p, dtype="float64"),
        "region_ids": [sorted(r) for r in region_ids],
    }
    for c in _LAYER_VALUE_COLUMNS:
        layers[c] = pd.Series([None] * n_genes, dtype=object)
    layers["start_codon"] = ["ATG"] * n_genes
    layers["flags"] = [["pseudogene"] if ftype[i] == "pseudogene" else [] for i in range(n_genes)]
    layers_df = pd.DataFrame(layers)[list(LAYER_COLUMNS)]

    # Peptides (C9): three targets and one decoy.
    pep_rows = [
        ("AAAAAAAK", 0.001, 0.01, 50.0, False, False, tag(min(4, max(n_genes - 1, 0)))),
        ("CCCCCCCR", 0.004, 0.02, 40.0, False, False, None),
        ("DDDDDDDK", 0.002, 0.01, 45.0, False, False, tag(7) if n_genes > 7 else None),
        ("KAAAAAAA", 0.30, 0.90, 5.0, True, False, "rev_" + tag(0)),
    ]
    peptides = pd.DataFrame(
        {
            "peptide": [r[0] for r in pep_rows],
            "q_value": pd.Series([r[1] for r in pep_rows], dtype="float64"),
            "pep": pd.Series([r[2] for r in pep_rows], dtype="float64"),
            "score": pd.Series([r[3] for r in pep_rows], dtype="float64"),
            "is_decoy": pd.Series([r[4] for r in pep_rows], dtype=bool),
            "is_contaminant": pd.Series([r[5] for r in pep_rows], dtype=bool),
            # Lead decision 1: native list of strings in Parquet (as M1-12 writes it).
            "protein_ids": [[r[6]] if r[6] is not None else [] for r in pep_rows],
            "modified_peptides": [_modified_peptides_value(r[0]) for r in pep_rows],
            "n_rows": pd.Series([2, 1, 1, 1], dtype="int64"),
            "n_runs": pd.Series([1, 1, 1, 1], dtype="int64"),
            "intensity": pd.Series([1.0e6, math.nan, 2.0e5, 3.0e3], dtype="float64"),
            "source": ["database"] * 4,
            "search_db": ["same"] * 4,
        }
    )[list(AGG_COLUMNS)]

    # C16 (M6-02): per-run detail of C9, one run per peptide (n_runs = 1, n_rows as in C9).
    peptide_runs = pd.DataFrame(
        {
            "peptide": [r[0] for r in pep_rows],
            "run": ["run1"] * len(pep_rows),
            "intensity": pd.Series([1.0e6, math.nan, 2.0e5, 3.0e3], dtype="float64"),
            "n_rows": pd.Series([2, 1, 1, 1], dtype="int64"),
        },
        columns=list(PEPTIDE_RUN_COLUMNS),
    )

    abundance = pd.DataFrame(
        {
            "protein_id": [tag(i) for i in range(min(n_genes, 5))],
            "value": pd.Series([20.0 + i for i in range(min(n_genes, 5))], dtype="float64"),
            # Merged M1-12 columns; the scale lives in preflight.json (lead decision 1, review).
            "raw_id": [tag(i) for i in range(min(n_genes, 5))],
        }
    )

    # Peptide hits (C11), internal 0-based half-open coordinates.
    hit_rows: list[dict[str, Any]] = []
    if n_genes > 4:
        g = 4
        hit_rows.append(
            _hit("AAAAAAAK", starts[g] + 30, strand_of[g], "cds_inframe", "unique", tag(g))
        )
    if n_genes > 7:
        g = 7
        hit_rows.append(
            _hit("DDDDDDDK", starts[g] + 60, strand_of[g], "pseudogene", "unique", tag(g))
        )
    hit_rows.append(
        _hit("CCCCCCCR", ends[-1] + 10 if n_genes else 10, "+", "intergenic", "unique", None)
    )
    peptide_hits = _hits_frame(hit_rows)

    prot_rows: list[dict[str, Any]] = []
    if n_genes > 5:
        g = 5
        prot_rows.append(
            {
                "protein_id": tag(g),
                "protein_len_aa": 99,
                "route": "miniprot",
                "contig": "chr1",
                "start": starts[g],
                "end": ends[g],
                "strand": strand_of[g],
                "rank": 1,
                "identity": 1.0,
                "query_coverage": 1.0,
                "reciprocal_overlap": 1.0,
                "class": "identical",
                "gene": tag(g),
                "flags": [],
            }
        )
    protein_hits = _typed_frame(prot_rows, PROTEIN_HITS_COLUMNS, PROTEIN_HITS_TYPES)

    unmapped = pd.DataFrame(
        [("EEEEEEEK", "EEEEEEEK", "peptide", "unmapped", "no_hit", None)],
        columns=list(UNMAPPED_COLUMNS),
    )
    id_map = pd.DataFrame(
        {
            "source_id": [f"orig_{i}" for i in range(n_genes)],
            "locus_tag": [tag(i) for i in range(n_genes)],
        }
    )
    # C15 (M5-04): an oriC in the free stretch after the last gene and one tRNA there too.
    free = OFFSET + GENE_STEP * max(n_genes, 1)
    features = pd.DataFrame(
        {
            "feature_id": ["GSIM_oric1", "GSIM_trna1"],
            "type": ["oric", "trna"],
            "contig": ["chr1", "chr1"],
            "start": pd.Series([free + 100, free + 500], dtype="int64"),
            "end": pd.Series([free + 332, free + 576], dtype="int64"),
            "strand": [".", "+"],
            "gene": [None, None],
            "product": [None, "tRNA-Leu"],
            "amino_acid": [None, "Leu"],
            "anticodon": [None, "CAG"],
        },
        columns=list(FEATURE_COLUMNS),
    )
    return {
        "replicons": replicons,
        "genes": genes,
        "xrefs": xrefs,
        "features": features,
        "annotation_source": "bakta_full",
        "abundance_scale": "log2",  # preflight.json routing.abundance_scale (M1-12)
        "tus": tus,
        "layers": layers_df,
        "regions": regions,
        "region_genes": region_genes,
        "peptides": peptides,
        "peptide_runs": peptide_runs,
        "abundance": abundance,
        "peptide_hits": peptide_hits,
        "protein_hits": protein_hits,
        "unmapped": unmapped,
        "id_map": id_map,
    }


PEPTIDE_HITS_TYPES = {
    "peptide_len": "int64", "best_q_input": "float64", "pep_input": "float64",
    "score_input": "float64", "start": "int64", "end": "int64", "frame": "int64",
    "crosses_origin": bool, "n_diff": "int64", "placement_q": "float64",
}  # fmt: skip
PROTEIN_HITS_TYPES = {
    "protein_len_aa": "int64", "start": "int64", "end": "int64", "rank": "int64",
    "identity": "float64", "query_coverage": "float64", "reciprocal_overlap": "float64",
}  # fmt: skip


def _hit(
    peptide: str, start: int, strand: str, locus_class: str, uniqueness: str, gene: str | None
) -> dict[str, Any]:
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "best_q_input": 0.001,
        "pep_input": 0.01,
        "score_input": 50.0,
        "contig": "chr1",
        "start": start,
        "end": start + 3 * len(peptide),
        "strand": strand,
        "frame": start % 3,
        "crosses_origin": False,
        "route": "exact",
        "n_diff": 0,
        "diff_types": [],
        "placement_q": 0.001,
        "length_class": "short",
        "locus_class": locus_class,
        "uniqueness": uniqueness,
        "gene": gene,
        "novel_candidate_id": None,
        "flags": [],
    }


def _typed_frame(
    rows: list[dict[str, Any]], columns: Sequence[str], types: dict[str, Any]
) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(columns))
    if not rows:
        df = pd.DataFrame({c: pd.Series([], dtype=object) for c in columns})
    return df.astype(types)


def _hits_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return _typed_frame(rows, PEPTIDE_HITS_COLUMNS, PEPTIDE_HITS_TYPES)


def empty_like(df: pd.DataFrame) -> pd.DataFrame:
    """Header-only copy (same columns and dtypes)."""
    return df.iloc[0:0].copy()


def make_inputs(frames: dict[str, Any]) -> Any:
    """``giap.io.work_tables.M4Inputs`` from ``build_frames`` output."""
    import dataclasses

    from giap.io.work_tables import M4Inputs

    # Lead decision 1 after review: the abundance scale (from preflight.json) is an M4Inputs field.
    names = {f.name for f in dataclasses.fields(M4Inputs)}
    extra = {"abundance_scale": frames["abundance_scale"]} if "abundance_scale" in names else {}
    return M4Inputs(
        replicons=frames["replicons"],
        genes=frames["genes"],
        xrefs=frames["xrefs"],
        annotation_source=frames["annotation_source"],
        tus=frames["tus"],
        layers=frames["layers"],
        regions=frames["regions"],
        region_genes=frames["region_genes"],
        peptides=frames["peptides"],
        abundance=frames["abundance"],
        **extra,
        peptide_hits=frames["peptide_hits"],
        protein_hits=frames["protein_hits"],
        unmapped=frames["unmapped"],
        id_map=frames["id_map"],
    )


def _tsv_rows(df: pd.DataFrame, one_based: bool) -> list[list[Any]]:
    out = []
    for rec in df.to_dict("records"):
        row = []
        for c in df.columns:
            v = rec[c]
            if one_based and c == "start" and v is not None:
                v = int(v) + 1  # 0-based half-open -> 1-based closed (end unchanged)
            if isinstance(v, float) and math.isnan(v):
                v = None
            row.append(v)
        out.append(row)
    return out


def write_work(frames: dict[str, Any], paths: dict[str, Path], skip: Sequence[str] = ()) -> None:
    """Write every contract file to ``paths`` (as ``default_paths`` returns), except ``skip``."""
    from giap.io.tsv import write_tsv

    for key, path in paths.items():
        if key == "preflight" and key not in skip:
            # M1-12 report subset the loader reads (lead decision 1 after review).
            path.parent.mkdir(parents=True, exist_ok=True)
            report = {"routing": {"abundance_scale": frames["abundance_scale"]}}
            path.write_text(json.dumps(report, sort_keys=True) + "\n", encoding="utf-8")
            continue
        if key in skip or key not in FILE_NAMES:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if key == "annotation_meta":
            meta = {"annotation_source": frames["annotation_source"], "flags": []}
            path.write_text(json.dumps(meta, sort_keys=True) + "\n", encoding="utf-8")
        elif key in TSV_KEYS:
            df = frames[key]
            write_tsv(path, list(df.columns), _tsv_rows(df, one_based=key != "id_map"))
        else:
            frames[key].to_parquet(path, engine="pyarrow", index=False)
