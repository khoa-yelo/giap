"""Builders for the M9-22 ``gippa compare`` tests (docs/specs/M9-22-compare-tables.md).

Two kinds of input:

- ``genes_frame(rows)``: a ``genes.tsv``-shaped frame as ``RunView.table("genes_tsv")`` returns it
  (every schema column, schema dtypes: ``Int64`` integers, float numbers, nullable ``boolean``,
  object strings; 0-based ``start``), for the pure functions of ``gippa.compare``;
- ``write_run(out, ...)``: a small ``OUT/`` tree (``tables/genes.tsv`` through the schema, the
  optional tables of ``compare_stability``, ``preflight/preflight.json`` and / or
  ``run_manifest.json`` carrying the genome SHA-256) for ``RunView`` and the CLI.

Interface assumed (test-author reading of the spec): the genome digest is the ``sha256`` of the
``inputs`` row with ``role = "genome"`` in ``run_manifest.json`` (when the file exists) or in
``preflight/preflight.json``. The manifest written here holds only ``schema_version`` and
``inputs``; the compare reader is expected to read the genome row only.

Nothing here touches ``$GIPPA_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from helpers.viz_out_tree import schema_fields, write_schema_table

SHA_A = "a" * 64
SHA_B = "b" * 64

#: compare table columns (spec "Tables")
DETECTION_COLUMNS = (
    "section", "locus_tag", "run", "detected", "detection_status", "intersection", "n_genes",
)  # fmt: skip
ABUNDANCE_COLUMNS = (
    "run_a", "run_b", "n_shared", "n_units", "spearman_rho", "ci_low", "ci_high", "ba_bias_log2",
    "ba_loa_low", "ba_loa_high", "method", "n_boot", "flags",
)  # fmt: skip
STABILITY_COLUMNS = ("family", "run_a", "run_b", "n_a", "n_b", "n_shared", "jaccard", "flags")
FAMILIES = ("islands", "operon_flags", "enriched_categories", "novel_orfs", "alt_sites")


def gene(
    locus_tag: str,
    *,
    start: int,
    detected: bool = True,
    log2: float | None = None,
    imputed: bool = False,
    tu: str | None = None,
    abundance_type: str = "maxlfq",
    abundance_scale: str = "relative",
    contig: str = "chr1",
    **extra: object,
) -> dict[str, object]:
    """One gene row (0-based ``start``; ``end = start + 300``; TU defaults to ``TU_<locus>``)."""
    return {
        "locus_tag": locus_tag,
        "contig": contig,
        "start": start,
        "end": start + 300,
        "strand": "+",
        "feature_type": "cds",
        "partial": False,
        "hypothetical": False,
        "crosses_origin": False,
        "tu_id": tu if tu is not None else f"TU_{locus_tag}",
        "detected": detected,
        "detection_status": "detected" if detected else "not_detected",
        "abundance_log2": log2,
        "abundance_raw": None if log2 is None else float(2.0**log2),
        "abundance_imputed": imputed,
        "abundance_type": abundance_type,
        "abundance_scale": abundance_scale,
        **extra,
    }


def genes_frame(rows: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    """Every ``genes`` schema column with the ``RunView`` dtypes; absent cells NA."""
    cols: dict[str, pd.Series] = {}
    for f in schema_fields("genes"):
        name, kind = f["name"], f.get("type", "string")
        values = [r.get(name) for r in rows]
        if kind == "integer":
            cols[name] = pd.Series(values, dtype="Int64")
        elif kind == "number":
            cols[name] = pd.Series(
                [np.nan if v is None else float(v) for v in values],  # type: ignore[arg-type]
                dtype=float,
            )
        elif kind == "boolean":
            cols[name] = pd.Series(values, dtype="boolean")
        else:
            cols[name] = pd.Series(values, dtype=object)
    return pd.DataFrame(cols)


def tus_of(frame: pd.DataFrame) -> pd.Series:
    """``tu_id`` indexed by ``locus_tag`` (the ``tus`` argument of ``abundance_pairs``)."""
    return pd.Series(frame["tu_id"].to_numpy(), index=frame["locus_tag"].to_numpy(), name="tu_id")


def _one_based(rows: Sequence[Mapping[str, object]], cols: Sequence[str]) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        row = dict(r)
        for c in cols:
            if row.get(c) is not None:
                row[c] = int(row[c]) + 1  # type: ignore[call-overload]
        out.append(row)
    return out


def write_genes(out: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    """``OUT/tables/genes.tsv`` (1-based ``start`` on disk)."""
    return write_schema_table(out / "tables" / "genes.tsv", "genes", _one_based(rows, ["start"]))


def write_preflight(out: Path, sha: str | None) -> Path:
    inputs = [] if sha is None else [{"role": "genome", "file": "genome.fna", "sha256": sha,
                                      "size": 1000}]  # fmt: skip
    doc = {"schema_version": "1.0", "inputs": inputs, "checks": [], "replicons": []}
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return path


def write_manifest(out: Path, sha: str) -> Path:
    doc = {
        "schema_version": "1.0",
        "inputs": [{"role": "genome", "file": "genome.fna", "path": "/x/genome.fna",
                    "sha256": sha, "size": 1000}],
    }  # fmt: skip
    path = out / "run_manifest.json"
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return path


def write_run(
    out: Path,
    genes: Sequence[Mapping[str, object]],
    *,
    preflight_sha: str | None = SHA_A,
    manifest_sha: str | None = None,
    tables: Mapping[str, Sequence[Mapping[str, object]]] | None = None,
) -> Path:
    """An ``OUT/`` tree; ``tables`` maps a schema stem (``islands``, ``operon_flags``,
    ``enrichment_detection``, ``enrichment_gsea``, ``novel_orfs``, ``variants``) to its rows."""
    (out / "tables").mkdir(parents=True, exist_ok=True)
    write_genes(out, genes)
    if preflight_sha is not None or manifest_sha is None:
        write_preflight(out, preflight_sha)
    if manifest_sha is not None:
        write_manifest(out, manifest_sha)
    for stem, rows in (tables or {}).items():
        write_table(out, stem, rows)
    return out


def write_table(out: Path, stem: str, rows: Sequence[Mapping[str, object]]) -> Path:
    starts = {"islands": ["start"], "novel_orfs": ["start_longest", "start_best_ostir"]}
    return write_schema_table(
        out / "tables" / f"{stem}.tsv", stem, _one_based(rows, starts.get(stem, []))
    )


# ------------------------------------------------------------- stability table rows


def island(genes: Sequence[str], p: float, island_id: str = "ISL_1") -> dict[str, object]:
    return {"island_id": island_id, "direction": "high", "contig": "chr1", "start": 0,
            "end": 1000, "n_genes": len(genes), "genes": ";".join(genes), "p_value": p,
            "q_value": p, "stat": "scan_max", "method": "scan_max_stat",
            "null_model": "circular_block_bootstrap", "n_perm": 999}  # fmt: skip


def operon_flag(a: str, b: str, q: float) -> dict[str, object]:
    return {"pair_id": f"{a}_{b}", "gene_a": a, "gene_b": b, "uniop_class": "OP",
            "abundance_pattern": "discordant", "flag": "unexplained_discordance",
            "q_value": q}  # fmt: skip


def ora(collection: str, term: str, direction: str, q: float) -> dict[str, object]:
    return {"collection": collection, "term": term, "name": term, "direction": direction,
            "p_value": q, "q_value": q, "method": "cpb_tu_ora",
            "null_model": "conditional_poisson_binomial_tu", "n_perm": 0}  # fmt: skip


def gsea(collection: str, term: str, nes: float, q: float) -> dict[str, object]:
    return {"collection": collection, "term": term, "name": term, "nes": nes, "es": nes / 2,
            "p_value": q, "q_value": q, "method": "fgsea_simple",
            "null_model": "tu_set_permutation", "n_perm": 1000}  # fmt: skip


def novel(orf_id: str, called: bool) -> dict[str, object]:
    return {"orf_id": orf_id, "class": "intergenic", "contig": "chr1", "strand": "+",
            "frame": 0, "start_longest": 10, "end": 100, "called": called}  # fmt: skip


def variant(codon_start: int, alt_aa: str, call: str = "ALT") -> dict[str, object]:
    return {"variant_id": f"V{codon_start}{alt_aa}", "contig": "chr1",
            "codon_start": codon_start, "strand": "+", "gene": "g1", "alt_aa": alt_aa,
            "ref_aa": "A", "alt_type": "substitution", "call": call,
            "confidence": "high"}  # fmt: skip


def flag_set(cell: object) -> set[str]:
    """The ``;``-joined ``flags`` cell as a set (NA / empty -> empty set)."""
    if cell is None or (isinstance(cell, float) and np.isnan(cell)) or cell is pd.NA:
        return set()
    text = str(cell)
    return {t for t in text.split(";") if t and t != "NA"}
