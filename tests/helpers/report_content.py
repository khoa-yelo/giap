# ruff: noqa: RUF003  (the report's wording is Unicode: –, ℓ̂, −∞, σ70)
"""Small synthetic ``OUT/`` trees for the M9-16 report-content tests (summary cards, caveats,
sensitivity, provenance, wording; docs/specs/M9-16-report-content.md).

Builds on tests/helpers/report_fixtures.py (``run_out``, ``parse``, ``cfg``) and
tests/helpers/viz_out_tree.py (``write_schema_table``, ``gene_row``). Every table goes through its
Frictionless schema where one exists on main; ``determinants_model.tsv`` (M8-10, #178, not yet on
main) is written with the column list of its branch schema (``DET_MODEL_COLUMNS``).

Interface pinned by these helpers (M9-16 test author, beyond the spec text):
- the report reads ``OUT/run_manifest.json`` (MVP-03a, ``giap.manifest``) when present;
- KEGG database entries of the manifest are those whose ``name`` starts with ``kegg``
  (``kegg_ko_pathway``, ``kegg_pathway_names`` in ``giap/data/db_sources.yaml``); the GO release is
  the ``version`` of ``go_basic``;
- run flags come from step events (``logs/steps/<step>.events.jsonl``, ``data.flag(s)``), as
  ``RunView.run_flags`` reads them.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from giap.io.tsv import write_tsv
from helpers.report_fixtures import Page, cfg, parse
from helpers.viz_out_tree import gene_row, write_events, write_schema_table

#: tables/determinants_model.tsv columns (origin/m8/determinants-step schema)
DET_MODEL_COLUMNS = (
    "model", "feature_or_block", "metric", "value", "ci_low", "ci_high", "n", "n_perm", "flags",
)  # fmt: skip

#: §7.8.2 "Caveats (printed in the report)", verbatim (``giap.determinants.CAVEATS`` on #178)
DETERMINANTS_CAVEATS = (
    "In a synthetic GFP library, mRNA folding near the RBS explained over half the variance; "
    "native genes add regulation and turnover, so expect a lower R².",
    "All features are predictions, and the data come from one condition.",
    "CAI, tAI and GC3 are collinear, so interpret blocks rather than single features.",
)
#: docs/methods/islands.md "Validity domain (ADR-0021 Amendment 1)", first sentence (no Markdown)
ISLANDS_VALIDITY = "The FWER is claimed for a dependence regime, not for a replicon size."

KEGG_DATE = "2026-10-15"
GO_RELEASE = "2026-05-19"
SEED = 424242
COMMIT = "0123456789abcdef0123456789abcdef01234567"
TOOL_DIGEST = "sha256:" + "f" * 64
MANIFEST_INPUTS = [
    {"role": "genome", "file": "genome.fa", "path": "genome.fa", "sha256": "d" * 64, "size": 1200},
    {"role": "peptides", "file": "peptides.tsv", "path": "peptides.tsv", "sha256": "e" * 64,
     "size": 300},
]  # fmt: skip
DATABASES = [
    {"name": "go_basic", "version": GO_RELEASE, "checksum": "sha256:" + "9" * 64,
     "retrieved": "2026-10-15", "size_bytes": 1000},
    {"name": "kegg_ko_pathway", "version": KEGG_DATE, "checksum": "sha256:" + "8" * 64,
     "retrieved": KEGG_DATE, "size_bytes": 2000},
    {"name": "kegg_pathway_names", "version": KEGG_DATE, "checksum": "sha256:" + "7" * 64,
     "retrieved": KEGG_DATE, "size_bytes": 300},
]  # fmt: skip
TOOLS = [
    {"name": "diamond", "version": "2.1.9", "container_digest": TOOL_DIGEST,
     "conda_lock_hash": None, "processes": ["MAP_PROTEINS"]},
]  # fmt: skip


def words(text: str) -> str:
    return " ".join(text.split())


def page(out: Path, *overrides: str) -> Page:
    return parse(html(out, *overrides))


def html(out: Path, *overrides: str) -> str:
    from giap.io.run_view import RunView
    from giap.report.build import build_report

    return build_report(RunView.from_outdir(out), {}, {}, cfg(*overrides))


def section_text(p: Page, sid: str) -> str:
    return words(p.section(sid).text())


# ---------------------------------------------------------------- tables


def write_genes(out: Path, n_detected: int, n_not: int = 0, *, tier: str = "decoy_picked",
                n_2pep: int = 0) -> None:  # fmt: skip
    """``n_detected`` detected genes (the first ``n_2pep`` with ``detected_2pep``) and ``n_not``
    undetected ones, all in ``tier``."""
    rows = []
    for i in range(n_detected + n_not):
        det = i < n_detected
        rows.append(gene_row(locus_tag=f"g{i:05d}", detected=det, detected_2pep=i < n_2pep,
                             detection_status="detected" if det else "not_detected",
                             detection_tier=tier))  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", rows)


def write_tier_fdr(out: Path, fdr_at: Mapping[str, float | None]) -> None:
    """``detection_summary`` tier rows: ``fdr`` at each q threshold (``group``)."""
    rows = [{"section": "tier", "group_type": "q_threshold", "group": g, "metric": "fdr",
             "value": v} for g, v in fdr_at.items()]  # fmt: skip
    write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", rows)


def write_variants(out: Path, rows: Sequence[tuple[str, str, str]]) -> None:
    """``(variant_id, confidence, class_error_rate)`` ALT rows."""
    body = [{"variant_id": v, "contig": "chr1", "codon_start": 100 + 3 * i, "strand": "+",
             "gene": "g00001", "aa_pos": 10 + i, "ref_codon": "GCT", "ref_aa": "A",
             "alt_aa": "V", "alt_type": "substitution", "call": "ALT", "confidence": c,
             "n_ref": 0, "n_alt": 3, "nt_changes": 1, "alt_codon_vcf": "GTT",
             "class_error_rate": e, "class_error_n": 50, "source": "pileup"}
            for i, (v, c, e) in enumerate(rows)]  # fmt: skip
    write_schema_table(out / "tables" / "variants.tsv", "variants", body)


def write_novel(out: Path, rows: Sequence[tuple[str, bool, str]]) -> None:
    body = [{"orf_id": i, "class": "intergenic", "called": k, "class_fdr": f} for i, k, f in rows]
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", body)


def write_det_model(out: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
    path = out / "tables" / "determinants_model.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, list(DET_MODEL_COLUMNS), [[r.get(c) for c in DET_MODEL_COLUMNS] for r in rows])
    return path


def det_rows(
    model: str, r2: float | None, sd: float | None, p: float | None, flags: str = ""
) -> list[dict[str, Any]]:
    """The whole-model rows of one model (``feature_or_block = all``)."""
    base = {"model": model, "feature_or_block": "all", "n": 300}
    return [
        {**base, "metric": "cv_r2_mean", "value": r2, "n_perm": 0},
        {**base, "metric": "cv_r2_sd", "value": sd, "n_perm": 0},
        {**base, "metric": "null_p", "value": p, "n_perm": 999, "flags": flags or None},
    ]


def write_mapping(
    out: Path,
    *,
    placed_q: Sequence[float],
    unmapped_peptides: int,
    proteins: Sequence[str],
    unmapped_proteins: int,
    strata: Mapping[str, float],
) -> None:
    """Peptide hits (one per q in ``placed_q``), protein hits, unmapped rows of both kinds and the
    ``mapping_fdr`` placement strata (``length_class -> fdr_stratum``) plus one de novo row."""
    m = out / "mapping"
    hits = [{"peptide": f"PEPTIDE{i}K", "peptide_len": 9, "contig": "chr1", "start": 101 + 3 * i,
             "end": 128 + 3 * i, "strand": "+", "placement_q": q, "length_class": "9-12"}
            for i, q in enumerate(placed_q)]  # fmt: skip
    write_schema_table(m / "peptide_hits.tsv", "peptide_hits", hits)
    prot = [{"protein_id": p, "protein_len_aa": 100, "route": "diamond", "contig": "chr1",
             "start": 101, "end": 400, "strand": "+", "rank": 1} for p in proteins]  # fmt: skip
    write_schema_table(m / "protein_hits.tsv", "protein_hits", prot)
    un = [{"sequence_id": f"u{i}", "sequence": "X", "kind": "peptide", "label": "unmapped",
           "reason": "no_hit"} for i in range(unmapped_peptides)]  # fmt: skip
    un += [{"sequence_id": f"up{i}", "sequence": "M", "kind": "protein", "label": "unmapped",
            "reason": "below_identity"} for i in range(unmapped_proteins)]  # fmt: skip
    write_schema_table(m / "unmapped.tsv", "unmapped", un)
    fdr = [{"section": "placement", "length_class": lc, "n_diff": 0, "targets": 100,
            "decoys": 1, "fdr_stratum": f, "n_accepted": 90, "decoy_type": "reverse"}
           for lc, f in strata.items()]  # fmt: skip
    fdr.append({"section": "denovo_precision", "score_bin": 1, "precision": 0.5,
                "fdr_stratum": 0.4321})  # fmt: skip
    write_schema_table(m / "mapping_fdr.tsv", "mapping_fdr", fdr)


def write_islands(out: Path, p_values: Sequence[float]) -> None:
    from helpers.report_fixtures import island_rows

    rows = island_rows(len(p_values), p_values)
    write_schema_table(out / "tables" / "islands.tsv", "islands", rows)


def write_operon_flags(out: Path, q_values: Sequence[float]) -> None:
    rows = [{"pair_id": f"PAIR_{i}", "gene_a": "g1", "gene_b": "g2", "uniop_class": "OP",
             "abundance_pattern": "discordant", "flag": "internal_promoter_supported",
             "promoter_between_yes": True, "terminator_between_yes": False, "lfdr": 0.5,
             "q_value": q, "evidence_note": "x"} for i, q in enumerate(q_values)]  # fmt: skip
    write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", rows)


def write_pairs(out: Path) -> None:
    write_schema_table(out / "tables" / "pairs.tsv", "pairs", [])


def write_category_summary(out: Path, coverage: Mapping[str, float]) -> None:
    """One row per collection with its ``coverage``."""
    rows = [{"collection": c, "category": f"{c}_1", "name": f"{c} one", "n_genes": 10,
             "n_detected": 5, "expected": 5.0, "o_minus_e": 0.0, "coverage": v}
            for c, v in coverage.items()]  # fmt: skip
    write_schema_table(out / "tables" / "category_summary.tsv", "category_summary", rows)


def write_gsea(out: Path, flags: Sequence[str | None]) -> None:
    rows = [{"collection": "cog", "term": f"T{i}", "name": f"term {i}", "size_tus": 20,
             "es": 0.4, "nes": 1.5, "coverage": 0.6, "p_value": 0.01, "q_value": 0.04,
             "method": "fgsea_simple", "null_model": "tu_set_permutation", "n_perm": 1000,
             "flags": f} for i, f in enumerate(flags)]  # fmt: skip
    write_schema_table(out / "tables" / "enrichment_gsea.tsv", "enrichment_gsea", rows)


def write_diagnostics(out: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
    """``context/island_diagnostics.parquet`` in the ``giap.context.diagnostics`` shape."""
    from giap.context.diagnostics import DIAGNOSTICS_COLUMNS

    frame = pd.DataFrame(
        {
            "contig": pd.Series([r["contig"] for r in rows], dtype=object),
            "scanned": pd.Series([r["scanned"] for r in rows], dtype=bool),
            "block_rule": pd.Series([r["block_rule"] for r in rows], dtype=object),
            "l_hat": pd.Series([r["l_hat"] for r in rows], dtype="Int64"),
            "block_len": pd.Series([r["block_len"] for r in rows], dtype="Int64"),
            "morans_i_lag1": pd.Series([r["morans_i_lag1"] for r in rows], dtype="float64"),
            "flags": pd.Series([list(r.get("flags", [])) for r in rows], dtype=object),
        },
        columns=list(DIAGNOSTICS_COLUMNS),
    )
    path = out / "context" / "island_diagnostics.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return path


def write_sensitivity(
    out: Path, rows: Sequence[tuple[str, str, int, int, int, float | None]]
) -> None:
    body = [{"family": f, "variant": v, "n_primary": a, "n_variant": b, "n_shared": s,
             "jaccard": j} for f, v, a, b, s, j in rows]  # fmt: skip
    write_schema_table(out / "tables" / "sensitivity.tsv", "sensitivity", body)


# ---------------------------------------------------------------- events


def write_flags(out: Path, step: str, flags: Sequence[str]) -> None:
    write_events(out / "logs" / "steps" / f"{step}.events.jsonl",
                 [{"event": f"{step}.run_flags", "data": {"flags": list(flags)}}])  # fmt: skip


def write_missed_event(out: Path, pi0: float = 0.5) -> None:
    """``operons.missed_operons`` as ``giap.steps.operons`` logs it: three distance bins, the
    outer bounds ``-Infinity`` / ``Infinity`` (strict JSON cannot hold them; Python's can)."""
    bins = [
        ("-Infinity", "0.0", 3, 1, 2),
        ("0.0", "50.0", 7, 4, 5),
        ("50.0", "Infinity", 11, 6, 9),
    ]
    rec = ", ".join(
        f'{{"bin": {i}, "lower_bp": {lo}, "upper_bp": {hi}, "N_NO": {a}, "N_cv": {b}, "N_sw": {c}}}'
        for i, (lo, hi, a, b, c) in enumerate(bins)
    )
    line = (f'{{"event": "operons.missed_operons", "level": "INFO", "step": "operons", '
            f'"data": {{"pi0": {pi0}, "bin_counts": [{rec}]}}}}')  # fmt: skip
    write_events(out / "logs" / "steps" / "operons.events.jsonl", [line])


# ---------------------------------------------------------------- manifest


def manifest(**over: Any) -> dict[str, Any]:
    """A schema-valid ``run_manifest.json`` (``helpers.manifest_data.build``) with a commit, tools,
    databases (KEGG retrieval date, GO release), seed ``SEED`` and relative-path inputs."""
    from giap.config import GiapConfig
    from helpers.manifest_data import build

    m = build(cfg=GiapConfig(seed=SEED), git_commit=COMMIT, git_dirty=False,
              command_line=["giap", "run", "--genome", "genome.fa", "--peptides",
                            "peptides.tsv", "--seed", str(SEED)])  # fmt: skip
    m["inputs"] = [dict(r) for r in MANIFEST_INPUTS]
    m["tools"] = [dict(t) for t in TOOLS]
    m["databases"] = [dict(d) for d in DATABASES]
    m.update(over)
    return m


def write_manifest(out: Path, doc: Mapping[str, Any] | None = None) -> Path:
    from giap.manifest import canonical_json

    path = out / "run_manifest.json"
    path.write_text(canonical_json(doc if doc is not None else manifest()), encoding="utf-8")
    return path


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True)
