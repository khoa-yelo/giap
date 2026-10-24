"""K-12 plot survey, batch 1 item 8: V01's legend lists only the categories it draws.

The legend listed every detection status, region type, island direction, tick track and marker
whatever the run held. Required: an entry only for a category drawn in this figure (one row of
that track / category in the data TSV). The ring palette itself is batch 2.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import pytest
from helpers.viz_out_tree import gene_row, make_out, write_schema_table

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.io.tsv import write_tsv
from giap.viz import base, style
from giap.viz.base import Instance, Skip

REGION_COLUMNS = ["region_id", "region_type", "contig", "start", "end", "source_tool", "score",
                  "genes", "notes"]  # fmt: skip
EVERY = {
    "detected", "not detected", "partial excluded", "single peptide unestimated",
    "prophage arc", "plasmid arc", "bgc arc", "crispr arc", "amr arc", "island high",
    "island low", "TU start", "promoter", "terminator", "ALT site (high)",
    "ALT site (candidate)", "novel ORF",
}  # fmt: skip


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _gene(tag: str, start: int, status: str) -> dict[str, object]:
    return gene_row(locus_tag=tag, contig="chr1", start=start, end=start + 899, strand="+",
                    feature_type="cds", detection_status=status, detected=status == "detected",
                    partial=False, abundance_log2=None, abundance_imputed=False,
                    abundance_type="none")  # fmt: skip


def _variant(vid: str, codon: int) -> dict[str, object]:
    return {"variant_id": vid, "contig": "chr1", "codon_start": codon, "strand": "+",
            "gene": "g1", "aa_pos": 10, "ref_codon": "GCT", "ref_aa": "A", "alt_aa": "V",
            "alt_type": "substitution", "call": "ALT", "confidence": "high", "n_ref": 0,
            "n_alt": 3, "nt_changes": 1}  # fmt: skip


def _render(tmp: Path, genes: Sequence[dict[str, object]], regions: list | None = None,
            variants: list | None = None) -> Instance:  # fmt: skip
    import giap.viz.circos  # noqa: F401  (registers V01)

    out = make_out(tmp)
    pf = out / "preflight" / "preflight.json"
    pf.parent.mkdir(parents=True, exist_ok=True)
    reps = [{"id": "chr1", "kind": "chromosome", "length": 100_000, "topology": "circular"}]
    pf.write_text(json.dumps({"replicons": reps,
                              "contig_map": [{"id": "chr1", "original_id": "chr1"}]}))  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if regions is not None:
        p = out / "annotation" / "regions.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(p, REGION_COLUMNS, regions)
    if variants is not None:
        write_schema_table(out / "tables" / "variants.tsv", "variants", variants)
    result = base.registered()["V01"].fn(RunView.from_outdir(out), resolve_config())
    assert not isinstance(result, Skip), result
    return result[0]


def _legend(inst: Instance) -> set[str]:
    legs = [*inst.figure.legends, *(ax.get_legend() for ax in inst.figure.axes)]
    return {t.get_text() for leg in legs if leg is not None for t in leg.get_texts()}


def test_genes_only_run_lists_only_its_status(tmp_path: Path) -> None:
    inst = _render(tmp_path, [_gene("g1", 1_001, "detected"), _gene("g2", 5_001, "detected")])
    labels = _legend(inst)
    assert "detected" in labels
    assert labels & EVERY == {"detected"}, sorted(labels & EVERY)


def test_drawn_categories_are_listed_and_absent_ones_are_not(tmp_path: Path) -> None:
    genes = [_gene("g1", 1_001, "detected"), _gene("g2", 5_001, "not_detected")]
    regions = [["prophage_1", "prophage", "chr1", 10_001, 20_000, "geNomad", None, "g1", ""]]
    inst = _render(tmp_path, genes, regions=regions, variants=[_variant("v1", 1_030)])
    labels = _legend(inst)
    drawn = {"detected", "not detected", "prophage arc", "ALT site (high)"}
    assert drawn <= labels, sorted(labels)
    assert labels & EVERY == drawn, sorted((labels & EVERY) - drawn)
