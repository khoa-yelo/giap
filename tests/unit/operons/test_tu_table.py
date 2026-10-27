"""``operons.tsv``: one row per C5 TU (issue M6-07 step 5; GIPPA_SPEC.md §7.4.7, §9.3;
docs/specs/M6-07-flags-and-operons-table.md example 11 and the degenerate cases).

Readings this file pins (the spec leaves them open):

- ``operons_table(tus, pairs, genes, flagged)``: ``pairs`` is the merged pair frame, ``flagged``
  the ``flag_pairs`` output (``n_flags`` needs the flags and the spec's three-argument form has no
  config to recompute them).
- ``detected_frac`` = detected eligible genes / eligible genes, eligible = not ``partial``
  (§7.3 "Eligible genes ... are the non-partial protein-coding rows").
- ``max_abs_delta_log2`` uses the TU's internal adjacent pairs only (consecutive genes of the
  TU's ``genes`` list, both detected): the NO pair joining two TUs never counts.
- Contig order ("C1 order", §6.5) is the order in which contigs first appear in ``tus``
  (``tus.parquet`` is written in C1 order); within a contig rows sort by ``start``.
- Coordinates stay 0-based half-open and unwrapped (``end > length`` across the origin); only
  the M6-08 writer converts them.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
import pytest
from helpers.flags_data import cfg, chain_pairs, gene_table, tu, tus

from gippa.operons import flags as fl
from gippa.operons import tu_table as tt

ROOT = Path(__file__).resolve().parents[3]
OPERON_COLUMNS = (
    "tu_id", "contig", "start", "end", "strand", "genes", "n_genes", "detected_frac",
    "median_abundance_log2", "max_abs_delta_log2", "n_flags",
)  # fmt: skip


def _table(tu_rows, tu_genes, abundance, partial=(), **chain):
    c = cfg()
    t = tus(*tu_rows)
    tags = [g for gs in tu_genes for g in gs]
    pairs = chain_pairs(tu_genes, abundance, **chain)
    genes = gene_table(abundance, tags, partial=partial)
    flagged = fl.flag_pairs(pairs, c)
    return tt.operons_table(t, pairs, genes, flagged), pairs, flagged


def _by_id(out: pd.DataFrame) -> pd.DataFrame:
    return out.set_index("tu_id")


def test_operon_columns_match_spec_9_3() -> None:
    line = next(
        ln
        for ln in (ROOT / "GIPPA_SPEC.md").read_text(encoding="utf-8").splitlines()
        if ln.startswith("**`tables/operons.tsv`:**")
    )
    cols = tuple(re.findall(r"`([^`]+)`", line))[1:]
    assert cols == OPERON_COLUMNS
    assert tuple(tt.OPERON_COLUMNS) == OPERON_COLUMNS


def test_example_11() -> None:
    out, _, _ = _table(
        [tu("T1", ["A", "B", "C"], 100, 1000)], [["A", "B", "C"]], {"A": 17.0, "B": 18.4}
    )
    assert tuple(out.columns) == OPERON_COLUMNS
    r = out.iloc[0]
    assert r["tu_id"] == "T1"
    assert list(r["genes"]) == ["A", "B", "C"]
    assert r["n_genes"] == 3
    assert r["detected_frac"] == pytest.approx(2 / 3)
    assert round(r["detected_frac"], 3) == 0.667
    assert r["median_abundance_log2"] == pytest.approx(17.7)
    assert r["max_abs_delta_log2"] == pytest.approx(1.4)
    assert r["n_flags"] == 0


def test_max_abs_delta_is_absolute_and_internal_only() -> None:
    tu_genes = [["A", "B", "C"], ["D", "E"]]
    ab = {"A": 10.0, "B": 8.0, "C": 9.0, "D": 20.0, "E": 20.5}
    out, _, _ = _table(
        [tu("T1", tu_genes[0], 0, 900), tu("T2", tu_genes[1], 1000, 1600)], tu_genes, ab
    )
    by = _by_id(out)
    assert by.loc["T1", "max_abs_delta_log2"] == pytest.approx(2.0)
    # C -> D (|Delta| = 11) joins the two TUs and counts in neither
    assert by.loc["T2", "max_abs_delta_log2"] == pytest.approx(0.5)


def test_single_gene_tu() -> None:
    out, _, _ = _table(
        [tu("T1", ["A"], 0, 300), tu("T2", ["B", "C"], 400, 1000)],
        [["A"], ["B", "C"]],
        {"A": 12.0, "B": 10.0, "C": 11.0},
        q={("A", "B"): 0.04},
    )
    by = _by_id(out)
    assert by.loc["T1", "n_genes"] == 1
    assert math.isnan(by.loc["T1", "max_abs_delta_log2"])
    assert by.loc["T1", "n_flags"] == 0
    assert by.loc["T1", "detected_frac"] == pytest.approx(1.0)
    assert by.loc["T1", "median_abundance_log2"] == pytest.approx(12.0)


def test_tu_with_no_detected_gene() -> None:
    out, _, _ = _table([tu("T1", ["A", "B"], 0, 800)], [["A", "B"]], {})
    r = out.iloc[0]
    assert r["detected_frac"] == 0
    assert math.isnan(r["median_abundance_log2"])
    assert math.isnan(r["max_abs_delta_log2"])
    assert r["n_flags"] == 0


def test_partial_genes_are_not_eligible() -> None:
    out, _, _ = _table(
        [tu("T1", ["A", "B", "C"], 0, 900)], [["A", "B", "C"]], {"B": 10.0}, partial=("A",)
    )
    assert out.iloc[0]["detected_frac"] == pytest.approx(0.5)


def test_n_flags_counts_pairs_with_both_genes_in_the_tu() -> None:
    tu_genes = [["A", "B", "C"], ["D", "E"]]
    ab = {"A": 10.0, "B": 11.4, "C": 11.5, "D": 11.5, "E": 11.6}
    q = {("A", "B"): 0.02, ("B", "C"): 0.02, ("C", "D"): 0.04, ("D", "E"): 0.5}
    out, _pairs, flagged = _table(
        [tu("T1", tu_genes[0], 0, 900), tu("T2", tu_genes[1], 1000, 1600)], tu_genes, ab, q=q
    )
    # A->B and B->C are discordant OP pairs; C->D a co-abundant NO pair between the TUs
    assert sorted(flagged["flag"]) == [
        "novel_operon_candidate", "unexplained_discordance", "unexplained_discordance",
    ]  # fmt: skip
    by = _by_id(out)
    assert by.loc["T1", "n_flags"] == 2
    assert by.loc["T2", "n_flags"] == 0
    n_op = int((flagged["uniop_class"] == "OP").sum())
    assert int(out["n_flags"].sum()) <= n_op


def test_minus_strand_genes_keep_transcription_order() -> None:
    out, _, _ = _table(
        [tu("T1", ["C", "B", "A"], 0, 900, strand="-")], [["C", "B", "A"]], {"C": 5.0}
    )
    r = out.iloc[0]
    assert list(r["genes"]) == ["C", "B", "A"]
    assert all(isinstance(g, str) for g in r["genes"])
    assert r["strand"] == "-"


def test_tu_spanning_the_origin_keeps_unwrapped_coordinates() -> None:
    out, _, _ = _table([tu("T1", ["A", "B"], 9800, 10300)], [["A", "B"]], {"A": 1.0, "B": 2.0})
    r = out.iloc[0]
    assert (r["start"], r["end"]) == (9800, 10300)
    assert r["contig"] == "c1"


def test_one_row_per_tu_sorted_by_contig_order_then_start() -> None:
    rows = [
        tu("T3", ["E"], 500, 800, contig="zeta"),
        tu("T1", ["A", "B"], 100, 400, contig="zeta"),
        tu("T4", ["F"], 900, 1200, contig="alpha"),
        tu("T2", ["C", "D"], 50, 300, contig="alpha"),
    ]
    tu_genes = [["E"], ["A", "B"], ["F"], ["C", "D"]]
    out, _, _ = _table(rows, tu_genes, {"A": 1.0})
    assert list(out["tu_id"]) == ["T1", "T3", "T2", "T4"]
    assert len(out) == 4


def test_tu_genes_absent_from_pairs_still_get_a_row() -> None:
    t = tus(tu("T1", ["A"], 0, 300), tu("T2", ["B"], 400, 700))
    genes = gene_table({"A": 3.0}, ["A", "B"])
    pairs = chain_pairs([], {})
    flagged = fl.flag_pairs(pairs, cfg())
    out = tt.operons_table(t, pairs, genes, flagged)
    assert list(out["tu_id"]) == ["T1", "T2"]
    assert list(out["n_flags"]) == [0, 0]


def test_integer_columns() -> None:
    out, _, _ = _table([tu("T1", ["A", "B"], 0, 800)], [["A", "B"]], {"A": 1.0})
    assert pd.api.types.is_integer_dtype(out["n_genes"])
    assert pd.api.types.is_integer_dtype(out["n_flags"])
    assert pd.api.types.is_integer_dtype(out["start"])
    assert pd.api.types.is_integer_dtype(out["end"])


def test_empty_tus() -> None:
    pairs = chain_pairs([], {})
    flagged = fl.flag_pairs(pairs, cfg())
    out = tt.operons_table(tus(), pairs, gene_table({}, []), flagged)
    assert len(out) == 0
    assert tuple(out.columns) == OPERON_COLUMNS


def test_does_not_mutate_and_is_deterministic() -> None:
    tu_genes = [["A", "B", "C"], ["D", "E"]]
    ab = {"A": 10.0, "B": 11.4, "C": 11.5, "E": 11.6}
    t = tus(tu("T1", tu_genes[0], 0, 900), tu("T2", tu_genes[1], 1000, 1600))
    pairs = chain_pairs(tu_genes, ab, q={("A", "B"): 0.02})
    genes = gene_table(ab, [g for gs in tu_genes for g in gs])
    flagged = fl.flag_pairs(pairs, cfg())
    copies = [x.copy(deep=True) for x in (t, pairs, genes, flagged)]
    a = tt.operons_table(t, pairs, genes, flagged)
    b = tt.operons_table(t, pairs, genes, flagged)
    for x, c in zip((t, pairs, genes, flagged), copies, strict=True):
        pd.testing.assert_frame_equal(x, c)
    pd.testing.assert_frame_equal(a, b)
