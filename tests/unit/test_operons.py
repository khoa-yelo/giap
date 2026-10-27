"""Operon pairs, distance-model fallback and TUs (issue M2-12; GIPPA_SPEC.md §7.1.3, §9.3, §6.5).

Examples 1-3 and 5-9 of docs/specs/M2-12-uniop-pairs-tus.md, its edge cases, the TSV/FASTA
writers and ``calibration.yaml`` (M2-index decision 5).

Conventions pinned here where the spec is silent (see the PR description):

* ``build_pairs`` / ``build_tus`` / ``adjacent_pairs`` take replicons as a mapping in FASTA order
  (contig order = mapping order); input gene order does not matter.
* A circular contig with two same-strand genes has both the forward pair and the wrap pair; a
  contig with one gene has no pair (no self pair).
* ``uniop=[]`` (an empty ``uniop.pred``) is a UniOP failure exactly like ``uniop=None``.
* "Nearest upstream gene": for a + TU with first start codon at ``S`` every other gene on the same
  contig that starts upstream of ``S`` (circularly on a circular contig) clips the region at
  ``min(gene end, S)``; genes starting at or after ``S`` never clip (mirror image on the minus
  strand). Genes on other contigs never clip.
* ``tu_upstream`` on a circular contig returns an unwrapped interval (``start >= 0``,
  ``end > length`` across the origin, §6.5).
* ``write_tu_regions`` writes ``<outdir>/<tu_id>.fa`` (one record, header starting with the TU ID,
  strand-aware sequence of the upstream interval) and returns the paths in TU order.
"""

from __future__ import annotations

import math
from importlib import resources
from pathlib import Path
from typing import Any

import pytest
import yaml

from gippa.annotate.operons import (
    Tu,
    adjacent_pairs,
    build_pairs,
    build_tus,
    distance_model_p,
    intergenic_distance,
    tu_upstream,
    write_operon_pairs_tsv,
    write_tu_regions,
    write_tus_tsv,
)
from gippa.annotate.uniop_adapter import PairProb
from gippa.io.tsv import read_tsv
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

D0, S = 50.0, 20.0
FLAG = "operon_source_distance_model"
CHR = Replicon("chr1", 100_000, "circular", "chromosome")
LIN = Replicon("lin1", 100_000, "linear", "unknown")
SMALL_LIN = Replicon("lin2", 5_000, "linear", "unknown")
PL = Replicon("p2", 8_000, "circular", "plasmid")
TUS_HEADER = [
    "tu_id", "contig", "start", "end", "strand", "genes", "upstream_start", "upstream_end",
]  # fmt: skip
PAIRS_HEADER = ["gene_a", "gene_b", "uniop_p", "source"]
_COMP = str.maketrans("ACGT", "TGCA")


def gene(
    tag: str, start: int, end: int, strand: str, contig: str = "chr1", ft: str = "cds"
) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type=ft,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def reps(*rs: Replicon) -> dict[str, Replicon]:
    return {r.id: r for r in rs}


def tags(pairs: Any) -> list[tuple[str, str]]:
    return [(a.locus_tag, b.locus_tag) for a, b in pairs]


def logistic(d: float) -> float:
    return 1.0 / (1.0 + math.exp((d - D0) / S))


# ---------------------------------------------------------------- examples 1-3: distance model


def test_example1_distance_and_p() -> None:
    a, b = gene("A", 100, 400, "+"), gene("B", 420, 700, "+")
    assert intergenic_distance(a, b, CHR) == 20
    assert distance_model_p(20, D0, S) == pytest.approx(0.817574, abs=1e-6)


def test_example2_overlap_is_negative_distance() -> None:
    a, b = gene("A", 100, 400, "+"), gene("B", 396, 700, "+")
    assert intergenic_distance(a, b, CHR) == -4
    assert distance_model_p(-4, D0, S) == pytest.approx(0.937027, abs=1e-6)


def test_example3_p_half_at_d0_and_linked_at_threshold() -> None:
    assert distance_model_p(50, D0, S) == pytest.approx(0.5, abs=1e-12)
    a, b = gene("A", 100, 400, "+", "lin1"), gene("B", 450, 700, "+", "lin1")
    pairs, _ = build_pairs([a, b], reps(LIN), None, D0, S)
    assert [(p.gene_a, p.gene_b) for p in pairs] == [("A", "B")]
    assert pairs[0].uniop_p == pytest.approx(0.5, abs=1e-12)
    tus = build_tus([a, b], [PairProb("A", "B", 0.5, "uniop")], 0.5, reps(LIN))
    assert [t.genes for t in tus] == [("A", "B")]  # p >= threshold links


def test_distance_model_extreme_distances_do_not_overflow() -> None:
    far = distance_model_p(10**7, D0, S)
    near = distance_model_p(-(10**7), D0, S)
    assert 0.0 <= far <= 1e-12
    assert 1.0 - 1e-12 <= near <= 1.0


def test_distance_model_uses_given_parameters() -> None:
    assert distance_model_p(100, 100.0, 5.0) == pytest.approx(0.5)
    assert distance_model_p(110, 100.0, 5.0) == pytest.approx(1 / (1 + math.e**2))


# ---------------------------------------------------------------- adjacent pairs


def test_adjacent_pairs_same_strand_consecutive_in_genome_order() -> None:
    g = [
        gene("D", 1250, 1500, "-", "lin1"),
        gene("A", 100, 400, "+", "lin1"),
        gene("C", 900, 1200, "-", "lin1"),
        gene("B", 420, 700, "+", "lin1"),
        gene("E", 1600, 1900, "+", "lin1"),
    ]
    assert tags(adjacent_pairs(g, reps(LIN))) == [("A", "B"), ("C", "D")]


def test_adjacent_pairs_never_cross_contigs() -> None:
    g = [gene("A", 100, 400, "+", "lin1"), gene("B", 100, 400, "+", "lin2")]
    assert adjacent_pairs(g, reps(LIN, SMALL_LIN)) == []


def test_adjacent_pairs_follow_fasta_contig_order() -> None:
    g = [
        gene("P1", 10, 300, "+", "p2"),
        gene("P2", 400, 700, "+", "p2"),
        gene("A", 100, 400, "+", "lin1"),
        gene("B", 420, 700, "+", "lin1"),
    ]
    assert tags(adjacent_pairs(g, reps(LIN, PL)))[:1] == [("A", "B")]
    assert tags(adjacent_pairs(g, reps(PL, LIN)))[:1] == [("P1", "P2")]


# ---------------------------------------------------------------- example 8: wrap pair


def test_example8_wrap_pair_on_circular_contig() -> None:
    first, mid, last = (
        gene("F", 30, 300, "+"),
        gene("M", 50_000, 50_300, "-"),
        gene("L", 99_800, 99_950, "+"),
    )
    assert tags(adjacent_pairs([first, mid, last], reps(CHR))) == [("L", "F")]
    assert intergenic_distance(last, first, CHR) == 80
    pairs, _ = build_pairs([first, mid, last], reps(CHR), None, D0, S)
    assert [(p.gene_a, p.gene_b) for p in pairs] == [("L", "F")]  # gene_a = last gene
    assert pairs[0].uniop_p == pytest.approx(logistic(80))


def test_example8_no_wrap_pair_on_linear_contig() -> None:
    g = [gene("F", 30, 300, "+", "lin1"), gene("L", 99_800, 99_950, "+", "lin1")]
    assert tags(adjacent_pairs(g, reps(LIN))) == [("F", "L")]
    g = [g[0], gene("M", 50_000, 50_300, "-", "lin1"), g[1]]
    assert adjacent_pairs(g, reps(LIN)) == []


def test_two_gene_circular_contig_has_forward_and_wrap_pair() -> None:
    a, b = gene("A", 100, 400, "+", "p2"), gene("B", 420, 700, "+", "p2")
    assert tags(adjacent_pairs([b, a], reps(PL))) == [("A", "B"), ("B", "A")]
    assert intergenic_distance(b, a, PL) == 100 + 8_000 - 700


def test_intergenic_distance_with_origin_crossing_gene() -> None:
    # an origin-crossing gene is stored unwrapped (end > length) and is the last gene by start
    x = gene("X", 99_500, 99_800, "+")
    o = gene("O", 99_900, 100_200, "+")
    b = gene("B", 300, 600, "+")
    assert intergenic_distance(x, o, CHR) == 100
    assert intergenic_distance(o, b, CHR) == 100  # 300 - (100_200 - 100_000)
    assert tags(adjacent_pairs([b, x, o], reps(CHR))) == [("B", "X"), ("X", "O"), ("O", "B")]


def test_single_gene_contig_has_no_pairs_and_one_tu() -> None:
    a = gene("A", 100, 400, "+")
    assert adjacent_pairs([a], reps(CHR)) == []
    pairs, _ = build_pairs([a], reps(CHR), None, D0, S)
    assert pairs == []
    tus = build_tus([a], pairs, 0.5, reps(CHR))
    assert [(t.tu_id, t.genes) for t in tus] == [("TU_1", ("A",))]


def test_no_genes_gives_empty_outputs() -> None:
    assert adjacent_pairs([], reps(CHR)) == []
    pairs, _ = build_pairs([], reps(CHR), None, D0, S)
    assert pairs == []
    assert build_tus([], [], 0.5, reps(CHR)) == []


# ---------------------------------------------------------------- example 9 and UniOP merge


def _abcd(contig: str = "chr1") -> list[Gene]:
    return [
        gene("A", 100, 400, "+", contig),
        gene("B", 420, 700, "+", contig),
        gene("C", 900, 1200, "-", contig),
        gene("D", 1250, 1500, "-", contig),
    ]


def test_example9_uniop_failure_all_distance_model_and_flag() -> None:
    g = _abcd()
    pairs, flags = build_pairs(g, reps(CHR), None, D0, S)
    assert flags == (FLAG,)
    assert [(p.gene_a, p.gene_b) for p in pairs] == [("A", "B"), ("C", "D")]
    assert all(p.source == "distance_model" for p in pairs)
    assert pairs[0].uniop_p == pytest.approx(logistic(20))
    assert pairs[1].uniop_p == pytest.approx(logistic(50))


def test_empty_uniop_pred_counts_as_failure() -> None:
    pairs, flags = build_pairs(_abcd(), reps(CHR), [], D0, S)
    assert flags == (FLAG,)
    assert all(p.source == "distance_model" for p in pairs)


def test_uniop_used_for_reported_pairs_distance_model_for_the_rest() -> None:
    g = _abcd()
    uni = [PairProb("A", "B", 0.12, "uniop")]
    pairs, flags = build_pairs(g, reps(CHR), uni, D0, S)
    assert flags == ()
    got = {(p.gene_a, p.gene_b): p for p in pairs}
    assert set(got) == {("A", "B"), ("C", "D")}
    assert got["A", "B"] == PairProb("A", "B", 0.12, "uniop")
    assert got["C", "D"].source == "distance_model"
    assert got["C", "D"].uniop_p == pytest.approx(logistic(50))


def test_uniop_never_covers_the_wrap_pair() -> None:
    first, last = gene("F", 30, 300, "+"), gene("L", 99_800, 99_950, "+")
    mid = gene("M", 400, 700, "+")
    uni = [PairProb("F", "M", 0.8, "uniop"), PairProb("M", "L", 0.7, "uniop")]
    pairs, flags = build_pairs([first, mid, last], reps(CHR), uni, D0, S)
    assert flags == ()
    assert [(p.gene_a, p.gene_b, p.source) for p in pairs] == [
        ("F", "M", "uniop"),
        ("M", "L", "uniop"),
        ("L", "F", "distance_model"),
    ]
    assert pairs[2].uniop_p == pytest.approx(logistic(80))


# ---------------------------------------------------------------- example 5: TU assembly


def test_example5_tus_split_by_strand_and_threshold() -> None:
    pairs = [PairProb("A", "B", 0.9, "uniop"), PairProb("C", "D", 0.3, "uniop")]
    tus = build_tus(_abcd(), pairs, 0.5, reps(CHR))
    assert [(t.tu_id, t.genes, t.strand) for t in tus] == [
        ("TU_1", ("A", "B"), "+"),
        ("TU_2", ("C",), "-"),
        ("TU_3", ("D",), "-"),
    ]
    assert tus[0].interval == Interval("chr1", 100, 700, "+")
    assert tus[1].interval == Interval("chr1", 900, 1200, "-")
    assert all(t.contig == "chr1" for t in tus)


def test_example5_minus_strand_tu_lists_genes_in_transcription_order() -> None:
    pairs = [PairProb("A", "B", 0.9, "uniop"), PairProb("C", "D", 0.8, "uniop")]
    tus = build_tus(_abcd(), pairs, 0.5, reps(CHR))
    assert [t.genes for t in tus] == [("A", "B"), ("D", "C")]
    assert tus[1].interval == Interval("chr1", 900, 1500, "-")


def test_tu_threshold_is_inclusive_and_configurable() -> None:
    pairs = [PairProb("A", "B", 0.7, "uniop")]
    g = _abcd()[:2]
    assert [t.genes for t in build_tus(g, pairs, 0.7, reps(CHR))] == [("A", "B")]
    assert [t.genes for t in build_tus(g, pairs, 0.71, reps(CHR))] == [("A",), ("B",)]


def test_tu_runs_are_maximal_chains() -> None:
    g = [gene(t, 100 + 400 * i, 400 + 400 * i, "+", "lin1") for i, t in enumerate("ABCDE")]
    pairs = [
        PairProb("A", "B", 0.9, "uniop"),
        PairProb("B", "C", 0.6, "uniop"),
        PairProb("C", "D", 0.1, "uniop"),
        PairProb("D", "E", 0.95, "uniop"),
    ]
    tus = build_tus(g, pairs, 0.5, reps(LIN))
    assert [t.genes for t in tus] == [("A", "B", "C"), ("D", "E")]


def test_tu_numbering_by_contig_then_start_then_strand() -> None:
    g = [
        gene("Q", 50, 300, "+", "p2"),
        gene("Y", 100, 350, "-", "lin1"),
        gene("X", 100, 400, "+", "lin1"),
        gene("Z", 20, 90, "-", "lin1"),
    ]
    tus = build_tus(g, [], 0.5, reps(LIN, PL))
    assert [(t.tu_id, t.genes) for t in tus] == [
        ("TU_1", ("Z",)),
        ("TU_2", ("X",)),  # same start as Y: + before -
        ("TU_3", ("Y",)),
        ("TU_4", ("Q",)),
    ]


def test_pseudogenes_and_sorfs_take_part_in_tus() -> None:
    g = [
        gene("A", 100, 400, "+", "lin1"),
        gene("S", 420, 480, "+", "lin1", ft="sorf"),
        gene("P", 500, 900, "+", "lin1", ft="pseudogene"),
    ]
    assert tags(adjacent_pairs(g, reps(LIN))) == [("A", "S"), ("S", "P")]
    pairs, _ = build_pairs(g, reps(LIN), None, D0, S)
    assert [t.genes for t in build_tus(g, pairs, 0.5, reps(LIN))] == [("A", "S", "P")]


def test_tu_spanning_the_origin_is_unwrapped() -> None:
    first, mid, last = (
        gene("F", 30, 300, "+"),
        gene("M", 50_000, 50_300, "-"),
        gene("L", 99_800, 99_950, "+"),
    )
    tus = build_tus([first, mid, last], [PairProb("L", "F", 0.9, "uniop")], 0.5, reps(CHR))
    assert [(t.tu_id, t.genes) for t in tus] == [("TU_1", ("M",)), ("TU_2", ("L", "F"))]
    assert tus[1].interval == Interval("chr1", 99_800, 100_300, "+")


def test_build_tus_uses_config_upstream_bp() -> None:
    a = gene("A", 1000, 1300, "+")
    (tu,) = build_tus([a], [], 0.5, reps(CHR), upstream_bp=(-100, 20))
    assert tu.upstream == Interval("chr1", 900, 1020, "+")


# ---------------------------------------------------------------- examples 6, 7: upstream regions


def test_example6_plus_upstream_clipped_at_nearest_upstream_gene() -> None:
    first, other = gene("A", 1000, 1300, "+"), gene("X", 500, 850, "-")
    assert tu_upstream(first, [other], CHR) == Interval("chr1", 850, 1050, "+")
    (tu, _) = build_tus([first, other], [], 0.5, reps(CHR))[::-1]
    assert tu.genes == ("A",) and tu.upstream == Interval("chr1", 850, 1050, "+")


def test_example7_minus_upstream_clipped_at_nearest_upstream_gene() -> None:
    first, other = gene("A", 5000, 5600, "-"), gene("Y", 5800, 6200, "+")
    assert tu_upstream(first, [other], CHR) == Interval("chr1", 5550, 5800, "-")


def test_upstream_unclipped_full_window() -> None:
    assert tu_upstream(gene("A", 1000, 1300, "+"), [], CHR) == Interval("chr1", 700, 1050, "+")
    assert tu_upstream(gene("A", 5000, 5600, "-"), [], CHR) == Interval("chr1", 5550, 5900, "-")


def test_upstream_clipped_at_the_nearest_of_several_genes() -> None:
    first = gene("A", 1000, 1300, "+")
    others = [gene("X", 200, 900, "+"), gene("W", 500, 850, "-"), gene("V", 100, 300, "+")]
    assert tu_upstream(first, others, CHR) == Interval("chr1", 900, 1050, "+")


def test_upstream_ignores_downstream_genes_and_other_contigs() -> None:
    first = gene("A", 1000, 1300, "+")
    others = [gene("D", 1020, 1500, "-"), gene("O", 800, 990, "+", contig="lin1")]
    assert tu_upstream(first, others, CHR) == Interval("chr1", 700, 1050, "+")


def test_q2_upstream_gene_overlapping_the_start_codon_plus() -> None:
    first = gene("A", 1000, 1300, "+")
    assert tu_upstream(first, [gene("X", 900, 1010, "+")], CHR) == Interval("chr1", 1000, 1050, "+")
    assert tu_upstream(first, [gene("Y", 600, 1000, "-")], CHR) == Interval("chr1", 1000, 1050, "+")
    assert tu_upstream(first, [gene("Z", 400, 2000, "-")], CHR) == Interval("chr1", 1000, 1050, "+")


def test_q2_upstream_gene_overlapping_the_start_codon_minus() -> None:
    first = gene("A", 5000, 5600, "-")
    assert tu_upstream(first, [gene("X", 5590, 5900, "-")], CHR) == Interval(
        "chr1", 5550, 5600, "-"
    )


def test_upstream_clipped_at_linear_contig_ends() -> None:
    assert tu_upstream(gene("A", 100, 400, "+", "lin2"), [], SMALL_LIN) == Interval(
        "lin2", 0, 150, "+"
    )
    assert tu_upstream(gene("B", 4700, 4900, "-", "lin2"), [], SMALL_LIN) == Interval(
        "lin2", 4850, 5000, "-"
    )


def test_upstream_wraps_on_circular_contig() -> None:
    assert tu_upstream(gene("A", 100, 400, "+"), [], CHR) == Interval("chr1", 99_800, 100_150, "+")
    assert tu_upstream(gene("B", 99_500, 99_900, "-"), [], CHR) == Interval(
        "chr1", 99_850, 100_200, "-"
    )


def test_upstream_wrap_clipped_by_gene_across_the_origin() -> None:
    first = gene("A", 100, 400, "+")
    assert tu_upstream(first, [gene("Z", 99_000, 99_950, "-")], CHR) == Interval(
        "chr1", 99_950, 100_150, "+"
    )
    first = gene("B", 99_500, 99_900, "-")
    assert tu_upstream(first, [gene("Y", 50, 400, "+")], CHR) == Interval(
        "chr1", 99_850, 100_050, "-"
    )


def test_upstream_custom_window() -> None:
    assert tu_upstream(gene("A", 1000, 1300, "+"), [], CHR, bp=(-100, 20)) == Interval(
        "chr1", 900, 1020, "+"
    )
    assert tu_upstream(gene("A", 5000, 5600, "-"), [], CHR, bp=(-100, 20)) == Interval(
        "chr1", 5580, 5700, "-"
    )


# ---------------------------------------------------------------- writers


def test_operon_pairs_tsv(tmp_path: Path) -> None:
    pairs = [PairProb("A", "B", 0.91, "uniop"), PairProb("C", "D", 0.5, "distance_model")]
    path = tmp_path / "operon_pairs.tsv"
    write_operon_pairs_tsv(path, pairs)
    header, rows = read_tsv(path)
    assert header == PAIRS_HEADER
    assert [(r[0], r[1], float(r[2]), r[3]) for r in rows] == [  # type: ignore[arg-type]
        ("A", "B", 0.91, "uniop"),
        ("C", "D", 0.5, "distance_model"),
    ]


def test_tus_tsv_is_one_based_closed_with_upstream(tmp_path: Path) -> None:
    g = [gene("A", 1000, 1300, "+"), gene("X", 500, 850, "-"), gene("B", 1320, 1600, "+")]
    g += [gene("C", 5000, 5600, "-"), gene("Y", 5800, 6200, "+")]
    pairs = [PairProb("A", "B", 0.9, "uniop")]
    tus = build_tus(g, pairs, 0.5, reps(CHR))
    path = tmp_path / "tus.tsv"
    write_tus_tsv(path, tus)
    header, rows = read_tsv(path)
    assert header == TUS_HEADER
    by_genes = {r[5]: r for r in rows}
    assert by_genes["A;B"] == ["TU_2", "chr1", "1001", "1600", "+", "A;B", "851", "1050"]
    assert by_genes["C"] == ["TU_3", "chr1", "5001", "5600", "-", "C", "5551", "5800"]


def test_tus_tsv_origin_spanning_tu_written_unwrapped(tmp_path: Path) -> None:
    g = [gene("F", 30, 300, "+"), gene("L", 99_800, 99_950, "+")]
    tus = build_tus(g, [PairProb("L", "F", 0.9, "uniop")], 0.5, reps(CHR))
    path = tmp_path / "tus.tsv"
    write_tus_tsv(path, tus)
    _, rows = read_tsv(path)
    assert [r[:6] for r in rows] == [["TU_1", "chr1", "99801", "100300", "+", "L;F"]]
    # upstream of L: [99_500, 99_850); F lies downstream of L across the origin and never clips
    assert rows[0][6:] == ["99501", "99850"]


def test_empty_tables_have_headers(tmp_path: Path) -> None:
    write_operon_pairs_tsv(tmp_path / "p.tsv", [])
    write_tus_tsv(tmp_path / "t.tsv", [])
    assert read_tsv(tmp_path / "p.tsv") == (PAIRS_HEADER, [])
    assert read_tsv(tmp_path / "t.tsv") == (TUS_HEADER, [])


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def test_write_tu_regions_one_fasta_per_tu(tmp_path: Path) -> None:
    import numpy as np

    rng = np.random.default_rng(12)
    seq = "".join(rng.choice(list("ACGT"), CHR.length))
    g = [
        gene("A", 1000, 1300, "+"),
        gene("X", 500, 850, "-"),
        gene("C", 5000, 5600, "-"),
        gene("Y", 5800, 6200, "+"),
        gene("W", 100, 400, "+"),
        gene("L", 99_000, 99_200, "-"),
    ]
    # W's region wraps the origin, L's region on the minus strand extends to 99_500
    tus = build_tus(g, [], 0.5, reps(CHR))
    out = tmp_path / "tu_upstream"
    out.mkdir()
    paths = write_tu_regions(out, tus, {"chr1": seq})
    assert paths == [out / f"{t.tu_id}.fa" for t in tus]
    double = seq + seq
    for t, p in zip(tus, paths, strict=True):
        lines = p.read_text().splitlines()
        assert lines[0].startswith(f">{t.tu_id}")
        got = "".join(lines[1:])
        iv = t.upstream
        want = double[iv.start : iv.end]
        assert got == (want if iv.strand == "+" else _rc(want)), t.tu_id
    by_gene = {t.genes[0]: t for t in tus}
    assert by_gene["W"].upstream == Interval("chr1", 99_800, 100_150, "+")


def test_tu_is_frozen_dataclass() -> None:
    tu = Tu("TU_1", "chr1", Interval("chr1", 0, 10, "+"), "+", ("A",), Interval("chr1", 0, 5, "+"))
    with pytest.raises(AttributeError):
        tu.tu_id = "TU_2"  # type: ignore[misc]


# ---------------------------------------------------------------- calibration.yaml


def _get(obj: Any, key: str) -> Any:
    return obj[key] if isinstance(obj, dict) else getattr(obj, key)


#: GIPPA_SPEC §7.1.3 / §15.4 placeholder values, binding only while an entry is uncalibrated.
PLACEHOLDERS = {"operon_distance": {"d0": 50, "s": 20}, "terminatornet_yes": {"value": 0.5}}
#: calibration.yaml entries mirrored by a defaults.yaml key (M6-12 decision 2: the mirror plus
#: this agreement test; defaults.yaml stays the value the code reads).
MIRRORS = {
    "terminatornet_yes": ("annotation", "terminators", "terminatornet_yes"),
    "uniop_pair_threshold": ("annotation", "uniop", "pair_threshold"),
}


def test_calibration_yaml_entries() -> None:
    """Placeholders while uncalibrated; once calibrated (ADR-0022) an entry names its
    provenance and keeps finite values in range."""
    text = resources.files("gippa.data").joinpath("calibration.yaml").read_text()
    data = yaml.safe_load(text)
    assert {"operon_distance", "terminatornet_yes"} <= set(data)
    for key, entry in data.items():
        assert isinstance(entry["calibrated"], bool), key
        # string metadata of a calibrated entry (ADR-0022)
        meta = ("calibrated", "provenance", "regulondb_version", "date", "rule")
        values = {k: v for k, v in entry.items() if k not in meta}
        assert values and all(math.isfinite(float(v)) for v in values.values()), key
        if not entry["calibrated"]:
            assert values == PLACEHOLDERS.get(key, values), key
        else:
            assert entry.get("provenance"), key
    od = data["operon_distance"]
    assert od["s"] > 0
    assert 0.0 <= data["terminatornet_yes"]["value"] <= 1.0
    assert 0.0 <= data["uniop_pair_threshold"]["value"] <= 1.0


def test_calibration_loader() -> None:
    from gippa.data import calibration

    cal = calibration.load()
    text = resources.files("gippa.data").joinpath("calibration.yaml").read_text()
    assert cal == yaml.safe_load(text)
    assert float(_get(_get(cal, "operon_distance"), "d0")) > 0


def test_calibration_and_defaults_agree() -> None:
    """A calibrated entry with a defaults.yaml mirror has the same value there (vacuous while
    every mirrored entry is a placeholder)."""
    from gippa.data import calibration

    text = resources.files("gippa.config").joinpath("defaults.yaml").read_text()
    defaults = yaml.safe_load(text)
    for key, entry in calibration.load().items():
        if key not in MIRRORS or not entry["calibrated"]:
            continue
        node: Any = defaults
        for part in MIRRORS[key]:
            node = node[part]
        assert float(node) == float(entry["value"]), key
