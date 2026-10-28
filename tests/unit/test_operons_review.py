"""M2-12 review mutants on the pure operon functions (lead decision 6 after PR #49 review).

Each test kills one surviving mutant: the first gene of a minus-strand multi-gene TU is the
transcription-first gene (highest coordinates, ``run[-1]`` in genome order), the ring-length cap on
a TU that closes the ring, and the first gene bounding its own upstream window on a small ring.
"""

from __future__ import annotations

from pathlib import Path

from gippa.annotate.operons import build_tus, tu_upstream, write_tus_tsv
from gippa.annotate.uniop_adapter import PairProb
from gippa.io.tsv import read_tsv
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

CHR = Replicon("chr1", 100_000, "circular", "chromosome")


def gene(tag: str, start: int, end: int, strand: str, contig: str = "chr1") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def test_minus_strand_tu_upstream_is_around_the_transcription_first_gene(tmp_path: Path) -> None:
    a, b = gene("A", 1000, 1500, "-"), gene("B", 1600, 2000, "-")
    x = gene("X", 2500, 3000, "+")
    pairs = [PairProb("A", "B", 0.9, "uniop")]
    tus = build_tus([a, b, x], pairs, 0.5, {"chr1": CHR})
    tu = next(t for t in tus if "A" in t.genes)
    assert tu.genes == ("B", "A")
    assert tu.interval == Interval("chr1", 1000, 2000, "-")
    # start codon of B is 2000 (minus strand): window [2000-50, 2000+300), X starts at 2500
    assert tu.upstream == Interval("chr1", 1950, 2300, "-")
    assert tu.upstream == tu_upstream(b, [a, x], CHR)
    write_tus_tsv(tmp_path / "tus.tsv", tus)
    header, rows = read_tsv(tmp_path / "tus.tsv")
    row = dict(zip(header, next(r for r in rows if "A" in str(r[5])), strict=True))
    assert (int(row["upstream_start"]), int(row["upstream_end"])) == (1951, 2300)  # type: ignore[arg-type]


def test_minus_strand_upstream_clipped_by_gene_beyond_the_first_gene() -> None:
    """The clip comes from genes upstream of B (the TU's first gene), not of A."""
    a, b = gene("A", 1000, 1500, "-"), gene("B", 1600, 2000, "-")
    y = gene("Y", 2100, 2600, "-")
    tus = build_tus([a, b, y], [PairProb("A", "B", 0.9, "uniop")], 0.5, {"chr1": CHR})
    tu = next(t for t in tus if "A" in t.genes)
    assert tu.upstream == Interval("chr1", 1950, 2100, "-")


def test_tu_closing_the_ring_is_capped_at_the_ring_length() -> None:
    ring = Replicon("r", 1000, "circular", "plasmid")
    a = gene("A", 50, 400, "+", "r")
    b = gene("B", 900, 1100, "+", "r")  # crosses the origin and overlaps A by 50 bp
    pairs = [PairProb("A", "B", 0.9, "uniop"), PairProb("B", "A", 0.9, "uniop")]
    (tu,) = build_tus([a, b], pairs, 0.5, {"r": ring})
    assert tu.genes == ("A", "B")
    assert tu.interval == Interval("r", 50, 1050, "+")
    assert tu.interval.end - tu.interval.start == ring.length


def test_ring_self_bound_plus() -> None:
    """On a 500-bp ring the -300 window would reach the gene's own 3' end across the origin."""
    ring = Replicon("r", 500, "circular", "plasmid")
    g = gene("G", 100, 400, "+", "r")
    assert tu_upstream(g, [], ring) == Interval("r", 400, 650, "+")
    (tu,) = build_tus([g], [], 0.5, {"r": ring})
    assert tu.upstream == Interval("r", 400, 650, "+")


def test_ring_self_bound_minus() -> None:
    ring = Replicon("r", 500, "circular", "plasmid")
    g = gene("G", 100, 400, "-", "r")
    assert tu_upstream(g, [], ring) == Interval("r", 350, 600, "-")
    (tu,) = build_tus([g], [], 0.5, {"r": ring})
    assert tu.upstream == Interval("r", 350, 600, "-")


def test_upstream_bp_argument_is_used_by_build_tus() -> None:
    g = gene("G", 1000, 1300, "+")
    (tu,) = build_tus([g], [], 0.5, {"chr1": CHR}, upstream_bp=(-100, 20))
    assert tu.upstream == Interval("chr1", 900, 1020, "+")
