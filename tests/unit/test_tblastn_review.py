"""tblastn route: lead decisions after code review (PR #53; docs/specs/M3-13-tblastn-route.md).

1. ``makeblastdb`` runs WITHOUT ``-parse_seqids`` (builder and main.nf); a real BLAST+ 2.16.0
   run on a contig named ``NC_000913.3`` (``tests/data/tblastn/k12_accession``) keeps the ID
   unchanged and chains without "unknown contig". The recorded ``-parse_seqids`` contrast run
   shows why: the subject becomes ``ref|NC_000913.3|``.
2. Chaining picks the BEST chain per protein x contig x strand (dynamic program maximising the
   summed bitscore; ties -> larger union query coverage, then lower start), not the nearest
   extension of a greedy seed. The reviewer's A/B/C case -> chain A+B, coverage 1.0, assigned.
3. Mutant killers: best-chain choice, a chain longer than a small circular contig (extension
   at the end and at the start), the ``b0 > a0`` genome order check, ``ValueError`` for hits
   outside the contig and for bad query ranges.
4. The parser rejects rows that do not have exactly 15 columns.

Ambiguity resolved here: "lower start" compares the chains' canonical interval starts (0-based).
"""

from __future__ import annotations

import itertools
import random
import re
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

from giap.config.model import ProteinsConfig
from giap.model import Interval, Replicon
from giap.model.coords import check_interval
from giap.model.genes import Gene, GeneModel

REPO = Path(__file__).resolve().parents[2]
K12 = REPO / "tests" / "data" / "tblastn" / "k12_accession"
MAIN_NF = REPO / "workflow" / "modules" / "local" / "tblastn" / "main.nf"
CFG = ProteinsConfig()

L = 50_000
CIRC = {"chr1": Replicon("chr1", L, "circular", "chromosome")}
LIN = {"chr1": Replicon("chr1", L, "linear", "unknown")}


def _tb():  # type: ignore[no-untyped-def]
    from giap.mapping import tblastn

    return tblastn


def _hsp(
    q: tuple[int, int],
    s: tuple[int, int],
    *,
    bitscore: float = 100.0,
    pid: str = "P1",
    contig: str = "chr1",
    qlen: int = 120,
):  # type: ignore[no-untyped-def]
    n = q[1] - q[0] + 1
    return _tb().Hsp(
        qseqid=pid,
        sseqid=contig,
        qstart=q[0],
        qend=q[1],
        sstart=s[0],
        send=s[1],
        length=n,
        nident=n,
        positive=n,
        bitscore=bitscore,
        qlen=qlen,
    )


def _merge(hsps, reps=CIRC, **kw):  # type: ignore[no-untyped-def]
    return _tb().merge_hsps(hsps, reps, **kw)


def _by_rank(alns, pid: str = "P1"):  # type: ignore[no-untyped-def]
    return sorted((a for a in alns if a.protein_id == pid), key=lambda a: a.rank)


def _gene(tag: str, start: int, end: int, contig: str = "chr1", strand: str = "+") -> Gene:
    n_aa = (end - start) // 3
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" + "A" * max(0, n_aa - 2),
        nt=None,
        start_codon="ATG",
        partial=False,
    )


def _assign(alns, gm):  # type: ignore[no-untyped-def]
    from giap.mapping.protein_assign import assign_alignments

    return assign_alignments(alns, gm, CFG)


# ---------------------------------------------------------------- 1. no -parse_seqids


def test_makeblastdb_cmd_has_no_parse_seqids() -> None:
    cmd = _tb().makeblastdb_cmd(Path("genome.fa"), Path("db/genome"))
    assert cmd == ["makeblastdb", "-in", "genome.fa", "-dbtype", "nucl", "-out", "db/genome"]


def test_main_nf_has_no_parse_seqids() -> None:
    text = re.sub(r"/\*.*?\*/", "", MAIN_NF.read_text(), flags=re.S)
    text = re.sub(r"(?m)^\s*//.*$", "", text)
    assert "-parse_seqids" not in text, "main.nf must build the BLAST db without -parse_seqids"
    assert re.search(r"makeblastdb\s+-in\s+\S+\s+-dbtype\s+nucl\s+-out\s+db/genome", text)


K12_REPS = {"NC_000913.3": Replicon("NC_000913.3", 3_000, "circular", "chromosome")}


def test_k12_fixture_parses_accession_id_unchanged() -> None:
    hsps = _tb().parse_tblastn(K12 / "tblastn.tsv")
    assert [(h.qseqid, h.sseqid) for h in hsps] == [
        ("P37005", "NC_000913.3"),
        ("P0AD86", "NC_000913.3"),
    ]
    p = hsps[0]
    assert (p.qstart, p.qend, p.sstart, p.send, p.qlen) == (1, 228, 942, 1625, 228)
    assert (p.length, p.nident, p.positive, p.bitscore) == (228, 228, 228, 466.0)


def test_k12_fixture_chains_and_assigns_without_unknown_contig() -> None:
    alns = _merge(_tb().parse_tblastn(K12 / "tblastn.tsv"), K12_REPS)
    got = {a.protein_id: a.interval for a in alns}
    assert got == {
        "P37005": Interval("NC_000913.3", 941, 1625, "+"),
        "P0AD86": Interval("NC_000913.3", 1841, 1904, "+"),
    }
    assert all(a.query_coverage == pytest.approx(1.0) for a in alns)
    gm = GeneModel(
        K12_REPS,
        (
            _gene("b0727x", 941, 1628, contig="NC_000913.3"),  # with stop codon
            _gene("b0728x", 1841, 1907, contig="NC_000913.3"),
        ),
        (),
    )
    df, left = _assign(alns, gm)
    assert dict(zip(df["protein_id"], df["gene"], strict=True)) == {
        "P0AD86": "b0728x",
        "P37005": "b0727x",
    }
    assert left == set()


def test_k12_parse_seqids_contrast_does_not_match_contig() -> None:
    """With ``-parse_seqids`` BLAST reports ``ref|NC_000913.3|``: why the flag is dropped."""
    hsps = _tb().parse_tblastn(K12 / "tblastn_parse_seqids.tsv")
    assert {h.sseqid for h in hsps} == {"ref|NC_000913.3|"}
    with pytest.raises(ValueError, match="unknown contig"):
        _merge(hsps, K12_REPS)


# ---------------------------------------------------------------- 2. best chain


def _abc():  # type: ignore[no-untyped-def]
    return [
        _hsp((1, 60), (1001, 1180), bitscore=120.0),  # A
        _hsp((61, 120), (1184, 1363), bitscore=120.0),  # B
        _hsp((61, 75), (1181, 1225), bitscore=20.0),  # C: nearer to A, but a worse chain
    ]


@pytest.mark.parametrize("order", list(itertools.permutations(range(3))))
def test_reviewer_abc_best_chain_is_a_plus_b(order: tuple[int, ...]) -> None:
    hsps = _abc()
    alns = _by_rank(_merge([hsps[i] for i in order]))
    assert len(alns) == 2
    best, rest = alns
    assert best.interval == Interval("chr1", 1000, 1363, "+")
    assert best.query_coverage == pytest.approx(1.0)
    assert (best.q_start, best.q_end) == (0, 120)
    assert rest.interval == Interval("chr1", 1180, 1225, "+")
    assert (rest.q_start, rest.q_end) == (60, 75)


def test_reviewer_abc_gene_assigned() -> None:
    gm = GeneModel(CIRC, (_gene("G1", 1000, 1363),), ())
    df, left = _assign(_merge(_abc()), gm)
    r1 = df[df["rank"] == 1]
    assert list(r1["gene"]) == ["G1"]
    assert list(r1["start"]) == [1000] and list(r1["end"]) == [1363]
    assert left == set()


def test_best_chain_beats_greedy_multi_hsp_extension() -> None:
    """A+B (180) beats A+C+D (160) although C is nearer and C+D extend A to the end."""
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=100.0),  # A
        _hsp((61, 120), (1184, 1363), bitscore=80.0),  # B
        _hsp((61, 75), (1181, 1225), bitscore=30.0),  # C
        _hsp((76, 120), (1226, 1360), bitscore=30.0),  # D
    ]
    best, rest = _by_rank(_merge(hsps))
    assert best.interval == Interval("chr1", 1000, 1363, "+")
    assert rest.interval == Interval("chr1", 1180, 1360, "+")  # C+D chain
    assert rest.query_coverage == pytest.approx(60 / 120)


def test_best_chain_minus_strand() -> None:
    """Mirror of A/B/C on the minus strand (N-terminus at the high coordinate)."""
    hsps = [
        _hsp((1, 60), (1363, 1184), bitscore=120.0),  # A
        _hsp((61, 120), (1180, 1001), bitscore=120.0),  # B
        _hsp((61, 75), (1183, 1139), bitscore=20.0),  # C, adjacent to A
    ]
    best, rest = _by_rank(_merge(hsps))
    assert best.interval == Interval("chr1", 1000, 1363, "-")
    assert best.query_coverage == pytest.approx(1.0)
    assert rest.interval == Interval("chr1", 1138, 1183, "-")


def test_tie_on_bitscore_larger_union_coverage_wins() -> None:
    """X+Y and X+Z both sum to 150; X+Y covers 120/120, X+Z 90/120 -> X+Y."""
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=100.0),  # X
        _hsp((61, 120), (1184, 1363), bitscore=50.0),  # Y
        _hsp((61, 90), (1184, 1273), bitscore=50.0),  # Z (same start as Y, shorter)
    ]
    for hs in (hsps, hsps[::-1]):
        best, rest = _by_rank(_merge(hs))
        assert best.interval == Interval("chr1", 1000, 1363, "+")
        assert best.query_coverage == pytest.approx(1.0)
        assert rest.interval == Interval("chr1", 1183, 1273, "+")


def test_tie_on_bitscore_and_coverage_lower_start_wins() -> None:
    """W1+Y and W2+Y: same bitscore and coverage; W1 starts lower -> W1+Y, W2 alone."""
    hsps = [
        _hsp((61, 120), (1184, 1363), bitscore=100.0),  # Y
        _hsp((1, 60), (1001, 1180), bitscore=50.0),  # W1
        _hsp((1, 60), (1004, 1183), bitscore=50.0),  # W2, adjacent to Y
    ]
    for hs in (hsps, hsps[::-1]):
        best, rest = _by_rank(_merge(hs))
        assert best.interval == Interval("chr1", 1000, 1363, "+")
        assert rest.interval == Interval("chr1", 1003, 1183, "+")


def test_top_hsp_alone_ranks_below_better_chain() -> None:
    """Rank is by summed bitscore of the chain, not by its best HSP."""
    hsps = [
        _hsp((1, 120), (5001, 5360), bitscore=150.0),  # H alone
        _hsp((1, 60), (1001, 1180), bitscore=100.0),
        _hsp((61, 120), (1184, 1363), bitscore=100.0),
    ]
    r1, r2 = _by_rank(_merge(hsps))
    assert r1.interval == Interval("chr1", 1000, 1363, "+")
    assert r2.interval == Interval("chr1", 5000, 5360, "+")


# Hypothesis: on a linear contig with genome-disjoint HSPs (so "small genome overlap" rules do
# not matter) the rank-1 chain is the brute-force maximum-bitscore chain. Bitscores are distinct
# powers of two, so the maximum is unique and no tie rule is involved.


def _compatible(a: tuple[int, int, int, int], b: tuple[int, int, int, int]) -> bool:
    """``b`` may follow ``a``; tuples are (t0, t1, q0, q1), 0-based half-open, a.t1 <= b.t0."""
    return b[0] - a[1] <= 300 and b[2] > a[2] and b[3] > a[3] and a[3] - b[2] <= 10


@st.composite
def _disjoint_hsps(draw):  # type: ignore[no-untyped-def]
    n = draw(st.integers(1, 6))
    pos, out = 1000, []
    for _ in range(n):
        pos += draw(st.integers(0, 350))
        aa = draw(st.integers(5, 30))
        q0 = draw(st.integers(0, min(90, 160 - aa)))
        out.append((pos, pos + 3 * aa, q0, q0 + aa))
        pos += 3 * aa
    bits = draw(st.permutations([float(2**k) for k in range(n)]))
    return out, bits


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_disjoint_hsps(), st.randoms(use_true_random=False))
@example(  # A/C/B: greedy takes the nearest C; the best chain is A+B
    ([(1000, 1180, 0, 60), (1180, 1225, 60, 75), (1228, 1408, 60, 120)], [4.0, 1.0, 2.0]),
    random.Random(0),
)
@example(  # found by Hypothesis against the greedy chainer
    ([(1000, 1015, 0, 5), (1015, 1030, 1, 6), (1030, 1045, 1, 6)], [4.0, 1.0, 2.0]),
    random.Random(0),
)
def test_rank1_is_brute_force_best_chain(data, rnd) -> None:  # type: ignore[no-untyped-def]
    pieces, bits = data
    best_sum, best = -1.0, None
    for r in range(1, len(pieces) + 1):
        for combo in itertools.combinations(range(len(pieces)), r):
            if all(_compatible(pieces[i], pieces[j]) for i, j in itertools.pairwise(combo)):
                total = sum(bits[i] for i in combo)
                if total > best_sum:
                    best_sum, best = total, combo
    assert best is not None
    hsps = [
        _hsp((q0 + 1, q1), (t0 + 1, t1), bitscore=b, qlen=160)
        for (t0, t1, q0, q1), b in zip(pieces, bits, strict=True)
    ]
    rnd.shuffle(hsps)
    alns = _by_rank(_merge(hsps, LIN))
    r1 = alns[0]
    exp_iv = Interval("chr1", pieces[best[0]][0], pieces[best[-1]][1], "+")
    assert r1.interval == exp_iv
    assert (r1.q_start, r1.q_end) == (
        min(pieces[i][2] for i in best),
        max(pieces[i][3] for i in best),
    )
    for a in alns:
        check_interval(a.interval, LIN["chr1"])


# ---------------------------------------------------------------- 3. mutant killers

SMALL = {"c": Replicon("c", 500, "circular", "plasmid")}


def _small_hsps(a: float, b: float, c: float):  # type: ignore[no-untyped-def]
    """A [100,400) q0-100; B [450,500) q100-117; C [50,200) q117-167 after the origin.

    A+B+C (and A+C) would span unwrapped [100, 700) = 600 nt > L = 500: never one chain.
    """
    return [
        _hsp((1, 100), (101, 400), bitscore=a, contig="c", qlen=200),
        _hsp((101, 117), (451, 500), bitscore=b, contig="c", qlen=200),
        _hsp((118, 167), (51, 200), bitscore=c, contig="c", qlen=200),
    ]


def test_chain_longer_than_small_circular_contig_forward() -> None:
    """Best A+B (110), C alone: C may not be appended after the origin (span > L)."""
    alns = _by_rank(_merge(_small_hsps(100.0, 10.0, 50.0), SMALL))
    assert [a.interval for a in alns] == [
        Interval("c", 100, 500, "+"),
        Interval("c", 50, 200, "+"),
    ]
    for a in alns:
        check_interval(a.interval, SMALL["c"])
        assert a.interval.end - a.interval.start <= 500


def test_chain_longer_than_small_circular_contig_backward() -> None:
    """Best B+C (110) across the origin, A alone: A may not be prepended (span > L)."""
    alns = _by_rank(_merge(_small_hsps(50.0, 10.0, 100.0), SMALL))
    assert [a.interval for a in alns] == [
        Interval("c", 450, 700, "+"),
        Interval("c", 100, 400, "+"),
    ]
    for a in alns:
        check_interval(a.interval, SMALL["c"])
        assert a.interval.end - a.interval.start <= 500


def test_chain_exactly_contig_length_allowed() -> None:
    hsps = [
        _hsp((1, 100), (1, 300), contig="c", qlen=166),
        _hsp((101, 166), (303, 500), contig="c", qlen=166),
    ]
    (a,) = _merge(hsps, SMALL)
    assert a.interval == Interval("c", 0, 500, "+")


def test_same_genome_start_not_chained() -> None:
    """``b0 > a0`` is strict: B starting where A starts (but ending later) does not follow A."""
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=100.0),
        _hsp((51, 120), (1001, 1210), bitscore=90.0),
    ]
    alns = _merge(hsps)
    assert len(alns) == 2
    assert sorted(a.interval.end for a in alns) == [1180, 1210]


def test_earlier_genome_start_not_chained() -> None:
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=100.0),
        _hsp((51, 120), (998, 1210), bitscore=90.0),
    ]
    assert len(_merge(hsps)) == 2


@pytest.mark.parametrize(
    "s",
    [(49901, 50001), (50001, 49901), (0, 150), (150, 0), (50001, 50100)],
    ids=["plus-past-end", "minus-past-end", "plus-zero", "minus-zero", "beyond"],
)
def test_hit_outside_contig_raises(s: tuple[int, int]) -> None:
    with pytest.raises(ValueError):
        _merge([_hsp((1, 34), s, qlen=120)])


def test_hit_ending_exactly_at_contig_end_ok() -> None:
    (a,) = _merge([_hsp((1, 50), (49851, 50000), qlen=50)])
    assert a.interval == Interval("chr1", 49850, 50000, "+")


@pytest.mark.parametrize(
    ("q", "qlen"),
    [((0, 60), 120), ((1, 121), 120), ((61, 60), 120), ((-5, 60), 120)],
    ids=["qstart-0", "qend-past-qlen", "qstart-after-qend", "negative"],
)
def test_bad_query_range_raises(q: tuple[int, int], qlen: int) -> None:
    h = _tb().Hsp(
        qseqid="P1",
        sseqid="chr1",
        qstart=q[0],
        qend=q[1],
        sstart=1001,
        send=1180,
        length=60,
        nident=60,
        positive=60,
        bitscore=100.0,
        qlen=qlen,
    )
    with pytest.raises(ValueError):
        _merge([h])


def test_query_range_full_length_ok() -> None:
    (a,) = _merge([_hsp((1, 120), (1001, 1360), qlen=120)])
    assert (a.q_start, a.q_end) == (0, 120)


# ---------------------------------------------------------------- 4. exactly 15 columns


def _real_line() -> str:
    return (K12 / "tblastn.tsv").read_text().splitlines()[0]


@pytest.mark.parametrize("n", [14, 16, 20])
def test_parser_rejects_wrong_column_count(tmp_path: Path, n: int) -> None:
    cols = _real_line().split("\t")
    assert len(cols) == 15
    cols = cols[:n] if n < 15 else cols + ["extra"] * (n - 15)
    p = tmp_path / "bad.tsv"
    p.write_text("\t".join(cols) + "\n")
    with pytest.raises(ValueError, match="line 1"):
        _tb().parse_tblastn(p)


def test_parser_rejects_trailing_tab(tmp_path: Path) -> None:
    p = tmp_path / "bad.tsv"
    p.write_text(_real_line() + "\t\n")
    with pytest.raises(ValueError):
        _tb().parse_tblastn(p)


def test_parser_accepts_exactly_15_columns(tmp_path: Path) -> None:
    p = tmp_path / "ok.tsv"
    p.write_text(_real_line() + "\n")
    (h,) = _tb().parse_tblastn(p)
    assert h.sseqid == "NC_000913.3"
