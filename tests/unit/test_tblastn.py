"""tblastn route: commands, outfmt-6 parser, collinear HSP chaining (issue M3-13).

GIAP_SPEC.md §7.2.3 (tblastn leftovers), §7.1.1 (``-db_gencode 11``), §6.5 (coordinates).
Examples 1-7 of docs/specs/M3-13-tblastn-route.md on hand-written rows (in the column layout of
the real tool) and on ``tests/data/tblastn/mini/tblastn.tsv``, a real BLAST+ 2.16.0 run.

Ambiguities resolved here (spec silent; lead decisions of M3-12/M3-13 applied):
- ``max_gap_nt = 300`` and ``max_q_overlap_aa = 10`` are keyword defaults (Khoa, 2026-10-14);
  both limits are inclusive. The genome gap is the number of nt strictly between two HSPs
  (0-based half-open: ``next.start - prev.end``); a small genome overlap (negative gap) between
  query-collinear HSPs is allowed (real frameshift HSPs overlap, see P8 of the fixture).
- Across the origin the gap is measured on unwrapped coordinates (end-of-contig part plus
  start-of-contig part), with the same 300-nt limit.
- Chains are normalised to canonical intervals before ``assign_alignments`` (M3-12 review
  decision 4): ``0 <= start < L``, unwrapped ``end <= start + L``; linear contigs never wrap.
- Rank is per protein, by summed bitscore, 1 = best.
- ``query_coverage`` is the ``ProteinAlignment`` span ``(q_end - q_start) / qlen``; the spec's
  "|union of query ranges|" equals it for every chain tested here (no query gaps), so the
  union-vs-span question is left open.
- tblastn cannot see a GTG/TTG start (no genome sequence in ``merge_hsps``): ``alt_start`` is
  False; ``frameshift`` and ``stop_codon`` are 0 (spec).
- An HSP on a contig absent from ``replicons`` raises ``ValueError`` (contract drift, as in M3-12).
"""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import pytest

from giap.config.model import ProteinsConfig
from giap.model import Interval, Replicon
from giap.model.coords import check_interval
from giap.model.genes import Gene, GeneModel

OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
    "qlen nident positive"
)
MINI = Path(__file__).resolve().parents[1] / "data" / "tblastn" / "mini"
CFG = ProteinsConfig()

L = 50_000
CIRC = {"chr1": Replicon("chr1", L, "circular", "chromosome")}
LIN = {"chr1": Replicon("chr1", L, "linear", "unknown")}
MINI_REPS = {
    "chr1": Replicon("chr1", 8_000, "circular", "chromosome"),
    "pl1": Replicon("pl1", 3_000, "circular", "plasmid"),
}


# ---------------------------------------------------------------- helpers


def _tb():  # type: ignore[no-untyped-def]
    from giap.mapping import tblastn

    return tblastn


def _hsp(
    q: tuple[int, int],
    s: tuple[int, int],
    *,
    pid: str = "P1",
    contig: str = "chr1",
    qlen: int = 120,
    nident: int | None = None,
    positive: int | None = None,
    bitscore: float = 100.0,
    length: int | None = None,
):  # type: ignore[no-untyped-def]
    n = q[1] - q[0] + 1 if length is None else length
    return _tb().Hsp(
        qseqid=pid,
        sseqid=contig,
        qstart=q[0],
        qend=q[1],
        sstart=s[0],
        send=s[1],
        length=n,
        nident=n if nident is None else nident,
        positive=(n if nident is None else nident) if positive is None else positive,
        bitscore=bitscore,
        qlen=qlen,
    )


def _row(h) -> str:  # type: ignore[no-untyped-def]
    """An outfmt-6 line in the OUTFMT column order (layout copied from the real tool)."""
    pident = 100.0 * h.nident / h.length
    return "\t".join(
        [
            h.qseqid, h.sseqid, f"{pident:.3f}", str(h.length), str(h.length - h.nident), "0",
            str(h.qstart), str(h.qend), str(h.sstart), str(h.send), "1.00e-30",
            f"{h.bitscore:g}", str(h.qlen), str(h.nident), str(h.positive),
        ]
    )  # fmt: skip


def _merge(hsps, reps=CIRC, **kw):  # type: ignore[no-untyped-def]
    return _tb().merge_hsps(hsps, reps, **kw)


def _ivs(alns) -> list[Interval]:  # type: ignore[no-untyped-def]
    return sorted((a.interval for a in alns), key=lambda iv: (iv.contig, iv.start, iv.strand))


def _by_rank(alns, pid: str = "P1"):  # type: ignore[no-untyped-def]
    return sorted((a for a in alns if a.protein_id == pid), key=lambda a: a.rank)


def _gene(
    tag: str, start: int, end: int, strand: str = "+", contig: str = "chr1", ftype: str = "cds"
) -> Gene:
    n_aa = (end - start) // 3
    return Gene(
        locus_tag=tag,
        feature_type=ftype,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" + "A" * max(0, n_aa - 2),
        nt=None,
        start_codon="ATG",
        partial=False,
    )


def _model(reps: dict[str, Replicon], *genes: Gene) -> GeneModel:
    return GeneModel(reps, tuple(genes), ())


def _assign(alns, gm):  # type: ignore[no-untyped-def]
    from giap.mapping.protein_assign import assign_alignments

    return assign_alignments(alns, gm, CFG)


def _canonical(alns, reps) -> None:  # type: ignore[no-untyped-def]
    for a in alns:
        rep = reps[a.interval.contig]
        check_interval(a.interval, rep)
        assert 0 <= a.interval.start < rep.length


# ---------------------------------------------------------------- example 1: commands


def test_outfmt_constant() -> None:
    assert _tb().OUTFMT == OUTFMT


def test_example1_tblastn_cmd_exact_argv() -> None:
    cmd = _tb().tblastn_cmd(Path("db"), Path("p"), 8, "-db_gencode 11 -seg no -evalue 1e-5")
    assert cmd == [
        "tblastn", "-db", "db", "-query", "p", "-num_threads", "8",
        "-db_gencode", "11", "-seg", "no", "-evalue", "1e-5", "-outfmt", OUTFMT,
    ]  # fmt: skip


def test_tblastn_cmd_default_config_args() -> None:
    """``mapping.proteins.tblastn_args`` (§7.2.3) is passed verbatim, split on whitespace."""
    assert CFG.tblastn_args == "-db_gencode 11 -seg no -evalue 1e-5"
    cmd = _tb().tblastn_cmd(Path("genome.db"), Path("leftovers.faa"), 2, CFG.tblastn_args)
    assert cmd == [
        "tblastn", "-db", "genome.db", "-query", "leftovers.faa", "-num_threads", "2",
        "-db_gencode", "11", "-seg", "no", "-evalue", "1e-5", "-outfmt", OUTFMT,
    ]  # fmt: skip
    assert all(isinstance(a, str) for a in cmd)


def test_tblastn_cmd_other_args_reach_argv() -> None:
    cmd = _tb().tblastn_cmd(Path("d"), Path("q"), 1, "-db_gencode 4 -evalue 1e-3")
    assert cmd[-2:] == ["-outfmt", OUTFMT]
    assert cmd[7:-2] == ["-db_gencode", "4", "-evalue", "1e-3"]
    assert "-seg" not in cmd, "no hard-coded default may leak past the configured args"


def test_makeblastdb_cmd() -> None:
    cmd = _tb().makeblastdb_cmd(Path("genome.fa"), Path("db/genome"))
    assert cmd[0] == "makeblastdb"
    assert all(isinstance(a, str) for a in cmd)
    pairs = {cmd[i]: cmd[i + 1] for i in range(1, len(cmd) - 1) if cmd[i].startswith("-")}
    assert pairs.get("-in") == "genome.fa"
    assert pairs.get("-dbtype") == "nucl"
    assert pairs.get("-out") == "db/genome"
    assert "-parse_seqids" in cmd
    assert sorted(a for a in cmd if a.startswith("-")) == [
        "-dbtype",
        "-in",
        "-out",
        "-parse_seqids",
    ]


# ---------------------------------------------------------------- parser (real tool output)


def test_parse_real_rows_plus_and_minus() -> None:
    hsps = _tb().parse_tblastn(MINI / "tblastn.tsv")
    assert len(hsps) == 14
    assert hsps[0] == _tb().Hsp(
        qseqid="P1", sseqid="chr1", qstart=1, qend=120, sstart=1001, send=1360,
        length=120, nident=120, positive=120, bitscore=248.0, qlen=120,
    )  # fmt: skip
    p2 = hsps[1]
    assert (p2.qseqid, p2.sstart, p2.send) == ("P2", 2360, 2001)  # minus: sstart > send
    p8 = [h for h in hsps if h.qseqid == "P8"]
    assert [(h.qstart, h.qend, h.length, h.nident, h.positive) for h in p8] == [
        (1, 63, 63, 60, 61),
        (61, 120, 60, 60, 60),
    ]


def test_parse_types_and_fractional_bitscore() -> None:
    hsps = _tb().parse_tblastn(MINI / "tblastn.tsv")
    p5 = [h for h in hsps if h.qseqid == "P5"]
    assert [h.sseqid for h in p5] == ["pl1", "pl1"]
    assert p5[1].bitscore == pytest.approx(80.9)
    for h in hsps:
        for name in ("qstart", "qend", "sstart", "send", "length", "nident", "positive", "qlen"):
            assert type(getattr(h, name)) is int, name
        assert isinstance(h.bitscore, float)


def test_parse_keeps_pipe_ids_verbatim() -> None:
    """Edge case: UniProt-style IDs with ``|`` come back exactly as written to the query FASTA."""
    hsps = _tb().parse_tblastn(MINI / "tblastn.tsv")
    assert "sp|Q12345|X_ECOLI" in {h.qseqid for h in hsps}


def test_parse_file_order(tmp_path: Path) -> None:
    rows = [
        _hsp((1, 60), (1001, 1180), pid="B"),
        _hsp((1, 60), (2001, 2180), pid="A"),
    ]
    p = tmp_path / "t.tsv"
    p.write_text("".join(_row(h) + "\n" for h in rows))
    assert _tb().parse_tblastn(p) == rows


def test_parse_empty_file(tmp_path: Path) -> None:
    p = tmp_path / "empty.tsv"
    p.write_text("")
    assert _tb().parse_tblastn(p) == []


def test_parse_too_few_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "bad.tsv"
    p.write_text("P1\tchr1\t100.0\t120\t0\t0\t1\t120\t1001\t1360\n")
    with pytest.raises(ValueError):
        _tb().parse_tblastn(p)


def test_parse_non_integer_raises(tmp_path: Path) -> None:
    line = _row(_hsp((1, 120), (1001, 1360))).split("\t")
    line[6] = "x"
    p = tmp_path / "bad.tsv"
    p.write_text("\t".join(line) + "\n")
    with pytest.raises(ValueError):
        _tb().parse_tblastn(p)


# ---------------------------------------------------------------- examples 2-5


def test_example2_single_hsp() -> None:
    (a,) = _merge([_hsp((1, 120), (1001, 1360), nident=110, positive=115)])
    assert a.interval == Interval("chr1", 1000, 1360, "+")
    assert a.identity == pytest.approx(110 / 120)
    assert round(a.identity, 3) == 0.917
    assert a.query_coverage == pytest.approx(1.0)
    assert (a.q_start, a.q_end) == (0, 120)
    assert a.protein_id == "P1" and a.protein_len_aa == 120
    assert a.route == "tblastn"
    assert a.rank == 1
    assert a.frameshift == 0 and a.stop_codon == 0
    assert a.alt_start is False
    assert a.positive is None or a.positive == pytest.approx(115 / 120)


def test_example3_frameshift_hsps_merge() -> None:
    hsps = [
        _hsp((1, 60), (1001, 1180), nident=58, bitscore=110.0),
        _hsp((61, 120), (1184, 1363), nident=57, bitscore=105.0),
    ]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 1000, 1363, "+")
    assert a.query_coverage == pytest.approx(1.0)
    assert (a.q_start, a.q_end) == (0, 120)
    assert a.identity == pytest.approx((58 + 57) / 120)
    assert a.frameshift == 0 and a.stop_codon == 0
    assert a.rank == 1


def test_example4_far_hsps_two_chains_not_assigned() -> None:
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=110.0),
        _hsp((61, 120), (9001, 9180), bitscore=105.0),
    ]
    alns = _merge(hsps)
    assert _ivs(alns) == [Interval("chr1", 1000, 1180, "+"), Interval("chr1", 9000, 9180, "+")]
    assert all(a.query_coverage == pytest.approx(0.5) for a in alns)
    r1, r2 = _by_rank(alns)
    assert (r1.rank, r2.rank) == (1, 2)
    assert r1.interval.start == 1000  # higher bitscore
    gm = _model(CIRC, _gene("G1", 1000, 1363), _gene("G2", 9000, 9363))
    df, left = _assign(alns, gm)
    assert df.empty
    assert left == {"P1"}


def test_example5_minus_strand() -> None:
    (a,) = _merge([_hsp((1, 120), (2360, 2001))])
    assert a.interval == Interval("chr1", 2000, 2360, "-")


def test_minus_strand_chain_in_transcription_order() -> None:
    """On minus the protein N-terminus is at the high coordinate."""
    hsps = [_hsp((1, 60), (1363, 1184)), _hsp((61, 120), (1180, 1001))]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 1000, 1363, "-")
    assert a.query_coverage == pytest.approx(1.0)


def test_minus_strand_wrong_order_not_chained() -> None:
    hsps = [_hsp((1, 60), (1180, 1001)), _hsp((61, 120), (1363, 1184))]
    assert len(_merge(hsps)) == 2


def test_plus_strand_wrong_order_not_chained() -> None:
    hsps = [_hsp((61, 120), (1001, 1180)), _hsp((1, 60), (1184, 1363))]
    assert len(_merge(hsps)) == 2


def test_different_strands_not_chained() -> None:
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((61, 120), (1363, 1184))]
    alns = _merge(hsps)
    assert sorted(a.interval.strand for a in alns) == ["+", "-"]


def test_different_contigs_not_chained() -> None:
    reps = {**CIRC, "pl1": Replicon("pl1", 5_000, "circular", "plasmid")}
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((61, 120), (1184, 1363), contig="pl1")]
    alns = _merge(hsps, reps)
    assert sorted(a.interval.contig for a in alns) == ["chr1", "pl1"]


def test_different_proteins_not_chained() -> None:
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((61, 120), (1184, 1363), pid="P2")]
    alns = _merge(hsps)
    assert sorted(a.protein_id for a in alns) == ["P1", "P2"]
    assert all(a.rank == 1 for a in alns)


# ---------------------------------------------------------------- chaining limits


@pytest.mark.parametrize(("gap", "chained"), [(0, True), (300, True), (301, False)])
def test_max_gap_nt_inclusive(gap: int, chained: bool) -> None:
    start2 = 1180 + gap + 1  # 1-based sstart; 0-based gap = (start2 - 1) - 1180
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((61, 120), (start2, start2 + 179))]
    alns = _merge(hsps)
    assert len(alns) == (1 if chained else 2)
    if chained:
        assert alns[0].interval == Interval("chr1", 1000, start2 + 179, "+")


@pytest.mark.parametrize(("overlap", "chained"), [(0, True), (10, True), (11, False)])
def test_max_q_overlap_aa_inclusive(overlap: int, chained: bool) -> None:
    q2 = 61 - overlap
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((q2, 120), (1184, 1184 + 3 * (120 - q2) + 2))]
    alns = _merge(hsps)
    assert len(alns) == (1 if chained else 2)
    if chained:
        assert alns[0].query_coverage == pytest.approx(1.0)
        assert (alns[0].q_start, alns[0].q_end) == (0, 120)
        assert alns[0].query_coverage <= 1.0


def test_limits_are_keyword_parameters() -> None:
    hsps = [_hsp((1, 60), (1001, 1180)), _hsp((61, 120), (1301, 1480))]  # gap 120 nt
    assert len(_merge(hsps)) == 1
    assert len(_merge(hsps, max_gap_nt=100)) == 2
    q = [_hsp((1, 60), (1001, 1180)), _hsp((56, 120), (1184, 1378))]  # overlap 5 aa
    assert len(_merge(q)) == 1
    assert len(_merge(q, max_q_overlap_aa=4)) == 2


def test_contained_query_range_not_chained() -> None:
    """The next HSP must start after the previous one on the protein (not inside it)."""
    hsps = [_hsp((1, 120), (1001, 1360), bitscore=240.0), _hsp((30, 50), (1401, 1463))]
    alns = _merge(hsps)
    assert len(alns) == 2
    assert _by_rank(alns)[0].interval == Interval("chr1", 1000, 1360, "+")


def test_three_hsps_chain() -> None:
    hsps = [
        _hsp((81, 120), (1271, 1390), bitscore=70.0),
        _hsp((1, 40), (1001, 1120), bitscore=70.0),
        _hsp((41, 80), (1131, 1250), bitscore=90.0),
    ]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 1000, 1390, "+")
    assert a.query_coverage == pytest.approx(1.0)


# ---------------------------------------------------------------- examples 6-7: origin


def test_example6_origin_crossing_chain_circular() -> None:
    hsps = [_hsp((1, 50), (49851, 50000), qlen=100), _hsp((51, 100), (1, 150), qlen=100)]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 49850, 50150, "+")
    assert a.interval.crosses_origin(CIRC["chr1"])
    assert a.query_coverage == pytest.approx(1.0)
    _canonical([a], CIRC)


def test_example6_input_order_does_not_change_canonical_interval() -> None:
    hsps = [_hsp((51, 100), (1, 150), qlen=100), _hsp((1, 50), (49851, 50000), qlen=100)]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 49850, 50150, "+")


def test_example7_same_on_linear_contig_two_chains() -> None:
    hsps = [_hsp((1, 50), (49851, 50000), qlen=100), _hsp((51, 100), (1, 150), qlen=100)]
    alns = _merge(hsps, LIN)
    assert _ivs(alns) == [Interval("chr1", 0, 150, "+"), Interval("chr1", 49850, 50000, "+")]
    _canonical(alns, LIN)


def test_origin_crossing_minus_strand() -> None:
    """Minus strand: N-terminus near 0 (high unwrapped end), C-terminus near the contig end."""
    hsps = [_hsp((1, 50), (150, 1), qlen=100), _hsp((51, 100), (50000, 49851), qlen=100)]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 49850, 50150, "-")
    _canonical([a], CIRC)
    assert len(_merge(hsps, LIN)) == 2


def test_origin_gap_measured_unwrapped() -> None:
    """100 nt before the end plus 100 nt after 0 = 200 nt <= 300: chained."""
    hsps = [_hsp((1, 50), (49751, 49900), qlen=100), _hsp((51, 100), (101, 250), qlen=100)]
    (a,) = _merge(hsps)
    assert a.interval == Interval("chr1", 49750, 50250, "+")
    _canonical([a], CIRC)


def test_origin_gap_over_limit_not_chained() -> None:
    """200 nt before the end plus 200 nt after 0 = 400 nt > 300: two chains."""
    hsps = [_hsp((1, 50), (49651, 49800), qlen=100), _hsp((51, 100), (201, 350), qlen=100)]
    alns = _merge(hsps)
    assert len(alns) == 2
    _canonical(alns, CIRC)


def test_far_from_origin_not_wrapped() -> None:
    """An HSP near 0 and one in the middle of the contig never chain across the origin."""
    hsps = [_hsp((1, 50), (24851, 25000), qlen=100), _hsp((51, 100), (1, 150), qlen=100)]
    assert len(_merge(hsps)) == 2


# ---------------------------------------------------------------- ranks, statistics


def test_rank_by_summed_bitscore() -> None:
    hsps = [
        _hsp((1, 60), (1001, 1180), bitscore=60.0),
        _hsp((61, 120), (1184, 1363), bitscore=60.0),  # chain total 120
        _hsp((1, 120), (20001, 20360), bitscore=100.0),  # single HSP 100
    ]
    r1, r2 = _by_rank(_merge(hsps))
    assert (r1.rank, r2.rank) == (1, 2)
    assert r1.interval == Interval("chr1", 1000, 1363, "+")
    assert r2.interval == Interval("chr1", 20000, 20360, "+")


def test_multi_locus_ranks_1_and_2() -> None:
    """Edge case: two chains on different CDSs of one protein -> both rows kept."""
    hsps = [_hsp((1, 80), (4001, 4240), qlen=80, bitscore=166.0)]
    hsps.append(_hsp((1, 80), (5001, 5240), qlen=80, bitscore=169.0))
    alns = _merge(hsps)
    r1, r2 = _by_rank(alns)
    assert r1.interval == Interval("chr1", 5000, 5240, "+")
    assert r2.interval == Interval("chr1", 4000, 4240, "+")
    gm = _model(CIRC, _gene("G4", 4000, 4240), _gene("G5", 5000, 5240))
    df, left = _assign(alns, gm)
    assert sorted(zip(df["rank"], df["gene"], strict=True)) == [(1, "G5"), (2, "G4")]
    assert left == set()


def test_merged_identity_is_sum_nident_over_sum_length() -> None:
    hsps = [
        _hsp((1, 60), (1001, 1180), nident=50, length=62),
        _hsp((61, 120), (1184, 1363), nident=55, length=60),
    ]
    (a,) = _merge(hsps)
    assert a.identity == pytest.approx(105 / 122)


def test_empty_input() -> None:
    assert _merge([]) == []


def test_unknown_contig_raises() -> None:
    with pytest.raises(ValueError):
        _merge([_hsp((1, 120), (1001, 1360), contig="nope")])


# ---------------------------------------------------------------- real run: merge + assign


def _mini_alns(reps=MINI_REPS):  # type: ignore[no-untyped-def]
    return _merge(_tb().parse_tblastn(MINI / "tblastn.tsv"), reps)


def test_real_run_merged_intervals() -> None:
    alns = _mini_alns()
    _canonical(alns, MINI_REPS)
    got = {(a.protein_id, a.rank): a.interval for a in alns}
    assert got == {
        ("P1", 1): Interval("chr1", 1000, 1360, "+"),
        ("P2", 1): Interval("chr1", 2000, 2360, "-"),
        ("P3", 1): Interval("chr1", 3000, 3367, "+"),  # engineered frameshift, 7-nt gap
        ("P4", 1): Interval("chr1", 7850, 8150, "+"),  # across the origin of chr1
        ("P5", 1): Interval("pl1", 2880, 3180, "-"),  # across the origin of pl1, minus
        ("P6", 1): Interval("chr1", 5000, 5240, "+"),  # bitscore 169
        ("P6", 2): Interval("chr1", 4000, 4240, "+"),  # bitscore 166
        ("P7", 1): Interval("chr1", 5500, 5950, "+"),
        ("P8", 1): Interval("chr1", 6500, 6861, "+"),  # natural frameshift, HSPs overlap 3 aa
        ("sp|Q12345|X_ECOLI", 1): Interval("chr1", 1000, 1150, "+"),
    }
    for a in alns:
        assert a.route == "tblastn"
        assert a.query_coverage == pytest.approx(1.0)
    p8 = next(a for a in alns if a.protein_id == "P8")
    assert p8.identity == pytest.approx(120 / 123)
    p7 = next(a for a in alns if a.protein_id == "P7")
    assert p7.identity == pytest.approx(131 / 150)


def test_real_run_linear_plasmid_splits_origin_chain() -> None:
    reps = {**MINI_REPS, "pl1": Replicon("pl1", 3_000, "linear", "plasmid")}
    p5 = [a for a in _mini_alns(reps) if a.protein_id == "P5"]
    assert _ivs(p5) == [Interval("pl1", 0, 180, "-"), Interval("pl1", 2880, 3000, "-")]
    assert sorted(a.rank for a in p5) == [1, 2]


def _mini_genes(p1_type: str = "cds") -> GeneModel:
    return _model(
        MINI_REPS,
        _gene("G1", 1000, 1360, ftype=p1_type),
        _gene("G2", 2000, 2360, "-"),
        _gene("G3", 3000, 3367),
        _gene("G4", 7850, 8150),
        _gene("G5", 2880, 3180, "-", contig="pl1"),
        _gene("G6a", 4000, 4240),
        _gene("G6b", 5000, 5240),
        _gene("G7", 5500, 5950),
        _gene("G8", 6500, 6861, ftype="sorf"),
    )


def test_real_run_assignment_via_m3_12() -> None:
    df, left = _assign(_mini_alns(), _mini_genes())
    cols = (df["protein_id"], df["rank"], df["gene"], df["class"])
    got = {(p, r): (None if pd.isna(g) else g, c) for p, r, g, c in zip(*cols, strict=True)}
    assert got[("P1", 1)] == ("G1", "identical")
    assert got[("P2", 1)] == ("G2", "identical")
    assert got[("P3", 1)][0] == "G3"
    assert got[("P4", 1)] == ("G4", "identical")
    assert got[("P5", 1)] == ("G5", "identical")
    assert got[("P6", 1)][0] == "G6b" and got[("P6", 2)][0] == "G6a"
    assert got[("P7", 1)] == ("G7", "divergent")
    assert got[("P8", 1)] == ("G8", "near")  # sORF target; identity 120/123
    # 50-aa UniProt fragment: identity/coverage pass, reciprocal overlap with G1 150/360 fails
    assert got[("sp|Q12345|X_ECOLI", 1)][0] is None
    assert left == {"sp|Q12345|X_ECOLI"}
    assert set(df["route"]) == {"tblastn"}
    assert set(df["flags"].dropna()) == set()  # frameshift/stop_codon are 0 for tblastn


def test_pseudogene_is_not_a_target() -> None:
    """M3-12 lead decision 1: CDS/sORF only; a chain over a pseudogene stays a leftover."""
    df, left = _assign(_mini_alns(), _mini_genes(p1_type="pseudogene"))
    assert "P1" in left
    p1 = df[df["protein_id"] == "P1"]
    assert len(p1) == 1 and p1["gene"].isna().all()


def test_origin_gene_assignment_on_circular_contig() -> None:
    hsps = [_hsp((1, 50), (49851, 50000), qlen=100), _hsp((51, 100), (1, 150), qlen=100)]
    gm = _model(CIRC, _gene("ORI", 49850, 50150))
    df, left = _assign(_merge(hsps), gm)
    assert list(df["gene"]) == ["ORI"]
    assert (int(df["start"].iloc[0]), int(df["end"].iloc[0])) == (49850, 50150)
    assert left == set()


# ---------------------------------------------------------------- determinism


def _key(a):  # type: ignore[no-untyped-def]
    return (a.protein_id, a.rank, a.interval.contig, a.interval.start, a.interval.end)


def test_input_order_invariance_real_run() -> None:
    hsps = _tb().parse_tblastn(MINI / "tblastn.tsv")
    ref = sorted(_merge(hsps, MINI_REPS), key=_key)
    rng = random.Random(0)
    for _ in range(20):
        shuffled = list(hsps)
        rng.shuffle(shuffled)
        assert sorted(_merge(shuffled, MINI_REPS), key=_key) == ref


def test_repeatable() -> None:
    hsps = _tb().parse_tblastn(MINI / "tblastn.tsv")
    assert _merge(hsps, MINI_REPS) == _merge(hsps, MINI_REPS)


def test_merge_accepts_tuple_and_does_not_mutate() -> None:
    hsps = tuple(_tb().parse_tblastn(MINI / "tblastn.tsv"))
    before = list(hsps)
    _merge(hsps, MINI_REPS)
    assert list(hsps) == before
