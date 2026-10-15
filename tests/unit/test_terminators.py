"""Terminators: TerminatorNet, TransTermHP, cross-check, BacTermFinder, ``terminator_between``
(M2-15; GIAP_SPEC.md §7.1.1 TerminatorNet/TransTermHP/BacTermFinder rows and gotchas, §7.1.5, G12).

Examples 1-10 and the edge cases of docs/specs/M2-15-terminators.md, plus the real formats of the
recordings copied to ``tests/data/recorded/terminators/`` (see PROVENANCE.md there).

Conventions pinned here where the spec is silent or its notation differs from the real tool:

* TerminatorNet's real output (commit 7218030) is TAB-separated without a header (the spec's
  ``name, TERM_i, …`` lists the columns); Example 2 is therefore fed as a TAB-separated line. A
  contig/strand without calls is the 2-column line ``<contig>\\tNO TERMINATORS FOUND``.
* Parsers validate every interval against ``replicons`` (unknown contig or out of range:
  ValueError).
* ``.tt`` gene lines (column 1) and sequence lines (indented, no ``TERM``) are ignored; the location
  code may have two letters (``GF``, ``RH``, ``Tf``) and the hairpin score may be an integer.
* ``terminators_bed`` keeps origin-crossing intervals unwrapped (``write_bed`` splits them);
  names are ``<source>_<n>`` numbered from 1 per source.
* Between/downstream windows are evaluated in unwrapped coordinates on circular replicons, so a
  window crossing the origin matches calls stored on either side of it.
* ``CrossRow`` exposes the Q2 columns as attributes ``contig``, ``n_terminatornet``,
  ``n_transtermhp``, ``n_matched``, ``share_matched``; the ``all`` row is last.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

import giapsim
from giap.annotate.terminatornet_cmd import terminatornet_command
from giap.annotate.terminators import (
    Terminator,
    crosscheck,
    parse_terminatornet,
    parse_transterm_tt,
    terminator_between,
    terminator_conf_downstream,
    terminators_bed,
    write_transterm_coords,
)
from giap.annotate.transtermhp_cmd import transterm_command
from giap.io.bed import BedRecord, bed_lines
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

DATA = Path(__file__).resolve().parents[1] / "data" / "recorded" / "terminators"
TN_SIM = DATA / "terminatornet" / "sim100k" / "terminatornet.tsv"
TT_SIM = DATA / "transtermhp" / "sim100k" / "transtermhp.tt"

CHR1 = Replicon("chr1", 100_000, "circular", "chromosome")
SIM_REPS = {"chr1": CHR1, "plasmid_1": Replicon("plasmid_1", 8_000, "circular", "plasmid")}
K12 = {"U00096.3": Replicon("U00096.3", 4_641_652, "circular", "chromosome")}
NISSLE_PL = {
    "CP171243.1": Replicon("CP171243.1", 4_523, "circular", "plasmid"),
    "CP171244.1": Replicon("CP171244.1", 3_173, "circular", "plasmid"),
    "CP171245.1": Replicon("CP171245.1", 5_514, "circular", "plasmid"),
}
COMP = str.maketrans("ACGTN", "TGCAN")


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


def term(
    start: int,
    end: int,
    strand: str,
    conf: float,
    source: str = "TerminatorNet",
    contig: str = "chr1",
) -> Terminator:
    return Terminator(Interval(contig, start, end, strand), conf, source)  # type: ignore[arg-type]


def tn_line(
    contig: str, i: int, start1: int, stop1: int, strand: str, pred: float, gh: float
) -> str:
    return "\t".join(
        [contig, f"TERM_{i}", "AAAAAAAA", "GCCGCCAAAAGGCGGC", "TTTTTTTT", "((((((....))))))",
         repr(gh), repr(pred), str(start1), str(stop1), strand]
    )  # fmt: skip


def _raw_tn(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln.count("\t") >= 10]


def _raw_tt(path: Path) -> list[tuple[str, int, int, str, int]]:
    """(contig, a, b, strand, conf) of every TERM line, independent of the parser under test."""
    out, contig = [], ""
    for ln in path.read_text().splitlines():
        if ln.startswith("SEQUENCE "):
            contig = ln.split()[1]
        elif ln.startswith("  TERM "):
            f = ln.split()
            out.append((contig, int(f[2]), int(f[4]), f[5], int(f[7])))
    return out


# ------------------------------------------------------------------- Example 1: TerminatorNet argv


def test_example1_terminatornet_command_exact_flags() -> None:
    home = Path("/opt/terminatornet")
    argv = terminatornet_command(Path("genome.fa"), Path("terms.tsv"), 8, home)
    assert argv == [
        "python", f"{home}/TerminatorNet.py", "-f", "genome.fa", "-s", "-n", "8", "-o", "terms.tsv",
    ]  # fmt: skip
    assert all(isinstance(a, str) for a in argv)


def test_terminatornet_command_always_passes_s_and_only_documented_flags() -> None:
    """G12 / gotcha: the parser matches flags by substring; pass only -f -s -n -o, always -s."""
    argv = terminatornet_command(Path("genome.fa"), Path("terms.tsv"), 1, Path("/h"))
    flags = [a for a in argv if a.startswith("-")]
    assert flags == ["-f", "-s", "-n", "-o"]


def test_terminatornet_command_matches_the_recorded_invocation() -> None:
    rec = (DATA / "terminatornet" / "sim100k" / "command.txt").read_text().splitlines()[1].split()
    home = Path(rec[1]).parent
    assert terminatornet_command(Path("genome.fa"), Path("terminatornet.tsv"), 17, home) == rec


@pytest.mark.parametrize("bad", ["my genome.fa", "dir with space/genome.fa"])
def test_terminatornet_command_rejects_names_with_spaces(bad: str) -> None:
    """Gotcha: staged file names without spaces (the tool's argument parser splits them)."""
    with pytest.raises(ValueError):
        terminatornet_command(Path(bad), Path("terms.tsv"), 2, Path("/h"))


def test_transterm_command_matches_the_recorded_invocation() -> None:
    rec = (DATA / "transtermhp" / "sim100k" / "command.txt").read_text().splitlines()[1].split()
    assert rec[:2] == ["transterm", "-p"]
    argv = transterm_command(Path("genome.fa"), Path("genes.coords"), Path(rec[2]))
    assert argv == ["transterm", "-p", rec[2], "genome.fa", "genes.coords"]
    assert argv == rec


# ------------------------------------------------------------------- Examples 2-3: TerminatorNet


def test_example2_minus_strand_line(tmp_path: Path) -> None:
    p = tmp_path / "terms.tsv"
    p.write_text(
        "chr1\tTERM_3\tAAAAAA\tGCCGCCAAAAGGCGGC\tTTTTTTT\t((((((....))))))\t-12.4\t0.87\t1201\t1240\t-\n"
    )
    got = parse_terminatornet(p, {"chr1": CHR1})
    assert got == [Terminator(Interval("chr1", 1200, 1240, "-"), 0.87, "TerminatorNet", -12.4)]


def test_example3_only_no_terminators_found(tmp_path: Path) -> None:
    p = tmp_path / "terms.tsv"
    p.write_text("NO TERMINATORS FOUND\n")
    assert parse_terminatornet(p, {"chr1": CHR1}) == []


def test_real_no_terminators_found_line_is_skipped(tmp_path: Path) -> None:
    """Real form (sim100k recording): ``plasmid_1<TAB>NO TERMINATORS FOUND``; also blank lines."""
    p = tmp_path / "terms.tsv"
    p.write_text(
        "\n" + tn_line("chr1", 1, 1888, 1907, "+", 0.5, -16.0) + "\n\n"
        "plasmid_1\tNO TERMINATORS FOUND\n\n"
    )
    got = parse_terminatornet(p, SIM_REPS)
    assert got == [Terminator(Interval("chr1", 1887, 1907, "+"), 0.5, "TerminatorNet", -16.0)]


def test_missing_strand_column_raises(tmp_path: Path) -> None:
    """Without ``-s`` TerminatorNet writes 10 columns; we always pass ``-s``: an error."""
    p = tmp_path / "terms.tsv"
    p.write_text(tn_line("chr1", 1, 100, 140, "+", 0.9, -10.0).rsplit("\t", 1)[0] + "\n")
    with pytest.raises(ValueError):
        parse_terminatornet(p, {"chr1": CHR1})


def test_unknown_contig_or_out_of_range_raises(tmp_path: Path) -> None:
    p = tmp_path / "a.tsv"
    p.write_text(tn_line("chrX", 1, 100, 140, "+", 0.9, -10.0) + "\n")
    with pytest.raises(ValueError):
        parse_terminatornet(p, {"chr1": CHR1})
    lin = {"c": Replicon("c", 1_000, "linear")}
    p.write_text(tn_line("c", 1, 990, 1010, "+", 0.9, -10.0) + "\n")
    with pytest.raises(ValueError):
        parse_terminatornet(p, lin)


def test_contig_with_safe_id_characters_survives(tmp_path: Path) -> None:
    """Edge case: a contig ``chr_1.a-b`` through a recorded-format line."""
    reps = {"chr_1.a-b": Replicon("chr_1.a-b", 50_000, "circular")}
    p = tmp_path / "t.tsv"
    p.write_text(tn_line("chr_1.a-b", 7, 30001, 30038, "-", 0.61, -14.5) + "\n")
    got = parse_terminatornet(p, reps)
    assert got == [
        Terminator(Interval("chr_1.a-b", 30000, 30038, "-"), 0.61, "TerminatorNet", -14.5)
    ]


def test_recorded_sim100k_terminatornet_parses_with_forward_coordinates() -> None:
    raw = _raw_tn(TN_SIM)
    got = parse_terminatornet(TN_SIM, SIM_REPS)
    assert len(got) == len(raw) == 133
    assert sum(t.interval.strand == "+" for t in got) == 60
    assert sum(t.interval.strand == "-" for t in got) == 73
    assert {t.interval.contig for t in got} == {"chr1", "plasmid_1"}
    assert all(t.source == "TerminatorNet" for t in got)
    first = Terminator(
        Interval("chr1", 1887, 1907, "+"), 0.8812014060630281, "TerminatorNet", -16.0
    )
    assert first in got
    expected = sorted(
        (f[0], int(f[8]) - 1, int(f[9]), f[10], float(f[7]), float(f[6])) for f in raw
    )
    seen = sorted(
        (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand, t.confidence,
         t.hairpin_dg)
        for t in got
    )  # fmt: skip
    assert seen == expected


@pytest.fixture(scope="module")
def sim_genome(tmp_path_factory: pytest.TempPathFactory) -> dict[str, str]:
    """The recorded sim100k runs used giapsim SIM100K's genome (same sequences, checked)."""
    res = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("s"))
    return {r.id: r.sequence.upper() for r in read_genome_fasta(res.genome_fasta)}


def test_recorded_hairpins_equal_the_genome_slice_on_both_strands(
    sim_genome: dict[str, str],
) -> None:
    """Module 1 test "a TerminatorNet -s and coordinate test": interval = hairpin on the genome."""
    raw = {(f[0], int(f[8]), int(f[9]), f[10]): f[3].upper() for f in _raw_tn(TN_SIM)}
    for t in parse_terminatornet(TN_SIM, SIM_REPS):
        iv = t.interval
        seq = sim_genome[iv.contig][iv.start : iv.end]
        if iv.strand == "-":
            seq = seq.translate(COMP)[::-1]
        assert seq == raw[(iv.contig, iv.start + 1, iv.end, iv.strand)], iv


def test_recorded_k12_excerpt_near_origin_and_end() -> None:
    got = parse_terminatornet(DATA / "terminatornet" / "k12_excerpt.tsv", K12)
    assert len(got) == 10
    assert got[0] == Terminator(
        Interval("U00096.3", 271, 309, "+"), 1.0, "TerminatorNet", -24.199999999999996
    )
    minus = [t for t in got if t.interval.strand == "-"]
    assert Interval("U00096.3", 273, 307, "-") in {t.interval for t in minus}
    assert max(t.interval.end for t in got) <= 4_641_652


def test_recorded_nissle_plasmids_parse_per_contig() -> None:
    got = parse_terminatornet(DATA / "terminatornet" / "nissle_plasmids.tsv", NISSLE_PL)
    by = {c: sum(t.interval.contig == c for t in got) for c in NISSLE_PL}
    assert by == {"CP171243.1": 9, "CP171244.1": 7, "CP171245.1": 10}
    assert Interval("CP171245.1", 3568, 3588, "-") in {t.interval for t in got}


def test_parse_is_deterministic() -> None:
    assert parse_terminatornet(TN_SIM, SIM_REPS) == parse_terminatornet(TN_SIM, SIM_REPS)
    assert parse_transterm_tt(TT_SIM, SIM_REPS) == parse_transterm_tt(TT_SIM, SIM_REPS)


# ------------------------------------------------------------------- Example 4: .coords writer


def test_example4_coords_lines(tmp_path: Path) -> None:
    p = tmp_path / "genes.coords"
    write_transterm_coords(
        p, [gene("GSIM_00005", 99, 400, "+"), gene("GSIM_00010", 1000, 1300, "-")]
    )
    assert p.read_text() == "GSIM_00005 100 400 chr1\nGSIM_00010 1300 1001 chr1\n"


def test_coords_include_cds_sorf_and_pseudogene(tmp_path: Path) -> None:
    p = tmp_path / "genes.coords"
    genes = [
        gene("A_1", 0, 300, "+"),
        gene("A_2", 400, 460, "-", ft="sorf"),
        gene("A_3", 500, 900, "+", contig="p2", ft="pseudogene"),
    ]
    write_transterm_coords(p, genes)
    assert p.read_text().splitlines() == ["A_1 1 300 chr1", "A_2 460 401 chr1", "A_3 501 900 p2"]


def _recorded_genes(path: Path) -> list[Gene]:
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == ["locus_tag", "feature_type", "contig", "start", "end", "strand"]
    out = []
    for ln in lines[1:]:
        tag, ft, contig, s, e, strand = ln.split("\t")
        out.append(gene(tag, int(s), int(e), strand, contig, ft))
    return out


def test_coords_writer_reproduces_the_recorded_transtermhp_input(tmp_path: Path) -> None:
    """The sim100k ``genes.coords`` TransTermHP ran on, rebuilt from the same Bakta genes."""
    genes = _recorded_genes(DATA / "transtermhp" / "sim100k" / "genes.tsv")
    p = tmp_path / "genes.coords"
    write_transterm_coords(p, genes)
    assert p.read_bytes() == (DATA / "transtermhp" / "sim100k" / "genes.coords").read_bytes()


# ------------------------------------------------------------------- Example 5: .tt parser


def test_example5_tt_line_under_sequence(tmp_path: Path) -> None:
    p = tmp_path / "t.tt"
    p.write_text(
        "TransTermHP v2.09\n\nSEQUENCE chr1  (length 100000)\n\n"
        "GSIM_0001        100 - 4300     + | \n"
        "  TERM 1   4342 - 4366   +   F   93 -11.5 -3.8 | \n"
        "  AAAAAAAAAAAAAAA      GCCGCC TTTT GGCGGC     TTTTTTTTTTTTTTT\n"
    )
    got = parse_transterm_tt(p, {"chr1": CHR1})
    assert len(got) == 1
    assert got[0].interval == Interval("chr1", 4341, 4366, "+")
    assert got[0].confidence == pytest.approx(0.93)
    assert got[0].source == "TransTermHP"


def test_recorded_sim100k_tt_parses() -> None:
    raw = _raw_tt(TT_SIM)
    got = parse_transterm_tt(TT_SIM, SIM_REPS)
    assert len(got) == len(raw) == 30
    assert all(t.source == "TransTermHP" and t.interval.contig == "chr1" for t in got)
    assert got[0].interval == Interval("chr1", 3909, 3927, "-")  # "TERM 1  3927 - 3910  -"
    assert got[0].confidence == pytest.approx(1.0)
    expected = sorted((c, min(a, b) - 1, max(a, b), s, conf / 100) for c, a, b, s, conf in raw)
    seen = sorted(
        (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand, t.confidence)
        for t in got
    )
    assert [r[:4] for r in seen] == [r[:4] for r in expected]
    assert [r[4] for r in seen] == pytest.approx([r[4] for r in expected])


def test_recorded_k12_tt_excerpt_loc_codes_and_integer_hp() -> None:
    got = parse_transterm_tt(DATA / "transtermhp" / "k12_excerpt.tt", K12)
    ivs = [(t.interval.start, t.interval.end, t.interval.strand) for t in got]
    assert len(got) == 9
    assert (279, 301, "+") in ivs  # TERM 1, hp -13 (integer)
    assert (270621, 270639, "+") in ivs  # TERM 93, loc GF, hp -7
    assert (363202, 363223, "-") in ivs  # TERM 120, loc RH
    assert (2347176, 2347220, "-") in ivs  # TERM 888, loc Tf
    assert (4639546, 4639564, "-") in ivs  # TERM 1847
    conf = {(t.interval.start, t.interval.strand): t.confidence for t in got}
    assert conf[(5043, "+")] == pytest.approx(0.91)
    assert conf[(270621, "+")] == pytest.approx(0.81)


def test_recorded_nissle_tt_plasmid_sections() -> None:
    reps = {"CP171242.1": Replicon("CP171242.1", 5_051_730, "circular"), **NISSLE_PL}
    got = parse_transterm_tt(DATA / "transtermhp" / "nissle_plasmids.tt", reps)
    by = {c: [t for t in got if t.interval.contig == c] for c in NISSLE_PL}
    assert [len(v) for v in by.values()] == [5, 0, 1]
    assert by["CP171245.1"][0].interval == Interval("CP171245.1", 3569, 3587, "-")
    assert Interval("CP171243.1", 1801, 1830, "-") in {t.interval for t in by["CP171243.1"]}


# ------------------------------------------------------------------- Example 6: cross-check


def _rows(rows: list) -> dict[str, tuple[int, int, int, float]]:  # type: ignore[type-arg]
    assert rows and rows[-1].contig == "all"
    return {
        r.contig: (r.n_terminatornet, r.n_transtermhp, r.n_matched, r.share_matched) for r in rows
    }


@pytest.mark.parametrize(("other_start", "matched"), [(1060, 1), (1061, 0)])
def test_example6_tolerance_boundary(other_start: int, matched: int) -> None:
    prim = [term(1000, 1040, "+", 0.9)]
    oth = [term(other_start, 1090, "+", 0.9, "TransTermHP")]
    rows = _rows(crosscheck(prim, oth))
    assert rows["chr1"][:3] == (1, 1, matched)
    assert rows["chr1"][3] == pytest.approx(float(matched))
    assert rows["all"][:3] == (1, 1, matched)


def test_crosscheck_overlap_strand_contig_and_upstream_side() -> None:
    prim = [
        term(1000, 1040, "+", 0.9),  # overlap → gap 0 → matched
        term(2000, 2040, "+", 0.9),  # other on the minus strand only → not matched
        term(3000, 3040, "-", 0.9),  # other ends 20 bp before → matched
        term(100, 140, "+", 0.9, contig="p2"),  # other on chr1 only → not matched
    ]
    oth = [
        term(1020, 1060, "+", 1.0, "TransTermHP"),
        term(2000, 2040, "-", 1.0, "TransTermHP"),
        term(2950, 2980, "-", 1.0, "TransTermHP"),
        term(100, 140, "+", 1.0, "TransTermHP"),
    ]
    rows = _rows(crosscheck(prim, oth))
    assert rows["chr1"][:3] == (3, 4, 2)
    assert rows["chr1"][3] == pytest.approx(2 / 3)
    assert rows["p2"][:3] == (1, 0, 0)
    assert rows["all"][:3] == (4, 4, 2)
    assert rows["all"][3] == pytest.approx(0.5)


def test_crosscheck_tol_parameter() -> None:
    prim, oth = [term(1000, 1040, "+", 0.9)], [term(1070, 1090, "+", 0.9, "TransTermHP")]
    assert _rows(crosscheck(prim, oth))["all"][2] == 0
    assert _rows(crosscheck(prim, oth, tol_bp=30))["all"][2] == 1


def test_crosscheck_on_the_recorded_sim100k_runs() -> None:
    prim = parse_terminatornet(TN_SIM, SIM_REPS)
    oth = parse_transterm_tt(TT_SIM, SIM_REPS)

    def gap(a: Interval, b: Interval) -> int:
        return max(0, max(a.start, b.start) - min(a.end, b.end))

    want = sum(
        any(
            o.interval.contig == p.interval.contig
            and o.interval.strand == p.interval.strand
            and gap(o.interval, p.interval) <= 20
            for o in oth
        )
        for p in prim
        if p.interval.contig == "chr1"
    )
    rows = _rows(crosscheck(prim, oth))
    assert rows["chr1"][:3] == (130, 30, want)
    assert rows["chr1"][3] == pytest.approx(want / 130)
    assert rows["plasmid_1"][:3] == (3, 0, 0)
    assert rows["all"][:3] == (133, 30, want)
    assert 0 < want < 130


# ------------------------------------------------------------------- Examples 7-8, 10: between


UP_P, DOWN_P = gene("U", 100, 400, "+"), gene("D", 420, 700, "+")


@pytest.mark.parametrize(
    ("t", "want"),
    [
        (term(360, 400, "+", 0.8), 0.8),
        (term(345, 400, "+", 0.8), None),  # starts before E - 50 = 350
        (term(360, 400, "-", 0.8), None),  # wrong strand
        (term(350, 420, "+", 0.8), 0.8),  # exactly the window [350, 420)
        (term(400, 421, "+", 0.8), None),  # ends past the downstream start
        (term(360, 400, "+", 0.8, contig="p2"), None),  # other contig
    ],
)
def test_example7_plus_pair(t: Terminator, want: float | None) -> None:
    assert terminator_between(UP_P, DOWN_P, [t], CHR1) == want


def test_example8_below_yes_is_none() -> None:
    assert terminator_between(UP_P, DOWN_P, [term(360, 400, "+", 0.4)], CHR1, yes=0.5) is None
    assert terminator_between(UP_P, DOWN_P, [term(360, 400, "+", 0.5)], CHR1, yes=0.5) == 0.5


def test_between_returns_the_maximum_confidence() -> None:
    ts = [term(360, 400, "+", 0.6), term(380, 410, "+", 0.9), term(355, 395, "+", 0.7)]
    assert terminator_between(UP_P, DOWN_P, ts, CHR1) == 0.9
    assert terminator_between(UP_P, DOWN_P, ts, CHR1, slack_bp=30) == 0.9  # [370, 420)
    assert terminator_between(UP_P, DOWN_P, ts[:1] + ts[2:], CHR1, slack_bp=30) is None


def test_between_minus_pair_mirrors_example7() -> None:
    """``-`` pair: up = higher coordinates ``[s_up, e_up)``, window ``[e_down, s_up + 50)``."""
    up, down = gene("U", 420, 700, "-"), gene("D", 100, 400, "-")  # window [400, 470)
    assert terminator_between(up, down, [term(420, 460, "-", 0.8)], CHR1) == 0.8
    assert terminator_between(up, down, [term(420, 475, "-", 0.8)], CHR1) is None
    assert terminator_between(up, down, [term(395, 460, "-", 0.8)], CHR1) is None
    assert terminator_between(up, down, [term(420, 460, "+", 0.8)], CHR1) is None


def test_overlapping_adjacent_genes_give_an_empty_window() -> None:
    """Edge case ``S < E - 50``: window empty → None."""
    up, down = gene("U", 100, 400, "+"), gene("D", 340, 700, "+")
    assert terminator_between(up, down, [term(340, 345, "+", 1.0)], CHR1) is None
    assert terminator_between(up, down, [term(350, 360, "+", 1.0)], CHR1) is None


def test_example10_bactermfinder_ignored_by_between_and_downstream() -> None:
    btf = term(360, 400, "+", 0.99, "BacTermFinder")
    assert terminator_between(UP_P, DOWN_P, [btf], CHR1) is None
    assert terminator_between(UP_P, DOWN_P, [btf, term(360, 400, "+", 0.6)], CHR1) == 0.6
    assert terminator_conf_downstream(UP_P, DOWN_P, [btf], CHR1) is None


def test_between_on_a_circular_replicon_across_the_origin() -> None:
    rep = Replicon("c", 10_000, "circular")
    up, down = gene("U", 9_500, 9_900, "+", "c"), gene("D", 100, 400, "+", "c")  # [9850, 10100)
    assert terminator_between(up, down, [term(9_950, 10_050, "+", 0.7, contig="c")], rep) == 0.7
    assert terminator_between(up, down, [term(20, 60, "+", 0.8, contig="c")], rep) == 0.8
    assert terminator_between(up, down, [term(90, 130, "+", 0.8, contig="c")], rep) is None
    mup, mdown = gene("U", 100, 400, "-", "c"), gene("D", 9_500, 9_900, "-", "c")  # [9900, 10150)
    assert terminator_between(mup, mdown, [term(9_950, 10_000, "-", 0.6, contig="c")], rep) == 0.6
    assert terminator_between(mup, mdown, [term(120, 145, "-", 0.65, contig="c")], rep) == 0.65
    assert terminator_between(mup, mdown, [term(9_880, 9_950, "-", 0.9, contig="c")], rep) is None


# ------------------------------------------------------------------- terminator_conf_downstream


def test_downstream_has_no_threshold_and_takes_the_max() -> None:
    ts = [term(360, 400, "+", 0.3), term(380, 415, "+", 0.2), term(500, 560, "+", 0.9)]
    assert terminator_conf_downstream(UP_P, DOWN_P, ts, CHR1) == 0.3
    assert terminator_conf_downstream(UP_P, DOWN_P, ts[1:], CHR1) == 0.2  # (500, 560) is past D
    assert terminator_conf_downstream(UP_P, DOWN_P, [], CHR1) is None
    assert terminator_conf_downstream(UP_P, DOWN_P, [term(345, 400, "+", 0.9)], CHR1) is None


def test_downstream_without_next_gene_uses_300_bp() -> None:
    """Q3: window up to 300 bp past the stop when there is no next same-strand gene."""
    g = gene("G", 100, 400, "+")
    assert terminator_conf_downstream(g, None, [term(600, 700, "+", 0.4)], CHR1) == 0.4
    assert terminator_conf_downstream(g, None, [term(650, 701, "+", 0.4)], CHR1) is None
    m = gene("M", 1_000, 1_300, "-")  # window [700, 1050)
    assert terminator_conf_downstream(m, None, [term(700, 760, "-", 0.45)], CHR1) == 0.45
    assert terminator_conf_downstream(m, None, [term(699, 760, "-", 0.45)], CHR1) is None


def test_downstream_window_clipped_at_a_linear_contig_end() -> None:
    lin = Replicon("L", 600, "linear")
    g = gene("G", 100, 400, "+", "L")
    assert terminator_conf_downstream(g, None, [term(500, 600, "+", 0.7, contig="L")], lin) == 0.7
    m = gene("M", 100, 400, "-", "L")
    assert terminator_conf_downstream(m, None, [term(0, 40, "-", 0.3, contig="L")], lin) == 0.3


def test_downstream_minus_gene_with_next_gene() -> None:
    g, nxt = gene("G", 420, 700, "-"), gene("N", 100, 400, "-")
    assert terminator_conf_downstream(g, nxt, [term(400, 470, "-", 0.25)], CHR1) == 0.25
    assert terminator_conf_downstream(g, nxt, [term(400, 471, "-", 0.25)], CHR1) is None


# ------------------------------------------------------------------- terminators.bed (Ex. 10)


def test_terminators_bed_names_scores_and_sources() -> None:
    ts = [
        term(1000, 1040, "+", 0.87),
        term(5000, 5030, "-", 0.8812014060630281),
        Terminator(Interval("chr1", 2000, 2030, "+"), 0.42, "BacTermFinder"),
    ]
    recs = terminators_bed(ts, {"chr1": CHR1})
    assert all(isinstance(r, BedRecord) for r in recs)
    assert all(r.blocks == () and r.thick is None for r in recs)  # BED6
    by = {r.name: r for r in recs}
    assert set(by) == {"TerminatorNet_1", "TerminatorNet_2", "BacTermFinder_1"}
    assert by["BacTermFinder_1"].interval == Interval("chr1", 2000, 2030, "+")
    assert by["BacTermFinder_1"].score == 420
    scores = sorted((r.interval.start, r.score) for r in recs if r.name.startswith("TerminatorNet"))
    assert scores == [(1000, 870), (5000, 881)]


def test_terminators_bed_is_independent_of_input_order() -> None:
    ts = parse_terminatornet(TN_SIM, SIM_REPS)
    first = terminators_bed(ts, SIM_REPS)
    shuffled = ts[:]
    random.Random(3).shuffle(shuffled)
    assert terminators_bed(shuffled, SIM_REPS) == first
    assert sorted(r.name for r in first) == sorted(f"TerminatorNet_{i}" for i in range(1, 134))


def test_origin_crossing_terminator_kept_unwrapped_and_split_by_bed() -> None:
    """Edge case: terminators crossing the origin keep ``end > length``; BED splits them."""
    t = term(99_980, 100_020, "+", 0.9)
    (rec,) = terminators_bed([t], {"chr1": CHR1})
    assert rec.interval == Interval("chr1", 99_980, 100_020, "+")
    lines = bed_lines(rec, CHR1)
    assert [ln.split("\t")[:3] for ln in lines] == [
        ["chr1", "99980", "100000"],
        ["chr1", "0", "20"],
    ]
