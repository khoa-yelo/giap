"""OSTIR input rows, writer, command and parser (M2-16 examples 1-8; GIAP_SPEC.md §7.1.1, §7.1.6).

Geometry (S = first start-codon base, E = interval end, 0-based): ``+`` → ``genome[S-50, S+50)``;
``-`` → ``revcomp(genome[E-50, E+50))``; circular replicons wrap, linear contigs clip and
``start = end = 1 + (bases kept upstream)``.

Ambiguities resolved here (test author, see the PR description):

* M2-16 Q1: the input header is ``id,seq,start,end,anti-Shine-Dalgarno`` (OSTIR 1.1.3 accepts it;
  the recorded runs used it); OSTIR's *output* id column is ``name``.
* ``ostir_command`` returns ``list[str]`` (paths and thread count stringified).
* ``parse_ostir`` fills ``RbsFeatures.start_codon`` from the row (it has no gene); the step
  ``rbs_layer`` replaces it by the gene's start codon (``tests/unit/test_ostir_steps.py``).
* Origin-crossing genes are stored unwrapped (``end > length``, ``giap.model.coords``); their
  windows wrap like any other window on a circular replicon.
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path

import pytest

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

ASD = "ACCTCCTTA"
N = 100_000
_COMP = str.maketrans("ACGT", "TGCA")
HEADER = "id,seq,start,end,anti-Shine-Dalgarno"
OSTIR_HEADER = (
    "name,start_codon,start_position,expression,RBS_distance_bp,dG_total,dG_rRNA:mRNA,dG_mRNA,"
    "dG_spacing,dG_standby,dG_start_codon"
)


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _genome() -> str:
    rnd = random.Random(1616)
    g = list("".join(rnd.choice("ACGT") for _ in range(N)))
    g[1000:1003] = "ATG"  # example 1 start codon (+)
    g[2897:2900] = revcomp("GTG")  # example 2 start codon (-): GTG on the minus strand
    g[20:23] = "TTG"  # examples 3/4
    g[99987:99990] = revcomp("ATG")  # example 5
    return "".join(g)


GENOME = _genome()
CIRC = Replicon("chr1", N, "circular", "chromosome")
LIN = Replicon("chr1", N, "linear", "chromosome")


def gene(start: int, end: int, strand: str, codon: str | None, tag: str = "GSIM_00001") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval("chr1", start, end, strand),  # type: ignore[arg-type]
        aa="",
        nt=None,
        start_codon=codon,
    )


# ---------------------------------------------------------------- examples 1-5: windows


def test_example_1_plus_strand_window() -> None:
    from giap.annotate.ostir import OstirRow, ostir_row

    row = ostir_row(gene(1000, 1300, "+", "ATG"), CIRC, GENOME, ASD)
    assert row == OstirRow("GSIM_00001", GENOME[950:1050], 51, 51, ASD)
    assert row.seq[50:53] == "ATG"


def test_example_2_minus_strand_window() -> None:
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(2000, 2900, "-", "GTG"), CIRC, GENOME, ASD)
    assert row.seq == revcomp(GENOME[2850:2950])
    assert (row.start, row.end) == (51, 51)
    assert row.seq[50:53] == revcomp(GENOME[2897:2900]) == "GTG"
    assert row.anti_sd == ASD and row.id == "GSIM_00001"


def test_example_3_linear_contig_clips_upstream() -> None:
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(20, 300, "+", "TTG"), LIN, GENOME, ASD)
    assert row.seq == GENOME[0:70]
    assert (row.start, row.end) == (21, 21)
    assert row.seq[row.start - 1 : row.start + 2] == "TTG"


def test_example_4_circular_contig_wraps_upstream() -> None:
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(20, 300, "+", "TTG"), CIRC, GENOME, ASD)
    assert row.seq == GENOME[99970:100000] + GENOME[0:70]
    assert (row.start, row.end) == (51, 51)
    assert row.seq[50:53] == "TTG"


def test_example_5_minus_gene_at_the_end_of_a_linear_contig() -> None:
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(99960, 99990, "-", "ATG"), LIN, GENOME, ASD)
    assert row.seq == revcomp(GENOME[99940:100000])
    assert len(row.seq) == 60  # 10 upstream kept + 50 downstream
    assert (row.start, row.end) == (11, 11)
    assert row.seq[10:13] == "ATG"


def test_minus_gene_at_the_end_of_a_circular_contig_wraps() -> None:
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(99960, 99990, "-", "ATG"), CIRC, GENOME, ASD)
    assert row.seq == revcomp(GENOME[99940:100000] + GENOME[0:40])
    assert (row.start, row.end, len(row.seq)) == (51, 51, 100)
    assert row.seq[50:53] == "ATG"


def test_minus_gene_at_the_start_of_a_linear_contig_clips_downstream_only() -> None:
    from giap.annotate.ostir import ostir_row

    g = GENOME[:10] + revcomp("ATG") + GENOME[13:]
    row = ostir_row(gene(5, 13, "-", "ATG"), LIN, g, ASD)  # E = 13: window [-37, 63) → [0, 63)
    assert row.seq == revcomp(g[0:63])
    assert (row.start, row.end) == (51, 51)  # all 50 upstream bases kept
    assert row.seq[50:53] == "ATG"


def test_plus_gene_near_the_end_of_a_linear_contig_clips_downstream_only() -> None:
    from giap.annotate.ostir import ostir_row

    g = GENOME[:99980] + "ATG" + GENOME[99983:]
    row = ostir_row(gene(99980, 99998, "+", "ATG"), LIN, g, ASD)
    assert row.seq == g[99930:100000]
    assert (row.start, row.end) == (51, 51)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_crossing_gene_unwrapped_end_wraps(strand: str) -> None:
    """Genes crossing the origin are stored unwrapped (end > length) and still get 100 nt."""
    from giap.annotate.ostir import ostir_row

    if strand == "+":
        g = GENOME[:99990] + "ATG" + GENOME[99993:]
        row = ostir_row(gene(99990, 100300, "+", "ATG"), CIRC, g, ASD)
        assert row.seq == g[99940:100000] + g[0:40]
    else:
        g = revcomp("ATG") + GENOME[3:]  # E = 100003 → last codon is g[0:3]
        row = ostir_row(gene(99700, 100003, "-", "ATG"), CIRC, g, ASD)
        assert row.seq == revcomp(g[99953:100000] + g[0:53])
    assert (row.start, row.end, len(row.seq)) == (51, 51, 100)
    assert row.seq[50:53] == "ATG"


def test_custom_window_sizes() -> None:
    """§5.2 ``ostir.upstream_nt`` / ``downstream_nt``: start = end = upstream + 1."""
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(1000, 1300, "+", "ATG"), CIRC, GENOME, ASD, up=30, down=20)
    assert row.seq == GENOME[970:1020]
    assert (row.start, row.end) == (31, 31)
    lin = ostir_row(gene(20, 300, "+", "TTG"), LIN, GENOME, ASD, up=30, down=20)
    assert lin.seq == GENOME[0:40] and (lin.start, lin.end) == (21, 21)


def test_lower_case_genome_gives_upper_case_rows() -> None:
    """The recorded inputs are upper case; soft-masked input must not leak into OSTIR."""
    from giap.annotate.ostir import ostir_row

    row = ostir_row(gene(1000, 1300, "+", "ATG"), CIRC, GENOME.lower(), ASD)
    assert row.seq == GENOME[950:1050]


def test_sorf_and_pseudogene_get_rows_too() -> None:
    """Edge cases: sORFs read +50 nt into the next gene; pseudogenes and partial genes get rows."""
    from giap.annotate.ostir import ostir_row

    sorf = dataclasses.replace(gene(1000, 1030, "+", "ATG"), feature_type="sorf")
    assert ostir_row(sorf, CIRC, GENOME, ASD).seq == GENOME[950:1050]
    pseudo = dataclasses.replace(gene(1000, 1300, "+", "ATG"), feature_type="pseudogene")
    assert ostir_row(pseudo, CIRC, GENOME, ASD).seq == GENOME[950:1050]
    partial = dataclasses.replace(gene(20, 300, "+", "TTG"), partial=True)
    assert ostir_row(partial, LIN, GENOME, ASD).start == 21


def test_ostir_row_is_frozen() -> None:
    from giap.annotate.ostir import OstirRow

    row = OstirRow("a", "ACGT", 1, 1, ASD)
    with pytest.raises(dataclasses.FrozenInstanceError):
        row.start = 2  # type: ignore[misc]


# ---------------------------------------------------------------- writer and command


def test_write_ostir_input_header_rows_and_bytes(tmp_path: Path) -> None:
    from giap.annotate.ostir import OstirRow, write_ostir_input

    rows = [
        OstirRow("GSIM_00002", "ACGTACGT", 5, 5, ASD),
        OstirRow("GSIM_00001", "TTTTATGCC", 5, 5, ASD),
    ]
    p = tmp_path / "ostir_input.csv"
    write_ostir_input(p, rows)
    assert (
        p.read_bytes()
        == (f"{HEADER}\nGSIM_00002,ACGTACGT,5,5,{ASD}\nGSIM_00001,TTTTATGCC,5,5,{ASD}\n").encode()
    )  # row order kept (§6.5 order comes from the caller), LF, trailing newline
    q = tmp_path / "again.csv"
    write_ostir_input(q, rows)
    assert q.read_bytes() == p.read_bytes()


def test_write_ostir_input_empty(tmp_path: Path) -> None:
    from giap.annotate.ostir import write_ostir_input

    p = tmp_path / "ostir_input.csv"
    write_ostir_input(p, [])
    assert p.read_text() == HEADER + "\n"


def test_ostir_command() -> None:
    from giap.annotate.ostir_cmd import ostir_command

    argv = ostir_command(Path("ostir_input.csv"), Path("ostir.csv"), 9)
    assert argv == ["ostir", "-i", "ostir_input.csv", "-o", "ostir.csv", "-j", "9"]
    assert all(isinstance(a, str) for a in argv)


# ---------------------------------------------------------------- examples 6-8: parser


def _csv(tmp_path: Path, *rows: str) -> Path:
    p = tmp_path / "ostir.csv"
    p.write_text("\n".join([OSTIR_HEADER, *rows]) + "\n")
    return p


ROW6 = "GSIM_00005,ATG,51,1523.4,7,-3.1,-9.8,-4.2,0.3,0.0,-1.19"


def test_example_6_row_to_features(tmp_path: Path) -> None:
    from giap.annotate.ostir import RbsFeatures, parse_ostir

    got = parse_ostir(_csv(tmp_path, ROW6), {"GSIM_00005": 51})
    assert got == {"GSIM_00005": RbsFeatures(1523.4, -9.8, -4.2, 0.3, 0.0, 7, "ATG")}
    f = got["GSIM_00005"]
    assert (f.tir_ostir, f.sd_dg, f.mrna_dg, f.spacing_dg, f.standby_dg) == (
        1523.4, -9.8, -4.2, 0.3, 0.0,
    )  # fmt: skip
    assert isinstance(f.spacer_bp, int) and f.spacer_bp == 7
    assert all(isinstance(x, float) for x in (f.tir_ostir, f.sd_dg, f.mrna_dg))


def test_example_7_only_the_requested_start_is_kept(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    p = _csv(
        tmp_path,
        "GSIM_00005,GTG,63,99.0,5,1.0,-2.0,-3.0,0.1,0.2,-0.07",  # other candidate first
        ROW6,
    )
    got = parse_ostir(p, {"GSIM_00005": 51})
    assert list(got) == ["GSIM_00005"]
    assert got["GSIM_00005"].tir_ostir == 1523.4 and got["GSIM_00005"].spacer_bp == 7


def test_example_7_requested_start_absent_drops_the_gene(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    p = _csv(tmp_path, "GSIM_00005,GTG,63,99.0,5,1.0,-2.0,-3.0,0.1,0.2,-0.07")
    assert parse_ostir(p, {"GSIM_00005": 51}) == {}


def test_example_8_missing_gene_is_absent(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    got = parse_ostir(_csv(tmp_path, ROW6), {"GSIM_00005": 51, "GSIM_00006": 51})
    assert set(got) == {"GSIM_00005"}
    assert "GSIM_00006" not in got


def test_unrequested_ids_are_ignored(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    assert parse_ostir(_csv(tmp_path, ROW6), {"GSIM_00009": 51}) == {}


def test_clipped_requested_start(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    p = _csv(tmp_path, ROW6.replace(",51,", ",21,"))
    assert set(parse_ostir(p, {"GSIM_00005": 21})) == {"GSIM_00005"}
    assert parse_ostir(p, {"GSIM_00005": 51}) == {}


def test_header_only_output(tmp_path: Path) -> None:
    from giap.annotate.ostir import parse_ostir

    assert parse_ostir(_csv(tmp_path), {"GSIM_00005": 51}) == {}
