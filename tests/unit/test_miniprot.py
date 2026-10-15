"""miniprot route: commands, GFF parser, CDS assignment, similarity class (issue M3-12).

GIAP_SPEC.md §7.2.3, §7.1.1 (miniprot gotcha), §9.3 (`protein_hits.tsv`), §6.5 (coordinates).
Examples 1-9 of docs/specs/M3-12-miniprot-route.md on hand-written GFF snippets; miniprot itself
is never run here.

Ambiguities resolved here (spec silent or open question):
- Module split: builders, ``ProteinAlignment`` and ``parse_miniprot_gff`` live in
  ``giap.mapping.miniprot``; ``reciprocal_overlap``, ``similarity_class`` and
  ``assign_alignments`` in ``giap.mapping.protein_assign``.
- The returned ``protein_hits`` frame is internal, so ``start``/``end`` are 0-based half-open
  (unwrapped); the 1-based conversion happens in the M3-16 writer.
- Thresholds are inclusive (``>=``) for identity, coverage and reciprocal overlap; the ``near``
  class means identity ``>= near_identity`` (not identical).
- Q1 (proposed default): an alignment passing identity and coverage where no CDS/sORF passes
  reciprocal overlap gives a row with ``gene`` NA and the protein stays a leftover.
- Q2 (proposed default): frameshifts are flags only, no identity penalty.
- "Largest overlap" is the overlap length in bp; ties go to the lower locus tag.
- Partial CDS are ordinary CDS for assignment; pseudogenes are not targets (spec: "CDS or sORF").
- Leftovers are protein IDs that occur in ``alns`` without a passing CDS assignment (proteins
  with no alignment at all are added by the caller), so an empty input gives an empty set.
- Output row order is a function of the alignment set (input order does not matter).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pandas as pd
import pytest

from giap.config.model import ProteinsConfig
from giap.model import Interval, Replicon
from giap.model.genes import Gene, GeneModel

L = 100_000
REPS = {"chr1": Replicon("chr1", L, "circular", "chromosome")}
LIN = {"ctg1": Replicon("ctg1", 20_000, "linear", "unknown")}
CFG = ProteinsConfig()

PROTEIN_HITS_COLUMNS = [
    "protein_id",
    "protein_len_aa",
    "route",
    "contig",
    "start",
    "end",
    "strand",
    "rank",
    "identity",
    "query_coverage",
    "reciprocal_overlap",
    "class",
    "gene",
    "flags",
]


# ---------------------------------------------------------------- helpers


def _mp():  # type: ignore[no-untyped-def]
    from giap.mapping import miniprot

    return miniprot


def _pa():  # type: ignore[no-untyped-def]
    from giap.mapping import protein_assign

    return protein_assign


def _gene(
    tag: str,
    start: int,
    end: int,
    strand: str = "+",
    contig: str = "chr1",
    ftype: str = "cds",
    partial: bool = False,
) -> Gene:
    n_aa = (end - start) // 3
    return Gene(
        locus_tag=tag,
        feature_type=ftype,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" + "A" * max(0, n_aa - 2),
        nt=None,
        start_codon="ATG",
        partial=partial,
    )


def _model(*genes: Gene, reps: dict[str, Replicon] = REPS) -> GeneModel:
    return GeneModel(reps, tuple(genes), ())


def _aln(
    pid: str = "P1",
    plen: int = 99,
    start: int = 100,
    end: int = 400,
    strand: str = "+",
    identity: float = 1.0,
    q_start: int = 0,
    q_end: int | None = None,
    rank: int = 1,
    contig: str = "chr1",
    frameshift: int = 0,
    stop_codon: int = 0,
    positive: float | None = 1.0,
):  # type: ignore[no-untyped-def]
    return _mp().ProteinAlignment(
        protein_id=pid,
        protein_len_aa=plen,
        route="miniprot",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        rank=rank,
        identity=identity,
        positive=positive,
        q_start=q_start,
        q_end=plen if q_end is None else q_end,
        frameshift=frameshift,
        stop_codon=stop_codon,
    )


def _mrna(
    start1: int,
    end1: int,
    target: str,
    *,
    strand: str = "+",
    contig: str = "chr1",
    attrs: str = "",
    mp_id: str = "MP000001",
    rank: int = 1,
    identity: str | None = "1.0000",
    positive: str | None = "1.0000",
) -> str:
    parts = [f"ID={mp_id}", f"Rank={rank}"]
    if identity is not None:
        parts.append(f"Identity={identity}")
    if positive is not None:
        parts.append(f"Positive={positive}")
    if attrs:
        parts.append(attrs)
    parts.append(f"Target={target}")
    cols = [contig, "miniprot", "mRNA", str(start1), str(end1), "500", strand, ".", ";".join(parts)]
    return "\t".join(cols)


def _cds_line(
    start1: int, end1: int, target: str, strand: str = "+", mp_id: str = "MP000001"
) -> str:
    attrs = f"Parent={mp_id};Rank=1;Identity=1.0000;Target={target}"
    return "\t".join(["chr1", "miniprot", "CDS", str(start1), str(end1), "500", strand, "0", attrs])


def _gff(tmp_path: Path, *lines: str, name: str = "aln.gff") -> Path:
    p = tmp_path / name
    p.write_text("##gff-version 3\n" + "".join(line + "\n" for line in lines))
    return p


def _flags(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)) or value is pd.NA:
        return []
    return [f for f in str(value).split(";") if f]


def _rows(df: pd.DataFrame, pid: str) -> pd.DataFrame:
    return df[df["protein_id"] == pid].reset_index(drop=True)


# ---------------------------------------------------------------- examples 1-2: commands


def test_example1_index_cmd_exact_argv() -> None:
    cmd = _mp().miniprot_index_cmd(Path("g"), Path("m"), 8)
    assert cmd == ["miniprot", "-t", "8", "-T", "11", "-d", "m", "g"]


def test_example2_align_cmd_exact_argv() -> None:
    cmd = _mp().miniprot_align_cmd(Path("m"), Path("p"), 8, "-S -N 5 --outs=0.95")
    assert cmd == ["miniprot", "-t", "8", "-S", "-N", "5", "--outs=0.95", "--gff", "m", "p"]


def test_align_cmd_default_config_args_and_no_translation_table() -> None:
    """§7.1.1 gotcha: ``-T 11`` is an indexing option, never passed when aligning."""
    cmd = _mp().miniprot_align_cmd(Path("genome.mpi"), Path("proteins.faa"), 4, CFG.miniprot_args)
    assert cmd == [
        "miniprot", "-t", "4", "-S", "-N", "5", "--outs=0.95", "--gff",
        "genome.mpi", "proteins.faa",
    ]  # fmt: skip
    assert "-T" not in cmd


def test_align_cmd_empty_args_and_paths_are_strings() -> None:
    cmd = _mp().miniprot_align_cmd(Path("/x/genome.mpi"), Path("/y/p.faa"), 1, "")
    assert cmd == ["miniprot", "-t", "1", "--gff", "/x/genome.mpi", "/y/p.faa"]
    assert all(isinstance(a, str) for a in cmd)


def test_align_cmd_args_split_on_any_whitespace() -> None:
    cmd = _mp().miniprot_align_cmd(Path("m"), Path("p"), 2, "  -S\t-N 5  ")
    assert cmd == ["miniprot", "-t", "2", "-S", "-N", "5", "--gff", "m", "p"]


def test_index_cmd_paths_are_strings() -> None:
    cmd = _mp().miniprot_index_cmd(Path("/a/genome.norm.fa"), Path("/b/genome.mpi"), 16)
    assert cmd == ["miniprot", "-t", "16", "-T", "11", "-d", "/b/genome.mpi", "/a/genome.norm.fa"]
    assert all(isinstance(a, str) for a in cmd)


# ---------------------------------------------------------------- parser


def test_example3_parse_mrna_line(tmp_path: Path) -> None:
    gff = _gff(
        tmp_path,
        "##PAF\tGSIM_00005\t99\t0\t99\t+\tchr1\t100000\t100\t400\t297\t297\t0\tAS:i:500",
        _mrna(101, 400, "GSIM_00005 1 99"),
        _cds_line(101, 400, "GSIM_00005 1 99"),
        "\t".join(
            ["chr1", "miniprot", "stop_codon", "401", "403", "0", "+", "0", "Parent=MP000001"]
        ),
    )
    alns = _mp().parse_miniprot_gff(gff, REPS, {"GSIM_00005": 99})
    assert len(alns) == 1, "one alignment per mRNA line; CDS/stop_codon/##PAF lines ignored"
    a = alns[0]
    assert a.protein_id == "GSIM_00005"
    assert a.protein_len_aa == 99
    assert a.route == "miniprot"
    assert a.interval == Interval("chr1", 100, 400, "+")
    assert a.rank == 1
    assert a.identity == pytest.approx(1.0)
    assert a.positive == pytest.approx(1.0)
    assert (a.q_start, a.q_end) == (0, 99)
    assert (a.frameshift, a.stop_codon) == (0, 0)
    assert a.query_coverage == pytest.approx(1.0)


def test_parse_minus_strand_from_column7(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(501, 800, "P2 11 90", strand="-"))
    (a,) = _mp().parse_miniprot_gff(gff, REPS, {"P2": 100})
    assert a.interval == Interval("chr1", 500, 800, "-")
    assert (a.q_start, a.q_end) == (10, 90)
    assert a.query_coverage == pytest.approx(0.8)


def test_parse_frameshift_stopcodon_counts(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(101, 400, "P1 1 99", attrs="Frameshift=1;StopCodon=2"))
    (a,) = _mp().parse_miniprot_gff(gff, REPS, {"P1": 99})
    assert (a.frameshift, a.stop_codon) == (1, 2)


def test_parse_positive_absent_is_none(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(101, 400, "P1 1 99", positive=None))
    (a,) = _mp().parse_miniprot_gff(gff, REPS, {"P1": 99})
    assert a.positive is None


def test_parse_identity_absent_raises_naming_the_line(tmp_path: Path) -> None:
    gff = _gff(
        tmp_path,
        _mrna(101, 400, "P1 1 99"),
        _mrna(1001, 1300, "P2 1 99", mp_id="MP000002", identity=None),
    )
    with pytest.raises(ValueError) as exc:
        _mp().parse_miniprot_gff(gff, REPS, {"P1": 99, "P2": 99})
    msg = str(exc.value)
    assert "Identity" in msg
    # the offending line is the 3rd line of the file (after ##gff-version)
    assert "3" in msg or "MP000002" in msg, msg


def test_parse_rank_and_multiple_alignments_in_file_order(tmp_path: Path) -> None:
    gff = _gff(
        tmp_path,
        _mrna(101, 400, "P1 1 99", rank=1),
        _mrna(5101, 5400, "P1 1 99", mp_id="MP000002", rank=2, identity="0.9800"),
    )
    alns = _mp().parse_miniprot_gff(gff, REPS, {"P1": 99})
    assert [a.rank for a in alns] == [1, 2]
    assert alns[1].interval == Interval("chr1", 5100, 5400, "+")
    assert alns[1].identity == pytest.approx(0.98)


def test_parse_empty_gff(tmp_path: Path) -> None:
    p = tmp_path / "empty.gff"
    p.write_text("")
    assert _mp().parse_miniprot_gff(p, REPS, {"P1": 99}) == []
    assert _mp().parse_miniprot_gff(_gff(tmp_path, name="hdr.gff"), REPS, {"P1": 99}) == []


def test_parse_unknown_contig_raises(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(101, 400, "P1 1 99", contig="nope"))
    with pytest.raises(ValueError):
        _mp().parse_miniprot_gff(gff, REPS, {"P1": 99})


def test_parse_unknown_protein_raises(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(101, 400, "PX 1 99"))
    with pytest.raises((KeyError, ValueError)):
        _mp().parse_miniprot_gff(gff, REPS, {"P1": 99})


def test_parse_is_deterministic(tmp_path: Path) -> None:
    gff = _gff(
        tmp_path,
        _mrna(101, 400, "P1 1 99"),
        _mrna(901, 1200, "P2 3 99", strand="-", mp_id="MP000002", identity="0.9100"),
    )
    a = _mp().parse_miniprot_gff(gff, REPS, {"P1": 99, "P2": 99})
    b = _mp().parse_miniprot_gff(gff, REPS, {"P1": 99, "P2": 99})
    assert a == b


def test_protein_alignment_is_frozen_with_coverage_property() -> None:
    a = _aln(plen=120, q_start=10, q_end=106)
    assert a.query_coverage == pytest.approx(96 / 120)
    with pytest.raises(dataclasses.FrozenInstanceError):
        a.identity = 0.5  # type: ignore[misc]


# ---------------------------------------------------------------- reciprocal overlap


def test_example3_reciprocal_overlap_value() -> None:
    ro = _pa().reciprocal_overlap(
        Interval("chr1", 100, 400, "+"), Interval("chr1", 100, 403, "+"), REPS["chr1"]
    )
    assert ro == pytest.approx(300 / 303)


def test_example7_reciprocal_overlap_half() -> None:
    ro = _pa().reciprocal_overlap(
        Interval("chr1", 100, 250, "+"), Interval("chr1", 100, 400, "+"), REPS["chr1"]
    )
    assert ro == pytest.approx(0.5)


def test_reciprocal_overlap_is_min_of_both_ratios() -> None:
    # alignment longer than the CDS: overlap / alignment span is the binding ratio
    ro = _pa().reciprocal_overlap(
        Interval("chr1", 100, 500, "+"), Interval("chr1", 100, 400, "+"), REPS["chr1"]
    )
    assert ro == pytest.approx(300 / 400)


def test_reciprocal_overlap_disjoint_is_zero() -> None:
    ro = _pa().reciprocal_overlap(
        Interval("chr1", 100, 400, "+"), Interval("chr1", 400, 700, "+"), REPS["chr1"]
    )
    assert ro == 0.0


def test_reciprocal_overlap_origin_crossing_cds_shifts_alignment() -> None:
    """CDS stored unwrapped ``[L-50, L+250)``; miniprot reports the linear ``[0, 250)``."""
    cds = Interval("chr1", L - 50, L + 250, "+")
    ro = _pa().reciprocal_overlap(Interval("chr1", 0, 250, "+"), cds, REPS["chr1"])
    assert ro == pytest.approx(min(250 / 300, 250 / 250))
    ro_pre = _pa().reciprocal_overlap(Interval("chr1", L - 50, L, "+"), cds, REPS["chr1"])
    assert ro_pre == pytest.approx(50 / 300)


# ---------------------------------------------------------------- similarity class


@pytest.mark.parametrize(
    ("identity", "coverage", "expected"),
    [
        (1.0, 1.0, "identical"),  # example 3
        (0.97, 0.92, "near"),  # example 4
        (0.85, 0.81, "divergent"),  # example 5
        (1.0, 0.95, "near"),  # identical needs coverage 1.0 too
        (0.95, 1.0, "near"),  # boundary inclusive
        (0.9499, 1.0, "divergent"),
        (0.80, 0.80, "divergent"),
    ],
)
def test_similarity_class(identity: float, coverage: float, expected: str) -> None:
    assert _pa().similarity_class(identity, coverage, 0.95) == expected


def test_similarity_class_uses_near_parameter() -> None:
    assert _pa().similarity_class(0.92, 1.0, 0.90) == "near"
    assert _pa().similarity_class(0.92, 1.0, 0.95) == "divergent"


# ---------------------------------------------------------------- assignment (examples 3-9)


def test_example3_identical_assigned_to_own_gene() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403))
    df, left = _pa().assign_alignments([_aln("GSIM_00005", 99, 100, 400)], gm, CFG)
    assert isinstance(df, pd.DataFrame)
    assert set(PROTEIN_HITS_COLUMNS) <= set(df.columns)
    assert len(df) == 1
    r = df.iloc[0]
    assert r["protein_id"] == "GSIM_00005"
    assert r["protein_len_aa"] == 99
    assert r["route"] == "miniprot"
    assert (r["contig"], r["start"], r["end"], r["strand"]) == ("chr1", 100, 400, "+")
    assert r["rank"] == 1
    assert r["identity"] == pytest.approx(1.0)
    assert r["query_coverage"] == pytest.approx(1.0)
    assert r["reciprocal_overlap"] == pytest.approx(300 / 303)
    assert round(float(r["reciprocal_overlap"]), 3) == 0.990
    assert r["class"] == "identical"
    assert r["gene"] == "GSIM_00005"
    assert _flags(r["flags"]) == []
    assert left == set()


def test_example4_near() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403))
    a = _aln("P4", 100, 100, 400, identity=0.97, q_start=4, q_end=96)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert df.iloc[0]["class"] == "near"
    assert df.iloc[0]["gene"] == "GSIM_00005"
    assert df.iloc[0]["query_coverage"] == pytest.approx(0.92)
    assert left == set()


def test_example5_divergent() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403))
    a = _aln("P5", 100, 100, 400, identity=0.85, q_start=10, q_end=91)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert df.iloc[0]["class"] == "divergent"
    assert df.iloc[0]["gene"] == "GSIM_00005"
    assert left == set()


def test_example6_identity_below_threshold_is_leftover() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403))
    a = _aln("P6", 99, 100, 400, identity=0.79)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == {"P6"}
    assert not (df["gene"] == "GSIM_00005").any()


def test_coverage_below_threshold_is_leftover() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403))
    a = _aln("P6b", 100, 100, 400, identity=1.0, q_start=0, q_end=79)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == {"P6b"}
    assert not (df["gene"] == "GSIM_00005").any()


def test_example7_low_reciprocal_overlap_is_leftover() -> None:
    gm = _model(_gene("GSIM_00007", 100, 400))
    a = _aln("P7", 50, 100, 250, identity=1.0)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == {"P7"}
    assert not (df["gene"] == "GSIM_00007").any()


def test_example8_rank2_second_copy_extra_row() -> None:
    gm = _model(_gene("GSIM_00005", 100, 403), _gene("GSIM_00020", 5100, 5403))
    alns = [
        _aln("P8", 99, 100, 400, rank=1),
        _aln("P8", 99, 5100, 5400, rank=2, identity=0.99),
    ]
    df, left = _pa().assign_alignments(alns, gm, CFG)
    rows = _rows(df, "P8").sort_values("rank").reset_index(drop=True)
    assert list(rows["rank"]) == [1, 2]
    assert list(rows["gene"]) == ["GSIM_00005", "GSIM_00020"]
    assert list(rows["start"]) == [100, 5100]
    assert list(rows["class"]) == ["identical", "near"]
    assert left == set()


def test_rank2_only_passing_alignment_still_assigns() -> None:
    """ "at least one alignment passes" — the passing one may be secondary."""
    gm = _model(_gene("GSIM_00005", 100, 403), _gene("GSIM_00020", 5100, 5403))
    alns = [
        _aln("P8b", 99, 100, 400, rank=1, identity=0.70),
        _aln("P8b", 99, 5100, 5400, rank=2, identity=0.90),
    ]
    df, left = _pa().assign_alignments(alns, gm, CFG)
    rows = _rows(df, "P8b")
    assert list(rows["gene"]) == ["GSIM_00020"]
    assert list(rows["rank"]) == [2]
    assert left == set()


def test_example9_flags_fs_st() -> None:
    gm = _model(_gene("GSIM_00009", 100, 403))
    a = _aln("P9", 99, 100, 400, identity=0.96, frameshift=1, stop_codon=1)
    df, left = _pa().assign_alignments([a], gm, CFG)
    r = df.iloc[0]
    assert r["flags"] == "fs;st"
    assert r["gene"] == "GSIM_00009"
    # Q2 default: no frameshift penalty on identity
    assert r["identity"] == pytest.approx(0.96)
    assert left == set()


@pytest.mark.parametrize(("fs", "st", "want"), [(2, 0, ["fs"]), (0, 3, ["st"]), (0, 0, [])])
def test_single_flags(fs: int, st: int, want: list[str]) -> None:
    gm = _model(_gene("G1", 100, 403))
    df, _ = _pa().assign_alignments(
        [_aln("P", 99, 100, 400, frameshift=fs, stop_codon=st)], gm, CFG
    )
    assert _flags(df.iloc[0]["flags"]) == want


def test_example9_parsed_end_to_end(tmp_path: Path) -> None:
    gff = _gff(
        tmp_path, _mrna(101, 400, "P9 1 99", identity="0.9600", attrs="Frameshift=1;StopCodon=1")
    )
    alns = _mp().parse_miniprot_gff(gff, REPS, {"P9": 99})
    df, left = _pa().assign_alignments(alns, _model(_gene("GSIM_00009", 100, 403)), CFG)
    assert df.iloc[0]["flags"] == "fs;st"
    assert df.iloc[0]["class"] == "near"
    assert left == set()


def test_example3_parsed_end_to_end(tmp_path: Path) -> None:
    gff = _gff(tmp_path, _mrna(101, 400, "GSIM_00005 1 99"))
    alns = _mp().parse_miniprot_gff(gff, REPS, {"GSIM_00005": 99})
    df, left = _pa().assign_alignments(alns, _model(_gene("GSIM_00005", 100, 403)), CFG)
    assert (df.iloc[0]["class"], df.iloc[0]["gene"]) == ("identical", "GSIM_00005")
    assert (df.iloc[0]["start"], df.iloc[0]["end"]) == (100, 400)  # 0-based internally
    assert left == set()


# ---------------------------------------------------------------- thresholds are inclusive


def test_thresholds_inclusive_at_080() -> None:
    gm = _model(_gene("G1", 100, 400))  # 300 bp
    # identity 0.80, coverage 80/100, reciprocal overlap 240/300 = 0.80
    a = _aln("P", 100, 100, 340, identity=0.80, q_start=0, q_end=80)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == set()
    assert list(df["gene"]) == ["G1"]
    assert df.iloc[0]["class"] == "divergent"


def test_thresholds_come_from_config() -> None:
    gm = _model(_gene("G1", 100, 403))
    a = _aln("P", 99, 100, 400, identity=0.85)
    strict = ProteinsConfig(min_identity=0.90)
    _, left = _pa().assign_alignments([a], gm, strict)
    assert left == {"P"}
    loose_cov = ProteinsConfig(min_coverage=0.5)
    b = _aln("Q", 100, 100, 400, identity=1.0, q_start=0, q_end=60)
    _, left = _pa().assign_alignments([b], gm, loose_cov)
    assert left == set()
    df, _ = _pa().assign_alignments(
        [_aln("R", 99, 100, 400, identity=0.92)], gm, ProteinsConfig(near_identity=0.90)
    )
    assert df.iloc[0]["class"] == "near"


# ---------------------------------------------------------------- strands, choice of CDS


def test_minus_strand_cds_assigned() -> None:
    gm = _model(_gene("GM", 1000, 1303, "-"))
    a = _aln("PM", 100, 1003, 1303, strand="-")
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert list(df["gene"]) == ["GM"]
    assert df.iloc[0]["strand"] == "-"
    assert left == set()


def test_opposite_strand_cds_not_assigned() -> None:
    gm = _model(_gene("GM", 1000, 1303, "-"))
    a = _aln("PX", 100, 1003, 1303, strand="+")
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == {"PX"}
    assert not (df["gene"] == "GM").any()


def test_opposite_strand_ignored_when_same_strand_cds_present() -> None:
    gm = _model(_gene("A_minus", 100, 403, "-"), _gene("B_plus", 100, 403, "+"))
    df, left = _pa().assign_alignments([_aln("P", 99, 100, 400, strand="+")], gm, CFG)
    assert list(df["gene"]) == ["B_plus"]
    assert left == set()


def test_largest_overlap_wins_over_lower_locus_tag() -> None:
    gm = _model(_gene("GSIM_00001", 130, 403), _gene("GSIM_00002", 100, 403))
    df, _ = _pa().assign_alignments([_aln("P", 99, 100, 400)], gm, CFG)
    assert list(df["gene"]) == ["GSIM_00002"]


def test_overlap_tie_goes_to_lower_locus_tag() -> None:
    gm = _model(_gene("GSIM_00009", 100, 403), _gene("GSIM_00003", 97, 400))
    df, _ = _pa().assign_alignments([_aln("P", 99, 100, 400)], gm, CFG)
    assert list(df["gene"]) == ["GSIM_00003"]


def test_sorf_is_a_target() -> None:
    gm = _model(_gene("S1", 2000, 2141, ftype="sorf"))  # 47 codons incl. stop
    df, left = _pa().assign_alignments([_aln("PS", 46, 2000, 2138)], gm, CFG)
    assert list(df["gene"]) == ["S1"]
    assert left == set()


def test_pseudogene_is_not_a_target() -> None:
    gm = _model(_gene("PSI", 3000, 3303, ftype="pseudogene"))
    df, left = _pa().assign_alignments([_aln("PP", 100, 3000, 3300)], gm, CFG)
    assert left == {"PP"}
    assert not (df["gene"] == "PSI").any()


def test_partial_gene_at_linear_contig_edge_assigned() -> None:
    gm = _model(_gene("PART", 19_700, 20_000, "+", contig="ctg1", partial=True), reps=LIN)
    a = _aln("PE", 100, 19_700, 20_000, contig="ctg1", q_start=0, q_end=100)
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert list(df["gene"]) == ["PART"]
    assert left == set()


def test_partial_alignment_on_full_gene_falls_below_thresholds() -> None:
    gm = _model(_gene("G1", 100, 403))
    a = _aln("PH", 99, 100, 250, q_start=0, q_end=50)  # half the protein, half the CDS
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert left == {"PH"}
    assert not (df["gene"] == "G1").any()


# ---------------------------------------------------------------- origin crossing


def test_origin_crossing_cds_assigned_via_unwrapped_overlap() -> None:
    gm = _model(_gene("ORI", L - 50, L + 250, "+"))
    a = _aln("PO", 99, 0, 250, q_start=16, q_end=99)  # coverage 83/99
    df, left = _pa().assign_alignments([a], gm, CFG)
    assert list(df["gene"]) == ["ORI"]
    assert df.iloc[0]["reciprocal_overlap"] == pytest.approx(250 / 300)
    assert left == set()


def test_origin_spanning_protein_split_in_two_is_leftover() -> None:
    """miniprot reports linear pieces; each covers < 80 % → tblastn (M3-13)."""
    gm = _model(_gene("ORI", L - 150, L + 150, "-"))
    alns = [
        _aln("PO2", 99, L - 150, L, strand="-", q_start=49, q_end=99, rank=1),
        _aln("PO2", 99, 0, 150, strand="-", q_start=0, q_end=50, rank=2),
    ]
    df, left = _pa().assign_alignments(alns, gm, CFG)
    assert left == {"PO2"}
    assert not (df["gene"] == "ORI").any()


# ---------------------------------------------------------------- Q1, empty, determinism


def test_q1_passing_alignment_without_cds_gives_na_row_and_leftover() -> None:
    gm = _model(_gene("FAR", 50_000, 50_303))
    df, left = _pa().assign_alignments([_aln("PN", 99, 100, 400)], gm, CFG)
    rows = _rows(df, "PN")
    assert len(rows) == 1
    assert pd.isna(rows.iloc[0]["gene"])
    assert (rows.iloc[0]["start"], rows.iloc[0]["end"]) == (100, 400)
    assert left == {"PN"}


def test_failing_alignment_without_cds_gives_no_row() -> None:
    gm = _model(_gene("FAR", 50_000, 50_303))
    df, left = _pa().assign_alignments([_aln("PF", 99, 100, 400, identity=0.5)], gm, CFG)
    assert _rows(df, "PF").empty
    assert left == {"PF"}


def test_empty_alignments() -> None:
    df, left = _pa().assign_alignments([], _model(_gene("G1", 100, 403)), CFG)
    assert df.empty
    assert set(PROTEIN_HITS_COLUMNS) <= set(df.columns)
    assert left == set()


def test_empty_gff_end_to_end(tmp_path: Path) -> None:
    alns = _mp().parse_miniprot_gff(_gff(tmp_path), REPS, {"P1": 99, "P2": 120})
    df, left = _pa().assign_alignments(alns, _model(_gene("G1", 100, 403)), CFG)
    assert df.empty and left == set()


def test_assigned_and_leftover_partition_the_aligned_proteins() -> None:
    gm = _model(_gene("G1", 100, 403), _gene("G2", 1000, 1303, "-"))
    alns = [
        _aln("A", 99, 100, 400),
        _aln("B", 99, 1000, 1300, strand="-", identity=0.6),
        _aln("C", 99, 1003, 1303, strand="-", identity=0.9),
        _aln("D", 99, 7000, 7300),
    ]
    df, left = _pa().assign_alignments(alns, gm, CFG)
    assigned = set(df.loc[df["gene"].notna(), "protein_id"])
    assert assigned == {"A", "C"}
    assert left == {"B", "D"}


def test_assignment_independent_of_input_order() -> None:
    gm = _model(_gene("G1", 100, 403), _gene("G2", 5100, 5403), _gene("G3", 1000, 1303, "-"))
    alns = [
        _aln("P1", 99, 100, 400),
        _aln("P1", 99, 5100, 5400, rank=2, identity=0.97),
        _aln("P2", 99, 1003, 1303, strand="-", identity=0.9),
        _aln("P3", 99, 9000, 9300),
    ]
    df1, left1 = _pa().assign_alignments(alns, gm, CFG)
    df2, left2 = _pa().assign_alignments(list(reversed(alns)), gm, CFG)
    pd.testing.assert_frame_equal(df1.reset_index(drop=True), df2.reset_index(drop=True))
    assert left1 == left2
