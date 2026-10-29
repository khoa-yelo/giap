"""TransTermHP ``genes.coords`` for a gene across the origin of a circular replicon (bug found by
the M. hyorhinis code-4 validation run, Slurm job 1903204; GIPPA_SPEC.md §6.5, §7.1.1, §7.1.5).

``write_transterm_coords`` wrote ``CAPMEI_0796`` (1-based 836688..838112 on an 836,897 bp
chromosome) unchanged; TransTermHP 2.09 exits 3 with ``Bad gene coordinates: 836688 .. 838112``
and ANNOTATE:TRANSTERMHP fails the run.

Lead decision (reversible), pinned here:

* ``write_transterm_coords(path, genes, replicons)`` takes ``replicons: Mapping[str, Replicon]``
  (the parsers' convention; called here by keyword). A gene whose unwrapped interval ends past the
  replicon length is written as two consecutive lines, both on the gene's strand, both with
  in-range 1-based coordinates: ``<locus_tag>`` from the gene start to the replicon end, then
  ``<locus_tag>_wrap`` from 1 to ``end - length``. Minus-strand pieces keep start > end.
* Every other gene is written byte-identically to the two-argument writer.
* ``transterm_input`` declares a ``replicons`` input (``replicons.parquet``) so the step can pass
  the lengths; ``workflow/tests/data/annotate/step_registry.json`` follows the registry.
* ``parse_transterm_tt`` keeps every terminator in range and never exposes a ``_wrap`` name:
  ``Terminator`` carries no gene name, so the pieces' names do not leak into the layer.
"""

from __future__ import annotations

import dataclasses
import importlib
from collections.abc import Iterator
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.annotate.terminators import parse_transterm_tt, write_transterm_coords
from gippa.config import GippaConfig
from gippa.model.coords import Interval, Replicon, check_interval
from gippa.model.genes import Gene

L_CAPMEI = 836_897
CHROM = "CP000000.1"
REPS = {CHROM: Replicon(CHROM, L_CAPMEI, "circular", "chromosome")}


def gene(tag: str, start: int, end: int, strand: str, contig: str = CHROM) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def _write(tmp: Path, genes: list[Gene], reps: dict[str, Replicon]) -> list[str]:
    p = tmp / "genes.coords"
    write_transterm_coords(p, genes, replicons=reps)
    return p.read_text().splitlines()


# ------------------------------------------------------------------- the CAPMEI_0796 case


def test_plus_gene_across_the_origin_is_written_as_two_in_range_pieces(tmp_path: Path) -> None:
    # 1-based 836688..838112 = 0-based unwrapped [836687, 838112); 838112 - 836897 = 1215
    lines = _write(tmp_path, [gene("CAPMEI_0796", 836_687, 838_112, "+")], REPS)
    assert lines == [
        f"CAPMEI_0796 836688 836897 {CHROM}",
        f"CAPMEI_0796_wrap 1 1215 {CHROM}",
    ]


def test_minus_gene_across_the_origin_keeps_start_after_end_in_each_piece(tmp_path: Path) -> None:
    lines = _write(tmp_path, [gene("CAPMEI_0796", 836_687, 838_112, "-")], REPS)
    assert lines == [
        f"CAPMEI_0796 836897 836688 {CHROM}",
        f"CAPMEI_0796_wrap 1215 1 {CHROM}",
    ]


def test_pieces_replace_the_gene_in_place_between_its_neighbours(tmp_path: Path) -> None:
    genes = [
        gene("CAPMEI_0001", 0, 300, "+"),
        gene("CAPMEI_0795", 835_000, 836_400, "-"),
        gene("CAPMEI_0796", 836_687, 838_112, "+"),
    ]
    lines = _write(tmp_path, genes, REPS)
    assert lines == [
        f"CAPMEI_0001 1 300 {CHROM}",
        f"CAPMEI_0795 836400 835001 {CHROM}",
        f"CAPMEI_0796 836688 836897 {CHROM}",
        f"CAPMEI_0796_wrap 1 1215 {CHROM}",
    ]


def test_one_base_past_the_origin(tmp_path: Path) -> None:
    lines = _write(tmp_path, [gene("G_1", 100, 201, "+", "p"), gene("G_2", 100, 201, "-", "p")],
                   {"p": Replicon("p", 200, "circular", "plasmid")})  # fmt: skip
    assert lines == ["G_1 101 200 p", "G_1_wrap 1 1 p", "G_2 200 101 p", "G_2_wrap 1 1 p"]


def test_only_the_crossing_contig_is_split(tmp_path: Path) -> None:
    reps = {**REPS, "pl": Replicon("pl", 5_000, "circular", "plasmid")}
    genes = [gene("A_1", 4_500, 5_400, "-", "pl"), gene("A_2", 4_500, 5_400, "+")]
    assert _write(tmp_path, genes, reps) == [
        "A_1 5000 4501 pl",
        "A_1_wrap 400 1 pl",
        f"A_2 4501 5400 {CHROM}",
    ]


# ------------------------------------------------------------------- non-crossing genes


def test_non_crossing_genes_are_byte_identical_to_the_two_argument_writer(tmp_path: Path) -> None:
    genes = [
        gene("B_1", 0, 300, "+"),
        gene("B_2", 400, 460, "-"),
        gene("B_3", L_CAPMEI - 300, L_CAPMEI, "+"),  # ends exactly at the replicon end
        gene("B_4", L_CAPMEI - 600, L_CAPMEI, "-"),
    ]
    old, new = tmp_path / "old.coords", tmp_path / "new.coords"
    write_transterm_coords(old, genes)
    write_transterm_coords(new, genes, replicons=REPS)
    assert new.read_bytes() == old.read_bytes()
    assert new.read_text() == (
        f"B_1 1 300 {CHROM}\nB_2 460 401 {CHROM}\n"
        f"B_3 836598 836897 {CHROM}\nB_4 836897 836298 {CHROM}\n"
    )


def test_recorded_sim100k_coords_unchanged_with_replicons(tmp_path: Path) -> None:
    data = Path(__file__).resolve().parents[2] / "data" / "recorded" / "terminators"
    d = data / "transtermhp" / "sim100k"
    rows = d.joinpath("genes.tsv").read_text().splitlines()[1:]
    genes = []
    for ln in rows:
        tag, _ft, contig, s, e, strand = ln.split("\t")
        genes.append(gene(tag, int(s), int(e), strand, contig))
    reps = {
        "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
        "plasmid_1": Replicon("plasmid_1", 8_000, "circular", "plasmid"),
    }
    p = tmp_path / "genes.coords"
    write_transterm_coords(p, genes, replicons=reps)
    assert p.read_bytes() == (d / "genes.coords").read_bytes()


# ------------------------------------------------------------------- properties


@st.composite
def genomes(draw: st.DrawFn) -> tuple[dict[str, Replicon], list[Gene]]:
    lengths = draw(st.lists(st.integers(10, 2_000), min_size=1, max_size=3))
    reps = {f"c{k}": Replicon(f"c{k}", n, "circular") for k, n in enumerate(lengths)}
    genes = []
    for i in range(draw(st.integers(0, 12))):
        c = draw(st.sampled_from(sorted(reps)))
        n = reps[c].length
        s = draw(st.integers(0, n - 1))
        e = draw(st.integers(s + 1, s + n))  # up to a full replicon length; may cross
        genes.append(gene(f"T_{i:04d}", s, e, draw(st.sampled_from("+-")), c))
    return reps, genes


@settings(max_examples=80, deadline=None)
@given(genomes())
def test_every_written_coordinate_is_in_range_and_names_are_unique(
    g: tuple[dict[str, Replicon], list[Gene]],
) -> None:
    reps, genes = g
    with TemporaryDirectory() as d:
        lines = _write(Path(d), genes, reps)
    names = [ln.split(" ")[0] for ln in lines]
    assert len(names) == len(set(names)), names
    n_cross = sum(x.interval.crosses_origin(reps[x.interval.contig]) for x in genes)
    assert len(lines) == len(genes) + n_cross
    for ln in lines:
        _tag, a, b, contig = ln.split(" ")
        assert 1 <= int(a) <= reps[contig].length and 1 <= int(b) <= reps[contig].length, ln


@settings(max_examples=80, deadline=None)
@given(genomes())
def test_pieces_cover_exactly_the_gene_on_its_strand(
    g: tuple[dict[str, Replicon], list[Gene]],
) -> None:
    reps, genes = g
    with TemporaryDirectory() as d:
        lines = _write(Path(d), genes, reps)
    by_name = {}
    for ln in lines:
        tag, a, b, contig = ln.split(" ")
        by_name[tag] = (int(a), int(b), contig)
    for x in genes:
        iv, n = x.interval, reps[x.interval.contig].length
        pieces = [by_name[x.locus_tag]]
        if iv.crosses_origin(reps[iv.contig]):
            pieces.append(by_name[f"{x.locus_tag}_wrap"])
        else:
            assert f"{x.locus_tag}_wrap" not in by_name
        covered: list[int] = []
        for a, b, contig in pieces:
            assert contig == iv.contig
            assert (a <= b) if iv.strand == "+" else (a >= b)
            covered += range(min(a, b) - 1, max(a, b))
        assert sorted(covered) == sorted(p % n for p in range(iv.start, iv.end))
        first = pieces[0]
        assert min(first[0], first[1]) == iv.start + 1  # the first piece starts at the gene start


# ------------------------------------------------------------------- .tt parser around the origin


def _tt(lines: list[str]) -> str:
    head = ["TransTermHP v2.09", "", f"SEQUENCE {CHROM}  (length {L_CAPMEI})", ""]
    return "\n".join(head + lines) + "\n"


def test_tt_with_a_wrap_piece_parses_in_range_without_wrap_names(tmp_path: Path) -> None:
    p = tmp_path / "t.tt"
    p.write_text(
        _tt(
            [
                "CAPMEI_0796_wrap        1 - 1215     + | ",
                "  TERM 1     1240 - 1268    +   F   96 -14.2 -4.1 | ",
                "  AAAAAAAAAAAAAAA      GCCGCC TTTT GGCGGC     TTTTTTTTTTTTTTT",
                "CAPMEI_0001          1301 - 1600     + | ",
                "CAPMEI_0795      836400 - 835001     - | ",
                "  TERM 2   834950 - 834975  -   R   88 -12.0 -3.5 | ",
                "  AAAAAAAAAAAAAAA      GCCGCC TTTT GGCGGC     TTTTTTTTTTTTTTT",
                "CAPMEI_0796      836688 - 836897     + | ",
            ]
        )
    )
    got = parse_transterm_tt(p, REPS)
    assert [t.interval for t in got] == [
        Interval(CHROM, 1239, 1268, "+"),
        Interval(CHROM, 834949, 834975, "-"),
    ]
    for t in got:
        check_interval(t.interval, REPS[CHROM])
        assert t.interval.end <= L_CAPMEI
        assert not any("_wrap" in str(v) for v in dataclasses.astuple(t))


def test_tt_term_past_the_replicon_end_is_rejected(tmp_path: Path) -> None:
    """TransTermHP scans each sequence linearly: a TERM ending past the length (or starting past
    it) is not a TransTermHP output and must not be accepted as an in-range call."""
    p = tmp_path / "t.tt"
    p.write_text(_tt(["  TERM 1   836900 - 836920  +   F   90 -11.0 -3.0 | "]))
    with pytest.raises(ValueError):
        parse_transterm_tt(p, REPS)


# ------------------------------------------------------------------- transterm_input step


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


def test_transterm_input_declares_replicons(steps: ModuleType) -> None:
    spec = steps.get("transterm_input")
    assert "replicons" in spec.inputs, spec.inputs


def _genes_parquet(path: Path, genes: list[Gene]) -> None:
    cols = {
        "locus_tag": [g.locus_tag for g in genes],
        "feature_type": [g.feature_type for g in genes],
        "contig": [g.interval.contig for g in genes],
        "start": pa.array([g.interval.start for g in genes], pa.int64()),
        "end": pa.array([g.interval.end for g in genes], pa.int64()),
        "strand": [g.interval.strand for g in genes],
        "aa": [g.aa for g in genes],
        "start_codon": [g.start_codon for g in genes],
    }
    pq.write_table(pa.table(cols), path)


def _replicons_parquet(path: Path, reps: dict[str, Replicon]) -> None:
    rs = list(reps.values())
    pq.write_table(
        pa.table(
            {
                "contig": [r.id for r in rs],
                "length": pa.array([r.length for r in rs], pa.int64()),
                "topology": [r.topology for r in rs],
                "kind": [r.kind for r in rs],
            }
        ),
        path,
    )


def test_transterm_input_step_splits_the_origin_crossing_gene(
    steps: ModuleType, tmp_path: Path
) -> None:
    genes = [gene("CAPMEI_0001", 0, 300, "+"), gene("CAPMEI_0796", 836_687, 838_112, "-")]
    gp, rp = tmp_path / "genes.parquet", tmp_path / "replicons.parquet"
    _genes_parquet(gp, genes)
    _replicons_parquet(rp, REPS)
    (key,) = steps.get("transterm_input").outputs
    out = tmp_path / "genes.coords"
    steps.run_step(
        "transterm_input", GippaConfig(seed=1), {"genes": gp, "replicons": rp}, {key: out}
    )
    assert out.read_text() == (
        f"CAPMEI_0001 1 300 {CHROM}\n"
        f"CAPMEI_0796 836897 836688 {CHROM}\n"
        f"CAPMEI_0796_wrap 1215 1 {CHROM}\n"
    )
