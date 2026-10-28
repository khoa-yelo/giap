"""geNOMAD command builder and summary parser (M2-18 examples 1-5; GIPPA_SPEC.md §7.1.1, §6.5).

Rules (spec M2-18): ``seq_name`` ``<contig>|provirus_<a>_<b>`` → ``Interval(contig, a-1, b, ".")``;
``seq_name`` without ``|provirus_`` → whole-contig virus (``interval=None``, Q1); ``NA`` → ``None``;
virus gene ``<seq_name>_<k>`` → contig = part before ``|`` (or before the last ``_<k>`` for a
whole-contig virus), strand ``1``/``-1`` → ``+``/``-``.

Ambiguities resolved here (test author):

* Q2 (gene coordinates relative to the provirus?) is decided by the recordings: geNOMAD 1.12.0
  writes **absolute** 1-based closed contig coordinates (K-12 gene ``…provirus_563848_584430_520``
  starts at 563848; ``tests/data/recorded/genomad/PROVENANCE.md``). Example 4 is therefore
  written with ``start 12201`` (absolute) for the expected ``Interval("chr1", 12200, 13100, "-")``;
  no ``a-1`` offset is applied.
* ``virus_hallmark`` is the 0/1 column of the same name; ``marker`` / ``annotation_description``
  are the ``marker`` / ``annotation_description`` columns.
* ``conjugation_genes`` / ``amr_genes`` are ``;``-separated, ``NA`` → ``()``.
* A virus ``seq_name`` whose contig is not in ``replicons`` → ``ValueError`` naming it (like the
  antiSMASH rule; needed to check the provirus against the contig length).
* Taxonomy is kept verbatim (including geNOMAD's trailing ``;;``).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval, Replicon, check_interval, to_one_based_closed

VS_HEADER = (
    "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\t"
    "n_hallmarks\tmarker_enrichment\ttaxonomy"
)
PS_HEADER = (
    "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\t"
    "marker_enrichment\tconjugation_genes\tamr_genes"
)
VG_HEADER = (
    "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\tevalue\t"
    "bitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\tannotation_conjscan\t"
    "annotation_amr\tannotation_accessions\tannotation_description"
)
TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"
REPS = {
    "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
    "p1": Replicon("p1", 8_000, "circular", "plasmid"),
    "plasmid_1": Replicon("plasmid_1", 8_000, "circular", "plasmid"),
}


def _g() -> ModuleType:
    from gippa.annotate import genomad

    return genomad


def _vs_row(
    seq_name: str,
    length: int,
    coords: str,
    score: str = "0.97",
    fdr: str = "NA",
    hallmarks: int = 9,
    tax: str = TAX,
) -> str:
    topo = "Provirus" if "|provirus_" in seq_name else "No terminal repeats"
    return (
        f"{seq_name}\t{length}\t{topo}\t{coords}\t42\t11\t{score}\t{fdr}\t{hallmarks}\t20.1\t{tax}"
    )


def _ps_row(
    seq_name: str,
    score: str = "0.93",
    fdr: str = "NA",
    hallmarks: int = 2,
    conj: str = "NA",
    amr: str = "NA",
) -> str:
    cols = [seq_name, "8000", "No terminal repeats", "8", "11", score, fdr, str(hallmarks), "5.3"]
    return "\t".join([*cols, conj, amr])


def _vg_row(
    gene: str,
    start: int,
    end: int,
    strand: int,
    marker: str = "NA",
    hallmark: int = 0,
    desc: str = "NA",
) -> str:
    return (
        f"{gene}\t{start}\t{end}\t{end - start + 1}\t{strand}\t0.5\t11\tNone\t{marker}\t"
        f"NA\tNA\t0\t0\t{hallmark}\t1\tNA\tNA\tNA\tNA\t{desc}"
    )


def _write(
    tmp: Path,
    virus: list[str] = (),
    plasmid: list[str] = (),
    genes: list[str] = (),
    prefix: str = "genome",
) -> Path:
    d = tmp / f"{prefix}_summary"
    d.mkdir(parents=True, exist_ok=True)
    for name, header, rows in (
        ("virus_summary", VS_HEADER, virus),
        ("plasmid_summary", PS_HEADER, plasmid),
        ("virus_genes", VG_HEADER, genes),
    ):
        (d / f"{prefix}_{name}.tsv").write_text("\n".join([header, *rows]) + "\n")
    return tmp


# ---------------------------------------------------------------- example 1: command


def test_example_1_command() -> None:
    from gippa.annotate.genomad_cmd import genomad_command

    got = genomad_command(Path("genome.fa"), Path("out"), Path("/db/genomad_db"), 16)
    assert got == [
        "genomad",
        "end-to-end",
        "--cleanup",
        "--splits",
        "8",
        "--threads",
        "16",
        "genome.fa",
        "out",
        "/db/genomad_db",
    ]
    assert all(isinstance(a, str) for a in got)


@pytest.mark.parametrize("threads", [1, 4, 128])
def test_command_threads_and_fixed_splits(threads: int) -> None:
    from gippa.annotate.genomad_cmd import genomad_command

    got = genomad_command(Path("g.fa"), Path("o"), Path("db/genomad_db"), threads)
    assert got[got.index("--threads") + 1] == str(threads)
    assert got[got.index("--splits") + 1] == "8"
    assert got[-3:] == ["g.fa", "o", "db/genomad_db"]


# ---------------------------------------------------------------- example 2: provirus


def test_example_2_provirus(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chr1|provirus_12001_43000", 31000, "12001-43000")])
    res = g.parse_genomad(out, "genome", REPS)
    assert res.proviruses == (
        g.Provirus("chr1", Interval("chr1", 12000, 43000, "."), 0.97, None, 9, TAX),
    )
    assert res.plasmids == () and res.virus_genes == ()


def test_provirus_fdr_value_is_float(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chr1|provirus_1_5000", 5000, "1-5000", fdr="0.0123")])
    (p,) = g.parse_genomad(out, "genome", REPS).proviruses
    assert p.fdr == pytest.approx(0.0123)
    assert p.interval == Interval("chr1", 0, 5000, ".")


def test_whole_contig_virus_has_no_interval_q1(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("p1", 8000, "NA", score="0.88", hallmarks=3)])
    (p,) = g.parse_genomad(out, "genome", REPS).proviruses
    assert p == g.Provirus("p1", None, 0.88, None, 3, TAX)


def test_provirus_ending_at_contig_end_is_valid(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chr1|provirus_90001_100000", 10000, "90001-100000")])
    (p,) = g.parse_genomad(out, "genome", REPS).proviruses
    assert p.interval == Interval("chr1", 90000, 100000, ".")


def test_provirus_beyond_contig_length_raises(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chr1|provirus_90001_100001", 10001, "90001-100001")])
    with pytest.raises(ValueError):
        g.parse_genomad(out, "genome", REPS)


def test_provirus_on_unknown_contig_raises_naming_it(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chrX|provirus_1_500", 500, "1-500")])
    with pytest.raises(ValueError, match="chrX"):
        g.parse_genomad(out, "genome", REPS)


# ---------------------------------------------------------------- example 3: plasmid


def test_example_3_plasmid(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, plasmid=[_ps_row("p1", conj="MOBP;T4SS_virB4", hallmarks=4)])
    (pc,) = g.parse_genomad(out, "genome", REPS).plasmids
    assert pc == g.PlasmidCall("p1", 0.93, None, 4, ("MOBP", "T4SS_virB4"), ())


def test_plasmid_amr_genes_and_fdr(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, plasmid=[_ps_row("p1", fdr="0.01", amr="AMR01;AMR02")])
    (pc,) = g.parse_genomad(out, "genome", REPS).plasmids
    assert pc.amr_genes == ("AMR01", "AMR02")
    assert pc.conjugation_genes == ()
    assert pc.fdr == pytest.approx(0.01)


# ---------------------------------------------------------------- example 4: virus genes


def test_example_4_virus_gene_absolute_coordinates_q2(tmp_path: Path) -> None:
    """Example 4 with the coordinate convention the recordings show (absolute, see module doc)."""
    g = _g()
    out = _write(
        tmp_path,
        virus=[_vs_row("chr1|provirus_12001_43000", 31000, "12001-43000")],
        genes=[
            _vg_row(
                "chr1|provirus_12001_43000_5", 12201, 13100, -1, "VV00123", 1, "phage integrase"
            )
        ],
    )
    (vg,) = g.parse_genomad(out, "genome", REPS).virus_genes
    assert vg == g.VirusGene(
        "chr1", Interval("chr1", 12200, 13100, "-"), "VV00123", "phage integrase", True
    )


def test_virus_gene_na_fields_and_plus_strand(tmp_path: Path) -> None:
    g = _g()
    out = _write(
        tmp_path,
        virus=[_vs_row("chr1|provirus_12001_43000", 31000, "12001-43000")],
        genes=[_vg_row("chr1|provirus_12001_43000_1", 12001, 12300, 1)],
    )
    (vg,) = g.parse_genomad(out, "genome", REPS).virus_genes
    assert vg == g.VirusGene("chr1", Interval("chr1", 12000, 12300, "+"), None, None, False)


def test_whole_contig_virus_genes_use_contig_with_underscores(tmp_path: Path) -> None:
    """``plasmid_1_3`` → contig ``plasmid_1`` (only the trailing ``_<k>`` is the gene number)."""
    g = _g()
    out = _write(
        tmp_path,
        virus=[_vs_row("plasmid_1", 8000, "NA")],
        genes=[_vg_row("plasmid_1_3", 10, 300, 1, "VV9", 1, "tail fiber")],
    )
    (vg,) = g.parse_genomad(out, "genome", REPS).virus_genes
    assert vg == g.VirusGene(
        "plasmid_1", Interval("plasmid_1", 9, 300, "+"), "VV9", "tail fiber", True
    )


# ---------------------------------------------------------------- example 5: empty


def test_example_5_headers_only(tmp_path: Path) -> None:
    g = _g()
    res = g.parse_genomad(_write(tmp_path), "genome", REPS)
    assert res == g.GenomadResult((), (), ())


def test_prefix_selects_the_summary_directory(tmp_path: Path) -> None:
    g = _g()
    out = _write(tmp_path, virus=[_vs_row("chr1|provirus_1_500", 500, "1-500")], prefix="sample")
    assert len(g.parse_genomad(out, "sample", REPS).proviruses) == 1
    with pytest.raises((FileNotFoundError, ValueError)):
        g.parse_genomad(out, "genome", REPS)


# ---------------------------------------------------------------- types and determinism


def test_result_types_are_frozen_slotted_tuples(tmp_path: Path) -> None:
    g = _g()
    out = _write(
        tmp_path,
        virus=[_vs_row("chr1|provirus_12001_43000", 31000, "12001-43000")],
        plasmid=[_ps_row("p1", conj="MOBP")],
        genes=[_vg_row("chr1|provirus_12001_43000_5", 12201, 13100, -1)],
    )
    res = g.parse_genomad(out, "genome", REPS)
    for obj in (res, res.proviruses[0], res.plasmids[0], res.virus_genes[0]):
        assert dataclasses.is_dataclass(obj)
        assert not hasattr(obj, "__dict__")  # slots
        field = dataclasses.fields(obj)[0].name
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(obj, field, getattr(obj, field))
    assert isinstance(res.proviruses, tuple) and isinstance(res.plasmids, tuple)
    assert isinstance(res.virus_genes, tuple)
    assert isinstance(res.plasmids[0].conjugation_genes, tuple)


def test_parsing_is_deterministic(tmp_path: Path) -> None:
    g = _g()
    out = _write(
        tmp_path,
        virus=[
            _vs_row("chr1|provirus_50001_60000", 10000, "50001-60000"),
            _vs_row("chr1|provirus_12001_43000", 31000, "12001-43000"),
            _vs_row("p1", 8000, "NA"),
        ],
        plasmid=[_ps_row("plasmid_1", conj="MOBP1")],
        genes=[
            _vg_row("chr1|provirus_12001_43000_5", 12201, 13100, -1),
            _vg_row("p1_1", 1, 300, 1),
        ],
    )
    assert g.parse_genomad(out, "genome", REPS) == g.parse_genomad(out, "genome", REPS)


# ---------------------------------------------------------------- Hypothesis properties


@st.composite
def _provirus(draw: st.DrawFn) -> tuple[int, int, int, int, int]:
    length = draw(st.integers(10, 200_000))
    a = draw(st.integers(1, length))
    b = draw(st.integers(a, length))
    s = draw(st.integers(a, b))
    e = draw(st.integers(s, b))
    return length, a, b, s, e


@settings(max_examples=60, deadline=None)
@given(_provirus(), st.sampled_from([1, -1]))
def test_property_one_based_closed_to_half_open(
    tmp_path_factory: pytest.TempPathFactory, pv: tuple[int, int, int, int, int], strand: int
) -> None:
    """Provirus and gene coordinates: 1-based closed → 0-based half-open, valid on the contig,
    the gene inside its provirus, round trip through ``to_one_based_closed``."""
    g = _g()
    length, a, b, s, e = pv
    reps = {"c_1": Replicon("c_1", length, "linear")}
    name = f"c_1|provirus_{a}_{b}"
    out = _write(
        tmp_path_factory.mktemp("gm"),
        virus=[_vs_row(name, b - a + 1, f"{a}-{b}")],
        genes=[_vg_row(f"{name}_7", s, e, strand)],
    )
    res = g.parse_genomad(out, "genome", reps)
    (p,) = res.proviruses
    (vg,) = res.virus_genes
    assert p.interval == Interval("c_1", a - 1, b, ".")
    assert to_one_based_closed(p.interval) == (a, b)
    check_interval(p.interval, reps["c_1"])
    assert vg.contig == "c_1"
    assert vg.interval == Interval("c_1", s - 1, e, "+" if strand == 1 else "-")
    assert p.interval.start <= vg.interval.start and vg.interval.end <= p.interval.end


@settings(max_examples=40, deadline=None)
@given(
    st.lists(st.from_regex(r"[A-Za-z0-9][A-Za-z0-9_.]{0,10}", fullmatch=True), max_size=5),
    st.lists(st.from_regex(r"[A-Za-z0-9][A-Za-z0-9_.]{0,10}", fullmatch=True), max_size=5),
)
def test_property_gene_lists_split_on_semicolon(
    tmp_path_factory: pytest.TempPathFactory, conj: list[str], amr: list[str]
) -> None:
    g = _g()
    row = _ps_row("p1", conj=";".join(conj) or "NA", amr=";".join(amr) or "NA")
    (pc,) = g.parse_genomad(
        _write(tmp_path_factory.mktemp("gm"), plasmid=[row]), "genome", REPS
    ).plasmids
    assert pc.conjugation_genes == tuple(conj)
    assert pc.amr_genes == tuple(amr)
