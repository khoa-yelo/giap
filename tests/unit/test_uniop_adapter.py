"""UniOP adapter: Prodigal-style GFF, ``uniop.pred`` parser, command (M2-12; GIPPA_SPEC.md §7.1.1,
§7.1.3; Module 1 test "adapter round-trip for UniOP (Bakta → Prodigal-style GFF → back to locus
tags)").

Examples 4 and 10 of docs/specs/M2-12-uniop-pairs-tus.md. UniOP has no licence (ADR-0008) and is
never imported or run here.

Conventions pinned here (lead decisions after code review, PR #49, items 1-2):

* ``uniop.pred`` is UniOP's real format (pin 3ed14dfc; excerpt in ``tests/data/uniop/``): a title
  line ``operon prediction scores``, a header ``Gene A<TAB>Gene B<TAB>Prediction``, then rows of
  1-based gene ROW NUMBERS ``k<TAB>k+1<TAB>p`` (p printed with a trailing space; blank-padded empty
  when the pair is on opposite strands or the gene is the last of its contig). Empty rows are
  skipped.
* Row numbers index UniOP's own gene order: every GFF gene, stable-sorted by start (1-based column
  4) across contigs. Row ``k`` is sorted gene ``k`` paired with the NEXT gene ON THE SAME CONTIG
  (the printed ``Gene B`` number is ``k+1`` and is not used for pairing).
* ``parse_uniop_pred(pred, gff, id_map)`` rebuilds that order from the exact GFF that
  ``write_prodigal_gff`` wrote and returns locus-tag pairs in row order. A row with a probability
  whose partner is missing or on the other strand, a row number outside the gene list, a GFF ID
  absent from ``id_map``, or a missing title/header raises ``ValueError``. An empty file, or
  title + header only, parses to ``[]``.
* ``uniop_map.tsv`` (step ``uniop_input``) has the header ``prodigal_id  locus_tag``.
* The GFF starts with a ``##gff-version`` line; each gene is one ``CDS`` line with 1-based closed
  coordinates and the attribute column starting with ``ID=<contig>_<k>``.
* ``uniop_command`` runs ``python <uniop_home>/src/UniOP`` (entry point ``UniOP`` or ``UniOP.py``,
  "confirm the entry-point name on the pinned commit").
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import gippasim
from gippa.annotate.bakta import parse_bakta_json
from gippa.annotate.operons import adjacent_pairs
from gippa.annotate.uniop_adapter import PairProb, parse_uniop_pred, write_prodigal_gff
from gippa.annotate.uniop_cmd import uniop_command
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

CHR = Replicon("chr1", 100_000, "circular", "chromosome")
PL = Replicon("p2", 8_000, "circular", "plasmid")


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


DATA = Path(__file__).resolve().parents[1] / "data" / "uniop"
TITLE = "operon prediction scores"
HEADER = "Gene A\tGene B\tPrediction"
BLANK = " " * 9  # UniOP's empty probability cell


def write_gff(path: Path, lines: list[tuple[str, int, int, str, str]]) -> Path:
    """(contig, start1, end1, strand, id) lines in the given order, Prodigal-style."""
    body = [f"{c}\tGIPPA\tCDS\t{a}\t{b}\t.\t{s}\t0\tID={i}" for c, a, b, s, i in lines]
    path.write_text("\n".join(["##gff-version 3", *body]) + "\n")
    return path


def write_pred(path: Path, rows: list[tuple[int, float | None]]) -> Path:
    """Real-format ``uniop.pred``: (row number k, p or None for UniOP's empty cell)."""
    body = [f"{k}\t{k + 1}\t" + (BLANK if p is None else f"{p:.6f} ") for k, p in rows]
    path.write_text("\n".join([TITLE, HEADER, *body]) + "\n")
    return path


def uniop_like_pred(gff: Path, path: Path, p: float) -> Path:
    """What UniOP 3ed14dfc writes for ``gff`` (probability ``p`` on every same-strand pair)."""
    rows = gff_rows(gff)
    order = sorted(range(len(rows)), key=lambda i: int(rows[i][3]))  # stable
    out: list[tuple[int, float | None]] = []
    for k, i in enumerate(order[:-1], start=1):
        nxt = next((rows[j] for j in order[k:] if rows[j][0] == rows[i][0]), None)
        out.append((k, p if nxt is not None and nxt[6] == rows[i][6] else None))
    return write_pred(path, out)


def fixture_map() -> dict[str, str]:
    lines = (DATA / "multicontig_map.tsv").read_text().splitlines()
    assert lines[0] == "prodigal_id\tlocus_tag"
    return dict(ln.split("\t") for ln in lines[1:])


def gff_rows(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


def gff_id(row: list[str]) -> str:
    m = re.match(r"ID=([^;]+)", row[8])
    assert m, f"attributes must start with ID=: {row[8]!r}"
    return m.group(1)


# ---------------------------------------------------------------- example 4 (real format)


def test_example4_pred_row_maps_to_locus_tags(tmp_path: Path) -> None:
    gff = write_gff(
        tmp_path / "in.gff", [("chr1", 101, 400, "+", "chr1_1"), ("chr1", 421, 700, "+", "chr1_2")]
    )
    pred = write_pred(tmp_path / "uniop.pred", [(1, 0.91)])
    got = parse_uniop_pred(pred, gff, {"chr1_1": "GSIM_00005", "chr1_2": "GSIM_00010"})
    assert got == [PairProb("GSIM_00005", "GSIM_00010", 0.91, "uniop")]


def test_real_uniop_excerpt_maps_rows_through_uniop_gene_order() -> None:
    """Real SIM100K excerpt: plasmid_1_1 (start 193) sorts before chr1_1 (288), so row 1 is the
    plasmid pair, row 3 (last plasmid gene) and row 4 (chr1_2 + / chr1_3 -) are empty."""
    pred = DATA / "multicontig_uniop.pred"
    lines = pred.read_text().splitlines()
    assert lines[:2] == [TITLE, HEADER]
    got = parse_uniop_pred(pred, DATA / "multicontig.gff", fixture_map())
    assert got == [
        PairProb("GSIM_00505", "GSIM_00510", 0.735206, "uniop"),
        PairProb("GSIM_00005", "GSIM_00010", 0.728063, "uniop"),
        PairProb("GSIM_00015", "GSIM_00020", 0.734771, "uniop"),
        PairProb("GSIM_00020", "GSIM_00025", 0.728503, "uniop"),
        PairProb("GSIM_00025", "GSIM_00030", 0.728210, "uniop"),
    ]


def test_real_excerpt_naive_row_mapping_would_be_wrong() -> None:
    """Guard on the fixture itself: reading row k as GFF line k / k+1 gives other pairs."""
    ids = [gff_id(r) for r in gff_rows(DATA / "multicontig.gff")]
    m = fixture_map()
    naive = (m[ids[0]], m[ids[1]])
    got = parse_uniop_pred(DATA / "multicontig_uniop.pred", DATA / "multicontig.gff", m)
    assert (got[0].gene_a, got[0].gene_b) != naive
    contig = {m[gff_id(r)]: r[0] for r in gff_rows(DATA / "multicontig.gff")}
    assert all(contig[p.gene_a] == contig[p.gene_b] for p in got)


def test_emulator_matches_the_real_excerpt_blank_pattern(tmp_path: Path) -> None:
    real = (DATA / "multicontig_uniop.pred").read_text().splitlines()
    fake = uniop_like_pred(DATA / "multicontig.gff", tmp_path / "f.pred", 0.5)
    fake_lines = fake.read_text().splitlines()
    assert len(fake_lines) == len(real)
    for r, f in zip(real[2:], fake_lines[2:], strict=True):
        assert r.split("\t")[:2] == f.split("\t")[:2]
        assert (r.split("\t")[2].strip() == "") == (f.split("\t")[2].strip() == "")


def test_stable_sort_by_start_across_contigs_with_ties(tmp_path: Path) -> None:
    """chrA_1 and chrB_1 share start 101: GFF file order breaks the tie (stable sort)."""
    gff = write_gff(
        tmp_path / "in.gff",
        [
            ("chrA", 101, 300, "+", "chrA_1"),
            ("chrA", 501, 800, "+", "chrA_2"),
            ("chrB", 101, 200, "+", "chrB_1"),
            ("chrB", 301, 400, "+", "chrB_2"),
        ],
    )
    # UniOP order: chrA_1, chrB_1, chrB_2, chrA_2
    pred = write_pred(tmp_path / "uniop.pred", [(1, 0.9), (2, 0.8), (3, None)])
    id_map = {"chrA_1": "A1", "chrA_2": "A2", "chrB_1": "B1", "chrB_2": "B2"}
    got = parse_uniop_pred(pred, gff, id_map)
    assert got == [PairProb("A1", "A2", 0.9, "uniop"), PairProb("B1", "B2", 0.8, "uniop")]


def test_pred_with_several_rows_keeps_row_order_and_values(tmp_path: Path) -> None:
    gff = write_gff(
        tmp_path / "in.gff",
        [
            ("chr1", 101, 400, "+", "chr1_1"),
            ("chr1", 421, 700, "+", "chr1_2"),
            ("chr1", 1001, 1300, "-", "chr1_3"),
            ("chr1", 1321, 1600, "-", "chr1_4"),
            ("p2", 2001, 2300, "+", "p2_1"),
            ("p2", 2401, 2600, "+", "p2_2"),
        ],
    )
    pred = write_pred(
        tmp_path / "uniop.pred", [(1, 0.91), (2, None), (3, 0.0), (4, None), (5, 1.0)]
    )
    id_map = {"chr1_1": "A", "chr1_2": "B", "chr1_3": "C", "chr1_4": "D", "p2_1": "E", "p2_2": "F"}
    got = parse_uniop_pred(pred, gff, id_map)
    assert [(p.gene_a, p.gene_b, p.uniop_p, p.source) for p in got] == [
        ("A", "B", 0.91, "uniop"),
        ("C", "D", 0.0, "uniop"),
        ("E", "F", 1.0, "uniop"),
    ]


@pytest.mark.parametrize(
    "content", ["", f"{TITLE}\n{HEADER}\n", f"{TITLE}\n{HEADER}\n1\t2\t{BLANK}\n"]
)
def test_pred_without_probabilities_parses_to_nothing(tmp_path: Path, content: str) -> None:
    gff = write_gff(
        tmp_path / "in.gff", [("chr1", 101, 400, "+", "chr1_1"), ("chr1", 421, 700, "-", "chr1_2")]
    )
    pred = tmp_path / "uniop.pred"
    pred.write_text(content)
    assert parse_uniop_pred(pred, gff, {"chr1_1": "A", "chr1_2": "B"}) == []


def _three_gene_gff(tmp_path: Path) -> Path:
    return write_gff(
        tmp_path / "in.gff",
        [
            ("chr1", 101, 400, "+", "chr1_1"),
            ("chr1", 421, 700, "-", "chr1_2"),
            ("p2", 201, 500, "+", "p2_1"),
        ],
    )  # UniOP order: chr1_1 (101), p2_1 (201), chr1_2 (421)


@pytest.mark.parametrize(
    "rows",
    [
        [(1, 0.5)],  # chr1_1 + -> chr1_2 -: opposite strands cannot carry a probability
        [(2, 0.5)],  # p2_1 is the last gene of p2: no partner
        [(3, 0.5)],  # row 3 would name a 4th gene
        [(7, 0.5)],
    ],
    ids=["opposite_strand", "no_partner_on_contig", "last_row_number", "out_of_range"],
)
def test_pred_row_without_valid_same_contig_partner_raises(
    tmp_path: Path, rows: list[tuple[int, float | None]]
) -> None:
    pred = write_pred(tmp_path / "uniop.pred", rows)
    with pytest.raises(ValueError):
        parse_uniop_pred(
            pred, _three_gene_gff(tmp_path), {"chr1_1": "A", "chr1_2": "B", "p2_1": "C"}
        )


def test_pred_gff_id_missing_from_map_raises(tmp_path: Path) -> None:
    gff = write_gff(
        tmp_path / "in.gff", [("chr1", 101, 400, "+", "chr1_1"), ("chr1", 421, 700, "+", "chr1_9")]
    )
    pred = write_pred(tmp_path / "uniop.pred", [(1, 0.5)])
    with pytest.raises(ValueError, match="chr1_9"):
        parse_uniop_pred(pred, gff, {"chr1_1": "A", "chr1_2": "B"})


@pytest.mark.parametrize(
    "content",
    ["gene_a\tgene_b\tp\nchr1_1\tchr1_2\t0.5\n", f"{HEADER}\n1\t2\t0.5 \n"],
    ids=["invented_3col", "no_title"],
)
def test_pred_not_in_uniop_format_raises(tmp_path: Path, content: str) -> None:
    gff = write_gff(
        tmp_path / "in.gff", [("chr1", 101, 400, "+", "chr1_1"), ("chr1", 421, 700, "+", "chr1_2")]
    )
    pred = tmp_path / "uniop.pred"
    pred.write_text(content)
    with pytest.raises(ValueError):
        parse_uniop_pred(pred, gff, {"chr1_1": "A", "chr1_2": "B"})


def test_pairprob_is_frozen() -> None:
    p = PairProb("A", "B", 0.5, "uniop")
    with pytest.raises(AttributeError):
        p.uniop_p = 0.1  # type: ignore[misc]


# ---------------------------------------------------------------- Prodigal-style GFF


def test_write_prodigal_gff_ids_in_genome_order(tmp_path: Path) -> None:
    genes = [
        gene("P_2", 400, 700, "-", "p2"),
        gene("C_2", 1000, 1300, "-"),
        gene("P_1", 10, 300, "+", "p2"),
        gene("C_1", 100, 400, "+"),
        gene("C_3", 1290, 1500, "+", ft="sorf"),
        gene("C_4", 2000, 2600, "+", ft="pseudogene"),
    ]
    path = tmp_path / "uniop_cds.gff"
    id_map = write_prodigal_gff(path, genes, [CHR, PL])
    assert id_map == [
        ("chr1_1", "C_1"),
        ("chr1_2", "C_2"),
        ("chr1_3", "C_3"),
        ("chr1_4", "C_4"),
        ("p2_1", "P_1"),
        ("p2_2", "P_2"),
    ]
    text = path.read_text()
    assert text.splitlines()[0].startswith("##gff-version")
    rows = gff_rows(path)
    assert all(len(r) == 9 for r in rows)
    assert [gff_id(r) for r in rows] == [pid for pid, _ in id_map]
    assert [(r[0], r[2], r[3], r[4], r[6]) for r in rows] == [
        ("chr1", "CDS", "101", "400", "+"),
        ("chr1", "CDS", "1001", "1300", "-"),
        ("chr1", "CDS", "1291", "1500", "+"),
        ("chr1", "CDS", "2001", "2600", "+"),
        ("p2", "CDS", "11", "300", "+"),
        ("p2", "CDS", "401", "700", "-"),
    ]


def test_write_prodigal_gff_contig_order_follows_replicons(tmp_path: Path) -> None:
    genes = [gene("C_1", 100, 400, "+"), gene("P_1", 10, 300, "+", "p2")]
    id_map = write_prodigal_gff(tmp_path / "a.gff", genes, [PL, CHR])
    assert id_map == [("p2_1", "P_1"), ("chr1_1", "C_1")]


def test_write_prodigal_gff_empty(tmp_path: Path) -> None:
    path = tmp_path / "e.gff"
    assert write_prodigal_gff(path, [], [CHR]) == []
    assert gff_rows(path) == []


def test_write_prodigal_gff_is_deterministic(tmp_path: Path) -> None:
    genes = [gene("C_2", 1000, 1300, "-"), gene("C_1", 100, 400, "+")]
    write_prodigal_gff(tmp_path / "a.gff", genes, [CHR])
    write_prodigal_gff(tmp_path / "b.gff", genes[::-1], [CHR])
    assert (tmp_path / "a.gff").read_bytes() == (tmp_path / "b.gff").read_bytes()


# ---------------------------------------------------------------- example 10: round trip


def _sim_genes(cfg: gippasim.SimConfig, seed: int, out: Path) -> tuple[list[Gene], list[Replicon]]:
    res = gippasim.simulate(cfg, seed, out)
    model = parse_bakta_json(res.bakta_json).model
    return list(model.genes), list(model.replicons.values())


def _non_wrap(genes: list[Gene], reps: list[Replicon]) -> set[tuple[str, str]]:
    pairs = adjacent_pairs(genes, {r.id: r for r in reps})
    by_contig: dict[str, list[Gene]] = {}
    for g in genes:
        by_contig.setdefault(g.interval.contig, []).append(g)
    wrap = set()
    for gs in by_contig.values():
        gs.sort(key=lambda g: (g.interval.start, g.interval.end))
        wrap.add((gs[-1].locus_tag, gs[0].locus_tag))
    return {(a.locus_tag, b.locus_tag) for a, b in pairs} - wrap


ROUND_TRIP = {
    "sim100k": (gippasim.SIM100K, gippasim.SIM100K_SEED),
    "origin_crossing": (
        gippasim.SimConfig(
            genome_length=40_000,
            n_plasmids=1,
            plasmid_length=(6_000, 9_000),
            origin_crossing_gene=True,
            n_sorfs=2,
            n_pseudogenes=1,
            linear_contig_length=5_000,
        ),
        7,
    ),
}


@pytest.mark.parametrize("name", sorted(ROUND_TRIP))
def test_example10_round_trip(tmp_path: Path, name: str) -> None:
    cfg, seed = ROUND_TRIP[name]
    genes, reps = _sim_genes(cfg, seed, tmp_path / "sim")
    gff = tmp_path / "uniop_cds.gff"
    id_map = write_prodigal_gff(gff, genes, reps)
    # every gene (CDS, sORF, pseudogene; Q1) appears exactly once
    assert sorted(t for _, t in id_map) == sorted(g.locus_tag for g in genes)
    assert len({p for p, _ in id_map}) == len(id_map)
    assert [gff_id(r) for r in gff_rows(gff)] == [p for p, _ in id_map]
    for pid, _ in id_map:
        contig, k = pid.rsplit("_", 1)
        assert contig in {r.id for r in reps} and int(k) >= 1
    pred = uniop_like_pred(gff, tmp_path / "uniop.pred", 0.5)  # UniOP never wraps
    got = parse_uniop_pred(pred, gff, dict(id_map))
    assert all(p.source == "uniop" and p.uniop_p == 0.5 for p in got)
    assert {(p.gene_a, p.gene_b) for p in got} == _non_wrap(genes, reps)
    assert len(got) == len(_non_wrap(genes, reps))


# ---------------------------------------------------------------- command builder


def test_uniop_command() -> None:
    home = Path("/opt/uniop")
    cmd = uniop_command(Path("work/uniop_cds.gff"), Path("work/out"), home)
    assert all(isinstance(a, str) for a in cmd)
    assert cmd[0] in {"python", "python3"}
    assert Path(cmd[1]).parent == home / "src"
    assert Path(cmd[1]).name in {"UniOP", "UniOP.py"}
    rest = cmd[2:]
    assert len(rest) == 6
    pairs = dict(zip(rest[::2], rest[1::2], strict=True))
    assert pairs == {
        "-f": "work/uniop_cds.gff",
        "--gff_source": "prodigal_gff",
        "-t": "work/out",
    }


def test_uniop_modules_never_import_uniop() -> None:
    import ast

    import gippa.annotate.operons as ops
    import gippa.annotate.uniop_adapter as uni
    import gippa.annotate.uniop_cmd as ucmd

    for mod in (uni, ucmd, ops):
        tree = ast.parse(Path(mod.__file__).read_text())  # type: ignore[arg-type]
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module]
            assert not [n for n in names if n.split(".")[0].lower() == "uniop"], mod.__name__
