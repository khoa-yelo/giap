"""UniOP adapter: Prodigal-style GFF, ``uniop.pred`` parser, command (M2-12; GIAP_SPEC.md §7.1.1,
§7.1.3; Module 1 test "adapter round-trip for UniOP (Bakta → Prodigal-style GFF → back to locus
tags)").

Examples 4 and 10 of docs/specs/M2-12-uniop-pairs-tus.md. UniOP has no licence (ADR-0008) and is
never imported or run here.

Conventions pinned here (spec Q3 and silences):

* Unit-test ``uniop.pred`` files are the doc-derived 3-column table of Q3: a tab-separated header
  ``gene_a  gene_b  p`` and one row per pair (the contract test covers the recorded header).
* ``uniop_map.tsv`` (step ``uniop_input``) has the header ``prodigal_id  locus_tag``.
* The GFF starts with a ``##gff-version`` line; each gene is one ``CDS`` line with 1-based closed
  coordinates and the attribute column starting with ``ID=<contig>_<k>``.
* A ``uniop.pred`` row naming a Prodigal ID absent from the map raises ``ValueError``.
* ``uniop_command`` runs ``python <uniop_home>/src/UniOP`` (entry point ``UniOP`` or ``UniOP.py``,
  "confirm the entry-point name on the pinned commit").
"""

from __future__ import annotations

import itertools
import re
from pathlib import Path

import pytest

import giapsim
from giap.annotate.bakta import parse_bakta_json
from giap.annotate.operons import adjacent_pairs
from giap.annotate.uniop import PairProb, parse_uniop_pred, write_prodigal_gff
from giap.annotate.uniop_cmd import uniop_command
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

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


def write_pred(path: Path, rows: list[tuple[str, str, float]]) -> Path:
    lines = ["gene_a\tgene_b\tp", *(f"{a}\t{b}\t{p}" for a, b, p in rows)]
    path.write_text("\n".join(lines) + "\n")
    return path


def gff_rows(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


def gff_id(row: list[str]) -> str:
    m = re.match(r"ID=([^;]+)", row[8])
    assert m, f"attributes must start with ID=: {row[8]!r}"
    return m.group(1)


# ---------------------------------------------------------------- example 4


def test_example4_pred_row_maps_to_locus_tags(tmp_path: Path) -> None:
    pred = write_pred(tmp_path / "uniop.pred", [("chr1_1", "chr1_2", 0.91)])
    got = parse_uniop_pred(pred, {"chr1_1": "GSIM_00005", "chr1_2": "GSIM_00010"})
    assert got == [PairProb("GSIM_00005", "GSIM_00010", 0.91, "uniop")]


def test_pred_with_several_rows_keeps_file_order_and_values(tmp_path: Path) -> None:
    pred = write_pred(
        tmp_path / "uniop.pred",
        [("chr1_1", "chr1_2", 0.91), ("chr1_3", "chr1_4", 0.0), ("p2_1", "p2_2", 1.0)],
    )
    id_map = {"chr1_1": "A", "chr1_2": "B", "chr1_3": "C", "chr1_4": "D", "p2_1": "E", "p2_2": "F"}
    got = parse_uniop_pred(pred, id_map)
    assert [(p.gene_a, p.gene_b, p.uniop_p, p.source) for p in got] == [
        ("A", "B", 0.91, "uniop"),
        ("C", "D", 0.0, "uniop"),
        ("E", "F", 1.0, "uniop"),
    ]


def test_pred_unknown_prodigal_id_raises(tmp_path: Path) -> None:
    pred = write_pred(tmp_path / "uniop.pred", [("chr1_1", "chr1_9", 0.5)])
    with pytest.raises(ValueError, match="chr1_9"):
        parse_uniop_pred(pred, {"chr1_1": "A", "chr1_2": "B"})


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


def _sim_genes(cfg: giapsim.SimConfig, seed: int, out: Path) -> tuple[list[Gene], list[Replicon]]:
    res = giapsim.simulate(cfg, seed, out)
    model = parse_bakta_json(res.bakta_json).model
    return list(model.genes), list(model.replicons.values())


def _pred_from_gff(gff: Path) -> list[tuple[str, str, float]]:
    """All adjacent same-strand pairs of consecutive GFF lines per contig (UniOP never wraps)."""
    rows = gff_rows(gff)
    out = []
    for r1, r2 in itertools.pairwise(rows):
        if r1[0] == r2[0] and r1[6] == r2[6]:
            out.append((gff_id(r1), gff_id(r2), 0.5))
    return out


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
    "sim100k": (giapsim.SIM100K, giapsim.SIM100K_SEED),
    "origin_crossing": (
        giapsim.SimConfig(
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
    pred = write_pred(tmp_path / "uniop.pred", _pred_from_gff(gff))
    got = parse_uniop_pred(pred, dict(id_map))
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

    import giap.annotate.operons as ops
    import giap.annotate.uniop as uni
    import giap.annotate.uniop_cmd as ucmd

    for mod in (uni, ucmd, ops):
        tree = ast.parse(Path(mod.__file__).read_text())  # type: ignore[arg-type]
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module]
            assert not [n for n in names if n.split(".")[0].lower() == "uniop"], mod.__name__
