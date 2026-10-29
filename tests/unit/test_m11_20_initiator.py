"""M11-20 initiator M at annotated non-ATG starts, exact route (docs/specs/
M11-20-initiator-methionine.md § Definition rules 1-7, § Examples 1-13, AC 1-9).

Interface pinned here (spec § Definition)::

    gippa.mapping.initiator.StartSite(contig, strand, codon_start, locus_tags, codon)
    gippa.mapping.initiator.StartCells(sites, target, decoys)
    annotated_starts(gm, genome, start_codons, *, table) -> (sites, {atg, not_a_start, aa_not_m,
                                                                       off_contig})
    start_frame_cells(index, sites) -> tuple[np.ndarray, ...]
    mirror_start_cells(kind, frames, cells) -> tuple[np.ndarray, ...]
    initiator_rows(rows, sites, lengths, circular) -> list[bool]
    match_exact(..., starts=None)            annotate_peptide_hits(..., initiator=<initiator_rows>)
    align_hit(..., starts=frozenset()) / align_hits(..., starts=frozenset())   OpKind "initiator"
    save_index(path, index, decoys, cells=None, starts=None) -> format 4
    load_start_cells(path, *, table) -> StartCells

Readings taken (report): ``annotate_peptide_hits(initiator=...)`` takes the per-row booleans of
``initiator_rows`` (as ``consumed=`` takes ``consumed_sites``); ``align_hit`` takes ``starts=`` as
``align_hits`` does; a pseudogene is not counted in any skip counter (as ``annotated_stops``);
the ``mapping.start_sites`` event's ``sites`` is a count; example 7's origin copy at p = 59 (it fits
the 180-nt extended frame, rule 2) is a cell whose window never fits.
"""

from __future__ import annotations

import importlib
import io
import json
import logging
from dataclasses import asdict, replace
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.initiator_toys import (
    STARTS,
    cds,
    row_set,
    toy_c1,
    toy_c2,
    toy_c3,
    toy_c4,
    toy_c5,
    with_decoys,
)
from helpers.readthrough_toys import model, write_annotation, write_genome

from gippa.mapping import loci
from gippa.mapping.exact import match_exact
from gippa.mapping.sixframe import decode, load_index, save_index
from gippa.model import Replicon

TABLES = (4, 11, 25)
ZERO = {"atg": 0, "not_a_start": 0, "aa_not_m": 0, "off_contig": 0}
FORMAT2 = Path(__file__).resolve().parents[1] / "data/mapping/readthrough/sixframe_format2.npz"
FORMAT1 = FORMAT2.with_name("sixframe_format1.npz")
WARN_EVENT = "mapping.index_without_start_cells"


def ini() -> ModuleType:
    """The module under test, imported per test so each test fails on its own."""
    return importlib.import_module("gippa.mapping.initiator")


def site(contig: str, strand: str, start: int, tags: tuple[str, ...], codon: str) -> Any:
    return ini().StartSite(contig, strand, start, tags, codon)


def starts_of(world: Any, table: int = 11, seed: int = 7) -> tuple[Any, Any, Any]:
    """(index, decoys, StartCells) as the step builds them."""
    genome, reps, gm = world
    m = ini()
    index, decoys = with_decoys(genome, reps, table, seed)
    sites, _ = m.annotated_starts(gm, genome, STARTS, table=table)
    target = m.start_frame_cells(index, sites)
    dc = {k: m.mirror_start_cells(k, fr, target) for k, fr in decoys.items()}
    return index, decoys, m.StartCells(sites, target, dc)


def exact(
    world: Any, queries: list[str], kind: str = "reverse", starts: Any = "built", table: int = 11
) -> pd.DataFrame:
    index, decoys, cells = starts_of(world, table)
    return match_exact(
        queries,
        index,
        decoys,
        dict.fromkeys(queries, kind),
        starts=cells if starts == "built" else starts,
    )


def lists(cells: tuple[np.ndarray, ...]) -> list[list[int]]:
    return [[int(x) for x in c] for c in cells]


def peptides(peps: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"peptide": peps, "search_db": "other", "protein_ids": [[] for _ in peps]})


def annotate(world: Any, hits: pd.DataFrame, relation: str = "related") -> pd.DataFrame:
    genome, reps, gm = world
    rows = hits[hits["search"] == "target"].reset_index(drop=True)
    sites, _ = ini().annotated_starts(gm, genome, STARTS, table=11)
    flags = ini().initiator_rows(
        rows,
        sites,
        {c: r.length for c, r in reps.items()},
        {c: r.topology == "circular" for c, r in reps.items()},
    )
    return loci.annotate_peptide_hits(
        rows, gm, genome, peptides(sorted(set(rows["peptide"]))), relation, initiator=flags
    )  # type: ignore[arg-type]


def flags_of(df: pd.DataFrame, start: int) -> set[str]:
    (f,) = df.loc[df["start"] == start, "flags"].tolist()
    return {x for x in str(f).split(";") if x and x != "None"}


# ---------------------------------------------------------------- AC 1: sites and skips


def test_example_1_site_and_no_skips() -> None:
    genome, _reps, gm = toy_c1()
    sites, skipped = ini().annotated_starts(gm, genome, STARTS, table=11)
    assert sites == (site("c1", "+", 12, ("G1",), "GTG"),)
    assert dict(skipped) == ZERO


def test_example_3_atg_start_is_skipped() -> None:
    genome, _reps, gm = toy_c1("ATG")
    assert ini().annotated_starts(gm, genome, STARTS, table=11) == ((), {**ZERO, "atg": 1})


def test_examples_6_7_7b_sites() -> None:
    m = ini()
    for world, want in [
        (toy_c2(), site("c2", "-", 27, ("G2",), "TTG")),
        (toy_c3(), site("c3", "+", 87, ("G3",), "TTG")),
        (toy_c4(), site("c4", "+", 89, ("G4",), "GTG")),
    ]:
        genome, _reps, gm = world
        assert m.annotated_starts(gm, genome, STARTS, table=11) == ((want,), ZERO)


def test_example_13_no_genes_and_only_atg() -> None:
    genome, reps, gm = toy_c1("ATG")
    m = ini()
    assert m.annotated_starts(model(reps, []), genome, STARTS, table=11) == ((), ZERO)
    assert m.annotated_starts(gm, genome, STARTS, table=11)[0] == ()


def test_skip_counters_each_rule() -> None:
    """Each feature fails exactly one rule: aa not M, codon not a start (CCC; and N), off the
    linear contig; a pseudogene with a GTG start is no site and not counted; sorf counts."""
    nt = list("C" * 120)
    nt[0:3], nt[30:36], nt[60:63], nt[90:93] = "GTG", "CCCGTG", "GNG", "TTG"
    genome = {"L": "".join(nt)}
    reps = {"L": Replicon("L", 120, "linear")}
    gm = model(
        reps,
        [
            cds("A", "L", 0, 27, "+", "VPPPPPPP", "GTG"),  # aa_not_m (5'-partial, Edge)
            cds("B", "L", 30, 57, "+", "MPPPPPPP", "CCC"),  # not_a_start
            cds("C", "L", 60, 87, "+", "MPPPPPPP", "GNG"),  # not_a_start (non-ACGT)
            cds("D", "L", 90, 117, "+", "MPPPPPPP", "TTG", "sorf"),  # site
            cds("E", "L", 102, 123, "-", "MPPPPPP", "GGG"),  # first codon [120, 123): off_contig
            cds("P", "L", 0, 27, "+", "MPPPPPPP", "GTG", "pseudogene"),
            cds("Z", "L", 33, 60, "+", "", None),  # empty aa (GTG at 33): aa_not_m
        ],
    )
    sites, skipped = ini().annotated_starts(gm, genome, STARTS, table=11)
    assert sites == (site("L", "+", 90, ("D",), "TTG"),)
    assert dict(skipped) == {"atg": 0, "not_a_start": 2, "aa_not_m": 2, "off_contig": 1}


def test_two_features_share_a_codon_and_codons_are_upper_cased() -> None:
    genome, reps, gm = toy_c1()
    gm2 = model(reps, [*gm.genes, cds("G0", "c1", 12, 39, "+", "MAELWHKG", "GTG", "sorf")])
    sites, _ = ini().annotated_starts(gm2, genome, ["atg", "gtg", "ttg"], table=11)
    assert sites == (site("c1", "+", 12, ("G0", "G1"), "GTG"),)


def test_user_set_without_gtg_gives_no_site() -> None:
    genome, _reps, gm = toy_c1()
    assert ini().annotated_starts(gm, genome, ["ATG", "TTG"], table=11) == (
        (),
        {**ZERO, "not_a_start": 1},
    )


# ---------------------------------------------------------------- AC 2: cells and mirrors


@pytest.mark.parametrize(
    ("world", "want"),
    [
        (toy_c1, [[4, 44], [], [], [], [], []]),
        (toy_c2, [[], [], [], [10], [], []]),
        (toy_c3, [[29, 59], [], [], [], [], []]),
        (toy_c4, [[], [], [29], [], [], []]),
    ],
)
def test_examples_cells(world: Any, want: list[list[int]]) -> None:
    index, _decoys, cells = starts_of(world())
    assert lists(cells.target) == want
    assert all(c.dtype == np.int64 for c in cells.target)
    n = [len(f.aa) for f in index.frames]
    assert lists(cells.decoys["reverse"]) == [
        sorted(x - 1 - p for p in w) for x, w in zip(n, want, strict=True)
    ]
    assert lists(cells.decoys["shuffle"]) == want


def test_example_8_reverse_cells_and_frame() -> None:
    _index, decoys, cells = starts_of(toy_c1())
    assert decode(decoys["reverse"][0].aa) == (
        "P" * 27 + "*GKHWLEAV" + "P" * 31 + "*GKHWLEAV" + "P" * 4
    )
    assert lists(cells.decoys["reverse"])[0] == [35, 75]


def test_start_frame_cells_raises_on_a_wrong_residue() -> None:
    genome, reps, _gm = toy_c1()
    index, _ = with_decoys(genome, reps)
    with pytest.raises(ValueError):
        ini().start_frame_cells(index, [site("c1", "+", 12, ("G1",), "TTG")])  # frame reads V
    with pytest.raises(ValueError):
        ini().start_frame_cells(index, [site("c1", "+", 15, ("G1",), "GTG")])  # reads A


def test_mirror_draws_no_randomness(monkeypatch: pytest.MonkeyPatch) -> None:
    _index, decoys, cells = starts_of(toy_c1())

    def boom(*_a: object, **_k: object) -> None:
        raise AssertionError("randomness drawn")

    monkeypatch.setattr(np.random, "default_rng", boom)
    for kind in ("reverse", "shuffle"):
        got = ini().mirror_start_cells(kind, decoys[kind], cells.target)
        assert lists(got) == lists(cells.decoys[kind])


# ---------------------------------------------------------------- AC 3: credit rule


def test_example_1_one_credited_target_row() -> None:
    df = exact(toy_c1(), ["MAELWHK"])
    tgt = df[df["search"] == "target"]
    assert len(tgt) == 1
    r = tgt.iloc[0]
    assert (r.frame_id, r.aa_start, r.aa_end, r.start, r.end, r.strand, r.frame) == (
        0,
        4,
        11,
        12,
        33,
        "+",
        1,
    )
    assert (r.route, r.n_diff, r.diff_types) == ("exact", 0, "")
    assert r.blocks is None or pd.isna(r.blocks)
    assert exact(toy_c1(), ["MAELWHK"], starts=None).empty  # today: no hit


def test_example_2_credited_row_and_chance_copy() -> None:
    world = toy_c1(copy_at_60=True)
    df = exact(world, ["MAELWHK"])
    tgt = df[df["search"] == "target"]
    assert sorted(zip(tgt["start"], tgt["end"], strict=True)) == [(12, 33), (60, 81)]
    ann = annotate(world, df)
    by = {int(s): r for s, r in zip(ann["start"], ann.to_dict("records"), strict=True)}
    assert (by[12]["locus_class"], by[12]["gene"]) == ("cds_inframe", "G1")
    assert by[60]["locus_class"] == "intergenic"
    assert {by[12]["uniqueness"], by[60]["uniqueness"]} == {"shared"}
    assert flags_of(ann, 12) == {"initiator_m", "multi"} and flags_of(ann, 60) == {"multi"}
    from gippa.novel.candidates import candidate_rows

    pep = pd.DataFrame({"peptide": ["MAELWHK"], "is_contaminant": [False], "q_value": [0.001]})
    assert candidate_rows(ann, pep, None).empty


@pytest.mark.parametrize("query", ["PPMAELWHK", "MAEMWHK"])
def test_example_4_m_at_position_ge_1_is_never_credited(query: str) -> None:
    assert exact(toy_c1(), [query]).empty


def test_example_5_v_initial_query_matches_as_today() -> None:
    world = toy_c1()
    new, old = exact(world, ["PPVAELWHK"]), exact(world, ["PPVAELWHK"], starts=None)
    assert row_set(new) == row_set(old)
    assert {
        (s, e)
        for s, e in zip(
            new.loc[new["search"] == "target", "start"],
            new.loc[new["search"] == "target", "end"],
            strict=True,
        )
    } == {(6, 33)}
    assert "initiator_m" not in flags_of(annotate(world, new), 6)


def test_example_6_minus_strand() -> None:
    world = toy_c2()
    df = exact(world, ["MKQWYFR"])
    (r,) = df[df["search"] == "target"].itertuples(index=False)
    assert (r.frame_id, r.aa_start, r.aa_end, r.start, r.end, r.strand, r.frame) == (
        3,
        10,
        17,
        9,
        30,
        "-",
        -1,
    )
    ann = annotate(world, df)
    assert (ann["locus_class"].iloc[0], ann["gene"].iloc[0]) == ("cds_inframe", "G2")
    assert flags_of(ann, 9) == {"initiator_m"}


@pytest.mark.parametrize(("world", "start", "end"), [(toy_c3, 87, 108), (toy_c4, 89, 110)])
def test_examples_7_and_7b_origin(world: Any, start: int, end: int) -> None:
    w = world()
    df = exact(w, ["MAELWHK"])
    tgt = df[df["search"] == "target"]
    assert list(zip(tgt["start"], tgt["end"], tgt["strand"], strict=True)) == [(start, end, "+")]
    assert flags_of(annotate(w, df), start) == {"initiator_m", "origin_spanning"}


def test_example_8_reverse_decoy_one_row() -> None:
    df = exact(toy_c1(), ["MKHWLEA", "MPPPPPP"])
    assert df[df["search"] == "target"].empty
    dec = df[df["search"] == "decoy"]
    assert len(dec) == 1
    r = dec.iloc[0]
    assert (r.peptide, r.decoy_type, r.frame_id, r.start, r.end, r.strand) == (
        "MKHWLEA",
        "reverse",
        0,
        12,
        33,
        "+",
    )
    assert (int(r.aa_start), int(r.aa_end)) in {(69, 76), (29, 36)}


def test_example_9_tolerant_route_unchanged() -> None:
    """A GTG gene whose 12-aa N-terminal peptide carries one more substitution (residue 6)."""
    from gippa.config.model import TolerantConfig
    from gippa.mapping.tolerant import match_tolerant

    aa = "MAELWHGNPQEK"
    body = "GTG" + "".join(
        {
            "A": "GCT",
            "E": "GAA",
            "L": "CTG",
            "W": "TGG",
            "H": "CAT",
            "G": "GGT",
            "N": "AAT",
            "P": "CCG",
            "Q": "CAG",
            "K": "AAA",
        }[x]
        for x in aa[1:]
    )
    nt = "C" * 12 + body + "TAA" + "C" * 63
    reps = {"t": Replicon("t", len(nt), "linear")}
    world = ({"t": nt}, reps, model(reps, [cds("T1", "t", 12, 51, "+", aa, "GTG")]))
    query = aa[:6] + "Y" + aa[7:]
    assert exact(world, [query]).empty
    index, decoys, _cells = starts_of(world)
    tol = match_tolerant([query], index, decoys, {query: "reverse"}, {query}, TolerantConfig())
    tgt = tol[tol["search"] == "target"]
    assert list(zip(tgt["start"], tgt["end"], tgt["n_diff"], strict=True)) == [(12, 48, 2)]


@pytest.mark.parametrize("table", (4, 25))
def test_example_11_codes_4_and_25(table: int) -> None:
    world = toy_c5(table)
    q = {4: "MAEWLHK", 25: "MAEGLHK"}[table]
    df = exact(world, [q], table=table)
    tgt = df[df["search"] == "target"]
    assert list(zip(tgt["start"], tgt["end"], strict=True)) == [(6, 27)]


def test_example_12_shuffle_decoy() -> None:
    world = toy_c1()
    index, decoys, cells = starts_of(world)
    assert lists(cells.decoys["shuffle"]) == lists(cells.target)
    df = match_exact(["MAELWHK"], index, decoys, {"MAELWHK": "shuffle"}, starts=cells)
    tgt, dec = df[df["search"] == "target"], df[df["search"] == "decoy"]
    assert list(zip(tgt["start"], tgt["end"], strict=True)) == [(12, 33)]
    sh = decode(decoys["shuffle"][0].aa)
    spelled = any(sh[p + 1 : p + 7] == "AELWHK" for p in (4, 44))
    plain = "MAELWHK" in sh
    assert (len(dec) > 0) == (spelled or plain)


def test_example_13_and_ac3_no_sites_is_byte_identical() -> None:
    world = toy_c1("ATG")
    index, decoys, cells = starts_of(world)
    assert cells.sites == () and all(len(c) == 0 for c in cells.target)
    qs = ["MAELWHK", "PPPPMAELWHK", "MKHWLEA", "AELWHKG"]
    of = dict.fromkeys(qs, "reverse")

    def pq(df: pd.DataFrame) -> bytes:
        buf = io.BytesIO()
        df.to_parquet(buf, index=False)
        return buf.getvalue()

    base = pq(match_exact(qs, index, decoys, of))
    assert pq(match_exact(qs, index, decoys, of, starts=None)) == base
    assert pq(match_exact(qs, index, decoys, of, starts=cells)) == base
    empty = ini().StartCells((), cells.target, cells.decoys)
    assert pq(match_exact(qs, index, decoys, of, starts=empty)) == base


# ---------------------------------------------------------------- AC 4-5: flag and schema


def test_initiator_rows_only_credited_target_rows() -> None:
    genome, _reps, gm = toy_c1()
    sites, _ = ini().annotated_starts(gm, genome, STARTS, table=11)
    base = {
        "query": "",
        "decoy_type": None,
        "contig": "c1",
        "frame_id": 0,
        "aa_start": 0,
        "aa_end": 7,
        "frame": 1,
        "n_diff": 0,
        "diff_types": "",
        "blocks": None,
    }
    rows = pd.DataFrame(
        [
            {
                **base,
                "peptide": "MAELWHK",
                "search": "target",
                "route": "exact",
                "start": 12,
                "end": 33,
                "strand": "+",
            },
            {
                **base,
                "peptide": "MAELWHK",
                "search": "decoy",
                "decoy_type": "reverse",
                "route": "exact",
                "start": 12,
                "end": 33,
                "strand": "+",
            },
            {
                **base,
                "peptide": "MAELWHKG",
                "search": "target",
                "route": "tolerant",
                "start": 12,
                "end": 36,
                "strand": "+",
                "n_diff": 1,
            },
            {
                **base,
                "peptide": "VAELWHK",
                "search": "target",
                "route": "exact",
                "start": 12,
                "end": 33,
                "strand": "+",
            },
            {
                **base,
                "peptide": "MAELWHK",
                "search": "target",
                "route": "exact",
                "start": 15,
                "end": 36,
                "strand": "+",
            },
            {
                **base,
                "peptide": "MAELWHK",
                "search": "target",
                "route": "exact",
                "start": 12,
                "end": 33,
                "strand": "-",
            },
        ]
    )
    got = ini().initiator_rows(rows, sites, {"c1": 120}, {"c1": True})
    assert list(got) == [True, False, False, False, False, False]


def test_initiator_rows_minus_strand_5prime_codon() -> None:
    genome, _reps, gm = toy_c2()
    sites, _ = ini().annotated_starts(gm, genome, STARTS, table=11)
    rows = pd.DataFrame(
        {
            "peptide": ["MKQWYFR", "MKQWYFR"],
            "search": "target",
            "route": "exact",
            "contig": "c2",
            "start": [9, 6],
            "end": [30, 27],
            "strand": "-",
        }
    )
    assert list(ini().initiator_rows(rows, sites, {"c2": 60}, {"c2": False})) == [True, False]


def test_atg_start_row_has_no_flag_and_paralogues_are_multi_identical() -> None:
    genome, reps, gm = toy_c1("ATG")
    ann = annotate((genome, reps, gm), exact((genome, reps, gm), ["MAELWHK"]))
    assert "initiator_m" not in flags_of(ann, 12)
    nt = list(toy_c1()[0]["c1"])
    nt[60:87] = "GTGGCTGAACTGTGGCATAAAGGCTAA"
    world = (
        {"c1": "".join(nt)},
        reps,
        model(
            reps,
            [
                cds("G1", "c1", 12, 39, "+", "MAELWHKG", "GTG"),
                cds("G2", "c1", 60, 87, "+", "MAELWHKG", "GTG"),
            ],
        ),
    )
    ann = annotate(world, exact(world, ["MAELWHK"]))
    assert sorted(ann["start"]) == [12, 60]
    assert set(ann["uniqueness"]) == {"multi_identical"}
    assert all("initiator_m" in flags_of(ann, s) for s in (12, 60))


def test_peptide_hits_schema_accepts_initiator_m_and_version_bumps() -> None:
    import re

    from gippa import schemas

    doc = json.loads(
        (Path(schemas.__file__).parent / "tables/peptide_hits.schema.json").read_text()
    )
    (pattern,) = [f["constraints"]["pattern"] for f in doc["fields"] if f["name"] == "flags"]
    assert re.fullmatch(pattern, "initiator_m;origin_spanning")
    assert re.fullmatch(pattern, "initiator_m")
    major, minor = (int(x) for x in schemas.OUTPUT_SCHEMA_VERSION.split("."))
    assert (major, minor) == (1, 2)  # one additive minor bump over 1.1 (M11-16)
    pkg = json.loads((Path(schemas.__file__).parent / "datapackage.json").read_text())
    assert pkg["output_schema_version"] == schemas.OUTPUT_SCHEMA_VERSION


# ---------------------------------------------------------------- AC 6: index format 4


@pytest.mark.parametrize("table", TABLES)
def test_format_4_round_trips(tmp_path: Path, table: int) -> None:
    from gippa.mapping.sixframe import load_start_cells, load_stop_cells

    index, decoys, cells = starts_of(toy_c1(), table)
    path = tmp_path / "six.npz"
    save_index(path, index, decoys, None, starts=cells)
    with np.load(path) as npz:
        meta = json.loads(str(npz["meta"]))
        keys = set(npz.keys())
    assert meta["format"] == 4 and meta["translation_table"] == table
    assert meta["start_sites"] == [
        asdict(s) | {"locus_tags": list(s.locus_tags)} for s in cells.sites
    ]
    assert {
        "target_init",
        "target_init_bounds",
        "reverse_init",
        "reverse_init_bounds",
        "shuffle_init",
        "shuffle_init_bounds",
    } <= keys
    got = load_start_cells(path, table=table)
    assert got.sites == cells.sites
    assert lists(got.target) == lists(cells.target)
    for k in ("reverse", "shuffle"):
        assert lists(got.decoys[k]) == lists(cells.decoys[k])
    assert load_stop_cells(path, table=table).sites == ()
    assert load_index(path, table=table)[0].translation_table == table


def test_format_4_without_starts_writes_empty_lists(tmp_path: Path) -> None:
    from gippa.mapping.sixframe import load_start_cells

    index, decoys, _ = starts_of(toy_c1())
    path = tmp_path / "six.npz"
    save_index(path, index, decoys)
    with np.load(path) as npz:
        meta = json.loads(str(npz["meta"]))
    assert meta["format"] == 4 and meta["start_sites"] == []
    got = load_start_cells(path, table=11)
    assert got.sites == () and all(len(c) == 0 for c in got.target)


def _as_format_3(src: Path, dst: Path) -> Path:
    with np.load(src) as npz:
        arrays = {k: npz[k] for k in npz.files if "_init" not in k}
    meta = json.loads(str(arrays["meta"]))
    meta["format"] = 3
    meta.pop("start_sites", None)
    arrays["meta"] = np.array(json.dumps(meta))
    with dst.open("wb") as fh:
        np.savez_compressed(fh, **arrays)  # type: ignore[arg-type]
    return dst


@pytest.mark.parametrize("which", ["format3", "format2"])
def test_old_formats_load_with_no_start_cells_and_a_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, which: str
) -> None:
    from gippa.mapping.sixframe import load_start_cells

    if which == "format3":
        index, decoys, cells = starts_of(toy_c1())
        save_index(tmp_path / "new.npz", index, decoys, None, starts=cells)
        path = _as_format_3(tmp_path / "new.npz", tmp_path / "old.npz")
    else:
        path = FORMAT2
    caplog.set_level(logging.INFO)
    got = load_start_cells(path, table=11)
    assert got.sites == ()
    assert sum(len(c) for c in got.target) == 0
    assert all(len(c) == 0 for d in got.decoys.values() for c in d)
    ev = [r for r in caplog.records if getattr(r, "event", None) == WARN_EVENT]
    assert len(ev) == 1 and ev[0].levelno == logging.WARNING
    assert "rebuild" in ev[0].getMessage()


def test_format_1_is_refused() -> None:
    from gippa.mapping.sixframe import load_start_cells

    with pytest.raises(ValueError, match="rebuild the six-frame index"):
        load_start_cells(FORMAT1, table=11)


# ---------------------------------------------------------------- AC 7, 9: the steps


def _step_inputs(tmp: Path, world: Any, peps: list[str]) -> dict[str, Path]:
    genome, _reps, gm = world
    ann = write_annotation(gm, tmp / "ann")
    pq = tmp / "peptides.parquet"
    pd.DataFrame(
        {
            "peptide": [*peps, "SAMPLEWLDK"],
            "source": [*["database"] * len(peps), "denovo"],
            "q_value": [0.001] * len(peps) + [np.nan],
            "pep": 0.01,
            "score": 0.9,
            "is_decoy": False,
            "search_db": "other",
        }
    ).to_parquet(pq, index=False)
    return {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": ann["replicons"],
        "genes": ann["genes"],
        "peptides": pq,
    }


def _run(tmp: Path, ins: dict[str, Path]) -> tuple[Path, Path]:
    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    cfg = resolve_config(overrides=["seed=7"])
    six, hx = tmp / "sixframe.npz", tmp / "hits_exact.parquet"
    reg.run_step(
        "sixframe_index",
        cfg,
        {k: ins[k] for k in reg.get("sixframe_index").inputs},
        {"sixframe": six},
    )
    reg.run_step(
        "map_exact", cfg, {"peptides": ins["peptides"], "sixframe": six}, {"hits_exact": hx}
    )
    return six, hx


def test_steps_event_format_4_and_map_exact_credit(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.INFO)
    six, hx = _run(tmp_path, _step_inputs(tmp_path, toy_c1(), ["MAELWHK"]))
    with np.load(six) as npz:
        assert json.loads(str(npz["meta"]))["format"] == 4
    ev = [r for r in caplog.records if getattr(r, "event", None) == "mapping.start_sites"]
    assert len(ev) == 1 and ev[0].levelno == logging.INFO
    data = ev[0].data  # type: ignore[attr-defined]
    want = {
        "sites": 1,
        "cells_target": 2,
        "cells_decoy": {"reverse": 2, "shuffle": 2},
        "skipped": ZERO,
        "translation_table": 11,
    }
    assert {k: data[k] for k in want} == want
    assert set(data) - set(want) <= {"usable_k7"}
    hits = pd.read_parquet(hx)
    tgt = hits[(hits["search"] == "target") & (hits["peptide"] == "MAELWHK")]
    assert list(zip(tgt["start"], tgt["end"], strict=True)) == [(12, 33)]


def test_ac9_atg_only_gippasim_hits_exact_unchanged(tmp_path: Path) -> None:
    """AC 9 at the mapping steps: ATG-only starts give no site and the ``starts=None`` table."""
    import gippasim
    from gippa.annotate.bakta import parse_bakta_json
    from gippa.mapping.sixframe import load_start_cells
    from gippa.steps.mapping import decoy_types, select_queries

    res = gippasim.simulate(
        gippasim.SimConfig(genome_length=20_000, start_codons=(("ATG", 1.0),)),
        seed=20,
        outdir=tmp_path / "sim",
    )
    ann = parse_bakta_json(res.bakta_json)
    gm = ann.model
    world = (ann.sequences, dict(gm.replicons), gm)
    peps = sorted({g.aa[:9] for g in gm.genes if len(g.aa) >= 9} | {"MAELWHK"})
    six, hx = _run(tmp_path, _step_inputs(tmp_path, world, peps))
    assert load_start_cells(six, table=11).sites == ()
    from gippa.config import resolve_config

    cfg = resolve_config(overrides=["seed=7"])
    q = select_queries(pd.read_parquet(tmp_path / "peptides.parquet"), cfg)
    index, decoys = load_index(six, table=11)
    old = match_exact(q["peptide"].tolist(), index, decoys, decoy_types(q, cfg))
    buf = io.BytesIO()
    old.to_parquet(buf, index=False)
    assert pd.read_parquet(hx).equals(pd.read_parquet(io.BytesIO(buf.getvalue())))


# ---------------------------------------------------------------- AC 8: variants


def _ph(rows: list[tuple[str, int, int]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "peptide": p,
                "contig": "c1",
                "start": s,
                "end": e,
                "strand": "+",
                "n_diff": 0,
                "diff_types": "",
            }
            for p, s, e in rows
        ]
    )


def test_example_10_align_and_pileup() -> None:
    from gippa.config.model import GippaConfig
    from gippa.variants import align, pileup

    genome, reps, gm = toy_c1()
    keys = frozenset({("c1", "+", 12)})
    aligns = align.align_hits(
        _ph([("MAELWHK", 12, 33), ("AELWHKG", 15, 36)]), None, genome, reps, starts=keys
    )
    first = aligns[0].ops[0]
    assert (first.kind, first.pep, first.codons, first.ref, first.alt) == (
        "initiator",
        (0, 1),
        (12,),
        "V",
        "M",
    )
    assert {o.kind for o in aligns[0].ops[1:]} == {"match"}
    cfg = GippaConfig.model_validate({"genome_relation": "same"})
    sup = {
        p: pileup.Support(p, "database", "same", 0.001, "unique") for p in ("MAELWHK", "AELWHKG")
    }
    pile = pileup.pileup(aligns, gm, genome, cfg, support=sup)
    assert pile.attrs["initiator"] == 1
    at12 = pile[pile["codon_start"] == 12]
    assert at12.empty or (int(at12["n_ref"].iloc[0]), int(at12["n_il"].iloc[0])) == (0, 0)
    assert at12.empty or at12["call"].iloc[0] == "no_call"
    for c in range(15, 33, 3):
        (call,) = pile.loc[pile["codon_start"] == c, "call"].tolist()
        assert call == "REF", c
    assert 12 not in set(pileup.sites(pile, cfg)["codon_start"])


def test_align_hit_raises_off_site_and_on_other_mismatches() -> None:
    from gippa.variants import align

    genome, reps, _ = toy_c1()
    rep = reps["c1"]
    keys = frozenset({("c1", "+", 12)})
    with pytest.raises(ValueError, match="not exact"):
        align.align_hit("MAELWHK", "c1", 12, 33, "+", 0, None, genome["c1"], rep)
    with pytest.raises(ValueError, match="not exact"):
        align.align_hit("MAEMWHK", "c1", 12, 33, "+", 0, None, genome["c1"], rep, starts=keys)
    got = align.align_hit("MAELWHK", "c1", 12, 33, "+", 0, None, genome["c1"], rep, starts=keys)
    assert got.ops[0].kind == "initiator"
    plain = align.align_hit("VAELWHK", "c1", 12, 33, "+", 0, None, genome["c1"], rep, starts=keys)
    assert {o.kind for o in plain.ops} == {"match"}
    assert "initiator" in align.OpKind.__args__  # type: ignore[attr-defined]


def test_tolerant_rows_keep_the_alt_m() -> None:
    """Edge case "tolerant hits": a tolerant row is never given the initiator op."""
    from gippa.variants import align

    genome, reps, _ = toy_c1()
    keys = frozenset({("c1", "+", 12)})
    got = align.align_hit(
        "MAELWHKG", "c1", 12, 36, "+", 1, None, genome["c1"], reps["c1"], starts=keys
    )
    assert got.ops[0].kind == "substitution" and got.ops[0].alt == "M"


def test_replace_keeps_site_dataclass_frozen() -> None:
    s = site("c1", "+", 12, ("G1",), "GTG")
    with pytest.raises(AttributeError):
        s.codon_start = 13  # type: ignore[misc]
    assert replace(s, codon="TTG").codon == "TTG"
