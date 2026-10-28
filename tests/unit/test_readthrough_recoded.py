"""Annotated recoded Sec/Pyl codons as stop sites (M7-16c; ADR-0058 item 6; spec "16c amendment
(2026-10-25)" rules 21-22, examples 21-28, edge cases; test plan C1-C5).

Interface pinned here (rule 22)::

    annotated_stops(gm, genome, recoded: bool = False) -> (sites, n_not_a_stop)   # unchanged type
    recoded_skips(gm, genome) -> {"length_mismatch": n, "codon_mismatch": n}       # pure
    event mapping.stop_sites data gains "recoded_skipped": recoded_skips(...)       # Q21 default

Rule 21 in brief: every CDS/sORF whose terminal codon is a stop (never a pseudogene, never a
``not_a_stop`` feature, Q24) with ``len(aa) == (e - s) / 3 - 1`` is scanned for upper-case ``U``
/ ``O`` at residue ``i >= 1``; residue ``i`` is the codon ``[s + 3i, s + 3i + 3)`` on ``+`` and
``[e - 3i - 3, e - 3i)`` on ``-`` (lowest coordinate wrapped mod n). ``U`` must sit over TGA and
``O`` over TAG, else skip and count ``codon_mismatch`` (never raise). Sites are merged by
``(contig, strand, codon_start)`` across all features (kind ``recoded`` wins, Q22) and sorted with
the terminal ones. ``recoded=False`` is exactly the 16a/16b output.

Toy contig ``c`` (examples 21-28): filler codons ``GCT`` around the named codons.
"""

from __future__ import annotations

import logging
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from helpers.readthrough_index import peptides_table
from helpers.readthrough_index import world as rt_world
from helpers.readthrough_toys import gene, model, rc, toy_c1, toy_c2, write_annotation, write_genome

from gippa.mapping.sixframe import STOP, build_sixframe
from gippa.model import Replicon
from gippa.model.genes import Gene, GeneModel

EVENT = "mapping.stop_sites"
FILL = "GCT"


def _rt():  # type: ignore[no-untyped-def]
    import gippa.mapping.readthrough as rt

    return rt


def g(
    tag: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    aa: str,
    ftype: str = "cds",
    partial: bool = False,
) -> Gene:
    return replace(gene(tag, contig, start, end, strand, ftype, partial=partial), aa=aa)


def site(contig: str, strand: str, start: int, tags: tuple[str, ...], codon: str, kind: str):  # type: ignore[no-untyped-def]
    return _rt().StopSite(contig, strand, start, tags, codon, kind)


def stops(gm: GeneModel, genome: dict[str, str], recoded: bool | None = None):  # type: ignore[no-untyped-def]
    rt = _rt()
    if recoded is None:
        return rt.annotated_stops(gm, genome)
    return rt.annotated_stops(gm, genome, recoded=recoded)


def skips(gm: GeneModel, genome: dict[str, str]) -> dict[str, int]:
    return dict(_rt().recoded_skips(gm, genome))


def put(seq: list[str], at: int, nt: str) -> None:
    n = len(seq)
    for i, c in enumerate(nt):
        seq[(at + i) % n] = c


def filler(n: int) -> list[str]:
    return list((FILL * (n // 3 + 1))[:n])


def lists(cells: tuple[np.ndarray, ...]) -> list[list[int]]:
    return [[int(x) for x in c] for c in cells]


ZERO = {"length_mismatch": 0, "codon_mismatch": 0}

# ---------------------------------------------------------------- the toys of examples 21-28

EX21_GENE = "ATGGCTGCTGCTTGAGCTGCTGCTGCTTAA"  # [3, 33), aa MAAAUAAAA


def ex21(aa: str = "MAAAUAAAA", codon4: str = "TGA") -> tuple[dict[str, str], GeneModel]:
    """Example 21: ``c`` linear 60 nt, ``G1`` cds ``+`` ``[3, 33)``, residue 4 at ``[15, 18)``."""
    seq = filler(60)
    put(seq, 3, EX21_GENE)
    put(seq, 15, codon4)
    reps = {"c": Replicon("c", 60, "linear")}
    return {"c": "".join(seq)}, model(reps, [g("G1", "c", 3, 33, "+", aa)])


def ex22() -> tuple[dict[str, str], GeneModel]:
    """Example 22: ``G2`` cds ``-`` ``[3, 33)``, the same gene-strand bases as example 21."""
    seq = filler(60)
    put(seq, 3, rc(EX21_GENE))
    reps = {"c": Replicon("c", 60, "linear")}
    return {"c": "".join(seq)}, model(reps, [g("G2", "c", 3, 33, "-", "MAAAUAAAA")])


def circular_gene(strand: str, residue: int, start: int = 270, n_codons: int = 20):  # type: ignore[no-untyped-def]
    """Circular ``c`` (300 nt, ``ext_nt`` 426) with one cds over unwrapped ``[start, start +
    3 * n_codons)``: ATG, sense codons, TGA at ``residue``, terminal TAA (gene strand)."""
    codons = ["ATG"] + [FILL] * (n_codons - 2) + ["TAA"]
    codons[residue] = "TGA"
    gnt = "".join(codons)
    seq = filler(300)
    put(seq, start, gnt if strand == "+" else rc(gnt))
    aa = "M" + "A" * (n_codons - 2)
    aa = aa[:residue] + "U" + aa[residue + 1 :]
    reps = {"c": Replicon("c", 300, "circular")}
    tag = "G3" if strand == "+" else "G4"
    gm = model(reps, [g(tag, "c", start, start + 3 * n_codons, strand, aa)])
    return {"c": "".join(seq)}, gm


def layout(
    parts: list[tuple[str, str, str, list[str], str]], contig: str = "c", gap: int = 3
) -> tuple[dict[str, str], GeneModel]:
    """Genes laid out left to right on a linear contig: (tag, strand, feature_type, gene-strand
    codons, aa); ``gap`` filler codons before each gene and after the last."""
    seq: list[str] = []
    genes = []
    for tag, strand, ftype, codons, aa in parts:
        seq += list(FILL * gap)
        nt = "".join(codons)
        s = len(seq)
        seq += list(nt if strand == "+" else rc(nt))
        genes.append(g(tag, contig, s, s + len(nt), strand, aa, ftype))
    seq += list(FILL * gap)
    reps = {contig: Replicon(contig, len(seq), "linear")}
    return {contig: "".join(seq)}, model(reps, genes)


# ---------------------------------------------------------------- C1: examples 21-24


def test_example_21_plus_linear() -> None:
    genome, gm = ex21()
    assert genome["c"][3:33] == EX21_GENE and genome["c"][15:18] == "TGA"
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 15, ("G1",), "TGA", "recoded"),
        site("c", "+", 30, ("G1",), "TAA", "terminal"),
    )
    assert n_not == 0
    assert skips(gm, genome) == ZERO
    assert stops(gm, genome, False) == ((site("c", "+", 30, ("G1",), "TAA", "terminal"),), 0)


def test_example_21_cells_round_trip() -> None:
    genome, gm = ex21()
    index = build_sixframe(genome, dict(gm.replicons))
    sites, _ = stops(gm, genome, True)
    cells = _rt().frame_cells(index, sites)
    assert lists(cells) == [[5, 10], [], [], [], [], []]
    iv = index.to_interval(0, 5, 6)
    assert (iv.contig, iv.start, iv.end, iv.strand) == ("c", 15, 18, "+")
    assert index.frames[0].aa[5] == STOP


def test_example_22_minus_linear() -> None:
    genome, gm = ex22()
    assert genome["c"][30:33] == "CAT" and genome["c"][18:21] == "TCA"
    assert genome["c"][3:6] == "TTA"
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "-", 3, ("G2",), "TAA", "terminal"),
        site("c", "-", 18, ("G2",), "TGA", "recoded"),
    )
    assert n_not == 0
    assert skips(gm, genome) == ZERO
    index = build_sixframe(genome, dict(gm.replicons))
    cells = _rt().frame_cells(index, sites)
    assert lists(cells) == [[], [], [], [13, 18], [], []]  # frame -0: (60 - 3 - x) / 3
    iv = index.to_interval(3, 13, 14)
    assert (iv.start, iv.end, iv.strand) == (18, 21, "-")


def test_example_23_plus_circular_wrap_two_cells() -> None:
    genome, gm = circular_gene("+", 12)
    assert genome["c"][6:9] == "TGA" and genome["c"][27:30] == "TAA"
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 6, ("G3",), "TGA", "recoded"),
        site("c", "+", 27, ("G3",), "TAA", "terminal"),
    )
    assert n_not == 0 and skips(gm, genome) == ZERO
    index = build_sixframe(genome, dict(gm.replicons))
    assert index.frames[0].ext_nt == 426
    rec = _rt().frame_cells(index, sites[:1])
    assert lists(rec) == [[2, 102], [], [], [], [], []]
    for p in (2, 102):  # the round trip of both copies, start wrapped
        iv = index.to_interval(0, p, p + 1)
        assert (iv.start % 300, iv.end - iv.start, iv.strand) == (6, 3, "+")
        assert index.frames[0].aa[p] == STOP


def test_example_24_minus_circular_wrap_two_cells() -> None:
    genome, gm = circular_gene("-", 2)
    assert rc(genome["c"][21:24]) == "TGA" and rc(genome["c"][270:273]) == "TAA"
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "-", 21, ("G4",), "TGA", "recoded"),
        site("c", "-", 270, ("G4",), "TAA", "terminal"),
    )
    assert n_not == 0 and skips(gm, genome) == ZERO
    index = build_sixframe(genome, dict(gm.replicons))
    rec = _rt().frame_cells(index, sites[:1])
    assert lists(rec) == [[], [], [], [34, 134], [], []]
    for p in (34, 134):
        iv = index.to_interval(3, p, p + 1)
        assert (iv.start % 300, iv.end - iv.start, iv.strand) == (21, 3, "-")


def test_every_recoded_cell_round_trips_and_mirrors() -> None:
    """Rule 2's round trip and rule 3's mirror for the recoded cells of examples 21-24."""
    from gippa.mapping.decoys import build_decoys
    from gippa.rng import for_step

    rt = _rt()
    for genome, gm in (ex21(), ex22(), circular_gene("+", 12), circular_gene("-", 2)):
        index = build_sixframe(genome, dict(gm.replicons))
        sites, _ = stops(gm, genome, True)
        cells = rt.frame_cells(index, sites)
        n = gm.replicons["c"].length
        keys = rt.site_keys(sites)
        for fid, c in enumerate(cells):
            for p in c:
                iv = index.to_interval(fid, int(p), int(p) + 1)
                assert ("c", iv.strand, iv.start % n) in keys
        for kind in ("reverse", "shuffle"):
            frames = build_decoys(index, kind, 50, for_step("sixframe_index", 1))
            dc = rt.mirror_cells(kind, frames, cells)
            assert [len(x) for x in dc] == [len(x) for x in cells]


# ---------------------------------------------------------------- C2: examples 25-26


def ex25() -> tuple[dict[str, str], GeneModel]:
    """Example 21 with ``TAA`` under G1's ``U``, plus ``H1`` (``[33, 60)``) with ``O`` over TGA."""
    genome, gm = ex21(codon4="TAA")
    seq = list(genome["c"])
    put(seq, 33, "ATGGCTTGAGCTGCTGCTGCTGCTTAA")
    reps = dict(gm.replicons)
    genes = [*gm.genes, g("H1", "c", 33, 60, "+", "MAOAAAAA")]
    return {"c": "".join(seq)}, model(reps, genes)


def test_example_25_codon_mismatch_is_skipped_and_counted() -> None:
    genome, gm = ex25()
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 30, ("G1",), "TAA", "terminal"),
        site("c", "+", 57, ("H1",), "TAA", "terminal"),
    )
    assert n_not == 0
    assert skips(gm, genome) == {"length_mismatch": 0, "codon_mismatch": 2}
    assert stops(gm, genome, False) == (sites, 0)


def test_example_26_length_mismatch_keeps_the_terminal_site() -> None:
    genome, gm = ex21(aa="MAAAUAAAA*")
    sites, n_not = stops(gm, genome, True)
    assert sites == (site("c", "+", 30, ("G1",), "TAA", "terminal"),)
    assert n_not == 0
    assert skips(gm, genome) == {"length_mismatch": 1, "codon_mismatch": 0}


def test_span_not_a_multiple_of_three_is_a_length_mismatch() -> None:
    """Rule 21 (b): span 31 with ``len(aa) == 31 // 3 - 1`` still fails (#346 review S3)."""
    genome, gm = ex21()
    gm = model(dict(gm.replicons), [g("G1", "c", 2, 33, "+", "MAAAUAAAA")])
    sites, _ = stops(gm, genome, True)
    assert [s.kind for s in sites] == ["terminal"]
    assert skips(gm, genome) == {"length_mismatch": 1, "codon_mismatch": 0}


def test_empty_protein_is_not_scanned_and_not_counted() -> None:
    """Rule 21 (b), lead decision (#346 review S1): an empty ``aa`` keeps the terminal site and
    adds to neither skip count."""
    genome, gm = ex21(aa="")
    sites, n_not = stops(gm, genome, True)
    assert sites == (site("c", "+", 30, ("G1",), "TAA", "terminal"),)
    assert n_not == 0
    assert skips(gm, genome) == ZERO


def test_length_mismatch_counted_once_per_feature() -> None:
    """Rule 21 (b): a too-short protein with two ``U`` letters counts once."""
    genome, gm = ex21(aa="MAUAU")
    sites, _ = stops(gm, genome, True)
    assert [s.kind for s in sites] == ["terminal"]
    assert skips(gm, genome) == {"length_mismatch": 1, "codon_mismatch": 0}


def test_codon_mismatch_counted_once_per_residue() -> None:
    """Rule 21 (e): two ``U`` over sense codons in one protein count 2."""
    genome, gm = ex21(aa="MUAAAAAUA")  # residues 1 and 7 sit over GCT
    sites, _ = stops(gm, genome, True)
    assert [s.kind for s in sites] == ["terminal"]
    assert skips(gm, genome) == {"length_mismatch": 0, "codon_mismatch": 2}


# ---------------------------------------------------------------- C3: example 27, ties


def ex27() -> tuple[dict[str, str], GeneModel]:
    genome, gm = ex21()
    genes = [*gm.genes, g("G5", "c", 3, 18, "+", "MAAA")]  # its terminal codon is G1's U codon
    return genome, model(dict(gm.replicons), genes)


def test_example_27_nested_merge_is_one_recoded_site_with_one_cell() -> None:
    genome, gm = ex27()
    sites, n_not = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 15, ("G1", "G5"), "TGA", "recoded"),
        site("c", "+", 30, ("G1",), "TAA", "terminal"),
    )
    assert n_not == 0
    assert skips(gm, genome) == ZERO
    index = build_sixframe(genome, dict(gm.replicons))
    cells = _rt().frame_cells(index, sites)
    assert lists(cells) == [[5, 10], [], [], [], [], []]  # p = 5 once, not twice
    assert stops(gm, genome, False) == (
        (
            site("c", "+", 15, ("G5",), "TGA", "terminal"),
            site("c", "+", 30, ("G1",), "TAA", "terminal"),
        ),
        0,
    )


def test_nested_terminal_scanned_after_the_recoding_gene_stays_recoded() -> None:
    """Q22, scan order reversed from example 27: ``G5`` starts after ``G1``, so its terminal-only
    write comes last and must not demote the merged site to ``terminal`` (#346 review S2)."""
    genome, gm = ex21()
    gm = model(dict(gm.replicons), [*gm.genes, g("G5", "c", 6, 18, "+", "AAA")])
    sites, _ = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 15, ("G1", "G5"), "TGA", "recoded"),
        site("c", "+", 30, ("G1",), "TAA", "terminal"),
    )
    assert skips(gm, genome) == ZERO


def test_two_features_recoding_one_codon_merge_their_tags() -> None:
    """Edge case "Ties": alternative starts recoding the same codon give one site."""
    genome, gm = ex21()
    genes = [*gm.genes, g("G1b", "c", 9, 33, "+", "MAUAAAA")]  # residue 2 = [15, 18)
    gm = model(dict(gm.replicons), genes)
    sites, _ = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 15, ("G1", "G1b"), "TGA", "recoded"),
        site("c", "+", 30, ("G1", "G1b"), "TAA", "terminal"),
    )
    assert skips(gm, genome) == ZERO


def _run_six(tmp: Path, genome: dict[str, str], gm: GeneModel) -> Path:
    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    paths = write_annotation(gm, tmp / "ann")
    ins = {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": paths["replicons"],
        "genes": paths["genes"],
        "peptides": peptides_table(tmp / "peptides.parquet"),
    }
    out = tmp / "sixframe.npz"
    spec = reg.get("sixframe_index")
    reg.run_step(
        "sixframe_index",
        resolve_config(overrides=["seed=7"]),
        {k: ins[k] for k in spec.inputs},
        {"sixframe": out},
    )
    return out


def _event(caplog: pytest.LogCaptureFixture) -> dict:  # type: ignore[type-arg]
    ev = [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    assert len(ev) == 1
    return ev[0].data  # type: ignore[attr-defined,no-any-return]


def test_example_27_event_counts_a_merged_site_once_as_recoded(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from gippa.mapping.sixframe import load_stop_cells

    caplog.set_level(logging.INFO)
    genome, gm = ex27()
    out = _run_six(tmp_path, genome, gm)
    assert _event(caplog) == {
        "terminal": 1,
        "recoded": 1,
        "not_a_stop": 0,
        "cells_target": 2,
        "cells_decoy": {"reverse": 2, "shuffle": 2},
        "recoded_skipped": ZERO,
    }
    cells = load_stop_cells(out)
    assert [(s.codon_start, s.kind, s.locus_tags) for s in cells.sites] == [
        (15, "recoded", ("G1", "G5")),
        (30, "terminal", ("G1",)),
    ]


def test_event_reports_the_skip_counts(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.INFO)
    genome, gm = ex25()
    _run_six(tmp_path, genome, gm)
    data = _event(caplog)
    assert data["recoded"] == 0 and data["terminal"] == 2
    assert data["recoded_skipped"] == {"length_mismatch": 0, "codon_mismatch": 2}


def test_event_on_example_21_counts_the_recoded_cell(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.INFO)
    genome, gm = ex21()
    _run_six(tmp_path, genome, gm)
    data = _event(caplog)
    assert (data["terminal"], data["recoded"], data["cells_target"]) == (1, 1, 2)
    assert data["cells_decoy"] == {"reverse": 2, "shuffle": 2}


# ---------------------------------------------------------------- C4: example 28, origin, sort


def ex28() -> tuple[dict[str, str], GeneModel]:
    return layout(
        [
            ("G6", "+", "cds", ["ATG", "GCT", "GCT", "TAG", "GCT", "TAA"], "MAAOA"),
            ("G7", "+", "pseudogene", ["ATG", "GCT", "TGA", "GCT", "TAA"], "MAUA"),
            ("G8", "+", "cds", ["TGA", "GCT", "GCT", "TAA"], "UAA"),
            ("G9", "+", "cds", ["ATG", "GCT", "GCT", "TGA", "GCT", "TAA"], "MAX*A"),
        ]
    )


def test_example_28_pyl_and_ignored_letters() -> None:
    genome, gm = ex28()
    by = {x.locus_tag: x.interval for x in gm.genes}
    sites, n_not = stops(gm, genome, True)
    g6, g8, g9 = by["G6"].start, by["G8"].start, by["G9"].start
    assert sites == (
        site("c", "+", g6 + 9, ("G6",), "TAG", "recoded"),
        site("c", "+", g6 + 15, ("G6",), "TAA", "terminal"),
        site("c", "+", g8 + 9, ("G8",), "TAA", "terminal"),
        site("c", "+", g9 + 15, ("G9",), "TAA", "terminal"),
    )
    assert n_not == 0
    assert skips(gm, genome) == ZERO  # G7, G8 and G9 are not counted


def test_o_over_tag_on_minus_reads_cta_on_plus() -> None:
    genome, gm = layout([("P1", "-", "sorf", ["ATG", "GCT", "TAG", "GCT", "TGA"], "MAOA")])
    iv = gm.genes[0].interval
    cs = iv.end - 9  # residue 2: [e - 9, e - 6)
    assert genome["c"][cs : cs + 3] == "CTA"
    sites, _ = stops(gm, genome, True)
    assert sites == (
        site("c", "-", iv.start, ("P1",), "TGA", "terminal"),
        site("c", "-", cs, ("P1",), "TAG", "recoded"),
    )


def test_lower_case_letters_are_not_recoded() -> None:
    genome, gm = ex21(aa="MAAAuAAAA")
    sites, _ = stops(gm, genome, True)
    assert [s.kind for s in sites] == ["terminal"]
    assert skips(gm, genome) == ZERO


def test_not_a_stop_feature_is_not_scanned() -> None:
    """Q24: a feature counted under ``not_a_stop`` gets no recoded site and no skip count."""
    seq = filler(60)
    put(seq, 3, "ATGGCTGCTGCTTGAGCTGCTGCTGCTGCT")
    reps = {"c": Replicon("c", 60, "linear")}
    gm = model(reps, [g("Q1", "c", 3, 33, "+", "MAAAUAAAA", partial=True)])
    genome = {"c": "".join(seq)}
    assert stops(gm, genome, True) == ((), 1)
    assert skips(gm, genome) == ZERO


def test_origin_crossing_recoded_codon() -> None:
    """Edge case: a recoded codon ``[n - 1, n + 2)``: ``codon_start = n - 1``, read wrapped, one
    cell (frame ``+2``, p = 99; 599 > ``ext_nt``)."""
    genome, gm = circular_gene("+", 9, start=272)
    assert genome["c"][299] + genome["c"][:2] == "TGA"
    sites, _ = stops(gm, genome, True)
    assert sites == (
        site("c", "+", 29, ("G3",), "TAA", "terminal"),
        site("c", "+", 299, ("G3",), "TGA", "recoded"),
    )
    index = build_sixframe(genome, dict(gm.replicons))
    cells = _rt().frame_cells(index, sites[1:])
    assert lists(cells) == [[], [], [99], [], [], []]
    iv = index.to_interval(2, 99, 100)
    assert (iv.start % 300, iv.end - iv.start, iv.strand) == (299, 3, "+")


def test_several_u_one_site_each_and_the_sort_order() -> None:
    """Rule 21 (g): one site per ``U``; sites of every kind sorted together by (contig order,
    strand, codon_start); contig order is the gene model's, not alphabetical."""
    gz, gmz = layout(
        [("Z1", "-", "cds", ["ATG", "TGA", "GCT", "TGA", "GCT", "TAA"], "MUAUA")], contig="z"
    )
    ga, gma = layout(
        [("A1", "+", "cds", ["ATG", "GCT", "TGA", "GCT", "TGA", "TAG"], "MAUAU")], contig="a"
    )
    reps = {"z": gmz.replicons["z"], "a": gma.replicons["a"]}
    gm = model(reps, [*gmz.genes, *gma.genes])
    genome = {"z": gz["z"], "a": ga["a"]}
    z, a = gmz.genes[0].interval, gma.genes[0].interval
    sites, _ = stops(gm, genome, True)
    assert sites == (
        site("z", "-", z.start, ("Z1",), "TAA", "terminal"),
        site("z", "-", z.end - 12, ("Z1",), "TGA", "recoded"),  # residue 3
        site("z", "-", z.end - 6, ("Z1",), "TGA", "recoded"),  # residue 1
        site("a", "+", a.start + 6, ("A1",), "TGA", "recoded"),
        site("a", "+", a.start + 12, ("A1",), "TGA", "recoded"),
        site("a", "+", a.start + 15, ("A1",), "TAG", "terminal"),
    )
    index = build_sixframe(genome, reps)
    assert sum(len(c) for c in _rt().frame_cells(index, sites)) == 6


# ---------------------------------------------------------------- C5: recoded=False, empty input


ALL_INPUTS = {
    "ex21": ex21,
    "ex22": ex22,
    "ex23": lambda: circular_gene("+", 12),
    "ex24": lambda: circular_gene("-", 2),
    "ex25": ex25,
    "ex26": lambda: ex21(aa="MAAAUAAAA*"),
    "ex27": ex27,
    "ex28": ex28,
    "origin": lambda: circular_gene("+", 9, start=272),
}


@pytest.mark.parametrize("name", sorted(ALL_INPUTS))
def test_recoded_false_is_the_16a_output(name: str) -> None:
    genome, gm = ALL_INPUTS[name]()
    default = stops(gm, genome)
    assert stops(gm, genome, False) == default
    assert all(s.kind == "terminal" for s in default[0])
    rec, n_not = stops(gm, genome, True)
    assert n_not == default[1]
    # the terminal part of the recoded output is the default output minus merged keys
    merged = {(s.contig, s.strand, s.codon_start) for s in rec if s.kind == "recoded"}
    want = [s for s in default[0] if (s.contig, s.strand, s.codon_start) not in merged]
    assert [s for s in rec if s.kind == "terminal"] == want


@pytest.mark.parametrize("toy", ["c1", "c2", "world"])
def test_16a_toys_without_u_or_o_give_the_same_sites(toy: str) -> None:
    if toy == "world":
        genome, _reps, gm = rt_world()
    else:
        genome, _reps, gm = toy_c1() if toy == "c1" else toy_c2()
    assert stops(gm, genome, True) == stops(gm, genome)


def test_proteins_without_u_or_o_give_the_same_sites_and_no_skips() -> None:
    genome, gm = ex21(aa="MAAAAAAAA")
    assert stops(gm, genome, True) == stops(gm, genome)
    assert skips(gm, genome) == ZERO


@pytest.mark.parametrize("recoded", [None, False, True])
def test_empty_gene_model(recoded: bool | None) -> None:
    genome, gm = ex21()
    empty = model(dict(gm.replicons), [])
    assert stops(empty, genome, recoded) == ((), 0)
    assert skips(empty, genome) == ZERO


def test_recoded_skips_keys_and_purity() -> None:
    genome, gm = ex25()
    before = dict(genome)
    a = skips(gm, genome)
    b = skips(gm, genome)
    assert a == b and set(a) == {"length_mismatch", "codon_mismatch"}
    assert all(type(v) is int for v in a.values())
    assert genome == before
