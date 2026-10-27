"""M2-19 regions: examples 1-10, rules and edge cases (docs/specs/M2-19-regions.md; GIPPA_SPEC.md
§5.2, §6.5, §7.1.7, §9.3).

The spec names this file ``tests/unit/test_regions.py``; M4-17 claims
``tests/unit/context/test_regions.py`` and test basenames must be unique (scripts/ci/
check_test_names.py), so it carries the ``regions_`` prefix instead.

Resolutions pinned here where the spec is silent (lead decision 6 of M2-index overrides Q2):

* Plasmid: every contig listed in geNOMAD's (already default-filtered) ``plasmid_summary`` is a
  plasmid call, whatever its ``plasmid_score`` — there is no 0.7 cut-off key.
* ``build_regions(..., cfg)`` takes the resolved ``GippaConfig``
  (``annotation.region_min_overlap_frac``, ``context.region_types``).
* CRISPR: a gene is a flanking ``cas`` gene when its gap to the array (``array.start - gene.end``
  upstream, ``gene.start - array.end`` downstream, 0 when overlapping) is ≤ 1,000 bp; the gap is
  measured to the array only (example 7: cas2 is 600 bp from cas1 but 2,200 bp from the array and
  is excluded). The region interval is the union; its genes are its ≥ 50 % members (Q3).
* AMR: ``expert`` containing ``amrfinder`` (Bakta's expert type, as parsed by M2-03).
* Prophage subsets: a Bakta gene takes the geNOMAD ``annotation_description`` of the virus gene with
  the identical interval; only that unambiguous case and the no-match case are pinned.
* ``Region.genes`` lists exactly the region's ``RegionGene`` locus tags; for a region that does not
  cross the origin they are in genome order.
* ``regions_bed`` returns one ``BedRecord`` per non-AMR region named by its ``region_id`` with the
  region's (unwrapped) interval; ``write_bed`` splits origin-crossing records (M1-04).
"""

from __future__ import annotations

import dataclasses
import importlib
import itertools
from pathlib import Path

import pytest
from helpers.regions_oracle import TSV_HEADER

from gippa.annotate.antismash import BgcRegion
from gippa.annotate.genomad import GenomadResult, PlasmidCall, Provirus, VirusGene
from gippa.config import GippaConfig, resolve_config
from gippa.config.model import AnnotationConfig
from gippa.io.bed import write_bed
from gippa.io.tsv import read_tsv
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, NcFeature

regions = importlib.import_module("gippa.annotate.regions")  # M2-19

TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"
CHR1 = Replicon("chr1", 100_000, "circular", "chromosome")
P1 = Replicon("p1", 10_000, "circular", "plasmid")
LIN = Replicon("lin", 10_000, "linear", "chromosome")
REPS = {"chr1": CHR1, "p1": P1}


def iv(start: int, end: int, contig: str = "chr1", strand: str = "+") -> Interval:
    return Interval(contig, start, end, strand)  # type: ignore[arg-type]


def gene(
    tag: str,
    start: int,
    end: int,
    contig: str = "chr1",
    *,
    product: str = "hypothetical protein",
    name: str | None = None,
    expert: tuple[str, ...] = (),
    strand: str = "+",
) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=iv(start, end, contig, strand),
        aa="M",
        nt=None,
        start_codon="ATG",
        gene=name,
        product=product,
        expert=expert,
    )


def crispr(fid: str, start: int, end: int, contig: str = "chr1") -> NcFeature:
    return NcFeature(id=fid, type="crispr", interval=iv(start, end, contig, "."))


def provirus(start: int | None, end: int | None, contig: str = "chr1", score: float = 0.97):
    span = None if start is None or end is None else iv(start, end, contig, ".")
    return Provirus(contig, span, score, None, 5, TAX)


def genomad(proviruses=(), plasmids=(), vgenes=()) -> GenomadResult:
    return GenomadResult(tuple(proviruses), tuple(plasmids), tuple(vgenes))


def plasmid_call(contig: str, score: float) -> PlasmidCall:
    return PlasmidCall(contig, score, None, 2, (), ())


def bgc(start: int, end: int, contig: str = "chr1", products=("NRPS",), core=()) -> BgcRegion:
    return BgcRegion(contig, iv(start, end, contig, "."), tuple(products), tuple(core))


def cfg(*overrides: str) -> GippaConfig:
    return resolve_config(overrides=list(overrides))


def build(genes=(), features=(), reps=REPS, gm=None, bgcs=None, config=None):
    return regions.build_regions(
        list(genes), list(features), reps, gm, bgcs, config if config is not None else cfg()
    )


def by_id(regs) -> dict:
    return {r.region_id: r for r in regs}


def rg_of(rgs, region_id: str) -> dict:
    return {x.locus_tag: x for x in rgs if x.region_id == region_id}


# ------------------------------------------------------------------ dataclasses and interface


def test_region_and_region_gene_are_frozen_slotted_dataclasses() -> None:
    r = regions.Region("amr_1", "amr", "chr1", None, "AMRFinderPlus", None, ("g1",), "")
    x = regions.RegionGene("amr_1", "g1", 1.0, None, False)
    for obj in (r, x):
        assert dataclasses.is_dataclass(obj)
        assert hasattr(type(obj), "__slots__")
        with pytest.raises(dataclasses.FrozenInstanceError):
            obj.region_id = "x"  # type: ignore[misc]
    assert [f.name for f in dataclasses.fields(regions.Region)] == [
        "region_id", "region_type", "contig", "interval", "source_tool", "score", "genes", "notes",
    ]  # fmt: skip
    assert [f.name for f in dataclasses.fields(regions.RegionGene)] == [
        "region_id", "locus_tag", "overlap_frac", "subset", "core_biosynthetic",
    ]  # fmt: skip


# ------------------------------------------------------------------ examples 1-3: overlap_frac


def test_example_1_exactly_half_is_a_member() -> None:
    g, r = iv(1000, 1300), iv(1150, 5000, strand=".")
    assert regions.overlap_frac(g, r, CHR1) == pytest.approx(0.5)
    assert regions.members([gene("g", 1000, 1300)], r, CHR1) == [("g", pytest.approx(0.5))]


def test_example_2_just_below_half_is_not_a_member() -> None:
    g, r = iv(1000, 1300), iv(1151, 5000, strand=".")
    assert regions.overlap_frac(g, r, CHR1) == pytest.approx(149 / 300)
    assert regions.members([gene("g", 1000, 1300)], r, CHR1) == []


def test_example_3_origin_crossing_gene() -> None:
    g, r = iv(99_900, 100_200), iv(0, 5000, strand=".")
    assert regions.overlap_frac(g, r, CHR1) == pytest.approx(200 / 300)
    assert regions.members([gene("g", 99_900, 100_200)], r, CHR1) == [
        ("g", pytest.approx(200 / 300))
    ]


@pytest.mark.parametrize(
    ("g", "r", "want"),
    [
        (iv(0, 300), iv(99_900, 100_200), 200 / 300),  # region crosses the origin
        (iv(99_950, 100_050), iv(99_900, 100_100), 1.0),  # both cross
        (iv(99_950, 100_050), iv(0, 100_000), 1.0),  # whole replicon
        (iv(10, 20), iv(30, 40), 0.0),  # disjoint
        (iv(10, 20), iv(20, 40), 0.0),  # half-open: touching is no overlap
        (iv(1000, 1301), iv(1150, 5000), 151 / 301),  # 301 bp gene, 151 bp inside
    ],
)
def test_overlap_frac_origin_aware_cases(g: Interval, r: Interval, want: float) -> None:
    assert regions.overlap_frac(g, r, CHR1) == pytest.approx(want)


def test_overlap_frac_on_a_linear_contig() -> None:
    g, r = iv(4900, 5100, "lin"), iv(0, 5000, "lin", ".")
    assert regions.overlap_frac(g, r, LIN) == pytest.approx(0.5)


def test_members_boundary_with_odd_length_gene_and_min_frac() -> None:
    genes = [gene("in", 1000, 1301), gene("out", 999, 1300)]  # 151/301 vs 150/301 inside
    r = iv(1150, 5000, strand=".")
    assert [t for t, _ in regions.members(genes, r, CHR1)] == ["in"]
    assert regions.members(genes, r, CHR1, min_frac=0.6) == []
    assert {t for t, _ in regions.members(genes, r, CHR1, min_frac=0.0)} == {"in", "out"}


def test_members_ignores_genes_on_other_contigs() -> None:
    genes = [gene("a", 1000, 1300, "chr1"), gene("b", 1000, 1300, "p1")]
    assert regions.members(genes, iv(0, 5000, strand="."), CHR1) == [("a", pytest.approx(1.0))]


# ------------------------------------------------------------------ examples 4-5: subsets


@pytest.mark.parametrize(
    ("product", "want"),
    [
        ("phage integrase", "integrase"),
        ("CI repressor", "repressor"),
        ("cI", "repressor"),
        ("tail fiber protein", "structural"),
        ("holin", "lysis"),
        ("DNA methylase", "cargo"),
        ("microcin immunity protein", "cargo"),
    ],
)
def test_example_4_subsets(product: str, want: str) -> None:
    assert regions.prophage_subset(product, None) == want


def test_example_5_rule_order() -> None:
    assert (
        regions.prophage_subset("tail fiber assembly integrase-like protein", None) == "integrase"
    )


@pytest.mark.parametrize(
    ("product", "marker", "want"),
    [
        ("site-specific recombinase", None, "integrase"),
        ("excisionase", None, "integrase"),
        ("INTEGRASE", None, "integrase"),  # case-insensitive
        ("antirepressor", None, "repressor"),
        ("Cro/CI family transcriptional regulator", None, "repressor"),  # Cro, case-sensitive
        ("cro", None, "cargo"),  # Cro is case-sensitive
        ("crossover junction endodeoxyribonuclease", None, "cargo"),  # no 'cro' substring match
        ("circular permuted protein", None, "cargo"),
        ("Phage portal protein", None, "structural"),
        ("phage tape measure protein", None, "structural"),
        ("baseplate assembly protein J", None, "structural"),
        ("terminase large subunit", None, "structural"),
        ("head-tail connector protein", None, "structural"),
        ("phage lysozyme", None, "lysis"),
        ("endolysin", None, "lysis"),
        ("spanin", None, "lysis"),
        ("hypothetical protein", "Phage major capsid protein E", "structural"),  # geNOMAD marker
        ("hypothetical protein", "Integrase", "integrase"),
        ("hypothetical protein", None, "cargo"),
        ("hypothetical protein", "", "cargo"),
        ("tailspike-like protein", None, "cargo"),  # whole words only
    ],
)
def test_subset_keyword_rules(product: str, marker: str | None, want: str) -> None:
    assert regions.prophage_subset(product, marker) == want


# ------------------------------------------------------------------ example 6: ids and numbering


def test_example_6_prophage_and_bgc_ids_by_contig_then_start() -> None:
    genes = [gene("c1", 50_100, 50_400), gene("p1g", 1100, 1400, "p1"), gene("b1", 20_100, 20_400)]
    gm = genomad(proviruses=[provirus(1000, 3000, "p1", 0.9), provirus(50_000, 62_000)])
    regs, _ = build(genes, gm=gm, bgcs=[bgc(20_000, 30_000)])
    ids = by_id(regs)
    assert ids["prophage_1"].contig == "chr1"
    assert ids["prophage_1"].interval == iv(50_000, 62_000, strand=".")
    assert ids["prophage_2"].contig == "p1"
    assert ids["bgc_1"].contig == "chr1"
    assert ids["bgc_1"].interval == iv(20_000, 30_000, strand=".")


def test_numbering_by_start_within_a_contig_regardless_of_input_order() -> None:
    # geNOMAD lists proviruses by score, not position (real K-12 / Nissle summaries)
    gm = genomad(proviruses=[provirus(70_000, 80_000, score=0.99), provirus(10_000, 20_000)])
    regs, _ = build(gm=gm, bgcs=[bgc(60_000, 61_000), bgc(5000, 6000), bgc(100, 200, "p1")])
    ids = by_id(regs)
    assert ids["prophage_1"].interval.start == 10_000
    assert ids["prophage_2"].interval.start == 70_000
    assert [ids[f"bgc_{k}"].interval.start for k in (1, 2, 3)] == [5000, 60_000, 100]
    assert ids["bgc_3"].contig == "p1"


# ------------------------------------------------------------------ prophage rules


def test_prophage_fields_members_and_subsets() -> None:
    genes = [
        gene("int", 50_000, 51_000, product="phage integrase"),
        gene("hyp", 51_000, 52_000),  # subset from the geNOMAD marker of the identical gene
        gene("half", 61_700, 62_300, product="holin"),  # exactly 300/600 inside
        gene("out", 61_701, 62_301, product="holin"),  # 299/600
        gene("far", 70_000, 71_000),
    ]
    vg = VirusGene("chr1", iv(51_000, 52_000), "GENOMAD.1.VV", "Phage major capsid protein", True)
    gm = genomad(proviruses=[provirus(50_000, 62_000, score=0.97)], vgenes=[vg])
    regs, rgs = build(genes, gm=gm)
    (r,) = [x for x in regs if x.region_type == "prophage"]
    assert (r.region_id, r.source_tool, r.score, r.notes) == ("prophage_1", "geNOMAD", 0.97, TAX)
    assert r.genes == ("int", "hyp", "half")
    members = rg_of(rgs, "prophage_1")
    assert {t: m.subset for t, m in members.items()} == {
        "int": "integrase",
        "hyp": "structural",
        "half": "lysis",
    }
    assert members["half"].overlap_frac == pytest.approx(0.5)
    assert members["int"].overlap_frac == pytest.approx(1.0)
    assert not any(m.core_biosynthetic for m in members.values())


def test_whole_contig_virus_becomes_a_whole_contig_prophage() -> None:
    gm = genomad(proviruses=[provirus(None, None, "p1", 0.88)])
    regs, _rgs = build([gene("pg", 100, 400, "p1")], gm=gm)
    (r,) = [x for x in regs if x.region_type == "prophage"]
    assert r.interval == Interval("p1", 0, 10_000, r.interval.strand)
    assert r.score == 0.88 and r.source_tool == "geNOMAD"
    assert "whole_contig" in r.notes
    assert r.genes == ("pg",)


def test_min_overlap_frac_comes_from_the_config() -> None:
    genes = [gene("g", 1000, 1300)]
    gm = genomad(proviruses=[provirus(1150, 5000)])
    assert by_id(build(genes, gm=gm)[0])["prophage_1"].genes == ("g",)
    strict = cfg("annotation.region_min_overlap_frac=0.6")
    assert by_id(build(genes, gm=gm, config=strict)[0])["prophage_1"].genes == ()


# ------------------------------------------------------------------ plasmid rules


def test_plasmid_from_bakta_only_when_genomad_is_absent() -> None:
    genes = [gene("c", 100, 400), gene("p_a", 100, 400, "p1"), gene("p_b", 9_900, 10_200, "p1")]
    regs, rgs = build(genes, gm=None)
    (r,) = [x for x in regs if x.region_type == "plasmid"]
    assert (r.region_id, r.contig, r.source_tool, r.score) == ("plasmid_1", "p1", "Bakta", None)
    assert r.interval == Interval("p1", 0, 10_000, r.interval.strand)
    assert set(r.genes) == {"p_a", "p_b"}  # all genes on it, origin-crossing one included
    assert all(m.overlap_frac == pytest.approx(1.0) for m in rg_of(rgs, "plasmid_1").values())


def test_plasmid_source_bakta_plus_genomad_and_genomad_only() -> None:
    reps = {"chr1": CHR1, "p1": P1, "c2": Replicon("c2", 5000, "linear", "unknown")}
    gm = genomad(plasmids=[plasmid_call("c2", 0.9932), plasmid_call("p1", 0.9951)])
    regs, _ = build([gene("x", 10, 400, "c2")], reps=reps, gm=gm)
    ids = by_id(regs)
    assert (ids["plasmid_1"].contig, ids["plasmid_1"].source_tool) == ("p1", "Bakta+geNOMAD")
    assert ids["plasmid_1"].score == pytest.approx(0.9951)
    assert (ids["plasmid_2"].contig, ids["plasmid_2"].source_tool) == ("c2", "geNOMAD")
    assert ids["plasmid_2"].score == pytest.approx(0.9932)
    assert ids["plasmid_2"].interval == Interval("c2", 0, 5000, ids["plasmid_2"].interval.strand)
    assert ids["plasmid_2"].genes == ("x",)
    assert all(";" not in r.source_tool for r in regs)


def test_genomad_plasmid_summary_rows_need_no_score_cutoff() -> None:
    """M2-index lead decision 6: no cut-off key; the summary is geNOMAD's filtered call set."""
    reps = {"chr1": CHR1, "c2": Replicon("c2", 5000, "linear", "unknown")}
    regs, _ = build(reps=reps, gm=genomad(plasmids=[plasmid_call("c2", 0.62)]))
    assert [r.contig for r in regs if r.region_type == "plasmid"] == ["c2"]
    assert "regions" not in AnnotationConfig.model_fields


# ------------------------------------------------------------------ BGC rules


def test_bgc_fields_and_core_flags() -> None:
    genes = [gene("k1", 20_000, 21_000), gene("k2", 21_000, 22_000), gene("edge", 29_800, 30_400)]
    b = bgc(20_000, 30_000, products=("NRP-metallophore", "NRPS"), core=("k2",))
    regs, rgs = build(genes, bgcs=[b])
    (r,) = [x for x in regs if x.region_type == "bgc"]
    assert (r.region_id, r.source_tool, r.score) == ("bgc_1", "antiSMASH", None)
    assert r.notes == "products=NRP-metallophore,NRPS"
    assert r.genes == ("k1", "k2")  # 'edge' has 200/600 inside
    members = rg_of(rgs, "bgc_1")
    assert {t: m.core_biosynthetic for t, m in members.items()} == {"k1": False, "k2": True}
    assert all(m.subset is None for m in members.values())


# ------------------------------------------------------------------ example 7: CRISPR


def test_example_7_crispr_with_flanking_cas_genes() -> None:
    genes = [
        gene("cas1", 61_500, 62_400, name="cas1"),
        gene("cas2", 63_000, 63_300, name="cas2"),
    ]
    regs, _rgs = build(genes, [crispr("arr", 60_000, 60_800)])
    (r,) = [x for x in regs if x.region_type == "crispr"]
    assert r.region_id == "crispr_1" and r.source_tool == "Bakta"
    assert r.genes == ("cas1",)
    assert (r.interval.start, r.interval.end) == (60_000, 62_400)
    assert r.score is None


def test_crispr_cas_gap_boundary_both_sides_and_product_rule() -> None:
    genes = [
        gene("up", 58_000, 59_000, product="CRISPR-associated endoribonuclease Cas2"),  # gap 1000
        gene("up_far", 56_000, 58_999, name="Cas9"),  # gap 1001
        gene("down", 61_800, 62_000, name="casA"),  # gap 1000
        gene("noncas", 60_850, 61_000, name="iap"),  # close but not cas
    ]
    regs, _ = build(genes, [crispr("arr", 60_000, 60_800)])
    (r,) = [x for x in regs if x.region_type == "crispr"]
    assert (r.interval.start, r.interval.end) == (58_000, 62_000)
    assert "up" in r.genes and "down" in r.genes and "up_far" not in r.genes


def test_crispr_without_cas_genes_is_the_array() -> None:
    regs, _ = build(
        [gene("x", 70_000, 71_000)], [crispr("a2", 62_387, 62_799), crispr("a1", 100, 500)]
    )
    ids = by_id(regs)
    assert (ids["crispr_1"].interval.start, ids["crispr_1"].interval.end) == (100, 500)
    assert (ids["crispr_2"].interval.start, ids["crispr_2"].interval.end) == (62_387, 62_799)
    assert ids["crispr_2"].genes == ()


# ------------------------------------------------------------------ example 8: AMR


def test_example_8_amr_one_gene_set_per_contig() -> None:
    genes = [
        gene("a1", 100, 400, expert=("amrfinder",)),
        gene("n1", 500, 800, expert=("expert_proteins",)),
        gene("a2", 5000, 5300, expert=("expert_proteins", "amrfinder")),
        gene("a3", 100, 400, "p1", expert=("amrfinder",)),
    ]
    regs, rgs = build(genes)
    ids = by_id(regs)
    assert (ids["amr_1"].contig, ids["amr_1"].interval, ids["amr_1"].genes) == (
        "chr1",
        None,
        ("a1", "a2"),
    )
    assert (ids["amr_2"].contig, ids["amr_2"].genes) == ("p1", ("a3",))
    assert ids["amr_1"].source_tool == "AMRFinderPlus" and ids["amr_1"].score is None
    assert set(rg_of(rgs, "amr_1")) == {"a1", "a2"}
    assert len([r for r in regs if r.region_type == "amr"]) == 2


# ------------------------------------------------------------------ examples 9-10: writers


def test_example_9_tsv_row(tmp_path: Path) -> None:
    genes = [gene("g1", 50_000, 51_000), gene("g2", 52_000, 53_000)]
    gm = genomad(proviruses=[provirus(1000, 3000, "p1", 0.9), provirus(50_000, 62_000, score=0.97)])
    regs, _ = build(
        genes, gm=gm, bgcs=[bgc(20_000, 30_000)], config=cfg("context.region_types=[prophage]")
    )
    out = tmp_path / "regions.tsv"
    regions.write_regions_tsv(out, regs)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == TSV_HEADER
    assert f"prophage_1\tprophage\tchr1\t50001\t62000\tgeNOMAD\t0.97\tg1;g2\t{TAX}" in lines
    assert len(lines) == 3


def test_example_10_origin_crossing_prophage_tsv_and_bed(tmp_path: Path) -> None:
    genes = [gene("ori", 99_900, 100_200), gene("g0", 100, 400)]
    gm = genomad(proviruses=[provirus(99_000, 101_000)])
    regs, _ = build(genes, gm=gm, config=cfg("context.region_types=[prophage]"))
    (r,) = regs
    assert r.interval == iv(99_000, 101_000, strand=r.interval.strand)
    assert set(r.genes) == {"ori", "g0"}
    tsv = tmp_path / "regions.tsv"
    regions.write_regions_tsv(tsv, regs)
    _, rows = read_tsv(tsv)
    assert rows[0][3:5] == ["99001", "101000"]  # end > length kept (§6.5)
    bed = tmp_path / "regions.bed"
    write_bed(bed, regions.regions_bed(regs, REPS), REPS)
    fields = [ln.split("\t") for ln in bed.read_text().splitlines()]
    assert sorted((f[0], f[1], f[2], f[3]) for f in fields) == [
        ("chr1", "0", "1000", "prophage_1_part2"),
        ("chr1", "99000", "100000", "prophage_1_part1"),
    ]


def test_tsv_amr_row_has_na_coordinates_and_bed_leaves_it_out(tmp_path: Path) -> None:
    genes = [
        gene("a1", 100, 400, expert=("amrfinder",)),
        gene("a2", 900, 1200, expert=("amrfinder",)),
    ]
    regs, _ = build(genes, [crispr("arr", 5000, 5400)], reps={"chr1": CHR1})
    tsv = tmp_path / "regions.tsv"
    regions.write_regions_tsv(tsv, regs)
    header, rows = read_tsv(tsv)
    amr = {r[0]: r for r in rows}["amr_1"]
    assert dict(zip(header, amr, strict=True)) == {
        "region_id": "amr_1", "region_type": "amr", "contig": "chr1", "start": None,
        "end": None, "source_tool": "AMRFinderPlus", "score": None, "genes": "a1;a2",
        "notes": amr[8],
    }  # fmt: skip
    beds = regions.regions_bed(regs, {"chr1": CHR1})
    assert sorted(b.name for b in beds) == ["crispr_1"]
    (b,) = beds
    assert (b.interval.contig, b.interval.start, b.interval.end) == ("chr1", 5000, 5400)


def test_tsv_bgc_and_plasmid_rows(tmp_path: Path) -> None:
    regs, _ = build(
        [gene("pg", 100, 400, "p1")],
        bgcs=[bgc(20_000, 30_000, products=("NRPS", "T1PKS"))],
        config=cfg("context.region_types=[plasmid,bgc]"),
    )
    tsv = tmp_path / "regions.tsv"
    regions.write_regions_tsv(tsv, regs)
    lines = tsv.read_text().splitlines()
    assert "bgc_1\tbgc\tchr1\t20001\t30000\tantiSMASH\tNA\tNA\tproducts=NRPS,T1PKS" in lines
    assert any(ln.startswith("plasmid_1\tplasmid\tp1\t1\t10000\tBakta\tNA\tpg\t") for ln in lines)


# ------------------------------------------------------------------ edge cases


def test_empty_everything_writes_a_header_only_tsv(tmp_path: Path) -> None:
    regs, rgs = build(reps={"chr1": CHR1})
    assert regs == [] and rgs == []
    out = tmp_path / "regions.tsv"
    regions.write_regions_tsv(out, [])
    assert out.read_text() == "\t".join(TSV_HEADER) + "\n"
    assert regions.regions_bed([], {"chr1": CHR1}) == []


def test_without_optional_layers_no_prophage_or_bgc() -> None:
    regs, _ = build([gene("g", 100, 400)], gm=None, bgcs=None)
    assert {r.region_type for r in regs} == {"plasmid"}  # p1 is a Bakta plasmid


@pytest.mark.parametrize("types", [["prophage"], ["bgc", "amr"], ["crispr"], []])
def test_region_types_restrict_the_built_types(types: list[str]) -> None:
    genes = [gene("a", 100, 400, expert=("amrfinder",)), gene("p", 100, 400, "p1")]
    gm = genomad(proviruses=[provirus(50_000, 60_000)], plasmids=[plasmid_call("p1", 0.99)])
    config = cfg(f"context.region_types=[{','.join(types)}]")
    regs, rgs = build(
        genes, [crispr("c", 70_000, 70_400)], gm=gm, bgcs=[bgc(1, 900)], config=config
    )
    assert {r.region_type for r in regs} == set(types)
    assert {x.region_id for x in rgs} <= {r.region_id for r in regs}


def test_gene_in_prophage_inside_a_plasmid_belongs_to_both() -> None:
    gm = genomad(proviruses=[provirus(1000, 3000, "p1")])
    _regs, rgs = build([gene("pg", 1500, 1800, "p1", product="phage integrase")], gm=gm)
    owners = sorted(x.region_id for x in rgs if x.locus_tag == "pg")
    assert owners == ["plasmid_1", "prophage_1"]
    assert rg_of(rgs, "prophage_1")["pg"].subset == "integrase"
    assert rg_of(rgs, "plasmid_1")["pg"].subset is None


def test_region_genes_agree_with_region_gene_lists() -> None:
    genes = [gene(f"g{k}", 1000 * k, 1000 * k + 600) for k in range(1, 60)]
    genes.append(gene("pg", 100, 400, "p1", expert=("amrfinder",)))
    gm = genomad(proviruses=[provirus(5000, 20_000)], plasmids=[plasmid_call("p1", 0.99)])
    regs, rgs = build(genes, [crispr("c", 40_100, 40_400)], gm=gm, bgcs=[bgc(30_000, 45_000)])
    for r in regs:
        assert tuple(x.locus_tag for x in rgs if x.region_id == r.region_id) == r.genes, r
    assert len({(x.region_id, x.locus_tag) for x in rgs}) == len(rgs)


def test_build_regions_is_deterministic_and_input_order_free() -> None:
    genes = [gene(f"g{k}", 1000 * k, 1000 * k + 600) for k in range(1, 60)]
    pvs = [provirus(5000, 20_000), provirus(50_000, 52_000), provirus(100, 900, "p1")]
    bgcs = [bgc(30_000, 45_000), bgc(10, 600), bgc(2000, 4000, "p1")]
    want = build(genes, gm=genomad(proviruses=pvs), bgcs=bgcs)
    assert build(genes, gm=genomad(proviruses=pvs), bgcs=bgcs) == want
    for p, b in itertools.islice(
        zip(itertools.permutations(pvs), itertools.permutations(bgcs), strict=False), 6
    ):
        got = build(list(reversed(genes)), gm=genomad(proviruses=p), bgcs=list(b))
        assert sorted(got[0], key=lambda r: r.region_id) == sorted(
            want[0], key=lambda r: r.region_id
        )


def test_ids_match_the_type_n_pattern_and_count_from_one() -> None:
    gm = genomad(proviruses=[provirus(5000, 20_000), provirus(100, 900, "p1")])
    regs, _ = build(gm=gm, bgcs=[bgc(30_000, 45_000)])
    by_type: dict[str, list[int]] = {}
    for r in regs:
        t, _, n = r.region_id.rpartition("_")
        assert t == r.region_type
        by_type.setdefault(t, []).append(int(n))
    assert all(sorted(v) == list(range(1, len(v) + 1)) for v in by_type.values())
