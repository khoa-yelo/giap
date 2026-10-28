"""Contract test: regions from the recorded K-12 MG1655 and Nissle 1917 runs (M2-19 acceptance:
"recorded Nissle Bakta + geNOMAD + antiSMASH: ≥ 1 prophage, ≥ 1 BGC, 2 plasmid regions (pMUT1,
pMUT2)").

Inputs are the merged producers' parses of real recordings: the committed Bakta 1.12.2 full-DB
JSON excerpts (``tests/data/recorded/regions/``, PROVENANCE.md), the complete committed geNOMAD
1.12.0 summaries and the trimmed antiSMASH 8.0.4 JSON (M2-18). Expected values come from the
recordings themselves and the independent oracles in ``helpers/regions_oracle.py``.

* Nissle: Bakta calls the three small contigs ``contig`` (kind ``unknown``); geNOMAD's plasmid
  summary lists CP171244.1 (3,173 bp, pMUT1) and CP171245.1 (5,514 bp, pMUT2) → exactly those two
  plasmid regions, ``source_tool = geNOMAD``; CP171243.1 is not a plasmid region.
* K-12: crispr_1 (array ``[2877883, 2878473)``) takes Cas2 (no gene name, product
  ``CRISPR-associated…``, gap 95) and cas1e (gap 381) but not cas6e (gap 1,314); crispr_2 has no
  cas gene within 1 kb.
* Full recordings under ``$GIPPA_ROOT/recordings`` (skipped when missing): the excerpt gives the
  same regions as the full Bakta JSON, and the ``region_layer`` step on the full Nissle run
  (``parse_bakta`` on the recorded genome) meets the acceptance criterion.
* Lead decision 4 after code review: Nissle has THREE plasmids (CP171243.1 4.5 kb, CP171244.1
  pMUT1, CP171245.1 pMUT2). The circular-replicon recording (``bakta_full_circ``, M2-10 lead
  decision 1) types all three ``plasmid`` → three plasmid regions, whole-contig, in FASTA order;
  CP171243.1 is Bakta-only, the other two Bakta+geNOMAD. The linear excerpt above calls
  CP171243.1 ``contig`` only for lack of topology tags.
"""

from __future__ import annotations

import hashlib
import importlib
import re
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from helpers.regions_oracle import (
    EXCERPTS,
    FULL,
    Sample,
    brute_members,
    contig_rank,
    exact_virus_gene,
    full_available,
    is_amr,
    load_excerpt,
    load_full,
    oracle_subset,
    overlapping_virus_genes,
)

from gippa.config import GippaConfig, resolve_config

regions = importlib.import_module("gippa.annotate.regions")  # M2-19


TAX = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"
K12_PROPHAGES = [  # 0-based half-open, numbered by start
    (563_847, 584_430, 0.8940),
    (1_196_866, 1_213_107, 0.9605),
    (1_411_999, 1_434_224, 0.9628),
    (1_627_516, 1_656_149, 0.8863),
    (2_463_011, 2_476_510, 0.9555),
]
NISSLE_PROPHAGES = [
    (1_694_617, 1_706_222, 0.9626),
    (2_171_336, 2_214_647, 0.9785),
    (2_869_018, 2_912_533, 0.9638),
]
SOURCE_SHA256 = {
    "k12_mg1655": "13c54a8961651053cda596ffa4cc0c611565588112275faf7f3ddc27e28bacec",
    "nissle1917": "49dbc4dc51e38443e822bf26481a7771c2d95f59315c5cb648cd23808809256f",
}


def _build(s: Sample, *overrides: str):
    cfg = resolve_config(overrides=list(overrides))
    return regions.build_regions(
        list(s.genes), list(s.features), s.replicons, s.genomad, s.bgc, cfg
    )


@pytest.fixture(scope="module")
def k12() -> Sample:
    return load_excerpt("k12_mg1655")


@pytest.fixture(scope="module")
def nissle() -> Sample:
    return load_excerpt("nissle1917")


def _of(regs, t: str) -> list:
    return sorted(
        (r for r in regs if r.region_type == t), key=lambda r: int(r.region_id[len(t) + 1 :])
    )


def _check_prophages(s: Sample, regs, rgs, want: list[tuple[int, int, float]]) -> None:
    pro = _of(regs, "prophage")
    assert [(r.interval.start, r.interval.end) for r in pro] == [(a, b) for a, b, _ in want]
    assert [r.score for r in pro] == [pytest.approx(x) for _, _, x in want]
    rep = next(iter(s.replicons.values()))
    order = {g.locus_tag: k for k, g in enumerate(s.genes)}
    for r in pro:
        assert (r.source_tool, r.notes, r.contig) == ("geNOMAD", TAX, rep.id)
        members = brute_members(s.genes, r.interval, s.replicons[r.contig])
        assert members, r.region_id
        assert list(r.genes) == sorted(members, key=order.__getitem__)
        rows = {x.locus_tag: x for x in rgs if x.region_id == r.region_id}
        assert set(rows) == set(members)
        for g in (g for g in s.genes if g.locus_tag in members):
            assert rows[g.locus_tag].overlap_frac == pytest.approx(members[g.locus_tag])
            v = exact_virus_gene(g, s.genomad.virus_genes)
            if v is not None:
                assert rows[g.locus_tag].subset == oracle_subset(g.product, v.description), g
            elif not overlapping_virus_genes(g, s.genomad.virus_genes):
                assert rows[g.locus_tag].subset == oracle_subset(g.product, None), g
    subsets = {x.subset for x in rgs if x.region_id.startswith("prophage_")}
    assert {"integrase", "structural", "cargo"} <= subsets


def _check_bgcs(s: Sample, regs, rgs) -> None:
    bgcs = _of(regs, "bgc")
    rank = contig_rank(s.replicons)
    areas = sorted(s.bgc, key=lambda b: (rank[b.contig], b.interval.start))
    assert len(bgcs) == len(areas)
    for r, a in zip(bgcs, areas, strict=True):
        assert (r.interval.start, r.interval.end) == (a.interval.start, a.interval.end)
        assert (r.source_tool, r.score) == ("antiSMASH", None)
        assert r.notes == "products=" + ",".join(a.products)
        members = brute_members(s.genes, r.interval, s.replicons[r.contig])
        assert set(r.genes) == set(members)
        rows = {x.locus_tag: x for x in rgs if x.region_id == r.region_id}
        assert {t for t, x in rows.items() if x.core_biosynthetic} == set(a.core_genes)
        assert set(a.core_genes) <= set(r.genes)


def test_k12_regions(k12: Sample) -> None:
    regs, rgs = _build(k12)
    assert {r.region_type for r in regs} == {"prophage", "bgc", "crispr", "amr"}
    _check_prophages(k12, regs, rgs, K12_PROPHAGES)
    _check_bgcs(k12, regs, rgs)
    assert [r.region_id for r in _of(regs, "bgc")] == ["bgc_1", "bgc_2", "bgc_3"]

    c1, c2 = _of(regs, "crispr")
    assert (c1.interval.start, c1.interval.end) == (2_877_883, 2_879_772)
    assert c1.genes == ("KOIAIB_02937", "KOIAIB_02938")
    assert (c2.interval.start, c2.interval.end) == (2_904_013, 2_904_417)
    assert c2.genes == ()
    assert {c1.source_tool, c2.source_tool} == {"Bakta"}

    (amr,) = _of(regs, "amr")
    want = [g.locus_tag for g in k12.genes if is_amr(g)]
    assert len(want) == 11
    assert (amr.region_id, amr.contig, amr.interval) == ("amr_1", "U00096.3", None)
    assert sorted(amr.genes) == sorted(want)
    assert amr.source_tool == "AMRFinderPlus"


def test_nissle_acceptance(nissle: Sample) -> None:
    regs, rgs = _build(nissle)
    assert len(_of(regs, "prophage")) >= 1 and len(_of(regs, "bgc")) >= 1
    _check_prophages(nissle, regs, rgs, NISSLE_PROPHAGES)
    _check_bgcs(nissle, regs, rgs)
    assert len(_of(regs, "bgc")) == 6

    plasmids = _of(regs, "plasmid")
    assert [(r.region_id, r.contig) for r in plasmids] == [
        ("plasmid_1", "CP171244.1"),  # pMUT1
        ("plasmid_2", "CP171245.1"),  # pMUT2
    ]
    for r, score in zip(plasmids, (0.9951, 0.9932), strict=True):
        assert (r.source_tool, r.score) == ("geNOMAD", pytest.approx(score))
        length = nissle.replicons[r.contig].length
        assert (r.interval.start, r.interval.end) == (0, length)
        on = [g.locus_tag for g in nissle.genes if g.interval.contig == r.contig]
        assert sorted(r.genes) == sorted(on) and on

    assert not _of(regs, "crispr")
    (amr,) = _of(regs, "amr")
    assert amr.contig == "CP171242.1" and len(amr.genes) == 37


def test_nissle_without_genomad_has_no_plasmid_or_prophage(nissle: Sample) -> None:
    """``modules`` without ``mobile``: Bakta alone calls no plasmid on the Nissle contigs."""
    cfg = resolve_config()
    regs, _ = regions.build_regions(
        list(nissle.genes), list(nissle.features), nissle.replicons, None, nissle.bgc, cfg
    )
    assert {r.region_type for r in regs} == {"bgc", "amr"}


def test_k12_tsv_and_bed_are_consistent(k12: Sample, tmp_path: Path) -> None:
    from gippa.io.bed import write_bed
    from gippa.io.tsv import read_tsv

    regs, _ = _build(k12)
    tsv = tmp_path / "regions.tsv"
    regions.write_regions_tsv(tsv, regs)
    _header, rows = read_tsv(tsv)
    by_id = {r.region_id: r for r in regs}
    assert {row[0] for row in rows} == set(by_id)
    for row in rows:
        r = by_id[row[0]]
        if r.interval is None:
            assert row[3:5] == [None, None]
        else:
            assert row[3:5] == [str(r.interval.start + 1), str(r.interval.end)]
        assert (row[7] or "") == ";".join(r.genes)
    bed = tmp_path / "regions.bed"
    write_bed(bed, regions.regions_bed(regs, k12.replicons), k12.replicons)
    names = [ln.split("\t")[3] for ln in bed.read_text().splitlines()]
    assert sorted(names) == sorted(r.region_id for r in regs if r.region_type != "amr")
    assert all(re.fullmatch(r"(prophage|plasmid|bgc|crispr)_\d+", n) for n in names)


def test_excerpt_provenance_names_the_source_hashes() -> None:
    text = (EXCERPTS / "PROVENANCE.md").read_text()
    for sha in SOURCE_SHA256.values():
        assert sha in text


# ----------------------------------------------------------------------------- full recordings


@pytest.mark.parametrize("sample", ["k12_mg1655", "nissle1917"])
def test_full_recordings_match_the_excerpts(sample: str) -> None:
    if not full_available(sample):
        pytest.skip(f"full recordings of {sample} not under {FULL}")
    src = FULL / "bakta_full" / sample / "bakta" / f"{sample}.json"
    assert hashlib.sha256(src.read_bytes()).hexdigest() == SOURCE_SHA256[sample]
    full = load_full(sample)
    got = _build(full)
    want = _build(load_excerpt(sample))
    assert got == want


CIRC = FULL / "bakta_full_circ" / "nissle1917"
CIRC_JSON = CIRC / "bakta" / "nissle1917.json"
CIRC_GENOME = CIRC / "preflight" / "work" / "preflight" / "genome.norm.fa"
NISSLE_PLASMIDS = [  # FASTA order; (contig, source_tool, geNOMAD plasmid_score)
    ("CP171243.1", "Bakta", None),
    ("CP171244.1", "Bakta+geNOMAD", 0.9951),  # pMUT1
    ("CP171245.1", "Bakta+geNOMAD", 0.9932),  # pMUT2
]


def _circ_available() -> bool:
    return full_available("nissle1917") and CIRC_JSON.exists() and CIRC_GENOME.exists()


def test_circular_nissle_recording_gives_three_plasmid_regions() -> None:
    if not _circ_available():
        pytest.skip(f"circular Nissle recording not under {CIRC}")
    from gippa.annotate.antismash import parse_antismash_json
    from gippa.annotate.bakta import parse_bakta_json
    from gippa.annotate.genomad import parse_genomad

    ann = parse_bakta_json(CIRC_JSON)
    reps = dict(ann.model.replicons)
    assert [(c, reps[c].kind) for c, _, _ in NISSLE_PLASMIDS] == [
        (c, "plasmid") for c, _, _ in NISSLE_PLASMIDS
    ]
    gm = parse_genomad(FULL / "genomad" / "nissle1917" / "genomad", "genome", reps)
    bgc = parse_antismash_json(FULL / "antismash" / "nissle1917" / "antismash" / "bakta.json", reps)
    genes = list(ann.model.genes)
    regs, _ = regions.build_regions(genes, list(ann.model.features), reps, gm, bgc,
                                    resolve_config())  # fmt: skip
    plasmids = _of(regs, "plasmid")
    assert [(r.region_id, r.contig) for r in plasmids] == [
        (f"plasmid_{k}", c) for k, (c, _, _) in enumerate(NISSLE_PLASMIDS, 1)
    ]
    for r, (c, tool, score) in zip(plasmids, NISSLE_PLASMIDS, strict=True):
        assert r.source_tool == tool, c
        assert r.score == (None if score is None else pytest.approx(score)), c
        assert (r.interval.start, r.interval.end) == (0, reps[c].length)
        on = [g.locus_tag for g in genes if g.interval.contig == c]
        assert on and sorted(r.genes) == sorted(on)
    assert len(_of(regs, "prophage")) >= 1 and len(_of(regs, "bgc")) >= 1


def test_region_layer_step_on_the_full_nissle_run(tmp_path: Path) -> None:
    genome = CIRC_GENOME
    if not _circ_available():
        pytest.skip(f"circular Nissle recording not under {CIRC}")
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    parse_out = {
        k: tmp_path / v
        for k, v in {
            "genes": "genes.parquet", "features": "features.parquet",
            "replicons": "replicons.parquet", "rrna16s": "rrna16s.fa",
            "trna_counts": "trna_counts.tsv", "xrefs": "xrefs.parquet",
            "annotation_meta": "annotation_meta.json",
        }.items()
    }  # fmt: skip
    bakta = CIRC_JSON
    reg.run_step("parse_bakta", GippaConfig(seed=1), {"bakta_json": bakta, "genome": genome},
                 parse_out)  # fmt: skip
    outs = {"regions": tmp_path / "regions.parquet", "region_genes": tmp_path / "rg.parquet"}
    inputs = {
        "genes": parse_out["genes"],
        "features": parse_out["features"],
        "replicons": parse_out["replicons"],
        "genomad": FULL / "genomad" / "nissle1917" / "genomad",
        "antismash": FULL / "antismash" / "nissle1917" / "antismash",
    }
    reg.run_step("region_layer", GippaConfig(seed=1), inputs, outs)
    rows = pq.read_table(outs["regions"]).to_pylist()
    types = [r["region_type"] for r in rows]
    assert types.count("prophage") >= 1 and types.count("bgc") >= 1
    plasmids = sorted(r["contig"] for r in rows if r["region_type"] == "plasmid")
    assert plasmids == [c for c, _, _ in NISSLE_PLASMIDS]
