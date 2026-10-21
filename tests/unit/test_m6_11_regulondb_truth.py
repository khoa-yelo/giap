"""K-12 pair and terminator truth from RegulonDB (M6-11; GIAP_SPEC.md §15.4 steps 1-2).

Nothing here reaches the network: the fixtures below are RegulonDB-*shaped* tables we wrote by
hand (never RegulonDB bytes, which are not committed — R6/R11), and the gene model is built
directly from :class:`giap.model.genes.Gene` objects.

The fixture genome, all on ``+`` unless said otherwise (1-based closed, as RegulonDB writes them):

==========  ==============  =============================================================
gene        1-based span    role
==========  ==============  =============================================================
lacZ        1001-1300       TU ``TU01`` = (lacZ, lacY, lacA)                   — example 1
lacY        1401-1700       "
lacA        1801-2100       "
nxtB        2201-2500       TU ``TU02``; no TU holds nxtB with lacA            — example 2
orfC        2601-2900       in no TU at all                                    — example 3
ghoX        3001-3300       Bakta only, absent from RegulonDB                  — example 4
selA        3601-3900       TU ``TU03`` with selB ...
rnaQ        4001-4100       ... RegulonDB-only between them: still positive    — decision 3
selB        4201-4500       TU ``TU03``
solo        5001-5300       TU ``TU04`` of one gene: contributes no positive   — degenerate
revX        5401-5700 (-)   opposite strand: the solo-revX pair is not same-strand
==========  ==============  =============================================================
"""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.operons import adjacent_pairs
from giap.model.coords import Interval, Replicon, Strand
from giap.model.genes import Gene, GeneModel, NcFeature

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
rt: Any = importlib.import_module("validation.regulondb_truth")

CONTIG = "c1"

# RegulonDB gene name -> (left, right, strand, bakta locus tag or None when Bakta misses it)
RDB_GENES: dict[str, tuple[int, int, str, str | None]] = {
    "lacZ": (1001, 1300, "forward", "T0001"),
    "lacY": (1401, 1700, "forward", "T0002"),
    "lacA": (1801, 2100, "forward", "T0003"),
    "nxtB": (2201, 2500, "forward", "T0004"),
    "orfC": (2601, 2900, "forward", "T0005"),
    "selA": (3601, 3900, "forward", "T0007"),
    "rnaQ": (4001, 4100, "forward", None),
    "selB": (4201, 4500, "forward", "T0008"),
    "solo": (5001, 5300, "forward", "T0009"),
    "revX": (5401, 5700, "reverse", "T0010"),
}
# Bakta-only: ghoX between orfC and selA (example 4).
BAKTA_ONLY: dict[str, tuple[int, int, Strand, str]] = {"T0006": (3001, 3300, "+", "ghoX")}

TU_ROWS = (
    ("TU01", "lacZYA", "OP01", "lacZYA", "lacZ;lacY;lacA", "TM01"),
    ("TU02", "nxtB", "OP02", "nxtB", "nxtB", ""),
    ("TU03", "selAB", "OP03", "selAB", "selA;selB", ""),
    ("TU04", "solo", "OP04", "solo", "solo", ""),
    # an overlapping TU holding lacY and lacA only: the degenerate "two TUs share a pair" row
    ("TU05", "lacYA", "OP01", "lacZYA", "lacY;lacA", ""),
)

GENE_SET = (
    "\n".join(
        [
            "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\t7)productId",
            *(
                f"RDB{name}\t{name}\tb{i:04d}\t{left}\t{right}\t{strand}\tPD{i:04d}"
                for i, (name, (left, right, strand, _)) in enumerate(RDB_GENES.items(), start=1)
            ),
            # a row RegulonDB writes without coordinates: dropped, never a crash
            "RDBnocoord\tphantom\tb9999\tNone\tNone\tNone\tPD9999",
        ]
    )
    + "\n"
)

TU_SET = (
    "\n".join(
        [
            "# Columns:",
            "# (1) Transcription Unit identifier assigned by RegulonDB",
            "1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes\t6)promoterId\t7)promoterName"
            "\t8)terminatorIds",
            *(
                f"{tu}\t{name}\t{op}\t{opname}\t{genes}\tPM{tu}\t{name}p\t{terms}"
                for tu, name, op, opname, genes, terms in TU_ROWS
            ),
        ]
    )
    + "\n"
)

TERMINATOR_SET = (
    "\n".join(
        [
            "1)id\t2)leftEndPos\t3)rightEndPos\t4)strand\t5)sequence\t6)relatedTus\t7)type",
            "TM01\t2150\t2180\tforward\tacgt\tlacZYA:lacZp|lacZYA:lacZp2\trho-independent",
            "TM02\t4600\t4640\treverse\ttgca\tselAB:selAp\trho-dependent",
            "TM03\tNone\tNone\tforward\t\t\t",
        ]
    )
    + "\n"
)


def _model() -> GeneModel:
    genes = []
    for name, (left, right, strand, locus) in RDB_GENES.items():
        if locus is None:
            continue
        genes.append(
            Gene(
                locus_tag=locus,
                feature_type="cds",
                interval=Interval(CONTIG, left - 1, right, "+" if strand == "forward" else "-"),
                aa="M" * ((right - left + 1) // 3),
                nt=None,
                start_codon=None,
                gene=name,
            )
        )
    for locus, (left, right, strand, name) in BAKTA_ONLY.items():
        genes.append(
            Gene(
                locus_tag=locus,
                feature_type="cds",
                interval=Interval(CONTIG, left - 1, right, strand),
                aa="M" * ((right - left + 1) // 3),
                nt=None,
                start_codon=None,
                gene=name,
            )
        )
    replicons = {CONTIG: Replicon(CONTIG, 6000, "circular", "chromosome")}
    return GeneModel(replicons=replicons, genes=tuple(genes), features=())


def _locus_of(mapping: Any, name: str) -> str | None:
    """The Bakta locus a RegulonDB gene mapped to, whether the mapping is keyed by name or id."""
    got = mapping.locus_by_name.get(f"RDB{name}")
    return got if got is not None else mapping.locus_by_name.get(name)


def _how_of(mapping: Any, name: str) -> str | None:
    got = mapping.how.get(f"RDB{name}")
    return got if got is not None else mapping.how.get(name)


def _unmapped_names(mapping: Any) -> tuple[str, ...]:
    """``mapping.unmapped`` as RegulonDB gene names (fixture gene ids are ``RDB<name>``)."""
    return tuple(sorted(u.removeprefix("RDB") for u in mapping.unmapped))


@pytest.fixture(scope="module")
def built() -> dict[str, Any]:
    model = _model()
    rdb_genes = rt.parse_gene_set(GENE_SET, CONTIG)
    tus = rt.parse_tu_set(TU_SET)
    mapping = rt.map_genes(rdb_genes, model, CONTIG)
    pairs = rt.pair_truth(model, tus, rdb_genes, mapping, CONTIG)
    return {
        "model": model,
        "rdb_genes": rdb_genes,
        "tus": tus,
        "mapping": mapping,
        "pairs": pairs,
        "by_pair": {(p.gene_a, p.gene_b): p for p in pairs},
    }


# ------------------------------------------------------------------ parsing


def test_parse_rif_strips_the_numbered_header_and_comments() -> None:
    header, rows = rt.parse_rif(TU_SET)
    assert header[:5] == ("id", "name", "operonId", "operonName", "tuGenes")
    assert len(rows) == len(TU_ROWS)
    assert rows[0]["tuGenes"] == "lacZ;lacY;lacA"


def test_parse_rif_rejects_a_row_longer_than_the_header() -> None:
    with pytest.raises(rt.RegulonDBFormatError, match="header"):
        rt.parse_rif("1)a\t2)b\nx\ty\tz\n")


def test_parse_gene_set_converts_to_half_open_and_drops_coordinateless_rows() -> None:
    genes = rt.parse_gene_set(GENE_SET, CONTIG)
    assert [g.name for g in genes] == list(RDB_GENES)  # "phantom" dropped
    lacz = next(g for g in genes if g.name == "lacZ")
    assert (lacz.interval.start, lacz.interval.end, lacz.interval.strand) == (1000, 1300, "+")
    revx = next(g for g in genes if g.name == "revX")
    assert revx.interval.strand == "-"


def test_parse_gene_set_rejects_an_unknown_strand() -> None:
    bad = (
        "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\n"
        "G\tg\tb\t1\t9\tup\n"
    )
    with pytest.raises(rt.RegulonDBFormatError, match="unknown strand"):
        rt.parse_gene_set(bad, CONTIG)


def test_parse_tu_set_splits_genes_and_terminators() -> None:
    tus = {t.tu_id: t for t in rt.parse_tu_set(TU_SET)}
    assert tus["TU01"].genes == ("lacZ", "lacY", "lacA")
    assert tus["TU01"].terminator_ids == ("TM01",)
    assert tus["TU02"].terminator_ids == ()


def test_parse_terminator_set_converts_coordinates_and_drops_empty_rows() -> None:
    terms = rt.parse_terminator_set(TERMINATOR_SET, CONTIG)
    assert [t.terminator_id for t in terms] == ["TM01", "TM02"]
    assert (terms[0].interval.start, terms[0].interval.end) == (2149, 2180)
    assert terms[1].interval.strand == "-"
    assert terms[0].related_tus == ("lacZYA:lacZp", "lacZYA:lacZp2")


def test_parse_rif_rejects_an_empty_table() -> None:
    with pytest.raises(rt.RegulonDBFormatError):
        rt.parse_rif("# only a comment\n")


# ------------------------------------------------------------------ mapping


def test_map_genes_matches_by_coordinates_and_reports_the_share(built: dict[str, Any]) -> None:
    mapping = built["mapping"]
    assert _locus_of(mapping, "lacZ") == "T0001"
    assert _how_of(mapping, "lacZ") == "coordinates"
    assert _unmapped_names(mapping) == ("rnaQ",)  # the RegulonDB-only gene
    assert mapping.n_mapped == len(RDB_GENES) - 1
    assert mapping.mapped_share == pytest.approx(9 / 10)


def test_map_genes_tolerates_a_shifted_start_but_not_a_distant_gene() -> None:
    model = _model()
    shifted = (
        "\n".join(
            [
                "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand",
                "RDB1\tlacZ\tb1\t1006\t1300\tforward",  # 5 bp in: still lacZ
                "RDB2\tfarAway\tb2\t1001\t1300\treverse",  # same span, other strand: no match
            ]
        )
        + "\n"
    )
    mapping = rt.map_genes(rt.parse_gene_set(shifted, CONTIG), model, CONTIG, tolerance_bp=10)
    assert list(mapping.locus_by_name.values()) == ["T0001"]
    assert _locus_of(mapping, "lacZ") == "T0001"
    assert _unmapped_names(mapping) == ("farAway",)


def test_map_genes_falls_back_to_a_unique_gene_name_that_overlaps_on_the_same_strand() -> None:
    model = _model()
    moved = (
        "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\n"
        # beyond the 30 bp tolerance at both ends, but overlapping T0001 (1001-1300) on +
        "RDBlacZ\tlacZ\tb1\t1101\t1500\tforward\n"
    )
    mapping = rt.map_genes(rt.parse_gene_set(moved, CONTIG), model, CONTIG)
    assert list(mapping.locus_by_name.values()) == ["T0001"]
    assert _how_of(mapping, "lacZ") == "name"


def test_map_genes_claims_each_bakta_gene_once() -> None:
    model = _model()
    twins = (
        "\n".join(
            [
                "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand",
                "RDB1\tlacZ\tb1\t1001\t1300\tforward",
                "RDB2\tlacZprime\tb2\t1002\t1300\tforward",
            ]
        )
        + "\n"
    )
    mapping = rt.map_genes(rt.parse_gene_set(twins, CONTIG), model, CONTIG)
    assert list(mapping.locus_by_name.values()) == ["T0001"]
    assert _locus_of(mapping, "lacZ") == "T0001"  # the closer one wins
    assert _unmapped_names(mapping) == ("lacZprime",)


def test_chromosome_contig_picks_the_chromosome() -> None:
    assert rt.chromosome_contig(_model()) == CONTIG


# ------------------------------------------------------- §15.4 step 2 labels


def test_pair_truth_skips_the_strand_change(built: dict[str, Any]) -> None:
    pairs = [(p.gene_a, p.gene_b) for p in built["pairs"]]
    assert ("T0009", "T0010") not in pairs
    assert ("T0010", "T0001") not in pairs  # the wrap pair of this ring changes strand
    assert pairs[:3] == [
        ("T0001", "T0002"),
        ("T0002", "T0003"),
        ("T0003", "T0004"),
    ]


def test_example_1_genes_of_one_tu_are_positive(built: dict[str, Any]) -> None:
    assert built["by_pair"][("T0001", "T0002")].label == "positive"
    assert built["by_pair"][("T0002", "T0003")].label == "positive"


def test_overlapping_tus_still_give_one_positive_with_both_ids(built: dict[str, Any]) -> None:
    row = built["by_pair"][("T0002", "T0003")]  # lacY-lacA: in TU01 and TU05
    assert row.label == "positive"
    assert row.tu_ids == ("TU01", "TU05")


def test_example_2_different_tus_and_no_tu_with_both_are_negative(built: dict[str, Any]) -> None:
    row = built["by_pair"][("T0003", "T0004")]  # lacA (TU01, TU05) and nxtB (TU02)
    assert row.label == "negative"
    assert row.tu_ids == ("TU01", "TU02", "TU05")


def test_example_3_a_gene_in_no_tu_is_ambiguous(built: dict[str, Any]) -> None:
    row = built["by_pair"][("T0004", "T0005")]  # nxtB and orfC, which is in no TU
    assert row.label == "ambiguous"


def test_example_4_every_pair_touching_a_bakta_only_gene_is_ambiguous(
    built: dict[str, Any],
) -> None:
    touching = [p for p in built["pairs"] if "T0006" in (p.gene_a, p.gene_b)]
    assert len(touching) == 2
    assert {p.label for p in touching} == {"ambiguous"}


def test_decision_3_one_tu_is_positive_even_across_a_regulondb_only_gene(
    built: dict[str, Any],
) -> None:
    # selA and selB share TU03 and rnaQ (RegulonDB only) lies between them. Lead decision 3
    # (2026-10-21) supersedes the degenerate row: §15.4 says "both genes lie in one TU" ->
    # positive (the real data had dnaK-dnaJ, ftsA-ftsZ and pdhR-aceE in exactly this shape).
    row = built["by_pair"][("T0007", "T0008")]
    assert row.label == "positive"
    assert row.tu_ids == ("TU03",)


def test_degenerate_a_tu_of_one_gene_contributes_no_positive(built: dict[str, Any]) -> None:
    solo_pairs = [p for p in built["pairs"] if "T0009" in (p.gene_a, p.gene_b)]
    assert all(p.label != "positive" for p in solo_pairs)


def test_label_pair_is_pure_and_total(built: dict[str, Any]) -> None:
    assert {p.label for p in built["pairs"]} <= {"positive", "negative", "ambiguous"}
    assert len(built["pairs"]) == len({(p.gene_a, p.gene_b) for p in built["pairs"]})


def test_tus_by_locus_keeps_every_tu_of_a_gene(built: dict[str, Any]) -> None:
    by_locus = rt.tus_by_locus(built["tus"], built["mapping"])
    assert by_locus["T0002"] == frozenset({"TU01", "TU05"})
    assert "T0006" not in by_locus  # the Bakta-only gene is in no RegulonDB TU


# ------------------------------------------------------------- the writers


def test_terminator_truth_is_sorted_and_half_open() -> None:
    rows = rt.terminator_truth(rt.parse_terminator_set(TERMINATOR_SET, CONTIG))
    assert [r.source_id for r in rows] == ["TM01", "TM02"]
    assert (rows[0].contig, rows[0].start, rows[0].end, rows[0].strand) == (CONTIG, 2149, 2180, "+")


def test_pair_tsv_has_the_spec_columns(built: dict[str, Any]) -> None:
    text = rt.pair_tsv(built["pairs"])
    lines = text.splitlines()
    assert lines[0].split("\t") == list(rt.PAIR_COLUMNS)
    assert len(lines) == len(built["pairs"]) + 1
    assert "T0001\tT0002\tpositive\tTU01" in text


def test_terminator_tsv_has_the_spec_columns() -> None:
    rows = rt.terminator_truth(rt.parse_terminator_set(TERMINATOR_SET, CONTIG))
    lines = rt.terminator_tsv(rows).splitlines()
    assert lines[0].split("\t") == list(rt.TERMINATOR_COLUMNS)
    assert lines[1].split("\t") == [CONTIG, "2149", "2180", "+", "TM01"]


def test_summarise_counts_and_shares(built: dict[str, Any]) -> None:
    s = rt.summarise(built["pairs"], built["mapping"], 2)
    assert s["pairs_total"] == len(built["pairs"])
    assert s["pairs_labelled"] == s["pairs_positive"] + s["pairs_negative"]
    assert s["pairs_total"] == s["pairs_labelled"] + s["pairs_ambiguous"]
    assert s["excluded_share"] == pytest.approx(s["pairs_ambiguous"] / s["pairs_total"])
    assert s["mapped_share"] == pytest.approx(0.9)
    assert s["mapped_share_ok"] is True
    assert s["terminators"] == 2


def test_summarise_flags_a_low_mapped_share() -> None:
    mapping = rt.GeneMapping(locus_by_name={"a": "T1"}, how={"a": "name"}, unmapped=("b", "c"))
    assert rt.summarise((), mapping, 0)["mapped_share_ok"] is False


def _bakta_json(path: Path) -> None:
    """A minimal Bakta JSON holding the fixture genes, so ``build`` can be exercised offline."""
    features = [
        {
            "type": "cds",
            "locus": locus,
            "sequence": CONTIG,
            "start": left,
            "stop": right,
            "strand": strand,
            "gene": name,
            "aa": "M" * ((right - left + 1) // 3),
        }
        for locus, (left, right, strand, name) in (
            {
                locus: (left, right, "+" if st == "forward" else "-", name)
                for name, (left, right, st, locus) in RDB_GENES.items()
                if locus is not None
            }
            | BAKTA_ONLY
        ).items()
    ]
    path.write_text(
        json.dumps(
            {
                "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
                "sequences": [
                    {
                        "id": CONTIG,
                        "length": 6000,
                        "type": "chromosome",
                        "topology": "circular",
                        "nt": "A" * 6000,
                    }
                ],
                "features": features,
            }
        ),
        encoding="utf-8",
    )


def _write_release(
    release: Path,
    gene_set: str = GENE_SET,
    tu_set: str = TU_SET,
    terminator_set: str = TERMINATOR_SET,
    version: str = "14.5.0",
) -> Path:
    """A release directory as ``fetch_regulondb.py`` writes it: three tables + provenance.json."""
    release.mkdir(parents=True)
    files = []
    for dataset, text in (
        ("TUSet", tu_set),
        ("TerminatorSet", terminator_set),
        ("GeneProductSet", gene_set),
    ):
        data = text.encode("utf-8")
        (release / f"{dataset}.txt").write_bytes(data)
        files.append(
            {
                "dataset": dataset,
                "file": f"{dataset}.txt",
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        )
    meta = {"source": "RegulonDB", "regulondb_version": version, "files": files}
    (release / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return release


def test_build_reads_a_release_directory_and_reproduces_the_labels(tmp_path: Path) -> None:
    """``build`` is the only file-reading entry point: paths in, values out, no network."""
    release = _write_release(tmp_path / "release")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)

    pairs, terminators, summary = rt.build(release, bakta)

    assert summary["contig"] == CONTIG
    labels = {(p.gene_a, p.gene_b): p.label for p in pairs}
    assert labels[("T0001", "T0002")] == "positive"
    assert labels[("T0003", "T0004")] == "negative"
    assert labels[("T0007", "T0008")] == "positive"  # lead decision 3
    assert [t.source_id for t in terminators] == ["TM01", "TM02"]
    assert summary["pairs_total"] == len(pairs)
    assert summary["mapped_share"] == pytest.approx(0.9)


def test_main_writes_the_three_outputs(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)
    out = tmp_path / "truth"

    code = rt.main(
        [
            "--regulondb-dir",
            str(release),
            "--bakta-json",
            str(bakta),
            "--out-dir",
            str(out),
        ]
    )

    assert code == 0
    assert (out / "pair_truth.tsv").read_text().startswith("\t".join(rt.PAIR_COLUMNS))
    assert (out / "terminator_truth.tsv").read_text().splitlines()[0] == "\t".join(
        rt.TERMINATOR_COLUMNS
    )
    summary = json.loads((out / "truth_summary.json").read_text())
    assert summary["mapped_share_ok"] is True


# =================================================================== review fixes (2026-10-21)
#
# Small RegulonDB-shaped worlds built per test. Each gene is
# ``(name, left, right, rdb_strand, bakta)`` in 1-based closed coordinates, where ``bakta`` is
# ``("cds", locus)`` for a Bakta gene, ``("ncrna" | "trna" | ..., feature id)`` for a Bakta
# non-coding feature, or ``None`` for a gene only RegulonDB has. Gene ids are ``RDB<name>``.

MiniGene = tuple[str, int, int, str, tuple[str, str] | None]

_GENE_HEADER = "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand"


def _gene_set(rows: list[MiniGene], ids: dict[str, str] | None = None) -> str:
    lines = [_GENE_HEADER]
    for k, (name, left, right, strand, _) in enumerate(rows):
        gid = (ids or {}).get(f"{name}#{k}", f"RDB{name}")
        lines.append(f"{gid}\t{name}\tb{k:04d}\t{left}\t{right}\t{strand}")
    return "\n".join(lines) + "\n"


def _tu_set(tus: dict[str, str]) -> str:
    lines = ["1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes"]
    lines += [f"{tu}\t{tu}\tOP{tu}\t{tu}\t{genes}" for tu, genes in tus.items()]
    return "\n".join(lines) + "\n"


def _mini_model(
    rows: list[MiniGene], length: int = 10_000, topology: str = "circular"
) -> GeneModel:
    genes, features = [], []
    for name, left, right, strand, bakta in rows:
        if bakta is None:
            continue
        kind, ident = bakta
        iv = Interval(CONTIG, left - 1, right, "+" if strand == "forward" else "-")
        if kind == "cds":
            genes.append(
                Gene(
                    locus_tag=ident,
                    feature_type="cds",
                    interval=iv,
                    aa="M" * ((right - left + 1) // 3),
                    nt=None,
                    start_codon=None,
                    gene=name,
                )
            )
        else:
            features.append(NcFeature(id=ident, type=kind, interval=iv, gene=name))  # type: ignore[arg-type]
    replicons = {CONTIG: Replicon(CONTIG, length, topology, "chromosome")}  # type: ignore[arg-type]
    return GeneModel(replicons=replicons, genes=tuple(genes), features=tuple(features))


def _mini(
    rows: list[MiniGene],
    tus: dict[str, str],
    length: int = 10_000,
    topology: str = "circular",
    bakta_rows: list[MiniGene] | None = None,
) -> dict[str, Any]:
    """Run the pure pipeline on a mini world; ``bakta_rows`` overrides the Bakta side."""
    model = _mini_model(bakta_rows if bakta_rows is not None else rows, length, topology)
    rdb_genes = rt.parse_gene_set(_gene_set(rows), CONTIG)
    mapping = rt.map_genes(rdb_genes, model, CONTIG)
    pairs = rt.pair_truth(model, rt.parse_tu_set(_tu_set(tus)), rdb_genes, mapping, CONTIG)
    return {
        "model": model,
        "mapping": mapping,
        "pairs": pairs,
        "labels": {(p.gene_a, p.gene_b): p.label for p in pairs},
    }


def _label(world: dict[str, Any], a: str, b: str) -> str:
    labels = world["labels"]
    assert (a, b) in labels, f"pair ({a}, {b}) missing from the pair truth: {sorted(labels)}"
    return str(labels[(a, b)])


# ----------------------------------------- lead decision 1: verify before use (ADR-0017)


def test_decision_1_build_refuses_a_tampered_file(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    # same length, different bytes: only the SHA-256 can catch it
    text = (release / "TUSet.txt").read_text(encoding="utf-8").replace("TU03", "TU33")
    (release / "TUSet.txt").write_text(text, encoding="utf-8")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)

    with pytest.raises(ValueError, match="TUSet"):
        rt.build(release, bakta)


def test_decision_1_build_refuses_a_file_whose_size_changed(tmp_path: Path) -> None:
    release = _write_release(tmp_path / "release")
    with (release / "GeneProductSet.txt").open("a", encoding="utf-8") as fh:
        fh.write("RDBextra\textra\tb7777\t5801\t5900\tforward\n")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)

    with pytest.raises(ValueError, match="GeneProductSet"):
        rt.build(release, bakta)


def test_decision_1_build_refuses_a_release_without_its_provenance_sidecar(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    (release / "provenance.json").unlink()
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)

    with pytest.raises((ValueError, OSError), match="provenance"):
        rt.build(release, bakta)


def test_decision_1_truth_summary_carries_the_release_and_the_three_hashes(
    tmp_path: Path,
) -> None:
    release = _write_release(tmp_path / "release")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)
    out = tmp_path / "truth"

    code = rt.main(
        ["--regulondb-dir", str(release), "--bakta-json", str(bakta), "--out-dir", str(out)]
    )

    assert code == 0
    summary = json.loads((out / "truth_summary.json").read_text(encoding="utf-8"))
    assert "14.5.0" in summary.values(), "the RegulonDB release is not in truth_summary.json"
    dumped = json.dumps(summary)
    for text in (GENE_SET, TU_SET, TERMINATOR_SET):
        assert hashlib.sha256(text.encode("utf-8")).hexdigest() in dumped


# --------------------------------------- lead decision 2: the pipeline's adjacency (M2-12)

# A 10 kb ring: A at the start, B in the middle, C at the end; all on +.
_RING: list[MiniGene] = [
    ("aaaA", 101, 400, "forward", ("cds", "T1")),
    ("bbbB", 3001, 3300, "forward", ("cds", "T2")),
    ("cccC", 9501, 9800, "forward", ("cds", "T3")),
]


def test_decision_2_the_circular_wrap_pair_is_in_the_pair_truth() -> None:
    world = _mini(_RING, {"TUw": "cccC;aaaA", "TUb": "bbbB"})
    assert _label(world, "T3", "T1") == "positive"  # an operon spanning the origin


def test_decision_2_a_linear_contig_has_no_wrap_pair() -> None:
    world = _mini(_RING, {"TUw": "cccC;aaaA", "TUb": "bbbB"}, topology="linear")
    assert ("T3", "T1") not in world["labels"]


def test_decision_2_pairs_are_exactly_the_pipelines_adjacent_pairs(built: dict[str, Any]) -> None:
    model = built["model"]
    expected = [(a.locus_tag, b.locus_tag) for a, b in adjacent_pairs(model.genes, model.replicons)]
    assert [(p.gene_a, p.gene_b) for p in built["pairs"]] == expected


@pytest.mark.parametrize(
    ("left", "right"),
    [(9901, 9990), (21, 80)],
    ids=["before-the-origin", "after-the-origin"],
)
def test_decision_2_the_regulondb_only_gap_check_works_across_the_origin(
    left: int, right: int
) -> None:
    tus = {"TUab": "aaaA;bbbB", "TUc": "cccC"}
    plain = _mini(_RING, tus)
    assert _label(plain, "T3", "T1") == "negative"  # different TUs, nothing in between

    blocked = _mini([*_RING, ("rnaW", left, right, "forward", None)], tus)
    assert _label(blocked, "T3", "T1") == "ambiguous"


def test_decision_2_a_regulondb_only_gene_elsewhere_does_not_block_the_wrap_pair() -> None:
    tus = {"TUab": "aaaA;bbbB", "TUc": "cccC"}
    world = _mini([*_RING, ("rnaW", 5001, 5100, "forward", None)], tus)
    assert _label(world, "T3", "T1") == "negative"
    assert _label(world, "T2", "T3") == "ambiguous"  # this is the pair it lies inside


@settings(max_examples=60, deadline=None)
@given(
    strands=st.lists(st.sampled_from(["forward", "reverse"]), min_size=2, max_size=12),
    circular=st.booleans(),
)
def test_decision_2_property_pair_truth_follows_adjacent_pairs(
    strands: list[str], circular: bool
) -> None:
    rows: list[MiniGene] = [
        (f"g{k:02d}", 101 + 500 * k, 400 + 500 * k, s, ("cds", f"L{k:02d}"))
        for k, s in enumerate(strands)
    ]
    world = _mini(
        rows,
        {"TUall": ";".join(r[0] for r in rows)},
        length=500 * len(rows) + 200,
        topology="circular" if circular else "linear",
    )
    model = world["model"]
    expected = [(a.locus_tag, b.locus_tag) for a, b in adjacent_pairs(model.genes, model.replicons)]
    assert [(p.gene_a, p.gene_b) for p in world["pairs"]] == expected
    assert {p.label for p in world["pairs"]} <= {"positive"}  # one TU holds every gene


# ------------------------------------------------- lead decision 3: §15.4 labels

_AB: list[MiniGene] = [
    ("aaaA", 101, 400, "forward", ("cds", "T1")),
    ("bbbB", 701, 1000, "forward", ("cds", "T2")),
]


def test_decision_3_one_tu_with_an_unmapped_gene_between_is_positive() -> None:
    world = _mini([*_AB, ("rnaX", 501, 600, "forward", None)], {"TU1": "aaaA;rnaX;bbbB"})
    assert _label(world, "T1", "T2") == "positive"


def test_decision_3_different_tus_with_an_unmapped_gene_between_is_ambiguous() -> None:
    tus = {"TU1": "aaaA", "TU2": "bbbB"}
    assert _label(_mini(_AB, tus), "T1", "T2") == "negative"  # the control
    world = _mini([*_AB, ("rnaX", 501, 600, "forward", None)], tus)
    assert _label(world, "T1", "T2") == "ambiguous"


def test_decision_3_a_minus_strand_pair_in_one_tu_is_positive() -> None:
    rows: list[MiniGene] = [
        ("aaaA", 101, 400, "reverse", ("cds", "T1")),
        ("bbbB", 701, 1000, "reverse", ("cds", "T2")),
        ("rnaX", 501, 600, "reverse", None),
    ]
    world = _mini(rows, {"TU1": "bbbB;rnaX;aaaA"})  # RegulonDB lists - strand TUs 5'->3'
    assert _label(world, "T1", "T2") == "positive"


def test_decision_3_an_opposite_strand_regulondb_only_gene_blocks_only_the_negative() -> None:
    blocker = ("revX", 501, 600, "reverse", None)
    negative = _mini([*_AB, blocker], {"TU1": "aaaA", "TU2": "bbbB"})
    assert _label(negative, "T1", "T2") == "ambiguous"
    positive = _mini([*_AB, blocker], {"TU1": "aaaA;bbbB"})
    assert _label(positive, "T1", "T2") == "positive"


# ------------------------------------- lead decision 4: RNA genes map onto Bakta features


@pytest.mark.parametrize("kind", ["ncrna", "trna", "rrna"])
def test_decision_4_a_regulondb_rna_gene_maps_onto_a_bakta_feature(kind: str) -> None:
    rows: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "forward", (kind, "F1"))]
    world = _mini(rows, {"TU1": "aaaA", "TU2": "bbbB"})
    mapping = world["mapping"]
    assert _locus_of(mapping, "rnaX") == "F1"
    assert _unmapped_names(mapping) == ()
    assert mapping.mapped_share == pytest.approx(1.0)
    # mapped, so it no longer makes the negative ambiguous
    assert _label(world, "T1", "T2") == "negative"


def test_decision_4_an_rna_feature_within_the_tolerance_maps() -> None:
    rdb: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "forward", None)]
    bakta: list[MiniGene] = [*_AB, ("rnaX", 531, 570, "forward", ("ncrna", "F1"))]
    world = _mini(rdb, {"TU1": "aaaA", "TU2": "bbbB"}, bakta_rows=bakta)
    assert _locus_of(world["mapping"], "rnaX") == "F1"


def test_decision_4_an_rna_feature_on_the_other_strand_does_not_map() -> None:
    rdb: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "forward", None)]
    bakta: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "reverse", ("ncrna", "F1"))]
    world = _mini(rdb, {"TU1": "aaaA", "TU2": "bbbB"}, bakta_rows=bakta)
    assert _unmapped_names(world["mapping"]) == ("rnaX",)
    assert _label(world, "T1", "T2") == "ambiguous"


def test_decision_4_a_crispr_array_is_not_an_rna_gene() -> None:
    rdb: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "forward", None)]
    bakta: list[MiniGene] = [*_AB, ("rnaX", 501, 600, "forward", ("crispr", "F1"))]
    world = _mini(rdb, {"TU1": "aaaA", "TU2": "bbbB"}, bakta_rows=bakta)
    assert _unmapped_names(world["mapping"]) == ("rnaX",)


def test_decision_4_the_summary_reports_the_cds_share_separately() -> None:
    world = _mini([*_AB, ("rnaX", 501, 600, "forward", ("ncrna", "F1"))], {"TU1": "aaaA"})
    summary = rt.summarise(world["pairs"], world["mapping"], 0)
    cds_keys = [k for k in summary if "cds" in k.lower() and "share" in k.lower()]
    assert cds_keys, f"no CDS mapped share in the summary: {sorted(summary)}"
    assert summary["mapped_share"] == pytest.approx(1.0)


# --------------------------------------------- lead decision 5: the name fallback


def test_decision_5_the_name_fallback_rejects_the_other_strand() -> None:
    rdb: list[MiniGene] = [("aaaA", 101, 500, "reverse", None)]  # overlaps T1, wrong strand
    bakta: list[MiniGene] = [("aaaA", 201, 400, "forward", ("cds", "T1"))]
    world = _mini(rdb, {"TU1": "aaaA"}, bakta_rows=bakta)
    assert _unmapped_names(world["mapping"]) == ("aaaA",)


def test_decision_5_the_name_fallback_rejects_a_gene_that_does_not_overlap() -> None:
    rdb: list[MiniGene] = [("aaaA", 5001, 5300, "forward", None)]
    bakta: list[MiniGene] = [("aaaA", 101, 400, "forward", ("cds", "T1"))]
    world = _mini(rdb, {"TU1": "aaaA"}, bakta_rows=bakta)
    assert _unmapped_names(world["mapping"]) == ("aaaA",)


def test_decision_5_the_name_fallback_accepts_a_one_bp_overlap_on_the_same_strand() -> None:
    rdb: list[MiniGene] = [("aaaA", 400, 900, "forward", None)]  # shares base 400 with T1
    bakta: list[MiniGene] = [("aaaA", 101, 400, "forward", ("cds", "T1"))]
    world = _mini(rdb, {"TU1": "aaaA"}, bakta_rows=bakta)
    assert _locus_of(world["mapping"], "aaaA") == "T1"
    assert _how_of(world["mapping"], "aaaA") == "name"


def test_decision_5_duplicate_regulondb_names_keep_every_copy() -> None:
    model = _mini_model(_AB)
    text = (
        _GENE_HEADER
        + "\n"
        + "RDBdup1\tdup\tb0001\t101\t400\tforward\n"
        + "RDBdup2\tdup\tb0002\t701\t1000\tforward\n"
    )
    mapping = rt.map_genes(rt.parse_gene_set(text, CONTIG), model, CONTIG)
    assert sorted(mapping.locus_by_name.values()) == ["T1", "T2"]
    assert mapping.n_mapped == 2
    assert mapping.unmapped == ()
    assert mapping.mapped_share == pytest.approx(1.0)


# --------------------------------------------------------- survived mutants


@pytest.mark.parametrize(("shift", "maps"), [(30, True), (31, False)])
@pytest.mark.parametrize("end", ["start", "end"])
def test_the_coordinate_tolerance_is_inclusive_at_30_bp(shift: int, maps: bool, end: str) -> None:
    left, right = (101 + shift, 400) if end == "start" else (101, 400 + shift)
    # a name Bakta does not carry, so the name fallback cannot rescue a 31 bp miss
    rdb: list[MiniGene] = [("alias", left, right, "forward", None)]
    bakta: list[MiniGene] = [("aaaA", 101, 400, "forward", ("cds", "T1"))]
    world = _mini(rdb, {}, bakta_rows=bakta)
    assert rt.DEFAULT_TOLERANCE_BP == 30
    assert (_locus_of(world["mapping"], "alias") == "T1") is maps


def test_a_zero_coordinate_is_rejected() -> None:
    with pytest.raises(rt.RegulonDBFormatError, match="invalid 1-based interval"):
        rt.parse_gene_set(_GENE_HEADER + "\nRDBz\tz\tb1\t0\t90\tforward\n", CONTIG)


def test_a_coordinate_of_one_is_accepted() -> None:
    (gene,) = rt.parse_gene_set(_GENE_HEADER + "\nRDBo\to\tb1\t1\t90\tforward\n", CONTIG)
    assert (gene.interval.start, gene.interval.end) == (0, 90)


def test_mapped_share_ok_holds_at_exactly_the_floor() -> None:
    mapping = rt.GeneMapping(
        locus_by_name={f"g{k}": f"T{k}" for k in range(4)},
        how={f"g{k}": "coordinates" for k in range(4)},
        unmapped=("g9",),
    )
    assert mapping.mapped_share == pytest.approx(0.8)
    assert rt.MIN_MAPPED_SHARE == 0.8
    assert rt.summarise((), mapping, 0)["mapped_share_ok"] is True


def test_main_returns_1_and_escalates_when_the_mapped_share_is_too_low(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    extra = "".join(
        f"RDBghost{k}\tghost{k}\tb8{k:03d}\t{5801 + 20 * k}\t{5810 + 20 * k}\tforward\n"
        for k in range(5)
    )
    release = _write_release(tmp_path / "release", gene_set=GENE_SET + extra)  # 9 / 15 = 0.6
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)
    out = tmp_path / "truth"

    code = rt.main(
        ["--regulondb-dir", str(release), "--bakta-json", str(bakta), "--out-dir", str(out)]
    )

    assert code == 1
    assert "escalate" in capsys.readouterr().err
    summary = json.loads((out / "truth_summary.json").read_text(encoding="utf-8"))
    assert summary["mapped_share_ok"] is False
