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
rnaQ        4001-4100       ... but RegulonDB-only, so it blocks selA-selB     — degenerate
selB        4201-4500       TU ``TU03``
solo        5001-5300       TU ``TU04`` of one gene: contributes no positive   — degenerate
revX        5401-5700 (-)   opposite strand: the solo-revX pair is not same-strand
==========  ==============  =============================================================
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest

from giap.model.coords import Interval, Replicon, Strand
from giap.model.genes import Gene, GeneModel

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
    assert mapping.locus_by_name["lacZ"] == "T0001"
    assert mapping.how["lacZ"] == "coordinates"
    assert mapping.unmapped == ("rnaQ",)  # the RegulonDB-only gene
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
    assert mapping.locus_by_name == {"lacZ": "T0001"}
    assert mapping.unmapped == ("farAway",)


def test_map_genes_falls_back_to_a_unique_gene_name() -> None:
    model = _model()
    moved = (
        "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\n"
        "RDB1\tlacZ\tb1\t4801\t5000\tforward\n"  # nowhere near T0001, but the name is unique
    )
    mapping = rt.map_genes(rt.parse_gene_set(moved, CONTIG), model, CONTIG)
    assert mapping.locus_by_name == {"lacZ": "T0001"}
    assert mapping.how == {"lacZ": "name"}


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
    assert mapping.locus_by_name == {"lacZ": "T0001"}  # the closer one wins
    assert mapping.unmapped == ("lacZprime",)


def test_chromosome_contig_picks_the_chromosome() -> None:
    assert rt.chromosome_contig(_model()) == CONTIG


# ------------------------------------------------------- §15.4 step 2 labels


def test_adjacent_same_strand_pairs_skips_the_strand_change(built: dict[str, Any]) -> None:
    pairs = rt.adjacent_same_strand_pairs(built["model"], CONTIG)
    assert ("T0009", "T0010") not in [(a.locus_tag, b.locus_tag) for a, b in pairs]
    assert [(a.locus_tag, b.locus_tag) for a, b in pairs][:3] == [
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


def test_degenerate_a_regulondb_only_gene_between_two_genes_is_ambiguous(
    built: dict[str, Any],
) -> None:
    # selA and selB share TU03, but rnaQ (RegulonDB only) lies between them: the degenerate
    # table of the spec states "ambiguous" without qualification, so it wins over the TU.
    row = built["by_pair"][("T0007", "T0008")]
    assert row.label == "ambiguous"
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


def test_build_reads_a_release_directory_and_reproduces_the_labels(tmp_path: Path) -> None:
    """``build`` is the only file-reading entry point: paths in, values out, no network."""
    release = tmp_path / "release"
    release.mkdir()
    (release / "GeneProductSet.txt").write_text(GENE_SET, encoding="utf-8")
    (release / "TUSet.txt").write_text(TU_SET, encoding="utf-8")
    (release / "TerminatorSet.txt").write_text(TERMINATOR_SET, encoding="utf-8")
    bakta = tmp_path / "k12.json"
    _bakta_json(bakta)

    pairs, terminators, summary = rt.build(release, bakta)

    assert summary["contig"] == CONTIG
    labels = {(p.gene_a, p.gene_b): p.label for p in pairs}
    assert labels[("T0001", "T0002")] == "positive"
    assert labels[("T0003", "T0004")] == "negative"
    assert labels[("T0007", "T0008")] == "ambiguous"
    assert [t.source_id for t in terminators] == ["TM01", "TM02"]
    assert summary["pairs_total"] == len(pairs)
    assert summary["mapped_share"] == pytest.approx(0.9)


def test_main_writes_the_three_outputs(tmp_path: Path) -> None:
    release = tmp_path / "release"
    release.mkdir()
    (release / "GeneProductSet.txt").write_text(GENE_SET, encoding="utf-8")
    (release / "TUSet.txt").write_text(TU_SET, encoding="utf-8")
    (release / "TerminatorSet.txt").write_text(TERMINATOR_SET, encoding="utf-8")
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
