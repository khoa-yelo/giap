"""Locus classes, uniqueness and gene assignment (issue M3-14; GIPPA_SPEC.md §7.2.6, §7.2.7, §9.4).

Examples 1-15 and the edge cases of docs/specs/M3-14-loci-and-assignment.md on a hand-made model
of the circular ``chr1`` (L = 10,000), plus flags, uniqueness, determinism and a timing guard.

Conventions these tests pin down where the spec is silent (lead decisions in M3-index.md:
``multi`` is a placement flag (Q1), repeats inside one gene are ``unique`` (Q2), pseudogene
placements are assigned to the pseudogene (Q3)):

* Example 10's G0 ``[9900, 10200)`` wraps onto ``[0, 200)`` and would overlap G1 ``[100, 400)``,
  which contradicts examples 2 and 5 (``[131,161) +`` is in frame with G0; ``[40,70) +`` overlaps
  G0 out of frame). The origin examples therefore use their own model (G0 alone); examples 1-9 and
  11-15 use the model without G0.
* ``flags`` is a ``;``-separated string (as M3-10 writes it); the empty string or NA mean "no
  flag"; flags are joined in alphabetical order (lead decision 1). Flags already present on the
  input (``same_db``, ``low_complexity``) are kept.
* ``accepted`` rows carry M3-04 ``HIT_COLUMNS`` plus ``source``, ``placement_q`` and M3-10's
  ``flags``; tolerant hits are recognised by ``route == "tolerant"``. ``peptides`` carries
  ``peptide``, ``protein_ids`` (``;``-separated) and ``search_db``.
* ``annotate_peptide_hits`` returns one row per input row (rows are matched here by peptide and
  interval); output rows are sorted by (peptide, contig, start, strand) (lead decision 1);
  ``gene`` is NA unless set; ``novel_candidate_id`` is NA.
* Frameshift ``blocks`` are ``start-end;start-end`` in genome coordinates (M3-08), listed by
  ascending coordinate.
* The `nterm_extension` stop scan uses table 11 stops (TAA, TAG, TGA) read on the gene's strand;
  no stop is looked for inside the hit itself.
"""

from __future__ import annotations

import time
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest

from gippa.mapping.hits import HIT_COLUMNS
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, GeneModel, NcFeature

L = 10_000
COMP = str.maketrans("ACGT", "TGCA")
# Stops in all three phases on both strands every 22 nt: TAA at 0, 4, 8; TTA (= TAA on -) at
# 11, 15, 19. Background made of this unit never lets a hit reach a CDS start by accident.
WALL = "TAAATAAATAATTATTTATTTA"
FILLER = "GCA"  # no stop codon in any of the six frames, whatever the phase


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def background(n: int) -> list[str]:
    return list((WALL * (n // len(WALL) + 1))[:n])


def fill(seq: list[str], start: int, end: int) -> None:
    """Stop-free filler on ``[start, end)``, wrapping the origin when ``end > len(seq)``."""
    n = len(seq)
    for k, i in enumerate(range(start, end)):
        seq[i % n] = FILLER[k % 3]


def put(seq: list[str], start: int, s: str) -> None:
    for k, c in enumerate(s):
        seq[(start + k) % len(seq)] = c


def _default_aa(tag: str, n_aa: int) -> str:
    """Distinct proteins per locus tag (identical `aa` would make loci `multi_identical`)."""
    body = "".join("ACDEFGHIKLMNPQRSTVWY"[ord(c) % 20] for c in tag)
    return ("M" + body + "A" * n_aa)[: max(n_aa, len(body) + 1)]


def cds(tag: str, contig: str, start: int, end: int, strand: str, aa: str | None = None,
        kind: str = "cds", partial: bool = False) -> Gene:  # fmt: skip
    n_aa = (end - start) // 3 - 1
    return Gene(
        locus_tag=tag,
        feature_type=kind,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa=aa if aa is not None else _default_aa(tag, n_aa),
        nt=None,
        start_codon="ATG",
        partial=partial,
    )


# --- the spec's example model (without G0, see the module docstring) ---------------------------

G1, G2, S1, P1, G5 = "GSIM_00005", "GSIM_00010", "GSIM_00015", "GSIM_00020", "GSIM_00025"
G3, G4, GA, GB, G0 = "GSIM_00030", "GSIM_00035", "GSIM_00040", "GSIM_00045", "GSIM_00001"
SAME_AA = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ" * 3  # G3 and G4 encode the same protein

REPS = {"chr1": Replicon("chr1", L, "circular", "chromosome")}
GENES = (
    cds(G1, "chr1", 100, 400, "+"),
    cds(G2, "chr1", 1000, 1300, "-"),
    cds(S1, "chr1", 2000, 2075, "+", kind="sorf"),
    cds(P1, "chr1", 3000, 3600, "+", kind="pseudogene"),
    cds(G5, "chr1", 5100, 5400, "+"),
    cds(GA, "chr1", 6600, 6900, "+"),
    cds(GB, "chr1", 6601, 6901, "+"),  # overlaps GA, other frame
    cds(G3, "chr1", 7000, 7300, "+", aa=SAME_AA),
    cds(G4, "chr1", 8000, 8300, "-", aa=SAME_AA),
)
FEATURES = (
    NcFeature("R1", "rrna", Interval("chr1", 5000, 6542, "+")),
    NcFeature("T1", "trna", Interval("chr1", 9000, 9076, "-")),
    NcFeature("C1", "crispr", Interval("chr1", 9300, 9500, ".")),
)
MODEL = GeneModel(replicons=REPS, genes=GENES, features=FEATURES)


def _genome_a() -> dict[str, str]:
    s = background(L)
    fill(s, 40, 100)  # upstream of G1: no stop (examples 4, 5 counterfactual)
    fill(s, 1300, 1360)  # upstream of G2 on the minus strand
    return {"chr1": "".join(s)}


GENOME = _genome_a()


def _genome_b() -> dict[str, str]:
    s = list(GENOME["chr1"])
    put(s, 76, "TAA")  # example 5: in frame with G1 ((76 - 100) % 3 == 0)
    put(s, 1306, "TTA")  # TAA on the minus strand, in frame with G2 ((1300 - 1309) % 3 == 0)
    return {"chr1": "".join(s)}


GENOME_B = _genome_b()


def iv(start: int, end: int, strand: str, contig: str = "chr1") -> Interval:
    return Interval(contig, start, end, strand)  # type: ignore[arg-type]


def classify(i: Interval, gm: GeneModel = MODEL, genome: Mapping[str, str] = GENOME) -> Any:
    from gippa.mapping.loci import classify_locus

    return classify_locus(i, gm, genome)


def test_fixture_sanity() -> None:
    """The test genome itself: stop-free where the examples need it, stops where they must."""
    stops = {"TAA", "TAG", "TGA"}
    a = GENOME["chr1"]
    assert all(a[i : i + 3] not in stops for i in range(40, 98))
    assert all(revcomp(a[i : i + 3]) not in stops for i in range(1300, 1358))
    assert GENOME_B["chr1"][76:79] == "TAA" and revcomp(GENOME_B["chr1"][1306:1309]) == "TAA"
    for ph in range(3):  # walls really contain stops in every phase on both strands
        assert any(a[i : i + 3] in stops for i in range(4000 + ph, 4100, 3))
        assert any(revcomp(a[i : i + 3]) in stops for i in range(4000 + ph, 4100, 3))


# --- examples 1-10: classify_locus -----------------------------------------------------------


@pytest.mark.parametrize(
    ("hit", "expected"),
    [
        pytest.param(iv(130, 160, "+"), ("cds_inframe", G1), id="ex1-inframe"),
        pytest.param(iv(131, 161, "+"), ("cds_outframe", None), id="ex2-outframe"),
        pytest.param(iv(130, 160, "-"), ("antisense", None), id="ex3-antisense"),
        pytest.param(iv(70, 100, "+"), ("nterm_extension", G1), id="ex4-nterm"),
        pytest.param(iv(1270, 1300, "-"), ("cds_inframe", G2), id="ex6-minus-inframe"),
        pytest.param(iv(1300, 1330, "-"), ("nterm_extension", G2), id="ex6-minus-nterm"),
        pytest.param(iv(5130, 5160, "+"), ("rna_overlap", None), id="ex7-rna-priority"),
        pytest.param(iv(2003, 2030, "+"), ("sorf", S1), id="ex8-sorf"),
        pytest.param(iv(3100, 3130, "+"), ("pseudogene", P1), id="ex9-pseudogene"),
    ],
)
def test_examples_1_to_9(hit: Interval, expected: tuple[str, str | None]) -> None:
    cls, gene = classify(hit)
    assert (cls, gene) == expected


def test_example_5_stop_separates_upstream_hit() -> None:
    assert tuple(classify(iv(40, 70, "+"), genome=GENOME_B)) == ("intergenic", None)
    # the same hit without the TAA is an N-terminal extension: the stop is what decides
    assert tuple(classify(iv(40, 70, "+"), genome=GENOME)) == ("nterm_extension", G1)


def test_minus_strand_stop_separates_upstream_hit() -> None:
    assert tuple(classify(iv(1330, 1360, "-"), genome=GENOME)) == ("nterm_extension", G2)
    assert tuple(classify(iv(1330, 1360, "-"), genome=GENOME_B)) == ("intergenic", None)


# --- example 10 and origin crossing ----------------------------------------------------------

ORIGIN_MODEL = GeneModel(replicons=REPS, genes=(cds(G0, "chr1", 9900, 10200, "+"),), features=())


def test_example_10_origin_gene() -> None:
    assert tuple(classify(iv(20, 50, "+"), ORIGIN_MODEL)) == ("cds_inframe", G0)
    assert tuple(classify(iv(9990, 10020, "+"), ORIGIN_MODEL)) == ("cds_inframe", G0)
    assert tuple(classify(iv(9951, 9981, "+"), ORIGIN_MODEL)) == ("cds_inframe", G0)
    # out of frame after the shift: (10021 - 9900) % 3 == 1
    assert tuple(classify(iv(21, 51, "+"), ORIGIN_MODEL)) == ("cds_outframe", None)
    assert tuple(classify(iv(20, 50, "-"), ORIGIN_MODEL)) == ("antisense", None)


# --- other class rules and priorities --------------------------------------------------------


@pytest.mark.parametrize(
    ("hit", "expected"),
    [
        pytest.param(iv(85, 115, "+"), ("nterm_extension", G1), id="nterm-reaches-into-cds"),
        pytest.param(iv(1285, 1315, "-"), ("nterm_extension", G2), id="nterm-minus-into-cds"),
        pytest.param(iv(2004, 2031, "+"), ("cds_outframe", None), id="sorf-outframe"),
        pytest.param(iv(2003, 2030, "-"), ("antisense", None), id="sorf-antisense"),
        pytest.param(iv(3101, 3131, "+"), ("pseudogene", P1), id="pseudogene-any-frame"),
        pytest.param(iv(2990, 3020, "+"), ("pseudogene", P1), id="pseudogene-partial-overlap"),
        pytest.param(iv(3100, 3130, "-"), ("intergenic", None), id="pseudogene-opposite-strand"),
        pytest.param(iv(5130, 5160, "-"), ("rna_overlap", None), id="rna-either-strand"),
        pytest.param(iv(4971, 5001, "+"), ("rna_overlap", None), id="rna-1bp-overlap"),
        pytest.param(iv(4970, 5000, "+"), ("intergenic", None), id="rna-adjacent-no-overlap"),
        pytest.param(iv(9030, 9060, "+"), ("rna_overlap", None), id="trna"),
        pytest.param(iv(9330, 9360, "+"), ("intergenic", None), id="crispr-is-not-rna"),
        pytest.param(iv(4000, 4030, "+"), ("intergenic", None), id="intergenic"),
        pytest.param(iv(4000, 4030, "-"), ("intergenic", None), id="intergenic-minus"),
        pytest.param(iv(6700, 6730, "+"), ("cds_inframe", GB), id="overlapping-cds-B"),
        pytest.param(iv(6699, 6729, "+"), ("cds_inframe", GA), id="overlapping-cds-A"),
        pytest.param(iv(7030, 7060, "+"), ("cds_inframe", G3), id="g3"),
        pytest.param(iv(8240, 8270, "-"), ("cds_inframe", G4), id="g4-minus"),
        pytest.param(iv(8241, 8271, "-"), ("cds_outframe", None), id="g4-minus-outframe"),
    ],
)
def test_class_rules(hit: Interval, expected: tuple[str, str | None]) -> None:
    assert tuple(classify(hit)) == expected


# a hit of 30 nt that overlaps R1 by one base on the minus strand
def test_rna_overlap_one_base_minus() -> None:
    assert tuple(classify(iv(6541, 6571, "-"))) == ("rna_overlap", None)


def test_partial_gene_classified_like_others() -> None:
    gm = GeneModel(REPS, (cds(G1, "chr1", 100, 400, "+", partial=True),), ())
    assert tuple(classify(iv(130, 160, "+"), gm)) == ("cds_inframe", G1)
    assert tuple(classify(iv(70, 100, "+"), gm)) == ("nterm_extension", G1)


def test_frameshift_hit_inside_cds_uses_prefix_frame() -> None:
    """A frameshift hit is classified by its span; the frame comes from the prefix block."""
    from gippa.mapping.loci import annotate_peptide_hits

    rows = [
        hit_row("FSPLUS", 130, 190, "+", route="tolerant", n_diff=1,
                diff_types="frameshift", blocks="130-160;161-190"),
        hit_row("FSMINUS", 1200, 1270, "-", route="tolerant", n_diff=1,
                diff_types="frameshift", blocks="1200-1239;1240-1270"),
    ]  # fmt: skip
    out = by_key(annotate_peptide_hits(acc(rows), MODEL, GENOME, peps(rows), "related"))
    assert out[("FSPLUS", 130, 190, "+")]["locus_class"] == "cds_inframe"
    assert out[("FSPLUS", 130, 190, "+")]["gene"] == G1
    assert out[("FSMINUS", 1200, 1270, "-")]["locus_class"] == "cds_inframe"
    assert out[("FSMINUS", 1200, 1270, "-")]["gene"] == G2


# --- nterm scan across the origin and at linear contig ends ---------------------------------


def _edge_model() -> tuple[GeneModel, dict[str, str]]:
    n = 3000
    reps = {
        "c2": Replicon("c2", n, "circular", "plasmid"),
        "c3": Replicon("c3", n, "circular", "plasmid"),
        "c4": Replicon("c4", 600, "circular", "plasmid"),
        "lin": Replicon("lin", 600, "linear", "unknown"),
    }
    genes = (
        cds("GX", "c2", 30, 330, "+"),
        cds("GY", "c3", 2700, 2970, "-"),
        cds("GZ", "c4", 30, 330, "+"),
        cds("GL", "lin", 30, 330, "+"),
        cds("GM", "lin", 400, 570, "-"),
    )
    seqs = {}
    for c, ln in (("c2", n), ("c3", n), ("c4", 600), ("lin", 600)):
        s = background(ln)
        if c == "c2":
            fill(s, 2900, 3030)
        elif c == "c3":
            fill(s, 2970, 3090)
        else:
            fill(s, 540, 630)  # [540, 600) and [0, 30)
        seqs[c] = "".join(s)
    return GeneModel(reps, genes, ()), seqs


def test_nterm_scan_crosses_origin_plus() -> None:
    gm, g = _edge_model()
    assert tuple(classify(iv(2940, 2970, "+", "c2"), gm, g)) == ("nterm_extension", "GX")
    s = list(g["c2"])
    put(s, 9, "TAA")  # (3009 - 3030) % 3 == 0: an in-frame stop just after the origin
    g2 = {**g, "c2": "".join(s)}
    assert tuple(classify(iv(2940, 2970, "+", "c2"), gm, g2)) == ("intergenic", None)


def test_nterm_scan_crosses_origin_minus() -> None:
    gm, g = _edge_model()
    # hit at [30, 60) on -: shifted to [3030, 3060); (2970 - 3060) % 3 == 0
    assert tuple(classify(iv(30, 60, "-", "c3"), gm, g)) == ("nterm_extension", "GY")
    s = list(g["c3"])
    put(s, 2985, "TTA")  # TAA on -, (2970 - 2988) % 3 == 0
    g2 = {**g, "c3": "".join(s)}
    assert tuple(classify(iv(30, 60, "-", "c3"), gm, g2)) == ("intergenic", None)


def test_hit_spanning_origin_upstream_of_gene() -> None:
    gm, g = _edge_model()
    assert tuple(classify(iv(2988, 3018, "+", "c2"), gm, g)) == ("nterm_extension", "GX")


def test_linear_contig_never_wraps() -> None:
    gm, g = _edge_model()
    # contig start, nothing to scan: fine on a linear contig
    assert tuple(classify(iv(0, 30, "+", "lin"), gm, g)) == ("nterm_extension", "GL")
    # minus gene near the contig end: the hit reaches the end of the contig
    assert tuple(classify(iv(570, 600, "-", "lin"), gm, g)) == ("nterm_extension", "GM")
    # the same layout on a circular contig reaches GZ through the origin; on the linear one the
    # hit at the contig end is not upstream of GL
    assert tuple(classify(iv(570, 600, "+", "c4"), gm, g)) == ("nterm_extension", "GZ")
    assert tuple(classify(iv(570, 600, "+", "lin"), gm, g)) == ("intergenic", None)


# --- annotate_peptide_hits ------------------------------------------------------------------


def frame_of(start: int, end: int, strand: str, length: int = L) -> int:
    if strand == "+":
        return 1 + start % 3
    return -(1 + ((length - end) % length) % 3)


def hit_row(peptide: str, start: int, end: int, strand: str, *, contig: str = "chr1",
            route: str = "exact", n_diff: int = 0, diff_types: str = "",
            blocks: str | None = None, flags: str = "", source: str = "database",
            length: int = L) -> dict[str, Any]:  # fmt: skip
    return {
        "peptide": peptide,
        "query": peptide.replace("I", "L"),
        "search": "target",
        "decoy_type": None,
        "contig": contig,
        "frame_id": 0,
        "aa_start": 0,
        "aa_end": (end - start) // 3,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": frame_of(start, end, strand, length),
        "route": route,
        "n_diff": n_diff,
        "diff_types": diff_types,
        "blocks": blocks,
        "source": source,
        "placement_q": 0.001,
        "flags": flags,
    }


def acc(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=[*HIT_COLUMNS, "source", "placement_q", "flags"])
    for c in ("frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff"):
        df[c] = df[c].astype("int64")
    return df


def peps(rows: Iterable[Mapping[str, Any]], protein_ids: Mapping[str, str | None] | None = None,
         search_db: str = "other") -> pd.DataFrame:  # fmt: skip
    names = sorted({r["peptide"] for r in rows})
    ids = protein_ids or {}
    return pd.DataFrame(
        {
            "peptide": names,
            "protein_ids": [ids.get(p) for p in names],
            "search_db": [search_db] * len(names),
            "source": ["database"] * len(names),
            "is_decoy": [False] * len(names),
            "is_contaminant": [False] * len(names),
        }
    )


Key = tuple[str, int, int, str]


def by_key(df: pd.DataFrame) -> dict[Key, dict[str, Any]]:
    out: dict[Key, dict[str, Any]] = {}
    for r in df.to_dict("records"):
        k = (str(r["peptide"]), int(r["start"]), int(r["end"]), str(r["strand"]))
        assert k not in out, f"duplicate output row {k}"
        out[k] = r
    return out


def flags_of(x: object) -> set[str]:
    if x is None or (not isinstance(x, str) and pd.isna(x)):
        return set()
    return {f for f in str(x).split(";") if f}


def na(x: object) -> bool:
    return x is None or bool(pd.isna(x))


def annotate(rows: Sequence[Mapping[str, Any]], relation: str = "related",
             genome: Mapping[str, str] = GENOME, gm: GeneModel = MODEL,
             peptides: pd.DataFrame | None = None) -> dict[Key, dict[str, Any]]:  # fmt: skip
    from gippa.mapping.loci import annotate_peptide_hits

    p = peps(rows) if peptides is None else peptides
    out = annotate_peptide_hits(acc(rows), gm, genome, p, relation)  # type: ignore[arg-type]
    assert len(out) == len(rows)
    return by_key(out)


def test_output_columns_and_assignment() -> None:
    rows = [
        hit_row("PEPA", 130, 160, "+"),
        hit_row("PEPB", 131, 161, "+"),
        hit_row("PEPC", 70, 100, "+"),
        hit_row("PEPD", 3100, 3130, "+"),
        hit_row("PEPE", 2003, 2030, "+"),
        hit_row("PEPF", 5130, 5160, "+"),
    ]
    from gippa.mapping.loci import annotate_peptide_hits

    df = annotate_peptide_hits(acc(rows), MODEL, GENOME, peps(rows), "related")
    for col in ("locus_class", "gene", "uniqueness", "flags", "novel_candidate_id"):
        assert col in df.columns
    for col in acc(rows).columns:  # input columns are kept
        assert col in df.columns
    assert df["novel_candidate_id"].isna().all()
    out = by_key(df)
    expect = {
        "PEPA": ("cds_inframe", G1, set()),
        "PEPB": ("cds_outframe", None, set()),
        "PEPC": ("nterm_extension", G1, {"nterm_extension_evidence"}),
        "PEPD": ("pseudogene", P1, set()),
        "PEPE": ("sorf", S1, set()),
        "PEPF": ("rna_overlap", None, set()),
    }
    for k, r in out.items():
        cls, gene, fl = expect[k[0]]
        assert r["locus_class"] == cls
        assert na(r["gene"]) if gene is None else r["gene"] == gene, (k, r["gene"])
        assert flags_of(r["flags"]) == fl, (k, r["flags"])
        assert r["uniqueness"] == "unique"


def test_example_10_origin_spanning_flag() -> None:
    rows = [hit_row("ORIA", 20, 50, "+"), hit_row("ORIB", 9990, 10020, "+")]
    out = annotate(rows, gm=ORIGIN_MODEL)
    a, b = out[("ORIA", 20, 50, "+")], out[("ORIB", 9990, 10020, "+")]
    assert (a["locus_class"], a["gene"]) == ("cds_inframe", G0)
    assert "origin_spanning" not in flags_of(a["flags"])
    assert (b["locus_class"], b["gene"]) == ("cds_inframe", G0)
    assert "origin_spanning" in flags_of(b["flags"])


def test_origin_spanning_on_intergenic_locus() -> None:
    out = annotate([hit_row("ORIC", 9995, 10025, "-")])
    r = out[("ORIC", 9995, 10025, "-")]
    assert "origin_spanning" in flags_of(r["flags"])


def test_example_11_multi_identical() -> None:
    rows = [hit_row("IDENT", 7030, 7060, "+"), hit_row("IDENT", 8240, 8270, "-")]
    out = annotate(rows)
    assert {r["gene"] for r in out.values()} == {G3, G4}
    for r in out.values():
        assert r["uniqueness"] == "multi_identical"
        assert "multi_identical" in flags_of(r["flags"])
        assert "multi" not in flags_of(r["flags"])  # same class on both loci


def test_different_genes_are_shared() -> None:
    rows = [hit_row("TWOG", 130, 160, "+"), hit_row("TWOG", 1270, 1300, "-")]
    for r in annotate(rows).values():
        assert r["uniqueness"] == "shared"
        assert not flags_of(r["flags"]) & {"multi", "multi_identical"}


def test_identical_aa_but_one_locus_unassigned_is_shared() -> None:
    rows = [hit_row("HALF", 7030, 7060, "+"), hit_row("HALF", 8241, 8271, "-")]  # G4 outframe
    out = annotate(rows)
    for r in out.values():
        assert r["uniqueness"] == "shared"
        assert "multi" in flags_of(r["flags"])
        assert "multi_identical" not in flags_of(r["flags"])


def test_example_12_shared_and_multi() -> None:
    rows = [hit_row("MIX", 130, 160, "+"), hit_row("MIX", 4000, 4030, "+")]
    out = annotate(rows)
    assert out[("MIX", 130, 160, "+")]["locus_class"] == "cds_inframe"
    assert out[("MIX", 4000, 4030, "+")]["locus_class"] == "intergenic"
    assert na(out[("MIX", 4000, 4030, "+")]["gene"])
    for r in out.values():
        assert r["uniqueness"] == "shared"
        assert "multi" in flags_of(r["flags"])


def test_q2_repeat_inside_one_gene_is_unique() -> None:
    rows = [hit_row("REP", 130, 160, "+"), hit_row("REP", 250, 280, "+")]
    for r in annotate(rows).values():
        assert (r["locus_class"], r["gene"], r["uniqueness"]) == ("cds_inframe", G1, "unique")
        assert not flags_of(r["flags"]) & {"multi", "multi_identical"}


def test_q2_same_gene_different_classes_unique_but_multi() -> None:
    rows = [hit_row("EXT", 70, 100, "+"), hit_row("EXT", 250, 280, "+")]
    out = annotate(rows)
    assert {r["locus_class"] for r in out.values()} == {"nterm_extension", "cds_inframe"}
    for r in out.values():
        assert r["gene"] == G1
        assert r["uniqueness"] == "unique"
        assert "multi" in flags_of(r["flags"])


def test_two_intergenic_loci_are_shared() -> None:
    rows = [hit_row("TWOI", 4000, 4030, "+"), hit_row("TWOI", 4400, 4430, "-")]
    for r in annotate(rows).values():
        assert r["uniqueness"] == "shared"
        assert "multi" not in flags_of(r["flags"])  # both intergenic


def test_example_13_strain_difference() -> None:
    rows = [
        hit_row("TOLER", 130, 160, "+", route="tolerant", n_diff=1, diff_types="substitution"),
        hit_row("EXACT", 250, 280, "+"),
    ]
    rel = annotate(rows, relation="related")
    same = annotate(rows, relation="same")
    assert "strain_difference" in flags_of(rel[("TOLER", 130, 160, "+")]["flags"])
    assert "strain_difference" not in flags_of(same[("TOLER", 130, 160, "+")]["flags"])
    assert "strain_difference" not in flags_of(rel[("EXACT", 250, 280, "+")]["flags"])
    # tolerant hits keep their class and gene under both relations
    for out in (rel, same):
        r = out[("TOLER", 130, 160, "+")]
        assert (r["locus_class"], r["gene"]) == ("cds_inframe", G1)


def test_example_14_engine_protein_mismatch() -> None:
    rows = [hit_row("MISM", 130, 160, "+")]
    out = annotate(rows, peptides=peps(rows, {"MISM": "GSIM_00010"}, search_db="same"))
    assert "engine_protein_mismatch" in flags_of(out[("MISM", 130, 160, "+")]["flags"])


@pytest.mark.parametrize(
    ("ids", "search_db", "flagged"),
    [
        ("GSIM_00005", "same", False),
        ("GSIM_00010;GSIM_00005", "same", False),
        ("rev_GSIM_00010;GSIM_00005", "same", False),
        ("rev_GSIM_00005", "same", True),  # a decoy entry never vouches for the gene
        ("CON__GSIM_00005;GSIM_00010", "same", True),
        ("GSIM_00010", "other", False),  # only with same-strain input
    ],
)
def test_engine_protein_mismatch_rules(ids: str, search_db: str, flagged: bool) -> None:
    rows = [hit_row("MISM", 130, 160, "+")]
    out = annotate(rows, peptides=peps(rows, {"MISM": ids}, search_db=search_db))
    assert ("engine_protein_mismatch" in flags_of(out[("MISM", 130, 160, "+")]["flags"])) == flagged


def test_engine_protein_mismatch_needs_a_gene() -> None:
    rows = [hit_row("NOGENE", 4000, 4030, "+")]
    out = annotate(rows, peptides=peps(rows, {"NOGENE": "GSIM_00010"}, search_db="same"))
    assert "engine_protein_mismatch" not in flags_of(out[("NOGENE", 4000, 4030, "+")]["flags"])


def test_engine_protein_mismatch_on_pseudogene_and_nterm() -> None:
    rows = [hit_row("PS", 3100, 3130, "+"), hit_row("NT", 70, 100, "+")]
    p = peps(rows, {"PS": "GSIM_00010", "NT": "GSIM_00005"}, search_db="same")
    out = annotate(rows, peptides=p)
    assert "engine_protein_mismatch" in flags_of(out[("PS", 3100, 3130, "+")]["flags"])
    assert "engine_protein_mismatch" not in flags_of(out[("NT", 70, 100, "+")]["flags"])


def test_same_db_flag_kept() -> None:
    rows = [hit_row("SAMEDB", 130, 160, "+", flags="same_db")]
    out = annotate(rows, relation="same", peptides=peps(rows, {"SAMEDB": G1}, search_db="same"))
    assert flags_of(out[("SAMEDB", 130, 160, "+")]["flags"]) == {"same_db"}


def test_low_complexity_and_multi_flags_combine() -> None:
    rows = [
        hit_row("KKKK", 130, 160, "+", flags="low_complexity"),
        hit_row("KKKK", 4000, 4030, "+", flags="low_complexity"),
    ]
    for r in annotate(rows).values():
        assert {"low_complexity", "multi"} <= flags_of(r["flags"])


def test_empty_input() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    out = annotate_peptide_hits(acc([]), MODEL, GENOME, peps([]), "related")
    assert len(out) == 0
    for col in ("locus_class", "gene", "uniqueness", "flags", "novel_candidate_id"):
        assert col in out.columns


def test_unknown_genome_relation_rejected() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    rows = [hit_row("PEPA", 130, 160, "+")]
    with pytest.raises(ValueError):
        annotate_peptide_hits(acc(rows), MODEL, GENOME, peps(rows), "cousin")  # type: ignore[arg-type]


# --- determinism ----------------------------------------------------------------------------

MANY = [
    hit_row("MIX", 130, 160, "+"),
    hit_row("MIX", 4000, 4030, "+"),
    hit_row("IDENT", 7030, 7060, "+"),
    hit_row("IDENT", 8240, 8270, "-"),
    hit_row("REP", 130, 160, "+"),
    hit_row("REP", 250, 280, "+"),
    hit_row("TOL", 1270, 1300, "-", route="tolerant", n_diff=2, diff_types="substitution;swap"),
    hit_row("ORI", 9995, 10025, "-"),
    hit_row("RNA", 5130, 5160, "+"),
    hit_row("NT", 70, 100, "+"),
]


SORT_KEY = ["peptide", "contig", "start", "strand"]


def _assert_canonical(df: pd.DataFrame) -> None:
    """Lead decision 1: rows sorted by (peptide, contig, start, strand), flags alphabetical."""
    keys = list(zip(*(df[c].tolist() for c in SORT_KEY), strict=True))
    assert keys == sorted(keys), keys
    for x in df["flags"]:
        fl = [] if x is None or (not isinstance(x, str) and pd.isna(x)) else str(x).split(";")
        fl = [f for f in fl if f]
        assert fl == sorted(fl), x
        assert len(fl) == len(set(fl)), x


def test_deterministic_and_row_order_invariant() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    p = peps(MANY, {"NT": G2}, search_db="same")
    a = annotate_peptide_hits(acc(MANY), MODEL, GENOME, p, "related")
    _assert_canonical(a)
    b = annotate_peptide_hits(acc(MANY), MODEL, GENOME, p, "related")
    pd.testing.assert_frame_equal(a, b)
    rng = np.random.default_rng(14)
    for _ in range(5):
        perm = [MANY[i] for i in rng.permutation(len(MANY))]
        c = annotate_peptide_hits(acc(perm), MODEL, GENOME, p.sample(frac=1.0, random_state=1),
                                  "related")  # fmt: skip
        pd.testing.assert_frame_equal(a.reset_index(drop=True), c.reset_index(drop=True))


def test_flags_joined_alphabetically() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    rows = [
        hit_row("ZZTOL", 130, 160, "+", route="tolerant", n_diff=1, diff_types="substitution",
                flags="low_complexity"),
        hit_row("ZZTOL", 4000, 4030, "+", route="tolerant", n_diff=1, diff_types="substitution",
                flags="low_complexity"),
        hit_row("AAORI", 9995, 10025, "-", route="tolerant", n_diff=1, diff_types="substitution"),
        hit_row("MMNT", 70, 100, "+", flags="same_db"),
    ]  # fmt: skip
    p = peps(rows, {"MMNT": G2}, search_db="same")
    p.loc[p["peptide"] != "MMNT", "search_db"] = "other"
    out = annotate_peptide_hits(acc(rows), MODEL, GENOME, p, "related")
    _assert_canonical(out)
    got = [(r["peptide"], int(r["start"]), r["flags"]) for r in out.to_dict("records")]
    assert got == [
        ("AAORI", 9995, "origin_spanning;strain_difference"),
        ("MMNT", 70, "engine_protein_mismatch;nterm_extension_evidence;same_db"),
        ("ZZTOL", 130, "low_complexity;multi;strain_difference"),
        ("ZZTOL", 4000, "low_complexity;multi;strain_difference"),
    ]


def test_rows_sorted_whatever_the_input_order() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    rows = [
        hit_row("B", 8240, 8270, "-"),
        hit_row("A", 4000, 4030, "-"),
        hit_row("A", 4000, 4030, "+"),
        hit_row("A", 130, 160, "+"),
    ]
    out = annotate_peptide_hits(acc(rows), MODEL, GENOME, peps(rows), "related")
    got = [(r["peptide"], int(r["start"]), r["strand"]) for r in out.to_dict("records")]
    assert got == [("A", 130, "+"), ("A", 4000, "+"), ("A", 4000, "-"), ("B", 8240, "-")]


def test_inputs_not_mutated() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    a, p = acc(MANY), peps(MANY)
    a0, p0 = a.copy(), p.copy()
    annotate_peptide_hits(a, MODEL, GENOME, p, "related")
    pd.testing.assert_frame_equal(a, a0)
    pd.testing.assert_frame_equal(p, p0)


def test_gene_is_set_exactly_for_assigning_classes() -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    out = annotate_peptide_hits(acc(MANY), MODEL, GENOME, peps(MANY), "related")
    assigning = {"cds_inframe", "sorf", "nterm_extension", "pseudogene"}
    for r in out.to_dict("records"):
        assert (not na(r["gene"])) == (r["locus_class"] in assigning), r
        if not na(r["gene"]):
            MODEL.by_locus(r["gene"])  # every gene is a locus tag of the model
        assert r["uniqueness"] in {"unique", "multi_identical", "shared"}


# --- proteins (example 15) ------------------------------------------------------------------


def prot_row(pid: str, gene: str | None, start: int = 100, end: int = 400, rank: int = 1,
             flags: str = "") -> dict[str, Any]:  # fmt: skip
    return {
        "protein_id": pid,
        "protein_len_aa": (end - start) // 3 - 1,
        "route": "miniprot",
        "contig": "chr1",
        "start": start,
        "end": end,
        "strand": "+",
        "rank": rank,
        "identity": 0.99,
        "query_coverage": 1.0,
        "reciprocal_overlap": 1.0,
        "class": "identical",
        "gene": gene,
        "flags": flags,
    }


def annotate_prot(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    from gippa.mapping.loci import annotate_protein_hits

    df = annotate_protein_hits(pd.DataFrame(list(rows)), MODEL)
    assert len(df) == len(rows)
    return sorted(df.to_dict("records"), key=lambda r: (r["protein_id"], int(r["rank"])))


def test_example_15_two_proteins_one_gene() -> None:
    out = annotate_prot([prot_row("P0A7V8", G1), prot_row(G1, G1)])
    fl = {r["protein_id"]: flags_of(r["flags"]) for r in out}
    assert fl["P0A7V8"] == {"two_proteins_one_gene", "cross_strain_protein"}
    assert fl[G1] == {"two_proteins_one_gene"}


def test_protein_flags_negative_cases() -> None:
    out = annotate_prot(
        [
            prot_row(G1, G1, flags="fs"),
            prot_row(G1, G1, rank=2),  # the same protein twice is one protein
            prot_row(G2, G2, start=1000, end=1300),
        ]
    )
    fl = [(r["protein_id"], int(r["rank"]), flags_of(r["flags"])) for r in out]
    assert fl == [(G1, 1, {"fs"}), (G1, 2, set()), (G2, 1, set())]


def test_cross_strain_protein_alone() -> None:
    out = annotate_prot([prot_row("WP_000001", G2, 1000, 1300), prot_row(G1, G1, flags="fs;st")])
    fl = {r["protein_id"]: flags_of(r["flags"]) for r in out}
    assert fl["WP_000001"] == {"cross_strain_protein"}
    assert fl[G1] == {"fs", "st"}


def test_two_proteins_ignores_unassigned_rows() -> None:
    out = annotate_prot([prot_row(G1, G1), prot_row(G2, None, 4000, 4300)])
    for r in out:
        assert "two_proteins_one_gene" not in flags_of(r["flags"])


# --- timing guard (the benchmark of the acceptance criteria) -------------------------------


def big_case() -> tuple[pd.DataFrame, GeneModel, dict[str, str], pd.DataFrame]:
    """500,000 hits on a 5-Mb, ~4,500-gene model (shared with benchmarks/test_loci.py)."""
    rng = np.random.default_rng(5)
    n_len = 5_000_000
    reps = {"big": Replicon("big", n_len, "circular", "chromosome")}
    genes, pos = [], 0
    for k in range(4_500):
        pos += int(rng.integers(50, 300))
        ln = 3 * int(rng.integers(50, 400))
        if pos + ln >= n_len:
            break
        genes.append(cds(f"B_{k:05d}", "big", pos, pos + ln, "+" if k % 2 else "-"))
        pos += ln
    codes = rng.integers(0, 4, size=n_len, dtype=np.uint8)
    genome = {"big": np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")}
    gm = GeneModel(reps, tuple(genes), ())
    n = 500_000
    starts = rng.integers(0, n_len - 60, size=n)
    strands = np.where(rng.random(n) < 0.5, "+", "-")
    pep_names = [f"P{i % 300_000:06d}" for i in range(n)]
    frames = np.where(strands == "+", 1 + starts % 3, -(1 + ((n_len - starts - 30) % n_len) % 3))
    df = pd.DataFrame(
        {
            "peptide": pep_names,
            "query": pep_names,
            "search": "target",
            "decoy_type": None,
            "contig": "big",
            "frame_id": 0,
            "aa_start": 0,
            "aa_end": 10,
            "start": starts.astype("int64"),
            "end": (starts + 30).astype("int64"),
            "strand": strands,
            "frame": frames.astype("int64"),
            "route": "exact",
            "n_diff": 0,
            "diff_types": "",
            "blocks": None,
            "source": "database",
            "placement_q": 0.001,
            "flags": "",
        }
    )
    p = pd.DataFrame({"peptide": sorted(set(pep_names))})
    p["protein_ids"], p["search_db"], p["source"] = None, "other", "database"
    p["is_decoy"], p["is_contaminant"] = False, False
    return df, gm, genome, p


def time_big_case(case: tuple[pd.DataFrame, GeneModel, dict[str, str], pd.DataFrame]) -> float:
    from gippa.mapping.loci import annotate_peptide_hits

    df, gm, genome, p = case
    t0 = time.perf_counter()
    out = annotate_peptide_hits(df, gm, genome, p, "related")
    secs = time.perf_counter() - t0
    assert len(out) == len(df)
    return secs


@pytest.mark.slow
def test_500k_hits_loose_guard() -> None:
    """Loose 10x guard; the 20 s budget itself lives in benchmarks/test_loci.py."""
    secs = time_big_case(big_case())
    assert secs <= 200.0, f"500k hits took {secs:.1f} s (10x the 20 s budget)"
