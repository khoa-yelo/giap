"""M7-02 pileup, calls and confidence (docs/specs/M7-02-pileup-and-calls.md; GIPPA_SPEC.md §7.5.1,
§7.5.2 steps 2-4, §6.5, §9.3).

Examples 1-12 and every degenerate row of the spec, on the hand-built genome of
``helpers.pileup_oracle``: circular ``c1`` (1,000 nt), G1 ``+`` [100, 400) with codon 112 = ``N``
(aa_pos 5) and an annotated GTG start, G2 ``-`` [500, 800) with codon 785 = ``K`` (aa_pos 5), G3
``+`` [950, 1100) across the origin. Alignments are built by hand from M7-01's dataclasses.

Interface reading (test author; see the report to the lead):

* ``eligible(aligns, support, cfg) -> list[HitAlignment]``; under ``genome_relation = same`` a
  ``database`` / ``search_db = other`` alignment keeps its ``match``/``il`` ops and loses every
  difference op and its frameshift event.
* ``pileup(aligns, gm, genome, cfg, support=support)``: the spec's ``sites`` row needs each ALT
  peptide's ``source`` and ``pileup`` has no other way to know it, so ``support`` is passed by
  keyword (it may be optional). Frameshift events reach ``sites(pile, cfg)`` through ``pile``
  in a way these tests do not constrain (extra columns or attributes), but ``pile`` holds exactly
  one row per covered codon.
* The §9.3 / rule-11 columns are present and in order; extra internal columns may follow.
* ``codon_start`` is 0-based (the writer of M7-03 converts); ``alt_counts`` empty is ``""`` or NA.
* Ties of rule 11 sort ``alt_type`` and ``alt_aa`` as strings.
* An I-peptide and an L-peptide that fold to one sequence read ``I or L`` at an ``L`` codon: one
  distinct peptide with IL support (degenerate row 5), not REF.
"""

from __future__ import annotations

import copy
import importlib
from collections.abc import Mapping, Sequence
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers.align_oracle import place
from helpers.fs_genomes import background
from helpers.pileup_oracle import (
    G1,
    G2,
    G3,
    SEQ,
    GeneSpec,
    L,
    codon_aa,
    codon_nt,
    gene_hit,
    gene_ops,
    make_hit,
    nt,
)

from gippa.config.model import GippaConfig
from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel

PILEUP_COLUMNS = [
    "contig", "codon_start", "strand", "gene", "aa_pos", "ref_aa", "n_ref", "n_il",
    "alt_counts", "call",
]  # fmt: skip
SITE_COLUMNS = [
    "variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon", "ref_aa",
    "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "source",
]  # fmt: skip

REP = Replicon("c1", L, "circular")


def pl() -> ModuleType:
    """The module under test (imported per test, so each test fails on its own)."""
    return importlib.import_module("gippa.variants.pileup")


def al() -> ModuleType:
    return importlib.import_module("gippa.variants.align")


def gene(spec: GeneSpec, ftype: str = "cds", start_codon: str | None = "ATG") -> Gene:
    return Gene(
        spec.tag, ftype, Interval("c1", spec.start, spec.end, spec.strand),  # type: ignore[arg-type]
        "M" + spec.protein[1:], None, start_codon,  # type: ignore[arg-type]
    )  # fmt: skip


GM = GeneModel({"c1": REP}, (gene(G1, start_codon="GTG"), gene(G2), gene(G3)), ())
GENOME = {"c1": SEQ}


def cfg(**over: Any) -> GippaConfig:
    return GippaConfig.model_validate(over)


CFG = cfg()


def sup(
    *peptides: str,
    source: str = "denovo",
    search_db: str = "same",
    q: float = 0.001,
    uniqueness: str = "unique",
) -> dict[str, Any]:
    return {p: pl().Support(p, source, search_db, q, uniqueness) for p in peptides}


def support_of(aligns: Sequence[Any], **kw: Any) -> dict[str, Any]:
    return sup(*sorted({a.peptide for a in aligns}), **kw)


def run(
    aligns: Sequence[Any],
    support: Mapping[str, Any] | None = None,
    config: GippaConfig = CFG,
    gm: GeneModel = GM,
    genome: Mapping[str, str] = GENOME,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = pl()
    support = support_of(aligns) if support is None else support
    el = m.eligible(list(aligns), support, config)
    pile = m.pileup(el, gm, genome, config, support=support)
    return pile, m.sites(pile, config)


def row(pile: pd.DataFrame, codon: int, strand: str = "+", contig: str = "c1") -> dict[str, Any]:
    sel = pile[
        (pile["contig"] == contig) & (pile["codon_start"] == codon) & (pile["strand"] == strand)
    ]
    assert len(sel) == 1, f"expected one pileup row at {contig}:{codon}{strand}, got {len(sel)}"
    return sel.iloc[0].to_dict()


def has_row(pile: pd.DataFrame, codon: int, strand: str = "+", contig: str = "c1") -> bool:
    return bool(
        (
            (pile["contig"] == contig) & (pile["codon_start"] == codon) & (pile["strand"] == strand)
        ).any()
    )


def site_rows(sites: pd.DataFrame, codon: int, strand: str = "+") -> list[dict[str, Any]]:
    sel = sites[(sites["codon_start"] == codon) & (sites["strand"] == strand)]
    return [r for _, r in sel.iterrows()]  # type: ignore[misc]


def empty_counts(v: object) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v)) or v is pd.NA or v == ""


def ordered_subset(cols: Sequence[str], want: Sequence[str]) -> list[str]:
    return [c for c in cols if c in want]


G1_112 = G1.codon(4)
assert G1_112 == 112
assert G2.codon(4) == 785


# --------------------------------------------------------------------------- columns


def test_pileup_columns_follow_section_9_3() -> None:
    pile, sites = run([gene_hit(al(), G1, 2, 8, {4: "D"})])
    assert set(PILEUP_COLUMNS) <= set(pile.columns)
    assert ordered_subset(list(pile.columns), PILEUP_COLUMNS) == PILEUP_COLUMNS
    assert set(SITE_COLUMNS) <= set(sites.columns)
    assert ordered_subset(list(sites.columns), SITE_COLUMNS) == SITE_COLUMNS


def test_no_accepted_peptides_gives_header_only_frames() -> None:
    pile, sites = run([], {})
    assert len(pile) == 0 and len(sites) == 0
    assert ordered_subset(list(pile.columns), PILEUP_COLUMNS) == PILEUP_COLUMNS
    assert ordered_subset(list(sites.columns), SITE_COLUMNS) == SITE_COLUMNS


def test_one_row_per_covered_codon_and_none_elsewhere() -> None:
    pile, _ = run([gene_hit(al(), G1, 2, 8)])
    assert sorted(pile["codon_start"].tolist()) == [106, 109, 112, 115, 118, 121]
    assert set(pile["strand"]) == {"+"} and set(pile["contig"]) == {"c1"}


# --------------------------------------------------------------------------- examples 1-12


def test_example1_ref_only() -> None:
    a = al()
    pile, sites = run([gene_hit(a, G1, 2, 8), gene_hit(a, G1, 3, 9)])
    r = row(pile, 112)
    assert (r["n_ref"], r["n_il"], r["call"]) == (2, 0, "REF")
    assert empty_counts(r["alt_counts"])
    assert (r["gene"], int(r["aa_pos"]), r["ref_aa"]) == ("G1", 5, "N")
    assert len(sites) == 0


def test_example2_single_alt_peptide_is_candidate() -> None:
    pile, sites = run([gene_hit(al(), G1, 2, 8, {4: "D"})])
    r = row(pile, 112)
    assert (r["n_ref"], r["alt_counts"], r["call"]) == (0, "D:1", "ALT")
    (s,) = site_rows(sites, 112)
    assert s["variant_id"] == "VAR_1"
    assert (s["contig"], s["gene"], int(s["aa_pos"])) == ("c1", "G1", 5)
    assert (s["ref_codon"], s["ref_aa"], s["alt_aa"]) == (codon_nt(SEQ, 112, "+"), "N", "D")
    assert s["ref_codon"] == "AAC"
    assert (s["alt_type"], s["call"], s["confidence"]) == ("substitution", "ALT", "candidate")
    assert (int(s["n_ref"]), int(s["n_alt"]), s["source"]) == (0, 1, "denovo")


def test_example3_two_distinct_alt_peptides_are_high() -> None:
    a = al()
    _, sites = run([gene_hit(a, G1, 2, 8, {4: "D"}), gene_hit(a, G1, 3, 9, {4: "D"})])
    (s,) = site_rows(sites, 112)
    assert (s["call"], s["confidence"], int(s["n_alt"])) == ("ALT", "high", 2)


def test_example4_charges_and_runs_collapse_to_one_peptide() -> None:
    """C9 is one row per stripped peptide; repeated alignments of it still count once."""
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    pile, sites = run([h, h, gene_hit(a, G1, 2, 8, {4: "D"})])
    assert row(pile, 112)["alt_counts"] == "D:1"
    (s,) = site_rows(sites, 112)
    assert (s["confidence"], int(s["n_alt"])) == ("candidate", 1)


def test_example5_mixed() -> None:
    a = al()
    hits = [gene_hit(a, G1, 2, 8), gene_hit(a, G1, 2, 9, {4: "D"}), gene_hit(a, G1, 3, 9, {4: "D"})]
    pile, sites = run(hits)
    r = row(pile, 112)
    assert (r["n_ref"], r["alt_counts"], r["call"]) == (1, "D:2", "mixed")
    (s,) = site_rows(sites, 112)
    assert (s["call"], s["confidence"], int(s["n_ref"]), int(s["n_alt"])) == ("mixed", "high", 1, 2)


def test_example6_il_ignored_when_ref_exists() -> None:
    a = al()
    assert G1.codon(10) == 130 and codon_aa(SEQ, 130, "+") == "L"
    hits = [gene_hit(a, G1, 8, 14, {10: "I"}), gene_hit(a, G1, 9, 14)]
    pile, sites = run(hits)
    r = row(pile, 130)
    assert (r["n_ref"], r["n_il"], r["call"]) == (1, 1, "REF")
    assert empty_counts(r["alt_counts"]) and len(sites) == 0


def test_example7_il_only_is_no_call() -> None:
    pile, sites = run([gene_hit(al(), G1, 8, 14, {10: "I"})])
    r = row(pile, 130)
    assert (r["n_ref"], r["n_il"], r["call"]) == (0, 1, "no_call")
    assert len(sites) == 0


def test_example8_i_and_l_at_a_v_codon_are_one_alt_j() -> None:
    a = al()
    assert G1.codon(11) == 133 and codon_aa(SEQ, 133, "+") == "V"
    hits = [gene_hit(a, G1, 9, 15, {11: "I"}), gene_hit(a, G1, 11, 15, {11: "L"})]
    pile, sites = run(hits)
    r = row(pile, 133)
    assert (r["alt_counts"], r["call"]) == ("J:2", "ALT")
    (s,) = site_rows(sites, 133)
    assert (s["alt_aa"], s["confidence"], int(s["n_alt"]), s["ref_aa"]) == ("J", "high", 2, "V")


def test_example9_two_alt_residues_two_sites() -> None:
    a = al()
    pile, sites = run([gene_hit(a, G1, 2, 8, {4: "D"}), gene_hit(a, G1, 3, 9, {4: "E"})])
    r = row(pile, 112)
    assert (r["alt_counts"], r["call"]) == ("D:1;E:1", "ALT")
    rows = site_rows(sites, 112)
    assert [s["alt_aa"] for s in rows] == ["D", "E"]
    assert all(s["confidence"] == "candidate" and s["call"] == "ALT" for s in rows)
    assert [s["variant_id"] for s in rows] == ["VAR_1", "VAR_2"]


def test_example10_block_at_first_codon_only() -> None:
    a = al()
    ops = [*gene_ops(G1, 2, 4), ("block", (112,), "N", "GG")]
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    r = row(pile, 112)
    assert (r["alt_counts"], r["call"], r["ref_aa"]) == ("N>GG:1", "ALT", "N")
    (s,) = site_rows(sites, 112)
    assert (s["alt_type"], s["alt_aa"], s["ref_aa"], s["confidence"]) == (
        "block", "GG", "N", "candidate",
    )  # fmt: skip
    assert not has_row(pile, 115)


def test_example11_other_strain_database_peptide_under_same() -> None:
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    pile, sites = run([h], sup(h.peptide, source="database", search_db="other"),
                      cfg(genome_relation="same"))  # fmt: skip
    assert not has_row(pile, 112)
    assert len(sites) == 0
    for c in (106, 109, 115, 118, 121):
        assert row(pile, c)["call"] == "REF"


def test_example12_minus_strand_aa_pos() -> None:
    pile, sites = run([gene_hit(al(), G2, 2, 8, {4: "R"})])
    r = row(pile, 785, "-")
    assert (r["gene"], int(r["aa_pos"]), r["ref_aa"]) == ("G2", 5, "K")
    (s,) = site_rows(sites, 785, "-")
    assert int(s["aa_pos"]) == (800 - 785 - 3) // 3 + 1 == 5
    assert s["ref_codon"] == codon_nt(SEQ, 785, "-") == "AAA"
    assert (s["ref_aa"], s["alt_aa"]) == ("K", "R")


# --------------------------------------------------------------------------- confidence


@pytest.mark.parametrize(("k", "n", "want"), [
    (2, 1, "candidate"), (2, 2, "high"), (2, 3, "high"),
    (3, 1, "candidate"), (3, 2, "candidate"), (3, 3, "high"),
])  # fmt: skip
def test_confidence_threshold(k: int, n: int, want: str) -> None:
    a = al()
    hits = [gene_hit(a, G1, 3 - i, 8, {4: "D"}) for i in range(n)]
    _, sites = run(hits, config=cfg(variants={"high_min_peptides": k}))
    (s,) = site_rows(sites, 112)
    assert (int(s["n_alt"]), s["confidence"]) == (n, want)


def test_mixed_alt_part_confidence_is_per_alt_key() -> None:
    a = al()
    hits = [gene_hit(a, G1, 2, 8), gene_hit(a, G1, 1, 8, {4: "D"}), gene_hit(a, G1, 2, 9, {4: "D"}),
            gene_hit(a, G1, 3, 8, {4: "E"})]  # fmt: skip
    pile, sites = run(hits)
    assert row(pile, 112)["alt_counts"] == "D:2;E:1"
    got = {
        s["alt_aa"]: (s["call"], s["confidence"], int(s["n_alt"])) for s in site_rows(sites, 112)
    }
    assert got == {"D": ("mixed", "high", 2), "E": ("mixed", "candidate", 1)}


def test_alt_counts_sorted_by_count_then_key() -> None:
    a = al()
    hits = [gene_hit(a, G1, 3, 8, {4: "E"})]
    hits += [gene_hit(a, G1, 3 - i, 9, {4: "Q"}) for i in range(2)]
    hits += [gene_hit(a, G1, 3 - i, 10, {4: "D"}) for i in range(2)]
    pile, _ = run(hits)
    assert row(pile, 112)["alt_counts"] == "D:2;Q:2;E:1"


# --------------------------------------------------------------------------- distinct peptides


def test_folded_peptides_at_an_alt_codon_are_one() -> None:
    """Rule 2: two peptides differing only by I/L elsewhere support D once."""
    a = al()
    h1 = gene_hit(a, G1, 2, 11, {4: "D", 10: "I"})
    h2 = gene_hit(a, G1, 2, 11, {4: "D"})
    assert h1.peptide != h2.peptide and h1.peptide.replace("I", "L") == h2.peptide
    pile, sites = run([h1, h2])
    assert row(pile, 112)["alt_counts"] == "D:1"
    (s,) = site_rows(sites, 112)
    assert s["confidence"] == "candidate"


def test_folded_i_and_l_peptides_at_an_l_codon_are_one_il() -> None:
    """Degenerate row: I-peptide and L-peptide at an L codon -> n_il 1, no_call."""
    a = al()
    h1 = gene_hit(a, G1, 8, 14, {10: "I"})
    h2 = gene_hit(a, G1, 8, 14)
    assert h1.peptide.replace("I", "L") == h2.peptide
    pile, _ = run([h1, h2])
    r = row(pile, 130)
    assert (r["n_ref"], r["n_il"], r["call"]) == (0, 1, "no_call")
    assert row(pile, 127)["n_ref"] == 1  # one distinct peptide everywhere else


def test_il_with_alt_is_alt() -> None:
    a = al()
    pile, sites = run([gene_hit(a, G1, 8, 14, {10: "I"}), gene_hit(a, G1, 9, 14, {10: "D"})])
    r = row(pile, 130)
    assert (r["n_il"], r["n_ref"], r["alt_counts"], r["call"]) == (1, 0, "D:1", "ALT")
    (s,) = site_rows(sites, 130)
    assert (int(s["n_ref"]), int(s["n_alt"])) == (0, 1)


def test_substitution_op_between_i_and_l_is_il_never_alt() -> None:
    """Rule 5: a ``substitution`` op with genome I/L and peptide L/I counts as IL."""
    a = al()
    ops = [*gene_ops(G1, 8, 10), ("substitution", (130,), "L", "I"), *gene_ops(G1, 11, 14)]
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    r = row(pile, 130)
    assert (r["n_il"], r["call"]) == (1, "no_call")
    assert empty_counts(r["alt_counts"]) and len(sites) == 0


def test_substitution_from_i_or_l_to_another_residue_is_not_j() -> None:
    pile, sites = run([gene_hit(al(), G1, 8, 14, {10: "V"})])
    assert row(pile, 130)["alt_counts"] == "V:1"
    (s,) = site_rows(sites, 130)
    assert (s["ref_aa"], s["alt_aa"]) == ("L", "V")


def test_multi_identical_counts_at_every_copy() -> None:
    """Decision 3: an identical-protein peptide supports both copies."""
    a = al()
    h1 = gene_hit(a, G1, 2, 8, {4: "D"})
    shifted = [(k, (c + 600,), r, x) for k, (c,), r, x in gene_ops(G1, 2, 8, {4: "D"})]
    copy_genome = SEQ[:700] + SEQ[100:130] + SEQ[730:]
    h2 = make_hit(a, "c1", "+", shifted)
    assert h1.peptide == h2.peptide
    pile, sites = run([h1, h2], sup(h1.peptide, uniqueness="multi_identical"),
                      genome={"c1": copy_genome})  # fmt: skip
    assert row(pile, 112)["alt_counts"] == "D:1"
    assert row(pile, 712)["alt_counts"] == "D:1"
    assert sorted(sites["codon_start"].tolist()) == [112, 712]


# --------------------------------------------------------------------------- eligibility


def test_shared_peptides_support_nothing() -> None:
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    pile, sites = run([h], sup(h.peptide, uniqueness="shared"))
    assert len(pile) == 0 and len(sites) == 0


def test_shared_peptide_beside_a_unique_one_is_not_counted() -> None:
    a = al()
    hs = gene_hit(a, G1, 2, 8, {4: "D"})
    hu = gene_hit(a, G1, 3, 9)
    support = {**sup(hs.peptide, uniqueness="shared"), **sup(hu.peptide)}
    pile, sites = run([hs, hu], support)
    r = row(pile, 112)
    assert (r["n_ref"], r["call"]) == (1, "REF") and empty_counts(r["alt_counts"])
    assert not has_row(pile, 106) and len(sites) == 0


@pytest.mark.parametrize(("alpha", "q", "kept"), [
    (0.01, 0.01, True), (0.01, 0.0100001, False),
    (0.001, 0.005, False), (0.001, 0.001, True),
])  # fmt: skip
def test_placement_q_threshold_is_alpha_variants(alpha: float, q: float, kept: bool) -> None:
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    config = cfg(stats={"alpha": {"variants": alpha}})
    pile, sites = run([h], sup(h.peptide, q=q), config)
    assert len(pile) == (6 if kept else 0)
    assert len(sites) == int(kept)
    assert len(pl().eligible([h], sup(h.peptide, q=q), config)) == int(kept)


def test_eligible_strips_differences_of_other_db_peptides_under_same() -> None:
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    s = sup(h.peptide, source="database", search_db="other")
    (e,) = pl().eligible([h], s, cfg(genome_relation="same"))
    assert e.peptide == h.peptide and e.frameshift is None
    assert [o for o in e.ops] == [o for o in h.ops if o.kind in ("match", "il")]


@pytest.mark.parametrize(("relation", "source", "search_db"), [
    ("related", "database", "other"),
    ("same", "denovo", "other"),
    ("same", "database", "same"),
    ("same", "denovo", "same"),
])  # fmt: skip
def test_every_other_combination_contributes_alt(
    relation: str, source: str, search_db: str
) -> None:
    a = al()
    h = gene_hit(a, G1, 2, 8, {4: "D"})
    pile, sites = run([h], sup(h.peptide, source=source, search_db=search_db),
                      cfg(genome_relation=relation))  # fmt: skip
    assert row(pile, 112)["call"] == "ALT"
    (s,) = site_rows(sites, 112)
    assert s["source"] == source


def test_other_db_alt_dropped_but_ref_kept_codon_may_be_ref() -> None:
    a = al()
    ho = gene_hit(a, G1, 2, 8, {4: "D"})
    hr = gene_hit(a, G1, 3, 9)
    support = {**sup(ho.peptide, source="database", search_db="other"),
               **sup(hr.peptide, source="database")}  # fmt: skip
    pile, sites = run([ho, hr], support, cfg(genome_relation="same"))
    r = row(pile, 112)
    assert (r["n_ref"], r["call"]) == (1, "REF") and empty_counts(r["alt_counts"])
    assert row(pile, 109)["n_ref"] == 2
    assert len(sites) == 0


def test_source_column_is_sorted_distinct_alt_sources() -> None:
    a = al()
    h1 = gene_hit(a, G1, 2, 8, {4: "D"})
    h2 = gene_hit(a, G1, 3, 9, {4: "D"})
    h3 = gene_hit(a, G1, 1, 9, {4: "D"})
    href = gene_hit(a, G1, 4, 9)
    support = {**sup(h1.peptide, h3.peptide), **sup(h2.peptide, href.peptide, source="database",
                                                   search_db="other")}  # fmt: skip
    _, sites = run([h1, h2, h3, href], support)
    (s,) = site_rows(sites, 112)
    assert (s["source"], int(s["n_alt"]), int(s["n_ref"])) == ("database;denovo", 3, 1)


# --------------------------------------------------------------------------- swaps and blocks


def test_swap_at_first_codon_nothing_at_second() -> None:
    a = al()
    ops = [("swap", (103, 106), "KT", "TK"), *gene_ops(G1, 3, 6)]
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    assert row(pile, 103)["alt_counts"] == "KT>TK:1"
    assert not has_row(pile, 106)
    (s,) = site_rows(sites, 103)
    assert (s["alt_type"], s["alt_aa"], s["ref_aa"]) == ("swap", "TK", "K")


def test_swap_second_codon_counts_only_other_peptides() -> None:
    a = al()
    ops = [("swap", (103, 106), "KT", "TK"), *gene_ops(G1, 3, 6)]
    pile, _ = run([make_hit(a, "c1", "+", ops), gene_hit(a, G1, 2, 6)])
    r = row(pile, 106)
    assert (r["n_ref"], r["call"]) == (1, "REF")
    assert row(pile, 103)["call"] == "ALT"


def test_block_one_residue_on_two_codons() -> None:
    """Block peptide only: ALT at the first codon; the second codon has no row."""
    a = al()
    assert (G1.codon(6), G1.codon(7)) == (118, 121)
    ops = [*gene_ops(G1, 3, 6), ("block", (118, 121), "GG", "N"), *gene_ops(G1, 8, 10)]
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    assert row(pile, 118)["alt_counts"] == "GG>N:1"
    assert not has_row(pile, 121)
    (s,) = site_rows(sites, 118)
    assert (s["alt_type"], s["alt_aa"]) == ("block", "N")


def test_near_block_is_written_as_block() -> None:
    a = al()
    assert G1.codon(13) == 139
    ops = [*gene_ops(G1, 11, 13), ("near_block", (139,), "W", "AD"), *gene_ops(G1, 14, 16)]
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    assert row(pile, 139)["alt_counts"] == "W>AD:1"
    (s,) = site_rows(sites, 139)
    assert (s["alt_type"], s["alt_aa"]) == ("block", "AD")


def test_minus_strand_block_at_codons0() -> None:
    a = al()
    assert (G2.codon(5), G2.codon(6)) == (782, 779)
    ops = [*gene_ops(G2, 2, 5), ("block", (782, 779), "GG", "N"), *gene_ops(G2, 7, 9)]
    pile, sites = run([make_hit(a, "c1", "-", ops)])
    assert row(pile, 782, "-")["alt_counts"] == "GG>N:1"
    assert not has_row(pile, 779, "-")
    (s,) = site_rows(sites, 782, "-")
    assert (s["gene"], int(s["aa_pos"]), s["alt_type"]) == ("G2", 6, "block")


def test_two_block_peptides_same_token_are_high() -> None:
    a = al()
    ops1 = [*gene_ops(G1, 2, 4), ("block", (112,), "N", "GG")]
    ops2 = [*gene_ops(G1, 3, 4), ("block", (112,), "N", "GG")]
    _, sites = run([make_hit(a, "c1", "+", ops1), make_hit(a, "c1", "+", ops2)])
    (s,) = site_rows(sites, 112)
    assert (s["confidence"], int(s["n_alt"])) == ("high", 2)


# --------------------------------------------------------------------------- genes, aa_pos, ref


def test_codon_outside_every_gene_has_na_gene_and_aa_pos() -> None:
    a = al()
    ops = [("match", (c,), codon_aa(SEQ, c, "+"), codon_aa(SEQ, c, "+")) for c in (430, 433)]
    ops.append(("substitution", (436,), codon_aa(SEQ, 436, "+"), "W"))
    pile, sites = run([make_hit(a, "c1", "+", ops)])
    r = row(pile, 436)
    assert pd.isna(r["gene"]) and pd.isna(r["aa_pos"])
    assert r["ref_aa"] == codon_aa(SEQ, 436, "+")
    (s,) = site_rows(sites, 436)
    assert pd.isna(s["gene"]) and pd.isna(s["aa_pos"])
    assert s["ref_codon"] == SEQ[436:439]


def test_gene_on_the_other_strand_or_out_of_frame_is_not_assigned() -> None:
    a = al()
    # minus-strand codons inside G1's span, and plus codons inside G1 in another frame
    minus = [("match", (c,), codon_aa(SEQ, c, "-"), codon_aa(SEQ, c, "-")) for c in (202, 205)]
    plus = [("match", (c,), codon_aa(SEQ, c, "+"), codon_aa(SEQ, c, "+")) for c in (203, 206)]
    pile, _ = run([make_hit(a, "c1", "-", minus), make_hit(a, "c1", "+", plus)])
    for c, s in ((202, "-"), (205, "-"), (203, "+"), (206, "+")):
        r = row(pile, c, s)
        assert pd.isna(r["gene"]) and pd.isna(r["aa_pos"])


def test_several_genes_smallest_locus_tag_and_its_aa_pos() -> None:
    a = al()
    b_gene = Gene("B_02", "cds", Interval("c1", 100, 400, "+"), "M", None, "ATG")
    a_gene = Gene("A_01", "sorf", Interval("c1", 106, 160, "+"), "M", None, "ATG")
    pseudo = Gene("C_03", "pseudogene", Interval("c1", 100, 400, "+"), "M", None, None)
    gm = GeneModel({"c1": REP}, (b_gene, pseudo, a_gene, gene(G2)), ())
    pile, sites = run([gene_hit(a, G1, 2, 8, {4: "D"})], gm=gm)
    r = row(pile, 112)
    assert (r["gene"], int(r["aa_pos"])) == ("A_01", (112 - 106) // 3 + 1)
    r = row(pile, 106)
    assert (r["gene"], int(r["aa_pos"])) == ("A_01", 1)
    (s,) = site_rows(sites, 112)
    assert (s["gene"], int(s["aa_pos"])) == ("A_01", 3)


@pytest.mark.parametrize("ftype", ["sorf", "pseudogene"])
def test_sorf_and_pseudogene_carry_the_gene(ftype: str) -> None:
    gm = GeneModel({"c1": REP}, (gene(G1, ftype=ftype), gene(G2)), ())
    pile, _ = run([gene_hit(al(), G1, 2, 8)], gm=gm)
    r = row(pile, 112)
    assert (r["gene"], int(r["aa_pos"])) == ("G1", 5)


def test_gtg_start_ref_is_v_and_m_peptides_give_alt_m() -> None:
    a = al()
    assert SEQ[100:103] == "GTG"
    hits = [gene_hit(a, G1, 0, 6, {0: "M"}), gene_hit(a, G1, 0, 7, {0: "M"})]
    pile, sites = run(hits)
    r = row(pile, 100)
    assert (r["ref_aa"], r["alt_counts"], int(r["aa_pos"])) == ("V", "M:2", 1)
    (s,) = site_rows(sites, 100)
    assert (s["ref_codon"], s["ref_aa"], s["alt_aa"], s["confidence"]) == ("GTG", "V", "M", "high")


# --------------------------------------------------------------------------- the origin


def test_codon_crossing_the_origin_plus() -> None:
    a = al()
    assert G3.codon(16) == 998
    pile, sites = run([gene_hit(a, G3, 14, 19, {16: "W"})])
    assert sorted(pile["codon_start"].tolist()) == [1, 4, 992, 995, 998]
    r = row(pile, 998)
    assert (r["gene"], int(r["aa_pos"]), r["call"]) == ("G3", 17, "ALT")
    assert r["ref_aa"] == codon_aa(SEQ, 998, "+") == "G"
    (s,) = site_rows(sites, 998)
    assert s["ref_codon"] == SEQ[998] + SEQ[999] + SEQ[0]
    r1 = row(pile, 1)
    assert (r1["gene"], int(r1["aa_pos"])) == ("G3", 18)


def test_same_peptide_unwrapped_and_plain_is_one_key() -> None:
    a = al()
    ops = gene_ops(G3, 17, 19)
    assert [o[1][0] for o in ops] == [1001, 1004]
    unwrapped = make_hit(a, "c1", "+", ops)
    plain = make_hit(a, "c1", "+", [(k, (c - L,), r, x) for k, (c,), r, x in ops])
    assert unwrapped.peptide == plain.peptide
    pile, _ = run([unwrapped, plain])
    assert sorted(pile["codon_start"].tolist()) == [1, 4]
    assert row(pile, 1)["n_ref"] == 1 and row(pile, 4)["n_ref"] == 1
    assert (row(pile, 1)["gene"], int(row(pile, 1)["aa_pos"])) == ("G3", 18)


def test_distinct_peptides_unwrapped_and_plain_add_up() -> None:
    a = al()
    unwrapped = gene_hit(a, G3, 16, 19)  # codons 998, 1001, 1004
    plain_ops = [(k, (c - L,), r, x) for k, (c,), r, x in gene_ops(G3, 17, 20)]
    pile, _ = run([unwrapped, make_hit(a, "c1", "+", plain_ops)])
    assert row(pile, 1)["n_ref"] == 2 and row(pile, 4)["n_ref"] == 2
    assert row(pile, 7)["n_ref"] == 1 and row(pile, 998)["n_ref"] == 1


def test_codon_crossing_the_origin_minus() -> None:
    a = al()
    g4 = GeneSpec("G4", 980, 1010, "-", "KDEGSRAWTQ")
    seq = place(SEQ, 980, nt(g4.protein), "-")
    gm = GeneModel({"c1": REP}, (gene(g4),), ())
    assert [g4.codon(i) for i in range(4)] == [1007, 1004, 1001, 998]
    ops = gene_ops(g4, 1, 5, {3: "W"}, seq=seq)
    pile, sites = run([make_hit(a, "c1", "-", ops)], gm=gm, genome={"c1": seq})
    assert sorted(pile["codon_start"].tolist()) == [1, 4, 995, 998]
    r = row(pile, 998, "-")
    assert (r["gene"], int(r["aa_pos"]), r["ref_aa"]) == ("G4", 4, "G")
    r1 = row(pile, 1, "-")
    assert (r1["gene"], int(r1["aa_pos"]), r1["ref_aa"]) == ("G4", 3, "E")
    (s,) = site_rows(sites, 998, "-")
    assert s["ref_codon"] == codon_nt(seq, 998, "-")
    assert (s["ref_aa"], s["alt_aa"]) == ("G", "W")


def test_linear_contig_no_wrap_codons_kept() -> None:
    a = al()
    lin = Replicon("c1", L, "linear")
    gm = GeneModel({"c1": lin}, (gene(G1), gene(G2)), ())
    pile, _ = run([gene_hit(a, G1, 2, 8)], gm=gm)
    assert sorted(pile["codon_start"].tolist()) == [106, 109, 112, 115, 118, 121]


# --------------------------------------------------------------------------- frameshifts


def fs_plus(a: ModuleType, first: int = 20, last_suffix: int = 3, shift: int = 1) -> Any:
    """G1 residues ``first``..22 (codons ..166), then ``last_suffix`` codons from 166+3+shift."""
    pre = gene_ops(G1, first, 23)
    s0 = 169 + shift
    suf = [("match", (c,), codon_aa(SEQ, c, "+"), codon_aa(SEQ, c, "+"))
           for c in range(s0, s0 + 3 * last_suffix, 3)]  # fmt: skip
    ev = a.FrameshiftEvent(len(pre), 166, s0, shift)
    return make_hit(a, "c1", "+", pre + suf, frameshift=ev)


def fs_rows(sites: pd.DataFrame) -> pd.DataFrame:
    return sites[sites["alt_type"] == "frameshift"]


def test_frameshift_event_is_a_site_not_a_pileup_row() -> None:
    a = al()
    h = fs_plus(a)
    pile, sites = run([h])
    assert sorted(pile["codon_start"].tolist()) == [160, 163, 166, 170, 173, 176]
    assert set(pile["call"]) == {"REF"}
    fs = fs_rows(sites)
    assert len(fs) == 1 and len(sites) == 1
    s = fs.iloc[0]
    assert (s["contig"], int(s["codon_start"]), s["strand"]) == ("c1", 166, "+")
    assert (s["gene"], int(s["aa_pos"])) == ("G1", 23)
    assert (s["call"], s["confidence"], int(s["n_alt"]), int(s["n_ref"])) == (
        "ALT", "candidate", 1, 0,
    )  # fmt: skip


def test_frameshift_n_ref_counts_peptides_spanning_both_codons() -> None:
    a = al()
    spanning = gene_hit(a, G1, 19, 26)  # covers 166 and 169
    ending = gene_hit(a, G1, 18, 23)  # covers 166 only
    starting = gene_hit(a, G1, 23, 27)  # covers 169 only
    _, sites = run([fs_plus(a), fs_plus(a, first=21), spanning, ending, starting])
    (s,) = [r for _, r in fs_rows(sites).iterrows()]
    assert (int(s["n_alt"]), int(s["n_ref"]), s["call"], s["confidence"]) == (
        2, 1, "mixed", "high",
    )  # fmt: skip


def test_frameshift_events_keyed_by_shift() -> None:
    a = al()
    _, sites = run([fs_plus(a), fs_plus(a, first=21, shift=-1)])
    fs = fs_rows(sites)
    assert len(fs) == 2
    assert set(fs["codon_start"].astype(int)) == {166}
    assert set(fs["n_alt"].astype(int)) == {1}


def test_frameshift_minus_strand() -> None:
    a = al()
    assert (G2.codon(20), G2.codon(22)) == (737, 731)
    pre = gene_ops(G2, 20, 23)
    suf = [("match", (c,), codon_aa(SEQ, c, "-"), codon_aa(SEQ, c, "-")) for c in (727, 724)]
    ev = a.FrameshiftEvent(3, 731, 727, 1)
    h = make_hit(a, "c1", "-", pre + suf, frameshift=ev)
    spanning = gene_hit(a, G2, 19, 26)  # covers 731 and 728
    _, sites = run([h, spanning])
    (s,) = [r for _, r in fs_rows(sites).iterrows()]
    assert (int(s["codon_start"]), s["strand"], s["gene"], int(s["aa_pos"])) == (731, "-", "G2", 23)
    assert (int(s["n_alt"]), int(s["n_ref"]), s["call"]) == (1, 1, "mixed")


def test_frameshift_from_other_db_under_same_is_dropped() -> None:
    a = al()
    h = fs_plus(a)
    _, sites = run([h], sup(h.peptide, source="database", search_db="other"),
                   cfg(genome_relation="same"))  # fmt: skip
    assert len(sites) == 0


# --------------------------------------------------------------------------- order and IDs


def two_contig_fixture() -> tuple[GeneModel, dict[str, str]]:
    zz, aa = background(300, seed=24), background(300, seed=22)
    for s in (zz, aa):
        assert "*" not in "".join(codon_aa(s, c, st) for c in range(24, 67, 3) for st in "+-")
    reps = {"zz": Replicon("zz", 300, "circular"), "aa": Replicon("aa", 300, "linear")}
    return GeneModel(reps, (), ()), {"zz": zz, "aa": aa}


def sub(a: ModuleType, genome: Mapping[str, str], contig: str, strand: str, codon: int,
        alt: str, pad: int = 1) -> Any:  # fmt: skip
    seq = genome[contig]
    step = 3 if strand == "+" else -3
    ops = []
    for k in range(-pad, 2):
        c = codon + k * step
        ref = codon_aa(seq, c, strand)
        ops.append(("substitution" if k == 0 else "match", (c,), ref, alt if k == 0 else ref))
    return make_hit(a, contig, strand, ops, length=300)


def test_sites_follow_section_6_5_order_and_ids() -> None:
    a = al()
    gm, genome = two_contig_fixture()
    seq = genome["zz"]
    ref30 = codon_aa(seq, 30, "+") + codon_aa(seq, 33, "+")
    swap = make_hit(a, "zz", "+", [("swap", (30, 33), ref30, ref30[::-1] if ref30[0] != ref30[1]
                                    else "WY")], length=300)  # fmt: skip
    block = make_hit(a, "zz", "+", [("block", (30,), codon_aa(seq, 30, "+"), "GG")], length=300)
    hits = [
        sub(a, genome, "aa", "+", 30, "W"),
        sub(a, genome, "zz", "+", 60, "W"),
        sub(a, genome, "zz", "-", 30, "W"),
        sub(a, genome, "zz", "+", 30, "Y", pad=2),
        sub(a, genome, "zz", "+", 30, "H"),
        swap,
        block,
    ]
    pile, sites = run(hits, gm=gm, genome=genome)
    got = list(zip(sites["contig"], sites["codon_start"].astype(int), sites["strand"],
                   sites["alt_type"], sites["alt_aa"], strict=True))  # fmt: skip
    assert got == [
        ("zz", 30, "+", "block", "GG"),
        ("zz", 30, "+", "substitution", "H"),
        ("zz", 30, "+", "substitution", "Y"),
        ("zz", 30, "+", "swap", swap.ops[0].alt),
        ("zz", 30, "-", "substitution", "W"),
        ("zz", 60, "+", "substitution", "W"),
        ("aa", 30, "+", "substitution", "W"),
    ]
    assert sites["variant_id"].tolist() == [f"VAR_{i}" for i in range(1, 8)]
    assert list(pile["contig"].drop_duplicates()) == ["zz", "aa"]
    keys = list(zip(pile["contig"].map({"zz": 0, "aa": 1}), pile["codon_start"],
                    pile["strand"].map({"+": 0, "-": 1}), strict=True))  # fmt: skip
    assert keys == sorted(keys)


def test_variant_id_is_a_pure_function_of_sorted_rows() -> None:
    a = al()
    base = [gene_hit(a, G1, 2, 8, {4: "D"}), gene_hit(a, G2, 2, 8, {4: "R"})]
    _, s1 = run(base)
    assert s1["variant_id"].tolist() == ["VAR_1", "VAR_2"]
    assert s1["codon_start"].astype(int).tolist() == [112, 785]
    _, s2 = run([*base, gene_hit(a, G1, 0, 3, {1: "R"})])  # earlier codon 103
    assert s2["codon_start"].astype(int).tolist() == [103, 112, 785]
    assert s2["variant_id"].tolist() == ["VAR_1", "VAR_2", "VAR_3"]
    _, s3 = run(list(reversed(base)))
    pd.testing.assert_frame_equal(s1.reset_index(drop=True), s3.reset_index(drop=True))


# --------------------------------------------------------------------------- purity


def test_deterministic_and_inputs_not_mutated() -> None:
    a = al()
    hits = [gene_hit(a, G1, 2, 8, {4: "D"}), gene_hit(a, G1, 3, 9), fs_plus(a),
            gene_hit(a, G2, 2, 8, {4: "R"})]  # fmt: skip
    support = support_of(hits)
    config = cfg()
    hits_c, support_c, genome_c = copy.deepcopy(hits), copy.deepcopy(support), dict(GENOME)
    dump = config.model_dump()
    m = pl()
    el = m.eligible(hits, support, config)
    pile = m.pileup(el, GM, GENOME, config, support=support)
    pile_c = pile.copy(deep=True)
    sites = m.sites(pile, config)
    assert hits == hits_c and support == support_c and genome_c == GENOME
    assert config.model_dump() == dump
    pd.testing.assert_frame_equal(pile, pile_c)
    pile2, sites2 = run(hits, support)
    assert pile.to_csv(index=False) == pile2.to_csv(index=False)
    assert sites.to_csv(index=False) == sites2.to_csv(index=False)


def test_support_dataclass_is_frozen() -> None:
    s = pl().Support("PEPK", "denovo", "same", 0.001, "unique")
    with pytest.raises((AttributeError, TypeError)):
        s.peptide = "X"  # type: ignore[misc]
