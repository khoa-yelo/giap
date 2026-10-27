"""M7-05 variant flags (docs/specs/M7-05-variant-flags.md rules 1-7, examples 1-13, every
degenerate row; GIPPA_SPEC.md §7.5.1, §7.5.3, §9.4, Appendix B.5).

Interface pinned here (test-author reading of the spec's "Public interface"):

* ``gippa.variants.flags`` (pure):

  - ``mass_equivalent(ref, alt, tol_da) -> bool``: ``ref``/``alt`` are residue strings (one
    residue for a substitution, one or two for a block or near-block; ``J`` weighs as I/L).
    Directional: true when ``|m(alt) - (m(ref) + Δ)| <= tol_da`` for Δ in {0} or one of the B.5
    modification deltas (negative deltas such as Amidated included).
  - ``nt_class(ref_codon, alt_aa) -> "snv" | "multi_nt"``.
  - ``site_flags(site, gm, cfg, alt_sources) -> list[str]``: ``site`` is one row of
    ``gippa.variants.pileup.sites`` (the ``SITE_COLUMNS`` keys, 0-based ``codon_start``);
    ``alt_sources`` is the ``(source, search_db)`` of each peptide supporting the site's ALT. The
    list holds the M7-05 flags only (never ``vcf_unrepresentable``, which M7-03 adds), each once,
    in the §9.4 order ``mass_equivalent, snv, multi_nt, strain_difference, database_derived,
    start_codon, stop_codon``.

* REF of a block or near-block is the table-11 translation of its ``ref_codon`` (``ref_aa`` of a
  ``sites`` row is the first codon only, e.g. ``G`` for ``GG>N``); ALT is ``alt_aa``.

Choices the test author made where the spec is ambiguous (reported to the lead):

1. Example 12 lists ``snv;strain_difference;start_codon``, but rule 1 also flags V -> M:
   m(M) - (m(V) + Dioxidation) = 131.04049 - 99.06841 - 31.989829 = -0.0177 Da, within 0.05.
   Rule 1 wins: example 12 (and the "ALT M at a GTG start / internal GTG" rows) expect
   ``mass_equivalent`` as well.
2. Example 13 lists ``stop_codon;snv;strain_difference``; rule 7 fixes the §9.4 order, so the
   expected list is ``snv, strain_difference, stop_codon``.
3. Degenerate row ``tol = 0``: with the B.5 masses N - GG and Q - AG are 0.00001 Da, not 0, so the
   formula of rule 1 does not flag them at tolerance 0. These tests pin only the unambiguous part:
   the exact equalities (GT = AS, AV = GL, AN = GQ, G + Methyl = A, S + Methyl = T; float noise
   below 1e-6 counts as 0, as ``gippa.mapping.masses`` rounds to 6 decimals) flag at 0 and K/Q does
   not. GG -> N at tolerance 0 is left unpinned.
4. ``nt_class``: ``snv`` when the minimal-change codon (G14, M7-03 ``nt_changes``) is at most one
   base away. It differs from a literal "exactly one" only where the codon already encodes the
   residue and no one-base neighbour does ((ATG, M) and (TGG, W)), which cannot be ALT sites.
   The 21 residues of the exhaustive table are the 20 standard ones (I and L separate) plus ``J``.
5. ``start_codon`` (rule 5, "a C2 CDS or sORF"): a pseudogene's first codon does not count; the
   first codon of a minus-strand gene is its lowest-coordinate codon ``end - 3``, and an
   origin-crossing gene's first codon is compared wrapped.
6. ``database_derived`` is tested under ``related`` only (rule 4); ``search_db = other`` alone
   (a de novo peptide) or ``search_db = same`` alone never flags.
"""

from __future__ import annotations

import importlib
from itertools import product
from types import ModuleType
from typing import Any

import pytest
from helpers.align_oracle import place
from helpers.pileup_oracle import G1, SEQ, gene_ops, make_hit

from gippa.config.model import GippaConfig
from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel

ORDER = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon",
)  # fmt: skip
DENOVO = frozenset({("denovo", "same")})
DB_OTHER = frozenset({("database", "other")})

# Appendix B.5, verbatim (independent of the package TSVs).
B5_RESIDUES: dict[str, float] = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}  # fmt: skip
B5_MODS: tuple[float, ...] = (
    15.994915, 31.989829, 0.984016, 14.015650, 28.031300, 42.046950, 42.010565, 57.021464,
    43.005814, 27.994915, 79.966331, -0.984016, -18.010565, -17.026549, 21.981943,
)  # fmt: skip

# The standard genetic code (table 11 translates identically), TCAG order.
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE: dict[str, str] = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for (i, a), (j, b), (k, c) in product(enumerate(_BASES), repeat=3)
}
RESIDUES21 = "ACDEFGHIKLMNPQRSTVWY" + "J"


def fl() -> ModuleType:
    """The module under test (imported per test, so each test fails on its own)."""
    return importlib.import_module("gippa.variants.flags")


def pl() -> ModuleType:
    return importlib.import_module("gippa.variants.pileup")


def al() -> ModuleType:
    return importlib.import_module("gippa.variants.align")


def cfg(relation: str = "related", tol: float | None = None) -> GippaConfig:
    over: dict[str, Any] = {"genome_relation": relation}
    if tol is not None:
        over["variants"] = {"mass_equivalent_tol_da": tol}
    return GippaConfig.model_validate(over)


# --------------------------------------------------------------------------- gene model

L = 1000
REP = Replicon("c1", L, "circular")


def _gene(tag: str, ftype: str, start: int, end: int, strand: str, sc: str) -> Gene:
    return Gene(tag, ftype, Interval("c1", start, end, strand), "M" * ((end - start) // 3),  # type: ignore[arg-type]
                None, sc)  # fmt: skip


# G1 + [100, 400) GTG start; P1 pseudogene + [420, 480); G2 - [500, 800) GTG start (first codon
# 797); S1 sORF + [600, 660) GTG start; G3 + [950, 1100) across the origin, TTG start at 950.
GM = GeneModel(
    {"c1": REP},
    (
        _gene("G1", "cds", 100, 400, "+", "GTG"),
        _gene("P1", "pseudogene", 420, 480, "+", "GTG"),
        _gene("G2", "cds", 500, 800, "-", "GTG"),
        _gene("S1", "sorf", 600, 660, "+", "GTG"),
        _gene("G3", "cds", 950, 1100, "+", "TTG"),
    ),
    (),
)


def site(
    ref_codon: str,
    alt_aa: str,
    alt_type: str = "substitution",
    *,
    codon: int = 103,
    strand: str = "+",
    ref_aa: str | None = None,
    call: str = "ALT",
    n_ref: int = 0,
    n_alt: int = 2,
) -> dict[str, Any]:
    """One ``sites`` row (0-based ``codon_start``); ``ref_aa`` is the first codon's residue."""
    return {
        "variant_id": "VAR_1", "contig": "c1", "codon_start": codon, "strand": strand,
        "gene": "G1", "aa_pos": 2, "ref_codon": ref_codon,
        "ref_aa": CODE[ref_codon[:3]] if ref_aa is None else ref_aa, "alt_aa": alt_aa,
        "alt_type": alt_type, "call": call, "confidence": "high" if n_alt >= 2 else "candidate",
        "n_ref": n_ref, "n_alt": n_alt, "source": "denovo",
    }  # fmt: skip


def flags_of(s: dict[str, Any], config: GippaConfig | None = None,
             sources: frozenset[tuple[str, str]] = DENOVO) -> list[str]:  # fmt: skip
    out = fl().site_flags(s, GM, cfg() if config is None else config, sources)
    assert isinstance(out, list)
    return out


# --------------------------------------------------------------------------- oracles


def mass_oracle(ref: str, alt: str, tol: float) -> bool | None:
    """Rule 1 by brute force over B.5; None when the pair sits within 1e-6 Da of ``tol``."""
    mr = sum(B5_RESIDUES[c] for c in ref.replace("J", "L"))
    ma = sum(B5_RESIDUES[c] for c in alt.replace("J", "L"))
    best = min(abs(ma - mr - d) for d in (0.0, *B5_MODS))
    best = 0.0 if best < 1e-9 else best
    if tol < best <= tol + 1e-6:
        return None
    return best <= tol


def nt_oracle(codon: str, aa: str) -> str:
    want = "IL" if aa == "J" else aa
    dist = min(sum(x != y for x, y in zip(codon, c, strict=True))
               for c, a in CODE.items() if a in want)  # fmt: skip
    return "snv" if dist <= 1 else "multi_nt"


def test_the_standard_code_oracle_is_table_11() -> None:
    assert len(CODE) == 64
    assert CODE["TGA"] == "*" and CODE["ATG"] == "M" and CODE["TGG"] == "W"
    assert CODE["GTG"] == "V" and CODE["TTG"] == "L" and CODE["AGA"] == "R"
    assert sorted(set(CODE.values()) - {"*"}) == sorted("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------- examples 1-13

EXAMPLES: list[tuple[str, dict[str, Any], str, frozenset[tuple[str, str]], list[str]]] = [
    ("1 GCT A->S", site("GCT", "S"), "related", DENOVO,
     ["mass_equivalent", "snv", "strain_difference"]),
    ("2 AAT N->D", site("AAT", "D"), "related", DENOVO,
     ["mass_equivalent", "snv", "strain_difference"]),
    ("3 CAG Q->E", site("CAG", "E"), "related", DENOVO,
     ["mass_equivalent", "snv", "strain_difference"]),
    ("4 AAA K->Q", site("AAA", "Q"), "related", DENOVO,
     ["mass_equivalent", "snv", "strain_difference"]),
    ("5 ATG M->F", site("ATG", "F"), "related", DENOVO,
     ["mass_equivalent", "multi_nt", "strain_difference"]),
    ("6 TTT F->M", site("TTT", "M"), "related", DENOVO, ["multi_nt", "strain_difference"]),
    ("7 TGG W->AD near-block", site("TGG", "AD", "block"), "related", DENOVO,
     ["mass_equivalent", "strain_difference"]),
    ("8 GGTGGC GG->N block", site("GGTGGC", "N", "block"), "related", DENOVO,
     ["mass_equivalent", "strain_difference"]),
    ("9 GCT A->W", site("GCT", "W"), "related", DENOVO, ["multi_nt", "strain_difference"]),
    ("10 example 1 under same", site("GCT", "S"), "same", DENOVO, ["mass_equivalent", "snv"]),
    ("11 database/other peptide", site("GCT", "S"), "related", DB_OTHER,
     ["mass_equivalent", "snv", "strain_difference", "database_derived"]),
    ("12 GTG start V->M", site("GTG", "M", codon=100), "related", DENOVO,
     ["mass_equivalent", "snv", "strain_difference", "start_codon"]),
    ("13 TGA *->W", site("TGA", "W", codon=160), "related", DENOVO,
     ["snv", "strain_difference", "stop_codon"]),
]  # fmt: skip


@pytest.mark.parametrize(
    ("s", "relation", "sources", "want"),
    [e[1:] for e in EXAMPLES],
    ids=[e[0] for e in EXAMPLES],
)
def test_examples(
    s: dict[str, Any], relation: str, sources: frozenset[tuple[str, str]], want: list[str]
) -> None:
    assert flags_of(s, cfg(relation), sources) == want


def test_example_12_from_n_terminal_peptides_on_a_hand_built_alignment() -> None:
    """Example 12 through M7-02: peptides over G1's first residues read M at the GTG start."""
    gm = GeneModel({"c1": REP}, (_gene("G1", "cds", 100, 400, "+", "GTG"),), ())
    aligns = [make_hit(al(), "c1", "+", gene_ops(G1, 0, b, {0: "M"})) for b in (8, 9)]
    rows = _sites(aligns, gm, {"c1": SEQ})
    (row,) = [r for r in rows if r["codon_start"] == 100]
    assert (row["ref_codon"], row["ref_aa"], row["alt_aa"]) == ("GTG", "V", "M")
    got = fl().site_flags(row, gm, cfg(), DENOVO)
    assert got == ["mass_equivalent", "snv", "strain_difference", "start_codon"]


def test_example_13_stop_codon_on_a_hand_built_alignment() -> None:
    """Rule 6 / decision 4: a hand-built alignment consumes a genome TGA (codon 160 of G1) as W."""
    seq = place(SEQ, 160, "TGA", "+")
    gm = GeneModel({"c1": REP}, (_gene("G1", "cds", 100, 400, "+", "GTG"),), ())
    ops = gene_ops(G1, 15, 25, {20: "W"}, seq=seq)
    assert ops[5] == ("substitution", (160,), "*", "W")
    rows = _sites([make_hit(al(), "c1", "+", ops)], gm, {"c1": seq})
    (row,) = [r for r in rows if r["codon_start"] == 160]
    assert (row["ref_codon"], row["ref_aa"]) == ("TGA", "*")
    got = fl().site_flags(row, gm, cfg(), DENOVO)
    assert got == ["snv", "strain_difference", "stop_codon"]


def _sites(aligns: list[Any], gm: GeneModel, genome: dict[str, str]) -> list[dict[str, Any]]:
    m, c = pl(), cfg()
    support = {a.peptide: m.Support(a.peptide, "denovo", "same", 0.001, "unique") for a in aligns}
    pile = m.pileup(m.eligible(aligns, support, c), gm, genome, c, support=support)
    return list(m.sites(pile, c).to_dict("records"))


# --------------------------------------------------------------------------- degenerate rows


def test_tol_zero_flags_exact_equalities_only() -> None:
    """Degenerate row ``tol = 0`` (see choice 3 in the module docstring)."""
    me = fl().mass_equivalent
    for ref, alt in [("GT", "AS"), ("AS", "GT"), ("AV", "GL"), ("AN", "GQ"), ("G", "A"),
                     ("S", "T")]:  # fmt: skip
        assert me(ref, alt, 0.0), (ref, alt)
    assert not me("K", "Q", 0.0)
    assert not me("M", "F", 0.0)
    zero = cfg(tol=0.0)
    assert "mass_equivalent" not in flags_of(site("AAA", "Q"), zero)
    assert "mass_equivalent" in flags_of(site("GGTACT", "AS", "block"), zero)
    assert "mass_equivalent" in flags_of(site("GCTGTT", "GL", "block"), zero)
    assert "mass_equivalent" in flags_of(site("GGT", "A"), zero)


def test_configured_tolerance_is_recomputed_in_process() -> None:
    """Degenerate row ``tol != 0.05`` / lead decision 2: the configured tolerance is honoured,
    not the default-tolerance list."""
    from gippa.mapping.masses import mass_equivalents

    tight = cfg(tol=0.02)
    assert "mass_equivalent" not in flags_of(site("AAA", "Q"), tight)  # Δ 0.0364
    assert "mass_equivalent" not in flags_of(site("ATG", "F"), tight)  # Δ 0.0329
    assert "mass_equivalent" in flags_of(site("TGG", "AD", "block"), tight)  # Δ 0.0153
    # G + Phospho vs H: Δ 0.0711, outside the default list, inside a 0.1 Da tolerance.
    default = {(e.left, e.right) for e in mass_equivalents(0.05)}
    assert ("G+Phospho", "H") not in default
    assert "mass_equivalent" not in flags_of(site("GGT", "H"))
    assert "mass_equivalent" in flags_of(site("GGT", "H"), cfg(tol=0.1))


def test_j_alt_at_a_v_codon() -> None:
    """Degenerate row: GTT V -> J is ``snv`` (ATT or CTT) and V + Methyl = L."""
    assert fl().nt_class("GTT", "J") == "snv"
    assert fl().mass_equivalent("V", "J", 0.05)
    assert flags_of(site("GTT", "J")) == ["mass_equivalent", "snv", "strain_difference"]


def test_alt_m_at_an_annotated_atg_start_is_impossible() -> None:
    """Degenerate row: at an ATG start a peptide reading M is REF support, so no site exists."""
    seq = place(SEQ, 100, "ATG", "+")
    gm = GeneModel({"c1": REP}, (_gene("G1", "cds", 100, 400, "+", "ATG"),), ())
    ops = gene_ops(G1, 0, 8, seq=seq)
    assert ops[0] == ("match", (100,), "M", "M")
    rows = _sites([make_hit(al(), "c1", "+", ops)], gm, {"c1": seq})
    assert rows == []
    assert fl().site_flags(site("ATG", "F", codon=100), gm, cfg(), DENOVO) == [
        "mass_equivalent", "multi_nt", "strain_difference",
    ]  # fmt: skip


def test_alt_m_at_a_gtg_start() -> None:
    """Degenerate row: ``start_codon`` and ``snv`` (GTG -> ATG); rule 1 adds ``mass_equivalent``."""
    assert flags_of(site("GTG", "M", codon=100)) == [
        "mass_equivalent", "snv", "strain_difference", "start_codon",
    ]  # fmt: skip


def test_alt_m_at_an_internal_gtg() -> None:
    assert flags_of(site("GTG", "M", codon=130)) == ["mass_equivalent", "snv", "strain_difference"]


@pytest.mark.parametrize(
    ("codon", "strand", "flagged"),
    [
        (100, "+", True),  # G1 first codon
        (797, "-", True),  # G2 (minus) first codon in transcription direction
        (500, "-", False),  # G2's last codon
        (600, "+", True),  # S1 sORF
        (420, "+", False),  # P1 pseudogene: not a CDS or sORF
        (950, "+", True),  # G3 across the origin
        (100, "-", False),  # no minus-strand gene starts there
        (103, "+", False),  # G1 second codon
    ],
)
def test_start_codon_is_the_first_codon_of_a_cds_or_sorf(
    codon: int, strand: str, flagged: bool
) -> None:
    got = flags_of(site("GTG", "M", codon=codon, strand=strand))
    assert ("start_codon" in got) == flagged


def test_start_codon_needs_alt_m() -> None:
    assert "start_codon" not in flags_of(site("GTG", "J", codon=100))
    assert "start_codon" not in flags_of(site("GTG", "A", codon=797, strand="-"))


@pytest.mark.parametrize(
    ("relation", "sources", "want"),
    [
        ("related", DENOVO, ["strain_difference"]),
        ("related", DB_OTHER, ["strain_difference", "database_derived"]),
        ("same", DENOVO, []),
    ],
)
def test_frameshift_carries_no_mass_nt_or_start_flags(
    relation: str, sources: frozenset[tuple[str, str]], want: list[str]
) -> None:
    """Degenerate row: a frameshift event (here at G1's GTG start) gets rules 3-4 only."""
    for alt in ("fs+1", "fs-1"):
        s = site("GTG", alt, "frameshift", codon=100)
        assert flags_of(s, cfg(relation), sources) == want


def test_mixed_sources_at_one_alt() -> None:
    """Degenerate row: ``database_derived`` if any database/other peptide supports the ALT."""
    s = site("GCT", "S", call="mixed", n_ref=1)
    both = frozenset({("denovo", "same"), ("database", "other")})
    assert flags_of(s, cfg(), both) == [
        "mass_equivalent", "snv", "strain_difference", "database_derived",
    ]  # fmt: skip


@pytest.mark.parametrize(
    "sources",
    [frozenset(), DENOVO, frozenset({("database", "same")}), frozenset({("denovo", "other")})],
)
def test_database_derived_needs_a_database_peptide_searched_against_another_strain(
    sources: frozenset[tuple[str, str]],
) -> None:
    assert "database_derived" not in flags_of(site("GCT", "S"), cfg(), sources)


def test_strain_difference_on_every_alt_under_related_only() -> None:
    for s in (site("GCT", "S"), site("GCT", "W", call="mixed", n_ref=3),
              site("TGG", "AD", "block"), site("AAAGCT", "AK", "swap")):  # fmt: skip
        assert "strain_difference" in flags_of(s, cfg("related"))
        assert "strain_difference" not in flags_of(s, cfg("same"))


def test_swaps_are_mass_equivalent_and_carry_no_nt_class() -> None:
    """ADR-0049 decision 2 (replaces M7-05 lead decision 1): a swap keeps composition, so it is
    always ``mass_equivalent``; lead decision 3 stands: swaps and blocks carry neither ``snv`` nor
    ``multi_nt``."""
    assert flags_of(site("AAAGCT", "AK", "swap")) == ["mass_equivalent", "strain_difference"]
    assert flags_of(site("GGTACT", "TG", "swap")) == ["mass_equivalent", "strain_difference"]
    assert flags_of(site("GGTACT", "TG", "swap"), cfg("same")) == ["mass_equivalent"]
    for s in (site("GGTGGC", "N", "block"), site("TGG", "AD", "block")):
        got = flags_of(s)
        assert "snv" not in got and "multi_nt" not in got


def test_block_ref_is_the_translation_of_ref_codon_not_ref_aa() -> None:
    """``GG>N``: the row's ``ref_aa`` is ``G`` (first codon), the REF composition is ``GG``."""
    s = site("GGTGGC", "N", "block", ref_aa="G")
    assert flags_of(s) == ["mass_equivalent", "strain_difference"]
    # A one-codon REF of N against a GG block ALT, and a two-codon REF of AG against Q.
    assert flags_of(site("AAT", "GG", "block")) == ["mass_equivalent", "strain_difference"]
    assert flags_of(site("GCTGGT", "Q", "block")) == ["mass_equivalent", "strain_difference"]
    assert flags_of(site("CGT", "GV", "block")) == ["mass_equivalent", "strain_difference"]


def test_vcf_unrepresentable_is_never_produced_here() -> None:
    for s in (site("GCT", "S", codon=997), site("GGTGGC", "N", "block", codon=997)):
        assert "vcf_unrepresentable" not in flags_of(s)


def test_site_flags_is_pure() -> None:
    s = site("GCT", "S", call="mixed", n_ref=1, n_alt=1)
    before = dict(s)
    first = flags_of(s, cfg(), DB_OTHER)
    assert flags_of(s, cfg(), DB_OTHER) == first
    assert s == before


# --------------------------------------------------------------------------- §7.5.3 bullets

# Every pair named in §7.5.3 "Expected hits include", in the direction(s) rule 1 flags them.
BULLET_PAIRS: list[tuple[str, str]] = [
    # N vs GG, Q vs AG, K vs Q (0.036 Da), W vs AD/EG, R vs GV: composition, either direction
    ("N", "GG"), ("GG", "N"), ("Q", "AG"), ("AG", "Q"), ("Q", "GA"), ("K", "Q"), ("Q", "K"),
    ("W", "AD"), ("AD", "W"), ("W", "DA"), ("W", "EG"), ("EG", "W"), ("R", "GV"), ("GV", "R"),
    # deamidation
    ("N", "D"), ("Q", "E"),
    # oxidation
    ("A", "S"), ("F", "Y"),
    # methylation
    ("G", "A"), ("S", "T"), ("V", "L"), ("V", "I"), ("V", "J"), ("D", "E"), ("N", "Q"),
    # M(oxidised) vs F
    ("M", "F"),
]  # fmt: skip


@pytest.mark.parametrize(("ref", "alt"), BULLET_PAIRS)
def test_section_7_5_3_pairs_are_mass_equivalent(ref: str, alt: str) -> None:
    assert fl().mass_equivalent(ref, alt, 0.05) is True


@pytest.mark.parametrize(
    ("ref", "alt"),
    [("S", "A"), ("Y", "F"), ("A", "G"), ("T", "S"), ("L", "V"), ("E", "D"), ("F", "M")],
)
def test_modification_pairs_are_directional(ref: str, alt: str) -> None:
    """Rule 1: the observed ALT is the genome residue plus a modification, not the reverse."""
    assert fl().mass_equivalent(ref, alt, 0.05) is False


def test_k_q_drops_out_below_its_delta() -> None:
    assert fl().mass_equivalent("K", "Q", 0.0364)
    assert not fl().mass_equivalent("K", "Q", 0.03)
    assert fl().mass_equivalent("M", "F", 0.034)  # Δ 0.033005
    assert not fl().mass_equivalent("M", "F", 0.03)


@pytest.mark.parametrize("tol", [0.0, 0.02, 0.05, 0.1])
def test_mass_equivalent_matches_a_brute_force_oracle(tol: float) -> None:
    me = fl().mass_equivalent
    checked = 0
    for ref in "ACDEFGHKLMNPQRSTVWY":
        for alt in "ACDEFGHIKLMNPQRSTVWYJ":
            if alt.replace("I", "L").replace("J", "L") == ref:
                continue
            want = mass_oracle(ref, alt, tol)
            if want is None:
                continue
            assert me(ref, alt, tol) is want, (ref, alt, tol)
            checked += 1
    assert checked > 300


@pytest.mark.parametrize("tol", [0.02, 0.05])
def test_mass_equivalent_agrees_with_the_generated_list(tol: float) -> None:
    """Rule 1 / decision 2: the flag is computed from ``mass_equivalents(tol)``: a substitution
    r -> a flags exactly when (r, a) is a composition entry or ``r+<mod>`` vs ``a`` is listed."""
    from gippa.mapping.masses import mass_equivalents

    entries = mass_equivalents(tol)
    pairs: set[tuple[str, str]] = set()
    for e in entries:
        if e.kind == "modification":
            pairs.add((e.left.split("+", 1)[0], e.right))
        else:
            pairs |= {(e.left, e.right), (e.right, e.left)}
    std = "ACDEFGHKLMNPQRSTVWY"
    for ref, alt in product(std, std):
        if ref != alt:
            assert fl().mass_equivalent(ref, alt, tol) is ((ref, alt) in pairs), (ref, alt)
    for e in entries:
        if e.kind in ("block", "near_block"):
            assert fl().mass_equivalent(e.left, e.right, tol), e
            assert fl().mass_equivalent(e.right, e.left, tol), e


# --------------------------------------------------------------------------- nt_class table


def test_nt_class_exhaustive_64_by_21_table() -> None:
    """§7.5.6: every genome codon (stops included) against every residue (and J)."""
    nt_class = fl().nt_class
    cells = 0
    for codon in sorted(CODE):
        for aa in RESIDUES21:
            assert nt_class(codon, aa) == nt_oracle(codon, aa), (codon, aa)
            cells += 1
    assert cells == 64 * 21


@pytest.mark.parametrize(
    ("codon", "aa", "want"),
    [
        ("GCT", "S", "snv"), ("AAT", "D", "snv"), ("CAG", "E", "snv"), ("AAA", "Q", "snv"),
        ("ATG", "F", "multi_nt"), ("TTT", "M", "multi_nt"), ("GCT", "W", "multi_nt"),
        ("GTG", "M", "snv"), ("TGA", "W", "snv"), ("GTT", "J", "snv"), ("TAA", "W", "multi_nt"),
        ("TGA", "C", "snv"), ("ATG", "J", "snv"), ("AAA", "J", "snv"), ("GGG", "J", "multi_nt"),
    ],
)  # fmt: skip
def test_nt_class_spot_values(codon: str, aa: str, want: str) -> None:
    assert fl().nt_class(codon, aa) == want


def test_snv_and_multi_nt_are_exclusive_on_substitutions() -> None:
    for codon in ("GCT", "TGG", "TAG"):
        for aa in "ACDEFGHKMNPQRSTVWYJ":
            if aa == CODE[codon] or (CODE[codon] in "IL" and aa == "J"):
                continue
            got = flags_of(site(codon, aa))
            assert sum(f in ("snv", "multi_nt") for f in got) == 1, (codon, aa, got)


def test_flags_follow_the_section_9_4_order_without_duplicates() -> None:
    cases = [s for _, s, *_ in EXAMPLES]
    for s in cases:
        for relation in ("related", "same"):
            got = flags_of(s, cfg(relation), DB_OTHER | DENOVO)
            assert got == sorted(set(got), key=ORDER.index), got
