"""M7-18 pileup on ambiguous and rare peptide residues
(docs/specs/M7-18-variants-ambiguous-residues.md: rules 3, 4, 6 and 10, examples 1-10, the
edge-case table, lead decisions 1-3).

Hand-built ``HitAlignment``s on the genomes of ``helpers.nonstandard_genome`` go through
``eligible`` -> ``pileup`` -> ``sites``.

Lead decision 1 supersedes the spec's readthrough rule: **every U and O is skipped as ``rare``**
(also U over TGA and O over TAG), ``attrs["nonstandard"]`` has no ``readthrough_sites`` key and
``eligible`` gets no keep-condition for ``ref == "*"``. Examples 7 and 8 are tested in that form.

Not statistical (counting rules): no calibration test is owed.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd
import pytest
from helpers.nonstandard_genome import (
    NONSTANDARD,
    Layout,
    P,
    hit,
    linear,
    minus,
    origin,
    with_x,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config.model import GiapConfig
from giap.variants import align as al
from giap.variants import pileup as pl
from giap.variants.vcf_records import vcf_records

LIN = linear()
KEYS = ["by_residue", "by_rule", "examples", "skipped_codons", "skipped_peptides", "skipped_sites"]
SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def cfg(relation: str = "same") -> GiapConfig:
    return GiapConfig.model_validate({"genome_relation": relation, "sample_name": "S1"})


def support(
    aligns: Sequence[Any], source: str = "denovo", search_db: str = "same"
) -> dict[str, Any]:
    return {
        a.peptide: pl.Support(a.peptide, source, search_db, 0.001, "unique")  # type: ignore[arg-type]
        for a in aligns
    }


def run(
    aligns: Sequence[Any],
    lay: Layout = LIN,
    sup: Mapping[str, Any] | None = None,
    relation: str = "same",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    c = cfg(relation)
    s = support(aligns) if sup is None else sup
    pile = pl.pileup(pl.eligible(list(aligns), s, c), lay.gm(), lay.genome(), c, support=s)
    return pile, pl.sites(pile, c)


def ns(pile: pd.DataFrame) -> dict[str, Any]:
    assert "nonstandard" in pile.attrs, "pileup() must set attrs['nonstandard'] (rule 6)"
    out = pile.attrs["nonstandard"]
    assert isinstance(out, dict)
    return out


def counts(
    sites: int = 0,
    codons: int = 0,
    peptides: int = 0,
    by_residue: Mapping[str, int] | None = None,
    ambiguous: int = 0,
    rare: int = 0,
    examples: Sequence[str] = (),
) -> dict[str, Any]:
    return {
        "by_residue": dict(by_residue or {}),
        "by_rule": {"ambiguous": ambiguous, "rare": rare},
        "examples": list(examples),
        "skipped_codons": codons,
        "skipped_peptides": peptides,
        "skipped_sites": sites,
    }


ZERO = counts()


def has_row(pile: pd.DataFrame, codon: int, strand: str = "+", contig: str = "c1") -> bool:
    sel = (pile["contig"] == contig) & (pile["codon_start"] == codon) & (pile["strand"] == strand)
    return bool(sel.any())


def row(pile: pd.DataFrame, codon: int, strand: str = "+", contig: str = "c1") -> dict[str, Any]:
    sel = pile[
        (pile["contig"] == contig) & (pile["codon_start"] == codon) & (pile["strand"] == strand)
    ]
    assert len(sel) == 1, f"expected one row at {contig}:{codon}{strand}, got {len(sel)}"
    return sel.iloc[0].to_dict()


def codons_of(pile: pd.DataFrame) -> list[int]:
    return sorted(int(c) for c in pile["codon_start"])


ALL12 = [3 * i for i in range(12)]


# --------------------------------------------------------------------------- examples


def test_ex1_x_codon_has_no_row_and_is_counted() -> None:
    pep = "AKXREFGHILWD"
    pile, sites = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    assert "readthrough_sites" not in ns(pile)  # lead decision 1
    assert codons_of(pile) == [c for c in ALL12 if c != 6]
    assert set(pile["call"]) == {"REF"} and set(pile["n_ref"].astype(int)) == {1}
    assert len(sites) == 0


def test_ex2_ref_peptide_alone_supports_the_skipped_codon() -> None:
    pile, sites = run([hit(al, LIN, P), hit(al, LIN, "AKXREFGHILWD")])
    assert ns(pile)["skipped_sites"] == 1
    r = row(pile, 6)
    assert (r["ref_aa"], int(r["n_ref"]), int(r["n_il"]), r["call"]) == ("C", 1, 0, "REF")
    assert pd.isna(r["alt_counts"]) or r["alt_counts"] == ""
    assert len(sites) == 0


def test_ex3_skipped_peptide_is_absent_from_the_alt_counts() -> None:
    p1, p2 = "AKXREFGHILWD", "AKDREFGHILWD"
    pile, sites = run([hit(al, LIN, p1), hit(al, LIN, p2)])
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[p1])
    r = row(pile, 6)
    assert (int(r["n_ref"]), int(r["n_il"]), r["alt_counts"], r["call"]) == (0, 0, "D:1", "ALT")
    (s,) = [x for _, x in sites.iterrows()]
    assert (int(s["codon_start"]), s["alt_aa"], int(s["n_alt"]), int(s["n_ref"])) == (6, "D", 1, 0)
    assert s["confidence"] == "candidate"
    alone, alone_sites = run([hit(al, LIN, p2)])
    assert row(alone, 6) == r
    assert alone_sites.drop(columns="variant_id").equals(sites.drop(columns="variant_id"))


def test_ex4_literal_j_over_an_i_codon_is_skipped() -> None:
    pep = "AKCREFGHJLWD"
    pile, sites = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(1, 1, 1, {"J": 1}, ambiguous=1, examples=[pep])
    assert not has_row(pile, 24), "a literal J is never an ALT (decision 3)"
    assert len(sites) == 0


def test_ex5_j_letter_and_j_token_are_separate() -> None:
    lay = linear("AKCREFGHVLWD")  # codon 24 reads GTT (V)
    p2, p3 = "AKCREFGHJLWD", "AKCREFGHILWD"
    pile, sites = run([hit(al, lay, p2), hit(al, lay, p3)], lay)
    assert ns(pile) == counts(1, 1, 1, {"J": 1}, ambiguous=1, examples=[p2])
    r = row(pile, 24)
    assert (r["ref_aa"], int(r["n_ref"]), r["alt_counts"], r["call"]) == ("V", 0, "J:1", "ALT")
    (s,) = [x for _, x in sites.iterrows()]
    assert (s["alt_aa"], int(s["n_alt"]), s["confidence"]) == ("J", 1, "candidate")
    assert pl.sites(pile, cfg()).attrs["alt_peptides"] == (frozenset({p3}),)


def test_ex6_u_over_a_sense_codon_is_rare() -> None:
    pep = "AKUREFGHILWD"
    pile, sites = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(1, 1, 1, {"U": 1}, rare=1, examples=[pep])
    assert not has_row(pile, 6) and len(sites) == 0


def _stop_case(lay: Layout, letter: str, i: int = 2) -> tuple[Any, str]:
    pep = with_x(lay.protein, i, letter)
    h = hit(al, lay, pep)
    assert h.ops[i].kind == "substitution" and h.ops[i].ref == "*", "fixture: op over a stop"
    return h, pep


@pytest.mark.parametrize(
    ("letter", "stop"),
    [("O", "TGA"), ("U", "TAG"), ("U", "TAA"), ("O", "TAA")],
)
def test_ex7_ex8_u_and_o_over_any_stop_are_rare_plus(letter: str, stop: str) -> None:
    """Lead decision 1, narrowed by M7-16 rule 11: U/O over a stop they do not recode is skipped
    as rare (U over TGA and O over TAG are ALTs: test_variants_readthrough.py)."""
    lay = linear(overrides={2: stop})
    assert lay.sequence()[6:9] == stop
    h, pep = _stop_case(lay, letter)
    pile, sites = run([h], lay)
    assert ns(pile) == counts(1, 1, 1, {letter: 1}, rare=1, examples=[pep])
    assert "readthrough_sites" not in ns(pile)
    assert not has_row(pile, 6) and len(sites) == 0
    assert codons_of(pile) == [c for c in ALL12 if c != 6]


@pytest.mark.parametrize("letter", ["U", "O"])
def test_ex7_u_over_tga_is_a_pileup_site_minus(letter: str) -> None:
    """M7-16 rule 11: U over TGA and O over TAG are not skipped; the codon is a pileup site."""
    stop = "TGA" if letter == "U" else "TAG"
    lay = minus(overrides={2: stop})
    assert lay.key(2) == 127
    h, _pep = _stop_case(lay, letter)
    pile, _ = run([h], lay)
    assert ns(pile) == counts()
    assert has_row(pile, 127, "-", "c2")
    assert has_row(pile, 130, "-", "c2") and has_row(pile, 124, "-", "c2")


def test_ex7_u_over_tga_across_the_origin_wrapped_and_unwrapped_count_once() -> None:
    """M7-16 rule 11: the site is a pileup row (not rare-skipped), counted once for both hits."""
    lay = origin(overrides={3: "TGA"})
    seq = lay.sequence()
    assert lay.codon(3) == 199 and seq[199] + seq[:2] == "TGA"
    pep = with_x(P, 3, "U")
    unwrapped, plain = hit(al, lay, pep), hit(al, lay, pep, shift=-lay.length)
    assert unwrapped.ops[3].ref == "*"
    pile, _ = run([unwrapped, plain], lay)
    assert ns(pile) == counts()
    r = row(pile, 199, "+", "c3")
    assert (r["ref_aa"], r["alt_counts"], r["call"]) == ("*", "U:1", "ALT")


def test_x_across_the_origin_wrapped_and_unwrapped_count_once() -> None:
    lay = origin()
    pep = with_x(P, 3)
    pile, _ = run([hit(al, lay, pep), hit(al, lay, pep, shift=-lay.length)], lay)
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    assert not has_row(pile, 199, "+", "c3")
    assert row(pile, 2, "+", "c3")["n_ref"] == 1  # residue 4 at [202, 205) -> key 2


def test_ex9_swap_then_x() -> None:
    pep = "KAXREFGHILWD"
    base = hit(al, LIN, pep)
    swap = al.ResidueOp("swap", (0, 2), (0, 3), "AK", "KA")
    h = al.HitAlignment(pep, "c1", "+", LIN.frame(), (swap, *base.ops[2:]), None)
    pile, sites = run([h])
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    r = row(pile, 0)
    assert (r["alt_counts"], r["call"]) == ("AK>KA:1", "ALT")
    assert not has_row(pile, 3), "the swap's second codon gets no support (M7-02 decision 4)"
    assert not has_row(pile, 6)
    (s,) = [x for _, x in sites.iterrows()]
    assert (int(s["codon_start"]), s["alt_type"], s["confidence"]) == (0, "swap", "candidate")
    assert s["ref_codon"] == "GCTAAA"


def test_ex10_minus_strand_x_codon_key_is_127() -> None:
    lay = minus()
    pep = "AKXREFGHILWD"
    pile, sites = run([hit(al, lay, pep)], lay)
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    assert not has_row(pile, 127, "-", "c2")
    assert codons_of(pile) == [lay.key(i) for i in range(11, -1, -1) if i != 2]
    assert set(pile["call"]) == {"REF"} and len(sites) == 0


# --------------------------------------------------------------------------- edge cases


@pytest.mark.parametrize(("i", "codon"), [(0, 0), (11, 33)])
def test_x_at_the_first_or_last_residue(i: int, codon: int) -> None:
    pep = with_x(P, i)
    pile, _ = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    assert codons_of(pile) == [c for c in ALL12 if c != codon]


def test_all_x_alignment_gives_no_rows() -> None:
    pep = "X" * 12
    pile, sites = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(12, 12, 1, {"X": 12}, ambiguous=12, examples=[pep])
    assert len(pile) == 0 and len(sites) == 0


def test_two_x_in_one_peptide() -> None:
    pep = "AKXXEFGHILWD"
    pile, _ = run([hit(al, LIN, pep)])
    assert ns(pile) == counts(2, 2, 1, {"X": 2}, ambiguous=2, examples=[pep])
    assert not has_row(pile, 6) and not has_row(pile, 9)


@pytest.mark.parametrize(("letter", "i"), [("B", 11), ("Z", 4), ("J", 9), ("O", 7)])
def test_b_z_j_o_over_their_look_alike_codons_are_skipped(letter: str, i: int) -> None:
    """B over D, Z over E, J over L: never REF (the mapper's NOMATCH); O over a sense codon."""
    pep = with_x(P, i, letter)
    pile, sites = run([hit(al, LIN, pep)])
    rule = {"rare": 1} if letter in "OU" else {"ambiguous": 1}
    assert ns(pile) == counts(1, 1, 1, {letter: 1}, examples=[pep], **rule)
    assert not has_row(pile, 3 * i) and len(sites) == 0


def test_mixed_letters_by_residue_and_by_rule() -> None:
    peps = [with_x(P, 0, "Z"), with_x(P, 1, "B"), with_x(P, 2, "X"), with_x(P, 3, "U")]
    pile, _ = run([hit(al, LIN, p) for p in peps])
    got = ns(pile)
    assert got == counts(4, 4, 4, {"B": 1, "U": 1, "X": 1, "Z": 1}, ambiguous=3, rare=1,
                         examples=sorted(peps))  # fmt: skip
    assert list(got["by_residue"]) == ["B", "U", "X", "Z"]
    assert list(got) == KEYS and list(got["by_rule"]) == ["ambiguous", "rare"]


def test_skip_resolved_by_another_hit_of_the_same_peptide_still_counts() -> None:
    """Q6: two hits of one folded peptide at one codon, one skipped and one match."""
    pep = "AKXREFGHILWD"
    skipped = hit(al, LIN, pep)
    matched = al.HitAlignment(pep, "c1", "+", LIN.frame(),
                              (al.ResidueOp("match", (2, 3), (6,), "C", "C"),), None)  # fmt: skip
    pile, _ = run([skipped, matched])
    assert ns(pile)["skipped_sites"] == 1
    r = row(pile, 6)
    assert (int(r["n_ref"]), r["call"]) == (1, "REF")


def test_examples_are_sorted_and_capped_at_ten() -> None:
    peps = [with_x(P, i) for i in range(12)]
    pile, _ = run([hit(al, LIN, p) for p in peps])
    got = ns(pile)
    assert got == counts(12, 12, 12, {"X": 12}, ambiguous=12, examples=sorted(peps)[:10])


def test_counts_are_json_native_with_sorted_keys() -> None:
    pile, _ = run([hit(al, LIN, "AKXREFGHILWD"), hit(al, LIN, "AKCREFGHJLWD")])
    got = ns(pile)
    assert json.loads(json.dumps(got)) == got
    for k in ("skipped_sites", "skipped_codons", "skipped_peptides"):
        assert type(got[k]) is int
    assert all(type(v) is int for v in got["by_residue"].values())
    assert all(type(v) is int for v in got["by_rule"].values())
    assert all(type(e) is str for e in got["examples"])
    assert list(got) == KEYS and list(got["by_residue"]) == ["J", "X"]


def test_empty_input_is_all_zeros() -> None:
    pile, sites = run([])
    assert ns(pile) == ZERO and len(pile) == 0 and len(sites) == 0


def test_all_standard_input_is_all_zeros() -> None:
    pile, _ = run([hit(al, LIN, P), hit(al, LIN, "AKDREFGHLLWD")])
    assert ns(pile) == ZERO


def test_defensive_swap_with_a_nonstandard_letter_is_skipped_not_raised() -> None:
    pep = "KXCREFGHILWD"
    base = hit(al, LIN, pep)
    swap = al.ResidueOp("swap", (0, 2), (0, 3), "AK", "KX")
    h = al.HitAlignment(pep, "c1", "+", LIN.frame(), (swap, *base.ops[2:]), None)
    pile, sites = run([h])
    assert ns(pile) == counts(1, 1, 1, {"KX": 1}, ambiguous=1, examples=[pep])
    assert not has_row(pile, 0) and len(sites) == 0


def test_literal_j_is_not_folded_into_l() -> None:
    """Q5: a J-peptide and an L-peptide stay two distinct peptides."""
    pj, pl_ = "AKCREFGHJLWD", "AKCREFGHLLWD"
    pile, _ = run([hit(al, LIN, pj), hit(al, LIN, pl_)])
    assert ns(pile)["skipped_peptides"] == 1
    assert int(row(pile, 27)["n_ref"]) == 2
    r = row(pile, 24)
    assert (int(r["n_ref"]), int(r["n_il"]), r["call"]) == (0, 1, "no_call")


@pytest.mark.parametrize("relation", ["same", "related"])
def test_relation_does_not_change_the_skip(relation: str) -> None:
    pep = "AKXREFGHILWD"
    pile, sites = run([hit(al, LIN, pep)], relation=relation)
    assert ns(pile) == counts(1, 1, 1, {"X": 1}, ambiguous=1, examples=[pep])
    assert not has_row(pile, 6) and len(sites) == 0


# --------------------------------------------------------------------------- eligible (rule 3)


def test_database_other_under_same_strips_an_x_op_and_does_not_count_it() -> None:
    h = hit(al, LIN, "AKXREFGHILWD")
    sup = support([h], source="database", search_db="other")
    pile, _ = run([h], sup=sup)
    assert ns(pile) == ZERO
    assert not has_row(pile, 6)


def test_database_other_under_same_strips_a_u_over_tga_op_too() -> None:
    """Lead decision 1: no keep-condition for ``ref == '*'``."""
    lay = linear(overrides={2: "TGA"})
    h, _ = _stop_case(lay, "U")
    el = pl.eligible([h], support([h], "database", "other"), cfg("same"))
    assert [o.kind for o in el[0].ops] == ["match"] * 11
    pile, _ = run([h], lay, sup=support([h], "database", "other"))
    assert ns(pile) == ZERO and not has_row(pile, 6)


def test_database_other_under_related_counts_the_skip() -> None:
    h = hit(al, LIN, "AKXREFGHILWD")
    pile, _ = run([h], sup=support([h], "database", "other"), relation="related")
    assert ns(pile)["skipped_sites"] == 1


# --------------------------------------------------------------------------- frameshift cover


@pytest.mark.parametrize("x_at", [5, 6])
def test_skipped_codon_still_counts_toward_frameshift_spanning_cover(x_at: int) -> None:
    """Rule 4 / Q4: the X peptide spans codons 15 and 18 contiguously, so it is a REF peptide of
    the +1 frameshift whose last prefix codon is 15 (the X sits at codon 15 or 18)."""
    seq = LIN.sequence()
    pre = [al.ResidueOp("match", (k, k + 1), (3 * k,), P[k], P[k]) for k in range(6)]
    suf_codons = (19, 22, 25)
    suf_aa = [al.translate_locus(seq, LIN.replicon(), c, c + 3, "+") for c in suf_codons]
    suf = [al.ResidueOp("match", (6 + j, 7 + j), (c,), a, a)
           for j, (c, a) in enumerate(zip(suf_codons, suf_aa, strict=True))]  # fmt: skip
    fs_pep = P[:6] + "".join(suf_aa)
    fs = al.HitAlignment(fs_pep, "c1", "+", 1, (*pre, *suf),
                         al.FrameshiftEvent(6, 15, 19, 1))  # fmt: skip
    span = with_x(P, x_at)
    pile, sites = run([fs, hit(al, LIN, span)])
    assert ns(pile)["skipped_sites"] == 1
    (s,) = [x for _, x in sites.iterrows()]
    assert s["alt_type"] == "frameshift"
    assert (int(s["codon_start"]), int(s["n_alt"]), int(s["n_ref"]), s["call"]) == (
        15, 1, 1, "mixed",
    )  # fmt: skip


# --------------------------------------------------------------------------- properties

WINDOW = st.tuples(st.integers(0, 11), st.integers(0, 11)).map(lambda t: (min(t), max(t) + 1))
STANDARD = "ACDEFGHIKLMNPQRSTVWY"


@st.composite
def standard_peptides(draw: st.DrawFn) -> list[tuple[int, str]]:
    """(lo, peptide) pairs over G1 of c1, with random standard substitutions."""
    out = []
    for _ in range(draw(st.integers(1, 6))):
        lo, hi = draw(WINDOW)
        pep = list(P[lo:hi])
        for k in range(len(pep)):
            if draw(st.booleans()) and draw(st.booleans()):
                pep[k] = draw(st.sampled_from(STANDARD))
        out.append((lo, "".join(pep)))
    return out


def _hits(pairs: Sequence[tuple[int, str]]) -> list[Any]:
    return [hit(al, LIN, pep, lo) for lo, pep in pairs]


@SETTINGS
@given(standard_peptides())
def test_property_standard_only_input_is_all_zeros(pairs: list[tuple[int, str]]) -> None:
    pile, _ = run(_hits(pairs))
    assert ns(pile) == ZERO


def _codon_view(pile: pd.DataFrame, sites: pd.DataFrame, codon: int) -> tuple[Any, ...]:
    r = pile[pile["codon_start"] == codon][list(pl.PILEUP_COLUMNS)]
    s = sites[sites["codon_start"] == codon].drop(columns="variant_id")
    return (r.to_dict("records"), s.to_dict("records"))


def _records_at(pile: pd.DataFrame, sites: pd.DataFrame, codon: int) -> list[tuple[Any, ...]]:
    recs, _ = vcf_records(sites, pile, LIN.genome(), {"c1": LIN.replicon()})
    return [(r.ref, r.alts, r.ad, r.dp) for r in recs if r.pos0 == codon]


@SETTINGS
@given(standard_peptides(), st.data())
def test_property_a_nonstandard_twin_changes_nothing_at_its_codon(
    pairs: list[tuple[int, str]], data: st.DataObject
) -> None:
    """Adding a peptide that differs from an existing one only by a nonstandard letter leaves
    the call, confidence, n_ref, n_il, n_alt, AD and DP at that letter's codon unchanged."""
    lo, pep = data.draw(st.sampled_from(pairs))
    k = data.draw(st.integers(0, len(pep) - 1))
    twin = with_x(pep, k, data.draw(st.sampled_from(NONSTANDARD)))
    codon = 3 * (lo + k)
    before = run(_hits(pairs))
    after = run([*_hits(pairs), hit(al, LIN, twin, lo)])
    assert ns(after[0])["skipped_sites"] == 1
    assert _codon_view(*after, codon) == _codon_view(*before, codon)
    assert _records_at(*after, codon) == _records_at(*before, codon)


@st.composite
def any_peptides(draw: st.DrawFn) -> list[tuple[int, str]]:
    out = []
    for _ in range(draw(st.integers(1, 6))):
        lo, hi = draw(WINDOW)
        pep = list(P[lo:hi])
        for k in range(len(pep)):
            if draw(st.booleans()):
                pep[k] = draw(st.sampled_from(STANDARD + NONSTANDARD + "JJJ"))
        out.append((lo, "".join(pep)))
    return out


@SETTINGS
@given(any_peptides(), st.sampled_from(["same", "related"]))
def test_property_no_vcf_record_has_alt_equal_to_ref(
    pairs: list[tuple[int, str]], relation: str
) -> None:
    pile, sites = run(_hits(pairs), relation=relation)
    recs, _ = vcf_records(sites, pile, LIN.genome(), {"c1": LIN.replicon()})
    for r in recs:
        assert r.ref not in r.alts, f"record at {r.pos0}: ALT == REF {r.ref}"
    assert all(set(str(x)) <= set(STANDARD) | {"J"} for x in sites["alt_aa"])
