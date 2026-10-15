"""Hypothesis properties of the M2-17 codon-usage functions (docs/specs/M2-17-codon-usage.md).

Every pure function equals the independent oracle (tests/helpers/codon_oracle.py) on arbitrary
input; invariants: shares sum to 1 per observed family, single-codon families are never rare,
``rare_codons`` is monotone in the cutoff, ``at_2_10`` depends only on nt 3..30, codon 1 and the
last codon never change a score, and ``codon_usage`` does not depend on gene order.
"""

from __future__ import annotations

import dataclasses
import importlib
import math
from types import ModuleType

from helpers import codon_oracle as O
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.model.coords import Interval
from giap.model.genes import Gene

ALL64 = sorted(O.CODE11)
codon = st.sampled_from(ALL64)
messy_codon = st.one_of(codon, st.sampled_from(["NNN", "CTN", "ANA", "GCR", "YAA"]))
cds_body = st.lists(st.sampled_from(O.SENSE), min_size=8, max_size=80).map("".join)
nt_any = st.text(alphabet="ACGTN", min_size=0, max_size=60)


def cu() -> ModuleType:
    return importlib.import_module("giap.annotate.codon_usage")


@settings(deadline=None)
@given(st.lists(st.lists(messy_codon, max_size=40).map("".join), max_size=6))
def test_shares_equal_oracle_and_sum_to_one(seqs: list[str]) -> None:
    got = cu().synonymous_shares(seqs)
    want = O.shares(seqs)
    assert set(got) == set(want)
    for c, v in want.items():
        assert math.isclose(got[c], v, rel_tol=1e-12, abs_tol=1e-15)
    for fam in O.FAMILY.values():
        if fam[0] in got:
            assert math.isclose(sum(got[c] for c in fam), 1.0, rel_tol=1e-12)


@settings(deadline=None)
@given(st.lists(st.lists(codon, max_size=40).map("".join), max_size=6))
def test_shares_do_not_depend_on_sequence_order(seqs: list[str]) -> None:
    assert cu().synonymous_shares(seqs) == cu().synonymous_shares(list(reversed(seqs)))


@settings(deadline=None)
@given(
    st.dictionaries(st.sampled_from(O.SENSE), st.floats(0, 1), max_size=61),
    st.floats(0, 1),
    st.floats(0, 1),
)
def test_rare_codons_monotone_strict_and_never_single_families(
    shares: dict[str, float], c1: float, c2: float
) -> None:
    for c in ("ATG", "TGG"):
        shares[c] = 1.0
    lo, hi = sorted((c1, c2))
    r_lo, r_hi = cu().rare_codons(shares, lo), cu().rare_codons(shares, hi)
    assert r_lo <= r_hi
    assert r_hi == O.rare(shares, hi)
    assert not {"ATG", "TGG"} & r_hi


@settings(deadline=None)
@given(nt_any)
def test_at_2_10_equals_oracle(nt: str) -> None:
    got, want = cu().at_2_10(nt), O.at_2_10(nt)
    assert (got is None) == (want is None)
    if want is not None:
        assert got is not None and math.isclose(got, want) and 0 <= got <= 1


@settings(deadline=None)
@given(
    nt_any.filter(lambda s: len(s) >= 30),
    st.text("ACGT", min_size=3, max_size=3),
    st.text("ACGT", max_size=30),
)
def test_at_2_10_depends_only_on_nt_3_to_30(nt: str, first: str, tail: str) -> None:  # fmt: skip
    m = cu()
    assert m.at_2_10(first + nt[3:30] + tail) == m.at_2_10(nt)


@settings(deadline=None)
@given(st.text("ACGT", min_size=30, max_size=60), st.frozensets(st.sampled_from(O.SENSE)))
def test_rare_2_10_equals_oracle(nt: str, rare: frozenset[str]) -> None:
    got = cu().rare_2_10(nt, rare)
    assert got == O.rare_2_10(nt, rare)
    assert got is not None and 0 <= got <= 9


@settings(deadline=None)
@given(st.lists(codon, min_size=1, max_size=60).map("".join))
def test_gc3_equals_oracle(seq: str) -> None:
    got = cu().gc3(seq)
    assert math.isclose(got, O.gc3(seq)) and 0 <= got <= 1


def _gene(locus: str, nt: str, k: int, **kw: object) -> Gene:
    return Gene(
        locus_tag=locus,
        feature_type=kw.pop("feature_type", "cds"),  # type: ignore[arg-type]
        interval=Interval("c1", 1000 * k, 1000 * k + len(nt), "+"),
        aa="M",
        nt=nt,
        start_codon=nt[:3],
        **kw,  # type: ignore[arg-type]
    )


@settings(deadline=None)
@given(
    cds_body,
    st.sampled_from(["ATG", "GTG", "TTG", "CTG"]),
    st.sampled_from(["TAA", "TAG", "TGA"]),
    st.booleans(),
    st.sampled_from(["cds", "sorf", "pseudogene"]),
)
def test_scored_sequence_equals_oracle(
    body: str, start: str, stop: str, partial: bool, kind: str
) -> None:
    nt = start + body + stop
    g = _gene("g", nt, 0, feature_type=kind, partial=partial)
    assert cu().scored_sequence(g) == O.scored(nt, kind, partial)


genes_st = st.lists(cds_body, min_size=3, max_size=8)


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    genes_st,
    st.integers(1, 3),
    st.dictionaries(
        st.sampled_from(["CAT", "CAG", "GCC", "TTC", "GTC"]), st.integers(1, 5), max_size=5
    ),
    st.randoms(use_true_random=False),
)
def test_codon_usage_oracle_ranges_and_order(
    bodies: list[str], n_rib: int, trna: dict[str, int], rnd: object
) -> None:  # fmt: skip
    genes = []
    for k, body in enumerate(bodies):
        rib = k < n_rib
        genes.append(
            _gene(
                f"L{k:02d}", "ATG" + body + "TAA", k,
                gene=f"rps{chr(65 + k)}" if rib else None,
                product="30S ribosomal protein" if rib else "hypothetical protein",
            )
        )  # fmt: skip
    feats, flags = cu().codon_usage(genes, trna)
    assert flags == ()
    ref = bodies[:n_rib]
    eligible = list(bodies)
    rare = O.rare(O.shares(eligible))
    for g, body in zip(genes, bodies, strict=True):
        f = feats[g.locus_tag]
        cai = O.cai(body, ref)
        assert cai is not None and f.cai is not None
        assert math.isclose(f.cai, cai, rel_tol=1e-9)
        assert math.isclose(f.fop, O.fop(body, ref), rel_tol=1e-9)
        assert math.isclose(f.gc3, O.gc3(body), rel_tol=1e-12)
        assert f.rare_codons_2_10 == O.rare_2_10("ATG" + body + "TAA", rare)
        t = O.tai(body, trna) if trna else None  # None also when no codon is readable
        assert (f.tai is None) == (t is None)
        assert t is None or (f.tai is not None and math.isclose(f.tai, t, rel_tol=1e-9))
        assert 20 <= f.enc <= 61 and math.isfinite(f.cpb)
    shuffled = list(genes)
    rnd.shuffle(shuffled)  # type: ignore[attr-defined]
    again, _ = cu().codon_usage(shuffled, trna)
    for lt, f in feats.items():
        for n in (x.name for x in dataclasses.fields(f)):
            a, b = getattr(f, n), getattr(again[lt], n)
            assert a == b or (a is not None and b is not None and math.isclose(a, b, rel_tol=1e-12))
