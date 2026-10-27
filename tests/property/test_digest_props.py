"""Hypothesis properties of the in-silico digest (issue M3-22; GIPPA_SPEC.md §11.2, §11.6).

- Acceptance: "concatenating missed = 0 products with the short fragments reconstructs the
  protein" (with unbounded lengths the missed = 0 products tile the protein; with the default
  7-30 bounds they are exactly the in-bounds tiles).
- The result equals a brute-force oracle for both enzymes, every ``missed`` and length bound.
- Monotone in ``missed``; trypsin cut sites are a subset of trypsin_p cut sites.
- Lead decision 1: the initiator M is kept; no start-1 product unless position 0 is a cut site.
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.mapping.digest import digest

# K, R and P over-represented so cut sites, KP/RP pairs and runs of sites are common
_ALPHABET = "ACDEFGHILMNQSTVWY" + "KKKRRRPPP"
proteins = st.text(alphabet=_ALPHABET, min_size=0, max_size=200)
enzymes = st.sampled_from(["trypsin_p", "trypsin"])


def _sites(p: str, enzyme: str) -> list[int]:
    """Product end positions: after K/R (trypsin: not before P), plus the protein end."""
    ends = [
        i + 1
        for i in range(len(p) - 1)
        if p[i] in "KR" and (enzyme == "trypsin_p" or p[i + 1] != "P")
    ]
    return [*ends, len(p)] if p else []


def oracle(p: str, enzyme: str, missed: int, lo: int, hi: int) -> list[tuple[int, int, str]]:
    ends = _sites(p, enzyme)
    starts = [0, *ends[:-1]]
    out = set()
    for i, s in enumerate(starts):
        for j in range(i, min(i + missed + 1, len(ends))):
            e = ends[j]
            if lo <= e - s <= hi:
                out.add((s, e, p[s:e]))
    return sorted(out)


def _sorted(products: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    return sorted((int(s), int(e), str(x)) for s, e, x in products)


@settings(max_examples=300, deadline=None)
@given(p=proteins, enzyme=enzymes)
def test_missed_0_unbounded_products_tile_the_protein(p: str, enzyme: str) -> None:
    got = _sorted(digest(p, enzyme=enzyme, missed=0, min_len=1, max_len=10**6))  # type: ignore[arg-type]
    assert "".join(x for _, _, x in got) == p
    pos = 0
    for s, e, _ in got:
        assert s == pos and e > s
        pos = e
    assert pos == len(p)


@settings(max_examples=300, deadline=None)
@given(p=proteins, enzyme=enzymes)
def test_missed_0_products_plus_short_fragments_reconstruct(p: str, enzyme: str) -> None:
    tiles = _sorted(digest(p, enzyme=enzyme, missed=0, min_len=1, max_len=10**6))  # type: ignore[arg-type]
    kept = _sorted(digest(p, enzyme=enzyme, missed=0))  # type: ignore[arg-type]
    assert kept == [t for t in tiles if 7 <= t[1] - t[0] <= 30]
    fragments = [t for t in tiles if not 7 <= t[1] - t[0] <= 30]
    assert "".join(x for _, _, x in sorted(kept + fragments)) == p


@settings(max_examples=400, deadline=None)
@given(
    p=proteins,
    enzyme=enzymes,
    missed=st.integers(0, 3),
    lo=st.integers(1, 12),
    span=st.integers(0, 40),
)
def test_equals_brute_force_oracle(p: str, enzyme: str, missed: int, lo: int, span: int) -> None:
    got = digest(p, enzyme=enzyme, missed=missed, min_len=lo, max_len=lo + span)  # type: ignore[arg-type]
    assert len(got) == len({(s, e) for s, e, _ in got}), "duplicate (start, end)"
    assert _sorted(got) == oracle(p, enzyme, missed, lo, lo + span)


@settings(max_examples=200, deadline=None)
@given(p=proteins, enzyme=enzymes, missed=st.integers(0, 3))
def test_products_are_slices_within_bounds(p: str, enzyme: str, missed: int) -> None:
    sites = set(_sites(p, enzyme))
    for s, e, x in digest(p, enzyme=enzyme, missed=missed):  # type: ignore[arg-type]
        assert 0 <= s < e <= len(p)
        assert p[s:e] == x
        assert 7 <= len(x) <= 30
        assert s == 0 or s in sites
        assert e in sites
        assert sum(1 for k in sites if s < k < e) <= missed


@settings(max_examples=200, deadline=None)
@given(p=proteins, enzyme=enzymes, missed=st.integers(0, 2))
def test_monotone_in_missed(p: str, enzyme: str, missed: int) -> None:
    a = set(_sorted(digest(p, enzyme=enzyme, missed=missed)))  # type: ignore[arg-type]
    b = set(_sorted(digest(p, enzyme=enzyme, missed=missed + 1)))  # type: ignore[arg-type]
    assert a <= b


@settings(max_examples=200, deadline=None)
@given(p=proteins)
def test_trypsin_ends_subset_of_trypsin_p_ends(p: str) -> None:
    unb = {"missed": 0, "min_len": 1, "max_len": 10**6}
    ends_t = {e for _, e, _ in digest(p, enzyme="trypsin", **unb)}  # type: ignore[arg-type]
    ends_tp = {e for _, e, _ in digest(p, enzyme="trypsin_p", **unb)}  # type: ignore[arg-type]
    assert ends_t <= ends_tp
    if "KP" not in p and "RP" not in p:
        assert ends_t == ends_tp


@settings(max_examples=200, deadline=None)
@given(rest=proteins, enzyme=enzymes, missed=st.integers(0, 3))
def test_initiator_met_kept_no_met_excised_products(rest: str, enzyme: str, missed: int) -> None:
    p = "M" + rest
    got = digest(p, enzyme=enzyme, missed=missed, min_len=1, max_len=10**6)  # type: ignore[arg-type]
    assert all(s != 1 for s, _, _ in got)  # M is never a cut site
    assert any(s == 0 and x.startswith("M") for s, _, x in got)
