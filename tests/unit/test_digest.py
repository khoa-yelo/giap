"""In-silico digestion (issue M3-22 examples 1-3; GIPPA_SPEC.md §11.6 TE1/TE8,
``abundance.enzyme``).

Interface (docs/specs/M3-22-digest-and-te-proxies.md "Public interface")::

    gippa.mapping.digest.digest(protein, enzyme="trypsin_p", missed=1, min_len=7, max_len=30)
        -> list[tuple[int, int, str]]          # (start, end, peptide), 0-based half-open

Semantics pinned here (test author):

* ``trypsin_p`` cuts after every K/R (also before P, the DIA-NN default); ``trypsin`` cuts after
  K/R except before P. The protein end is a cut site; no empty products.
* ``missed`` is the maximum number of internal cut sites a product may span (0..missed).
* Length bounds are inclusive: ``min_len <= end - start <= max_len``.
* ``protein[start:end] == peptide`` for every product; each (start, end) appears once; order is
  not pinned (compared as sorted lists).
* Lead decision 1 (2026-10-20): the protein is digested as Bakta writes it, with the initiator
  ``M``; no Met-excised (start 1) products are generated unless position 0 is itself a cut site.
* An unknown enzyme name is a ``ValueError``.
"""

from __future__ import annotations

import numpy as np
import pytest

from gippa.mapping.digest import digest

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _sorted(products: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    return sorted((int(s), int(e), str(p)) for s, e, p in products)


def _random_protein(seed: int, n: int) -> str:
    idx = np.random.default_rng(seed).integers(0, len(_AA), size=n - 1)
    return "M" + "".join(_AA[i] for i in idx)


# ------------------------------------------------------------------ example 1 (trypsin_p)

EX = "MKAAKPRGGRLLLK"  # cut sites after K1, K4, R6, R9 (and the protein end)


def test_example_1_trypsin_p_missed_1_min_len_3() -> None:
    got = _sorted(digest(EX, min_len=3))
    assert got == [
        (0, 5, "MKAAK"),
        (2, 5, "AAK"),
        (2, 7, "AAKPR"),  # missed cleavage at K4 (which trypsin_p cuts, before P)
        (5, 10, "PRGGR"),
        (7, 10, "GGR"),
        (7, 14, "GGRLLLK"),
        (10, 14, "LLLK"),
    ]


def test_example_1_explicit_arguments_equal_defaults() -> None:
    assert _sorted(digest(EX, "trypsin_p", 1, 3, 30)) == _sorted(digest(EX, min_len=3))
    assert _sorted(digest(EX)) == _sorted(
        digest(EX, enzyme="trypsin_p", missed=1, min_len=7, max_len=30)
    )


def test_example_1_defaults_keep_only_7_to_30() -> None:
    assert _sorted(digest(EX)) == [(7, 14, "GGRLLLK")]


def test_example_1_trypsin_p_cuts_before_proline() -> None:
    got = _sorted(digest(EX, missed=0, min_len=1))
    assert got == [(0, 2, "MK"), (2, 5, "AAK"), (5, 7, "PR"), (7, 10, "GGR"), (10, 14, "LLLK")]


# ------------------------------------------------------------------ example 2 (trypsin)


def test_example_2_trypsin_no_cut_before_proline() -> None:
    got = _sorted(digest(EX, enzyme="trypsin", min_len=3))
    assert got == [
        (0, 7, "MKAAKPR"),
        (2, 7, "AAKPR"),
        (2, 10, "AAKPRGGR"),
        (7, 10, "GGR"),
        (7, 14, "GGRLLLK"),
        (10, 14, "LLLK"),
    ]
    assert all(e != 5 for _, e, _ in got)  # K4 is followed by P: never a product end


def test_example_2_trypsin_missed_0() -> None:
    got = _sorted(digest(EX, enzyme="trypsin", missed=0, min_len=1))
    assert got == [(0, 2, "MK"), (2, 7, "AAKPR"), (7, 10, "GGR"), (10, 14, "LLLK")]


def test_trypsin_rp_is_not_cut_either() -> None:
    p = "AAAAAAARPAAAAAAAK"
    assert _sorted(digest(p, enzyme="trypsin", missed=0, min_len=1)) == [(0, 17, p)]
    assert _sorted(digest(p, enzyme="trypsin_p", missed=0, min_len=1)) == [
        (0, 8, "AAAAAAAR"),
        (8, 17, "PAAAAAAAK"),
    ]


# ------------------------------------------------------------------ example 3 (300 aa, missed 0)


def _cut_ends(protein: str, enzyme: str) -> list[int]:
    ends = [
        i + 1
        for i, a in enumerate(protein[:-1])
        if a in "KR" and (enzyme == "trypsin_p" or protein[i + 1] != "P")
    ]
    return [*ends, len(protein)]


@pytest.mark.parametrize("enzyme", ["trypsin_p", "trypsin"])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_example_3_missed_0_products_tile_between_cut_sites(enzyme: str, seed: int) -> None:
    p = _random_protein(seed, 300)
    ends = _cut_ends(p, enzyme)
    segments = list(zip([0, *ends[:-1]], ends, strict=True))
    got = _sorted(digest(p, enzyme=enzyme, missed=0))  # type: ignore[arg-type]
    assert got, "a 300-aa random protein yields some 7-30 aa products"
    for s, e, pep in got:
        assert 7 <= len(pep) <= 30
        assert (s, e) in segments
        assert p[s:e] == pep
    expected = [(s, e, p[s:e]) for s, e in segments if 7 <= e - s <= 30]
    assert got == expected


# ------------------------------------------------------------------ initiator Met (decision 1)


def test_initiator_met_is_kept_and_no_met_excised_product() -> None:
    p = "MAAAAAAAKWWWWWWWR"
    got = _sorted(digest(p, missed=0))
    assert got == [(0, 9, "MAAAAAAAK"), (9, 17, "WWWWWWWR")]
    assert all(s != 1 for s, _, _ in digest(p, missed=2, min_len=1))


def test_met_followed_by_cut_site_starts_next_product_at_1_only_via_cleavage() -> None:
    # position 0 is not K/R, so start 1 never appears; position 1 (K) gives start 2
    p = "MKAAAAAAAAR"
    starts = {s for s, _, _ in digest(p, missed=1, min_len=1)}
    assert starts == {0, 2}


# ------------------------------------------------------------------ edges


def test_terminal_cut_site_gives_no_empty_product() -> None:
    p = "AAAAAAAKK"
    got = _sorted(digest(p, missed=1, min_len=1))
    assert got == [(0, 8, "AAAAAAAK"), (0, 9, "AAAAAAAKK"), (8, 9, "K")]
    assert all(e > s for s, e, _ in got)


def test_consecutive_cut_sites_single_residue_fragments() -> None:
    got = _sorted(digest("KRK", missed=0, min_len=1))
    assert got == [(0, 1, "K"), (1, 2, "R"), (2, 3, "K")]


def test_no_cut_site_returns_whole_protein_within_bounds() -> None:
    assert _sorted(digest("ACDEFGHILMNPQ")) == [(0, 13, "ACDEFGHILMNPQ")]
    assert digest("ACDEFG") == []  # 6 aa < min_len 7
    assert digest("A" * 31) == []  # 31 aa > max_len 30


def test_length_bounds_are_inclusive() -> None:
    p = "AAAAAAK" + "C" * 29 + "R"  # 7 aa, then 30 aa
    got = _sorted(digest(p, missed=0))
    assert got == [(0, 7, "AAAAAAK"), (7, 37, "C" * 29 + "R")]


def test_empty_protein() -> None:
    assert digest("") == []


def test_missed_2_spans_two_internal_sites() -> None:
    p = "AAK" + "CCK" + "DDK" + "EER"
    got = {pep for _, _, pep in digest(p, missed=2, min_len=1)}
    assert "AAKCCKDDK" in got and "CCKDDKEER" in got
    assert p not in got  # three internal sites
    assert p in {pep for _, _, pep in digest(p, missed=3, min_len=1)}


def test_unknown_enzyme_is_rejected() -> None:
    with pytest.raises(ValueError):
        digest(EX, enzyme="chymotrypsin")  # type: ignore[arg-type]


@pytest.mark.parametrize("missed", [-1, -2])
def test_negative_missed_is_rejected(missed: int) -> None:
    with pytest.raises(ValueError):
        digest(EX, missed=missed)


def test_missed_zero_is_accepted() -> None:
    assert digest(EX, missed=0, min_len=1)
