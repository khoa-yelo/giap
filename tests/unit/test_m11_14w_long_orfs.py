"""M11-14w ``long_orf_cover``: examples 1-5, the brute-force Hypothesis property (AC 2), strand
symmetry, lowercase input and the non-ACGT break (docs/specs/M11-14w-wrong-genetic-code-warning.md
"Signal", "Examples", "Edge cases"; GIAP_SPEC.md §6.5).

Interface assumed (spec "Where it runs" item 1):
``giap.preflight.genetic_code.long_orf_cover(seq, stops, min_codons) -> tuple[int, int]``
returning ``(bp covered, number of runs)`` for one record over all six frames; ``stops`` is
passed as a ``frozenset`` of upper-case codons (``S11`` / ``S4`` below). The reference is
``helpers.genetic_code_fixture.reference_cover`` (codon-by-codon, pure Python).
"""

from __future__ import annotations

import pytest
from helpers.genetic_code_fixture import S4, S11, reference_cover, reference_runs, revcomp
from hypothesis import given, settings
from hypothesis import strategies as st


def cover(seq: str, stops: frozenset[str], min_codons: int) -> tuple[int, int]:
    from giap.preflight.genetic_code import long_orf_cover

    out = long_orf_cover(seq, stops, min_codons)
    assert isinstance(out, tuple) and len(out) == 2
    return int(out[0]), int(out[1])


def frame0_plus(seq: str, stops: frozenset[str], k: int) -> list[tuple[int, int]]:
    return [(a, b) for s, f, a, b in reference_runs(seq, stops, k) if s == "+" and f == 0]


EX1 = "ATG" + "GCC" * 8 + "TGA" + "GCC" * 8 + "TAA"
EX3 = "ATG" + "GCC" * 4 + "NNN" + "GCC" * 4 + "TAA"


# ---------------------------------------------------------------- examples 1-5


def test_example1_tga_read_through_only_under_s4() -> None:
    assert len(EX1) == 57
    # the reference reproduces the spec's named frame
    assert frame0_plus(EX1, S11, 1) == [(0, 27), (30, 54)]
    assert frame0_plus(EX1, S11, 10) == []
    assert frame0_plus(EX1, S4, 10) == [(0, 54)]
    for stops in (S11, S4):
        assert cover(EX1, stops, 10) == reference_cover(EX1, stops, 10)
    bp4, n4 = cover(EX1, S4, 10)
    assert bp4 >= 54 and n4 >= 1


def test_example2_minus_strand_maps_to_record_coordinates() -> None:
    rc = revcomp(EX1)
    minus = [(a, b) for s, f, a, b in reference_runs(rc, S4, 10) if s == "-"]
    assert (3, 57) in minus
    assert not [r for r in reference_runs(rc, S11, 10) if r[0] == "-" and r[2:] == (3, 57)]
    for stops in (S11, S4):
        assert cover(rc, stops, 10) == reference_cover(rc, stops, 10)
        assert cover(rc, stops, 10) == cover(EX1, stops, 10)


def test_example3_n_codon_breaks_and_is_not_covered() -> None:
    assert frame0_plus(EX3, S11, 4) == [(0, 15), (18, 30)]
    assert frame0_plus(EX3, S4, 4) == [(0, 15), (18, 30)]
    for stops in (S11, S4):
        assert cover(EX3, stops, 4) == reference_cover(EX3, stops, 4)


def test_example4_open_run_at_record_end_counts() -> None:
    # + frame 0: ATG [0,3) under both; + frame 1: TGA (a run under S4 only); revcomp TCAT gives
    # [1,4) and [0,3) on the minus strand. Hand-computed: S4 (4 bp, 4 runs), S11 (4 bp, 3 runs).
    assert frame0_plus("ATGA", S4, 1) == [(0, 3)]
    assert frame0_plus("ATGA", S11, 1) == [(0, 3)]
    assert cover("ATGA", S4, 1) == (4, 4)
    assert cover("ATGA", S11, 1) == (4, 3)


@pytest.mark.parametrize("stops", [S11, S4])
@pytest.mark.parametrize("k", [1, 300])
def test_example5_empty_record(stops: frozenset[str], k: int) -> None:
    assert cover("", stops, k) == (0, 0)


# ---------------------------------------------------------------- edge cases


def test_poly_a_has_six_open_runs() -> None:
    seq = "A" * 10_000
    assert cover(seq, S4, 300) == (10_000, 6)
    assert cover(seq, S11, 300) == (10_000, 6)


@pytest.mark.parametrize("seq", ["N" * 3000, "RYKM" * 900, "-" * 999], ids=["N", "iupac", "dash"])
def test_only_non_acgt_gives_nothing(seq: str) -> None:
    assert cover(seq, S4, 1) == (0, 0)
    assert cover(seq, S11, 1) == (0, 0)


def test_lowercase_equals_uppercase() -> None:
    seq = EX1 + "gccgcc" + EX3
    for stops in (S11, S4):
        for k in (1, 4, 10):
            assert cover(seq.lower(), stops, k) == cover(seq, stops, k)


def test_mixed_case_stop_is_a_stop() -> None:
    seq = "GCC" * 10 + "tGa" + "GCC" * 10
    assert cover(seq, S11, 11) == reference_cover(seq, S11, 11)
    assert cover(seq, S4, 21) == reference_cover(seq, S4, 21)


def test_min_codons_boundary_is_inclusive() -> None:
    seq = "GCC" * 300  # + frame 0: exactly 300 codons, open at both ends
    assert frame0_plus(seq, S4, 300) == [(0, 900)]
    assert cover(seq, S4, 300) == reference_cover(seq, S4, 300)
    assert cover(seq, S4, 300)[0] >= 900
    assert cover(seq, S4, 301) == reference_cover(seq, S4, 301)


# ---------------------------------------------------------------- AC 2: Hypothesis


SEQ = st.text(alphabet="ACGTN", min_size=0, max_size=600)
# stop-poor blocks so long runs (and overlaps across frames and strands) actually occur
BLOCKS = st.lists(
    st.sampled_from(["GCC", "GCCGCCGCC", "TGA", "TAA", "TAG", "A", "C", "N", "ATG", "ttg"]),
    max_size=200,
).map("".join)


@settings(max_examples=300, deadline=None)
@given(seq=SEQ, k=st.integers(1, 5), use_s4=st.booleans())
def test_agrees_with_brute_force_reference(seq: str, k: int, use_s4: bool) -> None:
    stops = S4 if use_s4 else S11
    assert cover(seq, stops, k) == reference_cover(seq, stops, k)


@settings(max_examples=200, deadline=None)
@given(seq=BLOCKS, k=st.integers(1, 5), use_s4=st.booleans())
def test_agrees_with_reference_on_structured_sequences(seq: str, k: int, use_s4: bool) -> None:
    stops = S4 if use_s4 else S11
    assert cover(seq, stops, k) == reference_cover(seq, stops, k)


@settings(max_examples=200, deadline=None)
@given(seq=st.one_of(SEQ, BLOCKS), k=st.integers(1, 5))
def test_s11_cover_inside_s4_cover(seq: str, k: int) -> None:
    """Every S11 run lies inside an S4 run, so bp11 <= bp4 (spec "Signal")."""
    bp11, _ = cover(seq, S11, k)
    bp4, _ = cover(seq, S4, k)
    assert 0 <= bp11 <= bp4 <= len(seq)


@settings(max_examples=200, deadline=None)
@given(seq=st.one_of(SEQ, BLOCKS), k=st.integers(1, 5), use_s4=st.booleans())
def test_reverse_complement_invariant(seq: str, k: int, use_s4: bool) -> None:
    stops = S4 if use_s4 else S11
    assert cover(revcomp(seq), stops, k) == cover(seq, stops, k)
