"""Hypothesis properties for decoy translations (issue M3-03; GIPPA_SPEC.md §7.2.9 "decoy
translations keep composition", §8.2, §7.6.4).

Frames are drawn as code arrays over ``ALPHABET`` with extra weight on ``*`` and ``X``. Windows
are counted over non-X residues of each stop-free segment (the issue's shuffle rule).
"""

from __future__ import annotations

from collections import Counter

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.mapping.decoys import build_decoys, mirror_aa, reverse_frame, shuffle_frame
from gippa.mapping.sixframe import AMBIG, STOP, build_sixframe
from gippa.model import Replicon
from gippa.rng import for_step

SETTINGS = settings(max_examples=60, deadline=None)

_code = st.one_of(
    st.integers(0, STOP - 1),
    st.integers(0, STOP - 1),
    st.integers(0, STOP - 1),
    st.just(STOP),
    st.just(AMBIG),
)
frames = st.lists(_code, min_size=0, max_size=400).map(lambda xs: np.array(xs, dtype=np.uint8))
windows = st.integers(1, 60)
seeds = st.integers(0, 2**32 - 1)


def _ms(a: np.ndarray) -> Counter[int]:
    return Counter(int(v) for v in a)


def _windows(aa: np.ndarray, window: int) -> list[list[int]]:
    """Index lists of every shuffle window (non-X positions, per stop-free segment)."""
    out: list[list[int]] = []
    seg: list[int] = []
    for i, c in enumerate([*aa.tolist(), STOP]):
        if c == STOP:
            out.extend(seg[k : k + window] for k in range(0, len(seg), window))
            seg = []
        elif c != AMBIG:
            seg.append(i)
    return out


@SETTINGS
@given(frames)
def test_reverse_keeps_composition_and_mirrors_stops(aa: np.ndarray) -> None:
    out = reverse_frame(aa)
    n = len(aa)
    assert out.dtype == np.uint8
    assert len(out) == n
    assert _ms(out) == _ms(aa)
    stops = set(np.flatnonzero(aa == STOP).tolist())
    assert set(np.flatnonzero(out == STOP).tolist()) == {n - 1 - i for i in stops}
    np.testing.assert_array_equal(reverse_frame(out), aa)  # reverse is an involution


@SETTINGS
@given(frames, windows, seeds)
def test_shuffle_keeps_composition_stops_x_and_window_multisets(
    aa: np.ndarray, window: int, seed: int
) -> None:
    out = shuffle_frame(aa, window, for_step("sixframe_index", seed))
    assert out.dtype == np.uint8
    assert len(out) == len(aa)
    assert _ms(out) == _ms(aa)
    fixed = (aa == STOP) | (aa == AMBIG)
    np.testing.assert_array_equal(out[fixed], aa[fixed])
    assert not np.any((out == STOP) & ~(aa == STOP))
    for idx in _windows(aa, window):
        assert _ms(out[idx]) == _ms(aa[idx]), idx


@SETTINGS
@given(frames, windows, seeds)
def test_shuffle_deterministic(aa: np.ndarray, window: int, seed: int) -> None:
    a = shuffle_frame(aa, window, for_step("sixframe_index", seed))
    b = shuffle_frame(aa, window, for_step("sixframe_index", seed))
    assert a.tobytes() == b.tobytes()


@SETTINGS
@given(
    st.integers(1, 10_000)
    .flatmap(lambda n: st.tuples(st.just(n), st.integers(0, n - 1)))
    .flatmap(lambda t: st.tuples(st.just(t[0]), st.just(t[1]), st.integers(t[1] + 1, t[0])))
)
def test_mirror_reverse_is_involution_inside_frame(t: tuple[int, int, int]) -> None:
    n, a, b = t
    ma, mb = mirror_aa("reverse", n, a, b)
    assert 0 <= ma < mb <= n
    assert mb - ma == b - a
    assert mirror_aa("reverse", n, ma, mb) == (a, b)
    assert mirror_aa("shuffle", n, a, b) == (a, b)


_seq = st.text(alphabet="ACGTN", min_size=1, max_size=600).filter(lambda s: s.count("N") < 30)


@settings(max_examples=25, deadline=None)
@given(_seq, st.sampled_from(["linear", "circular"]), windows, seeds)
def test_build_decoys_keep_composition_per_frame(
    seq: str, topology: str, window: int, seed: int
) -> None:
    idx = build_sixframe({"c": seq}, {"c": Replicon("c", len(seq), topology)}, wrap_nt=30)  # type: ignore[arg-type]
    rng = for_step("sixframe_index", seed)
    for kind in ("reverse", "shuffle"):
        dec = build_decoys(idx, kind, window, rng)  # type: ignore[arg-type]
        assert len(dec) == len(idx.frames)
        for t, d in zip(idx.frames, dec, strict=True):
            assert (d.contig, d.strand, d.offset, d.ext_nt) == (
                t.contig,
                t.strand,
                t.offset,
                t.ext_nt,
            )
            assert _ms(d.aa) == _ms(t.aa)
            if kind == "shuffle":
                np.testing.assert_array_equal(d.aa == STOP, t.aa == STOP)
                for w in _windows(t.aa, window):
                    assert _ms(d.aa[w]) == _ms(t.aa[w])
