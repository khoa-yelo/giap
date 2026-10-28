"""Decoy translations and the ``sixframe_index`` work file (issue M3-03; GIPPA_SPEC.md §7.2.1, §8.2,
§8.5 attack 5, §6.2, §6.6).

Conventions pinned here (from the issue text and the lead's defaults):

- ``*`` (``STOP``) and ``X`` (``AMBIG``) stay in place under ``shuffle``; windows of ``window``
  residues are counted over the *non-X* residues of each stop-free segment, from its start.
- ``reverse`` reverses the whole frame, so ``*`` positions are mirrored (``i -> n - 1 - i``).
- ``build_decoys`` processes frames in index order, drawing from one generator; the step seeds it
  with ``gippa.rng.for_step("sixframe_index", config.seed)``.
- Decoy frames keep ``contig``/``strand``/``offset``/``ext_nt`` of their target frame (decoy hits
  are mapped back through ``mirror_aa`` to a target position, §7.6.4).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from gippa.io.fasta import read_genome_fasta
from gippa.mapping.decoys import build_decoys, mirror_aa, reverse_frame, shuffle_frame
from gippa.mapping.sixframe import (
    AMBIG,
    STOP,
    Frame,
    SixFrameIndex,
    build_sixframe,
    decode,
    encode_target,
    load_index,
    save_index,
)
from gippa.model import Replicon
from gippa.rng import for_step

STEP = "sixframe_index"


def _rng(seed: int = 7) -> np.random.Generator:
    return for_step(STEP, seed)


def _aa(s: str) -> np.ndarray:
    return encode_target(s)


def _random_frame(n: int, seed: int, p_stop: float = 0.0, p_x: float = 0.0) -> np.ndarray:
    gen = np.random.default_rng(seed)
    aa = gen.integers(0, STOP, size=n).astype(np.uint8)
    u = gen.random(n)
    aa[u < p_stop] = STOP
    aa[(u >= p_stop) & (u < p_stop + p_x)] = AMBIG
    return aa


def _ms(a: np.ndarray) -> Counter[int]:
    return Counter(int(v) for v in a)


def _assert_frames_equal(a: tuple[Frame, ...], b: tuple[Frame, ...]) -> None:
    assert len(a) == len(b)
    for fa, fb in zip(a, b, strict=True):
        assert (fa.contig, fa.strand, fa.offset, fa.ext_nt) == (
            fb.contig,
            fb.strand,
            fb.offset,
            fb.ext_nt,
        )
        assert fa.aa.dtype == np.uint8
        assert fb.aa.dtype == np.uint8
        np.testing.assert_array_equal(fa.aa, fb.aa)


# --- example 1: reverse ----------------------------------------------------------------------


def test_ex1_reverse_frame() -> None:
    out = reverse_frame(_aa("MKV*AL"))
    assert decode(out) == "LA*VKM"
    assert out.dtype == np.uint8


def test_reverse_frame_does_not_mutate_input() -> None:
    aa = _aa("MKV*AL")
    before = aa.copy()
    reverse_frame(aa)
    np.testing.assert_array_equal(aa, before)


def test_reverse_frame_keeps_x_and_folds_i() -> None:
    # I is folded into L by the target encoding; decoys work on codes and never reintroduce I.
    assert decode(reverse_frame(_aa("MIXK*"))) == "*KXLM"


def test_reverse_frame_empty() -> None:
    out = reverse_frame(np.zeros(0, dtype=np.uint8))
    assert out.shape == (0,)
    assert out.dtype == np.uint8


# --- example 2: stops fixed, segments shuffled independently --------------------------------


def test_ex2_shuffle_keeps_stops_and_segments() -> None:
    aa = _aa("MKVAL*GG*")
    out = shuffle_frame(aa, 50, _rng())
    s = decode(out)
    assert len(s) == 9
    assert s[5] == "*"
    assert s[8] == "*"
    assert sorted(s[:5]) == sorted("MKVAL")
    assert s[6:8] == "GG"
    assert out.dtype == np.uint8


def test_shuffle_frame_does_not_mutate_input() -> None:
    aa = _random_frame(300, 1, p_stop=0.05, p_x=0.02)
    before = aa.copy()
    shuffle_frame(aa, 50, _rng())
    np.testing.assert_array_equal(aa, before)


# --- example 3: 50-aa windows --------------------------------------------------------------


def test_ex3_windows_keep_their_multisets() -> None:
    aa = _random_frame(120, 3)
    out = shuffle_frame(aa, 50, _rng())
    for lo, hi in ((0, 50), (50, 100), (100, 120)):
        assert _ms(out[lo:hi]) == _ms(aa[lo:hi]), (lo, hi)
    assert not np.array_equal(out, aa)  # 120 random residues: identity has negligible probability


def test_windows_restart_at_each_segment() -> None:
    # window 3: segment "ACDE" -> windows [0,3) and [3,4); segment "FGHK" after '*' -> [5,8), [8,9).
    aa = _aa("ACDE*FGHK")
    for seed in range(20):
        s = decode(shuffle_frame(aa, 3, _rng(seed)))
        assert sorted(s[:3]) == sorted("ACD")
        assert s[3] == "E"
        assert s[4] == "*"
        assert sorted(s[5:8]) == sorted("FGH")
        assert s[8] == "K"


def test_windows_count_non_x_residues() -> None:
    # X stays put and takes no window slot: non-X indices [0,1,3,4] -> windows {0,1}, {3,4}.
    aa = _aa("ACXDE")
    seen0 = set()
    for seed in range(30):
        s = decode(shuffle_frame(aa, 2, _rng(seed)))
        assert s[2] == "X"
        assert sorted(s[0:2]) == ["A", "C"]
        assert sorted(s[3:5]) == ["D", "E"]
        seen0.add(s[0])
    assert seen0 == {"A", "C"}  # the window is actually permuted


# --- example 4: X fixed ----------------------------------------------------------------------


def test_ex4_x_stays_in_place() -> None:
    aa = _aa("AXKL")
    seen = set()
    for seed in range(40):
        s = decode(shuffle_frame(aa, 50, _rng(seed)))
        assert s[1] == "X"
        assert sorted(s[0] + s[2] + s[3]) == ["A", "K", "L"]
        seen.add(s)
    assert len(seen) > 1  # {A,K,L} are really permuted over indices 0, 2, 3


# --- example 5: determinism ----------------------------------------------------------------


def test_ex5_same_seed_byte_identical_other_seed_differs() -> None:
    aa = _random_frame(1000, 5, p_stop=0.02, p_x=0.01)
    a = shuffle_frame(aa, 50, for_step(STEP, 42))
    b = shuffle_frame(aa, 50, for_step(STEP, 42))
    c = shuffle_frame(aa, 50, for_step(STEP, 43))
    assert a.tobytes() == b.tobytes()
    assert a.tobytes() != c.tobytes()


def test_build_decoys_deterministic_from_config_seed() -> None:
    idx = _circular_index()
    a = build_decoys(idx, "shuffle", 50, for_step(STEP, 2026))
    b = build_decoys(idx, "shuffle", 50, for_step(STEP, 2026))
    c = build_decoys(idx, "shuffle", 50, for_step(STEP, 2027))
    _assert_frames_equal(a, b)
    assert any(not np.array_equal(x.aa, y.aa) for x, y in zip(a, c, strict=True))


# --- examples 6, 7: mirror_aa ----------------------------------------------------------------


def test_ex6_mirror_reverse() -> None:
    assert mirror_aa("reverse", 100, 10, 18) == (82, 90)


def test_ex7_mirror_shuffle_identity() -> None:
    assert mirror_aa("shuffle", 100, 10, 18) == (10, 18)


def test_mirror_reverse_whole_frame_and_edges() -> None:
    assert mirror_aa("reverse", 100, 0, 100) == (0, 100)
    assert mirror_aa("reverse", 100, 0, 1) == (99, 100)
    assert mirror_aa("reverse", 100, 99, 100) == (0, 1)


def test_mirror_reverse_window_is_reversed_target_window() -> None:
    aa = _aa("MKVALGG*TRWYPPHE")
    rev = reverse_frame(aa)
    n = len(aa)
    for a in range(n):
        for b in range(a + 1, n + 1):
            ms, me = mirror_aa("reverse", n, a, b)
            assert decode(rev[a:b])[::-1] == decode(aa[ms:me])


# --- build_decoys -------------------------------------------------------------------------


G = "ATGGCTAAAATTCTGTGGTAACGTNNA"  # 27 nt, linear contig c1
R = "AAAGAACTGTAAGGCTCCATGCCGTTTACCAGGATCATGAAACGT"  # 45 nt, circular replicon r


def _circular_index() -> SixFrameIndex:
    return build_sixframe(
        {"c1": G, "r": R},
        {"c1": Replicon("c1", len(G), "linear"), "r": Replicon("r", len(R), "circular", "plasmid")},
        wrap_nt=12,
    )


@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_build_decoys_mirrors_index_frames(kind: str) -> None:
    idx = _circular_index()
    dec = build_decoys(idx, kind, 50, _rng())  # type: ignore[arg-type]
    assert isinstance(dec, tuple)
    assert len(dec) == len(idx.frames) == 12
    for t, d in zip(idx.frames, dec, strict=True):
        assert isinstance(d, Frame)
        assert (d.contig, d.strand, d.offset, d.ext_nt) == (t.contig, t.strand, t.offset, t.ext_nt)
        assert d.aa.dtype == np.uint8
        assert len(d.aa) == len(t.aa)  # circular extension decoyed with the rest of the frame
        assert _ms(d.aa) == _ms(t.aa)


def test_build_decoys_reverse_is_reverse_frame_on_every_frame() -> None:
    idx = _circular_index()
    dec = build_decoys(idx, "reverse", 50, _rng())
    for t, d in zip(idx.frames, dec, strict=True):
        np.testing.assert_array_equal(d.aa, reverse_frame(t.aa))
    # minus-strand frames are reversed as amino-acid strings too
    assert any(f.strand == "-" for f in dec)


def test_build_decoys_reverse_ignores_rng() -> None:
    idx = _circular_index()
    _assert_frames_equal(
        build_decoys(idx, "reverse", 50, _rng(1)), build_decoys(idx, "reverse", 50, _rng(2))
    )


def test_build_decoys_shuffle_processes_frames_in_index_order() -> None:
    idx = _circular_index()
    dec = build_decoys(idx, "shuffle", 5, _rng(11))
    rng = _rng(11)
    expected = [shuffle_frame(f.aa, 5, rng) for f in idx.frames]
    for d, e in zip(dec, expected, strict=True):
        np.testing.assert_array_equal(d.aa, e)


def test_build_decoys_does_not_mutate_index() -> None:
    idx = _circular_index()
    before = [f.aa.copy() for f in idx.frames]
    build_decoys(idx, "shuffle", 5, _rng())
    build_decoys(idx, "reverse", 5, _rng())
    for f, b in zip(idx.frames, before, strict=True):
        np.testing.assert_array_equal(f.aa, b)


# --- edge cases ----------------------------------------------------------------------------


def test_frame_without_stop_is_one_segment() -> None:
    aa = _aa("MKVALGGTRWYPPHE")  # 15 aa, no '*'; window 50 -> one window over the whole frame
    out = shuffle_frame(aa, 50, _rng())
    assert _ms(out) == _ms(aa)
    assert STOP not in out


def test_all_stops_and_x_unchanged() -> None:
    aa = _aa("**X*XX*")
    np.testing.assert_array_equal(shuffle_frame(aa, 50, _rng()), aa)


def test_shuffle_empty_frame() -> None:
    out = shuffle_frame(np.zeros(0, dtype=np.uint8), 50, _rng())
    assert out.shape == (0,)
    assert out.dtype == np.uint8


def test_build_decoys_with_empty_frames() -> None:
    # contig of 4 nt: offset 2 leaves 2 nt -> zero-length frames
    idx = build_sixframe({"t": "ATGA"}, {"t": Replicon("t", 4, "linear")})
    assert any(len(f.aa) == 0 for f in idx.frames)
    for kind in ("reverse", "shuffle"):
        dec = build_decoys(idx, kind, 50, _rng())  # type: ignore[arg-type]
        assert [len(f.aa) for f in dec] == [len(f.aa) for f in idx.frames]


def test_shuffle_window_one_is_identity() -> None:
    aa = _random_frame(200, 9, p_stop=0.05, p_x=0.05)
    np.testing.assert_array_equal(shuffle_frame(aa, 1, _rng()), aa)


def test_shuffle_is_uniform_over_window_permutations() -> None:
    # Calibration of the shuffle (§8.2 "decoys mimic targets"): every arrangement of a 3-residue
    # window is equally likely. 6 arrangements x 600 seeds -> expected 100 each; a chi-square
    # with 5 df above 30 has p < 2e-5, so this fixed-seed check is not flaky by construction.
    aa = _aa("ACD")
    counts = Counter(decode(shuffle_frame(aa, 50, for_step(STEP, s))) for s in range(600))
    assert set(counts) == {"ACD", "ADC", "CAD", "CDA", "DAC", "DCA"}
    chi2 = sum((c - 100) ** 2 / 100 for c in counts.values())
    assert chi2 < 30, counts


# --- save_index / load_index (example 8) --------------------------------------------------


def _sim_index(tmp_path: Path) -> SixFrameIndex:
    import gippasim

    res = gippasim.simulate(gippasim.SimConfig(genome_length=50_000), seed=1, outdir=tmp_path / "s")
    recs = read_genome_fasta(res.genome_fasta)
    genome = {r.id: r.sequence for r in recs}
    reps = {r.id: Replicon(r.id, len(r.sequence), "circular", "chromosome") for r in recs}
    return build_sixframe(genome, reps)


def _assert_index_equal(a: SixFrameIndex, b: SixFrameIndex) -> None:
    assert dict(a.replicons) == dict(b.replicons)
    assert list(a.replicons) == list(b.replicons)
    _assert_frames_equal(a.frames, b.frames)


def test_ex8_roundtrip_50kb_simulated_genome(tmp_path: Path) -> None:
    idx = _sim_index(tmp_path)
    decoys = {
        "reverse": build_decoys(idx, "reverse", 50, _rng()),
        "shuffle": build_decoys(idx, "shuffle", 50, _rng()),
    }
    path = tmp_path / "sixframe_index.npz"
    save_index(path, idx, decoys)  # type: ignore[arg-type]
    assert path.is_file()
    idx2, dec2 = load_index(path)
    _assert_index_equal(idx, idx2)
    assert set(dec2) == {"reverse", "shuffle"}
    for kind in ("reverse", "shuffle"):
        _assert_frames_equal(decoys[kind], dec2[kind])  # type: ignore[index]
    # the loaded index still maps coordinates like the original
    assert idx2.to_interval(0, 3, 10) == idx.to_interval(0, 3, 10)


@pytest.mark.parametrize("kinds", [(), ("reverse",), ("shuffle",)])
def test_roundtrip_subset_of_decoys(tmp_path: Path, kinds: tuple[str, ...]) -> None:
    idx = _circular_index()
    decoys = {k: build_decoys(idx, k, 50, _rng()) for k in kinds}  # type: ignore[arg-type]
    path = tmp_path / "idx.npz"
    save_index(path, idx, decoys)  # type: ignore[arg-type]
    idx2, dec2 = load_index(path)
    _assert_index_equal(idx, idx2)
    assert set(dec2) == set(kinds)
    for k in kinds:
        _assert_frames_equal(decoys[k], dec2[k])  # type: ignore[index]


def test_roundtrip_keeps_empty_frames_and_contig_order(tmp_path: Path) -> None:
    genome = {"z": "ATGAAATAGCCC", "t": "AC", "a": G}
    reps = {
        "z": Replicon("z", 12, "circular", "plasmid"),
        "t": Replicon("t", 2, "linear"),
        "a": Replicon("a", len(G), "linear", "chromosome"),
    }
    idx = build_sixframe(genome, reps, wrap_nt=6)
    assert any(len(f.aa) == 0 for f in idx.frames)
    decoys = {
        "reverse": build_decoys(idx, "reverse", 50, _rng()),
        "shuffle": build_decoys(idx, "shuffle", 50, _rng()),
    }
    path = tmp_path / "idx.npz"
    save_index(path, idx, decoys)  # type: ignore[arg-type]
    idx2, dec2 = load_index(path)
    _assert_index_equal(idx, idx2)
    assert [f.contig for f in idx2.frames] == [f.contig for f in idx.frames]
    for k in ("reverse", "shuffle"):
        _assert_frames_equal(decoys[k], dec2[k])  # type: ignore[index]
