"""Hypothesis properties for the six-frame index (issue M3-02; GIPPA_SPEC.md §7.2.1, §11.2).

Genomes are drawn as (seed, length, topology) and generated with a seeded RNG so that 5-kb
sequences stay cheap; about 1 % of positions are ``N`` to exercise ambiguous codons.
"""

from __future__ import annotations

import random

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.mapping.sixframe import (
    SixFrameIndex,
    build_sixframe,
    decode,
    dedupe_loci,
    encode_query,
    encode_target,
)
from gippa.model import Interval, Replicon, split_at_origin
from gippasim.codons import CODE_11

SETTINGS = settings(max_examples=40, deadline=None)
_RC = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _oracle(nt: str) -> str:
    """Table-11 translation (I -> L, non-ACGT codon -> X) via the simulator's codon table."""
    out = []
    for i in range(0, len(nt) - 2, 3):
        aa = CODE_11.get(nt[i : i + 3], "X")
        out.append("L" if aa == "I" else aa)
    return "".join(out)


def _genome(seed: int, length: int, n_frac: float = 0.01) -> str:
    rnd = random.Random(seed)
    return "".join("N" if rnd.random() < n_frac else rnd.choice("ACGT") for _ in range(length))


def _build(seq: str, topology: str, wrap_nt: int = 120) -> SixFrameIndex:
    return build_sixframe({"c": seq}, {"c": Replicon("c", len(seq), topology)}, wrap_nt=wrap_nt)


def _locus_nt(seq: str, rep: Replicon, iv: Interval) -> str:
    nt = "".join(seq[p.start : p.end] for p in split_at_origin(iv, rep))
    return _rc(nt) if iv.strand == "-" else nt


def _windows(idx: SixFrameIndex, k: int) -> dict[Interval, str]:
    """Every aa window of length ``k`` in every frame, mapped to its interval."""
    out: dict[Interval, str] = {}
    for fid, frame in enumerate(idx.frames):
        for a in range(len(frame.aa) - k + 1):
            out[idx.to_interval(fid, a, a + k)] = decode(frame.aa[a : a + k])
    return out


genomes = st.tuples(
    st.integers(0, 2**32 - 1),
    st.integers(1, 5000),
    st.sampled_from(["circular", "linear"]),
)


@SETTINGS
@given(genomes, st.integers(0, 2**32 - 1))
def test_interval_nucleotides_translate_to_frame_codes(g: tuple[int, int, str], wseed: int) -> None:
    seed, length, topology = g
    seq = _genome(seed, length)
    idx = _build(seq, topology)
    rep = idx.replicons["c"]
    rnd = random.Random(wseed)
    for fid, frame in enumerate(idx.frames):
        n = len(frame.aa)
        sign = 1 if frame.strand == "+" else -1
        assert frame.ext_nt == (length + min(120, length) if topology == "circular" else length)
        assert n == max(0, (frame.ext_nt - frame.offset) // 3)
        if n == 0:
            continue
        starts = [rnd.randrange(n) for _ in range(15)] + [n - 1, max(0, n - 40)]
        for a in starts:
            b = min(n, a + rnd.randint(1, 40))
            if 3 * (b - a) > length:
                continue
            iv = idx.to_interval(fid, a, b)
            assert iv.contig == "c" and iv.strand == frame.strand
            assert iv.end - iv.start == 3 * (b - a)
            assert 0 <= iv.start < length
            if topology == "linear":
                assert iv.end <= length
                assert idx.frame_label(iv) == sign * (frame.offset + 1)
            assert _oracle(_locus_nt(seq, rep, iv)) == decode(frame.aa[a:b])
            # frame label by its definition (TSV column `frame`)
            if iv.strand == "+":
                assert idx.frame_label(iv) == 1 + iv.start % 3
            else:
                assert idx.frame_label(iv) == -(1 + ((length - iv.end) % length) % 3)


@SETTINGS
@given(
    st.integers(0, 2**32 - 1),
    st.integers(3, 600),
    st.sampled_from(["circular", "linear"]),
    st.integers(1, 40),
)
def test_reverse_complement_swaps_the_strand(seed: int, length: int, topology: str, k: int) -> None:
    k = min(k, length // 3)
    seq = _genome(seed, length)
    fwd = _windows(_build(seq, topology), k)
    rev = _windows(_build(_rc(seq), topology), k)
    for iv, aa in fwd.items():
        if iv.strand != "+":
            continue
        start = (length - iv.end) % length
        mirrored = Interval("c", start, start + iv.end - iv.start, "-")
        assert rev.get(mirrored) == aa, (iv, mirrored)


@SETTINGS
@given(st.integers(0, 2**32 - 1), st.integers(3, 400), st.integers(1, 40))
def test_origin_windows_reachable_exactly_once_after_dedupe(seed: int, length: int, k: int) -> None:
    k = min(k, length // 3)
    w = 3 * k
    seq = _genome(seed, length, n_frac=0.0)
    idx = _build(seq, "circular")
    rep = idx.replicons["c"]
    all_ivs = [
        idx.to_interval(fid, a, a + k)
        for fid, frame in enumerate(idx.frames)
        for a in range(len(frame.aa) - k + 1)
    ]
    loci = dedupe_loci(all_ivs, idx.replicons)
    assert len(loci) == len(set(loci))
    assert set(loci) == set(all_ivs)  # dedupe only removes copies, never distinct loci
    assert all(0 <= iv.start < length for iv in loci)
    for strand in ("+", "-"):
        for s in range(max(0, length - w + 1), length):
            target = Interval("c", s, s + w, strand)
            assert target.crosses_origin(rep)
            assert loci.count(target) == 1, target
    assert dedupe_loci(all_ivs, idx.replicons) == loci  # deterministic


@SETTINGS
@given(
    st.integers(0, 2**32 - 1),
    st.integers(30, 2000),
    st.sampled_from(["circular", "linear"]),
    st.sampled_from(["+", "-"]),
    st.integers(0, 2**32 - 1),
    st.integers(5, 40),
)
def test_peptide_cut_from_genome_maps_back_to_origin(
    seed: int, length: int, topology: str, strand: str, cut: int, k: int
) -> None:
    k = min(k, length // 3)
    seq = _genome(seed, length, n_frac=0.0)
    rep = Replicon("c", length, topology)
    span = 3 * k
    start = cut % length if topology == "circular" else cut % (length - span + 1)
    origin = Interval("c", start, start + span, strand)
    peptide = _oracle(_locus_nt(seq, rep, origin))
    query = encode_query(peptide)
    assert np.array_equal(query, encode_target(peptide)) or "*" in peptide
    needle = encode_target(peptide).tobytes()
    idx = _build(seq, topology)
    hits = []
    for fid, frame in enumerate(idx.frames):
        hay = frame.aa.tobytes()
        pos = hay.find(needle)
        while pos != -1:
            hits.append(idx.to_interval(fid, pos, pos + k))
            pos = hay.find(needle, pos + 1)
    loci = dedupe_loci(hits, idx.replicons)
    assert loci.count(origin) == 1
