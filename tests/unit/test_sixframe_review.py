"""Review follow-ups for the six-frame index (issue M3-02b; docs/specs/M3-02-sixframe-index.md).

NOMATCH decodes to ``#``; ``to_interval`` rejects intervals longer than the replicon; lowercase and
``i`` queries fold like upper case; non-ASCII bytes keep coordinates; ``dedupe_loci`` reduces
unwrapped starts >= L; ``decode`` clamps out-of-range codes.
"""

from __future__ import annotations

import numpy as np
import pytest

from gippa.mapping.sixframe import (
    ALPHABET,
    AMBIG,
    NOMATCH,
    SixFrameIndex,
    build_sixframe,
    decode,
    dedupe_loci,
    encode_query,
    encode_target,
    translate11,
)
from gippa.model import Interval, Replicon

G = "ATGGCTAAAATTCTGTGGTAACGTNNA"  # 27 nt, linear
R = "AAAGAACTGTAAGGCTCCATGCCG"  # 24 nt, circular
P0, P1, P2, M0, M1, M2 = range(6)


def _circular(seq: str = R, wrap_nt: int = 12) -> SixFrameIndex:
    return build_sixframe({"r": seq}, {"r": Replicon("r", len(seq), "circular")}, wrap_nt=wrap_nt)


# --- NOMATCH decodes to '#' -----------------------------------------------------------------


def test_nomatch_decodes_to_hash_not_x() -> None:
    q = encode_query("SEUK")
    assert q[2] == NOMATCH
    s = decode(q)
    assert s == "SE#K"
    assert "X" not in s


@pytest.mark.parametrize("letter", list("BJOUXZ*"))
def test_nomatch_letters_decode_to_hash(letter: str) -> None:
    assert decode(encode_query(f"PE{letter}K")) == "PE#K"


def test_nomatch_never_equals_any_target_code() -> None:
    target = encode_target(ALPHABET)
    assert NOMATCH not in target.tolist()
    assert "#" not in ALPHABET
    assert decode(np.array([AMBIG], dtype=np.uint8)) == "X"


# --- decode clamp ---------------------------------------------------------------------------


@pytest.mark.parametrize("code", [NOMATCH, NOMATCH + 1, 100, 255])
def test_decode_clamps_out_of_range_codes_to_hash(code: int) -> None:
    out = decode(np.array([0, code, 0], dtype=np.uint8))
    assert out == "A#A"


def test_decode_preserves_length() -> None:
    codes = np.arange(256, dtype=np.uint8)
    out = decode(codes)
    assert len(out) == 256
    assert out[: len(ALPHABET)] == ALPHABET
    assert set(out[len(ALPHABET) :]) == {"#"}


# --- to_interval length guard ---------------------------------------------------------------


@pytest.mark.parametrize("fid", [P0, M0])
def test_to_interval_rejects_interval_longer_than_replicon(fid: int) -> None:
    idx = _circular()  # L = 24, frames of 12 aa over 36 nt
    assert len(idx.frames[fid].aa) == 12
    with pytest.raises(ValueError):
        idx.to_interval(fid, 0, 9)  # 27 nt > 24
    with pytest.raises(ValueError):
        idx.to_interval(fid, 0, 12)  # 36 nt > 24


def test_to_interval_accepts_interval_equal_to_replicon() -> None:
    idx = _circular()
    iv = idx.to_interval(P0, 0, 8)  # exactly 24 nt
    assert iv.length_bp == 24
    assert iv == Interval("r", 0, 24, "+")
    iv2 = idx.to_interval(P0, 1, 9)
    assert iv2 == Interval("r", 3, 27, "+")


# --- query case folding ---------------------------------------------------------------------


def test_lowercase_query_equals_uppercase() -> None:
    assert np.array_equal(encode_query("peptidek"), encode_query("PEPTIDEK"))
    assert np.array_equal(encode_query("PePtIdEk"), encode_query("PEPTLDEK"))


def test_lowercase_query_letters_outside_alphabet_are_nomatch() -> None:
    q = encode_query("seuk")
    assert q[2] == NOMATCH
    assert decode(q) == "SE#K"


@pytest.mark.parametrize("enc", [encode_query, encode_target])
def test_lowercase_i_folds_to_l(enc) -> None:  # type: ignore[no-untyped-def]
    assert enc("i").tolist() == enc("L").tolist() == enc("l").tolist() == enc("I").tolist()
    assert decode(enc("mike")) == "MLKE"


# --- non-ASCII input keeps coordinates ------------------------------------------------------


def test_non_ascii_query_keeps_length() -> None:
    q = encode_query("PEé K")
    assert len(q) == len("PEé K")
    assert q[2] == NOMATCH and q[3] == NOMATCH
    assert decode(q) == "PE##K"


def test_non_ascii_target_is_ambig_and_keeps_length() -> None:
    t = encode_target("MéK")
    assert t.tolist() == [encode_target("M")[0], AMBIG, encode_target("K")[0]]


def test_non_ascii_genome_byte_keeps_coordinates() -> None:
    seq = G[:3] + "é" + G[4:]  # replace one nt (code point) with a non-ASCII character
    assert len(seq) == len(G)
    idx = build_sixframe({"c": seq}, {"c": Replicon("c", len(seq), "linear")})
    ref = build_sixframe({"c": G}, {"c": Replicon("c", len(G), "linear")})
    for f, g in zip(idx.frames, ref.frames, strict=True):
        assert f.ext_nt == g.ext_nt == len(G)
        assert len(f.aa) == len(g.aa)
    # codon 1 of +0 (nt 3..5) is now ambiguous; every other codon is unchanged
    plus = decode(idx.frames[P0].aa)
    assert plus[1] == "X"
    assert plus[0] + plus[2:] == decode(ref.frames[P0].aa)[0] + decode(ref.frames[P0].aa)[2:]
    assert idx.to_interval(P0, 2, 5) == ref.to_interval(P0, 2, 5) == Interval("c", 6, 15, "+")
    assert decode(translate11(b"GC?")) == "X"


# --- dedupe_loci with unwrapped starts >= L -------------------------------------------------


def test_dedupe_reduces_unwrapped_start_modulo_length() -> None:
    reps = {"r": Replicon("r", 24, "circular")}
    ivs = [Interval("r", 24, 33, "+"), Interval("r", 0, 9, "+"), Interval("r", 30, 39, "-")]
    assert dedupe_loci(ivs, reps) == [Interval("r", 0, 9, "+"), Interval("r", 6, 15, "-")]


def test_dedupe_unwrapped_origin_crossing_kept_unwrapped() -> None:
    reps = {"r": Replicon("r", 24, "circular")}
    out = dedupe_loci([Interval("r", 42, 51, "+"), Interval("r", 18, 27, "+")], reps)
    assert out == [Interval("r", 18, 27, "+")]
    assert out[0].crosses_origin(reps["r"])
