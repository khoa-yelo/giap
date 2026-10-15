"""Six-frame translation and coordinate maps (issue M3-02; GIAP_SPEC.md §7.2.1, §6.5, §7.2.9).

Frame order per contig: ``+0, +1, +2, -0, -1, -2`` (frame ids 0..5 for a single contig).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from Bio.Seq import Seq

from giap.io.fasta import read_genome_fasta, write_fasta
from giap.mapping.sixframe import (
    ALPHABET,
    AMBIG,
    NOMATCH,
    STOP,
    WRAP_NT,
    Frame,
    SixFrameIndex,
    build_sixframe,
    decode,
    dedupe_loci,
    encode_query,
    encode_target,
    translate11,
)
from giap.model import Interval, Replicon, to_one_based_closed
from giapsim.codons import CODE_11

G = "ATGGCTAAAATTCTGTGGTAACGTNNA"  # 27 nt, linear contig c1
R = "AAAGAACTGTAAGGCTCCATGCCG"  # 24 nt, circular replicon r
P0, P1, P2, M0, M1, M2 = range(6)
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _linear(seq: str = G, cid: str = "c1") -> SixFrameIndex:
    return build_sixframe({cid: seq}, {cid: Replicon(cid, len(seq), "linear")})


def _circular(seq: str = R, cid: str = "r", wrap_nt: int = 12) -> SixFrameIndex:
    return build_sixframe({cid: seq}, {cid: Replicon(cid, len(seq), "circular")}, wrap_nt=wrap_nt)


def _frames(idx: SixFrameIndex) -> list[str]:
    return [decode(f.aa) for f in idx.frames]


# --- constants and encodings ---------------------------------------------------------------


def test_alphabet_and_codes() -> None:
    assert ALPHABET == "ACDEFGHKLMNPQRSTVWY*X"
    assert "I" not in ALPHABET
    assert (STOP, AMBIG, NOMATCH) == (19, 20, 21)
    # M3-07 review: the circular wrap covers a 40-aa peptide with k = 2 (3·(40 + 2) = 126 nt).
    assert WRAP_NT >= 126
    assert ALPHABET[STOP] == "*" and ALPHABET[AMBIG] == "X"


def test_encode_target_roundtrip_and_i_to_l() -> None:
    codes = encode_target("ACDEFGHKLMNPQRSTVWY*X")
    assert codes.dtype == np.uint8
    assert codes.tolist() == list(range(21))
    assert decode(codes) == ALPHABET
    assert encode_target("I").tolist() == encode_target("L").tolist()
    assert decode(encode_target("MIKE*")) == "MLKE*"


def test_example_10_query_i_equals_l_and_nomatch() -> None:
    a, b = encode_query("PEPTIDEK"), encode_query("PEPTLDEK")
    assert a.dtype == np.uint8
    assert np.array_equal(a, b)
    q = encode_query("SEUK")
    assert q[2] == NOMATCH
    assert q.tolist()[:2] == encode_target("SE").tolist()


@pytest.mark.parametrize("letter", list("BJOUXZ*"))
def test_query_letters_outside_alphabet_never_match(letter: str) -> None:
    code = int(encode_query(f"PE{letter}K")[2])
    assert code == NOMATCH
    assert code not in (STOP, AMBIG)


def test_query_codes_never_equal_stop_or_ambig() -> None:
    q = encode_query("ACDEFGHIKLMNPQRSTVWYBJOUXZ*")
    assert STOP not in q.tolist() and AMBIG not in q.tolist()


# --- table 11 ------------------------------------------------------------------------------


@pytest.mark.parametrize("codon", sorted(CODE_11))
def test_all_64_codons_table11(codon: str) -> None:
    expected = CODE_11[codon].replace("I", "L")
    oracle = str(Seq(codon).translate(table=11)).replace("I", "L")
    assert expected == oracle
    codes = translate11(codon.encode())
    assert codes.dtype == np.uint8
    assert decode(codes) == expected


def test_stop_codons_are_stop_code() -> None:
    assert translate11(b"TAATAGTGA").tolist() == [STOP, STOP, STOP]
    assert translate11(b"ATGTAAATG").tolist() == [
        int(encode_target("M")[0]),
        STOP,
        int(encode_target("M")[0]),
    ]


@pytest.mark.parametrize("codon", [b"NAA", b"ANA", b"AAN", b"RTG", b"GTN", b"TAY", b"NNN"])
def test_ambiguous_codons_are_x(codon: bytes) -> None:
    # Biopython resolves GTN -> V and TAY -> Y; the spec maps every non-ACGT codon to X.
    assert translate11(codon).tolist() == [AMBIG]


def test_translate11_drops_trailing_partial_codon() -> None:
    assert decode(translate11(b"ATGGC")) == "M"
    assert len(translate11(b"AT")) == 0


# --- examples 1-6 (linear) -------------------------------------------------------------------


def test_example_1_plus0() -> None:
    assert decode(_linear().frames[P0].aa) == "MAKLLW*RX"


def test_example_2_plus1_plus2() -> None:
    f = _frames(_linear())
    assert f[P1] == "WLKFCGNX"
    assert f[P2] == "G*NSVVTX"


def test_example_3_minus_frames() -> None:
    f = _frames(_linear())
    assert f[M0:] == ["XTLPQNFSH", "XRYHRLLA", "XVTTEF*P"]


def test_frame_order_and_fields() -> None:
    idx = _linear()
    assert len(idx.frames) == 6
    assert all(isinstance(f, Frame) for f in idx.frames)
    assert [(f.contig, f.strand, f.offset) for f in idx.frames] == [
        ("c1", "+", 0),
        ("c1", "+", 1),
        ("c1", "+", 2),
        ("c1", "-", 0),
        ("c1", "-", 1),
        ("c1", "-", 2),
    ]
    assert all(f.aa.dtype == np.uint8 for f in idx.frames)
    assert all(f.ext_nt == 27 for f in idx.frames)


def test_example_4_plus_interval() -> None:
    idx = _linear()
    iv = idx.to_interval(P0, 1, 4)
    assert iv == Interval("c1", 3, 12, "+")
    assert decode(idx.frames[P0].aa[1:4]) == "AKL"
    assert idx.frame_label(iv) == 1


def test_example_5_minus0_interval() -> None:
    idx = _linear()
    iv = idx.to_interval(M0, 1, 2)
    assert decode(idx.frames[M0].aa[1:2]) == "T"
    assert iv == Interval("c1", 21, 24, "-")
    assert idx.frame_label(iv) == -1
    assert _rc(G[21:24]) == "ACG"  # T


def test_example_6_minus1_interval() -> None:
    idx = _linear()
    iv = idx.to_interval(M1, 2, 4)
    assert decode(idx.frames[M1].aa[2:4]) == "YH"
    assert iv == Interval("c1", 14, 20, "-")
    assert idx.frame_label(iv) == -2


@pytest.mark.parametrize("fid", range(6))
def test_linear_frame_label_equals_offset_plus_one(fid: int) -> None:
    idx = _linear()
    frame = idx.frames[fid]
    sign = 1 if frame.strand == "+" else -1
    for a in range(len(frame.aa)):
        iv = idx.to_interval(fid, a, a + 1)
        assert idx.frame_label(iv) == sign * (frame.offset + 1)


def test_plus_and_minus_nt_translate_back_linear() -> None:
    idx = _linear()
    for fid, frame in enumerate(idx.frames):
        for a in range(len(frame.aa)):
            iv = idx.to_interval(fid, a, a + 1)
            nt = G[iv.start : iv.end]
            if iv.strand == "-":
                nt = _rc(nt)
            assert decode(translate11(nt.encode())) == decode(frame.aa[a : a + 1])


# --- examples 7-9 (circular) -----------------------------------------------------------------


def test_example_7_origin_crossing_plus() -> None:
    idx = _circular()
    frame = idx.frames[P0]
    assert decode(frame.aa) == "KEL*GSMPKEL*"
    assert frame.ext_nt == 36
    assert decode(frame.aa[6:11]) == "MPKEL"
    iv = idx.to_interval(P0, 6, 11)
    assert iv == Interval("r", 18, 33, "+")
    assert iv.crosses_origin(idx.replicons["r"])
    assert to_one_based_closed(iv) == (19, 33)
    assert idx.frame_label(iv) == 1


def test_example_7_mirrored_minus_strand_origin_crossing() -> None:
    idx = _circular(_rc(R))
    frame = idx.frames[M0]
    assert decode(frame.aa) == "GSMPKEL*GSMP"
    assert decode(frame.aa[2:7]) == "MPKEL"
    iv = idx.to_interval(M0, 2, 7)
    assert iv == Interval("r", 15, 30, "-")
    assert iv.crosses_origin(idx.replicons["r"])
    assert idx.frame_label(iv) == -1


def test_example_8_dedupe_origin_copy() -> None:
    idx = _circular()
    first, second = idx.to_interval(P0, 0, 3), idx.to_interval(P0, 8, 11)
    assert decode(idx.frames[P0].aa[8:11]) == "KEL"
    assert first == Interval("r", 0, 9, "+")
    assert second == Interval("r", 0, 9, "+")
    assert dedupe_loci([first, second], idx.replicons) == [Interval("r", 0, 9, "+")]


def test_dedupe_keeps_distinct_strands_and_loci() -> None:
    idx = _circular()
    ivs = [Interval("r", 0, 9, "+"), Interval("r", 0, 9, "-"), Interval("r", 3, 12, "+")]
    assert sorted(dedupe_loci(ivs + ivs, idx.replicons), key=repr) == sorted(ivs, key=repr)


def test_example_9_linear_never_wraps() -> None:
    idx = _linear(R, "r")
    assert decode(idx.frames[P0].aa) == "KEL*GSMP"
    for fid, frame in enumerate(idx.frames):
        assert frame.ext_nt == 24
        n = len(frame.aa)
        for a in range(n):
            for b in range(a + 1, n + 1):
                assert idx.to_interval(fid, a, b).end <= 24


# --- edge cases ------------------------------------------------------------------------------


@pytest.mark.parametrize("length", [1, 2, 3, 4, 5, 26, 28])
def test_frame_lengths_floor(length: int) -> None:
    seq = (G * 2)[:length]
    idx = _linear(seq)
    for frame in idx.frames:
        assert len(frame.aa) == max(0, (length - frame.offset) // 3)


def test_circular_shorter_than_wrap_appends_whole_replicon() -> None:
    seq = "ATGAAACCCGGGT"  # 13 nt < 120
    idx = build_sixframe({"p": seq}, {"p": Replicon("p", 13, "circular")})
    for frame in idx.frames:
        assert frame.ext_nt == 26
        assert len(frame.aa) == (26 - frame.offset) // 3
    assert decode(idx.frames[P0].aa) == decode(translate11((seq + seq).encode()))


def test_circular_default_wrap() -> None:
    # Default wrap = WRAP_NT (≥ 126 nt since the M3-07 review; was 120).
    rng = np.random.default_rng(7)
    seq = "".join(rng.choice(list("ACGT"), size=500))
    idx = build_sixframe({"c": seq}, {"c": Replicon("c", 500, "circular")})
    assert all(f.ext_nt == 500 + WRAP_NT for f in idx.frames)
    assert decode(idx.frames[P0].aa) == decode(translate11((seq + seq[:WRAP_NT]).encode()))
    assert decode(idx.frames[M0].aa) == decode(translate11(_rc(seq + seq[:WRAP_NT]).encode()))


def test_lower_case_is_upper_cased() -> None:
    assert _frames(_linear(G.lower())) == _frames(_linear())


def test_multi_contig_order_and_frame_ids() -> None:
    genome = {"z": R, "a": G}
    reps = {"a": Replicon("a", 27, "linear"), "z": Replicon("z", 24, "circular")}
    idx = build_sixframe(genome, reps, wrap_nt=12)
    assert [f.contig for f in idx.frames] == ["z"] * 6 + ["a"] * 6
    assert decode(idx.frames[6].aa) == "MAKLLW*RX"
    assert idx.to_interval(6, 1, 4) == Interval("a", 3, 12, "+")
    assert idx.to_interval(6 + M1, 2, 4) == Interval("a", 14, 20, "-")
    assert idx.to_interval(P0, 6, 11) == Interval("z", 18, 33, "+")


def test_build_from_genome_fasta(tmp_path: Path) -> None:
    path = tmp_path / "g.fa"
    write_fasta(path, [("c1", G), ("r", R)])
    recs = read_genome_fasta(path)
    reps = {"c1": Replicon("c1", 27, "linear"), "r": Replicon("r", 24, "circular")}
    idx = build_sixframe({r.id: r.sequence for r in recs}, reps, wrap_nt=12)
    assert decode(idx.frames[P0].aa) == "MAKLLW*RX"
    assert decode(idx.frames[6 + P0].aa) == "KEL*GSMPKEL*"


def test_deterministic() -> None:
    a, b = _circular(), _circular()
    assert len(a.frames) == len(b.frames)
    for fa, fb in zip(a.frames, b.frames, strict=True):
        assert (fa.contig, fa.strand, fa.offset, fa.ext_nt) == (
            fb.contig,
            fb.strand,
            fb.offset,
            fb.ext_nt,
        )
        assert np.array_equal(fa.aa, fb.aa)
    ivs = [a.to_interval(P0, 0, 3), a.to_interval(P0, 8, 11), a.to_interval(M0, 0, 2)]
    assert dedupe_loci(ivs, a.replicons) == dedupe_loci(ivs, b.replicons)
