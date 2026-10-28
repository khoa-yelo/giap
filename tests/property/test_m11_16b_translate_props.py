"""AC 6 of docs/specs/M11-16-translation-tables-4-25.md: ``translate_text``, ``translate`` and
``build_sixframe`` agree with a brute-force per-codon reference of each table, on Hypothesis
sequences of 0-600 nt over ``ACGTN``, both strands, linear and circular (wrap covered).

The reference reads the NCBI TCAG strings of spec item 1 (equal to ``codon_table(t)``, AC 5); a
codon with a non-ACGT base is ``X``; in index frames I is folded into L (§7.2.1).
"""

from __future__ import annotations

import itertools

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.mapping.sixframe import build_sixframe, decode, translate, translate_text
from gippa.model import Replicon

TCAG = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
CODE = {
    t: dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)), s, strict=True))
    for t, s in TCAG.items()
}
_RC = str.maketrans("ACGTN", "TGCAN")

seqs = st.text(alphabet="ACGTN", min_size=0, max_size=600)
replicon_seqs = st.text(alphabet="ACGTN", min_size=1, max_size=600)  # Replicon refuses length 0
tables = st.sampled_from((4, 11, 25))


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


def _ref(nt: str, table: int, offset: int = 0) -> str:
    return "".join(CODE[table].get(nt[i : i + 3], "X") for i in range(offset, len(nt) - 2, 3))


@settings(max_examples=150, deadline=None)
@given(seqs, tables, st.integers(0, 2))
def test_translate_text_equals_the_reference(nt: str, table: int, offset: int) -> None:
    assert translate_text(nt, table, offset) == _ref(nt, table, offset)


@settings(max_examples=100, deadline=None)
@given(seqs, tables)
def test_translate_codes_equal_the_reference(nt: str, table: int) -> None:
    assert decode(translate(nt.encode(), table)) == _ref(nt, table).replace("I", "L")


@settings(max_examples=100, deadline=None)
@given(replicon_seqs, tables, st.booleans(), st.integers(0, 60))
def test_build_sixframe_equals_the_reference(
    nt: str, table: int, circular: bool, wrap_nt: int
) -> None:
    rep = Replicon("c", len(nt), "circular" if circular else "linear")
    idx = build_sixframe({"c": nt}, {"c": rep}, wrap_nt=wrap_nt, table=table)
    assert idx.translation_table == table
    ext = nt + nt[: min(wrap_nt, len(nt))] if circular else nt
    want = [_ref(ext, table, o) for o in range(3)] + [_ref(_rc(ext), table, o) for o in range(3)]
    got = [decode(f.aa) for f in idx.frames]
    assert got == [w.replace("I", "L") for w in want]
    assert [(f.strand, f.offset, f.ext_nt) for f in idx.frames] == [
        (s, o, len(ext)) for s in "+-" for o in range(3)
    ]


@settings(max_examples=100, deadline=None)
@given(replicon_seqs, st.booleans())
def test_tables_differ_only_at_tga(nt: str, circular: bool) -> None:
    """The three tables give the same frames except where table 11 reads ``*`` over a TGA."""
    rep = Replicon("c", len(nt), "circular" if circular else "linear")
    frames = {
        t: [decode(f.aa) for f in build_sixframe({"c": nt}, {"c": rep}, table=t).frames]
        for t in (4, 11, 25)
    }
    for f11, f4, f25 in zip(frames[11], frames[4], frames[25], strict=True):
        for a, b, c in zip(f11, f4, f25, strict=True):
            assert (a, b, c) in {(a, a, a), ("*", "W", "G")}
