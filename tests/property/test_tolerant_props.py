"""Hypothesis properties of the banded tolerant verifier (issue M3-06; GIAP_SPEC.md §7.2.4, §11.2).

"Tolerant matcher = brute force": random (query, window) pairs are built by planting 0-2
operations of each type (substitution, swap, block, near_block; every residue order and both
orientations) between random matching chunks, optionally with a stop/X in the window, a query
NOMATCH and random flanks. The independent oracle (``oracle_tolerant.py``) is compared with
``verify``:

- the minimal cost and the reported diff types equal the banded oracle's (same tie-break);
- the reported window really has that cost and those types (oracle on that one window);
- zero false negatives: a clean plant of ``n <= k`` operations is always found with cost <= n;
- soundness: the unbanded, exhaustive oracle never finds more than ``verify`` reports missing
  (``verify`` never reports a cost below the true minimum);
- cost and diff types are symmetric under reversing both query and target.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

from hypothesis import assume, example, given, settings
from hypothesis import strategies as st
from oracle_tolerant import align_window, oracle_verify, window_cost

from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import decode, encode_query, encode_target
from giap.mapping.tolerant_dp import verify

TABLES = block_tables()
SETTINGS = settings(max_examples=300, deadline=None)
RESIDUES = "ACDEFGHKLMNPQRSTVWY"

# spec text, §7.2.4 (both orders of each pair)
_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _orders(s: str) -> list[str]:
    return sorted({"".join(p) for p in permutations(s)})


def _expand(raw: tuple[tuple[str, str], ...]) -> list[tuple[str, str]]:
    out = []
    for a, b in raw:
        for x in _orders(a):
            for y in _orders(b):
                out += [(x, y), (y, x)]
    return out


BLOCK_SEGS = _expand(_BLOCK_RAW)
NEAR_SEGS = _expand(_NEAR_RAW)


@dataclass(frozen=True)
class Case:
    query: str  # letters as given to encode_query (I allowed, B/J/X -> NOMATCH)
    frame: str  # letters as given to encode_target
    d0: int
    k: int
    planted: int  # number of planted operations
    clean: bool  # no stop/X/NOMATCH planted into the alignment


chunk = st.text(alphabet=RESIDUES + "I", min_size=0, max_size=4)


@st.composite
def op_piece(draw: st.DrawFn) -> tuple[str, str]:
    kind = draw(st.sampled_from(["substitution", "swap", "block", "near_block"]))
    if kind == "substitution":
        a = draw(st.sampled_from(RESIDUES))
        b = draw(st.sampled_from(RESIDUES.replace(a, "")))
        return a, b
    if kind == "swap":
        a = draw(st.sampled_from(RESIDUES))
        b = draw(st.sampled_from(RESIDUES.replace(a, "")))
        return a + b, b + a
    return draw(st.sampled_from(BLOCK_SEGS if kind == "block" else NEAR_SEGS))


@st.composite
def cases(draw: st.DrawFn) -> Case:
    n_ops = draw(st.integers(0, 4))
    first = draw(chunk)
    q, t = [first], [first]
    for _ in range(n_ops):
        qs, ts = draw(op_piece())
        c = draw(chunk)
        q += [qs, c]
        t += [ts, c]
    query, window = "".join(q), "".join(t)
    short = max(0, 8 - len(query))
    pad = draw(st.text(alphabet=RESIDUES, min_size=short, max_size=short + 3))
    query, window = query + pad, window + pad
    clean = True
    if window and draw(st.integers(0, 9)) == 0:  # a stop or X somewhere in the window
        pos = draw(st.integers(0, len(window) - 1))
        window = window[:pos] + draw(st.sampled_from("*X")) + window[pos + 1 :]
        clean = False
    if draw(st.integers(0, 9)) == 0:  # a query NOMATCH letter
        pos = draw(st.integers(0, len(query) - 1))
        query = query[:pos] + draw(st.sampled_from("BJXZ*U")) + query[pos + 1 :]
        clean = False
    flank = st.text(alphabet=RESIDUES + "*X", max_size=3)
    left, right = draw(flank), draw(flank)
    k = draw(st.integers(0, 2))
    return Case(query, left + window + right, len(left), k, n_ops, clean)


def _verify(c: Case) -> list:  # type: ignore[type-arg]
    return verify(encode_query(c.query), encode_target(c.frame), c.d0, c.k, TABLES)


def _letters(c: Case) -> tuple[str, str]:
    return decode(encode_query(c.query)), decode(encode_target(c.frame))


EXAMPLES = [
    Case("AKLMNPQRSTVW", "AKLMNPQRSTVF", 0, 1, 1, True),  # substitution
    Case("AKLMNPQRSTVW", "AKLMPNQRSTVW", 0, 1, 1, True),  # swap
    Case("PEPTNDEKAR", "PEPTGGDEKAR", 0, 1, 1, True),  # block 1 -> 2
    Case("AGLLKDEPTR", "QLLKDEPTR", 0, 1, 1, True),  # block 2 -> 1
    Case("VTEGTAPKLR", "VTEASAPKLR", 0, 1, 1, True),  # block 2 <-> 2
    Case("SWAPKELLTR", "SADAPKELLTR", 0, 1, 1, True),  # near_block 1 -> 2
    Case("MKGVLLAEQTR", "MKRLLAEQTR", 0, 1, 1, True),  # near_block 2 -> 1
    Case("AKNLMPQNSTVW", "AKGGLMPQGGSTVW", 0, 2, 2, True),  # two blocks, diagonal +2
    Case("MKWLHPCEGYHPK", "MKADLHPCWYHPK", 0, 2, 2, True),  # two near blocks, -/+
    Case("PEPTIDEKAR", "PEPTLDEKAR", 0, 0, 0, True),  # I = L
    Case("PEPTIDEKAR", "PEPT*DEKAR", 0, 2, 0, False),  # stop
    Case("PEPTBDEKAR", "PEPTXDEKAR", 0, 2, 0, False),  # NOMATCH vs X
]


def _with_examples(f):  # type: ignore[no-untyped-def]
    for ex in EXAMPLES:
        f = example(c=ex)(f)
    return f


@SETTINGS
@given(c=cases())
@_with_examples
def test_verify_equals_banded_brute_force(c: Case) -> None:
    q, f = _letters(c)
    got = _verify(c)
    want = oracle_verify(q, f, c.d0, c.k)
    if want is None:
        assert got == []
        return
    assert len(got) == 1, got
    aln = got[0]
    assert (aln.cost, aln.diff_types) == want[:2]
    assert aln.t_start == want[2]
    assert 0 <= aln.t_start < aln.t_end <= len(f)
    assert window_cost(q, f, aln.t_start, aln.t_end, c.d0, c.k) == (aln.cost, aln.diff_types)


@SETTINGS
@given(c=cases())
@_with_examples
def test_zero_false_negatives_for_planted_operations(c: Case) -> None:
    got = _verify(c)
    if c.clean and c.planted <= c.k:
        assert len(got) == 1, (c, got)
        assert got[0].cost <= c.planted


@SETTINGS
@given(c=cases())
@_with_examples
def test_never_below_the_unbanded_exhaustive_minimum(c: Case) -> None:
    q, f = _letters(c)
    got = _verify(c)
    full = oracle_verify(q, f, c.d0, c.k, band=False)
    if full is None:
        assert got == []
    elif got:
        assert got[0].cost >= full[0]
        # the reported window is a real alignment of that cost (no band needed to check it)
        best = align_window(q, f, got[0].t_start, got[0].t_end, None, c.k)
        assert best is not None and best[0] <= got[0].cost


@SETTINGS
@given(c=cases())
@_with_examples
def test_symmetric_under_reversal(c: Case) -> None:
    fwd = _verify(c)
    L = len(c.query)
    rev = Case(c.query[::-1], c.frame[::-1], len(c.frame) - L - c.d0, c.k, c.planted, c.clean)
    assume(rev.d0 >= -2)  # negative candidate diagonals beyond -2 are never seeded
    back = _verify(rev)
    assert [(a.cost, a.diff_types) for a in fwd] == [(a.cost, a.diff_types) for a in back]


@SETTINGS
@given(c=cases())
def test_deterministic(c: Case) -> None:
    assert _verify(c) == _verify(c)


# ---- the oracle itself agrees with the spec examples -----------------------------------------


def test_oracle_reproduces_spec_examples() -> None:
    rows = [
        ("AKLMNPQRSTVW", 2, "AKLMPNQRSTVW", (1, ("swap",), 0)),
        ("PEPTNDEKAR", 1, "PEPTGGDEKAR", (1, ("block",), 0)),
        ("VTEGTAPKLR", 1, "VTEASAPKLR", (1, ("block",), 0)),
        ("SWAPKELLTR", 1, "SADAPKELLTR", (1, ("near_block",), 0)),
        ("MKRLLAEQTR", 1, "MKGVLLAEQTR", (1, ("near_block",), 0)),
        ("PEPTLDEKAR", 1, "PEPTLDEKAR", (0, (), 0)),
        ("PEPTLDEKAR", 1, "PAPTLDEKGR", None),
        ("AKLMNPQRSTVW", 2, "AKLMPNQRSTVF", (2, ("substitution", "swap"), 0)),
        ("PEPTLDEKAR", 1, "PEPT*DEKAR", None),
        ("AGLLKDEPTR", 1, "QLLKDEPTR", (1, ("block",), 0)),
    ]
    for q, k, t, want in rows:
        assert oracle_verify(q, t, 0, k) == want, q
        assert oracle_verify(q, t, 0, k, band=False) == want, q
