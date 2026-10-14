"""Review pins for the banded tolerant verifier (issue M3-06, PR #39 code review).

docs/specs/M3-06-tolerant-verifier.md, "Lead decisions after code review", items 1-5:

1. ``verify`` raises ``ValueError`` for queries of ``L <= 2 * BAND`` (= 4 aa); ``L = 5`` is
   accepted.
2. A stop or ``X`` in each target slot, and a query NOMATCH (``#``) in each query slot, of every
   block shape (1<->1, 1<->2, 2<->1, 2<->2) is never consumed by a match, swap or (near-)block,
   including table-index aliasing cases such as ``P`` vs ``*G`` and ``GS`` vs ``*S``.
3. The kernel runs without numba: a subprocess with ``NUMBA_DISABLE_JIT=1`` reproduces the jitted
   results on a representative subset (examples, every block shape, substitutions, swaps).
4. ``verify`` normalises its inputs: int64 and strided (non-contiguous) arrays give the same result
   as C-contiguous uint8 ones (and do not compile a new kernel specialisation).
5. ``t_end`` tie-break: equal cost, equal diff types, equal ``t_start`` -> the smallest ``t_end``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest

from giap.mapping.diffops import BLOCKS, NEAR_BLOCKS, block_tables
from giap.mapping.sixframe import decode, encode_query, encode_target
from giap.mapping.tolerant_dp import BAND, Alignment, verify

TABLES = block_tables()
SRC = Path(__file__).resolve().parents[2] / "src"
RESIDUES = decode(np.arange(19, dtype=np.uint8))  # the 19 residue letters (I folded into L)
# flanks of residues that occur in no block, near-block or swap with the segment
LEFT, RIGHT = "MKPHC", "FYHPKM"


def run(query: str, frame: str, d0: int, k: int) -> list[Alignment]:
    return verify(encode_query(query), encode_target(frame), d0, k, TABLES)


# ---- 1. queries of <= 4 aa are rejected ------------------------------------------------------


def test_band_is_two() -> None:
    assert BAND == 2


@pytest.mark.parametrize("L", [1, 2, 3, 4])
def test_short_query_raises(L: int) -> None:
    query = "AKGAK"[:L]
    with pytest.raises(ValueError):
        verify(encode_query(query), encode_target("AKGAK"), 0, 0, TABLES)
    with pytest.raises(ValueError):
        verify(encode_query(query), encode_target("MMAKGAKMM"), 2, 1, TABLES)


def test_ak_in_akgak_raises() -> None:
    """The +-2 band around d0 = 0 holds two disjoint exact windows of ``AK`` in ``AKGAK``."""
    with pytest.raises(ValueError):
        verify(encode_query("AK"), encode_target("AKGAK"), 0, 0, TABLES)
    with pytest.raises(ValueError):
        verify(encode_query("AK"), encode_target("AKGAK"), 1, 0, TABLES)


def test_five_residues_accepted() -> None:
    assert run("AKGAK", "AKGAK", 0, 0) == [Alignment(0, 0, 5, ())]
    assert run("AKGAK", "MMAKGAKMM", 2, 1) == [Alignment(0, 2, 7, ())]
    assert run("AKNAK", "AKGGAK", 0, 1) == [Alignment(1, 0, 6, ("block",))]


# ---- 2. stop / X / NOMATCH never consumed, every block shape and slot -------------------------

# (query length, target length) of each operation shape
SHAPES = [(1, 1), (1, 2), (2, 1), (2, 2)]


def _segments(n: int) -> list[str]:
    if n == 1:
        return list(RESIDUES)
    return [a + b for a in RESIDUES for b in RESIDUES]


@pytest.mark.parametrize("special", ["*", "X"])
@pytest.mark.parametrize("shape", SHAPES, ids=[f"{a}to{b}" for a, b in SHAPES])
def test_target_stop_or_x_in_each_slot_never_consumed(shape: tuple[int, int], special: str) -> None:
    """Every query segment against every target segment with a stop/X in one slot: no alignment.

    The flanks force every window to span the stop/X, so any hit means it was consumed.
    """
    qa, tb = shape
    bad = []
    for qseg in _segments(qa):
        for slot in range(tb):
            for other in RESIDUES if tb == 2 else [""]:
                tseg = list(other * (tb - 1))
                tseg.insert(slot, special)
                tseg_s = "".join(tseg)
                out = run(LEFT + qseg + RIGHT, LEFT + tseg_s + RIGHT, 0, 3)
                if out:
                    bad.append((qseg, tseg_s, out))
    assert bad == [], bad[:10]


@pytest.mark.parametrize("shape", SHAPES, ids=[f"{a}to{b}" for a, b in SHAPES])
def test_query_nomatch_in_each_slot_only_by_substitution(shape: tuple[int, int]) -> None:
    """A ``#`` in any query slot costs a substitution; at k = 1 nothing else can be paid."""
    qa, tb = shape
    bad = []
    for slot in range(qa):
        for other in RESIDUES if qa == 2 else [""]:
            q = list(other * (qa - 1))
            q.insert(slot, "B")  # B encodes to NOMATCH ("#")
            qseg = "".join(q)
            for tseg in _segments(tb):
                query = LEFT + qseg + RIGHT
                assert decode(encode_query(query))[len(LEFT) + slot] == "#"
                out = run(query, LEFT + tseg + RIGHT, 0, 1)
                if out and out[0].diff_types != ("substitution",):
                    bad.append((qseg, tseg, out))
                if hit0 := run(query, LEFT + tseg + RIGHT, 0, 0):
                    bad.append((qseg, tseg, hit0))  # a k = 0 hit
    assert bad == [], bad[:10]


def _exchange_pairs() -> list[tuple[str, str]]:
    out = []
    for a, b in (*BLOCKS, *NEAR_BLOCKS):
        for x in sorted({"".join(p) for p in permutations(a)}):
            for y in sorted({"".join(p) for p in permutations(b)}):
                out += [(x, y), (y, x)]
    return out


@pytest.mark.parametrize("special", ["*", "X"])
def test_real_exchanges_with_a_stop_or_x_in_a_target_slot(special: str) -> None:
    for qseg, tseg in _exchange_pairs():
        for slot in range(len(tseg)):
            t = tseg[:slot] + special + tseg[slot + 1 :]
            assert run(LEFT + qseg + RIGHT, LEFT + t + RIGHT, 0, 3) == [], (qseg, t)


def test_real_exchanges_with_nomatch_in_a_query_slot() -> None:
    for qseg, tseg in _exchange_pairs():
        for slot in range(len(qseg)):
            q = qseg[:slot] + "B" + qseg[slot + 1 :]
            out = run(LEFT + q + RIGHT, LEFT + tseg + RIGHT, 0, 1)
            assert out == [] or out[0].diff_types == ("substitution",), (q, tseg, out)


@pytest.mark.parametrize(
    ("qseg", "tseg"),
    [("P", "*G"), ("P", "XG"), ("GS", "*S"), ("GS", "XS"), ("GT", "A*"),
     ("N", "G*"), ("N", "*G"), ("Q", "A*"), ("W", "*D"), ("AG", "*"), ("GG", "X")],
)  # fmt: skip
def test_named_aliasing_cases(qseg: str, tseg: str) -> None:
    """``P`` vs ``*G`` and ``GS`` vs ``*S`` are not blocks (a stop code past the residue range must
    not alias into another table entry)."""
    for k in (1, 2, 3):
        assert run(LEFT + qseg + RIGHT, LEFT + tseg + RIGHT, 0, k) == []


# ---- 3. the kernel without numba agrees with the jitted kernel --------------------------------


def _jit_cases() -> list[tuple[str, str, int, int]]:
    cases = [
        ("AKLMNPQRSTVW", "AKLMPNQRSTVW", 0, 2),  # swap
        ("PEPTNDEKAR", "PEPTGGDEKAR", 0, 1),  # block 1 -> 2
        ("VTEGTAPKLR", "VTEASAPKLR", 0, 1),  # block 2 <-> 2
        ("SWAPKELLTR", "SADAPKELLTR", 0, 1),  # near_block 1 -> 2
        ("MKRLLAEQTR", "MKGVLLAEQTR", 0, 1),  # near_block 1 -> 2
        ("PEPTIDEKAR", "PEPTLDEKAR", 0, 1),  # I = L
        ("PEPTIDEKAR", "PAPTIDEKGR", 0, 1),  # cost 2 > k
        ("PEPTIDEKAR", "PAPTIDEKGR", 0, 2),  # two substitutions
        ("AKLMNPQRSTVW", "AKLMPNQRSTVF", 0, 2),  # substitution + swap
        ("PEPTIDEKAR", "PEPT*DEKAR", 0, 1),  # stop
        ("AGLLKDEPTR", "QLLKDEPTR", 0, 1),  # block 2 -> 1
        ("AKNLMPQNSTVW", "CHAKGGLMPQGGSTVWHC", 2, 2),  # two blocks, diagonal +2
        ("AKGGLMPQGGSTVW", "CHAKNLMPQNSTVWHC", 2, 2),  # two blocks, diagonal -2
        ("MKWLHPCRYHPK", "MKADLHPCGVYHPK", 0, 2),  # two near blocks
        ("PEPTNDEKARW", "PEPTGGDEKARF", 0, 2),  # substitution + block
        ("PEPTBDEKAR", "PEPTLDEKAR", 0, 1),  # query NOMATCH
        ("PEPTJGDEKAR", "PEPTQGDEKAR", 0, 1),  # NOMATCH never in a block
        ("PEPTLDEKARN", "PEPTLDEKARGG", 0, 1),  # precedence tie
        ("GGGGGGGGGG", "GGGGGGGGGGGG", 1, 0),  # t_start tie
        ("PEPTLDEKMHAGGQ", "PEPTLDEKMHQGANCC", 0, 2),  # t_end tie
        ("PEPTLDEKAR", "MMMMMPEPTLDEKAR", 7, 0),  # band past the end
        ("PEPTLDEKAR", "MMMMMMMMPEPTLDEKARMMMMMMMM", 11, 0),  # outside the band
        ("AKGAK", "MMAKGAKMM", 2, 1),  # L = 5
    ]
    for qseg, tseg in _exchange_pairs()[::3]:
        cases.append((LEFT + qseg + RIGHT, "CHY" + LEFT + tseg + RIGHT + "YHC", 3, 1))
    rng = np.random.default_rng(39)
    for _ in range(20):
        frame = "".join(RESIDUES[i] for i in rng.integers(0, 19, size=30))
        s, L = int(rng.integers(0, 15)), int(rng.integers(6, 14))
        q = list(frame[s : s + L])
        for _ in range(int(rng.integers(0, 3))):
            q[int(rng.integers(L))] = RESIDUES[int(rng.integers(19))]
        cases.append(("".join(q), frame, s, 2))
    return cases


_SCRIPT = """
import json, sys
from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import encode_query, encode_target
from giap.mapping.tolerant_dp import verify
tables = block_tables()
out = []
for q, f, d0, k in json.load(sys.stdin):
    res = verify(encode_query(q), encode_target(f), d0, k, tables)
    out.append([[a.cost, a.t_start, a.t_end, list(a.diff_types)] for a in res])
print(json.dumps(out))
"""


def _as_lists(res: list[Alignment]) -> list[list[object]]:
    return [[a.cost, a.t_start, a.t_end, list(a.diff_types)] for a in res]


def test_kernel_without_numba_jit_agrees_with_jitted() -> None:
    cases = _jit_cases()
    assert len(cases) >= 20
    jitted = [_as_lists(run(*c)) for c in cases]
    env = dict(os.environ)
    env["NUMBA_DISABLE_JIT"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT],
        input=json.dumps(cases),
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    pure = json.loads(proc.stdout)
    diffs = [(c, j, p) for c, j, p in zip(cases, jitted, pure, strict=True) if j != p]
    assert diffs == [], diffs[:5]
    # the subset exercises every difference type and some misses
    seen = {t for res in jitted for a in res for t in a[3]}  # type: ignore[attr-defined]
    assert seen == {"substitution", "swap", "block", "near_block"}
    assert any(res == [] for res in jitted)


# ---- 4. dtype and stride normalisation ---------------------------------------------------------

DTYPE_CASES = [
    ("AKLMNPQRSTVW", "CHYAKLMPNQRSTVFYHC", 3, 2),
    ("PEPTNDEKAR", "CHYPEPTGGDEKARYHC", 3, 1),
    ("AGLLKDEPTR", "QLLKDEPTR", 0, 1),
    ("SWAPKELLTR", "SADAPKELLTR", 0, 1),
    ("PEPTBDEKAR", "PEPTLDEKAR", 0, 1),
    ("PEPTIDEKAR", "PEPT*DEKAR", 0, 1),
]


def _strided(a: np.ndarray) -> np.ndarray:
    buf = np.zeros(2 * a.shape[0], dtype=a.dtype)
    buf[::2] = a
    out = buf[::2]
    assert not out.flags.c_contiguous
    return out


def _variants(q: np.ndarray, f: np.ndarray) -> list[tuple[str, np.ndarray, np.ndarray]]:
    q64, f64 = q.astype(np.int64), f.astype(np.int64)
    return [
        ("int64", q64, f64),
        ("int64 query", q64, f),
        ("int64 frame", q, f64),
        ("strided", _strided(q), _strided(f)),
        ("strided int64", _strided(q64), _strided(f64)),
        ("reversed view", q[::-1].copy()[::-1], f[::-1].copy()[::-1]),
    ]


@pytest.mark.parametrize("case", DTYPE_CASES, ids=[c[0] for c in DTYPE_CASES])
def test_int64_and_strided_inputs_match_uint8(case: tuple[str, str, int, int]) -> None:
    query, frame, d0, k = case
    q, f = encode_query(query), encode_target(frame)
    assert q.dtype == np.uint8 and f.dtype == np.uint8
    want = verify(np.ascontiguousarray(q), np.ascontiguousarray(f), d0, k, TABLES)
    for name, qv, fv in _variants(q, f):
        assert verify(qv, fv, d0, k, TABLES) == want, name


_COMPILE_SCRIPT = """
import json
import numpy as np
from giap.mapping import tolerant_dp
from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import encode_query, encode_target
tables = block_tables()
kernels = [v for v in vars(tolerant_dp).values()
           if hasattr(v, "signatures") and hasattr(v, "py_func")]
q, f = encode_query("AKLMNPQRSTVW"), encode_target("CHYAKLMPNQRSTVFYHC")
tolerant_dp.verify(q, f, 3, 2, tables)  # the uint8 C-contiguous specialisation
before = sum(len(d.signatures) for d in kernels)
def strided(a):
    buf = np.zeros(2 * a.shape[0], dtype=a.dtype)
    buf[::2] = a
    return buf[::2]
q64, f64 = q.astype(np.int64), f.astype(np.int64)
for qv, fv in ((q64, f64), (q64, f), (q, f64), (strided(q), strided(f)),
               (strided(q64), strided(f64)), (q.astype(np.int32), f.astype(np.int16))):
    tolerant_dp.verify(qv, fv, 3, 2, tables)
after = sum(len(d.signatures) for d in kernels)
print(json.dumps([len(kernels), before, after]))
"""


def test_other_dtypes_and_strides_compile_no_new_kernel() -> None:
    """Fresh process (other tests may already have compiled other specialisations)."""
    env = dict(os.environ)
    env.pop("NUMBA_DISABLE_JIT", None)
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", _COMPILE_SCRIPT],
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    n_kernels, before, after = json.loads(proc.stdout)
    if n_kernels == 0:
        pytest.skip("no numba dispatcher found in tolerant_dp (kernel renamed or moved)")
    assert before >= 1
    assert after == before, f"{after - before} new kernel specialisation(s) compiled"


# ---- 5. t_end tie-break -----------------------------------------------------------------------


def test_t_end_tie_break_keeps_the_smallest_end() -> None:
    """Windows [0, 13) and [0, 14) both cost (substitution, block) from t_start 0."""
    query, frame = "PEPTLDEKMHAGGQ", "PEPTLDEKMHQGANCC"
    assert run(query, frame[:13], 0, 2) == [Alignment(2, 0, 13, ("substitution", "block"))]
    assert run(query, frame[:14] + "**", 0, 2)[0].diff_types == ("substitution", "block")
    assert run(query, frame, 0, 2) == [Alignment(2, 0, 13, ("substitution", "block"))]


@pytest.mark.parametrize("flank", ["", "CHY"])
def test_t_end_tie_break_more_cases(flank: str) -> None:
    for query, window in [
        ("PEPTLDEKMHAGGQ", "PEPTLDEKMHANAN"),
        ("PEPTLDEKMHAGAN", "PEPTLDEKMHAQGQ"),
        ("PEPTLDEKMHANAG", "PEPTLDEKMHAGQA"),
    ]:
        frame = flank + window + "CC"
        out = run(query, frame, len(flank), 2)
        assert len(out) == 1, out
        aln = out[0]
        assert (aln.cost, aln.diff_types) == (2, ("substitution", "block"))
        assert (aln.t_start, aln.t_end) == (len(flank), len(flank) + 13)
