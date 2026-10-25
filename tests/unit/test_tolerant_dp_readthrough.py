"""The banded verifier with a consumable stop (M7-16a rule 4; ADR-0058 item 1).

A marked cell carries ``READTHROUGH = 22`` in the search text. Only the (1, 1) move gains a case:
over code 22 it is a substitution (cost 1, rank 0) for any query code, ``NOMATCH`` included. A
marked cell never matches, and no swap, block or near-block reads it. Spec examples 1, 3, 4, 4b, 4c
on the toy contig ``c1`` (frame ``+0`` = ``MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*AAGLKEPSTR*``, the TGA at
p = 18 marked, p = 30 and 41 not), every block shape and a swap over or next to a marked cell, a
query ``U`` over a marked cell, the unmarked (code 19) control against the M3-06 oracle, and a
``NUMBA_DISABLE_JIT=1`` subset.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from helpers.readthrough_toys import C1_AA

from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import STOP, decode, encode_query, encode_target
from giap.mapping.tolerant_dp import Alignment, exchange_lookup, verify, verify_batch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "property"))
from oracle_tolerant import oracle_verify

TABLES = block_tables()
SRC = Path(__file__).resolve().parents[2] / "src"


def _rt_code() -> int:
    from giap.mapping.sixframe import READTHROUGH

    return int(READTHROUGH)


def frame_rt(letters: str, marked: tuple[int, ...]) -> np.ndarray:
    """``encode_target(letters)`` with the marked ``*`` cells set to ``READTHROUGH``."""
    t = encode_target(letters)
    for p in marked:
        assert t[p] == STOP, (letters, p)
        t[p] = _rt_code()
    return t


def run(query: str, letters: str, marked: tuple[int, ...], d0: int, k: int) -> list[Alignment]:
    return verify(encode_query(query), frame_rt(letters, marked), d0, k, TABLES)


# ---------------------------------------------------------------- spec examples (toy contig c1)

EXAMPLES = [
    # (example, query, d0, k, expected with p = 18 marked)
    ("1", "WQLNEFKWVPAGHTL", 11, 2, [Alignment(1, 11, 26, ("substitution",))]),
    ("3", "HTLDMRKQAAGLKE", 23, 2, []),
    ("4", "NEFKWVPAGHA", 14, 1, []),
    ("4b", "NEFKWVPAGHT", 14, 1, [Alignment(1, 14, 25, ("substitution",))]),
    ("4c", "NEFKWVPAGHAL", 14, 2, [Alignment(2, 14, 26, ("substitution", "substitution"))]),
]


@pytest.mark.parametrize(
    ("name", "query", "d0", "k", "want"), EXAMPLES, ids=[e[0] for e in EXAMPLES]
)
def test_spec_examples_at_verify_level(
    name: str, query: str, d0: int, k: int, want: list[Alignment]
) -> None:
    assert run(query, C1_AA, (18,), d0, k) == want


@pytest.mark.parametrize(
    ("name", "query", "d0", "k", "want"), EXAMPLES, ids=[e[0] for e in EXAMPLES]
)
def test_unmarked_stops_keep_the_m3_06_result(
    name: str, query: str, d0: int, k: int, want: list[Alignment]
) -> None:
    """Code 19 everywhere: the M3-06 behaviour (``*`` never consumed), equal to its oracle."""
    got = run(query, C1_AA, (), d0, k)
    assert oracle_verify(decode(encode_query(query)), C1_AA, d0, k) is None
    assert got == []  # every example window holds a '*'


@pytest.mark.parametrize("letter", ["U", "B", "X", "Z", "*", "O", "W"])
def test_any_query_code_over_a_marked_cell_costs_one(letter: str) -> None:
    """Rule 4: 'for any query code, NOMATCH included' (U, the Sec letter; example 9's X)."""
    assert run(f"NEFK{letter}VPAGHT", C1_AA, (18,), 14, 1) == [
        Alignment(1, 14, 25, ("substitution",))
    ]


def test_two_marked_cells_cost_two() -> None:
    """Edge case: two marked cells in one window need k = 2 and no other difference."""
    t = "CFHKMPY*CFHKM*PYCFH"
    q = "CFHKMPYWCFHKMWPYCFH"
    assert run(q, t, (7, 13), 0, 2) == [Alignment(2, 0, 19, ("substitution", "substitution"))]
    assert run(q, t, (7, 13), 0, 1) == []
    assert run(q, t, (7,), 0, 2) == []  # the unmarked one is never consumed


# ---------------------------------------------------------------- no block or swap reads a cell
# fillers from CFHKMPY take part in no BLOCKS / NEAR_BLOCKS pair

SHAPES = [
    # (shape, query, control target (the block holds), target with the cell, cell, k, want)
    ("1->2 G*", "CFHKMNPYCFH", "CFHKMGGPYCFH", "CFHKMG*PYCFH", 6, 1, []),
    ("1->2 *G", "CFHKMNPYCFH", "CFHKMGGPYCFH", "CFHKM*GPYCFH", 5, 1, []),
    ("2->1 *", "CFHKMAGPYCFH", "CFHKMQPYCFH", "CFHKM*PYCFH", 5, 2, []),
    ("2<->2 G*", "CFHKMASPYCFH", "CFHKMGTPYCFH", "CFHKMG*PYCFH", 6, 2,
     [Alignment(2, 0, 12, ("substitution", "substitution"))]),
    ("near 1->2 A*", "CFHKMWPYCFH", "CFHKMADPYCFH", "CFHKMA*PYCFH", 6, 1, []),
    ("near 2->1 *", "CFHKMGVPYCFH", "CFHKMRPYCFH", "CFHKM*PYCFH", 5, 2, []),
]  # fmt: skip


@pytest.mark.parametrize(
    ("shape", "query", "control", "target", "cell", "k", "want"), SHAPES, ids=[s[0] for s in SHAPES]
)
def test_a_marked_cell_inside_a_block_is_never_consumed(
    shape: str, query: str, control: str, target: str, cell: int, k: int, want: list[Alignment]
) -> None:
    ctl = verify(encode_query(query), encode_target(control), 0, k, TABLES)
    assert len(ctl) == 1 and ctl[0].cost == 1, ctl  # the construction is one block
    assert run(query, target, (cell,), 0, k) == want


def test_swap_and_near_block_next_to_a_marked_cell() -> None:
    assert run("CFHKMKWYPCFH", "CFHKM*YWPCFH", (5,), 0, 2) == [
        Alignment(2, 0, 12, ("substitution", "swap"))
    ]
    assert run("CFHKMKWPCFH", "CFHKM*ADPCFH", (5,), 0, 2) == [
        Alignment(2, 0, 12, ("substitution", "near_block"))
    ]


def test_verify_batch_rows_equal_verify_with_marked_cells() -> None:
    q = np.ascontiguousarray(encode_query("WQLNEFKWVPAGHTL"), dtype=np.uint8)
    t = np.ascontiguousarray(frame_rt(C1_AA, (18,)), dtype=np.uint8)
    out = verify_batch(q, t, np.array([9, 10, 11, 12, 13, 30], dtype=np.int64), 2,
                       exchange_lookup(TABLES))  # fmt: skip
    for d0, row in zip([9, 10, 11, 12, 13, 30], out, strict=True):
        got = verify(q, t, d0, 2, TABLES)
        if got:
            assert (int(row[1]), int(row[2])) == (got[0].t_start, got[0].t_end)
        else:
            assert list(row) == [-1, -1, -1]
    assert int(out[2, 1]) == 11 and int(out[2, 2]) == 26


# ---------------------------------------------------------------- without numba

_SCRIPT = """
import json, sys
from giap.mapping.diffops import block_tables
from giap.mapping.sixframe import READTHROUGH, encode_query, encode_target
from giap.mapping.tolerant_dp import verify
tables = block_tables()
out = []
for q, f, marked, d0, k in json.load(sys.stdin):
    t = encode_target(f)
    for p in marked:
        t[p] = READTHROUGH
    res = verify(encode_query(q), t, d0, k, tables)
    out.append([[a.cost, a.t_start, a.t_end, list(a.diff_types)] for a in res])
print(json.dumps(out))
"""


def test_kernel_without_numba_jit_agrees_on_marked_cells() -> None:
    cases = [[q, C1_AA, [18], d0, k] for _, q, d0, k, _ in EXAMPLES]
    cases += [[q, t, [c], 0, k] for _, q, _, t, c, k, _ in SHAPES]
    cases += [["NEFKUVPAGHT", C1_AA, [18], 14, 1], ["CFHKMKWYPCFH", "CFHKM*YWPCFH", [5], 0, 2]]
    jitted = [
        [[a.cost, a.t_start, a.t_end, list(a.diff_types)] for a in run(q, f, tuple(m), d0, k)]
        for q, f, m, d0, k in cases
    ]
    env = dict(os.environ)
    env["NUMBA_DISABLE_JIT"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT], input=json.dumps(cases), capture_output=True, text=True,
        env=env, timeout=600, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert json.loads(proc.stdout) == jitted


# ---------------------------------------------------------------- budget guard on the last residue
# Review round 1 item 1 (Q3, "window ending on the stop"): the (1, 1) move over a marked cell is
# taken only while the cost stays within k. Frame ``NEFKAVPAGHTLW*`` with the ``*`` marked, k = 1:
# ``...TLWQ`` costs 1 (the cell only); ``...TLAQ`` would cost 2 (W->A, then the cell on the last
# query residue) and must give no alignment.

BUDGET_FRAME = "NEFKAVPAGHTLW*"
BUDGET_CASES = [
    # (query, want with p = 13 marked, k = 1, d0 = 0)
    ("NEFKAVPAGHTLWQ", [Alignment(1, 0, 14, ("substitution",))]),
    ("NEFKAVPAGHTLAQ", []),
    ("NEFKAVPAGHTAAQ", []),
]


@pytest.mark.parametrize(("query", "want"), BUDGET_CASES, ids=[c[0] for c in BUDGET_CASES])
def test_a_marked_cell_on_the_last_residue_respects_the_budget(
    query: str, want: list[Alignment]
) -> None:
    assert run(query, BUDGET_FRAME, (13,), 0, 1) == want


def test_verify_batch_respects_the_budget_on_a_marked_last_cell() -> None:
    t = np.ascontiguousarray(frame_rt(BUDGET_FRAME, (13,)), dtype=np.uint8)
    d0s = np.array([0], dtype=np.int64)
    ex = exchange_lookup(TABLES)
    over = np.ascontiguousarray(encode_query("NEFKAVPAGHTLAQ"), dtype=np.uint8)
    assert [list(map(int, r)) for r in verify_batch(over, t, d0s, 1, ex)] == [[-1, -1, -1]]
    ok = np.ascontiguousarray(encode_query("NEFKAVPAGHTLWQ"), dtype=np.uint8)
    row = verify_batch(ok, t, d0s, 1, ex)[0]
    assert (int(row[1]), int(row[2])) == (0, 14)


def test_kernel_without_numba_jit_respects_the_budget_on_a_marked_last_cell() -> None:
    cases = [[q, BUDGET_FRAME, [13], 0, 1] for q, _ in BUDGET_CASES]
    want = [[[a.cost, a.t_start, a.t_end, list(a.diff_types)] for a in w] for _, w in BUDGET_CASES]
    env = dict(os.environ)
    env["NUMBA_DISABLE_JIT"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT], input=json.dumps(cases), capture_output=True, text=True,
        env=env, timeout=600, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert json.loads(proc.stdout) == want
