"""M2-15b: the pure plan and merge-order functions of the sharded TerminatorNet scan
(docs/specs/M2-15b-terminatornet-sharding.md items 1, 2, 4; examples 1-5, 7).

``workflow/bin/gippa_terminatornet_shard.py`` runs in the ``gippa-terminatornet`` env and is not
part of ``src/gippa``; its pure functions import nothing from upstream, so this file loads it by
path.
Interface pinned here (test author, from the spec):

* ``CORE = 250_000`` (windows per chunk core) and ``nwin(length) = max(length - 59, 1)``;
* ``chunk_id(seq_index, strand, k)`` -> ``c{seq_index:05d}{p|m}{k:04d}``;
* ``plan_rows(records, mode, core=CORE)``: ``records`` = ``[(seq_name, cleaned_length), ...]`` in
  the order upstream keeps them; returns one dict per chunk with the ``PLAN_COLUMNS`` keys (ints for
  ``seq_index``, ``seq_length``, ``core_start``, ``core_end``), ordered by
  (seq_index, + before -, k);
* ``output_blocks(rows)``: the plan rows grouped per (seq_index, strand), each group in k order,
  groups in upstream output order (every + block in FASTA order, then every - block);
* ``order_calls(calls_by_chunk, mode)``: the chunks' call lists concatenated in the given (k) order,
  stable-sorted by stop (``itemgetter(1)``) when ``mode == 'long'``, untouched when ``'short'``.
"""

from __future__ import annotations

import importlib.util
import sys
from operator import itemgetter
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "bin" / "gippa_terminatornet_shard.py"
COLUMNS = (
    "chunk_id",
    "seq_index",
    "seq_name",
    "strand",
    "seq_length",
    "core_start",
    "core_end",
    "mode",
)


def _shard() -> ModuleType:
    name = "gippa_terminatornet_shard"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    assert spec is not None and spec.loader is not None, f"cannot load {SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _cores(rows: list[dict], seq_index: int, strand: str) -> list[tuple[int, int]]:
    return [
        (r["core_start"], r["core_end"])
        for r in rows
        if (r["seq_index"], r["strand"]) == (seq_index, strand)
    ]


def _check_partition(rows: list[dict], lengths: list[int], core: int) -> None:
    sh = _shard()
    for i, length in enumerate(lengths):
        for strand in "+-":
            cores = _cores(rows, i, strand)
            assert cores, f"no chunk for record {i} strand {strand}"
            assert cores[0][0] == 0 and cores[-1][1] == sh.nwin(length), cores
            assert all(a < b for a, b in cores), f"empty core in {cores}"
            assert all(cores[k][1] == cores[k + 1][0] for k in range(len(cores) - 1)), cores
            assert all(b - a <= core for a, b in cores), cores
            assert all(a % core == 0 for a, _ in cores), "core k starts at k*CORE"


def test_constants_and_nwin() -> None:
    sh = _shard()
    assert sh.CORE == 250_000
    assert [sh.nwin(n) for n in (0, 2, 40, 59, 60, 61, 1_000_000)] == [1, 1, 1, 1, 1, 2, 999_941]


def test_plan_header_columns() -> None:
    assert tuple(_shard().PLAN_COLUMNS) == COLUMNS


def test_example_1_one_megabase_contig_four_cores_per_strand() -> None:
    rows = _shard().plan_rows([("chr1", 1_000_000)], "long")
    want = [(0, 250_000), (250_000, 500_000), (500_000, 750_000), (750_000, 999_941)]
    assert _cores(rows, 0, "+") == want and _cores(rows, 0, "-") == want
    assert [r["chunk_id"] for r in rows] == [f"c00000p{k:04d}" for k in range(4)] + [
        f"c00000m{k:04d}" for k in range(4)
    ]
    assert all(r["seq_length"] == 1_000_000 and r["seq_name"] == "chr1" for r in rows)
    assert all(set(r) == set(COLUMNS) for r in rows)


def test_example_2_plasmid_shorter_than_one_chunk() -> None:
    rows = _shard().plan_rows([("p1", 120_000)], "long")
    assert [(r["strand"], r["core_start"], r["core_end"]) for r in rows] == [
        ("+", 0, 119_941),
        ("-", 0, 119_941),
    ]


@pytest.mark.parametrize(("name", "length"), [("short40", 40), ("t", 0)])
def test_examples_3_and_5_shorter_than_a_window_or_cleaned_empty_one_window(
    name: str, length: int
) -> None:
    rows = _shard().plan_rows([(name, length)], "short")
    assert [(r["chunk_id"], r["core_start"], r["core_end"], r["seq_length"]) for r in rows] == [
        ("c00000p0000", 0, 1, length),
        ("c00000m0000", 0, 1, length),
    ]


def test_example_4_seq_index_counts_the_kept_records() -> None:
    # upstream already dropped the empty record between a and b: b is record 1
    rows = _shard().plan_rows([("a", 300), ("b", 120)], "short")
    assert [(r["chunk_id"], r["seq_index"], r["seq_name"]) for r in rows] == [
        ("c00000p0000", 0, "a"),
        ("c00000m0000", 0, "a"),
        ("c00001p0000", 1, "b"),
        ("c00001m0000", 1, "b"),
    ]


@pytest.mark.parametrize("mode", ["long", "short"])
def test_mode_is_copied_to_every_row(mode: str) -> None:
    rows = _shard().plan_rows([("a", 3000), ("b", 40)], mode, core=500)
    assert {r["mode"] for r in rows} == {mode}


@settings(max_examples=150, deadline=None)
@given(
    lengths=st.lists(st.integers(0, 3_000), min_size=1, max_size=4),
    core=st.sampled_from([1, 7, 61, 1000]),
    mode=st.sampled_from(["long", "short"]),
)
def test_cores_partition_windows_per_contig_and_strand(
    lengths: list[int], core: int, mode: str
) -> None:
    sh = _shard()
    records = [(f"s{i}", n) for i, n in enumerate(lengths)]
    rows = sh.plan_rows(records, mode, core=core)
    _check_partition(rows, lengths, core)
    ids = [r["chunk_id"] for r in rows]
    assert len(set(ids)) == len(ids), "chunk ids are not unique"
    order = [(r["seq_index"], "+-".index(r["strand"]), r["core_start"]) for r in rows]
    assert order == sorted(order), "rows not ordered by (seq_index, + before -, k)"
    for r in rows:
        k = r["core_start"] // core
        assert r["chunk_id"] == sh.chunk_id(r["seq_index"], r["strand"], k)
        assert r["chunk_id"] == f"c{r['seq_index']:05d}{'p' if r['strand'] == '+' else 'm'}{k:04d}"
        assert r["seq_length"] == lengths[r["seq_index"]] and r["seq_name"] == f"s{r['seq_index']}"
    for i in range(len(lengths)):
        for strand in "+-":
            block = [r["chunk_id"] for r in rows if (r["seq_index"], r["strand"]) == (i, strand)]
            assert block == sorted(block), "within a block, id order is k order"


@settings(max_examples=60, deadline=None)
@given(
    lengths=st.lists(st.integers(0, 3_000), min_size=1, max_size=4),
    core=st.sampled_from([7, 61, 1000]),
)
def test_output_blocks_plus_blocks_in_fasta_order_then_minus_blocks(
    lengths: list[int], core: int
) -> None:
    sh = _shard()
    rows = sh.plan_rows([(f"s{i}", n) for i, n in enumerate(lengths)], "long", core=core)
    blocks = sh.output_blocks(rows)
    keys = [(b[0]["seq_index"], b[0]["strand"]) for b in blocks]
    n = len(lengths)
    assert keys == [(i, "+") for i in range(n)] + [(i, "-") for i in range(n)]
    for b in blocks:
        assert len({(r["seq_index"], r["strand"]) for r in b}) == 1
        assert [r["core_start"] for r in b] == sorted(r["core_start"] for r in b)
    assert sorted(r["chunk_id"] for b in blocks for r in b) == sorted(r["chunk_id"] for r in rows)


def _call(start: int, stop: int, pred: float, tag: str = "") -> tuple:
    return (start, stop, "AAAAAAAA", "GCGC" + tag, "TTTTTTTT", "((..))", -10.0, pred)


def test_example_7_calls_cross_a_chunk_edge_in_core_order() -> None:
    a, b = _call(980, 1010, 0.93), _call(1000, 1030, 0.93)
    for mode in ("long", "short"):
        assert _shard().order_calls([[a], [b]], mode) == [a, b], (
            "later call B must follow A (pred >= best wins)"
        )


def test_long_mode_stable_sort_by_stop_short_mode_window_order() -> None:
    sh = _shard()
    c0 = [_call(10, 60, 0.7, "a"), _call(12, 40, 0.8, "b"), _call(14, 60, 0.6, "c")]
    c1 = [_call(30, 40, 0.9, "d"), _call(35, 50, 0.55, "e")]
    flat = c0 + c1
    assert sh.order_calls([c0, c1], "short") == flat
    got = sh.order_calls([c0, c1], "long")
    assert got == sorted(flat, key=itemgetter(1))
    assert [c[3] for c in got] == ["GCGCb", "GCGCd", "GCGCe", "GCGCa", "GCGCc"], (
        "sort must be stable"
    )


def test_order_calls_empty_chunks() -> None:
    assert _shard().order_calls([[], []], "long") == []
    assert _shard().order_calls([], "short") == []
