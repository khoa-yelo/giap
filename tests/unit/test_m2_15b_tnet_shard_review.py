"""M2-15b review (#455 S2): ``plan_problems`` - MERGE checks plan.tsv against the re-read genome
before writing, so a stale or edited plan fails (exit 1) instead of dropping or misplacing a
block."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "bin" / "gippa_terminatornet_shard.py"
RECORDS = [("chr", 2_000), ("pl", 100), ("tiny", 2)]


def _shard() -> ModuleType:
    name = "gippa_terminatornet_shard"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _rows() -> list[dict]:
    return _shard().plan_rows(RECORDS, "long", 500)


def test_a_fresh_plan_has_no_problems() -> None:
    assert _shard().plan_problems(_rows(), RECORDS) == []


def _drop(rows: list[dict]) -> list[dict]:
    return [
        r
        for r in rows
        if not (r["seq_index"] == 0 and r["strand"] == "-" and r["core_start"] == 500)
    ]


def _dup(rows: list[dict]) -> list[dict]:
    return [*rows, dict(rows[0])]


def _length(rows: list[dict]) -> list[dict]:
    return [{**r, "seq_length": 1_999} if r["seq_index"] == 0 else r for r in rows]


def _name(rows: list[dict]) -> list[dict]:
    return [{**r, "seq_name": "other"} if r["seq_index"] == 1 else r for r in rows]


def _index(rows: list[dict]) -> list[dict]:
    return [*rows, {**rows[-1], "chunk_id": "c00009p0000", "seq_index": 9}]


def _record(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r["seq_index"] != 2]


def _short(rows: list[dict]) -> list[dict]:
    last = max(r["core_end"] for r in rows if r["seq_index"] == 0 and r["strand"] == "+")
    return [
        {**r, "core_end": last - 1}
        if r["seq_index"] == 0 and r["strand"] == "+" and r["core_end"] == last
        else r
        for r in rows
    ]


@pytest.mark.parametrize(
    ("edit", "word"),
    [(_drop, "tile"), (_dup, "duplicate"), (_length, "genome has"), (_name, "genome has"),
     (_index, "not in the genome"), (_record, "tile"), (_short, "tile")],
)  # fmt: skip
def test_a_plan_that_does_not_match_the_genome_is_named(edit, word: str) -> None:
    problems = _shard().plan_problems(edit(_rows()), RECORDS)
    assert problems and any(word in p for p in problems), problems
