"""Slurm time-limit kills (exit 140) retry like 143 (#416 D3).

barley-01 killed a task at its Slurm time limit with exit 140, which base.config's errorStrategy
(137/143 only) terminated. 140 must take the 143 branch: retry while memory is below its cap, or
while a set time is below the time cap. The closure is evaluated here by a small Groovy-to-Python
translation of its one-line ``COND ? 'retry' : 'terminate'`` form (base.config keeps it on one
line for ``nextflow config -flat``), so the pin holds whatever list syntax carries 140.

barley-01 also joins the slow/failing-node exclusion for pipeline task submission
(``farmshare.config`` clusterOptions ``--exclude``).
"""

from __future__ import annotations

import itertools
import re
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

CONF = Path(__file__).resolve().parents[2] / "workflow" / "conf"
GB = 2**30
H = 3600_000  # ms


def _strip_comments(text: str) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)


def _closure() -> str:
    text = _strip_comments((CONF / "base.config").read_text(encoding="utf-8"))
    m = re.search(r"^\s*errorStrategy\s*=\s*\{(.*)\}\s*$", text, re.M)
    assert m, "no one-line errorStrategy closure in base.config"
    return m[1].strip()


def _safe(obj: Any, attr: str) -> Any:
    return None if obj is None else getattr(obj, attr)


def _compile(src: str) -> Any:
    m = re.fullmatch(r"(.*)\?\s*'retry'\s*:\s*'terminate'", src.replace("?.", "\0"), re.S)
    assert m, f"errorStrategy is not COND ? 'retry' : 'terminate': {src}"
    cond = m[1].replace("\0", "?.")
    cond = re.sub(r"\s+as\s+[\w.]+", "", cond)
    for _ in range(3):
        cond = re.sub(r"([\w.]+)\?\.(\w+)", r"_safe(\1, '\2')", cond)
    cond = re.sub(r"\bnull\b", "None", cond)
    cond = cond.replace("&&", " and ").replace("||", " or ")
    cond = re.sub(r"!(?!=)", " not ", cond)
    code = compile(cond, "errorStrategy", "eval")
    return lambda task: "retry" if eval(code, {"_safe": _safe}, {"task": task}) else "terminate"


STRATEGY = None


def _strategy(**kw: Any) -> str:
    global STRATEGY
    if STRATEGY is None:
        STRATEGY = _compile(_closure())
    limits = kw.pop("limits")
    return str(STRATEGY(SimpleNamespace(resourceLimits=limits, **kw)))


def _states() -> list[dict[str, Any]]:
    out = []
    for mem, time, lim_mem, lim_time, has_limits in itertools.product(
        (8 * GB, 64 * GB), (None, 12 * H, 48 * H), (64 * GB, None), (48 * H, None), (True, False)
    ):
        limits = SimpleNamespace(memory=lim_mem, time=lim_time) if has_limits else None
        out.append({"memory": mem, "time": time, "limits": limits})
    return out


def test_translation_reproduces_the_137_143_policy() -> None:
    at_cap = SimpleNamespace(memory=64 * GB, time=48 * H)
    assert _strategy(exitStatus=137, memory=8 * GB, time=None, limits=at_cap) == "retry"
    assert _strategy(exitStatus=137, memory=64 * GB, time=12 * H, limits=at_cap) == "terminate"
    assert _strategy(exitStatus=143, memory=64 * GB, time=12 * H, limits=at_cap) == "retry"
    assert _strategy(exitStatus=143, memory=64 * GB, time=48 * H, limits=at_cap) == "terminate"
    assert _strategy(exitStatus=1, memory=8 * GB, time=12 * H, limits=at_cap) == "terminate"


def test_exit_140_retries_on_time_below_the_cap() -> None:
    at_mem_cap = SimpleNamespace(memory=64 * GB, time=48 * H)
    assert _strategy(exitStatus=140, memory=64 * GB, time=12 * H, limits=at_mem_cap) == "retry"
    assert _strategy(exitStatus=140, memory=64 * GB, time=48 * H, limits=at_mem_cap) == "terminate"


@pytest.mark.parametrize("state", _states(), ids=lambda s: repr(s)[:60])
def test_exit_140_follows_the_143_branch(state: dict[str, Any]) -> None:
    assert _strategy(exitStatus=140, **dict(state)) == _strategy(exitStatus=143, **dict(state))


def _excluded_nodes(text: str) -> set[str]:
    text = _strip_comments(text)
    nodes: set[str] = set()
    for m in re.finditer(r"(?:--exclude[= ]|(?<!\S)-x\s*)([\w,\-\[\]]+)", text):
        nodes.update(n for n in m[1].split(",") if n)
    return nodes


def test_farmshare_excludes_barley_01_for_pipeline_tasks() -> None:
    text = (CONF / "farmshare.config").read_text(encoding="utf-8")
    opts = re.search(r"^\s*clusterOptions\s*=\s*(.+)$", _strip_comments(text), re.M)
    assert opts, "no clusterOptions in farmshare.config"
    assert "barley-01" in _excluded_nodes(opts[1]), opts[1]
    assert "--qos=normal" in opts[1], opts[1]
