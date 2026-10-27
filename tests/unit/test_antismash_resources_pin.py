"""antiSMASH resource request (#416 D2).

The module's ``memory { 32.GB * task.attempt }`` was dead: base.config's ``withLabel
'process_high'`` (16 CPUs, 48 GB) outranks a module directive. The deployed run peaked far below
that, so base.config sizes ANTISMASH itself (``withName: 'ANTISMASH'``: 4 CPUs, 8 GB x attempt)
and the module carries no memory directive. FarmShare keeps its max(cpus, ceil(memory / 4000 MB))
rule for the selector (base.config's withName CPUs outrank its label rule), as for Promotech.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

WF = Path(__file__).resolve().parents[2] / "workflow"
CONF = WF / "conf"
MODULE = WF / "modules" / "local" / "antismash" / "main.nf"


def _strip_comments(text: str) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)


def _block(text: str, selector: str) -> str:
    """Body of ``withName: '<selector>' { ... }`` with balanced braces (closures included)."""
    text = _strip_comments(text)
    m = re.search(rf"withName\s*:\s*(['\"]){re.escape(selector)}\1\s*\{{", text)
    assert m, f"no withName: '{selector}' block"
    depth, start = 0, m.end() - 1
    for i in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[start + 1 : i]
    raise AssertionError("unbalanced braces")


def test_base_config_sizes_antismash() -> None:
    body = _block((CONF / "base.config").read_text(encoding="utf-8"), "ANTISMASH")
    assert re.search(r"\bcpus\s*=\s*4\b", body), body
    assert re.search(r"memory\s*=\s*\{\s*8\.GB\s*\*\s*task\.attempt\s*\}", body), body


def test_module_has_no_dead_memory_directive() -> None:
    text = _strip_comments(MODULE.read_text(encoding="utf-8"))
    process = text.split("input:", 1)[0]  # the directives sit above input:
    assert not re.search(r"^\s*memory\b", process, re.M), process
    assert not re.search(r"^\s*cpus\b", process, re.M), process


def test_farmshare_derives_antismash_cpus_from_memory() -> None:
    body = _block((CONF / "farmshare.config").read_text(encoding="utf-8"), "ANTISMASH")
    m = re.search(
        r"Math\.max\(\s*(\d+)\s*,\s*Math\.ceil\(\s*task\.memory\.toBytes\(\)\s*/\s*"
        r"4000\.MB\.toBytes\(\)\s*\)",
        body,
    )
    assert m and int(m[1]) == 4, body
    # 8 / 16 / 24 GB at attempts 1 / 2 / 3 against MaxMemPerCPU = 4000 MB
    assert [max(4, math.ceil(8 * 1024 * a / 4000)) for a in (1, 2, 3)] == [4, 5, 7]
