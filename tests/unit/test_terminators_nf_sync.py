"""``terminatornet_cmd`` / ``transtermhp_cmd`` builders vs. their Nextflow processes (M2-15
acceptance "test_nf_cmd_sync covers TerminatorNet (-s present) and TransTermHP"; M2-07 sync rule,
M2-index decision 2: every ``giap.annotate.<tool>_cmd`` builder is kept in sync with
``workflow/modules/local/<tool>/main.nf``).

Skips until the processes exist (M2-07). G12: TerminatorNet's script block must pass ``-s``; it
passes only the documented flags (its parser matches flags by substring).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from giap.annotate.terminatornet_cmd import terminatornet_command
from giap.annotate.transtermhp_cmd import transterm_command

REPO = Path(__file__).resolve().parents[2]
FLAG = re.compile(r"(?<![\w-])(-{1,2}[A-Za-z][\w-]*)")


def _script(tool: str) -> str:
    nf = REPO / "workflow" / "modules" / "local" / tool / "main.nf"
    if not nf.is_file():
        pytest.skip(f"workflow/modules/local/{tool}/main.nf not merged yet (M2-07)")
    text = re.sub(r"/\*.*?\*/", "", nf.read_text(), flags=re.S)
    text = re.sub(r"(?m)^\s*//.*$", "", text)
    m = re.search(r"\bscript\s*:(.*?)(?:\bstub\s*:|\Z)", text, flags=re.S)
    assert m, f"{tool}: no script: block"
    return m.group(1)


def _call(script: str, entry: str) -> str:
    """The command line (with backslash continuations) that invokes ``entry``."""
    joined = re.sub(r"\\\s*\n", " ", script)
    lines = [ln for ln in joined.splitlines() if entry in ln]
    assert lines, f"no line invoking {entry}"
    return lines[0]


def test_terminatornet_process_passes_exactly_the_builder_flags() -> None:
    call = _call(_script("terminatornet"), "TerminatorNet.py")
    argv = terminatornet_command(Path("genome.fa"), Path("terminatornet.tsv"), 4, Path("/h"))
    want = [a for a in argv if a.startswith("-")]
    assert want == ["-f", "-s", "-n", "-o"]
    got = [f for f in FLAG.findall(call.split("TerminatorNet.py", 1)[1]) if not f.startswith("--")]
    assert got[: len(want)] == want or sorted(set(got)) == sorted(want), call
    assert re.search(r"(?<![\w-])-s(?![\w-])", call), "TerminatorNet must run with -s (G12)"


def test_transtermhp_process_matches_the_builder() -> None:
    call = _call(_script("transtermhp"), "transterm ")
    argv = transterm_command(Path("genome.fa"), Path("genes.coords"), Path("expterm.dat"))
    assert argv[:2] == ["transterm", "-p"]
    assert re.search(r"(?<![\w-])-p\s+\S*expterm\.dat", call), call
