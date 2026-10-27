"""``terminatornet_cmd`` / ``transtermhp_cmd`` builders vs. their Nextflow processes (M2-15
acceptance "test_nf_cmd_sync covers TerminatorNet (-s present) and TransTermHP"; M2-07 sync rule,
M2-index decision 2: every ``gippa.annotate.<tool>_cmd`` builder is kept in sync with
``workflow/modules/local/<tool>/main.nf``).

Skips until the processes exist (M2-07). G12: TerminatorNet's script block must pass ``-s``; it
passes only the documented flags (its parser matches flags by substring).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from gippa.annotate.terminatornet_cmd import terminatornet_command
from gippa.annotate.transtermhp_cmd import transterm_command

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


# ``-p`` takes the builder's expterm.dat path or a shell variable resolving it in the env (M2-15
# review decision 3; e.g. ``-p $EXPTERM`` / ``-p \$EXPTERM`` / ``-p ${EXPTERM}``).
_VAR = re.compile(r"^\\?\$(\{\w+\}|\w+)$")


def test_transtermhp_process_matches_the_builder() -> None:
    script = re.sub(r"\\\s*\n", " ", _script("transtermhp"))
    calls = [ln for ln in script.splitlines() if re.match(r"\s*transterm\s", ln)]
    assert calls, "no line invoking transterm"
    toks = calls[0].split(">", 1)[0].split()
    argv = transterm_command(Path("genome.fa"), Path("genes.coords"), Path("expterm.dat"))
    assert argv[:2] == ["transterm", "-p"]
    assert toks[:2] == ["transterm", "-p"], calls[0]
    assert len(toks) == len(argv), calls[0]
    exp = toks[2]
    assert exp.endswith("expterm.dat") or _VAR.match(exp), f"-p argument {exp!r}"
    # positional order genome.fa genes.coords, as the builder
    assert [Path(t).name for t in toks[3:]] == [Path(a).name for a in argv[3:]], calls[0]
    assert [Path(a).name for a in argv[3:]] == ["genome.fa", "genes.coords"]
