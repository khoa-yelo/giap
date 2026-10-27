"""versions.yml version capture keeps stderr out of the value (#416 D1).

The deployed antiSMASH run wrote ``antismash --version 2>&1`` into versions.yml; inside the
container matplotlib printed an MPLCONFIGDIR warning on stderr, the YAML became unparsable, and
``manifest_sources._yaml`` dropped the file, so ``run_manifest.tools`` lacked antismash. The
version commands of antismash, bakta, genomad, miniprot and tblastn must discard stderr.

Two pins per module: a static one (no ``2>&1`` in the version capture, stderr sent to
``/dev/null``) and a behavioural one: the module's ``versions.yml`` heredoc is rendered with bash
against a fake tool that prints its version on stdout and a matplotlib-style warning on stderr,
and the file must parse to exactly ``{task.process: {tool: version}}``. The antiSMASH file is
also read through ``manifest_sources.tool_reports`` (the run_manifest.tools path).
"""

from __future__ import annotations

import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest
import yaml
from helpers.manifest_data import trace_row

MODULES = Path(__file__).resolve().parents[2] / "workflow" / "modules" / "local"

WARNING = (
    "Matplotlib created a temporary cache directory at /tmp/matplotlib-x1y2z3 because the "
    "default path (/home/user/.config/matplotlib) is not a writable directory; it is highly "
    "recommended to set the MPLCONFIGDIR environment variable to a writable directory, in "
    "particular to speed up the import of Matplotlib and to better support multiprocessing."
)

# tool -> (executable, stdout of its version command, expected versions.yml value)
TOOLS: dict[str, tuple[str, str, str]] = {
    "antismash": ("antismash", "antiSMASH 8.0.4\n", "8.0.4"),
    "bakta": ("bakta", "bakta 1.12.2\n", "1.12.2"),
    "genomad": ("genomad", "geNomad, version 1.12.0\n", "1.12.0"),
    "miniprot": ("miniprot", "0.18-r281\n", "0.18-r281"),
    "tblastn": (
        "tblastn",
        "tblastn: 2.16.0+\n Package: blast 2.16.0, build Jun 25 2024 00:00:00\n",
        "2.16.0+",
    ),
}


def _script_section(text: str) -> str:
    """The ``script:`` block of a module (the real command, not ``stub:``)."""
    m = re.search(r"^\s*script:\s*$(.*?)^\s*stub:\s*$", text, re.M | re.S)
    assert m, "module has no script:/stub: sections"
    return m[1]


def _heredoc(module: str) -> list[str]:
    """Lines of the script's ``versions.yml`` heredoc (Groovy source, unrendered)."""
    body = _script_section((MODULES / module / "main.nf").read_text(encoding="utf-8"))
    m = re.search(r"cat\s+<<-?\s*END_VERSIONS\s*>\s*versions\.yml\n(.*?)^\s*END_VERSIONS\s*$",
                  body, re.M | re.S)  # fmt: skip
    assert m, f"{module}: no versions.yml heredoc in script:"
    return m[1].splitlines()


def _version_line(module: str, tool: str) -> str:
    lines = [ln for ln in _heredoc(module) if re.match(rf"\s*{tool}\s*:", ln)]
    assert len(lines) == 1, f"{module}: want one '{tool}:' line in versions.yml, got {lines}"
    return lines[0]


@pytest.mark.parametrize("module", sorted(TOOLS))
def test_version_capture_does_not_merge_stderr(module: str) -> None:
    line = _version_line(module, module)
    assert "2>&1" not in line and "|&" not in line, f"{module}: stderr merged into {line!r}"
    assert re.search(r"2>\s*/dev/null", line), f"{module}: stderr not discarded in {line!r}"


def _render(module: str, process: str) -> str:
    """The heredoc as Nextflow hands it to bash: Groovy escapes resolved, indentation stripped."""
    src = "\n".join(_heredoc(module)).replace("${task.process}", process)
    src = src.replace("\\$", "$").replace("\\\\", "\\")
    assert "${" not in src.replace("$(", ""), f"{module}: unresolved Groovy interpolation: {src}"
    return f"cat <<-END_VERSIONS > versions.yml\n{textwrap.dedent(src)}\nEND_VERSIONS\n"


def _fake_tool(bindir: Path, exe: str, stdout: str) -> None:
    out = stdout.replace("'", "'\\''")
    warn = WARNING.replace("'", "'\\''")
    script = f"#!/bin/sh\nprintf '%s\\n' '{warn}' >&2\nprintf '%s' '{out}'\n"
    path = bindir / exe
    path.write_text(script)
    path.chmod(0o755)


def _run_heredoc(tmp_path: Path, module: str, process: str) -> Path:
    exe, stdout, _ = TOOLS[module]
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    _fake_tool(bindir, exe, stdout)
    bash = shutil.which("bash")
    assert bash
    work = tmp_path / "task"
    work.mkdir(exist_ok=True)
    subprocess.run([bash, "-c", _render(module, process)], cwd=work, check=True,
                   env={"PATH": f"{bindir}:/usr/bin:/bin", "LC_ALL": "C"})  # fmt: skip
    return work


@pytest.mark.parametrize("module", sorted(TOOLS))
def test_versions_yml_parses_with_a_stderr_warning(tmp_path: Path, module: str) -> None:
    process = f"GIAP:ANNOTATE:{module.upper()}"
    work = _run_heredoc(tmp_path, module, process)
    text = (work / "versions.yml").read_text(encoding="utf-8")
    assert "Matplotlib" not in text, text
    doc = yaml.load(text, Loader=yaml.BaseLoader)
    assert doc == {process: {module: TOOLS[module][2]}}, text


def test_antismash_reaches_tool_reports(tmp_path: Path) -> None:
    from giap.io import manifest_sources

    work = _run_heredoc(tmp_path, "antismash", "GIAP:ANNOTATE:ANTISMASH")
    row = trace_row("GIAP:ANNOTATE:ANTISMASH (genome)", workdir=str(work), container="-")
    reps = manifest_sources.tool_reports([row], tmp_path / "work", tmp_path / "locks")
    assert [(r["name"], r["version"]) for r in reps] == [("antismash", "8.0.4")]
