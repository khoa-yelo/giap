"""Tests for ``scripts/ci/check_bash32.py``: shell scripts must run on bash 3.2 (macOS CI).

Lead decision (no separate spec file; PR #255 failed on macOS over ``declare -A``):

- ``check_bash32.py [paths...]``; default ``scripts/`` and ``workflow/bin`` (when it exists).
  Scans ``*.sh``, ``*.sbatch`` and files whose shebang names bash, recursively under directories.
- Flags ``declare/local/typeset -A``, ``mapfile``/``readarray``, case modification
  (``${x,,} ${x^^} ${x,} ${x^}``), ``|&``, ``;&`` and ``;;&``, ``coproc``, ``${x@Q}``-style
  transforms, negative array indices ``${a[-1]}``, ``shopt -s globstar``.
- Ignores comment lines and trailing comments.
- ``# bash32: exempt (<reason>)`` opts a file out; exempt files are printed, exit stays 0.
- Exit 1 with ``path:line: construct`` per finding; exit 0 when clean. Standard library only.

Interface assumed (the brief leaves the name open): ``scan_text(text: str) ->
list[tuple[int, str]]`` returns ``(1-based line, construct)`` findings for one file's text,
ignoring the exemption line (exemption is a per-file decision of the CLI).
"""

from __future__ import annotations

import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
TOOL = REPO / "scripts" / "ci" / "check_bash32.py"
EXEMPT_RE = re.compile(r"^\s*#\s*bash32:\s*exempt\s*\((?P<reason>[^)]*)\)")
FINDING_RE = re.compile(r"^(?P<path>.+):(?P<line>\d+): (?P<construct>\S.*)$")


def _load() -> ModuleType:
    assert TOOL.is_file(), f"missing {TOOL.relative_to(REPO)}"
    spec = importlib.util.spec_from_file_location("check_bash32", TOOL)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(
    *args: str, cwd: Path | None = None, tool: Path = TOOL
) -> subprocess.CompletedProcess[str]:
    assert tool.is_file(), f"missing {tool}"
    return subprocess.run(
        [sys.executable, str(tool), *args], capture_output=True, text=True, cwd=cwd
    )


def _findings(stdout: str) -> list[tuple[str, int, str]]:
    out = []
    for line in stdout.splitlines():
        m = FINDING_RE.match(line)
        if m:
            out.append((m["path"], int(m["line"]), m["construct"]))
    return out


def _script(tmp_path: Path, body: str, name: str = "x.sh") -> Path:
    path = tmp_path / name
    path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + body)
    return path


# --- each construct is flagged ---------------------------------------------------------------

BAD = {
    "declare -A": "declare -A MAP=([a]=1)\n",
    "local -A": "f() {\n  local -A seen\n}\n",
    "typeset -A": "typeset -A T\n",
    "declare -gA": "declare -gA G\n",
    "mapfile": "mapfile -t LINES < file.txt\n",
    "readarray": "readarray LINES < file.txt\n",
    "lower all": 'y="${x,,}"\n',
    "upper all": 'y="${x^^}"\n',
    "lower first": 'y="${x,}"\n',
    "upper first": 'y="${x^}"\n',
    "case mod on array": 'echo "${arr[0]^^}"\n',
    "pipe stderr": "make |& tee log\n",
    "case fallthrough": "case $x in\n  a) echo a ;&\n  b) echo b ;;\nesac\n",
    "case continue": "case $x in\n  a) echo a ;;&\n  b) echo b ;;\nesac\n",
    "coproc": "coproc cat\n",
    "transform Q": 'echo "${x@Q}"\n',
    "transform U": 'echo "${x@U}"\n',
    "transform a": 'echo "${x@a}"\n',
    "negative index": 'echo "${a[-1]}"\n',
    "negative index 2": "last=${PARTS[-2]}\n",
    "globstar": "shopt -s globstar\n",
    "globstar with others": "shopt -s nullglob globstar\n",
}


@pytest.mark.parametrize("label", sorted(BAD))
def test_each_construct_is_flagged_by_the_cli(tmp_path: Path, label: str) -> None:
    path = _script(tmp_path, BAD[label])
    proc = _run(str(path))
    assert proc.returncode == 1, (label, proc.stdout, proc.stderr)
    found = _findings(proc.stdout)
    assert found, proc.stdout
    assert all(p.endswith("x.sh") for p, _, _ in found)


@pytest.mark.parametrize("label", sorted(BAD))
def test_each_construct_is_flagged_on_its_own_line(label: str) -> None:
    text = "#!/usr/bin/env bash\n" + BAD[label]
    found = _load().scan_text(text)
    assert found, label
    bad_lines = {i for i, line in enumerate(text.splitlines(), 1) if i > 1 and line.strip()}
    assert {line for line, _ in found} <= bad_lines


# --- bash 3.2 constructs that look similar are not flagged ------------------------------------

GOOD = [
    'echo "${x:-1}"\n',
    'echo "${x: -1}"\n',
    'echo "${a[i-1]}"\n',
    'echo "${a[$((n-1))]}"\n',
    'echo "${#arr[@]}"\n',
    'for v in "${arr[@]}"; do echo "$v"; done\n',
    'echo "$@" "${@}" "${*}"\n',
    'list="${CSV//,/ }"\n',
    'echo "${x/^/caret}"\n',
    "cmd > log 2>&1\n",
    "echo err >&2\n",
    "a && b || c\n",
    "a | b\n",
    "case $x in\n  a) echo a ;;\n  *) echo other ;;\nesac\n",
    "declare -a ARR=(1 2 3)\n",
    "declare -r CONST=1\n",
    "local -a items\n",
    "shopt -s nullglob\n",
    'ARGS=()\nwhile IFS= read -r line; do ARGS+=("$line"); done < f\n',
    'echo "mapfile_dir=$HOME"\n',
    "my_coprocess=1\n",
    'email="a@b.org"\n',
]


@pytest.mark.parametrize("body", GOOD)
def test_bash32_constructs_are_not_flagged(tmp_path: Path, body: str) -> None:
    path = _script(tmp_path, body)
    proc = _run(str(path))
    assert proc.returncode == 0, (body, proc.stdout, proc.stderr)
    assert _findings(proc.stdout) == []


def test_a_clean_file_exits_0_with_no_findings(tmp_path: Path) -> None:
    path = _script(
        tmp_path,
        'for f in "$@"; do\n  [ -f "$f" ] && echo "$f" >&2\ndone\nexit 0\n',
        name="clean.sbatch",
    )
    proc = _run(str(path))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _findings(proc.stdout) == []


# --- comments ----------------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(BAD))
def test_full_line_comments_are_ignored(tmp_path: Path, label: str) -> None:
    body = "".join(f"# {line}\n" for line in BAD[label].splitlines())
    body += "".join(f"    #{line}\n" for line in BAD[label].splitlines())
    path = _script(tmp_path, body)
    proc = _run(str(path))
    assert proc.returncode == 0, (label, proc.stdout)


@pytest.mark.parametrize(
    "line",
    [
        "echo ok  # no declare -A here\n",
        "x=1 # was: mapfile -t L < f\n",
        "y=2 # ${x,,} needs bash 4\n",
        "z=3 # make |& tee\n",
        "true # shopt -s globstar\n",
    ],
)
def test_trailing_comments_are_ignored(tmp_path: Path, line: str) -> None:
    proc = _run(str(_script(tmp_path, line)))
    assert proc.returncode == 0, (line, proc.stdout)


@pytest.mark.parametrize(
    "line",
    [
        "n=${#a[@]}; declare -A m\n",
        'echo "${#x}" "${x,,}"\n',
        "echo 'a#b'; mapfile -t L < f\n",
    ],
)
def test_hash_that_is_not_a_comment_does_not_hide_a_construct(tmp_path: Path, line: str) -> None:
    proc = _run(str(_script(tmp_path, line)))
    assert proc.returncode == 1, (line, proc.stdout)


# --- line numbers, multiple findings, multiple files -------------------------------------------


def test_findings_report_path_and_line(tmp_path: Path) -> None:
    path = _script(tmp_path, "echo a\ndeclare -A M\necho b\nmapfile -t L < f\n")
    proc = _run(str(path))
    assert proc.returncode == 1
    found = _findings(proc.stdout)
    assert [(Path(p).resolve(), n) for p, n, _ in found] == [
        (path.resolve(), 4),
        (path.resolve(), 6),
    ]


def test_multiple_files_are_all_reported(tmp_path: Path) -> None:
    a = _script(tmp_path, "declare -A M\n", name="a.sh")
    b = _script(tmp_path, "coproc cat\n", name="b.sbatch")
    c = _script(tmp_path, "echo fine\n", name="c.sh")
    proc = _run(str(a), str(b), str(c))
    assert proc.returncode == 1, proc.stdout
    paths = {Path(p).name for p, _, _ in _findings(proc.stdout)}
    assert paths == {"a.sh", "b.sbatch"}


def test_directories_are_scanned_recursively_by_extension_and_shebang(tmp_path: Path) -> None:
    root = tmp_path / "scripts"
    (root / "slurm" / "deep").mkdir(parents=True)
    (root / "a.sh").write_text("declare -A M\n")  # extension alone suffices
    (root / "slurm" / "job.sbatch").write_text("#!/bin/bash\nmapfile -t L < f\n")
    (root / "slurm" / "deep" / "tool").write_text("#!/usr/bin/env bash\ncoproc cat\n")
    (root / "slurm" / "deep" / "tool2").write_text("#!/bin/bash -eu\necho ${x,,}\n")
    (root / "notbash.py").write_text('CMD = "declare -A M; mapfile x"\n')
    (root / "notes.txt").write_text("declare -A M\n")
    (root / "posix").write_text("#!/bin/sh\necho hi\n")
    proc = _run(str(root))
    assert proc.returncode == 1
    names = {Path(p).name for p, _, _ in _findings(proc.stdout)}
    assert names == {"a.sh", "job.sbatch", "tool", "tool2"}


# --- exemption ----------------------------------------------------------------------------------


def test_exempt_file_exits_0_and_is_printed(tmp_path: Path) -> None:
    path = _script(
        tmp_path,
        "# bash32: exempt (runs only on Linux Slurm nodes)\ndeclare -A M\nmapfile -t L < f\n",
        name="linux_only.sbatch",
    )
    proc = _run(str(path))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _findings(proc.stdout) == []
    assert "linux_only.sbatch" in proc.stdout + proc.stderr


def test_exemption_covers_only_its_own_file(tmp_path: Path) -> None:
    ex = _script(tmp_path, "# bash32: exempt (Slurm only)\ndeclare -A M\n", name="ex.sh")
    bad = _script(tmp_path, "declare -A M\n", name="bad.sh")
    proc = _run(str(ex), str(bad))
    assert proc.returncode == 1
    assert {Path(p).name for p, _, _ in _findings(proc.stdout)} == {"bad.sh"}


@pytest.mark.parametrize(
    "marker",
    [
        "# bash32: exempt\n",
        "# bash32: exempt ()\n",
        "# bash32: exempt (   )\n",
        "echo '# bash32: exempt (in a string)'\n",
    ],
)
def test_exemption_needs_a_comment_line_with_a_reason(tmp_path: Path, marker: str) -> None:
    proc = _run(str(_script(tmp_path, marker + "declare -A M\n")))
    assert proc.returncode == 1, (marker, proc.stdout)


# --- regression: PR #255 -----------------------------------------------------------------------


def test_regression_pr255_declare_A_shards_is_flagged(tmp_path: Path) -> None:
    slurm = tmp_path / "scripts" / "slurm"
    slurm.mkdir(parents=True)
    path = slurm / "calibration_full.sbatch"
    path.write_text(
        "#!/bin/bash\n#SBATCH --job-name=giap-calibration-full\nset -euo pipefail\n"
        "declare -A SHARDS=([13]=2)\n"
        'n="${SHARDS[$SLURM_ARRAY_TASK_ID]:-1}"\n'
    )
    proc = _run(str(path))
    assert proc.returncode == 1, proc.stdout
    found = _findings(proc.stdout)
    assert [(Path(p).name, n) for p, n, _ in found] == [("calibration_full.sbatch", 4)]


# --- default paths ------------------------------------------------------------------------------


def _tree_with_tool(tmp_path: Path) -> tuple[Path, Path]:
    """A copy of the tool at <tmp>/scripts/ci/ so repo-relative and cwd-relative defaults agree."""
    tool = tmp_path / "scripts" / "ci" / "check_bash32.py"
    tool.parent.mkdir(parents=True)
    assert TOOL.is_file(), f"missing {TOOL.relative_to(REPO)}"
    shutil.copy(TOOL, tool)
    return tmp_path, tool


def test_default_paths_are_scripts_and_workflow_bin(tmp_path: Path) -> None:
    root, tool = _tree_with_tool(tmp_path)
    (root / "scripts" / "slurm").mkdir()
    (root / "scripts" / "slurm" / "a.sbatch").write_text("#!/bin/bash\ndeclare -A M\n")
    (root / "workflow" / "bin").mkdir(parents=True)
    (root / "workflow" / "bin" / "helper").write_text("#!/usr/bin/env bash\nmapfile -t L < f\n")
    (root / "elsewhere").mkdir()
    (root / "elsewhere" / "c.sh").write_text("coproc cat\n")
    proc = _run(cwd=root, tool=tool)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    names = {Path(p).name for p, _, _ in _findings(proc.stdout)}
    assert names == {"a.sbatch", "helper"}


def test_default_paths_without_workflow_bin(tmp_path: Path) -> None:
    root, tool = _tree_with_tool(tmp_path)
    (root / "scripts" / "ok.sh").write_text("#!/bin/bash\necho ok\n")
    proc = _run(cwd=root, tool=tool)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (root / "scripts" / "bad.sh").write_text("#!/bin/bash\necho ${x^^}\n")
    proc = _run(cwd=root, tool=tool)
    assert proc.returncode == 1, proc.stdout + proc.stderr


# --- properties ---------------------------------------------------------------------------------

_bad_line = st.sampled_from(
    sorted({line for body in BAD.values() for line in body.splitlines() if line.strip()})
)
_good_line = st.sampled_from(
    sorted({line for body in GOOD for line in body.splitlines() if line.strip()})
)


@settings(max_examples=60, deadline=None)
@given(st.lists(st.one_of(_bad_line, _good_line), min_size=1, max_size=12))
def test_commenting_out_every_line_leaves_no_findings(lines: list[str]) -> None:
    text = "#!/bin/bash\n" + "".join(f"# {line}\n" for line in lines)
    assert _load().scan_text(text) == []


@settings(max_examples=60, deadline=None)
@given(st.lists(_good_line, max_size=8), st.lists(_good_line, max_size=8))
def test_clean_padding_only_shifts_line_numbers(before: list[str], after: list[str]) -> None:
    mod = _load()
    base = "declare -A M\n"
    shifted = "".join(f"{x}\n" for x in before) + base + "".join(f"{x}\n" for x in after)
    found = mod.scan_text(shifted)
    assert [line for line, _ in found] == [len(before) + 1]


# --- today's repository -------------------------------------------------------------------------


def test_repo_scripts_are_bash32_clean_or_exempt_with_a_reason() -> None:
    proc = _run(cwd=REPO)
    assert proc.returncode == 0, "bash-4 constructs in repo scripts:\n" + proc.stdout
    for path in sorted((REPO / "scripts").rglob("*")):
        if not path.is_file() or path.suffix not in {".sh", ".sbatch"}:
            continue
        for line in path.read_text(errors="replace").splitlines():
            m = EXEMPT_RE.match(line)
            if m:
                assert m["reason"].strip(), f"{path.relative_to(REPO)}: exemption without reason"
