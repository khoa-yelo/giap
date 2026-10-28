"""M11-19 item 1: the docker profile's ``docker.runOptions`` starts with ``--pull=never``.

docs/specs/M11-19-nextflow-docker-pull-never.md item 1, acceptance criteria 7 (second half: text
parse, no Nextflow) and 9. ``workflow/nextflow.config`` is read as text: the ``docker { }`` and
``apptainer { }`` bodies of ``profiles { }`` are cut out (comments stripped, quoted braces
skipped).

- The docker body has exactly one ``docker.runOptions`` statement, on one line; its expression
  starts with a string literal whose first whitespace-separated token is ``--pull=never``
  (equals form), and the line holds ``--pull=never`` exactly once and no other ``--pull``.
- The apptainer body and ``conf/container_common.config`` contain no ``--pull``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
NEXTFLOW = WF / "nextflow.config"
COMMON = WF / "conf" / "container_common.config"


def _strip_comments(text: str) -> str:
    out: list[str] = []
    i, quote = 0, ""
    while i < len(text):
        c = text[i]
        if quote:
            out.append(c)
            if c == quote and text[i - 1] != "\\":
                quote = ""
        elif c in "'\"":
            out.append(c)
            quote = c
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = len(text) if j < 0 else j
            continue
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = len(text) if j < 0 else j + 2
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def _body(text: str, start: int) -> str:
    """The text between the ``{`` at ``start`` and its matching ``}``."""
    depth, quote = 0, ""
    for i in range(start, len(text)):
        c = text[i]
        if quote:
            if c == quote and text[i - 1] != "\\":
                quote = ""
        elif c in "'\"":
            quote = c
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i]
    raise AssertionError("unbalanced braces")


def _profile(name: str) -> str:
    text = _strip_comments(NEXTFLOW.read_text(encoding="utf-8"))
    m = re.search(r"^profiles\s*\{", text, re.M)
    assert m, "no profiles block"
    profiles = _body(text, m.end() - 1)
    p = re.search(rf"(?m)^\s*{name}\s*\{{", profiles)
    assert p, f"no {name} profile"
    return _body(profiles, p.end() - 1)


def _run_options_lines() -> list[str]:
    return [ln for ln in _profile("docker").splitlines() if re.search(r"\brunOptions\s*=", ln)]


def test_one_run_options_statement_on_one_line() -> None:
    lines = _run_options_lines()
    assert len(lines) == 1, lines
    line = lines[0].rstrip()
    assert re.match(r"^\s*docker\.runOptions\s*=", line), line
    assert not line.endswith("+") and line.count("(") == line.count(")"), (
        "the expression must stay on one line (test_m11_14_profile_parity reads it from one line)"
    )


def test_pull_never_is_the_first_token() -> None:
    line = _run_options_lines()[0]
    rhs = line.split("=", 1)[1].strip()
    m = re.match(r"""^(['"])(.*?)\1""", rhs)
    assert m, f"runOptions must start with a string literal holding --pull=never: {rhs}"
    tokens = m.group(2).split()
    assert tokens and tokens[0] == "--pull=never", f"first token is not --pull=never: {rhs}"
    assert m.group(2).startswith("--pull=never "), "a space must separate it from what follows"


def test_pull_never_once_and_no_other_pull() -> None:
    line = _run_options_lines()[0]
    assert line.count("--pull=never") == 1, line
    assert line.count("--pull") == 1, f"another --pull in runOptions: {line}"
    assert not re.search(r"--pull\s+never", line), "use the equals form --pull=never"


def test_the_rest_of_run_options_is_kept() -> None:
    """The DB bind and ``-u`` stay (14p), after ``--pull=never``."""
    line = _run_options_lines()[0]
    i = line.index("--pull=never")
    assert line.index("${params.db_dir}:${params.db_dir}:ro") > i, line
    assert line.index("-u $(id -u):$(id -g)") > i, line


def test_docker_profile_has_no_other_pull() -> None:
    body = _profile("docker")
    assert body.count("--pull") == 1, "only docker.runOptions carries --pull"


def test_ac9_apptainer_and_common_config_have_no_pull() -> None:
    """Regression guard (AC 9): ``--pull`` is docker syntax only."""
    assert "--pull" not in _profile("apptainer")
    assert "--pull" not in _strip_comments(COMMON.read_text(encoding="utf-8"))
