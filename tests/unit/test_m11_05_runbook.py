"""``docs/release.md``, the release runbook (M11-05 acceptance 1-2, examples 1-2).

Spec: docs/specs/M11-05-release-runbook.md. Static: nothing is dispatched, no network. Names are
checked in context against the parsed ``release.yml``, ``pyproject.toml`` and ``CITATION.cff``.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
RUNBOOK = REPO / "docs" / "release.md"
RELEASE_YML = REPO / ".github" / "workflows" / "release.yml"
STEP = re.compile(r"^### (\d+)\. (.+)$", re.M)
HEADING = re.compile(r"^#{2,3} ", re.M)
ACTOR = re.compile(r"\*\*(agent|lead|Khoa)\*\*")
RULE = "Agents never push tags, and never dispatch `release.yml` with `dry_run=false`."
FORBIDDEN = re.compile(
    r"git tag|git push --tags|--follow-tags|git push origin v|gh release create"
    r"|dry_run\s*[=:]\s*false",
    re.I,
)
FIELDS = ("project name", "owner", "repository", "workflow name", "environment name")


@pytest.fixture(scope="module")
def text() -> str:
    return RUNBOOK.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def wf() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(RELEASE_YML.read_text(encoding="utf-8"))
    return data


def _steps(text: str) -> dict[int, tuple[str, str]]:
    """number -> (header, body up to the next ## or ### heading)."""
    out = {}
    for m in STEP.finditer(text):
        nxt = HEADING.search(text, m.end())
        out[int(m.group(1))] = (m.group(2), text[m.end() : nxt.start() if nxt else None])
    return out


def test_sections_a_to_e_and_steps_1_to_17(text: str) -> None:
    for s in "ABCDE":
        assert re.search(rf"^## {s}\. ", text, re.M), f"section {s} missing"
    assert [int(m.group(1)) for m in STEP.finditer(text)] == list(range(1, 18))


def test_every_step_names_exactly_one_actor(text: str) -> None:
    for n, (header, _) in _steps(text).items():
        assert len(ACTOR.findall(header)) == 1, f"step {n}: {header}"


def test_no_tag_or_publish_dispatch_outside_khoa_steps(text: str) -> None:
    assert RULE in text
    rest = text.replace(RULE, "")
    for header, body in _steps(text).values():
        if ACTOR.findall(header) == ["Khoa"]:
            rest = rest.replace(body, "")
    assert FORBIDDEN.findall(rest) == []


def test_step_1_publisher_fields_match_release_yml(text: str, wf: dict[str, Any]) -> None:
    body = _steps(text)[1][1]
    publish = wf["jobs"]["publish"]
    url = next(s for s in publish["steps"] if s.get("id") == "index")["run"]
    project = re.search(r"https://pypi\.org/pypi/([^/]+)/", url)
    assert (
        project
        and project.group(1)
        == tomllib.loads((REPO / "pyproject.toml").read_text())["project"]["name"]
    )
    cff = yaml.safe_load((REPO / "CITATION.cff").read_text(encoding="utf-8"))
    owner, repo = cff["repository-code"].rstrip("/").split("/")[-2:]
    targets = wf.get("on", wf.get(True))["workflow_dispatch"]["inputs"]["target"]["options"]
    assert publish["environment"] == "${{ inputs.target }}"
    want = {"pypi": "pypi", "test.pypi": "testpypi"}
    assert sorted(want.values()) == sorted(targets)
    for host, env in want.items():
        bullet = re.search(
            rf"^- \w*PyPI: <https://{re.escape(host)}\.org/[^\n]*(?:\n  [^\n]*)*", body, re.M
        )
        assert bullet, host
        got = {
            f: re.search(f.replace(" ", r"\s+") + r"\s+`([^`]+)`", bullet.group()) for f in FIELDS
        }
        assert {f: m and m.group(1) for f, m in got.items()} == dict(
            zip(FIELDS, (project.group(1), owner, repo, RELEASE_YML.name, env), strict=True)
        )


def test_dispatch_and_actor_names(text: str, wf: dict[str, Any]) -> None:
    runs = re.findall(r"gh workflow run (\S+)", text)
    assert runs and all((REPO / ".github" / "workflows" / r).is_file() for r in runs), runs
    for n in (9, 12):
        assert f"gh workflow run {RELEASE_YML.name} " in _steps(text)[n][1], n
    actor = re.search(r"github\.actor == '([^']+)'", wf["jobs"]["publish"]["if"])
    assert actor
    assert set(re.findall(r"github\.actor == '([^']+)'", text)) == {actor.group(1)}
    assert f"Required reviewers → `{actor.group(1)}`" in _steps(text)[2][1]
    for flag in wf.get("on", wf.get(True))["workflow_dispatch"]["inputs"]:
        assert f"-f {flag}=" in text, flag
