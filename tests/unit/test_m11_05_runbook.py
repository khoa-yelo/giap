"""``docs/release.md``, the release runbook (M11-05 acceptance 1-2, examples 1-2).

Spec: docs/specs/M11-05-release-runbook.md. Static: nothing is dispatched, no network.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
RUNBOOK = REPO / "docs" / "release.md"
STEP = re.compile(r"^### (\d+)\. (.+)$", re.M)
ACTOR = re.compile(r"\*\*(agent|lead|Khoa)\*\*")
FORBIDDEN = ("git tag", "git push --tags", "--follow-tags", "dry_run=false")


@pytest.fixture(scope="module")
def text() -> str:
    return RUNBOOK.read_text(encoding="utf-8")


def _steps(text: str) -> list[tuple[int, str, str]]:
    """(number, header, body up to the next ## or ### heading) per numbered step."""
    out = []
    for m in STEP.finditer(text):
        nxt = re.compile(r"^#{2,3} ", re.M).search(text, m.end())
        out.append((int(m.group(1)), m.group(2), text[m.end() : nxt.start() if nxt else None]))
    return out


def test_sections_a_to_e_and_steps_1_to_17(text: str) -> None:
    for s in "ABCDE":
        assert re.search(rf"^## {s}\. ", text, re.M), f"section {s} missing"
    assert [n for n, _, _ in _steps(text)] == list(range(1, 18))


def test_every_step_names_exactly_one_actor(text: str) -> None:
    for n, header, _ in _steps(text):
        assert len(ACTOR.findall(header)) == 1, f"step {n}: {header}"


def test_no_tag_or_publish_dispatch_outside_khoa_steps(text: str) -> None:
    for n, header, body in _steps(text):
        if ACTOR.findall(header) != ["Khoa"]:
            bad = [f for f in FORBIDDEN if f in header + body]
            assert bad == [], f"step {n} ({header}) contains {bad}"
    assert "never push tags" in text and "never dispatch `release.yml` with `dry_run=false`" in text


def test_names_match_release_yml(text: str) -> None:
    wf = yaml.safe_load((REPO / ".github/workflows/release.yml").read_text(encoding="utf-8"))
    inputs = wf.get("on", wf.get(True))["workflow_dispatch"]["inputs"]
    envs = inputs["target"]["options"]
    assert sorted(envs) == ["pypi", "testpypi"]
    assert wf["jobs"]["publish"]["environment"] == "${{ inputs.target }}"
    for name in ("release.yml", "`giap`", *(f"`{e}`" for e in envs)):
        assert name in text, name
    for flag in inputs:
        assert f"-f {flag}=" in text, flag
    assert "https://pypi.org/pypi/giap/" in wf["jobs"]["publish"]["steps"][1]["run"]
