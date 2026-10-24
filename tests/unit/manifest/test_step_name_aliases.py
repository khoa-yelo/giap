"""``giap.runner.nextflow._step_name`` on every ``GIAP_STEP`` alias in workflow/ (#323 review S1,
S2; GIAP_SPEC.md §6.4: the trace's step names equal the steps that ran).

``GIAP_STEP`` is tagged with its step (``tag "${step}"``, the first argument of ``stepIn``), so
the trace name of ``ALIAS(x.files, ...)`` with ``def x = stepIn('<step>', ...)`` is
``<WF>:<SUB>:<ALIAS> (<step>)``. The static test below reads every alias and its ``stepIn``
step from the .nf sources and checks that ``_step_name`` gives a registered step name (and the
step itself), so the next alias that is neither ``*_STEP`` nor its step's name is caught.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from giap.runner.nextflow import _step_name
from giap.steps import load_step_modules
from giap.steps import registry as reg

REPO = Path(__file__).resolve().parents[3]
WORKFLOW = REPO / "workflow"


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _aliases() -> list[tuple[str, str, str]]:
    """``(file, ALIAS, step)`` for every ``GIAP_STEP as ALIAS`` in workflow/*.nf."""
    found: list[tuple[str, str, str]] = []
    for nf in sorted(WORKFLOW.rglob("*.nf")):
        text = _strip(nf.read_text(encoding="utf-8"))
        for alias in re.findall(r"GIAP_STEP\s+as\s+(\w+)", text):
            calls = re.findall(rf"\b{alias}\s*\(\s*(\w+)\.files\b", text)
            assert calls, f"{nf.name}: {alias} is never called as {alias}(<x>.files, ...)"
            for var in calls:
                steps = re.findall(rf"\bdef\s+{var}\s*=\s*stepIn\(\s*'(\w+)'", text)
                assert len(steps) >= 1, f"{nf.name}: no `def {var} = stepIn('<step>', ...)`"
                found.extend((nf.name, alias, s) for s in steps)
    return found


ALIASES = _aliases()


def test_the_alias_scan_finds_the_known_aliases() -> None:
    names = {a for _, a, _ in ALIASES}
    assert {"REPORT_HTML", "PLOTS", "OPERONS_STEP", "VARIANTS_STEP", "CONTEXT_STATS",
            "PARSE_BAKTA", "FUNCTION_STEP"} <= names, sorted(names)  # fmt: skip
    assert len(ALIASES) >= 30, ALIASES


@pytest.mark.parametrize(("nf", "alias", "step"), ALIASES, ids=[a for _, a, _ in ALIASES])
def test_every_giap_step_alias_maps_to_its_registered_step(nf: str, alias: str, step: str) -> None:
    load_step_modules()
    registered = {s.name for s in reg.all_steps()}
    assert step in registered, f"{nf}: stepIn('{step}') is not a registered step"
    for name in (f"GIAP:WF:{alias} ({step})", f"WF:{alias} ({step})"):
        got = _step_name(name)
        assert got in registered, f"{nf}: {name!r} -> {got!r}, not a registered step"
        assert got == step, f"{nf}: {name!r} -> {got!r}, but the task runs {step!r}"


@pytest.mark.parametrize(
    ("name", "step"),
    [
        ("REPORT:REPORT_HTML (report)", "report"),
        ("GIAP:REPORT:REPORT_HTML (report)", "report"),
        ("GIAP:REPORT:PLOTS (plots)", "plots"),
        ("GIAP:CONTEXT:CONTEXT_STATS (context_stats)", "context_stats"),
    ],
)
def test_report_html_runs_report(name: str, step: str) -> None:
    assert _step_name(name) == step
