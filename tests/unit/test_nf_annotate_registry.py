"""ANNOTATE py step calls follow the step registry (M2-07 decision 1, PR #65 review).

Every ``GIAP_STEP`` alias call in ``workflow/subworkflows/local/annotate.nf`` whose step is
registered in ``giap.steps`` (``load_step_modules``) must declare exactly the registry's output
keys. Input keys depend on the channel wiring, so they are checked at run time by
``workflow/tests/annotate.nf.test`` (every ``.giap_inputs`` against
``workflow/tests/data/annotate/step_registry.json``); here that fixture must equal the live
registry, so the nf-test always checks the current contract. No Nextflow needed.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from giap.steps import load_step_modules
from giap.steps.registry import all_steps

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "annotate.nf"
FIXTURE = REPO / "workflow" / "tests" / "data" / "annotate" / "step_registry.json"

# Steps merged before PR #65 (M2-02, M2-12): the check below must cover at least these.
MERGED = {"parse_bakta", "uniop_input", "operon_pairs"}


def _registry() -> dict[str, dict[str, list[str]]]:
    load_step_modules()
    return {
        s.name: {
            "inputs": sorted(s.inputs),
            "optional_inputs": sorted(s.optional_inputs),
            "outputs": sorted(s.outputs),
        }
        for s in all_steps()
        if s.kind == "py"
    }


def _strip_comments(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _step_calls() -> dict[str, list[str]]:
    """step -> output keys of its call ``ALIAS(v.files, v.keys, [...] (+ scatterSpec)*)``."""
    text = _strip_comments(SUBWF.read_text(encoding="utf-8"))
    var_step = dict(re.findall(r"def\s+(\w+)\s*=\s*stepIn\(\s*'(\w+)'", text))
    out: dict[str, list[str]] = {}
    for m in re.finditer(r"\b([A-Z_]+)\(\s*(\w+)\.files\s*,\s*\2\.keys\s*,", text):
        alias, var = m.group(1), m.group(2)
        assert var in var_step, f"{alias}: {var} is not a stepIn(...) result"
        step = var_step[var]
        assert alias == step.upper(), f"alias {alias} runs step {step}"
        # the out spec runs to the closing parenthesis of the call
        depth, i = 1, m.end()
        while depth:
            depth += {"(": 1, ")": -1}.get(text[i], 0)
            i += 1
        spec = text[m.end() : i - 1]
        keys = re.findall(r"'(\w+)=", spec) + re.findall(r"scatterSpec\(\s*'(\w+)'", spec)
        assert step not in out, f"step {step} called twice"
        out[step] = keys
    return out


def test_step_registry_fixture_equals_the_live_registry() -> None:
    want = _registry()
    got = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert got == want, (
        "regenerate workflow/tests/data/annotate/step_registry.json from giap.steps "
        "(py steps: sorted inputs, optional_inputs, outputs)"
    )


def test_annotate_calls_every_merged_registered_step() -> None:
    called = set(_step_calls())
    assert set(_registry()) >= MERGED
    assert called >= MERGED, f"annotate.nf does not call {sorted(MERGED - called)}"


@pytest.mark.parametrize("step", sorted(_step_calls()))
def test_out_keys_equal_the_registry_outputs(step: str) -> None:
    reg = _registry()
    if step not in reg:
        pytest.skip(f"{step} not registered yet (its PR adds it)")
    got = _step_calls()[step]
    assert len(got) == len(set(got)), f"{step}: duplicate output key in {got}"
    assert sorted(got) == reg[step]["outputs"], f"{step}: annotate.nf outputs {sorted(got)}"
