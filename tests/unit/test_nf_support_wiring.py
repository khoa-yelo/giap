"""M9-17 follow-up: the plot support outputs under Nextflow (static parse; no Nextflow needed).

M9-17 gives the ``function``, ``determinants`` and ``context_stats`` steps an optional directory
output each (``function_support``, ``determinants_support``, ``context_stats_support``;
docs/specs/M9-17-plot-support-analysis.md lead decision 5). ``gippa.runner.local.publish`` merges
their ``support/`` subdirectories into ``OUT/plots/support`` (``SUPPORT_OUTPUTS``). GIPPA_SPEC.md
§6.4 executor equality asks the same of a Nextflow run, so:

1. the subworkflow that runs each step asks its GIPPA_STEP alias for the support key (an
   ``'<key>=...'`` out spec item; else the real step never writes it) and emits a channel named
   after the key, taken from that alias's outputs;
2. main.nf passes every support emission to the plots step under the ``plots_support`` key (one
   of the ``plots`` step's optional inputs, ``gippa.steps.plots._SUPPORT_INPUTS``);
3. the Nextflow support keys equal ``SUPPORT_OUTPUTS``, and both equal the registry's ``*_support``
   optional outputs (generated, so a fourth support step is caught on both sides).

The run-time halves (an empty support channel still runs, the files land in ``OUT/plots/support``)
are in workflow/tests/{function,determinants,context,main}_support.nf.test and
tests/integration/test_support_local_vs_nextflow.py.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from gippa.runner.local import SUPPORT_OUTPUTS
from gippa.steps import load_step_modules
from gippa.steps import registry as reg

REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "workflow" / "main.nf"
LOCAL = REPO / "workflow" / "subworkflows" / "local"
REPORT_NF = LOCAL / "report.nf"

#: the M9-17 support keys and the steps that write them (lead decision 5)
EXPECTED = {
    "function_support": "function",
    "determinants_support": "determinants",
    "context_stats_support": "context_stats",
}


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _support_keys() -> dict[str, str]:
    """Registry ``*_support`` optional outputs -> the step that writes each."""
    load_step_modules()
    return {
        k: s.name for s in reg.all_steps() for k in s.optional_outputs if k.endswith("_support")
    }


KEYS = sorted(EXPECTED)


def _owner(step: str) -> Path:
    """The subworkflow file that builds the step's GIPPA_STEP inputs (``stepIn('<step>'``)."""
    hits = [
        p
        for p in sorted(LOCAL.glob("*.nf"))
        if re.search(rf"stepIn\(\s*'{step}'", _strip(p.read_text(encoding="utf-8")))
    ]
    assert len(hits) == 1, f"step '{step}': stepIn('{step}', ...) in {[p.name for p in hits]}"
    return hits[0]


def _call(text: str, start: int) -> str:
    """The text of a call from just after its ``(`` to the matching ``)``."""
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _alias_call(text: str, step: str) -> tuple[str, str]:
    """(alias, argument text) of the GIPPA_STEP alias call fed by ``def v = stepIn('<step>'``."""
    m = re.search(rf"def\s+(\w+)\s*=\s*stepIn\(\s*'{step}'", text)
    assert m, f"no `def <v> = stepIn('{step}', ...)`"
    var = m.group(1)
    c = re.search(rf"\b([A-Z][A-Z0-9_]*)\(\s*{var}\.files\s*,\s*{var}\.keys\s*,", text)
    assert c, f"no <ALIAS>({var}.files, {var}.keys, [...]) call for step '{step}'"
    return c.group(1), _call(text, c.start() + len(c.group(1)) + 1)


def _emits(path: Path) -> dict[str, str]:
    """``emit:`` name -> its right-hand side."""
    text = _strip(path.read_text(encoding="utf-8"))
    m = re.search(r"\bemit:\s*\n(.*?)\n\}", text, re.S)
    assert m, f"{path.name}: no emit block"
    out = {}
    for line in m.group(1).splitlines():
        lm = re.match(r"\s*(\w+)\s*=\s*(.+)$", line)
        if lm:
            out[lm.group(1)] = lm.group(2)
    return out


def _subworkflow_name(path: Path) -> str:
    m = re.search(r"\bworkflow\s+([A-Z][A-Z0-9_]*)\s*\{", _strip(path.read_text(encoding="utf-8")))
    assert m, f"{path.name}: no named workflow"
    return m.group(1)


def _expanded_plots_support_args() -> str:
    """Every ``kv('plots_support', ...)`` / ``outs('plots_support', ...)`` argument in main.nf,
    with ``def x = ...`` variables substituted (two levels), joined."""
    text = _strip(MAIN.read_text(encoding="utf-8"))
    defs = dict(re.findall(r"(?m)^\s*def\s+(\w+)\s*=\s*(.+)$", text))
    args = []
    for m in re.finditer(r"\b(kv|outs)\(\s*'plots_support'\s*,", text):
        arg = _call(text, m.start() + len(m.group(1)) + 1)
        for _ in range(2):
            arg = re.sub(r"\b(\w+)\b", lambda w: f"({defs[w.group(1)]})"
                         if w.group(1) in defs else w.group(1), arg)  # fmt: skip
        args.append(arg)
    return "\n".join(args)


# ---------------------------------------------------------------- the key set (generated)


def test_registry_support_keys_are_the_m9_17_three() -> None:
    assert _support_keys() == EXPECTED


def test_local_support_outputs_equal_the_registry_support_keys() -> None:
    assert set(SUPPORT_OUTPUTS) == set(_support_keys())


def test_nextflow_support_keys_equal_the_local_support_outputs() -> None:
    """Every ``'<x>_support=`` out spec item in the subworkflows is a ``SUPPORT_OUTPUTS`` key and
    every ``SUPPORT_OUTPUTS`` key has one (executor equality, GIPPA_SPEC.md §6.4)."""
    found: set[str] = set()
    for nf in sorted(LOCAL.glob("*.nf")):
        found |= set(re.findall(r"['\"]([a-z_0-9]+_support)=", _strip(nf.read_text("utf-8"))))
    assert found == set(SUPPORT_OUTPUTS), (
        f"Nextflow out specs request {sorted(found)}; the local executor publishes "
        f"{sorted(SUPPORT_OUTPUTS)}"
    )


# ---------------------------------------------------------------- 1. each subworkflow


@pytest.mark.parametrize("key", KEYS)
def test_subworkflow_asks_the_step_for_its_support_output(key: str) -> None:
    step = EXPECTED[key]
    path = _owner(step)
    alias, args = _alias_call(_strip(path.read_text(encoding="utf-8")), step)
    assert re.search(rf"['\"]{key}=", args), (
        f"{path.name}: {alias}(...) gives no '{key}=...' out spec item, so the {step} step is "
        f"never asked to write its plot support under Nextflow; out specs: {args[-300:]!r}"
    )


@pytest.mark.parametrize("key", KEYS)
def test_subworkflow_emits_its_support_key_from_the_step(key: str) -> None:
    step = EXPECTED[key]
    path = _owner(step)
    alias, _ = _alias_call(_strip(path.read_text(encoding="utf-8")), step)
    emits = _emits(path)
    assert key in emits, f"{path.name}: no emit '{key}' (emits {sorted(emits)})"
    assert re.search(rf"\b{alias}\.out\b", emits[key]), (
        f"{path.name}: emit {key} = {emits[key]!r} does not come from {alias}'s outputs"
    )


# ---------------------------------------------------------------- 2. main.nf -> plots


def test_plots_step_declares_plots_support() -> None:
    load_step_modules()
    assert "plots_support" in reg.get("plots").optional_inputs


@pytest.mark.parametrize("key", KEYS)
def test_main_nf_passes_each_support_emission_as_plots_support(key: str) -> None:
    sub = _subworkflow_name(_owner(EXPECTED[key]))
    args = _expanded_plots_support_args()
    assert args, "main.nf has no kv('plots_support', ...) / outs('plots_support', ...) channel"
    assert re.search(rf"\b{sub}\.out\.{key}\b", args), (
        f"main.nf: {sub}.out.{key} is not in the channel given to the plots step as 'plots_support'"
    )


def test_main_nf_hands_plots_support_to_report() -> None:
    """The ``plots_support`` pairs are in the ``pub`` list REPORT receives (before ``REPORT(``),
    so report.nf's stepIn('plots', [pub_pairs], ...) stages them."""
    text = _strip(MAIN.read_text(encoding="utf-8"))
    end = text.index("REPORT(pub")
    assert re.search(r"\b(?:kv|outs)\(\s*'plots_support'\s*,", text[:end]), (
        "main.nf adds no 'plots_support' pair to `pub` before REPORT(pub ...)"
    )


def test_function_support_stays_behind_the_function_module_switch() -> None:
    """FUNCTION runs only under ``'function' in params.modules``; its support emission must not
    be read outside that branch (a run with function off would reference an unrun workflow)."""
    text = _strip(MAIN.read_text(encoding="utf-8"))
    m = re.search(r"if\s*\(\s*'function'\s+in\s+\(params\.modules[^)]*\)\s*\)\s*\{", text)
    assert m, "main.nf: no `if ('function' in (params.modules ...))` branch"
    depth, i = 1, m.end()
    while depth:
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        i += 1
    outside = text[: m.start()] + text[i:]
    assert "FUNCTION.out.function_support" not in outside, (
        "main.nf reads FUNCTION.out.function_support outside the function module branch"
    )
