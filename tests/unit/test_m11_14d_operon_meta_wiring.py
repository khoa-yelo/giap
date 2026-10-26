"""M11-14d, ADR-0061: a static guard on the Nextflow wiring of ``operon_meta``.

``operon_meta`` (operon_pairs' ``operon_meta.json``) is optional for its consumers, so a channel
that stops carrying it fails open: the ``tu_source_distance_model`` flags just vanish. This test
reads workflow/main.nf and workflow/subworkflows/local/context.nf as text: every registry step that
takes ``operon_meta`` must be one whose Nextflow call is checked here, and each check finds the
key in that call. The stub nf-test workflow/tests/main_operons_fallback.nf.test checks the staged
inputs at run time.
"""

from __future__ import annotations

import re
from pathlib import Path

from giap.steps import load_step_modules
from giap.steps.registry import all_steps

REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "workflow" / "main.nf"
CONTEXT = REPO / "workflow" / "subworkflows" / "local" / "context.nf"
KV = re.compile(r"kv\(\s*'operon_meta'\s*,\s*ANNOTATE\.out\.operon_meta\s*\)")


def _call(text: str, name: str) -> str:
    """The argument text of the one top-level ``NAME(...)`` call (balanced parentheses)."""
    hits = [m.end() for m in re.finditer(rf"(?m)^\s*{name}\(", text)]
    assert len(hits) == 1, f"{name}( called {len(hits)}x"
    depth, i = 1, hits[0]
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[hits[0] : i - 1]


def _consumers() -> set[str]:
    load_step_modules()
    return {s.name for s in all_steps() if "operon_meta" in (*s.optional_inputs, *s.inputs)}


def test_every_operon_meta_consumer_is_checked_here() -> None:
    assert _consumers() == {"operons", "context_stats"}


def test_main_nf_passes_operon_meta_to_operons() -> None:
    assert KV.search(_call(MAIN.read_text(encoding="utf-8"), "OPERONS"))


def test_main_nf_passes_operon_meta_to_context() -> None:
    assert KV.search(_call(MAIN.read_text(encoding="utf-8"), "CONTEXT"))


def test_context_nf_hands_operon_meta_to_context_stats() -> None:
    text = CONTEXT.read_text(encoding="utf-8")
    keep = r"layers_for_pairs\.filter\s*\{\s*k\s*,\s*f\s*->\s*k\s*==\s*'operon_meta'\s*\}"
    assert re.search(keep, text), "context.nf drops the operon_meta entry of layers_for_pairs"
    m = re.search(r"stepIn\('context_stats',\s*([^)]*)\)", text)
    assert m and "optional" in m.group(1), "context_stats' stepIn lost the optional list"
