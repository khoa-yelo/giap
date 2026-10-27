"""OPERONS wiring follows the step registry (M6-08 items 6 and 7; GIPPA_SPEC.md §6.2, §9.1).

These are static checks of ``workflow/subworkflows/local/operons.nf``, ``annotate.nf`` and
``main.nf`` against ``gippa.steps`` (no Nextflow needed), in the style of
``test_nf_annotate_registry.py``. ``workflow/tests/operons.nf.test`` checks the staged keys at
run time.

* operons.nf defines ``workflow OPERONS``. It calls the ``operons`` step once, through a
  ``GIPPA_STEP`` alias, as ``ALIAS(v.files, v.keys, [...])`` with ``v = stepIn('operons', ...)``.
  The alias name is free, because ``OPERONS`` is the subworkflow's own name.
* The out spec is exactly ``pairs=tables/pairs.tsv``, ``operons=tables/operons.tsv`` and
  ``operon_flags=tables/operon_flags.tsv``. Only tables/*.tsv are published, never work Parquet
  (item 7).
* ``stepIn``'s third argument (``need``) is exactly the registry's required inputs: "this issue
  passes the real required keys explicitly".
* Every statically resolvable input key is a declared input.
* ``emit:`` has ``pairs``, ``operons`` and ``operon_flags``.
* annotate.nf emits ``layers_for_pairs`` (lead decision 1). main.nf includes OPERONS, calls it
  with ``GENETABLE.out.gene_table`` and ``ANNOTATE.out.layers_for_pairs``, and makes no
  ``CONTEXT.out`` channel an input of OPERONS. Operons runs beside CONTEXT, not after it.
"""

from __future__ import annotations

import re
from pathlib import Path

from gippa.steps import load_step_modules
from gippa.steps import registry as reg

REPO = Path(__file__).resolve().parents[2]
LOCAL = REPO / "workflow" / "subworkflows" / "local"
SUBWF = LOCAL / "operons.nf"
ANNOTATE = LOCAL / "annotate.nf"
MAIN = REPO / "workflow" / "main.nf"
OUT_SPEC = {
    "pairs": "tables/pairs.tsv",
    "operons": "tables/operons.tsv",
    "operon_flags": "tables/operon_flags.tsv",
}


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _text(path: Path) -> str:
    assert path.is_file(), f"missing {path.relative_to(REPO)}"
    return _strip(path.read_text(encoding="utf-8"))


def _spec() -> reg.StepSpec:
    load_step_modules()
    return reg.get("operons")


def _body(text: str, start: int) -> str:
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1, "[": 0, "]": 0}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _step_in(text: str) -> tuple[str, str]:
    """(variable, argument text) of the ``stepIn('operons', ...)`` call."""
    m = re.search(r"def\s+(\w+)\s*=\s*stepIn\(\s*'operons'\s*,", text)
    assert m, "operons.nf has no `def v = stepIn('operons', ...)`"
    return m.group(1), _body(text, m.end())


def _split_top(args: str) -> list[str]:
    """Top-level comma-separated arguments."""
    out, depth, cur = [], 0, ""
    for ch in args:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return [a.strip() for a in out if a.strip()]


def _list_literal(text: str, expr: str) -> list[str]:
    """Strings of a list literal, or of the variable bound to one."""
    if not expr.startswith("["):
        m = re.search(rf"def\s+{re.escape(expr)}\s*=\s*(\[[^\]]*\])", text)
        assert m, f"`{expr}` is not bound to a list literal"
        expr = m.group(1)
    return re.findall(r"'(\w+)'", expr)


def test_operons_nf_defines_the_subworkflow() -> None:
    assert re.search(r"\bworkflow\s+OPERONS\s*\{", _text(SUBWF))


def test_operons_step_called_once_through_gippa_step() -> None:
    text = _text(SUBWF)
    var, _ = _step_in(text)
    calls = re.findall(rf"\b([A-Z_]+)\(\s*{var}\.files\s*,\s*{var}\.keys\s*,", text)
    assert len(calls) == 1, f"operons step called {len(calls)}x"
    alias = calls[0]
    inc = re.search(rf"GIPPA_STEP\s+as\s+{alias}\b", text)
    assert inc, f"{alias} is not a GIPPA_STEP alias"


def test_out_spec_is_the_three_tables_only() -> None:
    text = _text(SUBWF)
    var, _ = _step_in(text)
    m = re.search(rf"\b[A-Z_]+\(\s*{var}\.files\s*,\s*{var}\.keys\s*,", text)
    assert m
    spec = _body(text, m.end())
    items = dict(re.findall(r"'(\w+)=([^'\s]+)'", spec))
    assert items == OUT_SPEC, f"out spec {items}"
    assert set(items) == set(_spec().outputs)
    assert ".parquet" not in spec


def test_need_is_exactly_the_required_inputs() -> None:
    text = _text(SUBWF)
    _, args = _step_in(text)
    parts = _split_top(args)
    assert len(parts) == 2, "stepIn('operons', chans, need): pass the required keys explicitly"
    assert set(_list_literal(text, parts[1])) == set(_spec().inputs)


def test_static_input_keys_are_declared() -> None:
    text = _text(SUBWF)
    spec = _spec()
    allowed = set(spec.inputs) | set(spec.optional_inputs)
    keys = set(re.findall(r"kv\(\s*'(\w+)'", text))
    assert keys, "no kv('<key>', ...) in operons.nf"
    assert keys <= allowed, f"undeclared keys {sorted(keys - allowed)}"


def test_emit_names() -> None:
    text = _text(SUBWF)
    m = re.search(r"\bemit:\s*(.*)\}\s*$", text, flags=re.S)
    assert m, "operons.nf has no emit: block"
    names = set(re.findall(r"(?m)^\s*(\w+)\s*=", m.group(1)))
    assert {"pairs", "operons", "operon_flags"} <= names, names


def test_annotate_emits_layers_for_pairs() -> None:
    text = _text(ANNOTATE)
    m = re.search(r"\bemit:\s*(.*)$", text, flags=re.S)
    assert m
    assert re.search(r"(?m)^\s*layers_for_pairs\s*=", m.group(1)), (
        "annotate.nf: no layers_for_pairs emit"
    )


def test_main_calls_operons_beside_context() -> None:
    text = _text(MAIN)
    assert re.search(r"include\s*\{\s*OPERONS\s*\}\s*from\s*'./subworkflows/local/operons'", text)
    m = re.search(r"\bOPERONS\s*\(", text)
    assert m, "main.nf does not call OPERONS"
    args = _body(text, m.end())
    assert "GENETABLE.out.gene_table" in args, args
    assert "ANNOTATE.out.layers_for_pairs" in args, args
    assert "CONTEXT.out" not in args, "OPERONS runs beside CONTEXT, not on its outputs"
