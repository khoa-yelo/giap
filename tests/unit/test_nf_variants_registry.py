"""VARIANTS wiring follows the step registry (docs/specs/M7-04-variants-step-and-run.md rules 3
and 5, acceptance "a registry test that the step's declared keys match the subworkflow's"; lead
decision 3; GIPPA_SPEC.md §6.2, §6.4, §9.1).

Static checks of ``workflow/subworkflows/local/variants.nf``, ``mapping.nf`` and ``main.nf``
against ``gippa.steps`` (no Nextflow needed), in the style of ``test_nf_annotate_registry.py``.
``workflow/tests/variants.nf.test`` and ``main_variants.nf.test`` check the staged keys at run time.

* variants.nf defines ``workflow VARIANTS``. It calls the ``variants`` step once, through a
  ``GIPPA_STEP`` alias, as ``ALIAS(v.files, v.keys, [...])`` with ``v = stepIn('variants', ...)``.
  The spec writes ``GIPPA_STEP as VARIANTS``, but ``VARIANTS`` is the subworkflow's own name, so
  the alias name is left free (as OPERONS uses ``OPERONS_STEP``; the nf-tests find the task by its
  tag ``variants``).
* The out spec is exactly ``pileup=tables/pileup.tsv.gz``, ``variants=tables/variants.tsv`` and
  ``variants_vcf=tables/variants.vcf``. Only ``tables/*`` is published, never work Parquet. These
  names equal ``gippa.runner.local.PUBLISHED_FILES`` (lead decision 3).
* ``stepIn``'s third argument (``need``) is exactly the registry's required inputs (rule 5).
* Every statically resolvable input key is a declared input.
* ``take:`` names the six channels of rule 5, in order. ``emit:`` has ``pileup``, ``variants``
  and ``variants_vcf``.
* mapping.nf emits ``hits_tolerant`` (picked from MAP_TOLERANT's ``hits_tolerant.parquet``).
* main.nf includes VARIANTS and calls it inside the ``if (inputs.peptides)`` block, after
  MAPPING, with ``MAPPING.out.mapping`` and ``MAPPING.out.hits_tolerant``. No GENETABLE or
  CONTEXT output is an input of VARIANTS.
"""

from __future__ import annotations

import re
from pathlib import Path

from gippa.steps import load_step_modules
from gippa.steps import registry as reg

REPO = Path(__file__).resolve().parents[2]
LOCAL = REPO / "workflow" / "subworkflows" / "local"
SUBWF = LOCAL / "variants.nf"
MAPPING = LOCAL / "mapping.nf"
MAIN = REPO / "workflow" / "main.nf"
OUT_SPEC = {
    "pileup": "tables/pileup.tsv.gz",
    "variants": "tables/variants.tsv",
    "variants_vcf": "tables/variants.vcf",
    "variants_qc": "work/variants_qc.parquet",  # M7-14
}
TAKE = ["genome", "gene_model", "mapping", "hits_tolerant", "peptides", "reference_peptides",
        "peptide_spectra"]  # fmt: skip


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _text(path: Path) -> str:
    assert path.is_file(), f"missing {path.relative_to(REPO)}"
    return _strip(path.read_text(encoding="utf-8"))


def _spec() -> reg.StepSpec:
    load_step_modules()
    return reg.get("variants")


def _body(text: str, start: int, open_: str = "(", close: str = ")") -> str:
    depth, i = 1, start
    while depth:
        depth += {open_: 1, close: -1}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _step_in(text: str) -> tuple[str, str]:
    m = re.search(r"def\s+(\w+)\s*=\s*stepIn\(\s*'variants'\s*,", text)
    assert m, "variants.nf has no `def v = stepIn('variants', ...)`"
    return m.group(1), _body(text, m.end())


def _split_top(args: str) -> list[str]:
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
    if not expr.startswith("["):
        m = re.search(rf"def\s+{re.escape(expr)}\s*=\s*(\[[^\]]*\])", text)
        assert m, f"`{expr}` is not bound to a list literal"
        expr = m.group(1)
    return re.findall(r"'(\w+)'", expr)


def _call(text: str) -> tuple[str, str]:
    """(alias, out spec text) of the single ``ALIAS(v.files, v.keys, [...])`` call."""
    var, _ = _step_in(text)
    calls = list(re.finditer(rf"\b([A-Z_]+)\(\s*{var}\.files\s*,\s*{var}\.keys\s*,", text))
    assert len(calls) == 1, f"variants step called {len(calls)}x"
    return calls[0].group(1), _body(text, calls[0].end())


def _workflow_block(text: str, name: str) -> str:
    m = re.search(rf"\bworkflow\s+{name}\s*\{{", text)
    assert m, f"no `workflow {name} {{`"
    return _body(text, m.end(), "{", "}")


def _section(block: str, label: str) -> str:
    m = re.search(rf"\b{label}:\s*(.*?)(?:\n\s*(?:take|main|emit):|\Z)", block, flags=re.S)
    assert m, f"no {label}: section"
    return m.group(1)


def test_variants_nf_defines_the_subworkflow() -> None:
    _workflow_block(_text(SUBWF), "VARIANTS")


def test_variants_step_called_once_through_a_gippa_step_alias() -> None:
    text = _text(SUBWF)
    alias, _ = _call(text)
    assert re.search(rf"GIPPA_STEP\s+as\s+{alias}\b", text), f"{alias} is not a GIPPA_STEP alias"


def test_out_spec_is_the_three_tables_and_variants_qc() -> None:
    _, spec = _call(_text(SUBWF))
    items = dict(re.findall(r"'(\w+)=([^'\s]+)", spec))
    assert items == OUT_SPEC, f"out spec {items}"
    assert set(items) == set(_spec().outputs)
    # the only Parquet is the unpublished work intermediate (M7-14)
    assert [k for k, v in items.items() if v.endswith(".parquet")] == ["variants_qc"]


def test_out_spec_equals_the_local_published_files() -> None:
    """Lead decision 3: the local executor publishes the same names as Nextflow."""
    from gippa.runner.local import PUBLISHED_FILES

    for key, rel in OUT_SPEC.items():
        assert PUBLISHED_FILES.get(key) == rel, (key, PUBLISHED_FILES.get(key))


def test_need_is_exactly_the_required_inputs() -> None:
    text = _text(SUBWF)
    _, args = _step_in(text)
    parts = _split_top(args)
    assert len(parts) == 2, "stepIn('variants', chans, need): pass the required keys explicitly"
    assert set(_list_literal(text, parts[1])) == set(_spec().inputs)


def test_static_input_keys_are_declared() -> None:
    text = _text(SUBWF)
    spec = _spec()
    allowed = set(spec.inputs) | set(spec.optional_inputs)
    keys = set(re.findall(r"kv\(\s*'(\w+)'", text))
    assert keys, "no kv('<key>', ...) in variants.nf"
    assert keys <= allowed, f"undeclared keys {sorted(keys - allowed)}"
    assert set(spec.inputs) <= keys, f"required keys not wired: {sorted(set(spec.inputs) - keys)}"
    assert "hits_tolerant" in keys, "the optional hits_tolerant input is not wired"


def test_take_channels_in_rule_5_order() -> None:
    block = _workflow_block(_text(SUBWF), "VARIANTS")
    names = re.findall(r"(?m)^\s*(\w+)\s*$", _section(block, "take"))
    assert names == TAKE, names


def test_emit_names() -> None:
    block = _workflow_block(_text(SUBWF), "VARIANTS")
    names = set(re.findall(r"(?m)^\s*(\w+)\s*=", _section(block, "emit")))
    assert {"pileup", "variants", "variants_vcf"} <= names, names


def test_mapping_emits_hits_tolerant() -> None:
    block = _workflow_block(_text(MAPPING), "MAPPING")
    emit = _section(block, "emit")
    m = re.search(r"(?m)^\s*hits_tolerant\s*=\s*(.+)$", emit)
    assert m, f"mapping.nf: no hits_tolerant emit in {emit!r}"
    assert "hits_tolerant.parquet" in m.group(1), m.group(1)


def test_main_calls_variants_after_mapping_inside_the_peptide_block() -> None:
    text = _text(MAIN)
    assert re.search(r"include\s*\{\s*VARIANTS\s*\}\s*from\s*'./subworkflows/local/variants'", text)
    m = re.search(r"if\s*\(\s*inputs\.peptides\s*\)\s*\{", text)
    assert m, "main.nf has no `if (inputs.peptides) {` block"
    block = _body(text, m.end(), "{", "}")
    call = re.search(r"\bVARIANTS\s*\(", block)
    assert call, "VARIANTS is not called inside the peptide block"
    assert block.index("MAPPING(") < call.start(), "VARIANTS must follow MAPPING"
    args = _body(block, call.end())
    assert "MAPPING.out.mapping" in args, args
    assert "MAPPING.out.hits_tolerant" in args, args
    assert "GENETABLE.out" not in args and "CONTEXT.out" not in args, (
        "VARIANTS is independent of GENETABLE and CONTEXT"
    )
    outside = text.replace(block, "")
    assert not re.search(r"\bVARIANTS\s*\(", outside), "VARIANTS called outside the peptide block"


def test_nf_test_headers_equal_the_python_columns() -> None:
    """The headers that workflow/tests/variants.nf.test expects from the stub equal the §9.3
    columns, so the stub publishes what a real run would."""
    from gippa.variants.pileup import PILEUP_COLUMNS
    from gippa.variants.vcf_records import INFO_FIELDS, VARIANTS_COLUMNS

    text = (REPO / "workflow" / "tests" / "variants.nf.test").read_text(encoding="utf-8")
    for name, cols in (("PILEUP_HEADER", PILEUP_COLUMNS), ("VARIANTS_HEADER", VARIANTS_COLUMNS)):
        m = re.search(rf"def {name} = '([^']*)'", text)
        assert m, f"variants.nf.test pins no {name}"
        assert m.group(1).split(" ") == list(cols), name
    m = re.search(r"def INFO_IDS = \[([^\]]*)\]", text)
    assert m
    assert re.findall(r"'(\w+)'", m.group(1)) == [f.id for f in INFO_FIELDS]
