"""MAPPING py step calls follow the step registry (M3-17 lead decision 5: the subworkflow follows
M2-07, ``GIAP_STEP`` + ``stepIn``, no duplicate processes; decision 6: ``placement`` takes
genes/features/replicons; ``contaminants`` is wired into ``loci``).

Static checks of ``workflow/subworkflows/local/mapping.nf`` against ``giap.steps.mapping`` (no
Nextflow needed), like tests/unit/test_nf_annotate_registry.py does for ANNOTATE:

* every py step of ``giap.steps.mapping`` is called exactly once, through a ``GIAP_STEP`` alias
  named after the step (``ALIAS(v.files, v.keys, [out spec])`` with ``v = stepIn('<step>', ...)``),
  and its out-spec keys equal the registry's outputs;
* every statically resolvable ``stepIn`` input key (``kv('<key>', ...)`` literals, variables bound
  to them, and list variables of those) is a declared input of that step;
* the ``giap_mapping`` module (GIAP_MAPPING / GIAP_MAPPING_LOCI duplicates of GIAP_STEP) is gone;
* workflow/tests/data/mapping/step_registry.json (read by workflow/tests/mapping.nf.test, which
  checks every task's ``.giap_inputs`` at run time) equals the live registry;
* the §9.3 headers pinned in mapping.nf.test equal the Python column constants.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from giap.steps import load_step_modules
from giap.steps.registry import all_steps

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "mapping.nf"
NFTEST = REPO / "workflow" / "tests" / "mapping.nf.test"
FIXTURE = REPO / "workflow" / "tests" / "data" / "mapping" / "step_registry.json"
GIAP_MAPPING_MODULE = REPO / "workflow" / "modules" / "local" / "giap_mapping"

MERGED = {"sixframe_index", "map_exact", "map_tolerant", "placement", "protein_leftovers", "loci"}


def _registry() -> dict[str, dict[str, list[str]]]:
    load_step_modules()
    return {
        s.name: {
            "inputs": sorted(s.inputs),
            "optional_inputs": sorted(s.optional_inputs),
            "outputs": sorted(s.outputs),
        }
        for s in all_steps()
        if s.kind == "py" and getattr(s.func, "__module__", "") == "giap.steps.mapping"
    }


def _text() -> str:
    return re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))


def _call_body(text: str, start: int) -> str:
    """Text from ``start`` (just after an opening parenthesis) to its closing parenthesis."""
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _step_vars(text: str) -> dict[str, str]:
    """``def v = stepIn('<step>', ...)`` → {v: step}."""
    return dict(re.findall(r"def\s+(\w+)\s*=\s*stepIn\(\s*'(\w+)'", text))


def _step_calls() -> dict[str, tuple[str, list[str]]]:
    """step → (alias, out keys) of ``ALIAS(v.files, v.keys, [...])``."""
    text = _text()
    var_step = _step_vars(text)
    out: dict[str, tuple[str, list[str]]] = {}
    for m in re.finditer(r"\b([A-Z_]+)\(\s*(\w+)\.files\s*,\s*\2\.keys\s*,", text):
        alias, var = m.group(1), m.group(2)
        assert var in var_step, f"{alias}: {var} is not a stepIn(...) result"
        step = var_step[var]
        assert step not in out, f"step {step} called twice"
        spec = _call_body(text, m.end())
        keys = re.findall(r"'(\w+)=", spec) + re.findall(r"scatterSpec\(\s*'(\w+)'", spec)
        out[step] = (alias, keys)
    return out


def _gs_aliases(text: str) -> set[str]:
    """Aliases of ``GIAP_STEP`` included from modules/local/giap_step/main."""
    aliases: set[str] = set()
    for m in re.finditer(r"include\s*\{(.*?)\}\s*from\s*'([^']+)'", text, flags=re.S):
        if m.group(2).rstrip("/").endswith("giap_step/main"):
            aliases |= set(re.findall(r"GIAP_STEP\s+as\s+(\w+)", m.group(1)))
    return aliases


def _resolver(text: str) -> dict[str, set[str]]:
    """Variable → input keys it stands for (``kv('k', ..)`` vars and lists of them)."""
    keys: dict[str, set[str]] = {
        v: {k} for v, k in re.findall(r"def\s+(\w+)\s*=\s*kv\(\s*'(\w+)'", text)
    }
    lists = {
        m.group(1): text[m.end() : text.index("\n", m.end())]
        for m in re.finditer(r"def\s+(\w+)\s*=\s*\[", text)
    }
    for _ in range(5):  # nested list vars
        for v, body in lists.items():
            got = set(re.findall(r"kv\(\s*'(\w+)'", body))
            for tok in re.findall(r"\b(\w+)\b", body):
                got |= keys.get(tok, set())
            if got:
                keys[v] = keys.get(v, set()) | got
    return keys


def _step_in_keys(text: str, step: str) -> set[str]:
    m = re.search(rf"stepIn\(\s*'{step}'\s*,", text)
    assert m, f"no stepIn('{step}', ...) call in mapping.nf"
    body = _call_body(text, m.end())
    res = _resolver(text)
    keys = set(re.findall(r"kv\(\s*'(\w+)'", body))
    for tok in re.findall(r"\b(\w+)\b(?!\s*\[)", body):  # an indexed list var is indeterminate
        keys |= res.get(tok, set())
    return keys


# ---------------------------------------------------------------- registry fixture


def test_registry_has_the_merged_mapping_steps() -> None:
    assert set(_registry()) >= MERGED, sorted(_registry())


def test_step_registry_fixture_equals_the_live_registry() -> None:
    got = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert got == _registry(), (
        "regenerate workflow/tests/data/mapping/step_registry.json from giap.steps.mapping "
        "(py steps: sorted inputs, optional_inputs, outputs)"
    )


# ---------------------------------------------------------------- no duplicate processes


def test_giap_mapping_module_is_removed() -> None:
    assert not GIAP_MAPPING_MODULE.exists(), (
        "decision 5: no duplicate processes; py steps run through GIAP_STEP, so "
        "workflow/modules/local/giap_mapping/ must go"
    )
    text = _text()
    assert "giap_mapping" not in text, "mapping.nf still includes the giap_mapping module"
    assert not re.search(r"\bGIAP_MAPPING(_LOCI)?\b", text), "mapping.nf uses GIAP_MAPPING*"


def test_every_stepin_step_is_a_registered_mapping_step() -> None:
    text = _text()
    called = set(_step_vars(text).values())
    unknown = called - set(_registry())
    assert not unknown, f"stepIn steps not in giap.steps.mapping: {sorted(unknown)}"
    assert called == set(_registry()), (
        f"mapping.nf stepIn steps {sorted(called)} != registry {sorted(_registry())}"
    )


def test_every_mapping_step_runs_once_through_a_giap_step_alias() -> None:
    text = _text()
    calls = _step_calls()
    aliases = _gs_aliases(text)
    assert set(calls) == set(_registry()), f"called: {sorted(calls)}"
    for step, (alias, _) in calls.items():
        assert alias == step.upper(), f"alias {alias} runs step {step}"
        assert alias in aliases, f"{alias} is not `GIAP_STEP as {alias}` from giap_step/main"


@pytest.mark.parametrize("step", sorted(MERGED))
def test_out_keys_equal_the_registry_outputs(step: str) -> None:
    calls = _step_calls()
    assert step in calls, f"mapping.nf does not call {step} through GIAP_STEP"
    got = calls[step][1]  # a directory output may list one entry per stub file (scatterSpec)
    assert sorted(set(got)) == _registry()[step]["outputs"], f"{step}: mapping.nf outputs {got}"


@pytest.mark.parametrize("step", sorted(MERGED))
def test_static_stepin_keys_are_declared_inputs(step: str) -> None:
    reg = _registry()[step]
    keys = _step_in_keys(_text(), step)
    extra = keys - set(reg["inputs"]) - set(reg["optional_inputs"])
    assert not extra, f"stepIn('{step}') passes undeclared key(s) {sorted(extra)}"


# ---------------------------------------------------------------- decisions 5, 6


def test_placement_gets_genes_features_and_replicons() -> None:
    """Decision 6 / deviation b: placement resolves ``search_db: auto`` from the gene model."""
    keys = _step_in_keys(_text(), "placement")
    assert keys >= {"genes", "features", "replicons"}, f"placement inputs: {sorted(keys)}"


def test_loci_gets_the_contaminant_library() -> None:
    """Decision 5 (review finding 7): ``contaminants`` (optional) is wired into loci."""
    text = _text()
    assert "contaminants" in _step_in_keys(text, "loci"), "loci stepIn lacks kv('contaminants', ..)"
    take = re.search(r"take:(.*?)main:", text, flags=re.S)
    assert take and re.search(r"\bcontaminants\b", take.group(1)), (
        "MAPPING must take a `contaminants` channel (library FASTA, or channel.empty())"
    )


# ---------------------------------------------------------------- stub headers


def _nftest_header(name: str) -> list[str]:
    m = re.search(rf'def {name} = "([^"]+)"', NFTEST.read_text(encoding="utf-8"))
    assert m, f"mapping.nf.test does not pin {name}"
    return m.group(1).split("\\t")


def test_nf_test_headers_equal_the_python_columns() -> None:
    from giap.mapping.outputs import MAPPING_FDR_COLUMNS, PEPTIDE_HITS_COLUMNS, UNMAPPED_COLUMNS
    from giap.mapping.protein_assign import PROTEIN_HITS_COLUMNS

    assert _nftest_header("PEPTIDE_HITS") == list(PEPTIDE_HITS_COLUMNS)
    assert _nftest_header("PROTEIN_HITS") == list(PROTEIN_HITS_COLUMNS)
    assert _nftest_header("UNMAPPED") == list(UNMAPPED_COLUMNS)
    assert _nftest_header("MAPPING_FDR") == list(MAPPING_FDR_COLUMNS)
