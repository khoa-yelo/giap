"""NOVEL subworkflow and main.nf wiring against the step registry (issue M7-08 rules 7-8;
GIAP_SPEC.md §6.3, §9.1). Static checks, no Nextflow needed, as
tests/unit/test_mapping_nf_registry.py does for MAPPING:

* ``workflow/subworkflows/local/novel.nf`` takes ``genome, gene_model, contract, mapping,
  peptides, contaminants, gene_table, layers`` (in that order);
* every py step of ``giap.steps.novel`` runs exactly once through a ``GIAP_STEP`` alias named
  after the step, its out-spec keys equal the registry's outputs, and every static ``stepIn`` key
  is a declared input;
* the tools are ``PROMOTECH as PROMOTECH_NOVEL`` and ``OSTIR as OSTIR_NOVEL`` from the existing
  modules (rule 8: no new tool, env or image);
* ``novel_orfs`` publishes the four §9.1 files under ``tables/`` and the local executor's
  ``PUBLISHED_FILES`` gives the same names (executor parity);
* ``main.nf`` includes ``NOVEL`` and calls it inside the peptides block, after ``GENETABLE``,
  with ``GENETABLE.out.gene_table``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
NOVEL_NF = WF / "subworkflows" / "local" / "novel.nf"
MAIN_NF = WF / "main.nf"
TAKE = ["genome", "gene_model", "contract", "mapping", "peptides", "contaminants", "gene_table",
        "layers", "hits"]  # M7-11b: hits  # fmt: skip
PY_STEPS = {"novel_candidates", "novel_promoter_layer", "novel_rbs_layer", "novel_orfs"}
PUBLISHED = {
    "novel_orfs": "tables/novel_orfs.tsv",
    "novel_orfs_gff3": "tables/novel_orfs.gff3",
    "novel_orfs_faa": "tables/novel_orfs.faa",
    "annotation_patch": "tables/annotation_patch.gff3",
    "novel_class_fdr": "novel/novel_class_fdr.parquet",  # M7-11b side table
}


def _registry() -> dict[str, dict[str, list[str]]]:
    from giap.steps import load_step_modules
    from giap.steps.registry import all_steps

    load_step_modules()
    return {
        s.name: {"inputs": sorted(s.inputs), "optional_inputs": sorted(s.optional_inputs),
                 "outputs": sorted(s.outputs)}
        for s in all_steps()
        if s.kind == "py" and getattr(s.func, "__module__", "") == "giap.steps.novel"
    }  # fmt: skip


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _text() -> str:
    assert NOVEL_NF.is_file(), f"{NOVEL_NF.relative_to(REPO)} is missing"
    return _strip(NOVEL_NF.read_text(encoding="utf-8"))


def _call_body(text: str, start: int) -> str:
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _step_vars(text: str) -> dict[str, str]:
    return dict(re.findall(r"def\s+(\w+)\s*=\s*stepIn\(\s*'(\w+)'", text))


def _step_calls(text: str) -> dict[str, tuple[str, list[str], str]]:
    """step -> (alias, out keys, out-spec text) of ``ALIAS(v.files, v.keys, [...])``."""
    var_step = _step_vars(text)
    out: dict[str, tuple[str, list[str], str]] = {}
    for m in re.finditer(r"\b([A-Z_]+)\(\s*(\w+)\.files\s*,\s*\2\.keys\s*,", text):
        alias, var = m.group(1), m.group(2)
        assert var in var_step, f"{alias}: {var} is not a stepIn(...) result"
        step = var_step[var]
        assert step not in out, f"step {step} called twice"
        spec = _call_body(text, m.end())
        keys = re.findall(r"'(\w+)=", spec) + re.findall(r"scatterSpec\(\s*'(\w+)'", spec)
        out[step] = (alias, keys, spec)
    return out


def _includes(text: str) -> dict[str, set[str]]:
    """module path -> 'NAME as ALIAS' / 'NAME' entries included from it."""
    out: dict[str, set[str]] = {}
    for m in re.finditer(r"include\s*\{(.*?)\}\s*from\s*'([^']+)'", text, flags=re.S):
        items = {re.sub(r"\s+", " ", x.strip()) for x in m.group(1).split(";") if x.strip()}
        out.setdefault(m.group(2).rstrip("/"), set()).update(items)
    return out


def _step_in_keys(text: str, step: str) -> set[str]:
    m = re.search(rf"stepIn\(\s*'{step}'\s*,", text)
    assert m, f"no stepIn('{step}', ...) in novel.nf"
    body = _call_body(text, m.end())
    keys = set(re.findall(r"kv\(\s*'(\w+)'", body))
    var_keys = dict(re.findall(r"def\s+(\w+)\s*=\s*kv\(\s*'(\w+)'", text))
    for tok in re.findall(r"\b(\w+)\b", body):
        if tok in var_keys:
            keys.add(var_keys[tok])
    return keys


# ---------------------------------------------------------------- registry and steps


def test_registry_has_the_novel_py_steps() -> None:
    assert set(_registry()) == PY_STEPS


def test_take_block_in_order() -> None:
    m = re.search(r"workflow\s+NOVEL\s*\{\s*take:(.*?)main:", _text(), flags=re.S)
    assert m, "novel.nf defines no `workflow NOVEL { take: ... main: ... }`"
    names = re.findall(r"(?m)^\s*(\w+)\s*$", m.group(1))
    assert names == TAKE


def test_every_novel_step_runs_once_through_its_giap_step_alias() -> None:
    text = _text()
    calls = _step_calls(text)
    assert set(calls) == PY_STEPS, sorted(calls)
    gs = next((v for k, v in _includes(text).items() if k.endswith("giap_step/main")), set())
    for step, (alias, _, _) in calls.items():
        assert alias == step.upper(), f"alias {alias} runs {step}"
        assert f"GIAP_STEP as {alias}" in gs, f"{alias} is not `GIAP_STEP as {alias}`"
    assert set(_step_vars(text).values()) == PY_STEPS


@pytest.mark.parametrize("step", sorted(PY_STEPS))
def test_out_keys_equal_the_registry_outputs(step: str) -> None:
    calls = _step_calls(_text())
    assert step in calls, f"novel.nf does not call {step}"
    assert sorted(set(calls[step][1])) == _registry()[step]["outputs"]


@pytest.mark.parametrize("step", sorted(PY_STEPS))
def test_static_stepin_keys_are_declared_inputs(step: str) -> None:
    reg = _registry()[step]
    extra = _step_in_keys(_text(), step) - set(reg["inputs"]) - set(reg["optional_inputs"])
    assert not extra, f"stepIn('{step}') passes undeclared key(s) {sorted(extra)}"


def test_window_chunks_are_scattered_like_m2_13() -> None:
    """Rule 7: one PROMOTECH_NOVEL task per chunk (``scatterSpec`` on the window directory)."""
    spec = _step_calls(_text())["novel_candidates"][2]
    assert re.search(r"scatterSpec\(\s*'novel_windows'\s*,\s*'novel_windows/windows_\*\.fa'",
                     spec), spec  # fmt: skip


def test_tools_reuse_the_existing_modules() -> None:
    inc = _includes(_text())
    pt = next((v for k, v in inc.items() if k.endswith("modules/local/promotech/main")), set())
    os_ = next((v for k, v in inc.items() if k.endswith("modules/local/ostir/main")), set())
    assert "PROMOTECH as PROMOTECH_NOVEL" in pt, inc
    assert "OSTIR as OSTIR_NOVEL" in os_, inc
    mods = {p.name for p in (WF / "modules" / "local").iterdir()}
    assert not {m for m in mods if "novel" in m}, mods
    envs = {p.name for p in (WF / "envs").iterdir()}
    assert not {e for e in envs if "novel" in e}, envs


def test_steps_never_import_the_gpl_tools() -> None:
    src = (REPO / "src" / "giap" / "steps" / "novel.py").read_text(encoding="utf-8")
    assert not re.search(r"^\s*(import|from)\s+(ostir|promotech)\b", src, flags=re.M)


# ---------------------------------------------------------------- publishing (executor parity)


def test_novel_orfs_publishes_the_four_tables() -> None:
    spec = _step_calls(_text())["novel_orfs"][2]
    got = dict(re.findall(r"'([a-z_0-9]+)=([a-z_]+/[^' ]+)", spec))
    assert got == PUBLISHED


def test_local_executor_publishes_the_same_names() -> None:
    from giap.runner.local import PUBLISHED_FILES

    assert {k: PUBLISHED_FILES.get(k) for k in PUBLISHED} == PUBLISHED


def test_local_publish_copies_the_four_tables(tmp_path: Path) -> None:
    from giap.runner.local import publish

    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    for key in PUBLISHED:
        (work / key).write_text(key + "\n")
    publish(work, out)
    for key, rel in PUBLISHED.items():
        assert (out / rel).read_text() == key + "\n", rel


# ---------------------------------------------------------------- main.nf


def test_main_includes_and_calls_novel_after_genetable() -> None:
    text = _strip(MAIN_NF.read_text(encoding="utf-8"))
    inc = _includes(text)
    assert "NOVEL" in inc.get("./subworkflows/local/novel", set()), inc
    g, n = text.find("GENETABLE("), text.find("NOVEL(")
    assert g >= 0 and n > g, "main.nf must call NOVEL(...) after GENETABLE(...)"
    body = _call_body(text, n + len("NOVEL("))
    assert "GENETABLE.out.gene_table" in body, body
    pep = text.find("if (inputs.peptides)")
    assert 0 <= pep < n, "NOVEL runs only inside the peptides block"
