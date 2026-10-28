"""M11-18: docs/user-guide/outputs.md and the workflow publish config agree on annotation/tools/.

docs/specs/M11-18-publish-bakta-outputs.md acceptance 4 and 5. The tools/ line of the output tree
names native-output directories (``<tool>/``); each named tool must publish to annotation/tools:
through a ``publishDir`` on its module, a ``withName`` block in a config file, or (Bakta, lead
decision on #468) a ``tools/<tool>`` directory output that annotate.nf passes to CONTEXT_LAYERS
as the M2-20 ``tool_dirs`` input (``copy_native`` writes ``annotation/tools/<tool>/``). BAKTA
itself has no publishDir: CONTEXT_LAYERS publishes the whole ``annotation`` directory later and
Nextflow's default overwrite would delete annotation/tools/bakta first.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "workflow"
DOCS = ROOT / "docs"
ANNOTATE = WF / "subworkflows" / "local" / "annotate.nf"


def _tools_line() -> str:
    text = (DOCS / "user-guide" / "outputs.md").read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if re.search(r"[├└]── tools/", ln)]
    assert len(lines) == 1, f"expected one tools/ line in the outputs.md tree, got {lines}"
    return lines[0]


def _listed_tools() -> list[str]:
    line = _tools_line().split("tools/", 1)[1]
    return sorted(set(re.findall(r"\b([a-z][a-z0-9_]*)/", line)))


def _module(tool: str) -> str:
    module = WF / "modules" / "local" / tool / "main.nf"
    return module.read_text(encoding="utf-8") if module.is_file() else ""


def _process_block(proc: str) -> str:
    """The ``withName`` config blocks naming ``proc`` (every config file), concatenated."""
    out = []
    for cfg in [WF / "nextflow.config", *sorted((WF / "conf").glob("*.config"))]:
        if not cfg.is_file():
            continue
        text = cfg.read_text(encoding="utf-8")
        block = rf"withName\s*:\s*['\"]?[^'\"\s{{]*\b{proc}\b[^{{]*\{{([^}}]*)\}}"
        out += [m.group(1) for m in re.finditer(block, text)]
    return "\n".join(out)


def _tool_dir_emit(tool: str) -> str | None:
    """The emit name of the module's ``path 'tools/<tool>'`` output (M2-20 tool_dirs), or None."""
    m = re.search(
        rf"^\s*path\s*\(?\s*['\"]tools/{tool}/?['\"]\s*\)?\s*,\s*emit\s*:\s*(\w+)",
        _module(tool),
        re.M,
    )
    return m.group(1) if m else None


def _wired_as_tool_dirs(tool: str) -> bool:
    """annotate.nf passes ``<TOOL>.out.<emit>`` (the tools/<tool> dir) to a ``tool_dirs`` key."""
    emit = _tool_dir_emit(tool)
    if emit is None:
        return False
    src = ANNOTATE.read_text(encoding="utf-8")
    ref = rf"\b{tool.upper()}\.out\.{emit}\b"
    for m in re.finditer(r"kv\(\s*['\"]tool_dirs['\"]\s*,\s*([^\n]*)", src):
        arg = m.group(1)
        if re.search(ref, arg):
            return True
        var = re.match(r"\s*(\w+)", arg)
        if var and re.search(rf"\b{var.group(1)}\s*=[^\n]*{ref}", src):
            return True
    return False


def _publishes_to_tools(tool: str) -> bool:
    if re.search(r"^\s*publishDir\b[^\n]*annotation/tools", _module(tool), re.M):
        return True
    block = _process_block(tool.upper())
    if "publishDir" in block and "annotation/tools" in block:
        return True
    return _wired_as_tool_dirs(tool)


def test_outputs_md_tools_line_names_bakta() -> None:
    assert "bakta" in _listed_tools(), _tools_line()


def test_every_listed_tool_publishes_to_annotation_tools() -> None:
    """Acceptance 4: outputs.md lists only tools whose native outputs are published."""
    tools = _listed_tools()
    assert tools, _tools_line()
    missing = [t for t in tools if not _publishes_to_tools(t)]
    assert not missing, (
        f"outputs.md lists annotation/tools/{{{','.join(missing)}}}/ but neither a publishDir to "
        f"annotation/tools nor a tools/<tool> dir wired into CONTEXT_LAYERS tool_dirs exists for "
        f"them: {_tools_line()}"
    )


def test_bakta_has_no_publishdir() -> None:
    """Lead decision (#468): a BAKTA publishDir to annotation/tools is deleted again when
    CONTEXT_LAYERS publishes the whole annotation/ directory (Nextflow overwrite), so BAKTA
    publishes nothing itself, neither on the process nor through a withName block."""
    lines = [ln for ln in _module("bakta").splitlines() if re.match(r"\s*publishDir\b", ln)]
    assert not lines, f"BAKTA must not have a publishDir (tool_dirs instead): {lines}"
    assert "publishDir" not in _process_block("BAKTA"), "withName BAKTA sets a publishDir"


def test_bakta_emits_tools_bakta_with_only_json_and_gbff() -> None:
    """BAKTA declares a ``tools/bakta`` directory output holding hard links or copies of its
    declared ``bakta/<id>.json`` and ``.gbff`` only (never the other Bakta files), in the
    script and in the stub."""
    src = _module("bakta")
    assert _tool_dir_emit("bakta"), "BAKTA has no `path 'tools/bakta', emit: ...` output"
    script = src.split("script:", 1)[1].split("stub:", 1)[0]
    stub = src.split("stub:", 1)[1]
    for name, body in (("script", script), ("stub", stub)):
        assert "tools/bakta" in body, f"BAKTA {name} does not create tools/bakta"
        put = [ln for ln in body.splitlines() if "tools/bakta" in ln]
        for ext in ("json", "gbff"):
            assert any(re.search(rf"\b{ext}\b", ln) for ln in put), (
                f"BAKTA {name}: no link/copy of bakta/${{meta_id}}.{ext} into tools/bakta: {put}"
            )
        assert not [ln for ln in put if re.search(r"\bln\s+(-\w*s|--symbolic)", ln)], (
            f"BAKTA {name}: tools/bakta files must be hard links or copies, not symlinks: {put}"
        )
        assert not [
            ln for ln in put if re.search(r"bakta/\*|cp\s+-\w*[rR]|\.(gff3|faa|tsv)\b", ln)
        ], f"BAKTA {name}: tools/bakta gets more than the json and gbff: {put}"


def test_annotate_wires_bakta_tool_dir_into_context_layers() -> None:
    """annotate.nf passes BAKTA's tools/bakta dir to CONTEXT_LAYERS as ``tool_dirs`` (the
    context_layers stepIn), only inside the branch where BAKTA runs."""
    assert _wired_as_tool_dirs("bakta"), "annotate.nf: no kv('tool_dirs', BAKTA.out.<tools dir>)"
    src = ANNOTATE.read_text(encoding="utf-8")
    cl = re.search(r"stepIn\(\s*['\"]context_layers['\"][^\n]*", src)
    assert cl, "annotate.nf: no stepIn('context_layers', ...)"
    emit = _tool_dir_emit("bakta")
    first = src.index(f"BAKTA.out.{emit}")
    bakta_call = src.index("BAKTA(")
    else_at = src.rfind("} else {", 0, bakta_call)
    assert else_at != -1 and else_at < first, "BAKTA.out tools dir used outside the BAKTA branch"
    assert first < cl.start(), "the tools dir is referenced only after the context_layers stepIn"


def test_makedb_docs_point_at_published_bakta_json() -> None:
    """Acceptance 5: after a run, makedb --annotation takes annotation/tools/bakta/<sample>.json."""
    for rel in ("user-guide/makedb.md", "get-started/from-spectra.md"):
        text = (DOCS / rel).read_text(encoding="utf-8")
        assert re.search(r"--annotation\s+\S*annotation/tools/bakta/\S+\.json", text), rel
