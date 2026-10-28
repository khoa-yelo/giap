"""M11-18: docs/user-guide/outputs.md and the workflow publish config agree on annotation/tools/.

docs/specs/M11-18-publish-bakta-outputs.md acceptance 4 and 5. The tools/ line of the output tree
names native-output directories (``<tool>/``); each named tool must publish to annotation/tools,
either through a ``publishDir`` on its module or through a ``withName`` block in a config file.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "workflow"
DOCS = ROOT / "docs"


def _tools_line() -> str:
    text = (DOCS / "user-guide" / "outputs.md").read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if re.search(r"[├└]── tools/", ln)]
    assert len(lines) == 1, f"expected one tools/ line in the outputs.md tree, got {lines}"
    return lines[0]


def _listed_tools() -> list[str]:
    line = _tools_line().split("tools/", 1)[1]
    return sorted(set(re.findall(r"\b([a-z][a-z0-9_]*)/", line)))


def _publishes_to_tools(tool: str) -> bool:
    module = WF / "modules" / "local" / tool / "main.nf"
    if module.is_file():
        src = module.read_text(encoding="utf-8")
        if re.search(r"^\s*publishDir\b[^\n]*annotation/tools", src, re.M):
            return True
    proc = tool.upper()
    configs = [WF / "nextflow.config", *sorted((WF / "conf").glob("*.config"))]
    for cfg in configs:
        if not cfg.is_file():
            continue
        text = cfg.read_text(encoding="utf-8")
        block = rf"withName\s*:\s*['\"]?[^'\"\s{{]*\b{proc}\b[^{{]*\{{([^}}]*)\}}"
        for m in re.finditer(block, text):
            if "publishDir" in m.group(1) and "annotation/tools" in m.group(1):
                return True
    return False


def test_outputs_md_tools_line_names_bakta() -> None:
    assert "bakta" in _listed_tools(), _tools_line()


def test_every_listed_tool_publishes_to_annotation_tools() -> None:
    """Acceptance 4: outputs.md lists only tools whose native outputs are published."""
    tools = _listed_tools()
    assert tools, _tools_line()
    missing = [t for t in tools if not _publishes_to_tools(t)]
    assert not missing, (
        f"outputs.md lists annotation/tools/{{{','.join(missing)}}}/ but no publishDir to "
        f"annotation/tools exists for them: {_tools_line()}"
    )


    """Spec: mode copy (the GIPPA_STEP convention); a link into work/ dies with work/."""
    """Spec: mode copy (the GIPPA_STEP convention); a link into work/ breaks with work/."""
    src = (WF / "modules" / "local" / "bakta" / "main.nf").read_text(encoding="utf-8")
    lines = [ln for ln in src.splitlines() if re.match(r"\s*publishDir\b", ln)]
    if not lines:  # a withName config block is also allowed
        assert _publishes_to_tools("bakta"), "BAKTA has no publishDir to annotation/tools"
        return
    assert all("mode: 'copy'" in ln or 'mode: "copy"' in ln for ln in lines), lines


def test_makedb_docs_point_at_published_bakta_json() -> None:
    """Acceptance 5: after a run, makedb --annotation takes annotation/tools/bakta/<sample>.json."""
    for rel in ("user-guide/makedb.md", "get-started/from-spectra.md"):
        text = (DOCS / rel).read_text(encoding="utf-8")
        assert re.search(r"--annotation\s+\S*annotation/tools/bakta/\S+\.json", text), rel
