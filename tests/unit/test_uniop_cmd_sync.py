"""``uniop_cmd`` builder vs. the UNIOP Nextflow process (M2-12; M2-07 nf-cmd sync rule).

M2-index decision 2: every ``giap.annotate.<tool>_cmd`` builder is kept in sync with
``workflow/modules/local/<tool>/main.nf``. M2-07's ``tests/unit/test_nf_cmd_sync.py`` must gain a
``"uniop"`` entry in ``BUILDERS`` once both land; until the UNIOP process exists this check skips.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from giap.annotate.uniop_cmd import uniop_command

REPO = Path(__file__).resolve().parents[2]
MAIN_NF = REPO / "workflow" / "modules" / "local" / "uniop" / "main.nf"


def _script() -> str:
    if not MAIN_NF.is_file():
        pytest.skip("workflow/modules/local/uniop/main.nf not merged yet (M2-07)")
    text = re.sub(r"/\*.*?\*/", "", MAIN_NF.read_text(), flags=re.S)
    text = re.sub(r"(?m)^\s*//.*$", "", text)
    m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", text, flags=re.S)
    assert m, "uniop: no script: block followed by stub:"
    return m.group(1)


def test_builder_flags_appear_in_the_uniop_process() -> None:
    script = _script()
    argv = uniop_command(Path("uniop_cds.gff"), Path("out"), Path("/opt/uniop"))
    flags = {a for a in argv if re.fullmatch(r"--?[A-Za-z][\w-]*", a)}
    assert flags == {"-f", "--gff_source", "-t"}
    missing = [f for f in flags if not re.search(rf"(?<![\w-]){re.escape(f)}(?![\w-])", script)]
    assert not missing, f"builder flags missing from main.nf: {missing}"
    assert "prodigal_gff" in script
    assert re.search(r"src/UniOP", script)
