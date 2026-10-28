"""miniprot builders vs. ``workflow/modules/local/miniprot/main.nf`` (issue M3-12; M2-07 sync rule).

``gippa.mapping.miniprot`` is not under ``gippa.annotate``, so M2-07's ``test_nf_cmd_sync`` does not
discover it; this file applies the same rule: every flag the builders emit (index and align, minus
the configurable ``mapping.proteins.miniprot_args``, which reach the process via ``task.ext.args``)
appears in a ``script:`` block, the script runs ``miniprot`` and the ``stub:`` writes a GFF.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MAIN_NF = REPO / "workflow" / "modules" / "local" / "miniprot" / "main.nf"


def _text() -> str:
    assert MAIN_NF.is_file(), f"{MAIN_NF.relative_to(REPO)} missing"
    text = MAIN_NF.read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _scripts() -> str:
    blocks = re.findall(r"\bscript\s*:(.*?)\bstub\s*:", _text(), flags=re.S)
    assert blocks, "no script: block followed by stub:"
    return "\n".join(blocks)


def _has_flag(text: str, flag: str) -> bool:
    return re.search(rf"(?<![\w-]){re.escape(flag)}(?![\w-])", text) is not None


def test_builder_flags_appear_in_script() -> None:
    from gippa.mapping.miniprot import miniprot_align_cmd, miniprot_index_cmd

    argv = [
        a
        for table in (4, 11, 25)  # M11-16d: -T <annotation.translation_table>
        for a in miniprot_index_cmd(Path("g.fa"), Path("g.mpi"), 2, table=table)
    ] + miniprot_align_cmd(Path("g.mpi"), Path("p.faa"), 2, "")
    flags = {a for a in argv if re.fullmatch(r"--?[A-Za-z][\w-]*", a)}
    assert flags >= {"-t", "-T", "-d", "--gff"}
    script = _scripts()
    missing = sorted(f for f in flags if not _has_flag(script, f))
    assert not missing, f"builder flags missing from main.nf script: {missing}"
    assert re.search(r"(?m)^\s*miniprot\b", script), "script does not run `miniprot`"
    assert re.search(r"-T\s+\$\{[^}]*translation_table[^}]*\}", script), (
        "§7.1.1 / M11-16d: index with -T <annotation.translation_table>"
    )


def test_stub_writes_a_gff() -> None:
    stubs = re.findall(r"\bstub\s*:(.*?)(?:\n\s*\}\s*$|\Z)", _text(), flags=re.S | re.M)
    assert any(".gff" in s for s in stubs), "stub: must write a minimal GFF"
