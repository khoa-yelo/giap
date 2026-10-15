"""Stub outputs carry the real tool headers (M2-07 lead decisions 5, 6, 7 after the PR #65 review).

The committed real headers live in ``tests/data/recorded/headers/<tool>.<file>.txt`` (see
PROVENANCE.md there). Here: every committed header is the one the module writes (its text appears
in ``workflow/modules/local/<tool>/main.nf``, in the stub block for stubs), the committed copies
match the recordings under ``$GIAP_ROOT/recordings/<tool>/sim100k`` when those exist, and the
Promoter Calculator script writes the 27-column header when the tool wrote no CSV. The exact stub
output files are checked dynamically in ``workflow/tests/annotate.nf.test`` and
``workflow/tests/annotate_tools.nf.test``. No Nextflow needed.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
MODULES = REPO / "workflow" / "modules" / "local"
HEADERS = REPO / "tests" / "data" / "recorded" / "headers"

# committed header file -> (tool, regex for the output path in the stub, recording file glob)
CASES = {
    "genomad.virus_summary.tsv.txt": (
        "genomad",
        r"genome_virus_summary\.tsv",
        "*_virus_summary.tsv",
    ),
    "genomad.plasmid_summary.tsv.txt": (
        "genomad",
        r"genome_plasmid_summary\.tsv",
        "*_plasmid_summary.tsv",
    ),
    "genomad.virus_genes.tsv.txt": ("genomad", r"genome_virus_genes\.tsv", "*_virus_genes.tsv"),
    "ostir.ostir.csv.txt": ("ostir", r"ostir\.csv", "ostir.csv"),
    "promoter_calculator.csv.txt": ("promoter_calculator", r"\$\{meta_id\}\.csv", "*.csv"),
    "promotech.sequences_predictions.csv.txt": (
        "promotech",
        r"sequences_predictions\.csv",
        "sequences_predictions.csv",
    ),
    "uniop.uniop.pred.txt": ("uniop", r"uniop\.pred", "*.pred"),
}


def _header(name: str) -> list[str]:
    path = HEADERS / name
    assert path.is_file(), f"committed real header missing: {path.relative_to(REPO)}"
    return path.read_text(encoding="utf-8").splitlines()


def _main(tool: str) -> str:
    return (MODULES / tool / "main.nf").read_text(encoding="utf-8")


def _stub(tool: str) -> str:
    m = re.search(r"\bstub\s*:(.*)", _main(tool), flags=re.S)
    assert m, f"{tool}: no stub: block"
    return m.group(1)


def _unescape(literal: str, printf: bool) -> str:
    s = literal.replace("\\\\", "\\")  # Groovy GString: "\\t" reaches bash as "\t"
    if printf:
        s = s.replace("\\t", "\t").replace("\\n", "\n")
    return s


def _stub_writes(tool: str, target: str) -> list[str]:
    """Lines the stub's printf/echo writes to the file matching ``target`` (one statement)."""
    for line in _stub(tool).splitlines():
        if not re.search(r">\s*\S*" + target + r"\s*$", line.strip()):
            continue
        m = re.search(r"\b(printf|echo)\s+(['\"])(.*?)\2", line)
        assert m, f"{tool} stub: cannot read the literal written to {target}: {line.strip()}"
        text = _unescape(m.group(3), m.group(1) == "printf")
        return text.rstrip("\n").split("\n")
    raise AssertionError(f"{tool} stub: no printf/echo writes {target}")


@pytest.mark.parametrize("name", sorted(CASES))
def test_stub_writes_the_committed_real_header(name: str) -> None:
    tool, target, _ = CASES[name]
    want = _header(name)
    got = _stub_writes(tool, target)
    assert got == want, f"{tool} stub header differs from the real one ({name}):\n{got}\n{want}"


def test_genomad_virus_genes_header_has_20_columns() -> None:
    assert len(_header("genomad.virus_genes.tsv.txt")[0].split("\t")) == 20
    assert len(_stub_writes("genomad", r"genome_virus_genes\.tsv")[0].split("\t")) == 20


def test_ostir_stub_has_rbs_distance_in_the_real_order() -> None:
    cols = _stub_writes("ostir", r"ostir\.csv")[0].split(",")
    assert cols[4] == "RBS_distance_bp", cols


def test_promotech_stub_writes_sequences_predictions_tab_separated() -> None:
    stub = _stub("promotech")
    assert "sequences_predictions.csv" in stub
    assert "sequence_predictions.csv" not in stub, "real file name is sequences_predictions.csv"
    assert _stub_writes("promotech", r"sequences_predictions\.csv") == ["CHROM\tSEQ\tPRED"]


def test_promoter_calculator_script_writes_the_real_header_when_no_csv() -> None:
    """Decision 6: no promoter -> the tool writes no CSV -> the process writes the header."""
    want = _header("promoter_calculator.csv.txt")[0]
    assert len(want.split(",")) == 27
    m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", _main("promoter_calculator"), flags=re.S)
    assert m, "promoter_calculator: no script: block"
    assert want in m.group(1), "script block does not write the real header for an empty result"


def _recording(tool: str, glob: str) -> Path | None:
    root = os.environ.get("GIAP_ROOT")
    if not root:
        return None
    base = Path(root) / "recordings" / tool / "sim100k"
    hits = sorted(base.rglob(glob)) if base.is_dir() else []
    return hits[0] if hits else None


@pytest.mark.parametrize("name", sorted(CASES))
def test_committed_header_equals_the_recording(name: str) -> None:
    tool, _, glob = CASES[name]
    rec = _recording(tool, glob)
    if rec is None:
        pytest.skip(f"no recording $GIAP_ROOT/recordings/{tool}/sim100k/**/{glob}")
    want = _header(name)
    got = [ln.rstrip("\r") for ln in rec.read_text(encoding="utf-8").splitlines()[: len(want)]]
    assert got == want, f"{rec}: header differs from {name}"
