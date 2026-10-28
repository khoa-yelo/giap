"""tblastn builders vs. ``workflow/modules/local/tblastn/main.nf`` (issue M3-13; M2-07 sync rule).

Same rule as tests/unit/test_miniprot_nf.py: every flag the Python builders emit (minus the
configurable ``mapping.proteins.tblastn_args``, which reach the process via ``task.ext.args``)
appears in the ``script:`` block; the script runs ``makeblastdb`` then ``tblastn`` with the
15-column ``OUTFMT``; the ``stub:`` writes the tabular file; the process uses
``workflow/envs/blast.yml`` and has a pinned BLAST image; ``conf/base.config`` wires
``params.mapping.proteins.tblastn_args`` into ``ext.args``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MAIN_NF = REPO / "workflow" / "modules" / "local" / "tblastn" / "main.nf"
CONTAINERS = REPO / "workflow" / "conf" / "containers.config"
BASE = REPO / "workflow" / "conf" / "base.config"
BLAST_DIGEST = "9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c"


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _text() -> str:
    assert MAIN_NF.is_file(), f"{MAIN_NF.relative_to(REPO)} missing"
    return _strip_comments(MAIN_NF.read_text())


def _scripts() -> str:
    blocks = re.findall(r"\bscript\s*:(.*?)\bstub\s*:", _text(), flags=re.S)
    assert blocks, "no script: block followed by stub:"
    return "\n".join(blocks)


def _has_flag(text: str, flag: str) -> bool:
    return re.search(rf"(?<![\w-]){re.escape(flag)}(?![\w-])", text) is not None


def test_process_name_and_conda_env() -> None:
    text = _text()
    assert re.search(r"\bprocess\s+TBLASTN\s*\{", text), "process TBLASTN not defined"
    assert re.search(r"conda\s+[\"'].*envs/blast\.yml[\"']", text), "conda env blast.yml missing"


def test_builder_flags_appear_in_script() -> None:
    from gippa.mapping.tblastn import OUTFMT, makeblastdb_cmd, tblastn_cmd

    argv = makeblastdb_cmd(Path("g.fa"), Path("db/g")) + tblastn_cmd(
        Path("db/g"), Path("p.faa"), 2, ""
    )
    flags = {a for a in argv if re.fullmatch(r"-[A-Za-z][\w-]*", a)}
    assert flags >= {"-in", "-dbtype", "-out", "-db", "-query", "-num_threads"}
    script = _scripts()
    missing = sorted(f for f in flags if not _has_flag(script, f))
    assert not missing, f"builder flags missing from main.nf script: {missing}"
    mk = re.search(r"(?m)^\s*makeblastdb\b", script)
    tb = re.search(r"(?m)^\s*tblastn\b", script)
    assert mk and tb, "script must run makeblastdb and tblastn"
    assert mk.start() < tb.start(), "makeblastdb must run before tblastn"
    assert OUTFMT in script, "the 15-column OUTFMT must appear verbatim in the script"
    assert "task.ext.args" in script


def test_no_hard_coded_search_defaults() -> None:
    """``-db_gencode 11 -seg no -evalue 1e-5`` come from config via ext.args only."""
    script = _scripts()
    for flag in ("-db_gencode", "-seg", "-evalue"):
        assert not _has_flag(script, flag), f"{flag} hard-coded in main.nf"


def test_stub_writes_tsv() -> None:
    stubs = re.findall(r"\bstub\s*:(.*?)(?:\n\s*\}\s*$|\Z)", _text(), flags=re.S | re.M)
    assert any("tblastn.tsv" in s for s in stubs), "stub: must write tblastn.tsv"


def test_container_pinned_for_tblastn() -> None:
    conf = _strip_comments(CONTAINERS.read_text())
    sel = re.findall(
        r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{\s*container\s*=\s*['\"]([^'\"]+)", conf
    )
    matching = [img for pat, img in sel if re.fullmatch(pat, "TBLASTN")]
    assert matching, "no withName selector in containers.config matches the TBLASTN process"
    assert all(BLAST_DIGEST in img for img in matching), matching


def test_base_config_wires_tblastn_args() -> None:
    conf = _strip_comments(BASE.read_text())
    wired = [
        body
        for pat, body in re.findall(
            r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{(.*?)\}\s*$", conf, re.M
        )
        if re.fullmatch(pat, "TBLASTN")
    ]
    assert wired, "conf/base.config: no withName selector for TBLASTN"
    assert any("tblastn_args" in body and "ext.args" in body for body in wired), wired
