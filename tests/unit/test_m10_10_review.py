"""#364 review items S1 and S2 (M10-10a): the C7 production path and the ``.ssl`` integrity check.

Spec: docs/specs/M10-10-a3-run.md, step 5 (the three-way compare, rule C7) and the
``cascadia_convert.py`` ``summary`` section. S1: the 4-argument ``pilot_compare.sbatch`` passes
``--pair-cmp $(dirname OUT_A1)/cmp_a1_a2`` to ``run_evidence.py`` when that directory exists, and
the evidence records ``no_pair_cmp`` when it does not (N = 3 only). S2: ``summary`` fails when an
``.ssl`` file's sha256 differs from its ``ssl.sha256`` line. Fixtures: M10-04's compare-job
sandbox (``helpers.pilot_runs_m10_04``) and ``helpers.a3_run_m10_10``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_run_m10_10 import (
    ACQUISITION,
    RUNS,
    TWELVE,
    load_cascadia_convert,
    m10_04,
    sha256,
    write_a3_arm,
    write_ssl,
)

# ------------------------------------------------------------------ S1: C7 in pilot_compare.sbatch


@pytest.fixture
def job(tmp_path: Path):  # type: ignore[no-untyped-def]
    kit = m10_04()
    box = kit.compare_job_sandbox(tmp_path)
    pass1, pass2 = box.runs / "pass1", box.runs / "pass2"
    a1, a2 = kit.two_arms(pass1)
    a3 = write_a3_arm(pass2 / "out_a3", reference=kit.SHA_PEPTIDES)
    return box, a1, a2, a3, pass1, pass2


def _env(box) -> dict[str, str]:  # type: ignore[no-untyped-def]
    return {"GIAP_PYTHON": str(box.recording_python())}


def _evidence_calls(box, name: str) -> list[list[str]]:  # type: ignore[no-untyped-def]
    """The recorded ``run_evidence.py compare`` calls whose ``--out`` ends in ``name``."""
    calls = []
    for c in box.python_calls():
        if not any(a.endswith("run_evidence.py") for a in c) or "--out" not in c:
            continue
        if c[c.index("--out") + 1].endswith(name):
            calls.append(c)
    return calls


def test_four_arguments_pass_pass1_cmp_a1_a2_as_pair_cmp(job) -> None:  # type: ignore[no-untyped-def]
    box, a1, a2, a3, pass1, pass2 = job
    pair = pass1 / "cmp_a1_a2"
    res = box.run(str(a1), str(a2), str(pair), env=_env(box))  # pass 1: the real cmp_a1_a2
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    assert pair.is_dir()
    res = box.run(str(a1), str(a2), str(a3), str(pass2 / "cmp_a1_a2_a3"), env=_env(box))
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    calls = _evidence_calls(box, "compare3.json")
    assert len(calls) == 1, box.python_calls()
    c = calls[0]
    assert "--pair-cmp" in c, c
    assert c[c.index("--pair-cmp") + 1] == str(pair)
    doc = json.loads((pass2 / "evidence" / "compare3.json").read_text())
    assert doc["pair_cmp"] == str(pair)
    assert "no_pair_cmp" not in doc["warnings"], doc["warnings"]
    assert doc["failures"] == [], doc["failures"]


def test_four_arguments_without_pass1_cmp_warn_no_pair_cmp(job) -> None:  # type: ignore[no-untyped-def]
    box, a1, a2, a3, pass1, pass2 = job
    assert not (pass1 / "cmp_a1_a2").exists()
    res = box.run(str(a1), str(a2), str(a3), str(pass2 / "cmp_a1_a2_a3"), env=_env(box))
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    calls = _evidence_calls(box, "compare3.json")
    assert len(calls) == 1 and "--pair-cmp" not in calls[0], calls
    doc = json.loads((pass2 / "evidence" / "compare3.json").read_text())
    assert doc["pair_cmp"] is None
    assert doc["warnings"] == ["no_pair_cmp"], doc["warnings"]
    assert doc["failures"] == [], "no_pair_cmp is a warning, never a failure"


def test_three_arguments_never_pass_pair_cmp_nor_warn(job) -> None:  # type: ignore[no-untyped-def]
    """N = 2: C7 does not apply, even when a ``cmp_a1_a2`` sits next to OUT_A1."""
    box, a1, a2, _, pass1, _ = job
    (pass1 / "cmp_a1_a2" / "tables").mkdir(parents=True)
    cmp = pass1 / "cmp_again"
    res = box.run(str(a1), str(a2), str(cmp), env=_env(box))
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    calls = _evidence_calls(box, "compare.json")
    assert len(calls) == 1 and "--pair-cmp" not in calls[0], calls
    doc = json.loads((pass1 / "evidence" / "compare.json").read_text())
    assert doc.get("warnings", []) == [], doc.get("warnings")


# ------------------------------------------------------------------ S2: the .ssl sha256 check


def _summary(tmp: Path, *, tamper: bool) -> tuple[int, dict[str, Any], Path]:
    from click.testing import CliRunner

    from giap.cli import cli

    a3 = tmp / "A3"
    paths = [write_ssl(a3 / r / f"{r}.ssl", r, TWELVE[r]) for r in RUNS]
    sums = a3 / "ssl.sha256"
    sums.write_text("".join(f"{sha256(p)}  {p}\n" for p in paths))
    pep, log = tmp / "inputs" / "a3.tsv.gz", tmp / "inputs" / "a3.convert.log"
    pep.parent.mkdir()
    argv = ["convert", "--from", "cascadia", *map(str, paths), "--out", str(pep),
            "--source", "denovo"]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    assert res.exit_code == 0, res.output
    log.write_text("".join(ln + "\n" for ln in res.output.splitlines()
                           if ln.startswith("giap convert:")))  # fmt: skip
    victim = paths[1]
    if tamper:  # a blank line: same rows (rule 1 stays quiet), different bytes
        with victim.open("a") as fh:
            fh.write("\n")
    out = tmp / "a3.convert.json"
    rc = load_cascadia_convert().main(
        ["summary", "--ssl-sums", str(sums), "--peptides", str(pep), "--log", str(log),
         "--acquisition", str(ACQUISITION), "--out", str(out)]
    )  # fmt: skip
    doc = json.loads(out.read_text()) if out.is_file() else {}
    return rc, doc, victim


def test_summary_untampered_ssl_passes(tmp_path: Path) -> None:
    rc, doc, _ = _summary(tmp_path, tamper=False)
    assert rc == 0 and doc["failures"] == [], doc.get("failures")


def test_summary_fails_when_an_ssl_differs_from_ssl_sha256(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc, doc, victim = _summary(tmp_path, tamper=True)
    err = capsys.readouterr().err
    assert rc == 1, doc.get("failures")
    failures = [str(f) for f in doc["failures"]]
    assert len(failures) == 1, failures
    assert failures[0].startswith(f"ssl: sha256 of {victim} is {sha256(victim)}"), failures
    assert "ssl.sha256 lists" in failures[0]
    assert doc["ssl"][victim.stem]["sha256"] == sha256(victim), "the file's real hash is recorded"
    assert "cascadia_convert summary: ssl: sha256 of" in err
