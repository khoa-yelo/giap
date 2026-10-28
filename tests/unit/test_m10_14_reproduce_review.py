"""``reproduce_check.py`` and ``pilot_reproduce_check.sbatch``, review round 1 of M10-14a
(docs/specs/M10-14-pilot-reproducibility.md, "Lead decisions (2026-10-24), review round 1": S7,
S8, N2; examples 7 and 8).

The S-check job runs in a sandbox checkout (``helpers.slurm_sandbox.make_sandbox`` plus
``validation/``) with a stub ``gippa.cli`` on ``$GIPPA_REPO/src``: ``convert`` writes
``table:<report bytes>`` and ``abundance:<report bytes>``, or raises for the stems named in
``$FAKE_CONVERT_CRASH``. No Slurm, no network, no DIA-NN output. No statistical claim, so no
calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from helpers.m10_14_kit import (
    EVIDENCE_S5,
    OTHER_COMMIT,
    REPRO_SBATCH,
    STEP_S0,
    STEP_S7,
    F,
    evidence_s8,
    load_reproduce_check,
    run_cli,
    runbook,
    s5,
    sha256_bytes,
    synthetic_repo,
)
from helpers.slurm_sandbox import REPO, Sandbox, assert_honours_gippa_repo, make_sandbox

STEP_SCHECK = """### S-check reproducibility check (M10-14)

Job: none (session)

```bash
true
```

Reads: -
Writes: -
"""


@pytest.fixture(scope="module")
def rc() -> Any:
    return load_reproduce_check()


def _issues(rc: Any, text: str, repo: Path) -> list[str]:
    return [s for s in rc.static_check(rc.parse(text), repo) if not s.startswith("duplicate_hash:")]


# ---------------------------------------------------------------- S8 (a): Reads ids


def test_s8a_unknown_read_id_fails_static(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    step = STEP_S7.replace("Reads: genome", "Reads: genome, wt9_raw")
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), step]), repo)
    hits = [s for s in issues if s.startswith("unknown_read")]
    assert hits and all("wt9_raw" in s for s in hits), issues
    assert any("S7" in s for s in hits), hits


def test_s8a_outputs_ids_and_dashes_are_known_reads(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    out = f"| a1_table | $GIPPA_ROOT/pilot/inputs/a1.tsv.gz | {'d' * 64} | - | yes |"
    step = STEP_S7.replace("Reads: genome", "Reads: genome, a1_table")
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), step], outputs=[out]), repo)
    assert not [s for s in issues if s.startswith("unknown_read")], issues


# ---------------------------------------------------------------- S8 (b): step order


def test_s8b_numbered_steps_out_of_order_fail_static(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    issues = _issues(rc, runbook(steps=[STEP_S0, STEP_S7, s5()]), repo)
    assert any(s.startswith("step_order") for s in issues), issues


def test_s8b_word_steps_do_not_break_the_order(rc: Any, tmp_path: Path) -> None:
    """Only the numbered steps must ascend; an ``S-<word>`` step may stand between them."""
    repo = synthetic_repo(tmp_path)
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), STEP_SCHECK, STEP_S7]), repo)
    assert not [s for s in issues if s.startswith("step_order")], issues
    assert issues == [], issues


def test_s8b_numeric_not_lexical_order(rc: Any, tmp_path: Path) -> None:
    """S7 before S10 is ascending (10 > 7), although "S10" < "S7" as strings."""
    repo = synthetic_repo(tmp_path)
    s10 = STEP_S0.replace("### S0 checkout (M10-14)", "### S10 later (M10-14)")
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), STEP_S7, s10]), repo)
    assert not [s for s in issues if s.startswith("step_order")], issues


# ---------------------------------------------------------------- S8 (c): final needs a branch


def test_s8c_final_runbook_with_pending_a3_branch_fails(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path, verdict="fail")
    text = runbook(
        status="final",
        freeze=F,
        a3_branch="pending",
        evidence=[EVIDENCE_S5, evidence_s8()],
    )
    issues = _issues(rc, text, repo)
    assert "a3_pending_in_final" in issues, issues
    draft = runbook(a3_branch="pending")
    assert "a3_pending_in_final" not in _issues(rc, draft, repo)


# ---------------------------------------------------------------- N2: the CLI exits 2


@pytest.mark.parametrize("cmd", ["static", "verify", "protected"])
def test_n2_missing_runbook_exits_2(tmp_path: Path, cmd: str) -> None:
    missing = tmp_path / "no_such_reproduce.md"
    extra = ["--out", str(tmp_path / "v.json")] if cmd == "verify" else []
    res = run_cli(cmd, str(missing), *extra, env={"GIPPA_ROOT": str(tmp_path)})
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


# ---------------------------------------------------------------- S7: the S-check job

_STUB_CLI = r'''"""Stub gippa.cli for the S-check sandbox (test_m10_14_reproduce_review.py)."""
import json, os, sys


def cli():
    a = sys.argv[1:]
    with open(os.path.join(os.environ["FAKE_SLURM_DIR"], "convert.jsonl"), "a") as fh:
        fh.write(json.dumps(a) + "\n")
    assert a[0] == "convert" and a[1:3] == ["--from", "diann"], a
    report = a[3]
    out = a[a.index("--out") + 1]
    ab = a[a.index("--abundance-out") + 1]
    stem = os.path.basename(out).split(".")[0]
    if stem in os.environ.get("FAKE_CONVERT_CRASH", "").split(","):
        raise RuntimeError(f"fake converter crash for {stem}")
    data = open(report, "rb").read()
    open(out, "wb").write(b"table:" + data)
    open(ab, "wb").write(b"abundance:" + data)
'''

REPORT_A1 = b"A1 report bytes\n"
REPORT_A2 = b"A2 report bytes\n"


@dataclass
class ReproSandbox:
    sb: Sandbox

    @property
    def root(self) -> Path:
        return self.sb.root

    @property
    def pilot(self) -> Path:
        return self.sb.root / "pilot"

    def commit(self, text: str) -> None:
        (self.sb.repo / "validation" / "pilot" / "reproduce.md").write_text(text)
        for args in (["add", "-A"], ["commit", "-q", "-m", "runbook"]):
            subprocess.run(
                ["git", "-C", str(self.sb.repo), *args],
                check=True,
                capture_output=True,
                env=self.sb.env,
            )

    def run(self, **env: str) -> subprocess.CompletedProcess[str]:
        script = self.sb.repo / "scripts" / "slurm" / "pilot_reproduce_check.sbatch"
        assert script.is_file(), "scripts/slurm/pilot_reproduce_check.sbatch does not exist"
        return self.sb.run(["bash", str(script)], env=env, cwd=self.sb.tmp, timeout=120)

    def verify_json(self) -> dict[str, Any]:
        p = self.root / "pilot" / "reproduce" / f"verify-{self.sb.head[:7]}.json"
        assert p.is_file(), f"no {p}"
        doc = json.loads(p.read_text())
        assert isinstance(doc, dict), doc
        return doc

    def converts(self) -> list[list[str]]:
        return self.sb.jsonl("convert.jsonl")

    def repo_status(self) -> str:
        cmd = ["git", "-C", str(self.sb.repo), "status", "--porcelain"]
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()


def _scheck_runbook(prose: str = "", freeze: str = "pending") -> str:
    sha = sha256_bytes(REPORT_A1)
    row = (
        f"| a1_report | $GIPPA_ROOT/pilot/A1/A1_report.parquet | {sha} | "
        "$GIPPA_ROOT/pilot/A1/outputs.sha256 |"
    )
    return runbook(inputs=[row], steps=[STEP_S0], evidence=[], prose=prose, freeze=freeze)


@pytest.fixture
def rs(tmp_path: Path) -> ReproSandbox:
    assert REPRO_SBATCH.is_file(), "scripts/slurm/pilot_reproduce_check.sbatch does not exist"
    assert_honours_gippa_repo(REPRO_SBATCH.read_text(), "pilot_reproduce_check.sbatch")
    sb = make_sandbox(tmp_path.resolve())
    shutil.copytree(
        REPO / "validation",
        sb.repo / "validation",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    pkg = sb.repo / "src" / "gippa"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "cli.py").write_text(_STUB_CLI)
    sb.env.update({"GIPPA_PYTHON": sys.executable, "PYTHONDONTWRITEBYTECODE": "1"})
    r = ReproSandbox(sb)
    a1 = r.pilot / "A1"
    a1.mkdir(parents=True)
    (a1 / "A1_report.parquet").write_bytes(REPORT_A1)
    (a1 / "outputs.sha256").write_text(f"{sha256_bytes(REPORT_A1)}  A1_report.parquet\n")
    a2 = r.pilot / "A2"
    a2.mkdir(parents=True)
    (a2 / "A2_report.parquet").write_bytes(REPORT_A2)
    inputs = r.pilot / "inputs"
    inputs.mkdir(parents=True)
    (inputs / "a1.tsv.gz").write_bytes(b"table:" + REPORT_A1)
    (inputs / "a1_abundance.tsv.gz").write_bytes(b"abundance:" + REPORT_A1)
    (inputs / "a2.tsv.gz").write_bytes(b"an older table\n")
    (inputs / "a2_abundance.tsv.gz").write_bytes(b"abundance:" + REPORT_A2)
    r.commit(_scheck_runbook())
    return r


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return f"exit {res.returncode}\n--- stdout\n{res.stdout}\n--- stderr\n{res.stderr}"


def test_s7_examples_7_and_8_identical_and_changed(rs: ReproSandbox) -> None:
    """A1 re-converts to the same bytes (identical); A2 to other bytes (changed: a warning)."""
    before = {p: p.read_bytes() for p in (rs.pilot / "inputs").iterdir()}
    res = rs.run()
    assert res.returncode == 0, _out(res)
    doc = rs.verify_json()
    assert doc["reconvert"] == {"a1": "identical", "a2": "changed"}, doc
    assert doc["static_ok"] is True, doc
    assert doc["failures"] == [], doc
    assert doc["inputs"]["a1_report"]["status"] == "ok", doc
    assert doc["freeze_commit"] == "pending", doc
    assert "a2" in res.stderr.lower() or "A2" in res.stderr, "no warning for the changed table"
    assert {p: p.read_bytes() for p in (rs.pilot / "inputs").iterdir()} == before, (
        "the re-conversion wrote over $PILOT/inputs"
    )
    calls = {c[c.index("--out") + 1].rsplit("/", 1)[-1]: c for c in rs.converts()}
    assert calls["a1.tsv.gz"][calls["a1.tsv.gz"].index("--search-db") + 1] == "other"
    assert calls["a2.tsv.gz"][calls["a2.tsv.gz"].index("--search-db") + 1] == "same"
    assert rs.repo_status() == "", f"wrote into $GIPPA_REPO: {rs.repo_status()}"


def test_s7_absent_report(rs: ReproSandbox) -> None:
    (rs.pilot / "A2" / "A2_report.parquet").unlink()
    res = rs.run()
    assert res.returncode == 0, _out(res)
    assert rs.verify_json()["reconvert"] == {"a1": "identical", "a2": "absent"}


def test_s7_convert_crash_is_a_failure_recorded_in_the_json(rs: ReproSandbox) -> None:
    """Spec: a crash is a *failure* (exit 1), not a precondition (2); verify-<F7>.json must not
    be left with ``failures: []`` and no ``reconvert``, which only the clean-up's reconvert
    check would then catch."""
    res = rs.run(FAKE_CONVERT_CRASH="a1")
    assert res.returncode == 1, _out(res)
    doc = rs.verify_json()
    assert any("a1" in str(f).lower() for f in doc.get("failures", [])), doc
    assert doc.get("reconvert", {}).get("a1") not in ("identical", "changed", "absent"), doc


def test_s7_static_failure_exits_1_and_records_static_ok_false(rs: ReproSandbox) -> None:
    rs.commit(_scheck_runbook(prose=f"A stray hash {'ab12' * 16} in the prose."))
    res = rs.run()
    assert res.returncode == 1, _out(res)
    assert rs.verify_json()["static_ok"] is False


def test_s7_runbook_naming_another_commit_exits_2(rs: ReproSandbox) -> None:
    rs.commit(_scheck_runbook(freeze=OTHER_COMMIT))
    res = rs.run()
    assert res.returncode == 2, _out(res)
    assert rs.converts() == []


def test_s7_dirty_checkout_exits_2(rs: ReproSandbox) -> None:
    rb = rs.sb.repo / "validation" / "pilot" / "reproduce.md"
    rb.write_text(rb.read_text() + "\nuncommitted\n")
    res = rs.run()
    assert res.returncode == 2, _out(res)
    assert rs.converts() == []
