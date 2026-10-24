"""The path helpers and the assertion messages of the M10-03 contract test (review round 1).

Spec: docs/specs/M10-03-diann-convert.md, "The contract test" and job step 7; review findings
B1 (the job's report and arm reach the contract test through ``GIAP_M10_03_REPORT`` and
``GIAP_M10_03_ARM``), B2 (a report under a git-ignored, untracked directory of the checkout, the
FarmShare default ``GIAP_ROOT=$GIAP_REPO/local``, is allowed) and N5 (a failing check prints
counts and IDs only, never a peptide or protein string from the report; criterion 10 pastes the
output into the PR).

The contract module is imported by name, so pytest's assertion rewriting applies to it as it does
in the job; ``test_contract_module_is_assertion_rewritten`` checks that, so the N5 tests are not
vacuous. Report data here are invented two-row frames handed to the checks through monkeypatched
``_raw`` / ``_converted``; nothing reads a DIA-NN report.
"""

from __future__ import annotations

import importlib
import subprocess
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from giap.io.issues import Issue

CONTRACT_DIR = Path(__file__).resolve().parents[1] / "contract"
SECRET_PEP = "SECRETPEPTIDEK"
SECRET_PROT = "SECRETPROT1"
SECRET_CONT = "SECRETPROT2"
SECRETS = (SECRET_PEP, SECRET_PROT, SECRET_CONT, "SECRETPEPTIDE")


@pytest.fixture
def contract(monkeypatch: pytest.MonkeyPatch, request: pytest.FixtureRequest) -> ModuleType:
    """The contract module; assertion messages rendered as under ``-vv`` (full diffs), so a
    message that is clean only because pytest truncated it still fails the N5 tests."""
    monkeypatch.setattr(request.config.option, "verbose", 2)
    monkeypatch.syspath_prepend(str(CONTRACT_DIR))
    mod = importlib.import_module("test_m10_03_diann_contract")
    for name in ("_raw", "_converted", "_fasta_ids"):
        getattr(mod, name).cache_clear()
    return mod


def test_contract_module_is_assertion_rewritten(contract: ModuleType) -> None:
    """Guard for the N5 tests: the messages checked below are pytest's rewritten ones."""
    assert "@py_builtins" in vars(contract), "the contract module was imported without rewriting"


# ------------------------------------------------------------------ B1: the job's report and arm


def _git(repo: Path, *args: str) -> None:
    cmd = ["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", "-C", str(repo)]
    subprocess.run([*cmd, *args], check=True, capture_output=True)


def test_full_report_path_uses_the_job_report(
    contract: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep = tmp_path / "elsewhere" / "A1" / "A1_report.parquet"
    rep.parent.mkdir(parents=True)
    rep.write_bytes(b"PAR1")
    (tmp_path / "root" / "pilot" / "A1").mkdir(parents=True)
    (tmp_path / "root" / "pilot" / "A1" / "A1_report.parquet").write_bytes(b"PAR1-other")
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "root"))
    monkeypatch.setenv("GIAP_M10_03_REPORT", str(rep))
    monkeypatch.setenv("GIAP_M10_03_ARM", "A1")
    assert contract.full_report_path("A1") == rep


def test_full_report_path_never_checks_another_arm_against_the_job_report(
    contract: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The A1 job's report is never checked as A2 (with A2's ``--search-db``)."""
    rep = tmp_path / "A1_report.parquet"
    rep.write_bytes(b"PAR1")
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "root"))
    monkeypatch.setenv("GIAP_M10_03_REPORT", str(rep))
    monkeypatch.setenv("GIAP_M10_03_ARM", "A1")
    try:
        got = contract.full_report_path("A2")
    except pytest.skip.Exception:
        return
    assert got != rep


def test_full_report_path_without_job_env_reads_giap_root(
    contract: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Unchanged default: ``$GIAP_ROOT/pilot/<arm>/<arm>_report.parquet``."""
    rep = tmp_path / "root" / "pilot" / "A2" / "A2_report.parquet"
    rep.parent.mkdir(parents=True)
    rep.write_bytes(b"PAR1")
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "root"))
    monkeypatch.delenv("GIAP_M10_03_REPORT", raising=False)
    monkeypatch.delenv("GIAP_M10_03_ARM", raising=False)
    assert contract.full_report_path("A2") == rep


# ------------------------------------------------------------- B2: ignored dirs in the checkout


@pytest.fixture
def repo(tmp_path: Path, contract: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A temporary git checkout with ``/local/`` ignored, standing in for the contract's REPO."""
    r = tmp_path / "checkout"
    r.mkdir()
    _git(r, "init", "-q", "-b", "main")
    (r / ".gitignore").write_text("/local/\n")
    (r / "README.md").write_text("x\n")
    _git(r, "add", "-A")
    _git(r, "commit", "-q", "-m", "init")
    monkeypatch.setattr(contract, "REPO", r)
    return r


def _file(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"PAR1")
    return path


def test_outside_checkout_accepts_an_ignored_untracked_path(
    contract: ModuleType, repo: Path
) -> None:
    """Review B2: ``GIAP_ROOT=$GIAP_REPO/local`` (farmshare_env.sh's default)."""
    p = _file(repo / "local" / "pilot" / "A1" / "A1_report.parquet")
    assert contract._outside_checkout(p) == p


def test_sample_path_under_the_default_giap_root(
    contract: ModuleType, repo: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review B2, end to end: the sample under ``$GIAP_REPO/local`` is found, not refused."""
    p = _file(repo / "local" / "pilot/recorded/diann/2.3.1/A1/report.sample.parquet")
    monkeypatch.setenv("GIAP_ROOT", str(repo / "local"))
    assert contract.sample_path("A1") == p


def test_outside_checkout_accepts_a_path_outside(
    contract: ModuleType, repo: Path, tmp_path: Path
) -> None:
    p = _file(tmp_path / "root" / "pilot" / "A1" / "A1_report.parquet")
    assert contract._outside_checkout(p) == p


def test_outside_checkout_refuses_an_unignored_path(contract: ModuleType, repo: Path) -> None:
    p = _file(repo / "derived" / "A1_report.parquet")
    with pytest.raises(AssertionError):
        contract._outside_checkout(p)


def test_outside_checkout_refuses_a_tracked_file_under_an_ignored_dir(
    contract: ModuleType, repo: Path
) -> None:
    """``git add -f`` under ``local/``: tracked, so it is committed data, refused."""
    p = _file(repo / "local" / "pilot" / "A1" / "A1_report.parquet")
    _git(repo, "add", "-f", str(p))
    _git(repo, "commit", "-q", "-m", "forced")
    with pytest.raises(AssertionError):
        contract._outside_checkout(p)


# ------------------------------------------------------------------ N5: no report strings printed


def _clean(excinfo: pytest.ExceptionInfo[BaseException], *extra: str) -> None:
    text = str(excinfo.value) + "\n" + excinfo.exconly()
    leaked = [s for s in (*SECRETS, *extra) if s in text]
    assert not leaked, f"the assertion message prints report strings {leaked}:\n{text}"


def _raw_frame(**over: Any) -> pd.DataFrame:
    d: dict[str, Any] = {
        "Run": ["WT_1", "WT_2"],
        "Stripped.Sequence": [SECRET_PEP, SECRET_PEP],
        "Modified.Sequence": [SECRET_PEP, SECRET_PEP],
        "Precursor.Charge": [2, 2],
        "Decoy": [0, 0],
        "Global.Q.Value": [0.001, 0.002],
        "Protein.Ids": [SECRET_PROT, SECRET_PROT],
        "Protein.Group": [SECRET_PROT, SECRET_PROT],
        "PG.MaxLFQ": [5.0, 5.0],
    }
    d.update(over)
    return pd.DataFrame(d)


def _pep_frame(**over: Any) -> pd.DataFrame:
    d: dict[str, Any] = {
        "run": ["WT_1", "WT_2"],
        "peptide": [SECRET_PEP, SECRET_PEP],
        "modified_peptide": [SECRET_PEP, SECRET_PEP],
        "charge": pd.array([2, 2], dtype="Int64"),
        "is_decoy": pd.array([False, False], dtype="boolean"),
        "q_value": [0.001, 0.002],
        "protein_ids": [SECRET_PROT, SECRET_PROT],
    }
    d.update(over)
    return pd.DataFrame(d)


def _patch(
    contract: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    raw: pd.DataFrame,
    peptides: pd.DataFrame,
    abundance: pd.DataFrame | None = None,
    issues: list[Issue] | None = None,
) -> Path:
    res = SimpleNamespace(peptides=peptides, abundance=abundance, issues=issues or [])
    monkeypatch.setattr(contract, "_raw", lambda p: raw)
    monkeypatch.setattr(contract, "_converted", lambda p, s: res)
    path = Path("/nonexistent/A1/A1_report.parquet")
    monkeypatch.setattr(contract, "sample_path", lambda arm: path)
    return path


def test_f10_failure_prints_no_peptide(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _patch(contract, monkeypatch, _raw_frame(), _pep_frame(charge=pd.array([2, 3], "Int64")))
    with pytest.raises(AssertionError) as ei:
        contract.check_f10_rows(path, "A1")
    _clean(ei)


def test_f3_failure_prints_no_peptide(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    _patch(contract, monkeypatch, _raw_frame(), _pep_frame(q_value=[0.001, 0.5]))
    with pytest.raises(AssertionError) as ei:
        contract.test_f3_q_value_is_global_q("A1", lambda *a: None)
    _clean(ei)


@pytest.mark.parametrize("bad", ["strip", "shape"])
def test_f5_failure_prints_no_peptide(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    low = SECRET_PEP.lower()
    if bad == "strip":
        pep = _pep_frame(peptide=["SECRETPEPTIDE", SECRET_PEP])
    else:  # strips cleanly but is lower case, failing the shape check
        pep = _pep_frame(peptide=[low, SECRET_PEP], modified_peptide=[low, SECRET_PEP])
    _patch(contract, monkeypatch, _raw_frame(), pep)
    with pytest.raises(AssertionError) as ei:
        contract.test_f5_modifications("A1")
    _clean(ei, low)


@pytest.mark.parametrize("bad", ["value", "groups", "duplicate"])
def test_f9_failure_prints_no_protein(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    ab = {
        "value": pd.DataFrame({"protein_id": [SECRET_PROT], "maxlfq": [6.0]}),
        "groups": pd.DataFrame({"protein_id": [SECRET_CONT], "maxlfq": [5.0]}),
        "duplicate": pd.DataFrame({"protein_id": [SECRET_PROT] * 2, "maxlfq": [5.0, 5.0]}),
    }[bad]
    _patch(contract, monkeypatch, _raw_frame(), _pep_frame(), ab)
    with pytest.raises(AssertionError) as ei:
        contract.test_f9_maxlfq_matches_an_independent_oracle("A1")
    _clean(ei)


def test_f9_per_run_failure_prints_no_protein(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw = _raw_frame(Run=["WT_1", "WT_1"], **{"PG.MaxLFQ": [5.0, 6.0]})
    ab = pd.DataFrame({"protein_id": [SECRET_PROT], "maxlfq": [5.5]})
    _patch(contract, monkeypatch, raw, _pep_frame(), ab)
    with pytest.raises(AssertionError) as ei:
        contract.test_f9_maxlfq_matches_an_independent_oracle("A1")
    _clean(ei)


@pytest.mark.parametrize("bad", ["unknown", "dropped_prefix"])
def test_f7_failure_prints_no_protein(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    known = {"P1"} if bad == "unknown" else {"P1", f"Cont_{SECRET_CONT}"}
    prot = SECRET_PROT if bad == "unknown" else SECRET_CONT
    path = _patch(contract, monkeypatch, _raw_frame(), _pep_frame(protein_ids=[prot, prot]))
    monkeypatch.setattr(contract, "fasta_paths", lambda arm: [Path("/nonexistent.fasta")])
    monkeypatch.setattr(contract, "_fasta_ids", lambda paths: frozenset(known))
    with pytest.raises(AssertionError) as ei:
        contract.check_f7_protein_ids(path, "A1")
    _clean(ei)


def test_f2_failure_prints_no_report_value(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Issue messages may quote cell values; only codes and row counts are printed."""
    issue = Issue("bad_type:PEP", "warn", f"PEP: invalid value {SECRET_PEP!r}", 1)
    path = _patch(contract, monkeypatch, _raw_frame(), _pep_frame(), issues=[issue])
    with pytest.raises(AssertionError) as ei:
        contract.check_f2_types(path, "A1")
    _clean(ei)


def test_f4_failure_prints_no_protein(
    contract: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """100 decoys (the sample's count), one without ``protein_ids``; the others carry protein
    strings, which the failure must not print."""
    n = 100
    raw = pd.DataFrame({"Decoy": [1] * n, "Run": ["WT_1"] * n})
    pep = pd.DataFrame(
        {
            "is_decoy": pd.array([True] * n, dtype="boolean"),
            "protein_ids": [SECRET_PROT] * (n - 1) + [None],
            "peptide": [SECRET_PEP] * n,
        }
    )
    _patch(contract, monkeypatch, raw, pep)
    with pytest.raises(AssertionError) as ei:
        contract.test_f4_decoys("A1", lambda *a: None)
    _clean(ei)
