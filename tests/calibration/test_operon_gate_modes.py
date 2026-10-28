"""Harness modes of the M6-09 gate helpers (``_operon_gate``; PR #310 review B-1).

Under ``--shard`` (``GIPPA_CALIBRATION_REPLICATES``) ``calibrate*`` returns ``shard`` records and in
a ``--merge`` confirm run (``GIPPA_CALIBRATION_CONFIRM``) an unlisted check returns ``skipped``:
neither has a verdict, so the gate helpers must skip (after logging), never fail, and the
non-vacuity requirement must not fire on them. Toy count records, no simulation; every-commit tier.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _criteria as C
import _operon_gate as G
import pytest

SUITE = "missed_operons"


def _records(n_called: int, n_false: int) -> Any:
    def records(i: int) -> dict[str, Any]:
        rec = {"n_called": n_called, "n_false": n_false, "n_planted_called": 0, "n_planted": 0}
        return {f"q{a:g}": dict(rec) for a in G.levels()}

    return records


@pytest.fixture
def clean_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    for var in (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def test_shard_run_skips_and_logs(clean_env: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "0:2")
    recs = G.missed_fdr_gate(SUITE, "operon_gate_modes.shard", _records(0, 0), 4,
                             log_dir=clean_env)  # fmt: skip
    assert {r["status"] for r in recs} == {C.SHARD_STATUS}
    with pytest.raises(pytest.skip.Exception):
        G.assert_fdr_ok(recs, allow_vacuous=False, need_calls_at=0.1)
    logged = (clean_env / f"{SUITE}.jsonl").read_text().splitlines()
    assert len(logged) == len(G.levels())  # the shard records reach --merge


def test_confirm_run_unlisted_check_skips(clean_env: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.CONFIRM_ENV, json.dumps([]))
    recs = G.missed_fdr_gate(SUITE, "operon_gate_modes.confirm", _records(5, 0), 4,
                             log_dir=clean_env)  # fmt: skip
    assert {r["status"] for r in recs} == {C.SKIPPED_STATUS}
    with pytest.raises(pytest.skip.Exception):
        G.assert_fdr_ok(recs, allow_vacuous=True)


def test_unsharded_verdicts(clean_env: Path) -> None:
    """Unsharded: no calls is ``vacuous`` (accepted only when allowed, never at ``need_calls_at``),
    clean calls pass, all-false calls are a confirmed ``fail``."""
    vac = G.missed_fdr_gate(SUITE, "operon_gate_modes.vacuous", _records(0, 0), 4,
                            log_dir=clean_env)  # fmt: skip
    assert {r["status"] for r in vac} == {"vacuous"}
    G.assert_fdr_ok(vac, allow_vacuous=True)
    with pytest.raises(AssertionError):
        G.assert_fdr_ok(vac, allow_vacuous=True, need_calls_at=0.1)
    ok = G.missed_fdr_gate(SUITE, "operon_gate_modes.ok", _records(10, 0), 4, log_dir=clean_env)
    G.assert_fdr_ok(ok, allow_vacuous=False, need_calls_at=0.1)
    bad = G.missed_fdr_gate(SUITE, "operon_gate_modes.bad", _records(10, 10), 4, log_dir=clean_env)
    assert {r["status"] for r in bad} == {"fail"}
    with pytest.raises(AssertionError):
        G.assert_fdr_ok(bad, allow_vacuous=True)
