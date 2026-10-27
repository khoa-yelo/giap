"""Gate M5 evidence (GIPPA_SPEC.md §13 row M5: "island FWER <= 5 % under simulated nulls; power
curve recorded"; M5-11 acceptance "Gate evidence in PROGRESS.md" and "Benchmarks").

M5-11 makes no new statistical claim; the gate is the ledger of records the two M5 suites must
have produced at the nightly full replicate counts (``tolerances.yaml`` ``replicates_full``:
islands 500, autocorrelation 1,000), copied from tests/integration/test_m4_gate.py:

* ``islands`` -- criterion ``fwer`` at every tolerances level for ``islands.rho0``,
  ``islands.rho05``, ``islands.rho09`` (tests/calibration/test_island_fwer.py
  ``test_full_fwer_gate``: 4,500-gene chromosome + 60-gene plasmid per rho);
* ``autocorrelation`` -- criterion ``level`` at lags 1 and 20 and criterion ``fwer`` of the BH
  family over the lags, for every (statistic, null) of the label-null scenario and the TU-block
  null of the TU-null scenario at 4,500 genes (tests/calibration/
  test_autocorrelation_calibration.py, ``full`` mode).

The evidence check reads the JSON-lines logs (``$GIPPA_CALIBRATION_LOG_DIR`` or
``local/logs/calibration``) and skips unless a record at ``replicates_full`` exists (the fast suite
writes into the same directory, M4-20 decision 12) or ``GIPPA_M5_GATE=1`` demands the evidence.
Lead decision 39 / Q2: the gate run is one green full calibration run on main as a Slurm job.

The benchmark entry ``context_m5`` of ``benchmarks/baseline.json`` (benchmarks/
test_context_m5.py) must be recorded with wall time, peak memory and a budget <= 120 s on 4 CPUs.
"""

from __future__ import annotations

import json
import os
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
BASELINE = REPO / "benchmarks" / "baseline.json"
sys.path.insert(0, str(REPO / "tests" / "calibration"))
sys.path.insert(0, str(REPO / "tests" / "helpers"))

ISLAND_RHOS = {0.0: "islands.rho0", 0.5: "islands.rho05", 0.9: "islands.rho09"}
AUTOCORR_GENES = 4500
LAGS = (1, 20)
PASSED = frozenset({"pass", "unconfirmed"})
#: M5-11 acceptance: budget < 2 min on 4 CPUs
M5_BUDGET_S, M5_CPUS = 120.0, 4


def _tolerances() -> dict[str, Any]:
    with TOLERANCES.open() as fh:
        data: dict[str, Any] = yaml.safe_load(fh)
    return data


def _levels() -> list[float]:
    return [float(x) for x in _tolerances()["levels"]]


def _island_streams() -> list[tuple[str, str, float]]:
    return [(s, "fwer", a) for s in ISLAND_RHOS.values() for a in _levels()]


def _autocorr_streams() -> list[tuple[str, str, float]]:
    import autocorr_calib as acal
    import autocorr_data as ad

    label, tu = acal.label_null(AUTOCORR_GENES), acal.tu_null(AUTOCORR_GENES)
    combos = [(label, s, m) for s in ad.STATISTICS for m in ad.NULLS]
    combos += [(tu, s, "tu_block") for s in ad.STATISTICS]
    out: list[tuple[str, str, float]] = []
    for sc, stat, null in combos:
        base = f"{sc.stream}.{stat}.{null}"
        out += [(f"{base}.lag{k}", "level", a) for k in LAGS for a in _levels()]
        out += [(f"{base}.bh_family", "fwer", a) for a in _levels()]
    return out


def required() -> dict[str, list[tuple[str, str, float]]]:
    return {"islands": _island_streams(), "autocorrelation": _autocorr_streams()}


def log_dir() -> Path:
    return Path(os.environ.get("GIPPA_CALIBRATION_LOG_DIR", REPO / "local" / "logs" / "calibration"))


def _records(suite: str) -> list[dict[str, Any]]:
    path = log_dir() / f"{suite}.jsonl"
    if not path.is_file():
        return []
    with path.open(encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def _key(rec: dict[str, Any]) -> tuple[str, str, float]:
    return str(rec["stream"]), str(rec["criterion"]), round(float(rec["level"]), 10)


def demanded() -> bool:
    return os.environ.get("GIPPA_M5_GATE") == "1"


def has_full_run(records: Iterable[dict[str, Any]], full: int) -> bool:
    return any(int(r.get("R", 0)) >= full for r in records)


# --------------------------------------------------------------------------- ledger


def _constants(path: Path, *names: str) -> dict[str, Any]:
    """Module-level literal assignments of a calibration test file, read with ``ast`` (importing
    it outside tests/calibration fails on its conftest-registered markers)."""
    import ast

    out: dict[str, Any] = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign):
            targets = [t for t in node.targets if isinstance(t, ast.Name | ast.Tuple)]
            for t in targets:
                if isinstance(t, ast.Name) and t.id in names:
                    out[t.id] = ast.literal_eval(node.value)
                elif isinstance(t, ast.Tuple) and isinstance(node.value, ast.Tuple):
                    for el, val in zip(t.elts, node.value.elts, strict=True):
                        if isinstance(el, ast.Name) and el.id in names:
                            out[el.id] = ast.literal_eval(val)
    return out


def test_island_ledger_matches_the_gate_suite() -> None:
    """The island streams and rhos are the ones ``test_full_fwer_gate`` logs."""
    c = _constants(REPO / "tests" / "calibration" / "test_island_fwer.py",
                   "SUITE", "FULL_RHOS", "STREAM", "FULL_GENES")  # fmt: skip
    assert c["SUITE"] == "islands"
    assert tuple(c["FULL_RHOS"]) == tuple(ISLAND_RHOS)
    assert c["STREAM"] == ISLAND_RHOS
    assert c["FULL_GENES"] == AUTOCORR_GENES


def test_autocorrelation_ledger_matches_the_suite() -> None:
    import autocorr_calib as acal

    c = _constants(REPO / "tests" / "calibration" / "test_autocorrelation_calibration.py",
                   "FULL_GENES")  # fmt: skip
    assert acal.SUITE == "autocorrelation"
    assert c["FULL_GENES"] == AUTOCORR_GENES
    streams = {s for s, _c, _a in _autocorr_streams()}
    assert len(streams) == (4 + 2) * 3, sorted(streams)


# --------------------------------------------------------------------------- recorded evidence


@pytest.mark.parametrize("suite", ["islands", "autocorrelation"])
def test_m5_gate_evidence_is_recorded(suite: str) -> None:
    """Every required record exists at ``replicates_full`` and passed."""
    full = int(_tolerances()["replicates_full"][suite])
    records = _records(suite)
    if not demanded() and not has_full_run(records, full):
        pytest.skip(
            f"no nightly full record (R >= {full}) for suite {suite} in {log_dir()}: submit"
            " scripts/slurm/calibration_full.sbatch --array=4,5 on main and point"
            " GIPPA_CALIBRATION_LOG_DIR at its records, or set GIPPA_M5_GATE=1"
        )
    by_key: dict[tuple[str, str, float], list[dict[str, Any]]] = {}
    for rec in records:
        by_key.setdefault(_key(rec), []).append(rec)
    want = {(s, c, round(a, 10)) for s, c, a in required()[suite]}
    missing = sorted(want - set(by_key))
    assert not missing, f"{suite}: no record for {missing[:10]} ({len(missing)})"
    small, failed = [], []
    for key in sorted(want):
        runs = [r for r in by_key[key] if int(r.get("R", 0)) >= full]
        if not runs:
            small.append(key)
        elif not any(str(r.get("status")) in PASSED for r in runs):
            failed.append((key, sorted({str(r.get("status")) for r in runs})))
    assert not small, f"{suite}: fast-suite records only: {small[:10]}"
    assert not failed, f"{suite}: not passed: {failed[:10]}"


def test_the_gate_can_be_demanded_with_an_env_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GIPPA_M5_GATE", raising=False)
    assert not demanded()
    monkeypatch.setenv("GIPPA_M5_GATE", "1")
    assert demanded()
    tol = _tolerances()
    for suite in ("islands", "autocorrelation"):
        fast, full = int(tol["replicates_fast"][suite]), int(tol["replicates_full"][suite])
        assert not has_full_run([{"R": fast}], full) and has_full_run([{"R": full}], full)


# --------------------------------------------------------------------------- benchmark record


def test_context_m5_benchmark_is_recorded() -> None:
    """``benchmarks/baseline.json`` ``context_m5``: wall time and peak memory within budget."""
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    assert "context_m5" in base, "run benchmarks/test_context_m5.py on Slurm and record it"
    rec = base["context_m5"]
    assert rec.get("benchmark", "").startswith("benchmarks/test_context_m5.py"), rec
    assert float(rec["budget_s"]) <= M5_BUDGET_S
    assert 0.0 < float(rec["secs"]) < float(rec["budget_s"]), rec
    assert float(rec["peak_mb"]) > 0.0, rec
    assert int(rec["cpus"]) == M5_CPUS
    assert int(rec["n_genes"]) == 4500 and int(rec["bootstrap"]) == 1000
    for k in ("slurm_job", "commit", "date", "node"):
        assert rec.get(k), f"context_m5.{k} not recorded"
