"""Gate M4 evidence (GIPPA_SPEC.md §13 row M4: "circular-shift tests calibrated; gene-level FDR
calibrated (picked, MAYU)"; M4-20 acceptance "Gate evidence collected in PROGRESS.md").

M4-20 makes no new statistical claim (its methods note is "None new"): the gate is the ledger of
calibration records the two M4 suites must have produced at the nightly full replicate counts
(``tolerances.yaml`` ``replicates_full``: gene_fdr 300, regions 1,000).

``REQUIRED`` names, per suite, the ``_criteria`` streams and levels the gate needs:

* ``gene_fdr``: criterion ``fdr`` at q in {0.01, 0.05, 0.10} for picked target-decoy
  (``gene_fdr.picked.pi<pi>``), MAYU with PEP (``gene_fdr.mayu_pep.pi<pi>``) and MAYU without PEP
  (``gene_fdr.mayu_nopep.pi<pi>``), pi in {0.01, 0.05}
  (tests/calibration/test_detection_tiers_gene_fdr.py);
* ``regions``: criterion ``level`` at alpha in {0.01, 0.05, 0.10} for the pooled chromosome
  circular-shift check of each statistic in the full scenarios, the plasmid AMR set at
  {0.05, 0.10} (``PLASMID_LEVELS``), and criterion ``fdr`` at q in {0.01, 0.05, 0.10} for every
  (region type, statistic) family (tests/calibration/test_region_tests.py, lead decision 7).

``test_gate_ledger_matches_the_calibration_suites`` keeps that ledger from drifting away from the
stream names the suites actually build (it needs no logs and runs on every commit). The evidence
check reads the JSON-lines logs of a nightly full run -- ``$GIPPA_CALIBRATION_LOG_DIR`` or
``local/logs/calibration`` -- and fails when a required record is missing, was run at fewer than
``replicates_full`` replicates, or is a confirmed failure.

**Guard (lead decision 12).** The *fast* calibration suite writes into the same directory, so "the
directory exists" is not evidence of a nightly run: without the guard the check failed in every
full fast-suite run, CI included, demanding full-R records a fast run never writes. The evidence
check therefore skips unless the suite's log holds at least one record at ``replicates_full`` --
i.e. unless a nightly full run really happened -- or ``GIPPA_M4_GATE=1`` is set, which makes a
missing nightly a failure (that is how the gate is collected on purpose). The gate itself stays
the nightly full run (decision 6).
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
sys.path.insert(0, str(REPO / "tests" / "calibration"))
sys.path.insert(0, str(REPO / "tests" / "helpers"))

LEVELS = (0.01, 0.05, 0.10)
PLASMID_LEVELS = (0.05, 0.10)
PIS = (0.01, 0.05)
REGION_TYPES = ("prophage", "bgc", "crispr", "amr")
STATS = ("o_minus_e", "hl_shift")
#: full-mode region scenarios (test_region_tests.FULL / FULL_CHECKS, lead decision 7)
FULL_RHOS = (0.0, 0.5, 0.9)
#: statuses that are evidence of a passed check (``_gate.NOT_FAILED`` minus the no-verdict ones)
PASSED = frozenset({"pass", "unconfirmed"})


def _gene_fdr_streams() -> list[tuple[str, str, float]]:
    return [
        (f"gene_fdr.{variant}.pi{pi}", "fdr", q)
        for variant in ("picked", "mayu_pep", "mayu_nopep")
        for pi in PIS
        for q in LEVELS
    ]


def _region_streams() -> list[tuple[str, str, float]]:
    out: list[tuple[str, str, float]] = []
    hl = [(f"regions.region_tests.null_full_rho{r}", "hl_shift") for r in FULL_RHOS]
    for base, stat in [*hl, ("regions.region_tests.null_full_rho0.5", "o_minus_e")]:
        out += [(f"{base}.chr1_pooled.{stat}", "level", a) for a in LEVELS]
        out += [(f"{base}.amr@p1.{stat}", "level", a) for a in PLASMID_LEVELS]
    out += [
        (f"regions.region_tests.fdr_full.{typ}.{stat}", "fdr", q)
        for typ in REGION_TYPES
        for stat in STATS
        for q in LEVELS
    ]
    return out


#: suite -> required (stream, criterion, level) records of the nightly full run
REQUIRED: dict[str, list[tuple[str, str, float]]] = {
    "gene_fdr": _gene_fdr_streams(),
    "regions": _region_streams(),
}


def _tolerances() -> dict[str, Any]:
    with TOLERANCES.open() as fh:
        data: dict[str, Any] = yaml.safe_load(fh)
    return data


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


def _wanted(items: Iterable[tuple[str, str, float]]) -> set[tuple[str, str, float]]:
    return {(s, c, round(float(a), 10)) for s, c, a in items}


# --------------------------------------------------------------------------- ledger


def test_gate_ledger_matches_the_calibration_suites() -> None:
    """The ledger above names the streams the M4 suites build (no logs needed)."""
    import region_tests_calib as rc
    import region_tests_data as rd
    import test_detection_tiers_gene_fdr as tiers
    import test_region_tests as regions

    assert tiers.GATED == ("picked", "mayu_pep", "mayu_nopep")
    assert tuple(tiers.PIS) == PIS
    assert tuple(rd.REGION_TYPES) == REGION_TYPES and tuple(rd.STATS) == STATS
    assert tuple(rc.PLASMID_LEVELS) == PLASMID_LEVELS
    assert tuple(_tolerances()["levels"]) == LEVELS
    full = {sc.stream for sc in regions.FULL}
    assert full == {f"regions.region_tests.null_full_rho{r}" for r in FULL_RHOS}
    pooled = {f"{sc.stream}.chr1_pooled.{stat}" for sc, stat in regions.FULL_CHECKS}
    ledger = {s for s, c, _ in REQUIRED["regions"] if c == "level" and "chr1_pooled" in s}
    assert ledger == pooled


def test_every_gate_level_is_a_tolerances_level() -> None:
    levels = {round(float(x), 10) for x in _tolerances()["levels"]}
    used = {a for items in REQUIRED.values() for _s, _c, a in items}
    assert used <= levels, sorted(used - levels)


# --------------------------------------------------------------------------- recorded evidence


def required_replicates(suite: str) -> int:
    return int(_tolerances()["replicates_full"][suite])


def demanded() -> bool:
    """``GIPPA_M4_GATE=1``: a missing nightly is a failure, not a skip (gate collection)."""
    return os.environ.get("GIPPA_M4_GATE") == "1"


def has_full_run(records: Iterable[dict[str, Any]], full: int) -> bool:
    """True when some record was produced at the nightly full replicate count (decision 12:
    the fast suite writes into the same directory, so the directory alone proves nothing)."""
    return any(int(r.get("R", 0)) >= full for r in records)


@pytest.fixture(scope="module")
def evidence() -> dict[str, list[dict[str, Any]]]:
    return {suite: _records(suite) for suite in REQUIRED}


@pytest.mark.parametrize("suite", sorted(REQUIRED))
def test_m4_gate_evidence_is_recorded(
    evidence: dict[str, list[dict[str, Any]]], suite: str
) -> None:
    """Every required record exists, at ``replicates_full`` replicates, and none failed."""
    full = required_replicates(suite)
    records = evidence[suite]
    if not demanded() and not has_full_run(records, full):
        pytest.skip(
            f"no nightly full record (R >= {full}) for suite {suite} in {log_dir()}:"
            " run scripts/slurm/calibration_full.sbatch, point GIPPA_CALIBRATION_LOG_DIR at its"
            " logs, or set GIPPA_M4_GATE=1 to make the missing nightly a failure"
        )
    by_key: dict[tuple[str, str, float], list[dict[str, Any]]] = {}
    for rec in records:
        by_key.setdefault(_key(rec), []).append(rec)
    want = _wanted(REQUIRED[suite])
    missing = sorted(want - set(by_key))
    assert not missing, f"{suite}: no record in {log_dir()} for {missing[:10]} ({len(missing)})"
    small, failed = [], []
    for key in sorted(want):
        runs = [r for r in by_key[key] if int(r.get("R", 0)) >= full]
        if not runs:
            small.append((key, max(int(r.get("R", 0)) for r in by_key[key])))
            continue
        if not any(str(r.get("status")) in PASSED for r in runs):
            failed.append((key, sorted({str(r.get("status")) for r in runs})))
    assert not small, f"{suite}: fast-suite only (R < {full}): {small[:10]}"
    assert not failed, f"{suite}: not passed: {failed[:10]}"


@pytest.mark.parametrize("suite", sorted(REQUIRED))
def test_the_guard_does_not_accept_fast_suite_records(suite: str) -> None:
    """Decision 12: a log directory the fast suite created is not gate evidence."""
    tol = _tolerances()
    fast, full = int(tol["replicates_fast"][suite]), required_replicates(suite)
    assert fast < full, (suite, fast, full)
    assert not has_full_run([{"R": fast}, {"R": fast}], full)
    assert not has_full_run([], full)
    assert has_full_run([{"R": fast}, {"R": full}], full)


def test_the_gate_can_be_demanded_with_an_env_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GIPPA_M4_GATE", raising=False)
    assert not demanded()
    monkeypatch.setenv("GIPPA_M4_GATE", "1")
    assert demanded()
