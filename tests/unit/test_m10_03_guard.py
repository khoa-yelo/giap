"""Guard: no table derived from the DIA-NN reports is tracked (M10-03 criterion 9).

Spec: docs/specs/M10-03-diann-convert.md, "What is committed (and what is not)"; Khoa's rule of
2026-10-23 (nothing derived from the DIA-NN reports is committed), read by default (Q-K1) as
covering the converted tables. Extends M10-02's recorded-diann guard.

These tests read ``git ls-files`` and pass while nothing forbidden is tracked; they exist to fail
the commit that would add a converted table.
"""

from __future__ import annotations

import fnmatch
import subprocess
from pathlib import PurePosixPath

from helpers.slurm_sandbox import REPO

INPUTS_ALLOWED = {"pxd078459_files.tsv", "acquisition.json", "README.md"}
FORBIDDEN = ("a[1-3]*.tsv*", "*_abundance.tsv*", "*.convert.json", "*.convert.log")


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def test_pilot_inputs_hold_only_the_m10_01_files() -> None:
    names = {
        p.removeprefix("validation/pilot/inputs/") for p in _tracked("validation/pilot/inputs/")
    }
    assert names <= INPUTS_ALLOWED, (
        f"derived tables under validation/pilot/inputs/: {names - INPUTS_ALLOWED}"
    )


def test_no_converted_table_is_tracked_anywhere() -> None:
    bad = [
        p
        for p in _tracked()
        if any(fnmatch.fnmatchcase(PurePosixPath(p).name, pat) for pat in FORBIDDEN)
    ]
    assert not bad, f"converted DIA-NN tables are never committed (Khoa 2026-10-23): {bad}"


def test_recorded_diann_holds_only_the_provenance_note() -> None:
    files = _tracked("tests/data/recorded/diann/")
    extra = [p for p in files if p != "tests/data/recorded/diann/2.3.1/PROVENANCE.md"]
    assert not extra, f"only 2.3.1/PROVENANCE.md may be tracked there: {extra}"


def test_the_patterns_catch_the_job_outputs() -> None:
    """Self-check of the patterns against the names ``diann_convert.sbatch`` writes."""
    for name in (
        "a1.tsv.gz",
        "a2.tsv.gz",
        "a1_abundance.tsv.gz",
        "a2.convert.json",
        "a1.convert.log",
    ):
        assert any(fnmatch.fnmatchcase(name, pat) for pat in FORBIDDEN), name
