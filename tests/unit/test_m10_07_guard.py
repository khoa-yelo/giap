"""Guard and methods note (M10-07 acceptance criteria 10 and 11, as amended by the binding
"Statistician decisions", "Other required changes").

Spec: docs/specs/M10-07-free-entrapment.md. Data rule (Khoa, 2026-10-24): aggregates may be
committed, row-level DIA-NN content may not; ``entrapment_discoveries.tsv`` is row-level and is
never tracked. ``entrapment.json`` (aggregates) is the lead's decision at report time, so this
guard does not forbid it; the content guard (no discovered peptide string in it) is in
``test_m10_07_cli.py``. Extends M10-05's guard.

The ``git ls-files`` tests pass while nothing forbidden is tracked; they exist to fail the
commit that would add one. The methods-note tests fail until ``docs/methods/free-entrapment.md``
exists.
"""

from __future__ import annotations

import subprocess
from pathlib import PurePosixPath

from helpers.entrapment_m10_07 import METHODS_NOTE, REPO

DOI = "10.1038/s41592-025-02719-x"
EQ1 = "N_E (1 + 1/r) / (N_T + N_E)"
EQ2 = "N_E / (N_T + N_E)"


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def test_no_discoveries_table_is_tracked() -> None:
    bad = [p for p in _tracked() if PurePosixPath(p).name.startswith("entrapment_discoveries")]
    assert not bad, f"row-level DIA-NN content is never committed: {bad}"


def test_nothing_under_validation_pilot_metrics() -> None:
    assert not _tracked("validation/pilot/metrics/")


def test_the_entrapment_script_and_sbatch_are_tracked() -> None:
    files = set(_tracked("validation/pilot/", "scripts/slurm/"))
    assert "validation/pilot/entrapment.py" in files
    assert "scripts/slurm/entrapment.sbatch" in files


def test_sbatch_request() -> None:
    """Criterion 12's request: dev, 4 CPUs, 16 GB, 1 h, job name giap-entrapment (R14)."""
    text = (REPO / "scripts" / "slurm" / "entrapment.sbatch").read_text(encoding="utf-8")
    assert "--job-name=giap-entrapment" in text
    assert "--cpus-per-task=4" in text
    assert "--mem=16G" in text
    assert any(t in text for t in ("--time=01:00:00", "--time=1:00:00", "--time=60"))
    assert "--qos=dev" in text  # FarmShare: dev is a QOS, not a partition (GIAP_SPEC §3)


def test_methods_note_exists_and_cites_the_doi() -> None:
    assert METHODS_NOTE.is_file(), f"{METHODS_NOTE.relative_to(REPO)} does not exist"
    text = METHODS_NOTE.read_text(encoding="utf-8")
    assert DOI in text


def test_methods_note_has_equations_1_and_2() -> None:
    assert METHODS_NOTE.is_file(), f"{METHODS_NOTE.relative_to(REPO)} does not exist"
    text = " ".join(METHODS_NOTE.read_text(encoding="utf-8").split())
    assert EQ1 in text
    assert EQ2 in text


def test_methods_note_states_the_s3_bound_and_the_conditional_lower_bound() -> None:
    assert METHODS_NOTE.is_file(), f"{METHODS_NOTE.relative_to(REPO)} does not exist"
    low = " ".join(METHODS_NOTE.read_text(encoding="utf-8").lower().split())
    assert "ucb95" in low
    assert "equal chance" in low
    assert "realised" in low or "realized" in low
    assert "precursor" in low and ("not_assessed" in low or "not assessed" in low)
    assert "sample method" in low  # why (3) is not used


def test_methods_note_carries_no_pilot_numbers() -> None:
    """House style of target-decoy.md: no measured pilot number (job IDs, r, FDP values)."""
    assert METHODS_NOTE.is_file(), f"{METHODS_NOTE.relative_to(REPO)} does not exist"
    text = METHODS_NOTE.read_text(encoding="utf-8")
    assert "1879955" not in text  # the A1 search job
    assert "Measured" not in text
