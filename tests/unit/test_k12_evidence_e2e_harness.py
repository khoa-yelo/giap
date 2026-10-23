"""The e2e harness tests/e2e/test_k12_evidence.py against hand-built runs (MVP-04 part 2).

The e2e test itself only runs weekly on a real outdir; here it is run in a subprocess on the
all-pass fixture of tests/helpers/mvp04_summary_data.py (after ``scripts/dev/mvp04_summary.py``)
to show that it passes a passing run, fails a failing one, and skips without ``$GIAP_MVP04_OUT``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from helpers import mvp04_summary_data as S

REPO = Path(__file__).resolve().parents[2]
E2E = REPO / "tests" / "e2e" / "test_k12_evidence.py"
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"


def _pytest(env_out: Path | None) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env.pop("GIAP_MVP04_OUT", None)
    if env_out is not None:
        env["GIAP_MVP04_OUT"] = str(env_out)
    return subprocess.run(
        [sys.executable, "-m", "pytest", str(E2E), "-q", "-p", "no:cacheprovider", "-p",
         "no:randomly", "-n", "0", "-m", "e2e", "-W", "default"],
        cwd=REPO, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip


def _summarise(d: S.Dirs) -> None:
    subprocess.run([sys.executable, str(SCRIPT), str(d.out), str(d.inputs), "--inputs-again",
                    str(d.again)], check=True, capture_output=True)  # fmt: skip


def test_skips_without_the_env_var() -> None:
    res = _pytest(None)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "skipped" in res.stdout and "passed" not in res.stdout, res.stdout


def test_fails_when_the_summary_is_missing(tmp_path: Path) -> None:
    res = _pytest(tmp_path)
    assert res.returncode == 1, res.stdout + res.stderr


def test_passes_the_all_pass_run(tmp_path: Path) -> None:
    d = S.build(tmp_path)
    _summarise(d)
    res = _pytest(d.out)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "failed" not in res.stdout and "skipped" not in res.stdout, res.stdout


def test_fails_a_run_with_a_failing_gate(tmp_path: Path) -> None:
    """Example 4: the ``mayu`` tier fails gate 2."""
    d = S.build(tmp_path)

    def tier(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["detection_tier"] = "mayu"
        return rows

    S.edit_tsv(d.out / "tables" / "genes.tsv", tier)
    _summarise(d)
    res = _pytest(d.out)
    assert res.returncode == 1, res.stdout + res.stderr
    assert "test_gate2_tier_and_warnings" in res.stdout, res.stdout


def test_outside_expected_is_a_warning_not_a_failure(tmp_path: Path) -> None:
    d = S.build(tmp_path)
    S.set_tier_metric(d.out, "decoy_paired_share", 0.9)
    _summarise(d)
    res = _pytest(d.out)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "outside" in res.stdout, res.stdout
