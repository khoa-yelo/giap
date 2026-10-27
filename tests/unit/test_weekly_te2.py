"""M7-14 rule 5 and lead decision 2: the TE2 proxy entry in ``scripts/slurm/weekly.sbatch``
(docs/specs/M7-14-te2-proxy-and-id-mixups.md; GIPPA_SPEC.md §12, R14).

Static checks only: nothing is submitted and nothing runs (the script is parsed with ``bash -n``
and read as text). Conventions pinned here (test-author reading, following the TE8 entry):

* The Nissle task (``$S`` = ``nissle1917``) runs TE2 beside TE8, either as
  ``python validation/proxies/te_proxies.py te2 <nissle genome> <nissle Bakta JSON> --out ...``
  or as ``pytest tests/integration/test_te2_proxy.py``, on the committed Nissle genome and the
  recorded full-DB Bakta JSON (``tests/data/...``), never on a live download.
* Its result is reported through ``step`` with a summary naming ``TE2``, so a failing proxy turns
  the weekly line red (degenerate row "TE2 violation | the gate is red").
* A ``te2`` result JSON (CLI route) goes to ``$CI/validation/<date>_te2_<job>.json`` as TE8's does.
* R14 (decision 2): the job stays ``gippa-*``, ≤ 16 CPUs per task, no Slurm-executor profile.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WEEKLY = REPO / "scripts" / "slurm" / "weekly.sbatch"


def _code() -> str:
    keep = [ln for ln in WEEKLY.read_text().splitlines()
            if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")]  # fmt: skip
    return "\n".join(keep)


def _nissle_blocks(code: str) -> list[str]:
    """Bodies of ``if [ "$S" = nissle1917 ]; then ... fi`` blocks (one nesting level)."""
    out = []
    for m in re.finditer(r'if \[ "\$S" = "?nissle1917"? \]; then\n', code):
        body = code[m.end() :]
        depth, lines = 1, []
        for ln in body.splitlines():
            s = ln.strip()
            if re.match(r"^if\b", s):
                depth += 1
            if s == "fi" or s.startswith("fi;") or s.startswith("fi "):
                depth -= 1
                if depth == 0:
                    break
            lines.append(ln)
        out.append("\n".join(lines))
    return out


def _te2_block() -> str:
    blocks = [b for b in _nissle_blocks(_code()) if re.search(r"\bte2\b|test_te2_proxy", b)]
    assert blocks, 'weekly.sbatch has no Nissle-only ("$S" = nissle1917) TE2 entry (M7-14 rule 5)'
    return blocks[0]


def test_weekly_still_parses() -> None:
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash not on PATH")
    r = subprocess.run([bash, "-n", str(WEEKLY)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_te2_runs_on_the_nissle_task_only() -> None:
    block = _te2_block()
    assert "nissle1917.fna.gz" in block
    assert "nissle1917.json.gz" in block
    outside = _code()
    for b in _nissle_blocks(outside):
        outside = outside.replace(b, "")
    assert not re.search(r"te_proxies\.py te2\b|test_te2_proxy", outside), "TE2 outside Nissle"


def test_te2_uses_the_committed_inputs_not_a_download() -> None:
    block = _te2_block()
    assert not re.search(r"\b(curl|wget|fetch_)", block), block
    if "te_proxies.py te2" in block:
        assert re.search(r"tests/data|\$D/genomes", block), block


def test_te2_result_is_a_weekly_step() -> None:
    block = _te2_block()
    assert re.search(r'step\s+"?\$\{?PIPESTATUS\[0\]\}?"?\s+"TE2\b|step\s+\S+\s+"TE2\b', block), (
        "the TE2 exit code must reach `step` with a summary naming TE2"
    )


def test_te2_cli_route_writes_a_dated_json_under_ci_validation() -> None:
    block = _te2_block()
    if "te_proxies.py te2" not in block:
        pytest.skip("TE2 runs through pytest, not the CLI")
    assert "--out" in block
    assert re.search(r"\$CI/validation/\$\(date \+%F\)_te2_\$\{?JOB\}?\.json", block), block


def test_r14_caps_unchanged() -> None:
    text = WEEKLY.read_text()
    name = re.search(r"^#SBATCH\s+--job-name[= ](\S+)", text, re.M)
    cpus = re.search(r"^#SBATCH\s+--cpus-per-task[= ](\d+)", text, re.M)
    assert name and name.group(1).startswith("gippa-")
    assert cpus and int(cpus.group(1)) <= 16
    assert not re.search(r"te2[^\n]*-profile[^\n]*\bfarmshare\b(?!_local)", _te2_block())
