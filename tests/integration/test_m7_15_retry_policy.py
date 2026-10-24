"""The retry policy as Nextflow resolves it, per profile (M7-15; review #278 B3, lead decision 11).

``nextflow config <dir> -flat -profile <p>`` for base (no profile), farmshare and farmshare_local
(and their apptainer combinations): exactly one ``process.errorStrategy`` and one
``process.maxRetries`` (= 2), and no ``process.'withName:…'`` / ``process.'withLabel:…'`` key
sets ``errorStrategy`` or ``maxRetries``. A planted override, including one in a selector that
already holds a closure, makes the check fail (self-test on a copy of the configuration).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

WF = Path(__file__).resolve().parents[2] / "workflow"
PROFILES = [
    None,
    "farmshare",
    "farmshare_local",
    "farmshare,apptainer",
    "farmshare_local,apptainer",
]
SLURM_ENV = {"SLURM_CPUS_PER_TASK": "8", "SLURM_MEM_PER_NODE": "32768"}


def _flat(wf: Path, profile: str | None) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    env.update(NXF_ANSI_LOG="false", **SLURM_ENV)
    cmd = ["nextflow", "config", str(wf), "-flat"] + (["-profile", profile] if profile else [])
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    flat: dict[str, str] = {}
    for line in res.stdout.splitlines():
        m = re.match(r"^(\S+) = (.*)$", line)
        if m:
            flat[m.group(1)] = m.group(2).strip()
    return flat


def policy_problems(flat: dict[str, str]) -> list[str]:
    """Every way the resolved config departs from one process-wide retry policy."""
    problems = []
    for key in flat:
        if re.search(r"\.(errorStrategy|maxRetries)$", key) and key not in (
            "process.errorStrategy",
            "process.maxRetries",
        ):
            problems.append(f"{key} = {flat[key]} overrides the process-wide policy")
    if flat.get("process.maxRetries") != "2":
        problems.append(f"process.maxRetries = {flat.get('process.maxRetries')}, want 2")
    strategy = flat.get("process.errorStrategy", "")
    if not (strategy.startswith("{") and "137" in strategy and "143" in strategy):
        problems.append(f"process.errorStrategy is not the 137/143 closure: {strategy!r}")
    elif "'retry'" not in strategy or "'terminate'" not in strategy:
        problems.append(f"process.errorStrategy lacks retry/terminate: {strategy!r}")
    return problems


@pytest.mark.parametrize("profile", PROFILES, ids=lambda p: p or "base")
def test_one_process_wide_retry_policy(profile: str | None) -> None:
    flat = _flat(WF, profile)
    assert any(k.startswith("process.'withName:") for k in flat), "selectors not in -flat output"
    assert policy_problems(flat) == []


PLANTED = {
    # the review's probe: a selector body that already holds a closure
    "closure_selector_maxretries": (
        "base.config",
        "process { withName: 'VARIANTS_STEP' "
        "{ memory = { 8.GB * task.attempt }; maxRetries = 0 } }",
    ),
    "label_errorstrategy_closure": (
        "farmshare.config",
        "process { withLabel: 'process_high' { cpus = { 16 }; "
        "errorStrategy = { task.attempt < 5 ? 'retry' : 'ignore' } } }",
    ),
    "local_profile_maxretries": (
        "farmshare_local.config",
        "process { withName: 'BAKTA|GENOMAD' { maxRetries = 5 } }",
    ),
    "profile_wide_maxretries": ("farmshare.config", "process { maxRetries = 1 }"),
}


@pytest.mark.parametrize("case", sorted(PLANTED))
def test_the_check_catches_a_planted_override(tmp_path: Path, case: str) -> None:
    conf, extra = PLANTED[case]
    wf = tmp_path / "workflow"
    (wf / "conf").mkdir(parents=True)
    shutil.copy(WF / "nextflow.config", wf / "nextflow.config")
    for p in (WF / "conf").glob("*.config"):
        shutil.copy(p, wf / "conf" / p.name)
    with (wf / "conf" / conf).open("a", encoding="utf-8") as fh:
        fh.write("\n" + extra + "\n")
    profile = {"base.config": "farmshare", "farmshare.config": "farmshare"}.get(
        conf, "farmshare_local"
    )
    assert policy_problems(_flat(wf, profile)), f"{case} not caught"
