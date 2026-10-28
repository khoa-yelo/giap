"""A retry must change the request (M7-15; review #278 B1, lead decision 11; GIPPA_SPEC.md §6.3).

Mini Nextflow runs (as in test_nf_config.py example 7) with the repository's
``workflow/nextflow.config`` under each profile (base = no profile, farmshare,
farmshare_local) and an override that runs every task on the local executor; only ``exit N``
runs, nothing is submitted to Slurm. The trace records the requested memory and time per attempt.

* **143 (time):** under farmshare ``time`` grows with ``task.attempt``: 12 h x attempt on
  process_medium/high, 2 h x attempt on process_single/low.
* **137 (memory):** a retry happens only while the memory request still grows under the
  profile's cap (farmshare 64 GB, farmshare_local allocation - 2 GB = 30 GB here, base
  ``params.memory`` = 64 GB here); otherwise the task terminates.
* **Every profile:** two consecutive attempts of one task never ask for the same memory and
  time.
"""

from __future__ import annotations

import csv
import itertools
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

WF = Path(__file__).resolve().parents[2] / "workflow"
LABELS = {"single": 2, "low": 6, "medium": 16, "high": 48}  # base.config GB x attempt
CAP_GB = {"base": 64, "farmshare": 64, "farmshare_local": 30}
FARMSHARE_HOURS = {"single": 2, "low": 2, "medium": 12, "high": 12}  # x attempt

_OVERRIDE = """
process.executor = 'local'
process.scratch = false
executor {
    name = 'local'
    cpus = 64
    memory = '512 GB'
}
trace {
    enabled = true
    overwrite = true
    file = 'trace.tsv'
    fields = 'name,status,exit,attempt,cpus,memory,time'
}
"""


def _env(profile: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    env["NXF_ANSI_LOG"] = "false"
    if profile == "farmshare_local":
        env.update(SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE=str((30 + 2) * 1024))
    return env


def _process(name: str, label: str, script: str) -> str:
    return f"process {name} {{\n    label 'process_{label}'\n    script:\n    \"{script}\"\n}}\n"


def _run(tmp_path: Path, profile: str, procs: dict[str, tuple[str, str]]) -> tuple[int, list]:
    """``procs``: name -> (label, shell script). Returns the exit code and the trace rows."""
    body = "".join(_process(n, lab, s) for n, (lab, s) in procs.items())
    body += "workflow {\n" + "".join(f"    {n}()\n" for n in procs) + "}\n"
    (tmp_path / "main.nf").write_text(body)
    (tmp_path / "override.config").write_text(_OVERRIDE)
    params = {"threads": 64, "memory": f"{CAP_GB[profile]} GB"}
    (tmp_path / "params.json").write_text(json.dumps(params))
    cmd = ["nextflow", "-q", "run", "main.nf", "-c", str(WF / "nextflow.config")]
    cmd += ["-c", "override.config", "-params-file", "params.json", "-work-dir", "work"]
    if profile != "base":
        cmd += ["-profile", profile]
    res = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, env=_env(profile), timeout=480
    )
    trace = tmp_path / "trace.tsv"
    assert trace.is_file(), f"no trace:\n{res.stdout}\n{res.stderr}"
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return res.returncode, rows


def _gb(value: str) -> float:
    m = re.fullmatch(r"([\d.]+)\s*([KMGT]?B)", value.strip())
    assert m, f"not a memory value: {value!r}"
    scale = {"B": 2**-30, "KB": 2**-20, "MB": 2**-10, "GB": 1, "TB": 2**10}[m.group(2)]
    return float(m.group(1)) * scale


def _hours(value: str) -> float | None:
    """Trace ``time`` ('12h', '1d 12h', '-' when unset) in hours."""
    v = value.strip()
    if v in ("-", ""):
        return None
    units = {"d": 24.0, "h": 1.0, "m": 1 / 60, "s": 1 / 3600, "ms": 1 / 3_600_000}
    parts = re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", v)
    assert parts and "".join(f"{n}{u}" for n, u in parts) == v.replace(" ", ""), v
    return sum(float(n) * units[u] for n, u in parts)


def _attempts(rows: list[dict[str, str]], name: str) -> list[dict[str, str]]:
    return sorted((r for r in rows if r["name"] == name), key=lambda r: int(r["attempt"]))


def _request(row: dict[str, str]) -> tuple[float, float | None]:
    return _gb(row["memory"]), _hours(row["time"])


def _assert_each_retry_changes_the_request(name: str, mine: list[dict[str, str]]) -> None:
    for a, b in itertools.pairwise(mine):
        assert _request(a) != _request(b), (
            f"{name}: attempt {b['attempt']} repeats attempt {a['attempt']}: {_request(a)}"
        )


def _expected_gb(profile: str, label: str, attempt: int) -> float:
    return min(LABELS[label] * attempt, CAP_GB[profile])


# Tasks whose retries are all allowed: exit CODE on attempts 1 and 2, succeed on attempt 3.
# A 137 retry needs the memory to grow; a 143 retry under farmshare grows the time.
RETRIED = {
    "base": [("single", 137), ("low", 137), ("medium", 137), ("single", 143), ("low", 143),
             ("medium", 143)],
    "farmshare": [("single", 137), ("low", 137), ("medium", 137), ("single", 143),
                  ("low", 143), ("medium", 143), ("high", 143)],
    "farmshare_local": [("single", 137), ("low", 137), ("single", 143), ("low", 143)],
}  # fmt: skip


@pytest.mark.parametrize("profile", sorted(RETRIED))
def test_retries_grow_the_request(tmp_path: Path, profile: str) -> None:
    procs = {
        f"{label.upper()}_{code}": (label, f"exit ${{task.attempt < 3 ? {code} : 0}}")
        for label, code in RETRIED[profile]
    }
    code, rows = _run(tmp_path, profile, procs)
    assert code == 0, f"a retry that grows the request was refused: {rows}"
    for name in procs:
        label, exit_code = name.split("_")[0].lower(), name.split("_")[1]
        mine = _attempts(rows, name)
        assert [r["attempt"] for r in mine] == ["1", "2", "3"], (name, mine)
        assert [r["status"] for r in mine] == ["FAILED", "FAILED", "COMPLETED"], name
        assert [r["exit"] for r in mine] == [exit_code, exit_code, "0"], name
        gbs = [_gb(r["memory"]) for r in mine]
        assert gbs == [_expected_gb(profile, label, a) for a in (1, 2, 3)], (name, gbs)
        _assert_each_retry_changes_the_request(name, mine)
        if profile == "farmshare":
            hours = [_hours(r["time"]) for r in mine]
            assert hours == [FARMSHARE_HOURS[label] * a for a in (1, 2, 3)], (name, hours)


# 137 when the next memory request would not grow (already at the cap): terminate.
CAPPED_137 = [
    ("base", "high", ["1", "2"]),  # 48, then 64 = cap: no third attempt at 64
    ("farmshare", "high", ["1", "2"]),
    ("farmshare_local", "medium", ["1", "2"]),  # 16, then 30 = cap
    ("farmshare_local", "high", ["1"]),  # 48 is clamped to 30 already on attempt 1
]


@pytest.mark.parametrize(
    ("profile", "label", "want"), CAPPED_137, ids=[f"{p}-{lab}" for p, lab, _ in CAPPED_137]
)
def test_137_at_the_memory_cap_terminates(
    tmp_path: Path, profile: str, label: str, want: list[str]
) -> None:
    code, rows = _run(tmp_path, profile, {"OOM": (label, "exit 137")})
    assert code != 0
    mine = _attempts(rows, "OOM")
    assert [r["attempt"] for r in mine] == want, [(r["attempt"], r["memory"]) for r in mine]
    assert all((r["status"], r["exit"]) == ("FAILED", "137") for r in mine)
    gbs = [_gb(r["memory"]) for r in mine]
    assert gbs == [_expected_gb(profile, label, int(a)) for a in want]
    assert all(b > a for a, b in itertools.pairwise(gbs)), gbs


# 143 where the memory stops growing under the cap: a retry is allowed only if the time grows.
CAPPED_143 = [("base", "high"), ("farmshare", "high"), ("farmshare_local", "medium"),
              ("farmshare_local", "high")]  # fmt: skip


@pytest.mark.parametrize(
    ("profile", "label"), CAPPED_143, ids=[f"{p}-{lab}" for p, lab in CAPPED_143]
)
def test_143_never_repeats_an_identical_request(tmp_path: Path, profile: str, label: str) -> None:
    code, rows = _run(tmp_path, profile, {"TERM": (label, "exit 143")})
    assert code != 0
    mine = _attempts(rows, "TERM")
    assert 1 <= len(mine) <= 3
    assert all((r["status"], r["exit"]) == ("FAILED", "143") for r in mine)
    _assert_each_retry_changes_the_request("TERM", mine)
