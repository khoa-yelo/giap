"""#514 and #515 under real Nextflow, local executor, no container and no real tool.

The unit tests (tests/unit/test_threads_host_cap.py, tests/unit/test_antismash_image_databases.py)
pin the interface; these check what the real engine does with it:

* #514: a capped ``params.threads`` clamps ``process_high``, ``TERMINATORNET`` and the
  ``DETERMINANTS_STEP`` selector (a throw-away pipeline, as in tests/integration/test_nf_config.py;
  the trace records the CPUs each task got). ``resourceLimits`` follows ``params.threads`` under
  each local profile, and ``uses_local_executor`` agrees with the ``process.executor`` Nextflow
  resolves for every bundled profile.
* #515: ``workflow/tests/real/antismash.nf`` runs the real module with a fake ``antismash`` on
  ``PATH`` that records its argv, so ``.command.sh`` runs as in production. The ``-stub-run`` path
  is run as well.
"""

from __future__ import annotations

import csv
import json
import os
import re
import stat
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
REAL = WF / "tests" / "real"
THREADS = 2  # below any CI runner or allocation; the requests below ask for 16

FAKE_ANTISMASH = """#!/bin/bash
if [ "${1:-}" = --version ]; then echo 'antiSMASH 8.0.4'; exit 0; fi
printf '%s\\n' "$@" > "$FAKE_ARGV"
prev=
for a in "$@"; do
  if [ "$prev" = --databases ] && [ ! -d "$a/nrps_pks" ]; then
    echo "ERROR: No matching database in location $a" >&2
    exit 1
  fi
  prev="$a"
done
mkdir -p antismash
echo '{"records": []}' > antismash/genome.json
"""

MINI = """
process HIGH {
    label 'process_high'
    script: 'true'
}
process TERMINATORNET {
    script: 'true'
}
process DETERMINANTS_STEP {
    script: 'true'
}
workflow DETERMINANTS {
    DETERMINANTS_STEP()
}
workflow {
    HIGH()
    TERMINATORNET()
    DETERMINANTS()
}
"""

OVERRIDE = """
process.executor = 'local'
trace.enabled = true
trace.file = 'trace.tsv'
trace.fields = 'name,status,cpus'
trace.overwrite = true
"""


def _env(**extra: str) -> dict[str, str]:
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"}
    env = {k: v for k, v in os.environ.items() if k not in drop}
    env["NXF_ANSI_LOG"] = "false"
    env.update(extra)
    return env


def _flat(profile: str) -> dict[str, str]:
    cmd = ["nextflow", "config", str(WF), "-flat", "-profile", profile]
    res = subprocess.run(cmd, capture_output=True, text=True, env=_env(), timeout=300, check=False)
    assert res.returncode == 0, f"{' '.join(cmd)} failed:\n{res.stdout}\n{res.stderr}"
    pairs = (re.match(r"^([\w.$'\":-]+) = (.*)$", ln) for ln in res.stdout.splitlines())
    return {m[1]: m[2].strip() for m in pairs if m}


@pytest.mark.parametrize("profile", ["standard", "apptainer", "docker", "conda"])
def test_limit_follows_params_threads_under_local_profiles(profile: str) -> None:
    flat = _flat(profile)
    assert flat["process.resourceLimits.cpus"] == flat["params.threads"], flat


def test_every_bundled_profile_resolves_to_the_executor_the_cap_assumes() -> None:
    from gippa.runner.profiles import accepted_profiles, uses_local_executor

    for name in accepted_profiles():
        executor = _flat(name).get("process.executor", "'local'").strip("'\"")
        assert uses_local_executor(name) is (executor == "local"), f"{name}: {executor}"


def test_capped_threads_clamp_the_16_cpu_tasks(tmp_path: Path) -> None:
    """16 CPUs requested, ``params.threads`` 2: the tasks run, at 2 CPUs. Without the limit the
    local executor aborts as in #514 on a host below 16 CPUs."""
    (tmp_path / "main.nf").write_text(MINI, encoding="utf-8")
    (tmp_path / "override.config").write_text(OVERRIDE, encoding="utf-8")
    params = {"threads": THREADS, "memory": "4 GB"}
    (tmp_path / "params.json").write_text(json.dumps(params), encoding="utf-8")
    cmd = ["nextflow", "-q", "run", "main.nf", "-c", str(WF / "nextflow.config")]
    cmd += ["-c", "override.config", "-params-file", "params.json", "-work-dir", "work"]
    res = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, env=_env(),
                         timeout=600, check=False)  # fmt: skip
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    with (tmp_path / "trace.tsv").open(encoding="utf-8") as fh:
        rows = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert {n: r["cpus"] for n, r in rows.items()} == {
        "HIGH": str(THREADS),
        "TERMINATORNET": str(THREADS),
        "DETERMINANTS:DETERMINANTS_STEP": str(THREADS),
    }, rows


def _run_antismash(tmp_path: Path, db: Path, *extra: str) -> tuple[int, str, Path]:
    """The real module through ``workflow/tests/real/antismash.nf``: (exit, output, argv file)."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "antismash"
    exe.write_text(FAKE_ANTISMASH, encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    gbff = tmp_path / "genome.gbff"
    gbff.write_text("LOCUS x\n", encoding="utf-8")
    argv_file = tmp_path / "argv.txt"
    cmd = ["nextflow", "-q", "run", str(REAL / "antismash.nf"), "-c", str(REAL / "real.config")]
    cmd += ["--sample", "genome", "--gbff", str(gbff), "--db_dir", str(db)]
    cmd += ["--threads", str(THREADS), "-work-dir", "work", *extra]
    env = _env(PATH=f"{bindir}{os.pathsep}{os.environ['PATH']}", FAKE_ARGV=str(argv_file))
    res = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, env=env, timeout=600,
                         check=False)  # fmt: skip
    return res.returncode, f"{res.stdout}\n{res.stderr}", argv_file


def _argv(path: Path) -> list[str]:
    assert path.is_file(), "the fake antismash did not run"
    return path.read_text(encoding="utf-8").splitlines()


def test_antismash_fresh_setup_omits_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    db.mkdir()
    code, text, argv = _run_antismash(tmp_path, db)
    assert code == 0, text
    assert "--databases" not in _argv(argv)


def test_antismash_existing_install_keeps_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    (db / "antismash" / "nrps_pks" / "transATor").mkdir(parents=True)
    (db / "antismash" / ".giap-complete").touch()
    code, text, argv = _run_antismash(tmp_path, db)
    assert code == 0, text
    got = _argv(argv)
    assert got[got.index("--databases") + 1] == str(db / "antismash")


def test_antismash_stub_needs_neither_tool_nor_databases(tmp_path: Path) -> None:
    db = tmp_path / "db"
    db.mkdir()
    code, text, argv = _run_antismash(tmp_path, db, "-stub-run")
    assert code == 0, text
    assert not argv.exists(), "the stub ran the tool"
    assert json.loads((tmp_path / "out" / "antismash" / "genome.json").read_text()) == {
        "records": []
    }
