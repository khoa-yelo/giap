"""How Nextflow hands the pinned MINIPROT image to each engine (M11-09b examples 2-3, AC 3-4).

A one-process fixture pipeline in ``tmp_path`` (process ``MINIPROT``, a ``stub:`` that touches its
output) includes the checkout's real ``workflow/conf/containers.config`` and runs with
``-stub-run`` and the local executor. Neither Docker nor Apptainer is needed: a logging fake
``docker`` / ``apptainer`` is first on PATH and runs the task command itself.

- (a) ``docker.enabled``: the logged ``docker run`` argv holds the bare reference and no word
  starting with ``docker://`` (``docker run`` rejects that transport).
- (b) ``apptainer.enabled`` with a seeded ``apptainer.cacheDir``: no ``apptainer pull`` call, and
  ``.command.run`` runs ``apptainer exec … <cacheDir>/<simpleName>.img``.

Nextflow wraps ``apptainer exec`` in ``env - PATH="$PATH"``, so the fake apptainer has its log path
written into its body rather than read from the environment. The Nextflow version the run used is
recorded as the ``nextflow_version`` user property and printed.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
CONTAINERS = REPO / "workflow" / "conf" / "containers.config"
MINIPROT = (
    "quay.io/biocontainers/miniprot@sha256:"
    "2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f"
)
SIMPLE_NAME = MINIPROT.replace("/", "-").replace(":", "-") + ".img"

MAIN_NF = """\
process MINIPROT {
  output: path 'out.txt'
  script: "exit 1"
  stub: "touch out.txt"
}
workflow { MINIPROT() }
"""

FAKE_DOCKER = """\
#!/bin/bash
printf '%s\\0' "$@" >> {log}
printf '\\n' >> {log}
[ "$1" = run ] || exit 0
shift
wd=.
while [ "$#" -gt 0 ] && [ "$1" != /bin/bash ]; do
  [ "$1" = -w ] && wd="$2"
  shift
done
[ "$#" -gt 0 ] || exit 0
cd "$wd" && exec "$@"
"""

FAKE_APPTAINER = """\
#!/bin/bash
printf '%s\\0' "$@" >> {log}
printf '\\n' >> {log}
[ "$1" = exec ] || exit 0
shift
while [ "$#" -gt 0 ] && [[ "$1" != *.img ]]; do shift; done
shift
exec "$@"
"""


def _fake(bin_dir: Path, name: str, body: str, log: Path) -> None:
    bin_dir.mkdir(parents=True, exist_ok=True)
    path = bin_dir / name
    path.write_text(body.replace("{log}", shlex.quote(str(log))))
    path.chmod(0o755)
    log.write_text("")


def _calls(log: Path) -> list[list[str]]:
    """One argv per fake call (NUL-separated words, one call per line)."""
    out = []
    for chunk in log.read_text().split("\0\n"):
        if chunk:
            out.append(chunk.split("\0"))
    return out


def _version(env: dict[str, str]) -> str:
    res = subprocess.run(
        ["nextflow", "-version"], capture_output=True, text=True, env=env, timeout=120
    )
    m = re.search(r"version\s+(\S+)", res.stdout)
    return m.group(1) if m else "unknown"


def _run(tmp: Path, engine_cfg: str, record_property) -> tuple[subprocess.CompletedProcess, Path]:  # type: ignore[no-untyped-def]
    pipe = tmp / "pipe"
    pipe.mkdir()
    (pipe / "main.nf").write_text(MAIN_NF)
    (pipe / "nextflow.config").write_text(
        f"includeConfig '{CONTAINERS}'\nprocess.executor = 'local'\n{engine_cfg}\n"
    )
    env = {k: v for k, v in os.environ.items() if not k.startswith("NXF_APPTAINER")}
    env.update(
        PATH=f"{tmp / 'bin'}{os.pathsep}{os.environ['PATH']}",
        NXF_ANSI_LOG="false",
        NXF_DISABLE_CHECK_LATEST="true",
        NXF_OFFLINE="true",
    )
    version = _version(env)
    record_property("nextflow_version", version)
    print(
        f"nextflow version: {version}" + ("" if version.startswith("26.10.") else " (not 26.10.x)")
    )
    work = tmp / "work"
    res = subprocess.run(
        ["nextflow", "run", "main.nf", "-stub-run", "-w", str(work)],
        cwd=pipe, capture_output=True, text=True, env=env, timeout=300, check=False,
    )  # fmt: skip
    return res, work


def _task_dir(work: Path) -> Path:
    dirs = [p.parent for p in work.glob("*/*/.command.run")]
    assert len(dirs) == 1, dirs
    return dirs[0]


def test_example2_docker_run_gets_the_bare_reference(tmp_path: Path, record_property) -> None:  # type: ignore[no-untyped-def]
    log = tmp_path / "docker.log"
    _fake(tmp_path / "bin", "docker", FAKE_DOCKER, log)
    res, work = _run(tmp_path, "docker.enabled = true", record_property)
    runs = [c for c in _calls(log) if c[:1] == ["run"]]
    assert len(runs) == 1, f"docker calls: {_calls(log)}\n{res.stdout}\n{res.stderr}"
    argv = runs[0]
    assert [w for w in argv if w.startswith("docker://")] == [], argv
    assert MINIPROT in argv, argv
    assert res.returncode == 0, res.stdout + res.stderr
    assert (_task_dir(work) / "out.txt").is_file()


def test_example3_apptainer_uses_the_seeded_cache_without_pulling(
    tmp_path: Path,
    record_property,  # type: ignore[no-untyped-def]
) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    seeded = cache / SIMPLE_NAME
    seeded.write_bytes(b"seeded image\n")
    log = tmp_path / "apptainer.log"
    _fake(tmp_path / "bin", "apptainer", FAKE_APPTAINER, log)
    cfg = f"apptainer.enabled = true\napptainer.cacheDir = '{cache}'"
    res, work = _run(tmp_path, cfg, record_property)
    assert res.returncode == 0, res.stdout + res.stderr
    calls = _calls(log)
    assert not any(c[:1] == ["pull"] for c in calls), calls
    execs = [c for c in calls if c[:1] == ["exec"]]
    assert len(execs) == 1 and str(seeded) in execs[0], calls
    task = _task_dir(work)
    run_text = (task / ".command.run").read_text()
    assert re.search(rf"apptainer exec .*{re.escape(str(seeded))}", run_text), run_text
    assert (task / "out.txt").is_file()
    assert seeded.read_bytes() == b"seeded image\n", "the seeded image was replaced"
