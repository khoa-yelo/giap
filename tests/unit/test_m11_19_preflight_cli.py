"""M11-19 acceptance criterion 6: ``gippa preflight``/``gippa run --profile docker`` stop when a
pinned image the run needs is not local (runs never pull).

docs/specs/M11-19-nextflow-docker-pull-never.md item 2 "Surfaces": the reasons reach the
``environment.tools`` check (error, exit 2, ``value.missing``) through ``_tools_check``; ``gippa
run --profile docker`` exits 2 before Nextflow starts.

CliRunner with a fake ``$GIPPA_ROOT``/``$GIPPA_DB_DIR`` (tests/helpers/fake_env.py), a fake
``nextflow`` that records its argv (as in test_m11_15a_cli_auto.py) and a logging fake ``docker``
first on ``PATH`` whose ``image inspect`` succeeds only for the refs in a tags file. No network,
no real docker.

Under plain ``docker`` the conda-only tools are missing anyway (they would run on the host), so
those tests assert the Bakta reason itself; under ``docker,conda`` the fake root's prebuilt envs
satisfy the conda-only tools and the pinned images alone decide the exit code.
"""

from __future__ import annotations

import json
import os
import random
import stat
import sys
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, pinned_images, write_setup_manifest

from gippa.exitcodes import ExitCode

SETUP = "gippa setup --images --profile docker"
BAKTA = pinned_images()["bakta"]
ALL_REFS = sorted(set(pinned_images().values()))
ABSENT_BAKTA = (
    f"bakta: image {BAKTA} is not present locally and runs never pull it (docker image inspect "
    f"{BAKTA} failed); run: {SETUP}"
)

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""

FAKE_DOCKER = r"""#!{python} -I
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_DOCKER_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
if argv == ["--version"]:
    print("Docker version 27.3.1, build ce12230")
    sys.exit(0)
if argv[:2] == ["image", "inspect"]:
    refs = [a for a in argv[2:] if not a.startswith("-") and not a.startswith("{{")]
    with open(os.environ["FAKE_DOCKER_TAGS"]) as fh:
        tags = {{l.strip() for l in fh if l.strip()}}
    if refs and all(r in tags for r in refs):
        print("[]")
        sys.exit(0)
    sys.stderr.write("Error: No such image: " + refs[-1] + "\n")
    sys.exit(1)
sys.stderr.write("fake docker: unexpected call\n")
sys.exit(3)
"""


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    for var in ("NXF_APPTAINER_CACHEDIR", "GIPPA_CONTAINER_RUNTIME", "GIPPA_FGSEA_PREFIX",
                "SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    _exe(bindir / "nextflow", FAKE_NF)
    _exe(bindir / "docker", FAKE_DOCKER.format(python=sys.executable))
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    monkeypatch.setenv("FAKE_DOCKER_LOG", str(tmp_path / "docker.log"))
    monkeypatch.setenv("FAKE_DOCKER_TAGS", str(tmp_path / "tags.txt"))
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(seq[i:i + 60]
                      for i in range(0, 6000, 60)) + "\n")  # fmt: skip
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text("".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n"
                                for i in range(3)))  # fmt: skip
    return {"genome": genome, "proteins": proteins, "args": tmp_path / "nf_args.txt",
            "tags": tmp_path / "tags.txt", "log": tmp_path / "docker.log"}  # fmt: skip


def _tags(env: dict[str, Path], refs: list[str]) -> None:
    env["tags"].write_text("".join(r + "\n" for r in refs))


def _invoke(*argv: str) -> tuple[int, str]:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, list(argv))
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _preflight(env: dict[str, Path], out: Path, profile: str) -> tuple[int, str]:
    return _invoke("preflight", "--genome", str(env["genome"]), "--proteins",
                   str(env["proteins"]), "--outdir", str(out), "--profile", profile)  # fmt: skip


def _run(env: dict[str, Path], out: Path, profile: str) -> tuple[int, str]:
    return _invoke("run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"]),
                   "--outdir", str(out), "--profile", profile)  # fmt: skip


def _tools_check(out: Path) -> dict[str, Any]:
    doc = json.loads((out / "preflight" / "preflight.json").read_text())
    checks = [c for c in doc["checks"] if c["rule_id"] == "environment.tools"]
    assert len(checks) == 1, checks
    return dict(checks[0])


def _no_pulls(env: dict[str, Path]) -> None:
    log = env["log"]
    for line in log.read_text().splitlines() if log.exists() else []:
        argv = json.loads(line)
        assert argv == ["--version"] or argv[:2] == ["image", "inspect"], argv


# ---------------------------------------------------------------- gippa preflight


def test_preflight_docker_absent_bakta_exits_2_and_names_it(
    env: dict[str, Path], tmp_path: Path
) -> None:
    _tags(env, [r for r in ALL_REFS if r != BAKTA])
    out = tmp_path / "OUT"
    code, text = _preflight(env, out, "docker")
    assert code == ExitCode.ENVIRONMENT == 2, text
    check = _tools_check(out)
    assert check["level"] == "error" and "bakta" in check["value"]["missing"], check
    assert ABSENT_BAKTA in text, text
    assert ABSENT_BAKTA in (out / "preflight" / "preflight.md").read_text()
    _no_pulls(env)


def test_preflight_docker_all_refs_reports_no_pinned_tool(
    env: dict[str, Path], tmp_path: Path
) -> None:
    _tags(env, ALL_REFS)
    out = tmp_path / "OUT"
    _preflight(env, out, "docker")
    check = _tools_check(out)
    missing = set((check.get("value") or {}).get("missing") or [])
    assert not missing & set(pinned_images()), check
    log = env["log"].read_text() if env["log"].exists() else ""
    assert BAKTA in log, "the bakta ref was never checked with docker image inspect"
    _no_pulls(env)


def test_preflight_docker_conda_decided_by_images(env: dict[str, Path], tmp_path: Path) -> None:
    _tags(env, [r for r in ALL_REFS if r != BAKTA])
    code, text = _preflight(env, tmp_path / "OUT1", "docker,conda")
    assert code == 2, text
    assert set(_tools_check(tmp_path / "OUT1")["value"]["missing"]) == {"bakta"}, text
    _tags(env, ALL_REFS)
    code, text = _preflight(env, tmp_path / "OUT2", "docker,conda")
    assert _tools_check(tmp_path / "OUT2")["level"] != "error", text


# ---------------------------------------------------------------- gippa run


def test_run_docker_absent_bakta_exits_2_without_nextflow(
    env: dict[str, Path], tmp_path: Path
) -> None:
    _tags(env, [r for r in ALL_REFS if r != BAKTA])
    out = tmp_path / "OUT"
    code, text = _run(env, out, "docker")
    assert code == 2, text
    assert not env["args"].exists(), "Nextflow was started although preflight failed"
    assert ABSENT_BAKTA in text, text
    assert "bakta" in _tools_check(out)["value"]["missing"]
    _no_pulls(env)
