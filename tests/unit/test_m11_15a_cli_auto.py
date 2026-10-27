"""M11-15a: ``--profile auto`` in ``gippa run``/``preflight``/``setup``, the printed choice, the
manifest ``launch`` record, ``setup.runtime`` and the autouse probe fixture.

Spec: docs/specs/M11-15-portable-clusters.md item 1: "Printed" (stderr, before preflight),
"``--profile`` given → used as is ... no detection", R5 (exit 2 with the message), example 1.14 and
edge cases (a)-(d); acceptance criteria 2, 3 and 4.

Interface pinned here (with test_m11_15a_detect_profile.py):

- The CLI calls ``gippa.runner.profiles.probe_host(environ, db_dir)`` (as a module attribute) and
  ``detect_profile``; tests replace ``probe_host`` with a fixed ``Probe``. tests/conftest.py pins it
  for every test to a FarmShare host (``PINNED_PROBE``: ``sbatch`` + ``apptainer``, cluster name
  ``farmshare``), so ``auto`` resolves to today's default (edge case (d); lead decision).
- ``gippa run``: default ``--profile auto``; the detected profile goes to ``-profile``; stderr gets
  ``gippa run: profile <p> (auto: <reason>; pass --profile to choose)`` and each hint on its own
  line. R5: exit 2, ``no container runtime or conda found: install Apptainer (no root needed on
  most Linux systems) or Docker; see <https URL>``; Nextflow is not started.
- ``gippa preflight``: no ``--profile`` = ``auto`` (same detection, line prefixed ``gippa
  preflight:``); ``--profile none`` keeps today's input-only check, without probing.
- ``run_manifest.json`` gets a ``launch`` object: ``profile``, ``profile_source`` (``auto`` |
  ``flag``) and ``profile_reason`` (the reason text under ``auto``); the schema accepts it.
- ``gippa setup --images`` without ``--profile`` uses the detected runtime (FarmShare: apptainer).
  No runtime: the images are skipped with the edge-case (c) warning, the DBs still run, the exit
  code is the DB result, no ``setup.runtime``. Images pulled for a runtime record
  ``{"setup": {"runtime": <r>}}`` at the top of ``db_manifest.json``, schema version ``1.1``
  (schema minor), kept by later manifest writes.
"""

from __future__ import annotations

import json
import os
import random
import re
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit "${FAKE_NF_RC:-0}"
"""
R5 = (
    "no container runtime or conda found: install Apptainer (no root needed on most Linux "
    "systems) or Docker; see https://"
)
HINT_A = "sbatch on PATH but only docker found: tasks run on this host"
HINT_B = "the Nextflow head process runs on this host for the whole run"
IMAGES_SKIPPED = (
    "images: no container runtime found here; rerun gippa setup --images --profile "
    "apptainer|singularity|docker on a host that has one (same --db-dir)"
)
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models", "--no-uniop"]


def _mod() -> Any:
    try:
        from gippa.runner import profiles
    except ImportError as exc:
        pytest.fail(f"gippa.runner.profiles is not implemented (M11-15a): {exc}")
    return profiles


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """A fake ``$GIPPA_ROOT`` and ``$GIPPA_DB_DIR``, a fake ``nextflow`` and evidence inputs."""
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    for var in ("NXF_APPTAINER_CACHEDIR", "GIPPA_CONTAINER_RUNTIME", "FAKE_NF_RC",
                "SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
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
    return {"genome": genome, "proteins": proteins, "args": tmp_path / "nf_args.txt", "db": db}


@pytest.fixture
def set_probe(monkeypatch: pytest.MonkeyPatch) -> Any:
    """``set_probe(*names, **fields)`` makes ``probe_host`` return that ``Probe``; returns the
    list of ``(environ, db_dir)`` calls."""

    def install(*names: str, **kw: Any) -> list[Any]:
        mod = _mod()
        kw.setdefault("docker_ok", "docker" in names)
        probe = mod.Probe(on_path=frozenset(names), **kw)
        calls: list[Any] = []

        def fake(*a: Any, **k: Any) -> Any:
            calls.append((a, k))
            return probe

        monkeypatch.setattr(mod, "probe_host", fake)
        return calls

    return install


@pytest.fixture
def no_probe(monkeypatch: pytest.MonkeyPatch) -> list[Any]:
    """``probe_host`` records calls (none expected) and returns an empty probe. Before M11-15a
    there is nothing to probe, so the flag tests run against today's code as regression guards."""
    calls: list[Any] = []
    try:
        from gippa.runner import profiles as mod
    except ImportError:
        return calls

    def fake(*a: Any, **k: Any) -> Any:
        calls.append((a, k))
        return mod.Probe()

    monkeypatch.setattr(mod, "probe_host", fake)
    return calls


def _invoke(*argv: str) -> tuple[int, str]:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, list(argv))
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _run(env: dict[str, Path], out: Path, *extra: str) -> tuple[int, str]:
    return _invoke("run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"]),
                   "--outdir", str(out), *extra)  # fmt: skip


def _preflight(env: dict[str, Path], *extra: str) -> tuple[int, str]:
    return _invoke("preflight", "--genome", str(env["genome"]), "--proteins",
                   str(env["proteins"]), *extra)  # fmt: skip


def _nf_profile(env: dict[str, Path]) -> str:
    assert env["args"].is_file(), "gippa run did not launch nextflow"
    args = env["args"].read_text().splitlines()
    return args[args.index("-profile") + 1]


def _line(cmd: str, profile: str, reason: str) -> re.Pattern[str]:
    return re.compile(
        rf"gippa {cmd}: profile {re.escape(profile)} \(auto: [^\n]*{re.escape(reason)}[^\n]*; "
        r"pass --profile to choose\)"
    )


# ---------------------------------------------------------------- AC 4: the autouse fixture


def test_fixture_pins_probe_xdg_and_db_dir(tmp_path: Path) -> None:
    from conftest import PINNED_PROBE

    mod = _mod()
    assert getattr(mod.probe_host, "gippa_pinned", False)
    pinned = mod.probe_host(dict(os.environ), None)
    assert pinned == mod.Probe(
        on_path=frozenset({"sbatch", "apptainer"}), docker_ok=False, slurm_job_id=None,
        cluster_name="farmshare", setup_runtime=None,
    )  # fmt: skip
    assert pinned == mod.Probe(**PINNED_PROBE)
    assert mod.detect_profile(pinned).profile == "farmshare,apptainer", "today's default"
    assert Path(os.environ["XDG_CONFIG_HOME"]).is_relative_to(tmp_path)
    assert Path(os.environ["XDG_DATA_HOME"]).is_relative_to(tmp_path)
    assert os.environ["XDG_CONFIG_HOME"] != os.environ["XDG_DATA_HOME"]
    assert "GIPPA_DB_DIR" not in os.environ


def test_bare_preflight_sees_the_pinned_probe(env: dict[str, Path]) -> None:
    """A bare ``preflight`` (no ``--profile``) gets the pinned FarmShare probe, never the real
    host: it resolves to ``farmshare,apptainer`` and passes on the fake root."""
    _mod()
    code, text = _preflight(env)
    assert _line("preflight", "farmshare,apptainer", "FarmShare detected").search(text), text
    assert code == ExitCode.SUCCESS, text


def test_bare_run_sees_the_pinned_probe(env: dict[str, Path], tmp_path: Path) -> None:
    _mod()
    code, text = _run(env, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    assert _nf_profile(env) == "farmshare,apptainer"


@pytest.mark.real_host_probe
def test_real_host_probe_opts_out(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """The marker keeps the module's own ``probe_host`` and leaves the XDG variables alone."""
    mod = _mod()
    assert not getattr(mod.probe_host, "gippa_pinned", False)
    assert mod.probe_host.__module__ == mod.__name__
    for var in ("XDG_CONFIG_HOME", "XDG_DATA_HOME"):
        value = os.environ.get(var)
        assert value is None or not Path(value).is_relative_to(tmp_path)


def test_probe_override_wins_over_the_fixture(set_probe: Any) -> None:
    mod = _mod()
    set_probe("docker")
    assert mod.probe_host({}, None).on_path == frozenset({"docker"})


# ---------------------------------------------------------------- AC 2: gippa run


def test_run_auto_farmshare(env: dict[str, Path], set_probe: Any, tmp_path: Path) -> None:
    calls = set_probe("sbatch", "apptainer", cluster_name="farmshare", slurm_job_id="1")
    code, text = _run(env, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    assert calls, "no --profile: the host is probed"
    assert _nf_profile(env) == "farmshare,apptainer"
    assert _line("run", "farmshare,apptainer", "FarmShare detected").search(text), text


def test_run_auto_passes_the_db_dir_to_the_probe(
    env: dict[str, Path], set_probe: Any, tmp_path: Path
) -> None:
    """R0 reads ``setup.runtime`` of the run's DB directory (``--db-dir``, else
    ``$GIPPA_DB_DIR``)."""
    calls = set_probe("sbatch", "apptainer", cluster_name="farmshare")
    _run(env, tmp_path / "OUT")
    args, kwargs = calls[0]
    seen = [*args, *kwargs.values()]
    assert any(isinstance(v, Path) and v.resolve() == env["db"].resolve() for v in seen), seen


def test_run_explicit_auto_is_the_default(
    env: dict[str, Path], set_probe: Any, tmp_path: Path
) -> None:
    set_probe("sbatch", "apptainer", cluster_name="farmshare")
    code, text = _run(env, tmp_path / "OUT", "--profile", "auto")
    assert code == ExitCode.SUCCESS, text
    assert _nf_profile(env) == "farmshare,apptainer"


def test_run_auto_slurm_login_node(env: dict[str, Path], set_probe: Any, tmp_path: Path) -> None:
    """Example 1.13 through the CLI: ``slurm,apptainer`` on the command line and hint (b)."""
    set_probe("sbatch", "apptainer")
    code, text = _run(env, tmp_path / "OUT")
    assert _line("run", "slurm,apptainer", "sbatch on PATH").search(text), text
    assert HINT_B in text
    assert code == ExitCode.SUCCESS, text
    assert _nf_profile(env) == "slurm,apptainer"


def test_run_printed_before_preflight(env: dict[str, Path], set_probe: Any, tmp_path: Path) -> None:
    """The line and the hints precede the preflight output, even when preflight fails."""
    set_probe("sbatch", "docker", slurm_job_id="1")
    missing = tmp_path / "nodb"
    _code, text = _run(env, tmp_path / "OUT", "--db-dir", str(missing))
    assert re.search(r"gippa run: profile docker \(auto: ", text), text
    assert HINT_A in text  # example 1.12, edge case (a)


def test_run_r5_exits_2(env: dict[str, Path], set_probe: Any, tmp_path: Path) -> None:
    """Example 1.8: nothing on PATH."""
    set_probe()  # nothing found
    out = tmp_path / "OUT"
    code, text = _run(env, out)
    assert code == ExitCode.ENVIRONMENT, text
    assert R5 in text
    assert not env["args"].exists(), "Nextflow is not started"


def test_run_flag_skips_detection(
    env: dict[str, Path], no_probe: list[Any], tmp_path: Path
) -> None:
    code, text = _run(env, tmp_path / "OUT", "--profile", "farmshare,apptainer")
    assert code == ExitCode.SUCCESS, text
    assert no_probe == [], "--profile given: no detection"
    assert "(auto:" not in text
    assert _nf_profile(env) == "farmshare,apptainer"


def test_run_flag_docker_skips_detection(
    env: dict[str, Path], no_probe: list[Any], tmp_path: Path
) -> None:
    _run(env, tmp_path / "OUT", "--profile", "docker")
    assert no_probe == []


def test_run_accepts_slurm_profile(
    env: dict[str, Path], no_probe: list[Any], tmp_path: Path
) -> None:
    code, text = _run(env, tmp_path / "OUT", "--profile", "slurm,apptainer")
    assert code == ExitCode.SUCCESS, text
    assert _nf_profile(env) == "slurm,apptainer"


# ---------------------------------------------------------------- gippa preflight


def test_preflight_none_is_input_only(env: dict[str, Path], no_probe: list[Any]) -> None:
    code, text = _preflight(env, "--profile", "none")
    assert code == ExitCode.SUCCESS, text
    assert no_probe == []
    assert "(auto:" not in text


def test_preflight_auto_prints_and_checks_tools(env: dict[str, Path], set_probe: Any) -> None:
    calls = set_probe("sbatch", "apptainer", cluster_name="farmshare")
    code, text = _preflight(env)
    assert calls
    assert _line("preflight", "farmshare,apptainer", "FarmShare detected").search(text), text
    assert code == ExitCode.SUCCESS, text
    assert "environment.tools" in text, "auto: the tool check runs as under --profile"


def test_preflight_flag_skips_detection(env: dict[str, Path], no_probe: list[Any]) -> None:
    code, text = _preflight(env, "--profile", "farmshare,apptainer")
    assert code == ExitCode.SUCCESS, text
    assert no_probe == []


# ---------------------------------------------------------------- AC 3: manifest


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def test_manifest_launch_auto(env: dict[str, Path], set_probe: Any, tmp_path: Path) -> None:
    from helpers.manifest_data import validate

    set_probe("sbatch", "apptainer", cluster_name="farmshare")
    out = tmp_path / "OUT"
    code, text = _run(env, out)
    assert code == ExitCode.SUCCESS, text
    m = _manifest(out)
    launch = m.get("launch")
    assert isinstance(launch, dict), sorted(m)
    assert launch["profile"] == "farmshare,apptainer"
    assert launch["profile_source"] == "auto"
    assert "FarmShare detected" in launch["profile_reason"]
    assert m["profile"] == "farmshare,apptainer"
    validate(m)


def test_manifest_launch_flag(env: dict[str, Path], no_probe: list[Any], tmp_path: Path) -> None:
    from helpers.manifest_data import validate

    out = tmp_path / "OUT"
    code, text = _run(env, out, "--profile", "farmshare,apptainer")
    assert code == ExitCode.SUCCESS, text
    m = _manifest(out)
    assert m["launch"]["profile"] == "farmshare,apptainer"
    assert m["launch"]["profile_source"] == "flag"
    validate(m)


def test_launch_is_volatile() -> None:
    """``launch`` differs between equivalent runs (host, flag vs auto): D-C volatile."""
    from gippa.manifest import VOLATILE

    assert any(v == "/launch" or v.startswith("/launch/") for v in VOLATILE)


# ---------------------------------------------------------------- AC 3: setup.runtime


@pytest.fixture
def setup_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """``run_setup`` and the image pulls replaced by recorders; returns their call logs."""
    import gippa.setupdb
    from gippa.setupdb import images

    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "root"))
    for var in ("NXF_APPTAINER_CACHEDIR", "GIPPA_CONTAINER_RUNTIME"):
        monkeypatch.delenv(var, raising=False)
    log: dict[str, Any] = {"pulls": [], "dbs": [], "db_code": ExitCode.SUCCESS}

    def pull(runtime: str, tool: str, ref: str, version: str, cache: Path) -> Any:
        log["pulls"].append(runtime)
        return images.image_record(tool, ref, version, "ok", f"pulled {ref}")

    def run_setup(*a: Any, **k: Any) -> ExitCode:
        log["dbs"].append(a)
        return ExitCode(log["db_code"])

    monkeypatch.setattr(images, "_pull", pull)
    monkeypatch.setattr(gippa.setupdb, "run_setup", run_setup)
    log["db"] = tmp_path / "db"
    return log


def _setup(log: dict[str, Any], *extra: str) -> tuple[int, str]:
    return _invoke("setup", "--db-dir", str(log["db"]), *SKIP_DBS, *extra)


def _db_doc(db: Path) -> dict[str, Any]:
    path = db / "db_manifest.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def _validate_db(doc: dict[str, Any]) -> None:
    import jsonschema

    from gippa.schemas import __file__ as schemas_init

    schema = json.loads((Path(schemas_init).parent / "db_manifest.schema.json").read_text())
    jsonschema.Draft202012Validator(schema).validate(doc)


def test_1_14_setup_without_runtime(setup_env: dict[str, Any], set_probe: Any) -> None:
    """Edge case (c): images skipped with the warning, DBs run, exit by the DB result."""
    set_probe()  # nothing found
    code, text = _setup(setup_env, "--images")
    assert code == ExitCode.SUCCESS, text
    assert IMAGES_SKIPPED in text
    assert setup_env["dbs"], "the databases are still set up"
    assert setup_env["pulls"] == []
    assert "runtime" not in _db_doc(setup_env["db"]).get("setup", {})


def test_1_14_exit_is_the_db_result(setup_env: dict[str, Any], set_probe: Any) -> None:
    set_probe()  # nothing found
    setup_env["db_code"] = ExitCode.INVALID_INPUT
    code, text = _setup(setup_env, "--images")
    assert code == ExitCode.INVALID_INPUT, text
    assert IMAGES_SKIPPED in text


def test_setup_images_detected_runtime(setup_env: dict[str, Any], set_probe: Any) -> None:
    set_probe("apptainer")
    code, text = _setup(setup_env, "--images")
    assert code == ExitCode.SUCCESS, text
    assert setup_env["pulls"] and set(setup_env["pulls"]) == {"apptainer"}
    doc = _db_doc(setup_env["db"])
    assert doc["setup"]["runtime"] == "apptainer"
    assert doc["schema_version"] == "1.1"
    _validate_db(doc)


def test_setup_images_farmshare_is_apptainer(setup_env: dict[str, Any], set_probe: Any) -> None:
    set_probe("sbatch", "apptainer", "docker", cluster_name="farmshare")
    code, text = _setup(setup_env, "--images")
    assert code == ExitCode.SUCCESS, text
    assert set(setup_env["pulls"]) == {"apptainer"}
    assert _db_doc(setup_env["db"])["setup"]["runtime"] == "apptainer"


def test_setup_images_flag_skips_detection(setup_env: dict[str, Any], no_probe: list[Any]) -> None:
    code, text = _setup(setup_env, "--images", "--profile", "docker")
    assert code == ExitCode.SUCCESS, text
    assert no_probe == []
    assert set(setup_env["pulls"]) == {"docker"}
    doc = _db_doc(setup_env["db"])
    assert doc["setup"] == {"runtime": "docker"}, "only the runtime, never an executor"
    _validate_db(doc)


def test_setup_runtime_survives_later_writes(setup_env: dict[str, Any], set_probe: Any) -> None:
    from gippa.setupdb.items import DbRecord
    from gippa.setupdb.manifest import write_manifest

    set_probe("apptainer")
    _setup(setup_env, "--images")
    rec = DbRecord("genomad", "1.9", "https://x", "", 0, "", "", "absent", "")
    write_manifest(setup_env["db"], [rec])
    doc = _db_doc(setup_env["db"])
    assert doc["setup"]["runtime"] == "apptainer"
    _validate_db(doc)


def test_setup_without_images_records_no_runtime(setup_env: dict[str, Any], set_probe: Any) -> None:
    set_probe("apptainer")
    code, text = _setup(setup_env)
    assert code == ExitCode.SUCCESS, text
    assert setup_env["pulls"] == []
    assert "runtime" not in _db_doc(setup_env["db"]).get("setup", {})


def test_db_schema_accepts_setup_runtime() -> None:
    for runtime in ("apptainer", "singularity", "docker"):
        _validate_db({"schema_version": "1.1", "items": [], "setup": {"runtime": runtime}})
    import jsonschema

    with pytest.raises(jsonschema.ValidationError):
        _validate_db({"schema_version": "1.1", "items": [], "setup": {"runtime": "slurm"}})
