"""M11-14d review round (PR #387): console visibility of the UniOP fallback and the test gaps.

S1. ADR-0061 decision 5 and docs/specs/M11-14-container-install.md item 1 (c) "The fallback is
never silent" (Khoa 2026-10-26: "when uniop failed, need to flag the user in preflight"). When
preflight passes with a UniOP fallback warning, ``giap run`` prints that warning to stderr
before Nextflow starts, not only into ``OUT/preflight/preflight.{json,md}``. The cases:
``enabled: false``, ``auto`` resolving to false, a recorded failed setup and a failed probe.

S2. Guards for mutants that survived review:
- ``probe_uniop`` internals (micromamba argv, timeout 124, nothing startable 127), run against a
  fake env ``python`` (the global conftest stub is opted out with ``real_uniop_probe``);
- the ``test``/``stub`` no-probe gate of ``run_preflight``;
- re-resolving a recorded ``auto`` -> false;
- the local executor: ``uniop`` skipped under ``enabled: false``; ``giap run --executor local``
  resolves ``auto``;
- Nextflow wiring: OPERONS and CONTEXT_STATS receive ``operon_meta`` (static checks).

Fakes only: ``helpers.fake_env.fake_giap_root``, a fake ``nextflow`` on ``PATH``, a fake env
``python`` script. No network; never the real ``$GIAP_DB_DIR`` or ``$GIAP_ROOT``.
"""

from __future__ import annotations

import json
import os
import random
import re
import stat
import sys
import tomllib
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, write_setup_manifest

from giap.config import GiapConfig, resolve_config
from giap.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
PROFILE = "farmshare,apptainer"
QUIET = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
AUTO_WARNING = (
    "uniop not installed (unlicensed, not distributed); operons use the distance model; "
    f"run: {PKG} setup --uniop"
)
FALSE_WARNING = "annotation.uniop.enabled is false; operons use the distance model (see Methods)"
SETUP_MSG = "apptainer build failed (fakeroot unavailable or build error): FATAL: no fakeroot"
PROBE_TAIL = "ModuleNotFoundError: No module named 'sklearn'"
PROBE_STDERR = f'  File "UniOP", line 3, in <module>\n{PROBE_TAIL}\n'
#: written to the CliRunner's stderr when ``run_nextflow`` is entered
NF_MARKER = "<<FAKE-NEXTFLOW-LAUNCH>>"


def _set(enabled: bool | str) -> str:
    return f"annotation.uniop.enabled={str(enabled).lower()}"


def _cfg(enabled: bool | str, *extra: str) -> GiapConfig:
    return resolve_config(overrides=[*QUIET, _set(enabled), *extra])


def _env_mod() -> Any:
    from giap.preflight import environment

    return environment


def _probe(monkeypatch: pytest.MonkeyPatch, code: int, err: str) -> list[Any]:
    calls: list[Any] = []

    def fake(*args: Any, **kw: Any) -> tuple[int, str]:
        calls.append((args, kw))
        return code, err

    monkeypatch.setattr(_env_mod(), "probe_uniop", fake, raising=False)
    return calls


def _root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, uniop: bool) -> Path:
    root = fake_giap_root(tmp_path / "giap_root", skip=() if uniop else ("uniop",))
    monkeypatch.setenv("GIAP_ROOT", str(root))
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)
    return root


def _record_uniop(db: Path, status: str, message: str) -> Path:
    path = db / "db_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["items"].append({"name": "image:uniop", "version": "git:3ed14dfc", "checksum": None,
                         "status": status, "message": message})  # fmt: skip
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    return db


# ---------------------------------------------------------------- S1: giap run, stderr

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    rng = random.Random(3)
    path = tmp_path / "in" / "genome.fa"
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    path.write_text(f">chr1 [topology=circular]\n{seq}\n")
    return path


@pytest.fixture
def proteins(tmp_path: Path) -> Path:
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    path = tmp_path / "in" / "proteins.faa"
    path.parent.mkdir(parents=True, exist_ok=True)
    recs = [f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n" for i in range(3)]
    path.write_text("".join(recs))
    return path


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A fake ``nextflow`` on PATH; entering ``run_nextflow`` writes :data:`NF_MARKER` to the
    (CliRunner-captured) stderr first, so the tests can order the warning before the launch."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    from giap.runner import nextflow as nf

    real = nf.run_nextflow

    def marked(*a: Any, **kw: Any) -> Any:
        sys.stderr.write(f"{NF_MARKER}\n")
        sys.stderr.flush()
        return real(*a, **kw)

    monkeypatch.setattr(nf, "run_nextflow", marked)
    return args


def _giap_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path,
    *, uniop: bool, sets: tuple[str, ...] = (), record: str | None = None,
) -> tuple[int, str]:  # fmt: skip
    """``(exit code, stderr)`` of ``giap run`` under :data:`PROFILE`."""
    from giap.cli import cli

    _root(tmp_path, monkeypatch, uniop=uniop)
    db = write_setup_manifest(tmp_path / "db_env")
    if record is not None:
        _record_uniop(db, "failed", record)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    out = tmp_path / "OUT"
    argv = ["run", "--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out),
            "--profile", PROFILE]  # fmt: skip
    for s in (*QUIET, *sets):
        argv += ["--set", s]
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.stderr


def _before_launch(stderr: str, *wanted: str) -> None:
    assert NF_MARKER in stderr, f"nextflow was not launched:\n{stderr}"
    head = stderr.split(NF_MARKER, 1)[0]
    for text in wanted:
        assert text in head, (
            f"giap run did not print {text!r} to stderr before Nextflow started "
            f"(ADR-0061 decision 5); stderr:\n{stderr}"
        )


@pytest.mark.usefixtures("fake_nf")
def test_s1_false_warning_on_stderr_before_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, err = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=True,
                          sets=(_set(False),))  # fmt: skip
    assert code == ExitCode.SUCCESS, err
    _before_launch(err, FALSE_WARNING)


@pytest.mark.usefixtures("fake_nf")
def test_s1_auto_not_installed_warning_on_stderr_before_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, err = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False,
                          sets=(_set("auto"),))  # fmt: skip
    assert code == ExitCode.SUCCESS, err
    _before_launch(err, AUTO_WARNING)


@pytest.mark.usefixtures("fake_nf")
def test_s1_recorded_setup_failure_on_stderr_before_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, err = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False,
                          sets=(_set("auto"),), record=SETUP_MSG)  # fmt: skip
    assert code == ExitCode.SUCCESS, err
    _before_launch(err, AUTO_WARNING, f"(last setup: {SETUP_MSG})")


@pytest.mark.usefixtures("fake_nf")
def test_s1_failed_probe_on_stderr_before_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    _probe(monkeypatch, 1, PROBE_STDERR)
    code, err = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=True,
                          sets=(_set("auto"),))  # fmt: skip
    assert code == ExitCode.SUCCESS, err
    _before_launch(err, AUTO_WARNING, f"UniOP --help exited 1: {PROBE_TAIL}")


@pytest.mark.usefixtures("fake_nf")
def test_s1_no_uniop_warning_when_uniop_runs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    """Guard: the default ``true`` with UniOP installed prints no fallback warning."""
    code, err = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=True)
    assert code == ExitCode.SUCCESS, err
    assert NF_MARKER in err, err
    assert "distance model" not in err.split(NF_MARKER, 1)[0], err


# ---------------------------------------------------------------- S2.1: probe internals

FAKE_PY = """#!/bin/sh
printf '%s\\n' "$0" "$@" > "$FAKE_PY_ARGS"
{body}
"""


def _fake_env_python(env_root: Path, body: str, monkeypatch: pytest.MonkeyPatch) -> Path:
    """``<env_root>/<uniop env>/bin/python``: records its argv to ``$FAKE_PY_ARGS``."""
    from giap.annotate.tools import TOOLS

    exe = env_root / TOOLS["uniop"].env_name / "bin" / "python"
    exe.parent.mkdir(parents=True, exist_ok=True)
    exe.write_text(FAKE_PY.format(body=body))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    rec = env_root.parent / "py_args.txt"
    monkeypatch.setenv("FAKE_PY_ARGS", str(rec))
    monkeypatch.delenv("UNIOP_HOME", raising=False)
    return rec


def _real_probe(profile: str, tmp_path: Path, env_root: Path | None) -> tuple[int, str]:
    fn = _env_mod().probe_uniop
    code, err = fn(profile, image_cache=tmp_path / "cache", env_root=env_root, local_tools=None)
    return int(code), str(err)


@pytest.mark.real_uniop_probe
@pytest.mark.parametrize("profile", ["conda", "farmshare", PROFILE])
def test_probe_micromamba_runs_the_env_python_with_help(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, profile: str
) -> None:
    from giap.annotate.tools import TOOLS

    env_root = tmp_path / "envs"
    rec = _fake_env_python(env_root, "echo probe-ok >&2\nexit 0", monkeypatch)
    code, _ = _real_probe(profile, tmp_path, env_root)
    assert code == 0
    prefix = env_root / TOOLS["uniop"].env_name
    argv = rec.read_text(encoding="utf-8").splitlines()
    assert argv == [str(prefix / "bin" / "python"), f"{prefix}/share/uniop/src/UniOP", "--help"]


@pytest.mark.real_uniop_probe
def test_probe_returns_the_exit_code_and_stderr(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    env_root = tmp_path / "envs"
    _fake_env_python(env_root, f'echo "{PROBE_TAIL}" >&2\nexit 3', monkeypatch)
    code, err = _real_probe("conda", tmp_path, env_root)
    assert code == 3
    assert PROBE_TAIL in err


@pytest.mark.real_uniop_probe
def test_probe_timeout_gives_124(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    env_root = tmp_path / "envs"
    _fake_env_python(env_root, "exec sleep 10", monkeypatch)
    monkeypatch.setattr(_env_mod(), "PROBE_TIMEOUT_S", 0.5)
    code, err = _real_probe("conda", tmp_path, env_root)
    assert code == 124, (code, err)
    assert "timed out" in err


@pytest.mark.real_uniop_probe
@pytest.mark.parametrize(
    ("profile", "env_root"),
    [
        ("local", "envs"),  # no container or conda profile
        ("conda", None),  # no $GIAP_ROOT
        ("conda", "envs"),  # the env python does not exist (OSError)
    ],
)
def test_probe_nothing_startable_gives_127(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, profile: str, env_root: str | None
) -> None:
    monkeypatch.delenv("UNIOP_HOME", raising=False)
    root = None if env_root is None else tmp_path / env_root
    code, err = _real_probe(profile, tmp_path, root)
    assert code == 127, (code, err)
    assert err.strip(), "127 without a reason"


# ---------------------------------------------------------------- S2.2: test/stub no-probe gate

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture
def inputs(tmp_path: Path) -> Any:
    from giap.preflight.run import PreflightInputs

    rng = random.Random(0)
    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(parents=True, exist_ok=True)
    g.write_text(
        f">chr1 [topology=circular]\n{''.join(rng.choice('ACGT') for _ in range(120_000))}\n"
    )
    p = tmp_path / "in" / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return PreflightInputs(genome=g, proteins=p)


@pytest.mark.parametrize("profile", ["test", "stub", "test,docker"])
def test_no_probe_under_test_or_stub_profiles(
    inputs: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, profile: str
) -> None:
    from giap.preflight.run import exit_code, run_preflight

    _root(tmp_path, monkeypatch, uniop=True)
    calls = _probe(monkeypatch, 127, "no container or conda profile")
    db = write_setup_manifest(tmp_path / "db")
    report, _ = run_preflight(inputs, _cfg(True), None, "run", tmp_path / "OUT", db,
                              profile=profile, publication=False)  # fmt: skip
    assert calls == [], calls
    tools = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert all(c.level != "error" for c in tools), tools
    assert int(exit_code(report, False)) == ExitCode.SUCCESS


# ---------------------------------------------------------------- S2.3: recorded auto -> false


def test_re_resolving_a_recorded_auto_false_keeps_false_when_installed(tmp_path: Path) -> None:
    """A rerun or ``-resume`` re-runs preflight on ``resolved_config.yaml``: the recorded
    ``auto`` -> false keeps its outcome even if UniOP is installed now."""
    from giap.runner.nextflow import write_params

    resolve = _env_mod().resolve_uniop
    once = resolve(_cfg("auto"), installed=False, publication=False)
    assert once.run_uniop is False
    _, resolved = write_params(tmp_path / "OUT", once.cfg, {})
    back = resolve_config(resolved)
    for cfg in (once.cfg, back):
        again = resolve(cfg, installed=True, publication=False)
        assert again.run_uniop is False, again
        assert again.cfg.annotation.uniop.enabled is False
        assert again.cfg.model_dump(mode="json") == once.cfg.model_dump(mode="json")


# ---------------------------------------------------------------- S2.4: local executor


def test_local_plan_skips_uniop_when_disabled() -> None:
    from giap.runner import local
    from giap.steps import load_step_modules

    load_step_modules()
    plan = local.plan_local(_cfg(False), ["uniop", "operon_pairs"])
    assert "uniop" not in plan
    assert "operon_pairs" in plan
    with pytest.raises(local.LayerMissing):  # enabled: UniOP is external, local cannot run it
        local.plan_local(_cfg(True), ["uniop", "operon_pairs"])


def test_giap_run_executor_local_resolves_auto(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    from giap.cli import cli
    from giap.runner import local

    seen: list[GiapConfig] = []

    def plan(cfg: GiapConfig, *a: Any, **kw: Any) -> list[str]:
        seen.append(cfg)
        raise local.LayerMissing("stop here (test)")

    monkeypatch.setattr(local, "plan_local", plan)
    monkeypatch.delenv("GIAP_DB_DIR", raising=False)
    argv = ["run", "--executor", "local", "--genome", str(genome), "--proteins", str(proteins),
            "--outdir", str(tmp_path / "OUT"), "--set", _set("auto")]  # fmt: skip
    for s in QUIET:
        argv += ["--set", s]
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert len(seen) == 1, res.output
    uni = seen[0].annotation.uniop
    assert isinstance(uni.enabled, bool), f"--executor local planned with enabled={uni.enabled!r}"
    assert uni.enabled is True


# ---------------------------------------------------------------- S2.5: Nextflow wiring (static)

MAIN_NF = REPO / "workflow" / "main.nf"
CONTEXT_NF = REPO / "workflow" / "subworkflows" / "local" / "context.nf"


def _code_only(text: str) -> str:
    return "\n".join(re.sub(r"//.*$", "", ln) for ln in text.splitlines())


def _call(text: str, name: str) -> str:
    """The argument text of the first ``name(`` call (balanced parentheses)."""
    m = re.search(rf"(?<![\w.]){re.escape(name)}\s*\(", text)
    assert m is not None, f"no {name}( call"
    depth, i = 1, m.end()
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[m.end() : i - 1]


OPERON_META_KV = re.compile(r"kv\(\s*'operon_meta'\s*,\s*ANNOTATE\.out\.operon_meta\s*\)")


@pytest.mark.parametrize("process", ["CONTEXT", "OPERONS"])
def test_main_nf_passes_operon_meta(process: str) -> None:
    args = _call(_code_only(MAIN_NF.read_text(encoding="utf-8")), process)
    assert OPERON_META_KV.search(args), f"{process} does not receive operon_meta:\n{args}"


def test_context_nf_forwards_operon_meta_to_context_stats() -> None:
    text = _code_only(CONTEXT_NF.read_text(encoding="utf-8"))
    m = re.search(r"def optional = \[(.*?)\]\n", text, re.S)
    assert m is not None, "no `def optional = [...]` in context.nf"
    optional = m.group(1)
    assert re.search(r"layers_for_pairs\.filter\s*\{[^}]*'operon_meta'", optional), optional
    cs = _call(text, "stepIn")
    assert cs.startswith("'context_stats'") and "optional" in cs, cs
