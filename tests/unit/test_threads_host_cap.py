"""``gippa run`` caps the default threads to the host on local executors (#514).

``params.threads`` defaults to 16 and ``workflow/conf/base.config`` ties ``resourceLimits.cpus`` to
it, so under a local executor (plain ``apptainer``, ``docker``, ``standard``) a ``process_high``
task asks for 16 CPUs and Nextflow aborts: "Task requirement exceeds available CPUs -- req: 16;
avail: 9". Only ``farmshare_local`` capped by ``SLURM_CPUS_PER_TASK``.

Interface for the implementer (the pure functions decide; ``cli.run`` only wires them):

``gippa.runner.profiles.uses_local_executor(profile: str) -> bool``
    ``profile`` is the resolved, comma-joined value (never ``auto``). ``True`` when Nextflow's
    local executor runs the tasks: no executor profile at all (``apptainer``, ``docker``, ``conda``
    and their combinations), ``standard`` or ``farmshare_local``. ``False`` for every other
    executor profile: ``farmshare``, ``slurm`` and site profiles (M11-15 item 3 treats a site as a
    Slurm profile), whose tasks are separate jobs, so the head job's CPUs do not matter.

``gippa.runner.nextflow.available_cpus(environ: Mapping[str, str]) -> int | None``
    ``SLURM_CPUS_PER_TASK`` when it is a positive integer, else ``len(os.sched_getaffinity(0))``,
    else ``os.cpu_count()``; ``None`` when nothing is known. ``SLURM_CPUS_PER_TASK`` wins even over
    a different affinity: the issue's ``-c 8`` job showed Nextflow 9 CPUs. The two ``os`` functions
    are looked up at call time (the tests patch them; a platform may lack ``sched_getaffinity``).

``gippa.runner.nextflow.capped_threads(threads: int, *, explicit: bool, profile: str,
available: int | None) -> int``
    ``min(threads, available)`` when ``explicit`` is false, ``uses_local_executor(profile)`` and
    ``available`` is known; else ``threads`` unchanged. Pure.

``explicit`` is true when ``--threads``, the config file or ``--set threads=...`` gave the value,
including one equal to the default: comparing with 16 cannot tell ``threads: 16`` in a config file
from no setting at all. How ``cli.run`` learns this is free (e.g. a helper next to
``gippa.config.resolve_config``); the CLI tests pin the behaviour.

``gippa run`` (Nextflow executor, after ``auto`` is resolved) writes the capped value to
``cfg.threads``, so ``OUT/work/params.json`` and ``OUT/resolved_config.yaml`` both carry it, and
prints exactly one stderr line when the value changed, only then::

    gippa run: threads 8 (default 16 exceeds the 8 CPUs available; pass --threads to override)

The tests pin the prefix ``gippa run: threads <value> (``, and that the rest of the line names the
default and ``CPUs``. ``threads`` is in no config hash (tests/unit/test_config_slice.py), so a cap
does not invalidate a ``--resume``.

Real Nextflow propagation: tests/integration/test_fresh_run_nextflow.py.
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
import yaml
from click.testing import CliRunner, Result
from helpers.fake_env import fake_gippa_root, write_setup_manifest
from hypothesis import given
from hypothesis import strategies as st

from gippa.exitcodes import ExitCode

WF = Path(__file__).resolve().parents[2] / "workflow"
CONF = WF / "conf"
FAKE_NF = '#!/bin/sh\nprintf \'%s\\n\' "$@" > "$FAKE_NF_ARGS"\nexit 0\n'
CAP_LINE = re.compile(r"^gippa run: threads ", re.M)
DEFAULT = 16

LOCAL = ["apptainer", "docker", "conda", "docker,conda", "standard", "standard,apptainer",
         "farmshare_local", "farmshare_local,apptainer"]  # fmt: skip
SLURM = ["farmshare", "farmshare,apptainer", "slurm", "slurm,apptainer", "slurm,conda",
         "sherlock,apptainer"]  # fmt: skip


def _affinity(monkeypatch: pytest.MonkeyPatch, cpus: int | None) -> None:
    """``os.sched_getaffinity`` reports ``cpus`` CPUs (``None``: the platform lacks it)."""
    if cpus is None:
        monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    else:
        monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: set(range(cpus)), raising=False)


# ---------------------------------------------------------------- uses_local_executor


@pytest.mark.parametrize("profile", LOCAL)
def test_profiles_with_a_local_executor(profile: str) -> None:
    from gippa.runner.profiles import uses_local_executor

    assert uses_local_executor(profile) is True


@pytest.mark.parametrize("profile", SLURM)
def test_slurm_and_site_profiles_are_not_local(profile: str) -> None:
    from gippa.runner.profiles import uses_local_executor

    assert uses_local_executor(profile) is False


def test_every_bundled_profile_matches_its_config() -> None:
    """A new Slurm profile in nextflow.config cannot be silently treated as local."""
    from gippa.runner.profiles import accepted_profiles, uses_local_executor

    for name in accepted_profiles():
        path = CONF / f"{name}.config"
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        slurm = re.search(r"(?m)^\s*(process\.)?executor\s*=\s*'slurm'", text) is not None
        assert uses_local_executor(name) is not slurm, name


# ---------------------------------------------------------------- available_cpus


def test_slurm_cpus_per_task_wins_over_affinity(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.runner.nextflow import available_cpus

    _affinity(monkeypatch, 9)  # the issue: `-c 8`, and Nextflow saw avail: 9
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    assert available_cpus({"SLURM_CPUS_PER_TASK": "8"}) == 8


@pytest.mark.parametrize("bad", ["", "abc", "0", "-4", "8.5", "8(x2)"])
def test_invalid_slurm_value_falls_back_to_affinity(
    monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    from gippa.runner.nextflow import available_cpus

    _affinity(monkeypatch, 6)
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    assert available_cpus({"SLURM_CPUS_PER_TASK": bad}) == 6


def test_affinity_is_used_before_cpu_count(monkeypatch: pytest.MonkeyPatch) -> None:
    """A cgroup or taskset limit shows in the affinity mask, not in ``os.cpu_count``."""
    from gippa.runner.nextflow import available_cpus

    _affinity(monkeypatch, 4)
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    assert available_cpus({}) == 4


def test_cpu_count_without_affinity(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.runner.nextflow import available_cpus

    _affinity(monkeypatch, None)
    monkeypatch.setattr(os, "cpu_count", lambda: 12)
    assert available_cpus({}) == 12


def test_nothing_known_is_none(monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.runner.nextflow import available_cpus

    _affinity(monkeypatch, None)
    monkeypatch.setattr(os, "cpu_count", lambda: None)
    assert available_cpus({}) is None


@given(cpus=st.integers(min_value=1, max_value=4096))
def test_property_any_positive_slurm_value_is_taken(cpus: int) -> None:
    from gippa.runner.nextflow import available_cpus

    assert available_cpus({"SLURM_CPUS_PER_TASK": str(cpus)}) == cpus


# ---------------------------------------------------------------- capped_threads


@pytest.mark.parametrize(
    ("threads", "explicit", "profile", "available", "want"),
    [
        pytest.param(16, False, "apptainer", 8, 8, id="issue-apptainer-8"),
        pytest.param(16, False, "docker", 1, 1, id="docker-1"),
        pytest.param(16, False, "standard", 15, 15, id="standard-15"),
        pytest.param(16, False, "farmshare_local,apptainer", 8, 8, id="farmshare-local"),
        pytest.param(16, False, "apptainer", 16, 16, id="exactly-enough"),
        pytest.param(16, False, "apptainer", 64, 16, id="big-host-keeps-default"),
        pytest.param(16, False, "apptainer", None, 16, id="unknown-host"),
        pytest.param(16, True, "apptainer", 8, 16, id="explicit-default-value"),
        pytest.param(32, True, "docker", 8, 32, id="explicit-oversubscribes"),
        pytest.param(4, True, "apptainer", 8, 4, id="explicit-small"),
        pytest.param(16, False, "farmshare,apptainer", 2, 16, id="farmshare-head-job"),
        pytest.param(16, False, "slurm,apptainer", 2, 16, id="slurm-head-job"),
    ],
)
def test_capped_threads(
    threads: int, explicit: bool, profile: str, available: int | None, want: int
) -> None:
    from gippa.runner.nextflow import capped_threads

    got = capped_threads(threads, explicit=explicit, profile=profile, available=available)
    assert got == want


_threads = st.integers(min_value=1, max_value=512)
_avail = st.none() | st.integers(min_value=1, max_value=512)


@given(t=_threads, a=_avail, profile=st.sampled_from(LOCAL + SLURM))
def test_property_explicit_is_never_changed(t: int, a: int | None, profile: str) -> None:
    from gippa.runner.nextflow import capped_threads

    assert capped_threads(t, explicit=True, profile=profile, available=a) == t


@given(t=_threads, a=_avail, profile=st.sampled_from(SLURM))
def test_property_slurm_is_never_capped(t: int, a: int | None, profile: str) -> None:
    from gippa.runner.nextflow import capped_threads

    assert capped_threads(t, explicit=False, profile=profile, available=a) == t


@given(t=_threads, a=_avail, profile=st.sampled_from(LOCAL))
def test_property_local_default_fits_the_host(t: int, a: int | None, profile: str) -> None:
    """Never above the host (the abort of #514), never above the default, at least 1, and a second
    application changes nothing."""
    from gippa.runner.nextflow import capped_threads

    got = capped_threads(t, explicit=False, profile=profile, available=a)
    assert got == (t if a is None else min(t, a))
    assert 1 <= got <= t
    assert capped_threads(got, explicit=False, profile=profile, available=a) == got


# ---------------------------------------------------------------- gippa run, end to end


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """Fake ``nextflow``, ``$GIPPA_ROOT`` and ``$GIPPA_DB_DIR``, evidence inputs; an 8-CPU job on
    a 64-CPU node. The preflight tool matrix is not under test, so every profile may launch."""
    from gippa.preflight import run as pf_run

    monkeypatch.setattr(pf_run, "missing_tool_reasons", lambda *a, **k: {})
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    db = tmp_path / "db"
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    _affinity(monkeypatch, 64)
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + seq + "\n")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text(
        "".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n" for i in range(3))
    )
    return {"genome": genome, "proteins": proteins, "args": tmp_path / "nf_args.txt"}


def _run(env: dict[str, Path], out: Path, *extra: str) -> Result:
    from gippa.cli import cli

    argv = ["run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"])]
    res = CliRunner().invoke(cli, [*argv, "--outdir", str(out), *extra])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert env["args"].is_file(), "gippa run did not launch nextflow"
    return res


def _threads_written(out: Path) -> tuple[int, int]:
    """``threads`` of ``OUT/work/params.json`` and ``OUT/resolved_config.yaml``."""
    params = json.loads((out / "work" / "params.json").read_text(encoding="utf-8"))
    resolved = yaml.safe_load((out / "resolved_config.yaml").read_text(encoding="utf-8"))
    return params["threads"], resolved["threads"]


def _assert_capped(res: Result, out: Path, want: int) -> None:
    assert _threads_written(out) == (want, want), res.stderr
    lines = CAP_LINE.findall(res.stderr)
    assert len(lines) == 1, f"want one cap line on stderr, got:\n{res.stderr}"
    line = next(ln for ln in res.stderr.splitlines() if CAP_LINE.match(ln))
    assert line.startswith(f"gippa run: threads {want} ("), line
    assert str(DEFAULT) in line and "CPUs" in line, f"the line must say why: {line}"


def _assert_untouched(res: Result, out: Path, want: int) -> None:
    assert _threads_written(out) == (want, want), res.stderr
    assert not CAP_LINE.search(res.stderr), res.stderr


@pytest.mark.parametrize(
    "profile", ["apptainer", "docker", "standard", "conda", "farmshare_local,apptainer"]
)
def test_run_caps_the_default_under_a_local_executor(
    env: dict[str, Path], tmp_path: Path, profile: str
) -> None:
    """Example of #514: an 8-CPU Slurm job on a bigger node, no --threads."""
    out = tmp_path / "OUT"
    _assert_capped(_run(env, out, "--profile", profile), out, 8)


def test_run_caps_by_affinity_without_slurm(
    env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A workstation or container: no Slurm variable, the mask says 5."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK")
    _affinity(monkeypatch, 5)
    out = tmp_path / "OUT"
    _assert_capped(_run(env, out, "--profile", "docker"), out, 5)


def test_run_invalid_slurm_value_falls_back(
    env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "many")
    _affinity(monkeypatch, 3)
    out = tmp_path / "OUT"
    _assert_capped(_run(env, out, "--profile", "apptainer"), out, 3)


def test_run_big_host_prints_nothing(
    env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    out = tmp_path / "OUT"
    _assert_untouched(_run(env, out, "--profile", "apptainer"), out, DEFAULT)


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "slurm,apptainer"])
def test_run_slurm_executor_is_not_capped(
    env: dict[str, Path], tmp_path: Path, profile: str
) -> None:
    """Tasks are separate jobs: the head job's 8 CPUs do not matter."""
    out = tmp_path / "OUT"
    _assert_untouched(_run(env, out, "--profile", profile), out, DEFAULT)


@pytest.mark.parametrize("given_threads", [DEFAULT, 32, 4])
def test_run_explicit_flag_is_never_changed(
    env: dict[str, Path], tmp_path: Path, given_threads: int
) -> None:
    """Includes ``--threads 16``, the value of the default, on an 8-CPU host."""
    out = tmp_path / "OUT"
    res = _run(env, out, "--profile", "apptainer", "--threads", str(given_threads))
    _assert_untouched(res, out, given_threads)


@pytest.mark.parametrize("value", [DEFAULT, 12])
def test_run_config_file_threads_count_as_explicit(
    env: dict[str, Path], tmp_path: Path, value: int
) -> None:
    conf = tmp_path / "my.yaml"
    conf.write_text(f"threads: {value}\n", encoding="utf-8")
    out = tmp_path / "OUT"
    res = _run(env, out, "--profile", "apptainer", "--config", str(conf))
    _assert_untouched(res, out, value)


def test_run_config_file_without_threads_is_still_the_default(
    env: dict[str, Path], tmp_path: Path
) -> None:
    """A config file that sets other keys does not make the threads explicit."""
    conf = tmp_path / "my.yaml"
    conf.write_text("seed: 7\n", encoding="utf-8")
    out = tmp_path / "OUT"
    _assert_capped(_run(env, out, "--profile", "apptainer", "--config", str(conf)), out, 8)


def test_run_set_threads_counts_as_explicit(env: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    res = _run(env, out, "--profile", "apptainer", "--set", f"threads={DEFAULT}")
    _assert_untouched(res, out, DEFAULT)


def _probe(monkeypatch: pytest.MonkeyPatch, *names: str) -> None:
    from gippa.runner import profiles

    probe = profiles.Probe(on_path=frozenset(names))

    def fake(*_a: Any, **_k: Any) -> Any:
        return probe

    monkeypatch.setattr(profiles, "probe_host", fake)


def test_run_auto_is_decided_on_the_resolved_profile(
    env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``auto`` resolving to plain ``apptainer`` caps; the literal word ``auto`` decides nothing."""
    _probe(monkeypatch, "apptainer")
    out = tmp_path / "OUT"
    _assert_capped(_run(env, out), out, 8)


def test_run_auto_slurm_is_not_capped(
    env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``auto`` resolving to ``slurm,apptainer``, the choice on a cluster like Sherlock."""
    _probe(monkeypatch, "sbatch", "apptainer")
    out = tmp_path / "OUT"
    _assert_untouched(_run(env, out), out, DEFAULT)


# ---------------------------------------------------------------- the cap reaches the tasks


def _code(path: Path) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", path.read_text(encoding="utf-8"))


def test_base_config_ties_resource_limits_to_params_threads() -> None:
    """``resourceLimits`` clamps every request below, so ``params.threads`` is the one knob."""
    text = _code(CONF / "base.config")
    assert re.search(r"resourceLimits\s*=\s*\[\s*cpus\s*:\s*params\.threads\b", text), text


@pytest.mark.parametrize(
    "selector",
    [r"withLabel:\s*'process_high'", r"withName:\s*'TERMINATORNET'",
     r"withName:\s*'DETERMINANTS:DETERMINANTS_STEP'"],
)  # fmt: skip
def test_the_16_cpu_requests_of_the_issue_stay_in_base_config(selector: str) -> None:
    """The requests the limit must clamp; see the real run in the integration test."""
    assert re.search(selector + r"\s*\{\s*cpus\s*=\s*16\b", _code(CONF / "base.config"))


def test_local_profiles_do_not_replace_the_limit() -> None:
    """Under plain ``standard`` (and the runtime profiles) base.config's limit is the one in force;
    ``farmshare_local`` replaces it, but still with ``params.threads`` in the minimum."""
    assert "resourceLimits" not in _code(CONF / "standard.config")
    assert "resourceLimits" not in _code(CONF / "container_common.config")
    profiles = _code(WF / "nextflow.config")
    assert "resourceLimits" not in profiles
    local = _code(CONF / "farmshare_local.config")
    assert re.search(r"cpus:\s*Math\.min\([^\n]*params\.threads as int\)", local), local
