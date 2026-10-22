"""T3: ``giap.io.fgsea_r`` launcher and subprocess behaviour, no R (M8-17a §3, §4; A4).

A fake R front end sits in ``<prefix>/bin`` (``helpers.fgsea_r_data.fake_prefix``); the wrapper
starts the real ``workflow/bin/giap_fgsea.R``, whose ``env -S`` shebang finds the fake first on the
child ``PATH``. Covered:

- the CLI contract: argv (``--vanilla``, the absolute script path, exactly the nine spec options
  with the ``function.gsea`` values and the seed), cwd inside ``tmp_root``, the three input files
  as ``write_call`` writes them;
- the child environment (``PATH`` = ``<prefix>/bin`` first, R library variables neutralised,
  ``LC_ALL=C``);
- non-zero exit → ``FgseaEngineError`` with the exit code, job count and the last 20 stderr
  lines; timeout → error, no process left behind; missing env → ``FgseaEngineUnavailable``;
  wrong fgsea version, missing or malformed output → ``FgseaEngineError``;
- the temporary directory is removed on success and on every failure;
- launcher resolution (lead decision 9: conda-only, ``GIAP_FGSEA_IMAGE`` is not a launcher).
"""

from __future__ import annotations

import hashlib
import math
import os
import time
from pathlib import Path

import numpy as np
import pytest
from helpers.fgsea_r_data import SCRIPT, fake_log, fake_prefix, options

REPO = Path(__file__).resolve().parents[3]
ENV_FILE = REPO / "workflow" / "envs" / "fgsea.yml"
SPEC_OPTIONS = {
    "--stats",
    "--sets",
    "--jobs",
    "--sample-size",
    "--eps",
    "--n-perm-simple",
    "--seed",
    "--out",
    "--versions",
}


def _call(seed: int = 424242, scores: list[float] | None = None):
    from giap.config.model import GseaConfig
    from giap.function.gsea import EngineCall, GseaJob

    ids = tuple(f"TU_{i}" for i in (5, 2, 7, 1, 4, 6, 3, 8))
    return EngineCall(
        tu_ids=ids,
        scores=np.array(scores or [3.0, 2.5, 1.0, 0.5, -0.25, -1.0, -2.0, -3.5], dtype=np.float64),
        jobs=(
            GseaJob("main:cog", {"C1": ("TU_1", "TU_2", "TU_3"), "C2": ("TU_4", "TU_5")}, 2, 7),
            GseaJob("main:kegg", {"K1": ("TU_6", "TU_7", "TU_8")}, 2, 7),
        ),
        params=GseaConfig(),
        seed=seed,
    )


@pytest.fixture
def work(tmp_path: Path) -> Path:
    d = tmp_path / "work"
    d.mkdir()
    return d


def _engine(prefix: Path, work: Path, timeout_s: float = 60):
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    launcher = resolve_launcher({"GIAP_FGSEA_PREFIX": str(prefix)})
    return r_fgsea_engine(launcher, timeout_s=timeout_s, tmp_root=work)


# ---------------------------------------------------------------- success and the CLI contract


def test_success_parses_rows_and_removes_tmp(tmp_path: Path, work: Path) -> None:
    prefix = fake_prefix(tmp_path / "env")
    engine = _engine(prefix, work)
    assert engine.method == "fgsea_multilevel"
    result = engine(_call())
    assert [r.term for r in result["main:cog"]] == ["C1", "C2"]
    assert [r.term for r in result["main:kegg"]] == ["K1"]
    c1 = result["main:cog"][0]
    assert (c1.size, c1.es, c1.nes, c1.p, c1.log2err) == (3, 0.5, 1.5, 0.01, 0.25)
    assert c1.leading_edge == ("TU_1", "TU_2")
    assert list(work.iterdir()) == []


def test_cli_contract(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import write_call

    prefix = fake_prefix(tmp_path / "env")
    call = _call(seed=1999999999)
    _engine(prefix, work)(call)
    log = fake_log(prefix)
    argv = log["argv"]
    assert argv[0] == "--vanilla"
    assert Path(argv[1]).resolve() == SCRIPT.resolve()
    assert Path(argv[1]).is_absolute()
    opts = options(argv)
    assert set(opts) == SPEC_OPTIONS
    assert opts["--sample-size"] == "101"
    assert opts["--n-perm-simple"] == "1000"
    assert float(opts["--eps"]) == 0.0
    assert opts["--seed"] == "1999999999"
    cwd = Path(log["cwd"])
    assert cwd.parent == work
    # the inputs the script saw are exactly write_call's files
    ref = tmp_path / "ref"
    ref.mkdir()
    write_call(call, ref)
    seen = log["files"]
    for name in ("stats.tsv", "jobs.tsv", "sets.tsv"):
        on_disk = (cwd / opts[f"--{name.split('.')[0]}"]).name
        assert seen[on_disk] == (ref / name).read_text(encoding="utf-8"), name


def test_cli_passes_config_values(tmp_path: Path, work: Path) -> None:
    from giap.config.model import GseaConfig
    from giap.function.gsea import EngineCall

    prefix = fake_prefix(tmp_path / "env")
    base = _call()
    call = EngineCall(
        tu_ids=base.tu_ids,
        scores=base.scores,
        jobs=base.jobs,
        params=GseaConfig(sample_size=51, eps=1e-10, n_perm_simple=2000),
        seed=base.seed,
    )
    _engine(prefix, work)(call)
    opts = options(fake_log(prefix)["argv"])
    assert opts["--sample-size"] == "51"
    assert opts["--n-perm-simple"] == "2000"
    assert float(opts["--eps"]) == 1e-10


def test_child_environment(tmp_path: Path, work: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("R_LIBS_USER", "/home/someone/R/library")
    monkeypatch.setenv("R_LIBS", "/home/someone/R/other")
    monkeypatch.setenv("R_PROFILE_USER", "/home/someone/.Rprofile")
    monkeypatch.setenv("LC_ALL", "en_US.UTF-8")
    prefix = fake_prefix(tmp_path / "env")
    _engine(prefix, work)(_call())
    env = fake_log(prefix)["env"]
    for var in ("R_LIBS_USER", "R_LIBS", "R_PROFILE_USER"):
        assert env.get(var) in ("", "/dev/null"), var
    assert env["LC_ALL"] == "C"
    assert env["PATH"].startswith(f"{prefix / 'bin'}{os.pathsep}")
    assert env["PATH"].endswith(os.environ["PATH"])


def test_non_finite_score_never_starts_r(tmp_path: Path, work: Path) -> None:
    prefix = fake_prefix(tmp_path / "env")
    with pytest.raises(ValueError):
        _engine(prefix, work)(_call(scores=[3.0, math.inf, 1.0, 0.5, -0.25, -1.0, -2.0, -3.5]))
    assert not (prefix / "log.json").exists()
    assert list(work.iterdir()) == []


# ---------------------------------------------------------------- failures


def test_exit_1_carries_code_jobs_and_stderr_tail(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    lines = [f"line{i:02d}" for i in range(1, 26)] + ["Error in fgseaMultilevel(...) : boom"]
    prefix = fake_prefix(
        tmp_path / "env", {"mode": "fail", "code": 1, "stderr": "\n".join(lines) + "\n"}
    )
    with pytest.raises(FgseaEngineError, match=r"R fgsea failed \(exit 1, 2 jobs\)") as info:
        _engine(prefix, work)(_call())
    msg = str(info.value)
    assert "Error in fgseaMultilevel(...) : boom" in msg
    assert "line07" in msg and "line25" in msg  # the last 20 lines
    assert "line06" not in msg
    assert list(work.iterdir()) == []


@pytest.mark.parametrize("code", [2, 3])
def test_usage_and_input_exit_codes_fail(tmp_path: Path, work: Path, code: int) -> None:
    from giap.io.fgsea_r import FgseaEngineError, FgseaEngineUnavailable

    prefix = fake_prefix(tmp_path / "env", {"mode": "fail", "code": code, "stderr": "bad\n"})
    with pytest.raises(FgseaEngineError, match=rf"exit {code}") as info:
        _engine(prefix, work)(_call())
    assert not isinstance(info.value, FgseaEngineUnavailable)
    assert list(work.iterdir()) == []


def _alive(pid: int) -> bool:
    try:
        state = Path(f"/proc/{pid}/stat").read_text().split(")")[-1].split()[0]
    except (FileNotFoundError, ProcessLookupError, IndexError):
        return False
    return state not in ("Z", "X")


def test_timeout_kills_process_group(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"mode": "sleep", "seconds": 20})
    t0 = time.monotonic()
    with pytest.raises(FgseaEngineError, match=r"R fgsea timed out after 1 s"):
        _engine(prefix, work, timeout_s=1)(_call())
    assert time.monotonic() - t0 < 10, "the wrapper waited for the orphaned grandchild"
    pid = int((prefix / "grandchild.pid").read_text())
    deadline = time.monotonic() + 5
    while _alive(pid) and time.monotonic() < deadline:
        time.sleep(0.05)
    assert not _alive(pid), "a child of the R process survived the timeout"
    assert list(work.iterdir()) == []


@pytest.mark.parametrize("version", ["1.38.0", "1.36.1"])
def test_wrong_fgsea_version(tmp_path: Path, work: Path, version: str) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"fgsea": version})
    with pytest.raises(FgseaEngineError, match=version.replace(".", r"\.")):
        _engine(prefix, work)(_call())
    assert list(work.iterdir()) == []


def test_missing_output(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"mode": "no_output"})
    with pytest.raises(FgseaEngineError):
        _engine(prefix, work)(_call())
    assert list(work.iterdir()) == []


def test_size_mismatch(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"size_offset": 1})
    with pytest.raises(FgseaEngineError):
        _engine(prefix, work)(_call())
    assert list(work.iterdir()) == []


def test_malformed_results(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"results": "pathway\tpval\nC1\t0.1\n"})
    with pytest.raises(FgseaEngineError):
        _engine(prefix, work)(_call())
    assert list(work.iterdir()) == []


def test_never_retries(tmp_path: Path, work: Path) -> None:
    """Lead decision 5: a failed R call is not retried; the R front end starts exactly once."""
    from giap.io.fgsea_r import FgseaEngineError

    prefix = fake_prefix(tmp_path / "env", {"mode": "fail", "code": 1, "stderr": "x\n"})
    with pytest.raises(FgseaEngineError):
        _engine(prefix, work)(_call())
    assert len((prefix / "calls.txt").read_text().splitlines()) == 1


def test_each_call_gets_a_fresh_directory(tmp_path: Path, work: Path) -> None:
    prefix = fake_prefix(tmp_path / "env")
    engine = _engine(prefix, work)
    engine(_call())
    engine(_call())
    dirs = (prefix / "calls.txt").read_text().splitlines()
    assert len(dirs) == 2 and dirs[0] != dirs[1]
    assert list(work.iterdir()) == []


def test_error_classes() -> None:
    from giap.io.fgsea_r import FgseaEngineError, FgseaEngineUnavailable
    from giap.steps.registry import StepError

    assert issubclass(FgseaEngineUnavailable, FgseaEngineError)
    assert issubclass(FgseaEngineError, StepError)  # step exit code 4 (spec §4)
    exc = FgseaEngineError("R fgsea timed out after 1 s")
    assert "R fgsea timed out after 1 s" in str(exc)
    assert "x" in str(FgseaEngineUnavailable("x"))


# ---------------------------------------------------------------- launcher resolution


def _built_root(root: Path, marker: str | None) -> Path:
    env = root / "micromamba" / "envs" / "giap-fgsea"
    fake_prefix(env)
    if marker is not None:
        (env / ".giap-built").write_text(marker + "\n", encoding="utf-8")
    return env


def _yml_sha() -> str:
    return hashlib.sha256(ENV_FILE.read_bytes()).hexdigest()


def test_launcher_prefix_wins(tmp_path: Path) -> None:
    from giap.io.fgsea_r import resolve_launcher

    prefix = fake_prefix(tmp_path / "explicit")
    _built_root(tmp_path / "root", _yml_sha())
    launcher = resolve_launcher(
        {"GIAP_FGSEA_PREFIX": str(prefix), "GIAP_ROOT": str(tmp_path / "root")}
    )
    assert launcher.kind == "prefix"
    assert Path(launcher.path) == prefix


def test_launcher_built_env_under_giap_root(tmp_path: Path) -> None:
    from giap.io.fgsea_r import resolve_launcher

    env = _built_root(tmp_path / "root", _yml_sha())
    launcher = resolve_launcher({"GIAP_ROOT": str(tmp_path / "root")})
    assert launcher.kind == "prefix"
    assert Path(launcher.path) == env


def test_launcher_built_env_runs(tmp_path: Path, work: Path) -> None:
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    env = _built_root(tmp_path / "root", _yml_sha())
    engine = r_fgsea_engine(resolve_launcher({"GIAP_ROOT": str(tmp_path / "root")}), tmp_root=work)
    assert set(engine(_call())) == {"main:cog", "main:kegg"}
    assert fake_log(env)["env"]["PATH"].startswith(str(env / "bin"))


@pytest.mark.parametrize("marker", [None, "0" * 64, ""])
def test_launcher_unbuilt_env_is_unavailable(tmp_path: Path, marker: str | None) -> None:
    from giap.io.fgsea_r import FgseaEngineUnavailable, resolve_launcher

    _built_root(tmp_path / "root", marker)
    with pytest.raises(FgseaEngineUnavailable) as info:
        resolve_launcher({"GIAP_ROOT": str(tmp_path / "root")})
    msg = str(info.value)
    assert "giap-fgsea" in msg
    assert "sbatch scripts/slurm/build_envs.sbatch fgsea" in msg


def test_launcher_nothing_set_is_unavailable(tmp_path: Path) -> None:
    from giap.io.fgsea_r import FgseaEngineUnavailable, resolve_launcher

    with pytest.raises(FgseaEngineUnavailable):
        resolve_launcher({"GIAP_ROOT": str(tmp_path / "empty")})


def test_image_variable_is_not_a_launcher(tmp_path: Path) -> None:
    """Lead decision 9: conda-only; an image path never makes the engine available."""
    from giap.io.fgsea_r import FgseaEngineUnavailable, resolve_launcher

    img = tmp_path / "fgsea.img"
    img.write_bytes(b"\0")
    with pytest.raises(FgseaEngineUnavailable):
        resolve_launcher({"GIAP_FGSEA_IMAGE": str(img), "GIAP_ROOT": str(tmp_path / "empty")})
    prefix = fake_prefix(tmp_path / "env")
    launcher = resolve_launcher({"GIAP_FGSEA_IMAGE": str(img), "GIAP_FGSEA_PREFIX": str(prefix)})
    assert launcher.kind == "prefix" and Path(launcher.path) == prefix


def test_missing_prefix_is_unavailable(tmp_path: Path, work: Path) -> None:
    """A ``GIAP_FGSEA_PREFIX`` pointing nowhere: unavailable at resolution or at the call."""
    from giap.io.fgsea_r import FgseaEngineUnavailable, r_fgsea_engine, resolve_launcher

    with pytest.raises(FgseaEngineUnavailable):
        launcher = resolve_launcher({"GIAP_FGSEA_PREFIX": str(tmp_path / "nowhere")})
        r_fgsea_engine(launcher, tmp_root=work)(_call())
    assert list(work.iterdir()) == []


def test_prefix_without_r_fails(tmp_path: Path, work: Path) -> None:
    """An existing prefix with no R front end: the call fails, never silently succeeds."""
    from giap.io.fgsea_r import FgseaEngineError, r_fgsea_engine, resolve_launcher

    empty = tmp_path / "empty_env"
    (empty / "bin").mkdir(parents=True)
    with pytest.raises(FgseaEngineError):
        launcher = resolve_launcher({"GIAP_FGSEA_PREFIX": str(empty)})
        r_fgsea_engine(launcher, tmp_root=work)(_call())
    assert list(work.iterdir()) == []
