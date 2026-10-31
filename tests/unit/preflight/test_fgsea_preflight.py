"""M8-17b T10 and A9: fgsea in the tool registry, preflight and ``gippa run``
(docs/specs/M8-17-r-fgsea-engine.md §10; lead decisions 8 and 9; lead notes for M8-17b 1).

Interface pinned here (test-author reading of §10 with the lead decisions):

- ``gippa.annotate.tools.TOOLS["fgsea"]`` is the M8-17a ``FGSEA_TOOL`` entry: version 1.36.2, MIT,
  ``workflow/envs/fgsea.yml``, env ``gippa-fgsea``, module ``function``, an R ``version_cmd`` that
  never spells the R script front end. Lead decision 9: conda-only under every profile, so no
  image (``image`` absent or ``None``) and no ``containers.config`` pin.
- Lead decision 8: preflight requires fgsea only when ``function`` is in ``modules`` **and**
  abundance is given (``--abundance``). The check is the ostir one: ``environment.tools``, level
  ``error``, exit 2, ``value["missing"]`` names the tool, and ``gippa run`` exits 2 without
  launching Nextflow.
- Missing means ``$GIPPA_ROOT/micromamba/envs/gippa-fgsea/.gippa-built`` does not equal the sha256
of
  ``workflow/envs/fgsea.yml``, under farmshare and under apptainer profiles alike; the ``test`` and
  ``stub`` profiles require nothing.

Fakes only (``helpers.fake_env``); never the real $GIPPA_DB_DIR or $GIPPA_ROOT.
"""

from __future__ import annotations

import json
import os
import random
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import (
    CONTAINERS_CONFIG,
    MARKER,
    conda_only_tools,
    env_yml_sha256,
    fake_gippa_root,
    write_setup_manifest,
)

from gippa.config import resolve_config
from gippa.exitcodes import ExitCode

R_FRONT_END = "R" + "script"
PROFILES = [
    "farmshare",
    "farmshare_local",
    "farmshare,apptainer",
    "farmshare_local,apptainer",
    "apptainer",
]
NO_TOOL_PROFILES = ["test", "stub", "test,stub"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


# ---------------------------------------------------------------- the registry


def test_fgsea_is_a_registry_tool() -> None:
    from gippa.annotate.tools import TOOLS

    t = TOOLS["fgsea"]
    assert t.name == "fgsea"
    assert t.version == "1.36.2"
    assert t.licence == "MIT"
    assert t.env_file == "workflow/envs/fgsea.yml"
    assert t.env_name == "gippa-fgsea"
    assert t.module == "function"
    assert t.in_process is False
    assert t.version_cmd == (
        "R",
        "--vanilla",
        "--slave",
        "-e",
        "cat(as.character(packageVersion('fgsea')))",
    )
    assert not any(R_FRONT_END in w for w in t.version_cmd)


def test_registry_entry_is_the_wrapper_entry() -> None:
    from gippa.annotate.tools import TOOLS
    from gippa.io.fgsea_r import FGSEA_TOOL, FGSEA_VERSION

    assert TOOLS["fgsea"] == FGSEA_TOOL
    assert TOOLS["fgsea"].version == FGSEA_VERSION


def test_fgsea_is_conda_only() -> None:
    """Lead decision 9: no image pin; the BioContainers digest is deferred (backlog)."""
    from gippa.annotate.tools import TOOLS, read_containers_config

    assert getattr(TOOLS["fgsea"], "image", None) is None
    pinned = {p.lower() for p in read_containers_config(CONTAINERS_CONFIG)}
    assert "fgsea" not in pinned
    assert "fgsea" in conda_only_tools()


def _cfg(*over: str) -> Any:
    return resolve_config(overrides=list(over))


def _required(cfg: Any, abundance: bool) -> set[str]:
    """``required_tools`` under the lead decision 8 rule, whatever the abundance keyword is."""
    import inspect

    from gippa.annotate.tools import required_tools

    params = inspect.signature(required_tools).parameters
    extra = [p for p in params if p != "cfg"]
    assert extra, "required_tools needs to know whether abundance is given (lead decision 8)"
    return {t.name for t in required_tools(cfg, **{extra[0]: abundance})}


def test_required_with_function_and_abundance() -> None:
    assert "fgsea" in _required(_cfg(), abundance=True)


def test_not_required_without_abundance() -> None:
    assert "fgsea" not in _required(_cfg(), abundance=False)


def test_not_required_without_function() -> None:
    assert "fgsea" not in _required(_cfg("modules=[core,mobile]"), abundance=True)


def test_required_tools_without_abundance_is_unchanged() -> None:
    from gippa.annotate.tools import required_tools

    assert "fgsea" not in {t.name for t in required_tools(_cfg())}


# ---------------------------------------------------------------- missing_tools


def _missing(tools: list[Any], profile: str, root: Path) -> list[str]:
    """``missing_tools`` with a fake ``$GIPPA_ROOT``'s image cache and env directory (14f1)."""
    from gippa.preflight.environment import missing_tools

    return missing_tools(tools, profile, root / "containers", root / "micromamba" / "envs", None)


@pytest.mark.parametrize("profile", PROFILES)
def test_missing_when_env_not_built(tmp_path: Path, profile: str) -> None:
    from gippa.annotate.tools import TOOLS

    root = fake_gippa_root(tmp_path / "root", skip=["fgsea"])
    assert _missing([TOOLS["fgsea"]], profile, root) == ["fgsea"]


@pytest.mark.parametrize("profile", PROFILES)
def test_present_when_env_built(tmp_path: Path, profile: str) -> None:
    from gippa.annotate.tools import TOOLS

    root = fake_gippa_root(tmp_path / "root")
    marker = root / "micromamba" / "envs" / "gippa-fgsea" / MARKER
    assert marker.read_text().strip() == env_yml_sha256("fgsea")
    assert _missing([TOOLS["fgsea"]], profile, root) == []


@pytest.mark.parametrize("how", ["stale", "no_marker"])
def test_stale_or_unmarked_env_is_missing(tmp_path: Path, how: str) -> None:
    from gippa.annotate.tools import TOOLS

    root = fake_gippa_root(tmp_path / "root", **{how: ["fgsea"]})
    assert _missing([TOOLS["fgsea"]], "farmshare,apptainer", root) == ["fgsea"]


@pytest.mark.parametrize("profile", NO_TOOL_PROFILES)
def test_never_missing_under_test_or_stub(tmp_path: Path, profile: str) -> None:
    from gippa.annotate.tools import TOOLS

    empty = tmp_path / "empty"
    empty.mkdir()
    assert _missing([TOOLS["fgsea"]], profile, empty) == []


# ---------------------------------------------------------------- run_preflight and gippa run


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "MAMBA_ROOT_PREFIX", "NXF_APPTAINER_CACHEDIR"):
        empty = tmp_path / f"unset_{var.lower()}"
        empty.mkdir()
        monkeypatch.setenv(var, str(empty))
    # 14f1: a set cache directory wins over $GIPPA_ROOT/containers; unset, never the session's
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR")
    monkeypatch.delenv("GIPPA_FGSEA_PREFIX", raising=False)


def _acgt(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _files(tmp_path: Path) -> dict[str, Path]:
    d = tmp_path / "in"
    d.mkdir(parents=True, exist_ok=True)
    g = d / "g.fa"
    g.write_text(f">chr1 [topology=circular]\n{_acgt(120_000)}\n")
    p = d / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    a = d / "abundance.tsv"
    a.write_text("protein_id\tabundance\nP1\t1000\n")
    return {"genome": g, "proteins": p, "abundance": a}


def _tools_check(
    tmp_path: Path, root: Path, profile: str, *, abundance: bool, modules: str | None = None
) -> Any:
    from gippa.preflight.run import PreflightInputs, run_preflight

    f = _files(tmp_path)
    inp = PreflightInputs(
        genome=f["genome"], proteins=f["proteins"], abundance=f["abundance"] if abundance else None
    )
    over = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
    if modules is not None:
        over.append(f"modules={modules}")
    report, _ = run_preflight(
        inp,
        resolve_config(overrides=over),
        None,
        "run",
        tmp_path / "OUT",
        write_setup_manifest(tmp_path / "db"),
        gippa_root=root,
        profile=profile,
    )
    (check,) = [c for c in report.checks if c.rule_id == "environment.tools"]
    return check


@pytest.mark.parametrize("profile", ["farmshare", "farmshare,apptainer"])
def test_preflight_reports_a_missing_fgsea_env(tmp_path: Path, profile: str) -> None:
    """fgsea needs its prebuilt env under both profiles, with conda on (no image is used)."""
    fg = _tools_check(
        tmp_path / "fg", fake_gippa_root(tmp_path / "r1", skip=["fgsea"]), profile, abundance=True
    )
    assert fg.level == "error"
    assert fg.exit_code == 2
    assert fg.rule_id == "environment.tools"
    assert fg.value == {"missing": ["fgsea"]}
    assert fg.message == "missing: fgsea; fgsea: prebuilt env gippa-fgsea: not built"


@pytest.mark.parametrize("profile", ["farmshare", "farmshare,apptainer"])
def test_preflight_without_abundance_does_not_need_fgsea(tmp_path: Path, profile: str) -> None:
    check = _tools_check(
        tmp_path, fake_gippa_root(tmp_path / "root", skip=["fgsea"]), profile, abundance=False
    )
    assert check.level == "ok", check


def test_preflight_without_function_does_not_need_fgsea(tmp_path: Path) -> None:
    check = _tools_check(
        tmp_path,
        fake_gippa_root(tmp_path / "root", skip=["fgsea"]),
        "farmshare,apptainer",
        abundance=True,
        modules="[core]",
    )
    assert check.level == "ok", check


@pytest.mark.parametrize("profile", NO_TOOL_PROFILES)
def test_preflight_test_and_stub_never_need_fgsea(tmp_path: Path, profile: str) -> None:
    empty = tmp_path / "empty_root"
    empty.mkdir()
    check = _tools_check(tmp_path, empty, profile, abundance=True)
    assert check.level == "ok", check


def test_preflight_built_fgsea_is_ok(tmp_path: Path) -> None:
    check = _tools_check(
        tmp_path, fake_gippa_root(tmp_path / "root"), "farmshare,apptainer", abundance=True
    )
    assert check.level == "ok", check


FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    return args


def _gippa_run(tmp_path: Path, *extra: str, abundance: bool = True) -> tuple[int, str]:
    from gippa.cli import cli

    f = _files(tmp_path)
    argv = ["run", "--genome", str(f["genome"]), "--proteins", str(f["proteins"])]
    if abundance:
        argv += ["--abundance", str(f["abundance"])]
    argv += ["--outdir", str(tmp_path / "OUT"), *extra]
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


@pytest.mark.parametrize("profile", ["farmshare", "farmshare,apptainer"])
def test_gippa_run_refuses_without_fgsea(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path, profile: str
) -> None:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "root", skip=["fgsea"])))
    monkeypatch.setenv("GIPPA_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    code, text = _gippa_run(tmp_path, "--profile", profile)
    assert code == ExitCode.ENVIRONMENT, text
    assert "fgsea" in text
    assert not fake_nf.exists(), "nextflow launched despite a missing fgsea env"
    report = json.loads((tmp_path / "OUT" / "preflight" / "preflight.json").read_text())
    (tools,) = [c for c in report["checks"] if c["rule_id"] == "environment.tools"]
    assert tools["value"]["missing"] == ["fgsea"]


def test_gippa_run_without_abundance_launches_without_fgsea(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path
) -> None:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "root", skip=["fgsea"])))
    monkeypatch.setenv("GIPPA_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    code, text = _gippa_run(tmp_path, "--profile", "farmshare,apptainer", abundance=False)
    assert code == ExitCode.SUCCESS, text
    assert fake_nf.is_file()
