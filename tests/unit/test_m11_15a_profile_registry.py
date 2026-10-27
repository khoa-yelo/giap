"""M11-15a: the profile registry parsed from ``nextflow.config`` and the generic ``slurm`` profile.

Spec: docs/specs/M11-15-portable-clusters.md item 3, "Profile registry and combination rules
(15a)" and the ``slurm.config`` block (15a: the profile, reading the five ``slurm_*`` params with
``null`` defaults; the flags that set them are 15e); acceptance criteria 1 and 5 of item 3 (the
combination errors) and the CLI help (item 1: the default ``auto`` shown in ``--help``).

Interface pinned here:

- ``gippa.runner.profiles.parse_profiles(text) -> tuple[str, ...]``: the profile names of the
  ``profiles { }`` block of a ``nextflow.config`` text, in file order (comments ignored).
- ``gippa.runner.profiles.accepted_profiles() -> tuple[str, ...]``: the bundled
  ``workflow/nextflow.config`` names without ``test`` and ``stub``; the CLI accepts these
  (``_PROFILES`` is gone).
- Combination errors are usage errors (exit 2) of ``gippa run`` and ``gippa preflight``:
  ``more than one executor profile: slurm, farmshare`` (names in the order given) and
  ``more than one container runtime: apptainer, docker`` (Q-L16; the wording after the colon
  mirrors the executor message).

The ``slurm.config`` checks reuse the text-reader approach of test_m11_14_profile_parity.py.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
NEXTFLOW = WF / "nextflow.config"
SLURM = WF / "conf" / "slurm.config"
SLURM_PARAMS = ("slurm_queue", "slurm_account", "slurm_options", "slurm_max_jobs",
                "slurm_max_time")  # fmt: skip
OLD_PROFILES = ("farmshare", "farmshare_local", "apptainer", "conda", "docker", "standard")


def _mod() -> Any:
    try:
        from gippa.runner import profiles
    except ImportError as exc:
        pytest.fail(f"gippa.runner.profiles is not implemented (M11-15a): {exc}")
    return profiles


def _strip_comments(text: str) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)


def _block(text: str, header: str) -> str:
    m = re.search(rf"(?m)^\s*{header}\s*\{{", text)
    assert m, f"no block {header!r}"
    start = text.index("{", m.start())
    depth = 0
    for i in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[start + 1 : i]
    raise AssertionError(f"unbalanced block {header!r}")


def _cli(*argv: str) -> Any:
    from gippa.cli import cli

    return CliRunner().invoke(cli, list(argv))


def _out(res: Any) -> str:
    return res.output + (res.stderr if res.stderr_bytes is not None else "")


# ---------------------------------------------------------------- registry


def test_parse_profiles_fixture_text() -> None:
    text = (
        "params { x = 1 }\n"
        "profiles {\n"
        "    standard { includeConfig 'conf/standard.config' }\n"
        "    // commented { }\n"
        "    stub { conda.enabled = false }\n"
        "    acme { includeConfig 'conf/sites/acme.config' }\n"
        "    docker {\n        docker.enabled = true\n        process { withName: 'X' { a = 1 } }\n"
        "    }\n"
        "}\n"
    )
    assert _mod().parse_profiles(text) == ("standard", "stub", "acme", "docker")


def test_bundled_registry() -> None:
    names = _mod().accepted_profiles()
    assert set(OLD_PROFILES) <= set(names)
    assert "slurm" in names
    assert "test" not in names and "stub" not in names
    parsed = _mod().parse_profiles(NEXTFLOW.read_text(encoding="utf-8"))
    assert set(names) == set(parsed) - {"test", "stub"}, "every bundled profile, nothing else"


def test_cli_profiles_tuple_is_gone() -> None:
    import gippa.cli

    assert not hasattr(gippa.cli, "_PROFILES"), "the registry replaces the hard-coded tuple"


@pytest.mark.parametrize("cmd", ["run", "preflight"])
@pytest.mark.parametrize("value", ["test", "stub", "bogus", "slurm,test"])
def test_rejected_names(cmd: str, value: str, tmp_path: Path) -> None:
    g = tmp_path / "g.fa"
    g.write_text(">c\nACGT\n")
    argv = [cmd, "--genome", str(g), "--profile", value]
    res = _cli(*argv, *(["--outdir", str(tmp_path / "o")] if cmd == "run" else []))
    assert res.exit_code == 2, _out(res)
    assert "unknown profile" in _out(res)


@pytest.mark.parametrize("cmd", ["run", "preflight"])
@pytest.mark.parametrize(
    ("value", "message"),
    [
        ("slurm,farmshare", "more than one executor profile: slurm, farmshare"),
        ("farmshare,farmshare_local,apptainer",
         "more than one executor profile: farmshare, farmshare_local"),
        ("standard,slurm,apptainer", "more than one executor profile: standard, slurm"),
        ("apptainer,docker", "more than one container runtime: apptainer, docker"),
        ("slurm,docker,apptainer", "more than one container runtime: docker, apptainer"),
    ],
)  # fmt: skip
def test_combination_errors(cmd: str, value: str, message: str, tmp_path: Path) -> None:
    g = tmp_path / "g.fa"
    g.write_text(">c\nACGT\n")
    argv = [cmd, "--genome", str(g), "--profile", value]
    res = _cli(*argv, *(["--outdir", str(tmp_path / "o")] if cmd == "run" else []))
    assert res.exit_code == 2, _out(res)
    assert message in _out(res)
    assert not (tmp_path / "o").exists(), "a usage error before anything is written"


def test_run_help_shows_auto_and_slurm() -> None:
    res = _cli("run", "--help")
    assert res.exit_code == 0
    text = " ".join(res.output.split())
    assert re.search(r"--profile.*default: auto", text), text
    assert "slurm" in text


@pytest.mark.parametrize("cmd", ["preflight", "setup"])
def test_preflight_and_setup_help_list_slurm_or_auto(cmd: str) -> None:
    res = _cli(cmd, "--help")
    text = " ".join(res.output.split())
    assert "auto" in text
    if cmd == "preflight":
        assert "none" in text and "slurm" in text


# ---------------------------------------------------------------- slurm.config (static)


def test_slurm_profile_entry() -> None:
    body = _block(_strip_comments(NEXTFLOW.read_text(encoding="utf-8")), "profiles")
    assert re.search(r"(?m)^\s*slurm\s*\{\s*includeConfig\s+'conf/slurm.config'\s*\}", body)


def test_slurm_params_default_null() -> None:
    params = _block(_strip_comments(NEXTFLOW.read_text(encoding="utf-8")), "params")
    for name in SLURM_PARAMS:
        assert re.search(rf"(?m)^\s*{name}\s*=\s*null\s*$", params), f"params.{name} = null"


def test_slurm_config_generic() -> None:
    raw = SLURM.read_text(encoding="utf-8")
    text = _strip_comments(raw)
    process = _block(text, "process")
    executor = _block(text, "executor")
    assert re.search(r"executor\s*=\s*'slurm'", process)
    for word in ("qos", "exclude", "4000.MB", "scratch", "prebuilt_envs", "barley", "normal"):
        assert word not in text, f"nothing FarmShare-specific in slurm.config: {word}"
    assert "resourceLimits" in process
    assert re.search(r"jobName\s*=\s*\{\s*\"gippa-\"", executor)
    assert "submitRateLimit" in executor
    assert re.search(r"queueSize\s*=\s*\(\s*params\.slurm_max_jobs\s*\?:\s*20\s*\)", executor)


def test_slurm_config_reads_exactly_the_five_params() -> None:
    text = _strip_comments(SLURM.read_text(encoding="utf-8"))
    read = set(re.findall(r"params\.(slurm_[a-z_]+)", text))
    assert read == set(SLURM_PARAMS)
    assert "includeConfig" not in text
