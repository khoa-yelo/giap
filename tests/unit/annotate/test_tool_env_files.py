"""Conda env files and container pins (M2-01; GIPPA_SPEC.md §3.2, §4, §6.3, §7.1.1)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[3]
ENVS = REPO / "workflow" / "envs"
CONTAINERS = REPO / "workflow" / "conf" / "containers.config"
DIGEST = re.compile(r"^[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}$")
HEX64 = "a" * 64

EXPECTED_ENVS = {
    "bakta",
    "uniop",
    "promotech",
    "promoter_calculator",
    "terminatornet",
    "transtermhp",
    "bactermfinder",
    "ostir",
    "genomad",
    "antismash",
    "miniprot",
    "blast",
    "fgsea",  # M8-17a (ADR-0038)
}
# Q3 (lead decision 1): no upstream image yet -> no withName entry; conda even under apptainer.
CONDA_ONLY = {
    "uniop",
    "promotech",
    "promoter_calculator",
    "terminatornet",
    "bactermfinder",
    "ostir",
    "fgsea",  # M8-17 lead decision 9: conda-only under every profile
}


def _env_files() -> list[Path]:
    return sorted(ENVS.glob("*.yml"))


def _write(tmp_path: Path, stem: str, body: str) -> Path:
    path = tmp_path / f"{stem}.yml"
    path.write_text(body)
    return path


# ---------------------------------------------------------------- read_env_file / env_problems


def test_example1_bakta_env() -> None:
    from gippa.annotate.tools import EnvFile, env_problems, read_env_file

    env = read_env_file(ENVS / "bakta.yml")
    assert env == EnvFile("gippa-bakta", ("conda-forge", "bioconda"), ("bakta=1.12.2",), ())
    assert env_problems(env) == []


def test_example2_unpinned_conda(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [conda-forge, bioconda]\ndependencies:\n  - pandas\n",
    )
    assert env_problems(read_env_file(path)) == ["pandas: conda dependency without exact version"]


def test_example3_defaults_channel(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [defaults, bioconda]\ndependencies:\n  - samtools=1.21\n",
    )
    assert "channel 'defaults' is not allowed" in env_problems(read_env_file(path))


def test_example4_pip_range(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [conda-forge, bioconda]\ndependencies:\n"
        "  - python=3.11.9\n  - pip=24.2\n  - pip:\n    - ostir>=1.1\n",
    )
    env = read_env_file(path)
    assert env.pip == ("ostir>=1.1",)
    assert env_problems(env) == ["ostir>=1.1: pip dependency must be pinned with =="]


def test_pip_section_parsed(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [bioconda, conda-forge]\ndependencies:\n"
        "  - python=3.11.9\n  - pip=24.2\n  - pip:\n    - ostir==1.1.3\n",
    )
    env = read_env_file(path)
    assert env.conda == ("python=3.11.9", "pip=24.2")
    assert env.pip == ("ostir==1.1.3",)
    assert env_problems(env) == []  # channel order is free


def test_build_string_pin_allowed(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [conda-forge, bioconda]\ndependencies:\n"
        "  - python=3.6.10=h0371630_0\n",
    )
    assert env_problems(read_env_file(path)) == []


@pytest.mark.parametrize(
    "channels", ["[conda-forge, bioconda, r]", "[conda-forge]", "[conda-forge, bioconda, defaults]"]
)
def test_channel_set_must_be_exact(tmp_path: Path, channels: str) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        f"name: gippa-x\nchannels: {channels}\ndependencies:\n  - samtools=1.21\n",
    )
    assert env_problems(read_env_file(path))


def test_unpinned_pip_rejected(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-x\nchannels: [conda-forge, bioconda]\ndependencies:\n"
        "  - python=3.11.9\n  - pip:\n    - pandas\n",
    )
    assert env_problems(read_env_file(path))


def test_env_name_must_match_stem(tmp_path: Path) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    path = _write(
        tmp_path,
        "x",
        "name: gippa-y\nchannels: [conda-forge, bioconda]\ndependencies:\n  - samtools=1.21\n",
    )
    assert env_problems(read_env_file(path))


# ---------------------------------------------------------------- the committed env files


# M10-08: pilot-only envs (gippa.annotate.tools.PILOT_TOOLS), not pipeline tools; their own
# rules (conda-forge only, version ranges) are tested in tests/unit/test_m10_08_registry.py.
PILOT_ENVS = {"cascadia", "timsconvert"}


def test_all_env_files_present() -> None:
    assert {p.stem for p in _env_files()} - PILOT_ENVS == EXPECTED_ENVS


@pytest.mark.parametrize("stem", sorted(EXPECTED_ENVS))
def test_env_file_static(stem: str) -> None:
    path = ENVS / f"{stem}.yml"
    data = yaml.safe_load(path.read_text())
    assert data["name"] == f"gippa-{stem}"
    assert sorted(data["channels"]) == ["bioconda", "conda-forge"]
    assert "defaults" not in path.read_text().split("dependencies:")[0]


@pytest.mark.parametrize("stem", sorted(EXPECTED_ENVS))
def test_env_file_passes_env_problems(stem: str) -> None:
    from gippa.annotate.tools import env_problems, read_env_file

    assert env_problems(read_env_file(ENVS / f"{stem}.yml")) == []


def _deps(stem: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    from gippa.annotate.tools import read_env_file

    env = read_env_file(ENVS / f"{stem}.yml")
    return env.conda, env.pip


@pytest.mark.parametrize(
    ("stem", "package"),
    [
        ("bakta", "bakta"),
        ("transtermhp", "transtermhp"),
        ("genomad", "genomad"),
        ("antismash", "antismash"),
        ("miniprot", "miniprot"),
        ("blast", "blast"),
    ],
)
def test_conda_pin_matches_registry(stem: str, package: str) -> None:
    from gippa.annotate.tools import TOOLS

    conda, _ = _deps(stem)
    assert f"{package}={TOOLS[stem].version}" in {"=".join(d.split("=")[:2]) for d in conda}


@pytest.mark.parametrize(
    ("stem", "package"), [("ostir", "ostir"), ("promoter_calculator", "promotercalculator")]
)
def test_pip_pin_matches_registry(stem: str, package: str) -> None:
    from gippa.annotate.tools import TOOLS

    _, pip = _deps(stem)
    assert f"{package}=={TOOLS[stem].version}" in pip


def _has(stem: str, package: str, version_prefix: str) -> bool:
    conda, pip = _deps(stem)
    conda_hit = any(d.lower().startswith(f"{package}={version_prefix}") for d in conda)
    pip_hit = any(d.lower().startswith(f"{package}=={version_prefix}") for d in pip)
    return conda_hit or pip_hit


def test_ostir_env_has_viennarna() -> None:
    conda, _ = _deps("ostir")
    assert any(d.lower().startswith("viennarna=") for d in conda)


def test_promotech_frozen_env() -> None:
    assert _has("promotech", "python", "3.6.10")
    assert _has("promotech", "tensorflow", "2.1.0")
    assert _has("promotech", "scikit-learn", "0.23.1")


def test_uniop_env() -> None:
    assert _has("uniop", "python", "3.11")
    for pkg in ("pandas", "numpy", "scikit-learn"):
        assert _has("uniop", pkg, ""), pkg


def test_terminatornet_sklearn_pin_documented() -> None:
    assert _has("terminatornet", "scikit-learn", "")
    text = (ENVS / "terminatornet.yml").read_text()
    comments = "\n".join(line for line in text.splitlines() if line.lstrip().startswith("#"))
    assert "model.pkl" in comments


# ---------------------------------------------------------------- containers.config


def test_example8_read_containers_config(tmp_path: Path) -> None:
    from gippa.annotate.tools import read_containers_config

    path = tmp_path / "containers.config"
    ref = f"docker://oschwengers/bakta@sha256:{HEX64}"
    path.write_text(f"process {{\n  withName: 'BAKTA' {{ container = '{ref}' }}\n}}\n")
    assert read_containers_config(path) == {"BAKTA": ref}


def test_containers_config_digests_only() -> None:
    from gippa.annotate.tools import read_containers_config

    pins = read_containers_config(CONTAINERS)
    assert pins
    for process, ref in pins.items():
        assert DIGEST.match(ref), f"{process}: {ref}"
    # No container assignment anywhere in the file without a digest (no tag-only refs).
    for line in CONTAINERS.read_text().splitlines():
        if re.search(r"\bcontainer\s*=", line) and not line.lstrip().startswith("//"):
            assert "@sha256:" in line, line


def test_containers_config_named_images() -> None:
    from gippa.annotate.tools import read_containers_config

    refs = list(read_containers_config(CONTAINERS).values())
    assert any(r.startswith("oschwengers/bakta@") for r in refs)
    assert any(r.startswith("antismash/standalone@") for r in refs)


def _norm(s: str) -> str:
    return s.upper().replace("_", "")


def test_conda_only_tools_have_no_container_and_are_listed() -> None:
    from gippa.annotate.tools import read_containers_config

    keys = [_norm(k) for k in read_containers_config(CONTAINERS)]
    text = CONTAINERS.read_text()
    for tool in CONDA_ONLY:
        assert not any(_norm(tool) in k for k in keys), tool
        assert tool in text, f"{tool} not listed as conda-only in containers.config comments"


def test_every_tool_is_containerised_or_conda_only() -> None:
    from gippa.annotate.tools import TOOLS, read_containers_config

    keys = [_norm(k) for k in read_containers_config(CONTAINERS)]
    alias = {"transtermhp": "TRANSTERM"}
    for name, spec in TOOLS.items():
        if spec.in_process or name in CONDA_ONLY:
            continue
        needle = alias.get(name, _norm(name))
        assert any(needle in k for k in keys), f"{name}: no withName entry"
