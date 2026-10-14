"""Tool registry `giap.annotate.tools` (M2-01; GIAP_SPEC.md §7.1.1, R5)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
SHA40 = re.compile(r"^git:[0-9a-f]{40}$")

ALL_TOOLS = {
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
    "codon_bias",
}
GIT_TOOLS = {"uniop", "promotech", "terminatornet", "bactermfinder"}
# §7.1.1 exact versions.
VERSIONS = {
    "bakta": "1.12.2",
    "promoter_calculator": "1.2.4",
    "transtermhp": "2.09",
    "ostir": "1.1.3",
    "codon_bias": "0.5.0",
    "genomad": "1.12.0",
    "antismash": "8.0.4",
    "miniprot": "0.18",
}
# R5 / NOTICE.md (licences verified upstream at the pinned commits).
LICENCES = {
    "bakta": "GPL-3.0",
    "antismash": "AGPL-3.0",
    "genomad": "non-commercial",
    "miniprot": "MIT",
    "codon_bias": "MIT",
    "uniop": "unlicensed",
    "terminatornet": "MIT",
    "bactermfinder": "GPL-3.0",
}
GPL_FAMILY = {"promotech", "promoter_calculator", "transtermhp", "ostir"}
MODULES = {
    "bakta": "core",
    "uniop": "core",
    "promotech": "core",
    "promoter_calculator": "core",
    "terminatornet": "core",
    "transtermhp": "core",
    "ostir": "core",
    "codon_bias": "core",
    "bactermfinder": "opt_in",
    "genomad": "mobile",
    "antismash": "bgc",
    "miniprot": "mapping",
    "blast": "mapping",
}


@pytest.fixture(scope="module")
def tools():  # type: ignore[no-untyped-def]
    from giap.annotate.tools import TOOLS

    return TOOLS


def test_registry_has_every_7_1_1_tool(tools) -> None:  # type: ignore[no-untyped-def]
    assert set(tools) == ALL_TOOLS
    for key, spec in tools.items():
        assert spec.name == key


@pytest.mark.parametrize(("name", "version"), sorted(VERSIONS.items()))
def test_exact_versions(tools, name: str, version: str) -> None:  # type: ignore[no-untyped-def]
    assert tools[name].version == version


def test_blast_version_is_pinned(tools) -> None:  # type: ignore[no-untyped-def]
    assert re.fullmatch(r"\d+\.\d+\.\d+", tools["blast"].version)


@pytest.mark.parametrize("name", sorted(GIT_TOOLS))
def test_git_tools_pinned_to_full_sha(tools, name: str) -> None:  # type: ignore[no-untyped-def]
    spec = tools[name]
    assert SHA40.match(spec.version), spec.version
    assert spec.source_repo and spec.source_repo.startswith("https://")


def test_only_git_tools_have_git_versions(tools) -> None:  # type: ignore[no-untyped-def]
    for name, spec in tools.items():
        assert spec.version.startswith("git:") == (name in GIT_TOOLS), name


@pytest.mark.parametrize(("name", "licence"), sorted(LICENCES.items()))
def test_licences(tools, name: str, licence: str) -> None:  # type: ignore[no-untyped-def]
    assert tools[name].licence == licence


@pytest.mark.parametrize("name", sorted(GPL_FAMILY))
def test_gpl_tools(tools, name: str) -> None:  # type: ignore[no-untyped-def]
    assert tools[name].licence.startswith("GPL")


def test_every_licence_non_empty(tools) -> None:  # type: ignore[no-untyped-def]
    for spec in tools.values():
        assert spec.licence.strip()


@pytest.mark.parametrize(("name", "module"), sorted(MODULES.items()))
def test_modules(tools, name: str, module: str) -> None:  # type: ignore[no-untyped-def]
    assert tools[name].module == module


def test_env_paths_and_names(tools) -> None:  # type: ignore[no-untyped-def]
    for name, spec in tools.items():
        if spec.in_process:
            continue
        assert spec.env_file == f"workflow/envs/{name}.yml"
        assert spec.env_name == f"giap-{name}"
        assert (REPO / spec.env_file).is_file(), spec.env_file
        assert spec.version_cmd and all(isinstance(a, str) and a for a in spec.version_cmd)


def test_only_codon_bias_in_process(tools) -> None:  # type: ignore[no-untyped-def]
    assert {n for n, s in tools.items() if s.in_process} == {"codon_bias"}


def test_every_env_file_has_a_registry_entry(tools) -> None:  # type: ignore[no-untyped-def]
    stems = {p.stem for p in (REPO / "workflow" / "envs").glob("*.yml")}
    out_of_process = {n for n, s in tools.items() if not s.in_process}
    assert stems == out_of_process


def test_toolspec_is_frozen(tools) -> None:  # type: ignore[no-untyped-def]
    import dataclasses

    with pytest.raises(dataclasses.FrozenInstanceError):
        tools["bakta"].version = "0"


def test_bakta_prefetch_job_uses_the_pinned_container() -> None:
    # The prefetch job installs the DB with the image pinned for BAKTA (not a separate conda env).
    from giap.annotate.tools import read_containers_config

    pinned = read_containers_config(REPO / "workflow" / "conf" / "containers.config")["BAKTA"]
    text = (REPO / "scripts" / "slurm" / "fetch_bakta_db.sbatch").read_text()
    assert pinned.startswith("docker://oschwengers/bakta@sha256:")
    assert "withName: 'BAKTA'" in text
    assert "oschwengers-bakta@sha256-$DIGEST.img" in text


# ---------------------------------------------------------------- required_tools (examples 5-7)


def _names(overrides: list[str]) -> set[str]:
    from giap.annotate.tools import required_tools
    from giap.config import resolve_config

    return {t.name for t in required_tools(resolve_config(overrides=overrides))}


def test_required_tools_default_config() -> None:
    assert _names([]) == {
        "bakta",
        "uniop",
        "promotech",
        "promoter_calculator",
        "terminatornet",
        "transtermhp",
        "ostir",
        "genomad",
        "antismash",
        "miniprot",
        "blast",
    }


def test_required_tools_returns_toolspecs() -> None:
    from giap.annotate.tools import TOOLS, ToolSpec, required_tools
    from giap.config import GiapConfig

    out = required_tools(GiapConfig())
    assert isinstance(out, list)
    assert all(isinstance(t, ToolSpec) and TOOLS[t.name] is t for t in out)
    assert len({t.name for t in out}) == len(out)


def test_required_tools_bactermfinder_opt_in() -> None:
    assert "bactermfinder" in _names(["annotation.terminators.bactermfinder=true"])


def test_required_tools_precomputed_bakta_core_only() -> None:
    names = _names(["annotation.precomputed.bakta_json=x.json", "modules=[core]"])
    assert not names & {"bakta", "genomad", "antismash"}
    assert "uniop" in names


def test_required_tools_without_transtermhp_cross_check() -> None:
    names = _names(["annotation.terminators.cross_check_transtermhp=false"])
    assert "transtermhp" not in names
    assert "terminatornet" in names
