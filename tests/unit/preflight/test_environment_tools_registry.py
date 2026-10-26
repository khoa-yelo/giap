"""``environment.tools`` from the tool registry and the ``giap setup`` manifest (run-wiring lead
decision 4 a-c, docs/specs/M2-07-nextflow-tool-processes.md; GIAP_SPEC.md §7.0 Environment rows).

a. Binaries/environments come from the M2-01 registry (``giap.annotate.tools``), not from the DB
   manifest: a tool is present when its pinned image (``workflow/conf/containers.config``) exists
   under ``$GIAP_ROOT/containers`` (apptainer profiles) or, for a conda-only tool, its prebuilt env
   ``$GIAP_ROOT/micromamba/envs/giap-<tool>`` holds ``.giap-built`` equal to the sha256 of
   ``workflow/envs/<tool>.yml``. Under the ``test``/``stub`` profiles no tool is required.
b. Databases and models come from the ``giap setup`` manifest with item names mapped to
   requirement names: ``cog2020_def``+``cog2020_fun`` → ``cog``, ``go_basic``+``goslim_generic`` →
   ``go``, ``kegg_*`` → ``kegg``, ``promotech_models`` → ``promotech`` (model), ``contaminants`` →
   ``contaminants``; every part of a group must be ``ok``.
c. A manifest in the real ``giap setup`` format (tests/data/db_manifest/db_manifest.json, versions
   and paths faked) and fake container/env directories; never the real $GIAP_DB_DIR / $GIAP_ROOT.

Interface pinned here (implementer):
- ``giap.preflight.environment.missing_tools(tools, profile, image_cache, env_root, local_tools)
  -> list[str]`` (14f1; called here with a fake ``$GIAP_ROOT``'s two directories): sorted
  registry names of ``tools`` (``Sequence[ToolSpec]``) not found for the comma-joined ``profile``
  under the fake root; ``[]`` when ``profile`` names ``test`` or ``stub``.
- ``requirements(modules, mode)`` lists databases and models only (kinds ``db``/``model``).
- ``read_manifest(db_dir)`` maps setup item names to requirement names (groups all-``ok``).
- ``run_preflight(..., giap_root=<Path>, profile=<str>)`` keywords; with ``probe=None`` the
  ``environment.tools`` check reports missing tools and databases by name (exit 2).
- ``giap run`` reads ``$GIAP_ROOT`` and passes ``--profile`` to that check.
"""

from __future__ import annotations

import json
import os
import random
import stat
import tempfile
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import (
    MARKER,
    SETUP_MANIFEST,
    conda_only_tools,
    env_yml_sha256,
    fake_giap_root,
    image_file_name,
    pinned_images,
    setup_manifest,
    write_setup_manifest,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.annotate.tools import TOOLS, ToolSpec, required_tools
from giap.config import resolve_config
from giap.exitcodes import ExitCode
from giap.preflight.environment import (
    GIB,
    EnvProbe,
    environment_checks,
    read_manifest,
    requirements,
)

ALL_MODULES = ["core", "mobile", "bgc", "function"]
APPTAINER_PROFILES = ["farmshare,apptainer", "farmshare_local,apptainer", "apptainer"]
NO_TOOL_PROFILES = ["test", "stub", "test,stub"]
GROUPS = {
    "cog": ("cog2020_def", "cog2020_fun"),
    "go": ("go_basic", "goslim_generic"),
    "kegg": ("kegg_ko_pathway", "kegg_pathway_names"),
    "promotech": ("promotech_models",),
    "contaminants": ("contaminants",),
}
DB_REQS = {"bakta_db", "promotech", "genomad_db", "antismash_db", "go", "kegg", "cog"}


def _missing_tools() -> Any:
    from giap.preflight import environment

    fn = getattr(environment, "missing_tools", None)
    assert fn is not None, "giap.preflight.environment.missing_tools is not implemented"

    def under_root(tools: list[ToolSpec], profile: str, root: Path) -> list[str]:
        """Images in ``<root>/containers``, envs in ``<root>/micromamba/envs`` (14f1)."""
        return fn(tools, profile, root / "containers", root / "micromamba" / "envs", None)

    return under_root


def external_tools() -> list[ToolSpec]:
    """Every registry tool run as a separate process (in-process ``codon_bias`` excluded)."""
    return [t for t in TOOLS.values() if not t.in_process]


# ---------------------------------------------------------------- fixture sanity (c)


def test_fixture_has_the_real_setup_manifest_structure() -> None:
    data = setup_manifest()
    assert set(data) == {"schema_version", "items"}
    keys = {"name", "version", "url", "checksum", "size_bytes", "retrieved", "path", "status",
            "message"}  # fmt: skip
    names = [i["name"] for i in data["items"]]  # type: ignore[union-attr]
    assert len(names) == len(set(names))
    for item in data["items"]:  # type: ignore[union-attr]
        assert set(item) == keys, item
        assert item["status"] == "ok"
    for parts in GROUPS.values():
        assert set(parts) <= set(names)
    assert {"bakta_light", "genomad", "antismash"} <= set(names)


def test_fixture_helpers_match_the_repo_configs(tmp_path: Path) -> None:
    pinned = pinned_images()
    assert {"bakta", "miniprot", "blast", "genomad", "antismash", "transtermhp"} <= set(pinned)
    assert {"uniop", "promotech", "ostir"} <= set(conda_only_tools())
    ref = pinned["bakta"]
    assert image_file_name(ref).startswith("oschwengers-bakta@sha256-") and ref.endswith(
        image_file_name(ref).split("@sha256-")[1].removesuffix(".img")
    )
    root = fake_giap_root(tmp_path / "root")
    assert (root / "micromamba" / "envs" / "giap-uniop" / MARKER).read_text().strip() == (
        env_yml_sha256("uniop")
    )


# ---------------------------------------------------------------- a. tools from the registry


def _no_conda(profile: str) -> list[str]:
    """14f1: plain apptainer runs no conda, so conda-only tools are missing there (fgsea is
    started by the function step, not by a conda directive, so its env still counts)."""
    if "conda" in profile or profile.startswith("farmshare"):
        return []
    return sorted(t for t in conda_only_tools() if t != "fgsea")


@pytest.mark.parametrize("profile", APPTAINER_PROFILES)
def test_a_all_tools_present_in_fake_root(tmp_path: Path, profile: str) -> None:
    root = fake_giap_root(tmp_path / "root")
    assert _missing_tools()(external_tools(), profile, root) == _no_conda(profile)


@pytest.mark.parametrize("profile", APPTAINER_PROFILES)
def test_a_missing_image_reported_by_name(tmp_path: Path, profile: str) -> None:
    root = fake_giap_root(tmp_path / "root", skip=["bakta", "miniprot"])
    want = sorted(["bakta", "miniprot", *_no_conda(profile)])
    assert _missing_tools()(external_tools(), profile, root) == want


def test_a_image_must_match_the_pinned_digest(tmp_path: Path) -> None:
    """An image of the same tool under another digest is not the pinned one."""
    root = fake_giap_root(tmp_path / "root", skip=["bakta"])
    ref = pinned_images()["bakta"]
    other = ref.split("@sha256:")[0] + "@sha256:" + "f" * 64
    (root / "containers" / image_file_name(other)).write_bytes(b"other image\n")
    (root / "containers" / "oschwengers-bakta-v1.12.2.img").write_bytes(b"tagged image\n")
    assert _missing_tools()(external_tools(), "farmshare,apptainer", root) == ["bakta"]


def test_a_missing_conda_env_reported_by_name(tmp_path: Path) -> None:
    root = fake_giap_root(tmp_path / "root", skip=["uniop"])
    assert _missing_tools()(external_tools(), "farmshare,apptainer", root) == ["uniop"]


def test_a_env_with_stale_or_no_build_marker_is_missing(tmp_path: Path) -> None:
    root = fake_giap_root(tmp_path / "root", stale=["promotech"], no_marker=["ostir"])
    assert _missing_tools()(external_tools(), "farmshare,apptainer", root) == ["ostir", "promotech"]


def test_a_only_the_given_tools_are_checked(tmp_path: Path) -> None:
    """Tools a run does not need (e.g. opt-in BacTermFinder, antiSMASH with bgc off) never count."""
    root = fake_giap_root(tmp_path / "root", skip=["bactermfinder", "antismash", "genomad"])
    cfg = resolve_config(overrides=["modules=[core]"])
    tools = required_tools(cfg)
    names = {t.name for t in tools}
    assert "bactermfinder" not in names and "antismash" not in names, names
    assert _missing_tools()(tools, "farmshare,apptainer", root) == []


@pytest.mark.parametrize("profile", NO_TOOL_PROFILES)
def test_a_test_and_stub_profiles_require_no_tools(tmp_path: Path, profile: str) -> None:
    empty = tmp_path / "empty_root"
    empty.mkdir()
    assert _missing_tools()(external_tools(), profile, empty) == []


def test_a_empty_root_lists_every_tool(tmp_path: Path) -> None:
    empty = tmp_path / "empty_root"
    empty.mkdir()
    want = sorted(t.name for t in external_tools())
    assert _missing_tools()(external_tools(), "farmshare,apptainer", empty) == want


@settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(st.sets(st.sampled_from(sorted(t.name for t in TOOLS.values() if not t.in_process))))
def test_a_property_missing_is_exactly_the_skipped_tools(skipped: set[str]) -> None:
    with tempfile.TemporaryDirectory() as d:
        root = fake_giap_root(Path(d) / "root", skip=skipped)
        got = _missing_tools()(external_tools(), "farmshare,apptainer", root)
    assert got == sorted(skipped)


# ---------------------------------------------------------------- b. databases from the manifest


def test_b_requirements_are_databases_and_models_only() -> None:
    for mode in ("targeted", "genome"):
        reqs = requirements(ALL_MODULES, mode)
        kinds = {r.kind for r in reqs}
        assert kinds <= {"db", "model"}, f"binaries/envs still in requirements(): {reqs}"
        names = {r.name for r in reqs}
        assert names >= DB_REQS, names
        assert names <= DB_REQS | {"contaminants"}, names
        by_name = {r.name: r for r in reqs}
        assert by_name["promotech"].kind == "model"


def test_b_full_setup_manifest_satisfies_every_database(tmp_path: Path) -> None:
    db = write_setup_manifest(tmp_path / "db")
    found = read_manifest(str(db))
    assert found is not None
    for name in (*GROUPS, "bakta_db", "genomad_db", "antismash_db"):
        assert isinstance(found.get(name), str), (name, found)
    for mode in ("targeted", "genome"):
        reqs = requirements(ALL_MODULES, mode)
        probe = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
        tools = environment_checks(reqs, probe, 0, 0)[0]
        assert tools.level == "ok", tools


@pytest.mark.parametrize(("group", "item"), [(g, i) for g, parts in GROUPS.items() for i in parts])
@pytest.mark.parametrize("how", ["dropped", "absent", "mismatch"])
def test_b_group_with_one_part_not_ok_is_missing(
    tmp_path: Path, group: str, item: str, how: str
) -> None:
    db = tmp_path / "db"
    if how == "dropped":
        write_setup_manifest(db, drop=[item])
    else:
        write_setup_manifest(db, status={item: how})
    found = read_manifest(str(db))
    assert found is not None
    assert found.get(group) is None, (group, item, how, found)
    for other in GROUPS.keys() - {group}:
        assert isinstance(found.get(other), str), (other, found)


@pytest.mark.parametrize(("item", "req"), [("go_basic", "go"), ("cog2020_fun", "cog"),
                                           ("kegg_pathway_names", "kegg"),
                                           ("promotech_models", "promotech")])  # fmt: skip
def test_b_partial_group_fails_environment_tools_by_requirement_name(
    tmp_path: Path, item: str, req: str
) -> None:
    db = write_setup_manifest(tmp_path / "db", drop=[item])
    found = read_manifest(str(db))
    assert found is not None
    reqs = requirements(ALL_MODULES, "targeted")
    probe = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
    tools = environment_checks(reqs, probe, 0, 0)[0]
    assert tools.level == "error" and tools.exit_code == 2
    assert tools.value == {"missing": [req]}, tools.value


def test_b_bakta_light_or_full_satisfies_bakta_db(tmp_path: Path) -> None:
    for drop in (["bakta_full"], ["bakta_light"]):
        found = read_manifest(str(write_setup_manifest(tmp_path / drop[0], drop=drop)))
        assert found is not None and found.get("bakta_db") == "6.1", (drop, found)


def test_b_fixture_is_the_committed_file() -> None:
    assert SETUP_MANIFEST.is_file()
    json.loads(SETUP_MANIFEST.read_text())


# ---------------------------------------------------------------- run_preflight and giap run


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


def _acgt(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _inputs(tmp_path: Path) -> Any:
    from giap.preflight.run import PreflightInputs

    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(parents=True, exist_ok=True)
    g.write_text(f">chr1 [topology=circular]\n{_acgt(120_000)}\n")
    p = tmp_path / "in" / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return PreflightInputs(genome=g, proteins=p)


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Never the real $GIAP_DB_DIR / $GIAP_ROOT, even through a fallback."""
    for var in ("GIAP_DB_DIR", "GIAP_ROOT", "MAMBA_ROOT_PREFIX", "NXF_APPTAINER_CACHEDIR"):
        empty = tmp_path / f"unset_{var.lower()}"
        empty.mkdir()
        monkeypatch.setenv(var, str(empty))
    # 14f1: a set cache directory wins over $GIAP_ROOT/containers; unset, never the session's
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR")


def _run_preflight(tmp_path: Path, root: Path, db: Path, profile: str) -> Any:
    from giap.preflight.run import run_preflight

    low = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
    cfg = resolve_config(overrides=low)
    report, _ = run_preflight(
        _inputs(tmp_path), cfg, None, "run", tmp_path / "OUT", db, giap_root=root, profile=profile
    )
    return [c for c in report.checks if c.rule_id == "environment.tools"]


def test_run_preflight_ok_with_fake_root_and_setup_manifest(tmp_path: Path) -> None:
    root = fake_giap_root(tmp_path / "root")
    db = write_setup_manifest(tmp_path / "db")
    tools = _run_preflight(tmp_path, root, db, "farmshare,apptainer")
    assert [c.level for c in tools] == ["ok"], tools


def test_run_preflight_reports_missing_tools_and_databases_by_name(tmp_path: Path) -> None:
    root = fake_giap_root(tmp_path / "root", skip=["uniop", "bakta"])
    db = write_setup_manifest(tmp_path / "db", drop=["promotech_models"])
    tools = _run_preflight(tmp_path, root, db, "farmshare,apptainer")
    assert len(tools) == 1 and tools[0].level == "error" and tools[0].exit_code == 2, tools
    missing = tools[0].value["missing"]
    assert {"uniop", "bakta", "promotech"} <= set(missing), missing
    for name in ("uniop", "bakta", "promotech"):
        assert name in tools[0].message


@pytest.mark.parametrize("profile", NO_TOOL_PROFILES)
def test_run_preflight_stub_profiles_need_no_tools(tmp_path: Path, profile: str) -> None:
    empty = tmp_path / "empty_root"
    empty.mkdir()
    db = write_setup_manifest(tmp_path / "db")
    tools = _run_preflight(tmp_path, empty, db, profile)
    assert [c.level for c in tools] == ["ok"], tools


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


def _giap_run(tmp_path: Path, *extra: str) -> tuple[int, str]:
    from giap.cli import cli

    inp = _inputs(tmp_path)
    argv = ["run", "--genome", str(inp.genome), "--proteins", str(inp.proteins)]
    argv += ["--outdir", str(tmp_path / "OUT"), *extra]
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "farmshare_local,apptainer"])
def test_giap_run_launches_with_fake_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path, profile: str
) -> None:
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "root")))
    monkeypatch.setenv("GIAP_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    code, text = _giap_run(tmp_path, "--profile", profile)
    assert code == ExitCode.SUCCESS, text
    assert fake_nf.is_file(), "nextflow not launched"


def test_giap_run_missing_env_exits_2_naming_it_without_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_nf: Path
) -> None:
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "root", skip=["ostir"])))
    monkeypatch.setenv("GIAP_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    code, text = _giap_run(tmp_path)
    assert code == ExitCode.ENVIRONMENT, text
    assert "ostir" in text
    assert not fake_nf.exists(), "nextflow launched despite a missing tool"
    report = json.loads((tmp_path / "OUT" / "preflight" / "preflight.json").read_text())
    tools = [c for c in report["checks"] if c["rule_id"] == "environment.tools"]
    assert tools and "ostir" in tools[0]["value"]["missing"], tools
