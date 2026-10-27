"""M11-15c2: zero-flag ``gippa setup`` with everything on, ``--minimal``, the licence lines.

docs/specs/M11-15-portable-clusters.md item 6 ("New defaults (15c2)", "``--minimal``", the two
licence lines of the plan, the ``gippa_setup.sbatch`` static test), examples 6.1-6.3 and 6.7,
the acceptance criteria of 15c2; item 1 edge case (c) (no runtime: images skipped with a
warning); Khoa's decision Q-K35 ("turn everything on and just document the license stuff").

Interface pinned here (readings chosen where the spec leaves room):

- ``SetupOptions(db_dir)`` defaults: Bakta light, geNOMAD, antiSMASH, Promotech models, GO, KEGG,
  COG and contaminants on; UniOP off (it keeps its terms prompt); Bakta full stays
  ``--bakta full``.
- Click: every item has a ``--no-*`` form (``--no-antismash``, ``--no-go``, ``--no-kegg``,
  ``--no-cog``, ``--no-contaminants``, ``--no-images``; ``--no-genomad`` and
  ``--no-promotech-models`` exist already) and ``--minimal`` exists.
- Images: on by default for the runtime of ``--profile`` (and, with 15a, the detected one): the
  same tools as a bare ``--images`` today (``image_set(cfg, False)``), pulled through
  ``gippa.setupdb.images.install_images(tools, runtime, cache, db_dir)``; ``--no-images`` pulls
  nothing; with no runtime the images are skipped with the item 1 (c) warning and the exit code
  follows the DBs. With images on, the UniOP step runs (``setup_uniop(..., accepted=<flag>,
  explicit=False)``): consent stays the user's (``--accept-uniop-terms`` or the prompt).
- ``--minimal``: Bakta light, Promotech models, GO, KEGG, COG, contaminants; no geNOMAD, no
  antiSMASH; images for ``modules: [core, function]``; ``db_manifest.json`` gets
  ``{"setup": {"set": "minimal", ...}}`` (a later setup without ``--minimal`` no longer records
  ``minimal``). ``--images all`` with ``--minimal``: every image, minimal DBs.
- The plan prints the geNOMAD line when geNOMAD is planned and the KEGG line when KEGG is
  planned, exactly as below (stdout or stderr; ``CliRunner`` mixes them).

What is downloaded is observed at ``gippa.setupdb.run.install`` (faked: an ``ok`` record, nothing
fetched); images at ``install_images`` and UniOP at ``setup_uniop`` (both faked). No network, no
real tool, no real ``$GIPPA_DB_DIR``, ``$HOME`` or XDG directory.
"""

from __future__ import annotations

import dataclasses
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tomllib
from collections import namedtuple
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]

GENOMAD_LINE = (
    f"geNOMAD (module mobile) is for non-commercial use only; for commercial use: {PKG} setup "
    f"--no-genomad, {PKG} run --modules core,bgc,function"
)
KEGG_LINE = f"KEGG is free for academic use only; for other use: {PKG} setup --no-kegg"
NO_RUNTIME = (
    f"images: no container runtime found here; rerun {PKG} setup --images --profile "
    "apptainer|singularity|docker on a host that has one (same --db-dir)"
)

GROUPS: dict[str, tuple[str, ...]] = {
    "genomad": ("genomad",),
    "antismash": ("antismash",),
    "promotech_models": ("promotech_models",),
    "go": ("go_basic", "goslim_generic"),
    "kegg": ("kegg_ko_pathway", "kegg_pathway_names"),
    "cog": ("cog2020_def", "cog2020_fun"),
    "contaminants": ("contaminants",),
}
EVERYTHING = ["bakta_light", *(n for g in GROUPS.values() for n in g)]
MINIMAL = ["bakta_light", "promotech_models", "go_basic", "goslim_generic", "kegg_ko_pathway",
           "kegg_pathway_names", "cog2020_def", "cog2020_fun", "contaminants"]  # fmt: skip
NO_FLAG = {g: "--no-" + g.replace("_", "-") for g in GROUPS}

_Usage = namedtuple("_Usage", "total used free")
HUGE = _Usage(1 << 50, 0, 1 << 50)


@dataclasses.dataclass
class Calls:
    installed: list[str] = dataclasses.field(default_factory=list)
    images: list[tuple[list[str], str]] = dataclasses.field(default_factory=list)
    uniop: list[dict[str, Any]] = dataclasses.field(default_factory=list)
    opts: list[Any] = dataclasses.field(default_factory=list)
    blocked: list[str] = dataclasses.field(default_factory=list)
    printed_before_first_install: str | None = None


def _printed_so_far() -> str:
    """What ``CliRunner`` has captured on stdout and stderr up to now."""
    text = []
    for stream in (sys.stdout, sys.stderr):
        stream.flush()
        buf = getattr(stream, "buffer", None)
        if hasattr(buf, "getvalue"):
            text.append(buf.getvalue().decode("utf-8", "replace"))
    return "\n".join(text)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Calls:
    """Fakes for every side effect of ``gippa setup``; the host stays out of the test."""
    import gippa.setupdb as sdb
    from gippa.exitcodes import ExitCode
    from gippa.setupdb import images as img
    from gippa.setupdb import run as srun
    from gippa.setupdb import uniop_build
    from gippa.setupdb.items import DbRecord

    c = Calls()

    def fake_install(item: Any, db_dir: Path, _runner: Any, _fetcher: Any) -> DbRecord:
        if not c.installed:
            c.printed_before_first_install = _printed_so_far()
        c.installed.append(item.name)
        return DbRecord(name=item.name, version=item.version, url=item.url,
                        checksum=item.checksum or "sha256:" + "0" * 64, size_bytes=1,
                        retrieved="2026-10-27", path=item.path, status="ok",
                        message="downloaded")  # fmt: skip

    def fake_images(tools: list[str], runtime: str, *_a: Any, **_k: Any) -> ExitCode:
        c.images.append((list(tools), str(runtime)))
        return ExitCode.SUCCESS

    def fake_uniop(profile: str, *_a: Any, accepted: bool, explicit: bool, **_k: Any) -> ExitCode:
        c.uniop.append({"profile": profile, "accepted": accepted, "explicit": explicit})
        return ExitCode.SUCCESS

    real_run_setup = sdb.run_setup

    def recording_run_setup(opts: Any, *a: Any, **k: Any) -> Any:
        c.opts.append(opts)
        return real_run_setup(opts, *a, **k)

    def deny(what: str) -> Any:
        def _deny(*a: Any, **_k: Any) -> Any:
            c.blocked.append(f"{what}: {a!r}")
            raise FileNotFoundError(f"blocked in this test: {what}")

        return _deny

    def deny_net(*a: Any, **_k: Any) -> Any:
        raise AssertionError(f"gippa setup touched the network: {a!r}")

    monkeypatch.setattr(srun, "install", fake_install)
    monkeypatch.setattr(img, "install_images", fake_images)
    monkeypatch.setattr(uniop_build, "setup_uniop", fake_uniop)
    monkeypatch.setattr(sdb, "run_setup", recording_run_setup)
    monkeypatch.setattr(sdb.HttpsFetcher, "fetch", deny_net)
    monkeypatch.setattr(subprocess, "run", deny("subprocess.run"))  # quota, docker info, ...
    monkeypatch.setattr(subprocess, "Popen", deny("subprocess.Popen"))
    monkeypatch.setattr(socket, "create_connection", deny_net)
    monkeypatch.setattr(socket.socket, "connect", deny_net)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: HUGE)  # 15c1's free-space check passes
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR", "SLURM_CLUSTER_NAME"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")  # no "run --bakta full as a batch job" warning
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg_config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg_data"))
    return c


def setup(tmp_path: Path, *args: str, db: Path | None = None) -> Any:
    from gippa.cli import cli

    db = db or tmp_path / "db"
    argv = ["setup", "--db-dir", str(db), "--cache-dir", str(tmp_path / "cache"), *args]
    return CliRunner().invoke(cli, argv, catch_exceptions=False)


def empty_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """No container runtime, no Slurm, no conda on ``PATH`` (item 1 R5)."""
    empty = tmp_path / "empty_bin"
    empty.mkdir(exist_ok=True)
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)


def manifest(db: Path) -> dict[str, Any]:
    import json

    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def config_images(*overrides: str, all_images: bool = False) -> list[str]:
    from gippa.config import resolve_config
    from gippa.setupdb.images import image_set

    return image_set(resolve_config(overrides=list(overrides)), all_images)


# ---------------------------------------------------------------- SetupOptions() (the table)


def test_setup_options_defaults_turn_everything_on(tmp_path: Path) -> None:
    from gippa.setupdb import SetupOptions

    o = SetupOptions(db_dir=tmp_path)
    got = {k: getattr(o, k) for k in ("bakta", "genomad", "antismash", "promotech_models", "go",
                                      "kegg", "cog", "contaminants", "uniop", "check")}  # fmt: skip
    assert got == {"bakta": "light", "genomad": True, "antismash": True,
                   "promotech_models": True, "go": True, "kegg": True, "cog": True,
                   "contaminants": True, "uniop": False, "check": False}  # fmt: skip


def test_default_plan_is_everything_but_bakta_full_and_uniop(tmp_path: Path) -> None:
    from gippa.setupdb import SetupOptions
    from gippa.setupdb.items import plan

    names = [i.name for i in plan(SetupOptions(db_dir=tmp_path))]
    assert names == EVERYTHING
    assert "bakta_full" not in names and "uniop_source" not in names


@settings(max_examples=60, deadline=None)
@given(off=st.sets(st.sampled_from(sorted(GROUPS))))
def test_plan_is_the_defaults_minus_each_opt_out(off: set[str]) -> None:
    from gippa.setupdb import SetupOptions
    from gippa.setupdb.items import plan

    names = [i.name for i in plan(SetupOptions(db_dir=Path("/x"), **{g: False for g in off}))]
    dropped = {n for g in off for n in GROUPS[g]}
    assert names == [n for n in EVERYTHING if n not in dropped]


# ---------------------------------------------------------------- CLI defaults and --no-* forms


def test_help_lists_every_opt_out_and_minimal() -> None:
    from gippa.cli import cli

    text = CliRunner().invoke(cli, ["setup", "--help"]).output
    for flag in (*NO_FLAG.values(), "--no-images", "--minimal", "--accept-uniop-terms"):
        assert re.search(rf"(?<![\w-]){re.escape(flag)}(?![\w-])", text), flag


def test_zero_flag_setup_downloads_everything(tmp_path: Path, calls: Calls) -> None:
    """Example 6.1 (DB part): no item flag at all → every default item, Bakta light only."""
    r = setup(tmp_path, "--no-images")
    assert r.exit_code == 0, r.output
    assert calls.installed == EVERYTHING
    assert {i["name"] for i in manifest(tmp_path / "db")["items"]} >= set(EVERYTHING)


@pytest.mark.parametrize("group", sorted(GROUPS))
def test_each_no_flag_drops_only_its_items(tmp_path: Path, calls: Calls, group: str) -> None:
    r = setup(tmp_path, "--no-images", NO_FLAG[group])
    assert r.exit_code == 0, r.output
    assert calls.installed == [n for n in EVERYTHING if n not in GROUPS[group]]


def test_bakta_full_stays_opt_in(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--no-images", "--bakta", "full")
    assert r.exit_code == 0, r.output
    assert calls.installed == ["bakta_full", *EVERYTHING[1:]]


def test_check_verifies_the_default_set_and_downloads_nothing(tmp_path: Path, calls: Calls) -> None:
    """Example 6.7: ``--check`` alone verifies the new default set; nothing is fetched/pulled."""
    from gippa.setupdb.items import plan

    r = setup(tmp_path, "--check", "--profile", "docker")
    assert r.exit_code == 2, r.output  # nothing is installed in the empty directory
    assert calls.installed == [] and calls.images == [] and calls.uniop == []
    assert len(calls.opts) == 1 and calls.opts[0].check
    assert [i.name for i in plan(calls.opts[0])] == EVERYTHING


# ---------------------------------------------------------------- images for the runtime


@pytest.mark.parametrize("runtime", ["docker", "apptainer"])
def test_images_are_on_by_default_for_the_runtime(
    tmp_path: Path, calls: Calls, runtime: str
) -> None:
    r = setup(tmp_path, "--profile", runtime)
    assert r.exit_code == 0, r.output
    assert calls.images == [(config_images(), runtime)]


def test_no_images_pulls_nothing(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--profile", "docker", "--no-images")
    assert r.exit_code == 0, r.output
    assert calls.images == [] and calls.uniop == []
    assert calls.installed == EVERYTHING


def test_images_all_still_pulls_every_pin(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--profile", "docker", "--images", "all")
    assert r.exit_code == 0, r.output
    assert calls.images == [(config_images(all_images=True), "docker")]


def test_no_runtime_skips_images_with_a_warning_and_exits_by_the_dbs(
    tmp_path: Path, calls: Calls, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Item 1 edge case (c) / example 1.14: DBs downloaded, images skipped, exit 0."""
    empty_path(monkeypatch, tmp_path)
    r = setup(tmp_path)
    assert NO_RUNTIME in r.output
    assert calls.images == [] and calls.uniop == []
    assert calls.installed == EVERYTHING
    assert r.exit_code == 0, r.output
    assert "runtime" not in manifest(tmp_path / "db").get("setup", {})


def test_bactermfinder_stays_opt_in(tmp_path: Path, calls: Calls) -> None:
    """Q-K43 pending: the default image set has no BacTermFinder (``--add-bactermfinder``)."""
    r = setup(tmp_path, "--profile", "docker")
    assert r.exit_code == 0, r.output
    assert calls.images and all("bactermfinder" not in tools for tools, _ in calls.images)


# ---------------------------------------------------------------- UniOP keeps its consent


def test_default_images_run_the_uniop_step_without_consent(tmp_path: Path, calls: Calls) -> None:
    """Example 6.1: "asks for UniOP consent": the step runs, consent is not given for the user."""
    r = setup(tmp_path, "--profile", "docker")
    assert r.exit_code == 0, r.output
    assert calls.uniop == [{"profile": "docker", "accepted": False, "explicit": False}]
    assert "uniop_source" not in calls.installed


def test_accept_uniop_terms_is_passed_through(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--profile", "docker", "--accept-uniop-terms")
    assert r.exit_code == 0, r.output
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": False}]


def test_no_uniop_skips_the_uniop_step(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--profile", "docker", "--no-uniop")
    assert r.exit_code == 0, r.output
    assert calls.images and calls.uniop == []


# ---------------------------------------------------------------- licence lines (exact)


def test_zero_flag_plan_prints_both_licence_lines(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--no-images")
    assert r.exit_code == 0, r.output
    assert GENOMAD_LINE in r.output and KEGG_LINE in r.output
    assert r.output.count(GENOMAD_LINE) == 1 and r.output.count(KEGG_LINE) == 1


def test_licence_lines_come_before_any_download(tmp_path: Path, calls: Calls) -> None:
    """The lines are part of the plan, printed before the first item is installed."""
    r = setup(tmp_path, "--no-images")
    assert r.exit_code == 0, r.output
    assert calls.printed_before_first_install is not None
    before = calls.printed_before_first_install
    assert GENOMAD_LINE in before and KEGG_LINE in before, before


def test_no_genomad_no_kegg_prints_neither_line(tmp_path: Path, calls: Calls) -> None:
    """Example 6.2."""
    r = setup(tmp_path, "--no-images", "--no-genomad", "--no-kegg")
    assert r.exit_code == 0, r.output
    assert "genomad" not in calls.installed
    assert not {"kegg_ko_pathway", "kegg_pathway_names"} & set(calls.installed)
    assert "non-commercial" not in r.output and "academic use only" not in r.output


@settings(max_examples=12, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(no_genomad=st.booleans(), no_kegg=st.booleans())
def test_each_licence_line_is_printed_iff_its_item_is_planned(
    tmp_path_factory: pytest.TempPathFactory, calls: Calls, no_genomad: bool, no_kegg: bool
) -> None:
    tmp = tmp_path_factory.mktemp("lic")
    flags = ["--no-images"] + ["--no-genomad"] * no_genomad + ["--no-kegg"] * no_kegg
    r = setup(tmp, *flags)
    assert r.exit_code == 0, r.output
    assert (GENOMAD_LINE in r.output) is not no_genomad
    assert (KEGG_LINE in r.output) is not no_kegg


# ---------------------------------------------------------------- --minimal (example 6.3)


def test_minimal_downloads_the_small_set(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--minimal", "--no-images")
    assert r.exit_code == 0, r.output
    assert calls.installed == MINIMAL


def test_minimal_is_recorded_in_the_db_manifest(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--minimal", "--no-images")
    assert r.exit_code == 0, r.output
    assert manifest(tmp_path / "db")["setup"]["set"] == "minimal"


def test_full_setup_after_minimal_clears_the_minimal_record(tmp_path: Path, calls: Calls) -> None:
    """The run notice says "mobile and bgc need gippa setup without --minimal": doing so must
    lift the restriction."""
    assert setup(tmp_path, "--minimal", "--no-images").exit_code == 0
    r = setup(tmp_path, "--no-images")
    assert r.exit_code == 0, r.output
    assert manifest(tmp_path / "db").get("setup", {}).get("set") != "minimal"


def test_minimal_prints_only_the_kegg_line(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--minimal", "--no-images")
    assert r.exit_code == 0, r.output
    assert KEGG_LINE in r.output and GENOMAD_LINE not in r.output


def test_minimal_pulls_the_core_function_images(tmp_path: Path, calls: Calls) -> None:
    r = setup(tmp_path, "--minimal", "--profile", "docker")
    assert r.exit_code == 0, r.output
    want = config_images("modules=[core,function]")
    assert want != config_images()  # the minimal set is smaller (no geNOMAD, no antiSMASH)
    assert calls.images == [(want, "docker")]
    assert calls.installed == MINIMAL


def test_minimal_with_images_all_pulls_every_image(tmp_path: Path, calls: Calls) -> None:
    """Edge case: ``--images all`` with ``--minimal`` is allowed (images all, DBs minimal)."""
    r = setup(tmp_path, "--minimal", "--profile", "docker", "--images", "all")
    assert r.exit_code == 0, r.output
    assert calls.images == [(config_images(all_images=True), "docker")]
    assert calls.installed == MINIMAL


# ---------------------------------------------------------------- gippa_setup.sbatch (static)


def _sbatch_default_args() -> list[str]:
    text = (REPO / "scripts/slurm/gippa_setup.sbatch").read_text(encoding="utf-8")
    m = re.search(r"if \[ \$# -eq 0 \]; then set -- (.*?); fi", text)
    assert m, "gippa_setup.sbatch: default argument line not found"
    return shlex.split(m.group(1))


def test_gippa_setup_sbatch_installs_the_same_set_as_before(tmp_path: Path, calls: Calls) -> None:
    """Item 6: the sbatch keeps its explicit flags; the resulting DB set on FarmShare is the one
    it installed before 15c2 (Bakta full, geNOMAD, Promotech, GO, KEGG, COG, contaminants)."""
    r = setup(tmp_path, *_sbatch_default_args())
    assert r.exit_code == 0, r.output
    assert calls.installed == ["bakta_full", "genomad", "promotech_models", "go_basic",
                               "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names",
                               "cog2020_def", "cog2020_fun", "contaminants"]  # fmt: skip
