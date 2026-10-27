"""M11-15c2 review round 1 (PR #441): the B-items the first test set left open.

docs/specs/M11-15-portable-clusters.md item 6 ("New defaults (15c2)", "``--minimal``", the
licence lines, "Record of use"), Q-L25 ("auto with notice, recorded in the manifest"); the
unified ``db_manifest.json`` ``setup`` object of #439 + #441 (15a ``runtime``, 15c2 ``set``).

- B1: ``--images`` with ``--no-images`` is a usage error (exit 2), nothing installed.
- B2: ``setup --check`` prints no licence line and no no-runtime warning and never probes.
- B3: ``licence_notes`` has the KEGG note only when ``function`` ran with an ``ok`` KEGG item.
- B4: ``write_manifest`` merges ``setup`` (``None`` keeps it); ``1.0`` only when ``setup`` is
  empty; ``read_setup`` drops bad keys/values; the schema rejects an unknown key.
- B5: on a minimal DB dir, ``--set modules=...`` and a ``--config`` file's modules win (no
  override, no notice).
- B6 (lead decision): the minimal-run notice is also recorded in ``run_manifest.json``
  ``warnings`` (Q-L25 "recorded in the manifest").
- Nit: ``--minimal`` with an explicit ``--genomad``/``--antismash``: the spec is silent, so the
  reading pinned here is a usage error (exit 2) naming both flags, nothing installed.

No network, no real tool, no real ``$HOME``, XDG directory or DB directory.
"""

from __future__ import annotations

import json
import os
import random
import shutil
import socket
import stat
import subprocess
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
DB_SCHEMA = REPO / "src" / "gippa" / "schemas" / "db_manifest.schema.json"

GENOMAD_LINE = "geNOMAD (module mobile) is for non-commercial use only"
KEGG_LINE = "KEGG is free for academic use only"
NO_RUNTIME = "images: no container runtime found here"
# copied from gippa.cli.run (the notice of spec item 6 / Q-L25)
NOTICE = (
    "modules core,function (DB directory set up with --minimal; mobile and bgc need "
    f"{PKG} setup without --minimal)"
)
ALL = ["core", "mobile", "bgc", "function"]
GENOMAD_NOTE = "geNOMAD: non-commercial"
KEGG_NOTE = "KEGG: academic use"

_Usage = namedtuple("_Usage", "total used free")
HUGE = _Usage(1 << 50, 0, 1 << 50)


# ---------------------------------------------------------------- setup fakes


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, list[Any]]:
    """Every side effect of ``gippa setup`` faked; ``probe_host`` records its calls."""
    import gippa.setupdb as sdb
    from gippa.exitcodes import ExitCode
    from gippa.runner import profiles
    from gippa.setupdb import images as img
    from gippa.setupdb import run as srun
    from gippa.setupdb import uniop_build
    from gippa.setupdb.items import DbRecord

    log: dict[str, list[Any]] = {"installed": [], "images": [], "uniop": [], "probes": []}

    def fake_install(item: Any, _db: Path, _runner: Any, _fetcher: Any) -> DbRecord:
        log["installed"].append(item.name)
        return DbRecord(name=item.name, version=item.version, url=item.url,
                        checksum=item.checksum or "sha256:" + "0" * 64, size_bytes=1,
                        retrieved="2026-10-27", path=item.path, status="ok",
                        message="downloaded")  # fmt: skip

    def fake_images(tools: list[str], runtime: str, *_a: Any, **_k: Any) -> ExitCode:
        log["images"].append((list(tools), runtime))
        return ExitCode.SUCCESS

    def fake_uniop(*a: Any, **k: Any) -> ExitCode:
        log["uniop"].append((a, k))
        return ExitCode.SUCCESS

    pinned = profiles.probe_host  # the autouse FarmShare pin (tests/conftest.py)

    def probe(*a: Any, **k: Any) -> Any:
        log["probes"].append((a, k))
        return pinned(*a, **k)

    def deny(*a: Any, **_k: Any) -> Any:
        raise AssertionError(f"blocked in this test: {a!r}")

    monkeypatch.setattr(srun, "install", fake_install)
    monkeypatch.setattr(img, "install_images", fake_images)
    monkeypatch.setattr(uniop_build, "setup_uniop", fake_uniop)
    monkeypatch.setattr(profiles, "probe_host", probe)
    monkeypatch.setattr(sdb.HttpsFetcher, "fetch", deny)
    monkeypatch.setattr(subprocess, "run", deny)
    monkeypatch.setattr(subprocess, "Popen", deny)
    monkeypatch.setattr(socket, "create_connection", deny)
    monkeypatch.setattr(socket.socket, "connect", deny)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: HUGE)
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR", "SLURM_CLUSTER_NAME"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg_config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg_data"))
    return log


def setup(tmp_path: Path, *args: str) -> Any:
    from gippa.cli import cli

    argv = ["setup", "--db-dir", str(tmp_path / "db"), "--cache-dir", str(tmp_path / "cache"),
            *args]  # fmt: skip
    return CliRunner().invoke(cli, argv, catch_exceptions=False)


def no_runtime_host(monkeypatch: pytest.MonkeyPatch, log: dict[str, list[Any]]) -> None:
    """A host with no runtime that is not FarmShare (overrides the autouse pin)."""
    from gippa.runner import profiles

    probe = profiles.Probe(on_path=frozenset(), docker_ok=False, slurm_job_id=None,
                           cluster_name=None, setup_runtime=None)  # fmt: skip

    def fake(*a: Any, **k: Any) -> Any:
        log["probes"].append((a, k))
        return probe

    monkeypatch.setattr(profiles, "probe_host", fake)


# ---------------------------------------------------------------- B1


@pytest.mark.parametrize("images", [["--images"], ["--images", "all"]])
def test_b1_images_with_no_images_is_a_usage_error(
    tmp_path: Path, calls: dict[str, list[Any]], images: list[str]
) -> None:
    r = setup(tmp_path, "--profile", "docker", *images, "--no-images")
    assert r.exit_code == 2, r.output
    assert "--images" in r.output and "--no-images" in r.output
    assert calls["installed"] == [] and calls["images"] == [] and calls["uniop"] == []
    assert not (tmp_path / "db" / "db_manifest.json").exists()


# ---------------------------------------------------------------- B2


@pytest.mark.parametrize("profile", [[], ["--profile", "auto"], ["--profile", "docker"],
                                     ["--profile", "apptainer"]])  # fmt: skip
def test_b2_check_prints_no_licence_line_and_never_probes(
    tmp_path: Path, calls: dict[str, list[Any]], profile: list[str]
) -> None:
    r = setup(tmp_path, "--check", *profile)
    assert GENOMAD_LINE not in r.output and KEGG_LINE not in r.output, r.output
    assert "non-commercial" not in r.output and "academic use only" not in r.output
    assert NO_RUNTIME not in r.output, r.output
    assert calls["probes"] == []
    assert calls["installed"] == [] and calls["images"] == []


def test_b2_check_on_a_host_without_runtime_prints_no_warning(
    tmp_path: Path, calls: dict[str, list[Any]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The control: the same host without ``--check`` does print the warning."""
    no_runtime_host(monkeypatch, calls)
    r = setup(tmp_path, "--check")
    assert NO_RUNTIME not in r.output, r.output
    assert GENOMAD_LINE not in r.output and KEGG_LINE not in r.output
    assert calls["probes"] == []
    r2 = setup(tmp_path, "--no-genomad", "--no-kegg")
    assert NO_RUNTIME in r2.output, r2.output  # positive control: the warning exists


# ---------------------------------------------------------------- B3 (licence_notes)


def _kegg_items(status: str) -> list[dict[str, Any]]:
    return [{"name": n, "version": "2026-10-15", "url": "https://rest.kegg.jp/x",
             "checksum": "sha256:" + "8" * 64, "size_bytes": 2, "retrieved": "2026-10-15",
             "path": f"function/kegg/{n}.tsv", "status": status, "message": ""}
            for n in ("kegg_ko_pathway", "kegg_pathway_names")]  # fmt: skip


def _other_items() -> list[dict[str, Any]]:
    return [{"name": "go_basic", "version": "2026-10-01", "url": "https://purl.obolibrary.org/x",
             "checksum": "sha256:" + "2" * 64, "size_bytes": 2, "retrieved": "2026-10-01",
             "path": "function/go/go-basic.obo", "status": "ok", "message": ""}]  # fmt: skip


FUNCTION_RAN = ["GIPPA:FUNCTION:FUNCTION_STEP (function)"]
OTHER_RAN = ["GIPPA:CORE:GIPPA_STEP (gene_table)"]


def _notes(modules: list[str], processes: list[str], items: list[dict[str, Any]]) -> list[str]:
    from helpers.manifest_data import build, spec, trace_row

    from gippa.config import GippaConfig

    m = build(cfg=GippaConfig(seed=7, modules=modules),  # type: ignore[arg-type]
              registered=[spec("gene_table"), spec("function", module="function")],
              executor="nextflow", trace=[trace_row(p) for p in processes],
              db_manifest={"schema_version": "1.0", "items": items})  # fmt: skip
    notes = m["licence_notes"]
    assert isinstance(notes, list)
    return notes


def test_b3_positive_control_kegg_installed_and_function_ran() -> None:
    assert KEGG_NOTE in _notes(["core", "function"], FUNCTION_RAN, _kegg_items("ok"))


def test_b3_function_ran_without_kegg_in_the_db() -> None:
    assert KEGG_NOTE not in _notes(["core", "function"], FUNCTION_RAN, _other_items())


@pytest.mark.parametrize("status", ["failed", "absent", "mismatch"])
def test_b3_kegg_items_not_ok_give_no_note(status: str) -> None:
    assert KEGG_NOTE not in _notes(["core", "function"], FUNCTION_RAN, _kegg_items(status))


def test_b3_function_in_modules_but_not_run() -> None:
    assert KEGG_NOTE not in _notes(["core", "function"], OTHER_RAN, _kegg_items("ok"))


@settings(max_examples=40, deadline=None)
@given(function_on=st.booleans(), ran=st.booleans(),
       status=st.sampled_from(["ok", "failed", "absent", "mismatch", None]))  # fmt: skip
def test_b3_kegg_note_iff_function_on_ran_and_kegg_ok(
    function_on: bool, ran: bool, status: str | None
) -> None:
    modules = ["core", "function"] if function_on else ["core"]
    items = _other_items() + (_kegg_items(status) if status else [])
    got = _notes(modules, FUNCTION_RAN if ran else OTHER_RAN, items)
    assert (KEGG_NOTE in got) is (function_on and ran and status == "ok")


# ---------------------------------------------------------------- B4 (the setup object)


def _rec(name: str = "go_basic") -> Any:
    from gippa.setupdb.items import DbRecord

    return DbRecord(name, "1", "https://x", "", 0, "", "", "absent", "")


def _doc(db: Path) -> dict[str, Any]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _validate_db(doc: dict[str, Any]) -> None:
    import jsonschema

    jsonschema.validate(doc, json.loads(DB_SCHEMA.read_text(encoding="utf-8")))


def test_b4_setup_none_keeps_the_existing_setup(tmp_path: Path) -> None:
    from gippa.setupdb.manifest import write_manifest

    write_manifest(tmp_path, [_rec()], setup={"set": "minimal"})
    write_manifest(tmp_path, [_rec("cog2020_def")], setup=None)
    doc = _doc(tmp_path)
    assert doc["setup"] == {"set": "minimal"}
    assert doc["schema_version"] == "1.1"
    write_manifest(tmp_path, [_rec("cog2020_fun")])  # the default is None too
    assert _doc(tmp_path)["setup"] == {"set": "minimal"}
    _validate_db(_doc(tmp_path))


def test_b4_a_new_key_merges_with_the_existing_one(tmp_path: Path) -> None:
    from gippa.setupdb.manifest import read_setup, write_manifest

    write_manifest(tmp_path, [_rec()], setup={"set": "minimal"})
    write_manifest(tmp_path, [], setup={"runtime": "apptainer"})
    doc = _doc(tmp_path)
    assert doc["setup"] == {"set": "minimal", "runtime": "apptainer"}
    assert doc["schema_version"] == "1.1"
    assert [i["name"] for i in doc["items"]] == ["go_basic"]
    assert read_setup(tmp_path) == {"set": "minimal", "runtime": "apptainer"}
    _validate_db(doc)
    write_manifest(tmp_path, [], setup={"set": "default"})  # same key: the new value wins
    assert _doc(tmp_path)["setup"] == {"set": "default", "runtime": "apptainer"}


def test_b4_schema_1_0_only_without_setup(tmp_path: Path) -> None:
    from gippa.setupdb.manifest import write_manifest

    write_manifest(tmp_path, [_rec()])
    doc = _doc(tmp_path)
    assert doc["schema_version"] == "1.0" and "setup" not in doc
    write_manifest(tmp_path, [_rec()], setup={})
    doc = _doc(tmp_path)
    assert doc["schema_version"] == "1.0" and "setup" not in doc
    _validate_db(doc)


SETUP_VALUES = {"runtime": ["apptainer", "singularity", "docker"], "set": ["default", "minimal"]}


@settings(max_examples=40, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(
    writes=st.lists(
        st.dictionaries(
            st.sampled_from(sorted(SETUP_VALUES)),
            st.sampled_from(["apptainer", "singularity", "docker", "default", "minimal"]),
        ),
        max_size=5,
    ).map(lambda ws: [{k: v for k, v in w.items() if v in SETUP_VALUES[k]} for w in ws])
)
def test_b4_setup_is_the_union_of_every_write(
    tmp_path_factory: pytest.TempPathFactory, writes: list[dict[str, str]]
) -> None:
    """Property: after any sequence of writes, ``setup`` is the dict union (later wins), the
    version is 1.1 iff it is non-empty, and the document validates."""
    from gippa.setupdb.manifest import read_setup, write_manifest

    db = tmp_path_factory.mktemp("setup")
    want: dict[str, str] = {}
    for w in writes:
        write_manifest(db, [], setup=w)
        write_manifest(db, [], setup=None)
        want |= w
    if not writes:
        write_manifest(db, [])
    doc = _doc(db)
    assert doc.get("setup", {}) == want
    assert doc["schema_version"] == ("1.1" if want else "1.0")
    assert read_setup(db) == want
    _validate_db(doc)


@pytest.mark.parametrize(
    ("raw", "kept"),
    [
        ({"set": "minimal", "runtime": "slurm"}, {"set": "minimal"}),
        ({"set": "tiny", "runtime": "docker"}, {"runtime": "docker"}),
        ({"set": "minimal", "executor": "slurm"}, {"set": "minimal"}),
        ({"set": 1, "runtime": "apptainer"}, {"runtime": "apptainer"}),
        ({"set": ["minimal"], "runtime": None}, {}),
        ("minimal", {}),
    ],
)
def test_b4_read_setup_drops_bad_keys_and_values(tmp_path: Path, raw: Any, kept: Any) -> None:
    from gippa.setupdb.manifest import read_setup, write_manifest

    (tmp_path / "db_manifest.json").write_text(
        json.dumps({"schema_version": "1.1", "items": [], "setup": raw}), encoding="utf-8"
    )
    assert read_setup(tmp_path) == kept
    write_manifest(tmp_path, [])  # a later write keeps only the valid part
    doc = _doc(tmp_path)
    assert doc.get("setup", {}) == kept
    _validate_db(doc)


def test_b4_read_setup_without_manifest_is_empty(tmp_path: Path) -> None:
    from gippa.setupdb.manifest import read_setup

    assert read_setup(tmp_path) == {}
    (tmp_path / "db_manifest.json").write_text("{not json", encoding="utf-8")
    assert read_setup(tmp_path) == {}


@pytest.mark.parametrize("setup_obj", [{"executor": "slurm"}, {"set": "minimal", "x": "y"}])
def test_b4_schema_rejects_an_unknown_setup_key(setup_obj: dict[str, str]) -> None:
    import jsonschema

    with pytest.raises(jsonschema.ValidationError):
        _validate_db({"schema_version": "1.1", "items": [], "setup": setup_obj})


def test_b4_schema_accepts_both_keys() -> None:
    _validate_db({"schema_version": "1.1", "items": [],
                  "setup": {"set": "minimal", "runtime": "singularity"}})  # fmt: skip


# ---------------------------------------------------------------- B5 (run on a minimal DB dir)


@pytest.fixture
def seen(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[list[str]]:
    """The modules ``_run_preflight`` receives; no nextflow on PATH, so the run stops there."""
    import gippa.cli as cli_mod

    got: list[list[str]] = []

    def fake_preflight(cfg: Any, *_a: Any, **_k: Any) -> Any:
        got.append(list(cfg.modules))
        return {}, cfg

    monkeypatch.setattr(cli_mod, "_run_preflight", fake_preflight)
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "SLURM_JOB_ID", "GIPPA_CONTAINER_RUNTIME"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    return got


def _minimal_db(tmp_path: Path) -> Path:
    db = tmp_path / "db"
    db.mkdir()
    doc = {"schema_version": "1.1", "items": [], "setup": {"set": "minimal"}}
    (db / "db_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    return db


def _run_cli(tmp_path: Path, *args: str) -> Any:
    from gippa.cli import cli

    genome = tmp_path / "g.fa"
    genome.write_text(">c1\nACGTACGTACGT\n", encoding="utf-8")
    argv = ["run", "--genome", str(genome), "--outdir", str(tmp_path / "OUT"),
            "--profile", "docker", *args]  # fmt: skip
    return CliRunner().invoke(cli, argv, catch_exceptions=False)


@pytest.mark.parametrize("mods", [["core", "mobile"], ["core", "bgc"], ["core"]])
def test_b5_set_modules_wins_without_a_notice(
    tmp_path: Path, seen: list[list[str]], mods: list[str]
) -> None:
    db = _minimal_db(tmp_path)
    r = _run_cli(tmp_path, "--db-dir", str(db), "--set", f"modules=[{','.join(mods)}]")
    assert seen == [mods], r.output
    assert "set up with --minimal" not in r.output


def test_b5_config_file_modules_win_without_a_notice(tmp_path: Path, seen: list[list[str]]) -> None:
    db = _minimal_db(tmp_path)
    cfg = tmp_path / "c.yaml"
    cfg.write_text("modules: [core, mobile]\n", encoding="utf-8")
    r = _run_cli(tmp_path, "--db-dir", str(db), "--config", str(cfg))
    assert seen == [["core", "mobile"]], r.output
    assert "set up with --minimal" not in r.output


def test_b5_control_no_modules_anywhere_gets_the_override(
    tmp_path: Path, seen: list[list[str]]
) -> None:
    r = _run_cli(tmp_path, "--db-dir", str(_minimal_db(tmp_path)))
    assert seen == [["core", "function"]], r.output
    assert NOTICE in r.output


# ---------------------------------------------------------------- B6 (notice in the manifest)

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture
def run_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """A fake ``$GIPPA_ROOT``, a set-up DB dir, a fake ``nextflow`` and small inputs (the 15a
    ``env`` fixture); the run goes as far as the manifest."""
    from helpers.fake_env import fake_gippa_root, write_setup_manifest

    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    for var in ("NXF_APPTAINER_CACHEDIR", "GIPPA_CONTAINER_RUNTIME", "GIPPA_DB_DIR",
                "SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    db = write_setup_manifest(tmp_path / "db")
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
    return {"genome": genome, "proteins": proteins, "db": db}


def _mark_minimal(db: Path) -> None:
    doc = _doc(db)
    doc["schema_version"] = "1.1"
    doc["setup"] = {"set": "minimal"}
    (db / "db_manifest.json").write_text(json.dumps(doc, indent=2), encoding="utf-8")


def _full_run(env: dict[str, Path], out: Path, *extra: str) -> tuple[int, str, dict[str, Any]]:
    from gippa.cli import cli

    argv = ["run", "--genome", str(env["genome"]), "--proteins", str(env["proteins"]),
            "--outdir", str(out), "--db-dir", str(env["db"]), *extra]  # fmt: skip
    r = CliRunner().invoke(cli, argv)
    if r.exception is not None and not isinstance(r.exception, SystemExit):
        raise r.exception
    path = out / "run_manifest.json"
    assert path.is_file(), r.output
    return r.exit_code, r.output, json.loads(path.read_text(encoding="utf-8"))


def _minimal_warnings(m: dict[str, Any]) -> list[str]:
    return [w for w in m["warnings"] if "minimal" in w]


def test_b6_minimal_notice_is_recorded_in_the_manifest_warnings(
    run_env: dict[str, Path], tmp_path: Path
) -> None:
    from helpers.manifest_data import validate

    _mark_minimal(run_env["db"])
    code, text, m = _full_run(run_env, tmp_path / "OUT")
    assert code == 0, text
    assert NOTICE in text
    assert m["config"]["modules"] == ["core", "function"]
    assert _minimal_warnings(m), f"no minimal-set entry in warnings: {m['warnings']}"
    validate(m)


def test_b6_no_minimal_warning_on_a_default_db_dir(
    run_env: dict[str, Path], tmp_path: Path
) -> None:
    code, text, m = _full_run(run_env, tmp_path / "OUT")
    assert code == 0, text
    assert NOTICE not in text
    assert _minimal_warnings(m) == [], m["warnings"]


def test_b6_no_minimal_warning_when_modules_are_given(
    run_env: dict[str, Path], tmp_path: Path
) -> None:
    _mark_minimal(run_env["db"])
    code, text, m = _full_run(run_env, tmp_path / "OUT", "--modules", "core,function")
    assert code == 0, text
    assert _minimal_warnings(m) == [], m["warnings"]


# ---------------------------------------------------------------- nit: --minimal + explicit flag


@pytest.mark.parametrize("flag", ["--genomad", "--antismash"])
def test_minimal_with_an_explicit_dropped_item_is_a_usage_error(
    tmp_path: Path, calls: dict[str, list[Any]], flag: str
) -> None:
    """Spec item 6 is silent; ``--minimal`` must not silently drop an explicit flag."""
    r = setup(tmp_path, "--minimal", flag, "--no-images")
    assert r.exit_code == 2, r.output
    assert "--minimal" in r.output and flag in r.output, r.output
    assert calls["installed"] == []


@pytest.mark.parametrize("flag", ["--no-genomad", "--no-antismash", "--no-kegg"])
def test_minimal_with_a_consistent_opt_out_is_fine(
    tmp_path: Path, calls: dict[str, list[Any]], flag: str
) -> None:
    r = setup(tmp_path, "--minimal", flag, "--no-images")
    assert r.exit_code == 0, r.output
    assert "genomad" not in calls["installed"] and "antismash" not in calls["installed"]
