"""M11-14f2: ``setup --images`` and the docker path of ``ToolRunner``.

docs/specs/M11-14-container-install.md item 3 (the ``setup --images`` and ``ToolRunner`` bullets,
"Database pins per release"), AC 6, edge case "offline compute nodes"; ADR-0059 (file names);
ADR-0060. 14f2 stacks on 14a (``image_cache``, ``images.json``) and 14f1 (``missing_tools``).

Interface pinned here (readings chosen where the spec leaves room):

- New module ``gippa.setupdb.images``. It reads the pins from ``workflow_main().parent / "conf"``
  (``containers.config`` ``withName`` lines, process name lower-cased to the ``TOOLS`` name, and
  ``images.json`` entries, which add the process-less fgsea). It looks ``workflow_main`` up as the
  module attribute ``gippa.setupdb.images.workflow_main`` (the tests patch it with a fake tree).
- ``image_set(cfg, all_images: bool) -> list[str]``: sorted tool names to pull. ``False``: the
  ``required_tools(cfg, abundance=True)`` names that have a pin; ``True``: every pinned tool.
- ``gippa setup`` gains ``--images`` (bare: the set of the existing ``--config``/``--set``;
  ``--images all``: every pin), ``--profile apptainer|docker`` (required with ``--images``) and
  ``--cache-dir DIR`` (default ``image_cache(os.environ, db_dir)`` from 14a).
  - apptainer: ``apptainer pull <cache>/<simpleName>.img docker://<ref>`` (the last two
    arguments; flags may precede them). An image already in the cache is not pulled again.
  - docker: ``docker pull <ref>`` (the last argument); the pulled digest is the ``Digest: <d>``
    line of its output (what the docker CLI prints). ``<d>`` != the pin: item ``mismatch`` with
    ``pulled <d>, pinned <pin>``; setup exits 2.
  - A failed pull: item ``absent``, other images still pulled, exit 2.
  - ``db_manifest.json`` gets ``image:<tool>`` items with ``version`` = the ``images.json`` tag
    (GHCR) or the upstream tag the digest was resolved from (the ``//   <image>:<tag>`` comment
    lines of ``containers.config``), ``checksum`` = the pinned digest ``sha256:<hex>``,
    ``status`` ``ok``. Other manifest items are kept.
- ``ToolRunner(db_dir=...)`` (keyword; ``gippa setup`` passes its ``--db-dir``). ``command(tool,
  argv)`` looks up ``$MAMBA_ROOT_PREFIX``/``$GIPPA_ROOT/micromamba`` ``envs/gippa-<tool>`` → the
  ``.img`` in ``image_cache(os.environ, db_dir)`` → ``docker run -u <uid>:<gid> -v <p>:<p> …
  <bare ref> <argv>`` (every existing absolute path in ``argv`` bound at the same path). Neither
  ``$GIPPA_ROOT`` nor a checkout is required; nothing found (no docker on ``PATH``) raises
  ``ToolUnavailable``.
- ``setup --check``: an installed DB whose version differs from this release's pin is reported
  ``installed <v_old>, this release pins <v_new>`` (exit 2).

No test runs real docker/apptainer or touches a registry or the network: logging fakes sit first
on ``PATH``; ``$GIPPA_ROOT``, ``$GIPPA_DB_DIR`` and the caches are unset or temporary.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.annotate.tools import TOOLS, read_containers_config, required_tools
from gippa.config import resolve_config

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
CONTAINERS_CONFIG = WORKFLOW / "conf" / "containers.config"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
OWNER = "ghcr.io/khoa-yelo/"
OSTIR_REF = f"{OWNER}{PKG}-ostir@sha256:" + "a1" * 32
BTF_REF = f"{OWNER}{PKG}-bactermfinder@sha256:" + "b2" * 32
FGSEA_REF = f"{OWNER}{PKG}-fgsea@sha256:" + "c3" * 32
GHCR = {
    "ostir": (OSTIR_REF, "1.1.3-0123456789ab", ["OSTIR", "OSTIR_NOVEL"]),
    "bactermfinder": (BTF_REF, "git-deadbeef-0123456789ac", ["BACTERMFINDER"]),
    "fgsea": (FGSEA_REF, "1.36.2-0123456789ad", []),
}
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models", "--no-uniop"]
BAD = "sha256:" + "ee" * 32

FAKE_RUNTIME = r"""#!/bin/sh
printf '%s\n' "$(basename "$0") $*" >> "$FAKE_RUNTIME_LOG"
prog=$(basename "$0")
last=""; prev=""; out=""; grab=0
for a in "$@"; do
    prev=$last; last=$a
    if [ "$grab" = 1 ]; then out=$a; grab=0; fi
    [ "$a" = "--output" ] && grab=1
done
case "$prog $1" in
"docker pull")
    if [ -n "$FAKE_FAIL_FOR" ]; then
        case "$last" in *"$FAKE_FAIL_FOR"*)
            echo "Error response from daemon: manifest unknown" >&2; exit 1;;
        esac
    fi
    d=${last#*@}
    if [ -n "$FAKE_BAD_FOR" ]; then
        case "$last" in *"$FAKE_BAD_FOR"*) d=$FAKE_BAD_DIGEST;; esac
    fi
    echo "${last%@*}: Pulling from fake"
    echo "Digest: $d"
    echo "Status: Downloaded newer image for $last"
    exit 0;;
"apptainer pull")
    if [ -n "$FAKE_FAIL_FOR" ]; then
        case "$last" in *"$FAKE_FAIL_FOR"*)
            echo "FATAL: manifest unknown" >&2; exit 1;;
        esac
    fi
    printf 'fake sif %s\n' "$last" > "$prev"
    exit 0;;
"docker run")
    if [ -n "$out" ]; then
        mkdir -p "$out"
        printf '{"major": 6, "minor": 0, "type": "light"}\n' > "$out/version.json"
    fi
    exit 0;;
esac
echo "fake $prog: unexpected call: $*" >&2
exit 3
"""


# ---------------------------------------------------------------- helpers and fixtures


def _images() -> Any:
    try:
        from gippa.setupdb import images
    except ImportError as exc:  # pragma: no cover - until 14f2 exists
        pytest.fail(f"gippa.setupdb.images is not implemented ({exc})")
    return images


def simple_name(ref: str) -> str:
    """Nextflow's apptainer cache file name for ``ref`` (ADR-0059)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def upstream_tags(config: Path = CONTAINERS_CONFIG) -> dict[str, str]:
    """``<image>`` → ``<tag>`` from the ``//   <image>:<tag>`` comment lines."""
    out: dict[str, str] = {}
    for line in config.read_text(encoding="utf-8").splitlines():
        m = re.fullmatch(r"//\s+([a-z0-9./_-]+):([A-Za-z0-9_.-]+)\s*", line.strip())
        if m:
            out[m.group(1)] = m.group(2)
    return out


def tool_pins(config: Path) -> dict[str, str]:
    """Registry tool → ref from the ``withName`` lines (BLAST/TBLASTN share one image)."""
    return {p.lower(): r for p, r in read_containers_config(config).items() if p.lower() in TOOLS}


def install_fakes(bindir: Path, names: tuple[str, ...] = ("docker", "apptainer")) -> None:
    bindir.mkdir(parents=True, exist_ok=True)
    for exe in names:
        path = bindir / exe
        path.write_text(FAKE_RUNTIME)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Never the real $GIPPA_ROOT / $GIPPA_DB_DIR / caches / envs; fake runtimes first on PATH."""
    for var in (
        "GIPPA_DB_DIR",
        "GIPPA_ROOT",
        "GIPPA_REPO",
        "NXF_APPTAINER_CACHEDIR",
        "MAMBA_ROOT_PREFIX",
        "GIPPA_CONTAINER_RUNTIME",
        "FAKE_FAIL_FOR",
        "FAKE_BAD_FOR",
        "FAKE_BAD_DIGEST",
    ):
        monkeypatch.delenv(var, raising=False)
    # an empty micromamba root: a checkout's own local/micromamba envs never leak in
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no_mamba"))
    install_fakes(tmp_path / "fakebin")
    monkeypatch.setenv("PATH", f"{tmp_path / 'fakebin'}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_RUNTIME_LOG", str(tmp_path / "runtime.log"))


def calls(tmp_path: Path, prefix: str = "") -> list[list[str]]:
    log = tmp_path / "runtime.log"
    lines = log.read_text().splitlines() if log.exists() else []
    return [ln.split() for ln in lines if ln.startswith(prefix)]


def make_tree(root: Path) -> Path:
    """A workflow tree: the repo's containers.config plus GHCR pins, and an images.json with the
    GHCR tools (fgsea process-less). Returns ``main.nf``."""
    (root / "conf").mkdir(parents=True)
    (root / "main.nf").write_text("// fake\n")
    (root / "envs").symlink_to(WORKFLOW / "envs")
    text = CONTAINERS_CONFIG.read_text(encoding="utf-8").rstrip()
    assert text.endswith("}")
    extra = "".join(
        f"  withName: '{p}' {{ container = '{ref}' }}\n"
        for ref, _tag, procs in GHCR.values()
        for p in procs
    )
    (root / "conf" / "containers.config").write_text(text[:-1] + extra + "}\n", encoding="utf-8")
    entries = {
        tool: {
            "key": tag_[-12:],
            "tag": tag_,
            "ref": ref,
            "processes": sorted(procs),
            "built_from": "0" * 40,
            "licence": TOOLS[tool].licence,
        }
        for tool, (ref, tag_, procs) in GHCR.items()
    }
    (root / "conf" / "images.json").write_text(json.dumps(entries, indent=2) + "\n")
    return root / "main.nf"


@pytest.fixture
def tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, str]:
    """Patch the fake tree into ``gippa.setupdb.images``; return tool → pinned ref."""
    main = make_tree(tmp_path / "wf")
    monkeypatch.setattr(_images(), "workflow_main", lambda: main)
    pinned = tool_pins(main.parent / "conf" / "containers.config")
    pinned["fgsea"] = FGSEA_REF
    return pinned


def expected_version(tool: str, ref: str) -> str:
    if tool in GHCR:
        return GHCR[tool][1]
    return upstream_tags()[ref.split("@", 1)[0]]


def setup_cli(args: list[str]) -> Any:
    from gippa.cli import cli

    return CliRunner().invoke(cli, ["setup", *args], catch_exceptions=False)


def manifest_items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {item["name"]: item for item in data["items"]}


def default_set(pinned: dict[str, str], sets: list[str] | None = None) -> list[str]:
    cfg = resolve_config(None, sets or [])
    return sorted(t.name for t in required_tools(cfg, abundance=True) if t.name in pinned)


# ---------------------------------------------------------------- image_set


def test_image_set_default_config_is_required_tools_with_a_pin(tree: dict[str, str]) -> None:
    got = _images().image_set(resolve_config(None, []), False)
    assert got == default_set(tree)
    assert "fgsea" in got  # default modules contain `function` (process-less images.json pin)
    assert "ostir" in got and "bakta" in got
    assert "bactermfinder" not in got  # opt-in
    assert "uniop" not in got  # never pinned (14u1 builds it)


def test_image_set_all_is_every_pin(tree: dict[str, str]) -> None:
    got = _images().image_set(resolve_config(None, ['modules=["core"]']), True)
    assert got == sorted(tree)
    assert {"bactermfinder", "fgsea", "antismash", "genomad"} <= set(got)


def test_image_set_follows_modules(tree: dict[str, str]) -> None:
    got = _images().image_set(resolve_config(None, ['modules=["core"]']), False)
    assert not {"genomad", "antismash", "fgsea"} & set(got)
    assert {"bakta", "miniprot", "blast"} <= set(got)


def test_image_set_precomputed_layer_needs_no_image(tree: dict[str, str]) -> None:
    cfg = resolve_config(None, ['annotation.precomputed.bakta_json="/x/bakta.json"'])
    assert "bakta" not in _images().image_set(cfg, False)


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    mods=st.sets(st.sampled_from(["mobile", "bgc", "function"])),
    btf=st.booleans(),
    all_images=st.booleans(),
)
def test_image_set_property(
    tree: dict[str, str], mods: set[str], btf: bool, all_images: bool
) -> None:
    sets = [
        f"modules={json.dumps(sorted(mods))}",
        f"annotation.terminators.bactermfinder={str(btf).lower()}",
    ]
    cfg = resolve_config(None, sets)
    got = _images().image_set(cfg, all_images)
    assert got == sorted(got) and len(set(got)) == len(got)
    assert set(got) <= set(tree)
    if all_images:
        assert got == sorted(tree)
    else:
        assert got == default_set(tree, sets)
        assert ("fgsea" in got) == ("function" in mods)
        assert ("bactermfinder" in got) == btf


# ---------------------------------------------------------------- setup --images, apptainer


def test_apptainer_pulls_simple_names_into_the_db_dir_cache(
    tmp_path: Path, tree: dict[str, str]
) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"])
    assert res.exit_code == 0, res.output
    want = default_set(tree)
    pulls = calls(tmp_path, "apptainer pull")
    assert len(pulls) == len({tree[t] for t in want})  # one pull per distinct image
    cache = db / "containers"  # no NXF_APPTAINER_CACHEDIR, no GIPPA_ROOT (Q-L11)
    for tool in want:
        ref = tree[tool]
        dest = cache / simple_name(ref)
        assert dest.is_file(), f"{tool}: {dest} not pulled"
        assert [str(dest), f"docker://{ref}"] in [p[-2:] for p in pulls]
    assert not calls(tmp_path, "docker")
    assert {p.name for p in cache.iterdir()} == {simple_name(tree[t]) for t in want}


def test_apptainer_manifest_items(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"])
    assert res.exit_code == 0, res.output
    items = manifest_items(db)
    want = default_set(tree)
    assert {n for n in items if n.startswith("image:")} == {f"image:{t}" for t in want}
    for tool in want:
        item = items[f"image:{tool}"]
        assert item["status"] == "ok", item
        assert item["checksum"] == tree[tool].split("@", 1)[1]
        assert item["checksum"].startswith("sha256:")
        assert item["version"] == expected_version(tool, tree[tool]), tool


def test_upstream_versions_are_the_resolved_tags(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"])
    items = manifest_items(db)
    assert items["image:bakta"]["version"] == upstream_tags()[tree["bakta"].split("@")[0]]
    assert items["image:bakta"]["version"].startswith("v")  # oschwengers/bakta:v<x>, not TOOLS
    assert items["image:ostir"]["version"] == GHCR["ostir"][1]
    assert items["image:fgsea"]["version"] == GHCR["fgsea"][1]


@pytest.mark.parametrize("how", ["flag", "env"])
def test_apptainer_cache_dir_choice(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch, how: str
) -> None:
    db, cache = tmp_path / "db", tmp_path / "cache"
    extra: list[str] = []
    if how == "flag":
        extra = ["--cache-dir", str(cache)]
        monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(tmp_path / "not_this"))
    else:
        monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(cache))
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer", *extra])
    assert res.exit_code == 0, res.output
    assert (cache / simple_name(tree["bakta"])).is_file()
    assert not (db / "containers").exists()
    assert not (tmp_path / "not_this").exists()


def test_gippa_root_cache(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "root"
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    res = setup_cli(
        ["--db-dir", str(tmp_path / "db"), *SKIP_DBS, "--images", "--profile", "apptainer"]
    )
    assert res.exit_code == 0, res.output
    assert (root / "containers" / simple_name(tree["miniprot"])).is_file()


def test_apptainer_rerun_pulls_nothing(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    args = ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"]
    assert setup_cli(args).exit_code == 0
    first = len(calls(tmp_path, "apptainer pull"))
    res = setup_cli(args)
    assert res.exit_code == 0, res.output
    assert len(calls(tmp_path, "apptainer pull")) == first
    assert all(i["status"] == "ok" for n, i in manifest_items(db).items() if n.startswith("image:"))


def test_apptainer_failed_pull_is_absent_and_others_install(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_FAIL_FOR", "genomad")
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"])
    assert res.exit_code == 2, res.output
    items = manifest_items(db)
    assert items["image:genomad"]["status"] == "absent"
    assert items["image:genomad"]["message"]
    assert not (db / "containers" / simple_name(tree["genomad"])).exists()
    assert items["image:bakta"]["status"] == "ok"
    assert (db / "containers" / simple_name(tree["bakta"])).is_file()


def test_apptainer_then_preflight_finds_every_image(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The files setup writes are the ones 14f1's preflight (and Nextflow) look for."""
    from gippa.io.image_key import image_cache
    from gippa.preflight import environment

    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "all", "--profile", "apptainer"])
    assert res.exit_code == 0, res.output
    monkeypatch.setattr(environment, "workflow_main", _images().workflow_main)
    tools = [TOOLS[t] for t in tree if t != "fgsea"]
    cache = image_cache(os.environ, db)
    assert environment.missing_tools(tools, "apptainer", cache, None, None) == []


# ---------------------------------------------------------------- setup --images, docker


def test_docker_pulls_every_image_by_digest(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
    assert res.exit_code == 0, res.output
    want = default_set(tree)
    pulled = {p[-1] for p in calls(tmp_path, "docker pull")}
    assert pulled == {tree[t] for t in want}
    assert all("@sha256:" in ref and not ref.startswith("docker://") for ref in pulled)
    assert not calls(tmp_path, "apptainer")
    assert not (db / "containers").exists()  # docker keeps its own store
    items = manifest_items(db)
    for tool in want:
        assert items[f"image:{tool}"]["status"] == "ok"
        assert items[f"image:{tool}"]["checksum"] == tree[tool].split("@", 1)[1]
        assert items[f"image:{tool}"]["version"] == expected_version(tool, tree[tool])


def test_docker_digest_mismatch(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_BAD_FOR", "ostir")
    monkeypatch.setenv("FAKE_BAD_DIGEST", BAD)
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
    assert res.exit_code == 2, res.output
    pin = OSTIR_REF.split("@", 1)[1]
    items = manifest_items(db)
    assert items["image:ostir"]["status"] == "mismatch"
    assert f"pulled {BAD}, pinned {pin}" in items["image:ostir"]["message"]
    assert f"pulled {BAD}, pinned {pin}" in res.output
    assert items["image:bakta"]["status"] == "ok"


def test_docker_failed_pull(
    tmp_path: Path, tree: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_FAIL_FOR", "miniprot")
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
    assert res.exit_code == 2, res.output
    items = manifest_items(db)
    assert items["image:miniprot"]["status"] == "absent"
    assert items["image:blast"]["status"] == "ok"


def test_set_and_all_change_the_pulled_set(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    res = setup_cli(
        [
            "--db-dir",
            str(db),
            *SKIP_DBS,
            "--set",
            'modules=["core"]',
            "--images",
            "--profile",
            "docker",
        ]
    )
    assert res.exit_code == 0, res.output
    core = {p[-1] for p in calls(tmp_path, "docker pull")}
    assert tree["antismash"] not in core and tree["genomad"] not in core and FGSEA_REF not in core
    assert tree["bakta"] in core
    (tmp_path / "runtime.log").unlink()
    res = setup_cli(
        [
            "--db-dir",
            str(db),
            *SKIP_DBS,
            "--set",
            'modules=["core"]',
            "--images",
            "all",
            "--profile",
            "docker",
        ]
    )
    assert res.exit_code == 0, res.output
    assert {p[-1] for p in calls(tmp_path, "docker pull")} == set(tree.values())


def test_config_file_drives_the_set(tmp_path: Path, tree: dict[str, str]) -> None:
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("modules: [core, bgc]\n")
    res = setup_cli(
        [
            "--db-dir",
            str(tmp_path / "db"),
            *SKIP_DBS,
            "--config",
            str(cfg),
            "--images",
            "--profile",
            "docker",
        ]
    )
    assert res.exit_code == 0, res.output
    pulled = {p[-1] for p in calls(tmp_path, "docker pull")}
    assert tree["antismash"] in pulled and tree["genomad"] not in pulled


def test_existing_manifest_items_are_kept(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    db.mkdir()
    old = {"name": "genomad", "version": "1.9", "url": "https://zenodo.org/x", "checksum": "md5:0",
           "size_bytes": 1, "retrieved": "2026-01-01", "path": "genomad_db", "status": "ok",
           "message": ""}  # fmt: skip
    (db / "db_manifest.json").write_text(json.dumps({"schema_version": "1.0", "items": [old]}))
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
    assert res.exit_code == 0, res.output
    items = manifest_items(db)
    assert items["genomad"] == old
    assert "image:bakta" in items


def test_without_images_nothing_is_pulled(tmp_path: Path, tree: dict[str, str]) -> None:
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS])
    assert res.exit_code == 0, res.output
    assert not calls(tmp_path)


@pytest.mark.parametrize("profile", ["conda", "standard"])
def test_images_needs_a_container_profile(
    tmp_path: Path, tree: dict[str, str], profile: str
) -> None:
    args = ["--db-dir", str(tmp_path / "db"), *SKIP_DBS, "--images"]
    res = setup_cli([*args, "--profile", profile])
    assert res.exit_code != 0
    assert not calls(tmp_path)


def test_images_without_profile_use_the_detected_runtime(
    tmp_path: Path, tree: dict[str, str]
) -> None:
    """M11-15a: no ``--profile`` is ``auto``; the conftest's pinned FarmShare probe resolves to
    apptainer, so the images are pulled with apptainer (spec M11-15 item 1, D1)."""
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images"])
    assert res.exit_code == 0, res.output
    pulls = calls(tmp_path, "apptainer pull")
    assert len(pulls) == len({tree[t] for t in default_set(tree)})
    assert not calls(tmp_path, "docker")


# ---------------------------------------------------------------- ToolRunner


def _runner(db: Path) -> Any:
    from gippa.setupdb.run import ToolRunner

    return ToolRunner(db_dir=db)


def bakta_ref() -> str:
    return read_containers_config(CONTAINERS_CONFIG)["BAKTA"]


def check_docker_argv(cmd: list[str], ref: str, argv: list[str], bound: list[Path]) -> None:
    assert cmd[:2] == ["docker", "run"], cmd
    assert cmd[-len(argv) - 1 :] == [ref, *argv]
    head = cmd[: -len(argv) - 1]
    assert "-u" in head and head[head.index("-u") + 1] == f"{os.getuid()}:{os.getgid()}"
    vols = [head[i + 1] for i, a in enumerate(head) if a in ("-v", "--volume")]
    for p in bound:
        assert f"{p}:{p}" in vols or any(v.startswith(f"{p}:{p}:") for v in vols), (p, vols)


def test_toolrunner_docker_without_gippa_root(tmp_path: Path) -> None:
    db, out = tmp_path / "db", tmp_path / "db" / "bakta" / "db-light.partial"
    out.mkdir(parents=True)
    argv = ["bakta_db", "download", "--output", str(out), "--type", "light"]
    cmd = _runner(db).command("bakta", argv)
    check_docker_argv(cmd, bakta_ref(), argv, [out])


def test_toolrunner_cached_image_wins_over_docker(tmp_path: Path) -> None:
    db = tmp_path / "db"
    img = db / "containers" / simple_name(bakta_ref())
    img.parent.mkdir(parents=True)
    img.write_bytes(b"sif")
    cmd = _runner(db).command("bakta", ["bakta_db", "--help"])
    assert cmd[0] == "apptainer" and str(img) in cmd and cmd[-2:] == ["bakta_db", "--help"]


def test_toolrunner_env_wins(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    mamba = tmp_path / "mamba"
    (mamba / "envs" / "gippa-bakta").mkdir(parents=True)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(mamba))
    cmd = _runner(tmp_path / "db").command("bakta", ["bakta_db", "--help"])
    assert cmd[:4] == ["micromamba", "run", "-n", "gippa-bakta"]


def test_toolrunner_nothing_available(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from gippa.setupdb.run import ToolUnavailable

    empty = tmp_path / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    with pytest.raises(ToolUnavailable) as exc:
        _runner(tmp_path / "db").command("bakta", ["bakta_db", "--help"])
    assert "GIPPA_ROOT is not set" not in str(exc.value)


def test_setup_downloads_bakta_through_docker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``gippa setup`` needs no checkout env and no micromamba: the downloader runs in docker."""
    from collections import namedtuple

    from gippa.setupdb import run

    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(run.shutil, "disk_usage", lambda _p: usage(1 << 50, 0, 1 << 50))
    db = tmp_path / "db"
    res = setup_cli(
        ["--db-dir", str(db), "--bakta", "light", "--no-genomad", "--no-promotech-models"]
    )
    runs = calls(tmp_path, "docker run")
    assert len(runs) == 1, calls(tmp_path)
    assert bakta_ref() in runs[0] and "bakta_db" in runs[0]
    assert manifest_items(db)["bakta_light"]["status"] == "ok", res.output


# ---------------------------------------------------------------- an installed wheel (no checkout)


def wheel_site(tmp: Path) -> Path:
    """``src/gippa`` copied with ``workflow/`` bundled inside, as the wheel lays it out (M11-01)."""
    site = tmp / "a" / "b" / "site"
    shutil.copytree(
        REPO / "src" / "gippa", site / "gippa", ignore=shutil.ignore_patterns("__pycache__")
    )
    wf = site / "gippa" / "workflow"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// bundled\n")
    for name in ("containers.config", "images.json"):
        shutil.copy(WORKFLOW / "conf" / name, wf / "conf" / name)
    return site


def run_wheel(site: Path, cwd: Path, code: str, *args: str) -> subprocess.CompletedProcess[str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIPPA_", "MAMBA_", "NXF_"))}
    env["PYTHONPATH"] = str(site)
    return subprocess.run([sys.executable, "-c", code, *args], cwd=cwd, env=env, text=True,
                          capture_output=True, check=False, timeout=120)  # fmt: skip


def test_wheel_toolrunner_docker_no_checkout(tmp_path: Path) -> None:
    site = wheel_site(tmp_path)
    work = Path(tempfile.mkdtemp(dir=tmp_path))
    out = tmp_path / "db" / "out"
    out.mkdir(parents=True)
    code = (
        "import json, sys, gippa\n"
        "from pathlib import Path\n"
        "from gippa.setupdb.run import ToolRunner\n"
        "cmd = ToolRunner(db_dir=Path(sys.argv[1])).command('bakta', sys.argv[2:])\n"
        "print(json.dumps({'file': gippa.__file__, 'cmd': cmd}))\n"
    )
    argv = ["bakta_db", "download", "--output", str(out), "--type", "light"]
    proc = run_wheel(site, work, code, str(tmp_path / "db"), *argv)
    assert proc.returncode == 0, proc.stderr
    got = json.loads(proc.stdout.strip().splitlines()[-1])
    assert got["file"].startswith(str(site)), got["file"]  # the copy, not the checkout
    check_docker_argv(got["cmd"], bakta_ref(), argv, [out])


def test_wheel_setup_images_docker_no_checkout(tmp_path: Path) -> None:
    site = wheel_site(tmp_path)
    work = Path(tempfile.mkdtemp(dir=tmp_path))
    db = tmp_path / "db"
    code = "import sys\nfrom gippa.cli import cli\nsys.argv[0] = 'gippa'\ncli()\n"
    proc = run_wheel(site, work, code, "setup", "--db-dir", str(db), *SKIP_DBS, "--images",
                     "--profile", "docker")  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    pinned = tool_pins(CONTAINERS_CONFIG)
    want = {t for t in default_set(pinned)}
    assert {p[-1] for p in calls(tmp_path, "docker pull")} == {pinned[t] for t in want}
    items = manifest_items(db)
    assert {n for n in items if n.startswith("image:")} == {f"image:{t}" for t in want}


# ---------------------------------------------------------------- setup --check: DB pins


@pytest.mark.parametrize(
    ("item", "write", "old", "new", "args"),
    [
        ("genomad", lambda d: (d / "genomad_db").mkdir() or (d / "genomad_db" / "version.txt").write_text("1.8\n"),  # noqa: E501
         "1.8", "1.9", ["--bakta", "skip", "--no-promotech-models"]),
        ("bakta_light", lambda d: (d / "bakta" / "db-light").mkdir(parents=True) or (d / "bakta" / "db-light" / "version.json").write_text('{"major": 5, "minor": 1, "type": "light"}'),  # noqa: E501
         "5.1", "6.0", ["--no-genomad", "--no-promotech-models"]),
    ],
)  # fmt: skip
def test_check_reports_a_db_pin_change(
    tmp_path: Path, item: str, write: Any, old: str, new: str, args: list[str]
) -> None:
    from gippa.setupdb.items import SOURCES

    assert SOURCES[item].version == new  # the pin this release ships
    db = tmp_path / "db"
    db.mkdir()
    write(db)
    res = setup_cli(["--db-dir", str(db), *args, "--check"])
    assert res.exit_code == 2, res.output
    assert f"{item}: installed {old}, this release pins {new}" in res.output


# ---------------------------------------------------------------- review round (PR #384)


def test_toolrunner_docker_antismash_uses_the_downloader_as_entrypoint(tmp_path: Path) -> None:
    """The antiSMASH image's ENTRYPOINT is a wrapper: docker runs the downloader in its place."""
    db, out = tmp_path / "db", tmp_path / "db" / "antismash.partial"
    out.mkdir(parents=True)
    ref = read_containers_config(CONTAINERS_CONFIG)["ANTISMASH"]
    cmd = _runner(db).command(
        "antismash", ["download-antismash-databases", "--database-dir", str(out)]
    )
    assert cmd[:2] == ["docker", "run"], cmd
    assert cmd[-5:] == [
        "--entrypoint",
        "download-antismash-databases",
        ref,
        "--database-dir",
        str(out),
    ]
    head = cmd[:-5]
    assert "-u" in head and head[head.index("-u") + 1] == f"{os.getuid()}:{os.getgid()}"
    assert f"{out}:{out}" in [head[i + 1] for i, a in enumerate(head) if a in ("-v", "--volume")]


def test_docker_pulls_a_shared_ref_once(tmp_path: Path) -> None:
    """Two tools pinned to one image (as the BLAST and TBLASTN processes are): one pull, and an
    ``ok`` item for each tool under its own name."""
    blast = read_containers_config(CONTAINERS_CONFIG)["BLAST"]
    conf = tmp_path / "conf"
    conf.mkdir()
    (conf / "containers.config").write_text(
        "process {\n"
        f"  withName: 'BLAST' {{ container = '{blast}' }}\n"
        f"  withName: 'TBLASTN' {{ container = '{blast}' }}\n"
        f"  withName: 'MINIPROT' {{ container = '{blast}' }}\n"
        "}\n",
        encoding="utf-8",
    )
    db = tmp_path / "db"
    code = _images().install_images(["blast", "miniprot"], "docker", tmp_path / "c", db, conf)
    assert code == 0
    pulled = [p[-1] for p in calls(tmp_path, "docker pull")]
    assert pulled == [blast]
    items = manifest_items(db)
    for tool in ("blast", "miniprot"):
        assert items[f"image:{tool}"]["status"] == "ok", items
        assert items[f"image:{tool}"]["checksum"] == blast.split("@", 1)[1]


@pytest.mark.parametrize("profile", ["docker", "apptainer"])
def test_rerun_keeps_an_unchanged_ok_record(
    tmp_path: Path, tree: dict[str, str], profile: str
) -> None:
    db = tmp_path / "db"
    args = ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", profile]
    assert setup_cli(args).exit_code == 0
    path = db / "db_manifest.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    for item in data["items"]:
        if item["name"].startswith("image:"):
            item["retrieved"] = "2020-01-01"
    path.write_text(json.dumps(data), encoding="utf-8")
    before = manifest_items(db)
    res = setup_cli(args)
    assert res.exit_code == 0, res.output
    after = manifest_items(db)
    for name, item in before.items():
        if name.startswith("image:"):
            assert after[name] == item, name  # the installed date is not rewritten
            assert after[name]["retrieved"] == "2020-01-01"


def test_images_respects_the_setup_lock(tmp_path: Path, tree: dict[str, str]) -> None:
    import fcntl

    from gippa.setupdb.run import LOCK

    db = tmp_path / "db"
    db.mkdir()
    with (db / LOCK).open("a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "docker"])
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)
    assert res.exit_code == 2, res.output
    assert "another gippa setup is running" in res.output
    assert not calls(tmp_path)
    assert not (db / "db_manifest.json").exists()


def test_apptainer_pull_that_writes_no_image_is_absent(
    tmp_path: Path, tree: dict[str, str]
) -> None:
    fake = tmp_path / "fakebin" / "apptainer"
    fake.write_text('#!/bin/sh\nprintf \'%s\\n\' "apptainer $*" >> "$FAKE_RUNTIME_LOG"\nexit 0\n')
    db = tmp_path / "db"
    res = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer"])
    assert res.exit_code == 2, res.output
    assert calls(tmp_path, "apptainer pull")
    items = manifest_items(db)
    for tool in default_set(tree):
        item = items[f"image:{tool}"]
        assert item["status"] == "absent", item
        assert "wrote no" in item["message"] and simple_name(tree[tool]) in item["message"]


def test_check_bakta_other_minor_is_not_a_pin_change(tmp_path: Path) -> None:
    from gippa.setupdb.items import SOURCES

    major = SOURCES["bakta_light"].version.split(".")[0]
    other = f"{major}.{int(SOURCES['bakta_light'].version.split('.')[1]) + 1}"
    db = tmp_path / "db"
    (db / "bakta" / "db-light").mkdir(parents=True)
    # the wrong type keeps the copy from verifying, so only the minor differs from a pin change
    (db / "bakta" / "db-light" / "version.json").write_text(
        json.dumps({"major": int(major), "minor": int(other.split(".")[1]), "type": "full"})
    )
    res = setup_cli(["--db-dir", str(db), "--no-genomad", "--no-promotech-models", "--check"])
    assert res.exit_code == 2, res.output
    assert "this release pins" not in res.output
    assert f"found v{other} full" in res.output
