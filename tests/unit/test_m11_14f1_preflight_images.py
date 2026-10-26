"""M11-14f1: preflight treats pinned images as the install under the container profiles.

docs/specs/M11-14-container-install.md item 3 (the table), examples 6, 7 and 7b, AC 5 (14f1 part:
every cell of the item 3 table except the 14d ``enabled`` row and the 14c fgsea row), Q-L11.

Conda is enabled only by the ``conda`` profile and by ``prebuilt_envs.config`` (the
``farmshare*`` profiles). Under ``docker`` alone and under plain ``apptainer`` a process that has
only a ``conda`` directive would run on the host without its environment, so preflight reports it
missing with a reason instead of passing.

Interface pinned here (implementer, ``src/giap/preflight/environment.py``):

- ``missing_tools(tools, profile, image_cache, env_root, local_tools) -> list[str]``: sorted
  registry names of ``tools`` (``Sequence[ToolSpec]``) not usable under the comma-joined Nextflow
  ``profile``.
  - ``image_cache: Path``: the apptainer image cache (14a's ``image_cache(environ, db_dir)``:
    ``NXF_APPTAINER_CACHEDIR`` → ``$GIAP_ROOT/containers`` → ``<db_dir>/containers``); a pinned
    image is present under apptainer iff ``<image_cache>/<simpleName>.img`` exists (ADR-0059:
    scheme stripped, ``/`` and ``:`` → ``-``).
  - ``env_root: Path | None``: the directory holding the prebuilt envs
    (``$GIAP_ROOT/micromamba/envs``; ``None`` without ``$GIAP_ROOT``); an env is built iff
    ``<env_root>/giap-<tool>/.giap-built`` equals the sha256 of ``workflow/envs/<tool>.yml``.
  - ``local_tools: Path | None``: ``<db_dir>/local_tools.config`` written by ``setup --uniop``
    (item 1b), or ``None``. The file is one ``profiles { docker {...} apptainer {...} conda {...}
    }``; a block applies only under its own runtime ``container_runtime(profile)``. The docker
    block's ``withName: 'UNIOP' { container = '<tag>' }`` makes uniop present under docker iff
    ``docker image inspect <tag>`` exits 0; the apptainer block's ``container = '<absolute
    .sif>'`` makes it present under apptainer iff the file exists (a relative path is malformed:
    missing, never mapped to a cache file name). No block for the runtime: uniop is missing. Under
    a ``farmshare*`` list the file is ignored (``giap run`` never adds it there, example 19a).
  - Pins are the ``withName`` lines of ``workflow_main().parent / "conf/containers.config"``; the
    function looks ``workflow_main`` up as the module attribute
    ``giap.preflight.environment.workflow_main`` (the tests patch it with a fake workflow tree).
  - ``test``/``stub`` profiles: ``[]``; in-process tools are never missing.
- ``missing_tool_reasons(tools, profile, image_cache, env_root, local_tools) -> dict[str, str]``:
  one reason per missing tool (keys == ``set(missing_tools(...))``). Exact texts:
  - conda-only under ``docker`` / plain ``apptainer``:
    ``conda-only tool under a container profile without conda: <tool>``;
  - pinned image not in the apptainer cache: names ``<pkg> setup --images --profile apptainer``
    (``<pkg>`` = ``[project] name`` of pyproject.toml).
- ``container_runtime(profile) -> "docker" | "apptainer" | None``: the value ``giap run`` passes
  in ``GIAP_CONTAINER_RUNTIME`` (the export itself is 14c).
- ``run_preflight(..., profile=...)`` takes the cache from 14a's helper (no ``$GIAP_ROOT`` needed)
  and puts the reasons into the ``environment.tools`` message.

No test runs real docker/apptainer or touches a registry: a logging fake ``docker`` and
``apptainer`` sit first on ``PATH``. 14f1 stacks on 14a (``image_cache``, ``images.json``).
"""

from __future__ import annotations

import hashlib
import os
import random
import stat
import tempfile
import tomllib
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.annotate.tools import TOOLS, ToolSpec, read_containers_config, required_tools
from giap.config import resolve_config

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
CONTAINERS_CONFIG = WORKFLOW / "conf" / "containers.config"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
MARKER = ".giap-built"
OSTIR_REF = "ghcr.io/khoa-yelo/" + PKG + "-ostir@sha256:" + "a1" * 32
CONDA_REASON = "conda-only tool under a container profile without conda: {}"
SETUP_IMAGES = f"{PKG} setup --images --profile apptainer"

# Profiles without conda (item 3 table, columns 1 and 2) and with conda (column 3).
DOCKER = ["docker", "standard,docker"]
PLAIN_APPTAINER = ["apptainer"]
CONTAINER_NO_CONDA = ["docker", "apptainer"]
WITH_CONDA = ["conda", "farmshare,apptainer", "farmshare_local,apptainer", "conda,docker"]


# ---------------------------------------------------------------- helpers


def _env():  # type: ignore[no-untyped-def]
    from giap.preflight import environment

    return environment


def _fn(name: str) -> Any:
    fn = getattr(_env(), name, None)
    assert fn is not None, f"giap.preflight.environment.{name} is not implemented"
    return fn


def simple_name(ref: str) -> str:
    """Nextflow's apptainer cache file name for ``ref`` (ADR-0059)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def pins(config: Path = CONTAINERS_CONFIG) -> dict[str, str]:
    """Registry tool name → pinned ref from a containers.config."""
    return {
        p.lower(): ref for p, ref in read_containers_config(config).items() if p.lower() in TOOLS
    }


def external(*, exclude: Iterable[str] = ("fgsea",)) -> list[ToolSpec]:
    """Registry tools run as a process; fgsea's row (prefix → env → image) is 14c's."""
    skip = set(exclude)
    return [t for t in TOOLS.values() if not t.in_process and t.name not in skip]


def conda_only(pinned: dict[str, str]) -> list[str]:
    return sorted(t.name for t in external() if t.name not in pinned)


def fill_cache(cache: Path, tools: Iterable[str], pinned: dict[str, str]) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    for name in tools:
        (cache / simple_name(pinned[name])).write_bytes(b"fake sif\n")
    return cache


def build_envs(env_root: Path, tools: Iterable[str]) -> Path:
    """Prebuilt envs with a matching ``.giap-built`` marker (scripts/build_envs.sh)."""
    env_root.mkdir(parents=True, exist_ok=True)
    for name in tools:
        prefix = env_root / TOOLS[name].env_name
        (prefix / "bin").mkdir(parents=True, exist_ok=True)
        sha = hashlib.sha256((WORKFLOW / "envs" / f"{name}.yml").read_bytes()).hexdigest()
        (prefix / MARKER).write_text(sha + "\n")
    return env_root


FAKE_RUNTIME = """#!/bin/sh
printf '%s %s\\n' "$(basename "$0")" "$*" >> "$FAKE_RUNTIME_LOG"
if [ "$(basename "$0")" = docker ] && [ "$1" = image ] && [ "$2" = inspect ]; then
    grep -qxF "$3" "$FAKE_DOCKER_TAGS" 2>/dev/null && exit 0
    echo "Error: No such image: $3" >&2
    exit 1
fi
echo "fake $(basename "$0"): unexpected call" >&2
exit 3
"""


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Never the real $GIAP_ROOT / $GIAP_DB_DIR / caches; fake docker and apptainer on PATH."""
    for var in ("GIAP_DB_DIR", "GIAP_ROOT", "NXF_APPTAINER_CACHEDIR", "GIAP_FGSEA_PREFIX"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.delenv("GIAP_CONTAINER_RUNTIME", raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    for exe in ("docker", "apptainer", "singularity"):
        path = bindir / exe
        path.write_text(FAKE_RUNTIME)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_RUNTIME_LOG", str(tmp_path / "runtime.log"))
    monkeypatch.setenv("FAKE_DOCKER_TAGS", str(tmp_path / "docker_tags.txt"))


def runtime_calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "runtime.log"
    return log.read_text().splitlines() if log.exists() else []


@pytest.fixture
def ostir_pinned(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, str]:
    """A fake workflow tree whose containers.config also pins OSTIR (a GHCR image)."""
    wf = tmp_path / "wf"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// fake\n")
    (wf / "envs").symlink_to(WORKFLOW / "envs")
    text = CONTAINERS_CONFIG.read_text(encoding="utf-8").rstrip()
    assert text.endswith("}")
    text = text[:-1] + f"  withName: 'OSTIR' {{ container = '{OSTIR_REF}' }}\n}}\n"
    config = wf / "conf" / "containers.config"
    config.write_text(text, encoding="utf-8")
    monkeypatch.setattr(_env(), "workflow_main", lambda: wf / "main.nf")
    pinned = pins(config)
    assert pinned["ostir"] == OSTIR_REF
    return pinned


def call(
    name: str,
    tools: list[ToolSpec],
    profile: str,
    cache: Path,
    env_root: Path | None,
    local_tools: Path | None = None,
) -> Any:
    return _fn(name)(tools, profile, cache, env_root, local_tools)


def empty(tmp_path: Path, name: str) -> Path:
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------- signature


def test_signature_takes_cache_env_root_and_local_tools_by_keyword(tmp_path: Path) -> None:
    got = _fn("missing_tools")(
        [TOOLS["bakta"]],
        "docker",
        image_cache=empty(tmp_path, "cache"),
        env_root=None,
        local_tools=None,
    )
    assert got == []


# ---------------------------------------------------------------- examples 6, 7, 7b


def test_example_6_pinned_ostir_under_docker_needs_no_env(
    tmp_path: Path, ostir_pinned: dict[str, str]
) -> None:
    cache, envs = empty(tmp_path, "cache"), empty(tmp_path, "envs")
    assert call("missing_tools", [TOOLS["ostir"]], "docker", cache, envs) == []
    assert call("missing_tool_reasons", [TOOLS["ostir"]], "docker", cache, envs) == {}
    assert runtime_calls(tmp_path) == [], "docker was queried for a pinned image"


def test_example_7_pinned_ostir_under_apptainer_without_img(
    tmp_path: Path, ostir_pinned: dict[str, str]
) -> None:
    cache, envs = empty(tmp_path, "cache"), empty(tmp_path, "envs")
    assert call("missing_tools", [TOOLS["ostir"]], "apptainer", cache, envs) == ["ostir"]
    reasons = call("missing_tool_reasons", [TOOLS["ostir"]], "apptainer", cache, envs)
    assert set(reasons) == {"ostir"}
    assert SETUP_IMAGES in reasons["ostir"], reasons
    assert runtime_calls(tmp_path) == [], "preflight must check the cache, not pull"


def test_example_7_the_img_in_the_cache_makes_it_present(
    tmp_path: Path, ostir_pinned: dict[str, str]
) -> None:
    cache = fill_cache(tmp_path / "cache", ["ostir"], ostir_pinned)
    assert (cache / ("ghcr.io-khoa-yelo-" + PKG + "-ostir@sha256-" + "a1" * 32 + ".img")).is_file()
    assert call("missing_tools", [TOOLS["ostir"]], "apptainer", cache, None) == []


@pytest.mark.parametrize("profile", CONTAINER_NO_CONDA)
def test_example_7b_conda_only_tool_with_a_built_env_is_missing(
    tmp_path: Path, profile: str
) -> None:
    """A prebuilt env does not help where conda is off: the process would run on the host."""
    pinned = pins()
    tool = "promotech"
    assert tool not in pinned, "promotech is expected to be conda-only on this branch"
    envs = build_envs(tmp_path / "envs", [tool])
    cache = empty(tmp_path, "cache")
    assert call("missing_tools", [TOOLS[tool]], profile, cache, envs) == [tool]
    reasons = call("missing_tool_reasons", [TOOLS[tool]], profile, cache, envs)
    assert reasons == {tool: CONDA_REASON.format(tool)}


# ---------------------------------------------------------------- item 3 table, row by row


@pytest.mark.parametrize("profile", DOCKER)
def test_row_pinned_image_docker_is_present_with_nothing_installed(
    tmp_path: Path, profile: str
) -> None:
    pinned = pins()
    tools = [TOOLS[n] for n in sorted(pinned)]
    assert call("missing_tools", tools, profile, empty(tmp_path, "cache"), None) == []
    assert runtime_calls(tmp_path) == []


def test_row_pinned_image_apptainer_present_iff_img_in_cache(tmp_path: Path) -> None:
    pinned = pins()
    names = sorted(pinned)
    have, lack = names[::2], names[1::2]
    cache = fill_cache(tmp_path / "cache", have, pinned)
    tools = [TOOLS[n] for n in names]
    assert call("missing_tools", tools, "apptainer", cache, None) == sorted(lack)
    reasons = call("missing_tool_reasons", tools, "apptainer", cache, None)
    assert set(reasons) == set(lack)
    for name in lack:
        assert SETUP_IMAGES in reasons[name], reasons


def test_row_pinned_image_only_the_given_cache_counts(tmp_path: Path) -> None:
    """The image under the old ``$GIAP_ROOT/containers`` does not count when the cache is
    elsewhere: the cache is what Nextflow is told to use (``apptainer.cacheDir``)."""
    pinned = pins()
    root = tmp_path / "root"
    fill_cache(root / "containers", ["bakta"], pinned)
    envs = build_envs(root / "micromamba" / "envs", [])
    cache = empty(tmp_path, "nxf_cache")
    assert call("missing_tools", [TOOLS["bakta"]], "apptainer", cache, envs) == ["bakta"]


def test_row_pinned_image_digest_must_match(tmp_path: Path) -> None:
    pinned = pins()
    cache = empty(tmp_path, "cache")
    other = pinned["bakta"].split("@sha256:")[0] + "@sha256:" + "f" * 64
    (cache / simple_name(other)).write_bytes(b"other digest\n")
    assert call("missing_tools", [TOOLS["bakta"]], "apptainer", cache, None) == ["bakta"]


@pytest.mark.parametrize("profile", CONTAINER_NO_CONDA)
def test_row_pinned_image_and_env_the_env_is_ignored(
    tmp_path: Path, ostir_pinned: dict[str, str], profile: str
) -> None:
    """Pinned OSTIR plus an old prebuilt env: the image decides (edge case "old env, new
    image"). Under docker present without the env mattering; under apptainer missing without
    the ``.img`` although the env is built."""
    envs = build_envs(tmp_path / "envs", ["ostir"])
    cache = empty(tmp_path, "cache")
    want = [] if profile == "docker" else ["ostir"]
    assert call("missing_tools", [TOOLS["ostir"]], profile, cache, envs) == want
    if want:
        reasons = call("missing_tool_reasons", [TOOLS["ostir"]], profile, cache, envs)
        assert SETUP_IMAGES in reasons["ostir"]
        assert "conda-only" not in reasons["ostir"]


@pytest.mark.parametrize("profile", CONTAINER_NO_CONDA)
def test_row_conda_only_without_env_is_missing_with_the_conda_reason(
    tmp_path: Path, profile: str
) -> None:
    pinned = pins()
    names = conda_only(pinned)
    assert names, "no conda-only tool left to test"
    tools = [TOOLS[n] for n in names]
    cache = empty(tmp_path, "cache")
    assert call("missing_tools", tools, profile, cache, None) == names
    reasons = call("missing_tool_reasons", tools, profile, cache, None)
    assert reasons == {n: CONDA_REASON.format(n) for n in names}


@pytest.mark.parametrize("profile", WITH_CONDA)
def test_row_conda_only_with_conda_needs_the_prebuilt_env(tmp_path: Path, profile: str) -> None:
    """``conda`` or ``farmshare*``: as today, a conda-only tool is present iff its env is built."""
    pinned = pins()
    names = conda_only(pinned)
    built, unbuilt = names[::2], names[1::2]
    envs = build_envs(tmp_path / "envs", built)
    tools = [TOOLS[n] for n in names]
    got = call("missing_tools", tools, profile, empty(tmp_path, "cache"), envs)
    assert got == sorted(unbuilt)
    reasons = call("missing_tool_reasons", tools, profile, empty(tmp_path, "cache"), envs)
    assert set(reasons) == set(unbuilt)
    for name in unbuilt:
        assert "conda-only tool under a container profile" not in reasons[name]


@pytest.mark.parametrize("profile", WITH_CONDA)
def test_row_conda_only_with_conda_and_no_env_root(tmp_path: Path, profile: str) -> None:
    """No ``$GIAP_ROOT`` (``env_root=None``): every conda-only tool is missing."""
    names = conda_only(pins())
    tools = [TOOLS[n] for n in names]
    assert call("missing_tools", tools, profile, empty(tmp_path, "cache"), None) == names


def test_row_stale_marker_is_still_missing_under_conda(tmp_path: Path) -> None:
    envs = build_envs(tmp_path / "envs", ["promotech"])
    (envs / TOOLS["promotech"].env_name / MARKER).write_text("0" * 64 + "\n")
    got = call("missing_tools", [TOOLS["promotech"]], "conda", empty(tmp_path, "c"), envs)
    assert got == ["promotech"]


# ---------------------------------------------------------------- uniop with a local build


def _local_tools(db: Path, *, multiline: bool = False, **blocks: str) -> Path:
    """``<db>/local_tools.config`` in the item 1b format: one ``profiles { … }`` with one block
    per runtime (keyword order kept), each holding ``withName: 'UNIOP' { <body> }``."""
    db.mkdir(parents=True, exist_ok=True)
    lines = [f"// written by {PKG} setup --uniop", "profiles {"]
    for runtime, body in blocks.items():
        if multiline:
            stmts = [s.strip() for s in body.split(";") if s.strip()]
            lines += [f"    {runtime} {{", "        process {", "            withName: 'UNIOP' {"]
            lines += [f"                {s}" for s in stmts]
            lines += ["            }", "        }", "    }"]
        else:
            lines.append(f"    {runtime:<9} {{ process {{ withName: 'UNIOP' {{ {body} }} }} }}")
    lines.append("}")
    path = db / "local_tools.config"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _container_block(tag: str) -> str:
    return f"container = '{tag}'; conda = ''"


LOCAL_TAG = f"localhost/{PKG}-uniop:{TOOLS['uniop'].version}-0123456789ab"


def test_row_uniop_local_docker_tag_present(tmp_path: Path) -> None:
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    local = _local_tools(tmp_path / "db", docker=_container_block(LOCAL_TAG))
    got = call("missing_tools", [TOOLS["uniop"]], "docker", empty(tmp_path, "c"), None, local)
    assert got == []
    calls = runtime_calls(tmp_path)
    assert f"docker image inspect {LOCAL_TAG}" in calls, calls
    assert not any(" pull" in c or " build" in c or " push" in c for c in calls), calls


def test_row_uniop_local_docker_tag_absent(tmp_path: Path) -> None:
    (tmp_path / "docker_tags.txt").write_text("localhost/other:1\n")
    local = _local_tools(tmp_path / "db", docker=_container_block(LOCAL_TAG))
    cache = empty(tmp_path, "c")
    got = call("missing_tools", [TOOLS["uniop"]], "docker", cache, None, local)
    assert got == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], "docker", cache, None, local)
    assert set(reasons) == {"uniop"}
    assert LOCAL_TAG in reasons["uniop"], reasons
    assert f"docker image inspect {LOCAL_TAG}" in runtime_calls(tmp_path)


def test_row_uniop_local_sif_present_iff_file_exists(tmp_path: Path) -> None:
    sif = tmp_path / "cache" / f"localhost-{PKG}-uniop-{TOOLS['uniop'].version}-0123456789ab.sif"
    sif.parent.mkdir(parents=True)
    local = _local_tools(tmp_path / "db", apptainer=_container_block(str(sif)))
    cache = sif.parent
    got = call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local)
    assert got == ["uniop"]
    sif.write_bytes(b"fake sif\n")
    assert call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local) == []
    assert runtime_calls(tmp_path) == []


@pytest.mark.parametrize("profile", CONTAINER_NO_CONDA)
def test_row_uniop_local_micromamba_build_under_container_profile(
    tmp_path: Path, profile: str
) -> None:
    """A micromamba-built local UniOP works only with a conda-enabled profile (item 1b)."""
    prefix = tmp_path / "cache" / "envs" / "giap-uniop"
    (prefix / "bin").mkdir(parents=True)
    local = _local_tools(tmp_path / "db", conda=f"conda = '{prefix}'")
    cache = empty(tmp_path, "c")
    assert call("missing_tools", [TOOLS["uniop"]], profile, cache, None, local) == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], profile, cache, None, local)
    assert reasons == {"uniop": CONDA_REASON.format("uniop")}


def test_row_uniop_without_local_tools_under_docker_is_conda_only(tmp_path: Path) -> None:
    envs = build_envs(tmp_path / "envs", ["uniop"])
    cache = empty(tmp_path, "c")
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], "docker", cache, envs, None)
    assert reasons == {"uniop": CONDA_REASON.format("uniop")}


def test_row_local_tools_config_does_not_affect_other_tools(tmp_path: Path) -> None:
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    local = _local_tools(tmp_path / "db", docker=_container_block(LOCAL_TAG))
    names = [n for n in conda_only(pins()) if n != "uniop"]
    tools = [TOOLS[n] for n in names]
    got = call("missing_tools", tools, "docker", empty(tmp_path, "c"), None, local)
    assert got == names


def _two_blocks(tmp_path: Path, order: str, *, multiline: bool = False) -> tuple[Path, Path]:
    """A docker block (``LOCAL_TAG``) and an apptainer block (an absolute ``.sif``)."""
    sif = tmp_path / "builds" / "uniop.sif"
    sif.parent.mkdir(parents=True, exist_ok=True)
    blocks = {"docker": _container_block(LOCAL_TAG), "apptainer": _container_block(str(sif))}
    if order == "apptainer-first":
        blocks = dict(reversed(blocks.items()))
    return _local_tools(tmp_path / "db", multiline=multiline, **blocks), sif


ORDERS = ["docker-first", "apptainer-first"]


@pytest.mark.parametrize("multiline", [False, True], ids=["one-line", "multi-line"])
@pytest.mark.parametrize("order", ORDERS)
def test_b1_two_blocks_docker_uses_its_own_block(
    tmp_path: Path, order: str, multiline: bool
) -> None:
    """Item 1b: under docker only the docker block counts; its tag goes to docker image inspect."""
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    local, sif = _two_blocks(tmp_path, order, multiline=multiline)
    cache = empty(tmp_path, "c")
    assert call("missing_tool_reasons", [TOOLS["uniop"]], "docker", cache, None, local) == {}
    assert call("missing_tools", [TOOLS["uniop"]], "docker", cache, None, local) == []
    calls = runtime_calls(tmp_path)
    assert f"docker image inspect {LOCAL_TAG}" in calls, calls
    assert not any(str(sif) in c for c in calls), calls


@pytest.mark.parametrize("multiline", [False, True], ids=["one-line", "multi-line"])
@pytest.mark.parametrize("order", ORDERS)
def test_b1_two_blocks_apptainer_uses_its_own_block(
    tmp_path: Path, order: str, multiline: bool
) -> None:
    """Item 1b: under apptainer only the apptainer block's absolute ``.sif`` counts."""
    local, sif = _two_blocks(tmp_path, order, multiline=multiline)
    cache = empty(tmp_path, "c")
    assert call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local) == ["uniop"]
    sif.write_bytes(b"fake sif\n")
    assert call("missing_tool_reasons", [TOOLS["uniop"]], "apptainer", cache, None, local) == {}
    assert call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local) == []
    assert runtime_calls(tmp_path) == []


@pytest.mark.parametrize("order", ORDERS)
def test_b1_two_blocks_the_other_runtime_never_satisfies(tmp_path: Path, order: str) -> None:
    """The other block's artifact does not count: an accepted docker tag does not help apptainer,
    an existing ``.sif`` does not help docker."""
    local, sif = _two_blocks(tmp_path, order)
    sif.write_bytes(b"fake sif\n")
    (tmp_path / "docker_tags.txt").write_text("localhost/other:1\n")
    cache = empty(tmp_path, "c")
    assert call("missing_tools", [TOOLS["uniop"]], "docker", cache, None, local) == ["uniop"]
    sif.unlink()
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    got = call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local)
    assert got == ["uniop"]


def test_b1_docker_only_file_under_apptainer_is_missing(tmp_path: Path) -> None:
    """Item 1b: "a docker build followed by -profile apptainer finds no apptainer block", even if
    the cache holds a file named for the docker tag; apptainer never pulls ``docker://localhost``."""
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    local = _local_tools(tmp_path / "db", docker=_container_block(LOCAL_TAG))
    cache = empty(tmp_path, "c")
    (cache / simple_name(LOCAL_TAG)).write_bytes(b"fake sif\n")
    assert call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local) == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], "apptainer", cache, None, local)
    assert set(reasons) == {"uniop"}
    assert not any(c.startswith("apptainer ") and " pull" in c for c in runtime_calls(tmp_path))


def test_b1_apptainer_only_file_under_docker_is_missing(tmp_path: Path) -> None:
    sif = tmp_path / "builds" / "uniop.sif"
    sif.parent.mkdir(parents=True)
    sif.write_bytes(b"fake sif\n")
    (tmp_path / "docker_tags.txt").write_text(f"{sif}\n{LOCAL_TAG}\n")
    local = _local_tools(tmp_path / "db", apptainer=_container_block(str(sif)))
    cache = empty(tmp_path, "c")
    assert call("missing_tools", [TOOLS["uniop"]], "docker", cache, None, local) == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], "docker", cache, None, local)
    assert set(reasons) == {"uniop"}
    assert not any(str(sif) in c for c in runtime_calls(tmp_path))


@pytest.mark.parametrize("profile", CONTAINER_NO_CONDA)
def test_b1_conda_block_only_is_missing_under_a_container_profile(
    tmp_path: Path, profile: str
) -> None:
    """Only a conda block, no block for ``R(P)``: uniop stays missing under docker/apptainer."""
    prefix = tmp_path / "envs" / "giap-uniop"
    (prefix / "bin").mkdir(parents=True)
    local = _local_tools(tmp_path / "db", conda=f"conda = '{prefix}'")
    got = call("missing_tools", [TOOLS["uniop"]], profile, empty(tmp_path, "c"), None, local)
    assert got == ["uniop"]


@pytest.mark.parametrize("rel", ["uniop.sif", "builds/uniop.sif", LOCAL_TAG])
def test_b1_relative_apptainer_container_is_malformed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, rel: str
) -> None:
    """Item 1b: the apptainer block holds an absolute ``.sif`` path. A relative one is malformed:
    uniop is missing and the reason names the path. It is not resolved against the cwd or mapped
    to ``<cache>/<simpleName>.img`` (no basis in the spec)."""
    work = empty(tmp_path, "work")
    (work / rel).parent.mkdir(parents=True, exist_ok=True)
    (work / rel).write_bytes(b"fake sif\n")
    monkeypatch.chdir(work)
    cache = empty(tmp_path, "c")
    (cache / simple_name(rel)).write_bytes(b"fake sif\n")
    local = _local_tools(tmp_path / "db", apptainer=_container_block(rel))
    assert call("missing_tools", [TOOLS["uniop"]], "apptainer", cache, None, local) == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], "apptainer", cache, None, local)
    assert set(reasons) == {"uniop"}
    assert rel in reasons["uniop"], reasons
    assert simple_name(rel) not in reasons["uniop"], reasons


FARMSHARE_APPTAINER = ["farmshare,apptainer", "farmshare_local,apptainer"]


@pytest.mark.parametrize("profile", FARMSHARE_APPTAINER)
def test_b2_farmshare_ignores_local_tools_config(tmp_path: Path, profile: str) -> None:
    """AC 5, example 19a: under a ``farmshare*`` list ``giap run`` never adds local_tools.config,
    so a usable apptainer/docker block does not satisfy uniop without the prebuilt env."""
    local, sif = _two_blocks(tmp_path, "docker-first")
    sif.write_bytes(b"fake sif\n")
    (tmp_path / "docker_tags.txt").write_text(LOCAL_TAG + "\n")
    envs = build_envs(tmp_path / "envs", ["promotech"])
    cache = empty(tmp_path, "c")
    assert call("missing_tools", [TOOLS["uniop"]], profile, cache, envs, local) == ["uniop"]
    assert call("missing_tools", [TOOLS["uniop"]], profile, cache, None, local) == ["uniop"]
    reasons = call("missing_tool_reasons", [TOOLS["uniop"]], profile, cache, envs, local)
    assert set(reasons) == {"uniop"}
    assert str(sif) not in reasons["uniop"], reasons
    build_envs(envs, ["uniop"])
    assert call("missing_tools", [TOOLS["uniop"]], profile, cache, envs, local) == []


# ---------------------------------------------------------------- in-process, test/stub


@pytest.mark.parametrize("profile", [*CONTAINER_NO_CONDA, *WITH_CONDA])
def test_row_in_process_tools_are_never_missing(tmp_path: Path, profile: str) -> None:
    inproc = [t for t in TOOLS.values() if t.in_process]
    assert inproc, "codon_bias is expected to be in-process"
    assert call("missing_tools", inproc, profile, empty(tmp_path, "c"), None) == []


@pytest.mark.parametrize("profile", ["test", "stub", "test,docker", "stub,apptainer"])
def test_test_and_stub_profiles_need_no_tools(tmp_path: Path, profile: str) -> None:
    got = call("missing_tools", external(), profile, empty(tmp_path, "c"), None)
    assert got == []


# ---------------------------------------------------------------- properties


@settings(
    max_examples=30,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    cached=st.sets(st.sampled_from(sorted(pins()))),
    built=st.sets(st.sampled_from(conda_only(pins()))),
)
def test_property_plain_apptainer(cached: set[str], built: set[str]) -> None:
    """Plain apptainer: missing = pinned tools without an ``.img``, plus every conda-only tool,
    whatever envs are built."""
    pinned = pins()
    with tempfile.TemporaryDirectory() as d:
        cache = fill_cache(Path(d) / "cache", cached, pinned)
        envs = build_envs(Path(d) / "envs", built)
        got = call("missing_tools", external(), "apptainer", cache, envs)
        reasons = call("missing_tool_reasons", external(), "apptainer", cache, envs)
    want = sorted((set(pinned) - cached) | set(conda_only(pinned)))
    assert got == want
    assert sorted(reasons) == want
    for name in conda_only(pinned):
        assert reasons[name] == CONDA_REASON.format(name)


@settings(
    max_examples=30,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    cached=st.sets(st.sampled_from(sorted(pins()))),
    built=st.sets(st.sampled_from(conda_only(pins()))),
)
def test_property_docker_depends_on_nothing_installed(cached: set[str], built: set[str]) -> None:
    """Docker: missing is exactly the conda-only tools, independent of cache and envs."""
    pinned = pins()
    with tempfile.TemporaryDirectory() as d:
        cache = fill_cache(Path(d) / "cache", cached, pinned)
        envs = build_envs(Path(d) / "envs", built)
        got = call("missing_tools", external(), "docker", cache, envs)
    assert got == conda_only(pinned)


@settings(max_examples=30, deadline=None)
@given(
    built=st.sets(st.sampled_from(conda_only(pins()))),
    profile=st.sampled_from(WITH_CONDA),
)
def test_property_missing_and_reasons_agree(built: set[str], profile: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        cache = Path(d) / "cache"
        cache.mkdir()
        envs = build_envs(Path(d) / "envs", built)
        got = call("missing_tools", external(), profile, cache, envs)
        reasons = call("missing_tool_reasons", external(), profile, cache, envs)
    assert got == sorted(reasons)
    assert got == sorted(set(got))


# ---------------------------------------------------------------- GIAP_CONTAINER_RUNTIME


@pytest.mark.parametrize(
    ("profile", "want"),
    [
        ("docker", "docker"),
        ("test,docker", "docker"),
        ("apptainer", "apptainer"),
        ("farmshare,apptainer", "apptainer"),
        ("farmshare_local,apptainer", "apptainer"),
        ("conda", None),
        ("standard", None),
        ("test,stub", None),
    ],
)
def test_container_runtime_from_profile(profile: str, want: str | None) -> None:
    assert _fn("container_runtime")(profile) == want


# ---------------------------------------------------------------- run_preflight wiring


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


def _inputs(tmp_path: Path) -> Any:
    from giap.preflight.run import PreflightInputs

    rng = random.Random(0)
    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(parents=True, exist_ok=True)
    seq = "".join(rng.choice("ACGT") for _ in range(120_000))
    g.write_text(f">chr1 [topology=circular]\n{seq}\n")
    p = tmp_path / "in" / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return PreflightInputs(genome=g, proteins=p)


def _db(tmp_path: Path) -> Path:
    from helpers.fake_env import write_setup_manifest

    return write_setup_manifest(tmp_path / "db")


def _tools_check(tmp_path: Path, db: Path, profile: str) -> Any:
    from giap.preflight.run import run_preflight

    cfg = resolve_config(  # 14d2: the default is auto; pin the enabled: true path
        overrides=["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1",
                   "annotation.uniop.enabled=true"]
    )  # fmt: skip
    report, _ = run_preflight(_inputs(tmp_path), cfg, None, "run", tmp_path / "OUT", db,
                              profile=profile)  # fmt: skip
    checks = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert len(checks) == 1, checks
    return checks[0]


def _needed() -> tuple[list[str], list[str]]:
    """(pinned, conda-only) registry names a default run needs."""
    cfg = resolve_config(overrides=[])
    names = sorted(t.name for t in required_tools(cfg))
    pinned = pins()
    return [n for n in names if n in pinned], [n for n in names if n not in pinned]


def _assert_only_conda_only_missing(check: Any) -> None:
    need_pinned, need_conda = _needed()
    assert need_pinned and need_conda, (need_pinned, need_conda)
    assert check.level == "error" and check.exit_code == 2, check
    missing = set(check.value["missing"])
    assert set(need_conda) <= missing, missing
    assert not set(need_pinned) & missing, missing
    for name in need_conda:
        assert CONDA_REASON.format(name) in check.message, check.message


def test_run_preflight_apptainer_uses_nxf_cache_without_giap_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Q-L11: a pip user, no ``$GIAP_ROOT``; the images are in ``NXF_APPTAINER_CACHEDIR``."""
    need_pinned, _ = _needed()
    cache = fill_cache(tmp_path / "nxf", need_pinned, pins())
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(cache))
    _assert_only_conda_only_missing(_tools_check(tmp_path, _db(tmp_path), "apptainer"))


def test_run_preflight_apptainer_cache_falls_back_to_db_dir(tmp_path: Path) -> None:
    """No ``NXF_APPTAINER_CACHEDIR``, no ``$GIAP_ROOT``: ``<db_dir>/containers``."""
    need_pinned, _ = _needed()
    db = _db(tmp_path)
    fill_cache(db / "containers", need_pinned, pins())
    _assert_only_conda_only_missing(_tools_check(tmp_path, db, "apptainer"))


def test_run_preflight_apptainer_cache_from_giap_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    need_pinned, _ = _needed()
    root = tmp_path / "root"
    fill_cache(root / "containers", need_pinned, pins())
    monkeypatch.setenv("GIAP_ROOT", str(root))
    _assert_only_conda_only_missing(_tools_check(tmp_path, _db(tmp_path), "apptainer"))


def test_run_preflight_apptainer_missing_image_message_names_setup_images(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(empty(tmp_path, "nxf")))
    check = _tools_check(tmp_path, _db(tmp_path), "apptainer")
    need_pinned, _ = _needed()
    assert set(need_pinned) <= set(check.value["missing"]), check.value
    assert SETUP_IMAGES in check.message, check.message


def test_run_preflight_docker_with_nothing_installed(tmp_path: Path) -> None:
    """Under docker the pinned tools need nothing local; the conda-only ones are reported with
    the reason instead of running on the host."""
    _assert_only_conda_only_missing(_tools_check(tmp_path, _db(tmp_path), "docker"))
    assert not [c for c in runtime_calls(tmp_path) if " pull" in c], runtime_calls(tmp_path)


def test_run_preflight_farmshare_unchanged_with_full_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """FarmShare behaviour is unchanged: images in ``$GIAP_ROOT/containers`` and the prebuilt
    envs under ``$GIAP_ROOT/micromamba/envs`` pass."""
    need_pinned, need_conda = _needed()
    root = tmp_path / "root"
    fill_cache(root / "containers", need_pinned, pins())
    build_envs(root / "micromamba" / "envs", need_conda)
    monkeypatch.setenv("GIAP_ROOT", str(root))
    check = _tools_check(tmp_path, _db(tmp_path), "farmshare,apptainer")
    assert check.level == "ok", check
