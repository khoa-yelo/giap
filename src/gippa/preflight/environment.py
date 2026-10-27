"""Environment preflight rules (GIPPA_SPEC.md §7.0 Environment rows; ADR-0005 items 3-4).

``environment.tools``: every binary, environment, database and model of the enabled modules is
found at a usable version → else error (exit 2). Binaries and environments come from the M2-01 tool
registry (:func:`missing_tools`: pinned images, in the apptainer image cache under apptainer;
prebuilt micromamba envs with a matching build marker where conda is enabled); databases and
models from the ``gippa setup`` manifest (:func:`read_manifest`, run-wiring lead decision 4).
``annotation.uniop.enabled`` (M11-14d, ADR-0061 decision 5) is resolved against whether UniOP is
installed (:func:`resolve_uniop`, pure) after one ``--help`` run of UniOP in the profile's
runtime (:func:`probe_uniop`); the reason UniOP is missing includes a failed or absent
``image:uniop`` setup item (:func:`uniop_setup_record`).
``environment.resources``: the only warning is RAM below the 24 GiB that genome-wide Promotech
needs (the caller passes it as ``need_ram_bytes``, see :func:`resource_needs`). The check holds
the verdict only, never the raw RAM and disk numbers, so ``preflight.json`` is identical across
machines (§6.6); the caller logs the raw numbers (§6.7). Rules take an injected
:class:`EnvProbe`; the real sources are read by :func:`read_manifest`
(``$GIPPA_DB_DIR/db_manifest.json``) and :func:`probe_resources`.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from gippa.annotate.tools import TOOLS, ToolSpec, read_containers_config
from gippa.config import GippaConfig
from gippa.io.image_key import key as image_key
from gippa.preflight.checks import Check, JsonValue, Level
from gippa.runner.nextflow import workflow_main
from gippa.setupdb.items import RESOURCES, SOURCES

__all__ = [
    "GENOME_PROMOTECH_RAM_BYTES",
    "GIB",
    "LOCAL_TOOLS",
    "MANIFEST",
    "RESOURCES_RULE",
    "TOOLS_RULE",
    "EnvProbe",
    "ProbeError",
    "Requirement",
    "UniopResolution",
    "container_runtime",
    "environment_checks",
    "local_tools_applies",
    "missing_tool_reasons",
    "missing_tools",
    "probe_resources",
    "probe_uniop",
    "read_manifest",
    "requirements",
    "resolve_uniop",
    "resource_needs",
    "uniop_local_problem",
    "uniop_runtime",
    "uniop_runtime_problem",
    "uniop_setup_command",
    "uniop_setup_record",
]

GIB = 1024**3
GENOME_PROMOTECH_RAM_BYTES = 24 * GIB
PROMOTER_MODES = ("targeted", "genome")
TOOLS_RULE = "binaries, environments and DB versions of enabled modules found → else error (exit 2)"
RESOURCES_RULE = (
    "RAM ≥ 24 GiB for genome-wide Promotech → else warn; free disk and RAM logged (events.jsonl)"
)
MANIFEST = "db_manifest.json"
_CGROUP_LIMITS = (
    Path("/sys/fs/cgroup/memory.max"),  # cgroup v2
    Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),  # cgroup v1
)
_MB = 1024**2

Kind = Literal["binary", "env", "db", "model"]


@dataclass(frozen=True, slots=True)
class Requirement:
    """One tool, environment, database or model an enabled module needs."""

    name: str
    module: str
    kind: Kind
    min_version: str | None


@dataclass(frozen=True, slots=True)
class EnvProbe:
    """What the environment provides: version found per requirement (``None`` if missing)."""

    found: Mapping[str, str | None]
    free_disk_bytes: int
    total_ram_bytes: int


class ProbeError(Exception):
    """The environment probe could not read its sources."""


# ---------------------------------------------------------------- real probe


def read_manifest(db_dir: str | None) -> dict[str, str | None] | None:
    """Versions found per name in ``db_dir/db_manifest.json``; ``None`` if there is no manifest.

    The ``gippa setup`` manifest (``{"schema_version", "items": [...]}``, ``ok`` items only) or
    the provisional ``{"<name>": {"version": "<v>", ...}}`` / ``"<name>": "<v>"`` form. Raises
    :class:`ProbeError` when the manifest is unreadable or not a JSON object.
    """
    if not db_dir or not (Path(db_dir) / MANIFEST).exists():
        return None
    path = Path(db_dir) / MANIFEST
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:  # UnicodeDecodeError and JSONDecodeError included
        raise ProbeError(f"{MANIFEST} is unreadable: {type(exc).__name__}: {exc}") from exc
    if not isinstance(data, dict):
        raise ProbeError(f"{MANIFEST} is not a JSON object")
    if isinstance(data.get("items"), list):  # `gippa setup` manifest (M2-08)
        return _setup_manifest(data["items"])
    found: dict[str, str | None] = {}
    for name, entry in data.items():
        v = entry.get("version") if isinstance(entry, dict) else entry
        found[name] = v if isinstance(v, str) else None
    return found


# `gippa setup` item names → requirement names (M2-08 lead decision 4)
_SETUP_NAMES = {"bakta_light": "bakta_db", "bakta_full": "bakta_db", "genomad": "genomad_db",
                "antismash": "antismash_db"}  # fmt: skip
# Item groups → one requirement; every part must be ``ok`` (run-wiring lead decision 4b).
_SETUP_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    *RESOURCES,
    ("promotech", ("promotech_models",)),
)


def _setup_manifest(items: list[object]) -> dict[str, str | None]:
    """Versions of the ``status: "ok"`` items of a ``gippa setup`` manifest, by requirement name.

    Grouped items (:data:`_SETUP_GROUPS`) also keep their item name; the group's requirement
    counts only when every part is ``ok``, its version joining the distinct part versions with
    ``+``.
    """
    ok: dict[str, str | None] = {}
    for item in items:
        if isinstance(item, dict) and item.get("status") == "ok":
            version = item.get("version")
            ok[str(item.get("name"))] = version if isinstance(version, str) else None
    found: dict[str, str | None] = {}
    for name, version in ok.items():
        found[_SETUP_NAMES.get(name, name)] = version
    for req, parts in _SETUP_GROUPS:
        if all(part in ok for part in parts):
            versions = [v for v in dict.fromkeys(ok[part] for part in parts) if v is not None]
            found[req] = "+".join(versions) if versions else None
    return found


# ---------------------------------------------------------------- tools from the registry

MARKER = ".gippa-built"
LOCAL_TOOLS = "local_tools.config"  # written into db_dir by setup --uniop
_NO_TOOL_PROFILES = frozenset({"test", "stub"})


def _image_file_name(ref: str) -> str:
    """``a/b@sha256:<d>``, with or without ``docker://`` → ``a-b@sha256-<d>.img`` (Nextflow)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _env_built(tool: ToolSpec, env_root: Path | None, envs_dir: Path) -> bool:
    """The prebuilt env's ``.gippa-built`` marker equals the sha256 of its env file."""
    if env_root is None:
        return False
    marker = env_root / tool.env_name / MARKER
    want = _sha256(envs_dir / f"{tool.name}.yml")
    try:
        have = marker.read_text(encoding="utf-8").strip()
    except OSError:
        return False
    return want is not None and have == want


_PKG = __name__.split(".", 1)[0]
_CONDA_PROFILES = frozenset({"conda", "farmshare", "farmshare_local"})
# Started by a GIPPA step (gippa.io.fgsea_r), not by a Nextflow process: the prebuilt env counts under
# every profile, and under docker / plain apptainer its image from images.json (M11-14 14c).
_STEP_STARTED = frozenset({"fgsea"})
_WITH_NAME = re.compile(r"withName\s*:\s*(['\"]?)([A-Za-z0-9_]+)\1")
_CONTAINER = re.compile(r"(?:^|[\s;])container\s*=\s*(['\"])(.*?)\1")
_CONDA = re.compile(r"(?:^|[\s;])conda\s*=\s*(['\"])(.*?)\1")


def _strip_comments(text: str) -> str:
    """``text`` with ``//`` and ``/* */`` comments outside quotes blanked."""
    out: list[str] = []
    i, quote = 0, ""
    while i < len(text):
        c = text[i]
        if quote:
            out.append(c)
            quote = "" if c == quote and text[i - 1] != "\\" else quote
        elif c in "'\"":
            out.append(c)
            quote = c
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = len(text) if j < 0 else j
            continue
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = len(text) if j < 0 else j + 2
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def _blocks(text: str) -> tuple[list[tuple[str, str]], str]:
    """Top-level ``<header> { <body> }`` blocks of a comment-free Nextflow config fragment, in
    order, and the text outside them (its statements); braces inside quotes are skipped."""
    blocks: list[tuple[str, str]] = []
    rest: list[str] = []
    depth, quote, start, head = 0, "", 0, 0
    for i, c in enumerate(text):
        if quote:
            if c == quote and text[i - 1] != "\\":
                quote = ""
        elif c in "'\"":
            quote = c
        elif c == "{":
            if depth == 0:
                seg = text[head:i]
                cut = max(seg.rfind("\n"), seg.rfind(";"), seg.rfind("}"))
                rest.append(seg[: cut + 1])
                blocks.append((seg[cut + 1 :].strip(), ""))
                start = i + 1
            depth += 1
        elif c == "}" and depth > 0:
            depth -= 1
            if depth == 0:
                blocks[-1] = (blocks[-1][0], text[start:i])
                head = i + 1
    rest.append(text[head:] if depth == 0 else "")
    return blocks, "".join(rest)


def _local_containers(
    local_tools: Path | None, runtime: str | None, directive: str = "container"
) -> dict[str, str]:
    """Registry tool name → ``container`` (or ``conda``) of the ``withName`` selectors in the
    ``profiles { <runtime> { process { … } } }`` block of ``local_tools.config``; empty without
    a file or a block for ``runtime``."""
    pattern = _CONDA if directive == "conda" else _CONTAINER
    if local_tools is None or runtime is None:
        return {}
    try:
        text = _strip_comments(local_tools.read_text(encoding="utf-8"))
    except OSError:
        return {}
    found: dict[str, str] = {}
    for top, profiles in _blocks(text)[0]:
        if top != "profiles":
            continue
        for name, body in _blocks(profiles)[0]:
            if name != runtime:
                continue
            for proc, selectors in _blocks(body)[0]:
                if proc != "process":
                    continue
                for selector, inner in _blocks(selectors)[0]:
                    m = _WITH_NAME.fullmatch(selector)
                    c = pattern.search(_blocks(inner)[1])
                    if m and c:
                        found[m.group(2).lower()] = c.group(2)
    return found


def _profiles(profile: str) -> set[str]:
    return {n.strip() for n in profile.split(",") if n.strip()}


def container_runtime(profile: str) -> Literal["docker", "apptainer"] | None:
    """The container runtime the comma-joined Nextflow ``profile`` selects (apptainer first);
    ``None`` without one."""
    names = _profiles(profile)
    if "apptainer" in names:
        return "apptainer"
    return "docker" if "docker" in names else None


def _step_images(images_json: Path) -> dict[str, str]:
    """Step-started tool → its ``ref`` in ``images.json`` (``{}`` when unreadable)."""
    try:
        data = json.loads(images_json.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    out: dict[str, str] = {}
    for name in _STEP_STARTED:
        entry = data.get(name)
        ref = entry.get("ref") if isinstance(entry, dict) else None
        if isinstance(ref, str) and ref:
            out[name] = ref
    return out


def _docker_image_exists(tag: str) -> bool:
    docker = shutil.which("docker")
    if docker is None:
        return False
    try:
        proc = subprocess.run(
            [docker, "image", "inspect", tag], capture_output=True, check=False, timeout=60
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0


def missing_tool_reasons(
    tools: Sequence[ToolSpec],
    profile: str,
    image_cache: Path,
    env_root: Path | None,
    local_tools: Path | None,
) -> dict[str, str]:
    """Why each of ``tools`` is not usable under the comma-joined Nextflow ``profile``.

    A tool pinned in the bundled ``containers.config`` is satisfied by its image: under docker
    Nextflow pulls it by digest; under apptainer it needs ``<image_cache>/<simpleName>.img``
    (ADR-0059). A conda-only tool needs the prebuilt env ``<env_root>/gippa-<tool>`` with a
    ``.gippa-built`` marker equal to the sha256 of ``workflow/envs/<tool>.yml``; under a container
    profile without conda (``docker`` or plain ``apptainer``) it is missing, since it would run on
    the host. A step-started tool (fgsea) instead counts there through its ``images.json`` ref:
    under docker that entry suffices (``docker run`` pulls it), under plain apptainer it needs
    ``<image_cache>/<simpleName>.img``. A local UniOP build of ``local_tools`` counts under its
    own runtime (:func:`uniop_local_problem`). ``test``/``stub`` need nothing; in-process tools
    are never missing.
    """
    names = _profiles(profile)
    if names & _NO_TOOL_PROFILES:
        return {}
    workflow = workflow_main().parent
    pinned = {
        process.lower(): ref
        for process, ref in read_containers_config(workflow / "conf" / "containers.config").items()
    }
    conda_on = bool(names & _CONDA_PROFILES)
    runtime = container_runtime(profile)
    local_uniop = (
        uniop_local_problem(profile, local_tools) if any(t.name == "uniop" for t in tools) else None
    )
    step_images = {} if conda_on else _step_images(workflow / "conf" / "images.json")
    setup_images = f"run: {_PKG} setup --images --profile apptainer"
    reasons: dict[str, str] = {}
    for tool in tools:
        if tool.in_process:
            continue
        ref = pinned.get(tool.name)
        if tool.name == "uniop" and local_uniop is not None:
            if local_uniop:
                reasons[tool.name] = local_uniop
        elif ref is not None:
            image = image_cache / _image_file_name(ref)
            if runtime == "apptainer" and not image.is_file():
                where = (
                    f"not in the apptainer cache {image_cache}"
                    if image_cache.is_dir()
                    else f"not found: the apptainer cache {image_cache} does not exist (set "
                    "NXF_APPTAINER_CACHEDIR or $GIPPA_ROOT, or give a DB directory)"
                )
                reasons[tool.name] = f"image {ref} {where}; {setup_images}"
        elif tool.name in _STEP_STARTED and not conda_on and runtime is not None:
            if _env_built(tool, env_root, workflow / "envs"):
                continue
            step_ref = step_images.get(tool.name)
            if step_ref is None:
                where = "$GIPPA_ROOT is not set" if env_root is None else "not built"
                reasons[tool.name] = (
                    f"prebuilt env {tool.env_name}: {where}; "
                    "workflow/conf/images.json pins no image for it"
                )
            elif runtime == "apptainer":
                image = image_cache / _image_file_name(step_ref)
                if not image.is_file():
                    reasons[tool.name] = (
                        f"image {step_ref} not in the apptainer cache {image_cache}; {setup_images}"
                    )
        elif runtime is not None and not conda_on and tool.name not in _STEP_STARTED:
            reasons[tool.name] = (
                f"conda-only tool under a container profile without conda: {tool.name}"
            )
        elif not _env_built(tool, env_root, workflow / "envs"):
            where = "$GIPPA_ROOT is not set" if env_root is None else "not built"
            reasons[tool.name] = f"prebuilt env {tool.env_name}: {where}"
    return dict(sorted(reasons.items()))


_BUILD_BLOCK = {"docker": "docker", "apptainer": "apptainer", "micromamba": "conda"}


def _uniop_build(db_dir: Path) -> dict[str, str] | None:
    """``local_build`` of the ``ok`` ``image:uniop`` item of ``<db_dir>/db_manifest.json``."""
    try:
        data = json.loads((db_dir / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    items = data.get("items") if isinstance(data, dict) else None
    for item in items if isinstance(items, list) else []:
        if isinstance(item, dict) and item.get("name") == UNIOP_SETUP_ITEM:
            build = item.get("local_build") if item.get("status") == "ok" else None
            return {str(k): str(v) for k, v in build.items()} if isinstance(build, dict) else None
    return None


def _docker_image_id(tag: str) -> str | None:
    docker = shutil.which("docker")
    if docker is None:
        return None
    try:
        proc = subprocess.run([docker, "image", "inspect", "--format", "{{.Id}}", tag],
                              capture_output=True, text=True, check=False, timeout=60)  # fmt: skip
    except (OSError, subprocess.TimeoutExpired):
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def _uniop_commit_problem(prefix: Path) -> str:
    want = SOURCES["uniop_source"].version.removeprefix("git:")
    try:
        have = (prefix / "share" / "uniop" / "GIPPA_COMMIT").read_text(encoding="utf-8").strip()
    except OSError:
        have = ""
    return "" if have == want else f"{prefix}/share/uniop/GIPPA_COMMIT is not the pin {want[:8]}"


def _artifact_problem(runtime: str, ref: str, digest: str | None) -> str:
    """Why the built image/env ``ref`` is not usable (``digest`` ``None``: only its presence)."""
    if runtime == "micromamba":
        return _uniop_commit_problem(Path(ref))
    if runtime == "docker":
        if digest is None:
            ok = _docker_image_exists(ref)
        else:
            image_id = _docker_image_id(ref)
            ok = image_id is not None
            if ok and image_id != digest:
                return f"local image {ref} is {image_id}, not the recorded build {digest}"
        return "" if ok else f"local image {ref} not found (docker image inspect)"
    if not Path(ref).is_file():
        return f"local image {ref} not found"
    sha = _sha256(Path(ref))
    if digest is not None and f"sha256:{sha}" != digest:
        return f"local image {ref} differs from the recorded build (sha256 {digest})"
    return ""


def uniop_local_problem(profile: str, local_tools: Path | None) -> str | None:
    """The local UniOP build for ``profile`` (spec item 3): ``""`` usable, else why not; ``None``
    when none applies (no file, no block for ``R(P)``, or a ``farmshare*`` list, where
    ``prebuilt_envs.config`` owns UNIOP) and the generic rows decide. The block must match the
    ``local_build`` of ``db_manifest.json`` next to the file (runtime, ref, ``key("uniop")``,
    digest); without a manifest only the block's image or env is checked (14f1)."""
    runtime = uniop_runtime(profile)
    if local_tools is None or runtime is None or _profiles(profile) & _FARMSHARE_PROFILES:
        return None
    directive = "conda" if runtime == "micromamba" else "container"
    ref = _local_containers(local_tools, _BUILD_BLOCK[runtime], directive).get("uniop")
    if ref is None:
        return None
    if runtime != "docker" and not Path(ref).is_absolute():
        return f"{LOCAL_TOOLS} {_BUILD_BLOCK[runtime]} path '{ref}' is not absolute"
    if not (local_tools.parent / MANIFEST).is_file():
        return _artifact_problem(runtime, ref, None)
    build = _uniop_build(local_tools.parent)
    if build is None:
        return f"no local UniOP build recorded in {MANIFEST}"
    if build.get("runtime") != runtime:
        return f"the recorded local build is for {build.get('runtime')}, not {runtime}"
    if build.get("ref") != ref:
        return f"{LOCAL_TOOLS} names {ref}, not the recorded build {build.get('ref')}"
    want = image_key("uniop")
    if build.get("key") != want:
        return f"local build {build.get('key')} does not match this release ({want})"
    return _artifact_problem(runtime, ref, None if runtime == "micromamba" else build.get("digest"))


def missing_tools(
    tools: Sequence[ToolSpec],
    profile: str,
    image_cache: Path,
    env_root: Path | None,
    local_tools: Path | None,
) -> list[str]:
    """Sorted names of ``tools`` not usable under ``profile`` (:func:`missing_tool_reasons`)."""
    return sorted(missing_tool_reasons(tools, profile, image_cache, env_root, local_tools))


# ---------------------------------------------------------------- UniOP (M11-14d)

UniopRuntime = Literal["micromamba", "apptainer", "docker"]
UniopLevel = Literal["ok", "warn", "error"]
#: UniOP's entry point in an image (``modules/local/uniop/main.nf``; the image installs it under
#: ``/opt/conda/share/uniop``)
_UNIOP_HELP = 'python "${UNIOP_HOME:-${CONDA_PREFIX:-/opt/conda}/share/uniop}/src/UniOP" --help'
PROBE_TIMEOUT_S = 300
#: spec item 1 (c): the warning and error texts (``{cmd}``: :func:`uniop_setup_command`)
UNIOP_AUTO_WARNING = (
    "uniop not installed (unlicensed, not distributed); operons use the distance model; run: {cmd}"
)
UNIOP_AUTO_PUBLICATION_ERROR = (
    "uniop not installed; set annotation.uniop.enabled: false to publish with the "
    "distance-model operons (see Methods)"
)
UNIOP_FALSE_WARNING = (
    "annotation.uniop.enabled is false; operons use the distance model (see Methods)"
)
UNIOP_TRUE_ERROR = "uniop not installed; run: {cmd}"
UNIOP_SETUP_ITEM = "image:uniop"


@dataclass(frozen=True, slots=True)
class UniopResolution:
    """How ``annotation.uniop.enabled`` resolved (:func:`resolve_uniop`).

    ``cfg`` is the configuration ``resolved_config.yaml`` records: ``enabled`` is a bool, and
    ``resolved_from`` is ``"auto"`` when it came from ``auto``. ``message`` is the
    ``environment.tools`` warning or error text, ``""`` when ``level`` is ``ok``.
    """

    run_uniop: bool
    level: UniopLevel
    message: str
    cfg: GippaConfig


def uniop_setup_command(profile: str | None = None, db_dir: str | None = None) -> str:
    """The command that installs UniOP for ``profile`` into ``db_dir`` (spec item 3)."""
    return (
        f"{_PKG} setup --uniop --profile {profile or '<profile>'} "
        f"--db-dir {db_dir or '<db_dir>'} --accept-uniop-terms"
    )


def _with_detail(message: str, detail: str) -> str:
    return f"{message}; {detail}" if detail else message


def _set_enabled(
    cfg: GippaConfig, enabled: bool, resolved_from: Literal["auto"] | None
) -> GippaConfig:
    uniop = cfg.annotation.uniop.model_copy(
        update={"enabled": enabled, "resolved_from": resolved_from}
    )
    ann = cfg.annotation.model_copy(update={"uniop": uniop})
    return cfg.model_copy(update={"annotation": ann})


def resolve_uniop(
    cfg: GippaConfig,
    *,
    installed: bool,
    publication: bool,
    detail: str = "",
    setup_command: str | None = None,
) -> UniopResolution:
    """Resolve ``annotation.uniop.enabled`` (spec item 1 (c) table; ADR-0061 decision 5). Pure.

    ``true``: UniOP runs, or ``environment.tools`` errors (exit 2) when it is not installed.
    ``auto``: resolves to ``true`` when installed; else to ``false`` with a warning, or an error
    under ``publication``. ``false``: the distance model, always with a warning (accepted under
    ``publication``). ``detail`` (why UniOP is missing) is appended to the error or the warning
    when UniOP is not installed. A precomputed operon layer needs no UniOP: ``ok``, UniOP does
    not run. A configuration resolved before (``resolved_from: auto``) resolves the same way
    again.
    """
    uni = cfg.annotation.uniop
    cmd = setup_command or uniop_setup_command()
    if cfg.annotation.precomputed.operon_pairs is not None:
        on = uni.enabled is True
        return UniopResolution(False, "ok", "", _set_enabled(cfg, on, uni.resolved_from))
    auto = uni.enabled == "auto" or uni.resolved_from == "auto"
    if uni.enabled is False and not auto:
        return UniopResolution(False, "warn", UNIOP_FALSE_WARNING, cfg)
    if auto and uni.enabled is False:
        installed = False  # already resolved: keep the recorded outcome
    if installed:
        return UniopResolution(True, "ok", "", _set_enabled(cfg, True, "auto" if auto else None))
    if not auto:
        return UniopResolution(
            False, "error", _with_detail(UNIOP_TRUE_ERROR.format(cmd=cmd), detail), cfg
        )
    out = _set_enabled(cfg, False, "auto")
    if publication:
        return UniopResolution(
            False, "error", _with_detail(UNIOP_AUTO_PUBLICATION_ERROR, detail), out
        )
    return UniopResolution(
        False, "warn", _with_detail(UNIOP_AUTO_WARNING.format(cmd=cmd), detail), out
    )


_FARMSHARE_PROFILES = frozenset({"farmshare", "farmshare_local"})
_CONTAINER_PROFILES = frozenset({"docker", "apptainer"})


def _uniop_ambiguous(names: set[str]) -> bool:
    """Plain ``conda`` combined with a container profile (no ``farmshare*``): not supported."""
    return (
        "conda" in names and bool(names & _CONTAINER_PROFILES) and not names & _FARMSHARE_PROFILES
    )


def uniop_runtime(profile: str) -> UniopRuntime | None:
    """``R(P)`` (spec item 1b): micromamba when a profile enables conda, else apptainer, else
    docker; ``None`` without a container or conda profile, and for an ambiguous list (``conda``
    with ``docker`` or ``apptainer`` and no ``farmshare*``; :func:`uniop_runtime_problem`)."""
    names = _profiles(profile)
    if _uniop_ambiguous(names):
        return None
    if names & _CONDA_PROFILES:
        return "micromamba"
    if "apptainer" in names:
        return "apptainer"
    return "docker" if "docker" in names else None


def uniop_runtime_problem(profile: str) -> str:
    """Why ``R(P)`` is ``None`` (``""`` when it resolves): an ambiguous list, or no container or
    conda profile."""
    names = _profiles(profile)
    if _uniop_ambiguous(names):
        return f"ambiguous UniOP runtime for {profile}: use one of conda | docker | apptainer"
    return "" if uniop_runtime(profile) is not None else "no container or conda profile"


def local_tools_applies(profile: str) -> bool:
    """``gippa run`` adds ``-c <db_dir>/local_tools.config`` only when ``R(P)`` is docker,
    apptainer, or micromamba under plain ``conda``; never under ``farmshare*``, where
    ``prebuilt_envs.config`` owns UNIOP (spec item 1b)."""
    runtime = uniop_runtime(profile)
    if runtime == "micromamba":
        return not _profiles(profile) & _FARMSHARE_PROFILES
    return runtime is not None


def _uniop_argv(
    profile: str, image_cache: Path, env_root: Path | None, local_tools: Path | None
) -> list[str] | str:
    """The ``--help`` command line of UniOP under ``profile``, or why there is none."""
    runtime = uniop_runtime(profile)
    if runtime == "micromamba":
        local = "" if _profiles(profile) & _FARMSHARE_PROFILES else "conda"
        block = _local_containers(local_tools, local or None, "conda").get("uniop")
        if env_root is None and block is None:
            return "$GIPPA_ROOT is not set"
        prefix = Path(block) if block else Path(str(env_root)) / TOOLS["uniop"].env_name
        home = os.environ.get("UNIOP_HOME") or str(prefix / "share" / "uniop")
        return [str(prefix / "bin" / "python"), f"{home}/src/UniOP", "--help"]
    if runtime is None:
        return uniop_runtime_problem(profile)
    built = _local_containers(local_tools, runtime).get("uniop")
    if built is None:
        return f"no local UniOP image in {LOCAL_TOOLS}"
    if runtime == "docker":
        docker = shutil.which("docker")
        if docker is None:
            return "docker not found on PATH"
        return [docker, "run", "--rm", built, "sh", "-c", _UNIOP_HELP]
    exe = shutil.which("apptainer") or shutil.which("singularity")
    if exe is None:
        return "apptainer not found on PATH"
    sif = Path(built)
    if not sif.is_absolute():
        sif = image_cache / _image_file_name(built)
    return [exe, "exec", "--cleanenv", str(sif), "sh", "-c", _UNIOP_HELP]


def probe_uniop(
    profile: str, *, image_cache: Path, env_root: Path | None, local_tools: Path | None
) -> tuple[int, str]:
    """Run UniOP's entry point with ``--help`` in the runtime ``profile`` resolves to
    (:func:`uniop_runtime`): the prebuilt env, the local image of ``local_tools`` under docker,
    or its ``.sif`` under apptainer. Returns ``(returncode, stderr)``; 127 when it cannot be
    started, 124 on a timeout."""
    argv = _uniop_argv(profile, image_cache, env_root, local_tools)
    if isinstance(argv, str):
        return 127, argv
    try:
        proc = subprocess.run(
            argv, capture_output=True, text=True, check=False, timeout=PROBE_TIMEOUT_S
        )
    except subprocess.TimeoutExpired:
        return 124, f"UniOP --help timed out after {PROBE_TIMEOUT_S} s"
    except OSError as exc:
        return 127, f"{type(exc).__name__}: {exc}"
    return proc.returncode, proc.stderr or ""


def uniop_setup_record(db_dir: str | Path | None) -> str:
    """``(last setup: <message>)`` of a ``db_manifest.json`` ``image:uniop`` item whose
    ``status`` is not ``ok``; ``""`` without one (or without a readable manifest)."""
    if not db_dir:
        return ""
    try:
        data = json.loads((Path(db_dir) / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    items = data.get("items") if isinstance(data, dict) else None
    for item in items if isinstance(items, list) else []:
        if (
            isinstance(item, dict)
            and item.get("name") == UNIOP_SETUP_ITEM
            and item.get("status") != "ok"
        ):
            message = str(item.get("message") or item.get("status") or "").strip()
            return f"(last setup: {message})"
    return ""


def _slurm_mem_bytes() -> int | None:
    """Slurm job memory: ``SLURM_MEM_PER_NODE`` or ``SLURM_MEM_PER_CPU`` times CPUs (MB)."""
    per_node = os.environ.get("SLURM_MEM_PER_NODE", "")
    if per_node.isdigit() and int(per_node) > 0:
        return int(per_node) * _MB
    per_cpu = os.environ.get("SLURM_MEM_PER_CPU", "")
    cpus = os.environ.get("SLURM_CPUS_ON_NODE", "1")
    if per_cpu.isdigit() and int(per_cpu) > 0 and cpus.isdigit():
        return int(per_cpu) * max(int(cpus), 1) * _MB
    return None


def _cgroup_mem_bytes() -> int | None:
    for path in _CGROUP_LIMITS:
        try:
            text = path.read_text(encoding="ascii").strip()
        except (OSError, ValueError):
            continue
        if text.isdigit() and 0 < int(text) < 2**60:  # v1 reports "no limit" as a huge number
            return int(text)
    return None


def _total_ram_bytes() -> int:
    return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))


def probe_resources(outdir: Path) -> tuple[int, int]:
    """``(free_disk_bytes, ram_bytes)`` probed live: free disk where ``outdir`` is (or will be)
    and the Slurm or cgroup memory limit, else the machine's total RAM."""
    where = Path(outdir).absolute()
    while not where.exists() and where != where.parent:
        where = where.parent
    disk = shutil.disk_usage(where).free
    limits = [m for m in (_slurm_mem_bytes(), _cgroup_mem_bytes()) if m is not None]
    return int(disk), min(limits) if limits else _total_ram_bytes()


# ---------------------------------------------------------------- rules


# (name, kind, pinned version) per module: databases and models from the `gippa setup` manifest;
# binaries and environments come from the tool registry (:func:`missing_tools`, run-wiring
# decision 4a). A database matches the pin to its written precision (§3.2, §7.1.1).
_ROWS: dict[str, tuple[tuple[str, Kind, str | None], ...]] = {
    "core": (("bakta_db", "db", "6"), ("promotech", "model", None)),
    "mobile": (("genomad_db", "db", "1.9"),),
    "bgc": (("antismash_db", "db", None),),
    "function": (("go", "db", None), ("kegg", "db", None), ("cog", "db", None)),
}


def requirements(modules: Sequence[str], promoters_mode: str) -> list[Requirement]:
    """Requirements of ``modules`` (core always included), in module then table order."""
    if promoters_mode not in PROMOTER_MODES:
        raise ValueError(f"unknown promoters mode {promoters_mode!r}")
    unknown = sorted(set(modules) - set(_ROWS))
    if unknown:
        raise ValueError(f"unknown modules: {', '.join(unknown)}")
    return [
        Requirement(name, module, kind, version)
        for module in _ROWS
        if module == "core" or module in modules
        for name, kind, version in _ROWS[module]
    ]


def resource_needs(promoters_mode: str) -> tuple[int, int]:
    """``(need_disk_bytes, need_ram_bytes)``; only genome-wide Promotech states a need."""
    return 0, GENOME_PROMOTECH_RAM_BYTES if promoters_mode == "genome" else 0


_VERSION = re.compile(r"\s*v?(\d+(?:\.\d+)*)(.*)")


def _version(v: str) -> tuple[tuple[int, ...], bool] | None:
    """``(release numbers, is_final)``; any suffix but ``+local`` marks a pre-release."""
    m = _VERSION.fullmatch(v)
    if m is None:
        return None
    suffix = m.group(2).strip()
    return tuple(int(x) for x in m.group(1).split(".")), suffix == "" or suffix.startswith("+")


def _pad(t: tuple[int, ...], n: int = 3) -> tuple[int, ...]:
    return t + (0,) * (n - len(t))


def _version_problem(req: Requirement, found: str) -> str | None:
    """Why ``found`` does not satisfy the pinned ``req.min_version``, or ``None``.

    Tools: same major.minor as the pin, patch >= the pin's, a pre-release ranking below its
    release. Databases: equal to the pin over the components the pin is written with.
    """
    pin = req.min_version
    if pin is None:
        return None
    need, have = _version(pin), _version(found)
    if need is None:
        raise ValueError(f"unparseable pinned version {pin!r} for {req.name}")
    n = need[0]
    if req.kind == "db":
        ok = have is not None and have[0][: len(n)] == n
        return None if ok else f"{found} (need {pin}.x)"
    ok = (
        have is not None
        and _pad(have[0])[:2] == _pad(n)[:2]
        and (_pad(have[0]), have[1]) >= (_pad(n), need[1])
    )
    n2 = _pad(n)
    return None if ok else f"{found} (need {n2[0]}.{n2[1]}.x >= {pin})"


def _gib(n: int) -> int | float:
    g = float(f"{n / GIB:.6g}")
    return int(g) if g.is_integer() else g


def _reported_tools(reqs: Sequence[Requirement], probe: EnvProbe) -> list[Requirement]:
    """Registry tools of the enabled modules (those of ``reqs``, plus mapping) an injected probe
    reports (``None``: missing), pinned at their registry version (git pins carry no comparable
    version); databases come from ``reqs``."""
    names = {r.name for r in reqs}
    modules = {r.module for r in reqs} | {"core", "mapping"}
    out = []
    for name in probe.found:
        tool = TOOLS.get(name)
        if tool is None or tool.in_process or name in names or tool.module not in modules:
            continue
        pin = None if _version(tool.version) is None else tool.version
        out.append(Requirement(name, tool.module, "binary", pin))
    return out


def _tools_check(
    reqs: Sequence[Requirement],
    probe: EnvProbe,
    tools_missing: Sequence[str] = (),
    tool_reasons: Mapping[str, str] | None = None,
) -> Check:
    rule = TOOLS_RULE
    missing: list[str] = list(tools_missing)
    wrong: dict[str, JsonValue] = {}
    found: dict[str, JsonValue] = {}
    for r in [*reqs, *_reported_tools(reqs, probe)]:
        v = probe.found.get(r.name)
        if v is None:
            missing.append(r.name)
            continue
        found[r.name] = v
        problem = _version_problem(r, v)
        if problem is not None:
            wrong[r.name] = problem
    if not missing and not wrong:
        return Check("environment.tools", "environment", "ok", rule, {"found": found})
    value: dict[str, JsonValue] = {}
    parts = []
    if missing:
        value["missing"] = list[JsonValue](sorted(set(missing)))
        parts.append("missing: " + ", ".join(sorted(set(missing))))
        parts += [f"{k}: {v}" for k, v in sorted((tool_reasons or {}).items())]
    if wrong:
        value["wrong_version"] = dict(sorted(wrong.items()))
        parts.append("wrong version: " + ", ".join(f"{k} {v}" for k, v in sorted(wrong.items())))
    return Check("environment.tools", "environment", "error", rule, value, "; ".join(parts), 2)


def _resources_check(probe: EnvProbe, need_disk_bytes: int, need_ram_bytes: int) -> Check:
    """Verdict only: needs are recorded, the probed numbers are not (machine-independent)."""
    rule = RESOURCES_RULE
    value: dict[str, JsonValue] = {}
    if need_disk_bytes > 0:
        value["need_disk_gb"] = _gib(need_disk_bytes)
    if need_ram_bytes <= 0:
        return Check("environment.resources", "environment", "info", rule, value or None)
    value["need_ram_gb"] = _gib(need_ram_bytes)
    enough = probe.total_ram_bytes >= need_ram_bytes
    value["ram_sufficient"] = enough
    level: Level = "ok" if enough else "warn"
    message = "" if enough else f"RAM < {value['need_ram_gb']} GiB needed by genome-wide Promotech"
    return Check("environment.resources", "environment", level, rule, value, message)


def environment_checks(
    reqs: Sequence[Requirement],
    probe: EnvProbe,
    need_disk_bytes: int,
    need_ram_bytes: int,
    tools_missing: Sequence[str] = (),
    tool_reasons: Mapping[str, str] | None = None,
) -> list[Check]:
    """``environment.tools`` and ``environment.resources`` for ``reqs`` against ``probe``;
    ``tools_missing`` (from :func:`missing_tools`) are reported missing alongside databases,
    each with its ``tool_reasons`` entry (:func:`missing_tool_reasons`) in the message.

    Free disk never warns (ADR-0005 item 4): ``need_disk_bytes`` is only recorded.
    """
    return [
        _tools_check(reqs, probe, tools_missing, tool_reasons),
        _resources_check(probe, need_disk_bytes, need_ram_bytes),
    ]
