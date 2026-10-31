"""Nextflow profile registry, combination rules and auto-detection (spec M11-15 item 1, item 3).

``detect_profile`` is pure: it reads an injected :class:`Probe`. ``probe_host`` is the only part
that looks at the host (``shutil.which``, ``docker info``, ``scontrol show config`` with 10 s
timeouts, ``<db_dir>/db_manifest.json``).
"""

from __future__ import annotations

import functools
import json
import re
import shutil
import subprocess
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

__all__ = ["FARMSHARE_CLUSTER", "R5_MESSAGE", "Detected", "Probe", "Site", "accepted_profiles",
           "check_profile", "detect_profile", "parse_profiles", "probe_host",
           "read_sites", "uses_local_executor"]  # fmt: skip

#: D1: the Slurm cluster name of Stanford FarmShare (pinned; spec M11-15 item 1 AC 5).
FARMSHARE_CLUSTER = "farmshare"
INSTALL_DOCS = "https://khoa-yelo.github.io/gippa/get-started/install/"
R5_MESSAGE = (
    "no container runtime or conda found: install Apptainer (no root needed on most Linux "
    f"systems) or Docker; see {INSTALL_DOCS}"
)
HINT_DOCKER_ONLY = (
    "sbatch on PATH but only docker found: tasks run on this host; install Apptainer or pass "
    "--profile"
)
HINT_LOGIN_NODE = (
    "the Nextflow head process runs on this host for the whole run; submit gippa run as a job "
    "(see Install > Run on a cluster)"
)
IMAGES_SKIPPED = (
    "images: no container runtime found here; rerun gippa setup --images --profile "
    "apptainer|singularity|docker on a host that has one (same --db-dir)"
)
PROBED = ("sbatch", "apptainer", "singularity", "docker", "conda", "mamba", "micromamba")
CONTAINER_RUNTIMES = ("apptainer", "singularity", "docker")
_RUNTIMES = frozenset({*CONTAINER_RUNTIMES, "conda"})
_CONDA = frozenset({"conda", "mamba", "micromamba"})
_REJECTED = frozenset({"test", "stub"})
_LOCAL_EXECUTORS = frozenset({"standard", "farmshare_local"})
_RUNTIME_ORDER = ("apptainer", "singularity", "conda")  # R1, R2, R4 (R3 docker first off Slurm)
_TIMEOUT = 10


@dataclass(frozen=True)
class Probe:
    """What the launching host offers; ``Probe()`` is "nothing found"."""

    on_path: frozenset[str] = frozenset()
    docker_ok: bool = False
    slurm_job_id: str | None = None
    cluster_name: str | None = None
    setup_runtime: str | None = None


@dataclass(frozen=True)
class Detected:
    """The chosen profile (``None``: R5, no runtime), the reason printed after ``auto:``, the
    hints of edge cases (a)/(b) and warnings (an unusable recorded runtime)."""

    profile: str | None
    reason: str
    hints: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class Site:
    """A site profile's header: ``// site:``, ``// detect: cluster_name=``, ``// runtime:``."""

    name: str
    cluster_name: str | None
    runtime: str


def _usable(probe: Probe) -> set[str]:
    out: set[str] = {n for n in ("apptainer", "singularity") if n in probe.on_path}
    if "docker" in probe.on_path and probe.docker_ok:
        out.add("docker")
    if probe.on_path & _CONDA:
        out.add("conda")
    return out


def detect_profile(probe: Probe, sites: Iterable[Site] = ()) -> Detected:
    """Executor D1-D4 and runtime R0-R5 of spec M11-15 item 1; first match wins."""
    if probe.cluster_name == FARMSHARE_CLUSTER:
        return Detected("farmshare,apptainer", "FarmShare detected")
    for site in sites:
        if site.cluster_name is not None and site.cluster_name == probe.cluster_name:
            return Detected(f"{site.name},{site.runtime}", f"site {site.name} detected")
    slurm = "sbatch" in probe.on_path
    usable = _usable(probe)
    order = _RUNTIME_ORDER if slurm else ("docker", *_RUNTIME_ORDER)
    warnings: list[str] = []
    recorded = probe.setup_runtime
    runtime = None
    if recorded is not None:
        if recorded in usable and recorded in order:
            runtime = recorded
        elif recorded not in usable:
            warnings.append(f"setup used {recorded}, not available here")
    via_setup = runtime is not None
    runtime = runtime or next((r for r in order if r in usable), None)
    if runtime is None:
        if slurm and "docker" in usable:  # edge case (a): docker cannot run Slurm tasks
            return Detected("docker", "only docker usable", (HINT_DOCKER_ONLY,), tuple(warnings))
        return Detected(None, "", (), tuple(warnings))
    why = f"{runtime} as recorded by gippa setup" if via_setup else f"{runtime} found"
    if not slurm:
        return Detected(runtime, why, (), tuple(warnings))
    hints = (HINT_LOGIN_NODE,) if probe.slurm_job_id is None else ()
    reason = f"sbatch on PATH; {why}" if via_setup else "sbatch on PATH"
    return Detected(f"slurm,{runtime}", reason, hints, tuple(warnings))


def _run(argv: list[str]) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(argv, capture_output=True, text=True, timeout=_TIMEOUT,
                              stdin=subprocess.DEVNULL, check=False)  # fmt: skip
    except (OSError, subprocess.SubprocessError):
        return None


def _cluster_name(environ: Mapping[str, str]) -> str | None:
    if environ.get("SLURM_CLUSTER_NAME"):
        return environ["SLURM_CLUSTER_NAME"]
    scontrol = shutil.which("scontrol", path=environ.get("PATH"))
    res = _run([scontrol, "show", "config"]) if scontrol else None
    if res is None or res.returncode != 0:
        return None
    m = re.search(r"(?m)^\s*ClusterName\s*=\s*(\S+)", res.stdout or "")
    return m.group(1) if m else None


def _setup_runtime(db_dir: Path | None) -> str | None:
    if db_dir is None:
        return None
    try:
        doc = json.loads((db_dir / "db_manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    setup = doc.get("setup") if isinstance(doc, dict) else None
    runtime = setup.get("runtime") if isinstance(setup, dict) else None
    return runtime if runtime in CONTAINER_RUNTIMES else None


def probe_host(environ: Mapping[str, str], db_dir: Path | None) -> Probe:
    """The real :class:`Probe`: PATH lookups, ``docker info``, ``environ``'s Slurm variables, the
    cluster name and ``setup.runtime`` of ``db_dir``. A timeout or ``OSError`` is "not found"."""
    found = {n: p for n in PROBED if (p := shutil.which(n, path=environ.get("PATH")))}
    docker = found.get("docker")
    res = _run([docker, "info"]) if docker else None
    return Probe(
        on_path=frozenset(found),
        docker_ok=res is not None and res.returncode == 0,
        slurm_job_id=environ.get("SLURM_JOB_ID") or None,
        cluster_name=_cluster_name(environ),
        setup_runtime=_setup_runtime(db_dir),
    )


_HEADER = re.compile(r"(?m)^//\s*(site|runtime|detect):\s*(\S+)\s*$")


def read_sites(sites_dir: Path) -> tuple[Site, ...]:
    """The headers of ``sites_dir/*.config`` (``_``-prefixed files are not sites), by name."""
    if not sites_dir.is_dir():
        return ()
    out: list[Site] = []
    for path in sorted(sites_dir.glob("*.config")):
        if path.name.startswith("_"):
            continue
        head = dict(_HEADER.findall(path.read_text(encoding="utf-8")))
        detect = head.get("detect", "")
        cluster = (
            detect.removeprefix("cluster_name=") if detect.startswith("cluster_name=") else None
        )
        out.append(Site(head.get("site", path.stem), cluster or None, head.get("runtime", "")))
    return tuple(out)


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)


def parse_profiles(text: str) -> tuple[str, ...]:
    """The names of the ``profiles { }`` block of a ``nextflow.config`` text, in file order."""
    text = _strip_comments(text)
    m = re.search(r"(?m)^\s*profiles\s*\{", text)
    if m is None:
        return ()
    names: list[str] = []
    depth = 0
    for i in range(text.index("{", m.start()), len(text)):
        ch = text[i]
        if ch == "{":
            depth += 1
            if depth == 2 and (name := re.search(r"(\w+)\s*$", text[:i])):
                names.append(name.group(1))
        elif ch == "}":
            depth -= 1
            if depth == 0:
                break
    return tuple(names)


@functools.cache
def accepted_profiles() -> tuple[str, ...]:
    """The bundled ``workflow/nextflow.config`` profiles a user may pass (not ``test``/``stub``)."""
    from gippa.runner.workflow_path import workflow_main

    try:
        text = (workflow_main().parent / "nextflow.config").read_text(encoding="utf-8")
    except OSError:  # no bundled workflow: nothing to launch, so no profile is accepted
        return ()
    return tuple(n for n in parse_profiles(text) if n not in _REJECTED)


def uses_local_executor(profile: str) -> bool:
    """Nextflow's local executor runs the tasks of the resolved comma-joined ``profile``: every
    name is a runtime, ``standard`` or ``farmshare_local``. ``farmshare``, ``slurm`` and site or
    unknown profiles submit their tasks as separate jobs (#514)."""
    local = _RUNTIMES | _LOCAL_EXECUTORS
    return all(name.strip() in local for name in profile.split(","))


def check_profile(value: str) -> str:
    """``value`` if its comma-joined names are known, with at most one executor profile and at
    most one container runtime; else ``ValueError`` with the message."""
    names = [n.strip() for n in value.split(",")]
    known = accepted_profiles()
    bad = [n for n in names if n not in known]
    if bad:
        raise ValueError(
            f"unknown profile(s) {', '.join(repr(n) for n in bad)} in {value!r}; "
            f"choose comma-joined from: {', '.join(known)}"
        )
    executors = [n for n in names if n not in _RUNTIMES]
    if len(executors) > 1:
        raise ValueError(f"more than one executor profile: {', '.join(executors)}")
    runtimes = [n for n in names if n in CONTAINER_RUNTIMES]
    if len(runtimes) > 1:
        raise ValueError(f"more than one container runtime: {', '.join(runtimes)}")
    return ",".join(names)  # normalised: "slurm, apptainer" reaches Nextflow without the space
