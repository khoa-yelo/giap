"""Environment preflight rules (GIAP_SPEC.md §7.0 Environment rows; ADR-0005 items 3-4).

``environment.tools``: every binary, environment, database and model of the enabled modules is
found at a usable version → else error (exit 2). ``environment.resources``: the only warning is
RAM below the 24 GiB that genome-wide Promotech needs (the caller passes it as
``need_ram_bytes``, see :func:`resource_needs`). The check holds the verdict only, never the raw
RAM and disk numbers, so ``preflight.json`` is identical across machines (§6.6); the caller logs
the raw numbers (§6.7). Rules take an injected :class:`EnvProbe`; the real sources are
read by :func:`read_manifest` (``$GIAP_DB_DIR/db_manifest.json``) and :func:`probe_resources`.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from giap.preflight.checks import Check, JsonValue, Level

__all__ = [
    "GENOME_PROMOTECH_RAM_BYTES",
    "GIB",
    "MANIFEST",
    "RESOURCES_RULE",
    "TOOLS_RULE",
    "EnvProbe",
    "ProbeError",
    "Requirement",
    "environment_checks",
    "probe_resources",
    "read_manifest",
    "requirements",
    "resource_needs",
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

    Provisional format (pinned by M2's ``giap setup``): ``{"<name>": {"version": "<v>",
    "sha256": "<hex>", ...}}``; a bare ``"<name>": "<v>"`` is accepted too. Raises
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
    found: dict[str, str | None] = {}
    for name, entry in data.items():
        v = entry.get("version") if isinstance(entry, dict) else entry
        found[name] = v if isinstance(v, str) else None
    return found


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


# (name, kind, pinned version) per module; versions pinned in §3.2 and §7.1.1. A tool matches its
# pin's major.minor with patch >= the pin's; a database matches the pin to its written precision.
_ROWS: dict[str, tuple[tuple[str, Kind, str | None], ...]] = {
    "core": (
        ("bakta", "binary", "1.12.2"),
        ("bakta_db", "db", "6"),
        ("miniprot", "binary", "0.18"),
        ("tblastn", "binary", None),
        ("uniop", "env", None),
        ("promotech", "env", None),
        ("promotech_models", "model", None),
        ("promoter_calculator", "env", "1.2.4"),
        ("terminatornet", "env", None),
        ("transtermhp", "binary", "2.09"),
        ("ostir", "env", "1.1.3"),
    ),
    "mobile": (("genomad", "binary", "1.12.0"), ("genomad_db", "db", "1.9")),
    "bgc": (("antismash", "binary", "8.0.4"), ("antismash_db", "db", None)),
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


def _tools_check(reqs: Sequence[Requirement], probe: EnvProbe) -> Check:
    rule = TOOLS_RULE
    missing: list[str] = []
    wrong: dict[str, JsonValue] = {}
    found: dict[str, JsonValue] = {}
    for r in reqs:
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
        value["missing"] = list[JsonValue](sorted(missing))
        parts.append("missing: " + ", ".join(sorted(missing)))
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
    reqs: Sequence[Requirement], probe: EnvProbe, need_disk_bytes: int, need_ram_bytes: int
) -> list[Check]:
    """``environment.tools`` and ``environment.resources`` for ``reqs`` against ``probe``.

    Free disk never warns (ADR-0005 item 4): ``need_disk_bytes`` is only recorded.
    """
    return [_tools_check(reqs, probe), _resources_check(probe, need_disk_bytes, need_ram_bytes)]
