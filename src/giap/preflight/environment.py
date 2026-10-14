"""Environment preflight rules (GIAP_SPEC.md §7.0 Environment rows; ADR-0005 items 3-4).

``environment.tools``: every binary, environment, database and model of the enabled modules is
found at a usable version → else error (exit 2). ``environment.resources``: free disk and RAM are
recorded; the only warning is RAM below the 24 GiB that genome-wide Promotech needs (the caller
passes it as ``need_ram_bytes``, see :func:`resource_needs`). The probe is injected.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from giap.preflight.checks import Check, JsonValue, Level

__all__ = [
    "GENOME_PROMOTECH_RAM_BYTES",
    "GIB",
    "EnvProbe",
    "Requirement",
    "environment_checks",
    "requirements",
    "resource_needs",
]

GIB = 1024**3
GENOME_PROMOTECH_RAM_BYTES = 24 * GIB
PROMOTER_MODES = ("targeted", "genome")

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


# (name, kind, minimum version) per module; versions pinned in §3.2 and §7.1.1.
_ROWS: dict[str, tuple[tuple[str, Kind, str | None], ...]] = {
    "core": (
        ("bakta", "binary", "1.12"),
        ("bakta_db", "db", "6.0"),
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


def _version(v: str) -> tuple[int, ...]:
    return tuple(int(x) for x in re.findall(r"\d+", v))


def _version_problem(req: Requirement, found: str) -> str | None:
    """Why ``found`` does not satisfy ``req`` (databases need the same major version)."""
    if req.min_version is None:
        return None
    need, have = _version(req.min_version), _version(found)
    if req.kind == "db":
        return None if have[:1] == need[:1] else f"{found} (need {need[0]}.x)"
    return None if have >= need else f"{found} (need >= {req.min_version})"


def _gib(n: int) -> int | float:
    g = float(f"{n / GIB:.6g}")
    return int(g) if g.is_integer() else g


def _tools_check(reqs: Sequence[Requirement], probe: EnvProbe) -> Check:
    rule = "binaries, environments and DB versions of enabled modules found → else error (exit 2)"
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
    rule = "RAM ≥ 24 GiB for genome-wide Promotech → else warn; free disk and RAM recorded"
    value: dict[str, JsonValue] = {
        "ram_gb": _gib(probe.total_ram_bytes),
        "free_disk_gb": _gib(probe.free_disk_bytes),
    }
    if need_disk_bytes > 0:
        value["need_disk_gb"] = _gib(need_disk_bytes)
    if need_ram_bytes <= 0:
        return Check("environment.resources", "environment", "info", rule, value)
    value["need_ram_gb"] = _gib(need_ram_bytes)
    level: Level = "warn" if probe.total_ram_bytes < need_ram_bytes else "ok"
    message = (
        f"{value['ram_gb']} GiB RAM < {value['need_ram_gb']} GiB needed by genome-wide Promotech"
        if level == "warn"
        else ""
    )
    return Check("environment.resources", "environment", level, rule, value, message)


def environment_checks(
    reqs: Sequence[Requirement], probe: EnvProbe, need_disk_bytes: int, need_ram_bytes: int
) -> list[Check]:
    """``environment.tools`` and ``environment.resources`` for ``reqs`` against ``probe``.

    Free disk never warns (ADR-0005 item 4): ``need_disk_bytes`` is only recorded.
    """
    return [_tools_check(reqs, probe), _resources_check(probe, need_disk_bytes, need_ram_bytes)]
