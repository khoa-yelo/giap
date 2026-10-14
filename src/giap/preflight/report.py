"""``preflight.json`` and ``preflight.md`` (GIAP_SPEC.md §7.0 Outputs, §9.5, §6.6).

The JSON report follows ``giap/schemas/preflight.schema.json``; it holds file names only and no
timestamps, so the same inputs give byte-identical files. Floats keep 6 significant digits.
"""

from __future__ import annotations

import json
import math
import numbers
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from giap.exitcodes import ExitCode
from giap.io._open import open_text_write
from giap.model.coords import Replicon
from giap.preflight.checks import LEVEL_ORDER, Check, JsonValue

__all__ = [
    "SCHEMA_VERSION",
    "InputRef",
    "PreflightReport",
    "render_preflight_md",
    "to_json_dict",
    "write_preflight_json",
]

SCHEMA_VERSION = "1"
_MD_ORDER = tuple(reversed(LEVEL_ORDER))  # error, warn, info, ok
_MD_CELL_MAX = 200


@dataclass(frozen=True, slots=True)
class InputRef:
    """One input file: its role, file name (no directory), SHA-256 and size in bytes."""

    role: str
    file: str
    sha256: str
    size: int


@dataclass(frozen=True, slots=True)
class PreflightReport:
    giap_version: str
    inputs: list[InputRef]
    checks: list[Check]
    routing: Mapping[str, JsonValue]
    evidence_tier: str
    replicons: list[Replicon]
    contig_map: list[tuple[str, str]]

    @property
    def exit_code(self) -> int:
        """1 for an input error, else 2 for an environment error, else 0 (``--strict``: M1-12)."""
        errors = [c for c in self.checks if c.level == "error"]
        if any(c.input != "environment" for c in errors):
            return int(ExitCode.INVALID_INPUT)
        if errors:
            return int(ExitCode.ENVIRONMENT)
        return int(ExitCode.SUCCESS)


def _jsonable(v: object) -> JsonValue:
    """``v`` as plain JSON: numpy scalars converted, floats to 6 significant digits, NaN → null.

    Mapping keys become strings; keys that collide once converted raise ``ValueError``.
    """
    if v is None or isinstance(v, str):
        return v
    if isinstance(v, bool | np.bool_):
        return bool(v)
    if isinstance(v, numbers.Integral):
        return int(v)
    if isinstance(v, numbers.Real):
        f = float(v)
        return float(f"{f:.6g}") if math.isfinite(f) else None
    if isinstance(v, Mapping):
        out: dict[str, JsonValue] = {}
        for k, x in v.items():
            key = k if isinstance(k, str) else str(_jsonable(k))
            if key in out:
                raise ValueError(f"keys collide as {key!r} in preflight.json")
            out[key] = _jsonable(x)
        return dict(sorted(out.items()))
    if isinstance(v, list | tuple):
        return [_jsonable(x) for x in v]
    raise TypeError(f"cannot write {type(v).__name__} value {v!r} to preflight.json")


def _counts(r: PreflightReport) -> dict[str, int]:
    levels = [c.level for c in r.checks]
    return {f"n_{lv}": levels.count(lv) for lv in LEVEL_ORDER}


def to_json_dict(r: PreflightReport) -> dict[str, JsonValue]:
    """The ``preflight.json`` document as plain JSON values (numpy scalars converted throughout)."""
    checks: list[JsonValue] = [
        {
            "rule_id": c.rule_id,
            "input": c.input,
            "level": c.level,
            "rule": c.rule,
            "value": c.value,
            "message": c.message,
            "exit_code": c.exit_code,
        }
        for c in r.checks
    ]
    doc = {
        "schema_version": SCHEMA_VERSION,
        "giap_version": r.giap_version,
        "inputs": [
            {"role": i.role, "file": i.file, "sha256": i.sha256, "size": i.size} for i in r.inputs
        ],
        "checks": checks,
        "summary": dict(_counts(r)),
        "evidence_tier": r.evidence_tier,
        "routing": r.routing,
        "replicons": [
            {"id": x.id, "length": x.length, "topology": x.topology, "kind": x.kind}
            for x in r.replicons
        ],
        "contig_map": [{"id": a, "original_id": b} for a, b in r.contig_map],
    }
    out = _jsonable(doc)
    if not isinstance(out, dict):  # pragma: no cover - doc is a dict
        raise TypeError("preflight document is not an object")
    return out


def write_preflight_json(path: Path, r: PreflightReport) -> None:
    """Write ``preflight.json``: sorted keys, indent 2, UTF-8, one trailing newline."""
    text = json.dumps(
        to_json_dict(r), sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False
    )
    with open_text_write(path) as fh:
        fh.write(text + "\n")


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def _cell(v: JsonValue) -> str:
    text = v if isinstance(v, str) else json.dumps(v, sort_keys=True, ensure_ascii=False)
    text = " ".join(text.split())
    if len(text) > _MD_CELL_MAX:
        text = text[: _MD_CELL_MAX - 1] + "…"
    for ch in "\\|`":
        text = text.replace(ch, "\\" + ch)
    return text.replace("&", "&amp;").replace("<", "&lt;")


def render_preflight_md(r: PreflightReport) -> str:
    """A human-readable summary: counts, evidence tier, replicons, then every check by level."""
    n = _counts(r)
    lines = [
        "# GIAP preflight",
        "",
        f"GIAP {r.giap_version}. "
        f"{_plural(n['n_error'], 'error')}, {_plural(n['n_warn'], 'warning')}, "
        f"{n['n_info']} info, {n['n_ok']} ok. Exit code {r.exit_code}.",
        "",
        f"Evidence tier: {r.evidence_tier}.",
        "",
    ]
    if r.inputs:
        lines += ["| Role | File | Size (bytes) | SHA-256 |", "|---|---|---|---|"]
        lines += [
            f"| {_cell(i.role)} | {_cell(i.file)} | {i.size} | {_cell(i.sha256)} |"
            for i in r.inputs
        ]
        lines.append("")
    if r.replicons:
        lines += ["| Replicon | Length (bp) | Topology | Kind |", "|---|---|---|---|"]
        lines += [f"| {_cell(x.id)} | {x.length} | {x.topology} | {x.kind} |" for x in r.replicons]
        lines.append("")
    lines += ["| Level | Rule | Value | Message |", "|---|---|---|---|"]
    for c in sorted(r.checks, key=lambda c: _MD_ORDER.index(c.level)):
        value = _cell(_jsonable(c.value))
        lines.append(f"| {c.level} | {c.rule_id} | {value} | {_cell(c.message)} |")
    return "\n".join(lines) + "\n"
