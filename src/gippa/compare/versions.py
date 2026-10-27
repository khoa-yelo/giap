"""``compare_versions``: the compatibility check of ``gippa compare`` (M11-14h; ADR-0060 §6).

Pure. The parsed ``run_manifest.json`` of each run (``{}`` when a run has none), keyed by run
label in run order, in; the warning lines out, each starting with its event name:

- ``compare.cross_version``: the package ``major.minor`` differs (a patch alone gives nothing);
- ``compare.tool_versions``: one line per tool whose ``version`` or ``container_digest`` differs
  (a local build, ``local_build: true``, compares ``source_sha256``, ``key`` and ``runtime``
  instead of its host-specific digest);
- ``compare.db_versions``: one line per database whose ``version`` differs;
- ``compare.operon_source``: the operon source differs;
- ``compare.output_schema``: the ``output_schema_version`` major differs (``cross_major`` turns
  it into the exit-1 error unless ``--allow-cross-major``).

Every key is read with ``.get``: a value absent from a run (older manifests) is skipped, so a
key needs a value in two or more runs to warn.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

__all__ = ["CROSS_MAJOR", "compare_versions", "cross_major", "event_name"]

#: the exit-1 message of differing ``output_schema_version`` majors (``{majors}`` = ``1.x and 2.x``)
CROSS_MAJOR = "compare: output schema {majors} differ; rerun one with the same major"

_MAJOR_MINOR = re.compile(r"^\s*(\d+)\.(\d+)")
_LOCAL_KEYS = ("source_sha256", "key", "runtime")

Manifests = Mapping[str, Mapping[str, Any]]


def event_name(line: str) -> str:
    """The event of a warning line (its prefix before ``": "``)."""
    return line.split(": ", 1)[0]


def _join(values: Sequence[str]) -> str:
    return values[0] if len(values) == 1 else f"{', '.join(values[:-1])} and {values[-1]}"


def _distinct(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _str(doc: Mapping[str, Any], key: str) -> str | None:
    value = doc.get(key)
    return value if isinstance(value, str) and value else None


def _release(version: str) -> tuple[str, ...]:
    m = _MAJOR_MINOR.match(version)
    return (m.group(1), m.group(2)) if m else (version,)


def _major(version: str) -> str | None:
    m = _MAJOR_MINOR.match(version)
    return m.group(1) if m else None


def _cross_version(manifests: Manifests) -> list[str]:
    versions = [v for doc in manifests.values() if (v := _str(doc, "gippa_version"))]
    if len({_release(v) for v in versions}) < 2:
        return []
    return [f"compare.cross_version: runs from {_join(_distinct(versions))}"]


def _schema_majors(manifests: Manifests) -> list[str]:
    raw = [v for doc in manifests.values() if (v := _str(doc, "output_schema_version"))]
    return _distinct(m for v in raw if (m := _major(v)) is not None)


def cross_major(manifests: Manifests) -> str | None:
    """The exit-1 message when the runs' ``output_schema_version`` majors differ, else ``None``."""
    majors = _schema_majors(manifests)
    return CROSS_MAJOR.format(majors=_join([f"{m}.x" for m in majors])) if len(majors) > 1 else None


def _output_schema(manifests: Manifests) -> list[str]:
    majors = _schema_majors(manifests)
    if len(majors) < 2:
        return []
    return [f"compare.output_schema: output schema {_join([f'{m}.x' for m in majors])} differ"]


def _entries(doc: Mapping[str, Any], key: str) -> dict[str, Mapping[str, Any]]:
    items = doc.get(key)
    if not isinstance(items, list):
        return {}
    return {str(i["name"]): i for i in items if isinstance(i, Mapping) and i.get("name")}


def _short(digest: Any) -> str | None:
    if not isinstance(digest, str) or not digest:
        return None
    algo, sep, hexes = digest.partition(":")
    return f"{algo}:{hexes[:12]}" if sep else digest[:12]


def _tool_text(entry: Mapping[str, Any]) -> tuple[tuple[Any, ...], str]:
    """``(what is compared, what the line shows)`` of one tool entry."""
    version = entry.get("version")
    parts = [str(version) if version is not None else "NA"]
    if entry.get("local_build") is True:
        values = tuple(entry.get(k) for k in _LOCAL_KEYS)
        for k, v in zip(_LOCAL_KEYS, values, strict=True):
            if v is not None:
                parts.append(f"{k}={_short(v) if k == 'source_sha256' else v}")
        return (version, True, *values), " ".join(parts)
    digest = entry.get("container_digest")
    if (short := _short(digest)) is not None:
        parts.append(short)
    return (version, digest), " ".join(parts)


def _tool_versions(manifests: Manifests) -> list[str]:
    per_run = [_entries(doc, "tools") for doc in manifests.values()]
    names = _distinct(n for run in per_run for n in run)
    lines = []
    for name in names:
        seen = [_tool_text(run[name]) for run in per_run if name in run]
        if len(seen) > 1 and len({key for key, _ in seen}) > 1:
            lines.append(f"compare.tool_versions: {name} {' vs '.join(t for _, t in seen)}")
    return lines


def _db_versions(manifests: Manifests) -> list[str]:
    per_run = [_entries(doc, "databases") for doc in manifests.values()]
    names = _distinct(n for run in per_run for n in run)
    lines = []
    for name in names:
        seen = [str(v) for run in per_run if name in run and (v := run[name].get("version"))]
        if len(seen) > 1 and len(set(seen)) > 1:
            lines.append(f"compare.db_versions: {name} {' vs '.join(seen)}")
    return lines


def _operon_source(manifests: Manifests) -> list[str]:
    seen = [(lab, s) for lab, doc in manifests.items() if (s := _str(doc, "operon_source"))]
    if len(seen) < 2 or len({s for _, s in seen}) < 2:
        return []
    return ["compare.operon_source: " + ", ".join(f"run {lab} {s}" for lab, s in seen)]


def compare_versions(manifests: Manifests) -> list[str]:
    """The warning lines of ``manifests`` (run label → parsed manifest, in run order)."""
    return [
        *_cross_version(manifests),
        *_output_schema(manifests),
        *_tool_versions(manifests),
        *_db_versions(manifests),
        *_operon_source(manifests),
    ]
