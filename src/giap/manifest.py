"""``run_manifest.json``: the pure builder (GIAP_SPEC.md §6.7; ADR-0040;
docs/specs/MVP-03-manifest-finalize.md D-C, D-D, D-E).

``build_manifest`` turns already-parsed sources (:class:`ManifestSources`) into the manifest dict;
it reads no file. ``stable_view`` drops the :data:`VOLATILE` pointers so two runs of the same
inputs and seed compare equal across executors and ``-resume``; ``provenance_digest`` is the
SHA-256 of the canonical JSON of the provenance core (lead decision 6).
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import PurePosixPath
from typing import Any

import giap
from giap.config import GiapConfig
from giap.runner.nextflow import _step_name
from giap.steps.registry import StepSpec

#: JSON pointers (``*`` = every child) that differ between equivalent runs (D-C "Volatile").
VOLATILE: tuple[str, ...] = (
    "/command_line", "/executor", "/profile", "/resume", "/outdir", "/started_at",
    "/finished_at", "/failed_log", "/steps", "/inputs/*/path",
    "/config/annotation/precomputed/*",
)  # fmt: skip
_NAME_ONLY = "/config/annotation/precomputed/*"  # reduced to file names, not removed

#: Registered run-scope tokens (§9.4 Run row, ADR-0036 §2, ADR-0037); a trailing ``:`` is a prefix.
RUN_FLAGS: frozenset[str] = frozenset({
    "layer_missing:", "anti_sd_default", "bakta_db_light", "c_assumed",
    "operon_source_distance_model", "tu_recal_skipped", "too_few_ranked_tus",
    "no_eligible_genes", "p_detect_clipped", "null_degenerate", "tu_recal_reduced",
    # beyond §9.4, accepted by Khoa (ADR-0041)
    "input_prefiltered", "cai_reference_missing", "no_detection_limit_model", "uncalibrated:",
})  # fmt: skip
#: Run-scope tokens awaiting a Khoa decision; empty since ADR-0041 (MVP-03 "Needs Khoa" 1).
RUN_FLAGS_PENDING: frozenset[str] = frozenset()

_OPTIONAL_MODULES = ("bgc", "function", "mobile")
_COPIED = (
    "giap_version", "git_commit", "git_dirty", "executor", "profile", "resume", "outdir",
    "started_at", "finished_at", "exit_code", "failed_step", "failed_log", "publication",
)  # fmt: skip
_MODULE_LAYERS = {"genomad": "mobile", "antismash": "bgc"}  # decision 7
_LAYER = "layer_missing:"
_UNCAL = "uncalibrated:"
_DIAG = ("contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags")
_SENS = ("n_primary", "n_variant", "n_shared", "jaccard")
_DB = ("name", "version", "checksum", "retrieved", "size_bytes")
_TOOL = ("name", "version", "container_digest", "conda_lock_hash")
_STATUS = {"COMPLETED": "completed", "CACHED": "cached", "FAILED": "failed", "ABORTED": "failed"}
_DUR = re.compile(r"([0-9.]+)\s*(ms|d|h|m|s)")
_DUR_S = {"ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}
_RSS = re.compile(r"([0-9.]+)\s*([KMGT]?B)")
_RSS_B = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}


@dataclass(frozen=True, kw_only=True)
class ManifestSources:
    """Everything the builder reads, already parsed (tests/helpers/manifest_data.py)."""

    cfg: GiapConfig
    preflight: Mapping[str, Any]
    registered: Sequence[StepSpec]
    executor: str
    profile: str | None = None
    resume: bool = False
    outdir: str = ""
    command_line: Sequence[str] = ()
    started_at: str = ""
    finished_at: str = ""
    exit_code: int = 0
    failed_step: str | None = None
    failed_log: str | None = None
    publication: bool = False
    giap_version: str = giap.__version__
    git_commit: str | None = None
    git_dirty: bool | None = None
    input_paths: Mapping[str, str] = field(default_factory=dict)
    extra_inputs: Sequence[Mapping[str, Any]] = ()
    inputs_absent: Sequence[str] = ()
    events: Sequence[Mapping[str, Any]] = ()
    annotation_meta: Mapping[str, Any] | None = None
    trace: Sequence[Mapping[str, str]] | None = None
    local_steps: Sequence[Mapping[str, Any]] = ()
    tool_reports: Sequence[Mapping[str, Any]] = ()
    db_manifest: Mapping[str, Any] | None = None
    sensitivity: Sequence[Mapping[str, Any]] | None = None  # tables/sensitivity.tsv rows
    island_diagnostics: Sequence[Mapping[str, Any]] | None = None  # C1 order
    calibration: Mapping[str, Any] | None = None  # giap.data.calibration.load()


def canonical_json(obj: Any) -> str:
    """Sorted, indented, strict JSON plus a newline; NaN/±∞ raise ``ValueError``."""
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def provenance_digest(core: Mapping[str, Any]) -> str:
    """SHA-256 hex of ``canonical_json(core)``."""
    return hashlib.sha256(canonical_json(core).encode("utf-8")).hexdigest()


def _names_only(config: Mapping[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(dict(config))
    pre = out.get("annotation", {}).get("precomputed", {})
    for key, value in pre.items():
        if value is not None:
            pre[key] = PurePosixPath(str(value)).name
    return out


def provenance_core(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """The D-D core: everything known before launch, without paths."""
    keys = ("schema_version", "giap_version", "git_commit", "git_dirty", "publication")
    core: dict[str, Any] = {k: manifest[k] for k in keys}
    core["config"] = _names_only(manifest["config"])
    core["inputs"] = [{k: v for k, v in r.items() if k != "path"} for r in manifest["inputs"]]
    core["databases"] = copy.deepcopy(manifest["databases"])
    core["seeds"] = {"master": manifest["seeds"]["master"]}
    return core


def _drop(doc: Any, parts: list[str]) -> None:
    head, rest = parts[0], parts[1:]
    if isinstance(doc, dict):
        keys = list(doc) if head == "*" else [head] if head in doc else []
    elif isinstance(doc, list):
        keys = list(range(len(doc))) if head == "*" else [int(head)] if head.isdigit() else []
    else:
        return
    for k in reversed(keys):
        if rest:
            _drop(doc[k], rest)
        else:
            del doc[k]


def stable_view(manifest: Mapping[str, Any], volatile: Sequence[str] = VOLATILE) -> dict[str, Any]:
    """A copy of ``manifest`` without the ``volatile`` pointers (precomputed paths → names)."""
    out = copy.deepcopy(dict(manifest))
    for ptr in volatile:
        if ptr == _NAME_ONLY and isinstance(out.get("config"), dict):
            out["config"] = _names_only(out["config"])
        else:
            _drop(out, ptr.split("/")[1:])
    return out


def seconds(text: str | None) -> float | None:
    """Nextflow duration (``1m 2.5s``, ``746ms``) → seconds; None when absent."""
    parts = _DUR.findall(text or "")
    return sum(float(n) * _DUR_S[u] for n, u in parts) if parts else None


def rss_bytes(text: str | None) -> int | None:
    """``269.9 MB`` (binary units) → bytes; None when absent."""
    m = _RSS.fullmatch((text or "").strip())
    return round(float(m.group(1)) * _RSS_B[m.group(2)]) if m else None


def _percent(text: str | None) -> float | None:
    m = re.fullmatch(r"([0-9.]+)%?", (text or "").strip())
    return float(m.group(1)) if m else None


def _trace_steps(rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    """Trace rows grouped by step, last attempt of each task (lead decision 8)."""
    last: dict[str, tuple[int, Mapping[str, str]]] = {}
    tries: dict[str, int] = {}
    order: list[str] = []
    for row in rows:
        task = (row.get("name") or "").strip()
        step = _step_name(task)
        if (row.get("status") or "").strip() not in _STATUS or not step:
            continue
        if step not in order:
            order.append(step)
        att = (row.get("attempt") or "").strip()
        n = int(att) if att.isdigit() else 1
        tries[task] = tries.get(task, 0) + 1
        if task not in last or n >= last[task][0]:
            last[task] = (n, row)
    out = []
    for step in order:
        tasks = [t for t in last if _step_name(t) == step]
        rws = [last[t][1] for t in tasks]
        status = {_STATUS[(r.get("status") or "").strip()] for r in rws}
        entry: dict[str, Any] = {
            "name": step,
            "status": "failed" if "failed" in status else "cached" if status == {"cached"}
            else "completed",
            "n_tasks": len(tasks),
            "n_retries": sum(tries[t] - 1 for t in tasks),
        }  # fmt: skip
        wall = [w for r in rws if (w := seconds(r.get("realtime"))) is not None]
        cpu = [w * p / 100 for r in rws
               if (w := seconds(r.get("realtime"))) is not None
               and (p := _percent(r.get("%cpu"))) is not None]  # fmt: skip
        rss = [b for r in rws if (b := rss_bytes(r.get("peak_rss"))) is not None]
        for key, vals in (("wall_time_s", wall), ("cpu_time_s", cpu), ("peak_rss_bytes", rss)):
            if vals:  # a metric "-" or empty is left out, never 0
                entry[key] = max(vals) if key == "peak_rss_bytes" else sum(vals)
        out.append(entry)
    return out


def _local_steps(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    keep = ("name", "status", "wall_time_s", "cpu_time_s")
    return [{**{k: r[k] for k in keep if k in r}, "n_tasks": 1, "n_retries": 0} for r in rows]


def _skip_reason(cfg: GiapConfig, spec: StepSpec) -> str:
    if spec.module == "opt_in":
        return "opt_in_off"
    if spec.module not in cfg.modules:
        return f"module_off:{spec.module}"
    if spec.precomputed and getattr(cfg.annotation.precomputed, spec.precomputed) is not None:
        return f"precomputed:{spec.precomputed}"
    return "not_run"


def _flags(src: ManifestSources) -> list[str]:
    found: set[str] = set()
    for ev in src.events:
        data = ev.get("data")
        if not isinstance(data, Mapping):
            continue
        if isinstance(data.get("flag"), str):
            found.add(data["flag"])
        if isinstance(data.get("flags"), list):
            found |= {f for f in data["flags"] if isinstance(f, str)}
    meta_flags = (src.annotation_meta or {}).get("flags") or []
    found |= {f for f in meta_flags if isinstance(f, str)}
    found |= {_UNCAL + k for k, v in (src.calibration or {}).items()
              if isinstance(v, Mapping) and v.get("calibrated") is False}  # fmt: skip
    off = {_LAYER + layer for layer, mod in _MODULE_LAYERS.items() if mod not in src.cfg.modules}
    return sorted(found - off)


def finite(obj: Any) -> Any:
    """A copy of ``obj`` with every NaN/±∞ float replaced by None (strict JSON; M6-08)."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, Mapping):
        return {k: finite(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [finite(v) for v in obj]
    return obj


def _event_data(src: ManifestSources, name: str) -> list[Mapping[str, Any]]:
    """Unique ``data`` of the ``name`` events, ordered by last occurrence (a resume duplicate
    collapses; the last payload in file order comes last, lead decision 20)."""
    seen: dict[str, Mapping[str, Any]] = {}
    for ev in src.events:
        if ev.get("event") == name and isinstance(ev.get("data"), Mapping):
            key = json.dumps(finite(ev["data"]), sort_keys=True)
            seen.pop(key, None)
            seen[key] = ev["data"]
    return list(seen.values())


def _content(src: ManifestSources) -> dict[str, Any]:
    """The 03c blocks: sensitivity, context, operons, mapping, calibration (D-C rows c)."""
    fams: dict[str, dict[str, Any]] = {}
    for r in src.sensitivity or []:
        fams.setdefault(r["family"], {})[r["variant"]] = finite({k: r.get(k) for k in _SENS})
    oric = sorted((finite(d) for d in _event_data(src, "context_stats.oriter.multiple_oric")),
                  key=lambda d: (str(d.get("contig")), json.dumps(d, sort_keys=True)))  # fmt: skip
    missed = _event_data(src, "operons.missed_operons")
    search = [{"search_db": d.get("search_db"), "n_queries": d.get("n_queries")}
              for d in _event_data(src, "mapping.search_db_resolved")]  # fmt: skip
    search = [d for d in search if any(v is not None for v in d.values())]
    novel = sorted((finite(d) for d in _event_data(src, "novel_orfs.class_fdr")),
                   key=lambda d: (str(d.get("regime")), str(d.get("class"))))  # fmt: skip
    return {
        "sensitivity": {"level": src.cfg.output.sensitivity, "n_rows": len(src.sensitivity or []),
                        "families": fams},
        "context": {"replicons": [finite({k: r.get(k) for k in _DIAG})
                                  for r in src.island_diagnostics or []], "oric": oric},
        "operons": finite({"pi0": missed[-1].get("pi0"),
                           "bin_counts": missed[-1].get("bin_counts")}) if missed else None,
        "mapping": search[-1] if search else None,
        "novel_orfs": {"classes": novel} if novel else None,
        "calibration": {k: {"calibrated": bool(v.get("calibrated")),
                            "provenance": v.get("provenance")}
                        for k, v in (src.calibration or {}).items() if isinstance(v, Mapping)},
    }  # fmt: skip


def _databases(doc: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """``ok`` items of the ``giap setup`` form, or ``{name, version}`` of the provisional form."""
    if doc is None:
        return []
    if isinstance(doc.get("items"), list):
        rows = [{k: i.get(k) for k in _DB} for i in doc["items"] if i.get("status") == "ok"]
    else:
        rows = [{"name": n, "version": (e.get("version") if isinstance(e, Mapping) else e),
                 "checksum": None, "retrieved": None, "size_bytes": None}
                for n, e in doc.items()]  # fmt: skip
    return sorted(rows, key=lambda r: str(r["name"]))


def _tools(reports: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One row per (name, version) with its sorted processes; ``giap`` is ``giap_version``."""
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    for r in reports:
        if r["name"] == "giap":
            continue
        row = rows.setdefault((r["name"], r["version"]),
                              {**{k: r.get(k) for k in _TOOL}, "processes": []})  # fmt: skip
        for k in ("container_digest", "conda_lock_hash"):
            row[k] = row[k] or r.get(k)
        if r["process"] not in row["processes"]:
            row["processes"] = sorted([*row["processes"], r["process"]])
    return [rows[k] for k in sorted(rows)]


def _warnings(src: ManifestSources) -> list[str]:
    out = {f"preflight:{c['rule_id']}" for c in src.preflight.get("checks", [])
           if c.get("level") == "warn"}  # fmt: skip
    out |= {f"{e.get('step')}:{e.get('event')}" for e in src.events if e.get("level") == "WARNING"}
    return sorted(out)


def _inputs(src: ManifestSources) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for r in src.preflight.get("inputs", []):
        path = src.input_paths.get(r["role"], r["file"])
        rows.append({**{k: r[k] for k in ("role", "file", "sha256", "size")}, "path": path})
    seen = {r["role"] for r in rows}
    rows += [dict(r) for r in src.extra_inputs if r["role"] not in seen]
    return rows


def _annotation(meta: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if meta is None:
        return None
    return {"source": meta.get("annotation_source"), "bakta_version": meta.get("bakta_version"),
            "db_version": meta.get("db_version"), "db_type": meta.get("db_type")}  # fmt: skip


def build_manifest(src: ManifestSources) -> dict[str, Any]:
    """The manifest (D-C) of ``src``; pure and deterministic."""
    cfg = src.cfg
    nf = src.executor == "nextflow"
    steps = _trace_steps(src.trace or []) if nf else _local_steps(src.local_steps)
    executed = {s["name"] for s in steps}
    specs = {s.name: s for s in src.registered}
    seeds: dict[str, int] = {"master": cfg.seed}
    seeds |= {n: zlib.crc32(n.encode()) for n in executed if n in specs and specs[n].kind == "py"}
    flags = _flags(src)
    skipped = [{"name": n, "reason": _skip_reason(cfg, specs[n])}
               for n in sorted(specs) if n not in executed]  # fmt: skip
    m: dict[str, Any] = {k: getattr(src, k) for k in _COPIED}
    m |= {
        "schema_version": cfg.schema_version,
        "command_line": list(src.command_line),
        "status": "completed" if src.exit_code == 0 else "failed",
        "config": cfg.model_dump(mode="json"),
        "inputs": _inputs(src),
        "inputs_absent": sorted(src.inputs_absent),
        "tools": _tools(src.tool_reports),
        "databases": _databases(src.db_manifest),
        "annotation": _annotation(src.annotation_meta),
        "detection_tier": src.preflight.get("evidence_tier"),
        "seeds": seeds,
        "steps": steps,
        "skipped": skipped,
        "modules_off": sorted(set(_OPTIONAL_MODULES) - set(cfg.modules)),
        "flags": flags,
        "layers_missing": sorted({f.removeprefix(_LAYER) for f in flags if f.startswith(_LAYER)}),
        "warnings": _warnings(src),
        "volatile": list(VOLATILE),
        **_content(src),
    }
    m["provenance_digest"] = provenance_digest(provenance_core(m))
    return m
