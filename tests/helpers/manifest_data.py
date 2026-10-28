"""Builders for the MVP-03a manifest tests (docs/specs/MVP-03-manifest-finalize.md, part 03a).

Interface assumed (test-author reading of the spec; the spec names ``build_manifest(sources) ->
dict`` and leaves the sources open):

``gippa.manifest.ManifestSources`` is a keyword-only dataclass holding everything the pure builder
reads, already parsed (no file I/O in ``gippa.manifest``):

====================  =========================================================================
field                 meaning
====================  =========================================================================
``cfg``               the resolved ``GippaConfig``
``preflight``         the parsed ``preflight/preflight.json`` (``inputs``, ``checks``,
                      ``evidence_tier``, ...)
``registered``        every registered ``gippa.steps.registry.StepSpec``
``executor``          ``"local"`` or ``"nextflow"``
``profile``           the ``--profile`` string, or None
``resume``            the ``--resume`` flag
``outdir``            ``OUT`` as a string
``command_line``      ``sys.argv``
``started_at``        ISO-8601 string
``finished_at``       ISO-8601 string
``exit_code``         the run's exit code (0 → ``status completed``, else ``failed``)
``failed_step``       the failed step name, or None
``failed_log``        the failed step's log path, or None
``publication``       the ``--publication`` flag
``gippa_version``      ``gippa.__version__``
``git_commit``        ``str`` or None
``git_dirty``         ``bool`` or None
``input_paths``       preflight role → the path given on the command line (``inputs[].path``)
``extra_inputs``      rows the manifest hashed itself (``{role, file, path, sha256, size}``),
                      e.g. role ``precomputed:<key>``
``inputs_absent``     optional ``gippa run`` inputs that were not given
``events``            the parsed ``OUT/logs/events.jsonl`` lines, in file order
``annotation_meta``   the parsed ``annotation_meta.json``, or None
``trace``             the Nextflow ``trace.tsv`` rows (``csv.DictReader`` dicts, human-readable
                      units), or None (local executor, or no trace)
``local_steps``       local executor: ``{name, status, wall_time_s, cpu_time_s}`` per executed
                      step, in plan order
``tool_reports``      one ``{process, name, version, container_digest, conda_lock_hash}`` per tool
                      line of a task's ``versions.yml``
``db_manifest``       the parsed ``<db_dir>/db_manifest.json``, or None
====================  =========================================================================

Only ``cfg``, ``preflight``, ``registered`` and ``executor`` are required; every other field has
an empty default (``None``, ``()``, ``{}``, ``False``, ``0``, ``""``; ``gippa_version`` defaults to
``gippa.__version__``).
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "src" / "gippa" / "schemas" / "run_manifest.schema.json"

SHA_GENOME = "a" * 64
SHA_PEPTIDES = "b" * 64
SHA_LAYER = "c" * 64


def preflight_doc(
    *,
    inputs: Sequence[Mapping[str, Any]] | None = None,
    checks: Sequence[Mapping[str, Any]] = (),
    tier: str = "peptides",
) -> dict[str, Any]:
    """A ``preflight.json`` document with the fields the manifest reads."""
    rows = inputs if inputs is not None else [
        {"role": "genome", "file": "genome.fa", "sha256": SHA_GENOME, "size": 1200},
        {"role": "peptides", "file": "peptides.tsv", "sha256": SHA_PEPTIDES, "size": 300},
    ]  # fmt: skip
    return {
        "schema_version": "1.0",
        "gippa_version": "0.0.0",
        "inputs": [dict(r) for r in rows],
        "checks": [dict(c) for c in checks],
        "evidence_tier": tier,
        "summary": {},
        "routing": {},
        "replicons": [],
        "contig_map": [],
    }


def check(rule_id: str, level: str) -> dict[str, Any]:
    msg = f"{rule_id} message with /some/path"
    return {"rule_id": rule_id, "input": "genome", "level": level, "rule": rule_id,
            "value": None, "message": msg, "exit_code": 0}  # fmt: skip


def spec(
    name: str,
    *,
    kind: str = "py",
    module: str = "core",
    precomputed: str | None = None,
) -> Any:
    """A registry ``StepSpec`` (no function needed by the builder)."""
    from gippa.steps.registry import StepSpec

    return StepSpec(name=name, kind=kind, inputs=(), outputs=(f"{name}.out",),  # type: ignore[arg-type]
                    optional_inputs=(), module=module, precomputed=precomputed,  # type: ignore[arg-type]
                    func=None)  # fmt: skip


TS = "2026-10-22T00:00:00+00:00"


def event(
    step: str, name: str, data: Mapping[str, Any], level: str = "INFO", ts: str = TS
) -> dict[str, Any]:
    return {"ts": ts, "step": step, "level": level, "event": name, "data": dict(data)}


def layer_missing(step: str, layer: str) -> dict[str, Any]:
    return event(step, f"{step}.layer_missing", {"layer": layer, "flag": f"layer_missing:{layer}"},
                 level="WARNING")  # fmt: skip


def trace_row(name: str, status: str = "COMPLETED", **kw: str) -> dict[str, str]:
    row = {"task_id": "1", "hash": "ab/cdef12", "native_id": "1", "name": name, "tag": "-",
           "status": status, "exit": "0", "attempt": "1", "submit": "-", "duration": "-",
           "realtime": "-", "%cpu": "-", "peak_rss": "-", "peak_vmem": "-", "rchar": "-",
           "wchar": "-", "workdir": "/w/ab/cdef12", "hostname": "-", "container": "-"}  # fmt: skip
    row.update(kw)
    return row


def sources(**over: Any) -> Any:
    """``ManifestSources`` with a baseline (genome + peptides, two core py steps, local)."""
    from gippa.config import GippaConfig
    from gippa.manifest import ManifestSources

    base: dict[str, Any] = {
        "cfg": GippaConfig(seed=7),
        "preflight": preflight_doc(),
        "registered": [spec("detectability"), spec("gene_table")],
        "executor": "local",
        "profile": None,
        "resume": False,
        "outdir": "/runs/OUT",
        "command_line": ["gippa", "run", "--genome", "/data/genome.fa"],
        "started_at": "2026-10-22T10:00:00+00:00",
        "finished_at": "2026-10-22T10:05:00+00:00",
        "exit_code": 0,
        "input_paths": {"genome": "/data/genome.fa", "peptides": "/data/peptides.tsv"},
        "local_steps": [
            {
                "name": "detectability",
                "status": "completed",
                "wall_time_s": 1.5,
                "cpu_time_s": 1.25,
            },
            {"name": "gene_table", "status": "completed", "wall_time_s": 2.0, "cpu_time_s": 1.75},
        ],
    }
    base.update(over)
    return ManifestSources(**base)


def build(**over: Any) -> dict[str, Any]:
    from gippa.manifest import build_manifest

    out = build_manifest(sources(**over))
    assert isinstance(out, dict)
    return out


def schema() -> dict[str, Any]:
    data = json.loads(SCHEMA.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def validate(manifest: Mapping[str, Any]) -> None:
    """Raise ``jsonschema.ValidationError`` unless ``manifest`` validates (draft 2020-12)."""
    import jsonschema

    jsonschema.Draft202012Validator(schema()).validate(manifest)


def schema_errors(manifest: Mapping[str, Any]) -> list[str]:
    import jsonschema

    v = jsonschema.Draft202012Validator(schema())
    return [e.message for e in v.iter_errors(manifest)]
