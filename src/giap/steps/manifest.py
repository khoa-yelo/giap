"""The ``manifest`` step and the run tail (GIAP_SPEC.md §6.2, §6.7; MVP-03 D-A, D-B; ADR-0040).

``giap run`` gathers the launch-time sources (:func:`launch_sources`), writes the provenance core
before the graph (:func:`write_provenance`) and calls :func:`run_tail` after it, under both
executors and after a failed step. The registered step is the same builder on files only, for
``giap _step manifest``; the local runner never plans it (``TAIL_STEPS``).
"""

from __future__ import annotations

import sys
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import giap
from giap.config import GiapConfig
from giap.data import calibration
from giap.io import manifest_sources as ms
from giap.manifest import ManifestSources, build_manifest, provenance_core
from giap.runner.local import TAIL_STEPS
from giap.runner.nextflow import _step_name, _task_dir
from giap.steps.registry import StepContext, all_steps, step

_PREFLIGHT_ROLES = ("genome", "proteins", "peptides", "abundance", "annotation")
_OPTIONAL = ("proteins", "abundance", "annotation", "reference_peptides", "replicon_depth")
_SELF_HASHED = ("annotation", "reference_peptides", "replicon_depth")


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _source_root() -> Path:
    return Path(giap.__file__).resolve().parents[2]


def launch_sources(
    cfg: GiapConfig, outdir: Path, given: Mapping[str, Path | None], *, executor: str,
    profile: str | None, resume: bool, publication: bool, db_dir: Path | None,
    command_line: Sequence[str] = (),
) -> dict[str, Any]:  # fmt: skip
    """The :class:`ManifestSources` fields known before launch (D-C, D-D)."""
    from giap.steps import load_step_modules

    load_step_modules()  # the Nextflow path imports no step module itself
    preflight = ms.read_json(outdir / "preflight" / "preflight.json") or {}
    hashed = {r["role"] for r in preflight.get("inputs", [])}
    extra = [ms.input_row(r, Path(p).absolute()) for r in _SELF_HASHED
             if (p := given.get(r)) is not None and r not in hashed]  # fmt: skip
    pre = cfg.annotation.precomputed
    for key in sorted(type(pre).model_fields):
        layer = getattr(pre, key)
        if layer is not None and Path(layer).is_file():
            extra.append(ms.input_row(f"precomputed:{key}", Path(layer).absolute()))
    commit, dirty = ms.git_info(_source_root())
    return {
        "db_manifest": ms.read_db_manifest(db_dir),
        "calibration": calibration.load(),
        "cfg": cfg,
        "preflight": preflight,
        "registered": [s for s in all_steps() if s.name not in TAIL_STEPS],
        "executor": executor,
        "profile": profile,
        "resume": resume,
        "outdir": str(outdir),
        "command_line": list(command_line or sys.argv),
        "started_at": _now(),
        "publication": publication,
        "git_commit": commit,
        "git_dirty": dirty,
        "input_paths": {
            r: str(Path(p).absolute()) for r in _PREFLIGHT_ROLES if (p := given.get(r))
        },
        "extra_inputs": extra,
        "inputs_absent": [r for r in _OPTIONAL if given.get(r) is None],
    }


def write_provenance(outdir: Path, base: Mapping[str, Any]) -> None:
    """``OUT/work/provenance.json``: the provenance core, before the graph (D-B step 3)."""
    core = provenance_core(build_manifest(ManifestSources(**base)))
    ms.write_manifest(outdir / "work" / "provenance.json", core)


def _nf_annotation_meta(trace: Sequence[Mapping[str, str]], outdir: Path) -> Any:
    for row in trace:
        if _step_name(row.get("name") or "") == "parse_bakta" and row.get("status") in {
            "COMPLETED", "CACHED"}:  # fmt: skip
            return ms.read_json(Path(_task_dir(row, outdir / "work")) / "annotation_meta.json")
    return None


def _content_files(outdir: Path) -> dict[str, Any]:
    """The published 03c sources (``giap.runner.local.publish`` names; Nextflow publishDir)."""
    return {
        "sensitivity": ms.read_sensitivity(outdir / "tables" / "sensitivity.tsv"),
        "island_diagnostics": ms.read_island_diagnostics(
            outdir / "context" / "island_diagnostics.parquet"),
        "detection_tier": ms.read_detection_tier(outdir / "tables"),
    }  # fmt: skip


def run_tail(
    outdir: Path, base: Mapping[str, Any], exit_code: int, *, failed_step: str | None = None,
    failed_log: str | None = None, local_steps: Sequence[Mapping[str, Any]] = (),
) -> None:  # fmt: skip
    """Write ``OUT/run_manifest.json`` after the graph (D-B step 5); raises on any error."""
    events = ms.read_events(outdir / "logs" / "events.jsonl")
    trace = meta = None
    reports: list[dict[str, Any]] = []
    if base["executor"] == "nextflow":
        trace = ms.read_trace(outdir / "logs" / "trace.tsv")
        meta = _nf_annotation_meta(trace or [], outdir)
        locks = _source_root() / "workflow" / "envs" / "locks"
        reports = ms.tool_reports(trace or [], outdir / "work", locks)
    else:
        meta = ms.read_json(outdir / "work" / "annotation_meta")
    src = ManifestSources(**base, finished_at=_now(), exit_code=exit_code,
                          failed_step=failed_step, failed_log=failed_log, events=events,
                          annotation_meta=meta, trace=trace, local_steps=local_steps,
                          tool_reports=reports, **_content_files(outdir))  # fmt: skip
    ms.write_manifest(outdir / "run_manifest.json", build_manifest(src))


@step(
    "manifest",
    inputs=["preflight"],
    optional_inputs=["provenance", "events", "trace", "annotation_meta", "sensitivity",
                     "island_diagnostics"],
    outputs=["run_manifest"],
)  # fmt: skip
def manifest_step(ctx: StepContext) -> None:
    """``run_manifest.json`` from files only (the registry path; ``giap run`` uses the tail)."""
    ins = {k: Path(v[0] if isinstance(v, list) else v) for k, v in ctx.inputs.items()}
    prov = ms.read_json(ins.get("provenance")) or {}  # git and --publication of the run
    commit, dirty = prov.get("git_commit"), prov.get("git_dirty")
    trace = ms.read_trace(ins["trace"]) if "trace" in ins else None
    src = ManifestSources(
        cfg=ctx.cfg, preflight=ms.read_json(ins["preflight"]) or {},
        registered=[s for s in all_steps() if s.name not in TAIL_STEPS],
        executor="local" if trace is None else "nextflow", command_line=list(sys.argv),
        finished_at=_now(), publication=bool(prov.get("publication")), git_commit=commit,
        git_dirty=dirty, events=ms.read_events(ins["events"]) if "events" in ins else (),
        annotation_meta=ms.read_json(ins.get("annotation_meta")), trace=trace,
        sensitivity=ms.read_sensitivity(ins.get("sensitivity")),
        island_diagnostics=ms.read_island_diagnostics(ins.get("island_diagnostics")),
        calibration=calibration.load(),
    )  # fmt: skip
    ms.write_manifest(ctx.outputs["run_manifest"], build_manifest(src))
