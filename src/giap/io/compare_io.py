"""File I/O of ``giap compare`` (M9-22; ADR-0046): read and check the runs, write ``CMP/tables``.

Genome identity (ADR-0046 item 4, lead decision 1): the ``sha256`` of the ``role = genome`` input
of ``run_manifest.json`` when that file exists, else of ``preflight/preflight.json``.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from giap.io.run_view import NotAGiapRun, RunView
from giap.io.tsv import write_tsv

__all__ = ["CompareInputError", "LoadedRuns", "default_labels", "genome_digest", "load_runs",
           "read_manifests", "write_tables", "write_version_events"]  # fmt: skip


class CompareInputError(ValueError):
    """The runs cannot be compared (exit code 1)."""


@dataclass(frozen=True)
class LoadedRuns:
    labels: list[str]
    views: dict[str, RunView]
    genes: dict[str, pd.DataFrame]


def default_labels(paths: Sequence[Path]) -> list[str]:
    """Directory basenames; repeats suffixed ``_2``, ``_3``, ..."""
    seen: dict[str, int] = {}
    out = []
    for p in paths:
        name = p.resolve().name
        seen[name] = seen.get(name, 0) + 1
        out.append(name if seen[name] == 1 else f"{name}_{seen[name]}")
    return out


def _genome_sha(path: Path) -> str | None:
    doc = json.loads(path.read_text(encoding="utf-8"))
    for item in doc.get("inputs") or []:
        if isinstance(item, dict) and item.get("role") == "genome" and item.get("sha256"):
            return str(item["sha256"])
    return None


def genome_digest(out: Path) -> str | None:
    """The run's genome SHA-256 (manifest over preflight), or ``None``."""
    for rel in ("run_manifest.json", "preflight/preflight.json"):
        if (out / rel).is_file():
            return _genome_sha(out / rel)
    return None


def _labels(paths: Sequence[Path], labels: str | None) -> list[str]:
    if len(paths) < 2:
        raise CompareInputError(f"giap compare needs two or more runs, got {len(paths)}")
    out = default_labels(paths) if labels is None else [s.strip() for s in labels.split(",")]
    if len(out) != len(paths):
        raise CompareInputError(f"--labels has {len(out)} labels for {len(paths)} runs")
    if len(set(out)) != len(out) or not all(out):
        raise CompareInputError(f"--labels must be distinct and non-empty: {labels!r}")
    for lab in out:  # "&" and ";" join labels in the tables; "none" and "other" are rows there
        if "&" in lab or ";" in lab or lab in ("none", "other"):
            raise CompareInputError(f"run label {lab!r} must not hold '&' or ';' or be "
                                    "'none' or 'other' (set --labels)")  # fmt: skip
    return out


def _check_genome(paths: Sequence[Path], labels: Sequence[str]) -> None:
    digests = {lab: genome_digest(p) for lab, p in zip(labels, paths, strict=True)}
    for lab, sha in digests.items():
        if sha is None:
            raise CompareInputError(
                f"run {lab}: cannot verify the genome (no genome SHA-256 in run_manifest.json "
                "or preflight/preflight.json)"
            )
    first = labels[0]
    for lab in labels[1:]:
        if digests[lab] != digests[first]:
            raise CompareInputError(
                f"runs {first} and {lab} are on different genomes: "
                f"{first} {digests[first]}, {lab} {digests[lab]}"
            )


def load_runs(paths: Sequence[Path], labels: str | None) -> LoadedRuns:
    """Labels, views and genes tables of the runs, after every exit-1 check of the spec."""
    labs = _labels(paths, labels)
    _check_genome(paths, labs)
    views, genes = {}, {}
    for lab, path in zip(labs, paths, strict=True):
        try:
            views[lab] = RunView.from_outdir(path)
        except NotAGiapRun as exc:
            raise CompareInputError(str(exc)) from exc
        frame = views[lab].table("genes_tsv")
        if frame is None:
            raise CompareInputError(f"run {lab}: no tables/genes.tsv in {path}")
        genes[lab] = frame
    ref = set(genes[labs[0]]["locus_tag"])
    for lab in labs[1:]:
        other = set(genes[lab]["locus_tag"])
        if other != ref:
            raise CompareInputError(
                f"locus-tag sets differ: {labs[0]} has {len(ref)} genes, {lab} has {len(other)}; "
                f"{len(ref ^ other)} locus tags are in only one of them"
            )
    return LoadedRuns(labs, views, genes)


def write_tables(cmp: Path, tables: Mapping[str, pd.DataFrame]) -> list[Path]:
    """``CMP/tables/<name>.tsv`` for each frame (§9.2 conventions)."""
    done = []
    for name, frame in tables.items():
        path = cmp / "tables" / f"{name}.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(
            path, list(frame.columns), frame.astype(object).itertuples(index=False, name=None)
        )
        done.append(path)
    return done


def read_manifests(paths: Sequence[Path], labels: Sequence[str]) -> dict[str, dict[str, Any]]:
    """Run label → its parsed ``run_manifest.json`` (``{}`` when absent or not a JSON object)."""
    out: dict[str, dict[str, Any]] = {}
    for lab, path in zip(labels, paths, strict=True):
        file = path / "run_manifest.json"
        try:
            doc = json.loads(file.read_text(encoding="utf-8")) if file.is_file() else {}
        except (OSError, ValueError):
            doc = {}
        out[lab] = doc if isinstance(doc, dict) else {}
    return out


def write_version_events(cmp: Path, lines: Sequence[str]) -> Path:
    """``CMP/logs/events.jsonl``: one ``WARNING`` record (step ``compare``) per
    ``compare_versions`` line; the previous file is kept as ``events.previous.jsonl``."""
    from giap.compare.versions import event_name
    from giap.io.events import append_event, event_record
    from giap.io.manifest_sources import rotate_events

    logs = cmp / "logs"
    rotate_events(logs)
    logs.mkdir(parents=True, exist_ok=True)
    (logs / "events.jsonl").touch()
    for line in lines:
        data = {"check": "compare_versions", "message": line}
        append_event(logs, event_record("compare", event_name(line), data, "WARNING"))
    return logs / "events.jsonl"
