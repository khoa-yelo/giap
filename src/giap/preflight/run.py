"""``run_preflight``: read every input, evaluate all rules, normalise (GIAP_SPEC.md §5.1, §5.3,
§6.2, §7.0; ADR-0005 item 3).

Readers come from ``giap.io``, rules from ``giap.preflight``; this module wires them together,
probes the environment and writes ``OUT/preflight/`` plus, when there is no error, the normalised
work files under ``OUT/work/preflight/``. Outputs are atomic: a stale ``work/preflight/`` is moved
aside first, files are written into a temporary directory under ``OUT`` and moved into place with
``os.replace``, ``preflight/preflight.json`` last (the completion marker). A work file whose content
is unchanged keeps the earlier file (and its mtime), so Nextflow ``-resume`` caches (#78 review 3a).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import pandas as pd

import giap
from giap.annotate.tools import required_tools
from giap.config import GiapConfig
from giap.exitcodes import ExitCode
from giap.io._open import open_text_write
from giap.io.abundance import AbundanceTable, read_abundance, resolve_scale
from giap.io.fasta import (
    FastaFormatError,
    GenomeRecord,
    read_genome_fasta,
    sanitize_id,
    write_fasta,
)
from giap.io.issues import Issue
from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides
from giap.io.protein_fasta import ProteinRecord, Route, read_protein_fasta, route_by_length
from giap.io.tsv import write_tsv
from giap.model import Replicon
from giap.preflight.checks import Check, JsonValue
from giap.preflight.environment import (
    MANIFEST,
    RESOURCES_RULE,
    TOOLS_RULE,
    EnvProbe,
    ProbeError,
    environment_checks,
    missing_tools,
    probe_resources,
    read_manifest,
    requirements,
    resource_needs,
)
from giap.preflight.evidence import (
    abundance_checks,
    candidate_tier,
    inputs_checks,
    peptide_checks,
    protein_checks,
    tier_check,
)
from giap.preflight.genome import (
    classify_replicons,
    genome_checks,
    resolve_topology,
    summarize_genome,
)
from giap.preflight.report import (
    InputRef,
    PreflightReport,
    render_preflight_md,
    write_preflight_json,
)

__all__ = [
    "NormalisedInputs",
    "PreflightInputs",
    "exit_code",
    "run_preflight",
    "write_outputs",
]

Mode = Literal["preflight", "run"]
ANNOTATION_SUFFIXES = (".json", ".gff3", ".gff", ".gbff", ".gbk")
_REPLICON_TYPE = {"chromosome": "chromosome", "plasmid": "plasmid", "unknown": "contig"}


@dataclass(frozen=True, slots=True)
class PreflightInputs:
    genome: Path
    proteins: Path | None = None
    peptides: Path | None = None
    abundance: Path | None = None
    annotation: Path | None = None


@dataclass(frozen=True, slots=True, eq=False)
class NormalisedInputs:
    """In-memory normalised inputs; ``None`` where the input (or its rows) is absent."""

    genome: tuple[GenomeRecord, ...]
    replicons: tuple[Replicon, ...]
    proteins_long: tuple[ProteinRecord, ...] | None = None
    peptides: pd.DataFrame | None = None
    peptide_runs: pd.DataFrame | None = None
    abundance: pd.DataFrame | None = None
    resources: tuple[int, int] | None = None  # live (free_disk, ram) bytes; logged, not reported


# ---------------------------------------------------------------- environment


def _tools_error(message: str, value: JsonValue) -> Check:
    return Check("environment.tools", "environment", "error", TOOLS_RULE, value, message, 2)


def _environment(
    cfg: GiapConfig,
    probe: EnvProbe | None,
    mode: Mode,
    outdir: Path,
    db_dir: str | None,
    tools_missing: list[str],
) -> tuple[list[Check], tuple[int, int] | None]:
    """Environment checks and, for the real probe, the raw ``(free_disk, ram)`` bytes to log.

    ``probe=None`` selects the real probe: versions from ``<db_dir>/db_manifest.json``
    (missing: info under ``mode="preflight"``, error under ``mode="run"``; unreadable: one
    ``environment.tools`` error), ``tools_missing`` registry tools reported with the missing
    databases, and resources probed live where ``outdir`` is.
    """
    promoters = cfg.annotation.promoters.mode
    reqs = requirements(cfg.modules, promoters)
    needs = resource_needs(promoters)
    if probe is not None:
        return environment_checks(reqs, probe, *needs), None
    try:
        found = read_manifest(db_dir)
    except ProbeError as exc:
        return [_tools_error(f"environment probe failed: {exc}", {"probe_error": str(exc)})], None
    missing = f"{MANIFEST} not found in $GIAP_DB_DIR (written by giap setup)"
    if found is None and mode == "run":
        value: dict[str, JsonValue] = {MANIFEST: "missing"}
        if tools_missing:
            value["missing"] = list[JsonValue](tools_missing)
            missing += "; missing: " + ", ".join(tools_missing)
        return [_tools_error(missing, value)], None
    resources: tuple[int, int] | None = None
    try:
        resources = probe_resources(outdir)
    except OSError as exc:
        res = Check(
            "environment.resources", "environment", "info", RESOURCES_RULE, None,
            f"resources could not be probed: {type(exc).__name__}",
        )  # fmt: skip
    else:
        live = EnvProbe(found={}, free_disk_bytes=resources[0], total_ram_bytes=resources[1])
        res = environment_checks([], live, *needs)[1]
    if found is None:
        tools = Check(
            "environment.tools", "environment", "info", TOOLS_RULE, {MANIFEST: "missing"},
            missing + "; tool and database versions are checked by giap run",
        )  # fmt: skip
    else:
        tools = environment_checks(
            reqs, EnvProbe(found=found, free_disk_bytes=0, total_ram_bytes=0), 0, 0, tools_missing
        )[0]
    return [tools, res], resources


# ---------------------------------------------------------------- inputs


def _input_ref(role: str, path: Path) -> InputRef:
    with path.open("rb") as fh:
        digest = hashlib.file_digest(fh, "sha256").hexdigest()
    return InputRef(role=role, file=path.name, sha256=digest, size=path.stat().st_size)


def _scrub(text: object, path: Path) -> str:
    """``text`` with the input's directory removed (``preflight.json`` holds file names)."""
    return str(text).replace(str(path), path.name)


def _annotation_check(path: Path) -> Check:
    name = path.name.lower().removesuffix(".gz")
    ok = name.endswith(ANNOTATION_SUFFIXES)
    return Check(
        "inputs.annotation_kind",
        "inputs",
        "ok" if ok else "error",
        "bakta.json, GFF3 or GenBank by extension (optionally .gz) → else error; parsed in M2",
        path.name,
        "" if ok else f"unknown annotation kind {path.name!r}",
    )


def _short_protein_rows(
    records: Sequence[ProteinRecord], routes: Sequence[Route], cfg: GiapConfig
) -> pd.DataFrame:
    """Proteins of ``peptide_len`` aa as validated peptide rows naming their protein (§5.3)."""
    short = [r for r, route in zip(records, routes, strict=True) if route == "peptide"]
    raw = pd.DataFrame(
        {"peptide": [r.sequence for r in short], "protein_ids": [r.id for r in short]},
        dtype=pd.StringDtype(),
    )
    return validate_peptides(raw, cfg.input)[0]


def run_preflight(
    inp: PreflightInputs,
    cfg: GiapConfig,
    probe: EnvProbe | None,
    mode: Mode,
    outdir: Path | None = None,
    db_dir: Path | None = None,
    *,
    giap_root: Path | None = None,
    profile: str | None = None,
) -> tuple[PreflightReport, NormalisedInputs]:
    """Evaluate every preflight rule; ``probe=None`` selects the real environment probe, which
    measures free disk where ``outdir`` is (default: the working directory), reads the
    manifest from ``db_dir`` (default: ``$GIAP_DB_DIR``; the environment is never modified)
    and, when a Nextflow ``profile`` is given, checks the registry tools the run needs under
    ``giap_root`` (default: ``$GIAP_ROOT``; unset: every tool missing)."""
    lo, hi = cfg.input.peptide_len
    refs = [_input_ref("genome", inp.genome)]

    genome_error: str | None = None
    try:
        genome = read_genome_fasta(inp.genome)
    except FastaFormatError as exc:
        genome, genome_error = [], _scrub(exc, inp.genome)
    summary = summarize_genome(genome, cfg.preflight.n_run_min_bp)
    checks_genome = genome_checks(summary, cfg, genome_error)
    replicons = classify_replicons(
        summary,
        resolve_topology(summary, cfg.annotation.bakta.complete),
        small_bp=cfg.preflight.small_replicon_bp,
        chromosome_min_bp=cfg.preflight.min_largest_contig_warn,
    )

    proteins: list[ProteinRecord] | None = None
    checks_proteins: list[Check] = []
    if inp.proteins is not None:
        refs.append(_input_ref("proteins", inp.proteins))
        protein_error: str | None = None
        try:
            proteins = read_protein_fasta(inp.proteins)
        except FastaFormatError as exc:
            proteins, protein_error = [], _scrub(exc, inp.proteins)
        checks_proteins = protein_checks(proteins, protein_error, (lo, hi))

    user_peptides: pd.DataFrame | None = None
    checks_peptides: list[Check] = []
    if inp.peptides is not None:
        refs.append(_input_ref("peptides", inp.peptides))
        issues: list[Issue]
        try:
            raw, issues = read_peptide_table(inp.peptides)
        except ValueError as exc:
            raw, issues = (
                pd.DataFrame(),
                [Issue("malformed_table", "error", _scrub(exc, inp.peptides))],
            )
        user_peptides, more = validate_peptides(raw, cfg.input)
        checks_peptides = peptide_checks([*issues, *more], user_peptides, cfg.input)

    routes = [route_by_length(len(r.sequence), (lo, hi)) for r in proteins or []]
    frames = [f for f in (user_peptides,) if f is not None and len(f)]
    if proteins:
        short = _short_protein_rows(proteins, routes, cfg)
        if len(short):
            frames.append(short)
    peptides = runs = None
    if frames:
        combined = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
        peptides, runs = aggregate_peptides(combined)

    has_intensity = bool(user_peptides is not None and user_peptides["intensity"].notna().any())
    abundance: pd.DataFrame | None = None
    abundance_scale: str | None = None
    checks_abundance: list[Check] = []
    if inp.abundance is not None:
        refs.append(_input_ref("abundance", inp.abundance))
        table, ab_issues = read_abundance(inp.abundance)
        ab_issues = [replace(i, message=_scrub(i.message, inp.abundance)) for i in ab_issues]
        resolved = None if table is None else resolve_scale(table.values, cfg.abundance.scale)
        checks_abundance = abundance_checks(
            table,
            ab_issues,
            proteins if inp.proteins is not None and proteins else None,
            resolved or "linear",
            cfg.abundance.scale,
            has_intensity,
            cfg,
        )
        if table is not None:
            abundance, abundance_scale = _abundance_frame(table), resolved

    if inp.annotation is not None:
        refs.append(_input_ref("annotation", inp.annotation))

    env_db = os.environ.get("GIAP_DB_DIR")
    db = str(db_dir.absolute()) if db_dir is not None else env_db
    tools_missing: list[str] = []
    if probe is None and profile is not None:
        # unset $GIAP_ROOT: an empty root, so every tool that needs one is reported missing
        root = giap_root or Path(os.environ.get("GIAP_ROOT") or "/nonexistent/GIAP_ROOT")
        tools_missing = missing_tools(required_tools(cfg), profile, root)
    where = Path.cwd() if outdir is None else outdir
    env, resources = _environment(cfg, probe, mode, where, db, tools_missing)
    tier = candidate_tier(user_peptides, inp.proteins is not None)
    n_usable = None if user_peptides is None else len(user_peptides)
    checks = [
        *inputs_checks(
            inp.peptides is not None, inp.proteins is not None, n_usable_peptides=n_usable
        ),
        *([_annotation_check(inp.annotation)] if inp.annotation is not None else []),
        *checks_genome,
        *checks_proteins,
        *checks_peptides,
        tier_check(tier),
        *checks_abundance,
        *env,
    ]

    routing: dict[str, JsonValue] = {
        "n_proteins_long": routes.count("protein"),
        "n_proteins_peptide_path": routes.count("peptide"),
        "n_proteins_dropped": routes.count("dropped"),
        "n_peptides": 0 if peptides is None else len(peptides),
        "abundance_used": abundance is not None and not has_intensity,
        "abundance_scale": abundance_scale,
    }
    report = PreflightReport(
        giap_version=giap.__version__,
        inputs=refs,
        checks=checks,
        routing=routing,
        evidence_tier=tier,
        replicons=replicons,
        contig_map=[(r.id, r.original_id) for r in genome],
    )
    long = None
    if proteins is not None:
        long = tuple(r for r, route in zip(proteins, routes, strict=True) if route == "protein")
    norm = NormalisedInputs(
        genome=tuple(genome),
        replicons=tuple(replicons),
        proteins_long=long,
        peptides=peptides,
        peptide_runs=runs,
        abundance=abundance,
        resources=resources,
    )
    return report, norm


def _abundance_frame(table: AbundanceTable) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": pd.array(table.ids, dtype=pd.StringDtype()),
            "raw_id": pd.array(table.raw_ids, dtype=pd.StringDtype()),
            "value": pd.array(
                [float("nan") if v is None else v for v in table.values], dtype="float64"
            ),
        }
    )


# ---------------------------------------------------------------- exit code and outputs


def exit_code(report: PreflightReport, strict: bool) -> ExitCode:
    """1 for any input error, else 2 for any environment error, else 3 for a warning with
    ``strict``, else 0 (§5.1; ADR-0005 item 3)."""
    errors = [c.exit_code for c in report.checks if c.level == "error"]
    if ExitCode.INVALID_INPUT in errors:
        return ExitCode.INVALID_INPUT
    if errors:
        return ExitCode.ENVIRONMENT
    if strict and any(c.level == "warn" for c in report.checks):
        return ExitCode.WARNINGS_STRICT
    return ExitCode.SUCCESS


def _safe_protein_ids(records: Sequence[ProteinRecord]) -> list[str]:
    """``[A-Za-z0-9._-]`` IDs for ``records``; a collision gets a ``_2``, ``_3``, … suffix."""
    taken: set[str] = set()
    out: list[str] = []
    for r in records:
        base = sanitize_id(r.id)
        safe, n = base, 1
        while safe in taken:
            n += 1
            safe = f"{base}_{n}"
        taken.add(safe)
        out.append(safe)
    return out


def _log_resources(outdir: Path, resources: tuple[int, int]) -> None:
    """Append the raw probed resources to ``OUT/logs/events.jsonl`` (§6.7)."""
    disk, ram = resources
    event = {
        "ts": datetime.now(UTC).isoformat(timespec="seconds"),
        "step": "preflight",
        "level": "INFO",
        "event": "environment.resources",
        "data": {"free_disk_bytes": disk, "ram_bytes": ram},
    }
    logs = outdir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    with (logs / "events.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(event, sort_keys=True) + "\n")


def _write_work(work: Path, norm: NormalisedInputs) -> None:
    """The normalised work files (§6.2) into ``work``."""
    work.mkdir()
    write_fasta(work / "genome.norm.fa", ((r.id, r.sequence) for r in norm.genome))
    write_tsv(
        work / "contig_map.tsv",
        ["original_id", "contig"],
        [(r.original_id, r.id) for r in norm.genome],
    )
    with open_text_write(work / "replicons.tsv") as fh:
        for x in norm.replicons:
            topo = "c" if x.topology == "circular" else "l"
            fh.write(f"{x.id}\t{x.id}\t{_REPLICON_TYPE[x.kind]}\t{topo}\t-\n")
    if norm.proteins_long is not None:
        safe = _safe_protein_ids(norm.proteins_long)
        write_fasta(
            work / "proteins.long.faa",
            ((s, r.sequence) for s, r in zip(safe, norm.proteins_long, strict=True)),
        )
        write_tsv(
            work / "protein_map.tsv",
            ["protein_id", "original_id", "members"],
            [(s, r.id, ";".join(r.members)) for s, r in zip(safe, norm.proteins_long, strict=True)],
        )
    for name, frame in (
        ("peptides.norm.parquet", norm.peptides),
        ("peptide_runs.parquet", norm.peptide_runs),
        ("abundance.norm.parquet", norm.abundance),
    ):
        if frame is not None:
            frame.to_parquet(work / name, engine="pyarrow", index=False)


def _keep_unchanged(old: Path, new: Path) -> None:
    """Move each file of ``old`` whose bytes equal its counterpart in ``new`` over that
    counterpart, so an unchanged work file keeps its mtime (Nextflow hashes path, size, mtime)."""
    for prev in old.rglob("*"):
        cur = new / prev.relative_to(old)
        same_size = prev.is_file() and cur.is_file() and prev.stat().st_size == cur.stat().st_size
        if same_size and prev.read_bytes() == cur.read_bytes():
            os.replace(prev, cur)


def write_outputs(outdir: Path, report: PreflightReport, norm: NormalisedInputs) -> None:
    """``OUT/preflight/preflight.{json,md}`` always; ``OUT/work/preflight/`` files (§6.2) only
    when the report has no error.

    Atomic: the previous marker is removed and a stale ``work/preflight/`` is moved into a
    temporary directory under ``outdir`` first; everything is written into that directory and
    moved into place with ``os.replace``, ``preflight.json`` last. A crash leaves neither the
    marker nor normalised files behind. Work files with unchanged content are the earlier files
    (mtime kept), so a ``--resume`` run caches.
    """
    rep, work = outdir / "preflight", outdir / "work" / "preflight"
    rep.mkdir(parents=True, exist_ok=True)
    (rep / "preflight.json").unlink(missing_ok=True)
    if norm.resources is not None:
        _log_resources(outdir, norm.resources)
    tmp = Path(tempfile.mkdtemp(prefix=".preflight-", dir=outdir))
    try:
        if work.is_dir():
            os.replace(work, tmp / "old")
        else:
            shutil.rmtree(work, ignore_errors=True)
        if not any(c.level == "error" for c in report.checks):
            _write_work(tmp / "work", norm)
            if (tmp / "old").is_dir():
                _keep_unchanged(tmp / "old", tmp / "work")
        with open_text_write(tmp / "preflight.md") as fh:
            fh.write(render_preflight_md(report))
        write_preflight_json(tmp / "preflight.json", report)
        if (tmp / "work").exists():
            work.parent.mkdir(exist_ok=True)
            os.replace(tmp / "work", work)
        os.replace(tmp / "preflight.md", rep / "preflight.md")
        os.replace(tmp / "preflight.json", rep / "preflight.json")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
