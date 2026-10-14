"""``run_preflight``: read every input, evaluate all rules, normalise (GIAP_SPEC.md §5.1, §5.3,
§6.2, §7.0; ADR-0005 item 3).

Readers come from ``giap.io``, rules from ``giap.preflight``; this module wires them together,
probes the environment and writes ``OUT/preflight/`` plus, when there is no error, the normalised
work files under ``OUT/work/preflight/``.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

import pandas as pd

import giap
from giap.config import GiapConfig
from giap.exitcodes import ExitCode
from giap.io._open import open_text_write
from giap.io.abundance import AbundanceTable, read_abundance, resolve_scale
from giap.io.fasta import FastaFormatError, GenomeRecord, read_genome_fasta, write_fasta
from giap.io.issues import Issue
from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides
from giap.io.protein_fasta import ProteinRecord, read_protein_fasta, route_by_length
from giap.io.tsv import write_tsv
from giap.model import Replicon
from giap.preflight.checks import Check, JsonValue
from giap.preflight.environment import (
    GIB,
    RESOURCES_RULE,
    TOOLS_RULE,
    EnvProbe,
    environment_checks,
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
    "ProbeError",
    "exit_code",
    "probe_environment",
    "run_preflight",
    "write_outputs",
]

Mode = Literal["preflight", "run"]
MANIFEST = "db_manifest.json"
CAPABILITIES = Path(".giap") / "capabilities.json"
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


class ProbeError(Exception):
    """The environment probe could not read its sources."""


# ---------------------------------------------------------------- environment probe


def _read_object(path: Path) -> dict[str, object]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ProbeError(f"{path.name} not found") from exc
    except (OSError, ValueError) as exc:  # UnicodeDecodeError and JSONDecodeError included
        raise ProbeError(f"{path.name} is unreadable: {type(exc).__name__}: {exc}") from exc
    if not isinstance(data, dict):
        raise ProbeError(f"{path.name} is not a JSON object")
    return data


def _number(doc: dict[str, object], *keys: str) -> float | None:
    for key in keys:
        v = doc.get(key)
        if isinstance(v, int | float) and not isinstance(v, bool):
            return float(v)
    return None


def _resources(caps: dict[str, object]) -> tuple[int, int]:
    """``(free_disk_bytes, total_ram_bytes)`` from ``capabilities.json`` (§3.6 probe snapshot)."""
    alloc = caps.get("allocation")
    ram = _number(alloc, "mem_gb", "node_mem_gb") if isinstance(alloc, dict) else None
    disk = _number(caps, "giap_root_free_gb")
    if ram is None or disk is None:
        raise ProbeError("capabilities.json lacks allocation memory or giap_root_free_gb")
    return int(disk * 1e9), int(ram * GIB)


def probe_environment() -> EnvProbe | None:
    """Real probe: versions from ``$GIAP_DB_DIR/db_manifest.json``, resources from
    ``.giap/capabilities.json`` (under ``$GIAP_REPO``, else the working directory).

    Returns ``None`` when the manifest does not exist; raises :class:`ProbeError` when a source
    is unreadable. A manifest entry is ``name: "version"`` or ``name: {"version": ...}``.
    """
    db_dir = os.environ.get("GIAP_DB_DIR")
    if not db_dir or not (Path(db_dir) / MANIFEST).exists():
        return None
    manifest = _read_object(Path(db_dir) / MANIFEST)
    found: dict[str, str | None] = {}
    for name, entry in manifest.items():
        v = entry.get("version") if isinstance(entry, dict) else entry
        found[name] = v if isinstance(v, str) else None
    disk, ram = _resources(_read_object(_capabilities_path()))
    return EnvProbe(found=found, free_disk_bytes=disk, total_ram_bytes=ram)


def _capabilities_path() -> Path:
    return Path(os.environ.get("GIAP_REPO", ".")) / CAPABILITIES


def _tools_error(message: str, value: JsonValue) -> Check:
    return Check("environment.tools", "environment", "error", TOOLS_RULE, value, message, 2)


def _missing_manifest_checks(mode: Mode, needs: tuple[int, int]) -> list[Check]:
    msg = f"{MANIFEST} not found in $GIAP_DB_DIR (written by giap setup)"
    if mode == "run":
        return [_tools_error(msg, {MANIFEST: "missing"})]
    tools = Check(
        "environment.tools",
        "environment",
        "info",
        TOOLS_RULE,
        {MANIFEST: "missing"},
        msg + "; tool and database versions are checked by giap run",
    )
    try:
        disk, ram = _resources(_read_object(_capabilities_path()))
    except ProbeError as exc:
        res = Check("environment.resources", "environment", "info", RESOURCES_RULE, None, str(exc))
        return [tools, res]
    probe = EnvProbe(found={}, free_disk_bytes=disk, total_ram_bytes=ram)
    return [tools, environment_checks([], probe, *needs)[1]]


def _environment(cfg: GiapConfig, probe: EnvProbe | None, mode: Mode) -> list[Check]:
    promoters = cfg.annotation.promoters.mode
    reqs = requirements(cfg.modules, promoters)
    needs = resource_needs(promoters)
    if probe is None:
        try:
            probe = probe_environment()
        except ProbeError as exc:
            return [_tools_error(f"environment probe failed: {exc}", {"probe_error": str(exc)})]
        if probe is None:
            return _missing_manifest_checks(mode, needs)
    return environment_checks(reqs, probe, *needs)


# ---------------------------------------------------------------- inputs


def _input_ref(role: str, path: Path) -> InputRef:
    with path.open("rb") as fh:
        digest = hashlib.file_digest(fh, "sha256").hexdigest()
    return InputRef(role=role, file=path.name, sha256=digest, size=path.stat().st_size)


def _scrub(exc: Exception, path: Path) -> str:
    """The exception text with the input's directory removed (``preflight.json`` holds names)."""
    return str(exc).replace(str(path), path.name)


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


def _short_protein_rows(records: list[ProteinRecord], cfg: GiapConfig) -> pd.DataFrame:
    """Proteins of ``peptide_len`` aa as validated peptide rows naming their protein (§5.3)."""
    lo, hi = cfg.input.peptide_len
    short = [r for r in records if route_by_length(len(r.sequence), (lo, hi)) == "peptide"]
    raw = pd.DataFrame(
        {"peptide": [r.sequence for r in short], "protein_ids": [r.id for r in short]},
        dtype=pd.StringDtype(),
    )
    return validate_peptides(raw, cfg.input)[0]


def run_preflight(
    inp: PreflightInputs, cfg: GiapConfig, probe: EnvProbe | None, mode: Mode
) -> tuple[PreflightReport, NormalisedInputs]:
    """Evaluate every preflight rule; ``probe=None`` selects the real environment probe."""
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

    frames = [f for f in (user_peptides,) if f is not None and len(f)]
    if proteins:
        short = _short_protein_rows(proteins, cfg)
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
        where = str(inp.abundance)
        ab_issues = [
            replace(i, message=i.message.replace(where, inp.abundance.name)) for i in ab_issues
        ]
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
        *_environment(cfg, probe, mode),
    ]

    routes = [route_by_length(len(r.sequence), (lo, hi)) for r in proteins or []]
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
        long = tuple(r for r in proteins if route_by_length(len(r.sequence), (lo, hi)) == "protein")
    norm = NormalisedInputs(
        genome=tuple(genome),
        replicons=tuple(replicons),
        proteins_long=long,
        peptides=peptides,
        peptide_runs=runs,
        abundance=abundance,
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


def write_outputs(outdir: Path, report: PreflightReport, norm: NormalisedInputs) -> None:
    """``OUT/preflight/preflight.{json,md}`` always; ``OUT/work/preflight/`` files (§6.2) only
    when the report has no error."""
    rep = outdir / "preflight"
    rep.mkdir(parents=True, exist_ok=True)
    write_preflight_json(rep / "preflight.json", report)
    with open_text_write(rep / "preflight.md") as fh:
        fh.write(render_preflight_md(report))
    if any(c.level == "error" for c in report.checks):
        return

    work = outdir / "work" / "preflight"
    work.mkdir(parents=True, exist_ok=True)
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
        write_fasta(
            work / "proteins.long.faa",
            ((r.id, r.sequence) for r in norm.proteins_long),
            safe_ids=False,
        )
    for name, frame in (
        ("peptides.norm.parquet", norm.peptides),
        ("peptide_runs.parquet", norm.peptide_runs),
        ("abundance.norm.parquet", norm.abundance),
    ):
        if frame is not None:
            frame.to_parquet(work / name, engine="pyarrow", index=False)
