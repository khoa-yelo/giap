"""``giap`` command line (GIAP_SPEC.md §5.1).

Every subcommand resolves the configuration first (an unknown key is exit code 1). ``preflight``
(M1-12), ``_step``, ``run --executor local`` (M2-02) and ``run`` through Nextflow (M2-07) are
implemented; the other commands report that they are not implemented yet (exit code 5).
"""

from __future__ import annotations

import json
import logging
import os
import platform
import shutil
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, TypeVar, get_args

import click

import giap
from giap.config import GiapConfig, resolve_config
from giap.config import model as m
from giap.exitcodes import ExitCode

F = TypeVar("F", bound=Callable[..., Any])

LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")
REPORT_MESSAGE = (
    "giap: internal error (a bug). Please report it at https://github.com/khoa-yelo/giap/issues "
    "with the command line and the traceback above."
)

_existing_file = click.Path(exists=True, dir_okay=False, path_type=Path)
_existing_dir = click.Path(exists=True, file_okay=False, path_type=Path)
_out_path = click.Path(path_type=Path)


def _choice(literal: Any) -> click.Choice[str]:
    """Click choice built from a config ``Literal`` so flags and model never drift apart."""
    return click.Choice(list(get_args(literal)))


def shared_options(func: F) -> F:
    """``--config``, ``--set``, ``--log-level`` and ``--version``, shared by every subcommand."""
    func = click.option(
        "--log-level",
        type=click.Choice(LOG_LEVELS, case_sensitive=False),
        default="INFO",
        show_default=True,
        help="Logging verbosity.",
    )(func)
    func = click.option(
        "--set",
        "sets",
        multiple=True,
        metavar="KEY=VALUE",
        help="Override any config key (dotted, repeatable; VALUE is YAML).",
    )(func)
    func = click.option("--config", "config_file", type=_existing_file, help="YAML config file.")(
        func
    )
    func = click.version_option(giap.__version__, prog_name="giap")(func)
    return func


def _resolve(
    config_file: Path | None,
    sets: Sequence[str],
    log_level: str,
    flags: dict[str, Any] | None = None,
    publication: bool = False,
) -> GiapConfig:
    """Resolve the config: defaults ← file ← CLI flags ← ``--set``. Exit 1 on an invalid key."""
    logging.basicConfig(level=log_level.upper(), format="%(levelname)s %(name)s: %(message)s")
    flag_sets = [
        f"{key}={json.dumps(value)}" for key, value in (flags or {}).items() if value is not None
    ]
    try:
        return resolve_config(config_file, [*flag_sets, *sets], publication=publication)
    except ValueError as exc:
        click.echo(f"giap: invalid configuration: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)


def _not_implemented(command: str, milestone: str) -> None:
    click.echo(f"giap {command}: not implemented (milestone {milestone})", err=True)
    sys.exit(ExitCode.INTERNAL)


class _Cli(click.Group):
    """Keeps the raw argument list (the manifest's ``command_line``; ``sys.argv`` is not it when
    the CLI is invoked in-process)."""

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        ctx.meta["giap.argv"] = ["giap", *args]
        return super().parse_args(ctx, args)


@click.group(cls=_Cli, context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(giap.__version__, prog_name="giap")
def cli() -> None:
    """GIAP: genome-integrated analysis of proteomics data for bacterial genomes."""


# ---------------------------------------------------------------- setup


@cli.command()
@click.option("--db-dir", type=_out_path, required=True, help="Database directory.")
@click.option("--bakta", type=click.Choice(["light", "full", "skip"]), default="light")
@click.option("--genomad/--no-genomad", default=True)
@click.option("--antismash/--no-antismash", default=False)
@click.option("--promotech-models/--no-promotech-models", default=True)
@click.option("--go", is_flag=True, help="Gene Ontology and GO slims.")
@click.option("--kegg", is_flag=True, help="KEGG resources.")
@click.option("--cog", is_flag=True, help="COG categories.")
@click.option("--contaminants", is_flag=True, help="Contaminant protein sequences.")
@click.option("--check", is_flag=True, help="Only verify what is installed.")
@shared_options
def setup(config_file: Path | None, sets: tuple[str, ...], log_level: str, **opts: Any) -> None:
    """Download databases and models (writes DIR/db_manifest.json)."""
    from giap.setupdb import HttpsFetcher, SetupOptions, ToolRunner, run_setup

    _resolve(config_file, sets, log_level)
    code = run_setup(SetupOptions(**opts), ToolRunner(), HttpsFetcher())
    sys.exit(int(code))


# ---------------------------------------------------------------- preflight

# Nextflow profiles `giap run --profile` accepts, alone or comma-joined (#78 review 3c).
_PROFILES = ("farmshare", "farmshare_local", "apptainer", "conda", "docker", "standard")


def _check_profile(_ctx: click.Context, _param: click.Parameter, value: str) -> str:
    """``value`` if every comma-separated name is a known profile, else a usage error."""
    names = value.split(",")
    bad = [n for n in names if n not in _PROFILES]
    if bad:
        raise click.BadParameter(
            f"unknown profile(s) {', '.join(repr(n) for n in bad)} in {value!r}; "
            f"choose comma-joined from: {', '.join(_PROFILES)}"
        )
    return value


@cli.command()
@click.option("--genome", type=_existing_file, required=True, help="Genome FASTA.")
@click.option("--proteins", type=_existing_file, help="Protein FASTA.")
@click.option("--peptides", type=_existing_file, help="Peptide table (TSV).")
@click.option("--abundance", type=_existing_file, help="Protein abundance table (TSV).")
@click.option("--annotation", type=_existing_file, help="bakta.json, GFF3 or GenBank.")
@click.option("--outdir", type=_out_path, help="Output directory.")
@click.option("--strict", is_flag=True, help="Treat warnings as errors (exit code 3).")
@shared_options
def preflight(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    genome: Path,
    proteins: Path | None,
    peptides: Path | None,
    abundance: Path | None,
    annotation: Path | None,
    outdir: Path | None,
    strict: bool,
    **_: Any,
) -> None:
    """Validate inputs only (seconds)."""
    from giap.preflight.report import render_preflight_md
    from giap.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    cfg = _resolve(config_file, sets, log_level)
    if outdir is not None and outdir.exists() and not outdir.is_dir():
        click.echo(f"giap preflight: --outdir {outdir} exists and is not a directory", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    inp = PreflightInputs(genome, proteins, peptides, abundance, annotation)
    report, norm = run_preflight(inp, cfg, None, "preflight", outdir)
    if outdir is not None:
        write_outputs(outdir, report, norm)
    click.echo(render_preflight_md(report), nl=False)
    sys.exit(exit_code(report, strict))


# ---------------------------------------------------------------- run

_RUN_INPUTS = [
    "genome",
    "peptides",
    "proteins",
    "abundance",
    "annotation",
    "reference_peptides",
    "replicon_depth",
]
# params.inputs key -> (raw input, preflight work file): MAPPING and GENETABLE read the
# normalised files (the gene table's C10 abundance is the normalised Parquet, MVP-01)
_NORM_INPUTS = {
    "peptides": ("peptides", "peptides.norm.parquet"),
    "proteins_long": ("proteins", "proteins.long.faa"),
    "abundance": ("abundance", "abundance.norm.parquet"),
}


@cli.command()
@click.option("--genome", type=_existing_file, required=True, help="Genome FASTA.")
@click.option("--peptides", type=_existing_file, help="Peptide table (TSV).")
@click.option("--proteins", type=_existing_file, help="Protein FASTA.")
@click.option("--abundance", type=_existing_file, help="Protein abundance table (TSV).")
@click.option("--annotation", type=_existing_file, help="bakta.json, GFF3 or GenBank.")
@click.option("--reference-peptides", type=_existing_file, help="Reference peptide table.")
@click.option("--replicon-depth", type=_existing_file, help="Replicon depth table (TSV).")
@click.option("--outdir", type=_out_path, required=True, help="Output directory.")
@click.option("--db-dir", type=_out_path, help="Database directory.")
@click.option("--threads", type=click.IntRange(min=1), help="Threads (config: threads).")
@click.option("--memory", default="64.GB", show_default=True, help="Memory for the workflow.")
@click.option("--modules", help="Comma-separated modules: core,mobile,bgc,function.")
@click.option("--genome-relation", type=_choice(m.GenomeRelation))
@click.option("--abundance-missing", type=_choice(m.AbundanceMissing))
@click.option("--abundance-type", type=_choice(m.AbundanceType))
@click.option("--abundance-scale", type=_choice(m.AbundanceScale))
@click.option("--promoters", type=_choice(m.PromoterMode))
@click.option("--terminators", type=_choice(m.TerminatorTool))
@click.option("--add-bactermfinder", is_flag=True, default=None)
@click.option("--publication", is_flag=True, help="10,000 permutations and bootstrap draws.")
@click.option("--seed", type=int, help="Master seed (config: seed).")
@click.option("--parquet", is_flag=True, default=None, help="Also write Parquet tables.")
@click.option("--sensitivity", type=_choice(m.SensitivityLevel))
@click.option("--strict", is_flag=True, help="Treat warnings as errors (exit code 3).")
@click.option(
    "--profile",
    callback=_check_profile,
    default="farmshare,apptainer",
    show_default=True,
    help="Comma-joined Nextflow profiles: " + ", ".join(_PROFILES) + ".",
)
@click.option("--resume", is_flag=True, help="Resume a previous Nextflow run.")
@click.option(
    "--executor", type=click.Choice(["nextflow", "local"]), default="nextflow", show_default=True
)
@shared_options
def run(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    threads: int | None,
    modules: str | None,
    genome_relation: str | None,
    abundance_missing: str | None,
    abundance_type: str | None,
    abundance_scale: str | None,
    promoters: str | None,
    terminators: str | None,
    add_bactermfinder: bool | None,
    publication: bool,
    seed: int | None,
    parquet: bool | None,
    sensitivity: str | None,
    outdir: Path,
    executor: str,
    **_: Any,
) -> None:
    """Run the analysis."""
    flags: dict[str, Any] = {
        "threads": threads,
        "modules": [m.strip() for m in modules.split(",") if m.strip()] if modules else None,
        "genome_relation": genome_relation,
        "abundance.missing": abundance_missing,
        "abundance.type": abundance_type,
        "abundance.scale": abundance_scale,
        "annotation.promoters.mode": promoters,
        "annotation.terminators.tool": terminators,
        "annotation.terminators.bactermfinder": add_bactermfinder,
        "seed": seed,
        "output.parquet": parquet,
        "output.sensitivity": sensitivity,
        # §5.3: a Bakta JSON given as --annotation is used as is and Bakta is not run, so it is
        # the precomputed Bakta layer (it beats the config file; a --set still wins)
        "annotation.precomputed.bakta_json": _bakta_json(_.get("annotation")),
    }
    cfg = _resolve(config_file, sets, log_level, flags, publication)
    from giap.config import resolve_sample_name

    # §5.2 sample_name: auto is resolved once, before any step (M7-04 lead decision 1)
    sample = resolve_sample_name(cfg, _.get("peptides"), _.get("proteins"))
    cfg = cfg.model_copy(update={"sample_name": sample})
    from giap.annotate.promoters import check_promoter_mode

    try:  # before launching anything (M2-13 Q2, revised lead decision 4)
        check_promoter_mode(cfg.annotation.promoters.mode)
    except ValueError as exc:
        click.echo(f"giap run: error: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    from giap.io.manifest_sources import rotate_events

    rotate_events(outdir.absolute() / "logs")  # one events.jsonl per run (MVP-03 D-B step 1)
    env_db = os.environ.get("GIAP_DB_DIR")
    db_dir: Path | None = _.get("db_dir") or (Path(env_db) if env_db else None)
    given = {k: _.get(k) for k in _RUN_INPUTS}
    run_info = {"publication": publication, "db_dir": db_dir, "given": given}
    if executor == "local":
        keys = ("peptides", "proteins", "abundance", "annotation", "reference_peptides", "db_dir")
        _run_local(cfg, outdir, _["genome"], {k: _.get(k) for k in keys}, run_info)
    inputs = {k: _[k] for k in _RUN_INPUTS if _.get(k) is not None}
    outdir = outdir.absolute()
    inputs.update(_run_preflight(cfg, outdir, inputs, db_dir, bool(_["strict"]), str(_["profile"])))
    from giap.runner.local import contaminant_library

    lib = contaminant_library(db_dir)  # M3-15 library for MAPPING's loci
    if lib is not None:
        inputs["contaminants"] = lib
    if shutil.which("nextflow") is None:
        click.echo("giap run: nextflow not found on PATH (see AGENTS.md setup)", err=True)
        sys.exit(ExitCode.ENVIRONMENT)
    base = _launch_sources(cfg, outdir, run_info, "nextflow", str(_["profile"]), bool(_["resume"]))
    try:
        code = _run_nextflow(
            cfg, outdir, inputs, str(_["memory"]), str(_["profile"]), bool(_["resume"]), db_dir
        )
    except Exception:  # past the provenance core: exit 5, the tail still runs (review B2)
        code = _internal_error(outdir)
    from giap.runner.nextflow import failed_step

    trace = outdir / "logs" / "trace.tsv"
    hit = failed_step(trace) if code != ExitCode.SUCCESS and trace.is_file() else None
    failed = {"failed_step": hit[0], "failed_log": f"{hit[1]}/.command.log"} if hit else {}
    sys.exit(_run_tail(outdir, base, code, **failed))


def _launch_sources(
    cfg: GiapConfig, outdir: Path, info: dict[str, Any], executor: str, profile: str | None,
    resume: bool,
) -> dict[str, Any]:  # fmt: skip
    """Launch-time manifest sources; writes ``OUT/work/provenance.json`` (MVP-03 D-B step 3)."""
    from giap.steps.manifest import launch_sources, write_provenance

    argv = click.get_current_context().meta.get("giap.argv", sys.argv)
    base = launch_sources(cfg, outdir.absolute(), info["given"], executor=executor,
                          profile=profile, resume=resume, publication=info["publication"],
                          db_dir=info["db_dir"], command_line=argv)  # fmt: skip
    write_provenance(outdir.absolute(), base)
    return base


def _internal_error(outdir: Path) -> ExitCode:
    """Log the current traceback to ``OUT/logs/internal_error.log`` and stderr; exit 5."""
    import traceback

    log = outdir.absolute() / "logs" / "internal_error.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(tb := traceback.format_exc(), encoding="utf-8")
    click.echo(f"{tb}{REPORT_MESSAGE} (log: {log})", err=True)
    return ExitCode.INTERNAL


def _run_tail(outdir: Path, base: dict[str, Any], code: int, **kw: Any) -> int:
    """The run tail (MVP-03 D-B steps 5-6): the manifest, even after a failed step. A tail error
    after a successful graph is exit 4; after a failed graph it is only logged."""
    import traceback

    from giap.steps.manifest import run_tail

    outdir = outdir.absolute()
    try:
        run_tail(outdir, base, int(code), **kw)
    except Exception:
        log = outdir / "logs" / "manifest.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(traceback.format_exc(), encoding="utf-8")
        if code == ExitCode.SUCCESS:
            click.echo(f"giap run: step manifest failed; log: {log}", err=True)
            return ExitCode.STEP_FAILED
        click.echo(f"giap run: the run manifest was not written; log: {log}", err=True)
    return int(code)


def _bakta_json(annotation: Path | None) -> str | None:
    """The absolute path of ``--annotation`` when it is a Bakta JSON, else None (GFF3/GBFF go
    through Bakta ``--regions``, §5.3)."""
    if annotation is None or Path(annotation).suffix.lower() != ".json":
        return None
    return str(Path(annotation).absolute())


def _run_preflight(
    cfg: GiapConfig,
    outdir: Path,
    inputs: dict[str, Path],
    db_dir: Path | None,
    strict: bool,
    profile: str,
) -> dict[str, Path]:
    """Preflight in run mode (§5.1): exits on an error (or a warning with ``strict``) before
    Nextflow starts, else returns the normalised ``genome`` and ``replicons`` for params.json
    (plus ``peptides`` and ``proteins_long`` when given, for MAPPING).
    Tools are checked for the Nextflow ``profile`` under ``$GIAP_ROOT`` (run-wiring 4a)."""
    from giap.preflight.report import render_preflight_md
    from giap.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    if outdir.exists() and not outdir.is_dir():
        click.echo(f"giap run: --outdir {outdir} exists and is not a directory", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    inp = PreflightInputs(
        genome=inputs["genome"],
        proteins=inputs.get("proteins"),
        peptides=inputs.get("peptides"),
        abundance=inputs.get("abundance"),
        annotation=inputs.get("annotation"),
    )
    root = os.environ.get("GIAP_ROOT")
    report, norm = run_preflight(
        inp, cfg, None, "run", outdir, db_dir,
        giap_root=Path(root) if root else None, profile=profile,
    )  # fmt: skip
    outdir.mkdir(parents=True, exist_ok=True)
    write_outputs(outdir, report, norm)
    code = exit_code(report, strict)
    if code != ExitCode.SUCCESS:
        click.echo(render_preflight_md(report), err=True, nl=False)
        sys.exit(code)
    work = outdir / "work" / "preflight"
    roots = {"genome": work / "genome.norm.fa", "replicons": work / "replicons.tsv"}
    for key, (raw, name) in _NORM_INPUTS.items():
        if raw in inputs and (work / name).is_file():
            roots[key] = work / name
    report_json = outdir / "preflight" / "preflight.json"  # gene_table's abundance_scale (M1-12)
    if report_json.is_file():
        roots["preflight"] = report_json
    return roots


def _run_nextflow(
    cfg: GiapConfig,
    outdir: Path,
    inputs: dict[str, Path],
    memory: str,
    profile: str,
    resume: bool,
    db_dir: Path | None = None,
) -> ExitCode:
    """§6.3 "Launch": params.json (+ ``db_dir``) + resolved_config.yaml, then ``nextflow run``."""
    from giap.runner import nextflow as nf

    outdir = outdir.absolute()
    params, _resolved = nf.write_params(outdir, cfg, inputs, nf.nextflow_memory(memory), db_dir)
    cmd = nf.launch_command(nf.workflow_main(), params, profile, outdir, resume)
    code = nf.run_nextflow(cmd, outdir)
    if code == ExitCode.SUCCESS:
        from giap.runner.local import write_datapackage

        write_datapackage(outdir)
    return code


def _local_roots(
    cfg: GiapConfig, outdir: Path, genome: Path, inputs: dict[str, Path | None]
) -> dict[str, Path]:
    """Preflight-normalised root work files for the local executor (§6.2, §6.4).

    A preflight error (e.g. no evidence) exits with the preflight exit code; a GFF3/GenBank
    ``--annotation`` without ``annotation.precomputed.bakta_json`` exits 1 naming the key. The
    contaminant library (M3-15), the COG 2020 table (M4-15) and, with the ``function`` module,
    its installed resources (M8-06) come from ``--db-dir`` / ``$GIAP_DB_DIR``; preflight's
    ``abundance`` and ``preflight.json`` are staged when written.
    """
    import os

    from giap.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs
    from giap.runner.local import contaminant_library
    from giap.setupdb.installed import FUNCTION_RESOURCES, cog_table, resource

    inp = PreflightInputs(
        genome,
        inputs["proteins"],
        inputs["peptides"],
        inputs.get("abundance"),
        inputs["annotation"],
    )
    report, norm = run_preflight(inp, cfg, None, "preflight", outdir)
    write_outputs(outdir, report, norm)
    code = exit_code(report, False)
    if code != ExitCode.SUCCESS:
        for c in report.checks:
            if c.level == "error":
                click.echo(f"giap run: preflight error {c.rule_id}: {c.message}", err=True)
        sys.exit(code)
    work = outdir / "work" / "preflight"
    files = {"genome": "genome.norm.fa", "peptides": "peptides.norm.parquet"}
    roots = {k: work / f for k, f in files.items() if (work / f).is_file()}
    if (work / "proteins.long.faa").is_file():
        roots["proteins_long"] = work / "proteins.long.faa"
    ann, ref = inputs["annotation"], inputs["reference_peptides"]
    is_json = ann is not None and ann.name.lower().endswith(".json")
    bakta = ann if is_json else cfg.annotation.precomputed.bakta_json
    if ann is not None and bakta is None:  # a GFF3/GenBank annotation cannot replace Bakta JSON
        click.echo("giap run: error: annotation.precomputed.bakta_json is required with "
                   "--executor local (--annotation is not a Bakta JSON)", err=True)  # fmt: skip
        sys.exit(ExitCode.INVALID_INPUT)
    if bakta is not None:
        roots["bakta_json"] = Path(bakta)
    if ref is not None:
        roots["reference_peptides"] = ref
    db = inputs["db_dir"] or os.environ.get("GIAP_DB_DIR")
    db_path = Path(db) if db else None
    lib = contaminant_library(db_path)
    if lib is not None:
        roots["contaminants"] = lib
    cog = cog_table(db_path)  # M4-15 decision 15: cog_category under the local executor too
    if cog is not None:
        roots["cog2020_def"] = cog
    if "function" in cfg.modules:  # M8-06 item 7: only an installed resource, only when switched on
        for name in FUNCTION_RESOURCES:
            path = resource(db_path, name)
            if path is not None:
                roots[name] = path
    produced = {"abundance": work / "abundance.norm.parquet",
                "preflight": outdir / "preflight" / "preflight.json"}  # fmt: skip
    roots |= {k: p for k, p in produced.items() if p.is_file()}
    return roots


def _run_local(
    cfg: GiapConfig,
    outdir: Path,
    genome: Path,
    inputs: dict[str, Path | None],
    info: dict[str, Any],
) -> None:
    """``--executor local`` (§6.4): every registered step as target, then the run tail; exits
    when done. A layer or tool problem found while planning exits before the provenance core."""
    from giap.runner import local
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    workdir = outdir / "work"
    load_step_modules()
    targets = [s.name for s in reg.all_steps()]
    if not targets:
        click.echo("giap run --executor local: not implemented (no steps registered)", err=True)
        sys.exit(ExitCode.INTERNAL)
    roots = _local_roots(cfg, outdir, genome, inputs)
    try:
        local.plan_local(cfg, targets, set(roots))
    except local.LayerMissing as exc:
        click.echo(f"giap run: error: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    except local.ToolUnavailable as exc:
        click.echo(f"giap run: error: {exc}", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    base = _launch_sources(cfg, outdir, info, "local", None, False)
    timings: list[dict[str, Any]] = []
    code, failed = ExitCode.SUCCESS, {}
    # step events (e.g. the ``search_db: auto`` resolution) go to OUT/logs/events.jsonl (§6.7)
    restore = _record_step_events(outdir / "logs")
    try:
        local.run_local(cfg, workdir, targets, roots, timings)
        local.publish(workdir, outdir)
        local.write_datapackage(outdir)
    except reg.StepError as exc:
        log = exc.write_log(workdir.parent / "logs")
        click.echo(f"giap run: {exc} (see {log})", err=True)
        code, failed = ExitCode.STEP_FAILED, {"failed_step": exc.name, "failed_log": str(log)}
    except reg.WiringError as exc:
        click.echo(f"giap run: internal error: {exc}", err=True)
        code = ExitCode.INTERNAL
    except Exception:  # past the provenance core: exit 5, the tail still runs (review B2)
        code = _internal_error(outdir)
    finally:
        restore()
    sys.exit(_run_tail(outdir, base, code, local_steps=timings, **failed))


def _record_step_events(logs_dir: Path) -> Callable[[], None]:
    """Send ``giap.steps`` events to ``logs_dir/events.jsonl`` (§6.7); returns the undo."""
    from giap.io.events import EventsHandler

    steps_log = logging.getLogger("giap.steps")
    events = EventsHandler(logs_dir)
    steps_log.addHandler(events)
    prev_level = steps_log.level
    if steps_log.getEffectiveLevel() > logging.INFO:
        steps_log.setLevel(logging.INFO)

    def restore() -> None:
        steps_log.removeHandler(events)
        steps_log.setLevel(prev_level)

    return restore


# ---------------------------------------------------------------- makedb


@cli.command()
@click.option("--annotation", type=_existing_file, help="bakta.json.")
@click.option("--genome", type=_existing_file, help="Genome FASTA (default: the JSON sequences).")
@click.option("--outdir", type=_out_path, required=True, help="Output database directory.")
@click.option("--db-dir", type=_out_path, help="Database directory (default: $GIAP_DB_DIR).")
@click.option("--six-frame-orfs", is_flag=True, help="Add six-frame ORFs.")
@click.option("--min-orf-aa", type=int, default=20, show_default=True)
@click.option("--patch", type=_existing_file, help="annotation_patch.gff3 from a previous run.")
@click.option("--no-contaminants", is_flag=True, help="Do not append contaminant proteins.")
@shared_options
def makedb(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    annotation: Path | None,
    genome: Path | None,
    outdir: Path,
    db_dir: Path | None,
    six_frame_orfs: bool,
    min_orf_aa: int,
    patch: Path | None,
    no_contaminants: bool,
) -> None:
    """Build an optional strain search database (no FDR is built on it)."""
    import os

    from giap.makedb.build import MakedbOptions, build_strain_db

    _resolve(config_file, sets, log_level)
    if annotation is None:
        click.echo("giap makedb: no --annotation: Bakta wrapper required (M2)", err=True)
        sys.exit(ExitCode.ENVIRONMENT)
    if min_orf_aa < 1:
        click.echo(f"giap makedb: --min-orf-aa must be >= 1, got {min_orf_aa}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    contaminants: Path | None = None
    if not no_contaminants:
        env = os.environ.get("GIAP_DB_DIR")
        root = db_dir or (Path(env) if env else None)
        contaminants = root / "contaminants" / "contaminants.fasta" if root else None
        if contaminants is None or not contaminants.is_file():
            click.echo(
                f"giap makedb: contaminant library not found ({contaminants or 'no --db-dir'}); "
                "run giap setup --contaminants or pass --no-contaminants",
                err=True,
            )
            sys.exit(ExitCode.ENVIRONMENT)
    opts = MakedbOptions(
        annotation, genome, outdir, patch, six_frame_orfs, min_orf_aa, contaminants
    )
    try:
        result = build_strain_db(opts)
    except (ValueError, OSError) as exc:
        click.echo(f"giap makedb: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    click.echo(f"giap makedb: wrote {result.fasta} {dict(result.counts)}", err=True)
    sys.exit(ExitCode.SUCCESS)


# ---------------------------------------------------------------- convert


@cli.command()
@click.argument("inputs", nargs=-1, required=True, type=_existing_file)
@click.option(
    "--from",
    "engine",
    type=click.Choice(["diann", "casanovo", "instanovo", "cascadia", "generic"]),
    required=True,
)
@click.option("--out", type=_out_path, required=True, help="Output peptides.tsv.")
@click.option("--abundance-out", type=_out_path, help="Output abundance.tsv.")
@click.option("--source", type=click.Choice(["database", "denovo"]))
@click.option("--search-db", type=click.Choice(["same", "other"]))
@shared_options
def convert(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    inputs: tuple[Path, ...],
    engine: str,
    out: Path,
    abundance_out: Path | None,
    source: str | None,
    search_db: str | None,
    **_: Any,
) -> None:
    """Convert search-engine output to the generic peptide table."""
    from giap.io.converters import CONVERTERS
    from giap.io.issues import Issue
    from giap.io.tsv import write_tsv

    cfg = _resolve(config_file, sets, log_level)
    if engine not in CONVERTERS:
        click.echo(f"giap convert --from {engine}: not implemented (not in M1)", err=True)
        sys.exit(ExitCode.INTERNAL)
    if engine == "generic":
        res = CONVERTERS[engine](list(inputs), source=source, search_db=search_db, cfg=cfg.input)
    elif source == "denovo":
        click.echo(
            f"giap convert: error: --source denovo is invalid with --from {engine}", err=True
        )
        sys.exit(ExitCode.INVALID_INPUT)
    else:
        res = CONVERTERS[engine](list(inputs), search_db=search_db, cfg=cfg.input)
    issues = list(res.issues)
    if abundance_out is not None and res.abundance is None:
        issues.append(Issue("no_abundance", "error", f"--from {engine} gives no abundance table"))
    for issue in issues:
        examples = f" (e.g. rows {', '.join(issue.examples)})" if issue.examples else ""
        click.echo(
            f"giap convert: {issue.severity}: {issue.code}: {issue.message}{examples}", err=True
        )
    if any(i.severity == "error" for i in issues):
        sys.exit(ExitCode.INVALID_INPUT)
    peptides = res.peptides
    write_tsv(out, list(peptides.columns), peptides.itertuples(index=False, name=None))
    if abundance_out is not None and res.abundance is not None:
        ab = res.abundance
        write_tsv(abundance_out, list(ab.columns), ab.itertuples(index=False, name=None))
        click.echo(f"wrote {abundance_out}; use --abundance-type maxlfq with it")


# ---------------------------------------------------------------- compare


@cli.command()
@click.argument("runs", nargs=-1, required=True, type=_existing_dir)
@click.option("--outdir", type=_out_path, required=True, help="Comparison output directory.")
@click.option("--labels", help="Comma-separated run labels.")
@shared_options
def compare(config_file: Path | None, sets: tuple[str, ...], log_level: str, **_: Any) -> None:
    """Compare runs on the same genome."""
    _resolve(config_file, sets, log_level)
    _not_implemented("compare", "M9")


# ---------------------------------------------------------------- report


@cli.command()
@click.argument("out", type=_existing_dir)
@click.option("--no-interactive", is_flag=True, help="Static report only.")
@click.option("--plots", help="Comma-separated plot IDs (V01,V02,...).")
@shared_options
def report(config_file: Path | None, sets: tuple[str, ...], log_level: str, **_: Any) -> None:
    """Rebuild plots and report from existing tables."""
    _resolve(config_file, sets, log_level)
    _not_implemented("report", "M9")


# ---------------------------------------------------------------- hidden


@cli.command("_step", hidden=True)
@click.argument("name")
@click.option("--in", "inputs", multiple=True, metavar="KEY=PATH")
@click.option("--out", "outputs", multiple=True, metavar="KEY=PATH")
@shared_options
def step(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    name: str,
    inputs: tuple[str, ...],
    outputs: tuple[str, ...],
    **_: Any,
) -> None:
    """Run one pipeline step (used by Nextflow)."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    cfg = _resolve(config_file, sets, log_level)
    restore = _record_step_events(Path("logs"))  # e.g. placement's search_db event (§6.7)
    try:
        ins = _in_paths(name, inputs)
        outs = _key_paths(name, "--out", outputs)
        reg.run_step(name, cfg, ins, outs)
    except reg.StepError as exc:
        exc.write_log(Path("logs"))
        click.echo(f"giap _step: {exc} (see logs/{name}.log)", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    except reg.WiringError as exc:
        click.echo(f"giap _step: internal error: {exc}", err=True)
        sys.exit(ExitCode.INTERNAL)
    finally:
        restore()


def _key_paths(name: str, flag: str, items: Sequence[str]) -> dict[str, Path]:
    """Parse repeatable ``KEY=PATH`` options; a malformed or repeated key is a wiring error."""
    from giap.steps.registry import WiringError

    out: dict[str, Path] = {}
    for item in items:
        key, sep, path = item.partition("=")
        if not sep or not key or not path or key in out:
            raise WiringError(f"step '{name}': bad or repeated {flag} {item!r}")
        out[key] = Path(path)
    return out


def _in_paths(name: str, items: Sequence[str]) -> dict[str, Path | list[Path]]:
    """Parse ``--in KEY=PATH``; a repeated key collects a list of paths (in the given order)."""
    from giap.steps.registry import WiringError

    out: dict[str, Path | list[Path]] = {}
    for item in items:
        key, sep, path = item.partition("=")
        if not sep or not key or not path:
            raise WiringError(f"step '{name}': bad --in {item!r}")
        prev = out.get(key)
        if prev is None:
            out[key] = Path(path)
        elif isinstance(prev, list):
            prev.append(Path(path))
        else:
            out[key] = [prev, Path(path)]
    return out


@cli.command(hidden=True)
@click.option("--json", "as_json", is_flag=True, help="Machine-readable output.")
def version(as_json: bool) -> None:
    """Print version information."""
    info = {"giap": giap.__version__, "python": platform.python_version()}
    if as_json:
        click.echo(json.dumps(info, sort_keys=True))
    else:
        click.echo(f"giap {info['giap']} (python {info['python']})")


def main(argv: Sequence[str] | None = None) -> None:
    """Console-script entry point (§5.1 exit codes).

    Click usage errors (unknown option, missing file, bad choice) and aborts become exit code 1;
    exit code 2 is reserved for environment errors. Uncaught exceptions become exit code 5.
    """
    args = list(argv) if argv is not None else sys.argv[1:]
    if not args:
        args = ["--help"]
    try:
        # Without standalone mode click returns (not raises) the exit code of --help/--version.
        result = cli.main(args=args, prog_name="giap", standalone_mode=False)
    except click.exceptions.Exit as exc:
        sys.exit(exc.exit_code)
    except click.UsageError as exc:
        exc.show()
        sys.exit(ExitCode.INVALID_INPUT)
    except click.Abort:
        click.echo("Aborted!", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    except Exception:
        logging.getLogger("giap").exception("uncaught exception")
        click.echo(REPORT_MESSAGE, err=True)
        sys.exit(ExitCode.INTERNAL)
    sys.exit(result if isinstance(result, int) else ExitCode.SUCCESS)
