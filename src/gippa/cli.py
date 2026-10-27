"""``gippa`` command line (GIPPA_SPEC.md §5.1).

Every subcommand resolves the configuration first (an unknown key is exit code 1). ``preflight``
(M1-12), ``_step``, ``run --executor local`` (M2-02), ``run`` through Nextflow (M2-07),
``report`` (M9-02) and ``compare`` (M9-22) are implemented; an engine ``convert`` does not support
yet is exit code 5.
"""

from __future__ import annotations

import json
import logging
import os
import platform
import shutil
import sys
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TypeVar, get_args

import click

import gippa
from gippa.config import GippaConfig, resolve_config
from gippa.config import model as m
from gippa.exitcodes import ExitCode

F = TypeVar("F", bound=Callable[..., Any])

LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")
REPORT_MESSAGE = (
    "gippa: internal error (a bug). Please report it at https://github.com/khoa-yelo/gippa/issues "
    "with the command line and the traceback above."
)

_existing_file = click.Path(exists=True, dir_okay=False, path_type=Path)
_existing_dir = click.Path(exists=True, file_okay=False, path_type=Path)
_out_path = click.Path(path_type=Path)
_PKG = __name__.split(".", 1)[0]
_INSTALL_DOCS = "https://khoa-yelo.github.io/gippa/get-started/install/"


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
    func = click.version_option(gippa.__version__, prog_name="gippa")(func)
    return func


def _resolve(
    config_file: Path | None,
    sets: Sequence[str],
    log_level: str,
    flags: dict[str, Any] | None = None,
    publication: bool = False,
    base: Mapping[str, Any] | None = None,
) -> GippaConfig:
    """Resolve the config: defaults ← ``base`` ← file ← CLI flags ← ``--set``. Exit 1 on an
    invalid key."""
    logging.basicConfig(level=log_level.upper(), format="%(levelname)s %(name)s: %(message)s")
    flag_sets = [
        f"{key}={json.dumps(value)}" for key, value in (flags or {}).items() if value is not None
    ]
    try:
        return resolve_config(config_file, [*flag_sets, *sets], publication=publication, base=base)
    except ValueError as exc:
        click.echo(f"gippa: invalid configuration: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)


class _Cli(click.Group):
    """Keeps the raw argument list (the manifest's ``command_line``; ``sys.argv`` is not it when
    the CLI is invoked in-process)."""

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        ctx.meta["gippa.argv"] = ["gippa", *args]
        return super().parse_args(ctx, args)


@click.group(cls=_Cli, context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(gippa.__version__, prog_name="gippa")
def cli() -> None:
    """GIPPA: genome integrated prokaryotic proteomics analysis."""


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
@click.option(
    "--images",
    is_flag=False,
    flag_value="config",
    default=None,
    metavar="[all]",
    help="Pull the pinned container images the config needs (with 'all': every pinned image).",
)
@click.option(
    "--profile",
    default=None,
    metavar="apptainer|docker",
    help="Container runtime for --images (a comma-joined Nextflow profile is accepted).",
)
@click.option(
    "--cache-dir",
    type=_out_path,
    default=None,
    help="Apptainer image cache and UniOP build context [default: NXF_APPTAINER_CACHEDIR, "
    "$GIPPA_ROOT/containers, DIR/containers].",
)
@click.option("--uniop", is_flag=True, help="Fetch UniOP from its authors and build it locally.")
@click.option(
    "--accept-uniop-terms",
    is_flag=True,
    help="UniOP has no licence file: you decide that your use is permitted (needed without a TTY).",
)
@click.option("--no-uniop", is_flag=True, help="Do not build UniOP under --images.")
@shared_options
def setup(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    images: str | None,
    profile: str | None,
    cache_dir: Path | None,
    uniop: bool,
    accept_uniop_terms: bool,
    no_uniop: bool,
    **opts: Any,
) -> None:
    """Download databases and models (writes DIR/db_manifest.json)."""
    from gippa.annotate.tools import required_tools
    from gippa.io.image_key import image_cache
    from gippa.preflight.environment import container_runtime
    from gippa.setupdb import HttpsFetcher, SetupOptions, ToolRunner, run_setup
    from gippa.setupdb.images import image_set, install_images
    from gippa.setupdb.uniop_build import setup_uniop

    if uniop and no_uniop:
        raise click.UsageError("--uniop and --no-uniop exclude each other")
    cfg = _resolve(config_file, sets, log_level)
    runtime = container_runtime(profile or "")
    if images not in (None, "config", "all"):
        raise click.BadParameter(
            f"{images!r}: give --images alone or --images all", param_hint="--images"
        )
    if images is not None and runtime is None:
        raise click.UsageError("--images needs --profile apptainer or --profile docker")
    db_dir: Path = opts["db_dir"]
    check = bool(opts["check"])
    fetcher = HttpsFetcher()
    # under --check the archive is only verified; else the UniOP step fetches it (spec item 1b).
    # An explicit --uniop installs UniOP only (Q-K27): no database is planned.
    dbs = SetupOptions(db_dir, bakta="skip", genomad=False, promotech_models=False, check=check,
                       uniop=check) if uniop else SetupOptions(**opts)  # fmt: skip
    # one cache for the pull, the downloaders that run in an image, preflight and Nextflow
    cache = cache_dir or image_cache(os.environ, db_dir)
    code = ExitCode.SUCCESS
    if images is not None and runtime is not None and not check:  # before the downloaders
        tools = image_set(cfg, images == "all")
        if install_images(tools, runtime, cache, db_dir) != ExitCode.SUCCESS:
            code = ExitCode.ENVIRONMENT
    dbs_code = run_setup(dbs, ToolRunner(db_dir=db_dir, cache_dir=cache), fetcher)
    if dbs_code != ExitCode.SUCCESS:
        code = dbs_code
    needed = images is not None and not no_uniop and "uniop" in {
        t.name for t in required_tools(cfg)}  # fmt: skip
    if (uniop or needed) and not check:
        step = setup_uniop(profile or "", db_dir, cache, fetcher, accepted=accept_uniop_terms,
                           explicit=uniop)  # fmt: skip
        if step != ExitCode.SUCCESS:
            code = ExitCode.ENVIRONMENT
    sys.exit(int(code))


# ---------------------------------------------------------------- preflight

# Nextflow profiles `gippa run --profile` accepts, alone or comma-joined (#78 review 3c).
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


def _check_optional_profile(
    ctx: click.Context, param: click.Parameter, value: str | None
) -> str | None:
    """:func:`_check_profile` for an option that may be left out."""
    return None if value is None else _check_profile(ctx, param, value)


@cli.command()
@click.option("--genome", type=_existing_file, required=True, help="Genome FASTA.")
@click.option("--proteins", type=_existing_file, help="Protein FASTA.")
@click.option("--peptides", type=_existing_file, help="Peptide table (TSV).")
@click.option("--abundance", type=_existing_file, help="Protein abundance table (TSV).")
@click.option("--annotation", type=_existing_file, help="bakta.json(.gz), GFF3 or GenBank.")
@click.option("--outdir", type=_out_path, help="Output directory.")
@click.option("--strict", is_flag=True, help="Treat warnings as errors (exit code 3).")
@click.option(
    "--profile",
    callback=_check_optional_profile,
    default=None,
    help="Also check the tools a run under these profiles needs, as 'gippa run' does (DB "
    "directory: $GIPPA_DB_DIR). Comma-joined from: " + ", ".join(_PROFILES) + ".",
)
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
    profile: str | None,
    **_: Any,
) -> None:
    """Validate inputs only (seconds); with --profile also the installed tools."""
    from gippa.preflight.report import render_preflight_md
    from gippa.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    cfg = _resolve(config_file, sets, log_level)
    if outdir is not None and outdir.exists() and not outdir.is_dir():
        click.echo(f"gippa preflight: --outdir {outdir} exists and is not a directory", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    inp = PreflightInputs(genome, proteins, peptides, abundance, annotation)
    if profile is None:
        report, norm = run_preflight(inp, cfg, None, "preflight", outdir)
    else:  # the tool check of `gippa run --profile` (same DB directory and $GIPPA_ROOT)
        env_db, root = os.environ.get("GIPPA_DB_DIR"), os.environ.get("GIPPA_ROOT")
        report, norm = run_preflight(
            inp, cfg, None, "preflight", outdir, Path(env_db) if env_db else None,
            gippa_root=Path(root) if root else None, profile=profile,
        )  # fmt: skip
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
    "peptide_spectra": ("reference_peptides", "peptide_spectra.parquet"),  # C17 (M7-10)
}


@cli.command()
@click.option("--genome", type=_existing_file, required=True, help="Genome FASTA.")
@click.option("--peptides", type=_existing_file, help="Peptide table (TSV).")
@click.option("--proteins", type=_existing_file, help="Protein FASTA.")
@click.option("--abundance", type=_existing_file, help="Protein abundance table (TSV).")
@click.option("--annotation", type=_existing_file, help="bakta.json(.gz), GFF3 or GenBank.")
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
@click.option(
    "--publication",
    is_flag=True,
    help="10,000 permutations and bootstrap draws (also the determinants null); the run takes "
    "about 10x longer.",
)
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
    from gippa.config import resolve_sample_name

    # §5.2 sample_name: auto is resolved once, before any step (M7-04 lead decision 1)
    sample = resolve_sample_name(cfg, _.get("peptides"), _.get("proteins"))
    cfg = cfg.model_copy(update={"sample_name": sample})
    from gippa.annotate.promoters import check_promoter_mode

    try:  # before launching anything (M2-13 Q2, revised lead decision 4)
        check_promoter_mode(cfg.annotation.promoters.mode)
    except ValueError as exc:
        click.echo(f"gippa run: error: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    from gippa.io.manifest_sources import rotate_events

    rotate_events(outdir.absolute() / "logs")  # one events.jsonl per run (MVP-03 D-B step 1)
    env_db = os.environ.get("GIPPA_DB_DIR")
    db_dir: Path | None = _.get("db_dir") or (Path(env_db) if env_db else None)
    given = {k: _.get(k) for k in _RUN_INPUTS}
    run_info = {"publication": publication, "db_dir": db_dir, "given": given}
    if executor == "local":
        keys = ("peptides", "proteins", "abundance", "annotation", "reference_peptides", "db_dir")
        # no profile, so no container runtime: the function step must not see the shell's (14c)
        with _env_unset("GIPPA_CONTAINER_RUNTIME"):
            _run_local(cfg, outdir, _["genome"], {k: _.get(k) for k in keys}, run_info)
    inputs = {k: _[k] for k in _RUN_INPUTS if _.get(k) is not None}
    outdir = outdir.absolute()
    roots, cfg = _run_preflight(
        cfg, outdir, inputs, db_dir, bool(_["strict"]), str(_["profile"]), publication
    )
    inputs.update(roots)
    from gippa.runner.local import contaminant_library

    lib = contaminant_library(db_dir)  # M3-15 library for MAPPING's loci
    if lib is not None:
        inputs["contaminants"] = lib
    if shutil.which("nextflow") is None:
        click.echo(
            f"gippa run: nextflow not found on PATH; install Nextflow and Java ({_INSTALL_DOCS}), "
            f"then pull the images with '{_PKG} setup --images --profile <profile>'",
            err=True,
        )
        sys.exit(ExitCode.ENVIRONMENT)
    base = _launch_sources(cfg, outdir, run_info, "nextflow", str(_["profile"]), bool(_["resume"]))
    try:
        code = _run_nextflow(
            cfg, outdir, inputs, str(_["memory"]), str(_["profile"]), bool(_["resume"]), db_dir
        )
        if code == ExitCode.SUCCESS:
            _finalize(outdir)
    except Exception:  # past the provenance core: exit 5, the tail still runs (review B2)
        code = _internal_error(outdir)
    from gippa.runner.nextflow import failed_step

    trace = outdir / "logs" / "trace.tsv"
    hit = failed_step(trace) if code != ExitCode.SUCCESS and trace.is_file() else None
    failed = {"failed_step": hit[0], "failed_log": f"{hit[1]}/.command.log"} if hit else {}
    sys.exit(_run_tail(outdir, base, code, **failed))


def _launch_sources(
    cfg: GippaConfig, outdir: Path, info: dict[str, Any], executor: str, profile: str | None,
    resume: bool,
) -> dict[str, Any]:  # fmt: skip
    """Launch-time manifest sources; writes ``OUT/work/provenance.json`` (MVP-03 D-B step 3)."""
    from gippa.steps.manifest import launch_sources, write_provenance

    meta = click.get_current_context().meta
    argv = meta.get("gippa.argv", sys.argv)
    base = launch_sources(cfg, outdir.absolute(), info["given"], executor=executor,
                          profile=profile, resume=resume, publication=info["publication"],
                          db_dir=info["db_dir"], command_line=argv,
                          demo=bool(meta.get("gippa.demo")))  # fmt: skip
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


def _finalize(outdir: Path) -> None:
    """Edits of published files that need several steps' outputs, for either executor:
    ``novel_candidate_id`` in ``mapping/peptide_hits.tsv`` (M7-15 item 6)."""
    from gippa.steps.novel import backfill_candidate_ids

    backfill_candidate_ids(outdir)


def _run_tail(outdir: Path, base: dict[str, Any], code: int, **kw: Any) -> int:
    """The run tail (MVP-03 D-B steps 5-6): the manifest, even after a failed step (with one
    warning line when the operons came from the distance model, M11-14d), then the
    ``report.html`` re-render when the manifest says ``completed`` (M9-16b, ADR-0040 amendment
    1). A tail error after a successful graph is exit 4; after a failed graph it is only logged."""
    import traceback

    from gippa.steps.manifest import run_tail

    outdir = outdir.absolute()
    try:
        run_tail(outdir, base, int(code), **kw)
    except Exception:
        log = outdir / "logs" / "manifest.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(traceback.format_exc(), encoding="utf-8")
        if code == ExitCode.SUCCESS:
            click.echo(f"gippa run: step manifest failed; log: {log}", err=True)
            return ExitCode.STEP_FAILED
        click.echo(f"gippa run: the run manifest was not written; log: {log}", err=True)
        return int(code)
    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("operon_source") == "distance_model":  # never silent (ADR-0061 decision 5)
        why = manifest.get("operon_source_reason") or "reason not recorded"
        click.echo(f"gippa run: warning: operons used the distance model: {why}", err=True)
    if manifest.get("status") != "completed":
        return int(code)
    from gippa.steps import plots as P

    try:  # the page reads the manifest just written (M9-16b); nothing writes OUT/ after it
        P.rerender_report(outdir, base["cfg"])
    except Exception:
        log = outdir / "logs" / "report_tail.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(traceback.format_exc(), encoding="utf-8")
        click.echo(f"gippa run: report re-render failed; log: {log}", err=True)
        return ExitCode.STEP_FAILED
    return int(code)


def _is_bakta_json(annotation: Path) -> bool:
    """True when ``annotation`` is named like a Bakta JSON: ``.json`` or ``.json.gz``, any case."""
    return Path(annotation).name.lower().endswith((".json", ".json.gz"))


def _bakta_json(annotation: Path | None) -> str | None:
    """The absolute path of ``--annotation`` when it is a Bakta JSON (``.json`` or ``.json.gz``,
    passed as given), else None (GFF3/GBFF go through Bakta ``--regions``, §5.3)."""
    if annotation is None or not _is_bakta_json(annotation):
        return None
    return str(Path(annotation).absolute())


def _run_preflight(
    cfg: GippaConfig,
    outdir: Path,
    inputs: dict[str, Path],
    db_dir: Path | None,
    strict: bool,
    profile: str,
    publication: bool = False,
) -> tuple[dict[str, Path], GippaConfig]:
    """Preflight in run mode (§5.1): exits on an error (or a warning with ``strict``) before
    Nextflow starts, else returns the normalised ``genome`` and ``replicons`` for params.json
    (plus ``peptides`` and ``proteins_long`` when given, for MAPPING) and the configuration
    with ``annotation.uniop.enabled`` resolved (M11-14d).
    Tools are checked for the Nextflow ``profile`` under ``$GIPPA_ROOT`` (run-wiring 4a)."""
    from gippa.preflight.report import render_preflight_md
    from gippa.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    if outdir.exists() and not outdir.is_dir():
        click.echo(f"gippa run: --outdir {outdir} exists and is not a directory", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    inp = PreflightInputs(
        genome=inputs["genome"],
        proteins=inputs.get("proteins"),
        peptides=inputs.get("peptides"),
        abundance=inputs.get("abundance"),
        annotation=inputs.get("annotation"),
        reference_peptides=inputs.get("reference_peptides"),
    )
    root = os.environ.get("GIPPA_ROOT")
    report, norm = run_preflight(
        inp, cfg, None, "run", outdir, db_dir,
        gippa_root=Path(root) if root else None, profile=profile, publication=publication,
    )  # fmt: skip
    outdir.mkdir(parents=True, exist_ok=True)
    write_outputs(outdir, report, norm)
    code = exit_code(report, strict)
    if code != ExitCode.SUCCESS:
        click.echo(render_preflight_md(report), err=True, nl=False)
        sys.exit(code)
    for c in report.checks:  # the UniOP fallback is never silent (ADR-0061 decision 5)
        if c.rule_id == "environment.tools" and c.level == "warn":
            click.echo(f"gippa run: warning: environment.tools: {c.message}", err=True)
    work = outdir / "work" / "preflight"
    roots = {"genome": work / "genome.norm.fa", "replicons": work / "replicons.tsv"}
    for key, (raw, name) in _NORM_INPUTS.items():
        if raw in inputs and (work / name).is_file():
            roots[key] = work / name
    report_json = outdir / "preflight" / "preflight.json"  # gene_table's abundance_scale (M1-12)
    if report_json.is_file():
        roots["preflight"] = report_json
    return roots, norm.cfg or cfg


def _run_nextflow(
    cfg: GippaConfig,
    outdir: Path,
    inputs: dict[str, Path],
    memory: str,
    profile: str,
    resume: bool,
    db_dir: Path | None = None,
) -> ExitCode:
    """§6.3 "Launch": params.json (+ ``db_dir``) + resolved_config.yaml, then ``nextflow run``."""
    from gippa.runner import nextflow as nf

    outdir = outdir.absolute()
    from gippa.preflight.environment import LOCAL_TOOLS, local_tools_applies

    params, _resolved = nf.write_params(outdir, cfg, inputs, nf.nextflow_memory(memory), db_dir)
    from gippa.preflight.environment import container_runtime

    # a local UniOP build (M11-14); absolute, since Nextflow runs in outdir (#396 review B1)
    local = (db_dir / LOCAL_TOOLS).absolute() if db_dir is not None else None
    if local is None or not local.is_file() or not local_tools_applies(profile):
        local = None
    cmd = nf.launch_command(nf.workflow_main(), params, profile, outdir, resume, local)
    env = nf.launch_env(os.environ, container_runtime(profile), db_dir)
    code = nf.run_nextflow(cmd, outdir, env)
    if code == ExitCode.SUCCESS:
        from gippa.runner.local import write_datapackage

        write_datapackage(outdir)
    return code


def _local_roots(
    cfg: GippaConfig, outdir: Path, genome: Path, inputs: dict[str, Path | None]
) -> dict[str, Path]:
    """Preflight-normalised root work files for the local executor (§6.2, §6.4).

    A preflight error (e.g. no evidence) exits with the preflight exit code; a GFF3/GenBank
    ``--annotation`` without ``annotation.precomputed.bakta_json`` exits 1 naming the key. The
    contaminant library (M3-15), the COG 2020 table (M4-15) and, with the ``function`` module,
    its installed resources (M8-06) come from ``--db-dir`` / ``$GIPPA_DB_DIR``; preflight's
    ``abundance`` and ``preflight.json`` are staged when written.
    """
    import os

    from gippa.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs
    from gippa.runner.local import contaminant_library
    from gippa.setupdb.installed import FUNCTION_RESOURCES, cog_table, resource

    inp = PreflightInputs(
        genome,
        inputs["proteins"],
        inputs["peptides"],
        inputs.get("abundance"),
        inputs["annotation"],
        inputs.get("reference_peptides"),
    )
    report, norm = run_preflight(inp, cfg, None, "preflight", outdir)
    write_outputs(outdir, report, norm)
    code = exit_code(report, False)
    if code != ExitCode.SUCCESS:
        for c in report.checks:
            if c.level == "error":
                click.echo(f"gippa run: preflight error {c.rule_id}: {c.message}", err=True)
        sys.exit(code)
    work = outdir / "work" / "preflight"
    files = {"genome": "genome.norm.fa", "peptides": "peptides.norm.parquet",
             "peptide_spectra": "peptide_spectra.parquet"}  # fmt: skip
    roots = {k: work / f for k, f in files.items() if (work / f).is_file()}
    if (work / "proteins.long.faa").is_file():
        roots["proteins_long"] = work / "proteins.long.faa"
    ann, ref = inputs["annotation"], inputs["reference_peptides"]
    is_json = ann is not None and _is_bakta_json(ann)
    bakta = ann if is_json else cfg.annotation.precomputed.bakta_json
    if ann is not None and bakta is None:  # a GFF3/GenBank annotation cannot replace Bakta JSON
        click.echo("gippa run: error: annotation.precomputed.bakta_json is required with "
                   "--executor local (--annotation is not a Bakta JSON)", err=True)  # fmt: skip
        sys.exit(ExitCode.INVALID_INPUT)
    if bakta is not None:
        roots["bakta_json"] = Path(bakta)
    if ref is not None:
        roots["reference_peptides"] = ref
    db = inputs["db_dir"] or os.environ.get("GIPPA_DB_DIR")
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


@contextmanager
def _env_unset(name: str) -> Iterator[None]:
    """``os.environ`` without ``name`` for the block; the caller's value is restored after it
    (also on ``SystemExit``)."""
    saved = os.environ.pop(name, None)
    try:
        yield
    finally:
        if saved is not None:
            os.environ[name] = saved


def _run_local(
    cfg: GippaConfig,
    outdir: Path,
    genome: Path,
    inputs: dict[str, Path | None],
    info: dict[str, Any],
) -> None:
    """``--executor local`` (§6.4): every registered step as target, then the run tail; exits
    when done. A layer or tool problem found while planning exits before the provenance core."""
    from gippa.runner import local
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    workdir = outdir / "work"
    load_step_modules()
    targets = [s.name for s in reg.all_steps()]
    if not targets:
        click.echo("gippa run --executor local: not implemented (no steps registered)", err=True)
        sys.exit(ExitCode.INTERNAL)
    roots = _local_roots(cfg, outdir, genome, inputs)
    from gippa.preflight.environment import resolve_uniop

    # no tool runs here: ``auto`` resolves as preflight does without a profile (M11-14d)
    cfg = resolve_uniop(cfg, installed=True, publication=bool(info["publication"])).cfg
    try:
        local.plan_local(cfg, targets, set(roots))
    except local.LayerMissing as exc:
        click.echo(f"gippa run: error: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    except local.ToolUnavailable as exc:
        click.echo(f"gippa run: error: {exc}", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    base = _launch_sources(cfg, outdir, info, "local", None, False)
    timings: list[dict[str, Any]] = []
    code, failed = ExitCode.SUCCESS, {}
    # step events (e.g. the ``search_db: auto`` resolution) go to OUT/logs/events.jsonl (§6.7)
    restore = _record_step_events(outdir / "logs")
    try:
        # the live events log, so plots/report see this run's layer_missing flags (#205 B6)
        live = {"events": outdir / "logs" / "events.jsonl"}
        local.run_local(cfg, workdir, targets, roots, timings, live=live)
        local.publish(workdir, outdir)
        _finalize(outdir)
        local.write_datapackage(outdir)
    except reg.StepError as exc:
        log = exc.write_log(workdir.parent / "logs")
        click.echo(f"gippa run: {exc} (see {log})", err=True)
        code, failed = ExitCode.STEP_FAILED, {"failed_step": exc.name, "failed_log": str(log)}
    except reg.WiringError as exc:
        click.echo(f"gippa run: internal error: {exc}", err=True)
        code = ExitCode.INTERNAL
    except Exception:  # past the provenance core: exit 5, the tail still runs (review B2)
        code = _internal_error(outdir)
    finally:
        restore()
    sys.exit(_run_tail(outdir, base, code, local_steps=timings, **failed))


def _record_step_events(logs_dir: Path) -> Callable[[], None]:
    """Send ``gippa.steps`` events to ``logs_dir/events.jsonl`` (§6.7); returns the undo."""
    from gippa.io.events import EventsHandler

    steps_log = logging.getLogger("gippa.steps")
    events = EventsHandler(logs_dir)
    steps_log.addHandler(events)
    prev_level = steps_log.level
    if steps_log.getEffectiveLevel() > logging.INFO:
        steps_log.setLevel(logging.INFO)

    def restore() -> None:
        steps_log.removeHandler(events)
        steps_log.setLevel(prev_level)

    return restore


# ---------------------------------------------------------------- demo

_DEMO_MODULES = ("core", "mobile", "bgc")


@cli.command()
@click.option(
    "--outdir", type=_out_path, default="gippa-demo", show_default=True,
    help="Output directory; must not exist.",
)  # fmt: skip
@click.option(
    "--full", is_flag=True, default=False,
    help="After 'gippa setup': the same inputs through 'gippa run' with the real tools, the "
    "databases and Nextflow (checks the whole install).",
)  # fmt: skip
@click.option(
    "--open", "open", is_flag=True, default=False, help="Open the report in a web browser."
)
def demo(outdir: Path, full: bool, **kw: Any) -> None:
    """Run the bundled synthetic SIM100K demo (no Nextflow, containers, databases or network)."""
    import time

    if os.path.lexists(outdir):
        click.echo(f"gippa demo: {outdir} exists; use --outdir or remove it", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    created = _first_missing(outdir)  # what this call creates; removed if a --full run fails
    try:
        outdir.mkdir(parents=True)
    except OSError as exc:
        click.echo(f"gippa demo: cannot create {outdir} ({exc.strerror or exc}); choose a "
                   "writable directory with --outdir", err=True)  # fmt: skip
        sys.exit(ExitCode.INVALID_INPUT)
    ctx = click.get_current_context()
    ctx.meta["gippa.demo"] = True  # the manifest's ``demo: true`` and the report banner
    start = time.monotonic()
    code: int | str = 0
    try:
        (_demo_full if full else _demo_local)(ctx, outdir.absolute())
    except SystemExit as exc:
        code = exc.code if exc.code is not None else 0
    if code != 0:
        if full:  # typically preflight (missing setup): a kept outdir would refuse the retry;
            shutil.rmtree(created, ignore_errors=True)  # the local demo keeps its logs there
        sys.exit(code)
    report = (outdir / "report" / "report.html").resolve()
    click.echo(f"gippa demo: done in {round(time.monotonic() - start)} s; report: {report}")
    if kw["open"]:
        import webbrowser

        try:  # no browser (a headless host) is not an error: the path is printed above
            webbrowser.open(report.as_uri())
        except Exception as exc:
            click.echo(f"gippa demo: could not open a browser ({exc})", err=True)


def _first_missing(path: Path) -> Path:
    """The outermost directory ``path.mkdir(parents=True)`` would create (``path`` itself when
    its parent exists)."""
    path = path.absolute()
    while not os.path.lexists(path.parent) and path.parent != path:
        path = path.parent
    return path


def _demo_data() -> Path:
    """The bundled demo directory ``gippa/data/demo`` (package data, on disk)."""
    from importlib.resources import files

    return Path(str(files("gippa.data") / "demo"))


def _demo_local(ctx: click.Context, outdir: Path) -> None:
    """``gippa run --executor local`` on the bundled inputs and layers (ADR-0066): the layers are
    gunzipped into ``OUT/work/demo_layers``; no ``$GIPPA_DB_DIR``; numba caches under ``OUT``."""
    import gzip

    data, layers = _demo_data(), outdir / "work" / "demo_layers"
    layers.mkdir(parents=True)
    bakta = layers / "bakta.json"
    bakta.write_bytes(gzip.decompress((data / "layers" / "bakta.json.gz").read_bytes()))
    sets = []
    for key in m.PrecomputedConfig.model_fields:
        if key != "bakta_json":
            path = layers / f"{key}.tsv"
            path.write_bytes(gzip.decompress((data / "layers" / f"{key}.tsv.gz").read_bytes()))
            sets.append(f"annotation.precomputed.{key}={json.dumps(str(path))}")
    with _env_unset("GIPPA_DB_DIR"), _numba_cache(outdir / "work" / "numba_cache"):
        ctx.invoke(run, genome=data / "sim100k.fna.gz", peptides=data / "sim100k_peptides.tsv.gz",
                   annotation=bakta, modules=",".join(_DEMO_MODULES), outdir=outdir,
                   executor="local", sets=tuple(sets))  # fmt: skip


def _demo_full(ctx: click.Context, outdir: Path) -> None:
    """``gippa run`` (Nextflow, the profile it resolves itself) on the bundled genome and
    peptides with the real tools; adds ``function`` when its resources are installed."""
    from gippa.setupdb.installed import FUNCTION_RESOURCES, resource

    env_db = os.environ.get("GIPPA_DB_DIR")
    db = Path(env_db) if env_db else None
    modules = list(_DEMO_MODULES)
    if db is not None and all(resource(db, n) is not None for n in FUNCTION_RESOURCES):
        modules.append("function")
    data = _demo_data()
    ctx.invoke(run, genome=data / "sim100k.fna.gz", peptides=data / "sim100k_peptides.tsv.gz",
               modules=",".join(modules), outdir=outdir)  # fmt: skip


@contextmanager
def _numba_cache(path: Path) -> Iterator[None]:
    """``NUMBA_CACHE_DIR=path`` for the block unless the caller set it (numba's ``cache=True``
    must not write into a read-only site-packages); removed again after it."""
    if "NUMBA_CACHE_DIR" in os.environ:
        yield
        return
    os.environ["NUMBA_CACHE_DIR"] = str(path)
    try:
        yield
    finally:
        os.environ.pop("NUMBA_CACHE_DIR", None)


# ---------------------------------------------------------------- makedb


@cli.command()
@click.option("--annotation", type=_existing_file, help="bakta.json.")
@click.option("--genome", type=_existing_file, help="Genome FASTA (default: the JSON sequences).")
@click.option("--outdir", type=_out_path, required=True, help="Output database directory.")
@click.option("--db-dir", type=_out_path, help="Database directory (default: $GIPPA_DB_DIR).")
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

    from gippa.makedb.build import MakedbOptions, build_strain_db

    _resolve(config_file, sets, log_level)
    if annotation is None:
        click.echo("gippa makedb: no --annotation: Bakta wrapper required (M2)", err=True)
        sys.exit(ExitCode.ENVIRONMENT)
    if min_orf_aa < 1:
        click.echo(f"gippa makedb: --min-orf-aa must be >= 1, got {min_orf_aa}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    contaminants: Path | None = None
    if not no_contaminants:
        env = os.environ.get("GIPPA_DB_DIR")
        root = db_dir or (Path(env) if env else None)
        contaminants = root / "contaminants" / "contaminants.fasta" if root else None
        if contaminants is None or not contaminants.is_file():
            click.echo(
                f"gippa makedb: contaminant library not found ({contaminants or 'no --db-dir'}); "
                "run gippa setup --contaminants or pass --no-contaminants",
                err=True,
            )
            sys.exit(ExitCode.ENVIRONMENT)
    opts = MakedbOptions(
        annotation, genome, outdir, patch, six_frame_orfs, min_orf_aa, contaminants
    )
    try:
        result = build_strain_db(opts)
    except (ValueError, OSError) as exc:
        click.echo(f"gippa makedb: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    click.echo(f"gippa makedb: wrote {result.fasta} {dict(result.counts)}", err=True)
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
    from gippa.io.converters import CONVERTERS
    from gippa.io.issues import Issue
    from gippa.io.tsv import write_tsv

    cfg = _resolve(config_file, sets, log_level)
    if engine not in CONVERTERS:
        click.echo(f"gippa convert --from {engine}: not implemented (not in M1)", err=True)
        sys.exit(ExitCode.INTERNAL)
    if engine == "generic":
        res = CONVERTERS[engine](list(inputs), source=source, search_db=search_db, cfg=cfg.input)
    elif source is not None and source != ("denovo" if engine == "cascadia" else "database"):
        click.echo(
            f"gippa convert: error: --source {source} is invalid with --from {engine}", err=True
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
            f"gippa convert: {issue.severity}: {issue.code}: {issue.message}{examples}", err=True
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
@click.option(
    "--allow-cross-major",
    is_flag=True,
    help="Compare runs whose output_schema_version majors differ (kept as a warning).",
)
@shared_options
def compare(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    runs: tuple[Path, ...],
    outdir: Path,
    labels: str | None,
    allow_cross_major: bool,
    **_: Any,
) -> None:
    """Compare runs on the same genome (CMP/tables/compare_*.tsv, V16, report/compare.html)."""
    import pandas as pd

    from gippa.compare.abundance import abundance_pairs
    from gippa.compare.detection import detection_matrix, intersections
    from gippa.compare.stability import stability
    from gippa.compare.versions import compare_versions, cross_major
    from gippa.io.compare_io import (
        CompareInputError,
        load_runs,
        read_manifests,
        write_tables,
        write_version_events,
    )
    from gippa.rng import for_step

    cfg = _resolve(config_file, sets, log_level)
    try:
        loaded = load_runs(list(runs), labels)
    except CompareInputError as exc:
        click.echo(f"gippa compare: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    labs, genes = loaded.labels, loaded.genes
    manifests = read_manifests(list(runs), labs)
    major_error = cross_major(manifests)
    if major_error is not None and not allow_cross_major:
        click.echo(major_error, err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    versions = compare_versions(manifests)
    for line in versions:
        click.echo(line, err=True)
    write_version_events(outdir, versions)
    matrix = detection_matrix(genes)
    rng = for_step("compare", cfg.seed)
    pairs = []
    for i, la in enumerate(labs):  # decision 7: the TUs of the pair's first run
        tus = pd.Series(genes[la]["tu_id"].to_numpy(), index=genes[la]["locus_tag"].to_numpy())
        for lb in labs[i + 1 :]:
            pairs.append(abundance_pairs({la: genes[la], lb: genes[lb]}, tus,
                                         cfg.stats.bootstrap, rng))  # fmt: skip
    tables = {
        "compare_detection": pd.concat([matrix, intersections(matrix, labs)], ignore_index=True),
        "compare_abundance": pd.concat(pairs, ignore_index=True),
        "compare_stability": stability(loaded.views, labs, cfg),
    }
    written = write_tables(outdir, tables)
    for path in [*written, *_compare_report(outdir, list(runs), loaded, tables, cfg, versions)]:
        click.echo(str(path))


def _compare_report(
    cmp: Path,
    paths: list[Path],
    loaded: Any,
    tables: dict[str, Any],
    cfg: Any,
    versions: Sequence[str] = (),
) -> list[Path]:
    """V16 (``CMP/plots``) and ``CMP/report/compare.html`` (M9-23); returns the paths written."""
    import os

    import pandas as pd

    from gippa.io.compare_io import genome_digest
    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.report.compare import RunInfo, build_compare
    from gippa.viz import style
    from gippa.viz.compare import compare_plot
    from gippa.viz.interactive import figure_json

    style.apply()
    inst = compare_plot(tables["compare_detection"], tables["compare_abundance"], loaded.genes)
    plots = cmp / "plots"
    written = write_instance(inst, "V16", "compare", plots, cfg.report.formats, cfg.report.dpi)
    twin = None
    if cfg.report.interactive:
        data = pd.read_csv(plots / "data" / "V16_compare.tsv", sep="\t", na_values=["NA"],
                           keep_default_na=False)  # fmt: skip
        text = figure_json("V16", data, inst.title)
        twin = None if text is None else ("plotly-V16", text)
    runs = []
    for lab, path in zip(loaded.labels, paths, strict=True):
        tiers = sorted({str(t) for t in loaded.genes[lab]["detection_tier"].dropna()})
        runs.append(RunInfo(lab, os.path.relpath(path.resolve(), cmp.resolve()),
                            genome_digest(path) or "NA", ", ".join(tiers) or "NA"))  # fmt: skip
    page = cmp / "report" / "compare.html"
    page.parent.mkdir(parents=True, exist_ok=True)
    html = build_compare(runs, read_plot_files(plots).get("V16"), tables, twin, versions)
    page.write_text(html, encoding="utf-8")
    return [*written, page]


# ---------------------------------------------------------------- report


def _manifest_config(out: Path) -> Mapping[str, Any] | None:
    """The ``config`` that ``OUT/run_manifest.json`` records, when it is a valid configuration;
    ``None`` otherwise (no or unreadable manifest, a config this version cannot validate)."""
    try:
        recorded = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    config = recorded.get("config") if isinstance(recorded, dict) else None
    if not isinstance(config, dict):
        return None
    try:
        resolve_config(base=config)
    except ValueError:
        logging.getLogger(__name__).warning("report.manifest_config_invalid %s", out)
        return None
    return config


@cli.command()
@click.argument("out", type=_existing_dir)
@click.option("--no-interactive", is_flag=True, help="Static report only.")
@click.option("--plots", help="Comma-separated plot IDs (V01,V02,...).")
@click.option("--strict", is_flag=True, help="Exit 3 when the report exceeds its size budget.")
@shared_options
def report(
    config_file: Path | None,
    sets: tuple[str, ...],
    log_level: str,
    out: Path,
    no_interactive: bool,
    plots: str | None,
    strict: bool,
    **_: Any,
) -> None:
    """Rebuild plots and report from existing tables."""
    from gippa.io.run_view import NotAGippaRun, RunView
    from gippa.steps import plots as P
    from gippa.viz.base import normalise_ids

    base = None
    if config_file is None:  # the run's own thresholds, formats, dpi, modules (review B3)
        runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
        config_file = next((p for p in runs if p.is_file()), None)
        if config_file is None:  # a local run writes none; its manifest records the config
            base = _manifest_config(out)
    cfg = _resolve(config_file, sets, log_level, base=base)
    if no_interactive:
        cfg = cfg.model_copy(
            update={"report": cfg.report.model_copy(update={"interactive": False})}
        )
    P.load_plot_modules()
    try:
        view = RunView.from_outdir(out)
        only = None if plots is None else normalise_ids(plots)
    except (NotAGippaRun, ValueError) as exc:
        click.echo(f"gippa report: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    try:
        P.render_plots(view, cfg, out / "plots", only)
    except P.PlotFailed as exc:
        click.echo(f"gippa report: {exc}", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    try:
        page = P.rerender_report(out, cfg)
    except Exception as exc:  # unreadable tables are notes; anything else is exit 4 (B5)
        click.echo(f"gippa report: {type(exc).__name__}: {exc}", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    click.echo(str(page))
    if strict and page.stat().st_size > P.REPORT_BUDGET:  # report.size_exceeded (M9-14)
        sys.exit(ExitCode.WARNINGS_STRICT)


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
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    cfg = _resolve(config_file, sets, log_level)
    restore = _record_step_events(Path("logs"))  # e.g. placement's search_db event (§6.7)
    try:
        ins = _in_paths(name, inputs)
        outs = _key_paths(name, "--out", outputs)
        reg.run_step(name, cfg, ins, outs)
    except reg.StepError as exc:
        exc.write_log(Path("logs"))
        click.echo(f"gippa _step: {exc} (see logs/{name}.log)", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    except reg.WiringError as exc:
        click.echo(f"gippa _step: internal error: {exc}", err=True)
        sys.exit(ExitCode.INTERNAL)
    finally:
        restore()


def _key_paths(name: str, flag: str, items: Sequence[str]) -> dict[str, Path]:
    """Parse repeatable ``KEY=PATH`` options; a malformed or repeated key is a wiring error."""
    from gippa.steps.registry import WiringError

    out: dict[str, Path] = {}
    for item in items:
        key, sep, path = item.partition("=")
        if not sep or not key or not path or key in out:
            raise WiringError(f"step '{name}': bad or repeated {flag} {item!r}")
        out[key] = Path(path)
    return out


def _in_paths(name: str, items: Sequence[str]) -> dict[str, Path | list[Path]]:
    """Parse ``--in KEY=PATH``; a repeated key collects a list of paths (in the given order)."""
    from gippa.steps.registry import WiringError

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
    info = {"gippa": gippa.__version__, "python": platform.python_version()}
    if as_json:
        click.echo(json.dumps(info, sort_keys=True))
    else:
        click.echo(f"gippa {info['gippa']} (python {info['python']})")


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
        result = cli.main(args=args, prog_name="gippa", standalone_mode=False)
    except click.exceptions.Exit as exc:
        sys.exit(exc.exit_code)
    except click.UsageError as exc:
        exc.show()
        sys.exit(ExitCode.INVALID_INPUT)
    except click.Abort:
        click.echo("Aborted!", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    except Exception:
        logging.getLogger("gippa").exception("uncaught exception")
        click.echo(REPORT_MESSAGE, err=True)
        sys.exit(ExitCode.INTERNAL)
    sys.exit(result if isinstance(result, int) else ExitCode.SUCCESS)
