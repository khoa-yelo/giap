"""``giap`` command line (GIAP_SPEC.md §5.1).

Every subcommand resolves the configuration first (an unknown key is exit code 1). ``preflight``
(M1-12), ``_step`` and ``run --executor local`` (M2-02) are implemented; the other commands report
that they are not implemented yet (exit code 5).
"""

from __future__ import annotations

import json
import logging
import platform
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


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
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
def setup(config_file: Path | None, sets: tuple[str, ...], log_level: str, **_: Any) -> None:
    """Download databases and models (writes DIR/db_manifest.json)."""
    _resolve(config_file, sets, log_level)
    _not_implemented("setup", "M2")


# ---------------------------------------------------------------- preflight


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
    type=click.Choice(
        [
            "farmshare,apptainer",
            "farmshare,conda",
            "farmshare_local",
            "apptainer",
            "conda",
            "docker",
        ]
    ),
    default="farmshare,apptainer",
    show_default=True,
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
    }
    cfg = _resolve(config_file, sets, log_level, flags, publication)
    if executor == "local":
        _run_local(cfg, outdir / "work")
    _not_implemented("run", "M2")


def _run_local(cfg: GiapConfig, workdir: Path) -> None:
    """``--executor local`` (§6.4): every registered step as target; exits when done."""
    from giap.runner import local
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    targets = [s.name for s in reg.all_steps()]
    if not targets:
        click.echo("giap run --executor local: not implemented (no steps registered)", err=True)
        sys.exit(ExitCode.INTERNAL)
    try:
        local.plan_local(cfg, targets)
        local.run_local(cfg, workdir, targets)
    except local.LayerMissing as exc:
        click.echo(f"giap run: error: {exc}", err=True)
        sys.exit(ExitCode.INVALID_INPUT)
    except reg.StepError as exc:
        log = exc.write_log(workdir.parent / "logs")
        click.echo(f"giap run: {exc} (see {log})", err=True)
        sys.exit(ExitCode.STEP_FAILED)
    except reg.WiringError as exc:
        click.echo(f"giap run: internal error: {exc}", err=True)
        sys.exit(ExitCode.INTERNAL)
    sys.exit(ExitCode.SUCCESS)


# ---------------------------------------------------------------- makedb


@cli.command()
@click.option("--genome", type=_existing_file, required=True, help="Genome FASTA.")
@click.option("--annotation", type=_existing_file, help="bakta.json.")
@click.option("--outdir", type=_out_path, required=True, help="Output database directory.")
@click.option("--db-dir", type=_out_path, help="Database directory.")
@click.option("--six-frame-orfs", is_flag=True, help="Add six-frame ORFs.")
@click.option("--min-orf-aa", type=click.IntRange(min=1), default=20, show_default=True)
@click.option("--patch", type=_existing_file, help="annotation_patch.gff3 from a previous run.")
@click.option("--no-contaminants", is_flag=True, help="Do not append contaminant proteins.")
@shared_options
def makedb(config_file: Path | None, sets: tuple[str, ...], log_level: str, **_: Any) -> None:
    """Build an optional strain search database (no FDR is built on it)."""
    _resolve(config_file, sets, log_level)
    _not_implemented("makedb", "M3")


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
    from giap.steps import registry as reg

    cfg = _resolve(config_file, sets, log_level)
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
