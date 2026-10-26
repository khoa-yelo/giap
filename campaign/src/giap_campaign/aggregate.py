"""Pool the per-run GIAP findings over datasets and search arms (campaign SPEC §11).

For each job of the plan (§10) and each arm (``db`` → ``giap_db/``, ``denovo`` → ``giap_denovo/``)
whose ``run_manifest.json`` exists, the GIAP tables are read **as text** (``dtype=str``,
``keep_default_na=False``, no quoting) so per-run values such as ``true``, ``1.5e-05`` or
``unestimated`` pass through unchanged; GIAP's missing spelling ``NA`` becomes an empty cell (§1).
Rows are prefixed with the job's ``dataset``, ``tax_id``, ``organism``, ``genome_accession`` and the
``arm`` and sorted by dataset, numeric tax_id, arm (``db`` before ``denovo``), then the table's own
row order. Per-run q-values and FDRs are copied unchanged; nothing is pooled statistically (R4): the
only cross-run columns are counts of arms and datasets that report the same ORF or variant.
"""

from __future__ import annotations

import csv
import json
import math
import numbers
import re
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

__all__ = [
    "ARMS",
    "HYPOTHETICAL_COLUMNS",
    "OUTPUTS",
    "PREFIX_COLUMNS",
    "RUNS_COLUMNS",
    "aggregate",
]

#: Search arms in output order (``db`` before ``denovo``).
ARMS: tuple[str, ...] = ("db", "denovo")
PREFIX_COLUMNS: tuple[str, ...] = ("dataset", "tax_id", "organism", "genome_accession", "arm")
HYPOTHETICAL_DETAIL: tuple[str, ...] = (
    "locus_tag",
    "product",
    "contig",
    "start",
    "end",
    "strand",
    "length_aa",
    "detection_tier",
    "gene_q_value",
    "n_peptides_unique",
    "abundance_log2",
)
HYPOTHETICAL_COLUMNS: tuple[str, ...] = (*PREFIX_COLUMNS, *HYPOTHETICAL_DETAIL)
RUNS_COLUMNS: tuple[str, ...] = (
    *PREFIX_COLUMNS,
    "state",
    "acquisition",
    "giap_version",
    "n_genes_detected",
    "n_novel_called",
    "n_alt_high",
    "notes",
)

NOVEL_TSV = "campaign_novel_orfs.tsv"
NOVEL_FAA = "campaign_novel_orfs.faa"
VARIANTS_TSV = "campaign_variants.tsv"
OPERON_TSV = "campaign_operon_flags.tsv"
ISLANDS_TSV = "campaign_islands.tsv"
HYPOTHETICAL_TSV = "campaign_detected_hypothetical.tsv"
RUNS_TSV = "campaign_runs.tsv"
#: Every output, in the order of the returned row counts.
OUTPUTS: tuple[str, ...] = (
    NOVEL_TSV,
    NOVEL_FAA,
    VARIANTS_TSV,
    OPERON_TSV,
    ISLANDS_TSV,
    HYPOTHETICAL_TSV,
    RUNS_TSV,
)

_JOB_COLUMNS = ("dataset", "tax_id", "organism", "genome_accession", "workdir")
_GIAP_NA = "NA"
_FORBIDDEN = re.compile(r"[\t\r\n]")

# Source tables (below ``giap_<arm>/tables/``), the arms they are read for and the columns the
# aggregation needs (a table without them is noted as unreadable).
_NOVEL = "novel_orfs.tsv"
_FAA = "novel_orfs.faa"
_VARIANTS = "variants.tsv"
_OPERON = "operon_flags.tsv"
_ISLANDS = "islands.tsv"
_GENES = "genes.tsv"
_REQUIRED: dict[str, tuple[str, ...]] = {
    _NOVEL: ("orf_id", "called"),
    _VARIANTS: ("contig", "codon_start", "strand", "alt_aa", "call", "confidence"),
    _OPERON: (),
    _ISLANDS: (),
    _GENES: ("locus_tag", "hypothetical", "detected"),
}
_TABLES_BY_ARM: dict[str, tuple[str, ...]] = {
    "db": (_NOVEL, _VARIANTS, _OPERON, _ISLANDS, _GENES),
    "denovo": (_NOVEL, _VARIANTS, _GENES),
}
_VARIANT_KEY = ("contig", "codon_start", "strand", "alt_aa")

Row = dict[str, str]


@dataclass(frozen=True, slots=True)
class _Table:
    columns: tuple[str, ...]
    rows: tuple[Row, ...]


@dataclass(slots=True)
class _Run:
    """One job and arm: its prefix, sort key and what was read from its GIAP output."""

    prefix: Row
    sort_key: tuple[Any, ...]
    job_key: tuple[str, str, str]
    state: str
    acquisition: str
    giap_version: str = ""
    has_manifest: bool = False
    tables: dict[str, _Table] = field(default_factory=dict)
    fasta: list[str] | None = None
    notes: list[str] = field(default_factory=list)


# --- small pure helpers ---------------------------------------------------------------------------


def _text(value: object) -> str:
    """A job cell as text: missing → ``""``; integral floats lose their ``.0``."""
    if value is None or value is pd.NA or value is pd.NaT:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, numbers.Integral):
        return str(int(value))
    if isinstance(value, numbers.Real):
        f = float(value)
        if math.isnan(f):
            return ""
        return str(int(f)) if f.is_integer() else repr(f)
    return str(value)


def _tax_key(tax_id: str) -> tuple[int, int, str]:
    """Numeric tax_id order; non-numeric values after the numeric ones, by text."""
    try:
        return (0, int(tax_id), "")
    except ValueError:
        return (1, 0, tax_id)


def _arm_index(arm: str) -> int:
    return ARMS.index(arm)


def _is_true(cell: str) -> bool:
    return cell.strip().lower() == "true"


def _clean(cell: str) -> str:
    """A cell safe for a TSV without quoting (tabs and line breaks become spaces)."""
    return _FORBIDDEN.sub(" ", cell)


# --- reading --------------------------------------------------------------------------------------


def _read_table(path: Path) -> _Table:
    """A GIAP TSV read as text; ``NA`` → ``""``. Raises ``OSError``/``ValueError`` if unreadable."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        quoting=csv.QUOTE_NONE,
        encoding="utf-8",
        engine="c",
    )
    columns = tuple(str(c) for c in df.columns)
    rows: list[Row] = []
    for record in df.itertuples(index=False, name=None):
        row: Row = {}
        for col, value in zip(columns, record, strict=True):
            cell = "" if value is None or not isinstance(value, str) else value
            row[col] = "" if cell == _GIAP_NA else cell
        rows.append(row)
    return _Table(columns, tuple(rows))


def _read_fasta_lines(path: Path) -> list[str]:
    """The lines of a FASTA file (blank lines dropped); ``ValueError`` if it is not FASTA."""
    lines = [ln.rstrip("\r\n") for ln in path.read_text(encoding="utf-8").splitlines()]
    lines = [ln for ln in lines if ln.strip()]
    if lines and not lines[0].startswith(">"):
        raise ValueError("sequence before the first header")
    return lines


def _read_json_object(path: Path) -> dict[str, Any]:
    """A JSON object; ``OSError``/``ValueError`` when unreadable or not an object."""
    doc = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(doc, dict):
        raise ValueError("not a JSON object")
    return doc


def _first_failed_step(steps: object) -> str | None:
    if not isinstance(steps, Mapping):
        return None
    for name, info in steps.items():
        if isinstance(info, Mapping) and info.get("state") == "failed":
            return str(name)
    return None


def _step_failed(steps: object, name: str) -> bool:
    if not isinstance(steps, Mapping):
        return False
    info = steps.get(name)
    return isinstance(info, Mapping) and info.get("state") == "failed"


def _job_runs(job: Mapping[str, object]) -> list[_Run]:
    """Both arms of one job, with the tables of every arm that has a ``run_manifest.json``."""
    prefix_base = {
        "dataset": _text(job["dataset"]),
        "tax_id": _text(job["tax_id"]),
        "organism": _text(job["organism"]),
        "genome_accession": _text(job["genome_accession"]),
    }
    workdir_text = _text(job["workdir"])
    workdir = Path(workdir_text) if workdir_text else None
    job_key = (prefix_base["dataset"], prefix_base["tax_id"], prefix_base["genome_accession"])

    status: dict[str, Any] = {}
    status_note: str | None = None
    if workdir is not None and (workdir / "status.json").is_file():
        try:
            status = _read_json_object(workdir / "status.json")
        except (OSError, ValueError):
            status_note = "unreadable:status.json"
    steps = status.get("steps")
    first_failed = _first_failed_step(steps)
    acq = status.get("acquisition")
    acquisition = acq if isinstance(acq, str) else ""

    runs: list[_Run] = []
    for arm in ARMS:
        prefix = {**prefix_base, "arm": arm}
        sort_key = (
            prefix["dataset"],
            _tax_key(prefix["tax_id"]),
            prefix["genome_accession"],
            workdir_text,
            _arm_index(arm),
        )
        out = workdir / f"giap_{arm}" if workdir is not None else None
        manifest = out / "run_manifest.json" if out is not None else None
        has_manifest = manifest is not None and manifest.is_file()
        if _step_failed(steps, f"giap_{arm}"):
            state = f"failed:giap_{arm}"
        elif has_manifest:
            state = "done"
        elif first_failed is not None:
            state = f"failed:{first_failed}"
        else:
            state = "not_run"
        run = _Run(prefix, sort_key, job_key, state, acquisition, has_manifest=has_manifest)
        if status_note is not None:
            run.notes.append(status_note)
        if has_manifest and out is not None and manifest is not None:
            _load_arm(run, arm, out, manifest)
        runs.append(run)
    return runs


def _load_arm(run: _Run, arm: str, out: Path, manifest: Path) -> None:
    """Read the manifest version and the arm's tables into ``run`` (unreadable ones noted)."""
    try:
        version = _read_json_object(manifest).get("giap_version")
        run.giap_version = version if isinstance(version, str) else ""
    except (OSError, ValueError):
        run.notes.append("unreadable:run_manifest.json")
    tables = out / "tables"
    for name in _TABLES_BY_ARM[arm]:
        path = tables / name
        if not path.is_file():
            continue
        try:
            table = _read_table(path)
        except (OSError, ValueError):
            run.notes.append(f"unreadable:tables/{name}")
            continue
        missing = [c for c in _REQUIRED[name] if c not in table.columns]
        if missing:
            run.notes.append(f"missing_columns:tables/{name}:{','.join(missing)}")
            continue
        run.tables[name] = table
    faa = tables / _FAA
    if faa.is_file():
        try:
            run.fasta = _read_fasta_lines(faa)
        except (OSError, ValueError):
            run.notes.append(f"unreadable:tables/{_FAA}")


# --- pooling (pure) -------------------------------------------------------------------------------


def _merge_columns(tables: Iterable[_Table], extras: Sequence[str]) -> list[str]:
    """Prefix, the union of the source columns in first-seen order, then ``extras``."""
    reserved = {*PREFIX_COLUMNS, *extras}
    seen: dict[str, None] = {}
    for table in tables:
        for col in table.columns:
            if col not in reserved:
                seen.setdefault(col, None)
    return [*PREFIX_COLUMNS, *seen, *extras]


def _novel_rows(runs: Sequence[_Run]) -> tuple[list[str], list[Row]]:
    arms_called: dict[tuple[tuple[str, str, str], str], set[str]] = defaultdict(set)
    datasets_called: dict[tuple[str, str], set[str]] = defaultdict(set)
    for run in runs:
        table = run.tables.get(_NOVEL)
        if table is None:
            continue
        for row in table.rows:
            if _is_true(row["called"]):
                arms_called[(run.job_key, row["orf_id"])].add(run.prefix["arm"])
                key = (run.prefix["genome_accession"], row["orf_id"])
                datasets_called[key].add(run.prefix["dataset"])
    present = [r for r in runs if _NOVEL in r.tables]
    extras = ("arms_called", "n_datasets_called")
    columns = _merge_columns((r.tables[_NOVEL] for r in present), extras)
    out: list[Row] = []
    for run in present:
        for row in run.tables[_NOVEL].rows:
            arms = arms_called.get((run.job_key, row["orf_id"]), set())
            n = datasets_called.get((run.prefix["genome_accession"], row["orf_id"]), set())
            out.append(
                {
                    **row,
                    **run.prefix,
                    "arms_called": ";".join(sorted(arms, key=_arm_index)),
                    "n_datasets_called": str(len(n)),
                }
            )
    return columns, out


def _variant_rows(runs: Sequence[_Run]) -> tuple[list[str], list[Row]]:
    datasets: dict[tuple[str, ...], set[str]] = defaultdict(set)
    present = [r for r in runs if _VARIANTS in r.tables]
    for run in present:
        for row in run.tables[_VARIANTS].rows:
            key = (run.prefix["genome_accession"], *(row[c] for c in _VARIANT_KEY))
            datasets[key].add(run.prefix["dataset"])
    columns = _merge_columns((r.tables[_VARIANTS] for r in present), ("n_datasets",))
    out: list[Row] = []
    for run in present:
        for row in run.tables[_VARIANTS].rows:
            key = (run.prefix["genome_accession"], *(row[c] for c in _VARIANT_KEY))
            out.append({**row, **run.prefix, "n_datasets": str(len(datasets[key]))})
    return columns, out


def _db_rows(runs: Sequence[_Run], name: str) -> tuple[list[str], list[Row]]:
    present = [r for r in runs if r.prefix["arm"] == "db" and name in r.tables]
    columns = _merge_columns((r.tables[name] for r in present), ())
    rows = [{**row, **run.prefix} for run in present for row in run.tables[name].rows]
    return columns, rows


def _hypothetical_rows(runs: Sequence[_Run]) -> list[Row]:
    out: list[Row] = []
    for run in runs:
        table = run.tables.get(_GENES)
        if run.prefix["arm"] != "db" or table is None:
            continue
        for row in table.rows:
            if _is_true(row["hypothetical"]) and _is_true(row["detected"]):
                out.append({**run.prefix, **{c: row.get(c, "") for c in HYPOTHETICAL_DETAIL}})
    return out


def _count(run: _Run, name: str, predicate: Callable[[Row], bool]) -> str:
    table = run.tables.get(name)
    if table is None:
        return ""
    return str(sum(1 for row in table.rows if predicate(row)))


def _runs_rows(runs: Sequence[_Run]) -> list[Row]:
    out: list[Row] = []
    for run in runs:
        out.append(
            {
                **run.prefix,
                "state": run.state,
                "acquisition": run.acquisition,
                "giap_version": run.giap_version,
                "n_genes_detected": _count(run, _GENES, lambda r: _is_true(r["detected"])),
                "n_novel_called": _count(run, _NOVEL, lambda r: _is_true(r["called"])),
                "n_alt_high": _count(
                    run,
                    _VARIANTS,
                    lambda r: r["call"] == "ALT" and r["confidence"] == "high",
                ),
                "notes": ";".join(run.notes),
            }
        )
    return out


def _fasta_lines(runs: Sequence[_Run]) -> tuple[list[str], int]:
    lines: list[str] = []
    n = 0
    for run in runs:
        if run.fasta is None:
            continue
        p = run.prefix
        tag = f"{p['dataset']}|{p['arm']}|{p['genome_accession']}|"
        for line in run.fasta:
            if line.startswith(">"):
                lines.append(">" + tag + line[1:])
                n += 1
            else:
                lines.append(line)
    return lines, n


# --- writing --------------------------------------------------------------------------------------


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Row]) -> int:
    lines = ["\t".join(columns)]
    n = 0
    for row in rows:
        lines.append("\t".join(_clean(row.get(c, "")) for c in columns))
        n += 1
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return n


def aggregate(jobs: pd.DataFrame, out_dir: Path | str) -> dict[str, int]:
    """Write the campaign tables (§11) for ``jobs`` (the plan) into ``out_dir``.

    Returns the rows written per output (FASTA: records); every output is written, header-only
    when empty.
    """
    missing = [c for c in _JOB_COLUMNS if c not in jobs.columns]
    if missing:
        raise ValueError(f"jobs table lacks column(s): {', '.join(missing)}")
    records: list[dict[str, Any]] = [
        {str(k): v for k, v in rec.items()} for rec in jobs.to_dict("records")
    ]
    runs = [run for job in records for run in _job_runs(job)]
    runs.sort(key=lambda r: r.sort_key)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}

    columns, rows = _novel_rows(runs)
    counts[NOVEL_TSV] = _write_tsv(out / NOVEL_TSV, columns, rows)

    fasta, n_records = _fasta_lines(runs)
    with (out / NOVEL_FAA).open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("".join(line + "\n" for line in fasta))
    counts[NOVEL_FAA] = n_records

    columns, rows = _variant_rows(runs)
    counts[VARIANTS_TSV] = _write_tsv(out / VARIANTS_TSV, columns, rows)

    for output, source in ((OPERON_TSV, _OPERON), (ISLANDS_TSV, _ISLANDS)):
        columns, rows = _db_rows(runs, source)
        counts[output] = _write_tsv(out / output, columns, rows)

    counts[HYPOTHETICAL_TSV] = _write_tsv(
        out / HYPOTHETICAL_TSV, HYPOTHETICAL_COLUMNS, _hypothetical_rows(runs)
    )
    counts[RUNS_TSV] = _write_tsv(out / RUNS_TSV, RUNS_COLUMNS, _runs_rows(runs))
    return counts
