"""The R fgsea engine: ``workflow/bin/giap_fgsea.R`` as a subprocess (M8-17a; ADR-0038).

R (GPL-2) and fgsea (MIT) run only as a separate process in the ``giap-fgsea`` micromamba env
(R5); nothing here imports R bindings. One call writes ``stats.tsv``, ``jobs.tsv`` and ``sets.tsv``
to a fresh temporary directory, starts the script there, checks ``versions.txt`` (fgsea 1.36.2)
and parses ``results.tsv`` against what was sent. Every failure raises ``FgseaEngineError``
(``FgseaEngineUnavailable`` when the env cannot be started); the engine never retries and the
temporary directory is always removed. Conda-only under every profile (lead decision 9).
"""

from __future__ import annotations

import contextlib
import math
import os
import re
import shutil
import signal
import subprocess
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from giap.annotate.tools import TOOLS
from giap.function.gsea import EngineCall, EngineProvenance, EngineResult, EngineRow, GseaJob
from giap.preflight.environment import _env_built
from giap.runner.nextflow import workflow_main
from giap.steps.registry import StepError

__all__ = [
    "DATA_TABLE_VERSION",
    "ENV_NAME",
    "FGSEA_TOOL",
    "FGSEA_VERSION",
    "METHOD",
    "RESULTS_HEADER",
    "SCRIPT",
    "TIMEOUT_S",
    "FgseaEngineError",
    "FgseaEngineUnavailable",
    "Launcher",
    "RFgseaEngine",
    "child_env",
    "command",
    "r_fgsea_engine",
    "read_results",
    "read_versions",
    "resolve_launcher",
    "script_path",
    "write_call",
]

FGSEA_VERSION = "1.36.2"
SCRIPT = "giap_fgsea.R"  # resolved as workflow_main().parent / "bin" / SCRIPT
TIMEOUT_S = 3600.0  # per call (lead decision 6)
METHOD = "fgsea_multilevel"
ENV_NAME = "giap-fgsea"
RESULTS_HEADER = ("job_id", "pathway", "size", "ES", "NES", "pval", "log2err", "leadingEdge")
STDERR_TAIL = 20
DATA_TABLE_VERSION = "1.18.6.1"
FGSEA_TOOL = TOOLS["fgsea"]  # the env's registry entry, used for the ``.giap-built`` test
_NOT_BUILT = (
    f"the R fgsea environment {ENV_NAME} is not built: sbatch scripts/slurm/build_envs.sbatch fgsea"
)
_R_NUMBER = re.compile(r"^-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$")
_SIZE = re.compile(r"^\d+$")


class FgseaEngineError(StepError):
    """The R fgsea call failed; the step fails with exit code 4 (lead decision 5)."""

    def __init__(self, message: str) -> None:
        Exception.__init__(self, message)
        self.name = "function"
        self.cause: BaseException = self
        self.traceback = ""


class FgseaEngineUnavailable(FgseaEngineError):
    """The ``giap-fgsea`` env is missing or cannot be started."""


@dataclass(frozen=True, slots=True)
class Launcher:
    """How to start the script: a micromamba env prefix (conda-only, lead decision 9)."""

    kind: Literal["prefix"]
    path: Path


# ---------------------------------------------------------------- launcher


def script_path() -> Path:
    """Absolute path of ``workflow/bin/giap_fgsea.R`` in this checkout."""
    return (workflow_main().parent / "bin" / SCRIPT).resolve()


def _front_end() -> str:
    """The interpreter named by the script's ``#!/usr/bin/env -S <front end> ...`` shebang."""
    try:
        first = script_path().read_text(encoding="utf-8").split("\n", 1)[0]
    except OSError as exc:
        raise FgseaEngineUnavailable(f"cannot read {script_path()}: {exc}") from exc
    words = first.removeprefix("#!").split()
    if "-S" not in words or words.index("-S") + 1 >= len(words):
        raise FgseaEngineUnavailable(f"{script_path()}: unexpected shebang {first!r}")
    return words[words.index("-S") + 1]


def _has_r(prefix: Path) -> bool:
    return (prefix / "bin" / _front_end()).is_file()


def resolve_launcher(environ: Mapping[str, str]) -> Launcher:
    """``GIAP_FGSEA_PREFIX``, else ``$GIAP_ROOT/micromamba/envs/giap-fgsea`` when its
    ``.giap-built`` matches ``workflow/envs/fgsea.yml``; the prefix must hold R (lead
    decision 6). Otherwise ``FgseaEngineUnavailable``."""
    explicit = environ.get("GIAP_FGSEA_PREFIX", "")
    if explicit:
        prefix = Path(explicit)
        if not prefix.is_dir():
            raise FgseaEngineUnavailable(f"GIAP_FGSEA_PREFIX={explicit} is not a directory")
        if not _has_r(prefix):
            raise FgseaEngineUnavailable(f"GIAP_FGSEA_PREFIX={explicit} has no R in bin/")
        return Launcher("prefix", prefix)
    root = environ.get("GIAP_ROOT", "")
    if root:
        envs = Path(root) / "micromamba" / "envs"
        prefix = envs / ENV_NAME
        if _env_built(FGSEA_TOOL, envs, workflow_main().parent / "envs") and _has_r(prefix):
            return Launcher("prefix", prefix)
    raise FgseaEngineUnavailable(_NOT_BUILT)


def command(call: EngineCall) -> list[str]:
    """The script and its nine options (spec §3); file names are relative to the call dir."""
    p = call.params
    return [
        str(script_path()),
        "--stats", "stats.tsv",
        "--sets", "sets.tsv",
        "--jobs", "jobs.tsv",
        "--sample-size", str(int(p.sample_size)),
        "--eps", repr(float(p.eps)),
        "--n-perm-simple", str(int(p.n_perm_simple)),
        "--seed", str(int(call.seed)),
        "--out", "results.tsv",
        "--versions", "versions.txt",
    ]  # fmt: skip


def child_env(launcher: Launcher, environ: Mapping[str, str]) -> dict[str, str]:
    """``environ`` with the env's ``bin`` first on ``PATH``, user and site R libraries neutralised
    and ``LC_ALL=C``, so no other R library can shadow the pinned fgsea."""
    env = {k: v for k, v in environ.items() if k != "R_HOME"}
    path, bin_dir = environ.get("PATH", ""), str(launcher.path / "bin")
    env["PATH"] = f"{bin_dir}{os.pathsep}{path}" if path else bin_dir
    # R_LIBS_SITE must be non-empty: the env's Renviron reads ${R_LIBS_SITE:-'%S'}, and a site
    # library sits before .Library in the search path.
    nul = "/dev/null"
    env.update(R_LIBS_USER=nul, R_LIBS=nul, R_LIBS_SITE=nul, R_PROFILE_USER=nul, LC_ALL="C")
    return env


# ---------------------------------------------------------------- files


def write_call(call: EngineCall, d: Path) -> None:
    """``stats.tsv`` (``repr`` scores, call order), ``jobs.tsv`` and long-form ``sets.tsv``."""
    scores = [float(s) for s in call.scores]
    if not all(math.isfinite(s) for s in scores):
        raise ValueError("non-finite GSEA score; nothing is sent to R")
    stats = ["tu_id\tscore", *(f"{t}\t{s!r}" for t, s in zip(call.tu_ids, scores, strict=True))]
    jobs = ["job_id\tmin_size\tmax_size"]
    sets = ["job_id\tset_id\tmember"]
    for job in call.jobs:
        jobs.append(f"{job.job_id}\t{int(job.min_size)}\t{int(job.max_size)}")
        sets.extend(f"{job.job_id}\t{term}\t{m}" for term, ms in job.sets.items() for m in ms)
    for name, lines in (("stats.tsv", stats), ("jobs.tsv", jobs), ("sets.tsv", sets)):
        (d / name).write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


def _lines(path: Path, what: str) -> list[str]:
    try:
        text = path.read_bytes().decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise FgseaEngineError(f"R fgsea wrote no readable {what}: {exc}") from exc
    if not text.endswith("\n"):
        raise FgseaEngineError(f"R fgsea {what} does not end with a newline")
    return text[:-1].split("\n")


def read_versions(path: Path) -> EngineProvenance:
    """``versions.txt``; raises unless it names fgsea ``FGSEA_VERSION`` and data.table
    ``DATA_TABLE_VERSION`` and carries integer ``seed`` and ``jobs`` lines (§3)."""
    keys: dict[str, str] = {}
    warnings: list[tuple[str, str]] = []
    for line in _lines(path, "versions.txt"):
        key, _, rest = line.partition("\t")
        if key == "warning":
            job, _, msg = rest.partition("\t")
            if (job, msg) not in warnings:
                warnings.append((job, msg))
        elif key:
            keys.setdefault(key, rest)
    fgsea = keys.get("fgsea")
    if fgsea != FGSEA_VERSION:
        raise FgseaEngineError(f"R fgsea version {fgsea!r} is not the pinned {FGSEA_VERSION}")
    table = keys.get("data.table")
    if table != DATA_TABLE_VERSION:
        raise FgseaEngineError(f"R data.table {table!r} is not the pinned {DATA_TABLE_VERSION}")
    counts = {k: keys.get(k, "") for k in ("seed", "jobs")}
    for k, v in counts.items():
        if not _SIZE.match(v):
            raise FgseaEngineError(f"R fgsea versions.txt: {k} {v!r} is not an integer")
    return EngineProvenance(r=keys.get("R", ""), fgsea=fgsea, warnings=tuple(warnings),
                            seed=int(counts["seed"]), jobs=int(counts["jobs"]))  # fmt: skip


def _number(cell: str, what: str) -> float:
    if cell == "NA":
        return math.nan
    if not _R_NUMBER.match(cell):
        raise FgseaEngineError(f"R fgsea results: {what} {cell!r} is not a number")
    return float(cell)


def _row(fields: list[str], job: GseaJob, ranked: frozenset[str]) -> EngineRow:
    _, term, size, es, nes, pval, err, edge_cell = fields
    members = job.sets[term]
    distinct = set(members)
    if not _SIZE.match(size) or int(size) != len(distinct):
        raise FgseaEngineError(
            f"R fgsea results: {job.job_id} {term} size {size!r}, sent {len(distinct)} TUs"
        )
    edge = tuple(edge_cell.split(";")) if edge_cell else ()
    if len(set(edge)) != len(edge) or not set(edge) <= distinct or not set(edge) <= ranked:
        raise FgseaEngineError(
            f"R fgsea results: {job.job_id} {term} leading edge is not a subset of its set"
        )
    return EngineRow(
        term=term,
        size=int(size),
        es=_number(es, "ES"),
        nes=_number(nes, "NES"),
        p=_number(pval, "pval"),
        log2err=_number(err, "log2err"),
        leading_edge=edge,
    )


def read_results(path: Path, call: EngineCall) -> EngineResult:
    """``results.tsv`` → rows per sent job (every job is a key). Rows of a job must be in-range
    sets of that job in the sent order, sizes must equal the sent sizes and every leading-edge TU
    must be in its set; every in-range set must have its row (#204 review nit 1). Any violation
    raises ``FgseaEngineError``."""
    lines = _lines(path, "results.tsv")
    if tuple(lines[0].split("\t")) != RESULTS_HEADER:
        raise FgseaEngineError(f"R fgsea results: unexpected header {lines[0]!r}")
    jobs = {j.job_id: (i, j) for i, j in enumerate(call.jobs)}
    order = {j.job_id: {t: k for k, t in enumerate(j.sets)} for j in call.jobs}
    ranked = frozenset(call.tu_ids)
    rows: dict[str, list[EngineRow]] = {j.job_id: [] for j in call.jobs}
    last = (-1, -1)
    for n, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != len(RESULTS_HEADER):
            raise FgseaEngineError(f"R fgsea results line {n}: {len(fields)} fields")
        job_id, term = fields[0], fields[1]
        if job_id not in jobs:
            raise FgseaEngineError(f"R fgsea results line {n}: unknown job {job_id!r}")
        i, job = jobs[job_id]
        if term not in order[job_id]:
            raise FgseaEngineError(f"R fgsea results line {n}: {term!r} is not a set of {job_id}")
        k = len(set(job.sets[term]))
        if not job.min_size <= k <= job.max_size:
            raise FgseaEngineError(f"R fgsea results line {n}: {term} is outside the size range")
        pos = (i, order[job_id][term])
        if pos <= last:
            raise FgseaEngineError(f"R fgsea results line {n}: rows are not in the sent order")
        last = pos
        rows[job_id].append(_row(fields, job, ranked))
    for job in call.jobs:
        want = [t for t, ms in job.sets.items() if job.min_size <= len(set(ms)) <= job.max_size]
        if [r.term for r in rows[job.job_id]] != want:
            raise FgseaEngineError(
                f"R fgsea results: {job.job_id} has {len(rows[job.job_id])} of its "
                f"{len(want)} in-range sets"
            )
    return {j: tuple(r) for j, r in rows.items()}


# ---------------------------------------------------------------- the engine


def _tail(stderr: bytes) -> str:
    lines = stderr.decode("utf-8", errors="replace").rstrip("\n").split("\n")
    return "\n".join(lines[-STDERR_TAIL:])


def _kill_group(proc: subprocess.Popen[bytes]) -> None:
    with contextlib.suppress(ProcessLookupError):
        os.killpg(proc.pid, signal.SIGKILL)


class RFgseaEngine:
    """``GseaEngine`` backed by ``giap_fgsea.R``; ``provenance`` is that of the last call."""

    def __init__(self, launcher: Launcher, timeout_s: float, tmp_root: Path | None) -> None:
        self.method = METHOD
        self.launcher = launcher
        self.timeout_s = timeout_s
        self.tmp_root = tmp_root
        self.provenance: EngineProvenance | None = None

    def __call__(self, call: EngineCall) -> EngineResult:
        tmp = Path(tempfile.mkdtemp(prefix="giap-fgsea-", dir=self.tmp_root))
        try:
            write_call(call, tmp)
            self._run(call, tmp)
            prov = read_versions(tmp / "versions.txt")
            if (prov.seed, prov.jobs) != (call.seed, len(call.jobs)):
                raise FgseaEngineError(
                    f"R fgsea ran seed {prov.seed} with {prov.jobs} jobs; sent seed {call.seed} "
                    f"with {len(call.jobs)} jobs"
                )
            result = read_results(tmp / "results.tsv", call)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        self.provenance = prov
        return result

    def _run(self, call: EngineCall, cwd: Path) -> None:
        try:
            proc = subprocess.Popen(
                command(call),
                cwd=cwd,
                env=child_env(self.launcher, os.environ),
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                start_new_session=True,  # own process group: a timeout kills R and its children
            )
        except (FileNotFoundError, PermissionError) as exc:
            raise FgseaEngineUnavailable(f"cannot start R fgsea: {exc}") from exc
        try:
            _, stderr = proc.communicate(timeout=self.timeout_s)
        except subprocess.TimeoutExpired:
            _kill_group(proc)
            proc.communicate()
            raise FgseaEngineError(f"R fgsea timed out after {self.timeout_s:g} s") from None
        except BaseException:
            _kill_group(proc)
            proc.communicate()
            raise
        if proc.returncode != 0:
            raise FgseaEngineError(
                f"R fgsea failed (exit {proc.returncode}, {len(call.jobs)} jobs): {_tail(stderr)}"
            )


def r_fgsea_engine(
    launcher: Launcher, *, timeout_s: float = TIMEOUT_S, tmp_root: Path | None = None
) -> RFgseaEngine:
    """The R engine for ``gsea_table`` (``method = "fgsea_multilevel"``)."""
    return RFgseaEngine(launcher, timeout_s, tmp_root)
