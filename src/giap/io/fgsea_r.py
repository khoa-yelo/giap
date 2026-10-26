"""The R fgsea engine: ``workflow/bin/giap_fgsea.R`` as a subprocess (M8-17a; ADR-0038).

R (GPL-2) and fgsea (MIT) run only as a separate process (R5): in the ``giap-fgsea`` micromamba
env, or in the ``giap-fgsea`` image through ``docker run`` / ``apptainer exec`` (M11-14 14c);
nothing here imports R bindings or a container SDK. One call writes ``stats.tsv``, ``jobs.tsv``
and ``sets.tsv`` to a fresh temporary directory, starts the script there, checks
``versions.txt`` (fgsea 1.36.2) and parses ``results.tsv`` against what was sent. Every failure
raises ``FgseaEngineError`` (``FgseaEngineUnavailable`` when R cannot be started); the engine
never retries and the temporary directory is always removed. The launcher is
``GIAP_FGSEA_PREFIX`` (never falls through), else the prebuilt env, else the image under
``GIAP_CONTAINER_RUNTIME`` (set by ``giap run``).
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
from giap.io.image_key import image_cache
from giap.preflight.environment import _env_built, _image_file_name, _step_images
from giap.runner.nextflow import workflow_main
from giap.steps.registry import StepError

__all__ = [
    "BAKED_SCRIPT",
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
_PKG = __name__.split(".", 1)[0]
BAKED_SCRIPT = f"/opt/{_PKG}/{SCRIPT}"  # in the image (workflow/containers/fgsea/Dockerfile)
RUNTIME_ENV = "GIAP_CONTAINER_RUNTIME"
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
    """How to start the script: a micromamba env ``prefix`` (``path`` is the env), or the image
    under ``docker`` (``image`` is the ref, pulled by ``docker run``) or ``apptainer`` (``path``
    is the cached ``.img``, ``image`` its ref)."""

    kind: Literal["prefix", "docker", "apptainer"]
    path: Path | None
    image: str = ""


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


def _image_ref() -> str | None:
    """``images.json["fgsea"]["ref"]`` of the workflow tree, or ``None`` without an entry."""
    return _step_images(workflow_main().parent / "conf" / "images.json").get("fgsea")


def _image_launcher(environ: Mapping[str, str]) -> Launcher:
    """The image under ``GIAP_CONTAINER_RUNTIME`` (``docker`` | ``apptainer``), else
    ``FgseaEngineUnavailable`` with the unbuilt-env message."""
    runtime = environ.get(RUNTIME_ENV)
    if runtime is None:
        raise FgseaEngineUnavailable(_NOT_BUILT)
    if runtime not in ("docker", "apptainer"):
        raise FgseaEngineUnavailable(f"{_NOT_BUILT}; {RUNTIME_ENV}={runtime!r} is not a runtime")
    ref = _image_ref()
    if ref is None:
        raise FgseaEngineUnavailable(f"{_NOT_BUILT}; workflow/conf/images.json pins no fgsea image")
    if runtime == "docker":
        return Launcher("docker", None, ref)
    keys = ("NXF_APPTAINER_CACHEDIR", "GIAP_ROOT", "GIAP_DB_DIR")
    if not any(environ.get(k) for k in keys):
        raise FgseaEngineUnavailable(f"{_NOT_BUILT}; no apptainer image cache ({', '.join(keys)})")
    img = image_cache(environ, Path(environ.get("GIAP_DB_DIR", ""))) / _image_file_name(ref)
    if not img.is_file():
        raise FgseaEngineUnavailable(f"{_NOT_BUILT}; image {ref} not in the apptainer cache {img}")
    return Launcher("apptainer", img, ref)


def resolve_launcher(environ: Mapping[str, str]) -> Launcher:
    """In order (M11-14 item 1): ``GIAP_FGSEA_PREFIX``, which never falls through (it must be a
    directory holding R); ``$GIAP_ROOT/micromamba/envs/giap-fgsea`` when its ``.giap-built``
    matches ``workflow/envs/fgsea.yml``; the image under ``GIAP_CONTAINER_RUNTIME``. Otherwise
    ``FgseaEngineUnavailable``."""
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
    return _image_launcher(environ)


def _options(call: EngineCall) -> list[str]:
    """The nine options (spec §3); file names are relative to the call dir."""
    p = call.params
    return [
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


def command(
    call: EngineCall, launcher: Launcher | None = None, call_dir: Path | None = None
) -> list[str]:
    """The script and its nine options for a prefix (or no) launcher. For an image launcher, the
    ``docker run`` / ``apptainer exec`` argv of M11-14 item 1 binding only ``call_dir``, which
    runs ``Rscript`` on the baked-in script (the shebang is not read)."""
    if launcher is None or launcher.kind == "prefix":
        return [str(script_path()), *_options(call)]
    if call_dir is None:
        raise ValueError(f"the {launcher.kind} fgsea launcher needs the call directory")
    d = str(call_dir)
    tail = ["Rscript", BAKED_SCRIPT, *_options(call)]
    if launcher.kind == "docker":
        user = f"{os.getuid()}:{os.getgid()}"
        return ["docker", "run", "--rm", "-u", user, "-v", f"{d}:{d}", "-w", d,
                "-e", "LC_ALL=C", launcher.image, *tail]  # fmt: skip
    return ["apptainer", "exec", "--cleanenv", "-B", d, "--pwd", d, str(launcher.path), *tail]


def child_env(launcher: Launcher, environ: Mapping[str, str]) -> dict[str, str]:
    """``environ`` with the env's ``bin`` first on ``PATH``, user and site R libraries neutralised
    and ``LC_ALL=C``, so no other R library can shadow the pinned fgsea. An image launcher gets
    ``environ`` unchanged: the container does not see it (``-e LC_ALL=C``, ``--cleanenv``)."""
    if launcher.kind != "prefix" or launcher.path is None:
        return dict(environ)
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
                command(call, self.launcher, cwd),
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
