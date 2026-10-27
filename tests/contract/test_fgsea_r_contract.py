"""T4 and T5: ``gippa_fgsea.R`` on the real ``gippa-fgsea`` env (M8-17a; ADR-0038 item 7; A1, A2).

Slow and ``requires_r``: runs only on Slurm (``scripts/slurm/test_branch.sbatch -m requires_r``,
at most 2 CPUs) and is skipped, with the reason, where the env or #192's golden data are absent.

- **A1 golden contract.** #192's frozen ``inputs/gsea_stats.tsv.gz`` / ``gsea_sets.tsv.gz`` as one
  job ``main:golden`` ([10, 500], GOLDEN_SEED, 101 / 0 / 1000) give ``results.tsv`` columns
  ``pathway, pval, log2err, ES, NES, size, leadingEdge`` equal to ``outputs/fgsea.tsv.gz`` as
  strings, row for row; through the wrapper the parsed numbers equal them exactly.
  The golden sets name 9 TUs absent from the ranking (``TU_ABSENT_*``); the spec's input rule
  (every member is a ranked TU, exit 3 otherwise) means the call carries each set restricted to
  ranked TUs, in file order, which is what fgsea itself matches against (the golden ``size`` is
  the matched size).
- ``set.seed`` before every job: the golden job's rows do not change when another job precedes it.
- **A2** two identical calls give byte-identical ``results.tsv``.
- **T5 tie order**: with 40 % tied scores, the leading edge follows the input order of ties.
- The script's own exit codes (2 usage, 3 invalid input) and its total behaviour (no job → header
  only), with no partial ``results.tsv`` on failure.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml
from helpers.fgsea_r_data import GOLDEN, RESULTS_HEADER, SCRIPT, real_prefix

REPO = Path(__file__).resolve().parents[2]
TOLERANCES = REPO / "tests" / "calibration" / "tolerances.yaml"
PREFIX = real_prefix()

pytestmark = [
    pytest.mark.slow,
    pytest.mark.requires_r,
    pytest.mark.skipif(
        PREFIX is None,
        reason="the gippa-fgsea env is not built (sbatch scripts/slurm/build_envs.sbatch fgsea); "
        "run on Slurm with scripts/slurm/test_branch.sbatch -m requires_r",
    ),
]
needs_golden = pytest.mark.skipif(
    not (GOLDEN / "outputs" / "fgsea.tsv.gz").is_file(),
    reason="#192's golden R data (tests/data/golden/r/) are not on this branch yet",
)
GOLDEN_COLUMNS = ("pathway", "pval", "log2err", "ES", "NES", "size", "leadingEdge")


def _golden_seed() -> int:
    return int(yaml.safe_load(TOLERANCES.read_text())["seeds"]["GOLDEN_SEED"])


def _gz_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def _golden_job(job_id: str = "main:golden"):
    from gippa.function.gsea import GseaJob

    ranked = {r["tu_id"] for r in _gz_rows(GOLDEN / "inputs" / "gsea_stats.tsv.gz")}
    sets: dict[str, list[str]] = {}
    for r in _gz_rows(GOLDEN / "inputs" / "gsea_sets.tsv.gz"):
        members = sets.setdefault(r["set_id"], [])
        if r["member"] in ranked:
            members.append(r["member"])
    return GseaJob(job_id, {k: tuple(v) for k, v in sets.items() if v}, 10, 500)


def _golden_call(jobs: tuple | None = None):
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall

    stats = _gz_rows(GOLDEN / "inputs" / "gsea_stats.tsv.gz")
    return EngineCall(
        tu_ids=tuple(r["tu_id"] for r in stats),
        scores=np.array([float(r["score"]) for r in stats], dtype=np.float64),
        jobs=jobs or (_golden_job(),),
        params=GseaConfig(),
        seed=_golden_seed(),
    )


def _child_env() -> dict[str, str]:
    assert PREFIX is not None
    env = dict(os.environ)
    env.update(
        PATH=f"{PREFIX / 'bin'}{os.pathsep}{os.environ['PATH']}",
        R_LIBS_USER="",
        R_LIBS="",
        R_PROFILE_USER="/dev/null",
        LC_ALL="C",
    )
    return env


def _run(d: Path, call, *extra: str, args: list[str] | None = None) -> subprocess.CompletedProcess:
    """Run the script directly on ``write_call(call)`` in ``d`` (string-level comparisons)."""
    from gippa.io.fgsea_r import write_call

    if call is not None:
        write_call(call, d)
    argv = (
        args
        if args is not None
        else [
            "--stats",
            "stats.tsv",
            "--sets",
            "sets.tsv",
            "--jobs",
            "jobs.tsv",
            "--sample-size",
            str(call.params.sample_size),
            "--eps",
            "0",
            "--n-perm-simple",
            str(call.params.n_perm_simple),
            "--seed",
            str(call.seed),
            "--out",
            "results.tsv",
            "--versions",
            "versions.txt",
        ]
    )
    return subprocess.run(
        [str(SCRIPT), *argv, *extra],
        cwd=d,
        env=_child_env(),
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )


def _results(d: Path) -> list[dict[str, str]]:
    text = (d / "results.tsv").read_text(encoding="utf-8")
    assert text.split("\n", 1)[0] == RESULTS_HEADER
    return list(csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE))


# ---------------------------------------------------------------- A1: golden contract


@needs_golden
def test_golden_contract_strings(tmp_path: Path) -> None:
    proc = _run(tmp_path, _golden_call())
    assert proc.returncode == 0, proc.stderr
    got = [{c: r[c] for c in GOLDEN_COLUMNS} for r in _results(tmp_path)]
    want = [
        {c: r[c] for c in GOLDEN_COLUMNS} for r in _gz_rows(GOLDEN / "outputs" / "fgsea.tsv.gz")
    ]
    assert len(got) == len(want) == 119
    for i, (g, w) in enumerate(zip(got, want, strict=True)):
        assert g == w, f"row {i}"
    assert {r["job_id"] for r in _results(tmp_path)} == {"main:golden"}
    versions = (tmp_path / "versions.txt").read_text(encoding="utf-8").splitlines()
    assert "fgsea\t1.36.2" in versions
    assert f"seed\t{_golden_seed()}" in versions
    assert "jobs\t1" in versions
    assert any(v.startswith("R\tR version 4.5.3") for v in versions)
    assert not (tmp_path / "results.tsv.tmp").exists()


@needs_golden
def test_golden_contract_through_wrapper(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import r_fgsea_engine, resolve_launcher

    engine = r_fgsea_engine(resolve_launcher({"GIPPA_FGSEA_PREFIX": str(PREFIX)}), tmp_root=tmp_path)
    rows = engine(_golden_call())["main:golden"]
    want = _gz_rows(GOLDEN / "outputs" / "fgsea.tsv.gz")
    assert [r.term for r in rows] == [w["pathway"] for w in want]
    for r, w in zip(rows, want, strict=True):
        assert r.size == int(w["size"])
        assert (r.es, r.nes, r.p, r.log2err) == (
            float(w["ES"]),
            float(w["NES"]),
            float(w["pval"]),
            float(w["log2err"]),
        )
        assert r.leading_edge == tuple(w["leadingEdge"].split(";"))
    assert list(tmp_path.iterdir()) == []


@needs_golden
def test_seed_is_set_before_every_job(tmp_path: Path) -> None:
    from gippa.function.gsea import GseaJob

    golden = _golden_job()
    first = GseaJob("main:other", dict(list(golden.sets.items())[:30]), 10, 500)
    proc = _run(tmp_path, _golden_call((first, golden)))
    assert proc.returncode == 0, proc.stderr
    got = [r for r in _results(tmp_path) if r["job_id"] == "main:golden"]
    want = _gz_rows(GOLDEN / "outputs" / "fgsea.tsv.gz")
    assert [{c: r[c] for c in GOLDEN_COLUMNS} for r in got] == [
        {c: r[c] for c in GOLDEN_COLUMNS} for r in want
    ]


# ---------------------------------------------------------------- A2 and T5: determinism, ties


def _ties_call(order: str, seed: int = 20261022):
    """50 TUs; 20 (40 %) tie at the top, 10 of them in set S, interleaved with 10 others."""
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall, GseaJob

    tied = [f"T{i:02d}" for i in range(20)]  # even = members of S
    if order == "reversed":
        tied = tied[::-1]
    rest = [f"U{i:02d}" for i in range(30)]
    ids = tuple(tied + rest)
    scores = [2.0] * 20 + [1.5 - 0.1 * i for i in range(30)]
    members = tuple(f"T{i:02d}" for i in range(0, 20, 2))
    other = tuple(f"U{i:02d}" for i in range(0, 30, 3))
    jobs = (GseaJob("main:ties", {"S": members, "O": other}, 5, 49),)
    return EngineCall(ids, np.array(scores), jobs, GseaConfig(), seed)


def test_ties_follow_input_order(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import r_fgsea_engine, resolve_launcher

    engine = r_fgsea_engine(resolve_launcher({"GIPPA_FGSEA_PREFIX": str(PREFIX)}), tmp_root=tmp_path)
    edges = {}
    for order in ("forward", "reversed"):
        call = _ties_call(order)
        row = next(r for r in engine(call)["main:ties"] if r.term == "S")
        assert row.es > 0 and row.leading_edge
        in_order = [t for t in call.tu_ids if t in set(call.jobs[0].sets["S"])]
        assert list(row.leading_edge) == in_order[: len(row.leading_edge)], order
        edges[order] = row.leading_edge
    assert edges["forward"] != edges["reversed"]


def test_identical_calls_are_byte_identical(tmp_path: Path) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    for d in (a, b):
        proc = _run(d, _ties_call("forward"))
        assert proc.returncode == 0, proc.stderr
    assert (a / "results.tsv").read_bytes() == (b / "results.tsv").read_bytes()
    assert (a / "versions.txt").read_bytes() == (b / "versions.txt").read_bytes()


# ---------------------------------------------------------------- the script's exit codes


def _write(d: Path, stats: str, jobs: str, sets: str) -> None:
    (d / "stats.tsv").write_text(stats)
    (d / "jobs.tsv").write_text(jobs)
    (d / "sets.tsv").write_text(sets)


GOOD_STATS = "tu_id\tscore\n" + "".join(f"T{i}\t{5 - i}.5\n" for i in range(10))
GOOD_JOBS = "job_id\tmin_size\tmax_size\nmain:x\t1\t9\n"
GOOD_SETS = "job_id\tset_id\tmember\n" + "".join(f"main:x\tS\tT{i}\n" for i in (0, 2, 4))
BAD_INPUTS = {
    "NA score": (GOOD_STATS.replace("T3\t2.5", "T3\tNA"), GOOD_JOBS, GOOD_SETS),
    "Inf score": (GOOD_STATS.replace("T3\t2.5", "T3\tInf"), GOOD_JOBS, GOOD_SETS),
    "duplicate tu_id": (GOOD_STATS.replace("T3\t", "T2\t"), GOOD_JOBS, GOOD_SETS),
    "member not ranked": (GOOD_STATS, GOOD_JOBS, GOOD_SETS + "main:x\tS\tT99\n"),
    "job not in jobs.tsv": (GOOD_STATS, GOOD_JOBS, GOOD_SETS + "main:y\tS\tT1\n"),
    "min_size < 1": (GOOD_STATS, GOOD_JOBS.replace("\t1\t9", "\t0\t9"), GOOD_SETS),
    "max_size < min_size": (GOOD_STATS, GOOD_JOBS.replace("\t1\t9", "\t5\t4"), GOOD_SETS),
}


def _args(**over: str) -> list[str]:
    base = {
        "--stats": "stats.tsv",
        "--sets": "sets.tsv",
        "--jobs": "jobs.tsv",
        "--sample-size": "101",
        "--eps": "0",
        "--n-perm-simple": "1000",
        "--seed": "7",
        "--out": "results.tsv",
        "--versions": "versions.txt",
    }
    base.update(over)
    return [x for kv in base.items() for x in kv if x != ""]


@pytest.mark.parametrize("name", sorted(BAD_INPUTS))
def test_invalid_input_exits_3(tmp_path: Path, name: str) -> None:
    _write(tmp_path, *BAD_INPUTS[name])
    proc = _run(tmp_path, None, args=_args())
    assert proc.returncode == 3, (name, proc.stderr)
    assert proc.stderr.strip()
    assert not (tmp_path / "results.tsv").exists()


@pytest.mark.parametrize(
    "args",
    [
        _args(**{"--seed": "1.5"}),
        _args(**{"--seed": "abc"}),
        [a for a in _args() if a not in ("--out", "results.tsv")],
        [*_args(), "--bogus", "1"],
    ],
    ids=["non-integer seed", "text seed", "missing --out", "unknown argument"],
)
def test_usage_error_exits_2(tmp_path: Path, args: list[str]) -> None:
    _write(tmp_path, GOOD_STATS, GOOD_JOBS, GOOD_SETS)
    proc = _run(tmp_path, None, args=args)
    assert proc.returncode == 2, proc.stderr
    assert not (tmp_path / "results.tsv").exists()


def test_no_job_writes_header_only(tmp_path: Path) -> None:
    _write(tmp_path, GOOD_STATS, "job_id\tmin_size\tmax_size\n", "job_id\tset_id\tmember\n")
    proc = _run(tmp_path, None, args=_args())
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "results.tsv").read_text() == RESULTS_HEADER + "\n"
    assert "jobs\t0" in (tmp_path / "versions.txt").read_text().splitlines()


def test_out_of_range_job_adds_no_rows(tmp_path: Path) -> None:
    _write(tmp_path, GOOD_STATS, GOOD_JOBS.replace("\t1\t9", "\t5\t9"), GOOD_SETS)
    proc = _run(tmp_path, None, args=_args())
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "results.tsv").read_text() == RESULTS_HEADER + "\n"
