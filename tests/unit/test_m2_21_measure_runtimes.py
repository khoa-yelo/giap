"""``scripts/measure_runtimes.py`` (docs/specs/M2-21-integration-nightly-runtimes.md example 7).

``measure_runtimes.py TRACE.tsv [--sacct JOBID]`` prints markdown rows
``step | genome | wall | peak RSS | node class (job ID)`` for ``PROGRESS.md`` > Measured numbers.

Inputs are real: the Nextflow trace of the M2-20 sim100k run (Slurm job 1851578 on barley-01) and
that job's real ``sacct`` rows (``tests/data/recorded/runtimes/PROVENANCE.md``). ``sacct`` is a
fake executable first on ``PATH`` that answers whatever ``--format`` field list it is asked for
from the recorded rows (``-P``/``--parsable2`` → ``|``-separated, ``-n``/``--noheader`` honoured)
and logs its argv.

Conventions settled by the test author (the spec names the columns, not these details):

* The trace does not name the genome (the tag is ``genome``), so the genome label is passed with
  ``--genome NAME``.
* One row per process (the last ``:`` component of the trace ``name``, before `` (tag)``),
  lowercased; the step cell starts with that name (``bakta (full)`` allowed). Several tasks of
  one process (16 ``PROMOTER_CALCULATOR`` tasks) give one row whose peak RSS is the largest task's.
* ``wall`` is ``H:MM:SS`` from the task's ``realtime`` (or ``duration``) column; peak RSS keeps a
  unit (``KB``/``MB``/``GB``, numerically equal to the trace value).
* Node class is the node-name prefix (``wheat``/``barley``/``rye``/``oat``) of the ``sacct``
  ``NodeList``, written ``<class> (<job ID>)``.
* ``CACHED`` and ``FAILED`` tasks are not measurements: no row comes from them.
* A markdown header row and ``|---|`` separator may precede the rows; outer ``|`` optional.
"""

from __future__ import annotations

import csv
import os
import re
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "measure_runtimes.py"
DATA = REPO / "tests" / "data" / "recorded" / "runtimes"
TRACE = DATA / "m2_21_trace_sim100k.tsv"
SACCT = DATA / "m2_21_sacct_1851578.txt"
JOB = "1851578"
TOOLS = ("bakta", "genomad", "antismash", "terminatornet", "transtermhp", "ostir", "uniop",
         "promoter_calculator", "promotech")  # fmt: skip

FAKE_SACCT = r"""#!{python}
import os, sys
rows = [ln.rstrip("\n").split("|") for ln in open({data!r}) if ln.strip()]
head, body = rows[0], rows[1:]
nodelist = os.environ.get("M2_21_FAKE_NODELIST")
argv = sys.argv[1:]
with open(os.environ["M2_21_FAKE_SACCT_LOG"], "a") as fh:
    fh.write(" ".join(argv) + "\n")
fields = ["JobID", "JobName", "Partition", "Account", "AllocCPUS", "State", "ExitCode"]
parsable = noheader = trailing = False
i = 0
while i < len(argv):
    a = argv[i]
    if a in ("--format", "-o"):
        fields = argv[i + 1].split(","); i += 1
    elif a.startswith("--format="):
        fields = a.split("=", 1)[1].split(",")
    elif a.startswith("-o") and len(a) > 2:
        fields = a[2:].split(",")
    elif a in ("-P", "--parsable2"):
        parsable = True
    elif a in ("-p", "--parsable"):
        parsable = trailing = True
    elif a in ("-n", "--noheader"):
        noheader = True
    elif a.startswith("-") and not a.startswith("--") and set(a[1:]) <= set("nPpX"):
        parsable = parsable or "P" in a or "p" in a
        trailing = trailing or "p" in a
        noheader = noheader or "n" in a
    i += 1
fields = [f.split("%")[0] for f in fields]
extra = {"Partition": "normal", "AllocCPUS": "16", "ExitCode": "0:0", "Account": "giap",
         "Start": "2026-10-16T14:01:00", "End": "2026-10-16T14:17:59", "ReqMem": "64G",
         "NNodes": "1", "CPUTime": "04:31:44", "TotalCPU": "01:12:00"}
out = []
for r in body:
    rec = dict(zip(head, r))
    if nodelist:
        rec["NodeList"] = nodelist
    out.append([rec.get(f, extra.get(f, "")) for f in fields])
sep = "|" if parsable else "  "
if not noheader:
    print(sep.join(fields))
for r in out:
    print(sep.join(r) + ("|" if trailing else ""))
"""


@pytest.fixture
def fake_sacct(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "sacct"
    exe.write_text(
        FAKE_SACCT.replace("{python}", sys.executable).replace("{data!r}", repr(str(SACCT)))
    )
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    log = tmp_path / "sacct.argv"
    log.write_text("")
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("M2_21_FAKE_SACCT_LOG", str(log))
    return log


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} does not exist (M2-21)"
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        timeout=120,
        env=os.environ.copy(),
        check=False,
    )


def _rows(stdout: str) -> list[list[str]]:
    rows = []
    for ln in stdout.splitlines():
        if "|" not in ln:
            continue
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if not cells or cells[0].lower() == "step" or all(set(c) <= set("-: ") for c in cells):
            continue
        rows.append(cells)
    return rows


def _proc(name: str) -> str:
    return name.split(" (")[0].split(":")[-1].lower()


def _trace() -> list[dict[str, str]]:
    with TRACE.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


_UNIT = {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1.0, "TB": 1e3}


def _gb(text: str) -> float:
    m = re.fullmatch(r"([0-9.]+)\s*([KMGT]?B)", text.strip())
    assert m, f"peak RSS without a unit: {text!r}"
    return float(m.group(1)) * _UNIT[m.group(2)]


def _hms(nf: str) -> str:
    """Nextflow ``1m 8s`` / ``50.4s`` / ``13m 15s`` → ``H:MM:SS`` (seconds truncated)."""
    secs = 0.0
    for num, unit in re.findall(r"([0-9.]+)(ms|h|m|s|d)", nf):
        secs += float(num) * {"ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    s = int(secs)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


@pytest.fixture
def table(fake_sacct: Path) -> dict[str, list[str]]:
    res = _run(str(TRACE), "--sacct", JOB, "--genome", "sim100k")
    assert res.returncode == 0, res.stderr
    rows = _rows(res.stdout)
    by_step: dict[str, list[str]] = {}
    for r in rows:
        assert len(r) == 5, f"expected 5 columns: {r}"
        key = r[0].split(" (")[0].strip().lower()
        assert key not in by_step, f"two rows for {key}: {rows}"
        by_step[key] = r
    return by_step


def test_example_7_one_row_per_tool(table: dict[str, list[str]]) -> None:
    missing = [t for t in TOOLS if t not in table]
    assert not missing, f"no row for {missing}: {sorted(table)}"
    for row in table.values():
        assert row[1] == "sim100k", row
        assert re.fullmatch(r"\d+:\d\d:\d\d", row[2]), row
        _gb(row[3])
        assert row[4] == f"barley ({JOB})", row


def test_example_7_values_come_from_the_trace(table: dict[str, list[str]]) -> None:
    rows = {_proc(r["name"]): r for r in _trace() if _proc(r["name"]) in ("bakta", "genomad")}
    assert table["bakta"][2] == "0:13:15"
    assert table["genomad"][2] in {_hms(rows["genomad"]["realtime"]),
                                   _hms(rows["genomad"]["duration"])}  # fmt: skip
    assert table["bakta"][3] and abs(_gb(table["bakta"][3]) - 11.0) < 0.05
    assert abs(_gb(table["genomad"][3]) - 4.2) < 0.05
    assert abs(_gb(table["terminatornet"][3]) - 0.6014) < 0.001


def test_parallel_tasks_give_one_row_with_the_largest_peak_rss(
    table: dict[str, list[str]],
) -> None:
    tasks = [r for r in _trace() if _proc(r["name"]) == "promoter_calculator"]
    assert len(tasks) == 16
    biggest = max(_gb(r["peak_rss"]) for r in tasks)
    assert abs(_gb(table["promoter_calculator"][3]) - biggest) < 1e-3


def test_sacct_is_asked_about_the_given_job(fake_sacct: Path, table: dict[str, list[str]]) -> None:
    calls = fake_sacct.read_text().splitlines()
    assert calls, "sacct was never called"
    assert any(re.search(rf"(^|\s)(-j|--jobs)(\s+|=){JOB}(\s|$)", c) for c in calls), calls


@pytest.mark.parametrize(
    ("node", "klass"),
    [("wheat-05", "wheat"), ("barley-03", "barley"), ("rye-02", "rye"), ("oat-11", "oat")],
)
def test_node_class_from_node_name(
    fake_sacct: Path, monkeypatch: pytest.MonkeyPatch, node: str, klass: str
) -> None:
    monkeypatch.setenv("M2_21_FAKE_NODELIST", node)
    res = _run(str(TRACE), "--sacct", JOB, "--genome", "sim100k")
    assert res.returncode == 0, res.stderr
    rows = _rows(res.stdout)
    assert rows and all(r[4] == f"{klass} ({JOB})" for r in rows), rows


def test_cached_and_failed_tasks_are_not_measurements(fake_sacct: Path, tmp_path: Path) -> None:
    lines = TRACE.read_text().splitlines(keepends=True)
    head = lines[0].rstrip("\n").split("\t")
    st, rt = head.index("status"), head.index("realtime")
    out = [lines[0]]
    for ln in lines[1:]:
        f = ln.rstrip("\n").split("\t")
        if _proc(f[head.index("name")]) == "genomad":
            f[st], f[rt] = "CACHED", "-"
        if _proc(f[head.index("name")]) == "antismash":
            f[st] = "FAILED"
        out.append("\t".join(f) + "\n")
    trace = tmp_path / "trace.tsv"
    trace.write_text("".join(out))
    res = _run(str(trace), "--sacct", JOB, "--genome", "sim100k")
    assert res.returncode == 0, res.stderr
    steps = {r[0].split(" (")[0].strip().lower() for r in _rows(res.stdout)}
    assert "genomad" not in steps and "antismash" not in steps, steps
    assert "bakta" in steps


def test_missing_trace_is_an_error(fake_sacct: Path, tmp_path: Path) -> None:
    res = _run(str(tmp_path / "nope.tsv"), "--sacct", JOB, "--genome", "sim100k")
    assert res.returncode != 0
    assert "nope.tsv" in res.stderr
