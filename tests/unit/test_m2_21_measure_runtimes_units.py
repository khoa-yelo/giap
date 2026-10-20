"""``scripts/measure_runtimes.py``: units, aggregation, node class per task (M2-21 dec. 16, 17).

Kills the reviewer's surviving hand mutants of PR #107 (longest vs shortest wall, ``realtime`` vs
``duration``, the ``ms``/``d``/``KB`` factors, the ``ms`` alternative before ``m`` in the unit
regex, the minutes field of ``H:MM:SS``, ``ABORTED`` rows, ignored ``sacct`` failures) with small
synthetic traces in the format of ``tests/data/recorded/runtimes/m2_21_trace_sim100k.tsv`` and
with the recorded trace itself; and pins decision 16 (the node class of a row is the node its own
task ran on, from the trace ``hostname`` column, not the head job's node).

Conventions are those of ``test_m2_21_measure_runtimes.py``; in addition (test author):

* Wall times are truncated to whole seconds (``61.5 s`` → ``0:01:01``; ``746ms`` → ``0:00:00``).
* With several tasks the peak-RSS cell is the trace text of the largest task (units compared as
  bytes: 1 KB = 1024 B, …), the wall cell the longest task's ``realtime``.
* With a ``hostname`` trace column the row's node cell starts with that host's class; a process
  whose tasks ran on different classes is not required to be split, only its cell must name a
  class its tasks ran on. ``--sacct`` stays the fallback for traces without ``hostname``.
* ``workflow/nextflow.config`` lists ``hostname`` in ``trace.fields`` so new traces carry it.
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.slurm_sandbox import write_exe

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "measure_runtimes.py"
TRACE = REPO / "tests" / "data" / "recorded" / "runtimes" / "m2_21_trace_sim100k.tsv"
HEAD = [
    "task_id",
    "hash",
    "native_id",
    "name",
    "tag",
    "status",
    "exit",
    "attempt",
    "submit",
    "duration",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
    "workdir",
]

FAKE_SACCT = r"""#!{python}
import os, sys
with open(os.environ["M2_21_SACCT_LOG"], "a") as fh:
    fh.write(" ".join(sys.argv[1:]) + "\n")
rc = int(os.environ.get("M2_21_SACCT_RC", "0"))
if rc:
    print("sacct: error: Unable to contact slurm controller (fake)", file=sys.stderr)
    sys.exit(rc)
node = os.environ.get("M2_21_SACCT_NODE", "barley-01")
job = sys.argv[sys.argv.index("-j") + 1] if "-j" in sys.argv else "1"
if os.environ.get("M2_21_SACCT_UNASSIGNED"):  # a pending step row before the allocation rows
    print(f"{job}|None assigned")
print(f"{job}|{node}")
print(f"{job}.batch|{node}")
"""


@pytest.fixture
def sacct(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    write_exe(bindir / "sacct", FAKE_SACCT.replace("{python}", sys.executable))
    log = tmp_path / "sacct.log"
    log.write_text("")
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("M2_21_SACCT_LOG", str(log))
    return log


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True,
                          timeout=120, env=os.environ.copy(), check=False)  # fmt: skip


def _row(name: str, realtime: str, duration: str | None = None, rss: str = "1 GB",
         status: str = "COMPLETED", **extra: str) -> dict[str, str]:  # fmt: skip
    r = dict.fromkeys(HEAD, "-")
    r.update(
        {
            "task_id": "1",
            "hash": "aa/bbbbbb",
            "native_id": "1",
            "name": name,
            "tag": "genome",
            "status": status,
            "exit": "0",
            "attempt": "1",
            "duration": duration or realtime,
            "realtime": realtime,
            "peak_rss": rss,
            "workdir": "/nonexistent/work/aa/bbbbbb",
        }
    )
    r.update(extra)
    return r


def _write(tmp: Path, rows: list[dict[str, str]], head: list[str] = HEAD) -> Path:
    p = tmp / "trace.tsv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=head, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return p


def _table(stdout: str) -> dict[str, list[str]]:
    out = {}
    for ln in stdout.splitlines():
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) < 5 or cells[0].lower() == "step" or set(cells[0]) <= set("-: "):
            continue
        out[cells[0].split(" (")[0].strip().lower()] = cells
    return out


def _measure(tmp: Path, rows: list[dict[str, str]], *args: str,
             head: list[str] = HEAD) -> dict[str, list[str]]:  # fmt: skip
    res = _run(str(_write(tmp, rows, head)), "--genome", "g", *args)
    assert res.returncode == 0, res.stderr
    return _table(res.stdout)


@pytest.mark.parametrize(
    ("realtime", "wall"),
    [
        ("746ms", "0:00:00"),  # `m` before `ms` in the regex reads 746 minutes
        ("61500ms", "0:01:01"),  # ms factor 1e-3
        ("50.4s", "0:00:50"),
        ("13m 15s", "0:13:15"),
        ("1h 5m 7s", "1:05:07"),  # minutes field is minutes within the hour
        ("2h 59m 59s", "2:59:59"),
        ("1d 2h 3m 4s", "26:03:04"),  # d factor 86400
        ("1d", "24:00:00"),
    ],
)
def test_wall_units(tmp_path: Path, sacct: Path, realtime: str, wall: str) -> None:
    table = _measure(tmp_path, [_row("ANNOTATE:GENOMAD (genome)", realtime)])
    assert table["genomad"][2] == wall, table


def test_wall_is_realtime_not_duration(tmp_path: Path, sacct: Path) -> None:
    rows = [_row("ANNOTATE:GENOMAD (genome)", realtime="1m 0s", duration="2m 30s")]
    assert _measure(tmp_path, rows)["genomad"][2] == "0:01:00"


def test_duration_is_the_fallback_without_realtime(tmp_path: Path, sacct: Path) -> None:
    rows = [_row("ANNOTATE:GENOMAD (genome)", realtime="-", duration="2m 30s")]
    assert _measure(tmp_path, rows)["genomad"][2] == "0:02:30"


def test_several_tasks_report_the_longest_wall(tmp_path: Path, sacct: Path) -> None:
    rows = [_row("ANNOTATE:PROMOTECH (w1)", "5s"), _row("ANNOTATE:PROMOTECH (w2)", "1m 40s"),
            _row("ANNOTATE:PROMOTECH (w3)", "20s")]  # fmt: skip
    assert _measure(tmp_path, rows)["promotech"][2] == "0:01:40"


def test_recorded_promoter_calculator_wall_is_the_longest_of_16(sacct: Path) -> None:
    with TRACE.open(newline="") as fh:
        pc = [r for r in csv.DictReader(fh, delimiter="\t") if "PROMOTER_CALCULATOR" in r["name"]]
    assert len(pc) == 16
    res = _run(str(TRACE), "--sacct", "1851578", "--genome", "sim100k")
    assert res.returncode == 0, res.stderr
    # realtime 8.8 s is the longest (duration 9 s; the shortest realtime is 4.1 s)
    assert _table(res.stdout)["promoter_calculator"][2] == "0:00:08"


@pytest.mark.parametrize(
    ("sizes", "largest"),
    [
        (["900 KB", "1 MB"], "1 MB"),  # KB factor 1024 (a MB-sized KB factor picks 900 KB)
        (["2048 KB", "1.5 MB"], "2048 KB"),
        (["1030 KB", "1 MB"], "1030 KB"),  # 1 MB = 1024 KB, not 1000 KB
        (["1030 MB", "1 GB"], "1030 MB"),
        (["1023 MB", "1 GB"], "1 GB"),
        (["1.5 GB", "1600 MB"], "1600 MB"),
        (["512 B", "1 KB"], "1 KB"),
        (["0.5 TB", "600 GB"], "600 GB"),
    ],
)
def test_peak_rss_units(tmp_path: Path, sacct: Path, sizes: list[str], largest: str) -> None:
    rows = [_row(f"ANNOTATE:UNIOP (t{i})", "1s", rss=s) for i, s in enumerate(sizes)]
    rows += [_row(f"ANNOTATE:OSTIR (t{i})", "1s", rss=s) for i, s in enumerate(reversed(sizes))]
    table = _measure(tmp_path, rows)
    assert table["uniop"][3] == largest, table
    assert table["ostir"][3] == largest, table


@pytest.mark.parametrize("status", ["ABORTED", "FAILED", "CACHED"])
def test_unfinished_tasks_are_not_measured(tmp_path: Path, sacct: Path, status: str) -> None:
    rows = [_row("ANNOTATE:GENOMAD (genome)", "1s"),
            _row("ANNOTATE:BAKTA (genome)", "9h", status=status, rss="60 GB"),
            _row("ANNOTATE:UNIOP (a)", "10s", rss="1 GB"),
            _row("ANNOTATE:UNIOP (b)", "3h", status=status, rss="50 GB")]  # fmt: skip
    table = _measure(tmp_path, rows)
    assert "bakta" not in table, table
    assert table["uniop"][2:4] == ["0:00:10", "1 GB"], table


def test_sacct_failure_is_an_error(sacct: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("M2_21_SACCT_RC", "1")
    res = _run(str(TRACE), "--sacct", "1851578", "--genome", "sim100k")
    assert res.returncode != 0, res.stdout
    assert "sacct" in res.stderr, res.stderr
    assert not _table(res.stdout), "rows printed although the node class is unknown"


# --- decision 16: node class from each task's own node -------------------------------------------


def test_trace_config_records_the_task_hostname() -> None:
    conf = (REPO / "workflow" / "nextflow.config").read_text()
    m = re.search(r"trace\s*\{[^}]*fields\s*=\s*'([^']*)'", conf, re.S)
    assert m, "no trace.fields in workflow/nextflow.config"
    assert "hostname" in m.group(1).split(","), m.group(1)


def _host_rows() -> list[dict[str, str]]:
    return [_row("ANNOTATE:BAKTA (genome)", "13m 15s", hostname="wheat-05"),
            _row("ANNOTATE:GENOMAD (genome)", "1m 8s", hostname="oat-11"),
            _row("ANNOTATE:UNIOP (genome)", "4s", hostname="rye-02.stanford.edu")]  # fmt: skip


@pytest.mark.parametrize("with_sacct", [True, False])
def test_node_class_from_the_task_hostname(tmp_path: Path, sacct: Path, with_sacct: bool) -> None:
    """Under -profile farmshare the tasks run on other nodes than the head (barley-01 here)."""
    args = ("--sacct", "1851578") if with_sacct else ()
    table = _measure(tmp_path, _host_rows(), *args, head=[*HEAD, "hostname"])
    got = {k: table[k][4].split()[0] for k in ("bakta", "genomad", "uniop")}
    assert got == {"bakta": "wheat", "genomad": "oat", "uniop": "rye"}, table


def test_sacct_head_node_is_the_fallback_without_hostname(
    tmp_path: Path, sacct: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("M2_21_SACCT_NODE", "barley-03")
    rows = [_row("ANNOTATE:BAKTA (genome)", "13m 15s")]
    table = _measure(tmp_path, rows, "--sacct", "1851578")
    assert table["bakta"][4] == "barley (1851578)", table


def test_hostname_dash_falls_back_to_the_head_node(
    tmp_path: Path, sacct: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("M2_21_SACCT_NODE", "barley-03")
    rows = [_row("ANNOTATE:BAKTA (genome)", "1s", hostname="wheat-01"),
            _row("ANNOTATE:GENOMAD (genome)", "1s", hostname="-")]  # fmt: skip
    table = _measure(tmp_path, rows, "--sacct", "1851578", head=[*HEAD, "hostname"])
    assert table["bakta"][4].split()[0] == "wheat", table
    assert table["genomad"][4].split()[0] == "barley", table


# --- second review of PR #107: decision 26 -------------------------------------------------------


def test_one_task_has_no_task_count_two_tasks_do(tmp_path: Path, sacct: Path) -> None:
    rows = [_row("ANNOTATE:GENOMAD (genome)", "1s"),
            _row("ANNOTATE:UNIOP (a)", "1s"), _row("ANNOTATE:UNIOP (b)", "2s")]  # fmt: skip
    res = _run(str(_write(tmp_path, rows)), "--genome", "g", "--sacct", "1851578")
    assert res.returncode == 0, res.stderr
    steps = [ln.strip().strip("|").split("|")[0].strip() for ln in res.stdout.splitlines()]
    assert "genomad" in steps and "uniop (2 tasks)" in steps, steps


def test_none_assigned_rows_are_skipped(
    tmp_path: Path, sacct: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("M2_21_SACCT_UNASSIGNED", "1")
    monkeypatch.setenv("M2_21_SACCT_NODE", "rye-02")
    table = _measure(tmp_path, [_row("ANNOTATE:BAKTA (genome)", "1s")], "--sacct", "1851578")
    assert table["bakta"][4] == "rye (1851578)", table


def test_no_sacct_call_when_every_task_has_a_hostname(tmp_path: Path, sacct: Path) -> None:
    _measure(tmp_path, _host_rows(), "--sacct", "1851578", head=[*HEAD, "hostname"])
    assert not sacct.read_text().strip(), f"sacct called needlessly: {sacct.read_text()}"
