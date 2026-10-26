"""M11-14 14s2: the weekly K-12 arm on the distance-model fallback (acceptance criterion 10).

Spec (docs/specs/M11-14-container-install.md, 14s2): "Add one arm to `weekly.sbatch`: K-12
standard with `annotation.uniop.enabled=false`. Record the TU agreement with RegulonDB next to the
UniOP arm in `weekly.tsv`. The arm is reported, not gated, and adds no bound." AC 10: "the weekly
arm writes its `weekly.tsv` row". ADR-0022 addendum (14s2): "a `weekly.tsv` row, reported, not
gated".

``scripts/slurm/weekly.sbatch`` runs in the M2-21 sandbox (``helpers.slurm_sandbox``: fake Slurm,
``gh``, ``python -m pytest``, a sandbox ``GIAP_REPO`` with a local ``origin``, a temporary
``GIAP_ROOT``; no network, no real Slurm, no ``$GIAP_DB_DIR``). This file replaces the logging
``giap`` fake with one that also writes ``<outdir>/run_manifest.json`` (``operon_source`` as the
real run records it, M11-14d) and records the ``--genome`` header, and puts two one-record
genomes under the sandbox ``$GIAP_ROOT/data/genomes`` (K-12 ``U00096.3``, Nissle ``CP171242``) so
the K-12 run can be told apart.

Conventions settled by the test author (the spec leaves the shape open):

* **Placement-agnostic.** The arm may live inside the K-12 task or be its own array task: every
  task of ``#SBATCH --array`` runs with one job ID and all ``weekly.tsv`` rows / ``giap`` calls are
  pooled.
* **The fallback row** is the one row whose label (``summary`` up to the first ``:``, the rule of
  ``scripts/release/readiness.py``) is neither ``k12_mg1655`` nor ``nissle1917``; the label names
  K-12 (``k12``). A distinct label keeps readiness R5 (the newest ``k12_mg1655`` row must carry
  TE8) reading the UniOP K-12 row. The row has the five ``parse_log`` columns.
* **Operon source** is reported as ``operon_source=<value>`` with the value read from the arm's
  ``run_manifest.json`` (never a literal): ``distance_model`` when the run recorded it, the other
  value when it did not, and never ``distance_model`` without a manifest.
* **TU agreement** with RegulonDB: the row's text names RegulonDB and, besides the operon-source
  token, the UniOP arm (the value "next to the UniOP arm"); the agreement reads the RegulonDB
  truth under ``$GIAP_DATA_DIR`` (static check). The metric itself is not pinned here.
* **Not gated** (reported only): a failing fallback run changes neither the exit column of the
  other rows nor the exit code of a task that writes one of them, the fallback row's exit column
  stays ``0`` (a red exit column would gate), and its text reports the run's exit.
"""

from __future__ import annotations

import gzip
import importlib.util
import re
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import pytest
from helpers.slurm_sandbox import REPO, Sandbox, assert_honours_giap_repo, make_sandbox, write_exe

WEEKLY = REPO / "scripts" / "slurm" / "weekly.sbatch"
READINESS = REPO / "scripts" / "release" / "readiness.py"
JOB = "4400"
FLAG = "annotation.uniop.enabled=false"
MAIN_LABELS = {"k12_mg1655", "nissle1917"}
K12_ID, NISSLE_ID = "U00096.3", "CP171242.1"

# A ``giap`` fake: logs like the sandbox logger (``calls.log``), records argv and the genome header
# (``giap_runs.jsonl``) and writes ``<outdir>/run_manifest.json`` for ``giap run``. With
# ``annotation.uniop.enabled=false`` the manifest's operon_source is $FAKE_FALLBACK_SOURCE
# (default distance_model), no manifest is written under FAKE_FALLBACK_NO_MANIFEST=1, and the exit
# is $FAKE_FALLBACK_RC (default 0); every other run records uniop and exits 0.
_GIAP = r"""#!{python}
import json, os, sys
state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "calls.log"), "a") as fh:
    fh.write("giap " + " ".join(argv) + "\n")
def opt(name):
    val = None
    for i, a in enumerate(argv):
        if a == name and i + 1 < len(argv):
            val = argv[i + 1]
        elif a.startswith(name + "="):
            val = a.split("=", 1)[1]
    return val
sets = [argv[i + 1] for i, a in enumerate(argv) if a == "--set" and i + 1 < len(argv)]
sets += [a.split("=", 1)[1] for a in argv if a.startswith("--set=")]
fallback = any(s.replace(" ", "").lower() == "annotation.uniop.enabled=false" for s in sets)
genome, outdir, header = opt("--genome"), opt("--outdir"), None
if genome and os.path.isfile(genome):
    with open(genome) as fh:
        header = fh.readline().strip()
rec = {"argv": argv, "sets": sets, "fallback": fallback, "genome_header": header,
       "outdir": outdir, "cwd": os.getcwd()}
with open(os.path.join(state, "giap_runs.jsonl"), "a") as fh:
    fh.write(json.dumps(rec) + "\n")
if argv[:1] == ["run"] and outdir:
    os.makedirs(outdir, exist_ok=True)
    source = os.environ.get("FAKE_FALLBACK_SOURCE", "distance_model") if fallback else "uniop"
    if not (fallback and os.environ.get("FAKE_FALLBACK_NO_MANIFEST")):
        reason = "annotation.uniop.enabled=false" if source == "distance_model" else None
        with open(os.path.join(outdir, "run_manifest.json"), "w") as fh:
            json.dump({"operon_source": source, "operon_source_reason": reason}, fh)
sys.exit(int(os.environ.get("FAKE_FALLBACK_RC", "0")) if fallback else 0)
"""


@dataclass
class Week:
    sb: Sandbox
    rows: list[list[str]]
    rcs: dict[int, int]
    task_rows: dict[int, list[list[str]]]
    runs: list[dict[str, object]]
    out: str


def _array_tasks() -> list[int]:
    m = re.search(r"^#SBATCH\s+--array[= ](\S+)", WEEKLY.read_text(), re.M)
    assert m, "weekly.sbatch is a job array"
    ids: set[int] = set()
    for part in m.group(1).split("%")[0].split(","):
        lo, _, hi = part.partition("-")
        ids |= set(range(int(lo), int(hi or lo) + 1))
    return sorted(ids)


def _genome(root: Path, sample: str, rid: str) -> None:
    d = root / "data" / "genomes" / sample
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / f"{rid}.fna.gz", "wt") as fh:
        fh.write(f">{rid} sandbox stand-in\nACGTACGTACGT\n")


def _week(tmp: Path, **env: str) -> Week:
    assert_honours_giap_repo(WEEKLY.read_text(), "weekly.sbatch")
    sb = make_sandbox(tmp)
    write_exe(sb.bin / "giap", _GIAP.replace("{python}", sys.executable))
    _genome(sb.root, "k12_mg1655", K12_ID)
    _genome(sb.root, "nissle1917", NISSLE_ID)
    tasks = _array_tasks()
    rcs: dict[int, int] = {}
    task_rows: dict[int, list[list[str]]] = {}
    out = []
    seen = 0
    for t in tasks:
        base = {
            "SLURM_ARRAY_TASK_ID": str(t),
            "SLURM_ARRAY_JOB_ID": JOB,
            "SLURM_JOB_ID": str(6100 + t),
            "SLURM_ARRAY_TASK_COUNT": str(len(tasks)),
            "SLURM_CPUS_PER_TASK": "2",
            "SLURM_JOB_NAME": "giap-weekly",
        }
        res = sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "weekly.sbatch")],
                     env={**base, **env})  # fmt: skip
        rcs[t] = res.returncode
        rows = sb.tsv_rows("weekly.tsv")
        task_rows[t] = rows[seen:]
        seen = len(rows)
        out.append(f"--- task {t} rc={res.returncode}\n{res.stdout}\n{res.stderr}")
    return Week(sb, sb.tsv_rows("weekly.tsv"), rcs, task_rows, sb.jsonl("giap_runs.jsonl"),
                "\n".join(out))  # fmt: skip


def _label(row: list[str]) -> str:
    return row[4].split(":", 1)[0].strip() if len(row) >= 5 else ""


def _fallback_rows(w: Week) -> list[list[str]]:
    return [r for r in w.rows if _label(r) not in MAIN_LABELS]


def _fallback_row(w: Week) -> list[str]:
    rows = _fallback_rows(w)
    assert len(rows) == 1, (
        f"expected exactly one weekly.tsv row for the distance-model arm (label other than "
        f"{sorted(MAIN_LABELS)}), got {rows}\nall rows: {w.rows}\n{w.out[-4000:]}"
    )
    return rows[0]


def _code() -> str:
    keep = [ln for ln in WEEKLY.read_text().splitlines()
            if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")]  # fmt: skip
    return "\n".join(keep)


@pytest.fixture(scope="module")
def week(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Week]:
    yield _week(tmp_path_factory.mktemp("week"))


# ------------------------------------------------------------------ static


def test_script_runs_k12_with_uniop_disabled() -> None:
    code = _code()
    assert re.search(r"\bgiap\s+run\b", code)
    assert re.search(r"--set\s+\"?annotation\.uniop\.enabled=false\b", code), (
        "weekly.sbatch has no `giap run ... --set annotation.uniop.enabled=false` arm (M11-14 14s2)"
    )


def test_agreement_reads_the_regulondb_truth() -> None:
    """The arm's TU agreement is against RegulonDB, read from $GIAP_DATA_DIR (real data)."""
    code = _code()
    assert re.search(r"regulondb", code, re.I), "no RegulonDB TU agreement in weekly.sbatch"
    assert re.search(r"\$\{?GIAP_DATA_DIR\}?/regulondb", code), (
        "the RegulonDB truth is read from $GIAP_DATA_DIR/regulondb"
    )


# ------------------------------------------------------------------ the run


def test_exactly_one_fallback_run_on_k12(week: Week) -> None:
    fb = [r for r in week.runs if r["argv"][:1] == ["run"] and r["fallback"]]  # type: ignore[index]
    assert len(fb) == 1, f"fallback giap runs: {fb}\n{week.out[-4000:]}"
    header = str(fb[0]["genome_header"] or "")
    assert header.startswith(f">{K12_ID}"), f"the fallback arm runs on K-12, not {header!r}"


def test_uniop_arm_unchanged_and_separate(week: Week) -> None:
    runs = [r for r in week.runs if r["argv"][:1] == ["run"]]  # type: ignore[index]
    k12 = [r for r in runs if str(r["genome_header"] or "").startswith(f">{K12_ID}")]
    uniop = [r for r in k12 if not r["fallback"]]
    assert uniop, f"the K-12 UniOP arm(s) no longer run: {k12}"
    for r in uniop:
        assert not any("uniop.enabled" in s for s in r["sets"]), r  # type: ignore[attr-defined]
    fb = [r for r in k12 if r["fallback"]]
    assert fb and fb[0]["outdir"] not in {r["outdir"] for r in uniop}, (
        "the fallback arm writes into its own --outdir"
    )


def test_fallback_row_written_once_with_a_k12_label(week: Week) -> None:
    row = _fallback_row(week)
    assert len(row) >= 5, row
    assert re.fullmatch(r"-?\d+", row[3].strip()) and len(row[2].strip()) >= 7, (
        f"not a readiness parse_log row (date, job, commit, exit, summary): {row}"
    )
    assert row[1].startswith(JOB), row
    label = _label(row)
    assert "k12" in label.lower(), f"the fallback row's label names K-12: {label!r}"


def test_fallback_row_reports_operon_source_from_manifest(week: Week) -> None:
    row = _fallback_row(week)
    assert re.search(r"\boperon_source=distance_model\b", row[4]), row[4]


@pytest.mark.parametrize("source", ["uniop"])
def test_operon_source_is_read_not_hard_coded(tmp_path: Path, source: str) -> None:
    """A manifest that recorded another source is reported as such."""
    w = _week(tmp_path, FAKE_FALLBACK_SOURCE=source)
    row = _fallback_row(w)
    assert re.search(rf"\boperon_source={source}\b", row[4]), row[4]
    assert "operon_source=distance_model" not in row[4], row[4]


def test_missing_manifest_never_reports_distance_model(tmp_path: Path) -> None:
    w = _week(tmp_path, FAKE_FALLBACK_NO_MANIFEST="1")
    row = _fallback_row(w)
    assert "operon_source=distance_model" not in row[4], row[4]


def test_fallback_row_reports_regulondb_agreement_next_to_uniop(week: Week) -> None:
    row = _fallback_row(week)
    text = row[4].split(":", 1)[1] if ":" in row[4] else row[4]
    assert re.search(r"regulondb", text, re.I), f"no RegulonDB TU agreement: {row[4]}"
    rest = re.sub(r"\boperon_source=\S+", "", text)
    assert re.search(r"\buniop\b", rest, re.I), (
        f"the UniOP arm's agreement is not reported next to the fallback's: {row[4]}"
    )


def test_readiness_still_reads_the_uniop_k12_row(week: Week) -> None:
    """R5 takes the newest `k12_mg1655` row for TE8; the fallback row must not shadow it."""
    spec = importlib.util.spec_from_file_location("giap_readiness_14s2", READINESS)
    assert spec is not None and spec.loader is not None
    mod: ModuleType = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    text = (week.sb.root / "logs" / "weekly.tsv").read_text()
    rows, bad = mod.parse_log(text)
    assert bad == 0, text
    _fallback_row(week)
    hit = mod.newest(rows, lambda c: True, lambda r: r.label == "k12_mg1655")
    assert hit is not None and "TE8" in hit.summary, (hit, text)


# ------------------------------------------------------------------ reported, not gated


def test_failing_fallback_run_is_reported_not_gated(tmp_path: Path) -> None:
    ok = _week(tmp_path / "ok", FAKE_FALLBACK_RC="0")
    bad = _week(tmp_path / "bad", FAKE_FALLBACK_RC="1")

    def main(w: Week) -> tuple[list[tuple[str, str]], dict[int, int]]:
        rows = sorted((_label(r), r[3]) for r in w.rows if _label(r) in MAIN_LABELS)
        rcs = {t: rc for t, rc in w.rcs.items()
               if any(_label(r) in MAIN_LABELS for r in w.task_rows[t])}  # fmt: skip
        return rows, rcs

    assert main(ok)[0], f"no K-12/Nissle rows at all:\n{ok.out[-4000:]}"
    assert main(bad) == main(ok), (
        "a failing distance-model run changed the gated rows or task exits "
        f"(reported, not gated):\nok  {main(ok)}\nbad {main(bad)}"
    )
    for w in (ok, bad):
        row = _fallback_row(w)
        assert row[3] == "0", f"the fallback row's exit column gates: {row}"
    assert re.search(r"\bexit 1\b", _fallback_row(bad)[4]), (
        f"the fallback run's failure is not reported: {_fallback_row(bad)}"
    )
