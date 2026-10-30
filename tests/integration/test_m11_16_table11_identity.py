"""AC 1: table 11 is unchanged by M11-16 (docs/specs/M11-16-translation-tables-4-25.md, AC 1;
test plan row ``tests/integration/test_m11_16_table11_identity.py``; ``slow``, Slurm via
``scripts/slurm/test_branch.sbatch``).

One set of gippasim pilot inputs (``tests/integration/test_m10_05_sim_pilot.py``: SIM100K,
``pi_wrong = 0``, the test ``$GIPPA_DB_DIR``, ``SEED``), made once by this checkout's gippasim,
and two ``gippa run --executor local`` runs on them in subprocesses: one with ``src`` +
``workflow`` of ``HEAD_ARM`` (main after the last M11-16 PR, 16f) and one with those of
``REFERENCE`` (main before any M11-16 code: the 16c merge, whose diff is ``tests/sim`` only), both
extracted with ``git archive``. Each subprocess checks that ``gippa`` is imported from its own
tree. AC 1 is a property of the M11-16 change, so both arms are pinned: later changes to the
outputs (M11-20, #490: initiator methionine, sixframe meta format 4) have their own tests.

- Every file under ``tables/`` and ``mapping/`` is byte-identical between the two runs (same set
  of files, same bytes).
- The documented changes, and only as documented: ``sixframe.npz`` meta (format 2 -> 3,
  ``translation_table`` 11), ``run_manifest.json`` ``annotation.translation_table`` (11) and
  ``output_schema_version`` (one minor version up), ``preflight.json`` (``genome.genetic_code``
  ``value.translation_table`` 11) and the report's Methods sentence (no genetic-code note at 11).
  The Bakta argv is unchanged at 11 (D2): these runs use the precomputed Bakta JSON.
- Every other published file outside ``work/`` is byte-identical, except ``DIFFER``: the AC 1
  files above plus ``BEYOND_AC1``, which AC 1 does not list but which follow from other parts of
  the spec (16f finding, for the lead): the new ``config.tool_genetic_code`` preflight row (2.3)
  and the new ``config.translation_table`` rule text, so ``preflight.md`` and the V17 QC plot
  count one more ``ok``; ``schemas/datapackage.json``'s ``output_schema_version``; and the
  manifest's config snapshot, whose ``mapping.proteins.tblastn_args`` default lost
  ``-db_gencode 11`` (D1, ADR-0067), and so its ``provenance_digest``.

Skips (never passes) when ``git`` or either pinned commit is unavailable (a shallow clone).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m4_real as mr

pytestmark = pytest.mark.slow

REPO = Path(__file__).resolve().parents[2]
#: main at the fork of the M11-16 code PRs: 0980bf7c (#474, M11-16c, tests/sim only)
REFERENCE = "0980bf7caabd0ed7cd6606606ccb8a22de146b35"
#: main after the last M11-16 PR: 7326f545 (#481, M11-16f), the parent of M11-20 (#490)
HEAD_ARM = "7326f5452c748949b116bc1fec2a4c96d1560e81"
SEED = 20261023
IDENTICAL = ("tables", "mapping")
METHODS = "Genes were called and all six frames translated with NCBI genetic code 11"
#: published files AC 1 names as changed (``logs/events.jsonl`` carries timestamps)
AC1_FILES = frozenset({"run_manifest.json", "preflight/preflight.json", "report/report.html",
                       "logs/events.jsonl"})  # fmt: skip
#: published files that change at 11 beyond AC 1's list (module docstring)
BEYOND_AC1 = frozenset({"preflight/preflight.md", "schemas/datapackage.json",
                        "plots/V17_qc.pdf", "plots/V17_qc.png", "plots/V17_qc.svg",
                        "plots/data/V17_qc.tsv"})  # fmt: skip
DIFFER = AC1_FILES | BEYOND_AC1
TBLASTN = ("-seg no -evalue 1e-5", "-db_gencode 11 -seg no -evalue 1e-5")  # (head, reference)

_RUN = """
import sys
from pathlib import Path
import gippa
want = Path(sys.argv[1]).resolve()
assert Path(gippa.__file__).resolve().is_relative_to(want), (gippa.__file__, want)
from gippa.cli import cli
cli(sys.argv[2:], standalone_mode=True)
"""


def _git(*args: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=False)


def _tree(commit: str, dest: Path) -> Path:
    if shutil.which("git") is None or _git("cat-file", "-e", f"{commit}^{{commit}}").returncode:
        pytest.skip(f"commit {commit[:8]} not available (shallow clone?)")
    archive = _git("archive", "--format=tar", commit, "src", "workflow")
    assert archive.returncode == 0, archive.stderr.decode()
    dest.mkdir(parents=True)
    subprocess.run(["tar", "-x", "-C", str(dest)], input=archive.stdout, check=True)
    return dest


def _run(src: Path, sim: Any, out: Path, db: Path) -> Path:
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml),
            "--outdir", str(out), "--seed", str(SEED), "--executor", "local"]  # fmt: skip
    env = {**os.environ, "GIPPA_DB_DIR": str(db), "PYTHONPATH": str(src)}
    r = subprocess.run([sys.executable, "-c", _RUN, str(src / "gippa"), *args], env=env,
                       cwd=out.parent, capture_output=True, text=True, check=False)  # fmt: skip
    assert r.returncode == 0, (src, r.stdout[-4000:], r.stderr[-4000:])
    return out


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """The two runs, made once per pytest session: under xdist the workers share one directory
    under the session's base temp, and the first worker builds it under an ``flock`` while the
    others wait and reuse it (a ``done`` marker; a failed build leaves none, so the next worker
    retries and fails the same way)."""
    import fcntl

    base = tmp_path_factory.getbasetemp()
    shared = base.parent if os.environ.get("PYTEST_XDIST_WORKER") else base
    root = shared / "m11_16_identity"
    out = {"head": root / "out_head", "ref": root / "out_ref"}
    with (shared / "m11_16_identity.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not (root / "done").exists():
            shutil.rmtree(root, ignore_errors=True)
            root.mkdir()
            ref = _tree(REFERENCE, root / "reference")
            head = _tree(HEAD_ARM, root / "head_arm")
            db = mc.test_db(root / "db")
            sim = mr.sim100k(root / "sim", pi_wrong=0.0)
            _run(head / "src", sim, out["head"], db)
            _run(ref / "src", sim, out["ref"], db)
            (root / "done").touch()
    return out


def _files(out: Path) -> dict[str, Path]:
    return {
        p.relative_to(out).as_posix(): p
        for d in IDENTICAL
        for p in sorted((out / d).rglob("*"))
        if p.is_file()
    }


def test_tables_and_mapping_are_byte_identical(runs: dict[str, Path]) -> None:
    head, ref = _files(runs["head"]), _files(runs["ref"])
    assert ref, "the reference run wrote no tables/ or mapping/ files"
    assert {k.split("/")[0] for k in ref} == set(IDENTICAL)
    assert sorted(head) == sorted(ref)
    differ = [k for k in ref if head[k].read_bytes() != ref[k].read_bytes()]
    assert not differ, differ


def _published(out: Path) -> dict[str, bytes]:
    return {
        p.relative_to(out).as_posix(): p.read_bytes()
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.relative_to(out).parts[0] not in ("work", *IDENTICAL)
    }


def test_other_published_files_differ_only_as_listed(runs: dict[str, Path]) -> None:
    head, ref = _published(runs["head"]), _published(runs["ref"])
    assert sorted(head) == sorted(ref)
    differ = {k for k in ref if head[k] != ref[k]}
    assert differ <= DIFFER, sorted(differ - DIFFER)
    assert {"run_manifest.json", "preflight/preflight.json", "report/report.html"} <= differ


def _npz_meta(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as npz:
        meta: dict[str, Any] = json.loads(str(npz["meta"]))
    return meta


def test_sixframe_meta_format_3_at_table_11(runs: dict[str, Path]) -> None:
    head = _npz_meta(runs["head"] / "work" / "sixframe")
    ref = _npz_meta(runs["ref"] / "work" / "sixframe")
    assert (ref["format"], "translation_table" in ref) == (2, False)
    assert (head["format"], head["translation_table"]) == (3, 11)
    rest = {k: v for k, v in head.items() if k not in ("format", "translation_table")}
    assert rest == {k: v for k, v in ref.items() if k != "format"}


def _json(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


def test_manifest_records_table_11_and_bumps_the_schema_minor(runs: dict[str, Path]) -> None:
    head = _json(runs["head"] / "run_manifest.json")
    ref = _json(runs["ref"] / "run_manifest.json")
    assert head["annotation"]["translation_table"] == 11
    assert "translation_table" not in (ref.get("annotation") or {})
    h_major, h_minor = map(int, str(head["output_schema_version"]).split("."))
    r_major, r_minor = map(int, str(ref["output_schema_version"]).split("."))
    assert (h_major, h_minor) == (r_major, r_minor + 1)
    assert head["config"]["annotation"]["translation_table"] == 11
    h_cfg, r_cfg = head["config"], ref["config"]
    tb = (
        h_cfg["mapping"]["proteins"]["tblastn_args"],
        r_cfg["mapping"]["proteins"]["tblastn_args"],
    )
    assert tb == TBLASTN  # D1: -db_gencode now comes from the table
    h_cfg["mapping"]["proteins"].pop("tblastn_args")
    r_cfg["mapping"]["proteins"].pop("tblastn_args")
    assert h_cfg == r_cfg  # otherwise the config snapshot is unchanged
    skip = {
        "annotation",
        "output_schema_version",
        "config",
        "provenance_digest",
        "outdir",
        "command_line",
        "git_commit",
        "git_dirty",
        "steps",
        "started_at",
        "finished_at",
    }
    assert {k: v for k, v in head.items() if k not in skip} == {
        k: v for k, v in ref.items() if k not in skip
    }
    assert {k: v for k, v in head["annotation"].items() if k != "translation_table"} == (
        ref.get("annotation") or {}
    )
    assert [s["name"] for s in head["steps"]] == [s["name"] for s in ref["steps"]]


def _checks(out: Path) -> dict[str, dict[str, Any]]:
    doc = _json(out / "preflight" / "preflight.json")
    return {c["rule_id"]: c for c in doc["checks"]}


def test_preflight_differs_only_in_the_table_rows(runs: dict[str, Path]) -> None:
    head, ref = _checks(runs["head"]), _checks(runs["ref"])
    assert set(head) - set(ref) == {"config.tool_genetic_code"}  # 2.3, new at every table
    assert head["config.tool_genetic_code"]["level"] == "ok"
    assert set(ref) <= set(head)
    h, r = head["genome.genetic_code"], ref["genome.genetic_code"]
    assert (h["level"], h["exit_code"]) == (r["level"], r["exit_code"])
    assert h["value"] == {**r["value"], "translation_table": 11}
    assert h["rule"] != r["rule"]
    tt_h, tt_r = head["config.translation_table"], ref["config.translation_table"]
    assert (tt_h["level"], tt_h["value"]) == (tt_r["level"], tt_r["value"]) == ("ok", 11)
    same = set(ref) - {"genome.genetic_code", "config.translation_table"}
    assert {k: head[k] for k in same} == {k: ref[k] for k in same}


def test_report_methods_sentence_and_no_note_at_11(runs: dict[str, Path]) -> None:
    head = (runs["head"] / "report" / "report.html").read_text(encoding="utf-8")
    ref = (runs["ref"] / "report" / "report.html").read_text(encoding="utf-8")
    assert METHODS in head and METHODS not in ref
    assert "data-genetic-code" not in head
