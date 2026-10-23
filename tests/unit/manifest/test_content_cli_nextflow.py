"""MVP-03c content blocks through ``giap run`` with the Nextflow executor (a fake ``nextflow``)
(docs/specs/MVP-03-manifest-finalize.md D-C rows marked c, D-B step 5, acceptance 03c "the
oriter event reaches events.jsonl under both executors", "``context.replicons`` equals
``island_diagnostics.parquet`` row for row"; lead decisions 9 and 10; ADR-0017 / M9-16 KEGG date).

The fake ``nextflow`` (as ``test_manifest_cli_nextflow.py``) exits ``$FAKE_NF_RC``, copies
``$FAKE_NF_TRACE`` to the ``-with-trace`` path, publishes ``$FAKE_NF_STEP_EVENTS`` as
``OUT/logs/steps/context_stats.events.jsonl`` (what the CONTEXT_STATS ``GIAP_STEP`` task
publishes; ``append_step_events`` copies it into ``OUT/logs/events.jsonl``) and copies the tree
``$FAKE_NF_FILES`` into ``OUT`` (the published ``context/island_diagnostics.parquet`` and
``tables/sensitivity.tsv``, the ``GIAP_STEP`` publishDir names). The step events file is the
lines a fixed ``giap _step context_stats`` would write: the real ``multiple_oric`` event,
captured through ``giap.cli._record_step_events``, plus an ``operons.missed_operons`` line with
the ``-Infinity`` token a pre-fix task wrote, and a ``mapping.search_db_resolved`` line. Nothing
real is launched; ``$GIAP_ROOT`` and ``$GIAP_DB_DIR`` are fakes from ``helpers.fake_env``.
"""

from __future__ import annotations

import json
import math
import os
import random
import shutil
import stat
import sys
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, setup_manifest
from helpers.manifest_data import validate

from giap.exitcodes import ExitCode

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import oriter_data as od

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
if [ -n "$FAKE_NF_TRACE" ]; then
  prev=""
  for a in "$@"; do
    if [ "$prev" = "-with-trace" ]; then mkdir -p "$(dirname "$a")"; cp "$FAKE_NF_TRACE" "$a"; fi
    prev="$a"
  done
fi
if [ -n "$FAKE_NF_STEP_EVENTS" ]; then
  mkdir -p logs/steps
  cp "$FAKE_NF_STEP_EVENTS" logs/steps/context_stats.events.jsonl
fi
if [ -n "$FAKE_NF_FILES" ]; then cp -R "$FAKE_NF_FILES"/. .; fi
exit "${FAKE_NF_RC:-0}"
"""
COLS = ["task_id", "hash", "native_id", "name", "tag", "status", "exit", "attempt", "realtime",
        "%cpu", "peak_rss", "workdir", "container"]  # fmt: skip
SENS_TSV = ("family\tvariant\tn_primary\tn_variant\tn_shared\tjaccard\tagreement\n"
            "islands\tmissing_min\t5\t6\t4\t0.571429\tNA\n")  # fmt: skip
DIAG = [("chr1", True, "auto", 2, 5, 0.33, []),
        ("pA", False, "not_scanned", None, None, math.nan, ["few_values"])]  # fmt: skip
KEGG_DATES = {"kegg_ko_pathway": "2026-10-01", "kegg_pathway_names": "2026-10-02"}
MISSED_LINE = ('{"data": {"bin_counts": [{"N_NO": 4, "N_cv": 3, "N_sw": 1, "bin": 0, '
               '"lower_bp": -Infinity, "upper_bp": 0.0}, {"N_NO": 7, "N_cv": 5, "N_sw": 2, '
               '"bin": 1, "lower_bp": 0.0, "upper_bp": Infinity}], "pi0": 1.0}, '
               '"event": "operons.missed_operons", "level": "INFO", "step": "operons", '
               '"ts": "2026-10-22T00:00:00+00:00"}')  # fmt: skip
SEARCH_LINE = json.dumps({"ts": "2026-10-22T00:00:00+00:00", "step": "placement",
                          "level": "INFO", "event": "mapping.search_db_resolved",
                          "data": {"key": "search_db", "search_db": "other",
                                   "n_queries": 31}}, sort_keys=True)  # fmt: skip


def _db(root: Path) -> Path:
    items = [dict(i) for i in setup_manifest()["items"]]  # type: ignore[attr-defined]
    for n, item in enumerate(sorted(items, key=lambda i: i["name"])):
        item["retrieved"] = KEGG_DATES.get(item["name"], f"2026-08-{n + 1:02d}")
        item["checksum"] = "sha256:" + format(n % 16, "x") * 64
    root.mkdir(parents=True, exist_ok=True)
    (root / "db_manifest.json").write_text(json.dumps({"schema_version": "1.0", "items": items}))
    return root


def _write_diag(path: Path) -> None:
    schema = pa.schema([("contig", pa.string()), ("scanned", pa.bool_()),
                        ("block_rule", pa.string()), ("l_hat", pa.int64()),
                        ("block_len", pa.int64()), ("morans_i_lag1", pa.float64()),
                        ("flags", pa.list_(pa.string()))])  # fmt: skip
    cols = list(zip(*DIAG, strict=True))
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_arrays([pa.array(list(c), type=f.type) for c, f in
                                         zip(cols, schema, strict=True)], schema=schema),
                   path)  # fmt: skip


def _oric_lines(tmp_path: Path) -> list[str]:
    """What ``giap _step context_stats`` writes to the task's ``logs/events.jsonl``."""
    from giap.cli import _record_step_events
    from giap.context import oriter

    logs = tmp_path / "task_logs"
    restore = _record_step_events(logs)
    try:
        oriter.oric_position(od.oric_at("chr1", 1_000, 200_000), "chr1", 4_641_652)
    finally:
        restore()
    path = logs / "events.jsonl"
    return [ln for ln in path.read_text().splitlines() if ln.strip()] if path.is_file() else []


@pytest.fixture
def nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "giap_root")))
    monkeypatch.setenv("GIAP_DB_DIR", str(_db(tmp_path / "db")))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    for var in ("FAKE_NF_RC", "FAKE_NF_TRACE", "FAKE_NF_STEP_EVENTS", "FAKE_NF_FILES"):
        monkeypatch.delenv(var, raising=False)
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(seq[i:i + 60]
                      for i in range(0, 6000, 60)) + "\n")  # fmt: skip
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text("".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n"
                                for i in range(3)))  # fmt: skip
    files = tmp_path / "published"
    _write_diag(files / "context" / "island_diagnostics.parquet")
    (files / "tables").mkdir(parents=True)
    (files / "tables" / "sensitivity.tsv").write_text(SENS_TSV)
    task = tmp_path / "tasks" / "cs"
    task.mkdir(parents=True)
    (task / "versions.yml").write_text('"CONTEXT:CONTEXT_STATS":\n    giap: 0.1.0\n')
    trace = tmp_path / "trace_src.tsv"
    row = {c: "-" for c in COLS} | {"task_id": "1", "attempt": "1", "exit": "0",
                                    "name": "CONTEXT:CONTEXT_STATS (context_stats)",
                                    "status": "COMPLETED", "realtime": "3s", "%cpu": "90.0%",
                                    "peak_rss": "100 MB", "workdir": str(task)}  # fmt: skip
    trace.write_text("\t".join(COLS) + "\n" + "\t".join(row[c] for c in COLS) + "\n")
    events = tmp_path / "step_events.jsonl"
    events.write_text("\n".join([*_oric_lines(tmp_path), MISSED_LINE, SEARCH_LINE]) + "\n")
    monkeypatch.setenv("FAKE_NF_TRACE", str(trace))
    monkeypatch.setenv("FAKE_NF_STEP_EVENTS", str(events))
    monkeypatch.setenv("FAKE_NF_FILES", str(files))
    return {"genome": genome, "proteins": proteins, "events": events, "files": files}


def _run(nf: dict[str, Path], out: Path, *extra: str) -> tuple[int, str]:
    from giap.cli import cli

    argv = ["run", "--genome", str(nf["genome"]), "--proteins", str(nf["proteins"]),
            "--outdir", str(out), *extra]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}; OUT has {sorted(p.name for p in out.glob('*'))}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


@pytest.fixture
def ok(nf: dict[str, Path], tmp_path: Path) -> Path:
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    return out


# ----------------------------------------------------------------------------- tests


def test_the_step_events_file_holds_the_oriter_event(nf: dict[str, Path]) -> None:
    """The Nextflow half of "the oriter event reaches events.jsonl": the task-local handler that
    ``giap _step`` installs records it."""
    evs = [json.loads(ln) for ln in nf["events"].read_text().splitlines()[:-2]]
    assert [e["event"] for e in evs] == ["context_stats.oriter.multiple_oric"], evs
    assert evs[0]["step"] == "context_stats"


def test_manifest_validates_and_is_strict(ok: Path) -> None:
    m = _manifest(ok)
    validate(m)
    assert m["executor"] == "nextflow"
    assert "Infinity" not in (ok / "run_manifest.json").read_text(encoding="utf-8")


def test_oriter_event_in_run_events_and_manifest(ok: Path) -> None:
    evs = [json.loads(ln) for ln in (ok / "logs" / "events.jsonl").read_text().splitlines()
           if ln.strip() and "Infinity" not in ln]  # fmt: skip
    ori = [e for e in evs if e["event"] == "context_stats.oriter.multiple_oric"]
    assert len(ori) == 1
    (o,) = _manifest(ok)["context"]["oric"]
    assert o == ori[0]["data"]


def test_context_replicons_equal_the_published_parquet(ok: Path) -> None:
    p = ok / "context" / "island_diagnostics.parquet"
    assert p.is_file()
    want = [{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()}
            for r in pq.read_table(p).to_pylist()]  # fmt: skip
    assert _manifest(ok)["context"]["replicons"] == want


def test_operons_minus_infinity_from_a_task_becomes_null(ok: Path) -> None:
    ops = _manifest(ok)["operons"]
    assert ops["pi0"] == 1.0
    assert ops["bin_counts"][0]["lower_bp"] is None
    assert ops["bin_counts"][1]["upper_bp"] is None
    assert ops["bin_counts"][1]["N_sw"] == 2


def test_mapping_and_sensitivity_blocks(ok: Path) -> None:
    m = _manifest(ok)
    assert m["mapping"] == {"search_db": "other", "n_queries": 31}
    assert m["sensitivity"]["n_rows"] == 1
    cell = m["sensitivity"]["families"]["islands"]["missing_min"]
    assert (cell["n_primary"], cell["n_variant"], cell["n_shared"]) == (5, 6, 4)
    assert cell["jaccard"] == pytest.approx(0.571429)


def test_calibration_block_and_flags(ok: Path) -> None:
    m = _manifest(ok)
    assert set(m["calibration"]) == {"operon_distance", "terminatornet_yes"}
    assert {"uncalibrated:operon_distance", "uncalibrated:terminatornet_yes"} <= set(m["flags"])


def test_databases_with_kegg_dates(ok: Path) -> None:
    m = _manifest(ok)
    names = [d["name"] for d in m["databases"]]
    assert names and names == sorted(names)
    kegg = {d["name"]: d["retrieved"] for d in m["databases"] if d["name"].startswith("kegg")}
    assert kegg == KEGG_DATES
    core = json.loads((ok / "work" / "provenance.json").read_text(encoding="utf-8"))
    assert core["databases"] == m["databases"]


def test_giap_only_versions_give_no_tools(ok: Path) -> None:
    assert _manifest(ok)["tools"] == []


def test_failed_run_still_has_the_content_known_so_far(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Edge case "A failed run": the manifest has the blocks known so far."""
    monkeypatch.setenv("FAKE_NF_RC", "1")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    validate(m)
    assert m["status"] == "failed"
    assert m["mapping"] == {"search_db": "other", "n_queries": 31}
    assert set(m["calibration"]) == {"operon_distance", "terminatornet_yes"}


def test_without_published_content_the_blocks_are_empty(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("FAKE_NF_FILES")
    monkeypatch.delenv("FAKE_NF_STEP_EVENTS")
    shutil.rmtree(nf["files"])
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    m = _manifest(out)
    assert m["context"] == {"replicons": [], "oric": []}
    assert m["operons"] is None and m["mapping"] is None
    assert m["sensitivity"] == {"level": "basic", "n_rows": 0, "families": {}}
    validate(m)
