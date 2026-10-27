"""MVP-03c content blocks through ``gippa run --executor local`` and ``gippa _step manifest``
(docs/specs/MVP-03-manifest-finalize.md D-B step 5, D-C rows marked c, edge cases "No --db-dir",
"A step skipped", acceptance 03c; lead decision 10; ADR-0017 KEGG date read by M9-16).

A toy registry (``tests/unit/conftest.py::registry`` plus a reload of ``gippa.steps.manifest``):
``ctx`` writes the work keys ``island_diagnostics`` and ``sensitivity`` (published by
``gippa.runner.local.publish`` to ``OUT/context/island_diagnostics.parquet`` and
``OUT/tables/sensitivity.tsv``), logs an ``operons.missed_operons`` and a
``mapping.search_db_resolved`` event and calls ``gippa.context.oriter.oric_position`` on two oriC
calls (the real ``multiple_oric`` event). Preflight and the CLI are real. The database directory
is a tmp ``db_manifest.json`` in the ``gippa setup`` item form with one ``retrieved`` date per
item (never the real ``$GIPPA_DB_DIR``).

The run tail is expected to read those files and pass ``gippa.data.calibration.load()`` and the
db manifest to the builder; ``databases`` is also in ``OUT/work/provenance.json`` (D-D).
"""

from __future__ import annotations

import importlib
import json
import math
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from click.testing import CliRunner
from helpers.manifest_data import validate

from gippa.exitcodes import ExitCode

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import oriter_data as od

SENS_TSV = ("family\tvariant\tn_primary\tn_variant\tn_shared\tjaccard\tagreement\n"
            "regions\tmissing_min\t3\t4\t2\t0.4\tNA\n"
            "regions\ttop3\t3\t3\t3\t1\tNA\n")  # fmt: skip
BINS = [
    {"bin": 0, "lower_bp": None, "upper_bp": 0.0, "N_NO": 4, "N_cv": 3, "N_sw": 1},
    {"bin": 1, "lower_bp": 0.0, "upper_bp": None, "N_NO": 7, "N_cv": 5, "N_sw": 2},
]
DIAG = [("chr1", True, "auto", 3, 7, 0.21, []),
        ("lin1", False, "not_scanned", None, None, math.nan, ["few_values"])]  # fmt: skip
#: retrieval dates pinned per item (M9-16 prints the ``kegg*`` ones)
KEGG_DATES = {"kegg_ko_pathway": "2026-10-01", "kegg_pathway_names": "2026-10-02"}


def _items(kegg: dict[str, str] | None = None) -> list[dict[str, Any]]:
    """The real-format fixture (every resource preflight needs), one date and checksum each."""
    from helpers.fake_env import setup_manifest

    dates = {**KEGG_DATES, **(kegg or {})}
    items = [dict(i) for i in setup_manifest()["items"]]  # type: ignore[attr-defined]
    for n, item in enumerate(sorted(items, key=lambda i: i["name"])):
        item["retrieved"] = dates.get(item["name"], f"2026-09-{n + 1:02d}")
        item["checksum"] = "sha256:" + format(n % 16, "x") * 64
        item["size_bytes"] = 10 * (n + 1)
    return items


def _db(root: Path, kegg: dict[str, str] | None = None) -> Path:
    items = _items(kegg)
    items.append({"name": "zz_failed", "version": "1.0", "url": "https://example.invalid/z",
                  "checksum": "", "size_bytes": 0, "retrieved": "2026-09-03", "path": "z",
                  "status": "failed", "message": "download failed"})  # fmt: skip
    root.mkdir(parents=True, exist_ok=True)
    (root / "db_manifest.json").write_text(
        json.dumps({"schema_version": "1.0", "items": items}, indent=2), encoding="utf-8"
    )
    return root


def _write_diag(path: Path) -> None:
    schema = pa.schema([("contig", pa.string()), ("scanned", pa.bool_()),
                        ("block_rule", pa.string()), ("l_hat", pa.int64()),
                        ("block_len", pa.int64()), ("morans_i_lag1", pa.float64()),
                        ("flags", pa.list_(pa.string()))])  # fmt: skip
    cols = list(zip(*DIAG, strict=True))
    pq.write_table(pa.Table.from_arrays([pa.array(list(c), type=f.type) for c, f in
                                         zip(cols, schema, strict=True)], schema=schema),
                   path)  # fmt: skip


def _toy(reg: ModuleType, *, with_content: bool = True) -> None:
    def x(ctx: Any) -> None:
        Path(ctx.outputs["b.txt"]).write_text("b\n")

    def content(ctx: Any) -> None:
        from gippa.context import oriter

        _write_diag(Path(ctx.outputs["island_diagnostics"]))
        Path(ctx.outputs["sensitivity"]).write_text(SENS_TSV)
        ctx.log.info("bins", extra={"event": "operons.missed_operons",
                                    "data": {"pi0": 0.8, "bin_counts": BINS}})  # fmt: skip
        ctx.log.info("search_db", extra={"event": "mapping.search_db_resolved", "data": {
            "key": "search_db", "search_db": "same", "n_queries": 12}})  # fmt: skip
        oriter.oric_position(od.oric_at("chr1", 1_000, 200_000), "chr1", 4_641_652)

    reg.step("x", inputs=[], outputs=["b.txt"])(x)
    if with_content:
        reg.step("ctx", inputs=[], outputs=["island_diagnostics", "sensitivity"])(content)
    importlib.reload(importlib.import_module("gippa.steps.manifest"))


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    genome = tmp_path / "in" / "g.fa"
    genome.parent.mkdir()
    genome.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    pep = tmp_path / "in" / "p.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    from helpers.fake_env import fake_gippa_root

    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    db = _db(tmp_path / "test_db")
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    return {"genome": genome, "peptides": pep, "db": db}


def _argv(env: dict[str, Path], out: Path, *extra: str) -> list[str]:
    return ["run", "--genome", str(env["genome"]), "--peptides", str(env["peptides"]),
            "--outdir", str(out), "--executor", "local", "--modules", "core", "--seed", "5",
            *extra]  # fmt: skip


def _invoke(argv: list[str]) -> tuple[int, str]:
    from gippa.cli import cli

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
def ok_run(registry: ModuleType, env: dict[str, Path], tmp_path: Path) -> Iterator[Path]:
    _toy(registry)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    yield out


# ----------------------------------------------------------------------------- blocks


def _calibration() -> dict[str, Any]:
    from gippa.data.calibration import load

    return load()


def _uncalibrated_flags() -> set[str]:
    """The ``uncalibrated:<key>`` run flags calibration.yaml implies (M6-12 C0)."""
    return {f"uncalibrated:{k}" for k, e in _calibration().items() if not e["calibrated"]}


def test_content_manifest_validates(ok_run: Path) -> None:
    m = _manifest(ok_run)
    validate(m)
    assert {"sensitivity", "context", "operons", "mapping", "calibration"} <= set(m)
    assert "Infinity" not in (ok_run / "run_manifest.json").read_text(encoding="utf-8")


def test_sensitivity_block_from_the_published_table(ok_run: Path) -> None:
    assert (ok_run / "tables" / "sensitivity.tsv").is_file()
    s = _manifest(ok_run)["sensitivity"]
    assert s["level"] == "basic" and s["n_rows"] == 2
    assert s["families"] == {"regions": {
        "missing_min": {"n_primary": 3, "n_variant": 4, "n_shared": 2, "jaccard": 0.4},
        "top3": {"n_primary": 3, "n_variant": 3, "n_shared": 3, "jaccard": 1.0}}}  # fmt: skip


def test_context_replicons_equal_the_published_parquet(ok_run: Path) -> None:
    (diag,) = list(ok_run.rglob("island_diagnostics.parquet"))[:1]
    want = [{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()}
            for r in pq.read_table(diag).to_pylist()]  # fmt: skip
    assert _manifest(ok_run)["context"]["replicons"] == want
    assert [r["contig"] for r in want] == ["chr1", "lin1"]


def test_context_oric_from_the_real_oriter_event(ok_run: Path) -> None:
    evs = [json.loads(ln) for ln in (ok_run / "logs" / "events.jsonl").read_text().splitlines()
           if ln.strip()]  # fmt: skip
    ori = [e for e in evs if e["event"] == "context_stats.oriter.multiple_oric"]
    assert len(ori) == 1, "the oriter event did not reach logs/events.jsonl (local executor)"
    assert ori[0]["step"] == "context_stats"
    (o,) = _manifest(ok_run)["context"]["oric"]
    assert o == ori[0]["data"]


def test_operons_and_mapping_blocks(ok_run: Path) -> None:
    m = _manifest(ok_run)
    assert m["operons"] == {"pi0": 0.8, "bin_counts": BINS}
    assert m["mapping"] == {"search_db": "same", "n_queries": 12}


def test_calibration_block_and_flags(ok_run: Path) -> None:
    from gippa.data.calibration import load

    m = _manifest(ok_run)
    cal = load()
    assert set(m["calibration"]) == set(cal)
    for key, entry in cal.items():
        assert m["calibration"][key]["calibrated"] is entry["calibrated"]
        assert m["calibration"][key]["provenance"] == entry.get("provenance")
    assert _uncalibrated_flags() == {f for f in m["flags"] if f.startswith("uncalibrated:")}
    assert m["flags"] == sorted(set(m["flags"]))


# ----------------------------------------------------------------------------- databases


def test_databases_from_the_setup_manifest_with_retrieved(ok_run: Path) -> None:
    dbs = _manifest(ok_run)["databases"]
    want = sorted(_items(), key=lambda i: i["name"])  # every ok item; the failed one left out
    assert [d["name"] for d in dbs] == [i["name"] for i in want]
    for d, item in zip(dbs, want, strict=True):
        assert d == {k: item[k] for k in ("name", "version", "checksum", "retrieved",
                                          "size_bytes")}  # fmt: skip


def test_kegg_dates_in_manifest_and_provenance_core(ok_run: Path) -> None:
    """M9-16 reads the KEGG date from ``databases[]`` of the manifest or of
    ``work/provenance.json`` (the core it stages into ``report``)."""
    m = _manifest(ok_run)
    core = json.loads((ok_run / "work" / "provenance.json").read_text(encoding="utf-8"))
    assert core["databases"] == m["databases"]
    kegg = {d["name"]: d["retrieved"] for d in m["databases"] if d["name"].startswith("kegg")}
    assert kegg == KEGG_DATES


def test_db_dir_flag_overrides_the_environment(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    other = _db(tmp_path / "other_db", kegg={"kegg_ko_pathway": "2027-01-02"})
    _toy(registry, with_content=False)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out, "--db-dir", str(other)))
    assert code == ExitCode.SUCCESS, log
    kegg = {d["name"]: d["retrieved"] for d in _manifest(out)["databases"]
            if d["name"].startswith("kegg")}  # fmt: skip
    assert kegg == {"kegg_ko_pathway": "2027-01-02", "kegg_pathway_names": "2026-10-02"}


def test_no_db_dir_gives_no_databases(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Edge case "No --db-dir and no $GIPPA_DB_DIR": ``databases == []`` (M9-16: not recorded)."""
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    _toy(registry, with_content=False)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    m = _manifest(out)
    assert m["databases"] == []
    validate(m)


# ----------------------------------------------------------------------------- skipped content


def test_without_content_steps_the_blocks_are_empty(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    """Edge case "A step skipped": no table, no event → empty blocks, no flag invented; the
    calibration block is there all the same."""
    _toy(registry, with_content=False)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    m = _manifest(out)
    assert m["sensitivity"] == {"level": "basic", "n_rows": 0, "families": {}}
    assert m["context"] == {"replicons": [], "oric": []}
    assert m["operons"] is None and m["mapping"] is None
    assert set(m["calibration"]) == set(_calibration())
    validate(m)


def test_sensitivity_off_recorded(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    """Backlog M4-19 d7 / MVP-01 pass 2: ``output.sensitivity = off`` is in the manifest."""
    _toy(registry, with_content=False)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out, "--set", 'output.sensitivity="off"'))
    assert code == ExitCode.SUCCESS, log
    assert _manifest(out)["sensitivity"] == {"level": "off", "n_rows": 0, "families": {}}


def test_two_runs_have_equal_content_blocks(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    from gippa.manifest import stable_view

    _toy(registry)
    a, b = tmp_path / "A" / "OUT", tmp_path / "B" / "OUT"
    for out in (a, b):
        code, log = _invoke(_argv(env, out))
        assert code == ExitCode.SUCCESS, log
    va, vb = stable_view(_manifest(a)), stable_view(_manifest(b))
    for key in ("sensitivity", "context", "operons", "mapping", "calibration", "databases"):
        assert va[key] == vb[key], key
    assert va == vb


# ----------------------------------------------------------------------------- registry path


def test_step_manifest_reads_the_optional_content_inputs(
    ok_run: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``gippa _step manifest`` with ``sensitivity`` and ``island_diagnostics`` (its declared
    optional inputs) fills the same blocks as the run tail, calibration included."""
    task = tmp_path / "task"
    task.mkdir()
    monkeypatch.chdir(task)
    target = task / "run_manifest.json"
    (diag,) = list(ok_run.rglob("island_diagnostics.parquet"))[:1]
    code, log = _invoke([
        "_step", "manifest",
        "--in", f"preflight={ok_run / 'preflight' / 'preflight.json'}",
        "--in", f"events={ok_run / 'logs' / 'events.jsonl'}",
        "--in", f"sensitivity={ok_run / 'tables' / 'sensitivity.tsv'}",
        "--in", f"island_diagnostics={diag}",
        "--out", f"run_manifest={target}",
    ])  # fmt: skip
    assert code == ExitCode.SUCCESS, log
    m = json.loads(target.read_text(encoding="utf-8"))
    validate(m)
    tail = _manifest(ok_run)
    for key in ("sensitivity", "context", "operons", "mapping", "calibration"):
        assert m[key] == tail[key], key
    assert _uncalibrated_flags() <= set(m["flags"])
