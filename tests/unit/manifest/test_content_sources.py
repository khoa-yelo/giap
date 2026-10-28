"""MVP-03c readers and the two producer fixes (docs/specs/MVP-03-manifest-finalize.md "The gap"
bullets on ``oriter`` and ``-Infinity``, D-C rows ``sensitivity`` / ``context.replicons`` /
``context.oric`` / ``operons`` / ``databases``, edge cases "No --db-dir", "db_manifest.json in the
provisional form"; acceptance 03c "the oriter event reaches events.jsonl under both executors";
backlog M6-06/M6-08 "must-do before the manifest writer").

Interface assumed (test-author names, ``gippa.io.manifest_sources``):

* ``read_sensitivity(path: Path | None) -> list[dict] | None``: the rows of
  ``tables/sensitivity.tsv`` (counts as ``int``, ``jaccard``/``agreement`` as ``float`` or None
  for ``NA``), ``[]`` for a header-only file, None when ``path`` is None or absent.
* ``read_island_diagnostics(path: Path | None) -> list[dict] | None``: the rows of
  ``context/island_diagnostics.parquet`` in file (C1) order with the seven
  ``DIAGNOSTICS_COLUMNS`` keys (``flags`` a list), None when absent.
* ``read_db_manifest(db_dir: Path | None) -> dict | None``: the parsed
  ``<db_dir>/db_manifest.json`` (either shape), None when ``db_dir`` is None or the file is absent.

Producer fixes pinned:

* ``gippa.context.oriter``'s logger is a child of ``gippa.steps`` whose leaf is ``context_stats``,
  so the ``gippa.steps`` events handler (``gippa run`` local and ``gippa _step`` inside a Nextflow
  task both install it, ``gippa.cli._record_step_events``) writes the
  ``context_stats.oriter.multiple_oric`` event with ``step == "context_stats"``.
* ``gippa.steps.operons._plain`` maps ±∞ to None (the event line is strict JSON) and keeps every
  other value (numpy scalars → Python).
"""

from __future__ import annotations

import json
import logging
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from helpers.manifest_data import build, event

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import oriter_data as od

DIAG_COLS = ("contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags")
SENS_HEADER = "family\tvariant\tn_primary\tn_variant\tn_shared\tjaccard\tagreement\n"
ORIC_EVENT = "context_stats.oriter.multiple_oric"


def _ms():  # type: ignore[no-untyped-def]
    from gippa.io import manifest_sources

    return manifest_sources


def write_diagnostics(path: Path, rows: list[tuple[Any, ...]]) -> None:
    """The ``_write_diagnostics`` Arrow schema (``gippa.steps.context_stats``)."""
    schema = pa.schema([("contig", pa.string()), ("scanned", pa.bool_()),
                        ("block_rule", pa.string()), ("l_hat", pa.int64()),
                        ("block_len", pa.int64()), ("morans_i_lag1", pa.float64()),
                        ("flags", pa.list_(pa.string()))])  # fmt: skip
    cols = list(zip(*rows, strict=True)) if rows else [[] for _ in DIAG_COLS]
    table = pa.Table.from_arrays([pa.array(list(c), type=f.type) for c, f in
                                  zip(cols, schema, strict=True)], schema=schema)  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path)


DIAG_ROWS = [
    ("chr1", True, "auto", 4, 9, 0.31, []),
    ("pA", True, "fixed", None, 5, 0.02, []),
    ("lin1", False, "not_scanned", None, None, math.nan, ["few_values"]),
]


# ----------------------------------------------------------------------------- sensitivity


def test_read_sensitivity_rows(tmp_path: Path) -> None:
    p = tmp_path / "sensitivity.tsv"
    p.write_text(SENS_HEADER + "regions\tmissing_min\t3\t4\t2\t0.4\tNA\n"
                 "islands\ttop3\t0\t0\t0\tNA\tNA\n")  # fmt: skip
    rows = _ms().read_sensitivity(p)
    assert rows == [
        {"family": "regions", "variant": "missing_min", "n_primary": 3, "n_variant": 4,
         "n_shared": 2, "jaccard": 0.4, "agreement": None},
        {"family": "islands", "variant": "top3", "n_primary": 0, "n_variant": 0,
         "n_shared": 0, "jaccard": None, "agreement": None},
    ]  # fmt: skip
    assert all(isinstance(r["n_primary"], int) for r in rows)


def test_read_sensitivity_header_only_is_empty(tmp_path: Path) -> None:
    p = tmp_path / "sensitivity.tsv"
    p.write_text(SENS_HEADER)
    assert _ms().read_sensitivity(p) == []


def test_read_sensitivity_absent_is_none(tmp_path: Path) -> None:
    assert _ms().read_sensitivity(tmp_path / "nope.tsv") is None
    assert _ms().read_sensitivity(None) is None


def test_example_8_off_from_files_header_only_and_absent(tmp_path: Path) -> None:
    """Example 8 through the reader: both cases give the same block."""
    from gippa.config import GippaConfig

    cfg = GippaConfig.model_validate({"seed": 7, "output": {"sensitivity": "off"}})
    p = tmp_path / "sensitivity.tsv"
    p.write_text(SENS_HEADER)
    a = build(cfg=cfg, sensitivity=_ms().read_sensitivity(p))
    b = build(cfg=cfg, sensitivity=_ms().read_sensitivity(tmp_path / "absent.tsv"))
    assert a["sensitivity"] == b["sensitivity"] == {"level": "off", "n_rows": 0, "families": {}}


def test_read_sensitivity_written_by_the_real_writer(tmp_path: Path) -> None:
    import pandas as pd

    from gippa.context.sensitivity import COLUMNS
    from gippa.io.sensitivity_tsv import write_sensitivity_tsv

    df = pd.DataFrame([["regions", "top3", 5, 6, 4, 4 / 7, math.nan]], columns=list(COLUMNS))
    p = tmp_path / "sensitivity.tsv"
    write_sensitivity_tsv(p, df)
    (row,) = _ms().read_sensitivity(p)
    assert row["jaccard"] == pytest.approx(4 / 7, rel=1e-5)  # 6 significant digits
    assert row["agreement"] is None
    m = build(sensitivity=[row])
    assert m["sensitivity"]["families"]["regions"]["top3"]["n_shared"] == 4


# ----------------------------------------------------------------------------- diagnostics


def test_read_island_diagnostics_rows_in_file_order(tmp_path: Path) -> None:
    p = tmp_path / "context" / "island_diagnostics.parquet"
    write_diagnostics(p, DIAG_ROWS)
    rows = _ms().read_island_diagnostics(p)
    assert rows is not None
    assert [r["contig"] for r in rows] == ["chr1", "pA", "lin1"]
    for r, src in zip(rows, DIAG_ROWS, strict=True):
        assert set(r) == set(DIAG_COLS)
        assert list(r["flags"]) == src[6]
        assert r["l_hat"] == src[3] and r["block_len"] == src[4]


def test_read_island_diagnostics_absent_is_none(tmp_path: Path) -> None:
    assert _ms().read_island_diagnostics(tmp_path / "none.parquet") is None
    assert _ms().read_island_diagnostics(None) is None


def test_context_replicons_equal_the_parquet_row_for_row(tmp_path: Path) -> None:
    """Acceptance 03c: ``context.replicons`` equals ``island_diagnostics.parquet`` row for row
    (NaN → null, the one change strict JSON forces)."""
    p = tmp_path / "island_diagnostics.parquet"
    write_diagnostics(p, DIAG_ROWS)
    m = build(island_diagnostics=_ms().read_island_diagnostics(p))
    want = [
        {k: (None if isinstance(v, float) and math.isnan(v) else v)
         for k, v in r.items()}
        for r in pq.read_table(p).to_pylist()
    ]  # fmt: skip
    assert m["context"]["replicons"] == want
    json.dumps(m, allow_nan=False)


def test_context_replicons_of_the_real_writer(tmp_path: Path) -> None:
    """The frame ``context_stats`` persists (``island_diagnostics`` → ``_write_diagnostics``)."""
    import pandas as pd

    from gippa.context.diagnostics import DIAGNOSTICS_COLUMNS
    from gippa.steps.context_stats import _write_diagnostics

    df = pd.DataFrame(
        {"contig": ["chr1", "lin1"], "scanned": [True, False],
         "block_rule": ["auto", "not_scanned"],
         "l_hat": pd.Series([2, None], dtype="Int64"),
         "block_len": pd.Series([5, None], dtype="Int64"),
         "morans_i_lag1": [0.4, math.nan], "flags": [(), ("few_values",)]},
        columns=list(DIAGNOSTICS_COLUMNS),
    )  # fmt: skip
    p = tmp_path / "island_diagnostics.parquet"
    _write_diagnostics(p, df)
    reps = build(island_diagnostics=_ms().read_island_diagnostics(p))["context"]["replicons"]
    assert reps == [
        {"contig": "chr1", "scanned": True, "block_rule": "auto", "l_hat": 2, "block_len": 5,
         "morans_i_lag1": 0.4, "flags": []},
        {"contig": "lin1", "scanned": False, "block_rule": "not_scanned", "l_hat": None,
         "block_len": None, "morans_i_lag1": None, "flags": ["few_values"]},
    ]  # fmt: skip


def test_empty_diagnostics_file_gives_no_replicons(tmp_path: Path) -> None:
    p = tmp_path / "island_diagnostics.parquet"
    write_diagnostics(p, [])
    assert _ms().read_island_diagnostics(p) == []
    assert build(island_diagnostics=[])["context"]["replicons"] == []


# ----------------------------------------------------------------------------- db manifest


def test_read_db_manifest_setup_form(tmp_path: Path) -> None:
    doc = {"schema_version": "1.0", "items": [
        {"name": "kegg_ko_pathway", "version": "2026-09-30", "url": "https://example.invalid/k",
         "checksum": "sha256:" + "1" * 64, "size_bytes": 10, "retrieved": "2026-10-01",
         "path": "kegg", "status": "ok", "message": ""}]}  # fmt: skip
    (tmp_path / "db_manifest.json").write_text(json.dumps(doc))
    got = _ms().read_db_manifest(tmp_path)
    assert got == doc
    (d,) = build(db_manifest=got)["databases"]
    assert d["retrieved"] == "2026-10-01"


def test_read_db_manifest_no_dir_or_no_file(tmp_path: Path) -> None:
    assert _ms().read_db_manifest(None) is None
    assert _ms().read_db_manifest(tmp_path) is None  # directory without db_manifest.json
    assert build(db_manifest=_ms().read_db_manifest(None))["databases"] == []


def test_read_db_manifest_provisional_form(tmp_path: Path) -> None:
    from helpers.mapping_steps import write_test_db

    db = write_test_db(tmp_path / "db")
    got = _ms().read_db_manifest(db)
    assert isinstance(got, dict) and "items" not in got
    dbs = build(db_manifest=got)["databases"]
    assert dbs and [d["name"] for d in dbs] == sorted(got)
    assert all(d["checksum"] is None and d["retrieved"] is None for d in dbs)


# ----------------------------------------------------------------------------- -Infinity


def test_events_file_with_minus_infinity_reaches_the_manifest_as_null(tmp_path: Path) -> None:
    """An events line written before the ``_plain`` fix (``-Infinity`` token): the manifest
    half of backlog M6-06/M6-08 still writes strict JSON."""
    p = tmp_path / "events.jsonl"
    ev = event(
        "operons",
        "operons.missed_operons",
        {
            "pi0": 1.0,
            "bin_counts": [
                {"bin": 0, "lower_bp": -math.inf, "upper_bp": 0.0, "N_NO": 1, "N_cv": 1, "N_sw": 0},
                {"bin": 1, "lower_bp": 0.0, "upper_bp": math.inf, "N_NO": 2, "N_cv": 0, "N_sw": 1},
            ],
        },
    )
    p.write_text(json.dumps(ev) + "\n")
    assert "-Infinity" in p.read_text()
    m = build(events=_ms().read_events(p))
    assert m["operons"]["bin_counts"][0]["lower_bp"] is None
    assert m["operons"]["bin_counts"][1]["upper_bp"] is None
    out = tmp_path / "run_manifest.json"
    _ms().write_manifest(out, m)
    assert "Infinity" not in out.read_text(encoding="utf-8")


@pytest.mark.parametrize("value", [-math.inf, math.inf, np.float64(-np.inf), np.float64(np.inf)])
def test_operons_plain_maps_infinity_to_none(value: Any) -> None:
    from gippa.steps.operons import _plain

    assert _plain(value) is None


@pytest.mark.parametrize(
    ("value", "want"),
    [(np.int64(3), 3), (np.float64(2.5), 2.5), (0.0, 0.0), (-50.0, -50.0), (7, 7), ("x", "x")],
)
def test_operons_plain_keeps_finite_values(value: Any, want: Any) -> None:
    from gippa.steps.operons import _plain

    got = _plain(value)
    assert got == want and type(got) is type(want)


def test_operons_event_line_is_strict_json(tmp_path: Path) -> None:
    """The event the operons step builds from ``missed_operons().bin_counts`` serialises with
    ``allow_nan=False`` (its first bin edge is ``-inf`` today)."""
    import pandas as pd

    from gippa.steps.operons import _plain

    bins = pd.DataFrame({"bin": np.arange(2, dtype=np.int64),
                         "lower_bp": np.array([-np.inf, 0.0]),
                         "upper_bp": np.array([0.0, np.inf]),
                         "N_NO": np.array([1, 2], dtype=np.int64),
                         "N_cv": np.array([0, 3], dtype=np.int64),
                         "N_sw": np.array([4, 0], dtype=np.int64)})  # fmt: skip
    data = {
        "pi0": 1.0,
        "bin_counts": [{k: _plain(v) for k, v in r.items()} for r in bins.to_dict("records")],
    }
    text = json.dumps(data, allow_nan=False)
    assert json.loads(text)["bin_counts"][0]["lower_bp"] is None


# ----------------------------------------------------------------------------- oriter logger


def _two_calls() -> Any:
    return od.oric_at("chr1", 1_000, 200_000)  # within 0.05 L of K-12 length: averaged


def test_oriter_logger_is_under_gippa_steps() -> None:
    from gippa.context import oriter

    assert oriter.LOGGER.name.startswith("gippa.steps.")
    assert oriter.LOGGER.name.rsplit(".", 1)[-1] == "context_stats"


def test_oriter_event_reaches_the_step_events_file(tmp_path: Path) -> None:
    """The handler ``gippa run`` (local) and ``gippa _step`` (each Nextflow task) install."""
    from gippa.cli import _record_step_events
    from gippa.context import oriter

    logs = tmp_path / "logs"
    restore = _record_step_events(logs)
    try:
        oriter.oric_position(_two_calls(), "chr1", 4_641_652)
    finally:
        restore()
    path = logs / "events.jsonl"
    assert path.is_file(), "the multiple_oric event was lost (logger outside gippa.steps)"
    evs = [json.loads(ln) for ln in path.read_text().splitlines() if ln.strip()]
    (ev,) = [e for e in evs if e["event"] == ORIC_EVENT]
    assert ev["step"] == "context_stats"
    assert ev["level"] == "INFO"
    assert sorted(ev["data"]["midpoints"]) == pytest.approx([1_000.0, 200_000.0], abs=1.0)
    assert ev["data"]["contig"] == "chr1"


def test_refused_oric_event_also_reaches_the_file(tmp_path: Path) -> None:
    from gippa.cli import _record_step_events
    from gippa.context import oriter

    logs = tmp_path / "logs"
    restore = _record_step_events(logs)
    try:
        assert oriter.oric_position(od.oric_at("chr1", 1_000, 300_000), "chr1", 4_641_652) is None
    finally:
        restore()
    evs = [json.loads(ln) for ln in (logs / "events.jsonl").read_text().splitlines()]
    (ev,) = [e for e in evs if e["event"] == ORIC_EVENT]
    assert ev["data"]["accepted"] is False and ev["data"]["oric"] is None


def test_oriter_event_from_file_becomes_context_oric(tmp_path: Path) -> None:
    from gippa.cli import _record_step_events
    from gippa.context import oriter

    logs = tmp_path / "logs"
    restore = _record_step_events(logs)
    try:
        oriter.oric_position(_two_calls(), "chr1", 4_641_652)
    finally:
        restore()
    m = build(events=_ms().read_events(logs / "events.jsonl"))
    (o,) = m["context"]["oric"]
    assert o["contig"] == "chr1" and o["accepted"] is True
    assert o["oric"] == pytest.approx(100_500.0, abs=1.0)


def test_one_oric_call_logs_nothing(caplog: pytest.LogCaptureFixture) -> None:
    from gippa.context import oriter

    with caplog.at_level(logging.INFO, logger="gippa"):
        oriter.oric_position(od.oric_at("chr1", 5_000), "chr1", 4_641_652)
    assert not [r for r in caplog.records if getattr(r, "event", None) == ORIC_EVENT]
