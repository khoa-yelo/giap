"""M4-16 ``write_detection_summary`` and the ``detection_summary`` table schema (§9.2, §9.3, §9.5).

Resolved ambiguities (pinned here):
- The schema file is ``src/gippa/schemas/tables/detection_summary.json`` (spec name); the M3-16
  naming ``detection_summary.schema.json`` is accepted too (exactly one must exist).
- Missing values are written as ``NA``; ``flags`` is ``;``-joined (empty when none).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_summary_fixture as fx
from detection_summary_fixture import gene, qc_hit, qc_peptide

TABLE_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"


def schema_path() -> Path:
    found = [p for p in (TABLE_DIR / "detection_summary.json",
                         TABLE_DIR / "detection_summary.schema.json") if p.exists()]  # fmt: skip
    assert len(found) == 1, f"expected one detection_summary schema in {TABLE_DIR}, got {found}"
    return found[0]


def _descriptor() -> dict[str, Any]:
    d = json.loads(schema_path().read_text(encoding="utf-8"))
    return {k: v for k, v in d.items() if k != "coordinates"}


def _table() -> pd.DataFrame:
    from gippa.context import summary
    from gippa.rng import for_step

    rows = [gene(f"G{i}", f"T{i // 2}", i % 3 != 0, 0.1 + 0.05 * i,
                 cog_category=["KT", None, "J"][i % 3]) for i in range(14)]  # fmt: skip
    basic = summary.basic_rows(fx.genes_frame(rows), 100, for_step("context_stats", 2))
    qc = summary.qc_rows(fx.peptides_frame([qc_peptide("AAAAAAAK")]),
                         fx.hits_frame([qc_hit("AAAAAAAK")]), None)  # fmt: skip
    det = fx.detectability_rows([i % 2 == 0 for i in range(30)], [0.5] * 30)
    return summary.detection_summary([basic, qc, det, fx.tier_rows(9, 9)])  # type: ignore[no-any-return]


def _write(tmp_path: Path, df: pd.DataFrame) -> Path:
    from gippa.io.summary_tsv import write_detection_summary

    path = tmp_path / "detection_summary.tsv"
    write_detection_summary(path, df)
    return path


def test_schema_fields_match_summary_columns() -> None:
    fields = _descriptor()["fields"]
    assert tuple(f["name"] for f in fields) == fx.SUMMARY_COLUMNS
    sec = next(f for f in fields if f["name"] == "section")
    assert sec["constraints"]["enum"] == list(fx.SECTIONS)
    assert "NA" in _descriptor()["missingValues"]


def test_schema_is_valid_frictionless() -> None:
    frictionless = pytest.importorskip("frictionless")
    frictionless.Schema.from_descriptor(_descriptor())


def test_written_table_header_and_na(tmp_path: Path) -> None:
    from gippa.io.tsv import read_tsv

    df = _table()
    path = _write(tmp_path, df)
    header, rows = read_tsv(path)
    assert tuple(header) == fx.SUMMARY_COLUMNS
    assert len(rows) == len(df)
    cells = [c for line in path.read_text(encoding="utf-8").splitlines() for c in line.split("\t")]
    assert not {c.lower() for c in cells} & {"nan", "<na>"}  # "none" is a COG group label
    # contaminant counts are NA without C13
    qc = [r for r in rows if r[0] == "qc" and r[3] == "contaminant_peptides"]
    assert len(qc) == 1 and qc[0][4] is None
    sections = [r[0] for r in rows]
    assert sections == sorted(sections, key=fx.SECTIONS.index)


def test_written_values_round_trip(tmp_path: Path) -> None:
    from gippa.io.tsv import read_tsv

    df = _table()
    _, rows = read_tsv(_write(tmp_path, df))
    for r, (_, want) in zip(rows, df.iterrows(), strict=True):
        assert r[:4] == [str(want[c]) for c in ("section", "group_type", "group", "metric")]
        for j, c in ((4, "value"), (5, "ci_low"), (6, "ci_high")):
            if pd.isna(want[c]):
                assert r[j] is None, (c, r)
            else:
                assert math.isclose(float(r[j] or "nan"), float(want[c]), rel_tol=1e-9), (c, r)
        if pd.isna(want["n"]):
            assert r[7] is None
        else:
            assert int(r[7] or -1) == int(want["n"])
        assert fx.flag_set(r[8]) == fx.flag_set(want["flags"])


def test_written_table_validates_with_frictionless(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    path = _write(tmp_path, _table())
    resource = frictionless.Resource(
        path=path.name, basepath=str(tmp_path),
        schema=frictionless.Schema.from_descriptor(_descriptor()),
    )  # fmt: skip
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])
