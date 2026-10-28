"""PR #467 review: pins three gaps left by the M11-17 tests (surviving mutants).

1. A directory named ``<stem>.gbff`` is not a GBFF (decision Q3; ``is_file``, not ``exists``).
2. ``preflight.md``: the Message cell is never cut, the Value cell keeps its 200-character cap.
3. Only a final ``.json[.gz]`` is replaced; ``.json`` inside a directory name stays.
"""

from __future__ import annotations

import json
from pathlib import Path

from unit.preflight.test_report_review import check_rows, chk, report

from gippa.config import resolve_config
from gippa.preflight.gbff import bgc_gbff_check, expected_gbff
from gippa.preflight.report import render_preflight_md


def test_a_directory_named_gbff_does_not_count(tmp_path: Path) -> None:
    js = tmp_path / "x.json"
    js.write_text("{}", encoding="utf-8")
    (tmp_path / "x.gbff").mkdir()
    (tmp_path / "x.gbff.gz").mkdir()
    cfg = resolve_config(overrides=[f"annotation.precomputed.bakta_json={json.dumps(str(js))}"])
    c = bgc_gbff_check(cfg, None)
    assert c is not None
    assert c.level == "warn"


def test_message_cell_is_full_and_value_cell_is_capped() -> None:
    message = "m" * 500
    value = "v" * 500
    md = render_preflight_md(report([chk("inputs.long_rule", "warn", value, message)]))
    (row,) = (r for r in check_rows(md) if r[1] == "inputs.long_rule")
    assert row[3] == message
    assert len(row[2]) == 200
    assert row[2].endswith("…")


def test_json_inside_a_directory_name_is_kept() -> None:
    assert expected_gbff(Path("/d/a.json.d/x.json")) == (
        Path("/d/a.json.d/x.gbff"),
        Path("/d/a.json.d/x.gbff.gz"),
    )
    assert expected_gbff(Path("/d/a.json/x.JSON.gz"))[0] == Path("/d/a.json/x.gbff")
