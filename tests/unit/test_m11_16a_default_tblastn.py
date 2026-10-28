"""Lead pin (M11-16a): a run at table 4 or 25 with the shipped ``tblastn_args`` passes preflight.

Until 16d drops ``-db_gencode 11`` from the default (D1, Khoa 2026-10-28), the shipped default is
exempt from ``config.tool_genetic_code``; after 16d it carries no ``-db_gencode`` at all. Either
way a default run must not be an input error.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from unit.test_m11_16a_preflight import SMALL, TOOL_RULE, TT_RULE, _only, _run

from gippa.exitcodes import ExitCode


@pytest.mark.parametrize("table", [4, 11, 25])
def test_default_tblastn_args_pass_at_every_table(tmp_path: Path, table: int) -> None:
    report, _ = _run(tmp_path, SMALL, f"annotation.translation_table={table}")
    assert _only(report, TT_RULE).level == "ok"
    assert _only(report, TOOL_RULE).level == "ok"
    assert report.exit_code == int(ExitCode.SUCCESS)


def test_a_user_db_gencode_11_at_table_4_is_still_an_error(tmp_path: Path) -> None:
    report, _ = _run(
        tmp_path,
        SMALL,
        "annotation.translation_table=4",
        "mapping.proteins.tblastn_args=-db_gencode 11 -evalue 1e-3",
    )
    assert _only(report, TOOL_RULE).level == "error"
