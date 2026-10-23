"""The report size budget (M9-14 behaviour 5; degenerate row "a 30-MB report (huge run): written;
warning event; exit 3 only with ``--strict``"; lead decision 2 "the 25 MB budget is measured with
both" the static SVGs and the twins; GIAP_SPEC.md §10.3 "The size target is < 25 MB").

The oversized pages are made of a hand-written ``plots/`` directory (a padded SVG and a data TSV)
and a test-only twin padding its layout, so no huge run is needed. Sizes are chosen far from the
25 MB line on both readings (25 000 000 and 25 * 2**20 bytes), so the tests do not pin which.
Never touches the real ``$GIAP_DB_DIR``.

Interface assumed beyond the spec text (M9-14 test author):
- the event is logged by the ``report`` step (``ctx.log`` or a module logger) as a record with
  ``event == "report.size_exceeded"``, level WARNING, ``data == {"bytes": <size of report.html>}``;
- ``giap report OUT`` gains ``--strict``: an oversized report is still written and the command
  exits 3 (``WARNINGS_STRICT``) with ``--strict``, 0 without.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import ModuleType

import pytest
from click.testing import CliRunner, Result
from helpers.interactive_fixtures import BANNER, register_twin, report_step, swap_twins
from helpers.report_fixtures import run_out, swap_registry

from giap.exitcodes import ExitCode

MB = 1_000_000
EVENT = "report.size_exceeded"


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


@pytest.fixture
def twins(monkeypatch: pytest.MonkeyPatch, fake_registry: ModuleType) -> ModuleType:
    return swap_twins(monkeypatch)


def _plots(root: Path, svg_bytes: int) -> Path:
    """``plots/V17_qc.svg`` padded to about ``svg_bytes`` and its data TSV."""
    plots = root / "plots"
    (plots / "data").mkdir(parents=True, exist_ok=True)
    pad = "p" * svg_bytes
    (plots / "V17_qc.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
        f'<desc id="pad">{pad}</desc><rect id="r" width="1" height="1"/></svg>\n',
        encoding="utf-8",
    )
    (plots / "data" / "V17_qc.tsv").write_text("x\ty\n0\t0\n1\t1\n", encoding="utf-8")
    return plots


def _events(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [r for r in caplog.records if getattr(r, "event", None) == EVENT]


def _cli(out: Path, *args: str) -> Result:
    from giap.cli import cli

    fast = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]
    res = CliRunner().invoke(cli, ["report", str(out), *fast, *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def test_small_report_logs_no_size_event(
    twins: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    out = run_out(tmp_path)
    register_twin("V17")
    with caplog.at_level(logging.INFO):
        html = report_step(tmp_path, plots=_plots(tmp_path, 1000), out=out)
    assert BANNER in html.read_text(encoding="utf-8")
    assert _events(caplog) == []


def test_oversized_report_is_written_with_a_warning_event(
    fake_registry: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    out = run_out(tmp_path)
    with caplog.at_level(logging.INFO):
        html = report_step(tmp_path, *["report.interactive=false"],
                           plots=_plots(tmp_path, 27 * MB), out=out)  # fmt: skip
    size = html.stat().st_size
    assert size > 27 * MB
    (rec,) = _events(caplog)
    assert rec.levelno == logging.WARNING
    assert dict(rec.data) == {"bytes": size}  # type: ignore[attr-defined]


def test_budget_counts_static_and_interactive_together(
    twins: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Lead decision 2: a 12 MB static SVG plus a 12 MB twin plus plotly.js (~4.8 MB) is over
    budget; the same page static-only is not."""
    out = run_out(tmp_path)
    register_twin("V17", meta_bytes=12 * MB)
    plots = _plots(tmp_path, 12 * MB)
    with caplog.at_level(logging.INFO):
        static = report_step(tmp_path, "report.interactive=false", plots=plots, out=out)
    assert static.stat().st_size < 13 * MB
    assert _events(caplog) == []
    caplog.clear()
    with caplog.at_level(logging.INFO):
        both = report_step(tmp_path, plots=plots, out=out)
    size = both.stat().st_size
    assert size > 28 * MB
    (rec,) = _events(caplog)
    assert dict(rec.data) == {"bytes": size}  # type: ignore[attr-defined]


def test_cli_oversized_exit_0_without_strict_3_with(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    _plots(out, 27 * MB)
    res = _cli(out)
    assert res.exit_code == ExitCode.SUCCESS, res.output[-2000:]
    page = out / "report" / "report.html"
    assert page.stat().st_size > 27 * MB
    page.unlink()
    res = _cli(out, "--strict")
    assert res.exit_code == ExitCode.WARNINGS_STRICT, res.output[-2000:]
    assert page.is_file() and page.stat().st_size > 27 * MB, "the page must still be written"


def test_cli_strict_small_report_exits_0(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _plots(out, 1000)
    res = _cli(out, "--strict")
    assert res.exit_code == ExitCode.SUCCESS, res.output[-2000:]
