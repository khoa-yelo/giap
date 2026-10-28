"""M9-16b behaviour 4: ``gippa report OUT`` keeps its plot rendering, then writes the page through
the same ``rerender_report(out, cfg)`` as the ``gippa run`` tail; ``write_report`` stays for the
in-graph ``report`` step only (docs/specs/M9-16b-report-after-manifest.md, "Public interface";
acceptance "tests/unit/test_cli_report.py: gippa report goes through rerender_report").

Interface pinned here (M9-16b test author): the command looks ``render_plots``,
``rerender_report`` and ``write_report`` up on ``gippa.steps.plots`` at call time (``from
gippa.steps import plots as P``, as today), and prints the returned path on stdout. Fake plots
(tests/helpers/report_fixtures.py); never touches the real ``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from click.testing import CliRunner, Result
from helpers.report_fixtures import register_fake, run_out, swap_registry

from gippa.exitcodes import ExitCode

FAST = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, tuple[Any, ...]]]:
    """Record ``render_plots``, ``rerender_report`` and ``write_report`` calls, in order."""
    from gippa.steps import plots as P

    seen: list[tuple[str, tuple[Any, ...]]] = []
    for name in ("render_plots", "rerender_report", "write_report"):
        real = getattr(P, name, None)

        def wrapper(*a: Any, _n: str = name, _r: Any = real, **k: Any) -> Any:
            seen.append((_n, a))
            assert _r is not None, f"gippa.steps.plots.{_n} does not exist"
            return _r(*a, **k)

        monkeypatch.setattr(P, name, wrapper, raising=False)
    return seen


def _report(out: Path, *args: str) -> Result:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["report", str(out), *FAST, *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def test_gippa_report_renders_plots_then_calls_rerender_report(
    fake_registry: ModuleType, tmp_path: Path, calls: list[tuple[str, tuple[Any, ...]]]
) -> None:
    register_fake("V07")
    out = run_out(tmp_path, islands=3)
    res = _report(out)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    names = [n for n, _ in calls]
    assert names == ["render_plots", "rerender_report"], names
    _name, args = calls[1]
    assert Path(args[0]) == out
    assert args[1].report.dpi == 50, "the command's resolved cfg"
    assert str(out / "report" / "report.html") in res.stdout


def test_no_interactive_reaches_rerender_report(
    fake_registry: ModuleType, tmp_path: Path, calls: list[tuple[str, tuple[Any, ...]]]
) -> None:
    out = run_out(tmp_path)
    res = _report(out, "--no-interactive")
    assert res.exit_code == ExitCode.SUCCESS, res.output
    rer = [a for n, a in calls if n == "rerender_report"]
    assert len(rer) == 1 and rer[0][1].report.interactive is False


def test_a_failing_rerender_is_exit_4(
    fake_registry: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.steps import plots as P

    def boom(*_a: Any, **_k: Any) -> Path:
        raise RuntimeError("re-render failed on purpose")

    monkeypatch.setattr(P, "rerender_report", boom, raising=False)
    monkeypatch.setattr(P, "write_report", boom)
    res = _report(run_out(tmp_path))
    assert res.exit_code == ExitCode.STEP_FAILED, res.output
    assert "re-render failed on purpose" in res.output
