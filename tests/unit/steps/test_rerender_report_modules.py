"""``rerender_report(outdir, cfg)`` takes the modules off from ``cfg.modules`` (#323 review S1):
the page's skip reasons follow the config it is given, also when ``OUT/`` has no
``resolved_config.yaml`` (``RunView.with_modules(cfg.modules)``, as the in-graph step's view).

V14 (``enrichment``) is registered with ``module="function"``; with ``function`` off it skips as
``module_off:function``, shown as "not run: the function module is off".
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers.report_content import write_manifest
from helpers.report_fixtures import cfg, preload, run_out

OFF = "modules=[core,mobile,bgc]"
OFF_WORDS = "the function module is off"


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))


def _skipped(out: Path, c: Any, monkeypatch: pytest.MonkeyPatch) -> dict[str, str]:
    """The ``skipped`` mapping ``rerender_report`` hands to ``build_report``."""
    import giap.report.build as B
    from giap.steps.plots import rerender_report

    seen: list[dict[str, str]] = []
    real = B.build_report

    def spy(view: Any, plots: Any, skipped: Any, c: Any, *rest: Any, **kw: Any) -> str:
        seen.append(dict(skipped))
        return real(view, plots, skipped, c, *rest, **kw)

    monkeypatch.setattr(B, "build_report", spy)
    rerender_report(out, c)
    assert len(seen) == 1
    return seen[0]


def _out(tmp_path: Path, resolved_modules: list[str] | None) -> Path:
    preload()
    out = run_out(tmp_path)
    write_manifest(out)
    if resolved_modules is not None:
        body = "modules:\n" + "".join(f"  - {m}\n" for m in resolved_modules)
        (out / "resolved_config.yaml").write_text(body, encoding="utf-8")
    return out


def test_no_resolved_config_and_function_off_gives_v14_module_off(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _out(tmp_path, None)
    assert not (out / "resolved_config.yaml").exists()
    c = cfg(OFF)
    assert "function" not in c.modules
    skipped = _skipped(out, c, monkeypatch)
    assert skipped.get("V14") == "module_off:function", skipped
    page = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert OFF_WORDS in " ".join(page.split())


def test_cfg_modules_win_over_resolved_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _out(tmp_path, ["core", "mobile", "bgc", "function"])
    skipped = _skipped(out, cfg(OFF), monkeypatch)
    assert skipped.get("V14") == "module_off:function", skipped


def test_function_on_in_cfg_is_not_module_off(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _out(tmp_path, ["core", "mobile", "bgc"])
    c = cfg()
    assert "function" in c.modules
    skipped = _skipped(out, c, monkeypatch)
    assert skipped.get("V14") != "module_off:function", skipped
