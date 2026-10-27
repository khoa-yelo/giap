"""``gippa report OUT`` on a local run (no ``resolved_config.yaml`` in ``OUT`` or ``OUT/work``)
takes the run's config, modules included, from ``run_manifest.json`` ``config``, so the page
it writes equals the ``gippa run`` tail's (M9-16b example 3; the nightly regression after
``rerender_report`` started honouring ``cfg.modules``). New guard test authorised by the lead.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.report_fixtures import run_out, swap_registry

from gippa.config import resolve_config
from gippa.exitcodes import ExitCode


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


def _cfg_seen(out: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    from gippa.cli import cli
    from gippa.steps import plots as P

    swap_registry(monkeypatch)
    seen: list[Any] = []
    monkeypatch.setattr(P, "render_plots", lambda *a, **k: None)
    monkeypatch.setattr(P, "rerender_report", lambda o, c: seen.append(c) or o / "r.html")
    res = CliRunner().invoke(cli, ["report", str(out)])
    assert res.exit_code == ExitCode.SUCCESS, res.output
    return seen[0]


def test_no_resolved_config_takes_the_manifest_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = run_out(tmp_path)
    for p in (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml"):
        assert not p.exists()
    run_cfg = resolve_config(overrides=["modules=[core]", "report.dpi=72"])
    manifest = {"config": run_cfg.model_dump(mode="json")}
    (out / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    cfg = _cfg_seen(out, monkeypatch)
    assert list(cfg.modules) == ["core"]
    assert cfg.report.dpi == 72


def test_invalid_manifest_config_falls_back_to_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = run_out(tmp_path)
    manifest = {"config": {"no_such_key": 1}}
    (out / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    assert _cfg_seen(out, monkeypatch) == resolve_config()
