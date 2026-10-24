"""M9-16b ``rerender_report(outdir, cfg)`` on its own (docs/specs/M9-16b-report-after-manifest.md
behaviour 3, examples 4 and 10, the degenerate rows "``OUT/plots/`` absent or empty" and
"``OUT/report/report.html`` absent", and example 8's temp-file check).

Interface pinned here (M9-16b test author, beyond the spec text):

* ``giap.steps.plots.rerender_report(outdir: Path, cfg: GiapConfig) -> Path`` returns
  ``outdir / "report" / "report.html"``;
* the page is byte-identical to what ``write_report(RunView.from_outdir(outdir), outdir /
  "plots" or None, <dir>, cfg)`` writes (``build_report`` on ``read_plot_files``, with the
  interactive twins once M9-14 is in the tree): the same page ``giap report`` builds;
* the atomic write goes through ``os.replace(report/.report.html.tmp, report/report.html)``
  (``Path.replace`` calls ``os.replace`` too), and the target still holds its old bytes when
  ``os.replace`` is called;
* an exception from ``build_report`` propagates (the CLI turns it into exit 4) and leaves no
  ``.report.html.tmp`` behind.

Fake plots (tests/helpers/report_fixtures.py) stand in for the real ones. Never touches the real
``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.report_content import sha256_file, write_manifest
from helpers.report_fixtures import cfg, register_fake, run_out, swap_registry

TMP_NAME = ".report.html.tmp"


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _rerender() -> Callable[[Path, Any], Path]:
    from giap.steps import plots as P

    fn = getattr(P, "rerender_report", None)
    assert fn is not None, "giap.steps.plots.rerender_report is not implemented (M9-16b)"
    return fn  # type: ignore[no-any-return]


def _tree(root: Path) -> dict[str, str]:
    """Relative path -> SHA-256 of every file under ``root`` (symlinks as ``link:<target>``)."""
    out: dict[str, str] = {}
    if not root.exists():
        return out
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root).as_posix()
        if p.is_symlink():
            out[rel] = f"link:{os.readlink(p)}"
        elif p.is_file():
            out[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def _rendered_out(tmp_path: Path) -> Path:
    """``OUT/`` with tables, a preflight, a manifest and ``plots/`` drawn by two fake plots."""
    from giap.io.run_view import RunView
    from giap.steps import plots as P

    register_fake("V07")
    register_fake("V17", suffixes=(None, "b"))
    out = run_out(tmp_path, islands=3)
    write_manifest(out)
    P.render_plots(RunView.from_outdir(out), cfg(), out / "plots", None)
    assert (out / "plots" / "V07_obs_vs_exp.svg").is_file()
    return out


def _expected_page(out: Path, c: Any, scratch: Path) -> str:
    """What ``write_report`` (the in-graph step's writer) gives on ``RunView.from_outdir(out)``
    and ``out/plots``: the same page ``giap report`` builds."""
    from giap.io.run_view import RunView
    from giap.steps import plots as P

    plots_dir = out / "plots" if (out / "plots").is_dir() else None
    page = P.write_report(RunView.from_outdir(out), plots_dir, scratch, c)
    return page.read_text(encoding="utf-8")


def _forbid_drawing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every way of drawing or writing under ``plots/`` raises."""
    from giap.steps import plots as P
    from giap.viz import base

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("rerender_report must not draw plots or write under plots/")

    monkeypatch.setattr(P, "render_plots", boom)
    monkeypatch.setattr(P, "write_plot_support", boom)
    monkeypatch.setattr(base, "_REGISTRY", {})
    register_fake("V07", raises=True)
    register_fake("V17", suffixes=(None, "b"), raises=True)


# ---------------------------------------------------------------- behaviour 3


def test_returns_the_report_path_and_writes_the_giap_report_page(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _rendered_out(tmp_path)
    c = cfg()
    path = _rerender()(out, c)
    assert Path(path) == out / "report" / "report.html"
    assert path.read_text(encoding="utf-8") == _expected_page(out, c, tmp_path / "expected")


def test_the_page_reads_the_manifest_on_disk(fake_registry: ModuleType, tmp_path: Path) -> None:
    from helpers.report_fixtures import parse

    out = _rendered_out(tmp_path)
    page = parse(_rerender()(out, cfg()).read_text(encoding="utf-8"))
    text = " ".join(page.section("provenance").text().split())
    assert sha256_file(out / "run_manifest.json") in text, text
    assert "run_manifest.json not written by this GIAP version" not in text


def test_example_4_twice_identical_and_plots_untouched(
    fake_registry: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _rendered_out(tmp_path)
    before = _tree(out / "plots")
    assert before
    c = cfg()
    _forbid_drawing(monkeypatch)
    first = _rerender()(out, c).read_bytes()
    second = _rerender()(out, c).read_bytes()
    assert first == second
    assert _tree(out / "plots") == before
    assert b"V07_obs_vs_exp.svg" in first, "the existing plot files are embedded or linked"


def test_the_run_cfg_reaches_build_report(
    fake_registry: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import giap.report.build as B

    out = _rendered_out(tmp_path)
    seen: list[Any] = []
    real = B.build_report

    def spy(view: Any, plots: Any, skipped: Any, c: Any, *rest: Any, **kw: Any) -> str:
        seen.append(c)
        return real(view, plots, skipped, c, *rest, **kw)

    monkeypatch.setattr(B, "build_report", spy)
    c = cfg("report.dpi=51")
    _rerender()(out, c)
    assert len(seen) == 1 and seen[0] == c


# ---------------------------------------------------------------- degenerate rows


def test_plots_dir_absent(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    out = run_out(tmp_path, islands=3)
    write_manifest(out)
    assert not (out / "plots").exists()
    c = cfg()
    path = _rerender()(out, c)
    assert path.read_text(encoding="utf-8") == _expected_page(out, c, tmp_path / "expected")
    assert not (out / "plots").exists(), "rerender_report never creates plots/"


def test_plots_dir_empty(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = run_out(tmp_path)
    (out / "plots").mkdir()
    c = cfg()
    path = _rerender()(out, c)
    assert path.read_text(encoding="utf-8") == _expected_page(out, c, tmp_path / "expected")
    assert _tree(out / "plots") == {}


def test_report_dir_absent_is_written_fresh(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = _rendered_out(tmp_path)
    assert not (out / "report").exists()
    path = _rerender()(out, cfg())
    assert path.is_file()
    assert sorted(p.name for p in (out / "report").iterdir()) == ["report.html"]


# ---------------------------------------------------------------- atomic write


def test_atomic_replace_from_a_temp_file_in_the_same_directory(
    fake_registry: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _rendered_out(tmp_path)
    page = out / "report" / "report.html"
    page.parent.mkdir()
    page.write_bytes(b"<!-- in-graph page -->")
    calls: list[tuple[Path, Path, bytes]] = []
    real = os.replace

    def spy(src: Any, dst: Any, *a: Any, **k: Any) -> None:
        calls.append((Path(src), Path(dst), Path(dst).read_bytes() if Path(dst).exists() else b""))
        real(src, dst, *a, **k)

    monkeypatch.setattr(os, "replace", spy)
    _rerender()(out, cfg())
    mine = [c for c in calls if c[1] == page]
    assert len(mine) == 1, f"os.replace calls: {calls}"
    src, _dst, old = mine[0]
    assert src == out / "report" / TMP_NAME
    assert old == b"<!-- in-graph page -->", "the target is untouched until the replace"
    assert not src.exists()
    assert sorted(p.name for p in page.parent.iterdir()) == ["report.html"]


def test_build_report_raising_keeps_the_old_page_and_leaves_no_temp(
    fake_registry: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import giap.report.build as B
    from giap.steps import plots as P

    out = _rendered_out(tmp_path)
    page = out / "report" / "report.html"
    page.parent.mkdir()
    page.write_bytes(b"<!-- in-graph page -->")
    manifest = (out / "run_manifest.json").read_bytes()

    def boom(*_a: Any, **_k: Any) -> str:
        raise RuntimeError("build_report failed on purpose")

    rerender = _rerender()
    monkeypatch.setattr(B, "build_report", boom)
    if hasattr(P, "build_report"):  # a module-level import in plots.py
        monkeypatch.setattr(P, "build_report", boom)
    with pytest.raises(RuntimeError, match="on purpose"):
        rerender(out, cfg())
    assert page.read_bytes() == b"<!-- in-graph page -->"
    assert not (out / "report" / TMP_NAME).exists()
    assert (out / "run_manifest.json").read_bytes() == manifest


def test_example_10_a_symlinked_page_is_replaced_not_written_through(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = _rendered_out(tmp_path)
    scratch = tmp_path / "scratch" / "report.html"
    scratch.parent.mkdir()
    scratch.write_bytes(b"<!-- work copy -->")
    page = out / "report" / "report.html"
    page.parent.mkdir()
    page.symlink_to(scratch)
    _rerender()(out, cfg())
    assert page.is_file() and not page.is_symlink()
    assert scratch.read_bytes() == b"<!-- work copy -->"
    assert b"<html" in page.read_bytes().lower()
