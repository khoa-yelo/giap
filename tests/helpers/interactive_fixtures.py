"""Test-only Plotly twins, a local external-load checker and step drivers for the M9-14 tests
(the interactive report: Plotly.js inlined once, twins, sortable/filterable tables, size budget).

M9-14 lead decision 3: this branch registers no twin itself (twins are M9-15), so every twin here
is registered into a swapped ``gippa.viz.interactive._REGISTRY`` (``swap_twins``), the way
``helpers.report_fixtures.swap_registry`` swaps the plot registry.

``external_loads`` is a local stand-in for M9-13's ``gippa.report.checks.external_loads`` (not on
main when these tests were written); M9-13 adds the shared checker, and these tests may switch to
it once it lands.

Interface pinned here (M9-14 test author, beyond the spec text):
- ``gippa.viz.interactive`` keeps its twins in a module-level dict ``_REGISTRY`` (plot ID ->
  ``InteractiveFn``), read at call time, so a test can swap it with ``monkeypatch``;
- ``gippa.steps.plots.write_report(view, plots_dir, report_dir, cfg)`` builds the twins from the
  data TSVs under ``plots_dir/data/`` (the step and ``gippa report`` both go through it).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

from helpers.report_fixtures import Node, Page, cfg, parse, preload, run_out

#: the banner at the top of ``plotly.offline.get_plotlyjs()`` (plotly 7.1.0)
BANNER = "plotly.js v4.1.1"

EXTERNAL = re.compile(r"^\s*(https?:|//|file:|ftp:)", re.I)
_CSS_URL = re.compile(r"url\(\s*['\"]?([^)'\"\s]+)")
_IMPORT = re.compile(r"@import\s+(?:url\()?\s*['\"]?([^)'\";\s]+)", re.I)
#: (tag, attribute) pairs whose value the browser fetches when the page is viewed
_LOAD_ATTRS = {
    ("script", "src"),
    ("img", "src"),
    ("img", "srcset"),
    ("iframe", "src"),
    ("embed", "src"),
    ("object", "data"),
    ("source", "src"),
    ("source", "srcset"),
    ("audio", "src"),
    ("video", "src"),
    ("video", "poster"),
    ("input", "src"),
    ("link", "href"),
    ("image", "href"),
    ("image", "xlink:href"),
    ("use", "href"),
    ("use", "xlink:href"),
    ("feimage", "href"),
    ("feimage", "xlink:href"),
}


def external_loads(html: str) -> list[str]:
    """Every resource the page would fetch from outside itself at view time: any
    ``<script src>`` (scripts must be inline), any ``<link>`` href, and every external URL in a
    loading attribute, a ``style`` attribute or a ``<style>`` element (``url(...)``,
    ``@import``). Text inside ``<script>`` (plotly.js's own URL strings) is not a load
    (M9-14 example 3)."""
    page = parse(html)
    found: list[str] = []
    for n in page.root.iter():
        for tag, attr in sorted(_LOAD_ATTRS):
            if n.tag != tag or attr not in n.attrs:
                continue
            if tag in ("script", "link") or EXTERNAL.match(n.attrs[attr]):
                found.append(n.attrs[attr])
        if "style" in n.attrs:
            found += [u for u in _CSS_URL.findall(n.attrs["style"]) if EXTERNAL.match(u)]
    for n in page.root.by_tag("style"):
        text = n.text()
        found += [u for u in _CSS_URL.findall(text) if EXTERNAL.match(u)]
        found += _IMPORT.findall(text)
    return found


# ---------------------------------------------------------------- twins


def swap_twins(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Give ``gippa.viz.interactive`` an empty twin registry for one test."""
    preload()
    import gippa.viz.interactive as interactive

    monkeypatch.setattr(interactive, "_REGISTRY", {})
    return interactive


def scatter_twin(data: pd.DataFrame, title: str, *, meta_bytes: int = 0) -> Any:
    """A Plotly scatter of the data TSV's ``x`` and ``y`` columns (plain lists, so the JSON holds
    numbers, not typed arrays); ``meta_bytes`` pads the layout with a string of that length."""
    import plotly.graph_objects as go

    fig = go.Figure(
        go.Scatter(x=data["x"].tolist(), y=data["y"].tolist(), mode="lines+markers"),
        layout={"title": {"text": title}},
    )
    if meta_bytes:
        fig.update_layout(meta="m" * meta_bytes)
    return fig


def register_twin(
    pid: str,
    *,
    calls: list[tuple[str, list[str], str]] | None = None,
    raises: bool = False,
    meta_bytes: int = 0,
) -> Callable[..., Any]:
    """Register a test-only twin for ``pid`` through ``gippa.viz.interactive.interactive``;
    ``calls`` collects ``(pid, columns, title)`` per call."""
    from gippa.viz.interactive import interactive

    @interactive(pid)
    def _twin(data: pd.DataFrame, title: str) -> Any:
        if calls is not None:
            calls.append((pid, [str(c) for c in data.columns], title))
        if raises:
            raise RuntimeError(f"twin {pid} failed")
        return scatter_twin(data, title, meta_bytes=meta_bytes)

    return _twin


# ---------------------------------------------------------------- driving the steps


def tables(out: Path) -> dict[str, Path]:
    """Every published file of ``out`` that ``RunView`` knows, by key (the steps' inputs)."""
    from gippa.io.run_view import PATHS
    from gippa.steps.plots import ALL_PUBLISHED_KEYS

    return {k: out / PATHS[k] for k in ALL_PUBLISHED_KEYS
            if k in PATHS and (out / PATHS[k]).is_file()}  # fmt: skip


def run_steps(tmp_path: Path, *overrides: str, out: Path | None = None) -> Path:
    """The ``plots`` then ``report`` steps on a small ``OUT/``; returns ``report.html``."""
    from gippa.steps import registry as reg

    out = out or run_out(tmp_path)
    work = tmp_path / "work"
    pre = out / "preflight" / "preflight.json"
    c = cfg(*overrides)
    reg.run_step("plots", c, {"preflight": pre, **tables(out)}, {"plots": work / "plots"})
    return report_step(tmp_path, *overrides, plots=work / "plots", out=out)


def report_step(tmp_path: Path, *overrides: str, plots: Path, out: Path) -> Path:
    """Only the ``report`` step, on an existing ``plots`` directory."""
    from gippa.steps import registry as reg

    pre = out / "preflight" / "preflight.json"
    report = tmp_path / "work" / "report"
    inputs = {"plots": plots, "preflight": pre, **tables(out)}
    reg.run_step("report", cfg(*overrides), inputs, {"report": report})
    html = report / "report.html"
    assert html.is_file()
    return html


def page_of(path: Path) -> Page:
    return parse(path.read_text(encoding="utf-8"))


def head(page: Page) -> Node:
    found = page.root.by_tag("head")
    assert len(found) == 1, f"{len(found)} <head> elements"
    return found[0]


def twin_divs(node: Node) -> list[Node]:
    """``<div id="plotly-...">`` elements under ``node``."""
    return [n for n in node.by_tag("div") if n.attrs.get("id", "").startswith("plotly-")]


def plotly_scripts(page: Page) -> list[Node]:
    """``<script>`` elements whose text holds the plotly.js banner."""
    return [s for s in page.root.by_tag("script") if BANNER in s.text()]
