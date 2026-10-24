"""K-12 plot survey, batch 1 item 9: the checklist tool (``giap.viz.checklist``).

- The neutral greys the plots draw as chart furniture (gene-arrow edges ``#333333``, V10's
  gene-label envelope ``#CCCCCC``, V01's sector and ring backgrounds ``#DDDDDD`` / ``#F4F4F4``)
  failed "palette as §10.1". They are allowed through an allow-list of neutral greys kept as
  ``giap.viz.style`` constants (one place for the palette), not inline in the checklist.
- "Moran's I" (V10) and "candidates" (V13) are unitless axis labels.

Interface assumed: the allow-list is some public ``giap.viz.style`` constant (any name) whose
values (a mapping's values, or the members of a tuple / frozenset / set) include the four greys.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure

from giap.io.plot_files import write_instance
from giap.viz import checklist as ck
from giap.viz import style
from giap.viz.base import Instance

GREYS = ("#333333", "#CCCCCC", "#DDDDDD", "#F4F4F4")
PALETTE, UNITS = "palette as §10.1", "axes carry units"
TITLE = "Prophage genes: 12 % detected vs 31 % expected"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _check(fig: Figure, tmp: Path) -> dict[str, ck.CheckResult]:
    data = pd.DataFrame({"category": ["a", "b"], "value": [1.0, 2.0]})
    inst = Instance(suffix=None, figure=fig, data=data, title=TITLE)
    files = write_instance(inst, "V90", "fake", tmp / "plots", ["svg", "png"], 100)
    return {r.item: r for r in ck.check_instance(inst, "V90", files)}


def _bars(colours: list[str], edge: str | None = None) -> Figure:
    fig = Figure(figsize=(style.SINGLE_IN, 2.5))
    ax = fig.add_subplot()
    kw = {} if edge is None else {"edgecolor": edge, "linewidth": 0.5}
    ax.bar(range(len(colours)), [1.0] * len(colours), color=colours, **kw)
    ax.set_xticks([])
    ax.set_ylabel("q")
    ax.set_title(TITLE)
    return fig


def _style_colours() -> set[str]:
    """Every hex string held by a public ``giap.viz.style`` constant."""
    found: set[str] = set()
    for name in dir(style):
        if name.startswith("_"):
            continue
        value = getattr(style, name)
        if isinstance(value, Mapping):
            items = list(value.values())
        elif isinstance(value, tuple | frozenset | set | list):
            items = list(value)
        else:
            continue
        for item in items:
            for v in item if isinstance(item, tuple) else (item,):
                if isinstance(v, str) and v.startswith("#"):
                    found.add(v.upper())
    return found


def test_the_neutral_greys_are_style_constants() -> None:
    missing = [g for g in GREYS if g not in _style_colours()]
    assert missing == [], f"{missing} not held by any giap.viz.style constant"


@pytest.mark.parametrize("grey", GREYS)
def test_a_neutral_grey_fill_passes_the_palette_item(tmp_path: Path, grey: str) -> None:
    res = _check(_bars([style.DETECTION["detected"], grey]), tmp_path)[PALETTE]
    assert res.passed is True, res.detail


def test_neutral_grey_edges_pass_the_palette_item(tmp_path: Path) -> None:
    """Gene arrows: a viridis / detection fill with a ``#333333`` edge (V03, V13)."""
    fig = _bars([style.DETECTION["detected"], style.DETECTION["not_detected"]], edge="#333333")
    res = _check(fig, tmp_path)[PALETTE]
    assert res.passed is True, res.detail


def test_all_four_greys_in_one_figure_pass(tmp_path: Path) -> None:
    res = _check(_bars(list(GREYS)), tmp_path)[PALETTE]
    assert res.passed is True, res.detail


def test_a_non_grey_off_palette_colour_still_fails(tmp_path: Path) -> None:
    """Guard: the allow-list is the neutral greys, not a blanket pass."""
    res = _check(_bars(["#333333", "#123456"]), tmp_path)[PALETTE]
    assert res.passed is False
    assert "#123456" in res.detail.upper()


@pytest.mark.parametrize("label", ["Moran's I", "candidates"])
def test_unitless_labels_from_the_survey(tmp_path: Path, label: str) -> None:
    fig = _bars([style.DETECTION["detected"]])
    fig.axes[0].set_xlabel(label)
    res = _check(fig, tmp_path)[UNITS]
    assert res.passed is True, res.detail


@pytest.mark.parametrize("label", ["Moran's I", "candidates"])
def test_unitless_entries_are_listed(label: str) -> None:
    assert label.lower() in {u.lower() for u in ck.UNITLESS}


def test_a_quantity_label_still_needs_a_unit(tmp_path: Path) -> None:
    """Guard: adding the two labels does not exempt every label."""
    fig = _bars([style.DETECTION["detected"]])
    fig.axes[0].set_xlabel("distance")
    assert _check(fig, tmp_path)[UNITS].passed is False
