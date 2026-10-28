"""``gippa.io.plot_files.write_instance``: SVG ids do not depend on float layout (R7 byte
identity; M9-01 behaviour 7; M9-26 work item 4).

matplotlib's SVG backend names clip paths ``p`` + 10 hex digits of a sha256 over the
full-precision clip rectangle (``backend_svg._get_clip_attrs`` -> ``_make_id``), and markers and
hatches ``m``/``h`` + 10 hex digits likewise. A 1-ulp layout difference (seen on macOS arm64 CI for
V18_stat_diagnostics.svg) then changes only the ids, while PDF and PNG stay identical.

Assumed (the planned fix): for SVG, ``write_instance`` renames every hashed id
(``id="(p|m|h)[0-9a-f]{10}"``) and each ``url(#...)`` / ``xlink:href="#..."`` reference to it to
sequential ids in order of first appearance, so the id names no longer carry float bits; the
drawn coordinates themselves are rounded by matplotlib's ``short_float_fmt`` and stay equal.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gippa.io.plot_files import write_instance
from gippa.viz import style
from gippa.viz.base import Instance

RECT = (0.125, 0.11, 0.775, 0.77)
HASHED_ID = re.compile(r'id="([pmh][0-9a-f]{10})"')
ID = re.compile(r'\bid="([^"]+)"')
URL_REF = re.compile(r"url\(#([^)]+)\)")
HREF_REF = re.compile(r'xlink:href="#([^"]+)"')


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _instance(nudge: bool) -> Instance:
    """A figure with clip paths (bars, a line), markers and a hatch; ``nudge`` moves the axes'
    left edge by one ulp (``np.nextafter``)."""
    fig = style.new_figure("single", 2.0)
    x0, y0, w, h = RECT
    ax = fig.add_axes((float(np.nextafter(x0, 1.0)) if nudge else x0, y0, w, h))
    ax.bar(np.arange(5), [3, 1, 4, 1, 5], hatch="//")
    ax.plot([0, 2, 4], [1.0, 1.0 / 3.0, 2.0], marker="o", label="trend")
    ax.scatter([1, 3], [2.0, 4.0], marker="s")
    ax.set_title("Fake: ulp")
    data = pd.DataFrame({"start": pd.array([0, 10], dtype="int64"), "value": [1.0, 2.0]})
    return Instance(suffix=None, figure=fig, data=data, title="Fake: ulp", coords=("start",))


def _svg(tmp: Path, nudge: bool) -> bytes:
    paths = write_instance(_instance(nudge), "V90", "fake", tmp, ["svg"], 100)
    (svg,) = [p for p in paths if p.suffix == ".svg"]
    return svg.read_bytes()


def test_one_ulp_axes_nudge_gives_byte_identical_svg(tmp_path: Path) -> None:
    a = _svg(tmp_path / "a", nudge=False)
    b = _svg(tmp_path / "b", nudge=True)
    ids_a = sorted(set(HASHED_ID.findall(a.decode())))
    ids_b = sorted(set(HASHED_ID.findall(b.decode())))
    assert a == b, f"SVG differs after a 1-ulp axes nudge; hashed ids a={ids_a} b={ids_b}"


@pytest.mark.parametrize("nudge", [False, True])
def test_every_svg_reference_resolves_to_an_id(tmp_path: Path, nudge: bool) -> None:
    text = _svg(tmp_path, nudge).decode("utf-8")
    ids = ID.findall(text)
    assert len(ids) == len(set(ids)), "duplicate SVG ids"
    refs = URL_REF.findall(text) + HREF_REF.findall(text)
    assert refs, "the figure should carry clip-path and marker references"
    assert sorted(set(refs) - set(ids)) == []


def test_no_hashed_svg_id_remains(tmp_path: Path) -> None:
    text = _svg(tmp_path, nudge=False).decode("utf-8")
    assert HASHED_ID.findall(text) == []
    assert re.findall(r"#[pmh][0-9a-f]{10}\b", text) == []


def test_ids_are_stable_across_two_identical_renders(tmp_path: Path) -> None:
    assert _svg(tmp_path / "a", nudge=False) == _svg(tmp_path / "b", nudge=False)
