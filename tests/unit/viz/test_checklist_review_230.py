"""PR #230 review B1: the colour-blind item's values are pinned — the CIELAB conversion, the
``min ΔE before X.X → after Y.Y`` detail on a known palette (lead decision 8), and the written
``deuteranopia_*.png`` equal to the simulation of the original PNG."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_hex

from giap.io.plot_files import read_png, write_instance
from giap.viz import checklist as ck
from giap.viz import style
from giap.viz.base import Instance


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def test_lab_of_red_and_white() -> None:
    lab = ck._lab(np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 1.0]]))
    np.testing.assert_allclose(lab[0], [53.24, 80.09, 67.20], atol=0.02)
    np.testing.assert_allclose(lab[1], [100.0, 0.0, 0.0], atol=0.01)


def _written(tmp: Path) -> tuple[Instance, list[Path]]:
    fig = plt.figure(figsize=(style.SINGLE_IN, 2.0))
    ax = fig.add_subplot()
    colours = list(style.DETECTION.values())
    ax.bar(range(len(colours)), [1, 2, 3, 4][: len(colours)], color=colours)
    ax.set_xlabel("category")
    ax.set_ylabel("genes")
    inst = Instance(suffix=None, figure=fig, data=pd.DataFrame({"x": [1]}), title="3 genes")
    files = write_instance(inst, "V99", "fake", tmp / "plots", ["png"], 100)
    return inst, files


def test_detail_on_the_detection_palette(tmp_path: Path) -> None:
    inst, files = _written(tmp_path)
    colours = {to_hex(c).upper() for c in style.DETECTION.values()}
    assert ck._cvd(inst, files, "V99_fake", colours) == "min ΔE before 51.1 → after 55.3"


def test_simulated_png_is_the_simulation_of_the_original(tmp_path: Path) -> None:
    inst, files = _written(tmp_path)
    ck._cvd(inst, files, "V99_fake", {"#0072B2", "#E69F00"})
    png = next(p for p in files if p.suffix == ".png")
    orig, sim = read_png(png), read_png(png.parent / "deuteranopia_V99_fake.png")
    want = ck.simulate_deuteranopia(orig[..., :3])
    assert sim.shape == orig.shape
    assert np.abs(sim[..., :3] - want).max() <= 1 / 255 + 1e-9
    assert np.abs(sim[..., :3] - orig[..., :3]).max() > 0.05  # it was simulated
