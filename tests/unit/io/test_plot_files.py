"""``giap.io.plot_files``: image and data-TSV writers (M9-01 behaviour 7, example 9, index
decision 6). Render and structure checks only: no ``--mpl`` baselines (Khoa approves them at
M9-25).

Assumed: ``write_instance`` returns every path it wrote (images and the data TSV);
``remove_plot_files(plots_dir, pid)`` removes that ID's images and data TSVs, all instances, and
nothing else.
"""

from __future__ import annotations

import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from giap.io.plot_files import remove_plot_files, write_instance
from giap.viz import style
from giap.viz.base import Instance

ALL = ["svg", "png", "pdf"]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _instance(suffix: str | None = None, title: str = "Fake: 2 of 3 detected") -> Instance:
    fig = style.new_figure("single", 2.0)
    ax = fig.add_subplot()
    ax.plot([0, 1, 2], [0.0, 1.0 / 3.0, 2.0], label="trend")
    ax.set_title(title)
    ax.set_xlabel("position (bp)")
    data = pd.DataFrame(
        {
            "start": pd.array([1, 2501], dtype="int64"),
            "value": [1.0 / 3.0, np.nan],
            "detected": [True, False],
            "genes": [["g1", "g2"], ["g3"]],
            "label": ["a", None],
        }
    )
    return Instance(suffix=suffix, figure=fig, data=data, title=title)


def _png_dpi(path: Path) -> float:
    blob = path.read_bytes()
    i = blob.index(b"pHYs")
    ppm_x, _ppm_y, unit = struct.unpack(">IIB", blob[i + 4 : i + 13])
    assert unit == 1  # metre
    return ppm_x * 0.0254


def test_names_and_returned_paths(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    paths = write_instance(_instance(), "V90", "fake", plots, ALL, 100)
    expected = {plots / f"V90_fake.{f}" for f in ALL} | {plots / "data" / "V90_fake.tsv"}
    assert set(paths) == expected
    assert all(p.is_file() and p.stat().st_size > 0 for p in expected)


def test_multi_instance_suffix(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance("PRO0001"), "V03", "locus", plots, ["svg"], 100)
    assert (plots / "V03_locus_PRO0001.svg").is_file()
    assert (plots / "data" / "V03_locus_PRO0001.tsv").is_file()


@pytest.mark.parametrize(
    ("suffix", "clean"),
    [
        ("NC_000913.3:prophage 1/2", "NC_000913.3_prophage_1_2"),
        ("a b", "a_b"),
        ("x/y", "x_y"),
        ("ok-1.2_x", "ok-1.2_x"),
    ],
)
def test_suffix_is_sanitised(tmp_path: Path, suffix: str, clean: str) -> None:
    plots = tmp_path / "plots"
    paths = write_instance(_instance(suffix), "V03", "locus", plots, ["svg"], 100)
    assert plots / f"V03_locus_{clean}.svg" in paths
    assert plots / "data" / f"V03_locus_{clean}.tsv" in paths
    assert all(p.parent in (plots, plots / "data") for p in paths)


def test_formats_are_honoured(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance(), "V90", "fake", plots, ["svg"], 100)
    assert sorted(p.name for p in plots.iterdir() if p.is_file()) == ["V90_fake.svg"]
    assert (plots / "data" / "V90_fake.tsv").is_file()


def test_png_written_at_the_configured_dpi(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance(), "V90", "fake", plots, ["png"], 150)
    assert _png_dpi(plots / "V90_fake.png") == pytest.approx(150, abs=0.5)


def test_data_tsv_encoding_follows_adr_0006(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance(), "V90", "fake", plots, ["svg"], 100)
    text = (plots / "data" / "V90_fake.tsv").read_bytes().decode("utf-8")
    assert text == (
        "start\tvalue\tdetected\tgenes\tlabel\n1\t0.333333\ttrue\tg1;g2\ta\n2501\tNA\tfalse\tg3\tNA\n"
    )


def test_header_only_data_tsv(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    inst = _instance()
    empty = Instance(suffix=None, figure=inst.figure, data=inst.data.iloc[0:0], title="none")
    write_instance(empty, "V90", "fake", plots, ["svg"], 100)
    assert (plots / "data" / "V90_fake.tsv").read_text() == "start\tvalue\tdetected\tgenes\tlabel\n"


def test_no_dates_or_software_in_the_images(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance(), "V90", "fake", plots, ALL, 100)
    svg = (plots / "V90_fake.svg").read_text(encoding="utf-8")
    assert "<dc:date>" not in svg
    assert "Fake: 2 of 3 detected" in svg  # svg.fonttype none: text stays text
    pdf = (plots / "V90_fake.pdf").read_bytes()
    assert b"/CreationDate" not in pdf
    assert b"/ModDate" not in pdf
    png = (plots / "V90_fake.png").read_bytes()
    assert b"Software" not in png


def test_example_9_same_instance_twice_is_byte_identical(tmp_path: Path) -> None:
    inst = _instance("PRO0001")
    a = write_instance(inst, "V03", "locus", tmp_path / "a" / "plots", ALL, 100)
    b = write_instance(inst, "V03", "locus", tmp_path / "b" / "plots", ALL, 100)
    assert len(a) == len(b) == 4
    for pa in a:
        pb = tmp_path / "b" / pa.relative_to(tmp_path / "a")
        assert pb in b
        assert pa.read_bytes() == pb.read_bytes(), pa.name


def test_two_renders_of_the_same_data_are_byte_identical(tmp_path: Path) -> None:
    a = write_instance(_instance(), "V90", "fake", tmp_path / "a", ALL, 100)
    b = write_instance(_instance(), "V90", "fake", tmp_path / "b", ALL, 100)
    for pa, pb in zip(sorted(a), sorted(b), strict=True):
        assert pa.name == pb.name
        assert pa.read_bytes() == pb.read_bytes(), pa.name


def test_figures_are_closed_after_writing(tmp_path: Path) -> None:
    inst = _instance()
    write_instance(inst, "V90", "fake", tmp_path / "plots", ["svg"], 100)
    assert getattr(inst.figure, "number", None) not in plt.get_fignums()


def test_remove_plot_files_removes_only_that_id(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    write_instance(_instance("A"), "V03", "locus", plots, ALL, 100)
    write_instance(_instance("B"), "V03", "locus", plots, ["svg"], 100)
    write_instance(_instance(), "V02", "abundance_landscape", plots, ["svg"], 100)
    support = plots / "support" / "gc_skew.tsv"
    support.parent.mkdir(parents=True)
    support.write_text("contig\n")
    remove_plot_files(plots, "V03")
    left = sorted(str(p.relative_to(plots)) for p in plots.rglob("*") if p.is_file())
    assert left == [
        "V02_abundance_landscape.svg",
        "data/V02_abundance_landscape.tsv",
        "support/gc_skew.tsv",
    ]


def test_remove_plot_files_on_a_missing_directory_is_a_no_op(tmp_path: Path) -> None:
    remove_plot_files(tmp_path / "plots", "V01")
