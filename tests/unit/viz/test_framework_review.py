"""PR #191 review round (M9-01): ``binned`` integer edges (blocker B2) and two nits.

Pinned here (lead decisions on the review):
- ``binned`` returns integer edges (an integer NumPy dtype) of ``n = min(max_elements, length)``
  bins over ``[0, length)``; every bin is at least 1 bp, and bin widths differ by at most 1 bp.
  Written to a data TSV, the edges carry no decimal point or exponent (§9.2 ``%.6g`` floats lose
  precision on a 5 Mb replicon).
- ``title_pct`` rounds half up: ``title_pct(1, 8) == "13 %"``.
- ``style.COG_SUPER`` agrees with the bundled ``giap/data/cog_categories.tsv`` ``super_category``
  column for every letter.
"""

from __future__ import annotations

import re
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.io.plot_files import write_instance
from giap.viz import style
from giap.viz.base import Instance, binned, title_pct


def _intervals(n: int, length: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, length, size=n)
    ends = np.minimum(starts + rng.integers(1, 3000, size=n), length)
    return starts, ends


# ---------------------------------------------------------------- B2 binned


def test_binned_edges_are_integers() -> None:
    starts, ends = _intervals(25_000, 1_000_000, 3)
    edges = binned(starts, ends, 1_000_000)
    assert edges is not None
    assert np.issubdtype(edges.dtype, np.integer)
    assert len(edges) == 20_001
    assert edges[0] == 0 and edges[-1] == 1_000_000
    assert set(np.diff(edges).tolist()) == {50}


def test_binned_short_replicon_gives_one_bp_bins() -> None:
    starts, ends = _intervals(1_001, 100, 5)
    edges = binned(starts, ends, 100, max_elements=1_000)
    assert edges is not None
    assert np.issubdtype(edges.dtype, np.integer)
    assert edges.tolist() == list(range(101))  # n = min(1000, 100) = 100 bins of 1 bp


def test_binned_5_mb_replicon_data_tsv_has_plain_integer_edges(tmp_path: Path) -> None:
    style.apply()
    length = 4_999_999
    starts, ends = _intervals(25_001, length, 7)
    edges = binned(starts, ends, length)
    assert edges is not None
    data = pd.DataFrame({"bin_start": edges[:-1], "bin_end": edges[1:]})
    fig = style.new_figure("single", 2.0)
    fig.add_subplot().plot([0, 1], [0, 1])
    paths = write_instance(Instance(None, fig, data, "bins"), "V99", "bins", tmp_path, ["svg"], 72)
    tsv = next(p for p in paths if p.suffix == ".tsv")
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "bin_start\tbin_end"
    assert len(lines) == 20_001
    bad = [line for line in lines[1:] if not re.fullmatch(r"\d+\t\d+", line)]
    assert bad == [], bad[:3]
    assert lines[-1].split("\t")[1] == str(length)


@settings(max_examples=60, deadline=None)
@given(
    max_elements=st.integers(1, 5_000),
    extra=st.integers(1, 50),
    length=st.integers(1, 10_000_000),
    seed=st.integers(0, 2**32 - 1),
)
def test_property_binned_integer_edges_at_least_one_bp(
    max_elements: int, extra: int, length: int, seed: int
) -> None:
    starts, ends = _intervals(max_elements + extra, length, seed)
    edges = binned(starts, ends, length, max_elements=max_elements)
    assert edges is not None
    assert np.issubdtype(edges.dtype, np.integer)
    assert len(edges) - 1 == min(max_elements, length)
    assert edges[0] == 0 and edges[-1] == length
    widths = np.diff(edges)
    assert widths.min() >= 1
    assert widths.max() - widths.min() <= 1


# ---------------------------------------------------------------- nits


@pytest.mark.parametrize(("n", "d", "text"), [(1, 8, "13 %"), (1, 3, "33 %"), (2, 3, "67 %")])
def test_title_pct_rounds_half_up(n: int, d: int, text: str) -> None:
    assert title_pct(n, d) == text


def test_cog_super_agrees_with_bundled_table() -> None:
    text = resources.files("giap.data").joinpath("cog_categories.tsv").read_text("utf-8")
    rows = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
    header, body = rows[0], rows[1:]
    i_letter, i_super = header.index("letter"), header.index("super_category")
    table = {r[i_letter]: r[i_super] for r in body}
    assert len(table) == 26
    assert dict(style.COG_SUPER) == table
