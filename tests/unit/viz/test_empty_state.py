"""``empty_state`` and ``title_pct`` (M9-01 behaviour 9, index decision 4(b)).

Assumed: ``empty_state``'s ``width`` takes the ``new_figure`` literals ``"single"``/``"double"``;
``title_pct`` rounds to a whole percent.
"""

from __future__ import annotations

import matplotlib.text
import pytest

from giap.viz import style
from giap.viz.base import Instance, empty_state, title_pct

MESSAGE = "No islands at FWER ≤ 0.05 in this run"
COLUMNS = ["island_id", "contig", "start", "end", "q_value"]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def test_empty_state_returns_header_only_data_and_the_message_as_title() -> None:
    inst = empty_state("V02", "abundance_landscape", MESSAGE, COLUMNS, "double")
    assert isinstance(inst, Instance)
    assert inst.suffix is None
    assert inst.title == MESSAGE
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0


def test_empty_state_draws_the_message_without_axes() -> None:
    inst = empty_state("V02", "abundance_landscape", MESSAGE, COLUMNS, "single")
    texts = [t.get_text() for t in inst.figure.findobj(matplotlib.text.Text)]
    assert MESSAGE in texts
    assert all(not ax.axison for ax in inst.figure.axes)
    assert round(float(inst.figure.get_size_inches()[0]), 6) == 3.5


def test_empty_state_double_width() -> None:
    inst = empty_state("V01", "circos", "No genes in this run", ["locus_tag"], "double")
    assert round(float(inst.figure.get_size_inches()[0]), 6) == 7.2


@pytest.mark.parametrize(
    ("n", "d", "text"), [(12, 100, "12 %"), (0, 7, "0 %"), (5, 5, "100 %"), (0, 0, "NA")]
)
def test_title_pct(n: int, d: int, text: str) -> None:
    assert title_pct(n, d) == text
