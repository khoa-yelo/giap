# ruff: noqa: RUF001  (the en dash is an accepted NA marker)
"""K-12 plot survey, batch 1 items 2 and 3: V12 panels (a) and (b).

(a) stacked per-flag counts, so a site with 3 flags counted 3 times and the bars rose above the
site total. Required: one "sites carrying this flag" bar per (call, confidence, flag), not stacked
(each starts at 0) and not hidden behind another (bars of one group do not overlap: dodged); the
printed total equals the number of sites; no bar exceeds it; every flag present gets a bar.

(b) printed Python ``None`` ("None aa None") in the header of a site with no gene. Required: an
NA-safe header (the plots' "NA" convention, ``base.title_pct``; "intergenic" or a dash also
accepted), never the text "None".
"""

from __future__ import annotations

import itertools
from pathlib import Path

import matplotlib.colors as mcolors
import pytest
from helpers import viz_variants as vv
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _render(out: Path) -> Instance:
    import giap.viz.variants  # noqa: F401  (registers V12)

    result = base.registered()["V12"].fn(RunView.from_outdir(out), resolve_config())
    assert not isinstance(result, Skip), result
    return result[0]


def _ax(inst: Instance, label: str) -> Axes:
    (ax,) = [a for a in inst.figure.axes if a.get_label() == label]
    return ax


def _flag_bars(ax: Axes) -> list[Rectangle]:
    from giap.viz.variants import FLAG_COLOUR, NO_FLAG_COLOUR

    colours = {c.upper() for c in [*FLAG_COLOUR.values(), NO_FLAG_COLOUR]}
    return [p for p in ax.patches if isinstance(p, Rectangle) and p.get_height() > 0
            and mcolors.to_hex(p.get_facecolor()[:3]).upper() in colours]  # fmt: skip


#: ALT high: 2 sites, one carrying 3 flags; mixed candidate: 1 site without a flag
SITES = [
    vv.variant("v1", aa_pos=20, flags="snv;start_codon;database_derived"),
    vv.variant("v2", aa_pos=40, flags="snv"),
    vv.variant("v3", aa_pos=60, call="mixed", confidence="candidate", flags=None),
]


def test_a_flag_bars_start_at_zero(tmp_path: Path) -> None:
    bars = _flag_bars(_ax(_render(vv.build_out(tmp_path, SITES, [])), "a"))
    assert bars
    assert [round(b.get_y(), 9) for b in bars] == [0.0] * len(bars)


def test_a_no_bar_rises_above_the_site_total(tmp_path: Path) -> None:
    bars = _flag_bars(_ax(_render(vv.build_out(tmp_path, SITES, [])), "a"))
    assert max(b.get_y() + b.get_height() for b in bars) == pytest.approx(2.0)


def test_a_one_bar_per_flag_with_its_site_count(tmp_path: Path) -> None:
    from giap.viz.variants import FLAG_COLOUR, NO_FLAG_COLOUR

    bars = _flag_bars(_ax(_render(vv.build_out(tmp_path, SITES, [])), "a"))
    got = sorted((mcolors.to_hex(b.get_facecolor()[:3]).upper(), round(b.get_height(), 9))
                 for b in bars)  # fmt: skip
    want = sorted([(FLAG_COLOUR["snv"].upper(), 2.0), (FLAG_COLOUR["start_codon"].upper(), 1.0),
                   (FLAG_COLOUR["database_derived"].upper(), 1.0),
                   (NO_FLAG_COLOUR.upper(), 1.0)])  # fmt: skip
    assert got == want


def test_a_bars_of_one_group_are_dodged_not_hidden(tmp_path: Path) -> None:
    bars = _flag_bars(_ax(_render(vv.build_out(tmp_path, SITES, [])), "a"))
    spans = sorted((b.get_x(), b.get_x() + b.get_width()) for b in bars)
    overlaps = [(a, b) for a, b in itertools.pairwise(spans) if b[0] < a[1] - 1e-9]
    assert overlaps == [], spans


def test_a_total_labels_equal_the_site_counts(tmp_path: Path) -> None:
    ax = _ax(_render(vv.build_out(tmp_path, SITES, [])), "a")
    bars = _flag_bars(ax)
    labels = [t for t in ax.texts if t.get_text().strip().isdigit()]
    by_value = sorted(int(t.get_text()) for t in labels)
    assert by_value == [1, 2], [t.get_text() for t in ax.texts]
    # each total sits over its group: above every bar whose x range holds the label's x
    for t in labels:
        x, y = t.get_position()
        under = [b for b in bars if b.get_x() - 0.5 <= x <= b.get_x() + b.get_width() + 0.5]
        assert under, (t.get_text(), x)
        assert y >= max(b.get_height() for b in under) - 1e-9


def test_a_every_drawn_flag_is_in_the_legend(tmp_path: Path) -> None:
    ax = _ax(_render(vv.build_out(tmp_path, SITES, [])), "a")
    leg = ax.get_legend() or (ax.figure.legends[0] if ax.figure.legends else None)
    assert leg is not None
    names = {t.get_text() for t in leg.get_texts()}
    assert {"snv", "start_codon", "database_derived", "no flag"} <= names, names


# --------------------------------------------------------------------- (b) no gene


def _no_gene(vid: str) -> dict[str, object]:
    row = vv.variant(vid, aa_pos=30, codon_start=7_001)
    row.update(gene=None, aa_pos=None)
    return row


def test_b_header_of_a_site_without_a_gene_has_no_python_none(tmp_path: Path) -> None:
    hit = vv.plus_hit("MKTAYIAK", 6_995, 1)
    inst = _render(vv.build_out(tmp_path, [vv.variant("VAR0001"), _no_gene("VAR0002")], [hit]))
    texts = [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]
    assert not [t for t in texts if "None" in t], [t for t in texts if "None" in t]
    joined = "".join(t.get_text() for t in _ax(inst, "b").texts)
    assert any(m in joined for m in ("NA aa NA", "intergenic", "–", "NA")), joined


# --------------------------------------------------------------------- #322 review N2


def test_a_all_label_sits_on_the_tallest_bar_of_its_group(tmp_path: Path) -> None:
    """Dodged bars: the total is printed at the tallest bar (2), not the stacked sum (4)."""
    ax = _ax(_render(vv.build_out(tmp_path, SITES, [])), "a")
    bars = _flag_bars(ax)
    for t in [t for t in ax.texts if t.get_text().strip().isdigit()]:
        x, y = t.get_position()
        group = [b for b in bars if abs(b.get_x() + b.get_width() / 2 - x) <= 0.4 + 1e-9]
        assert group, (t.get_text(), x)
        assert y == pytest.approx(max(b.get_height() for b in group)), (t.get_text(), y)


def test_b_gene_present_with_aa_pos_na_is_na_safe(tmp_path: Path) -> None:
    """A site whose gene is known but whose ``aa_pos`` is NA: no "None", the position reads NA,
    and the site is listed in the data."""
    row = vv.variant("VAR0002", aa_pos=30, codon_start=7_001)
    row.update(aa_pos=None)
    inst = _render(vv.build_out(tmp_path, [vv.variant("VAR0001"), row], []))
    texts = [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]
    assert not [t for t in texts if "None" in t], [t for t in texts if "None" in t]
    b = inst.data[inst.data["panel"] == "b"]
    assert "VAR0002" in set(b["variant_id"])
    joined = "".join(t.get_text() for t in _ax(inst, "b").texts)
    assert "NA" in joined, joined
