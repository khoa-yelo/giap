"""K-12 plot survey, batch 1 item 10: the V10 join-count scale (statistician's verdict,
docs/specs/stat-m5-joincount-scale.md; ADR-0020 Amendment 1).

In ``autocorrelation.tsv`` a ``join_count_bb`` row's ``value`` is the raw count BB_k, while
``null_mean`` / ``null_lo`` / ``null_hi`` and ``effect`` are on the Pearson-residual statistic
S_k; S_k = ``effect + null_mean``. V10 panel (b) and its Plotly twin plotted the raw ``value``
against the S-scale band. Required: both plot S_k against ``null_lo`` / ``null_hi``; the (b) y
label names the Pearson-residual join statistic; the twin's hover shows S_k and the raw BB_k.

Interface assumed (test-author choice, reported to the lead): the V10 data TSV
(``giap.viz.correlogram.COLUMNS``, which the twin reads) gains the table's ``effect`` column so
the twin can form S_k. Panel (a) (Moran's I) is unchanged: it plots ``value``.

Fixture, in the K-12 proportions: lag 1 raw BB 1,973, effect 1,117, null mean 0, band ±130
(S_k = 1,117, far above the band); lag 2 raw BB 1,900, effect 20 (S_k = 20, null-like, inside
the band). Every other join-count row has ``value`` far from its band, ``effect`` = 0.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers import twins as tw
from helpers.viz_landscape import correlogram_rows, gene, landscape_out, markers, plot_spec

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import style
from giap.viz.base import Instance, Skip

CHR_A = ("chrA", 100_000, "circular", "chromosome")
RAW = {1: 1_973.0, 2: 1_900.0}
EFFECT = {1: 1_117.0, 2: 20.0}
HALF = 130.0


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


def _rows() -> list[dict[str, object]]:
    rows = correlogram_rows("chrA", lags=tuple(range(1, 21)))
    for r in rows:
        if r["statistic"] != "join_count_bb":
            continue
        lag = int(str(r["lag"]))
        r["value"] = RAW.get(lag, 1_800.0 - lag)
        r["null_mean"], r["null_lo"], r["null_hi"] = 0.0, -HALF, HALF
        r["effect"] = EFFECT.get(lag, 0.0)
        r["q_value"] = r["p_value"] = 0.001 if lag == 1 else 0.6
    return rows


def _render(tmp: Path, cfg: GiapConfig) -> Instance:
    out = landscape_out(tmp, replicons=(CHR_A,), genes=[gene("g1", 10_000, 1.0)],
                        autocorr=_rows())  # fmt: skip
    result = plot_spec("V10").fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    return result[0]


def _ax(inst: Instance, label: str):  # type: ignore[no-untyped-def]
    (ax,) = [a for a in inst.figure.axes if a.get_label() == label]
    return ax


def _ys_at(inst: Instance, label: str, lag: int) -> list[float]:
    return [y for x, y, _, _ in markers(_ax(inst, label)) if x == float(lag)]


# --------------------------------------------------------------------- the static figure


def test_b_plots_effect_plus_null_mean_not_the_raw_count(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(tmp_path, cfg)
    ys = _ys_at(inst, "b:chrA", 1)
    assert ys == [pytest.approx(1_117.0)], ys
    assert all(abs(y - 1_973.0) > 1.0 for x, y, _, _ in markers(_ax(inst, "b:chrA")))


def test_a_null_like_row_is_drawn_inside_the_band(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(tmp_path, cfg)
    ys = _ys_at(inst, "b:chrA", 2)
    assert ys == [pytest.approx(20.0)], ys
    assert -HALF < ys[0] < HALF
    for lag in range(3, 21):  # effect 0: S_k = null mean, inside the band
        assert _ys_at(inst, "b:chrA", lag) == [pytest.approx(0.0)], lag


def test_b_y_limits_frame_the_band_and_the_statistic(tmp_path: Path, cfg: GiapConfig) -> None:
    lo, hi = _ax(_render(tmp_path, cfg), "b:chrA").get_ylim()
    assert lo <= -HALF and 1_117.0 <= hi < 1_973.0, (lo, hi)


def test_b_y_label_names_the_pearson_residual_statistic(tmp_path: Path, cfg: GiapConfig) -> None:
    label = _ax(_render(tmp_path, cfg), "b:chrA").get_ylabel()
    assert "pearson" in label.lower(), label
    assert "join" in label.lower(), label


def test_a_still_plots_morans_i_value(tmp_path: Path, cfg: GiapConfig) -> None:
    """Guard: Moran's I rows are on one scale; panel (a) keeps ``value``."""
    ys = _ys_at(_render(tmp_path, cfg), "a:chrA", 1)
    assert ys == [pytest.approx(0.16)], ys


def test_data_tsv_carries_effect(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(tmp_path, cfg)
    assert "effect" in inst.data.columns, list(inst.data.columns)
    d = inst.data
    row = d[(d["statistic"] == "join_count_bb") & (d["null_model"] == "tu_block")
            & (d["lag"] == 1)].iloc[0]  # fmt: skip
    assert float(row["effect"]) + float(row["null_mean"]) == pytest.approx(1_117.0)
    assert float(row["value"]) == pytest.approx(1_973.0)


# --------------------------------------------------------------------- the twin


def _twin_frame():  # type: ignore[no-untyped-def]
    from giap.viz import correlogram

    cols = tuple(correlogram.COLUMNS)
    assert "effect" in cols, f"V10 data columns {cols} carry no effect: the twin cannot form S_k"
    rows = []
    for r in _rows():
        if r["statistic"] not in ("morans_i", "join_count_bb"):
            continue
        q = r["q_value"]
        full = {**r, "significant": q is not None and float(str(q)) <= 0.05}
        rows.append(tuple(full.get(c) for c in cols))
    return tw.frame(cols, rows)


def _twin_bb_points(s: dict) -> list[tuple[float, float]]:  # type: ignore[type-arg]
    ys = []
    for t in tw.traces(s):
        if t.get("fill", "none") != "none" or "markers" not in tw.mode(t):
            continue
        ys += tw.point_pairs([t])
    return ys


def test_twin_plots_effect_plus_null_mean() -> None:
    s = tw.spec("V10", _twin_frame())
    pts = _twin_bb_points(s)
    assert (1.0, pytest.approx(1_117.0)) in [(x, y) for x, y in pts], pts
    assert (2.0, pytest.approx(20.0)) in [(x, y) for x, y in pts], pts
    assert not any(abs(y - 1_973.0) < 1.0 for _, y in pts), pts


def _numbers(text: str) -> list[float]:
    return [float(m.replace(",", "")) for m in re.findall(r"\d[\d,]*\.?\d*", text)]


def test_twin_hover_shows_the_statistic_and_the_raw_count() -> None:
    s = tw.spec("V10", _twin_frame())
    hovers = [h for t in tw.traces(s) for h in tw.hovers(t) if h.startswith("lag 1")]
    both = [h for h in hovers
            if any(n == pytest.approx(1_117.0, rel=0.03) for n in _numbers(h))
            and any(n == pytest.approx(1_973.0, rel=0.03) for n in _numbers(h))]  # fmt: skip
    assert both, hovers


def test_twin_y_title_names_the_pearson_residual_statistic() -> None:
    titles = tw.axis_titles(tw.spec("V10", _twin_frame()))
    assert any("pearson" in t.lower() and "join" in t.lower() for t in titles), titles
