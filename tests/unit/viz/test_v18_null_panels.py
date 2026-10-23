"""V18 panel (c): null distributions against the observed statistic (M9-21; GIAP_SPEC.md §10.2
V18, §7.3.5(b, d), §7.8.2 step 4; ADR-0018, ADR-0034; docs/specs/M9-21-v18-null-panels.md and its
lead decisions 1-2).

The null draws come from the M9-17 support files, written here from the M9-17 contract
(tests/helpers/support_analysis_data.py and M9-17 lead decision 5), not from M9-17 code:

* ``null_regions``: ``region_id``, ``contig``, ``stat``, ``shift_k``, ``null_value``,
  ``observed`` (the exhaustive form), or ``percentile`` (1..99) in place of ``shift_k`` (the Q2
  grid form). ``null_value`` / ``observed`` are the studentised statistic the test compares
  (ADR-0018), not the reported raw O - E.
* ``null_islands``: per scanned contig, ``direction = abs``, ``replicate`` 1..B, ``max_scan_s``,
  ``observed_max``.
* ``null_determinants``: per model, ``permutation`` 1..B, ``cv_r2``, ``observed_cv_r2_mean``.

They sit at ``OUT/plots/support/<name>.tsv`` (where M9-17's local ``publish`` puts them and
``RunView.support`` reads them).

Interface pinned here where the spec leaves it open (test author; reported to the lead):

* **Axes.** Each small multiple is an axes labelled ``"c:<table> / <family>"`` (the panel-b
  pattern), ``table`` in {``regions_stats``, ``islands``, ``determinants_model``} and ``family``
  ``"<region_id> / o_minus_e"``, ``"<contig> / <direction>"`` or ``"<model>"``. With no small
  multiple at all, one axes labelled ``"c"`` carries "null draws not available". No other axes
  label starts with ``"c"``.
* **Regions.** Only ``stat = o_minus_e`` (the spec's selection), ranked by the ``subset = all``
  ``q_value`` of ``regions_stats`` (ties: ``p_value``, then ``region_id``), among regions that have
  null rows and a finite q; at most 12. Without ``regions_stats`` no region is drawn (no q to rank
  by). x label "test statistic (studentised O <U+2212> E)" (lead decision 1). Caption holds
  ``"p <p>, q <q>"``.
* **Islands.** One small multiple per (contig, direction) group of ``null_islands``, in file
  order; the caption holds ``"FWER p <p>"``, the smallest ``islands`` ``p_value`` on that contig.
* **Determinants.** One per model in file order; caption ``"null p <p> (uncalibrated,
  ADR-0034)"`` when the model's ``null_p`` row in ``determinants_model`` carries the
  ``null_uncalibrated`` flag, ``"null p <p>"`` and no "uncalibrated" otherwise. A model whose
  ``determinants_model`` rows are flagged ``too_few_rows`` is not drawn.
* **Histogram data rows** (every table): ``panel = c``, 30 rows with ``metric`` NA, ``bin`` =
  the left edge (ascending), ``n`` = the count, ``x``/``y``/``value`` NA; the bins cover every
  finite null value (``bin_0 <= min``), value v in bin k when ``bin_k <= v < bin_{k+1}``, the
  last bin takes the rest; NaN / NA null values are not counted. Then one row ``metric =
  observed``, ``value`` = the file's observed statistic. The observed statistic is drawn as a
  vertical line (``axvline``).
* **Grid form** (``percentile`` column): a step CDF (a ``Line2D`` with a ``steps*`` drawstyle,
  no bars), a text in the axes containing "percentile"; data rows: 99 rows ``x`` = null value,
  ``y`` = percentile / 100, ``bin``/``n``/``metric`` NA, then the ``observed`` row.
* **Order.** Panel c rows follow panels a and b; tables in the order regions_stats, islands,
  determinants_model; regions in ascending q.
* **Plots step.** ``plots`` declares the M9-17 directories ``context_stats_support`` and
  ``determinants_support`` as optional inputs; the plots it renders see their files through
  ``view.support(name)``, also when the step writes its own support tables, and whether a file
  sits directly in the directory or under ``support/`` (the two layouts M9-17 allows).

No ``--mpl`` baselines (M9-25 owns them, as for pass 1). Not statistical: nothing is computed
beyond binning, so there is no calibration test.
"""

from __future__ import annotations

import math
import random
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers import support_analysis_data as sa
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.io.tsv import write_tsv
from giap.viz import base, style
from giap.viz.base import Instance, Skip

COLUMNS = ["panel", "table", "family", "bin", "x", "y", "n", "metric", "value"]
N_BINS = 30
MAX_REGIONS = 12
XLABEL_REGIONS = "test statistic (studentised O " + chr(0x2212) + " E)"
NOT_AVAILABLE = "null draws not available"
TABLE_ORDER = ["regions_stats", "islands", "determinants_model"]


def _with_contig(cols: Sequence[str]) -> tuple[str, ...]:
    """M9-17 lead decision 5: ``null_regions`` gains ``contig`` (after ``region_id``)."""
    cols = list(cols)
    if "contig" not in cols:
        cols.insert(1, "contig")
    return tuple(cols)


NR_COLUMNS = _with_contig(sa.COLUMNS["null_regions"])
NR_GRID_COLUMNS = _with_contig(sa.GRID_COLUMNS)
NI_COLUMNS = tuple(sa.COLUMNS["null_islands"])
ND_COLUMNS = tuple(sa.COLUMNS["null_determinants"])


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def _cell(v: object) -> object:
    """Floats as their ``repr`` (exact through the file), None / NaN as ``NA``."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    if isinstance(v, float):
        return repr(v)
    return v


def _write(path: Path, cols: Sequence[str], rows: Iterable[Sequence[object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, list(cols), [[_cell(c) for c in r] for r in rows])
    return path


def _support(out: Path, name: str, cols: Sequence[str], rows: Iterable[Sequence[object]]) -> Path:
    return _write(out / "plots" / "support" / f"{name}.tsv", cols, rows)


def _rs(
    rid: str,
    p: float | None,
    q: float | None,
    *,
    stat: str = "o_minus_e",
    subset: str = "all",
    effect: float = 3.0,
    n_perm: int = 99,
    contig: str = "chr1",
    i: int = 0,
) -> dict[str, object]:
    return {
        "region_id": rid,
        "region_type": "prophage",
        "subset": subset,
        "contig": contig,
        "start": 101 + 1000 * i,
        "end": 900 + 1000 * i,
        "n_genes": 5,
        "n_eligible": 5,
        "observed": 4,
        "expected": 1.0,
        "o_minus_e": effect if stat == "o_minus_e" else None,
        "stat": stat,
        "effect": effect,
        "p_value": p,
        "q_value": q,
        "method": "x",
        "null_model": "circular_shift" if n_perm else "none",
        "n_perm": n_perm,
    }


def _nr(
    rid: str,
    values: Sequence[float],
    observed: float,
    stat: str = "o_minus_e",
    contig: str = "chr1",
) -> list[tuple[object, ...]]:
    return [(rid, contig, stat, k, v, observed) for k, v in enumerate(values, start=1)]


def _nr_grid(
    rid: str,
    values: Sequence[float],
    observed: float,
    stat: str = "o_minus_e",
    contig: str = "chr1",
) -> list[tuple[object, ...]]:
    assert len(values) == 99
    return [(rid, contig, stat, k, v, observed) for k, v in enumerate(values, start=1)]


def _island(iid: str, contig: str, p: float, direction: str = "high", i: int = 0) -> dict:
    return {
        "island_id": iid,
        "direction": direction,
        "contig": contig,
        "start": 1001 + 10_000 * i,
        "end": 5400 + 10_000 * i,
        "tu_first": "TU0001",
        "tu_last": "TU0003",
        "n_tus": 3,
        "n_genes": 5,
        "genes": "g1;g2;g3;g4;g5",
        "mean_abundance_log2": 21.5,
        "detected_frac": 0.8,
        "scan_s": 4.25,
        "block_len": 3,
        "stat": "scan_max",
        "effect": 1.5,
        "p_value": p,
        "method": "scan_max_stat",
        "null_model": "circular_block_bootstrap",
        "n_perm": 999,
    }


def _ni(contig: str, values: Sequence[float], observed: float) -> list[tuple[object, ...]]:
    return [(contig, "abs", b, v, observed) for b, v in enumerate(values, start=1)]


def _nd(model: str, values: Sequence[float], observed: float) -> list[tuple[object, ...]]:
    return [(model, b, v, observed) for b, v in enumerate(values, start=1)]


def _dm(model: str, null_p: float | None, flags: str | None, n_perm: int = 999) -> list[dict]:
    return [
        {"model": model, "metric": "cv_r2_mean", "value": 0.21, "n": 500, "flags": flags},
        {"model": model, "metric": "null_p", "value": null_p, "n_perm": n_perm, "flags": flags},
    ]


def _vals(n: int, seed: int, scale: float = 1.0) -> list[float]:
    rng = random.Random(seed)
    return [round(rng.gauss(0.0, scale), 4) for _ in range(n)]


def _tables(out: Path, **tables: Sequence[Mapping[str, object]] | None) -> Path:
    for name, rows in tables.items():
        if rows is not None:
            write_schema_table(out / "tables" / f"{name}.tsv", name, rows)
    return out


#: a region set with distinct q (R00 smallest), 20 null values each, observed = 2.0 + i / 10
def _regions(out: Path, n: int, *, order: Sequence[int] | None = None) -> dict[str, float]:
    """``n`` regions with q = 0.001 * (i + 1) and p = q / 2; rows written in ``order``. Returns
    region_id -> observed (the null file's studentised statistic)."""
    idx = list(order) if order is not None else list(range(n))
    rows, nulls, obs = [], [], {}
    for i in idx:
        rid = f"R{i:02d}"
        q = 0.001 * (i + 1)
        rows.append(_rs(rid, q / 2, q, i=i))
        rows.append(_rs(rid, 0.5, 0.6, stat="hl_shift", i=i))
        obs[rid] = 2.0 + i / 10
        nulls += _nr(rid, _vals(20, i), obs[rid])
        nulls += _nr(rid, _vals(20, 100 + i), 0.5, stat="hl_shift")
    _tables(out, regions_stats=rows)
    _support(out, "null_regions", NR_COLUMNS, nulls)
    return obs


def _full(tmp: Path) -> Path:
    """Every null panel: 3 regions, chr1 islands, both determinants models (flagged)."""
    out = make_out(tmp)
    _regions(out, 3)
    _tables(out, islands=[_island("ISL1", "chr1", 0.03, i=0), _island("ISL2", "chr1", 0.2, i=1)])
    _support(out, "null_islands", NI_COLUMNS, _ni("chr1", _vals(50, 7, 1.5), 4.2))
    _tables(
        out,
        determinants_model=[
            *_dm("lightgbm", 0.001, "null_uncalibrated"),
            *_dm("elastic_net", 0.002, "null_uncalibrated"),
        ],
    )
    _support(
        out,
        "null_determinants",
        ND_COLUMNS,
        _nd("lightgbm", _vals(999, 11, 0.05), 0.21)
        + _nd("elastic_net", _vals(999, 12, 0.05), 0.15),
    )
    return out


# --- calling the plot -------------------------------------------------------------------------


def _v18() -> base.PlotSpec:
    import giap.viz.diagnostics  # noqa: F401  (registers V18)

    return base.registered()["V18"]


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _v18().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result[0]


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _num(value: object) -> float:
    return float(str(value))


def _c_axes(inst: Instance) -> dict[str, Axes]:
    """``"<table> / <family>"`` -> small-multiple axes of panel c."""
    return {ax.get_label()[2:]: ax for ax in inst.figure.axes if ax.get_label().startswith("c:")}


def _ax(inst: Instance, key: str) -> Axes:
    axes = _c_axes(inst)
    assert key in axes, f"no panel-c axes {key!r}: {[a.get_label() for a in inst.figure.axes]}"
    return axes[key]


def _texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _bars(ax: Axes) -> list[Rectangle]:
    return [
        p
        for p in ax.patches
        if isinstance(p, Rectangle) and p.get_width() > 0 and p.get_height() > 0
    ]


def _vlines(ax: Axes) -> list[float]:
    """x of every vertical line (constant x data)."""
    out = []
    for line in ax.lines:
        xs = [float(x) for x in line.get_xdata()]
        if len(xs) >= 2 and all(x == xs[0] for x in xs):
            out.append(xs[0])
    return out


def _families(inst: Instance, table: str | None = None) -> list[tuple[str, str]]:
    c = _rows(inst, "c")
    pairs = list(dict.fromkeys(zip(c["table"].astype(str), c["family"].astype(str), strict=True)))
    return [p for p in pairs if table is None or p[0] == table]


def _family(inst: Instance, table: str, family: str) -> pd.DataFrame:
    c = _rows(inst, "c")
    return c[(c["table"] == table) & (c["family"] == family)].reset_index(drop=True)


def _bins(inst: Instance, table: str, family: str) -> pd.DataFrame:
    f = _family(inst, table, family)
    return f[f["metric"].map(_isna)].reset_index(drop=True)


def _observed(inst: Instance, table: str, family: str) -> float:
    f = _family(inst, table, family)
    hit = f[f["metric"].astype(str) == "observed"]
    assert len(hit) == 1, f"{table} / {family}: {len(hit)} observed rows\n{f}"
    return _num(hit.iloc[0]["value"])


def _check_histogram(inst: Instance, table: str, family: str, values: Sequence[float]) -> None:
    """The pinned histogram rows: 30 left edges ascending, counts by interval, NaN dropped."""
    b = _bins(inst, table, family)
    assert len(b) == N_BINS, f"{table} / {family}: {len(b)} bin rows"
    edges = [_num(e) for e in b["bin"]]
    counts = [int(_num(n)) for n in b["n"]]
    assert edges == sorted(edges) and len(set(edges)) == N_BINS
    finite = [v for v in values if math.isfinite(v)]
    assert sum(counts) == len(finite)
    assert edges[0] <= min(finite)
    want = [
        sum(1 for v in finite if edges[k] <= v and (k == N_BINS - 1 or v < edges[k + 1]))
        for k in range(N_BINS)
    ]
    assert counts == want
    for col in ("x", "y", "value"):
        assert all(_isna(v) for v in b[col]), col


# --- registration ------------------------------------------------------------------------------


def test_v18_declares_the_null_tables_optional() -> None:
    spec = _v18()
    assert spec.needs == ()
    assert set(spec.optional) >= {"islands", "determinants_model", "regions_stats"}


# --- regions -------------------------------------------------------------------------------------


def test_example1_fifteen_regions_draw_the_twelve_smallest_q(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    out = make_out(tmp_path)
    order = [7, 14, 0, 3, 11, 9, 1, 13, 5, 2, 12, 4, 10, 8, 6]  # rows not in q order
    _regions(out, 15, order=order)
    inst = _render(out, cfg)
    want = [("regions_stats", f"R{i:02d} / o_minus_e") for i in range(MAX_REGIONS)]
    assert _families(inst, "regions_stats") == want
    assert set(_c_axes(inst)) >= {f"{t} / {f}" for t, f in want}
    assert not [k for k in _c_axes(inst) if k.startswith("regions_stats") and "hl_shift" in k]
    assert not any(ax.get_label() == "c" for ax in inst.figure.axes)


def test_region_histogram_rows_and_observed_row(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    obs = _regions(out, 3)
    inst = _render(out, cfg)
    for i in range(3):
        rid = f"R{i:02d}"
        fam = f"{rid} / o_minus_e"
        _check_histogram(inst, "regions_stats", fam, _vals(20, i))
        assert _observed(inst, "regions_stats", fam) == pytest.approx(obs[rid], abs=1e-12)
        ax = _ax(inst, f"regions_stats / {fam}")
        heights = sum(r.get_height() for r in _bars(ax))
        assert heights == pytest.approx(20)
        assert obs[rid] in [pytest.approx(x, abs=1e-12) for x in _vlines(ax)]
    assert list(inst.data.columns) == COLUMNS


def test_observed_line_is_the_studentised_statistic_not_raw_o_minus_e(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    """Acceptance: the line is the null file's ``observed`` (ADR-0018), not regions_stats'
    ``effect`` / ``o_minus_e`` (the reported raw O - E)."""
    out = make_out(tmp_path)
    _tables(out, regions_stats=[_rs("R01", 0.01, 0.02, effect=7.5)])
    _support(out, "null_regions", NR_COLUMNS, _nr("R01", _vals(30, 1), 2.25))
    inst = _render(out, cfg)
    fam = "R01 / o_minus_e"
    assert _observed(inst, "regions_stats", fam) == 2.25
    lines = _vlines(_ax(inst, f"regions_stats / {fam}"))
    assert 2.25 in lines
    assert 7.5 not in lines


def test_region_x_label_names_the_studentised_statistic(tmp_path: Path, cfg: GiapConfig) -> None:
    """Lead decision 1."""
    inst = _render(_full(tmp_path), cfg)
    regions = [ax for k, ax in _c_axes(inst).items() if k.startswith("regions_stats / ")]
    assert len(regions) == 3
    assert {ax.get_xlabel() for ax in regions} == {XLABEL_REGIONS}


def test_region_caption_gives_p_and_q(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(out, regions_stats=[_rs("R01", 0.012, 0.034)])
    _support(out, "null_regions", NR_COLUMNS, _nr("R01", _vals(30, 1), 1.0))
    inst = _render(out, cfg)
    texts = _texts(_ax(inst, "regions_stats / R01 / o_minus_e"))
    assert any("p 0.012, q 0.034" in t for t in texts), texts


def test_q_ties_break_on_p_then_region_id(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    rows = [_rs(f"R{i:02d}", 0.02, 0.05, i=i) for i in range(14)]
    rows[13] = _rs("R13", 0.01, 0.05, i=13)  # same q, smaller p: first
    _tables(out, regions_stats=rows)
    nulls = [r for i in range(14) for r in _nr(f"R{i:02d}", _vals(20, i), 1.0)]
    _support(out, "null_regions", NR_COLUMNS, nulls)
    got = [f for _, f in _families(_render(out, cfg), "regions_stats")]
    assert got == ["R13 / o_minus_e", *(f"R{i:02d} / o_minus_e" for i in range(11))]


def test_constant_null_draws_one_bar_and_the_observed_line(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(out, regions_stats=[_rs("R01", 1.0, 1.0)])
    _support(out, "null_regions", NR_COLUMNS, _nr("R01", [1.5] * 40, 1.5))
    inst = _render(out, cfg)
    fam = "R01 / o_minus_e"
    _check_histogram(inst, "regions_stats", fam, [1.5] * 40)
    nonzero = [n for n in _bins(inst, "regions_stats", fam)["n"] if _num(n) > 0]
    assert [int(_num(n)) for n in nonzero] == [40]
    ax = _ax(inst, f"regions_stats / {fam}")
    assert len(_bars(ax)) == 1
    assert 1.5 in _vlines(ax)


def test_k_zero_region_is_not_drawn(tmp_path: Path, cfg: GiapConfig) -> None:
    """M9-17: a region with K = 0 shifts has no null rows (n_perm 0, p and q NA)."""
    out = make_out(tmp_path)
    _tables(
        out,
        regions_stats=[
            _rs("R00", None, None, n_perm=0, i=0),
            _rs("R01", 0.02, 0.04, i=1),
            _rs("R02", 0.03, 0.05, i=2),
        ],
    )
    nulls = _nr("R01", _vals(20, 1), 1.0) + _nr("R02", _vals(20, 2), 1.0)
    _support(out, "null_regions", NR_COLUMNS, nulls)
    inst = _render(out, cfg)
    assert _families(inst, "regions_stats") == [
        ("regions_stats", "R01 / o_minus_e"),
        ("regions_stats", "R02 / o_minus_e"),
    ]
    assert not [k for k in _c_axes(inst) if "R00" in k]


def test_na_null_values_are_not_counted(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    values = [*_vals(25, 3), math.nan, math.nan]
    _tables(out, regions_stats=[_rs("R01", 0.02, 0.04)])
    _support(out, "null_regions", NR_COLUMNS, _nr("R01", values, 0.3))
    inst = _render(out, cfg)
    _check_histogram(inst, "regions_stats", "R01 / o_minus_e", values)


def test_percentile_grid_is_a_labelled_step_cdf(tmp_path: Path, cfg: GiapConfig) -> None:
    """M9-17 Q2 fallback: ``percentile`` column -> a step CDF, labelled, not a histogram."""
    out = make_out(tmp_path)
    grid = [round(-2.0 + 4.0 * k / 100, 4) for k in range(1, 100)]
    _tables(out, regions_stats=[_rs("R01", 0.02, 0.04), _rs("R02", 0.03, 0.05, i=1)])
    rows = _nr_grid("R01", grid, 1.25) + _nr_grid("R02", grid, -0.5)
    rows += _nr_grid("R01", grid, 0.1, stat="hl_shift")
    _support(out, "null_regions", NR_GRID_COLUMNS, rows)
    inst = _render(out, cfg)
    assert _families(inst, "regions_stats") == [
        ("regions_stats", "R01 / o_minus_e"),
        ("regions_stats", "R02 / o_minus_e"),
    ]
    for rid, obs in (("R01", 1.25), ("R02", -0.5)):
        fam = f"{rid} / o_minus_e"
        ax = _ax(inst, f"regions_stats / {fam}")
        assert _bars(ax) == []
        steps = [ln for ln in ax.lines if str(ln.get_drawstyle()).startswith("steps")]
        assert steps, [ln.get_drawstyle() for ln in ax.lines]
        xs = {round(float(x), 9) for ln in steps for x in ln.get_xdata()}
        assert {round(v, 9) for v in grid} <= xs
        assert any("percentile" in t.lower() for t in _texts(ax)), _texts(ax)
        assert ax.get_xlabel() == XLABEL_REGIONS
        assert obs in _vlines(ax)
        f = _family(inst, "regions_stats", fam)
        pts = f[f["metric"].map(_isna)]
        assert len(pts) == 99
        assert [_num(x) for x in pts["x"]] == pytest.approx(grid)
        assert [_num(y) for y in pts["y"]] == pytest.approx([k / 100 for k in range(1, 100)])
        assert all(_isna(v) for v in pts["bin"]) and all(_isna(v) for v in pts["n"])
        assert _observed(inst, "regions_stats", fam) == obs


def test_regions_need_regions_stats_to_rank(tmp_path: Path, cfg: GiapConfig) -> None:
    """Proposed default: without ``regions_stats`` there is no q, so no region is drawn."""
    out = make_out(tmp_path)
    _support(out, "null_regions", NR_COLUMNS, _nr("R01", _vals(20, 1), 1.0))
    inst = _render(out, cfg)
    assert _families(inst) == []
    assert any(ax.get_label() == "c" for ax in inst.figure.axes)


# --- islands -----------------------------------------------------------------------------------


def test_islands_histogram_observed_max_and_fwer_p(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(
        out,
        islands=[
            _island("ISL1", "chr1", 0.2, i=0),
            _island("ISL2", "chr1", 0.03, "low", i=1),
            _island("ISL3", "plasmid_1", 0.01, i=2),
        ],
    )
    null = _vals(50, 7, 1.5)
    _support(out, "null_islands", NI_COLUMNS, _ni("chr1", null, 4.2))
    inst = _render(out, cfg)
    assert _families(inst, "islands") == [("islands", "chr1 / abs")]
    _check_histogram(inst, "islands", "chr1 / abs", null)
    assert _observed(inst, "islands", "chr1 / abs") == 4.2
    ax = _ax(inst, "islands / chr1 / abs")
    assert 4.2 in _vlines(ax)
    texts = _texts(ax)
    assert any("FWER p 0.03" in t for t in texts), texts
    assert not any("0.01" in t for t in texts)


def test_islands_one_panel_per_scanned_contig(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(out, islands=[_island("ISL1", "chr1", 0.03), _island("ISL2", "chr2", 0.4, i=1)])
    rows = _ni("chr2", _vals(40, 1), 3.0) + _ni("chr1", _vals(40, 2), 5.0)
    _support(out, "null_islands", NI_COLUMNS, rows)
    inst = _render(out, cfg)
    assert _families(inst, "islands") == [("islands", "chr2 / abs"), ("islands", "chr1 / abs")]
    assert _observed(inst, "islands", "chr2 / abs") == 3.0
    assert _observed(inst, "islands", "chr1 / abs") == 5.0


@pytest.mark.parametrize("form", ["absent", "header_only"])
def test_no_scanned_contig_no_islands_panel(tmp_path: Path, cfg: GiapConfig, form: str) -> None:
    out = make_out(tmp_path)
    _regions(out, 2)
    _tables(out, islands=[])
    if form == "header_only":
        _support(out, "null_islands", NI_COLUMNS, [])
    inst = _render(out, cfg)
    assert _families(inst, "islands") == []
    assert len(_families(inst, "regions_stats")) == 2
    assert not [k for k in _c_axes(inst) if k.startswith("islands")]


# --- determinants -------------------------------------------------------------------------------


def test_example2_determinants_one_histogram_per_model(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(out, determinants_model=[*_dm("lightgbm", 0.001, None), *_dm("elastic_net", 0.3, None)])
    lgb, enet = _vals(999, 11, 0.05), _vals(999, 12, 0.05)
    _support(
        out,
        "null_determinants",
        ND_COLUMNS,
        _nd("lightgbm", lgb, 0.21) + _nd("elastic_net", enet, 0.21),
    )
    inst = _render(out, cfg)
    assert _families(inst, "determinants_model") == [
        ("determinants_model", "lightgbm"),
        ("determinants_model", "elastic_net"),
    ]
    for model, null in (("lightgbm", lgb), ("elastic_net", enet)):
        _check_histogram(inst, "determinants_model", model, null)
        assert _observed(inst, "determinants_model", model) == 0.21
        assert 0.21 in _vlines(_ax(inst, f"determinants_model / {model}"))


def test_uncalibrated_caption_while_flagged(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(
        out,
        determinants_model=[
            *_dm("lightgbm", 0.001, "null_uncalibrated"),
            *_dm("elastic_net", 0.001, None),
        ],
    )
    rows = _nd("lightgbm", _vals(999, 1, 0.05), 0.21) + _nd("elastic_net", _vals(999, 2, 0.05), 0.2)
    _support(out, "null_determinants", ND_COLUMNS, rows)
    inst = _render(out, cfg)
    flagged = _texts(_ax(inst, "determinants_model / lightgbm"))
    assert any("null p 0.001 (uncalibrated, ADR-0034)" in t for t in flagged), flagged
    plain = _texts(_ax(inst, "determinants_model / elastic_net"))
    assert any("null p 0.001" in t for t in plain), plain
    assert not any("uncalibrated" in t for t in plain), plain


@pytest.mark.parametrize("with_file", [False, True])
def test_too_few_rows_no_determinants_panel(
    tmp_path: Path, cfg: GiapConfig, with_file: bool
) -> None:
    out = make_out(tmp_path)
    _regions(out, 2)
    _tables(
        out,
        determinants_model=[
            {"model": m, "metric": "cv_r2_mean", "value": None, "flags": "too_few_rows"}
            for m in ("lightgbm", "elastic_net")
        ],
    )
    if with_file:
        _support(out, "null_determinants", ND_COLUMNS, _nd("lightgbm", _vals(9, 1), 0.1))
    inst = _render(out, cfg)
    assert _families(inst, "determinants_model") == []
    assert not [k for k in _c_axes(inst) if k.startswith("determinants_model")]
    assert len(_families(inst, "regions_stats")) == 2


# --- absent files, panel order, panels a and b ------------------------------------------------


def test_example3_no_support_draws_panel_c_text_only(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    _tables(out, regions_stats=[_rs("R01", 0.01, 0.02)])
    inst = _render(out, cfg)
    c = [ax for ax in inst.figure.axes if ax.get_label() == "c"]
    assert len(c) == 1
    assert NOT_AVAILABLE in _texts(c[0])
    assert _c_axes(inst) == {}
    assert len(_rows(inst, "c")) == 0


def test_panels_a_and_b_unchanged_by_the_null_support(tmp_path: Path, cfg: GiapConfig) -> None:
    full = _full(tmp_path / "with")
    bare = make_out(tmp_path / "without")
    for name in ("regions_stats", "islands", "determinants_model"):
        (bare / "tables" / f"{name}.tsv").write_bytes(
            (full / "tables" / f"{name}.tsv").read_bytes()
        )
    with_c, without_c = _render(full, cfg), _render(bare, cfg)
    assert len(_rows(with_c, "c")) > 0
    assert any(ax.get_label() == "c" for ax in without_c.figure.axes)
    ab = with_c.data[with_c.data["panel"].isin(["a", "b"])].reset_index(drop=True)
    pd.testing.assert_frame_equal(
        ab, without_c.data[without_c.data["panel"].isin(["a", "b"])].reset_index(drop=True)
    )
    assert with_c.title == without_c.title


@pytest.mark.parametrize(
    "present",
    [
        ("null_regions",),
        ("null_islands",),
        ("null_determinants",),
        ("null_regions", "null_determinants"),
    ],
)
def test_each_support_file_alone(tmp_path: Path, cfg: GiapConfig, present: tuple[str, ...]) -> None:
    out = _full(tmp_path)
    for name in ("null_regions", "null_islands", "null_determinants"):
        if name not in present:
            (out / "plots" / "support" / f"{name}.tsv").unlink()
    inst = _render(out, cfg)
    table_of = {
        "null_regions": "regions_stats",
        "null_islands": "islands",
        "null_determinants": "determinants_model",
    }
    assert {t for t, _ in _families(inst)} == {table_of[p] for p in present}
    assert not any(ax.get_label() == "c" for ax in inst.figure.axes)


def test_panel_c_rows_follow_b_in_table_order(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_full(tmp_path), cfg)
    panels = list(dict.fromkeys(inst.data["panel"]))
    assert panels == ["b", "c"]  # no detection_summary: no panel a rows
    tables = list(dict.fromkeys(_rows(inst, "c")["table"]))
    assert tables == TABLE_ORDER
    assert [f for _, f in _families(inst, "regions_stats")] == [
        f"R{i:02d} / o_minus_e" for i in range(3)
    ]
    assert len(_families(inst)) == 3 + 1 + 2
    assert len(_c_axes(inst)) == 6


# --- properties ---------------------------------------------------------------------------------


_finite = st.floats(min_value=-50, max_value=50, allow_nan=False, allow_infinity=False, width=32)


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    values=st.lists(st.one_of(_finite, st.just(math.nan)), min_size=1, max_size=200).filter(
        lambda v: any(math.isfinite(x) for x in v)
    ),
    observed=_finite,
)
def test_property_histogram_rows_bin_every_finite_value(
    values: list[float], observed: float
) -> None:
    values = [float(v) for v in values]
    with tempfile.TemporaryDirectory() as tmp:
        out = make_out(Path(tmp))
        _tables(out, regions_stats=[_rs("R01", 0.5, 0.5)])
        _support(out, "null_regions", NR_COLUMNS, _nr("R01", values, float(observed)))
        inst = _render(out, resolve_config())
    _check_histogram(inst, "regions_stats", "R01 / o_minus_e", values)
    assert _observed(inst, "regions_stats", "R01 / o_minus_e") == float(observed)


@settings(max_examples=8, deadline=None)
@given(
    qs=st.lists(st.integers(min_value=1, max_value=10_000), min_size=1, max_size=20, unique=True),
    seed=st.integers(0, 2**16),
)
def test_property_at_most_twelve_smallest_q(qs: list[int], seed: int) -> None:
    rng = random.Random(seed)
    order = list(range(len(qs)))
    rng.shuffle(order)
    with tempfile.TemporaryDirectory() as tmp:
        out = make_out(Path(tmp))
        rows = [_rs(f"R{i:02d}", qs[i] / 20_000, qs[i] / 10_000, i=i) for i in order]
        _tables(out, regions_stats=rows)
        nulls = [r for i in order for r in _nr(f"R{i:02d}", _vals(5, i), 0.0)]
        _support(out, "null_regions", NR_COLUMNS, nulls)
        inst = _render(out, resolve_config())
    want = sorted(range(len(qs)), key=lambda i: qs[i])[:MAX_REGIONS]
    got = [f for _, f in _families(inst, "regions_stats")]
    assert got == [f"R{i:02d} / o_minus_e" for i in want]


# --- determinism, files, style checklist ---------------------------------------------------------


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _full(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V18", "stat_diagnostics", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert written[0] == written[1]
    tsv = written[0][Path("data/V18_stat_diagnostics.tsv")].decode("utf-8").splitlines()
    assert "c\tdeterminants_model\tlightgbm\tNA\tNA\tNA\tNA\tobserved\t0.21" in tsv
    assert "c\tislands\tchr1 / abs\tNA\tNA\tNA\tNA\tobserved\t4.2" in tsv
    assert "c\tregions_stats\tR00 / o_minus_e\tNA\tNA\tNA\tNA\tobserved\t2" in tsv


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GiapConfig) -> None:
    view = RunView.from_outdir(_full(tmp_path))
    first, second = _v18().fn(view, cfg), _v18().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    assert len(_rows(first[0], "c")) > 0
    pd.testing.assert_frame_equal(first[0].data, second[0].data)


@pytest.mark.parametrize("case", ["full", "grid", "no_support"])
def test_style_checklist_passes_every_gated_item(
    tmp_path: Path, cfg: GiapConfig, case: str
) -> None:
    from giap.viz import checklist

    if case == "full":
        out = _full(tmp_path / "run")
        _regions(out, 12)
    elif case == "grid":
        out = make_out(tmp_path / "run")
        grid = [round(-2.0 + 4.0 * k / 100, 4) for k in range(1, 100)]
        _tables(out, regions_stats=[_rs("R01", 0.02, 0.04)])
        _support(out, "null_regions", NR_GRID_COLUMNS, _nr_grid("R01", grid, 1.25))
    else:
        out = make_out(tmp_path / "run")
        _tables(out, regions_stats=[_rs("R01", 0.01, 0.02)])
    inst = _render(out, cfg)
    if case != "no_support":
        assert len(_rows(inst, "c")) > 0
    else:
        assert any(ax.get_label() == "c" for ax in inst.figure.axes)
    files = write_instance(inst, "V18", "stat_diagnostics", tmp_path / "plots", ["svg", "png"], 100)
    bad = [f"{r.item}: {r.detail}" for r in checklist.check_instance(inst, "V18", files)
           if r.passed is False]  # fmt: skip
    assert bad == []


# --- the plots step hands the M9-17 support directories to the plots ----------------------------


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    from helpers.report_fixtures import swap_registry

    return swap_registry(monkeypatch)


def test_plots_step_declares_the_m9_17_support_inputs() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("plots")
    assert {sa.SUPPORT_KEYS["context_stats"], sa.SUPPORT_KEYS["determinants"]} <= set(
        spec.optional_inputs
    )


@pytest.mark.parametrize("layout", ["flat", "support_subdir"])
def test_plots_step_view_reads_the_null_files(
    fake_registry: ModuleType, tmp_path: Path, layout: str
) -> None:
    """The plots rendered by the step see ``null_regions`` / ``null_islands`` /
    ``null_determinants`` from the M9-17 directories, next to the step's own
    ``replicon_landmarks`` (written because ``replicons`` is given)."""
    from helpers.report_fixtures import cfg as small_cfg
    from helpers.report_fixtures import figure, run_out

    from giap.steps import registry as reg

    seen: dict[str, Any] = {}

    @base.plot("V98", "probe")
    def _probe(view: Any, cfg: Any) -> Any:
        for name in ("null_regions", "null_islands", "null_determinants", "replicon_landmarks"):
            seen[name] = view.support(name)
        return [Instance(suffix=None, figure=figure(), data=pd.DataFrame({"a": [1]}), title="1")]

    sub = "support" if layout == "support_subdir" else ""
    ctx_dir = tmp_path / "work" / "context_stats_support"
    det_dir = tmp_path / "work" / "determinants_support"
    _write(ctx_dir / sub / "null_regions.tsv", NR_COLUMNS, _nr("R01", [0.5, 1.5], 1.0))
    _write(ctx_dir / sub / "null_islands.tsv", NI_COLUMNS, _ni("chr1", [2.0, 3.0], 4.0))
    _write(det_dir / sub / "null_determinants.tsv", ND_COLUMNS, _nd("lightgbm", [0.1], 0.2))
    reps = tmp_path / "inputs" / "replicons.parquet"
    reps.parent.mkdir(parents=True)
    pd.DataFrame(
        [("chr1", 25_000, "circular", "chromosome")],
        columns=["contig", "length", "topology", "kind"],
    ).astype({"length": "int64"}).to_parquet(reps, index=False)

    out = run_out(tmp_path)
    inputs = {
        "preflight": out / "preflight" / "preflight.json",
        "replicons": reps,
        sa.SUPPORT_KEYS["context_stats"]: ctx_dir,
        sa.SUPPORT_KEYS["determinants"]: det_dir,
    }
    reg.run_step("plots", small_cfg(), inputs, {"plots": tmp_path / "work" / "plots"})
    assert seen["replicon_landmarks"] is not None
    for name, n in (("null_regions", 2), ("null_islands", 2), ("null_determinants", 1)):
        frame = seen[name]
        assert frame is not None, name
        assert len(frame) == n, name
    assert list(seen["null_regions"].columns) == list(NR_COLUMNS)
