# ruff: noqa: RUF001, RUF003  (Δ, ², ± and – are the spec's own characters)
"""V15 ``determinants``: what explains abundance? (M9-20; GIAP_SPEC.md §10.2 V15, §7.8.2, §9.3
``determinants_model.tsv`` / ``determinants_attribution.tsv``; ADR-0031..0034; M9-17 support
files, its lead decision 5; M9-20 lead decisions 1-3 of 2026-10-23).

Covers examples 1-6, the four panels, the title rule (the elastic net's ΔR² leads, decision 1),
the feature naming ("promoter score" = ``promotech_max``, decision 2), the footnotes, the
degenerate cases (``too_few_rows``, support files absent, a block without a row, ΔR² CI NA, NA
predictions, > 20,000 genes), the §14.4 style checklist (M9-24), determinism, the registry entry,
the ``plots`` step rendering and the Plotly twin. No ``--mpl`` baselines: Khoa approves those at
M9-25 (as M9-18 and M9-23 left theirs). V15 presents the M8 statistics; it makes no new
statistical claim, so there is no calibration test here (the ΔR² t-interval, the null p and the
PD spread are calibrated or flagged by M8 / M9-17).

Fixtures are written from the M9-17 contract (``helpers.support_analysis_data``), not from M9-17
code (``helpers.viz_determinants``).

Interface assumed where the spec leaves it open (M9-20 test author):

- The plot registers on ``giap.steps.plots.load_plot_modules()`` as ``V15`` (name
  ``determinants``, module ``core``), ``needs`` ⊇ {``determinants_model``}, ``optional`` ⊇
  {``determinants_attribution``}. It returns one ``Instance`` (``suffix=None``, ``coords=()``),
  double width (7.2 in), ``Instance.title`` drawn in the figure.
- Data columns exactly ``panel, model, block, feature, x, y, ci_low, ci_high, locus_tag,
  category, flags`` (spec order), ``panel`` ∈ {a, b, c, d}, rows sorted by panel.
  (a) one row per ``delta_r2`` row of the model table: ``model`` verbatim (``elastic_net``,
  ``lightgbm``), ``block``, ``x`` = ΔR², ``ci_low`` / ``ci_high`` copied (NA stays NA).
  (b) ``model = lightgbm``; the rows with a ``y`` are the PD grid: ``feature``, ``x`` =
  ``grid_value``, ``y`` = ``pd_mean``, ``ci_low`` / ``ci_high`` = ``pd_low`` / ``pd_high``
  (rug rows, if any, have ``y`` NA). (c) ``model = lightgbm``, one row per drawn
  (``locus_tag``, ``feature``) point, ``x`` = SHAP, features in drawing order (mean |SHAP|
  descending), ``x`` ascending within a feature. (d) one row per drawn gene: ``x`` =
  ``predicted_log2``, ``y`` = ``abundance_log2``, ``category`` ``unexplained`` /
  ``ribosomal_reference`` for the flagged rows, ``flags`` the attribution flags.
- Axes (``Axes.get_label()``): ``"a"``; the PD sub-panels ``"b_<feature>"`` (``b_cai``,
  ``b_tai``, ``b_tir_ostir``, ``b_promotech_max``) or one ``"b"`` when the PD file is absent;
  ``"c"``; ``"d"``. Other axes (a colour bar) may exist. Footnotes are ``Text`` artists anywhere
  in the figure.
- The ``too_few_rows`` title counts model rows (§7.8.2 rule 1: detected, non-imputed abundance,
  not ribosomal) from ``determinants_features.tsv`` and their distinct ``tu_id``.
- A missing PD feature's sub-panel reads "not available (layer_missing:<layer>)"; with the run
  flag ``layer_missing:ostir`` the ``tir_ostir`` sub-panel names it.
"""

from __future__ import annotations

import json
import math
import re
import tempfile
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import support_analysis_data as sa
from helpers import viz_determinants as vd
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import Collection
from matplotlib.lines import Line2D
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

COLUMNS = ["panel", "model", "block", "feature", "x", "y", "ci_low", "ci_high", "locus_tag",
           "category", "flags"]  # fmt: skip
VERMILLION, ORANGE = "#D55E00", "#E69F00"
EX1_TITLE = (
    "Elastic net explains CV R² 0.21 ± 0.03; largest block ΔR²: transcription 0.08 "
    "(95 % CI 0.05–0.11)"
)
NULL_FOOT = "null p 0.02 (uncalibrated: valid only for exchangeable TUs, ADR-0034)"
UNI_FOOT = "univariate CIs uncalibrated (ADR-0032)"
EN_LEGEND, LGB_LEGEND = "elastic net (primary, ADR-0033)", "LightGBM (secondary check)"
PD_LABEL = "spread over 25 CV models (not a CI)"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- calling the plot -------------------------------------------------------------------------


def _v15() -> base.PlotSpec:
    from giap.steps.plots import load_plot_modules

    load_plot_modules()
    return base.registered()["V15"]


def _result(out: Path, cfg: GiapConfig) -> list[Instance] | Skip:
    return _v15().fn(RunView.from_outdir(out), cfg)


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _result(out, cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _main(tmp: Path, cfg: GiapConfig, **kw: object) -> Instance:
    return _render(vd.build_out(tmp, **kw), cfg)  # type: ignore[arg-type]


def _isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _f(v: object) -> float:
    return math.nan if _isna(v) else float(str(v))


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"].astype(str) == panel].reset_index(drop=True)


def _axes(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _b_axes(inst: Instance) -> list[Axes]:
    found = [ax for ax in inst.figure.axes if ax.get_label() == "b"
             or ax.get_label().startswith("b_")]  # fmt: skip
    assert found, [a.get_label() for a in inst.figure.axes]
    return found


def _ax_texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _texts(inst: Instance) -> list[str]:
    return [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]


def _alltext(inst: Instance) -> str:
    return "\n".join(_texts(inst))


def _colours(ax: Axes) -> set[str]:
    """Hex face / edge / line colours of the data artists of ``ax``."""
    rgba: list[np.ndarray] = []
    for c in ax.collections:
        for get in (c.get_facecolor, c.get_edgecolor):
            arr = np.asarray(get())
            if arr.size:
                rgba += list(arr.reshape(-1, 4))
    for ln in ax.lines:
        for col in (ln.get_color(), ln.get_markerfacecolor(), ln.get_markeredgecolor()):
            try:
                rgba.append(np.asarray(mcolors.to_rgba(col)))
            except ValueError:
                continue
    return {mcolors.to_hex(tuple(float(x) for x in c[:3])).upper() for c in rgba if c[3] > 0}


# --- registration, plots step, twin -----------------------------------------------------------


def test_v15_is_registered_as_determinants() -> None:
    spec = _v15()
    assert spec.pid == "V15"
    assert spec.name == "determinants"
    assert spec.module == "core"
    assert "determinants_model" in spec.needs
    assert "determinants_attribution" in spec.optional
    assert "determinants_attribution" not in spec.needs
    assert base.normalise_ids("v15") == ["V15"]


def test_model_table_absent_skips_input_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vd.build_out(tmp_path, no_model=True)
    assert _result(out, cfg) == Skip("input_absent:determinants_model")


def test_plots_step_renders_v15_files(tmp_path: Path, cfg: GiapConfig) -> None:
    from giap.steps.plots import render_plots

    out = vd.build_out(tmp_path)
    plots = tmp_path / "plots"
    events = render_plots(RunView.from_outdir(out), cfg, plots, ["V15"])
    assert [e.event for e in events] == ["plots.rendered"], events
    files = set(events[0].data["files"])
    assert {"V15_determinants.svg", "V15_determinants.png", "V15_determinants.pdf",
            "data/V15_determinants.tsv"} <= files  # fmt: skip
    back = pd.read_csv(plots / "data" / "V15_determinants.tsv", sep="\t", na_values=["NA"],
                       keep_default_na=False)  # fmt: skip
    assert list(back.columns) == COLUMNS


def test_plotly_twin_registered_and_built_from_the_data_tsv(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    from giap.steps.plots import TWIN_ID_DTYPES
    from giap.viz.interactive import figure_json, has_twin

    inst = _main(tmp_path, cfg)
    assert has_twin("V15")
    write_instance(inst, "V15", "determinants", tmp_path / "plots", ["svg"], 50)
    data = pd.read_csv(tmp_path / "plots" / "data" / "V15_determinants.tsv", sep="\t",
                       na_values=["NA"], keep_default_na=False, dtype=TWIN_ID_DTYPES)  # fmt: skip
    text = figure_json("V15", data, inst.title, None)
    assert text is not None
    fig = json.loads(text)
    assert fig["data"], "the twin draws no trace"
    assert inst.title in json.dumps(fig["layout"], ensure_ascii=False)
    blob = json.dumps(fig, ensure_ascii=False)
    assert "G_0007" in blob  # an unexplained gene is in the (d) hover text


# --- data frame shape -------------------------------------------------------------------------


def test_data_columns_suffix_coords_panel_order_width(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert list(inst.data.columns) == COLUMNS
    assert inst.suffix is None
    assert inst.coords == ()
    panels = [str(p) for p in inst.data["panel"]]
    assert set(panels) == {"a", "b", "c", "d"}
    assert panels == sorted(panels)
    assert round(float(inst.figure.get_size_inches()[0]), 6) == 7.2


# --- title ------------------------------------------------------------------------------------


def test_example_1_title(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert inst.title == EX1_TITLE
    assert EX1_TITLE in _texts(inst)


def test_title_is_in_the_svg(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    write_instance(inst, "V15", "determinants", tmp_path / "p", ["svg"], 50)
    assert inst.title in (tmp_path / "p" / "V15_determinants.svg").read_text(encoding="utf-8")


def test_title_uses_the_elastic_net_not_lightgbm(tmp_path: Path, cfg: GiapConfig) -> None:
    """Decision 1: LightGBM's largest block (codon 0.12) and CV R² (0.25) never lead."""
    inst = _main(tmp_path, cfg)
    assert "codon" not in inst.title
    assert "0.25" not in inst.title and "0.12" not in inst.title


def test_example_2_no_block_ci_excludes_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vd.build_out(tmp_path, model=vd.model_table(vd.EN_DELTAS_NULL))
    t = _render(out, cfg).title
    assert t.startswith("Elastic net explains CV R² 0.21 ± 0.03"), t
    assert t.endswith("no block's ΔR² CI excludes 0"), t
    assert "largest block" not in t


def test_example_6_too_few_rows_title(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vd.build_out(
        tmp_path,
        model=vd.too_few_model_table(),
        attribution=vd.too_few_attribution(60),
        features=vd.feature_rows(60, tu=3, extra_non_model=3),
        support=False,
    )
    t = _render(out, cfg).title
    assert t == "Too few genes for the determinants model (60 rows, 20 TUs)"


def test_title_never_says_expressed(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert "express" not in inst.title.lower()


# --- panel (a) --------------------------------------------------------------------------------


def test_panel_a_rows_are_the_delta_r2_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    a = _rows(_main(tmp_path, cfg), "a")
    got = {(str(m), str(b)): (_f(x), _f(lo), _f(hi))
           for m, b, x, lo, hi in zip(a["model"], a["block"], a["x"], a["ci_low"], a["ci_high"],
                                      strict=True)}  # fmt: skip
    assert len(got) == len(a)
    want = {("elastic_net", b): v for b, v in vd.EN_DELTAS.items()}
    want |= {("lightgbm", b): v for b, v in vd.LGB_DELTAS.items()}
    assert set(got) == set(want)
    for k, (x, lo, hi) in want.items():
        gx, glo, ghi = got[k]
        assert gx == pytest.approx(x), k
        assert (_isna(lo) and math.isnan(glo)) or glo == pytest.approx(lo), k
        assert (_isna(hi) and math.isnan(ghi)) or ghi == pytest.approx(hi), k


def test_block_without_a_row_is_not_drawn(tmp_path: Path, cfg: GiapConfig) -> None:
    a = _rows(_main(tmp_path, cfg), "a")
    lgb = a[a["model"].astype(str) == "lightgbm"]
    assert "context" not in set(lgb["block"].astype(str))


def test_panel_a_legend_primary_and_secondary(tmp_path: Path, cfg: GiapConfig) -> None:
    texts = _ax_texts(_axes(_main(tmp_path, cfg), "a"))
    assert EN_LEGEND in texts
    assert LGB_LEGEND in texts


def test_panel_a_vertical_line_at_zero(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _axes(_main(tmp_path, cfg), "a")
    xs = [np.asarray(ln.get_xdata(), float) for ln in ax.lines]
    zero = [x for x in xs if len(x) and np.allclose(x, 0.0)]
    assert zero, "no vertical line at ΔR² = 0"


def test_delta_r2_ci_na_reads_ci_unestimated(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    a = _rows(inst, "a")
    ctx = a[(a["model"].astype(str) == "elastic_net") & (a["block"].astype(str) == "context")]
    assert len(ctx) == 1
    assert _isna(ctx.iloc[0]["ci_low"]) and _isna(ctx.iloc[0]["ci_high"])
    assert any("CI unestimated" in t for t in _ax_texts(_axes(inst, "a")))


# --- panel (b) --------------------------------------------------------------------------------


def test_panel_b_rows_are_the_pd_grid(tmp_path: Path, cfg: GiapConfig) -> None:
    b = _rows(_main(tmp_path, cfg), "b")
    assert set(b["model"].astype(str)) == {"lightgbm"}
    grid = b[[not _isna(v) for v in b["y"]]].reset_index(drop=True)
    want = vd.pd_frame()
    assert list(grid["feature"].astype(str)) == list(want["feature"])
    assert np.allclose([_f(v) for v in grid["x"]], want["grid_value"])
    assert np.allclose([_f(v) for v in grid["y"]], want["pd_mean"])
    assert np.allclose([_f(v) for v in grid["ci_low"]], want["pd_low"])
    assert np.allclose([_f(v) for v in grid["ci_high"]], want["pd_high"])


def test_panel_b_features_and_promoter_score_name(tmp_path: Path, cfg: GiapConfig) -> None:
    """Decision 2: "promoter score" is ``promotech_max``, the feature name in the data."""
    inst = _main(tmp_path, cfg)
    feats = list(dict.fromkeys(_rows(inst, "b")["feature"].astype(str)))
    assert feats == list(sa.PD_FEATURES)
    assert "promoter_score" not in set(inst.data["feature"].astype(str))
    for f in sa.PD_FEATURES:
        _axes(inst, f"b_{f}")
    assert "promoter score" in _axes(inst, "b_promotech_max").get_xlabel().lower()


def test_panel_b_pd_is_on_each_features_own_scale(tmp_path: Path, cfg: GiapConfig) -> None:
    """M9-17 decision 5: ``tir_ostir`` grid values are raw (not ln) TIR."""
    b = _rows(_main(tmp_path, cfg), "b")
    tir = [_f(x) for x, f, y in zip(b["x"], b["feature"], b["y"], strict=True)
           if str(f) == "tir_ostir" and not _isna(y)]  # fmt: skip
    assert min(tir) == pytest.approx(50.0) and max(tir) == pytest.approx(5000.0)


def test_panel_b_spread_label_is_not_a_ci(tmp_path: Path, cfg: GiapConfig) -> None:
    """25 = ``determinants.outer_folds`` x ``repeats`` (5 x 5 by default)."""
    assert PD_LABEL in _alltext(_main(tmp_path, cfg))


def test_panel_b_spread_is_shaded_and_mean_is_a_line(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _axes(_main(tmp_path, cfg), "b_cai")
    want = vd.pd_frame()
    want = want[want["feature"] == "cai"]
    means = [ln for ln in ax.lines if len(ln.get_ydata()) == sa.PD_GRID_POINTS
             and np.allclose(np.asarray(ln.get_ydata(), float), want["pd_mean"])]  # fmt: skip
    assert means, "no pd_mean line"
    assert any(
        type(c).__name__ in ("PolyCollection", "FillBetweenPolyCollection") for c in ax.collections
    ), "no shaded pd_low-pd_high band"


def test_missing_pd_feature_is_not_available(tmp_path: Path, cfg: GiapConfig) -> None:
    pdf = vd.pd_frame(("cai", "tai", "promotech_max"))
    inst = _main(tmp_path, cfg, pd_=pdf, flags=["layer_missing:ostir"])
    feats = set(_rows(inst, "b")["feature"].astype(str))
    assert "tir_ostir" not in feats
    text = " ".join(_ax_texts(_axes(inst, "b_tir_ostir")))
    assert "not available" in text
    assert "layer_missing:ostir" in text
    assert "not available" not in " ".join(_ax_texts(_axes(inst, "b_cai")))


# --- panel (c) --------------------------------------------------------------------------------


def _mean_abs(shap: pd.DataFrame) -> list[str]:
    m = shap.assign(a=shap["shap"].abs()).groupby("feature")["a"].mean()
    return sorted(m.index, key=lambda f: (-m[f], f))


def test_panel_c_points_are_the_shap_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    c = _rows(_main(tmp_path, cfg), "c")
    assert set(c["model"].astype(str)) == {"lightgbm"}
    shap = vd.shap_frame()
    got = {(str(t), str(f)): _f(x)
           for t, f, x in zip(c["locus_tag"], c["feature"], c["x"], strict=True)}  # fmt: skip
    cols = (shap["locus_tag"], shap["feature"], shap["shap"])
    want = {(str(t), str(f)): float(s) for t, f, s in zip(*cols, strict=True)}
    assert len(got) == len(c)
    assert set(got) == set(want)
    assert all(got[k] == pytest.approx(want[k], rel=1e-5) for k in want)  # %.6g file


def test_panel_c_feature_order_and_x_sorted(tmp_path: Path, cfg: GiapConfig) -> None:
    c = _rows(_main(tmp_path, cfg), "c")
    feats = list(dict.fromkeys(c["feature"].astype(str)))
    assert feats == _mean_abs(vd.shap_frame())
    for f, grp in c.groupby(c["feature"].astype(str), sort=False):
        xs = [_f(x) for x in grp["x"]]
        assert xs == sorted(xs), f


def test_panel_c_caps_at_top_20_features(tmp_path: Path, cfg: GiapConfig) -> None:
    shap = vd.shap_frame(features=vd.SHAP_FEATURES)  # 22 features
    c = _rows(_main(tmp_path, cfg, shap=shap), "c")
    feats = list(dict.fromkeys(c["feature"].astype(str)))
    assert feats == list(vd.SHAP_FEATURES[:20])
    assert len(feats) == sa.SHAP_TOP


def test_panel_c_top_feature_at_the_top(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _axes(_main(tmp_path, cfg), "c")
    pairs = zip(ax.get_yticks(), ax.get_yticklabels(), strict=True)
    ticks = [(float(y), lab.get_text()) for y, lab in pairs if lab.get_text()]
    assert len(ticks) == 20
    top = max(ticks)[1] if not ax.yaxis_inverted() else min(ticks)[1]
    assert "tir_ostir" in top or "TIR" in top or "OSTIR" in top.upper(), ticks


def test_panel_c_axis_header_and_viridis(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    ax = _axes(inst, "c")
    assert ax.get_xlabel() == "SHAP (log2 abundance)"
    assert any("LightGBM (secondary)" in t for t in _ax_texts(ax))
    maps = {c.get_cmap().name for c in ax.collections if c.get_array() is not None}
    assert "viridis" in maps, maps


def test_shap_row_order_does_not_change_the_plot(tmp_path: Path, cfg: GiapConfig) -> None:
    """Deterministic jitter (by rank within the feature): input row order is irrelevant."""
    shap = vd.shap_frame()
    a = _main(tmp_path / "a", cfg, shap=shap)
    b = _main(tmp_path / "b", cfg, shap=shap.sample(frac=1.0, random_state=1))
    pd.testing.assert_frame_equal(a.data, b.data)


# --- panel (d) --------------------------------------------------------------------------------


def test_panel_d_rows_and_na_predictions(tmp_path: Path, cfg: GiapConfig) -> None:
    d = _rows(_main(tmp_path, cfg), "d")
    att = {str(r["locus_tag"]): r for r in vd.attribution_rows()}
    tags = list(d["locus_tag"].astype(str))
    assert len(tags) == len(set(tags))
    assert set(tags) == {t for t, r in att.items() if r["predicted_log2"] is not None}
    assert vd.NA_PRED_TAG not in tags
    for t, x, y in zip(tags, d["x"], d["y"], strict=True):
        assert _f(x) == pytest.approx(float(att[t]["predicted_log2"]))  # type: ignore[arg-type]
        assert _f(y) == pytest.approx(float(att[t]["abundance_log2"]))  # type: ignore[arg-type]


def test_example_4_unexplained_category_and_colour(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    d = _rows(inst, "d")
    cat = dict(zip(d["locus_tag"].astype(str), d["category"].astype(str), strict=True))
    for t in vd.UNEXPLAINED_TAGS:
        assert cat[t] == "unexplained", t
    assert sum(c == "unexplained" for c in cat.values()) == len(vd.UNEXPLAINED_TAGS)
    assert VERMILLION in _colours(_axes(inst, "d"))


def test_example_5_ribosomal_reference_category_and_colour(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    d = _rows(inst, "d")
    cat = dict(zip(d["locus_tag"].astype(str), d["category"].astype(str), strict=True))
    for t in vd.RIBOSOMAL:
        assert cat[t] == "ribosomal_reference", t
    assert ORANGE in _colours(_axes(inst, "d"))


def test_example_5_r2_text_is_the_elastic_net_cv_r2(tmp_path: Path, cfg: GiapConfig) -> None:
    """The R² text quotes the model table (model rows only), not a refit with ribosomal genes."""
    texts = _ax_texts(_axes(_main(tmp_path, cfg), "d"))
    assert any("CV R² 0.21 ± 0.03 (elastic net)" in t for t in texts), texts


def test_panel_d_identity_line(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _axes(_main(tmp_path, cfg), "d")
    ident = [ln for ln in ax.lines if len(ln.get_xdata()) >= 2
             and np.allclose(np.asarray(ln.get_xdata(), float),
                             np.asarray(ln.get_ydata(), float))]  # fmt: skip
    assert ident or any(type(ln).__name__ == "AxLine" for ln in ax.get_children()), "no y = x"


def test_more_than_20000_genes_bin_panel_d(tmp_path: Path, cfg: GiapConfig) -> None:
    n = 20_500
    out = vd.build_out(tmp_path, attribution=vd.many_attribution(n),
                       features=vd.many_features(n), support=False)  # fmt: skip
    inst = _render(out, cfg)
    ax = _axes(inst, "d")
    assert any(c.get_array() is not None for c in ax.collections), "no binned (hexbin) counts"
    for c in inst.figure.findobj(Collection):
        assert len(c.get_offsets()) <= 20_000, c
    for ln in inst.figure.findobj(Line2D):
        assert len(ln.get_xdata()) <= 20_000, ln


# --- footnotes --------------------------------------------------------------------------------


def test_example_3_null_p_footnote(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _alltext(_main(tmp_path, cfg))
    assert NULL_FOOT in text
    assert "(uncalibrated" in text


def test_univariate_footnote_only_with_univariate_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    assert UNI_FOOT in _alltext(_main(tmp_path / "a", cfg))
    out = vd.build_out(tmp_path / "b", model=vd.model_table(univariate=False))
    assert UNI_FOOT not in _alltext(_render(out, cfg))


def test_native_genes_caveat(tmp_path: Path, cfg: GiapConfig) -> None:
    from giap.determinants import CAVEATS

    assert "expect a lower R²" in CAVEATS[0]
    assert "expect a lower R²" in _alltext(_main(tmp_path, cfg))


# --- degenerate cases -------------------------------------------------------------------------


def _too_few(tmp: Path) -> Path:
    return vd.build_out(
        tmp,
        model=vd.too_few_model_table(),
        attribution=vd.too_few_attribution(60),
        features=vd.feature_rows(60, tu=3, extra_non_model=3),
        support=False,
    )


def test_too_few_rows_empties_a_to_c(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_too_few(tmp_path), cfg)
    panels = set(inst.data["panel"].astype(str))
    assert not panels & {"a", "b", "c"}, panels
    assert all(not _isna(x) for x in _rows(inst, "d")["x"])  # NA predictions never drawn


def test_support_files_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg, support=False)
    panels = set(inst.data["panel"].astype(str))
    assert {"a", "d"} <= panels
    assert not panels & {"b", "c"}, panels
    b_text = " ".join(t for ax in _b_axes(inst) for t in _ax_texts(ax))
    assert "not available" in b_text
    assert "not available" in " ".join(_ax_texts(_axes(inst, "c")))
    assert inst.title == EX1_TITLE


# --- style checklist (M9-24) and determinism ---------------------------------------------------


def _gated_failures(inst: Instance, tmp: Path) -> list[str]:
    from giap.viz import checklist

    files = write_instance(inst, "V15", "determinants", tmp, ["svg", "png"], 100)
    return [f"{r.item}: {r.detail}" for r in checklist.check_instance(inst, "V15", files)
            if r.passed is False]  # fmt: skip


def test_checklist_main_fixture(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _gated_failures(_main(tmp_path, cfg), tmp_path / "plots") == []


def test_checklist_without_support(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _gated_failures(_main(tmp_path, cfg, support=False), tmp_path / "plots") == []


def test_checklist_too_few_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _gated_failures(_render(_too_few(tmp_path), cfg), tmp_path / "plots") == []


def test_two_renders_give_equal_data_and_identical_svg(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vd.build_out(tmp_path)
    a, b = _render(out, cfg), _render(out, cfg)
    pd.testing.assert_frame_equal(a.data, b.data)
    pa = write_instance(a, "V15", "determinants", tmp_path / "pa", ["svg"], 100)
    pb = write_instance(b, "V15", "determinants", tmp_path / "pb", ["svg"], 100)
    svg_a = next(p for p in pa if p.suffix == ".svg")
    svg_b = next(p for p in pb if p.suffix == ".svg")
    assert svg_a.read_bytes() == svg_b.read_bytes()


# --- properties -------------------------------------------------------------------------------

#: per block: (value in hundredths, half-width in hundredths); the CI is value ± half-width
_BLOCK = st.tuples(st.integers(-5, 20), st.integers(1, 10))


@settings(
    max_examples=15,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(st.dictionaries(st.sampled_from(vd.BLOCKS), _BLOCK, min_size=1, max_size=5))
def test_property_title_rule_and_panel_a_rows(draw: dict[str, tuple[int, int]]) -> None:
    """The alternative ending iff no elastic-net ΔR² CI excludes 0; when the numerically largest
    block's CI excludes 0 it is named, with its value and CI; one (a) row per input row."""
    cfg = resolve_config()
    en = {b: (v / 100, (v - h) / 100, (v + h) / 100) for b, (v, h) in draw.items()}
    # a CI wholly below 0 is left out: whether it "excludes 0" for the title is open (gap)
    assume(all(hi > 0 for _, _, hi in en.values()))  # type: ignore[operator]
    excl = {b for b, (_, lo, _) in en.items() if lo > 0}  # type: ignore[operator]
    with tempfile.TemporaryDirectory() as tmp:
        out = vd.build_out(Path(tmp), model=vd.model_table(en), support=False)
        inst = _render(out, cfg)
        t = inst.title
        assert t.startswith("Elastic net explains CV R² 0.21 ± 0.03; "), t
        if not excl:
            assert t.endswith("no block's ΔR² CI excludes 0"), (t, en)
        else:
            assert "largest block ΔR²: " in t, (t, en)
            named = re.search(r"largest block ΔR²: (\w+) ", t)
            assert named and named.group(1) in excl, (t, en)
            top = max(en, key=lambda b: (en[b][0], b))
            if top in excl and sum(en[b][0] == en[top][0] for b in en) == 1:
                assert named.group(1) == top, (t, en)
        a = _rows(inst, "a")
        en_rows = a[a["model"].astype(str) == "elastic_net"]
        assert sorted(en_rows["block"].astype(str)) == sorted(en)
        inst.figure.clear()
