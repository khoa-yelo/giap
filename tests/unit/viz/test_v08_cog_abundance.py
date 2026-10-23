"""V8 ``cog_abundance``: log2 abundance per COG category (M9-06; GIAP_SPEC.md §10.2 V8).

Tests the plotted data (``Instance.data`` = ``plots/data/V08_cog_abundance.tsv``), category order,
the strip/violin switch at > 50 measured values, that imputed and partial genes are never drawn,
the detected-fraction labels, colours, titles, degenerate cases and determinism. No ``--mpl``
baselines (M9-25).

Interface assumed where the spec leaves it open:

- ``giap.viz.detection`` registers ``V08`` (name ``cog_abundance``, ``needs=("genes_tsv",)``,
  ``module="core"``, no ``layers``: ``layer_missing:cog`` draws the "none" category, it does not
  skip). One ``Instance`` with ``suffix=None``.
- Data columns, exactly: ``category, super_category, locus_tag, abundance_log2, detected,
  median_log2, n_genes, n_detected, n_imputed``. Gene rows: one per (non-partial gene, letter);
  ``abundance_log2`` is the drawn value (detected, ``abundance_imputed = false``, not NA) and NA
  for every other gene, so an imputed value never appears. Summary rows: ``locus_tag`` NA, one per
  category; ``n_imputed`` counts the category's non-partial genes with ``abundance_imputed``.
  Categories appear in drawing order; gene rows within a category by ``locus_tag``.
- With abundance, the plot is the axes whose y label is "log2 abundance"; x tick labels start
  with the category (letter or "none"), left to right in drawing order. Each drawn value is one
  marker (scatter or line marker); a violin is a ``PolyCollection`` from ``violinplot`` and its
  values are not also drawn as markers; nothing else in that axes is a marker.
- "``n_detected``/``n_genes``" is a text of the figure. "none" is drawn in ``#999999``.
"""

from __future__ import annotations

import math
import os
import statistics
import tempfile
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
import pytest
from helpers.viz_artists import markers, polys, rects, texts, texts_left_to_right
from helpers.viz_out_tree import gene_row, make_out, write_events, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "m4"
VAL_OUT = Path(os.environ.get("GIAP_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = [
    "category",
    "super_category",
    "locus_tag",
    "abundance_log2",
    "detected",
    "median_log2",
    "n_genes",
    "n_detected",
    "n_imputed",
]
YLABEL = "log2 abundance"
NO_ABUNDANCE_TITLE = "No measured abundance; detected fraction per COG category"
GREY = "#999999"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def g(
    tag: str,
    cog: str | None,
    log2: float | None,
    *,
    detected: bool = True,
    imputed: bool = False,
    partial: bool = False,
) -> dict[str, object]:
    return gene_row(
        locus_tag=tag,
        feature_type="cds",
        cog_category=cog,
        abundance_log2=log2,
        abundance_imputed=imputed,
        abundance_type="ibaq" if log2 is not None else "none",
        detected=detected,
        detection_status=(
            "partial_excluded" if partial else "detected" if detected else "not_detected"
        ),
        partial=partial,
    )


#: example 4, plus: an undetected J gene, an imputed (detected) C gene, a partial J gene
EXAMPLE4 = [
    g("j1", "J", 20.0),
    g("j2", "J", 22.0),
    g("j3", "J", 24.0),
    g("j4", "J", None, detected=False),
    g("c1", "C", 18.0),
    g("c2", "C", 19.0),
    g("c3", "C", 30.0, imputed=True),
    g("kl", "KL", 21.0),
    g("jp", "J", 99.0, partial=True),
]
DRAWN = [20.0, 22.0, 24.0, 18.0, 19.0, 21.0, 21.0]


def _out(tmp: Path, genes: Sequence[dict[str, object]] | None, flags: Sequence[str] = ()) -> Path:
    out = make_out(tmp)
    if genes is not None:
        write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if flags:
        records = [{"event": "layer_missing", "data": {"flag": f}} for f in flags]
        write_events(out / "logs" / "events.jsonl", records)
    return out


def _v08() -> base.PlotSpec:
    import giap.viz.detection  # noqa: F401  (registers V07 and V08)

    return base.registered()["V08"]


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _v08().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _ax(inst: Instance) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_ylabel() == YLABEL]
    assert len(found) == 1, [ax.get_ylabel() for ax in inst.figure.axes]
    return found[0]


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _num(value: object) -> float:
    return math.nan if _isna(value) else float(str(value))


def _summary(inst: Instance) -> pd.DataFrame:
    return inst.data[inst.data["locus_tag"].isna()].set_index("category")


def _genes(inst: Instance, category: str) -> pd.DataFrame:
    d = inst.data
    return d[(d["category"].astype(str) == category) & d["locus_tag"].notna()]


def _order(inst: Instance) -> list[str]:
    seen: list[str] = []
    for c in inst.data["category"].astype(str):
        if c not in seen:
            seen.append(c)
    return seen


def _drawn_values(inst: Instance) -> list[float]:
    return sorted(round(m.y, 9) for m in markers(_ax(inst)))


# --- registration and skips -------------------------------------------------------------------


def test_v08_is_registered_as_cog_abundance() -> None:
    spec = _v08()
    assert spec.pid == "V08"
    assert spec.name == "cog_abundance"
    assert spec.needs == ("genes_tsv",)
    assert spec.module == "core"
    assert spec.layers == ()


def test_no_genes_tsv_is_skipped(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _v08().fn(RunView.from_outdir(_out(tmp_path, None)), cfg) == Skip(
        "input_absent:genes_tsv"
    )


# --- data TSV (example 4) ---------------------------------------------------------------------


def test_columns(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    assert inst.suffix is None
    assert list(inst.data.columns) == COLUMNS


def test_example4_category_order(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _order(_render(_out(tmp_path, EXAMPLE4), cfg)) == ["J", "K", "L", "C"]


def test_example4_summary_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    s = _summary(_render(_out(tmp_path, EXAMPLE4), cfg))
    assert sorted(s.index.astype(str)) == ["C", "J", "K", "L"]
    want = {"J": (22.0, 4, 3, 0), "K": (21.0, 1, 1, 0), "L": (21.0, 1, 1, 0), "C": (18.5, 3, 3, 1)}
    for cat, (median, n_genes, n_det, n_imp) in want.items():
        row = s.loc[cat]
        assert _num(row["median_log2"]) == pytest.approx(median), cat
        assert (_num(row["n_genes"]), _num(row["n_detected"]), _num(row["n_imputed"])) == (
            n_genes,
            n_det,
            n_imp,
        ), cat


def test_example4_multi_letter_gene_under_each_letter(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    for cat in ("K", "L"):
        rows = _genes(inst, cat)
        assert list(rows["locus_tag"].astype(str)) == ["kl"]
        assert _num(rows.iloc[0]["abundance_log2"]) == 21.0


def test_gene_rows_sorted_by_locus_tag_with_drawn_values(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    j = _genes(inst, "J")
    assert list(j["locus_tag"].astype(str)) == ["j1", "j2", "j3", "j4"]
    assert [_num(v) for v in j["abundance_log2"]][:3] == [20.0, 22.0, 24.0]
    assert _isna(j.iloc[3]["abundance_log2"])
    c = _genes(inst, "C")
    assert list(c["locus_tag"].astype(str)) == ["c1", "c2", "c3"]
    assert _isna(c.iloc[2]["abundance_log2"])  # imputed: never a drawn value


def test_partial_genes_and_imputed_values_are_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    assert "jp" not in set(inst.data["locus_tag"].dropna().astype(str))
    values = {_num(v) for v in inst.data["abundance_log2"] if not _isna(v)}
    assert values.isdisjoint({30.0, 99.0})


def test_super_category_column(tmp_path: Path, cfg: GiapConfig) -> None:
    s = _summary(_render(_out(tmp_path, EXAMPLE4), cfg))
    assert str(s.loc["J", "super_category"]) == "information"
    assert str(s.loc["C", "super_category"]) == "metabolism"


# --- drawing ----------------------------------------------------------------------------------


def test_example4_drawn_points_are_the_measured_values(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    assert _drawn_values(inst) == sorted(DRAWN)
    assert polys(_ax(inst)) == []  # every category has <= 50 values: strips only


def test_example4_x_ticks_left_to_right(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, EXAMPLE4), cfg)
    labels = texts_left_to_right(inst.figure, _ax(inst), set())
    assert [t.split()[0] for t in labels] == ["J", "K", "L", "C"]


def test_detected_fraction_labels(tmp_path: Path, cfg: GiapConfig) -> None:
    found = texts(_render(_out(tmp_path, EXAMPLE4), cfg).figure)
    assert {"3/4", "1/1", "3/3"} <= set(found)
    assert found.count("1/1") == 2  # K and L


def test_points_coloured_by_super_category(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _ax(_render(_out(tmp_path, EXAMPLE4), cfg))
    by_y = {round(m.y, 9): m for m in markers(ax)}
    assert style.COG_SUPER_COLOUR["information"] in by_y[20.0].colours()
    assert style.COG_SUPER_COLOUR["metabolism"] in by_y[18.0].colours()


def test_jitter_is_deterministic_and_within_the_category(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _out(tmp_path / "a", EXAMPLE4)
    first, second = _render(out, cfg), _render(out, cfg)
    pts = [sorted((round(m.x, 9), round(m.y, 9)) for m in markers(_ax(i))) for i in (first, second)]
    assert pts[0] == pts[1]
    j = [m for m in markers(_ax(first)) if m.y in (20.0, 22.0, 24.0)]
    assert len({round(m.x, 9) for m in j}) == 3  # distinct offsets
    assert max(m.x for m in j) - min(m.x for m in j) < 1.0


def test_input_row_order_does_not_change_the_plot(tmp_path: Path, cfg: GiapConfig) -> None:
    a = _render(_out(tmp_path / "a", EXAMPLE4), cfg)
    b = _render(_out(tmp_path / "b", list(reversed(EXAMPLE4))), cfg)
    pd.testing.assert_frame_equal(a.data.reset_index(drop=True), b.data.reset_index(drop=True))
    pa = sorted((round(m.x, 9), round(m.y, 9)) for m in markers(_ax(a)))
    pb = sorted((round(m.x, 9), round(m.y, 9)) for m in markers(_ax(b)))
    assert pa == pb


@pytest.mark.parametrize(("n", "violin"), [(50, False), (51, True), (60, True)])
def test_strip_violin_switch_above_50(
    tmp_path: Path, cfg: GiapConfig, n: int, violin: bool
) -> None:
    genes = [g(f"j{i:03d}", "J", 15.0 + i / 10) for i in range(n)]
    genes += [g("c1", "C", 10.0), g("c2", "C", 11.0)]
    inst = _render(_out(tmp_path, genes), cfg)
    ax = _ax(inst)
    j_values = {round(15.0 + i / 10, 9) for i in range(n)}
    j_markers = [m for m in markers(ax) if round(m.y, 9) in j_values]
    if violin:
        assert len(polys(ax)) >= 1
        assert j_markers == []
    else:
        assert polys(ax) == []
        assert len(j_markers) == n
    assert {10.0, 11.0} <= {round(m.y, 9) for m in markers(ax)}  # C stays a strip


def test_single_measured_value_is_one_point_no_violin(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, [g("k1", "K", 21.0), g("k2", "K", None, detected=False)]), cfg)
    ax = _ax(inst)
    assert _drawn_values(inst) == [21.0]
    assert polys(ax) == []


def test_none_and_unmeasured_categories_last(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [
        g("n1", None, 30.0),  # highest value, still last
        g("e1", "E", None, detected=False),
        g("c1", "C", 18.0),
        g("j1", "J", 22.0),
    ]
    inst = _render(_out(tmp_path, genes), cfg)
    assert _order(inst) == ["J", "C", "E", "none"]
    ax = _ax(inst)
    none_marker = [m for m in markers(ax) if m.y == 30.0]
    assert none_marker and GREY in none_marker[0].colours()
    assert _isna(_summary(inst).loc["E", "median_log2"])


def test_median_ties_broken_by_letter(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [g("m1", "M", 20.0), g("e1", "E", 20.0), g("a1", "A", 20.0)]
    assert _order(_render(_out(tmp_path, genes), cfg)) == ["A", "E", "M"]


# --- titles -----------------------------------------------------------------------------------


def test_example4_title(tmp_path: Path, cfg: GiapConfig) -> None:
    title = _render(_out(tmp_path, EXAMPLE4), cfg).title
    head = "J (Translation, ribosomal structure and biogenesis) has the highest median abundance ("
    assert title.startswith(head)
    assert "log2" in title
    assert title.endswith(" genes)")


def test_no_measured_abundance_draws_bars(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [
        g("j1", "J", None),
        g("j2", "J", None, detected=False),
        g("c1", "C", None),
        g("c2", "C", 25.0, imputed=True),  # imputed only: still "no measured abundance"
    ]
    inst = _render(_out(tmp_path, genes), cfg)
    assert inst.title == NO_ABUNDANCE_TITLE
    heights = sorted(round(r.get_height(), 9) for ax in inst.figure.axes for r in rects(ax))
    assert heights in ([0.5, 1.0], [50.0, 100.0])
    assert all(_isna(v) for v in inst.data["abundance_log2"])
    assert _order(inst) == ["C", "J"]  # no measured value anywhere: by letter


def test_undetected_abundance_is_not_measured(tmp_path: Path, cfg: GiapConfig) -> None:
    """The validation run: abundance on every gene but none detected."""
    genes = [g(f"g{i}", "J", -0.1 * i, detected=False) for i in range(4)]
    inst = _render(_out(tmp_path, genes), cfg)
    assert inst.title == NO_ABUNDANCE_TITLE
    assert _num(_summary(inst).loc["J", "n_detected"]) == 0


def test_layer_missing_cog_one_none_category(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [g("a", None, 20.0), g("b", None, 22.0), g("c", None, None, detected=False)]
    inst = _render(_out(tmp_path, genes, flags=["layer_missing:cog"]), cfg)
    assert _order(inst) == ["none"]
    s = _summary(inst)
    assert (_num(s.loc["none", "n_genes"]), _num(s.loc["none", "n_detected"])) == (3, 2)
    assert "no COG annotation (layer_missing:cog)" in inst.title


# --- empty state and files --------------------------------------------------------------------


@pytest.mark.parametrize("case", ["header_only", "partial_only"])
def test_no_eligible_genes_is_an_empty_state(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    genes = [] if case == "header_only" else [g("p1", "J", 20.0, partial=True)]
    inst = _render(_out(tmp_path, genes), cfg)
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0
    assert all(not ax.axison for ax in inst.figure.axes)


def test_example6_two_runs_identical_tsv(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _out(tmp_path / "run", EXAMPLE4)
    blobs = []
    for name in ("x", "y"):
        paths = write_instance(
            _render(out, cfg), "V08", "cog_abundance", tmp_path / name, ["svg", "png"], 100
        )
        blobs.append({p.name: p.read_bytes() for p in paths})
    assert blobs[0] == blobs[1]
    tsv = blobs[0]["V08_cog_abundance.tsv"].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)


# --- properties -------------------------------------------------------------------------------

_GENES = st.lists(
    st.tuples(
        st.sampled_from([None, "J", "C", "K", "JC", "KJ"]),
        st.one_of(st.none(), st.sampled_from([10.0, 12.5, 15.0, 20.0])),
        st.booleans(),  # detected
        st.booleans(),  # imputed
        st.integers(0, 4).map(lambda v: v == 0),  # partial (1 in 5)
    ),
    min_size=1,
    max_size=15,
).filter(lambda rows: any(not r[4] for r in rows))


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_GENES)
def test_property_summary_order_and_drawn_values(
    rows: list[tuple[str | None, float | None, bool, bool, bool]],
) -> None:
    style.apply()
    genes = [
        g(f"g{i:02d}", cog, v, detected=det, imputed=imp and v is not None, partial=part)
        for i, (cog, v, det, imp, part) in enumerate(rows)
    ]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), genes), resolve_config())

    cats: dict[str, dict[str, list]] = {}
    for cog, v, det, imp, part in rows:
        if part:
            continue
        for letter in sorted(set(cog)) if cog else ["none"]:
            c = cats.setdefault(letter, {"all": [], "det": [], "vals": [], "imp": []})
            c["all"].append(1)
            c["det"].append(det)
            c["imp"].append(imp and v is not None)
            if det and v is not None and not (imp and v is not None):
                c["vals"].append(v)

    def key(letter: str) -> tuple[int, int, float, str]:
        vals = cats[letter]["vals"]
        return (letter == "none", not vals, -statistics.median(vals) if vals else 0.0, letter)

    assert _order(inst) == sorted(cats, key=key)
    s = _summary(inst)
    for letter, c in cats.items():
        assert _num(s.loc[letter, "n_genes"]) == len(c["all"])
        assert _num(s.loc[letter, "n_detected"]) == sum(c["det"])
        assert _num(s.loc[letter, "n_imputed"]) == sum(c["imp"])
        if c["vals"]:
            assert _num(s.loc[letter, "median_log2"]) == pytest.approx(statistics.median(c["vals"]))
        else:
            assert _isna(s.loc[letter, "median_log2"])
    drawn = Counter(round(v, 9) for c in cats.values() for v in c["vals"])
    if drawn:
        assert Counter(_drawn_values(inst)) == drawn
    else:
        assert inst.title.startswith("No measured abundance")


# --- golden fixture and the validation run ----------------------------------------------------


def test_renders_on_the_sim50k_golden_genes(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    (out / "tables" / "genes.tsv").write_bytes((GOLDEN / "genes.tsv").read_bytes())
    inst = _render(out, cfg)
    assert _order(inst) == ["none"]  # the golden genes carry no COG letters
    s = _summary(inst)
    assert _num(s.loc["none", "n_detected"]) > 0
    write_instance(inst, "V08", "cog_abundance", tmp_path / "plots", ["svg", "png"], 100)


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_renders_on_the_validation_run(tmp_path: Path, cfg: GiapConfig) -> None:
    """2026-10-22 validation run: no COG letters, nothing detected."""
    inst = _render(VAL_OUT, cfg)
    assert _order(inst) == ["none"]
    assert _num(_summary(inst).loc["none", "n_detected"]) == 0
    assert inst.title.startswith("No measured abundance")
    paths = write_instance(inst, "V08", "cog_abundance", tmp_path / "plots", ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)
