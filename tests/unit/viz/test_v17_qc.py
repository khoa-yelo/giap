"""V17 ``qc``: are the inputs sane? (M9-04; GIAP_SPEC.md §10.2 V17, §7.0, §9.3).

Tests the plotted data (the ``Instance`` data frame = ``plots/data/V17_qc.tsv``), the axes
structure, the title, the degenerate cases, determinism and the colours the spec names. No
``--mpl`` baselines (Khoa approves those at M9-25).

Interface assumed where the spec leaves it open:

- The plot lives in ``giap.viz.qc``; importing it registers ``V17`` (name ``qc``,
  ``needs=["preflight"]``, ``optional`` ⊇ {``genes_tsv``, ``peptide_hits``, ``unmapped``}). The
  tests call it through ``giap.viz.base.registered()["V17"].fn(view, cfg)``, which returns a list
  of exactly one ``Instance`` with ``suffix=None`` and ``coords=()``.
- The five panel axes carry ``Axes.get_label()`` ``"a"`` … ``"e"``; the sixth cell (the warning
  list) is the axes labelled ``"warnings"``. Panels (c), (d) are vertical histograms (one
  ``Rectangle`` per drawn bin and layer); panels (b), (e) are horizontal bars (one ``Rectangle``
  per row); (b)'s x scale is ``"log"``.
- Data frame columns are exactly ``panel, group, key, value, n``. Panel (a): count rows have
  ``group = "level"``, ``key`` = the level, ``value`` = count, one row per level (zero counts
  included); message rows have ``group = "message"``, ``key = rule_id``, ``value`` NA, one row per
  ``warn`` / ``error`` check (not truncated: only the drawn list stops at 10 lines). Panel (b):
  ``value`` = length, ``n`` NA. Panel (c): ``value`` = number of distinct sequences of that length,
  only lengths with a non-zero count. Panel (d): ``value`` = genes in the bin, ``n`` = genes in the
  layer (``measured`` or ``imputed``) on every row of that layer. Panel (e): ``abundance_raw`` is in
  ``group = "abundance"``; the feature groups are named ``promoters and transcription``,
  ``translation initiation``, ``codon usage``, ``protein``.
- Rows are sorted by panel, then group, then key; numeric keys (panels c, d) sort numerically.
- "measured abundance" = ``abundance_imputed = false`` and ``abundance_log2`` not NA; the title
  percentage is over the ``partial = false`` genes (example 4). The "no abundance input" title is
  used when no non-partial gene has a measured or imputed ``abundance_log2``.
- The degenerate-case texts ("not available", "no peptide input", "no warnings") are drawn inside
  the corresponding panel axes.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_out_tree import gene_row, make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PathCollection
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
VAL_PREFLIGHT = REPO / "tests" / "data" / "viz" / "val_mvp_preflight.json"
VAL_OUT = Path(os.environ.get("GIAP_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = ["panel", "group", "key", "value", "n"]
LEVELS = ("ok", "info", "warn", "error")
FEATURE_GROUPS: dict[str, tuple[str, ...]] = {
    "promoters and transcription": (
        "promoter_score",
        "promoter_score_tu",
        "tx_rate_max_tu",
        "tx_rate_sum_tu",
        "terminator_conf_downstream",
    ),
    "translation initiation": ("tir_ostir", "sd_dg", "mrna_dg", "spacer_bp", "start_codon"),
    "codon usage": ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10"),
    "protein": ("mw_da", "pi", "gravy", "tm_segments", "n_observable_peptides"),
    "abundance": ("abundance_raw",),
}
TWO_REPLICONS = [
    {"id": "chr1", "kind": "chromosome", "length": 100_000, "topology": "circular"},
    {"id": "plasmid_1", "kind": "plasmid", "length": 8_000, "topology": "circular"},
]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def _check(rule_id: str, level: str, message: str = "") -> dict[str, Any]:
    return {
        "rule_id": rule_id,
        "input": "genome",
        "level": level,
        "rule": f"rule text of {rule_id}",
        "value": None,
        "message": message,
        "exit_code": 1 if level == "error" else 0,
    }


def _preflight(
    checks: Sequence[Mapping[str, Any]], replicons: Sequence[Mapping[str, Any]] = TWO_REPLICONS
) -> dict[str, Any]:
    counts = {lvl: sum(c["level"] == lvl for c in checks) for lvl in LEVELS}
    return {
        "schema_version": "1.0",
        "giap_version": "0.0.1",
        "inputs": [],
        "checks": list(checks),
        "summary": {f"n_{lvl}": counts[lvl] for lvl in LEVELS},
        "evidence_tier": "decoy_picked",
        "routing": {
            "n_proteins_long": 0,
            "n_proteins_peptide_path": 0,
            "n_proteins_dropped": 0,
            "n_peptides": 0,
            "abundance_used": True,
            "abundance_scale": "log2",
        },
        "replicons": [dict(r) for r in replicons],
        "contig_map": [{"id": r["id"], "original_id": r["id"]} for r in replicons],
    }


OK_CHECKS = [_check("genome.parse", "ok"), _check("genome.gc", "ok")]


def _abundance_gene(tag: str, log2: float | None, *, imputed: bool = False, **kw: object) -> dict:
    row = gene_row(
        locus_tag=tag,
        feature_type="cds",
        abundance_log2=log2,
        abundance_raw=None if log2 is None or imputed else 2.0**log2,
        abundance_imputed=imputed,
        abundance_type="ibaq" if log2 is not None else "none",
    )
    row.update(kw)
    return row


#: example 4: three measured, one imputed, one partial gene without abundance
EXAMPLE4_GENES = [
    _abundance_gene("g1", 20.0),
    _abundance_gene("g2", 21.0),
    _abundance_gene("g3", 22.0),
    _abundance_gene("g4", 20.5, imputed=True),
    _abundance_gene("g5", None, partial=True),
]


def _peptide_hit(peptide: str, start: int, **kw: object) -> dict[str, object]:
    row: dict[str, object] = {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "contig": "chr1",
        "start": start,
        "end": start + 3 * len(peptide) - 1,
        "strand": "+",
        "frame": 1,
        "crosses_origin": False,
        "route": "exact",
        "n_diff": 0,
        "uniqueness": "unique",
    }
    row.update(kw)
    return row


def _unmapped(seq: str, kind: str = "peptide", label: str = "unmapped") -> dict[str, object]:
    return {"sequence_id": seq, "sequence": seq, "kind": kind, "label": label, "reason": "no_hit"}


EXAMPLE3_HITS = [
    _peptide_hit("AAAAAAAK", 101),
    _peptide_hit("AAAAAAAK", 5001),
    _peptide_hit("CCCCCCCCCK", 9001),
]
EXAMPLE3_UNMAPPED = [_unmapped("DDDDDDDR"), _unmapped("MKVLAAGIVGLLLAAPLAAQAEPK", "protein")]


def _out(
    tmp: Path,
    *,
    preflight: Mapping[str, Any] | None = None,
    genes: Sequence[Mapping[str, object]] | None = None,
    hits: Sequence[Mapping[str, object]] | None = None,
    unmapped: Sequence[Mapping[str, object]] | None = None,
) -> Path:
    out = make_out(tmp)
    if preflight is not None:
        path = out / "preflight" / "preflight.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(preflight, sort_keys=True), encoding="utf-8")
    if genes is not None:
        write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if hits is not None:
        write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits", hits)
    if unmapped is not None:
        write_schema_table(out / "mapping" / "unmapped.tsv", "unmapped", unmapped)
    return out


def _default_out(tmp: Path, **kw: Any) -> Path:
    args: dict[str, Any] = {
        "preflight": _preflight([*OK_CHECKS, _check("genome.contiguity", "warn", "small")]),
        "genes": EXAMPLE4_GENES,
        "hits": EXAMPLE3_HITS,
        "unmapped": EXAMPLE3_UNMAPPED,
    }
    args.update(kw)
    return _out(tmp, **args)


# --- calling the plot -------------------------------------------------------------------------


def _v17() -> base.PlotSpec:
    import giap.viz.qc  # noqa: F401  (registers V17)

    return base.registered()["V17"]


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _v17().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _num(value: object) -> float:
    return float(str(value))


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _cell(frame: pd.DataFrame, group: str, key: object) -> pd.Series:
    hit = frame[(frame["group"].astype(str) == group) & (frame["key"].astype(str) == str(key))]
    assert len(hit) == 1, f"expected one row {group}/{key}, got {len(hit)}:\n{frame}"
    return hit.iloc[0]


def _int_key_counts(frame: pd.DataFrame, group: str) -> dict[int, int]:
    sub = frame[frame["group"].astype(str) == group]
    return {int(_num(k)): int(_num(v)) for k, v in zip(sub["key"], sub["value"], strict=True)}


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _rects(ax: Axes) -> list[Rectangle]:
    return [p for p in ax.patches if isinstance(p, Rectangle)]


def _drawn(ax: Axes) -> list[Rectangle]:
    return [r for r in _rects(ax) if r.get_width() > 0 and r.get_height() > 0]


def _hex(rect: Rectangle) -> str:
    return mcolors.to_hex(rect.get_facecolor()).upper()


def _marker_points(ax: Axes) -> int:
    n = 0
    for line in ax.get_lines():
        if line.get_marker() not in (None, "", "None", " ", "none"):
            n += len(line.get_xdata())
    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            n += len(coll.get_offsets())
    return n


def _vlines(ax: Axes, linestyle: str) -> list[float]:
    xs = []
    for line in ax.get_lines():
        x = [float(v) for v in line.get_xdata()]
        if line.get_linestyle() == linestyle and len(set(x)) == 1:
            xs.append(x[0])
    return sorted(xs)


# --- registration and skips -------------------------------------------------------------------


def test_v17_is_registered_as_qc_needing_preflight() -> None:
    spec = _v17()
    assert spec.pid == "V17"
    assert spec.name == "qc"
    assert spec.needs == ("preflight",)
    assert {"genes_tsv", "peptide_hits", "unmapped"} <= set(spec.optional)
    assert spec.module == "core"


def test_no_preflight_is_skipped_as_input_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _out(tmp_path, genes=EXAMPLE4_GENES)
    assert _v17().fn(RunView.from_outdir(out), cfg) == Skip("input_absent:preflight")


# --- structure, title, style ------------------------------------------------------------------


def test_one_instance_double_width_with_six_labelled_cells(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    assert round(float(inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    labels = {ax.get_label() for ax in inst.figure.axes}
    assert {"a", "b", "c", "d", "e", "warnings"} <= labels
    assert list(inst.data.columns) == COLUMNS


def test_axis_labels_per_panel(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    assert _ax(inst, "b").get_xlabel() == "length (bp)"
    assert _ax(inst, "b").get_xscale() == "log"
    assert _ax(inst, "c").get_xlabel() == "length (aa)"
    assert _ax(inst, "d").get_xlabel() == "log2 abundance"
    assert _ax(inst, "e").get_xlabel() == "missing (%)"


def test_title_states_counts_and_measured_share(tmp_path: Path, cfg: GiapConfig) -> None:
    """Example 4: 3 of 4 non-partial genes measured; 1 warning, 0 errors."""
    inst = _render(_default_out(tmp_path), cfg)
    expected = "1 preflight warnings, 0 errors; 75 % of genes have a measured abundance"
    assert inst.title == expected
    assert expected in [t.get_text() for t in inst.figure.findobj(Text)]
    assert "express" not in inst.title.lower()


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


def test_rows_are_sorted_by_panel_group_key(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    data = inst.data
    panels = list(data["panel"])
    assert panels == sorted(panels)
    assert set(panels) == {"a", "b", "c", "d", "e"}
    for panel, frame in data.groupby("panel", sort=False):
        groups = [str(g) for g in frame["group"]]
        assert groups == sorted(groups), panel
        for group, sub in frame.groupby("group", sort=False):
            keys = list(sub["key"])
            if panel in ("c", "d"):
                assert [_num(k) for k in keys] == sorted(_num(k) for k in keys), (panel, group)
            else:
                assert [str(k) for k in keys] == sorted(str(k) for k in keys), (panel, group)


# --- panel (a) --------------------------------------------------------------------------------


def test_example1_validation_preflight_counts_and_title(tmp_path: Path, cfg: GiapConfig) -> None:
    pf = json.loads(VAL_PREFLIGHT.read_text(encoding="utf-8"))
    inst = _render(_out(tmp_path, preflight=pf, genes=EXAMPLE4_GENES), cfg)
    a = _rows(inst, "a")
    for level, count in {"ok": 14, "info": 6, "warn": 1, "error": 0}.items():
        row = _cell(a, "level", level)
        assert _num(row["value"]) == count
    msg = _cell(a, "message", "genome.contiguity")
    assert _isna(msg["value"])
    assert len(a[a["group"] == "message"]) == 1
    assert inst.title.startswith("1 preflight warnings, 0 errors")


def test_example1_validation_replicons_are_panel_b(tmp_path: Path, cfg: GiapConfig) -> None:
    pf = json.loads(VAL_PREFLIGHT.read_text(encoding="utf-8"))
    inst = _render(_out(tmp_path, preflight=pf, genes=EXAMPLE4_GENES), cfg)
    b = _rows(inst, "b")
    assert len(b) == 2
    assert _num(_cell(b, "chromosome", "chr1")["value"]) == 100_000
    assert _num(_cell(b, "plasmid", "plasmid_1")["value"]) == 8_000


def test_warning_cell_lists_warn_and_error_messages(tmp_path: Path, cfg: GiapConfig) -> None:
    checks = [
        *OK_CHECKS,
        _check("genome.contiguity", "warn", "largest contig 100000 bp"),
        _check("peptides.content", "error", "no usable peptide"),
        _check("genome.topology", "info", "info is not listed"),
    ]
    inst = _render(_default_out(tmp_path, preflight=_preflight(checks)), cfg)
    text = " ".join(" ".join(_texts(_ax(inst, "warnings"))).split())  # wrapped lines rejoined
    assert "genome.contiguity: largest contig 100000 bp" in text
    assert "peptides.content: no usable peptide" in text
    assert "genome.topology" not in text
    assert "genome.parse" not in text
    a = _rows(inst, "a")
    assert sorted(a[a["group"] == "message"]["key"]) == ["genome.contiguity", "peptides.content"]
    assert inst.title.startswith("1 preflight warnings, 1 errors;")


def test_warning_cell_truncates_to_10_lines(tmp_path: Path, cfg: GiapConfig) -> None:
    warns = [_check(f"genome.rule_{chr(97 + i)}", "warn", f"message {i}") for i in range(13)]
    inst = _render(_default_out(tmp_path, preflight=_preflight([*OK_CHECKS, *warns])), cfg)
    lines = "\n".join(_texts(_ax(inst, "warnings"))).splitlines()
    assert sum("genome.rule_" in line for line in lines) == 10
    assert any(line.strip() == "+ 3 more" for line in lines)
    a = _rows(inst, "a")
    assert len(a[a["group"] == "message"]) == 13  # the data TSV keeps every message


def test_every_check_ok_draws_one_bar_and_no_warnings(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path, preflight=_preflight(OK_CHECKS)), cfg)
    assert len(_drawn(_ax(inst, "a"))) == 1
    assert "no warnings" in _texts(_ax(inst, "warnings"))
    a = _rows(inst, "a")
    assert _num(_cell(a, "level", "ok")["value"]) == 2
    assert (a["group"] == "message").sum() == 0


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(levels=st.lists(st.sampled_from(LEVELS), min_size=1, max_size=25))
def test_property_panel_a_counts_every_check(levels: list[str]) -> None:
    style.apply()
    checks = [
        _check(f"genome.rule_{i:c}", lvl) for i, lvl in zip(range(97, 200), levels, strict=False)
    ]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_default_out(Path(tmp), preflight=_preflight(checks)), resolve_config())
    a = _rows(inst, "a")
    for lvl in LEVELS:
        assert _num(_cell(a, "level", lvl)["value"]) == levels.count(lvl)
    n_msg = levels.count("warn") + levels.count("error")
    assert (a["group"] == "message").sum() == n_msg
    assert len(a) == 4 + n_msg
    assert inst.title.startswith(
        f"{levels.count('warn')} preflight warnings, {levels.count('error')} errors;"
    )


# --- panel (b) --------------------------------------------------------------------------------


def test_example2_contig_bars_coloured_by_kind(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    b = _rows(inst, "b")
    assert sorted(zip(b["group"], b["key"], strict=True)) == [
        ("chromosome", "chr1"),
        ("plasmid", "plasmid_1"),
    ]
    assert all(_isna(v) for v in b["n"])
    rects = _rects(_ax(inst, "b"))
    assert len(rects) == 2
    by_end = {round(r.get_x() + r.get_width()): _hex(r) for r in rects}
    assert by_end == {100_000: "#0072B2", 8_000: "#009E73"}


def test_unknown_kind_is_grey_and_one_contig_one_bar(tmp_path: Path, cfg: GiapConfig) -> None:
    rep = [{"id": "ctg7", "kind": "unknown", "length": 54_321, "topology": "linear"}]
    pf = _preflight(OK_CHECKS, replicons=rep)
    inst = _render(_default_out(tmp_path, preflight=pf), cfg)
    b = _rows(inst, "b")
    assert len(b) == 1
    assert _num(_cell(b, "unknown", "ctg7")["value"]) == 54_321
    rects = _rects(_ax(inst, "b"))
    assert [_hex(r) for r in rects] == ["#999999"]


def test_circular_replicons_carry_a_marker(tmp_path: Path, cfg: GiapConfig) -> None:
    rep = [
        {"id": "chr1", "kind": "chromosome", "length": 2_000_000, "topology": "circular"},
        {"id": "ctg2", "kind": "unknown", "length": 30_000, "topology": "linear"},
        {"id": "p1", "kind": "plasmid", "length": 9_000, "topology": "circular"},
    ]
    inst = _render(_default_out(tmp_path, preflight=_preflight(OK_CHECKS, replicons=rep)), cfg)
    assert _marker_points(_ax(inst, "b")) == 2


def test_dashed_line_at_min_largest_contig_warn(tmp_path: Path) -> None:
    out = _default_out(tmp_path)
    assert _vlines(_ax(_render(out, resolve_config()), "b"), "--") == [1_000_000.0]
    cfg = resolve_config(overrides=["preflight.min_largest_contig_warn=50000"])
    assert _vlines(_ax(_render(out, cfg), "b"), "--") == [50_000.0]


def test_contig_bars_follow_fasta_order(tmp_path: Path, cfg: GiapConfig) -> None:
    rep = [
        {"id": "zeta", "kind": "chromosome", "length": 500_000, "topology": "linear"},
        {"id": "alpha", "kind": "plasmid", "length": 7_000, "topology": "linear"},
        {"id": "mid", "kind": "unknown", "length": 60_000, "topology": "linear"},
    ]
    inst = _render(_default_out(tmp_path, preflight=_preflight(OK_CHECKS, replicons=rep)), cfg)
    ax = _ax(inst, "b")
    ticks = {
        t.get_text(): float(y) for t, y in zip(ax.get_yticklabels(), ax.get_yticks(), strict=True)
    }
    assert set(ticks) == {"zeta", "alpha", "mid"}
    # top to bottom on screen: the largest y first, unless the y axis is inverted
    top_down = sorted(ticks, key=ticks.__getitem__, reverse=not ax.yaxis_inverted())
    assert top_down == ["zeta", "alpha", "mid"]


# --- panel (c) --------------------------------------------------------------------------------


def test_example3_distinct_peptide_lengths(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    c = _rows(inst, "c")
    assert _int_key_counts(c, "placed") == {8: 1, 10: 1}
    assert _int_key_counts(c, "unmapped") == {8: 1}  # the protein row is not a peptide
    ax = _ax(inst, "c")
    drawn = _drawn(ax)
    assert sum(r.get_height() for r in drawn) == pytest.approx(3)
    assert {_hex(r) for r in drawn} == {"#0072B2", "#BBBBBB"}


def test_peptide_length_bounds_are_dotted_lines(tmp_path: Path) -> None:
    out = _default_out(tmp_path)
    assert _vlines(_ax(_render(out, resolve_config()), "c"), ":") == [7.0, 40.0]
    cfg = resolve_config(overrides=["input.peptide_len=[8, 30]"])
    assert _vlines(_ax(_render(out, cfg), "c"), ":") == [8.0, 30.0]


@pytest.mark.parametrize("protein_only", ["absent", "no_peptide_rows"])
def test_no_peptides_draws_no_peptide_input(
    tmp_path: Path, cfg: GiapConfig, protein_only: str
) -> None:
    if protein_only == "absent":
        out = _default_out(tmp_path, hits=None, unmapped=None)
    else:
        out = _default_out(tmp_path, hits=[], unmapped=[_unmapped("MKLV" * 12, "protein")])
    inst = _render(out, cfg)
    assert len(_rows(inst, "c")) == 0
    ax = _ax(inst, "c")
    assert "no peptide input" in _texts(ax)
    assert _drawn(ax) == []


def test_many_peptides_draw_at_most_34_lengths(tmp_path: Path, cfg: GiapConfig) -> None:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    peptides = []
    for i in range(20_400):
        length = 7 + i % 34
        digits, x = [], i
        for _ in range(5):
            digits.append(alphabet[x % 20])
            x //= 20
        peptides.append(("".join(digits) + "G" * 40)[: length - 1] + "K")
    assert len(set(peptides)) == 20_400
    hits = [_peptide_hit(p, 1 + 3 * j) for j, p in enumerate(peptides[:20_000])]
    unmapped = [_unmapped(p) for p in peptides[20_000:]]
    inst = _render(_default_out(tmp_path, hits=hits, unmapped=unmapped), cfg)
    c = _rows(inst, "c")
    keys = {int(_num(k)) for k in c["key"]}
    assert keys <= set(range(7, 41))
    assert sum(_num(v) for v in c["value"]) == 20_400
    assert len(_drawn(_ax(inst, "c"))) <= 2 * 34


# --- panel (d) --------------------------------------------------------------------------------


def test_example4_measured_and_imputed_layers(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    d = _rows(inst, "d")
    assert set(d["group"]) == {"measured", "imputed"}
    for group, total in (("measured", 3), ("imputed", 1)):
        sub = d[d["group"] == group]
        assert sum(_num(v) for v in sub["value"]) == total
        assert {int(_num(n)) for n in sub["n"]} == {total}
    edges = sorted({_num(k) for k in d["key"]})
    assert len(edges) <= 40
    assert edges[0] == pytest.approx(20.0)
    assert edges[-1] < 22.0
    drawn = _drawn(_ax(inst, "d"))
    assert sum(r.get_height() for r in drawn) == pytest.approx(4)
    assert "#E69F00" in {_hex(r) for r in drawn}


def test_imputed_layer_is_orange_only(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [_abundance_gene("g1", 20.0, imputed=True), _abundance_gene("g2", 24.0, imputed=True)]
    inst = _render(_default_out(tmp_path, genes=genes), cfg)
    assert {_hex(r) for r in _drawn(_ax(inst, "d"))} == {"#E69F00"}


@pytest.mark.parametrize("case", ["abundance_type_none", "all_na"])
def test_no_abundance_empty_panel_and_title(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    kind = "none" if case == "abundance_type_none" else "ibaq"
    genes = [
        gene_row(
            locus_tag=f"g{i}", feature_type="cds", abundance_type=kind, abundance_imputed=False
        )
        for i in range(3)
    ]
    inst = _render(_default_out(tmp_path, genes=genes), cfg)
    assert inst.title == "1 preflight warnings; no abundance input"
    assert len(_rows(inst, "d")) == 0
    ax = _ax(inst, "d")
    assert _drawn(ax) == []
    assert _texts(ax)  # the empty-state text
    assert len(_rows(inst, "e")) > 0  # missingness is still drawn


# --- panel (e) --------------------------------------------------------------------------------


def test_example5_cai_missing_share(tmp_path: Path, cfg: GiapConfig) -> None:
    genes = [_abundance_gene(f"g{i}", 20.0 + i, cai=None if i == 0 else 0.7) for i in range(4)]
    genes.append(_abundance_gene("gp", None, partial=True, cai=None))
    inst = _render(_default_out(tmp_path, genes=genes), cfg)
    e = _rows(inst, "e")
    row = _cell(e, "codon usage", "cai")
    assert _num(row["value"]) == pytest.approx(0.25)
    assert int(_num(row["n"])) == 4


def test_missingness_rows_are_the_93_feature_columns(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path), cfg)
    e = _rows(inst, "e")
    expected = {(g, c) for g, cols in FEATURE_GROUPS.items() for c in cols}
    assert set(zip(e["group"].astype(str), e["key"].astype(str), strict=True)) == expected
    assert len(e) == len(expected)
    assert all(0.0 <= _num(v) <= 1.0 for v in e["value"])
    # example 4: abundance_raw is NA on the imputed gene only (1 of 4 non-partial)
    assert _num(_cell(e, "abundance", "abundance_raw")["value"]) == pytest.approx(0.25)
    rects = _rects(_ax(inst, "e"))
    assert len(rects) == len(e)
    widths = sorted(round(r.get_width(), 6) for r in rects)
    assert widths == sorted(round(100 * _num(v), 6) for v in e["value"])


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    genes=st.lists(
        st.tuples(st.booleans(), st.booleans(), st.booleans()),  # partial, cai NA, measured
        min_size=1,
        max_size=15,
    ).filter(lambda g: any(not p for p, _, _ in g))
)
def test_property_missingness_and_measured_share(genes: list[tuple[bool, bool, bool]]) -> None:
    style.apply()
    rows = [
        _abundance_gene(
            f"g{i}",
            20.0 + i if measured else None,
            partial=partial,
            cai=None if cai_na else 0.5,
        )
        for i, (partial, cai_na, measured) in enumerate(genes)
    ]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_default_out(Path(tmp), genes=rows), resolve_config())
    whole = [g for g in genes if not g[0]]
    row = _cell(_rows(inst, "e"), "codon usage", "cai")
    assert int(_num(row["n"])) == len(whole)
    assert _num(row["value"]) == pytest.approx(sum(g[1] for g in whole) / len(whole))
    n_measured = sum(g[2] for g in whole)
    if n_measured:
        pct = int(100 * n_measured / len(whole) + 0.5)
        assert inst.title.endswith(f"; {pct} % of genes have a measured abundance")
    else:
        assert inst.title.endswith("; no abundance input")


# --- degenerate: no genes_tsv -----------------------------------------------------------------


def test_example6_no_genes_tsv_marks_d_and_e_not_available(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_default_out(tmp_path, genes=None), cfg)
    assert set(inst.data["panel"]) == {"a", "b", "c"}
    for label in ("d", "e"):
        ax = _ax(inst, label)
        assert any("not available" in t for t in _texts(ax)), label
        assert _drawn(ax) == []
    assert len(_drawn(_ax(inst, "c"))) > 0
    assert inst.title.startswith("1 preflight warnings")


# --- files and determinism --------------------------------------------------------------------


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _default_out(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(inst, "V17", "qc", tmp_path / name, ["svg", "png", "pdf"], 100)
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V17_qc.svg"),
        Path("V17_qc.png"),
        Path("V17_qc.pdf"),
        Path("data/V17_qc.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V17_qc.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert "e\tcodon usage\tcai\t1\t4" in tsv  # cai is NA on every example-4 gene


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _default_out(tmp_path)
    view = RunView.from_outdir(out)
    first = _v17().fn(view, cfg)
    second = _v17().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


# --- the validation run (local/, slow) --------------------------------------------------------


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_renders_on_the_validation_run(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(VAL_OUT, cfg)
    assert set(inst.data["panel"]) >= {"a", "b", "c"}
    assert inst.title.startswith("1 preflight warnings")
    paths = write_instance(inst, "V17", "qc", tmp_path / "plots", ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)


def test_long_warnings_stay_in_the_figure_and_keep_the_more_line(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    """Review B1 (#197): wrapped messages fit the cell; "+ n more" is drawn, not cut off."""
    long = "largest contig 54321 bp is below the threshold; N50 12345 bp; WWWW mmmm " * 2
    warns = [_check(f"genome.rule_{chr(97 + i)}", "warn", long) for i in range(12)]
    inst = _render(_default_out(tmp_path, preflight=_preflight([*OK_CHECKS, *warns])), cfg)
    cell = _ax(inst, "warnings")
    lines = "\n".join(_texts(cell)).splitlines()
    shown = sum(line.startswith("genome.rule_") for line in lines)
    assert any(line == f"+ {12 - shown} more" for line in lines)
    assert len(lines) <= 20
    fig = inst.figure
    FigureCanvasAgg(fig).draw()  # the writer's renderer; the bare canvas mis-measures text
    box = fig.bbox
    for text in cell.findobj(Text):
        if text.get_text():
            ext = text.get_window_extent()
            assert ext.y0 >= box.y0 - 1 and ext.x1 <= box.x1 + 1


def test_missingness_bars_have_one_colour_and_no_legend(tmp_path: Path, cfg: GiapConfig) -> None:
    """Review B2 (#197): no legend over the bars; the feature groups live in the data TSV."""
    inst = _render(_default_out(tmp_path), cfg)
    ax = _ax(inst, "e")
    assert ax.get_legend() is None
    assert len({_hex(r) for r in _rects(ax) if r.get_width() > 0}) <= 1
