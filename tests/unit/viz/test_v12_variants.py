"""V12 ``variants``: where do peptides support ALT residues? (M9-18; GIAP_SPEC.md §10.2 V12, §7.5,
§9.3 ``variants.tsv`` / ``peptide_hits.tsv``, §9.4 flag order; lead decisions 1-3 of 2026-10-22).

Tests the plotted data (the ``Instance`` data frame = ``plots/data/V12_variants.tsv``), the two
panels, the title, the degenerate cases, strand and origin handling, truncation, determinism and
the colours the spec names. No ``--mpl`` baselines (Khoa approves those at M9-25).

Interface assumed where the spec leaves it open (M9-18 test author):

- The plot lives in ``giap.viz.variants``; importing it registers ``V12`` (name ``variants``,
  module ``core``, ``needs`` ⊇ {``variants``}, ``optional`` ⊇ {``peptide_hits``}). It returns one
  ``Instance`` with ``suffix=None`` and ``coords=()``. ``FLAG_COLOUR`` maps every §9.4 variant
  flag (``giap.variants.flags.FLAG_ORDER``) to a distinct Okabe-Ito colour, assigned in §9.4
  order (colour indices in ``style.OKABE_ITO`` increase along the flag order).
- The two panel axes carry ``Axes.get_label()`` ``"a"`` and ``"b"``. Panel (a) draws its bars as
  ``Rectangle`` patches whose height is the count (stacked: a flag's segments sum to its ``n``)
  and prints each "all" total with ``ax.text`` / ``bar_label`` (``ax.texts``). Panel (b) is drawn
  with ``Text`` artists in ``ax.texts``, in reading order, so ``"".join`` of their strings holds
  each line; the ALT residue of each aligned peptide is a separate ``Text`` in bold ``#D55E00``.
- Data columns (set; ``panel`` first): ``panel, call, confidence, flag, n, rank, variant_id,
  gene, aa_pos, ref_aa, alt_aa, alt_type, n_alt, n_ref, class_error_rate, ref_window,
  peptides``. Panel (a) ``flag`` is a §9.4 flag, ``"all"`` or ``"no flag"``; rows with ``n = 0``
  may be present. Panel (b) ``peptides`` is the ``;``-joined distinct covering peptides;
  ``class_error_rate`` is a number, or NA / ``"unestimated"``; ``ref_window`` is NA without
  ``plots/support/variant_windows.tsv``.
- "Accepted ``peptide_hits`` rows": every row of ``peptide_hits.tsv`` (the mapping step writes
  accepted placements only). A placement covers a site when it is on the site's contig and
  strand and the whole codon lies in ``[start, end]``; an origin-crossing placement (``end`` >
  replicon length) covers the wrapped codon past the origin (positions are protein positions).
- Title counts ``h`` / ``c`` and the class-error range are over ``call = ALT`` sites; the range
  is printed with each rate's shortest ``%g`` form and an en dash.
"""

from __future__ import annotations

import math
import re
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers import viz_variants as vv
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

COLUMNS = {
    "panel", "call", "confidence", "flag", "n", "rank", "variant_id", "gene", "aa_pos",
    "ref_aa", "alt_aa", "alt_type", "n_alt", "n_ref", "class_error_rate", "ref_window",
    "peptides",
}  # fmt: skip
VERMILLION = "#D55E00"
NOTE = "length-changing: shown unaligned"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- calling the plot -------------------------------------------------------------------------


def _v12() -> base.PlotSpec:
    import giap.viz.variants  # noqa: F401  (registers V12)

    return base.registered()["V12"]


def _result(out: Path, cfg: GiapConfig) -> list[Instance] | Skip:
    return _v12().fn(RunView.from_outdir(out), cfg)


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _result(out, cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _main(tmp: Path, cfg: GiapConfig, *, windows: bool = True, hits: bool = True) -> Instance:
    sites = vv.main_variants()
    alt = [s for s in sites if s["call"] == "ALT"]
    out = vv.build_out(
        tmp,
        sites,
        vv.main_hits() if hits else None,
        vv.windows_frame(alt) if windows else None,
    )
    return _render(out, cfg)


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _int(value: object) -> int:
    return int(float(str(value)))


def _rate(value: object) -> float | None:
    if _isna(value) or str(value) == "unestimated":
        return None
    return float(str(value))


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _counts(inst: Instance) -> dict[tuple[str, str, str], int]:
    a = _rows(inst, "a")
    keys = list(zip(a["call"], a["confidence"], a["flag"], strict=True))
    assert len(keys) == len(set(keys)), a
    return {(str(c), str(k), str(f)): _int(n) for (c, k, f), n in zip(keys, a["n"], strict=True)}


def _nonzero(inst: Instance) -> dict[tuple[str, str, str], int]:
    return {k: v for k, v in _counts(inst).items() if v}


def _order(inst: Instance) -> list[str]:
    b = _rows(inst, "b")
    assert [_int(r) for r in b["rank"]] == list(range(1, len(b) + 1)), b
    return [str(v) for v in b["variant_id"]]


def _b_row(inst: Instance, vid: str) -> pd.Series:
    b = _rows(inst, "b")
    sel = b[b["variant_id"] == vid]
    assert len(sel) == 1, (vid, b)
    return sel.iloc[0]


def _peptides(value: object) -> set[str]:
    return set() if _isna(value) or str(value) == "" else set(str(value).split(";"))


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _joined(ax: Axes) -> str:
    return "".join(t.get_text() for t in ax.texts)


def _all_texts(inst: Instance) -> list[str]:
    return [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]


def _bold(t: Text) -> bool:
    w = t.get_fontweight()
    return w in ("bold", "heavy", "extra bold", "black", "semibold", "demibold", "demi") or (
        isinstance(w, int | float) and w >= 600
    )


def _bold_vermillion(ax: Axes) -> list[str]:
    return [
        t.get_text()
        for t in ax.texts
        if t.get_text() and _bold(t) and mcolors.to_hex(t.get_color()).upper() == VERMILLION
    ]


def _site(vid: str, **kw: object) -> dict[str, object]:
    """A GSIM_00105 ALT site; distinct ``aa_pos`` per call keeps codons apart."""
    return vv.variant(vid, **kw)  # type: ignore[arg-type]


# --- registration, skips, empty states ------------------------------------------------------


def test_v12_is_registered_as_variants_needing_variants() -> None:
    spec = _v12()
    assert spec.pid == "V12"
    assert spec.name == "variants"
    assert spec.module == "core"
    assert "variants" in spec.needs
    assert "peptide_hits" in spec.optional
    assert "peptide_hits" not in spec.needs


def test_variants_absent_skips_input_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vv.build_out(tmp_path, None, vv.main_hits())
    assert _result(out, cfg) == Skip("input_absent:variants")


def test_example_5_header_only_is_the_empty_state(tmp_path: Path, cfg: GiapConfig) -> None:
    out = vv.build_out(tmp_path, [], [])
    inst = _render(out, cfg)
    assert len(inst.data) == 0
    assert set(inst.data.columns) == COLUMNS
    assert inst.data.columns[0] == "panel"
    assert inst.title == "No ALT site supported by peptides"
    assert inst.suffix is None


def test_genome_relation_same_empty_state_title(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["genome_relation=same"])
    out = vv.build_out(tmp_path, [], [])
    inst = _render(out, cfg)
    assert len(inst.data) == 0
    assert inst.title in (
        "No ALT calls (genome relation same)",
        "No ALT calls (genome relation `same`)",
    )


def test_renders_without_peptide_hits(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg, hits=False)
    assert _order(inst) == ["VAR0001", "VAR0002", "VAR0003", "VAR0004"]
    assert all(not _peptides(v) for v in _rows(inst, "b")["peptides"])


# --- data frame shape -------------------------------------------------------------------------


def test_data_columns_suffix_coords_and_panel_order(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert set(inst.data.columns) == COLUMNS
    assert inst.data.columns[0] == "panel"
    assert inst.suffix is None
    assert inst.coords == ()
    panels = list(inst.data["panel"])
    assert set(panels) == {"a", "b"}
    assert panels == sorted(panels)  # panel a rows first, then panel b


def test_panel_b_carries_the_site_columns(tmp_path: Path, cfg: GiapConfig) -> None:
    r = _b_row(_main(tmp_path, cfg), "VAR0001")
    assert str(r["gene"]) == "GSIM_00105"
    assert _int(r["aa_pos"]) == 214
    assert (str(r["ref_aa"]), str(r["alt_aa"])) == ("Q", "K")
    assert str(r["alt_type"]) == "substitution"
    assert str(r["confidence"]) == "high"
    assert (_int(r["n_alt"]), _int(r["n_ref"])) == (3, 0)
    assert _rate(r["class_error_rate"]) == pytest.approx(0.004)


# --- panel (a) --------------------------------------------------------------------------------


def _example1() -> list[dict[str, object]]:
    return [
        _site("v1", aa_pos=20, flags="snv"),
        _site("v2", aa_pos=40, flags="snv;start_codon"),
        _site("v3", aa_pos=60, call="mixed", confidence="candidate", flags=None),
    ]


def test_example_1_counts_by_call_confidence_flag(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(vv.build_out(tmp_path, _example1(), []), cfg)
    assert _nonzero(inst) == {
        ("ALT", "high", "snv"): 2,
        ("ALT", "high", "start_codon"): 1,
        ("ALT", "high", "all"): 2,
        ("mixed", "candidate", "no flag"): 1,
        ("mixed", "candidate", "all"): 1,
    }


def test_example_1_stacked_segments_sum_to_each_flag_count(tmp_path: Path, cfg: GiapConfig) -> None:
    from giap.viz.variants import FLAG_COLOUR

    inst = _render(vv.build_out(tmp_path, _example1(), []), cfg)
    ax = _ax(inst, "a")
    heights: dict[str, float] = {}
    for p in ax.patches:
        if isinstance(p, Rectangle) and p.get_facecolor()[3] > 0:
            key = mcolors.to_hex(p.get_facecolor()[:3]).upper()
            heights[key] = heights.get(key, 0.0) + abs(float(p.get_height()))
    assert heights.get(FLAG_COLOUR["snv"].upper(), 0.0) == pytest.approx(2.0)
    assert heights.get(FLAG_COLOUR["start_codon"].upper(), 0.0) == pytest.approx(1.0)
    for f in ("mass_equivalent", "multi_nt", "stop_codon"):
        assert heights.get(FLAG_COLOUR[f].upper(), 0.0) == pytest.approx(0.0)


def test_example_1_all_totals_printed_above_the_bars(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(vv.build_out(tmp_path, _example1(), []), cfg)
    tokens = re.findall(r"\d+", " ".join(t.get_text() for t in _ax(inst, "a").texts))
    assert "2" in tokens and "1" in tokens


def test_flag_colours_are_okabe_ito_in_flag_order() -> None:
    from giap.variants.flags import FLAG_ORDER
    from giap.viz.variants import FLAG_COLOUR

    assert tuple(FLAG_ORDER) == vv.FLAG_ORDER
    palette = [c.upper() for c in style.OKABE_ITO.values()]
    colours = [FLAG_COLOUR[f].upper() for f in FLAG_ORDER]
    assert len(set(colours)) == len(colours)
    assert all(c in palette for c in colours)
    idx = [palette.index(c) for c in colours]
    assert idx == sorted(idx)


def test_panel_a_rows_sorted_call_confidence_then_flag_order(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _main(tmp_path, cfg)
    a = _rows(inst, "a")
    calls = {"ALT": 0, "mixed": 1}
    confs = {"high": 0, "candidate": 1}
    keys = [(calls[str(c)], confs[str(k)]) for c, k in zip(a["call"], a["confidence"], strict=True)]
    assert keys == sorted(keys)
    for (c, k), grp in a.groupby(["call", "confidence"], sort=False):
        flags = [str(f) for f in grp["flag"] if str(f) in vv.FLAG_ORDER]
        assert flags == sorted(flags, key=vv.FLAG_ORDER.index), (c, k, flags)


def test_main_fixture_panel_a_counts(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _nonzero(_main(tmp_path, cfg)) == {
        ("ALT", "high", "snv"): 2,
        ("ALT", "high", "start_codon"): 1,
        ("ALT", "high", "all"): 2,
        ("ALT", "candidate", "database_derived"): 1,
        ("ALT", "candidate", "no flag"): 1,
        ("ALT", "candidate", "all"): 2,
        ("mixed", "candidate", "no flag"): 1,
        ("mixed", "candidate", "all"): 1,
    }


# --- panel (b): order ---------------------------------------------------------------------------


def test_example_2_order_confidence_then_n_alt(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = [
        _site("v1", aa_pos=20, n_alt=3),
        _site("v2", aa_pos=40, n_alt=5),
        _site("v4", aa_pos=60, confidence="candidate", n_alt=9),
    ]
    assert _order(_render(vv.build_out(tmp_path, sites, []), cfg)) == ["v2", "v1", "v4"]


def test_order_then_class_error_ascending_na_last_then_id(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = [
        _site("VAR0001", aa_pos=10, n_alt=2, class_error_rate="unestimated"),
        _site("VAR0002", aa_pos=20, n_alt=2, class_error_rate="0.01"),
        _site("VAR0004", aa_pos=30, n_alt=2, class_error_rate="0.002"),
        _site("VAR0003", aa_pos=40, n_alt=2, class_error_rate="0.002"),
    ]
    inst = _render(vv.build_out(tmp_path, sites, []), cfg)
    assert _order(inst) == ["VAR0003", "VAR0004", "VAR0002", "VAR0001"]


def test_mixed_sites_are_not_in_panel_b(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _order(_main(tmp_path, cfg)) == ["VAR0001", "VAR0002", "VAR0003", "VAR0004"]


def test_more_than_30_alt_sites_draw_30_and_count_all(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = [
        _site(f"VAR{i:04d}", gene=f"GA{i:02d}", aa_pos=5 + 3 * i, n_alt=100 - i, flags="snv")
        for i in range(1, 36)
    ]
    inst = _render(vv.build_out(tmp_path, sites, []), cfg)
    assert _order(inst) == [f"VAR{i:04d}" for i in range(1, 31)]
    counts = _nonzero(inst)
    assert counts[("ALT", "high", "all")] == 35
    assert counts[("ALT", "high", "snv")] == 35
    text = _joined(_ax(inst, "b"))
    for i in range(1, 31):
        assert f"GA{i:02d}" in text
    for i in range(31, 36):
        assert f"GA{i:02d}" not in text


# --- panel (b): peptides, strands, origin, truncation ------------------------------------------


def test_covering_peptides_are_distinct_same_contig_strand_whole_codon(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _main(tmp_path, cfg)
    for vid, peps in vv.covering().items():
        assert _peptides(_b_row(inst, vid)["peptides"]) == peps, vid


def test_ref_window_comes_from_the_support_table(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert str(_b_row(inst, "VAR0001")["ref_window"]) == vv.P105[203:224]
    assert str(_b_row(inst, "VAR0002")["ref_window"]) == vv.P200[39:60]  # - strand: as is
    assert str(_b_row(inst, "VAR0004")["ref_window"]) == vv.P900[24:45]  # origin gene


def test_header_line_of_a_site(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _joined(_ax(_main(tmp_path, cfg), "b"))
    for token in (
        "GSIM_00105", "aa 214", "Q→K", "substitution", "high", "n_alt 3 / n_ref 0",
        "class error 0.004",
    ):  # fmt: skip
        assert token in text, token


def test_ref_line_and_aligned_peptides(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _joined(_ax(_main(tmp_path, cfg), "b"))
    assert vv.P105[203:224] in text  # the ref window, residues 204-224
    assert vv.PEP1 in text
    assert vv.PEP2 in text


def test_long_peptide_is_truncated_to_21_residues_around_the_site(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    text = _joined(_ax(_main(tmp_path, cfg), "b"))
    shown = vv.P105[203:213] + "K" + vv.P105[214:224]  # PEP_LONG residues 204-224
    assert shown in text
    assert vv.PEP_LONG not in text
    assert vv.P105[202:213] + "K" not in text  # residue 203 is cut
    assert "K" + vv.P105[214:225] not in text  # residue 225 is cut


def test_no_residue_line_longer_than_21(tmp_path: Path, cfg: GiapConfig) -> None:
    for t in _ax(_main(tmp_path, cfg), "b").texts:
        for run in re.findall(r"[A-Z]+", t.get_text()):
            assert len(run) <= 21, t.get_text()


def test_alt_residue_is_bold_vermillion_on_both_strands(tmp_path: Path, cfg: GiapConfig) -> None:
    bold = _bold_vermillion(_ax(_main(tmp_path, cfg), "b"))
    assert set(bold) == {"K", "N", "T"}  # VAR0001 (+), VAR0002 (-), VAR0004 (origin)


def test_minus_strand_site_peptide_and_window(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _joined(_ax(_main(tmp_path, cfg), "b"))
    assert vv.PEP_MINUS in text
    assert vv.P200[39:60] in text
    assert "LVAEKTRGLW" not in text  # the + strand placement at the - site


def test_origin_crossing_gene_site(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _joined(_ax(_main(tmp_path, cfg), "b"))
    assert vv.PEP_ORIGIN in text
    assert vv.P900[24:45] in text


def test_example_4_frameshift_row_present_peptides_unaligned(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _main(tmp_path, cfg)
    r = _b_row(inst, "VAR0003")
    assert _int(r["rank"]) == 3
    assert str(r["alt_type"]) == "frameshift"
    assert _peptides(r["peptides"]) == {vv.PEP_FS}
    text = _joined(_ax(inst, "b"))
    assert NOTE in text
    assert vv.PEP_FS in text


def test_block_site_is_shown_unaligned(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = [_site("VAR0001", aa_pos=245, ref_aa="QW", alt_aa="GT", alt_type="block")]
    inst = _render(vv.build_out(tmp_path, sites, [vv.plus_hit(vv.PEP_FS, 1001, 245)]), cfg)
    ax = _ax(inst, "b")
    assert NOTE in _joined(ax)
    assert vv.PEP_FS in _joined(ax)
    assert not _bold_vermillion(ax)


def test_without_support_panel_b_shows_peptides_only(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg, windows=False)
    assert all(_isna(v) for v in _rows(inst, "b")["ref_window"])
    text = _joined(_ax(inst, "b"))
    assert vv.P105[203:224] not in text
    assert vv.PEP1 in text
    assert _peptides(_b_row(inst, "VAR0001")["peptides"]) == vv.covering()["VAR0001"]


# --- title and wording ------------------------------------------------------------------------


def test_title_counts_and_class_error_range(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    expected = "2 high-confidence and 2 candidate ALT sites (class error 0.004\u20130.02)"
    assert inst.title == expected
    assert expected in _all_texts(inst)


def test_title_when_every_class_error_is_unestimated(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = [
        _site("v1", aa_pos=20, class_error_rate="unestimated"),
        _site("v2", aa_pos=40, confidence="candidate", class_error_rate="unestimated"),
        _site("v3", aa_pos=60, confidence="candidate", class_error_rate="unestimated"),
    ]
    inst = _render(vv.build_out(tmp_path, sites, []), cfg)
    assert inst.title == "1 high-confidence and 2 candidate ALT sites (class error unestimated)"


def test_never_says_mutation_or_modification(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    words = " ".join([inst.title, *_all_texts(inst)]).lower()
    assert "mutation" not in words
    assert "modification" not in words


# --- determinism ------------------------------------------------------------------------------


def test_two_renders_give_equal_data_and_identical_svg(tmp_path: Path, cfg: GiapConfig) -> None:
    sites = vv.main_variants()
    out = vv.build_out(
        tmp_path, sites, vv.main_hits(), vv.windows_frame([s for s in sites if s["call"] == "ALT"])
    )
    a, b = _render(out, cfg), _render(out, cfg)
    pd.testing.assert_frame_equal(a.data, b.data)
    pa = write_instance(a, "V12", "variants", tmp_path / "pa", ["svg"], 100)
    pb = write_instance(b, "V12", "variants", tmp_path / "pb", ["svg"], 100)
    svg_a = next(p for p in pa if p.suffix == ".svg")
    svg_b = next(p for p in pb if p.suffix == ".svg")
    assert svg_a.read_bytes() == svg_b.read_bytes()


# --- properties -------------------------------------------------------------------------------

_SITE = st.tuples(
    st.sampled_from(["ALT", "ALT", "mixed"]),
    st.sampled_from(["high", "candidate"]),
    st.integers(1, 12),
    st.one_of(st.none(), st.sampled_from([0.001, 0.004, 0.01, 0.02, 0.2])),
    st.lists(st.sampled_from(vv.FLAG_ORDER[:7]), unique=True, max_size=3),
)


def _key(s: Mapping[str, object]) -> tuple[int, int, int, float, str]:
    rate = _rate(s["class_error_rate"])
    return (
        0 if s["confidence"] == "high" else 1,
        -_int(s["n_alt"]),
        1 if rate is None else 0,
        0.0 if rate is None else rate,
        str(s["variant_id"]),
    )


def _sites(draw: Sequence[tuple]) -> list[dict[str, object]]:
    return [
        _site(
            f"VAR{i + 1:04d}",
            aa_pos=3 * i + 2,
            call=c,
            confidence=k,
            n_alt=n,
            class_error_rate="unestimated" if r is None else f"{r:g}",
            flags=";".join(sorted(fl, key=vv.FLAG_ORDER.index)) or None,
        )
        for i, (c, k, n, r, fl) in enumerate(draw)
    ]


@settings(
    max_examples=15,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(st.lists(_SITE, min_size=1, max_size=40))
def test_property_counts_and_top_site_order(draw: list[tuple]) -> None:
    cfg = resolve_config()
    sites = _sites(draw)
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(vv.build_out(Path(tmp), sites, []), cfg)
        counts = _counts(inst)
        for call in ("ALT", "mixed"):
            for conf in ("high", "candidate"):
                grp = [s for s in sites if s["call"] == call and s["confidence"] == conf]
                if not grp:
                    assert counts.get((call, conf, "all"), 0) == 0
                    continue
                assert counts[(call, conf, "all")] == len(grp)
                unflagged = sum(1 for s in grp if s["flags"] is None)
                assert counts.get((call, conf, "no flag"), 0) == unflagged
                for f in vv.FLAG_ORDER[:7]:
                    n = sum(1 for s in grp if s["flags"] and f in str(s["flags"]).split(";"))
                    assert counts.get((call, conf, f), 0) == n
        alt = sorted((s for s in sites if s["call"] == "ALT"), key=_key)
        assert _order(inst) == [str(s["variant_id"]) for s in alt[:30]]
        inst.figure.clear()
