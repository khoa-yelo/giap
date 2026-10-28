"""V13 ``novel_orfs``: class counts by call status and the top-20 locus views (M9-12, examples
1-4, degenerate cases, lead decisions 1-4; GIPPA_SPEC.md §10.2 V13, §7.6, §9.3; M9-index
decision 19).

Interface assumed where the spec leaves it open (reported to the lead):

- ``gippa.viz.novel_orfs`` registers ``V13`` as ``"novel_orfs"`` with ``needs = ("novel_orfs",)``
  (genes, peptide hits, BEDs optional). The plot returns: the panel (a) instance (``suffix =
  None``; or the empty state), one locus instance per rendered ORF (``suffix = "locus_<orf_id>"``,
  raw ID; files sanitised by ``gippa.io.plot_files.stem`` as for V3, decision 4) in ``rank``
  order, and the index instance (``suffix = "locus"``, columns exactly ``rank, orf_id,
  rendered``) whenever any locus view is rendered.
- Panel (a) is the axes labelled ``"counts"``; its data columns are exactly ``class, status, n,
  class_fdr_max_called``; rows with ``n = 0`` may be present or absent; the ``called`` row of a
  class carries the maximum ``class_fdr`` of its called rows. Status: ``called = true`` ->
  called; otherwise ``class_fdr`` NA / ``unestimated`` or flag ``class_fdr_unestimated`` ->
  unestimated; otherwise not_called. Bars are ``Rectangle`` patches on that axes coloured
  ``#0072B2`` / ``#BBBBBB`` / ``#E69F00``; class names are the category tick labels.
- Locus views are built by ``gippa.viz.locus.render_locus`` with ``Locus(kind="novel_orf")``
  (spied on both ``gippa.viz.locus`` and, if it imports the name, ``gippa.viz.novel_orfs``);
  ``Locus.start0/end0`` is the 0-based half-open span. The data frame has the V3 columns plus
  one ``track = novel_orf`` row (``id = orf_id``, the span, the strand). The ORF arrow is a patch
  with edge colour ``#CC79A7``; the ORF's peptides are highlighted on the ``"pileup"`` axes with
  an artist coloured ``#CC79A7``; the ``start_best_ostir`` tick is a ``Line2D`` on the
  ``"genes"`` axes at the codon's first base (``+``) within 1 bp.
- Origin-crossing rows: on ``+`` ``end < start_longest``, on ``-`` ``start_longest < end``; the
  span is unwrapped past the length.
- Circular windows: the ``"genes"`` x major formatter labels a position ``x`` as ``x mod length``
  in kb (the real-data fixes); linear windows clip to ``[0, length]``.
"""

from __future__ import annotations

import importlib
import inspect
import itertools
import random
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text

import gippa.viz.locus as locus_mod
from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

CLASSES = (
    "nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
    "pseudogene", "sorf",
)  # fmt: skip
STATUSES = ("called", "not_called", "unestimated")
COUNT_COLUMNS = ["class", "status", "n", "class_fdr_max_called"]
INDEX_COLUMNS = ["rank", "orf_id", "rendered"]
LOCUS_COLUMNS = ["track", "id", "contig", "start", "end", "strand", "value", "category"]
CALLED, NOT_CALLED, UNEST = "#0072B2", "#BBBBBB", "#E69F00"
ORF_EDGE = "#CC79A7"
THIN, NNBSP, MINUS = chr(0x2009), chr(0x202F), chr(0x2212)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def nov() -> ModuleType:
    """The module under test (imported per test so each test fails on its own)."""
    return importlib.import_module("gippa.viz.novel_orfs")


# --- fixtures and helpers -----------------------------------------------------------------------


def orf(
    orf_id: str,
    *,
    cls: str = "intergenic",
    contig: str = "chr1",
    strand: str = "+",
    start_longest: int | None = 5_001,
    end: int = 5_300,
    start_best_ostir: int | None = None,
    called: bool = True,
    rank: int | None = 1,
    class_fdr: str | None = "0.004",
    n_peptides: int = 2,
    peptides: str | None = "PEPTWDEK;SAMPQER",
    best_q: float | None = 0.0001,
    **kw: object,
) -> dict[str, object]:
    """One ``novel_orfs.tsv`` row; coordinates 1-based, in transcription direction."""
    row: dict[str, object] = {
        "orf_id": orf_id,
        "class": cls,
        "contig": contig,
        "strand": strand,
        "frame": 1 if strand == "+" else -1,
        "start_longest": start_longest,
        "start_best_ostir": start_longest if start_best_ostir is None else start_best_ostir,
        "end": end,
        "length_aa_longest": 99,
        "length_aa_best": 99,
        "start_codon_best": "ATG",
        "n_peptides": n_peptides,
        "peptides": peptides,
        "best_peptide_q": best_q,
        "ostir_tir": 1000.0,
        "promoter_score": 0.9,
        "promoter_supported": True,
        "terminator_downstream": True,
        "terminator_conf": 80.0,
        "terminator_distance_nt": 30,
        "rbs_supported": True,
        "class_fdr": class_fdr,
        "called": called,
        "rank": rank,
        "overlapping_gene": None,
        "flags": None,
    }
    row.update(kw)
    return row


def _out(tmp: Path, rows: Sequence[Mapping[str, object]], **kw: Any) -> Path:
    return vl.build_out(tmp, novel_orfs=rows, **kw)


def _run(out: Path, cfg: GippaConfig, nov: ModuleType) -> list[Instance]:
    assert nov is not None
    result = base.registered()["V13"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    return list(result)


def _main(insts: Sequence[Instance]) -> Instance:
    (m,) = [i for i in insts if i.suffix is None]
    return m


def _views(insts: Sequence[Instance]) -> list[Instance]:
    return [i for i in insts if i.suffix is not None and i.suffix.startswith("locus_")]


def _index(insts: Sequence[Instance]) -> Instance:
    (ix,) = [i for i in insts if i.suffix == "locus"]
    return ix


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _agg(fig: Figure) -> Any:
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    fig.draw(renderer)
    return renderer


def _collapse(s: str) -> str:
    return " ".join(s.split())


def _hex(c: Any) -> str | None:
    try:
        return mcolors.to_hex(c).upper()
    except (ValueError, TypeError):
        return None


def _texts(fig: Figure) -> list[Text]:
    return [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]


def _counts(inst: Instance) -> list[tuple[str, str, int]]:
    d = inst.data
    return [(str(r["class"]), str(r["status"]), int(r["n"])) for _, r in d.iterrows() if r["n"]]


def _bar_lengths(ax: Axes, colour: str) -> list[float]:
    """Lengths (along the count axis) of the non-empty bars of ``colour`` on ``ax``; the
    orientation is the one in which every bar has the same thickness."""
    status = {CALLED, NOT_CALLED, UNEST}
    rects = [
        p for p in ax.patches
        if isinstance(p, Rectangle) and _hex(p.get_facecolor()) in status
        and p.get_width() > 0 and p.get_height() > 0
    ]  # fmt: skip
    vertical = len({round(p.get_width(), 6) for p in rects}) == 1 and (
        len({round(p.get_height(), 6) for p in rects}) > 1
        or len({round(p.get_x(), 6) for p in rects}) > 1
    )
    mine = [p for p in rects if _hex(p.get_facecolor()) == colour]
    return sorted(round(p.get_height() if vertical else p.get_width(), 6) for p in mine)


def _example1_rows() -> list[dict[str, object]]:
    return [
        orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, class_fdr="0.004", rank=1),
        orf("NOVEL_chr1_3001_p", start_longest=3_001, end=3_300, class_fdr="0.008", rank=2),
        orf("NOVEL_chr1_7001_p", start_longest=7_001, end=7_300, class_fdr="0.2", called=False,
            rank=3),
        orf("NOVEL_chr1_9001_m", cls="antisense", strand="-", start_longest=9_300, end=9_001,
            class_fdr="unestimated", called=False, rank=4, flags="class_fdr_unestimated"),
        orf("NOVEL_chr1_12001_m", cls="antisense", strand="-", start_longest=12_300,
            end=12_001, class_fdr="unestimated", called=False, rank=5),
    ]  # fmt: skip


# --- registration -------------------------------------------------------------------------------


def test_v13_is_registered_with_novel_orfs_as_its_only_hard_need(nov: ModuleType) -> None:
    spec = base.registered()["V13"]
    assert spec.name == "novel_orfs"
    assert spec.needs == ("novel_orfs",)


def test_absent_table_is_a_declared_skip(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    out = vl.build_out(tmp_path)
    result = base.registered()["V13"].fn(RunView.from_outdir(out), cfg)
    assert result == Skip("input_absent:novel_orfs")


def test_module_is_pure_and_uses_no_pyplot(nov: ModuleType) -> None:
    assert "pyplot" not in inspect.getsource(nov)


# --- example 1: panel (a) -----------------------------------------------------------------------


def test_example1_counts(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    main = _main(_run(_out(tmp_path, _example1_rows()), cfg, nov))
    assert list(main.data.columns) == COUNT_COLUMNS
    assert _counts(main) == [
        ("intergenic", "called", 2),
        ("intergenic", "not_called", 1),
        ("antisense", "unestimated", 2),
    ]


def test_example1_class_fdr_max_called(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    d = _main(_run(_out(tmp_path, _example1_rows()), cfg, nov)).data
    called = d[(d["class"] == "intergenic") & (d["status"] == "called")]
    assert float(called["class_fdr_max_called"].iloc[0]) == pytest.approx(0.008)
    anti = d[(d["class"] == "antisense") & (d["status"] == "unestimated")]
    assert pd.isna(anti["class_fdr_max_called"].iloc[0])


def test_example1_title(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    main = _main(_run(_out(tmp_path, _example1_rows()), cfg, nov))
    want = "2 novel ORFs called (class FDR ≤ 0.008), 2 unestimated"
    assert _collapse(main.title) == want
    assert any(want in _collapse(t.get_text()) for t in _texts(main.figure))


def test_example1_bars_are_stacked_in_the_status_colours(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    ax = _axes(_main(_run(_out(tmp_path, _example1_rows()), cfg, nov)).figure, "counts")
    assert _bar_lengths(ax, CALLED) == [2.0]
    assert _bar_lengths(ax, NOT_CALLED) == [1.0]
    assert _bar_lengths(ax, UNEST) == [2.0]


def test_example1_bar_labels_give_the_fdr_range_or_unestimated(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    ax = _axes(_main(_run(_out(tmp_path, _example1_rows()), cfg, nov)).figure, "counts")
    texts = [_collapse(t.get_text()) for t in ax.texts]
    assert "FDR ≤ 0.008" in texts
    assert "unestimated" in texts


def test_rows_sorted_by_class_enum_then_status(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rng = random.Random(7)
    rows = []
    for i in range(30):
        cls = CLASSES[i % len(CLASSES)]
        fdr = rng.choice(["0.01", "0.3", "unestimated", None])
        called = fdr == "0.01" and i % 2 == 0
        rows.append(orf(f"NOVEL_chr1_{1001 + 1000 * i}_p", cls=cls, start_longest=1001 + 1000 * i,
                        end=1300 + 1000 * i, class_fdr=fdr, called=called, rank=i + 1))  # fmt: skip
    rng.shuffle(rows)
    d = _main(_run(_out(tmp_path, rows), cfg, nov)).data
    keys = [
        (CLASSES.index(c), STATUSES.index(s)) for c, s in zip(d["class"], d["status"], strict=True)
    ]
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys)
    assert int(d["n"].sum()) == 30


def test_category_ticks_follow_the_class_enum_order(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", cls="sorf", start_longest=1001, end=1100),
        orf("NOVEL_chr1_2001_p", cls="nterm_extension", start_longest=2001, end=2300),
        orf("NOVEL_chr1_4001_p", cls="antisense", start_longest=4001, end=4300),
    ]
    fig = _main(_run(_out(tmp_path, rows), cfg, nov)).figure
    ax = _axes(fig, "counts")
    r = _agg(fig)
    xs = [(t.get_window_extent(r).x0, t.get_text()) for t in ax.get_xticklabels()]
    ys = [(-t.get_window_extent(r).y0, t.get_text()) for t in ax.get_yticklabels()]
    for labels in (xs, ys):
        shown = [txt for _, txt in sorted(labels) if txt in CLASSES]
        if shown:
            assert {"sorf", "nterm_extension", "antisense"} <= set(shown)
            assert shown == sorted(shown, key=CLASSES.index)
            return
    pytest.fail("no class tick labels on the counts axes")


# --- the literal `unestimated` and the flag ------------------------------------------------------


@pytest.mark.parametrize(
    ("fdr", "flags"),
    [("unestimated", None), (None, None), ("0.05", "class_fdr_unestimated"),
     ("unestimated", "partial;class_fdr_unestimated")],
)  # fmt: skip
def test_unestimated_class_fdr_is_status_unestimated(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, fdr: str | None, flags: str | None
) -> None:
    rows = [orf("NOVEL_chr1_1001_p", class_fdr=fdr, called=False, flags=flags)]
    insts = _run(_out(tmp_path, rows), cfg, nov)
    assert _counts(_main(insts)) == [("intergenic", "unestimated", 1)]
    assert _views(insts) == []


def test_literal_unestimated_is_read_as_na(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", class_fdr="0.004"),
        orf("NOVEL_chr1_3001_p", start_longest=3001, end=3300, class_fdr="unestimated",
            called=False, rank=2),
    ]  # fmt: skip
    d = _main(_run(_out(tmp_path, rows), cfg, nov)).data
    vals = pd.to_numeric(d["class_fdr_max_called"], errors="raise")
    assert vals.dtype.kind == "f"
    assert sorted(vals.dropna().tolist()) == [pytest.approx(0.004)]


def test_called_row_with_unestimated_fdr_titles_class_fdr_unestimated(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [orf("NOVEL_chr1_5001_p", class_fdr="unestimated", n_peptides=3, best_q=0.002)]
    (view,) = _views(_run(_out(tmp_path, rows), cfg, nov))
    assert _collapse(view.title) == (
        "NOVEL_chr1_5001_p (intergenic): 3 peptides, best q 0.002, class FDR unestimated"
    )


# --- degenerate cases ---------------------------------------------------------------------------


def test_example4_header_only_is_the_empty_state(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    out = _out(tmp_path / "run", [])
    insts = _run(out, cfg, nov)
    main = _main(insts)
    assert main.title == "No novel ORF candidates in this run"
    assert list(main.data.columns) == COUNT_COLUMNS
    assert len(main.data) == 0
    assert _views(insts) == []
    for inst in insts:
        assert inst.suffix in (None, "locus")
        if inst.suffix == "locus":
            assert len(inst.data) == 0
        paths = write_instance(inst, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100)
        assert not [p for p in paths if p.name.startswith("V13_novel_orfs_locus_")]


def test_no_called_but_unestimated_is_panel_a_only(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", class_fdr="unestimated", called=False, rank=1),
        orf("NOVEL_chr1_3001_p", start_longest=3001, end=3300, class_fdr=None, called=False,
            rank=2),
    ]  # fmt: skip
    insts = _run(_out(tmp_path, rows), cfg, nov)
    main = _main(insts)
    assert _collapse(main.title) == "No novel ORF called; 2 candidates unestimated"
    assert _counts(main) == [("intergenic", "unestimated", 2)]
    assert _views(insts) == []
    assert all(len(i.data) == 0 for i in insts if i.suffix == "locus")


def test_fewer_than_20_called_gives_that_many_views(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf(f"NOVEL_chr1_{1001 + 5000 * i}_p", start_longest=1001 + 5000 * i, end=1300 + 5000 * i,
            rank=i + 1)
        for i in range(3)
    ]  # fmt: skip
    insts = _run(_out(tmp_path, rows), cfg, nov)
    assert [v.suffix for v in _views(insts)] == [f"locus_{r['orf_id']}" for r in rows]
    ix = _index(insts)
    assert list(ix.data.columns) == INDEX_COLUMNS
    assert ix.data["rank"].tolist() == [1, 2, 3]


def test_evidence_columns_na_leave_panel_a_unchanged_and_views_render(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    """PR #180 not merged: the M7-09 evidence columns are all NA."""
    evidence = (
        "promoter_score", "promoter_supported", "terminator_downstream", "terminator_conf",
        "terminator_distance_nt", "rbs_supported", "ostir_tir",
    )  # fmt: skip
    full = _example1_rows()
    bare = [{**r, **dict.fromkeys(evidence)} for r in full]
    a = _run(_out(tmp_path / "a", full), cfg, nov)
    b = _run(_out(tmp_path / "b", bare), cfg, nov)
    pd.testing.assert_frame_equal(_main(a).data, _main(b).data)
    assert _main(a).title == _main(b).title
    assert [v.suffix for v in _views(a)] == [v.suffix for v in _views(b)]
    assert len(_views(b)) == 2


# --- example 2 and decision 1: top 20 by rank among called rows -----------------------------------


def _many(n_called: int, *, uncalled: bool = False) -> list[dict[str, object]]:
    rows = []
    for i in range(n_called):
        s = 1001 + 3000 * i
        rank = 2 * (i + 1) if uncalled else i + 1
        rows.append(orf(f"NOVEL_chr1_{s}_p", start_longest=s, end=s + 299, rank=rank))
        if uncalled:
            u = s + 1500
            rows.append(orf(f"NOVEL_chr1_{u}_p", start_longest=u, end=u + 299, rank=rank - 1,
                            called=False, class_fdr="0.3"))  # fmt: skip
    random.Random(3).shuffle(rows)
    return rows


def test_example2_25_called_gives_20_views_ranks_1_to_20(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = _many(25)
    insts = _run(_out(tmp_path, rows), cfg, nov)
    by_rank = {int(str(r["rank"])): str(r["orf_id"]) for r in rows}
    assert [v.suffix for v in _views(insts)] == [f"locus_{by_rank[k]}" for k in range(1, 21)]
    ix = _index(insts)
    assert list(ix.data.columns) == INDEX_COLUMNS
    assert ix.data["rank"].tolist() == list(range(1, 21))
    assert ix.data["orf_id"].tolist() == [by_rank[k] for k in range(1, 21)]
    assert ix.data["rendered"].astype(bool).tolist() == [True] * 20


def test_uncalled_rows_never_get_a_view_even_with_smaller_rank(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = _many(22, uncalled=True)
    insts = _run(_out(tmp_path, rows), cfg, nov)
    called = {str(r["orf_id"]): int(str(r["rank"])) for r in rows if r["called"]}
    ids = [v.suffix.removeprefix("locus_") for v in _views(insts) if v.suffix]
    assert len(ids) == 20
    assert set(ids) <= set(called)
    assert [called[i] for i in ids] == list(range(2, 42, 2))
    assert _index(insts).data["rank"].tolist() == list(range(2, 42, 2))


def test_example2_index_file_name(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    ix = _index(_run(_out(tmp_path / "run", _many(2)), cfg, nov))
    paths = write_instance(ix, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100)
    names = sorted(p.relative_to(tmp_path / "plots").as_posix() for p in paths)
    assert names == ["V13_novel_orfs_locus.svg", "data/V13_novel_orfs_locus.tsv"]


# --- the V3 renderer (no second implementation) --------------------------------------------------


class _Spy:
    def __init__(self) -> None:
        self.loci: list[locus_mod.Locus] = []
        self.real = locus_mod.render_locus

    def __call__(self, lo: locus_mod.Locus, *args: Any, **kw: Any) -> Instance:
        self.loci.append(lo)
        return self.real(lo, *args, **kw)


@pytest.fixture
def spy(nov: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Iterator[_Spy]:
    s = _Spy()
    monkeypatch.setattr(locus_mod, "render_locus", s)
    if hasattr(nov, "render_locus"):
        monkeypatch.setattr(nov, "render_locus", s)
    yield s


def test_every_view_is_drawn_by_render_locus(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, spy: _Spy
) -> None:
    insts = _run(_out(tmp_path, _many(4)), cfg, nov)
    assert len(spy.loci) == len(_views(insts)) == 4
    assert all(isinstance(lo, locus_mod.Locus) and lo.kind == "novel_orf" for lo in spy.loci)
    assert [f"locus_{lo.locus_id}" for lo in spy.loci] == [v.suffix for v in _views(insts)]


def test_views_have_the_v3_layout_plus_a_novel_orf_track(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [orf("NOVEL_chr1_5001_p", start_longest=5_001, end=5_300)]
    (view,) = _views(_run(_out(tmp_path, rows), cfg, nov))
    assert list(view.data.columns) == LOCUS_COLUMNS
    assert "start" in view.coords
    assert {"pileup", "genes"} <= {ax.get_label() for ax in view.figure.axes}
    nv = view.data[view.data["track"] == "novel_orf"]
    assert [(r.id, r.contig, int(r.start), int(r.end), r.strand) for r in nv.itertuples()] == [
        ("NOVEL_chr1_5001_p", "chr1", 5_000, 5_300, "+")
    ]


def test_locus_title(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    rows = [orf("NOVEL_chr1_5001_p", n_peptides=3, best_q=0.0001, class_fdr="0.004")]
    (view,) = _views(_run(_out(tmp_path, rows), cfg, nov))
    want = "NOVEL_chr1_5001_p (intergenic): 3 peptides, best q 1e-04, class FDR 0.004"
    assert _collapse(view.title) == want
    assert want in _collapse(" ".join(t.get_text() for t in _texts(view.figure)))


# --- example 3: minus strand, transcription direction, circular window ---------------------------


def _example3_out(tmp: Path) -> Path:
    rows = [
        orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500, end=1_201,
            start_best_ostir=1_500, peptides="PEPTWDEK;SAMPQER")
    ]  # fmt: skip
    hits = [
        vl.hit("PEPTWDEK", "chr1", 1_301, 1_324, "-"),
        vl.hit("SAMPQER", "chr1", 1_401, 1_421, "-"),
        vl.hit("GGGGGGGK", "chr1", 2_501, 2_524, "+"),
    ]
    genes = [vl.gene("g1", "chr1", 1_001, 1_700, strand="+", log2=20.0)]
    return _out(tmp, rows, hits=hits, genes=genes)


def test_example3_span_is_zero_based_half_open(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, spy: _Spy
) -> None:
    (view,) = _views(_run(_example3_out(tmp_path), cfg, nov))
    ((lo),) = spy.loci
    assert (lo.contig, lo.start0, lo.end0) == ("chr1", 1_200, 1_500)
    nv = view.data[view.data["track"] == "novel_orf"]
    assert [(int(r.start), int(r.end), r.strand) for r in nv.itertuples()] == [(1_200, 1_500, "-")]


def test_example3_data_tsv_is_one_based(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    (view,) = _views(_run(_example3_out(tmp_path / "run"), cfg, nov))
    paths = write_instance(view, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100)
    (tsv,) = [p for p in paths if p.suffix == ".tsv"]
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == LOCUS_COLUMNS
    nv = [ln.split("\t") for ln in lines[1:] if ln.startswith("novel_orf\t")]
    assert [(r[1], r[3], r[4], r[5]) for r in nv] == [("NOVEL_chr1_1201_m", "1201", "1500", "-")]


def test_example3_window_is_unwrapped_below_zero(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    (view,) = _views(_run(_example3_out(tmp_path), cfg, nov))
    assert _axes(view.figure, "genes").get_xlim() == pytest.approx((-800, 3_500))


def _tick_int(label: str) -> int:
    """The bp position of a kb tick label."""
    s = label.strip().replace(",", "").replace(THIN, "").replace(NNBSP, "")
    return round(float(s.replace(" ", "")) * 1000)


def test_example3_ticks_are_labelled_by_wrapped_position(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    (view,) = _views(_run(_example3_out(tmp_path), cfg, nov))
    ax = _axes(view.figure, "genes")
    fmt = ax.xaxis.get_major_formatter()
    assert _tick_int(fmt(-500.0, 0)) == 99_500
    assert _tick_int(fmt(1_000.0, 1)) == 1_000
    _agg(view.figure)
    for t in ax.get_xticklabels():
        txt = t.get_text().strip()
        if txt:
            assert not txt.startswith(("-", MINUS)), txt
            assert 0 <= _tick_int(txt) < 100_000


def test_wrapped_ticks_past_the_length(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    """A ``+`` ORF crossing the origin (``end < start_longest``) is unwrapped past 100 kb."""
    rows = [orf("NOVEL_chr1_99901_p", start_longest=99_901, end=200)]
    (view,) = _views(_run(_out(tmp_path, rows), cfg, nov))
    ax = _axes(view.figure, "genes")
    assert ax.get_xlim() == pytest.approx((97_900, 102_200))
    nv = view.data[view.data["track"] == "novel_orf"]
    assert [(int(r.start), int(r.end)) for r in nv.itertuples()] == [(99_900, 100_200)]
    assert _tick_int(ax.xaxis.get_major_formatter()(100_500.0, 0)) == 500


def test_minus_orf_crossing_the_origin(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, spy: _Spy
) -> None:
    """On ``-`` the row crosses the origin when ``start_longest < end``."""
    rows = [orf("NOVEL_chr1_99851_m", strand="-", start_longest=150, end=99_851)]
    _run(_out(tmp_path, rows), cfg, nov)
    ((lo),) = spy.loci
    assert (lo.start0, lo.end0) == (99_850, 100_150)


@pytest.mark.parametrize(
    ("start", "end", "xlim"),
    [(301, 600, (0, 2_600)), (29_701, 29_900, (27_700, 30_000))],
)
def test_linear_window_clips(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, start: int, end: int, xlim: tuple[int, int]
) -> None:
    rows = [orf(f"NOVEL_lin1_{start}_p", contig="lin1", start_longest=start, end=end)]
    (view,) = _views(_run(_out(tmp_path, rows), cfg, nov))
    assert _axes(view.figure, "genes").get_xlim() == pytest.approx(xlim)


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    strand=st.sampled_from(["+", "-"]),
    lo=st.integers(min_value=3_000, max_value=90_000),
    codons=st.integers(min_value=10, max_value=400),
)
def test_property_span_is_min_minus_1_to_max(
    cfg: GippaConfig, nov: ModuleType, strand: str, lo: int, codons: int
) -> None:
    hi = lo + 3 * codons - 1
    first, last = (lo, hi) if strand == "+" else (hi, lo)
    rows = [orf(f"NOVEL_chr1_{lo}_{'p' if strand == '+' else 'm'}", strand=strand,
                start_longest=first, end=last)]  # fmt: skip
    with tempfile.TemporaryDirectory() as tmp, pytest.MonkeyPatch.context() as mp:
        s = _Spy()
        mp.setattr(locus_mod, "render_locus", s)
        if hasattr(nov, "render_locus"):
            mp.setattr(nov, "render_locus", s)
        (view,) = _views(_run(_out(Path(tmp), rows), cfg, nov))
    ((locus),) = s.loci
    assert (locus.start0, locus.end0) == (lo - 1, hi)
    nv = view.data[view.data["track"] == "novel_orf"]
    assert [(int(r.start), int(r.end), r.strand) for r in nv.itertuples()] == [(lo - 1, hi, strand)]


# --- start_longest NA: the span starts at the first peptide ---------------------------------------


@pytest.mark.parametrize(
    ("strand", "end", "hits", "span"),
    [
        ("+", 5_300, [(4_801, 4_824), (5_001, 5_021)], (4_800, 5_300)),
        ("-", 1_201, [(1_401, 1_424), (1_601, 1_621)], (1_200, 1_621)),
    ],
)
def test_no_start_codon_span_starts_at_the_first_peptide(
    tmp_path: Path,
    cfg: GippaConfig,
    nov: ModuleType,
    spy: _Spy,
    strand: str,
    end: int,
    hits: list[tuple[int, int]],
    span: tuple[int, int],
) -> None:
    oid = f"NOVEL_chr1_{span[0] + 1}_{'p' if strand == '+' else 'm'}"
    rows = [orf(oid, strand=strand, start_longest=None, end=end, peptides="PEPTWDEK;SAMPQER",
                flags="no_start_codon")]  # fmt: skip
    peps = ["PEPTWDEK", "SAMPQER"]
    placed = [vl.hit(p, "chr1", s, e, strand) for p, (s, e) in zip(peps, hits, strict=True)]
    placed.append(vl.hit("PEPTWDEK", "lin1", 101, 124, strand))  # another contig: ignored
    _run(_out(tmp_path, rows, hits=placed), cfg, nov)
    ((lo),) = spy.loci
    assert (lo.contig, lo.start0, lo.end0) == ("chr1", *span)


# --- drawing: the ORF arrow, the OSTIR tick, the peptides -----------------------------------------


def _patches(fig: Figure) -> list[tuple[Axes, Patch]]:
    return [(ax, p) for ax in fig.axes for p in ax.patches]


def test_example3_orf_arrow_outline_points_left_above_the_genes(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    fig = _views(_run(_example3_out(tmp_path), cfg, nov))[0].figure
    genes_ax = _axes(fig, "genes")
    _agg(fig)
    orfs = [(ax, p) for ax, p in _patches(fig) if _hex(p.get_edgecolor()) == ORF_EDGE]
    assert len(orfs) == 1
    _, p = orfs[0]
    disp = p.get_transform().transform(p.get_path().vertices)
    data_x = genes_ax.transData.inverted().transform(disp)[:, 0]
    width = 4_300
    assert data_x.min() == pytest.approx(1_200, abs=0.002 * width)
    assert data_x.max() == pytest.approx(1_500, abs=0.002 * width)
    at_min = disp[np.abs(data_x - data_x.min()) <= 1e-6 * width]
    at_max = disp[np.abs(data_x - data_x.max()) <= 1e-6 * width]
    assert np.ptp(at_min[:, 1]) < np.ptp(at_max[:, 1])  # the head is at the low end on '-'
    gene_tops = [
        q.get_transform().transform(q.get_path().vertices)[:, 1].max()
        for q in genes_ax.patches
        if q is not p and q.get_facecolor()[3] > 0.5 and _hex(q.get_edgecolor()) != ORF_EDGE
    ]
    assert gene_tops
    assert disp[:, 1].min() >= max(gene_tops) - 0.5


def test_ostir_start_differing_gets_a_tick(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [orf("NOVEL_chr1_5001_p", start_longest=5_001, start_best_ostir=5_031, end=5_300)]
    fig = _views(_run(_out(tmp_path, rows), cfg, nov))[0].figure
    ax = _axes(fig, "genes")

    def at(x0: float) -> list[Line2D]:
        out = []
        for ln in ax.get_lines():
            x = np.asarray(ln.get_xdata(), dtype=float)
            if len(x) and np.all(np.abs(x - x0) <= 1.0):
                out.append(ln)
        return out

    assert at(5_030)


def test_orf_peptides_are_highlighted_in_the_pileup(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    def highlighted(peptides: str | None) -> bool:
        rows = [orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500,
                    end=1_201, peptides=peptides)]  # fmt: skip
        hits = [vl.hit("PEPTWDEK", "chr1", 1_301, 1_324, "-")]
        sub = tmp_path / str(peptides)
        fig = _views(_run(_out(sub, rows, hits=hits), cfg, nov))[0].figure
        pile = _axes(fig, "pileup")
        cols: list[Any] = [ln.get_color() for ln in pile.get_lines()]
        cols += [c for p in pile.patches for c in (p.get_facecolor(), p.get_edgecolor())]
        for c in pile.collections:
            if isinstance(c, Collection):
                cols += list(c.get_facecolors()) + list(c.get_edgecolors())
        return ORF_EDGE in {_hex(c) for c in cols}

    assert highlighted("PEPTWDEK")
    assert not highlighted(None)


# --- file names, text, layout, determinism ------------------------------------------------------


def test_locus_files_are_named_by_sanitised_orf_id(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [orf("NOVEL_chr/1 x_5001_p")]
    (view,) = _views(_run(_out(tmp_path / "run", rows), cfg, nov))
    assert view.suffix == "locus_NOVEL_chr/1 x_5001_p"
    paths = write_instance(
        view, "V13", "novel_orfs", tmp_path / "plots", ["svg", "png", "pdf"], 100
    )
    assert sorted(p.relative_to(tmp_path / "plots").as_posix() for p in paths) == [
        "V13_novel_orfs_locus_NOVEL_chr_1_x_5001_p.pdf",
        "V13_novel_orfs_locus_NOVEL_chr_1_x_5001_p.png",
        "V13_novel_orfs_locus_NOVEL_chr_1_x_5001_p.svg",
        "data/V13_novel_orfs_locus_NOVEL_chr_1_x_5001_p.tsv",
    ]


def test_main_data_tsv_name_and_header(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    main = _main(_run(_out(tmp_path / "run", _example1_rows()), cfg, nov))
    paths = write_instance(main, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100)
    tsv = tmp_path / "plots" / "data" / "V13_novel_orfs.tsv"
    assert tsv in paths
    assert tsv.read_text(encoding="utf-8").splitlines()[0].split("\t") == COUNT_COLUMNS


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    insts = _run(_example3_out(tmp_path / "a"), cfg, nov)
    insts += _run(_out(tmp_path / "b", _example1_rows()), cfg, nov)
    for inst in insts:
        _agg(inst.figure)
        sizes = [t.get_fontsize() for t in _texts(inst.figure)]
        assert min(sizes, default=6.0) >= 6, inst.suffix
        assert sizes or inst.suffix == "locus", inst.suffix


def test_main_texts_inside_the_figure_and_not_overlapping(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    fig = _main(_run(_out(tmp_path, _example1_rows()), cfg, nov)).figure
    r = _agg(fig)
    box = fig.bbox
    texts = _texts(fig)
    exts = [t.get_window_extent(r) for t in texts]
    for t, e in zip(texts, exts, strict=True):
        assert e.x0 >= box.x0 - 1 and e.x1 <= box.x1 + 1, t.get_text()
        assert e.y0 >= box.y0 - 1 and e.y1 <= box.y1 + 1, t.get_text()
    ax = _axes(fig, "counts")
    labels = [t for t in ax.texts if t.get_text().strip()]
    lab_ext = [t.get_window_extent(r) for t in labels]
    for i, j in itertools.combinations(range(len(labels)), 2):
        assert not lab_ext[i].overlaps(lab_ext[j]), (labels[i].get_text(), labels[j].get_text())
    title = [
        e
        for t, e in zip(texts, exts, strict=True)
        if "novel ORFs called" in _collapse(t.get_text())
    ]
    assert title
    assert not any(title[0].overlaps(e) for e in lab_ext)


def test_files_are_byte_identical_across_renders(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    out = _example3_out(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        files: dict[Path, bytes] = {}
        for inst in _run(out, cfg, nov):
            for p in write_instance(
                inst, "V13", "novel_orfs", tmp_path / name, ["svg", "png"], 100
            ):
                files[p.relative_to(tmp_path / name)] = p.read_bytes()
        written.append(files)
    assert Path("V13_novel_orfs.svg") in written[0]
    assert Path("V13_novel_orfs_locus_NOVEL_chr1_1201_m.svg") in written[0]
    assert Path("data/V13_novel_orfs_locus.tsv") in written[0]
    assert written[0] == written[1]


def test_render_leaves_no_pyplot_figure(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    insts = _run(_example3_out(tmp_path), cfg, nov)
    assert all(isinstance(i.figure, Figure) for i in insts)
    assert plt.get_fignums() == before


def test_titles_never_say_express(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    for inst in _run(_example3_out(tmp_path), cfg, nov):
        assert "express" not in inst.title.lower()
