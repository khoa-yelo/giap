"""V13 review round (PR #219; M9-12 lead decisions 11-14).

- 11: the novel ORF arrow has its own band: after an Agg draw its patch's window extent meets
  no gene-label ``Text`` and no ``+`` promoter / terminator mark ``Line2D``.
- 12: the ``start_best_ostir`` tick is pinned exactly (< 0.5 bp): ``+`` at the codon's first base,
  ``-`` with the add-back, ``-`` across the origin shifted by the length, none when OSTIR equals
  the drawn start.
- 13: every locus view carries a one-entry key for the arrow, the tick and the purple pileup.
- 14: "none called", the singular title, the title's FDR bound over every class, NA rank last.

Fixtures are 1-based ``novel_orfs.tsv`` rows on the 100 kb circular ``chr1`` of
``helpers.viz_locus``; BED rows are 0-based half-open.
"""

from __future__ import annotations

import importlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pytest
from helpers import viz_locus as vl
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

ORF_EDGE = "#CC79A7"
PROMOTER, TERMINATOR = "#009E73", "#D55E00"
EXACT = 0.5  # bp


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def nov() -> ModuleType:
    return importlib.import_module("gippa.viz.novel_orfs")


def orf(
    orf_id: str,
    *,
    cls: str = "intergenic",
    strand: str = "+",
    start_longest: int | None = 5_001,
    end: int = 5_300,
    start_best_ostir: int | None = None,
    called: bool = True,
    rank: int | None = 1,
    class_fdr: str | None = "0.004",
    peptides: str | None = "PEPTWDEK;SAMPQER",
    **kw: object,
) -> dict[str, object]:
    """One ``novel_orfs.tsv`` row (the ``class`` column keyed as in the schema)."""
    row: dict[str, object] = {
        "orf_id": orf_id,
        "class": cls,
        "contig": "chr1",
        "strand": strand,
        "frame": 1 if strand == "+" else -1,
        "start_longest": start_longest,
        "start_best_ostir": start_longest if start_best_ostir is None else start_best_ostir,
        "end": end,
        "length_aa_longest": 99,
        "length_aa_best": 99,
        "start_codon_best": "ATG",
        "n_peptides": 2,
        "peptides": peptides,
        "best_peptide_q": 0.0001,
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


def _run(out: Path, cfg: GippaConfig, nov: ModuleType) -> list[Instance]:
    assert nov is not None
    result = base.registered()["V13"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    return list(result)


def _rows_run(
    tmp: Path, cfg: GippaConfig, nov: ModuleType, rows: Sequence[Mapping[str, object]], **kw: Any
) -> list[Instance]:
    return _run(vl.build_out(tmp, novel_orfs=rows, **kw), cfg, nov)


def _main(insts: Sequence[Instance]) -> Instance:
    (m,) = [i for i in insts if i.suffix is None]
    return m


def _views(insts: Sequence[Instance]) -> list[Instance]:
    return [i for i in insts if i.suffix is not None and i.suffix.startswith("locus_")]


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _agg(fig: Figure) -> Any:
    canvas = FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    fig.draw(renderer)
    return renderer


def _hex(c: Any) -> str | None:
    try:
        return mcolors.to_hex(c).upper()
    except (ValueError, TypeError):
        return None


def _collapse(s: str) -> str:
    return " ".join(s.split())


def _orf_patches(fig: Figure) -> list[Patch]:
    return [p for ax in fig.axes for p in ax.patches if _hex(p.get_edgecolor()) == ORF_EDGE]


# --- decision 11: the ORF band is clear of gene labels and + marks ------------------------------

GENE_TAGS = ("gLeft", "gUnder", "gRight")


def _band_out(tmp: Path, strand: str) -> Path:
    """An ORF over the middle of ``gUnder`` (its label at the gene centre, 40,000), with ``+``
    promoter and terminator marks inside the ORF span and flanking labelled genes."""
    if strand == "+":
        row = orf("NOVEL_chr1_39501_p", cls="cds_outframe", start_longest=39_501, end=40_499,
                  start_best_ostir=39_531)  # fmt: skip
    else:
        row = orf("NOVEL_chr1_39501_m", cls="antisense", strand="-", start_longest=40_499,
                  end=39_501, start_best_ostir=40_399)  # fmt: skip
    genes = [
        vl.gene("gLeft", "chr1", 38_001, 38_900, strand="+", log2=21.0),
        vl.gene("gUnder", "chr1", 39_001, 41_000, strand="-", log2=22.0),
        vl.gene("gRight", "chr1", 41_101, 41_900, strand="+", log2=23.0),
    ]
    promoters = [("chr1", 39_700, 39_750, "P1", 0.9, "+"), ("chr1", 40_300, 40_350, "P2", 0.9, "+")]
    terminators = [("chr1", 39_900, 39_940, "T1", 80.0, "+"), ("chr1", 40_100, 40_140, "T2", 80.0,
                                                                  "+")]  # fmt: skip
    hits = [vl.hit("PEPTWDEK", "chr1", 39_801, 39_824, strand)]
    return vl.build_out(tmp, novel_orfs=[row], genes=genes, promoters=promoters,
                        terminators=terminators, hits=hits)  # fmt: skip


@pytest.mark.parametrize("strand", ["+", "-"])
def test_orf_arrow_meets_no_gene_label(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, strand: str
) -> None:
    (view,) = _views(_run(_band_out(tmp_path, strand), cfg, nov))
    fig = view.figure
    r = _agg(fig)
    (arrow,) = _orf_patches(fig)
    box = arrow.get_window_extent(r)
    labels = [t for ax in fig.axes for t in ax.texts if t.get_text().strip() in GENE_TAGS]
    assert {t.get_text().strip() for t in labels} == set(GENE_TAGS)
    hit = [t.get_text() for t in labels if t.get_window_extent(r).overlaps(box)]
    assert hit == [], f"ORF arrow covers gene labels {hit}"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_orf_arrow_meets_no_plus_strand_mark(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, strand: str
) -> None:
    (view,) = _views(_run(_band_out(tmp_path, strand), cfg, nov))
    fig = view.figure
    r = _agg(fig)
    (arrow,) = _orf_patches(fig)
    box = arrow.get_window_extent(r)
    marks = [
        ln for ax in fig.axes for ln in ax.get_lines() if _hex(ln.get_color()) in (PROMOTER,
                                                                                   TERMINATOR)
    ]  # fmt: skip
    assert len({_hex(ln.get_color()) for ln in marks}) == 2  # both kinds drawn
    hit = [
        (_hex(ln.get_color()), tuple(ln.get_xdata()))
        for ln in marks
        if ln.get_window_extent(r).overlaps(box)
    ]
    assert hit == [], f"ORF arrow covers marks {hit}"


# --- decision 12: exact OSTIR tick positions ---------------------------------------------------


def _ticks(fig: Figure) -> list[float]:
    """x of every vertical ``Line2D`` in the ORF colour on the ``"genes"`` axes."""
    out = []
    for ln in _axes(fig, "genes").get_lines():
        x = np.asarray(ln.get_xdata(), dtype=float)
        if len(x) >= 2 and _hex(ln.get_color()) == ORF_EDGE and np.ptp(x) == 0:
            out.append(float(x[0]))
    return out


@pytest.mark.parametrize(
    ("row", "want"),
    [
        # + : the codon's first base, 1-based 5,031 -> 0-based 5,030
        (orf("NOVEL_chr1_5001_p", start_longest=5_001, start_best_ostir=5_031, end=5_300), 5_030),
        # - : 1-based OSTIR 1,401 on example 3 -> 0-based 1,400, +1 add-back -> 1,401
        (orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500,
             start_best_ostir=1_401, end=1_201), 1_401),
        # - crossing the origin: span [99,850, 100,150); OSTIR 90 -> 89 + 1 + 100,000
        (orf("NOVEL_chr1_99851_m", cls="antisense", strand="-", start_longest=150,
             start_best_ostir=90, end=99_851), 100_090),
    ],
    ids=["plus", "minus", "minus_crossing"],
)  # fmt: skip
def test_ostir_tick_is_exact(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, row: dict[str, object], want: int
) -> None:
    (view,) = _views(_rows_run(tmp_path, cfg, nov, [row]))
    ticks = _ticks(view.figure)
    assert len(ticks) == 1, ticks
    assert abs(ticks[0] - want) < EXACT, (ticks[0], want)


@pytest.mark.parametrize(
    "row",
    [
        orf("NOVEL_chr1_5001_p", start_longest=5_001, start_best_ostir=5_001, end=5_300),
        orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500,
            start_best_ostir=1_500, end=1_201),
        orf("NOVEL_chr1_99851_m", cls="antisense", strand="-", start_longest=150,
            start_best_ostir=150, end=99_851),
    ],
    ids=["plus", "minus", "minus_crossing"],
)  # fmt: skip
def test_no_tick_when_ostir_equals_the_drawn_start(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType, row: dict[str, object]
) -> None:
    (view,) = _views(_rows_run(tmp_path, cfg, nov, [row]))
    assert _ticks(view.figure) == []


# --- decision 13: a one-entry key for the novel ORF marks -----------------------------------------


def _legends(fig: Figure) -> list[Legend]:
    out = list(fig.legends)
    out += [lg for ax in fig.axes if (lg := ax.get_legend()) is not None]
    return out


def test_every_locus_view_has_a_novel_orf_key(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_5001_p", start_longest=5_001, start_best_ostir=5_031, end=5_300, rank=1),
        orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500, end=1_201,
            rank=2),
    ]  # fmt: skip
    views = _views(_rows_run(tmp_path, cfg, nov, rows))
    assert len(views) == 2
    for view in views:
        _agg(view.figure)
        entries = [
            (lg, i, _collapse(t.get_text()))
            for lg in _legends(view.figure)
            for i, t in enumerate(lg.get_texts())
            if "novel ORF" in t.get_text()
        ]
        assert len(entries) == 1, (view.suffix, [e[2] for e in entries])
        lg, i, text = entries[0]
        assert "OSTIR" in text and "tick" in text, text
        assert "pileup" in text and "peptides" in text, text
        handle = lg.legend_handles[i]
        cols = set()
        if isinstance(handle, Patch):
            cols = {_hex(handle.get_edgecolor()), _hex(handle.get_facecolor())}
        elif isinstance(handle, Line2D):
            cols = {_hex(handle.get_color()), _hex(handle.get_markerfacecolor())}
        assert ORF_EDGE in cols, (view.suffix, cols)


# --- decision 14: bar label, title, rank NA ------------------------------------------------------


def test_class_with_only_uncalled_finite_rows_is_labelled_none_called(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, class_fdr="0.004", rank=1),
        orf("NOVEL_chr1_3001_p", cls="sorf", start_longest=3_001, end=3_090, class_fdr="0.3",
            called=False, rank=2),
        orf("NOVEL_chr1_4001_p", cls="sorf", start_longest=4_001, end=4_090, class_fdr="0.4",
            called=False, rank=3),
    ]  # fmt: skip
    ax = _axes(_main(_rows_run(tmp_path, cfg, nov, rows)).figure, "counts")
    texts = sorted(_collapse(t.get_text()) for t in ax.texts if t.get_text().strip())
    assert texts == ["FDR ≤ 0.004", "none called"]


def test_singular_title(tmp_path: Path, cfg: GippaConfig, nov: ModuleType) -> None:
    rows = [orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, class_fdr="0.004")]
    main = _main(_rows_run(tmp_path, cfg, nov, rows))
    assert _collapse(main.title) == "1 novel ORF called (class FDR ≤ 0.004), 0 unestimated"


def test_title_fdr_bound_is_the_max_over_called_rows_of_every_class(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, class_fdr="0.004", rank=1),
        orf("NOVEL_chr1_2001_p", start_longest=2_001, end=2_300, class_fdr="0.001", rank=2),
        orf("NOVEL_chr1_3300_m", cls="antisense", strand="-", start_longest=3_300, end=3_001,
            class_fdr="0.02", rank=3),
        orf("NOVEL_chr1_4300_m", cls="antisense", strand="-", start_longest=4_300, end=4_001,
            class_fdr="0.005", rank=6),
        orf("NOVEL_chr1_5001_p", cls="sorf", start_longest=5_001, end=5_090, class_fdr="0.01",
            rank=4),
        orf("NOVEL_chr1_7001_p", cls="pseudogene", start_longest=7_001, end=7_300,
            class_fdr="0.5", called=False, rank=5),
    ]  # fmt: skip
    main = _main(_rows_run(tmp_path, cfg, nov, rows))
    assert _collapse(main.title) == "5 novel ORFs called (class FDR ≤ 0.02), 0 unestimated"


def test_called_row_with_na_rank_sorts_after_every_ranked_row(
    tmp_path: Path, cfg: GippaConfig, nov: ModuleType
) -> None:
    rows = [
        orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, rank=None),  # smallest ID
        orf("NOVEL_chr1_5001_p", start_longest=5_001, end=5_300, rank=2),
        orf("NOVEL_chr1_9001_p", start_longest=9_001, end=9_300, rank=1),
    ]
    views = _views(_rows_run(tmp_path, cfg, nov, rows))
    assert [v.suffix for v in views] == [
        "locus_NOVEL_chr1_9001_p",
        "locus_NOVEL_chr1_5001_p",
        "locus_NOVEL_chr1_1001_p",
    ]
