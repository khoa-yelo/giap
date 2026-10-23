"""V1 ``circos``, review round of PR #212 (M9-08 lead decisions 5-19; review items B1-B8).

Builders and readers come from ``tests/unit/viz/test_v01_circos.py``; nothing here touches
``$GIAP_DB_DIR`` or the network. Angles are read from the drawn artists (polar data coordinates,
radians) and converted to replicon positions with a per-layout calibration: a ``+`` gene
``1..half`` on one replicon is one arc whose angular extent is exactly ``[x(0), x(half))``.

Interface choices (where the spec and the decisions leave it open):

- ``novel_orf`` rows span the ORF's low to high coordinate after unwrapping with the strand
  (decision 14): on ``+`` the start column is the low end, on ``-`` the high end; an ORF whose stop
  lies before its start in transcription direction gets ``end > length``. The position column is
  ``start_best_ostir``, else ``start_longest``, else a 1-bp row at ``end`` (decision 15).
- ``tu_start`` rows are 1 bp: ``start == end ==`` the 1-based transcription start (the low
  coordinate on ``+``, the high one on ``-``, modulo the length; decision 6, 14). The tick is drawn
  within 1.5 bp of that position. Promoter ticks sit at the BED ``start`` on ``+`` and the BED
  ``end`` on ``-`` (the transcription start).
- Linear-contig gap marks are artists with ``gid = "sector"`` near both ends of the sector.
- Binned tick tracks (decision 19) write ``binned = true`` rows whose ``value`` is the number of
  elements whose tick position (1-based) lies in ``[start, end]``; the bins of every binned track
  total at most 20,000 over all sectors (per-sector share of the cap). Binned gene rows carry the
  share of ``detected`` genes among all genes whose start lies in the bin (every status counts in
  the denominator); binned abundance rows the median of the measured values in the bin (``NA``
  without one).
- Text layout (decision 17) is measured on every visible non-blank ``Text`` of the figure after a
  ``FigureCanvasAgg`` draw (window extents, 0.5-px tolerance). A sector name may be truncated with
  "…" (a prefix of the replicon id followed by "…").
- The missing-layers note (decision 18) is judged on the texts a flagged run adds to an unflagged
  run of the same data.
- With support files present (``plots/support/replicon_landmarks.tsv``, ``gc_skew.tsv``, M9-03
  columns), ``oric``/``ter`` rows are 1 bp at the 1-based landmark and ``gc_skew`` rows keep the
  1-based window coordinates.
"""

from __future__ import annotations

import io
import itertools
import math
import re
import struct
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.viz_tracks import elements, tagged
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st
from matplotlib.artist import Artist
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import Collection, LineCollection, PathCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from unit.viz.test_v01_circos import (
    FULL_ISLANDS,
    FULL_NOVEL,
    FULL_REGIONS,
    FULL_VARIANTS,
    TWO,
    full_run,
    gene,
    make_run,
    novel,
    pep,
    region,
    render,
    rows,
    spans,
    texts,
    tsv,
)

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.tsv import write_tsv
from giap.viz import style
from giap.viz.base import Instance

L = 100_000
BASE_TITLE = "3 of 5 genes detected (60 %)"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- geometry helpers -----------------------------------------------------------------------


def _xy(artist: Artist) -> list[np.ndarray]:
    """Per drawn element, its (theta, r) vertices in polar data coordinates."""
    if not artist.get_visible():
        return []
    if isinstance(artist, LineCollection):
        return [np.asarray(s, dtype=float).reshape(-1, 2) for s in artist.get_segments()]
    if isinstance(artist, PathCollection):
        return [np.asarray(o, dtype=float).reshape(1, 2) for o in artist.get_offsets()]
    if isinstance(artist, Collection):
        return [np.asarray(p.vertices, dtype=float) for p in artist.get_paths()]
    if isinstance(artist, Patch):
        path = artist.get_patch_transform().transform_path(artist.get_path())
        return [np.asarray(path.vertices, dtype=float)]
    if isinstance(artist, Line2D):
        xy = np.column_stack([artist.get_xdata(), artist.get_ydata()]).astype(float)
        if artist.get_marker() in (None, "", "None", " ", "none"):
            return [xy]
        return [row.reshape(1, 2) for row in xy]
    return []


def shapes(fig: Figure, track: str) -> list[np.ndarray]:
    found: list[np.ndarray] = []
    for a in tagged(fig).get(track, []):
        found += _xy(a)
    return found


class Axis:
    """theta <-> position on one replicon of one layout (calibrated from a gene arc)."""

    def __init__(self, theta0: float, per_bp: float, length: int) -> None:
        self.theta0, self.k, self.length = theta0, per_bp, length

    def theta(self, x: float) -> float:
        return self.theta0 + self.k * x

    def pos(self, theta: float) -> float:
        return (theta - self.theta0) / self.k

    def to_origin(self, theta: float) -> float:
        """Distance in bp from ``theta`` to the replicon's origin (either side)."""
        turn = 2 * math.pi / self.k
        x = ((theta - self.theta0) % (2 * math.pi)) / self.k
        return min(x, abs(x - self.length), turn - x)


def calibrate(
    tmp: Path, cfg: GiapConfig, replicons: Sequence[tuple[str, int, str]], contig: str
) -> Axis:
    length = next(n for c, n, _ in replicons if c == contig)
    half = length // 2
    g = [gene("cal", contig, 1, half, "+", "detected", length=length)]
    inst = render(make_run(tmp, replicons=replicons, genes=g, name="CAL_" + contig), cfg)
    (arc,) = shapes(inst.figure, "gene_plus")
    lo, hi = float(arc[:, 0].min()), float(arc[:, 0].max())
    return Axis(lo, (hi - lo) / half, length)


def text_boxes(fig: Figure) -> tuple[list[tuple[str, Any]], float, float]:
    canvas = fig.canvas if isinstance(fig.canvas, FigureCanvasAgg) else FigureCanvasAgg(fig)
    canvas.draw()
    r = canvas.get_renderer()
    boxes = [
        (t.get_text(), t.get_window_extent(r))
        for t in fig.findobj(Text)
        if t.get_visible() and t.get_text().strip()
    ]
    return boxes, fig.bbox.width, fig.bbox.height


def png_size(data: bytes) -> tuple[int, int]:
    w, h = struct.unpack(">II", data[16:24])
    return int(w), int(h)


def orf(
    oid: str,
    strand: str,
    end: int,
    longest: int | None,
    ostir: int | None,
    contig: str = "chr1",
) -> dict[str, object]:
    row = novel(oid, contig, strand, 1, end, True)
    row["start_longest"], row["start_best_ostir"] = longest, ostir
    return row


# --- B1 / B2 (decision 14): origin-crossing novel ORFs and '-' TU starts --------------------


def test_b1_origin_crossing_novel_orfs_rows_unwrap_with_the_strand(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    orfs = [
        orf("nplus", "+", 150, 99_901, 99_901),  # 99,901 .. 100,150 across the origin
        orf("nminus", "-", 99_950, 100, 100),  # high end 100,100 wrapped to 100
    ]
    data = tsv(render(make_run(tmp_path, novel_orfs=orfs), cfg), tmp_path)
    assert spans(data, "novel_orf") == [("chr1", 99_901, 100_150), ("chr1", 99_950, 100_100)]


def test_b1_origin_crossing_novel_orfs_are_drawn_at_the_origin(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    orfs = [orf("nplus", "+", 150, 99_901, 99_901), orf("nminus", "-", 99_950, 100, 100)]
    fig = render(make_run(tmp_path, novel_orfs=orfs), cfg).figure
    pts = shapes(fig, "novel_orf")
    assert pts
    far = [round(ax.to_origin(float(t))) for s in pts for t in s[:, 0]]
    assert max(far) <= 300, far  # every drawn piece lies on the short arc across the origin


def test_b2_minus_tu_across_the_origin_wraps_its_start(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    tus = [["TU_x", "chr1", 99_001, 101_000, "-", ["g1"], 101_001, 101_200]]
    inst = render(make_run(tmp_path, tus=tus), cfg)
    assert spans(tsv(inst, tmp_path), "tu_start") == [("chr1", 1_000, 1_000)]
    (tick,) = shapes(inst.figure, "tu_start")
    got = [ax.pos(float(t)) for t in tick[:, 0]]
    assert all(abs(x - 1_000) <= 1.5 for x in got), got


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(start=st.integers(1, L), width=st.integers(3, 3_000), minus=st.booleans())
@example(start=99_901, width=250, minus=False)
@example(start=99_950, width=151, minus=True)
def test_b1_novel_orf_rows_are_positions_modulo_length_never_clamped(
    tmp_path_factory: pytest.TempPathFactory, start: int, width: int, minus: bool
) -> None:
    tmp = tmp_path_factory.mktemp("prop")
    lo, hi = start, start + width - 1  # unwrapped, 1-based; hi may exceed L
    wrap = lambda x: (x - 1) % L + 1  # noqa: E731
    row = orf("n", "-", wrap(lo), wrap(hi), wrap(hi)) if minus else orf("n", "+", wrap(hi), lo, lo)
    data = tsv(render(make_run(tmp, novel_orfs=[row]), resolve_config()), tmp)
    assert spans(data, "novel_orf") == [("chr1", lo, hi)]


# --- decision 15: novel-ORF position fallback -----------------------------------------------


def test_decision15_novel_orf_position_fallback(tmp_path: Path, cfg: GiapConfig) -> None:
    orfs = [
        orf("a", "+", 80_300, 80_001, 80_101),  # OSTIR start wins
        orf("b", "+", 82_300, 82_001, None),  # then the longest start
        orf("c", "+", 84_300, None, None),  # then the stop codon, 1 bp
        orf("d", "-", 86_001, 86_400, 86_300),  # '-': OSTIR start is the high end
        orf("e", "-", 88_001, None, None),
    ]
    data = tsv(render(make_run(tmp_path, novel_orfs=orfs), cfg), tmp_path)
    assert spans(data, "novel_orf") == [
        ("chr1", 80_101, 80_300),
        ("chr1", 82_001, 82_300),
        ("chr1", 84_300, 84_300),
        ("chr1", 86_001, 86_300),
        ("chr1", 88_001, 88_001),
    ]


# --- B3 (decision 16): rcParams after importing pyCirclize ----------------------------------


def test_b3_style_rc_pins_savefig_bbox_and_padding() -> None:
    assert "savefig.bbox" in style.RC and style.RC["savefig.bbox"] is None
    assert style.RC["savefig.pad_inches"] == pytest.approx(0.1)


def test_b3_other_figures_keep_their_pixel_size_after_importing_circos(tmp_path: Path) -> None:
    import giap.viz.circos  # noqa: F401  (imports pyCirclize)

    style.apply()
    fig = style.new_figure("single", 2.0)
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    ax.plot([0, 1], [0, 1])
    fig.text(0.99, 0.5, "x" * 40)  # runs past the right edge: a tight bbox would grow
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    assert png_size(buf.getvalue()) == (350, 200)
    inst = Instance(None, fig, pd.DataFrame({"start": [0]}), "t", coords=("start",))
    (png,) = [
        p
        for p in write_instance(inst, "V99", "probe", tmp_path, ["png"], 100)
        if p.suffix == ".png"
    ]
    assert png_size(png.read_bytes()) == (350, 200)


def test_b3_v01_saves_at_double_width(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = render(full_run(tmp_path), cfg)
    h = float(inst.figure.get_size_inches()[1])
    (png,) = [
        p
        for p in write_instance(inst, "V01", "circos", tmp_path / "p", ["png"], 100)
        if p.suffix == ".png"
    ]
    assert png_size(png.read_bytes()) == (720, round(100 * h))


# --- B4 / B5 (decision 17): text inside the canvas, no overlaps -----------------------------

NCBI = [
    ("NZ_CP000001.1", 4_600_000, "circular"),
    ("NZ_CP000002.1_plasmid_unnamed1", 60_000, "circular"),
    ("NZ_CP000003.1_plasmid_unnamed2", 20_000, "circular"),
    ("contig_00042_linear", 300_000, "linear"),
    ("NZ_CP000004.1", 1_100_000, "circular"),
]
SKEW_BP = 10_000


def _support(out: Path, replicons: Sequence[tuple[str, int, str]]) -> None:
    sup = out / "plots" / "support"
    sup.mkdir(parents=True, exist_ok=True)
    marks = [
        [c, n, t, 1, n // 2 + 1, "called"] if t == "circular" else [c, n, t, None, None, "linear"]
        for c, n, t in replicons
    ]
    write_tsv(
        sup / "replicon_landmarks.tsv",
        ["contig", "length_bp", "topology", "oric_position", "ter_position", "status"],
        marks,
    )
    skew = []
    for c, n, _ in replicons:
        cum = 0.0
        for s in range(0, n, SKEW_BP):
            g = 0.05 * math.sin(2 * math.pi * s / n)
            cum += g
            skew.append([c, s + 1, min(s + SKEW_BP, n), round(g, 6), round(cum, 6)])
    write_tsv(
        sup / "gc_skew.tsv", ["contig", "start", "end", "gc_skew", "cumulative_gc_skew"], skew
    )


@pytest.fixture
def ncbi_inst(tmp_path: Path, cfg: GiapConfig) -> Instance:
    genes = [
        gene("a1", "NZ_CP000001.1", 1_001, 2_000, "+", "detected", 20.0, length=4_600_000),
        gene("b1", "NZ_CP000002.1_plasmid_unnamed1", 101, 900, "-", "not_detected", length=60_000),
        gene(
            "c1", "NZ_CP000003.1_plasmid_unnamed2", 101, 900, "+", "detected", 21.0, length=20_000
        ),
    ]
    regs = [region("plasmid_1", "plasmid", "NZ_CP000002.1_plasmid_unnamed1", 1, 60_000)]
    out = make_run(
        tmp_path,
        replicons=NCBI,
        genes=genes,
        regions=regs,
        flags=["layer_missing:genomad", "layer_missing:antismash"],
    )
    _support(out, NCBI)
    return render(out, cfg)


def _inside(fig: Figure) -> list[str]:
    boxes, w, h = text_boxes(fig)
    tol = 0.5
    return [t for t, b in boxes if b.x0 < -tol or b.y0 < -tol or b.x1 > w + tol or b.y1 > h + tol]


def _overlaps(fig: Figure) -> list[tuple[str, str]]:
    boxes, _, _ = text_boxes(fig)
    tol = 0.5
    return [
        (a[:30], b[:30])
        for (a, ba), (b, bb) in itertools.combinations(boxes, 2)
        if min(ba.x1, bb.x1) - max(ba.x0, bb.x0) > tol
        and min(ba.y1, bb.y1) - max(ba.y0, bb.y0) > tol
    ]


def test_b4_every_text_inside_the_canvas_5_replicons(ncbi_inst: Instance) -> None:
    assert round(float(ncbi_inst.figure.get_size_inches()[0]), 6) == style.DOUBLE_IN
    assert _inside(ncbi_inst.figure) == []


def test_b4_every_text_inside_the_canvas_full_run(tmp_path: Path, cfg: GiapConfig) -> None:
    out = full_run(tmp_path, peptides=None, flags=["layer_missing:genomad"])
    assert _inside(render(out, cfg).figure) == []


def test_b5_no_two_texts_overlap_5_replicons(ncbi_inst: Instance) -> None:
    assert _overlaps(ncbi_inst.figure) == []


def test_b5_every_sector_is_named_whole_or_truncated(ncbi_inst: Instance) -> None:
    shown = texts(ncbi_inst.figure)
    for name, *_ in NCBI:
        assert any(
            t == name or (t.endswith("…") and len(t) > 1 and name.startswith(t[:-1])) for t in shown
        ), name


def test_b5_oric_and_ter_labelled_on_every_circular_sector(ncbi_inst: Instance) -> None:
    shown = Counter(texts(ncbi_inst.figure))
    circular = sum(t == "circular" for *_, t in NCBI)
    assert shown["oriC"] == circular
    assert shown["ter"] == circular


def test_example7_landmark_and_gc_skew_rows(ncbi_inst: Instance, tmp_path: Path) -> None:
    data = tsv(ncbi_inst, tmp_path)
    circ = [(c, n) for c, n, t in NCBI if t == "circular"]
    assert spans(data, "oric") == [(c, 1, 1) for c, _ in circ]
    assert spans(data, "ter") == [(c, n // 2 + 1, n // 2 + 1) for c, n in circ]
    g = data[(data["track"] == "gc_skew") & (data["contig"] == "NZ_CP000003.1_plasmid_unnamed2")]
    assert list(zip(g["start"].astype(int), g["end"].astype(int), strict=True)) == [
        (1, 10_000),
        (10_001, 20_000),
    ]


# --- B6 (decision 18): the missing-layers note ----------------------------------------------

NO_RING = ("cog", "go", "kegg", "tx_rates", "features", "hits_tolerant")


def test_b6_missing_layers_note_names_only_v1_rings(tmp_path: Path, cfg: GiapConfig) -> None:
    base_texts = set(texts(render(full_run(tmp_path, name="A"), cfg).figure))
    flags = [f"layer_missing:{x}" for x in (*NO_RING, "genomad", "antismash")]
    flagged = set(texts(render(full_run(tmp_path, name="B", flags=flags), cfg).figure))
    extra = " ".join(sorted(flagged - base_texts)).lower()
    assert "genomad" in extra and "antismash" in extra
    for word in NO_RING:
        assert not re.search(rf"\b{re.escape(word)}\b", extra), (word, extra)


def test_b6_only_ringless_layers_missing_adds_no_layer_note(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    base_texts = set(texts(render(full_run(tmp_path, name="A"), cfg).figure))
    flags = [f"layer_missing:{x}" for x in NO_RING]
    flagged = set(texts(render(full_run(tmp_path, name="B", flags=flags), cfg).figure))
    extra = " ".join(sorted(flagged - base_texts)).lower()
    for word in NO_RING:
        assert not re.search(rf"\b{re.escape(word)}\b", extra), (word, extra)


# --- B7 (decision 19) and binned gene tracks ------------------------------------------------

BIG = [("chr1", 5_000_000, "circular"), ("p1", 50_000, "circular")]
STATUSES = (
    "detected",
    "single_peptide_unestimated",
    "not_detected",
    "detected",
    "partial_excluded",
)


def _big_genes() -> list[dict[str, object]]:
    out = []
    for c, n, step, count in (("chr1", 5_000_000, 200, 25_000), ("p1", 50_000, 400, 125)):
        for i in range(count):
            status = STATUSES[i % len(STATUSES)]
            log2 = None if status == "not_detected" else float((3 * i) % 11)
            s = 1 + step * i
            out.append(gene(f"{c}_{i}", c, s, s + 149, "+", status, log2, length=n))
    return out


def _tick_positions() -> dict[str, list[tuple[str, int]]]:
    """1-based tick positions per tick track: 20,001 elements each (19,901 + 100)."""
    pos: dict[str, list[tuple[str, int]]] = {}
    for track, off in (("promoter", 11), ("terminator", 101), ("tu_start", 51)):
        pos[track] = [("chr1", off + 250 * i) for i in range(19_901)]
        pos[track] += [("p1", off + 500 * i) for i in range(100)]
    return pos


@pytest.fixture(scope="module")
def big(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    style.apply()
    tmp = tmp_path_factory.mktemp("big")
    pos = _tick_positions()
    promoters = [(c, p - 1, p + 49, f"P{i}", 1, "+") for i, (c, p) in enumerate(pos["promoter"])]
    terminators = [
        (c, p - 1, p + 29, f"T{i}", 1, "+") for i, (c, p) in enumerate(pos["terminator"])
    ]
    tus = [
        [f"TU_{i}", c, p, p + 99, "+", [f"g{i}"], max(1, p - 200), max(1, p - 1)]
        for i, (c, p) in enumerate(pos["tu_start"])
    ]
    out = make_run(
        tmp,
        replicons=BIG,
        genes=_big_genes(),
        promoters=promoters,
        terminators=terminators,
        tus=tus,
    )
    inst = render(out, resolve_config())
    return inst, tsv(inst, tmp)


def _bins_ok(d: pd.DataFrame) -> None:
    length = dict((c, n) for c, n, _ in BIG)
    for c, g in d.groupby("contig"):
        s, e = g["start"].astype(int).to_numpy(), g["end"].astype(int).to_numpy()
        assert (s <= e).all() and (s[1:] > e[:-1]).all(), c  # sorted, disjoint
        assert s.min() >= 1 and e.max() <= length[str(c)], c


@pytest.mark.parametrize("track", ["promoter", "terminator", "tu_start"])
def test_b7_tick_tracks_binned_above_20000(big: tuple[Instance, pd.DataFrame], track: str) -> None:
    inst, data = big
    d = data[data["track"] == track]
    assert 0 < len(d) <= 20_000
    assert set(d["binned"]) == {"true"}
    assert len(elements(inst.figure, track)) <= 20_000
    _bins_ok(d)
    assert d["value"].astype(float).sum() == 20_001
    pos = _tick_positions()[track]
    for c in ("chr1", "p1"):
        x = np.sort(np.array([p for cc, p in pos if cc == c]))
        g = d[d["contig"] == c]
        s, e = g["start"].astype(int).to_numpy(), g["end"].astype(int).to_numpy()
        want = np.searchsorted(x, e, side="right") - np.searchsorted(x, s, side="left")
        assert g["value"].astype(float).astype(int).tolist() == want.tolist(), c


def test_b8_binned_gene_bins_are_per_sector_shares_of_the_cap(
    big: tuple[Instance, pd.DataFrame],
) -> None:
    _, data = big
    d = data[data["track"] == "gene_plus"]
    assert set(d["binned"]) == {"true"}
    assert len(d) <= 20_000  # total over both sectors, not 20,000 per sector
    assert set(d["contig"]) == {"chr1", "p1"}
    _bins_ok(d)


def _in_bins(d: pd.DataFrame, genes: pd.DataFrame) -> list[pd.DataFrame]:
    out = []
    for c, s, e in zip(d["contig"], d["start"].astype(int), d["end"].astype(int), strict=True):
        out.append(genes[(genes["contig"] == c) & (genes["start"] >= s) & (genes["start"] <= e)])
    return out


def test_b8_detected_share_is_detected_over_all_genes_in_the_bin(
    big: tuple[Instance, pd.DataFrame],
) -> None:
    _, data = big
    genes = pd.DataFrame(_big_genes())
    d = data[data["track"] == "gene_plus"]
    members = _in_bins(d, genes)
    assert sum(len(m) for m in members) == len(genes)  # every gene in exactly one written bin
    want = [float((m["detection_status"] == "detected").mean()) for m in members]
    assert d["value"].astype(float).tolist() == pytest.approx(want)


def test_b8_binned_abundance_is_the_bin_median(big: tuple[Instance, pd.DataFrame]) -> None:
    _, data = big
    genes = pd.DataFrame(_big_genes())
    d = data[data["track"] == "abundance"]
    assert set(d["binned"]) == {"true"}
    assert len(d) <= 20_000
    got, want = [], []
    for v, m in zip(d["value"], _in_bins(d, genes), strict=True):
        vals = m["abundance_log2"].dropna().astype(float)
        want.append(None if vals.empty else float(np.median(vals)))
        got.append(None if v == "NA" else float(v))
    assert got == want
    assert any(w is not None and w != max(w2 for w2 in want if w2 is not None) for w in want)


def test_b8_exactly_20000_genes_are_not_binned(tmp_path: Path, cfg: GiapConfig) -> None:
    n = 5_000_000
    genes = [
        gene(f"g{i}", "chr1", 1 + 200 * i, 150 + 200 * i, "+", length=n) for i in range(20_000)
    ]
    data = tsv(
        render(make_run(tmp_path, replicons=[("chr1", n, "circular")], genes=genes), cfg), tmp_path
    )
    d = data[data["track"] == "gene_plus"]
    assert set(d["binned"]) == {"false"}
    assert len(d) == 20_000


# --- B8: lead decisions 5, 6, 8, 10, 12, 13a ------------------------------------------------


def test_decision5_amr_without_coordinates_not_drawn_and_noted(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    with_xy = render(full_run(tmp_path, name="A"), cfg)
    regs = [*FULL_REGIONS, region("amr_2", "amr", "chr1", None, None)]
    inst = render(full_run(tmp_path, name="B", regions=regs), cfg)
    data = tsv(inst, tmp_path)
    assert [cat for *_, cat in rows(data, "region")].count("amr") == 1  # the one with coordinates
    assert len(elements(inst.figure, "region")) == 5
    extra = set(texts(inst.figure)) - set(texts(with_xy.figure))
    assert any("amr" in t.lower() for t in extra), extra


def test_decision6_minus_tu_start_is_the_high_coordinate(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    tus = [
        ["TU_1", "chr1", 1_001, 2_000, "+", ["g2"], 801, 1_000],
        ["TU_2", "chr1", 3_001, 6_000, "-", ["g4"], 6_001, 6_200],
    ]
    inst = render(make_run(tmp_path, tus=tus), cfg)
    assert spans(tsv(inst, tmp_path), "tu_start") == [
        ("chr1", 1_001, 1_001),
        ("chr1", 6_000, 6_000),
    ]
    got = sorted(
        float(np.mean([ax.pos(t) for t in s[:, 0]])) for s in shapes(inst.figure, "tu_start")
    )
    assert len(got) == 2
    assert abs(got[0] - 1_000) <= 1.5 and abs(got[1] - 6_000) <= 1.5, got


def test_promoter_ticks_at_the_transcription_start(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    beds = [("chr1", 2_000, 2_050, "p_plus", 500, "+"), ("chr1", 900, 950, "p_minus", 500, "-")]
    inst = render(make_run(tmp_path, promoters=beds), cfg)
    got = sorted(
        float(np.mean([ax.pos(t) for t in s[:, 0]])) for s in shapes(inst.figure, "promoter")
    )
    assert len(got) == 2
    assert abs(got[0] - 950) <= 1.5 and abs(got[1] - 2_000) <= 1.5, got


def test_decision8_singular_counts_in_the_title(tmp_path: Path, cfg: GiapConfig) -> None:
    out = full_run(
        tmp_path, islands=FULL_ISLANDS[:1], variants=FULL_VARIANTS[:1], novel_orfs=FULL_NOVEL[:1]
    )
    assert render(out, cfg).title == (
        f"{BASE_TITLE}; 1 island at FWER ≤ 0.05; 1 ALT site, 1 novel ORF"
    )


def test_decision8_suffix_lists_only_present_tables_zero_is_plural(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    a = render(full_run(tmp_path, name="A", variants=None), cfg)
    assert a.title == f"{BASE_TITLE}; 2 islands at FWER ≤ 0.05; 2 novel ORFs"
    b = render(full_run(tmp_path, name="B", variants=[]), cfg)
    assert b.title == f"{BASE_TITLE}; 2 islands at FWER ≤ 0.05; 0 ALT sites, 2 novel ORFs"


def test_decision10_every_density_bin_has_a_row(tmp_path: Path, cfg: GiapConfig) -> None:
    data = tsv(render(full_run(tmp_path), cfg), tmp_path)
    d = data[data["track"] == "peptide_density"]
    want = [("chr1", 1 + 10_000 * k, 10_000 * (k + 1)) for k in range(10)]
    assert spans(d, "peptide_density") == [*want, ("plasmid_1", 1, 8_000)]
    assert d["value"].astype(float).tolist() == [2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]


def test_density_bars_use_a_per_sector_scale(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    peps = [pep("chr1", s) for s in (5, 9_999, 10_001)] + [pep("plasmid_1", 101)]
    fig = render(make_run(tmp_path, peptides=peps), cfg).figure
    bars: dict[str, list[float]] = {"chr1": [], "plasmid_1": []}
    for s in shapes(fig, "peptide_density"):
        height = float(s[:, 1].max() - s[:, 1].min())
        if height <= 1e-9:
            continue
        mid = ax.pos(float(np.mean(s[:, 0])))
        bars["chr1" if 0 <= mid <= L else "plasmid_1"].append(height)
    (p,) = bars["plasmid_1"]
    tall, short = sorted(bars["chr1"], reverse=True)
    assert tall == pytest.approx(2 * short, rel=1e-6)  # counts 2 and 1 on chr1
    assert p == pytest.approx(tall, rel=1e-6)  # the plasmid's 1 is its own maximum


def test_decision12_islands_absent_drop_the_island_clause(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = render(full_run(tmp_path, islands=None), cfg)
    assert inst.title == f"{BASE_TITLE}; 2 ALT sites, 2 novel ORFs"


def test_decision13a_zero_genes_draws_the_circle(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_run(tmp_path, regions=FULL_REGIONS)
    inst = render(out, cfg)
    assert inst.title == "0 of 0 genes detected (NA)"
    data = tsv(inst, tmp_path)
    s = data[data["track"] == "sector"]
    assert s["value"].astype(int).tolist() == [100_000, 8_000]
    assert len(elements(inst.figure, "region")) == 5


def test_rows_with_equal_coordinates_sort_by_category(tmp_path: Path, cfg: GiapConfig) -> None:
    regs = [
        region("x_plasmid", "plasmid", "chr1", 10_001, 20_000),
        region("x_crispr", "crispr", "chr1", 10_001, 20_000),
    ]
    data = tsv(render(make_run(tmp_path, regions=regs), cfg), tmp_path)
    assert [cat for *_, cat in rows(data, "region")] == ["crispr", "plasmid"]


# --- B8: wrap geometry, linear clipping and gap marks ---------------------------------------


def test_second_piece_of_a_wrapped_gene_spans_0_to_end_minus_length(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    ax = calibrate(tmp_path, cfg, TWO, "chr1")
    g = [gene("g1", "chr1", 99_001, 101_000, "+", "detected")]
    fig = render(make_run(tmp_path, genes=g), cfg).figure
    got = sorted(
        (round(ax.pos(float(s[:, 0].min())), 2), round(ax.pos(float(s[:, 0].max())), 2))
        for s in shapes(fig, "gene_plus")
    )
    assert len(got) == 2
    (a0, a1), (b0, b1) = got
    assert abs(a0) <= 0.25 and abs(a1 - 1_000) <= 0.25, got  # [0, end - length)
    assert abs(b0 - 99_000) <= 0.25 and abs(b1 - L) <= 0.25, got  # [start, length)


LINEAR = [("chr1", L, "linear")]


def test_linear_contig_clips_a_feature_past_its_end(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, LINEAR, "chr1")
    g = [gene("g1", "chr1", 99_501, 100_400, "+", "detected")]
    fig = render(make_run(tmp_path, replicons=LINEAR, genes=g), cfg).figure
    (arc,) = shapes(fig, "gene_plus")
    lo, hi = ax.pos(float(arc[:, 0].min())), ax.pos(float(arc[:, 0].max()))
    assert abs(lo - 99_500) <= 0.25 and abs(hi - L) <= 0.25, (lo, hi)


def test_linear_contig_has_gap_marks_at_both_ends(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = calibrate(tmp_path, cfg, LINEAR, "chr1")
    fig = render(make_run(tmp_path, replicons=LINEAR), cfg).figure
    marks = [ax.pos(float(np.mean(s[:, 0]))) for s in shapes(fig, "sector")]
    assert len(marks) >= 2
    assert any(abs(m) <= 300 for m in marks), marks
    assert any(abs(m - L) <= 300 for m in marks), marks
