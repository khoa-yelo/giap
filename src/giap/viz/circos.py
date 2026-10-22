"""V1 ``circos``: where in the genome is the proteome? (GIAP_SPEC.md §10.2 V1; M9-08).

One pyCirclize sector per replicon (FASTA order, angle ∝ length) with ticks; outer to inner:
feature arcs, genes on ``+`` and ``-`` coloured by ``detection_status``, an abundance heat track,
peptide density per 10-kb bin, TU starts with promoter / terminator ticks, claimed islands, the
cumulative GC skew and the ALT-site / novel-ORF markers; oriC / ter as radial lines. pyCirclize
lays out the sectors and ticks; every data track is drawn as its own tagged collection
(``gid`` = the track's data-TSV enum value) on the polar axes. Data rows are 0-based half-open
(``start`` written ``+ 1``); circular intervals with ``end > length`` keep that end and are drawn
as two pieces. A track with more than 20,000 elements is binned (index decision 12).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.artist import Artist
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.patches import ArcRectangle

from giap.viz import base, style

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from pycirclize.sector import Sector

    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = ("track", "contig", "start", "end", "value", "category", "binned")
TRACKS = (
    "sector", "region", "gene_plus", "gene_minus", "abundance", "peptide_density", "tu_start",
    "promoter", "terminator", "island", "gc_skew", "alt_site", "novel_orf", "oric", "ter",
)  # fmt: skip
#: track -> (inner, outer) radius (pyCirclize units, 0-100)
RADII: dict[str, tuple[float, float]] = {
    "sector": (96.0, 98.0),
    "region": (91.0, 95.0),
    "gene_plus": (85.0, 90.0),
    "gene_minus": (79.0, 84.0),
    "abundance": (73.0, 78.0),
    "peptide_density": (63.0, 72.0),
    "ticks": (56.0, 62.0),
    "island": (50.0, 55.0),
    "gc_skew": (38.0, 49.0),
    "alt_site": (34.0, 34.0),
    "novel_orf": (30.0, 30.0),
}
MAX_ELEMENTS = 20_000
DENSITY_BP = 10_000
PARTIAL_HATCH = "////"
WHITE = "#FFFFFF"
BLACK = style.OKABE_ITO["black"]
PROMOTER = style.OKABE_ITO["bluish_green"]
TERMINATOR = style.OKABE_ITO["vermillion"]
NOVEL = style.OKABE_ITO["reddish_purple"]
Row = tuple[str, str, int, int, object, object, bool]


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _num(v: object) -> float | None:
    return None if _na(v) else float(str(v))


def pieces(start: int, end: int, length: int, circular: bool) -> list[tuple[int, int]]:
    """0-based half-open ``[start, end)`` as drawn: split at the origin on circular replicons,
    clipped at the end on linear ones."""
    start = max(0, start)
    if end <= length:
        return [(start, end)]
    if not circular:
        return [(start, length)]
    return [(start, length), (0, min(end - length, length))]


def tick_step(length: int) -> int:
    """The largest 1, 2 or 5 x 10^n step giving at least 4 (so at most 10) ticks."""
    step = 1
    for n in range(12):
        for m in (1, 2, 5):
            if length // (m * 10**n) >= 4:
                step = m * 10**n
    return step


def tick_label(x: float, step: int) -> str:
    if x == 0:
        return ""
    if step >= 1_000_000:
        return f"{x / 1e6:g} Mb"
    return f"{x / 1e3:g} kb" if step >= 1_000 else f"{x:g} bp"


def bin_edges(length: int, n: int) -> np.ndarray:
    """``n + 1`` integer edges of ``n`` equal-width bins over ``[0, length)``, as in
    ``base.binned``."""
    n = max(1, min(n, length))
    return (np.arange(n + 1, dtype=np.int64) * np.int64(length)) // np.int64(n)


def _groups(
    starts: pd.Series, length: int, total: int
) -> tuple[np.ndarray, list[tuple[int, np.ndarray]]]:
    """This sector's share of ``MAX_ELEMENTS`` bins and the row indices per non-empty bin."""
    edges = bin_edges(length, MAX_ELEMENTS * length // total)
    pos = starts.astype(int).to_numpy()
    if len(pos) == 0:
        return edges, []
    idx = np.clip(np.searchsorted(edges, pos, side="right") - 1, 0, len(edges) - 2)
    order = np.argsort(idx, kind="stable")
    keys, first = np.unique(idx[order], return_index=True)
    return edges, list(zip(keys.tolist(), np.split(order, first[1:]), strict=True))


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


class _Plot:
    """Collects the data rows and drawn pieces of every track for one run."""

    def __init__(self, replicons: Sequence[dict[str, Any]]) -> None:
        self.reps = {str(r["id"]): r for r in replicons}
        self.order = {c: i for i, c in enumerate(self.reps)}
        self.rows: list[Row] = []
        self.notes: list[str] = []

    def length(self, contig: str) -> int:
        return int(self.reps[contig]["length"])

    def circular(self, contig: str) -> bool:
        return str(self.reps[contig].get("topology", "circular")) == "circular"

    def known(self, frame: pd.DataFrame, column: str = "contig") -> pd.DataFrame:
        return frame[frame[column].astype(str).isin(list(self.reps))]

    def add(
        self, track: str, contig: str, s: int, e: int, v: object, cat: object, b: bool = False
    ) -> None:
        self.rows.append((track, contig, int(s), int(e), None if _na(v) else v, cat, b))

    def frame(self) -> pd.DataFrame:
        def key(r: Row) -> tuple[int, int, int, int, str]:
            return (
                TRACKS.index(r[0]),
                self.order[r[1]],
                r[2],
                r[3],
                "" if r[5] is None else str(r[5]),
            )

        return pd.DataFrame(sorted(self.rows, key=key), columns=list(COLUMNS), dtype=object)


class _Canvas:
    """Draws pieces on the polar axes, one tagged artist per (track, style)."""

    def __init__(self, ax: Axes, sectors: dict[str, Sector]) -> None:
        self.ax, self.sectors = ax, sectors

    def rad(self, contig: str, x: float) -> float:
        sector = self.sectors[contig]
        return float(sector.x_to_rad(min(max(x, 0.0), sector.size)))

    def arc(self, contig: str, a: float, b: float, r0: float, r1: float) -> ArcRectangle:
        ra, rb = self.rad(contig, a), self.rad(contig, b)
        return ArcRectangle((ra, r0), max(rb - ra, 1e-9), r1 - r0)

    def rects(self, gid: str, arcs: list[tuple[ArcRectangle, str]], **kw: Any) -> None:
        if not arcs:
            return
        coll = PatchCollection(
            [a for a, _ in arcs], facecolors=[c for _, c in arcs], gid=gid, zorder=2, **kw
        )
        coll.set_clip_on(False)
        self.ax.add_collection(coll)

    def ticks(self, gid: str, segs: list[tuple[str, float, float, float]], colour: str) -> None:
        if not segs:
            return
        lines = [[(self.rad(c, x), r0), (self.rad(c, x), r1)] for c, x, r0, r1 in segs]
        coll = LineCollection(lines, colors=colour, linewidths=0.5, gid=gid, zorder=3)
        coll.set_clip_on(False)
        self.ax.add_collection(coll)

    def markers(self, gid: str, pts: list[tuple[str, float]], r: float, **kw: Any) -> None:
        if pts:
            rads = [self.rad(c, x) for c, x in pts]
            self.ax.scatter(rads, [r] * len(pts), s=14, linewidths=0.6, gid=gid, zorder=4, **kw)


def _intervals(
    plot: _Plot, canvas: _Canvas, track: str, frame: pd.DataFrame, colour: Callable[[Any], str]
) -> list[tuple[ArcRectangle, str]]:
    """Data rows for ``frame`` (``contig``, ``start``, ``end``, ``category``, ``value``) and its
    arc pieces on ``track``."""
    r0, r1 = RADII[track]
    arcs = []
    for r in frame.itertuples(index=False):
        c = str(r.contig)
        plot.add(track, c, int(r.start), int(r.end), r.value, r.category)
        for a, b in pieces(int(r.start), int(r.end), plot.length(c), plot.circular(c)):
            arcs.append((canvas.arc(c, a, b, r0, r1), colour(r)))
    return arcs


def _fixed(colour: str) -> Callable[[Any], str]:
    return lambda _row: colour


def _abundance_colours(values: Sequence[float | None]) -> list[str]:
    finite = np.array([v for v in values if v is not None], dtype=float)
    if len(finite) == 0:
        return [WHITE] * len(values)
    lo, hi = (float(x) for x in np.percentile(finite, [2, 98]))
    norm = Normalize(lo, hi if hi > lo else lo + 1, clip=True)
    cmap = matplotlib.colormaps[style.CONTINUOUS]
    mid = matplotlib.colors.to_hex(cmap(0.5))
    hexes = [WHITE if v is None else matplotlib.colors.to_hex(cmap(norm(v))) for v in values]
    return hexes if hi > lo else [WHITE if v is None else mid for v in values]


def _genes(plot: _Plot, canvas: _Canvas, genes: pd.DataFrame) -> bool:
    """Gene and abundance tracks; ``True`` when binned."""
    genes = plot.known(genes).copy()
    imputed = genes["abundance_imputed"].astype("boolean").fillna(False).astype(bool)
    genes["measured"] = [
        None if i or _na(v) else float(v)
        for v, i in zip(genes["abundance_log2"], imputed, strict=True)
    ]
    total_len = sum(plot.length(c) for c in plot.reps)
    binned = len(genes) > MAX_ELEMENTS
    for track, strand in (("gene_plus", "+"), ("gene_minus", "-")):
        sub = genes[genes["strand"] == strand]
        r0, r1 = RADII[track]
        if not binned:
            normal = sub.assign(category=sub["detection_status"].astype(str), value=None)
            part = normal["category"] == "partial_excluded"
            for frame, kw in ((normal[~part], {}), (normal[part], {"hatch": PARTIAL_HATCH})):
                arcs = _intervals(
                    plot, canvas, track, frame, lambda r: style.STATUS_COLOUR[r.category]
                )
                canvas.rects(
                    track, arcs, edgecolors="#555555" if kw else "none", linewidths=0, **kw
                )
            continue
        arcs = []
        for c in plot.reps:
            g = sub[sub["contig"] == c]
            edges, groups = _groups(g["start"], plot.length(c), total_len)
            det = (g["detection_status"] == "detected").to_numpy()
            for k, ix in groups:
                share = float(det[ix].mean())
                plot.add(track, c, int(edges[k]), int(edges[k + 1]), share, "detected_share", True)
                hexc = matplotlib.colors.to_hex(matplotlib.colormaps[style.CONTINUOUS](share))
                arcs.append((canvas.arc(c, edges[k], edges[k + 1], r0, r1), hexc))
        canvas.rects(track, arcs, edgecolors="none", linewidths=0)
    _abundance(plot, canvas, genes, binned, total_len)
    return binned


def _abundance(plot: _Plot, canvas: _Canvas, genes: pd.DataFrame, binned: bool, total: int) -> None:
    r0, r1 = RADII["abundance"]
    cells: list[tuple[str, int, int, float | None]] = []
    for c in plot.reps:
        g = genes[genes["contig"] == c]
        if not binned:
            cells += [
                (c, int(s), int(e), _num(v))
                for s, e, v in zip(g["start"], g["end"], g["measured"], strict=True)
            ]
            continue
        edges, groups = _groups(g["start"], plot.length(c), total)
        vals = [_num(v) for v in g["measured"]]
        for k, ix in groups:
            m = [vals[i] for i in ix.tolist() if vals[i] is not None]
            cells.append((c, int(edges[k]), int(edges[k + 1]), float(np.median(m)) if m else None))
    colours = _abundance_colours([v for *_, v in cells])
    arcs = []
    for (c, s, e, v), hexc in zip(cells, colours, strict=True):
        plot.add("abundance", c, s, e, v, "bin_median" if binned else "abundance_log2", binned)
        for a, b in pieces(s, e, plot.length(c), plot.circular(c)):
            arcs.append((canvas.arc(c, a, b, r0, r1), hexc))
    canvas.rects("abundance", arcs, edgecolors="none", linewidths=0)


def _density(plot: _Plot, canvas: _Canvas, hits: pd.DataFrame) -> None:
    r0, r1 = RADII["peptide_density"]
    hits = plot.known(hits)
    arcs = []
    for c in plot.reps:
        n = plot.length(c)
        k = hits.loc[hits["contig"] == c, "start"].astype(int).to_numpy() // DENSITY_BP
        counts = np.bincount(k[(k >= 0) & (k * DENSITY_BP < n)], minlength=-(-n // DENSITY_BP))
        top = max(int(counts.max()) if len(counts) else 0, 1)
        for i, v in enumerate(counts.tolist()):
            a, b = i * DENSITY_BP, min((i + 1) * DENSITY_BP, n)
            plot.add("peptide_density", c, a, b, float(v), "placements")
            if v > 0:
                arcs.append(
                    (canvas.arc(c, a, b, r0, r0 + (r1 - r0) * v / top), style.OKABE_ITO["grey"])
                )
    canvas.rects("peptide_density", arcs, edgecolors="none", linewidths=0)


def _ticks(plot: _Plot, canvas: _Canvas, view: RunView) -> None:
    r0, r1 = RADII["ticks"]
    mid = (r0 + r1) / 2
    tus = view.table("tus_tsv")
    if tus is not None:
        segs = []
        for r in plot.known(tus).itertuples(index=False):
            s, e, c = int(r.start), int(r.end), str(r.contig)
            x = s if r.strand == "+" else e - 1
            plot.add("tu_start", c, x, x + 1, None, r.tu_id)
            segs.append((c, float(x if r.strand == "+" else x + 1), r0, r1))
        canvas.ticks("tu_start", segs, BLACK)
    for key, track, colour, (lo, hi) in (
        ("promoters_bed", "promoter", PROMOTER, (mid, r1)),
        ("terminators_bed", "terminator", TERMINATOR, (r0, mid)),
    ):
        bed = view.table(key)
        if bed is None:
            continue
        segs = []
        for r in plot.known(bed, "chrom").itertuples(index=False):
            c, strand = str(r.chrom), getattr(r, "strand", "+")
            plot.add(
                track,
                c,
                int(r.start),
                int(r.end),
                _num(getattr(r, "score", None)),
                getattr(r, "name", None),
            )
            segs.append((c, float(r.start if strand != "-" else r.end), lo, hi))
        canvas.ticks(track, segs, colour)


def _support(plot: _Plot, canvas: _Canvas, view: RunView) -> None:
    skew = view.support("gc_skew")
    if skew is None:
        plot.notes.append("GC skew omitted (no plots/support/gc_skew.tsv)")
    else:
        r0, r1 = RADII["gc_skew"]
        for c in plot.reps:
            w = skew[(skew["contig"] == c) & skew["cumulative_gc_skew"].notna()]
            if w.empty:
                continue
            y = w["cumulative_gc_skew"].astype(float).to_numpy()
            top = float(np.abs(y).max()) or 1.0
            x = []
            for s, e, v in zip(w["start"].astype(int), w["end"].astype(int), y, strict=True):
                plot.add("gc_skew", c, s - 1, e, float(v), "cumulative")
                x.append(canvas.rad(c, (s - 1 + e) / 2))
            rs = (r0 + r1) / 2 + (r1 - r0) / 2 * y / top
            canvas.ax.plot(x, rs, color=BLACK, linewidth=0.5, gid="gc_skew", zorder=3)
    marks = view.support("replicon_landmarks")
    if marks is None:
        plot.notes.append("oriC / ter omitted (no plots/support/replicon_landmarks.tsv)")
        return
    called = plot.known(marks)
    called = called[(called["status"] == "called") & (called["topology"] == "circular")]
    for r in called.itertuples(index=False):
        for track, pos, ls in (("oric", r.oric_position, "-"), ("ter", r.ter_position, ":")):
            if _na(pos):
                continue
            c, p = str(r.contig), int(pos)
            plot.add(track, c, p - 1, p, None, None)
            rad = canvas.rad(c, p - 1)
            canvas.ax.plot(
                [rad, rad],
                [RADII["novel_orf"][0] - 4, 100],
                color=BLACK,
                linestyle=ls,
                linewidth=0.6,
                gid=track,
            )
            canvas.ax.text(
                rad,
                RADII["novel_orf"][0] - 8,
                "oriC" if track == "oric" else "ter",
                ha="center",
                va="center",
                fontsize=6,
            )


def _markers(plot: _Plot, canvas: _Canvas, view: RunView) -> list[str]:
    counts = []
    variants = view.table("variants")
    if variants is not None:
        v = plot.known(variants[variants["call"].isin(["ALT", "mixed"])])
        counts.append(_plural(len(v), "ALT site"))
        pts: dict[str, list[tuple[str, float]]] = {"high": [], "candidate": []}
        for r in v.itertuples(index=False):
            s = int(r.codon_start) - 1
            plot.add("alt_site", str(r.contig), s, s + 3, None, r.confidence)
            pts["high" if r.confidence == "high" else "candidate"].append((str(r.contig), s))
        r0 = RADII["alt_site"][0]
        canvas.markers("alt_site", pts["high"], r0, marker="^", facecolors=BLACK, edgecolors=BLACK)
        canvas.markers(
            "alt_site", pts["candidate"], r0, marker="^", facecolors="none", edgecolors=BLACK
        )
    orfs = view.table("novel_orfs")
    if orfs is not None:
        o = plot.known(orfs[orfs["called"].astype("boolean").fillna(False).astype(bool)])
        counts.append(_plural(len(o), "novel ORF"))
        pts_o = []
        for r in o.itertuples(index=False):
            first = r.start_longest if not _na(r.start_longest) else r.start_best_ostir
            ends = [int(r.end)] + ([] if _na(first) else [int(first) + 1])  # 1-based
            lo, hi = min(ends), max(ends)
            plot.add("novel_orf", str(r.contig), lo - 1, hi, _num(r.n_peptides), r[1])
            pts_o.append((str(r.contig), (lo - 1 + hi) / 2))
        canvas.markers(
            "novel_orf",
            pts_o,
            RADII["novel_orf"][0],
            marker="o",
            facecolors=NOVEL,
            edgecolors=NOVEL,
        )
    return counts


def _legend(fig: Any, binned: bool, notes: Sequence[str]) -> None:
    def box(colour: str, label: str, hatch: str | None = None) -> tuple[Patch, str]:
        return Patch(
            facecolor=colour, edgecolor="#555555" if hatch else colour, hatch=hatch, linewidth=0
        ), label

    items = (
        []
        if binned
        else [
            box(c, k.replace("_", " "), PARTIAL_HATCH if k == "partial_excluded" else None)
            for k, c in style.STATUS_COLOUR.items()
        ]
    )
    items += [box(c, f"{k} arc") for k, c in style.REGION.items()]
    items += [box(c, f"island {k}", h) for k, (c, h) in style.ISLAND.items()]
    lines = [(BLACK, "TU start"), (PROMOTER, "promoter"), (TERMINATOR, "terminator")]
    handles: list[Artist] = [h for h, _ in items]
    handles += [Line2D([], [], color=c) for c, _ in lines]
    labels = [lab for _, lab in items] + [lab for _, lab in lines]
    for marker, face, edge, lab in (
        ("^", BLACK, BLACK, "ALT site (high)"),
        ("^", "none", BLACK, "ALT site (candidate)"),
        ("o", NOVEL, NOVEL, "novel ORF"),
    ):
        handles.append(
            Line2D(
                [], [], marker=marker, linestyle="none", markerfacecolor=face, markeredgecolor=edge
            )
        )
        labels.append(lab)
    fig.legend(
        handles, labels, loc="upper left", bbox_to_anchor=(0.70, 0.92), frameon=False, fontsize=6
    )
    key = "gene colour: detected share (viridis)" if binned else "gene colour: detection status"
    text = [key, "abundance: log2, viridis 2nd-98th pct; white = NA or imputed", *notes]
    fig.text(0.70, 0.04, "\n".join(text), ha="left", va="bottom", fontsize=6, color="#555555")


@base.plot(
    "V01",
    "circos",
    needs=("preflight", "genes_tsv"),
    optional=(
        "peptide_hits", "regions_tsv", "tus_tsv", "promoters_bed", "terminators_bed", "islands",
        "variants", "novel_orfs", "events",
    ),
)  # fmt: skip
def circos(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V1: the genome as a circle with every per-locus track of the run."""
    genes = view.table("genes_tsv")
    pre = view.preflight()
    assert genes is not None and pre is not None  # guaranteed by ``needs``
    plot = _Plot([r for r in pre["replicons"] if int(r["length"]) > 0])
    one = len(plot.reps) == 1 and plot.circular(next(iter(plot.reps)))
    layout = Circos({c: plot.length(c) for c in plot.reps}, space=0 if one else 2)
    step = tick_step(max(plot.length(c) for c in plot.reps))
    for sector in layout.sectors:
        c = sector.name
        plot.add("sector", c, 0, plot.length(c), plot.length(c), str(plot.reps[c].get("topology")))
        track = sector.add_track(RADII["sector"])
        track.axis(fc="#DDDDDD", ec="none")
        track.xticks_by_interval(
            step, label_size=6, tick_length=1.5, label_formatter=lambda x: tick_label(x, step)
        )
        sector.text(c, r=112, size=7)
        for name in ("region", "gene_plus", "gene_minus", "abundance"):
            sector.add_track(RADII[name]).axis(fc="#F4F4F4", ec="none")
    fig = style.new_figure("double", 6.0)
    fig.set_label("V01_circos")
    ax = fig.add_axes((0.02, 0.04, 0.64, 0.86), projection="polar")
    layout.plotfig(ax=ax)
    canvas = _Canvas(ax, {s.name: s for s in layout.sectors})
    linear = [
        (c, float(x), 60.0, 100.0)
        for c in plot.reps
        if not plot.circular(c)
        for x in (0, plot.length(c))
    ]
    canvas.ticks("sector", linear, BLACK)

    regions = view.table("regions_tsv")
    if regions is not None:
        regions = plot.known(regions)
        no_xy = regions["start"].isna() | regions["end"].isna()
        if (no_xy & (regions["region_type"] == "amr")).any():
            plot.notes.append("AMR sets without coordinates not drawn")
        reg = regions[~no_xy]
        reg = reg.assign(category=reg["region_type"], value=None)
        canvas.rects(
            "region",
            _intervals(
                plot, canvas, "region", reg, lambda r: style.REGION.get(r.category, "#999999")
            ),
            edgecolors="none",
            linewidths=0,
        )
    missing = sorted(f.split(":", 1)[1] for f in view.run_flags() if f.startswith("layer_missing:"))
    if missing:
        plot.notes.append(f"layers missing (their arcs absent): {', '.join(missing)}")
    binned = _genes(plot, canvas, genes)
    hits = view.table("peptide_hits")
    if hits is None:
        plot.notes.append("peptide density omitted (no peptide_hits.tsv)")
    else:
        _density(plot, canvas, hits)
    _ticks(plot, canvas, view)

    alpha = float(cfg.stats.alpha.context)
    islands = view.table("islands")
    n_isl = None
    if islands is not None:
        isl = plot.known(islands[islands["p_value"].astype(float) <= alpha])
        n_isl = len(isl)
        for direction, (colour, hatch) in style.ISLAND.items():
            sub = isl[isl["direction"] == direction]
            sub = sub.assign(category=direction, value=sub["p_value"])
            arcs = _intervals(plot, canvas, "island", sub, _fixed(colour + "40"))
            canvas.rects("island", arcs, hatch=hatch, edgecolors=colour, linewidths=0.4)
    _support(plot, canvas, view)
    markers = _markers(plot, canvas, view)

    d = int((genes["detection_status"] == "detected").sum())
    parts = [f"{d} of {len(genes)} genes detected ({base.title_pct(d, len(genes))})"]
    if n_isl is not None:
        parts.append(f"{_plural(n_isl, 'island')} at FWER ≤ {alpha:g}")
    if markers:
        parts.append(", ".join(markers))
    title = "; ".join(parts)
    fig.suptitle(title, y=0.98)
    _legend(fig, binned, plot.notes)
    return [base.Instance(None, fig, plot.frame(), title, coords=("start",))]
