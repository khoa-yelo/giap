"""Draw the GIPPA process-level diagram for the docs (light and dark, SVG and 300-dpi PNG).

Usage: python scripts/docs/process_figure.py [--outdir docs/assets] [--check]

Writes ``gippa-processes-{light,dark}.{svg,png}``: every Nextflow process of workflow/main.nf as a
card (the process name and its main outputs), grouped by subworkflow, with the channels between
them as edges, in the palette and themes of the overview figure (pipeline_figure.py).
Matplotlib only, with the bundled DejaVu fonts, text converted to paths, no SVG
date and a fixed ``svg.hashsalt``, so a re-run with the same Matplotlib version is byte-identical.
``--check`` writes nothing and exits 1 when a committed file differs, when a text label overlaps
another label or leaves its card, or when an edge runs through a card or a label.

The graph is the one ``nextflow run workflow/main.nf -preview -with-dag dag.mmd`` reports with
every module on (49 processes); the three TERMINATORNET_SHARDED processes are drawn as one node.
Edges follow the ``stepIn`` / ``pick`` wiring in workflow/subworkflows/local/*.nf. Inputs that
nearly every process reads (the genome, Bakta's gene model and the peptide table) are badges, not
edges.
"""

from __future__ import annotations

import argparse
import io
import itertools
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.path import Path as MPath
from matplotlib.text import Text

# the overview figure's palette and themes (the script directory is on sys.path when run)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_figure import DARK, FONT, HUE, LIGHT, Theme, _mix

NOVEL_LEFT = ["NOVEL_CANDIDATES", "PROMOTECH_NOVEL", "NOVEL_PROMOTER_LAYER", "OSTIR_NOVEL"]
REPO = Path(__file__).resolve().parents[2]
NAME = "gippa-processes"

MONO = "DejaVu Sans Mono"
PT = 1 / 72  # one point in inches; the axes use inches as data units
F_HEAD, F_NAME, F_OUT, F_SMALL = 9.2, 7.8, 6.3, 6.1

# badges: inputs read by most processes, drawn as a lettered dot instead of an edge
BADGE = {
    "G": ("io", "genome.norm.fa"),
    "A": ("ann", "Bakta's gene model or ANNOTATE's tables"),
    "P": ("map", "peptides.norm.parquet"),
}


# ---------------------------------------------------------------------------------------------
# layout grid: an input column, five process columns and an output column
W = 13.0
G = 0.37  # gutter between columns: room for the vertical channels
XI, WI = 0.12, 1.34
CW = 1.58
COLS = [XI + WI + G + i * (CW + G) for i in range(5)]  # left edges of c1..c5
XO = COLS[4] + CW + G
WO = W - 0.12 - XO


def cx(c: int) -> float:
    return COLS[c - 1] + CW / 2


def gut(c: int) -> float:
    """x of the gutter right of column ``c`` (0: right of the inputs column)."""
    return (XI + WI + G / 2) if c == 0 else COLS[c - 1] + CW + G / 2


def card_h(n: int) -> float:
    return 0.235 + 0.11 * n


@dataclass
class Node:
    name: str
    col: int
    y: float
    hue: str
    outs: list[str]
    tool: bool = False  # an external tool's process; else a `gippa _step` (GIPPA_STEP alias)
    opt: bool = False  # runs only with an optional module, input or switch
    badges: str = ""
    note: str = ""  # a small line under the name (a collapsed subworkflow)
    badges_left: bool = False  # badges at the top-left corner (an edge enters top-right)
    x: float = 0.0
    w: float = CW
    h: float = 0.0

    def __post_init__(self) -> None:
        self.x = COLS[self.col - 1]
        self.h = card_h(len(self.outs) + (1 if self.note else 0))

    @property
    def box(self) -> tuple[float, float, float, float]:
        return (self.x, self.y - self.h / 2, self.w, self.h)

    def l(self, dy: float = 0.0) -> tuple[float, float]:  # noqa: E743
        return (self.x, self.y + dy)

    def r(self, dy: float = 0.0) -> tuple[float, float]:
        return (self.x + self.w, self.y + dy)

    def t(self, dx: float = 0.0) -> tuple[float, float]:
        return (self.x + self.w / 2 + dx, self.y + self.h / 2)

    def b(self, dx: float = 0.0) -> tuple[float, float]:
        return (self.x + self.w / 2 + dx, self.y - self.h / 2)


@dataclass
class Edge:
    pts: list[tuple[float, float]]
    hue: str
    arrow: bool = True
    dashed: bool = False
    prio: int = 0  # at a crossing the lower priority line breaks
    ends: tuple[str, ...] = ()  # the nodes it may touch
    gaps: list[tuple[int, float, float]] = field(default_factory=list)


def path(*pts: tuple[float, float] | float, start_h: bool = True) -> list[tuple[float, float]]:
    """An orthogonal polyline: points, or bare floats that move alternately in x and in y."""
    out: list[tuple[float, float]] = []
    horiz = start_h
    for p in pts:
        if isinstance(p, tuple):
            out.append(p)
            continue
        x, y = out[-1]
        out.append((p, y) if horiz else (x, p))
        horiz = not horiz
    return out


# The processes: (name, column, row, outputs "a | b", flags). A fractional row sits between two
# rows. Flags: t = an external tool's process, o = optional, < = badges at the top-left corner;
# G, A, P = the badges (BADGE). Outputs are each process's main files (out_spec / output:).
ANNOTATE_NODES = [
    ("BAKTA", 1, 0, "<id>.json, <id>.gbff | tools/bakta/", "t"),
    (
        "PARSE_BAKTA",
        2,
        0,
        "genes.parquet, features.parquet | replicons.parquet, xrefs.parquet | trna_counts.tsv …",
        "G",
    ),
    ("CODON_USAGE", 3, 0, "codon.parquet", "A"),
    ("ANTISMASH", 1, 1, "antismash/", "to"),
    ("GENOMAD", 1, 2, "genomad/", "toG"),
    ("REGION_LAYER", 3, 1.5, "regions.parquet | region_genes.parquet", "A"),
    ("UNIOP_INPUT", 1, 3, "uniop_cds.gff | uniop_map.tsv", "A"),
    ("UNIOP", 2, 3, "uniop.pred", "t"),
    ("OPERON_PAIRS", 3, 3, "pairs_raw.parquet, tus.parquet | tu_upstream/, operon_meta.json", "GA"),
    ("PROMOTER_CALCULATOR", 3, 4, "*.csv (100 TUs a task)", "t"),
    ("TX_LAYER", 4, 4, "tx_rates.parquet", "A"),
    ("PROMOTER_WINDOWS", 1, 5, "windows/windows_*.fa | windows_map.parquet", "GA"),
    ("PROMOTECH", 2, 5, "promotech_NNN/", "t"),
    ("PROMOTER_LAYER", 3, 5, "promoter_windows.parquet", ""),
    ("TERMINATORNET_SHARDED", 2, 6, "terminatornet.tsv", "tG"),
    ("TRANSTERM_INPUT", 1, 7, "genes.coords", "A"),
    ("TRANSTERMHP", 2, 7, "transtermhp.tt", "tG"),
    ("TERMINATOR_LAYER", 3, 7, "terminators.parquet | terminators_crosscheck.tsv", "A"),
    ("BACTERMFINDER", 2, 8, "bactermfinder/", "toG"),
    ("OSTIR_INPUT", 1, 9, "ostir_input.csv", "GA"),
    ("OSTIR", 2, 9, "ostir.csv", "t"),
    ("RBS_LAYER", 3, 9, "rbs.parquet", ""),
    ("CONTEXT_LAYERS", 5, 4.5, "layers.parquet | annotation/", "GA"),
]
MAPPING_NODES = [
    ("SIXFRAME_INDEX", 1, 0, "sixframe.npz", "GA"),
    ("MAP_EXACT", 2, 0, "hits_exact.parquet", "P"),
    ("MAP_TOLERANT", 3, 0, "hits_tolerant.parquet", "P"),
    ("PLACEMENT", 4, 0, "placement.parquet | mapping_fdr.parquet", "PA"),
    ("MINIPROT", 1, 1, "miniprot.gff", "toG"),
    ("PROTEIN_LEFTOVERS", 2, 1, "proteins.leftover.faa", "oA"),
    ("TBLASTN", 3, 1, "tblastn.tsv", "toG"),
    ("LOCI", 5, 0.5, "mapping/peptide_hits.tsv | protein_hits.tsv, unmapped.tsv", "GPA"),
]
ANALYSIS_NODES = [
    ("GENE_TABLE", 1, 0, "genes_pre.parquet | tier.parquet", "PA"),
    ("DETECTABILITY", 2, 0, "genetable/genes.parquet | tu_folds.parquet", ""),
    ("VARIANTS_STEP", 1, 1, "tables/variants.tsv, .vcf | variants_qc.parquet", "GPA"),
    ("CONTEXT_STATS", 2, 1, "tables/genes.tsv, islands.tsv … | gene_context.parquet", "PA"),
    ("DETERMINANTS_STEP", 3, 1, "tables/determinants_*.tsv | determinants_support/", "A"),
    ("SENSITIVITY", 4, 1, "tables/sensitivity.tsv", "A<"),
    ("OPERONS_STEP", 2, 2, "tables/pairs.tsv, operons.tsv | operon_flags.tsv", "A"),
    ("FUNCTION_STEP", 3, 2, "tables/category_summary.tsv | enrichment_*.tsv", "oA"),
    (
        "NOVEL_CANDIDATES",
        1,
        3,
        "novel_candidates.parquet | novel_windows/*.fa | ostir_input_novel.csv",
        "GPA",
    ),
    ("PROMOTECH_NOVEL", 2, 3, "promotech_NNN/", "t"),
    ("NOVEL_PROMOTER_LAYER", 3, 3, "novel_promoter_windows.parquet", ""),
    (
        "NOVEL_ORFS",
        4,
        3.5,
        "tables/novel_orfs.tsv, .gff3 | novel_orfs.faa | annotation_patch.gff3",
        "GPA",
    ),
    ("OSTIR_NOVEL", 2, 4, "ostir.csv", "t"),
    ("NOVEL_RBS_LAYER", 3, 4, "rbs_novel.parquet", ""),
    ("PLOTS", 5, 2, "plots/", ""),
    ("REPORT_HTML", 5, 3, "report/report.html", ""),
]
NOTES = {"TERMINATORNET_SHARDED": "PLAN → n CHUNKs → MERGE"}


# ---------------------------------------------------------------------------------------------
class Fig:
    def __init__(self, theme: Theme, h: float) -> None:
        self.t = theme
        self.H = h
        # 300 dpi: the QA measures the text at the PNG resolution (no pixel-snapped font sizes)
        self.fig = plt.figure(figsize=(W, h), dpi=300, facecolor=theme.bg)
        ax = self.fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(0, W)
        ax.set_ylim(0, h)
        ax.set_axis_off()
        self.ax = ax
        self.labels: list[tuple[Text, tuple[float, float, float, float] | None]] = []
        self.boxes: dict[str, tuple[float, float, float, float]] = {}
        self.edges: list[Edge] = []
        self.discs: list[tuple[float, float, float, float]] = []  # badge dots, for the QA

    def fill(self, hue: str) -> str:
        return _mix(HUE[hue], self.t.bg, self.t.tint)

    def edge(self, hue: str) -> str:
        return _mix("#FFFFFF", HUE[hue], self.t.edge_mix)

    def text(
        self,
        x: float,
        y: float,
        s: str,
        *,
        size: float = F_OUT,
        weight: str = "normal",
        color: str | None = None,
        ha: str = "center",
        va: str = "center",
        family: str = FONT,
        box: tuple[float, float, float, float] | None = None,
        style: str = "normal",
        z: float = 6,
    ) -> Text:
        t = self.ax.text(
            x,
            y,
            s,
            fontsize=size,
            fontweight=weight,
            fontfamily=family,
            fontstyle=style,
            color=color or self.t.text,
            ha=ha,
            va=va,
            linespacing=1.15,
            zorder=z,
        )
        self.labels.append((t, box))
        return t

    def rbox(
        self,
        b: tuple[float, float, float, float],
        hue: str,
        *,
        filled: bool = True,
        lw: float = 0.8,
        r: float = 0.05,
        dashed: bool = False,
        z: float = 2,
    ) -> None:
        x, y, w, h = b
        self.ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle=f"round,pad=0,rounding_size={r}",
                facecolor=self.fill(hue) if filled else "none",
                edgecolor=self.edge(hue),
                linewidth=lw,
                linestyle=(0, (3, 2)) if dashed else "solid",
                zorder=z,
            )
        )

    def doc(
        self, b: tuple[float, float, float, float], hue: str = "io", *, dashed: bool = False
    ) -> None:
        """A file: a rectangle with a folded top-right corner."""
        x, y, w, h = b
        c = min(0.07, h * 0.3)
        pts = [(x, y), (x + w, y), (x + w, y + h - c), (x + w - c, y + h), (x, y + h)]
        ls = (0, (3, 2)) if dashed else "solid"
        kw = dict(
            edgecolor=self.edge(hue), linewidth=0.8, zorder=2, joinstyle="round", linestyle=ls
        )
        self.ax.add_patch(Polygon(pts, closed=True, facecolor=self.fill(hue), **kw))
        fold = [(x + w - c, y + h), (x + w - c, y + h - c), (x + w, y + h - c)]
        kw["linestyle"] = "solid"
        self.ax.add_patch(Polygon(fold, closed=False, facecolor="none", **kw))

    # cards ----------------------------------------------------------------------------------
    def card(self, n: Node) -> None:
        t = self.t
        b = n.box
        self.boxes[n.name] = b
        self.ax.add_patch(
            FancyBboxPatch(
                (b[0], b[1]),
                b[2],
                b[3],
                boxstyle="round,pad=0,rounding_size=0.06",
                facecolor=self.fill(n.hue) if n.tool else t.bg,
                edgecolor=self.edge(n.hue),
                linewidth=1.1 if n.tool else 0.85,
                linestyle=(0, (3, 1.6)) if n.opt else "solid",
                zorder=3,
            )
        )
        top = b[1] + b[3] - 0.065
        self.text(
            n.x + n.w / 2,
            top - F_NAME * PT * 0.86,
            n.name,
            size=F_NAME,
            weight="bold",
            family=MONO,
            va="baseline",
            box=b,
        )
        lines = ([n.note] if n.note else []) + n.outs
        for k, s in enumerate(lines):
            y = top - 0.172 - k * 0.11 - F_OUT * PT * 0.5
            italic = bool(n.note) and k == 0
            self.text(
                n.x + n.w / 2,
                y,
                s,
                size=F_OUT,
                color=t.sub,
                va="baseline",
                style="italic" if italic else "normal",
                box=b,
            )
        for k, ch in enumerate(reversed(n.badges)):
            bx = n.x + 0.11 + k * 0.15 if n.badges_left else n.x + n.w - 0.11 - k * 0.15
            self.badge(bx, b[1] + b[3] + 0.012, ch)

    def badge(self, x: float, y: float, ch: str) -> None:
        hue = BADGE[ch][0]
        self.ax.add_patch(
            Circle((x, y), 0.068, facecolor=self.edge(hue), edgecolor=self.t.bg, lw=0.9, zorder=7)
        )
        self.ax.text(
            x,
            y - 0.002,
            ch,
            fontsize=5.3,
            fontweight="bold",
            fontfamily=FONT,
            color=self.t.bg,
            ha="center",
            va="center",
            zorder=8,
        )
        self.discs.append((x - 0.068, y - 0.068, x + 0.068, y + 0.068))

    # edges ----------------------------------------------------------------------------------
    def add(self, e: Edge) -> None:
        self.edges.append(e)

    def _crossings(self) -> None:
        """Give the lower-priority line of every interior crossing a small gap (a line hop)."""
        segs = []
        for i, e in enumerate(self.edges):
            for k in range(len(e.pts) - 1):
                (x0, y0), (x1, y1) = e.pts[k], e.pts[k + 1]
                segs.append((i, k, x0, y0, x1, y1))
        eps = 0.03
        for a in range(len(segs)):
            i, k, ax0, ay0, ax1, ay1 = segs[a]
            for bb in range(a + 1, len(segs)):
                j, m, bx0, by0, bx1, by1 = segs[bb]
                if i == j:
                    continue
                a_h, b_h = ay0 == ay1, by0 == by1
                if a_h == b_h:
                    continue
                (hx0, hx1, hy), (vx, vy0, vy1) = (
                    ((min(ax0, ax1), max(ax0, ax1), ay0), (bx0, min(by0, by1), max(by0, by1)))
                    if a_h
                    else ((min(bx0, bx1), max(bx0, bx1), by0), (ax0, min(ay0, ay1), max(ay0, ay1)))
                )
                if not (hx0 + eps < vx < hx1 - eps and vy0 + eps < hy < vy1 - eps):
                    continue
                ei, ej = self.edges[i], self.edges[j]
                # the lower priority breaks; on a tie the vertical one
                if ei.prio != ej.prio:
                    brk, seg = (i, k) if ei.prio < ej.prio else (j, m)
                else:
                    brk, seg = (j, m) if a_h else (i, k)
                be = self.edges[brk]
                p, q = be.pts[seg], be.pts[seg + 1]
                along = vy0 if p[0] == q[0] else hx0
                pos = hy if p[0] == q[0] else vx
                be.gaps.append((seg, pos, along))

    def draw_edges(self) -> None:
        self._crossings()
        gap = 0.055
        for e in self.edges:
            col = self.edge(e.hue)
            ls = (0, (3, 2)) if e.dashed else "solid"
            pieces: list[list[tuple[float, float]]] = [[e.pts[0]]]
            for k in range(len(e.pts) - 1):
                p, q = e.pts[k], e.pts[k + 1]
                vert = p[0] == q[0]
                cuts = sorted(
                    (pos for seg, pos, _ in e.gaps if seg == k),
                    reverse=(q[1] < p[1]) if vert else (q[0] < p[0]),
                )
                sgn = (1 if q[1] > p[1] else -1) if vert else (1 if q[0] > p[0] else -1)
                for c in cuts:
                    if vert:
                        pieces[-1].append((p[0], c - sgn * gap))
                        pieces.append([(p[0], c + sgn * gap)])
                    else:
                        pieces[-1].append((c - sgn * gap, p[1]))
                        pieces.append([(c + sgn * gap, p[1])])
                pieces[-1].append(q)
            for n, pc in enumerate(pieces):
                last = n == len(pieces) - 1
                if last and e.arrow:
                    self.ax.add_patch(
                        FancyArrowPatch(
                            path=MPath(pc),
                            arrowstyle="-|>,head_length=3.0,head_width=1.8",
                            color=col,
                            linewidth=0.85,
                            linestyle=ls,
                            shrinkA=0,
                            shrinkB=0,
                            joinstyle="miter",
                            capstyle="butt",
                            zorder=4,
                        )
                    )
                else:
                    xs, ys = zip(*pc, strict=True)
                    self.ax.plot(
                        xs,
                        ys,
                        color=col,
                        lw=0.85,
                        ls=ls,
                        zorder=4,
                        solid_joinstyle="miter",
                        solid_capstyle="butt",
                        dash_capstyle="butt",
                    )

    def dot(self, x: float, y: float, hue: str) -> None:
        self.ax.add_patch(Circle((x, y), 0.028, facecolor=self.edge(hue), lw=0, zorder=5))

    # QA -------------------------------------------------------------------------------------
    def problems(self) -> list[str]:
        """Labels that overlap or leave their box; edges through a card or a label."""
        r = self.fig.canvas.get_renderer()
        inv = self.ax.transData.inverted()
        boxes = []
        for t, owner in self.labels:
            e = t.get_window_extent(r)
            (x0, y0), (x1, y1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
            boxes.append((t.get_text(), (x0, y0, x1, y1), owner))
        eps = 0.01
        out = []
        pad = 0.03
        for s, b in [(s, b) for s, b, _ in boxes]:
            for d in self.discs:
                if b[0] < d[2] and d[0] < b[2] and b[1] < d[3] and d[1] < b[3]:
                    out.append(f"badge over {s!r}")
        for s, (x0, y0, x1, y1), owner in boxes:
            if owner is not None:
                bx, by, bw, bh = owner
                if x0 < bx + pad or x1 > bx + bw - pad or y0 < by + 0.015 or y1 > by + bh:
                    out.append(f"outside its box: {s!r}")
            if x0 < 0.02 or x1 > W - 0.02 or y0 < 0.02 or y1 > self.H - 0.02:
                out.append(f"outside the figure: {s!r}")
        for i, (s, a, _) in enumerate(boxes):
            for s2, b, _ in boxes[i + 1 :]:
                if a[0] < b[2] and b[0] < a[2] and a[1] + eps < b[3] and b[1] + eps < a[3]:
                    out.append(f"overlap: {s!r} / {s2!r}")
        # edges: no segment may enter a card (other than at its own ends) or cross a label
        shrink = 0.012
        rects = [(n, (x, y, x + w, y + h)) for n, (x, y, w, h) in self.boxes.items()]
        rects += [(s, bb) for s, bb, _ in boxes]
        for e in self.edges:
            for k in range(len(e.pts) - 1):
                (x0, y0), (x1, y1) = e.pts[k], e.pts[k + 1]
                sx0, sx1, sy0, sy1 = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
                for name, (rx0, ry0, rx1, ry1) in rects:
                    if name in e.ends:
                        continue
                    if (
                        sx0 < rx1 - shrink
                        and rx0 + shrink < sx1 + 1e-9
                        and sy0 < ry1 - shrink
                        and ry0 + shrink < sy1 + 1e-9
                    ):
                        out.append(f"edge {e.ends} segment {k} crosses {name!r}")
        return out


# ---------------------------------------------------------------------------------------------
def rows(top: float, heights: list[float], gap: float) -> list[float]:
    """Centres of stacked rows of the given heights, ``gap`` apart, from ``top`` down."""
    ys, y = [], top
    for h in heights:
        ys.append(y - h / 2)
        y -= h + gap
    return ys


def draw(theme: Theme) -> Fig:
    H = 14.9
    f = Fig(theme, H)
    t = theme
    nodes: dict[str, Node] = {}

    def place(table: list[tuple[str, int, float, str, str]], ys: list[float], hue: str) -> None:
        for name, col, row, outs, fl in table:
            i = int(row)
            y = ys[i] if row == i else (ys[i] + ys[i + 1]) / 2
            hu = "io" if name in ("PLOTS", "REPORT_HTML") else hue
            if name in ("GENE_TABLE", "DETECTABILITY"):
                hu = "gt"
            nodes[name] = Node(
                name, col, y, hu, [o.strip() for o in outs.split("|")],
                tool="t" in fl, opt="o" in fl, badges="".join(c for c in fl if c in BADGE),
                note=NOTES.get(name, ""), badges_left="<" in fl,
            )  # fmt: skip

    def e(
        pts: list[tuple[float, float]],
        hue: str,
        a: str,
        b: str = "",
        **kw: object,
    ) -> None:
        f.add(Edge(pts, hue, ends=tuple(x for x in (a, b) if x), **kw))  # type: ignore[arg-type]

    # column headers ---------------------------------------------------------------------------
    y_head = H - 0.22
    for x, w, s in [(XI, WI, "Pipeline inputs"), (XO, WO, "Published outputs")]:
        f.text(x + w / 2, y_head, s, size=F_HEAD, weight="bold", va="baseline")
        f.ax.plot([x, x + w], [y_head - 0.08] * 2, color=t.sub, lw=0.8, zorder=3)
    f.text(
        (COLS[0] + COLS[4] + CW) / 2,
        y_head,
        "Nextflow processes, grouped by subworkflow",
        size=F_HEAD,
        weight="bold",
        va="baseline",
    )
    f.ax.plot([COLS[0], COLS[4] + CW], [y_head - 0.08] * 2, color=t.sub, lw=0.8, zorder=3)

    # ANNOTATE ---------------------------------------------------------------------------------
    gap_a = 0.17
    ha = [card_h(2), card_h(1), card_h(1), card_h(2), card_h(1), card_h(2), card_h(2)]
    ha += [card_h(2), card_h(1), card_h(1)]
    ha[0] = card_h(3)
    ann_top = H - 0.75
    a = rows(ann_top, ha, gap_a)
    A = "ann"
    place(ANNOTATE_NODES, a, A)
    ann_bot = a[9] - ha[9] / 2

    # pipeline inputs of ANNOTATE
    di = 0.40
    gx = gut(0)
    docs = [
        ("genome.norm.fa", "inputs.genome", a[0], False),
        ("replicons.tsv", "inputs.replicons", a[1] - 0.06, False),
        ("regions.gff3", "inputs.regions_gff3", a[2] - 0.14, True),
    ]
    for title, sub, y, dashed in docs:
        b = (XI, y - di / 2, WI, di)
        f.doc(b, dashed=dashed)
        f.text(XI + WI / 2, y + 0.06, title, size=F_NAME, weight="bold", family=MONO, box=b)
        f.text(XI + WI / 2, y - 0.08, sub, size=F_SMALL, color=t.sub, box=b)
    bk = nodes["BAKTA"]
    e([(XI + WI, a[0]), bk.l()], "io", "BAKTA")
    e(path((XI + WI, docs[1][2]), gx - 0.07, a[0] - 0.09, bk.x), "io", "BAKTA")
    e(path((XI + WI, docs[2][2]), gx + 0.07, a[0] - 0.15, bk.x), "io", "BAKTA", dashed=True)

    n = nodes
    # BAKTA -> PARSE_BAKTA, ANTISMASH; PARSE_BAKTA -> CODON_USAGE
    e([bk.r(), n["PARSE_BAKTA"].l()], A, "BAKTA", "PARSE_BAKTA")
    e([bk.b(), n["ANTISMASH"].t()], A, "BAKTA", "ANTISMASH", dashed=True)
    e([n["PARSE_BAKTA"].r(), n["CODON_USAGE"].l()], A, "PARSE_BAKTA", "CODON_USAGE")
    rl = n["REGION_LAYER"]
    xg2 = gut(2)
    e(path(n["ANTISMASH"].r(), xg2, rl.y + 0.07, rl.x), A, "ANTISMASH", "REGION_LAYER", dashed=True)
    e(path(n["GENOMAD"].r(), xg2, rl.y - 0.07, rl.x), A, "GENOMAD", "REGION_LAYER", dashed=True)
    # operons
    ui, un, op = n["UNIOP_INPUT"], n["UNIOP"], n["OPERON_PAIRS"]
    e([ui.r(), un.l()], A, "UNIOP_INPUT", "UNIOP")
    e([un.r(), op.l()], A, "UNIOP", "OPERON_PAIRS")
    yo = op.y + op.h / 2 + gap_a / 2
    e(path(ui.t(), yo, op.t()[0], op.t()[1], start_h=False), A, "UNIOP_INPUT", "OPERON_PAIRS")
    pc, tx = n["PROMOTER_CALCULATOR"], n["TX_LAYER"]
    e([op.b(), pc.t()], A, "OPERON_PAIRS", "PROMOTER_CALCULATOR")
    e([pc.r(), tx.l()], A, "PROMOTER_CALCULATOR", "TX_LAYER")
    # promoters
    pw, pt, pl = n["PROMOTER_WINDOWS"], n["PROMOTECH"], n["PROMOTER_LAYER"]
    e([pw.r(), pt.l()], A, "PROMOTER_WINDOWS", "PROMOTECH")
    e([pt.r(), pl.l()], A, "PROMOTECH", "PROMOTER_LAYER")
    yw = pw.y + pw.h / 2 + gap_a / 2
    e(
        path(pw.t(), yw, pl.t(-0.35)[0], pl.t()[1], start_h=False),
        A,
        "PROMOTER_WINDOWS",
        "PROMOTER_LAYER",
    )
    # terminators
    tn, ti, th, tl, bf = (
        n["TERMINATORNET_SHARDED"],
        n["TRANSTERM_INPUT"],
        n["TRANSTERMHP"],
        n["TERMINATOR_LAYER"],
        n["BACTERMFINDER"],
    )
    e([ti.r(), th.l()], A, "TRANSTERM_INPUT", "TRANSTERMHP")
    e([th.r(), tl.l()], A, "TRANSTERMHP", "TERMINATOR_LAYER")
    e(path(tn.r(), xg2, tl.y + 0.1, tl.x), A, "TERMINATORNET_SHARDED", "TERMINATOR_LAYER")
    e(path(bf.r(), xg2, tl.y - 0.1, tl.x), A, "BACTERMFINDER", "TERMINATOR_LAYER", dashed=True)
    # RBS
    oi, os_, rb = n["OSTIR_INPUT"], n["OSTIR"], n["RBS_LAYER"]
    e([oi.r(), os_.l()], A, "OSTIR_INPUT", "OSTIR")
    e([os_.r(), rb.l()], A, "OSTIR", "RBS_LAYER")
    yr = oi.y + oi.h / 2 + gap_a / 2
    e(path(oi.t(), yr, rb.t()[0], rb.t()[1], start_h=False), A, "OSTIR_INPUT", "RBS_LAYER")

    # the layers bus into CONTEXT_LAYERS
    cl = n["CONTEXT_LAYERS"]
    xb = gut(4)
    feeders = ["CODON_USAGE", "REGION_LAYER", "OPERON_PAIRS", "TX_LAYER", "PROMOTER_LAYER"]
    feeders += ["TERMINATOR_LAYER", "RBS_LAYER"]
    ys = [n[k].y for k in feeders]
    # BAKTA's tools/bakta/ goes to CONTEXT_LAYERS too (tool_dirs -> annotation/tools/bakta/):
    # over the top of row 0, into the head of the bus
    y_over = max(nodes[k].box[1] + nodes[k].box[3] for k in ("BAKTA", "PARSE_BAKTA")) + 0.14
    e(path(bk.t(), y_over, xb, start_h=False), A, "BAKTA", arrow=False)
    f.add(Edge([(xb, y_over), (xb, min(ys))], A, arrow=False, prio=2, ends=()))
    for k in feeders:
        e([n[k].r(), (xb, n[k].y)], A, k, arrow=False)
        f.dot(xb, n[k].y, A)
    e([(xb, cl.y), cl.l()], A, "CONTEXT_LAYERS", prio=1)
    f.dot(xb, cl.y, A)
    f.text(
        xb + 0.07,
        max(ys) + 0.02,
        "layers",
        size=F_SMALL,
        style="italic",
        color=f.edge(A),
        ha="left",
        va="bottom",
    )

    # ANNOTATE cluster
    clusters: list[tuple[str, list[str], str, bool]] = []
    clusters.append(("ANNOTATE", [k for k, v in nodes.items() if v.hue == A], A, False))

    # MAPPING ----------------------------------------------------------------------------------
    M = "map"
    gap_m = 0.20
    m_top = ann_bot - 0.62
    m = rows(m_top, [card_h(2), card_h(1)], gap_m)
    place(MAPPING_NODES, m, M)
    for title, sub, y, dashed in [
        ("peptides.norm.parquet", "inputs.peptides", m[0], False),
        ("proteins.long.faa", "inputs.proteins_long", m[1], True),
    ]:
        b = (XI, y - di / 2, WI, di)
        f.doc(b, dashed=dashed)
        f.text(XI + WI / 2, y + 0.06, title, size=F_SMALL + 0.4, weight="bold", family=MONO, box=b)
        f.text(XI + WI / 2, y - 0.08, sub, size=F_SMALL, color=t.sub, box=b)
    chain = ["SIXFRAME_INDEX", "MAP_EXACT", "MAP_TOLERANT", "PLACEMENT"]
    e([(XI + WI, m[0]), n["SIXFRAME_INDEX"].l()], "io", "SIXFRAME_INDEX")
    for p, q in itertools.pairwise(chain):
        e([n[p].r(), n[q].l()], M, p, q)
    lo = n["LOCI"]
    e(path(n["PLACEMENT"].r(), gut(4), lo.y + 0.08, lo.x), M, "PLACEMENT", "LOCI")
    e([(XI + WI, m[1]), n["MINIPROT"].l()], "io", "MINIPROT", dashed=True)
    prot = ["MINIPROT", "PROTEIN_LEFTOVERS", "TBLASTN"]
    for p, q in itertools.pairwise(prot):
        e([n[p].r(), n[q].l()], M, p, q, dashed=True)
    e(path(n["TBLASTN"].r(), gut(4), lo.y - 0.08, lo.x), M, "TBLASTN", "LOCI", dashed=True)
    # miniprot's hits go to LOCI themselves (only its failures go on to tblastn): under the row
    y_mb = m[1] - card_h(1) / 2 - 0.13
    mp_ = n["MINIPROT"]
    e(
        path(mp_.b(), y_mb, lo.x + 0.32, lo.b()[1], start_h=False),
        M,
        "MINIPROT",
        "LOCI",
        dashed=True,
    )
    clusters.append(("MAPPING", [k for k, v in nodes.items() if v.hue == M], M, False))
    map_bot = y_mb - 0.04

    # outputs of ANNOTATE and MAPPING
    do = 0.44

    def out_doc(y: float, title: str, sub: str, src: Node | None) -> None:
        b = (XO, y - do / 2, WO, do)
        f.doc(b)
        f.text(XO + WO / 2, y + 0.07, title, size=F_NAME, weight="bold", family=MONO, box=b)
        f.text(XO + WO / 2, y - 0.08, sub, size=F_SMALL, color=t.sub, box=b)
        if src is not None:
            e([src.r(), (XO, y)], src.hue, src.name)

    out_doc(cl.y, "annotation/", "genes.gff3, regions, layers", cl)
    out_doc(lo.y, "mapping/", "peptide_hits.tsv …", lo)

    # GENETABLE and the analyses ---------------------------------------------------------------
    gap_d = 0.50
    d_top = map_bot - 0.95
    hd = [card_h(2), card_h(2), card_h(2), card_h(3), card_h(1)]
    d = rows(d_top, hd, gap_d)
    GT, AN = "gt", "ana"
    place(ANALYSIS_NODES, d, AN)

    yl = d_top + 0.48
    x_or = gut(0)
    gt_, va, cs, opn = n["GENE_TABLE"], n["VARIANTS_STEP"], n["CONTEXT_STATS"], n["OPERONS_STEP"]
    nc, de, se, fu, no = (
        n["NOVEL_CANDIDATES"],
        n["DETECTABILITY"],
        n["SENSITIVITY"],
        n["FUNCTION_STEP"],
        n["NOVEL_ORFS"],
    )
    dm = n["DETERMINANTS_STEP"]
    # LOCI's mapping/: down from LOCI, left above GENETABLE, then down the inputs side to
    # GENE_TABLE, VARIANTS, OPERONS and NOVEL, and down the c1|c2 gutter to CONTEXT
    f.add(Edge([lo.b(), (lo.b()[0], yl), (x_or, yl), (x_or, nc.y)], M, arrow=False, prio=3))
    for tgt in (gt_, va, opn, nc):
        e([(x_or, tgt.y), tgt.l()], M, tgt.name, prio=1)
        f.dot(x_or, tgt.y, M)
    x_cm = gut(1)
    e(path((x_cm, yl), cs.y + 0.08, cs.x, start_h=False), M, "CONTEXT_STATS")
    f.dot(x_cm, yl, M)
    f.text(
        lo.b()[0] - 0.12,
        yl + 0.03,
        "mapping/ (LOCI) → GENE_TABLE, VARIANTS, CONTEXT, OPERONS, NOVEL",
        size=F_SMALL,
        style="italic",
        color=f.edge(M),
        va="bottom",
        ha="right",
    )

    # GENE_TABLE -> DETECTABILITY; the gene table to every analysis
    e([gt_.r(), de.l()], GT, "GENE_TABLE", "DETECTABILITY", prio=4)
    e([de.b(), cs.t()], GT, "DETECTABILITY", "CONTEXT_STATS", prio=1)
    xs_ = se.t(0.6)[0]  # right of the SENSITIVITY_RERUNS label
    f.add(Edge([de.r(), (xs_, de.y)], GT, arrow=False, prio=2, ends=("DETECTABILITY",)))
    xd = dm.t(0.3)[0]
    e([(xd, de.y), (xd, dm.t()[1])], GT, "DETERMINANTS_STEP", prio=1)
    f.dot(xd, de.y, GT)
    e([(xs_, de.y), se.t(0.6)], GT, "SENSITIVITY", prio=1)
    # the c2|c3 trunk: OPERONS (left), FUNCTION (right), then over NOVEL
    x_gt = gut(2)
    yn = fu.y - fu.h / 2 - gap_d / 2  # above NOVEL
    f.add(Edge([(x_gt, de.y), (x_gt, yn)], GT, arrow=False, prio=2))
    f.dot(x_gt, de.y, GT)
    e([(x_gt, opn.y), opn.r()], GT, "OPERONS_STEP", prio=1)
    e([(x_gt, fu.y), fu.l()], GT, "FUNCTION_STEP", prio=1)
    f.dot(x_gt, opn.y, GT)
    e(path((x_gt, yn), nc.t()[0], nc.t()[1], start_h=True), GT, "NOVEL_CANDIDATES", prio=1)
    e(path((x_gt, yn), no.t()[0], no.t()[1], start_h=True), GT, "NOVEL_ORFS", prio=1)
    f.text(
        (xd + xs_) / 2,
        de.y + 0.04,
        "genetable/genes.parquet",
        size=F_SMALL,
        style="italic",
        color=f.edge(GT),
        va="bottom",
    )

    # the analyses' own channels
    e([va.r(-0.08), cs.l(-0.08)], AN, "VARIANTS_STEP", "CONTEXT_STATS", prio=3)
    e([cs.r(), dm.l()], AN, "CONTEXT_STATS", "DETERMINANTS_STEP", prio=3)
    e([dm.r(), se.l()], AN, "DETERMINANTS_STEP", "SENSITIVITY", prio=3)
    # NOVEL: NOVEL_CANDIDATES feeds both tools, both layer steps and NOVEL_ORFS
    pn, npl, on, nrl = (
        n["PROMOTECH_NOVEL"],
        n["NOVEL_PROMOTER_LAYER"],
        n["OSTIR_NOVEL"],
        n["NOVEL_RBS_LAYER"],
    )
    e([nc.r(0.08), pn.l(0.08)], AN, "NOVEL_CANDIDATES", "PROMOTECH_NOVEL")
    e([pn.r(), npl.l()], AN, "PROMOTECH_NOVEL", "NOVEL_PROMOTER_LAYER")
    x1 = gut(1)
    e(path(nc.r(-0.12), x1, on.y, on.x), AN, "NOVEL_CANDIDATES", "OSTIR_NOVEL")
    e([on.r(), nrl.l()], AN, "OSTIR_NOVEL", "NOVEL_RBS_LAYER")
    e([(x1, no.y), no.l()], AN, "NOVEL_ORFS", prio=1)
    f.dot(x1, no.y, AN)
    xn = npl.b(-0.35)[0]
    e([(xn, no.y), npl.b(-0.35)], AN, "NOVEL_PROMOTER_LAYER", prio=1)
    e([(xn, no.y), nrl.t(-0.35)], AN, "NOVEL_RBS_LAYER", prio=1)
    f.dot(xn, no.y, AN)
    e(path(npl.r(), gut(3), no.y + 0.12, no.x), AN, "NOVEL_PROMOTER_LAYER", "NOVEL_ORFS")
    e(path(nrl.r(), gut(3), no.y - 0.12, no.x), AN, "NOVEL_RBS_LAYER", "NOVEL_ORFS")
    # REPORT
    pp, rh = n["PLOTS"], n["REPORT_HTML"]
    e([pp.b(), rh.t()], "io", "PLOTS", "REPORT_HTML")
    xp = gut(4)
    f.add(Edge([(xp, pp.y + 0.35), (xp, rh.y + 0.06)], "io", arrow=False, prio=2))
    e([(xp, pp.y), pp.l()], "io", "PLOTS", prio=1)
    e([(xp, rh.y + 0.06), rh.l(0.06)], "io", "REPORT_HTML", prio=1)
    f.dot(xp, pp.y, "io")
    f.text(
        xp - 0.05,
        pp.y + 0.35,
        "pub: every published\n[key, path] pair",
        size=F_SMALL,
        style="italic",
        color=t.sub,
        ha="right",
        va="top",
    )

    clusters += [
        ("GENETABLE", ["GENE_TABLE", "DETECTABILITY"], GT, False),
        ("VARIANTS", ["VARIANTS_STEP"], AN, False),
        ("CONTEXT", ["CONTEXT_STATS"], AN, False),
        ("DETERMINANTS", ["DETERMINANTS_STEP"], AN, False),
        ("SENSITIVITY_RERUNS", ["SENSITIVITY"], AN, False),
        ("OPERONS", ["OPERONS_STEP"], AN, False),
        ("FUNCTION", ["FUNCTION_STEP"], AN, True),
        (
            "NOVEL",
            [*NOVEL_LEFT, "NOVEL_RBS_LAYER", "NOVEL_ORFS"],
            AN,
            False,
        ),
        ("REPORT", ["PLOTS", "REPORT_HTML"], "io", False),
    ]

    # optional inputs of the analyses (no edges: each lists its readers)
    opt_in = [
        ("reference_peptides", "PLACEMENT, VARIANTS"),
        ("peptide_spectra.parquet", "VARIANTS_STEP"),
        ("contaminants", "LOCI, NOVEL_CANDIDATES"),
        ("abundance.norm.parquet", "GENE_TABLE, CONTEXT"),
        ("replicon_depth", "CONTEXT_STATS"),
        ("preflight.json", "GENE_TABLE, REPORT"),
    ]
    oy = rows(d[0] + hd[0] / 2, [0.36] * len(opt_in), 0.10)
    for (title, sub), y in zip(opt_in, oy, strict=True):
        b = (XI, y - 0.18, WI - 0.12, 0.36)
        f.doc(b, dashed=title != "preflight.json")
        f.text(b[0] + b[2] / 2, y + 0.06, title, size=F_SMALL, weight="bold", family=MONO, box=b)
        f.text(b[0] + b[2] / 2, y - 0.08, "→ " + sub, size=F_SMALL - 0.4, color=t.sub, box=b)

    # outputs of the analyses
    out_doc(
        (d[1] + d[2]) / 2 + 0.1,
        "tables/",
        "the tables/… outputs\nof the analysis cards",
        None,
    )
    out_doc(pp.y, "plots/", "SVG, PNG, PDF + data", pp)
    out_doc(rh.y, "report/", "report.html", rh)

    # draw -------------------------------------------------------------------------------------
    for name, members, hue, dashed in clusters:
        bs = [nodes[k].box for k in members]
        x0 = min(b[0] for b in bs) - 0.09
        y0 = min(b[1] for b in bs) - 0.09
        x1 = max(b[0] + b[2] for b in bs) + 0.09
        y1 = max(b[1] + b[3] for b in bs) + (0.30 if name in ("ANNOTATE", "MAPPING") else 0.15)
        if name == "ANNOTATE":
            x1 = COLS[4] + CW + 0.09
        if name == "MAPPING":
            x1 = COLS[4] + CW + 0.09
            y0 = map_bot - 0.06
        f.rbox((x0, y0, x1 - x0, y1 - y0), hue, filled=False, lw=0.9, r=0.08, dashed=dashed, z=1)
        tt = f.text(
            x0 + 0.12,
            y1,
            f" {name} ",
            size=F_SMALL + 1.0,
            weight="bold",
            color=f.edge(hue),
            ha="left",
            family=MONO,
        )
        tt.set_bbox(dict(facecolor=t.bg, edgecolor="none", pad=0.4))
    for nd in nodes.values():
        f.card(nd)
    f.draw_edges()

    # legend and notes -------------------------------------------------------------------------
    ly = min(b[1] for k, b in f.boxes.items() if k.startswith("NOVEL")) - 0.09 - 0.40
    lx = XI
    sw = 0.42  # swatch width; row 1 the card styles, row 2 the badges below their swatches
    xs = [lx, lx + 3.0, lx + 6.4]
    f.rbox((xs[0], ly - 0.13, sw, 0.26), "ann", r=0.05, lw=1.1, z=3)
    f.ax.add_patch(
        FancyBboxPatch(
            (xs[1], ly - 0.13),
            sw,
            0.26,
            boxstyle="round,pad=0,rounding_size=0.05",
            facecolor=t.bg,
            edgecolor=f.edge("ann"),
            lw=0.85,
            zorder=3,
        )
    )
    f.rbox((xs[2], ly - 0.13, sw, 0.26), "io", filled=False, dashed=True, z=3)
    styles = [
        "an external tool's process (container)",
        "a gippa _step process (GIPPA_STEP alias)",
        "optional: an input, a module (all on by default) or an off-by-default switch",
    ]
    for x, s in zip(xs, styles, strict=True):
        f.text(x + sw + 0.1, ly, s, size=F_SMALL, ha="left")
    yb = ly - 0.36
    for x, (ch, (_hue, label)) in zip(xs, BADGE.items(), strict=True):
        f.badge(x + sw / 2, yb, ch)
        f.text(x + sw + 0.1, yb, "reads " + label, size=F_SMALL, ha="left")
    note = (
        "From nextflow run workflow/main.nf -preview -with-dag with every module on: 49 processes; "
        "TERMINATORNET_SHARDED's PLAN, CHUNK and MERGE are one card "
        "(params.terminatornet_sharded = false runs one TERMINATORNET).\n"
        "Edges are each process's main input channels. Not drawn: database paths, step event "
        "logs and some second readers (e.g. of sixframe.npz and the hits tables). Without "
        "peptides only ANNOTATE and REPORT run."
    )
    f.text(XI, yb - 0.27, note, size=F_SMALL, color=t.sub, ha="left", va="top")
    return f


def render(theme: Theme) -> tuple[bytes, bytes, list[str]]:
    plt.rcParams["svg.hashsalt"] = "gippa-processes"
    plt.rcParams["svg.fonttype"] = "path"
    plt.rcParams["font.family"] = FONT
    f = draw(theme)
    f.fig.canvas.draw()
    bad = f.problems()
    svg, png = io.BytesIO(), io.BytesIO()
    f.fig.savefig(svg, format="svg", facecolor=theme.bg, metadata={"Date": None})
    f.fig.savefig(png, format="png", dpi=300, facecolor=theme.bg, metadata={"Software": None})
    plt.close(f.fig)
    return svg.getvalue(), png.getvalue(), bad


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--outdir", type=Path, default=REPO / "docs" / "assets")
    ap.add_argument("--check", action="store_true", help="compare with the files; write nothing")
    args = ap.parse_args(argv)
    status = 0
    for theme in (LIGHT, DARK):
        svg, png, bad = render(theme)
        for msg in bad:
            print(f"{theme.name}: {msg}", file=sys.stderr)
            status = 1
        for ext, data in (("svg", svg), ("png", png)):
            path_ = args.outdir / f"{NAME}-{theme.name}.{ext}"
            if args.check:
                if not path_.is_file() or path_.read_bytes() != data:
                    print(f"out of date: {path_}", file=sys.stderr)
                    status = 1
            else:
                path_.parent.mkdir(parents=True, exist_ok=True)
                path_.write_bytes(data)
                print(path_)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
