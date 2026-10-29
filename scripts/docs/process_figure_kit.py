"""Building blocks of the GIPPA process-level diagram (scripts/docs/process_figure.py).

The grid, the ``Node`` / ``Edge`` model, the node tables of every subworkflow and ``Fig``, which
draws cards, edges, badges and labels and records them for the ``--check`` QA. process_figure.py
lays the nodes out and writes the files.
"""

from __future__ import annotations

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
from pipeline_figure import FONT, HUE, Theme, _mix

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
