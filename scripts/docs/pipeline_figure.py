"""Draw the GIPPA pipeline schematic for the README (light and dark, SVG and 300-dpi PNG).

Usage: python scripts/docs/pipeline_figure.py [--outdir docs/assets] [--check]

Writes ``gippa-pipeline-{light,dark}.{svg,png}``. Matplotlib only (patches, arrows, text), with
the bundled DejaVu Sans, text converted to paths, no SVG date and a fixed ``svg.hashsalt``, so a
re-run with the same Matplotlib version is byte-identical. ``--check`` writes nothing and exits 1
when a committed file differs, or when a text label overlaps another label or leaves its box.

The boxes follow workflow/main.nf: ANNOTATE and MAPPING after preflight, GENETABLE joining them
on the Bakta locus tags, CONTEXT, OPERONS, FUNCTION and DETERMINANTS reading the gene table,
VARIANTS and NOVEL reading the placements, then REPORT on every published table.
"""

from __future__ import annotations

import argparse
import io
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.text import Text

REPO = Path(__file__).resolve().parents[2]
NAME = "gippa-pipeline"

W, H = 8.0, 4.55  # inches; the axes use inches as data units
PT = 1 / 72  # one point in inches
FONT = "DejaVu Sans"
F_HEAD, F_MAIN, F_SUB = 8.2, 7.5, 6.4  # stage headers, box labels, secondary text

# Okabe-Ito hues (colour-blind safe): annotation, mapping, gene table, analysis, files
HUE = {
    "ann": "#0072B2",
    "map": "#D55E00",
    "gt": "#009E73",
    "ana": "#CC79A7",
    "io": "#6E7781",
}


@dataclass(frozen=True)
class Theme:
    name: str
    bg: str
    text: str
    sub: str
    tint: float  # share of the hue in a box fill
    edge_mix: float  # share of white mixed into the hue for edges and headers
    arrow: str


LIGHT = Theme("light", "#FFFFFF", "#1F2328", "#57606A", 0.11, 0.0, "#57606A")
DARK = Theme("dark", "#0D1117", "#E6EDF3", "#9DA7B3", 0.22, 0.35, "#8B949E")


def _mix(a: str, b: str, t: float) -> str:
    """``t`` of colour ``a`` over colour ``b``."""
    ra, rb = to_rgb(a), to_rgb(b)
    return to_hex(tuple(t * x + (1 - t) * y for x, y in zip(ra, rb, strict=True)))


class Fig:
    def __init__(self, theme: Theme) -> None:
        self.t = theme
        self.fig = plt.figure(figsize=(W, H), facecolor=theme.bg)
        ax = self.fig.add_axes((0, 0, 1, 1))
        ax.set_xlim(0, W)
        ax.set_ylim(0, H)
        ax.set_axis_off()
        self.ax = ax
        self.labels: list[tuple[Text, tuple[float, float, float, float] | None]] = []

    # colours -------------------------------------------------------------------------------
    def fill(self, hue: str) -> str:
        return _mix(HUE[hue], self.t.bg, self.t.tint)

    def edge(self, hue: str) -> str:
        return _mix("#FFFFFF", HUE[hue], self.t.edge_mix)

    # primitives ----------------------------------------------------------------------------
    def text(
        self,
        x: float,
        y: float,
        s: str,
        *,
        size: float = F_MAIN,
        weight: str = "normal",
        color: str | None = None,
        ha: str = "center",
        va: str = "center",
        box: tuple[float, float, float, float] | None = None,
        style: str = "normal",
    ) -> Text:
        t = self.ax.text(
            x,
            y,
            s,
            fontsize=size,
            fontweight=weight,
            fontfamily=FONT,
            fontstyle=style,
            color=color or self.t.text,
            ha=ha,
            va=va,
            linespacing=1.15,
            zorder=5,
        )
        self.labels.append((t, box))
        return t

    def rbox(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        hue: str,
        *,
        filled: bool = True,
        lw: float = 0.8,
        r: float = 0.05,
        dashed: bool = False,
        z: float = 2,
    ) -> tuple[float, float, float, float]:
        self.ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle=f"round,pad=0,rounding_size={r}",
                facecolor=self.fill(hue) if filled else self.t.bg,
                edgecolor=self.edge(hue),
                linewidth=lw,
                linestyle=(0, (3, 2)) if dashed else "solid",
                zorder=z,
            )
        )
        return (x, y, w, h)

    def doc(
        self, x: float, y: float, w: float, h: float, hue: str = "io", *, dashed: bool = False
    ) -> tuple[float, float, float, float]:
        """A file: a rectangle with a folded top-right corner."""
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
        return (x, y, w, h)

    def arrow(
        self,
        p: tuple[float, float],
        q: tuple[float, float],
        *,
        color: str | None = None,
        conn: str = "arc3",
        lw: float = 0.9,
    ) -> None:
        self.ax.add_patch(
            FancyArrowPatch(
                p,
                q,
                arrowstyle="-|>,head_length=3.2,head_width=1.9",
                connectionstyle=conn,
                color=color or self.t.arrow,
                linewidth=lw,
                shrinkA=0,
                shrinkB=0,
                joinstyle="miter",
                capstyle="butt",
                zorder=4,
            )
        )

    def labelled(
        self,
        b: tuple[float, float, float, float],
        title: str,
        sub: str = "",
        *,
        weight: str = "bold",
        size: float = F_MAIN,
    ) -> None:
        """A title and optional secondary text centred in box ``b``."""
        x, y, w, h = b
        cx, cy = x + w / 2, y + h / 2
        if not sub:
            self.text(cx, cy, title, size=size, weight=weight, box=b)
            return
        subs = sub.split("\n")
        lh_main, lh_sub, sep = size * PT * 1.2, F_SUB * PT * 1.2, 0.02
        top = cy + (lh_main + sep + len(subs) * lh_sub) / 2
        # every line on its own baseline, so titles in neighbouring boxes line up exactly
        self.text(cx, top - size * PT * 0.95, title, size=size, weight=weight, va="baseline", box=b)
        for k, line in enumerate(subs):
            y0 = top - lh_main - sep - k * lh_sub - F_SUB * PT * 0.95
            self.text(cx, y0, line, size=F_SUB, color=self.t.sub, va="baseline", box=b)

    # QA -------------------------------------------------------------------------------------
    def problems(self) -> list[str]:
        """Labels that overlap each other or stick out of their box."""
        r = self.fig.canvas.get_renderer()
        inv = self.ax.transData.inverted()
        boxes = []
        for t, owner in self.labels:
            e = t.get_window_extent(r)
            (x0, y0), (x1, y1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
            boxes.append((t.get_text(), (x0, y0, x1, y1), owner))
        eps = 0.012  # line boxes of stacked lines touch by their leading
        out = []
        pad = 0.025
        for s, (x0, y0, x1, y1), owner in boxes:
            if owner is not None:
                bx, by, bw, bh = owner
                if x0 < bx + pad or x1 > bx + bw - pad or y0 < by + pad * 0.6 or y1 > by + bh:
                    out.append(f"outside its box: {s!r}")
            if x0 < 0.02 or x1 > W - 0.02 or y0 < 0.02 or y1 > H - 0.02:
                out.append(f"outside the figure: {s!r}")
        for i, (s, a, _) in enumerate(boxes):
            for s2, b, _ in boxes[i + 1 :]:
                if a[0] < b[2] and b[0] < a[2] and a[1] + eps < b[3] and b[1] + eps < a[3]:
                    out.append(f"overlap: {s!r} / {s2!r}")
        return out


def draw(theme: Theme) -> Fig:
    f = Fig(theme)
    t = theme
    gap = 0.16  # between columns: room for an arrow

    # columns (x, width)
    xi, wi = 0.10, 1.00  # inputs
    xp, wp = xi + wi + gap, 0.62  # preflight
    xt, wt = xp + wp + gap, 2.32  # annotation and mapping tracks
    xg, wg = xt + wt + gap, 0.98  # gene table
    xa, wa = xg + wg + gap, 1.10  # analysis modules
    xo = xa + wa + gap  # outputs
    wo = W - 0.10 - xo

    y_head = H - 0.16
    top = H - 0.40  # top of the content area
    band_y, band_h = 0.10, 0.46  # error-control band
    bot = band_y + band_h + 0.16  # bottom of the content area

    for x, w, s in [
        (xi, wi, "Inputs"),
        (xp, wp, "Check"),
        (xt, wt, "Annotate and map"),
        (xg, wg, "Gene table"),
        (xa, wa, "Analyse"),
        (xo, wo, "Outputs"),
    ]:
        f.text(x + w / 2, y_head - 0.05, s, size=F_HEAD, weight="bold", color=t.text, va="baseline")
        f.ax.plot([x, x + w], [y_head - 0.12] * 2, color=t.sub, lw=0.8, zorder=3)

    # inputs ---------------------------------------------------------------------------------
    span = top - bot
    dh, cvh = 0.44, 0.22
    dg = (span - 5 * dh - cvh) / 5  # five documents and the convert pill, evenly spaced
    g = f.doc(xi, top - dh, wi, dh)
    f.labelled(g, "Genome", "FASTA; optional\nannotation")
    s = f.doc(xi, g[1] - dg - dh, wi, dh)
    f.labelled(s, "Search results", "DIA-NN, Cascadia\nor generic")
    cv = f.rbox(xi + 0.10, s[1] - dg - cvh, wi - 0.20, cvh, "io", r=0.11)
    f.text(xi + wi / 2, cv[1] + cv[3] / 2, "gippa convert", size=F_SUB, weight="bold", box=cv)
    p = f.doc(xi, cv[1] - dg - dh, wi, dh)
    f.labelled(p, "Peptide table", "peptides.tsv")
    o = f.doc(xi, p[1] - dg - dh, wi, dh, dashed=True)
    f.labelled(o, "Optional", "proteins.faa,\nabundance.tsv")
    c = f.doc(xi, bot, wi, dh)
    f.labelled(c, "Config", "config.yaml")
    f.arrow((xi + wi / 2, s[1]), (xi + wi / 2, cv[1] + cv[3]))
    f.arrow((xi + wi / 2, cv[1]), (xi + wi / 2, p[1] + p[3]))

    # preflight ------------------------------------------------------------------------------
    pf = f.rbox(xp, bot, wp, span, "io")
    f.labelled(pf, "Preflight", "validate,\nnormalise,\nfail early")
    for b in (g, p, o, c):
        y = b[1] + b[3] / 2
        f.arrow((xi + wi, y), (xp, y))

    # annotation track -----------------------------------------------------------------------
    pad, hdr, cg = 0.08, 0.27, 0.06
    cw = (wt - 2 * pad - cg) / 2
    chh = 0.40
    ah = hdr + 3 * chh + 2 * cg + pad
    ann = f.rbox(xt, top - ah, wt, ah, "ann", lw=0.9)
    f.text(
        xt + pad, top - 0.135, "Annotation", weight="bold", color=f.edge("ann"), ha="left", box=ann
    )
    f.text(
        xt + wt - pad, top - 0.135, "Bakta locus tags", size=F_SUB, color=t.sub, ha="right", box=ann
    )
    chips = [
        ("Gene model", "Bakta", False),
        ("Operons, TUs", "UniOP", False),
        ("Promoters", "Promotech,\nPromoter Calculator", False),
        ("Terminators", "TerminatorNet,\nTransTermHP", False),
        ("RBS, codons", "OSTIR, codon-bias", False),
        ("Regions", "Bakta, geNOMAD,\nantiSMASH", False),
    ]
    for i, (a, b, opt) in enumerate(chips):
        col, row = i % 2, i // 2
        cx = xt + pad + col * (cw + cg)
        cy = top - hdr - (row + 1) * chh - row * cg
        ch = f.rbox(cx, cy, cw, chh, "ann", filled=False, lw=0.6, r=0.03, z=3, dashed=opt)
        f.labelled(ch, a, b)

    # mapping track --------------------------------------------------------------------------
    my = bot + 0.22  # a lane below carries placements to variants and novel ORFs
    mtop = top - ah - 0.20
    mp = f.rbox(xt, my, wt, mtop - my, "map", lw=0.9)
    f.text(xt + pad, mtop - 0.135, "Mapping", weight="bold", color=f.edge("map"), ha="left", box=mp)
    f.text(
        xt + wt - pad, mtop - 0.135, "six-frame index", size=F_SUB, color=t.sub, ha="right", box=mp
    )
    mg = 0.18  # room for the arrows between the two chip columns
    mw = (wt - 2 * pad - mg) / 2
    rh = (mtop - hdr - my - pad - 0.15) / 2
    x1, x2 = xt + pad, xt + pad + mw + mg
    y1, y2 = mtop - hdr - rh, my + pad
    kw = dict(filled=False, lw=0.6, r=0.03, z=3)
    pe = f.rbox(x1, y1, mw, rh, "map", **kw)
    f.labelled(pe, "Peptides", "exact, tolerant")
    pl = f.rbox(x2, y1, mw, rh, "map", lw=0.6, r=0.03, z=3)
    f.labelled(pl, "Placement FDR", "genome decoys,\nq per stratum")
    pr = f.rbox(x1, y2, mw, rh, "map", dashed=True, **kw)
    f.labelled(pr, "Proteins", "miniprot, then\ntblastn")
    lo = f.rbox(x2, y2, mw, rh, "map", **kw)
    f.labelled(lo, "Loci", "genes; unmapped,\ncontaminants")
    ac = f.edge("map")
    f.arrow((x1 + mw, y1 + rh / 2), (x2, y1 + rh / 2), color=ac, lw=0.7)
    f.arrow((x1 + mw, y2 + rh / 2), (x2, y2 + rh / 2), color=ac, lw=0.7)
    f.arrow((x2 + mw / 2, y1), (x2 + mw / 2, y2 + rh), color=ac, lw=0.7)
    gx = xt + wt / 2  # the gene model feeds locus assignment
    f.arrow((gx, ann[1]), (gx, mtop), color=f.edge("ann"))
    for b in (ann, mp):
        y = b[1] + b[3] / 2
        f.arrow((xp + wp, y), (xt, y))

    # gene table -----------------------------------------------------------------------------
    gt = f.rbox(xg, my, wg, top - my, "gt", lw=0.9)
    f.text(xg + wg / 2, top - 0.135, "one row per gene", size=F_SUB, color=t.sub, box=gt)
    gg = 0.08
    gh = (top - hdr - my - pad - 3 * gg) / 4
    for i, (a, b) in enumerate(
        [
            ("Context", "layers per\nlocus tag"),
            ("Detection", "tiers and\ngene q-value"),
            ("Abundance", "peptide\nintensities"),
            ("Detectability", "cross-fitted\np_detect"),
        ]
    ):
        y = top - hdr - (i + 1) * gh - i * gg
        ch = f.rbox(xg + pad, y, wg - 2 * pad, gh, "gt", filled=False, lw=0.6, r=0.03, z=3)
        f.labelled(ch, a, b)
    for b in (ann, mp):
        y = b[1] + b[3] / 2
        f.arrow((xt + wt, y), (xg, y))

    # analysis modules -----------------------------------------------------------------------
    mods = [
        ("Islands, regions", "autocorrelation", False),
        ("Operons", "discordant pairs,\nmissed operons", False),
        ("Function", "ORA, GSEA", True),
        ("Determinants", "promoter, RBS,\ncodon usage", False),
        ("Variants", "REF/ALT per codon", False),
        ("Novel ORFs", "and start sites", False),
    ]
    mg, split = 0.07, 0.16
    mh = (top - my - 4 * mg - split) / 6
    ys = []
    for i, (a, b, opt) in enumerate(mods):
        y = top - (i + 1) * mh - i * mg - (split - mg if i >= 4 else 0)
        ch = f.rbox(xa, y, wa, mh, "ana", lw=0.8, r=0.04, dashed=opt)
        f.labelled(ch, a, b)
        ys.append(y + mh / 2)
    bx = xg + wg + gap / 2  # gene table to the four context modules
    f.ax.plot([bx, bx], [ys[3], ys[0]], color=t.arrow, lw=0.9, zorder=4, solid_capstyle="butt")
    f.ax.plot([xg + wg, bx], [(ys[1] + ys[2]) / 2] * 2, color=t.arrow, lw=0.9, zorder=4)
    for y in ys[:4]:
        f.arrow((bx, y), (xa, y))
    lane = (bot + my) / 2  # placements to variants and novel ORFs, under the gene table
    lx = x2 + mw / 2
    f.ax.plot(
        [lx, lx, bx, bx],
        [my, lane, lane, ys[4]],
        color=f.edge("map"),
        lw=0.9,
        zorder=4,
        solid_joinstyle="miter",
    )
    for y in ys[4:]:
        f.arrow((bx, y), (xa, y), color=f.edge("map"))

    # outputs --------------------------------------------------------------------------------
    outs = [
        ("Tables", "TSV, GFF3,\nBED, VCF"),
        ("Plots", "SVG, PNG,\nPDF + data"),
        ("report.html", "self-contained"),
        ("Manifest", "versions,\nhashes, seeds"),
    ]
    oh = mh  # aligned with the first and last module rows
    og = (top - my - 4 * oh) / 3
    oys = []
    for i, (a, b) in enumerate(outs):
        y = top - (i + 1) * oh - i * og
        f.labelled(f.doc(xo, y, wo, oh), a, b)
        oys.append(y + oh / 2)
    ox = xa + wa + gap / 2  # every module's tables collect on one bus to the outputs
    for y in ys:
        f.ax.plot([xa + wa, ox], [y, y], color=t.arrow, lw=0.9, zorder=4)
    f.ax.plot(
        [ox, ox],
        [min(ys[-1], oys[-1]), max(ys[0], oys[0])],
        color=t.arrow,
        lw=0.9,
        zorder=4,
        solid_capstyle="projecting",
    )
    for y in oys:
        f.arrow((ox, y), (xo, y))

    # error control band ---------------------------------------------------------------------
    f.rbox(xi, band_y, xa + wa - xi, band_h, "io", filled=False, dashed=True, lw=0.7)
    lab = (xi, band_y, wi + gap + wp, band_h)
    f.labelled(lab, "Error control", "every claim carries\nan error estimate")
    for x, w, a, b in [
        (xt, wt, "Placement q", "target-decoy per stratum:\nexact, mismatch, de novo"),
        (xg, wg, "Gene q", "picked decoys\nor MAYU-style"),
        (xa, wa, "Module nulls", "BH, FWER, decoy\nor class FDR"),
    ]:
        f.labelled((x, band_y, w, band_h), a, b)
    # legend: dashed outline = optional module
    lg = f.rbox(
        xo + 0.02,
        band_y + 0.10,
        0.20,
        band_h - 0.20,
        "io",
        filled=False,
        lw=0.7,
        dashed=True,
        r=0.03,
    )
    f.text(
        lg[0] + lg[2] + 0.05, band_y + band_h / 2, "optional", size=F_SUB, color=t.sub, ha="left"
    )
    return f


def render(theme: Theme) -> tuple[bytes, bytes, list[str]]:
    plt.rcParams["svg.hashsalt"] = "gippa-pipeline"
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
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
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
            path = args.outdir / f"{NAME}-{theme.name}.{ext}"
            if args.check:
                if not path.is_file() or path.read_bytes() != data:
                    print(f"out of date: {path}", file=sys.stderr)
                    status = 1
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
                print(path)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
