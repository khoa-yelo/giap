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
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# the script directory is on sys.path when run
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_figure import DARK, FONT, LIGHT, Theme
from process_figure_kit import (
    ANALYSIS_NODES,
    ANNOTATE_NODES,
    BADGE,
    COLS,
    CW,
    F_HEAD,
    F_NAME,
    F_SMALL,
    MAPPING_NODES,
    MONO,
    NAME,
    NOTES,
    NOVEL_LEFT,
    REPO,
    WI,
    WO,
    XI,
    XO,
    Edge,
    Fig,
    Node,
    card_h,
    gut,
    path,
)


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
