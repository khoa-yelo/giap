"""V13 ``novel_orfs``: which new genes or upstream extensions are supported? (GIAP_SPEC.md
§10.2 V13, §7.6, §9.3; M9-12 and its lead decisions 1-10).

Panel (a) counts the candidates per ``class`` (schema enum order), stacked by status: ``called``
is authoritative (decision 5); an uncalled row whose ``class_fdr`` is NA, the literal
``unestimated`` or flagged ``class_fdr_unestimated`` is ``unestimated``; the rest are
``not_called``. Panel (b) draws the ``TOP_N`` called rows with the smallest ``rank`` with the V3
renderer (``giap.viz.locus.render_locus``, ``Locus(kind="novel_orf")``), one instance each, plus
an index instance. ``start_longest`` and ``end`` are in transcription direction; ``RunView``
already subtracted 1 from ``start_longest`` / ``start_best_ostir``, which on ``-`` are the high
ends, so 1 is added back there (decision 6). A row crosses the origin when ``end <
start_longest`` on ``+`` or ``start_longest < end`` on ``-``; its span is unwrapped past the
replicon length (decision 7).
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from giap.viz import base, locus, style
from giap.viz.locus import Locus, OrfMark, render_locus

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

CLASSES: tuple[str, ...] = (
    "nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
    "pseudogene", "sorf",
)  # fmt: skip
STATUSES: tuple[str, ...] = ("called", "not_called", "unestimated")
STATUS_COLOUR = {
    "called": style.DETECTION["detected"],
    "not_called": style.DETECTION["not_detected"],
    "unestimated": style.DETECTION["unestimated"],
}
COUNT_COLUMNS = ("class", "status", "n", "class_fdr_max_called")
INDEX_COLUMNS = ("rank", "orf_id", "rendered")
TOP_N = 20
EMPTY = "No novel ORF candidates in this run"
_OPTIONAL = ("genes_tsv", "peptide_hits", "promoters_bed", "terminators_bed", "regions_tsv",
             "islands")  # fmt: skip


def _fdr(v: object) -> float | None:
    """``class_fdr`` as a finite float; NA and the literal ``unestimated`` are ``None``."""
    try:
        return locus._num(v)
    except ValueError:
        return None


def _status(row: Any) -> str:
    if row.called is not pd.NA and bool(row.called):
        return "called"
    flags = str(row.flags).split(";") if isinstance(row.flags, str) else []
    if _fdr(row.class_fdr) is None or "class_fdr_unestimated" in flags:
        return "unestimated"
    return "not_called"


def _fmt(v: float | None) -> str:
    return "unestimated" if v is None else base.sig(v)


def _counts(orfs: pd.DataFrame) -> pd.DataFrame:
    status = [_status(r) for r in orfs.itertuples(index=False)]
    fdr = [_fdr(v) for v in orfs["class_fdr"]]
    rows = []
    for c in CLASSES:
        for s in STATUSES:
            hit = [i for i, k in enumerate(orfs["cls"]) if k == c and status[i] == s]
            if not hit:
                continue
            vals = [f for i in hit if (f := fdr[i]) is not None] if s == "called" else []
            rows.append([c, s, len(hit), max(vals) if vals else math.nan])
    data = pd.DataFrame(rows, columns=list(COUNT_COLUMNS))
    return data.astype({"n": "int64", "class_fdr_max_called": float})


def _bar_label(d: pd.DataFrame) -> str:
    called = d[d["status"] == "called"]
    if called.empty:
        return "unestimated" if set(d["status"]) == {"unestimated"} else "none called"
    fdr = called["class_fdr_max_called"].dropna()
    return f"FDR ≤ {fdr.max():g}" if len(fdr) else "unestimated"


def _title(data: pd.DataFrame) -> str:
    n = {s: int(data.loc[data["status"] == s, "n"].sum()) for s in STATUSES}
    if not n["called"]:
        return f"No novel ORF called; {n['unestimated']} candidates unestimated"
    fdr = data["class_fdr_max_called"].dropna()
    bound = f"class FDR ≤ {fdr.max():g}" if len(fdr) else "class FDR unestimated"
    orfs = "ORF" if n["called"] == 1 else "ORFs"
    return f"{n['called']} novel {orfs} called ({bound}), {n['unestimated']} unestimated"


def _panel_a(orfs: pd.DataFrame) -> base.Instance:
    data = _counts(orfs)
    title = _title(data)
    fig = style.new_figure("single", 2.4, layout="constrained")
    ax = fig.add_subplot(label="counts")
    left = np.zeros(len(CLASSES))
    for s in STATUSES:
        sub = data[data["status"] == s]
        n = [int(sub.loc[sub["class"] == c, "n"].sum()) for c in CLASSES]
        ax.barh(range(len(CLASSES)), n, left=left, color=STATUS_COLOUR[s], height=0.7,
                label=s.replace("_", " "))  # fmt: skip
        left += n
    for i, c in enumerate(CLASSES):
        if left[i]:
            label = _bar_label(data[data["class"] == c])
            ax.text(left[i] + 0.02 * left.max(), i, label, va="center", fontsize=6)
    ax.set_yticks(range(len(CLASSES)), list(CLASSES))
    ax.invert_yaxis()
    ax.set_xlim(0, math.ceil(1.6 * left.max()))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("candidates")
    # below the axes: inside, the key met the bar labels (sim100k: "unestimated" twice)
    fig.legend(loc="outside lower center", ncols=len(STATUSES), frameon=False)
    fig.suptitle(title, fontsize=7)
    return base.Instance(None, fig, data, title)


def _int(v: object) -> int | None:
    f = locus._num(v)
    return None if f is None else int(f)


def _span(row: Any, view: RunView, length: int | None) -> tuple[int, int]:
    """The 0-based half-open span (decisions 6, 7); ``start_longest`` NA -> first peptide."""
    minus, end, first = str(row.strand) == "-", int(row.end), _int(row.start_longest)
    if first is None:
        peps = (
            {locus.fold_il(p) for p in row.peptides.split(";")}
            if isinstance(row.peptides, str)
            else set()
        )
        hits = view.table("peptide_hits")
        if hits is not None and len(hits):
            keys = hits["peptide"].astype(str).map(locus.fold_il)
            hits = hits[(hits["contig"] == row.contig) & keys.isin(peps)]
        if hits is None or hits.empty:
            return (end - 1, end + 2) if minus else (end - 3, end)
        return (end - 1, int(hits["end"].max())) if minus else (int(hits["start"].min()), end)
    if minus:
        lo, hi = end - 1, first + 1
    else:
        lo, hi = first, end
    if hi <= lo:  # crosses the origin
        hi += length or 0
    return (lo, hi) if hi > lo else (hi - 1, lo + 1)


def _tick(row: Any, span: tuple[int, int], length: int | None) -> int | None:
    best, first = _int(row.start_best_ostir), _int(row.start_longest)
    if best is None or best == first:
        return None
    x = best + 1 if str(row.strand) == "-" else best
    for k in (0, 1, -1):
        if span[0] <= x + k * (length or 0) <= span[1]:
            return x + k * (length or 0)
    return None


def _view(row: Any, view: RunView, cfg: GiapConfig, scale: tuple[float, float]) -> base.Instance:
    length, _ = locus._replicon(view, str(row.contig))
    span = _span(row, view, length)
    q, fdr = locus._num(row.best_peptide_q), _fdr(row.class_fdr)
    n = _int(row.n_peptides)
    label = (f"{row.orf_id} ({row.cls}): {'NA' if n is None else n} "
             f"peptides, best q {base.sig(q)}, class FDR {_fmt(fdr)}")  # fmt: skip
    lo = Locus(str(row.orf_id), "novel_orf", str(row.contig), *span, fdr, label)
    peps = frozenset(str(row.peptides).split(";")) if isinstance(row.peptides, str) else frozenset()
    mark = OrfMark(str(row.strand), _tick(row, span, length), peps)
    inst = render_locus(lo, view, cfg, scale, orf=mark)
    return dataclasses.replace(inst, suffix=f"locus_{row.orf_id}")


def _index(ranks: Sequence[int], ids: Sequence[str]) -> base.Instance:
    data = pd.DataFrame({"rank": list(ranks), "orf_id": list(ids), "rendered": True},
                        columns=list(INDEX_COLUMNS))  # fmt: skip
    fig = style.new_figure("single", 0.4 + 0.13 * len(ids))
    ax = fig.add_axes((0.02, 0.02, 0.96, 0.96))
    ax.set_axis_off()
    for i, (r, oid) in enumerate(zip(ranks, ids, strict=True)):
        ax.text(0.0, 1.0 - (i + 0.5) / len(ids), f"{r}. {oid}", fontsize=6, va="center",
                transform=ax.transAxes)  # fmt: skip
    return base.Instance("locus", fig, data, f"V13 novel ORF loci: top {len(ids)} by rank")


@base.plot("V13", "novel_orfs", needs=["novel_orfs"], optional=list(_OPTIONAL))
def novel_orfs(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """Panel (a), then one locus view per top-ranked called ORF, then the index."""
    orfs = view.table("novel_orfs")
    if orfs is None or orfs.empty:
        return [base.empty_state("V13", "novel_orfs", EMPTY, COUNT_COLUMNS, "single")]
    orfs = orfs.rename(columns={"class": "cls"})  # ``class`` is no attribute name
    out = [_panel_a(orfs)]
    called = orfs[orfs["called"].fillna(False).astype(bool)]
    called = called.assign(_r=called["rank"].astype("Float64").fillna(math.inf))
    top = called.sort_values(["_r", "orf_id"], kind="stable").head(TOP_N)
    if top.empty:
        return out
    scale = locus.run_scale(view.table("genes_tsv"))
    out += [_view(r, view, cfg, scale) for r in top.itertuples(index=False)]
    ranks = [-1 if (k := _int(r)) is None else k for r in top["rank"]]
    return [*out, _index(ranks, [str(i) for i in top["orf_id"]])]
