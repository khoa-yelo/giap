"""V12 ``variants``: where do peptides support ALT residues? (GIAP_SPEC.md §10.2 V12, §7.5, §9.3
``variants.tsv`` / ``peptide_hits.tsv``; M9-18 and its lead decisions 1-10).

(a) Grouped bars of the ``variants.tsv`` sites: x = ``call`` (ALT, mixed), one bar per
``confidence``, stacked by §7.5.3 flag (a site with several flags counts once in each; sites
without flags count under "no flag"); the "all" total is printed above each bar. (b) Monospace
text alignments of the ``giap.viz.support.top_alt_sites`` (at most 30): a header line, the
reference window from ``plots/support/variant_windows.tsv`` (absent: no ``ref`` line) and one
peptide line, the covering peptide that spans most of the ±10 aa display window, with the ALT
residue bold vermillion, then "(N peptides)", the distinct covering count (lead decision 8).

A ``peptide_hits`` row covers a site when it is on the site's contig and strand, the whole codon
lies in ``[start, end)`` and the codon is in frame with the placement (lead decision 7); on a
circular replicon the codon is also tried one replicon length further (origin-crossing
placements). ``alt_type`` block / frameshift sites list a peptide unaligned with ``NOTE``.
Panel (b) residue positions are protein positions, so ``-`` strand and origin-crossing genes need
no special case. Pure: tables come from the ``RunView``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.transforms import offset_copy

from giap.variants.flags import FLAG_ORDER
from giap.viz import base, style
from giap.viz.support import VARIANT_FLANK_AA, class_error_value, top_alt_sites

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = (
    "panel", "call", "confidence", "flag", "n", "rank", "variant_id", "gene", "aa_pos",
    "ref_aa", "alt_aa", "alt_type", "n_alt", "n_ref", "class_error_rate", "ref_window",
    "peptides",
)  # fmt: skip
CALLS = ("ALT", "mixed")
CONFIDENCES = ("high", "candidate")
NO_FLAG, ALL = "no flag", "all"
#: §9.4 flag -> Okabe-Ito colour, assigned in flag order (grey is kept for "no flag")
FLAG_COLOUR: Mapping[str, str] = dict(
    zip(FLAG_ORDER, [c for k, c in style.OKABE_ITO.items() if k != "grey"], strict=False)
)
NO_FLAG_COLOUR = style.OKABE_ITO["grey"]
ALT_COLOUR = style.OKABE_ITO["vermillion"]
EN_DASH = "\N{EN DASH}"
NOTE = "length-changing: shown unaligned"
UNALIGNED = frozenset({"block", "frameshift"})
EMPTY = "No ALT site supported by peptides"
EMPTY_SAME = "No ALT calls (genome relation same)"
MONO = "DejaVu Sans Mono"
FONT_PT = 6.0
CHAR_PT = FONT_PT * 1233 / 2048  # DejaVu Sans Mono advance width
LINE_PT = 7.5
LABEL = 6  # "ref   " / "pep   "
WIDTH = 2 * VARIANT_FLANK_AA + 1


def _na(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _flags(value: object) -> list[str]:
    return [] if _na(value) else [f for f in str(value).split(";") if f]


def _flag_rank(flag: str) -> tuple[int, str]:
    return (FLAG_ORDER.index(flag), "") if flag in FLAG_ORDER else (len(FLAG_ORDER), flag)


def panel_a(variants: pd.DataFrame) -> list[dict[str, object]]:
    """Counts by call x confidence x flag (n > 0 only), sorted by call, confidence, §9.4 flag
    order, then "no flag" and "all"."""
    rows: list[dict[str, object]] = []
    for call in CALLS:
        for conf in CONFIDENCES:
            grp = variants[(variants["call"] == call) & (variants["confidence"] == conf)]
            if not len(grp):
                continue
            counts: dict[str, int] = {}
            for value in grp["flags"]:
                for f in _flags(value) or [NO_FLAG]:
                    counts[f] = counts.get(f, 0) + 1
            named = sorted((f for f in counts if f != NO_FLAG), key=_flag_rank)
            tail = [NO_FLAG] if NO_FLAG in counts else []
            for f, n in [*((f, counts[f]) for f in [*named, *tail]), (ALL, len(grp))]:
                rows.append({"panel": "a", "call": call, "confidence": conf, "flag": f, "n": n})
    return rows


def _lengths(view: RunView) -> dict[str, int]:
    """Circular replicon -> length, from ``preflight.json``."""
    pf = view.preflight() or {}
    return {
        str(r["id"]): int(r["length"])
        for r in pf.get("replicons", [])
        if r.get("topology") == "circular"
    }


def covering(
    site: pd.Series, hits: pd.DataFrame | None, lengths: Mapping[str, int]
) -> dict[str, int]:
    """Distinct covering peptide -> 0-based index of the site's residue in it."""
    if hits is None or not len(hits):
        return {}
    sel = hits[(hits["contig"] == site["contig"]) & (hits["strand"] == site["strand"])]
    c0 = int(site["codon_start"]) - 1
    codons = [c0] + ([c0 + lengths[str(site["contig"])]] if str(site["contig"]) in lengths else [])
    plus = str(site["strand"]) == "+"
    out: dict[str, int] = {}
    for pep, s0, e in zip(sel["peptide"], sel["start"], sel["end"], strict=True):
        for c in codons:
            if not (int(s0) <= c and c + 3 <= int(e)):
                continue
            off = c - int(s0) if plus else int(e) - c - 3
            if off % 3 == 0 and off // 3 < len(str(pep)):
                k = off // 3
                out[str(pep)] = min(k, out.get(str(pep), k))
    return out


def _windows(view: RunView) -> dict[str, tuple[int, str]]:
    sup = view.support("variant_windows")
    if sup is None:
        return {}
    return {
        str(v): (int(s), str(w))
        for v, s, w in zip(
            sup["variant_id"], sup["window_start_aa"], sup["ref_window"], strict=True
        )
        if not _na(s) and not _na(w)
    }


def _py(value: object) -> object:
    if _na(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    return float(value) if isinstance(value, np.floating) else value


def panel_b(
    view: RunView, variants: pd.DataFrame
) -> tuple[list[dict[str, object]], list[dict[str, int]]]:
    """Panel (b) rows in top-site order, and each site's covering peptides."""
    hits, lengths, windows = view.table("peptide_hits"), _lengths(view), _windows(view)
    rows: list[dict[str, object]] = []
    peps: list[dict[str, int]] = []
    for rank, (_, s) in enumerate(top_alt_sites(variants).iterrows(), start=1):
        cov = covering(s, hits, lengths)
        win = windows.get(str(s["variant_id"]))
        row = {"panel": "b", "rank": rank, "confidence": s["confidence"]}
        for col in ("variant_id", "gene", "aa_pos", "ref_aa", "alt_aa", "alt_type", "n_alt",
                    "n_ref", "class_error_rate", "contig", "codon_start", "strand"):  # fmt: skip
            row[col] = _py(s[col])
        row["ref_window"] = None if win is None else win[1]
        row["peptides"] = ";".join(sorted(cov)) or None
        rows.append(row)
        peps.append(cov)
    return rows, peps


def data_frame(rows: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    """The plotted data (``V12_variants.tsv``): ``COLUMNS``, object cells, NA as ``None``."""
    cells = {c: pd.Series([r.get(c) for r in rows], dtype=object) for c in COLUMNS}
    return pd.DataFrame(cells, columns=list(COLUMNS))


def title(variants: pd.DataFrame) -> str:
    """Counts and class-error range over the ``call = ALT`` sites (lead decision 4)."""
    alt = variants[variants["call"] == "ALT"]
    h, c = int((alt["confidence"] == "high").sum()), int((alt["confidence"] == "candidate").sum())
    rates = [r for r in (class_error_value(v) for v in alt["class_error_rate"]) if r is not None]
    lo, hi = (base.sig(min(rates)), base.sig(max(rates))) if rates else ("", "")
    rng = "unestimated" if not rates else lo if lo == hi else f"{lo}{EN_DASH}{hi}"
    return f"{h} high-confidence and {c} candidate ALT sites (class error {rng})"


def _draw_a(ax: Axes, rows: Sequence[Mapping[str, object]]) -> None:
    ticks: list[float] = []
    labels: list[str] = []
    used: list[str] = []
    for i, call in enumerate(CALLS):
        for j, conf in enumerate(CONFIDENCES):
            x = i * 2.6 + j
            ticks.append(x)
            labels.append(f"{conf}\n{call}" if j else conf)
            grp = [r for r in rows if r["call"] == call and r["confidence"] == conf]
            bars = [r for r in grp if r["flag"] != ALL]
            width = 0.8 / max(len(bars), 1)  # dodged: a site with k flags is in k bars
            top = 0.0
            for k, r in enumerate(bars):
                f, n = str(r["flag"]), float(str(r["n"]))
                colour = FLAG_COLOUR.get(f, NO_FLAG_COLOUR)
                ax.bar(x - 0.4 + (k + 0.5) * width, n, width, color=colour, edgecolor="white",
                       linewidth=0.4)  # fmt: skip
                top = max(top, n)
                if f not in used:
                    used.append(f)
            for r in grp:
                if r["flag"] == ALL:
                    ax.text(x, top, str(r["n"]), ha="center", va="bottom", fontsize=6)
    ax.set_xticks(ticks, labels)
    ax.set_ylim(0, max([1.0, *ax.get_ylim()[1:]]) * 1.15)
    ax.set_xlim(-0.7, 2.6 + 1.7)
    ax.set_ylabel("ALT sites (per flag)")
    ax.set_title("(a) sites by call, confidence and flag", loc="left")
    order = sorted(used, key=lambda f: (f == NO_FLAG, _flag_rank(f)))
    handles = [Patch(color=FLAG_COLOUR.get(f, NO_FLAG_COLOUR), label=f) for f in order]
    if handles:
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)


class _Lines:
    """Monospace text runs placed by character column and line number from the top-left."""

    def __init__(self, ax: Axes) -> None:
        self.ax, self.line = ax, 0.0

    def put(self, col: int, text: str, **kw: object) -> None:
        if not text:
            return
        fig = self.ax.get_figure(root=True)
        assert isinstance(fig, Figure)
        tr = offset_copy(self.ax.transAxes, fig=fig, x=col * CHAR_PT, y=-(self.line + 1) * LINE_PT,
                         units="points")  # fmt: skip
        self.ax.text(0.0, 1.0, text, transform=tr, fontfamily=MONO, fontsize=FONT_PT,
                     va="baseline", ha="left", **kw)  # type: ignore[arg-type]  # fmt: skip

    def next(self, gap: float = 1.0) -> None:
        self.line += gap


def _header(r: Mapping[str, object]) -> str:
    """The site line; a site with no gene shows its locus, ``contig:pos (strand)``, and "no gene"
    (NA cells read "NA", never Python ``None``)."""

    def cell(key: str) -> str:
        v = r.get(key)
        return "NA" if _na(v) else str(v)

    rate = class_error_value(r["class_error_rate"])
    shown = cell("class_error_rate")
    if _na(r["gene"]) or _na(r["aa_pos"]):
        where = f"{cell('contig')}:{cell('codon_start')} ({cell('strand')})  no gene (aa NA)"
    else:
        where = f"{r['gene']} aa {r['aa_pos']}"
    return (f"{where}  {cell('ref_aa')}→{cell('alt_aa')}  {cell('alt_type')}  "
            f"{cell('confidence')}  n_alt {cell('n_alt')} / n_ref {cell('n_ref')}  "
            f"class error {shown if rate is None else base.sig(rate)}")  # fmt: skip


def _best(cov: Mapping[str, int]) -> str:
    """The peptide spanning most of the display window (then the longer, then by name)."""

    def key(p: str) -> tuple[int, int, str]:
        k = cov[p]
        return (-(min(k, VARIANT_FLANK_AA) + min(len(p) - k, VARIANT_FLANK_AA + 1)), -len(p), p)

    return min(cov, key=key)


def _draw_b(ax: Axes, rows: Sequence[Mapping[str, object]], peps: Sequence[Mapping[str, int]],
            windows: Mapping[str, tuple[int, str]]) -> None:  # fmt: skip
    ax.set_axis_off()
    ax.set_title("(b) peptides at the top ALT sites", loc="left")
    lines = _Lines(ax)
    for r, cov in zip(rows, peps, strict=True):
        lines.put(0, _header(r))
        lines.next()
        genic = not pd.isna(r["aa_pos"])  # no gene: shown unaligned (lead decision 11)
        pos = int(str(r["aa_pos"])) if genic else 0
        unaligned = str(r["alt_type"]) in UNALIGNED or not genic
        win = windows.get(str(r["variant_id"]))
        if win is not None and not unaligned:
            lines.put(0, "ref   " + " " * (win[0] - (pos - VARIANT_FLANK_AA)) + win[1])
            lines.next()
        n = len(cov)
        count = f"   ({n} peptide{'' if n == 1 else 's'})"
        if not cov:
            lines.put(0, "pep   (no covering peptide)")
        elif unaligned:
            pep = sorted(cov)[0]
            shown = pep[:WIDTH] + ("…" if len(pep) > WIDTH else "")
            note = f"  {NOTE}" if genic else "  (no gene: not aligned)"
            lines.put(0, f"pep   {shown}{count}{note}")
        else:
            pep = _best(cov)
            k, la = cov[pep], max(1, len(str(r["alt_aa"])))
            left = pep[max(0, k - VARIANT_FLANK_AA) : k]
            right = pep[k + la : k + VARIANT_FLANK_AA + 1]
            col = LABEL + VARIANT_FLANK_AA - len(left)
            lines.put(0, "pep   " + " " * (col - LABEL) + left)
            lines.put(col + len(left), pep[k : k + la], color=ALT_COLOUR, fontweight="bold")
            lines.put(col + len(left) + la, right + count)
        lines.next(1.4)


@base.plot("V12", "variants", needs=("variants",), optional=("peptide_hits",))
def variants(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V12: ALT-site counts by call x confidence x flag; peptide alignments at the top sites."""
    table = view.table("variants")
    assert table is not None  # guaranteed by ``needs``
    if not len(table):
        msg = EMPTY_SAME if cfg.genome_relation == "same" else EMPTY
        return [base.empty_state("V12", "variants", msg, COLUMNS, "double")]
    rows_a = panel_a(table)
    rows_b, peps = panel_b(view, table)
    refs = [r["ref_window"] is not None and r["alt_type"] not in UNALIGNED for r in rows_b]
    h_a, h_b = 1.9, (sum(2.4 + ref for ref in refs) + 1.5) * LINE_PT / 72 + 0.15
    height = h_a + h_b + 0.35
    fig = style.new_figure("double", height)
    fig.set_label("V12_variants")
    ax_a = fig.add_axes((0.08, (h_b + 0.6) / height, 0.7, (h_a - 0.9) / height), label="a")
    ax_b = fig.add_axes((0.03, 0.1 / height, 0.94, h_b / height), label="b")
    _draw_a(ax_a, rows_a)
    _draw_b(ax_b, rows_b, peps, _windows(view))
    if not rows_b:
        ax_b.text(0.0, 0.9, "no ALT site", transform=ax_b.transAxes, color="#555555")
    text = title(table)
    fig.suptitle(text, y=1 - 0.12 / height)
    return [base.Instance(None, fig, data_frame([*rows_a, *rows_b]), text)]
