"""V5 ``mapping_sankey``: where did my sequences go? (GIAP_SPEC.md §10.2 V5; M9-05).

A custom alluvial drawn with matplotlib patches: input class (peptide, protein) → route (exact,
tolerant, miniprot, tblastn, unmapped) → outcome (placed, contaminant, rejected by FDR,
unexplained). Units are distinct sequences: peptides by sequence (``peptide_hits.peptide``,
``unmapped.sequence``), proteins by ID (``protein_hits.protein_id``, ``unmapped.sequence_id``).
A sequence both placed and unmapped counts once, as placed, and is counted in the ``conflicts``
row (M9-05 table, index decision 14). Pure: tables come from the ``RunView``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

from giap.viz import base, style

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = ("input_class", "route", "outcome", "n")
CLASSES = ("peptide", "protein")
ROUTES = ("exact", "tolerant", "miniprot", "tblastn", "unmapped")
OUTCOMES = ("placed", "contaminant", "rejected by FDR", "unexplained")
_LOST = OUTCOMES[1:]
#: the ten valid flows in node order (the data TSV grid)
GRID: tuple[tuple[str, str, str], ...] = (
    *(("peptide", r, "placed") for r in ("exact", "tolerant")),
    *(("peptide", "unmapped", o) for o in _LOST),
    *(("protein", r, "placed") for r in ("miniprot", "tblastn")),
    *(("protein", "unmapped", o) for o in _LOST),
)
CONFLICTS = ("all", "all", "conflicts")
#: placed is blue; the other outcomes take the next Okabe-Ito category colours in palette order
OUTCOME_COLOUR: Mapping[str, str] = {
    "placed": style.OKABE_ITO["blue"],
    "contaminant": style.OKABE_ITO["orange"],
    "rejected by FDR": style.OKABE_ITO["sky_blue"],
    "unexplained": style.OKABE_ITO["bluish_green"],
}
NODE_GREY = style.OKABE_ITO["grey"]
_X = (0.2, 0.47, 0.72)  # left edges of the class, route and outcome node columns
_W = 0.025  # node width (axes fraction)
_GAP = 0.04  # vertical gap between nodes (axes fraction)
_LABEL_SEP = (0.12, 0.06, 0.12)  # minimum label spacing per column (two-line labels at 0, 2)

Flow = tuple[str, str, str]


def _distinct(values: pd.Series) -> set[str]:
    return {str(v) for v in values.dropna()}


def _unmapped_outcome(label: object, reason: object) -> str:
    if label == "contaminant" or reason == "contaminant_id":
        return "contaminant"
    return "rejected by FDR" if reason == "rejected_placement_fdr" else "unexplained"


def flows(
    hits: pd.DataFrame | None, proteins: pd.DataFrame | None, unmapped: pd.DataFrame | None
) -> tuple[dict[Flow, int], int]:
    """Counts of distinct sequences per grid flow, and the number of placed-and-unmapped
    conflicts."""
    counts: dict[Flow, int] = dict.fromkeys(GRID, 0)
    placed: dict[str, set[str]] = {"peptide": set(), "protein": set()}
    if hits is not None and len(hits):
        exact = _distinct(hits.loc[hits["route"] == "exact", "peptide"])
        every = _distinct(hits["peptide"])
        counts["peptide", "exact", "placed"] = len(exact)
        counts["peptide", "tolerant", "placed"] = len(every - exact)
        placed["peptide"] = every
    if proteins is not None and "gene" in proteins.columns:
        proteins = proteins[proteins["gene"].notna()]  # placed = assigned to a gene, as in M3
    if proteins is not None and len(proteins):
        best = proteins.sort_values(["rank", "route"], kind="mergesort").drop_duplicates(
            "protein_id"
        )
        for route in ("miniprot", "tblastn"):
            counts["protein", route, "placed"] = int((best["route"] == route).sum())
        placed["protein"] = _distinct(best["protein_id"])
    n_conflicts = 0
    if unmapped is not None and len(unmapped):
        seen: set[tuple[str, str]] = set()
        for row in unmapped.itertuples(index=False):
            kind = str(row.kind)
            key = str(row.sequence if kind == "peptide" else row.sequence_id)
            if (kind, key) in seen:
                continue
            seen.add((kind, key))
            if key in placed[kind]:
                n_conflicts += 1
                continue
            counts[kind, "unmapped", _unmapped_outcome(row.label, row.reason)] += 1
    return counts, n_conflicts


def flow_table(counts: Mapping[Flow, int], n_conflicts: int) -> pd.DataFrame:
    rows = [(*key, int(counts[key])) for key in GRID] + [(*CONFLICTS, int(n_conflicts))]
    return pd.DataFrame(rows, columns=list(COLUMNS))


def title(counts: Mapping[Flow, int]) -> str:
    parts = []
    for cls, noun in (("peptide", "input peptides"), ("protein", "proteins")):
        total = sum(n for (c, _, _), n in counts.items() if c == cls)
        if total == 0:
            continue
        n_placed = sum(n for (c, _, o), n in counts.items() if c == cls and o == "placed")
        parts.append(f"{base.title_pct(n_placed, total)} of {total} {noun}")
    text = "; ".join(parts)
    head, sep, rest = text.partition(";")
    return f"{head} placed in the genome{sep}{rest}"


def _ribbon(x0: float, y0: tuple[float, float], x1: float, y1: tuple[float, float]) -> Path:
    """A band from ``x0`` (bottom, top ``y0``) to ``x1`` (``y1``) with cubic Bézier edges."""
    xm = (x0 + x1) / 2
    verts = [
        (x0, y0[1]),
        (xm, y0[1]),
        (xm, y1[1]),
        (x1, y1[1]),
        (x1, y1[0]),
        (xm, y1[0]),
        (xm, y0[0]),
        (x0, y0[0]),
        (x0, y0[1]),
    ]
    codes = [Path.MOVETO, *[Path.CURVE4] * 3, Path.LINETO, *[Path.CURVE4] * 3, Path.CLOSEPOLY]
    return Path(verts, codes)


def _stack(order: Sequence[str], size: Mapping[str, int], scale: float) -> dict[str, float]:
    """Top y of each node with flow, stacked downward from 1 in ``order``."""
    tops: dict[str, float] = {}
    y = 1.0
    for node in order:
        if size.get(node, 0) > 0:
            tops[node] = y
            y -= size[node] * scale + _GAP
    return tops


def _draw(ax: Axes, counts: Mapping[Flow, int]) -> None:
    live = [(key, n) for key, n in counts.items() if n > 0]
    columns = (CLASSES, ROUTES, OUTCOMES)
    sizes: list[dict[str, int]] = [{}, {}, {}]
    for key, n in live:
        for i in range(3):
            sizes[i][key[i]] = sizes[i].get(key[i], 0) + n
    total = sum(n for _, n in live)
    k_max = max(len(s) for s in sizes)
    scale = (1.0 - _GAP * (k_max - 1)) / total
    tops = [_stack(order, s, scale) for order, s in zip(columns, sizes, strict=True)]
    # stacking offsets per node side: left of routes by (class, outcome), right by (outcome, class)
    out_cursor = [dict(t) for t in tops]
    in_cursor = [dict(t) for t in tops]
    ordered = sorted(live, key=lambda kn: (CLASSES.index(kn[0][0]), _key(kn[0])))
    spans: dict[tuple[Flow, int], tuple[float, float]] = {}
    for key, n in ordered:  # class out-side and route in-side
        h = n * scale
        for col, cursor in ((0, out_cursor), (1, in_cursor)):
            top = cursor[col][key[col]]
            spans[key, col] = (top - h, top)
            cursor[col][key[col]] = top - h
    for key, n in sorted(live, key=lambda kn: (OUTCOMES.index(kn[0][2]), _key(kn[0]))):
        h = n * scale
        for col, cursor, tag in ((1, out_cursor, 3), (2, in_cursor, 2)):
            top = cursor[col][key[col]]
            spans[key, tag] = (top - h, top)
            cursor[col][key[col]] = top - h
    for key, _ in live:
        colour = OUTCOME_COLOUR[key[2]]
        for (xa, sa), (xb, sb) in (((_X[0], 0), (_X[1], 1)), ((_X[1], 3), (_X[2], 2))):
            path = _ribbon(xa + _W, spans[key, sa], xb, spans[key, sb])
            ax.add_patch(PathPatch(path, facecolor=colour, edgecolor="none", alpha=0.45))
    for col, (order, size, col_tops) in enumerate(zip(columns, sizes, tops, strict=True)):
        nodes = [node for node in order if node in col_tops]
        centres = []
        for node in nodes:
            h = size[node] * scale
            colour = OUTCOME_COLOUR[node] if col == 2 else NODE_GREY
            ax.add_patch(Rectangle((_X[col], col_tops[node] - h), _W, h, facecolor=colour, lw=0))
            centres.append(col_tops[node] - h / 2)
        ha, x = ("right", _X[col] - 0.01) if col == 0 else ("left", _X[col] + _W + 0.01)
        sep = "\n" if col != 1 else " "
        for node, y in zip(nodes, spread(centres, _LABEL_SEP[col]), strict=True):
            ax.text(x, y, f"{node}{sep}({size[node]:,})", ha=ha, va="center", fontsize=6)


def spread(ys: Sequence[float], sep: float, floor: float = 0.0) -> list[float]:
    """Label positions for descending ``ys``: each at least ``sep`` below the previous one,
    shifted up together when the last would fall below ``floor``."""
    out: list[float] = []
    for y in ys:
        out.append(y if not out else min(y, out[-1] - sep))
    shift = max(0.0, floor - out[-1]) if out else 0.0
    return [y + shift for y in out]


def _key(flow: Flow) -> tuple[int, int, int]:
    return CLASSES.index(flow[0]), ROUTES.index(flow[1]), OUTCOMES.index(flow[2])


@base.plot("V05", "mapping_sankey", needs=("peptide_hits", "unmapped"), optional=("protein_hits",))
def mapping_sankey(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V5: one alluvial of distinct sequences, input class → route → outcome."""
    counts, n_conflicts = flows(
        view.table("peptide_hits"), view.table("protein_hits"), view.table("unmapped")
    )
    if sum(counts.values()) == 0:
        msg = "No sequences in the mapping tables"
        return [base.empty_state("V05", "mapping_sankey", msg, COLUMNS, "single")]
    both = all(sum(n for k, n in counts.items() if k[0] == c) for c in CLASSES)
    fig = style.new_figure("double" if both else "single", 2.8)
    fig.set_label("V05_mapping_sankey")
    ax = fig.add_axes((0.02, 0.03, 0.96, 0.82))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.02)
    ax.set_axis_off()
    _draw(ax, counts)
    text = title(counts)
    fig.suptitle(text, y=0.97)
    return [base.Instance(None, fig, flow_table(counts, n_conflicts), text)]
