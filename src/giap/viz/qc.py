"""V17 ``qc``: are the inputs sane? (M9-04; GIAP_SPEC.md §10.2 V17, §7.0, §9.3).

Five panels and a warning list on a double-width 2 x 3 grid: (a) preflight check counts by level,
(b) contig lengths, (c) distinct peptide lengths (placed / unmapped), (d) the measured and imputed
``abundance_log2`` histogram, (e) the NA share of the §9.3 feature columns over ``partial = false``
genes. The data frame is the long ``panel, group, key, value, n`` table of
``plots/data/V17_qc.tsv``. Pure: every table comes from the ``RunView``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from giap.viz.base import Instance, Skip, plot, title_pct
from giap.viz.style import DETECTION, OKABE_ITO, REGION, new_figure

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = ("panel", "group", "key", "value", "n")
LEVELS = ("ok", "info", "warn", "error")
LEVEL_COLOUR: Mapping[str, str] = {
    "ok": OKABE_ITO["blue"],
    "info": OKABE_ITO["sky_blue"],
    "warn": OKABE_ITO["orange"],
    "error": OKABE_ITO["vermillion"],
}
KIND_COLOUR: Mapping[str, str] = {
    "chromosome": OKABE_ITO["blue"],
    "plasmid": REGION["plasmid"],
    "unknown": OKABE_ITO["grey"],
}
PLACED, UNMAPPED, IMPUTED = DETECTION["detected"], DETECTION["not_detected"], OKABE_ITO["orange"]
#: §9.3 feature groups of ``genes.tsv`` plus ``abundance_raw`` (panel e)
FEATURE_GROUPS: Mapping[str, tuple[str, ...]] = {
    "abundance": ("abundance_raw",),
    "codon usage": ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10"),
    "promoters and transcription": (
        "promoter_score",
        "promoter_score_tu",
        "tx_rate_max_tu",
        "tx_rate_sum_tu",
        "terminator_conf_downstream",
    ),
    "protein": ("mw_da", "pi", "gravy", "tm_segments", "n_observable_peptides"),
    "translation initiation": ("tir_ostir", "sd_dg", "mrna_dg", "spacer_bp", "start_codon"),
}
MAX_WARN_LINES = 10
N_BINS = 40
_MUTED = "#555555"
#: characters kept per warning; matplotlib wraps the line at the figure edge when drawing
_LINE_CHARS = 110
#: panel (e) bar colour per feature group
GROUP_COLOUR: Mapping[str, str] = {
    "abundance": OKABE_ITO["orange"],
    "codon usage": OKABE_ITO["bluish_green"],
    "promoters and transcription": OKABE_ITO["blue"],
    "protein": OKABE_ITO["reddish_purple"],
    "translation initiation": OKABE_ITO["sky_blue"],
}

Row = tuple[str, str, Any, Any, Any]


def _note(ax: Axes, text: str) -> None:
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=_MUTED)
    ax.set_xticks([])
    ax.set_yticks([])


def _warning_line(check: Mapping[str, Any]) -> str:
    rule_id = str(check.get("rule_id", ""))
    message = str(check.get("message") or "").strip() or str(check.get("rule") or "").strip()
    line = f"{rule_id}: {message}" if message else rule_id
    line = " ".join(line.split())  # one logical line: no embedded newlines
    return line if len(line) <= _LINE_CHARS else line[: _LINE_CHARS - 1] + "\u2026"


def _panel_a(ax: Axes, cell: Axes, checks: Sequence[Mapping[str, Any]]) -> list[Row]:
    counts = {lvl: sum(str(c.get("level")) == lvl for c in checks) for lvl in LEVELS}
    ax.bar(range(4), [counts[lvl] for lvl in LEVELS], color=[LEVEL_COLOUR[x] for x in LEVELS])
    ax.set_xticks(range(4), list(LEVELS))
    ax.set_ylabel("checks")
    ax.set_title("(a) preflight checks", loc="left")
    flagged = [c for c in checks if str(c.get("level")) in ("warn", "error")]
    lines = [_warning_line(c) for c in flagged[:MAX_WARN_LINES]]
    if len(flagged) > MAX_WARN_LINES:
        lines.append(f"+ {len(flagged) - MAX_WARN_LINES} more")
    cell.set_axis_off()
    cell.set_title("warnings and errors", loc="left")
    text = "\n".join(lines) if lines else "no warnings"
    note = cell.text(
        0.0, 1.0, text, ha="left", va="top", transform=cell.transAxes, fontsize=6, wrap=True
    )
    note.set_in_layout(False)
    rows: list[Row] = [("a", "level", lvl, counts[lvl], None) for lvl in sorted(LEVELS)]
    msgs = sorted(str(c.get("rule_id", "")) for c in flagged)
    return rows + [("a", "message", rule_id, None, None) for rule_id in msgs]


def _panel_b(ax: Axes, replicons: Sequence[Mapping[str, Any]], threshold: int) -> list[Row]:
    ax.set_title("(b) contig lengths", loc="left")
    ax.set_xscale("log")
    ax.set_xlabel("length (bp)")
    ys = list(range(len(replicons)))
    lengths = [int(r["length"]) for r in replicons]
    kinds = [str(r.get("kind") or "unknown") for r in replicons]
    colours = [KIND_COLOUR.get(k, KIND_COLOUR["unknown"]) for k in kinds]
    ax.barh(ys, lengths, left=0, color=colours, height=0.7)
    circ = [
        (n, y)
        for n, y, r in zip(lengths, ys, replicons, strict=True)
        if r.get("topology") == "circular"
    ]
    if circ:
        ax.plot(
            [n for n, _ in circ],
            [y for _, y in circ],
            ls="none",
            marker="o",
            ms=3,
            color="black",
            label="circular",
        )
    ax.axvline(threshold, ls="--", lw=0.8, color=_MUTED)
    label = ax.text(
        threshold,
        0.5,
        "genome.contiguity ",
        transform=ax.get_xaxis_transform(),
        rotation=90,
        ha="right",
        va="center",
        fontsize=6,
        color=_MUTED,
    )
    label.set_in_layout(False)
    ax.set_yticks(ys, [str(r["id"]) for r in replicons])
    ax.invert_yaxis()  # FASTA order top to bottom
    span = [*lengths, threshold]
    ax.set_xlim(10 ** math.floor(math.log10(max(min(span), 1) / 2)), max(span) * 3)
    data = [
        ("b", k, str(r["id"]), n, None) for k, r, n in zip(kinds, replicons, lengths, strict=True)
    ]
    return sorted(data, key=lambda row: (row[1], row[2]))


def _panel_c(
    ax: Axes, hits: pd.DataFrame | None, unmapped: pd.DataFrame | None, bounds: Sequence[int]
) -> list[Row]:
    ax.set_title("(c) peptide lengths", loc="left")
    ax.set_xlabel("length (aa)")
    placed = set() if hits is None else {str(p) for p in hits["peptide"].dropna()}
    lost: set[str] = set()
    if unmapped is not None:
        sel = unmapped[unmapped["kind"].astype(str) == "peptide"]["sequence"].dropna()
        lost = {str(s) for s in sel}
    if not placed and not lost:
        _note(ax, "no peptide input")
        return []
    layers = {
        "placed": pd.Series([len(p) for p in placed], dtype=int).value_counts(),
        "unmapped": pd.Series([len(p) for p in lost], dtype=int).value_counts(),
    }
    bottom: dict[int, int] = {}
    for group, colour in (("placed", PLACED), ("unmapped", UNMAPPED)):
        counts = layers[group].sort_index()
        xs = [int(x) for x in counts.index]
        under = [bottom.get(x, 0) for x in xs]
        ax.bar(xs, counts.to_numpy(), width=1.0, bottom=under, color=colour, label=group)
        for x, v in zip(xs, counts.to_numpy(), strict=True):
            bottom[x] = bottom.get(x, 0) + int(v)
    for bound in bounds:
        ax.axvline(bound, ls=":", lw=0.8, color=_MUTED)
    ax.set_ylabel("distinct peptides")
    ax.legend(frameon=False, loc="best")
    return [
        ("c", g, int(k), int(v), None)
        for g in ("placed", "unmapped")
        for k, v in layers[g].sort_index().items()
    ]


def _whole(genes: pd.DataFrame) -> pd.DataFrame:
    return genes[~genes["partial"].fillna(False).astype(bool)]


def _panel_d(ax: Axes, whole: pd.DataFrame) -> list[Row]:
    ax.set_title("(d) abundance", loc="left")
    ax.set_xlabel("log2 abundance")
    log2 = pd.to_numeric(whole["abundance_log2"], errors="coerce")
    imputed = whole["abundance_imputed"].fillna(False).astype(bool)
    layers = {
        "measured": log2[~imputed].dropna().to_numpy(float),
        "imputed": log2[imputed].dropna().to_numpy(float),
    }
    pooled = np.concatenate(list(layers.values()))
    if pooled.size == 0:
        _note(ax, "no abundance input")
        return []
    lo, hi = float(pooled.min()), float(pooled.max())
    if lo == hi:
        lo, hi = lo - 0.5, hi + 0.5
    edges = np.linspace(lo, hi, N_BINS + 1)
    width = float(edges[1] - edges[0])
    rows: list[Row] = []
    bottom = np.zeros(N_BINS, dtype=int)
    for group, colour in (("measured", PLACED), ("imputed", IMPUTED)):
        counts, _ = np.histogram(layers[group], bins=edges)
        keep = counts > 0
        ax.bar(
            edges[:-1][keep],
            counts[keep],
            width=width,
            bottom=bottom[keep],
            align="edge",
            color=colour,
            label=group,
        )
        bottom += counts
        n = int(layers[group].size)
        rows += [
            ("d", group, float(e), int(c), n)
            for e, c in zip(edges[:-1][keep], counts[keep], strict=True)
        ]
    ax.set_ylabel("genes")
    ax.legend(frameon=False, loc="best")
    return sorted(rows, key=lambda row: (row[1], row[2]))


def _panel_e(ax: Axes, whole: pd.DataFrame) -> list[Row]:
    ax.set_title("(e) missingness (non-partial genes)", loc="left")
    ax.set_xlabel("missing (%)")
    n = len(whole)
    rows: list[Row] = []
    for group in sorted(FEATURE_GROUPS):
        for col in sorted(FEATURE_GROUPS[group]):
            share = float(whole[col].isna().sum()) / n if n else None
            rows.append(("e", group, col, share, n))
    ys = list(range(len(rows)))
    widths = [0.0 if r[3] is None else 100 * r[3] for r in rows]
    for group in sorted(FEATURE_GROUPS):
        sel = [i for i, r in enumerate(rows) if r[1] == group]
        ax.barh(sel, [widths[i] for i in sel], color=GROUP_COLOUR[group], height=0.7, label=group)
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    ax.set_yticks(ys, [str(r[2]) for r in rows], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    return rows


def _title(n_warn: int, n_err: int, whole: pd.DataFrame | None) -> str:
    head = f"{n_warn} preflight warnings"
    if whole is None:
        return f"{head}, {n_err} errors; abundance not available"
    log2 = pd.to_numeric(whole["abundance_log2"], errors="coerce")
    if log2.notna().sum() == 0:
        return f"{head}; no abundance input"
    measured = int((log2.notna() & ~whole["abundance_imputed"].fillna(False).astype(bool)).sum())
    share = title_pct(measured, len(whole))
    return f"{head}, {n_err} errors; {share} of genes have a measured abundance"


@plot("V17", "qc", needs=["preflight"], optional=["genes_tsv", "peptide_hits", "unmapped"])
def qc(view: RunView, cfg: GiapConfig) -> list[Instance] | Skip:
    """The V17 figure on one run (see the module docstring)."""
    pf = view.preflight()
    if pf is None:  # the registry wrapper skips first; kept for direct calls
        return Skip("input_absent:preflight")
    checks = list(pf.get("checks") or [])
    fig = new_figure("double", 6.0, layout="constrained")
    fig.set_label("V17_qc")
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 1.4))
    axes = {
        lab: fig.add_subplot(grid[i // 3, i % 3], label=lab)
        for i, lab in enumerate(("a", "b", "c", "d", "e", "warnings"))
    }
    rows = _panel_a(axes["a"], axes["warnings"], checks)
    rows += _panel_b(
        axes["b"], list(pf.get("replicons") or []), int(cfg.preflight.min_largest_contig_warn)
    )
    bounds = [int(b) for b in cfg.input.peptide_len]
    rows += _panel_c(axes["c"], view.table("peptide_hits"), view.table("unmapped"), bounds)
    genes = view.table("genes_tsv")
    whole = None if genes is None else _whole(genes)
    if whole is None:
        for lab in ("d", "e"):
            axes[lab].set_title(
                f"({lab}) {'abundance' if lab == 'd' else 'missingness'}", loc="left"
            )
            _note(axes[lab], "not available (no genes.tsv)")
    else:
        rows += _panel_d(axes["d"], whole)
        rows += _panel_e(axes["e"], whole)
    n_warn = sum(str(c.get("level")) == "warn" for c in checks)
    n_err = sum(str(c.get("level")) == "error" for c in checks)
    title = _title(n_warn, n_err, whole)
    fig.suptitle(title)
    data = pd.DataFrame(rows, columns=list(COLUMNS), dtype=object)
    return [Instance(suffix=None, figure=fig, data=data, title=title)]
