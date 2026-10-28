"""Row builders for the M9-09 plots (V2 ``abundance_landscape``, V10 ``correlogram``).

Genes are placed by their 0-based midpoint ``mid`` (``(start0 + end) / 2``): the published
1-based closed row is ``start = mid - half + 1``, ``end = mid + half``. Tables are written with
``helpers.viz_out_tree.write_schema_table`` under the ``genes``, ``islands`` and
``autocorrelation`` schemas; ``preflight/preflight.json`` carries the replicons in FASTA order.
Nothing here touches ``$GIPPA_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import numpy as np

from helpers.viz_out_tree import gene_row, make_out, write_schema_table

#: (id, length, topology, kind)
Replicon = tuple[str, int, str, str]
CHR1: Replicon = ("chr1", 100_000, "circular", "chromosome")


def write_preflight(out: Path, replicons: Sequence[Replicon]) -> Path:
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema_version": 1,
        "replicons": [
            {"id": r[0], "length": r[1], "topology": r[2], "kind": r[3]} for r in replicons
        ],
        "contig_map": [{"id": r[0], "original_id": r[0]} for r in replicons],
    }
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return path


def gene(
    tag: str,
    mid: int,
    log2: float | None = None,
    *,
    contig: str = "chr1",
    status: str = "detected",
    imputed: bool = False,
    partial: bool = False,
    half: int = 500,
    crosses_origin: bool = False,
    abundance_type: str = "ibaq",
) -> dict[str, object]:
    """A ``genes.tsv`` row centred on the 0-based midpoint ``mid`` (``end`` may exceed the
    length for an origin-crossing gene; pass ``mid`` unwrapped then)."""
    return gene_row(
        locus_tag=tag,
        contig=contig,
        start=mid - half + 1,
        end=mid + half,
        feature_type="cds",
        detection_status=status,
        detected=status == "detected",
        partial=partial,
        crosses_origin=crosses_origin,
        abundance_log2=log2,
        abundance_raw=None if log2 is None else 2.0**log2,
        abundance_imputed=imputed,
        abundance_type=abundance_type,
        replicon_type="chromosome",
    )


def island(
    island_id: str,
    direction: str,
    start: int,
    end: int,
    p_value: float | None,
    contig: str = "chr1",
) -> dict[str, object]:
    """An ``islands.tsv`` row; ``start``/``end`` are the published 1-based closed coordinates."""
    return {
        "island_id": island_id,
        "direction": direction,
        "contig": contig,
        "start": start,
        "end": end,
        "tu_first": "TU0001",
        "tu_last": "TU0004",
        "n_tus": 4,
        "n_genes": 6,
        "genes": "g1;g2",
        "mean_abundance_log2": 2.5,
        "detected_frac": 0.9,
        "scan_s": 4.0,
        "block_len": 3,
        "stat": "scan_max",
        "effect": 1.25 if direction == "high" else -1.25,
        "p_value": p_value,
        "method": "scan_max_stat",
        "null_model": "circular_block_bootstrap",
        "n_perm": 999,
        "flags": "ci_unestimated",
    }


def ac_row(
    contig: str,
    statistic: str,
    null_model: str,
    lag: int | None,
    value: float | None,
    null_mean: float | None,
    null_lo: float | None,
    null_hi: float | None,
    q_value: float | None,
    flags: str | None = "ci_unestimated",
) -> dict[str, object]:
    """An ``autocorrelation.tsv`` row as M5-05/M5-09 write it (ADR-0020)."""
    finite = value is not None and null_mean is not None
    return {
        "contig": contig,
        "statistic": statistic,
        "lag": lag,
        "null_model": null_model,
        "value": value,
        "null_mean": null_mean,
        "null_lo": null_lo,
        "null_hi": null_hi,
        "stat": statistic,
        "effect": (value - null_mean) if finite else None,  # type: ignore[operator]
        "p_value": q_value,
        "q_value": q_value,
        "method": "moran_i_permutation" if statistic == "morans_i" else "join_count_permutation",
        "n_perm": 999 if finite else 0,
        "flags": flags,
    }


#: (value, null_mean, null_lo, null_hi, q) for one (statistic, null, lag)
Cell = tuple[float | None, float | None, float | None, float | None, float | None]


def default_cell(statistic: str, null_model: str, lag: int) -> Cell:
    """A smooth, writer-exact (≤ 4 significant digits) correlogram cell; nothing significant."""
    if statistic == "morans_i":
        value = round(0.2 * 0.8**lag, 4)
        mean = -0.002 if null_model == "gene_label" else -0.004
        half = 0.05 if null_model == "gene_label" else 0.09
        return value, mean, round(mean - half, 4), round(mean + half, 4), 0.5
    value = float(40 - lag)
    mean = 38.0 if null_model == "gene_label" else 37.5
    half = 3.0 if null_model == "gene_label" else 5.0
    return value, mean, mean - half, mean + half, 0.6


def correlogram_rows(
    contig: str = "chr1",
    lags: Sequence[int] = tuple(range(1, 21)),
    cell: Callable[[str, str, int], Cell] = default_cell,
    overrides: Mapping[tuple[str, str, int], Cell] | None = None,
    oriter: bool = True,
) -> list[dict[str, object]]:
    """Every morans_i and join_count_bb row (both nulls, every lag) plus the oriter_slope
    rows M5-09 writes; ``overrides[(statistic, null_model, lag)]`` replaces one cell. The
    observed ``value`` of a (statistic, lag) is shared by both nulls unless overridden."""
    overrides = overrides or {}
    rows = []
    for statistic in ("morans_i", "join_count_bb"):
        for null_model in ("gene_label", "tu_block"):
            for lag in lags:
                c = overrides.get((statistic, null_model, lag), cell(statistic, null_model, lag))
                rows.append(ac_row(contig, statistic, null_model, lag, *c))
    if oriter:
        rows.append(
            {
                **ac_row(contig, "oriter_slope", "bootstrap", None, 0.12, None, None, None, None),
                "flags": "descriptive",
            }
        )
    return rows


def na_correlogram_rows(contig: str = "chr1") -> list[dict[str, object]]:
    """The validation run's rows, flagged ``few_values`` / ``not_estimable``: Moran's I NA, the
    join count written as 0 with zero nulls and q NA, as the producer writes it (#218 review B1)."""
    rows = []
    for statistic, flags, v in (
        ("morans_i", "ci_unestimated;few_values;not_estimable", None),
        ("join_count_bb", "ci_unestimated;not_estimable", 0.0),
    ):
        for null_model in ("gene_label", "tu_block"):
            for lag in range(1, 21):
                rows.append(ac_row(contig, statistic, null_model, lag, v, v, v, v, None, flags))
    rows.append(
        {
            **ac_row(contig, "oriter_slope", "bootstrap", None, None, None, None, None, None),
            "flags": "descriptive;not_estimable",
        }
    )
    return rows


def landscape_out(
    tmp: Path,
    *,
    replicons: Sequence[Replicon] = (CHR1,),
    genes: Sequence[Mapping[str, object]] | None = None,
    islands: Sequence[Mapping[str, object]] | None = None,
    autocorr: Sequence[Mapping[str, object]] | None = None,
    preflight: bool = True,
) -> Path:
    """An ``OUT/`` with the given tables; ``None`` leaves a table out, ``[]`` writes it
    header-only."""
    out = make_out(tmp)
    if preflight:
        write_preflight(out, replicons)
    for name, key, rows in (
        ("genes", "genes", genes),
        ("islands", "islands", islands),
        ("autocorrelation", "autocorrelation", autocorr),
    ):
        if rows is not None:
            write_schema_table(out / "tables" / f"{name}.tsv", key, rows)
    return out


def lowess_oracle(
    y: Sequence[float], pos: Sequence[float], frac: float, circular: bool, period: float
) -> np.ndarray:
    """statsmodels ``lowess(it=0, delta=0)`` with M5-01's span ``k = min(n, max(3, ⌈frac·n⌉))``;
    a circular series is padded with its last and first ``k`` points shifted by ``∓period`` and
    the middle ``n`` fits are returned. Independent of ``gippa.stats.smooth``."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    x = np.asarray(pos, dtype=float)
    v = np.asarray(y, dtype=float)
    n = x.size
    k = int(min(n, max(3, math.ceil(frac * n - 1e-9))))
    if circular:
        x = np.concatenate([x[n - k :] - period, x, x[:k] + period])
        v = np.concatenate([v[n - k :], v, v[:k]])
    fit = lowess(v, x, frac=k / x.size, it=0, delta=0.0, return_sorted=False)
    return np.asarray(fit[k : k + n] if circular else fit, dtype=float)


# --- figure inspection -------------------------------------------------------------------------


def plot_spec(pid: str):  # type: ignore[no-untyped-def]
    """The registered ``PlotSpec`` of ``pid`` after importing every ``gippa.viz`` module (the
    ``plots`` step discovers plots the same way), so the tests do not fix a module name."""
    import importlib
    import pkgutil

    import gippa.viz
    from gippa.viz import base

    for info in pkgutil.walk_packages(gippa.viz.__path__, prefix="gippa.viz."):
        importlib.import_module(info.name)
    return base.registered()[pid]


def markers(ax) -> list[tuple[float, float, bool, str]]:  # type: ignore[no-untyped-def]
    """(x, y, hollow, face-or-edge colour hex) of every marker in ``ax``: Line2D markers and
    scatter offsets."""
    import matplotlib.colors as mcolors
    from matplotlib.collections import PathCollection

    out = []
    for line in ax.get_lines():
        if line.get_marker() in (None, "", "None", " ", "none", "|", "_"):
            continue
        face = line.get_markerfacecolor()
        hollow = (isinstance(face, str) and face.lower() == "none") or mcolors.to_rgba(face)[3] == 0
        colour = line.get_markeredgecolor() if hollow else face
        for x, y in zip(line.get_xdata(), line.get_ydata(), strict=True):
            out.append((float(x), float(y), hollow, mcolors.to_hex(colour).upper()))
    for coll in ax.collections:
        if not isinstance(coll, PathCollection):
            continue
        faces = coll.get_facecolors()
        edges = coll.get_edgecolors()
        for i, (x, y) in enumerate(coll.get_offsets()):
            hollow = len(faces) == 0 or faces[i % len(faces)][3] == 0
            src = edges if hollow else faces
            colour = mcolors.to_hex(src[i % len(src)]).upper() if len(src) else "NONE"
            out.append((float(x), float(y), hollow, colour))
    return out


def text_overlaps(fig) -> list[tuple[str, str]]:  # type: ignore[no-untyped-def]
    """Pairs of visible, non-empty texts whose rendered boxes overlap by more than 1 px, measured
    with ``FigureCanvasAgg`` after a draw (``FigureCanvasBase`` gave false passes)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.text import Text

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    tick_labels: set[int] = set()
    drawn: list[Text] = []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                tick_labels.update((id(tick.label1), id(tick.label2)))
            if not (ax.get_visible() and ax.axison and axis.get_visible()):
                continue  # an axes switched off draws no tick labels
            # the ticks inside the view interval are the ones drawn
            for tick in axis._update_ticks():
                drawn += [lab for lab in (tick.label1, tick.label2) if lab.get_visible()]
    texts = [
        t
        for t in fig.findobj(Text)
        if id(t) not in tick_labels and t.get_visible() and t.get_figure() is not None
    ]
    texts = [t for t in [*texts, *drawn] if t.get_text().strip()]
    boxes = [(t.get_text(), t.get_window_extent(renderer).expanded(1.0, 1.0)) for t in texts]
    hits = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (s1, b1), (s2, b2) = boxes[i], boxes[j]
            w = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
            h = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
            if w > 1.0 and h > 1.0:
                hits.append((s1, s2))
    return hits
