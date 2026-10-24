"""Data-TSV frames and figure-JSON readers for the M9-15 twin tests (interactive Plotly twins of
V2, V5, V6, V7, V8, V10, V11, V17, V18; GIAP_SPEC.md §10.1, §10.2 "Interactive in report").

A twin sees only its plot's data TSV, read back the way ``giap.steps.plots.twins`` reads it
(``read_csv(sep="\\t", na_values=["NA"], keep_default_na=False)``). ``frame`` therefore writes
the rows with the TSV writer's own cell formatter (``giap.io.tsv.format_value``) and reads them
back with those flags, so a frame here is byte-for-byte what a twin gets from a file. Coordinates
in these frames are the TSV's 1-based closed ones.

``spec`` parses ``figure_json``; ``values`` decodes Plotly's typed arrays (``{"dtype", "bdata"}``)
so tests do not depend on whether a twin passes lists or numpy arrays. Nothing here touches the
network or ``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import base64
import io
import json
import re
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

#: §10.2 rows whose "Interactive in report" column is Plotly (V16: the concordance view)
PLOTLY_PLOTS = ("V02", "V05", "V06", "V07", "V08", "V10", "V11", "V14", "V15", "V16", "V17", "V18")
#: twins landing with their plots (M9-19, M9-20, M9-23), not in M9-15
LATER = ("V14", "V15", "V16")
M915 = tuple(p for p in PLOTLY_PLOTS if p not in LATER)
STATIC_ONLY = ("V01", "V03", "V12", "V13")

V02_COLUMNS = ("contig", "kind", "locus_tag", "product", "position_bp", "start", "end", "value",
               "category")  # fmt: skip
V05_COLUMNS = ("input_class", "route", "outcome", "n")
V06_COLUMNS = ("panel", "length_class", "n_diff", "score_bin", "metric", "value", "ci_low",
               "ci_high", "flags")  # fmt: skip
V07_COLUMNS = ("block", "group", "region_type", "n", "observed", "expected", "o_minus_e",
               "ci_low", "ci_high", "q_value", "flags")  # fmt: skip
V08_COLUMNS = ("category", "super_category", "locus_tag", "abundance_log2", "detected",
               "median_log2", "n_genes", "n_detected", "n_imputed")  # fmt: skip
V10_COLUMNS = ("contig", "statistic", "null_model", "lag", "value", "null_mean", "null_lo",
               "null_hi", "effect", "q_value", "significant", "flags")  # fmt: skip
V11_COLUMNS = ("panel", "pair_id", "class", "abs_delta_log2", "lfdr", "upstream_log2",
               "downstream_log2", "flag", "q_value", "bin", "density", "n")  # fmt: skip
V17_COLUMNS = ("panel", "group", "key", "value", "n")
V18_COLUMNS = ("panel", "table", "family", "bin", "x", "y", "n", "metric", "value")

Row = Sequence[object]


# ---------------------------------------------------------------- frames as the step reads them


def tsv_text(columns: Sequence[str], rows: Iterable[Row]) -> str:
    from giap.io.tsv import format_value

    lines = ["\t".join(columns)]
    for r in rows:
        assert len(r) == len(columns), (r, columns)
        lines.append("\t".join(format_value(v) for v in r))
    return "\n".join(lines) + "\n"


def read(text: str) -> pd.DataFrame:
    """A data TSV's text as ``giap.steps.plots.twins`` reads the file."""
    from giap.steps.plots import TWIN_ID_DTYPES

    return pd.read_csv(io.StringIO(text), sep="\t", na_values=["NA"], keep_default_na=False,
                       dtype=TWIN_ID_DTYPES)  # fmt: skip


def frame(columns: Sequence[str], rows: Iterable[Row]) -> pd.DataFrame:
    return read(tsv_text(columns, rows))


def instance_frame(data: pd.DataFrame, coords: Sequence[str] = ()) -> pd.DataFrame:
    """A plot ``Instance``'s data as its written data TSV reads back (``coords`` + 1)."""
    out = data.copy()
    for c in coords:
        out[c] = out[c] + 1
    return frame([str(c) for c in out.columns], out.itertuples(index=False, name=None))


# ---------------------------------------------------------------- example frames per plot


def v02_rows(contig: str = "chr1") -> list[Row]:
    """chr1 (100 kb): two rug genes, four points (one with product NA), a trend, two islands."""
    return [
        (contig, "point", "g1", "DNA polymerase III subunit beta", 10_001, None, None, 1.5,
         "detected"),
        (contig, "point", "g2", None, 30_001, None, None, 2.25, "detected"),
        (contig, "point", "g3", "ribosomal protein L2", 50_001, None, None, 3.0, "detected"),
        (contig, "point", "g4", "hypothetical protein", 70_001, None, None, 2.0, "detected"),
        (contig, "rug", "g5", "transposase", 20_001, None, None, None, "not_detected"),
        (contig, "rug", "g6", None, 60_001, None, None, None, "single_peptide_unestimated"),
        (contig, "trend", "g1", None, 10_001, None, None, 1.75, None),
        (contig, "trend", "g2", None, 30_001, None, None, 2.25, None),
        (contig, "trend", "g3", None, 50_001, None, None, 2.5, None),
        (contig, "trend", "g4", None, 70_001, None, None, 2.25, None),
        (contig, "island", None, None, 45_001, 40_001, 50_000, None, "high"),
        (contig, "island", None, None, 85_001, 80_001, 90_000, None, "low"),
    ]  # fmt: skip


def v02() -> pd.DataFrame:
    return frame(V02_COLUMNS, v02_rows())


V05_EX1 = {
    ("peptide", "exact", "placed"): 1,
    ("peptide", "tolerant", "placed"): 1,
    ("peptide", "unmapped", "contaminant"): 1,
    ("peptide", "unmapped", "rejected by FDR"): 1,
    ("peptide", "unmapped", "unexplained"): 1,
}
V05_GRID = (
    ("peptide", "exact", "placed"),
    ("peptide", "tolerant", "placed"),
    ("peptide", "unmapped", "contaminant"),
    ("peptide", "unmapped", "rejected by FDR"),
    ("peptide", "unmapped", "unexplained"),
    ("protein", "miniprot", "placed"),
    ("protein", "tblastn", "placed"),
    ("protein", "unmapped", "contaminant"),
    ("protein", "unmapped", "rejected by FDR"),
    ("protein", "unmapped", "unexplained"),
)


def v05(counts: Mapping[tuple[str, str, str], int] = V05_EX1, conflicts: int = 0) -> pd.DataFrame:
    """The V05 data TSV (M9-05 grid order, then the conflicts row); default: M9-05 example 1."""
    rows = [(*k, int(counts.get(k, 0))) for k in V05_GRID] + [
        ("all", "all", "conflicts", conflicts)
    ]
    return frame(V05_COLUMNS, rows)


def v06(with_c: bool = True) -> pd.DataFrame:
    lc = "8"
    rows: list[Row] = [
        ("a", lc, 0, None, "decoys", 2, None, None, None),
        ("a", lc, 0, None, "targets", 120, None, None, None),
        ("a", lc, 1, None, "decoys", 1, None, None, None),
        ("a", lc, 1, None, "targets", 30, None, None, None),
        ("a", lc, 2, None, "decoys", None, None, None, None),
        ("a", lc, 2, None, "targets", None, None, None, None),
        ("b", "all", 0, None, "fdr_stratum", 0.004, None, None, None),
        ("b", "all", 0, None, "n_accepted", 118, None, None, None),
        ("b", "all", 1, None, "fdr_stratum", 0.03, None, None, None),
        ("b", "all", 1, None, "n_accepted", 0, None, None, None),
        ("b", "all", 2, None, "fdr_stratum", None, None, None, None),
        ("b", "all", 2, None, "n_accepted", 0, None, None, None),
    ]
    if with_c:
        rows += [
            ("c", None, None, 1, "precision", 0.9, 0.85, 0.94, None),
            ("c", None, None, 2, "precision", 0.95, 0.9, 0.98, "c_assumed"),
        ]
    return frame(V06_COLUMNS, rows)


#: M9-06 example 1 (regions) and example 3 (basic replicon row), as the V07 TSV writes them
V07_ROWS: list[Row] = [
    ("replicon", "chr1", None, 100, 52.0, 49.9, 2.1, -1.0, 5.2, None, None),
    ("feature_class", "cds", None, 90, 50.0, 50.0, 0.0, None, None, None, "ci_unestimated"),
    ("region", "prophage_1", "prophage", 25, 3.0, 7.75, -4.75, -7.1, -2.2, 0.003, None),
    ("region", "prophage_2", "prophage", 20, 9.0, 8.0, 1.0, -1.5, 3.5, 0.4, None),
]


def v07(rows: Sequence[Row] = V07_ROWS) -> pd.DataFrame:
    return frame(V07_COLUMNS, rows)


def v08_rows(measured: bool = True) -> list[Row]:
    """J: 60 measured genes (a violin), C: 5 measured (a strip), none: 1 of 2 measured."""
    j = [(f"j{i:03d}", 8.0 + 0.05 * i) for i in range(60)]
    c = [(f"c{i}", 4.0 + 0.5 * i) for i in range(5)]
    cats = [("J", "information", j), ("C", "metabolism", c), ("none", "none", [("n1", 1.0)])]
    rows: list[Row] = []
    for letter, sup, genes in cats:
        vals = sorted(v for _, v in genes)
        extra = [("n2", None)] if letter == "none" else []
        everything = sorted([*genes, *extra])
        med = float(np.median(vals)) if measured else None
        n_det = len(genes)
        rows.append((letter, sup, None, None, None, med, len(everything), n_det, 0))
        for tag, v in everything:
            det = v is not None
            rows.append((letter, sup, tag, v if measured else None, det, None, None, None, None))
    return rows


def v08(measured: bool = True) -> pd.DataFrame:
    return frame(V08_COLUMNS, v08_rows(measured))


def v10_rows(contig: str = "chr1", lags: Sequence[int] = (1, 2, 3, 4, 5)) -> list[Row]:
    rows: list[Row] = []
    for stat in ("morans_i", "join_count_bb"):
        for null in ("gene_label", "tu_block"):
            for lag in lags:
                if stat == "morans_i":
                    v = round(0.2 * 0.8**lag, 4)
                    mean = -0.002 if null == "gene_label" else -0.004
                    half = 0.05 if null == "gene_label" else 0.09
                else:
                    v, mean, half = float(40 - lag), 38.0 if null == "gene_label" else 37.5, 3.0
                q = 0.001 if (stat, null, lag) == ("morans_i", "tu_block", 1) else 0.5
                # effect = value - null_mean, so the plotted statistic (effect + null_mean) is value
                rows.append((contig, stat, null, lag, v, mean, round(mean - half, 4),
                             round(mean + half, 4), round(v - mean, 4), q, q <= 0.05,
                             "ci_unestimated"))  # fmt: skip
    return rows


def v10(rows: Sequence[Row] | None = None) -> pd.DataFrame:
    return frame(V10_COLUMNS, v10_rows() if rows is None else rows)


def v10_not_estimable() -> pd.DataFrame:
    """Every row flagged ``not_estimable``: Moran's I NA, the join count written as 0."""
    rows: list[Row] = []
    for stat, v in (("morans_i", None), ("join_count_bb", 0.0)):
        for null in ("gene_label", "tu_block"):
            for lag in (1, 2, 3):
                rows.append(("chr1", stat, null, lag, v, v, v, v, v, None, False,
                             "ci_unestimated;not_estimable"))  # fmt: skip
    return frame(V10_COLUMNS, rows)


V11_OP_DENSITY = [round(1.0 - 0.03 * i, 4) for i in range(30)]
V11_NO_DENSITY = [round(0.1 + 0.02 * i, 4) for i in range(30)]


def v11() -> pd.DataFrame:
    n = None
    rows: list[Row] = []
    for cls, dens in (("OP", V11_OP_DENSITY), ("NO", V11_NO_DENSITY)):
        rows += [("a", n, cls, n, n, n, n, n, n, i, d, n) for i, d in enumerate(dens)]
        rows.append(("a", n, cls, n, n, n, n, n, n, n, n, 4))
    for cls in ("CV", "SW", "DV"):
        rows.append(("a", n, cls, n, n, n, n, n, n, n, n, 0))
    rows += [
        ("b", "P0001", "OP", 0.25, 0.05, n, n, n, n, n, n, n),
        ("b", "P0002", "OP", 0.5, 0.1, n, n, "internal_promoter_supported", 0.001, n, n, n),
        ("b", "P0003", "NO", 1.5, 0.6, n, n, n, n, n, n, n),
        ("b", "P0004", "NO", 3.0, 0.95, n, n, "novel_operon_candidate", 0.002, n, n, n),
        ("c", "P0001", "OP", 0.25, 0.05, 5.0, 5.25, n, n, n, n, n),
        ("c", "P0002", "OP", 0.5, 0.1, 6.0, 5.5, "internal_promoter_supported", 0.001, n, n, n),
        ("c", "P0003", "NO", 1.5, 0.6, 2.0, 3.5, n, n, n, n, n),
        ("c", "P0004", "NO", 3.0, 0.95, 7.0, 4.0, "novel_operon_candidate", 0.002, n, n, n),
    ]
    return frame(V11_COLUMNS, rows)


def v17() -> pd.DataFrame:
    n = None
    rows: list[Row] = [
        ("a", "level", "error", 1, n),
        ("a", "level", "info", 2, n),
        ("a", "level", "ok", 7, n),
        ("a", "level", "warn", 3, n),
        ("a", "message", "PF-E01", n, n),
        ("a", "message", "PF-W03", n, n),
        ("a", "message", "PF-W07", n, n),
        ("a", "message", "PF-W11", n, n),
        ("b", "chromosome", "chr1", 4_600_000, n),
        ("b", "plasmid", "p1", 50_000, n),
        ("c", "placed", 8, 3, n),
        ("c", "placed", 9, 5, n),
        ("c", "unmapped", 9, 2, n),
        ("d", "imputed", -2.0, 1, 1),
        ("d", "measured", 1.5, 4, 10),
        ("d", "measured", 2.0, 6, 10),
        ("e", "abundance", "abundance_raw", 0.1, 100),
        ("e", "codon usage", "cai", 0.0, 100),
        ("e", "protein", "mw_da", 0.25, 100),
    ]
    return frame(V17_COLUMNS, rows)


V18_FAMILIES = {
    ("regions_stats", "o_minus_e"): [3, 1, 0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 0, 1],
    ("autocorrelation", "morans_i / tu_block"): [7] + [0] * 18 + [1],
}


def v18_rows(families: Mapping[tuple[str, str], Sequence[int]] = V18_FAMILIES) -> list[Row]:
    n = None
    rows: list[Row] = [
        ("a", "detection_summary", n, 0, 0.05, 0.02, 40, n, n),
        ("a", "detection_summary", n, 1, 0.3, 0.25, 25, n, n),
        ("a", "detection_summary", n, 2, 0.7, 0.75, 30, n, n),
        ("a", "detection_summary", n, n, n, n, n, "brier", 0.12),
        ("a", "detection_summary", n, n, n, n, n, "calibration_slope", 0.95),
    ]
    for (table, family), counts in families.items():
        rows += [("b", table, family, k / 20, n, n, c, n, n) for k, c in enumerate(counts)]
        if (table, family) == ("regions_stats", "o_minus_e"):
            rows.append(("b", table, family, n, n, n, n, "n_outside_0_1", 2))
    return rows


def v18(families: Mapping[tuple[str, str], Sequence[int]] = V18_FAMILIES) -> pd.DataFrame:
    return frame(V18_COLUMNS, v18_rows(families))


EXAMPLES: Mapping[str, Callable[[], pd.DataFrame]] = {
    "V02": v02,
    "V05": v05,
    "V06": v06,
    "V07": v07,
    "V08": v08,
    "V10": v10,
    "V11": v11,
    "V17": v17,
    "V18": v18,
}
HEADERS: Mapping[str, Sequence[str]] = {
    "V02": V02_COLUMNS,
    "V05": V05_COLUMNS,
    "V06": V06_COLUMNS,
    "V07": V07_COLUMNS,
    "V08": V08_COLUMNS,
    "V10": V10_COLUMNS,
    "V11": V11_COLUMNS,
    "V17": V17_COLUMNS,
    "V18": V18_COLUMNS,
}


# ---------------------------------------------------------------- building and reading twins


def load() -> Any:
    """``giap.viz.interactive`` after every ``giap.viz`` module is imported (the ``report`` step
    discovers plots and twins the same way)."""
    from giap.steps.plots import load_plot_modules

    load_plot_modules()
    import giap.viz.interactive as interactive

    return interactive


def figure_text(pid: str, data: pd.DataFrame, title: str | None = None) -> str | None:
    text: str | None = load().figure_json(pid, data, title or f"{pid} twin")
    return text


def spec(pid: str, data: pd.DataFrame, title: str | None = None) -> dict[str, Any]:
    text = figure_text(pid, data, title)
    assert isinstance(text, str), f"{pid}: no twin JSON for a non-empty data TSV"
    out: dict[str, Any] = json.loads(text)
    return out


def values(v: Any) -> list[Any]:
    """A Plotly data array as a list: typed arrays decoded, scalars wrapped, ``None`` -> []."""
    if v is None:
        return []
    if isinstance(v, dict) and "bdata" in v:
        arr = np.frombuffer(base64.b64decode(v["bdata"]), dtype=np.dtype(v["dtype"]))
        if "shape" in v:
            arr = arr.reshape([int(s) for s in str(v["shape"]).split(",")])
        return arr.tolist()  # type: ignore[no-any-return]
    if isinstance(v, list):
        return v
    return [v]


def floats(v: Any) -> list[float]:
    return [float("nan") if x is None else float(x) for x in values(v)]


def traces(s: Mapping[str, Any], *types: str) -> list[dict[str, Any]]:
    """The figure's traces (not the template's), optionally of the given types."""
    out = [t for t in s.get("data", []) if not types or t.get("type", "scatter") in types]
    return out


def ttype(t: Mapping[str, Any]) -> str:
    return str(t.get("type", "scatter"))


MINUS = "\N{MINUS SIGN}"
_SPACE = re.compile("[" + "".join(map(chr, (0x00A0, 0x2009, 0x202F, 0x2007))) + "]")


def norm(text: str) -> str:
    """Hover text compared on ASCII minus and plain spaces (U+2212 -> ``-``; thin/nbsp -> space)."""
    return _SPACE.sub(" ", text.replace(MINUS, "-"))


def hovers(t: Mapping[str, Any]) -> list[str]:
    """Per-point hover strings: ``hovertext`` else ``text`` (scalar strings broadcast)."""
    raw = t.get("hovertext", t.get("text"))
    if raw is None:
        return []
    vals = values(raw)
    n = max(len(values(t.get("x"))), len(values(t.get("y"))), 1)
    if isinstance(raw, str):
        vals = [raw] * n
    return [norm(str(v)) for v in vals]


def all_hovers(s: Mapping[str, Any]) -> list[str]:
    return [h for t in traces(s) for h in hovers(t)]


def axis_key(t: Mapping[str, Any], axis: str) -> str:
    """``"x"`` -> ``"xaxis"``, ``"x3"`` -> ``"xaxis3"`` for the trace's ``axis`` (x or y)."""
    ref = str(t.get(f"{axis}axis", axis))
    return f"{axis}axis{ref[1:]}"


def axis(s: Mapping[str, Any], key: str) -> dict[str, Any]:
    out: dict[str, Any] = s.get("layout", {}).get(key, {})
    return out


def title_of(obj: Mapping[str, Any]) -> str:
    t = obj.get("title")
    if isinstance(t, dict):
        return str(t.get("text", ""))
    return "" if t is None else str(t)


def axis_title(s: Mapping[str, Any], t: Mapping[str, Any], ax: str) -> str:
    return title_of(axis(s, axis_key(t, ax)))


def axis_titles(s: Mapping[str, Any]) -> set[str]:
    lay = s.get("layout", {})
    return {title_of(v) for k, v in lay.items() if re.match(r"^[xy]axis\d*$", k)} - {""}


def on_axes(s: Mapping[str, Any], xtitle: str, *types: str) -> list[dict[str, Any]]:
    """The traces drawn on the x axis titled ``xtitle``."""
    return [t for t in traces(s, *types) if axis_title(s, t, "x") == xtitle]


def mode(t: Mapping[str, Any]) -> str:
    return str(t.get("mode", "lines+markers" if ttype(t) == "scatter" else ""))


def point_pairs(ts: Iterable[Mapping[str, Any]]) -> list[tuple[float, float]]:
    out = []
    for t in ts:
        out += list(zip(floats(t.get("x")), floats(t.get("y")), strict=False))
    return out


def annotations(s: Mapping[str, Any]) -> list[str]:
    return [norm(str(a.get("text", ""))) for a in s.get("layout", {}).get("annotations", [])]


# ---------------------------------------------------------------- colours


_HEX = re.compile(r"^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$")
_RGB = re.compile(r"^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*[\d.]+\s*)?\)$")


def to_hex(colour: str) -> str | None:
    """A colour string as upper-case ``#RRGGBB`` (alpha dropped); ``None`` when not a colour."""
    import matplotlib.colors as mcolors

    c = colour.strip()
    m = _RGB.match(c)
    if m:
        r, g, b = (round(float(x)) for x in m.groups())
        return f"#{r:02X}{g:02X}{b:02X}"
    if _HEX.match(c):
        c = c[:7] if len(c) == 9 else c
        return str(mcolors.to_hex(c)).upper()
    if c.lower() in mcolors.CSS4_COLORS:
        return str(mcolors.to_hex(c.lower())).upper()
    return None


def colours(obj: Any, key: str = "") -> Iterator[str]:
    """Every colour string in a figure JSON: values of ``*color*`` keys (scalars or lists),
    ``colorway`` and explicit ``colorscale`` stops, plus any hex/rgb string elsewhere."""
    if isinstance(obj, dict):
        if "bdata" in obj:
            return
        for k, v in obj.items():
            yield from colours(v, k)
    elif isinstance(obj, list):
        for v in obj:
            yield from colours(v, key)
    elif isinstance(obj, str):
        c = obj.strip()
        if _HEX.match(c) or _RGB.match(c) or ("color" in key.lower() and to_hex(c) is not None):
            yield c


def allowed_colours() -> set[str]:
    """§10.1: the category table and Okabe-Ito, viridis / cividis (256 samples), plus achromatic
    neutrals (black, white, greys: text, rules, the null bands)."""
    import matplotlib

    from giap.viz import style

    out = {c.upper() for c in style.OKABE_ITO.values()}
    out |= {c.upper() for c in style.DETECTION.values()}
    out |= {c.upper() for c in style.REGION.values()}
    out |= {c.upper() for c, _ in style.ISLAND.values()}
    out |= {c.upper() for c in style.COG_SUPER_COLOUR.values()}
    for name in ("viridis", "cividis"):
        cmap = matplotlib.colormaps[name]
        out |= {matplotlib.colors.to_hex(cmap(i / 255)).upper() for i in range(256)}
    import plotly.colors

    for scale in (plotly.colors.sequential.Viridis, plotly.colors.sequential.Cividis):
        out |= {h for c in scale if (h := to_hex(c)) is not None}
    return out


def achromatic(hex_colour: str) -> bool:
    return hex_colour[1:3] == hex_colour[3:5] == hex_colour[5:7]


def disallowed(s: Mapping[str, Any]) -> set[str]:
    ok = allowed_colours()
    bad = set()
    for c in colours(s):
        h = to_hex(c)
        if h is not None and h not in ok and not achromatic(h):
            bad.add(c)
    return bad
