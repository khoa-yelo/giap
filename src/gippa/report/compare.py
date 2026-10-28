"""``CMP/report/compare.html`` (GIPPA_SPEC.md §7.9.3, §10.3; M9-23).

Pure: the run labels, directories (relative to ``CMP/``), genome SHA-256 and detection tiers,
the V16 plot files, the three ``compare_*`` frames and the optional V16 Plotly twin in, one HTML
string out. Sections ``compare_versions`` (only when the runs' versions differ, M11-14h),
``runs``, ``detection`` (V16 inlined, the intersection rows), ``abundance`` (the concordance
twin, ``compare_abundance`` with its Bland-Altman columns), ``stability`` and ``provenance``;
the run report's CSS, table script and Plotly embedding (M9-02, M9-14). Every link is relative
and names a file under ``CMP/``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pandas as pd

import gippa
from gippa.io.plot_files import PlotFiles
from gippa.io.tsv import format_value
from gippa.report.build import _BANNER, _ENV, Twin, _plot_slot, _plotlyjs

#: intersection rows shown (the V16 bars); the full table is linked
TOP_INTERSECTIONS = 30


@dataclass(frozen=True)
class RunInfo:
    label: str
    directory: str
    genome_sha256: str
    detection_tier: str


def _table(key: str, frame: pd.DataFrame, limit: int | None = None) -> dict[str, Any]:
    shown = frame if limit is None else frame.head(limit)
    rows = [[format_value(v) for v in r] for r in shown.astype(object).itertuples(index=False)]
    return {"key": key, "name": f"{key}.tsv", "header": [str(c) for c in frame.columns],
            "rows": rows, "href": f"../tables/{key}.tsv"}  # fmt: skip


def build_compare(
    runs: Sequence[RunInfo],
    v16: PlotFiles | None,
    tables: Mapping[str, pd.DataFrame],
    twin: Twin | None = None,
    versions: Sequence[str] = (),
) -> str:
    """The comparison page; deterministic, offline, no absolute path. ``versions`` are the
    ``compare_versions`` warning lines (M11-14h), shown in their own section when any."""
    det = tables["compare_detection"]
    inter = det[det["section"] == "intersection"][["intersection", "n_genes"]]
    twins = {} if twin is None or v16 is None else {i.stem: twin for i in v16.instances}
    slot = _plot_slot("V16", v16, {}, {})
    live = _plot_slot("V16", v16, {}, twins) if twins else None
    plotlyjs = _plotlyjs() if twins else None
    version = None if plotlyjs is None else _BANNER.search(plotlyjs)
    return _ENV.get_template("compare.html.j2").render(
        runs=runs, slot=slot, live=live, version=gippa.__version__,
        detection=_table("compare_detection", inter, TOP_INTERSECTIONS),
        abundance=_table("compare_abundance", tables["compare_abundance"]),
        stability=_table("compare_stability", tables["compare_stability"]),
        plotlyjs=plotlyjs, plotly_version=version and version.group(1), versions=list(versions),
    )  # fmt: skip
