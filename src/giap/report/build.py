"""The static ``report/report.html`` (GIAP_SPEC.md §10.3 bullets 1, 2, 4; M9-02 behaviour 2-6).

Pure: a ``RunView``, the plot files read by ``giap.io.plot_files.read_plot_files`` and the skip
reasons in, one HTML string out. The 13 ``SECTIONS`` always exist; each shows its plots (SVG text
inlined, ids prefixed with the instance stem so no two SVGs share an id), then the first
``TOP_ROWS`` rows of its tables as written, with a link to the full TSV. The template inlines its
CSS; every link is relative (``../plots/...``, ``../tables/...``) and names a file that exists.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import jinja2
import pandas as pd
from markupsafe import Markup

import giap
from giap.config import GiapConfig
from giap.io.plot_files import PlotFiles
from giap.io.run_view import PATHS, RunView
from giap.io.tsv import format_value

TOP_ROWS = 20
UNESTIMATED = "unestimated"


@dataclass(frozen=True)
class Section:
    id: str
    heading: str
    plots: tuple[str, ...] = ()
    tables: tuple[str, ...] = ()


SECTIONS: tuple[Section, ...] = (
    Section("summary", "Summary"),
    Section("inputs-qc", "Inputs and quality control", ("V17",)),
    Section("mapping", "Peptide mapping", ("V05", "V06")),
    Section("detection", "Gene detection", ("V07", "V18"), ("genes_tsv", "detection_summary")),
    Section(
        "genome-context",
        "Genome context",
        ("V01", "V02", "V10", "V03"),
        ("islands", "regions_stats", "autocorrelation"),
    ),
    Section("operons", "Operons", ("V11",), ("pairs", "operon_flags")),
    Section("variants", "Variants", ("V12",)),
    Section("novel-orfs", "Novel ORFs", ("V13",), ("novel_orfs",)),
    Section("function", "Function", ("V08", "V14")),
    Section("determinants", "Determinants of detection", ("V15",)),
    Section("sensitivity", "Sensitivity", (), ("sensitivity",)),
    Section("methods", "Methods"),
    Section("provenance", "Provenance"),
)

_ENV = jinja2.Environment(
    loader=jinja2.PackageLoader("giap.report", "templates"),
    autoescape=jinja2.select_autoescape(["html", "j2"]),
    keep_trailing_newline=True,
    trim_blocks=True,
    lstrip_blocks=True,
)
_ID = re.compile(r'(\sid=")([^"]*)"')
_HREF = re.compile(r'(href=")#([^"]*)"')
_URL = re.compile(r"url\(\s*#([^)\s]+)\s*\)")


def embed_svg(text: str, prefix: str) -> str:
    """``text`` from its ``<svg`` on (XML prolog and doctype dropped), every id and local
    reference (``href="#x"``, ``url(#x)``) prefixed with ``<prefix>-``."""
    start = text.find("<svg")
    body = text[start:] if start >= 0 else ""
    body = _ID.sub(lambda m: f'{m.group(1)}{prefix}-{m.group(2)}"', body)
    body = _HREF.sub(lambda m: f'{m.group(1)}#{prefix}-{m.group(2)}"', body)
    return _URL.sub(lambda m: f"url(#{prefix}-{m.group(1)})", body)


def reason_words(reason: str) -> str:
    """``module_off:function`` -> "the function module is off" (behaviour 2)."""
    kind, _, what = reason.partition(":")
    words = {
        "module_off": f"the {what} module is off",
        "input_absent": f"input {what} absent",
        "layer_missing": f"layer {what} missing",
    }
    return words.get(kind, reason)


def _plot_slot(pid: str, files: PlotFiles | None, skipped: Mapping[str, str]) -> dict[str, Any]:
    if files is None or not files.instances:
        note = f"not run: {reason_words(skipped[pid])}" if pid in skipped else "plot file absent"
        return {"pid": pid, "note": note, "instances": []}
    insts = []
    for inst in files.instances:
        links = [(rel.rsplit(".", 1)[1].upper(), f"../plots/{rel}") for rel in inst.images]
        if inst.data is not None:
            links.append(("data TSV", f"../plots/{inst.data}"))
        svg = None if inst.svg is None else Markup(embed_svg(inst.svg, inst.stem))
        insts.append({"stem": inst.stem, "svg": svg, "links": links})
    absent = None if any(i["svg"] is not None for i in insts) else "plot file absent"
    return {"pid": pid, "note": absent, "instances": insts}


def _table(view: RunView, key: str) -> dict[str, Any]:
    head = view.head(key, TOP_ROWS)
    name = PATHS[key].rsplit("/", 1)[1]
    if head is None:
        return {"key": key, "name": name, "absent": True}
    header, rows = head
    return {"key": key, "name": name, "header": header, "rows": rows, "href": f"../{PATHS[key]}"}


def _count(frame: pd.DataFrame, column: str, alpha: float | None = None) -> int:
    col = frame[column]
    if alpha is None:
        return int(col.fillna(False).astype(bool).sum())
    return int((pd.to_numeric(col, errors="coerce") <= alpha).sum())


def _tier_fdr(view: RunView, alpha: float) -> str | None:
    ds = view.table("detection_summary")
    if ds is None:
        return None
    hit = ds[(ds["section"] == "tier") & (ds["group_type"] == "q_threshold")
             & (ds["group"] == format_value(alpha)) & (ds["metric"] == "fdr")]  # fmt: skip
    vals = hit["value"].dropna()
    return None if vals.empty else format_value(float(vals.iloc[0]))


def cards(view: RunView, cfg: GiapConfig) -> list[dict[str, str]]:
    """The pass-1 summary cards (behaviour 5): a number and its error estimate or
    ``unestimated``; a card whose table is absent is left out."""
    a = cfg.stats.alpha
    out: list[dict[str, str]] = []

    def add(name: str, label: str, value: int | str, estimate: str | None) -> None:
        out.append({"name": name, "label": label, "value": str(value),
                    "estimate": estimate or UNESTIMATED})  # fmt: skip

    genes = view.table("genes_tsv")
    if genes is not None:
        tiers = ", ".join(sorted({str(t) for t in genes["detection_tier"].dropna()}))
        fdr = _tier_fdr(view, a.identification)
        est = None if fdr is None else f"tier FDR {fdr} at q ≤ {format_value(a.identification)}"
        add("genes_detected", f"genes detected of {len(genes)} ({tiers or 'no tier'})",
            _count(genes, "detected"), est)  # fmt: skip
    hits = view.table("peptide_hits")
    if hits is not None:
        placed = hits[pd.to_numeric(hits["placement_q"], errors="coerce") <= a.placement]
        unmapped = view.table("unmapped")
        n_in = hits["peptide"].nunique() + (0 if unmapped is None else len(unmapped))
        add("peptides_placed", f"peptides placed of {n_in} input", placed["peptide"].nunique(),
            f"mapping FDR controlled at placement q ≤ {format_value(a.placement)}")  # fmt: skip
    islands = view.table("islands")
    if islands is not None:
        add("islands", "abundance islands", _count(islands, "p_value", a.context),
            f"FWER-controlled at p ≤ {format_value(a.context)}")  # fmt: skip
    flags = view.table("operon_flags")
    if flags is not None:
        add("operon_flags", "flagged operon pairs", _count(flags, "q_value", a.operon_flags),
            f"FDR-controlled at q ≤ {format_value(a.operon_flags)}")  # fmt: skip
    novel = view.table("novel_orfs")
    if novel is not None:
        fdrs = pd.to_numeric(novel["class_fdr"], errors="coerce").dropna()
        est = None if fdrs.empty else f"class FDR ≤ {format_value(float(fdrs.max()))}"
        add("novel_orfs", "novel ORFs called", _count(novel, "called"), est)
    return out


def _provenance(view: RunView) -> dict[str, Any]:
    pre = view.preflight() or {}
    inputs = [{"file": str(i.get("file", "")).replace("\\", "/").rsplit("/", 1)[-1],
               "role": str(i.get("role", "")), "sha256": str(i.get("sha256", ""))}
              for i in pre.get("inputs", []) if isinstance(i, Mapping)]  # fmt: skip
    return {"version": giap.__version__, "inputs": inputs, "config": view.digest("resolved_config")}


def build_report(
    view: RunView, plots: Mapping[str, PlotFiles], skipped: Mapping[str, str], cfg: GiapConfig
) -> str:
    """The report page (behaviour 2-6); deterministic, offline, no absolute path."""
    sections = [
        {"id": s.id, "heading": s.heading,
         "plots": [_plot_slot(p, plots.get(p), skipped) for p in s.plots],
         "tables": [_table(view, k) for k in s.tables]}
        for s in SECTIONS
    ]  # fmt: skip
    return _ENV.get_template("report.html.j2").render(
        sections=sections, cards=cards(view, cfg), provenance=_provenance(view),
        version=giap.__version__,
    )  # fmt: skip
