# ruff: noqa: RUF001  (the report's wording is Unicode)
"""The static ``report/report.html`` (GIAP_SPEC.md §10.3 bullets 1, 2, 4; M9-02 behaviour 2-6).

Pure: a ``RunView``, the plot files read by ``giap.io.plot_files.read_plot_files`` and the skip
reasons in, one HTML string out. The 13 ``SECTIONS`` always exist; each shows its plots (SVG text
inlined, ids prefixed with the instance stem so no two SVGs share an id), then the first
``TOP_ROWS`` rows of its tables as written, with a link to the full TSV. The template inlines its
CSS and ``tables.js``; ``twins`` (instance stem -> ``(div id, figure JSON)``, M9-14) adds one Plotly
div per instance beside its static SVG and inlines plotly.js once in ``<head>`` (only when a twin
exists). Every link is relative (``../plots/...``, ``../tables/...``) and names a file that exists.
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
from giap.report import content

TOP_ROWS = 20
UNESTIMATED = "unestimated"
#: a Plotly twin: ``(div id, figure JSON)``
Twin = tuple[str, str]
_BANNER = re.compile(r"plotly\.js v(\S+)")


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


def _script_json(text: str) -> str:
    """JSON that cannot end its ``<script>`` element (``<`` written as ``\\u003c``)."""
    return text.replace("<", "\\u003c")


def reason_words(reason: str) -> str:
    """``module_off:function`` -> "the function module is off" (behaviour 2)."""
    kind, _, what = reason.partition(":")
    words = {
        "module_off": f"the {what} module is off",
        "input_absent": f"input {what} absent",
        "layer_missing": f"layer {what} missing",
        "not_registered": "this plot is not available in this GIAP version",
    }
    return words.get(kind, reason)


def _plot_slot(
    pid: str, files: PlotFiles | None, skipped: Mapping[str, str], twins: Mapping[str, Twin]
) -> dict[str, Any]:
    if files is None or not files.instances:
        note = f"not run: {reason_words(skipped[pid])}" if pid in skipped else "plot file absent"
        return {"pid": pid, "note": note, "instances": []}
    insts = []
    for inst in files.instances:
        links = [(rel.rsplit(".", 1)[1].upper(), f"../plots/{rel}") for rel in inst.images]
        if inst.data is not None:
            links.append(("data TSV", f"../plots/{inst.data}"))
        svg = None if inst.svg is None else Markup(embed_svg(inst.svg, inst.stem))
        twin = twins.get(inst.stem)
        spec = None if twin is None else {"id": twin[0], "spec": Markup(_script_json(twin[1]))}
        insts.append({"stem": inst.stem, "svg": svg, "links": links, "twin": spec})
    absent = None if any(i["svg"] is not None for i in insts) else "plot file absent"
    return {"pid": pid, "note": absent, "instances": insts}


def _table(view: RunView, key: str) -> dict[str, Any]:
    name = PATHS[key].rsplit("/", 1)[1]
    try:
        head = view.head(key, TOP_ROWS)
        whole = view.table(key) if key == "sensitivity" and head is not None else None
    except Exception as exc:  # an unreadable file is a note, never a crash (review B5)
        return {"key": key, "name": name, "error": _unreadable(key, exc)}
    if head is None:
        return {"key": key, "name": name, "absent": True}
    header, rows = head
    found: dict[str, Any] = {
        "key": key,
        "name": name,
        "header": header,
        "rows": rows,
        "href": f"../{PATHS[key]}",
    }
    if whole is not None and len(whole):  # rendered whole, by family (M9-16 behaviour 5)
        found |= {"header": list(whole.columns), "groups": content.sensitivity_groups(whole)}
    return found


def _unreadable(key: str, exc: BaseException) -> str:
    return f"{PATHS[key].rsplit('/', 1)[1]}: could not be read: {type(exc).__name__}: {exc}"


def _count(frame: pd.DataFrame, column: str, alpha: float | None = None) -> int:
    col = frame[column]
    if alpha is None:
        return int(col.fillna(False).astype(bool).sum())
    return int((pd.to_numeric(col, errors="coerce") <= alpha).sum())


def _tier_fdr(ds: pd.DataFrame | None, alpha: float) -> str | None:
    if ds is None:
        return None
    hit = ds[(ds["section"] == "tier") & (ds["group_type"] == "q_threshold")
             & (ds["group"] == format_value(alpha)) & (ds["metric"] == "fdr")]  # fmt: skip
    vals = hit["value"].dropna()
    return None if vals.empty else format_value(float(vals.iloc[0]))


def _novel_fdr(novel: pd.DataFrame) -> str | None:
    """The max class FDR over CALLED rows; ``None`` (unestimated) if any called row is (B1).
    M7-17: a suffix counts the called rows flagged ``denovo_only_unsupported`` (ADR-0050)."""
    called = novel[novel["called"].fillna(False).astype(bool)]
    fdrs = pd.to_numeric(called["class_fdr"], errors="coerce")
    flags = called["flags"].fillna("").astype(str).str.split(";")
    k = int(flags.map(lambda f: "denovo_only_unsupported" in f).sum())
    suffix = (f"; {k} of {len(called)} calls de novo-only without database or RBS support "
              "(exact de novo errors are not covered by the class FDR)") if k else ""  # fmt: skip
    if called.empty or bool(fdrs.isna().any()):
        return UNESTIMATED + suffix if k else None
    return f"class FDR ≤ {format_value(float(fdrs.max()))}{suffix}"


_TIER_FDR = {"decoy_picked": "picked target–decoy FDR", "mayu": "Mayu FDR"}
TWO_PEPTIDE = "(two-peptide rule; FDR unestimated)"


def _alt_estimate(alt: pd.DataFrame) -> str | None:
    """The ``class_error_rate`` range over the numeric rows, "(N unestimated)" appended (lead
    decision 5); ``None`` when no row has a number."""
    rates = pd.to_numeric(alt["class_error_rate"], errors="coerce")
    known, n_un = rates.dropna(), int(rates.isna().sum())
    if known.empty:
        return None
    lo, hi = format_value(float(known.min())), format_value(float(known.max()))
    tail = f" ({n_un} unestimated)" if n_un else ""
    return f"class error rate {lo if lo == hi else f'{lo}–{hi}'} (pooled per class){tail}"


def _metric(en: pd.DataFrame, metric: str) -> tuple[float | None, list[str]]:
    rows = en[en["metric"] == metric]
    vals = pd.to_numeric(rows["value"], errors="coerce").dropna()
    flags = [f for v in rows["flags"] if isinstance(v, str) for f in v.split(";")]
    return (None if vals.empty else float(vals.iloc[0])), flags


def _strata(fdr: pd.DataFrame | None) -> list[str]:
    """The descriptive placement rows of ``mapping_fdr``, each with its full key (decoy type,
    differences, length class); they are not q-values (review B1 on #228, lead decision 15)."""
    if fdr is None:
        return []
    cols = ["decoy_type", "n_diff", "length_class", "fdr_stratum"]
    rows = fdr[fdr["section"] == "placement"].reindex(columns=cols)
    return [f"{d} {format_value(n)} diff {lc}: {format_value(v)}"
            for d, n, lc, v in rows.itertuples(index=False) if not pd.isna(v)]  # fmt: skip


def cards(view: RunView, cfg: GiapConfig) -> tuple[list[dict[str, Any]], list[str]]:
    """The summary cards (M9-02 behaviour 5, M9-16 behaviour 1): a number and its error estimate
    or ``unestimated``, and under ``sources`` the ``<path>:<column>`` entries read (decision 3);
    a card whose table is absent is left out. Also the notes naming each table that could not be
    read (its card left out, the other cards kept; review B5)."""
    a = cfg.stats.alpha
    out: list[dict[str, Any]] = []
    notes: dict[str, str] = {}

    def table(key: str) -> pd.DataFrame | None:
        try:
            return view.table(key)
        except Exception as exc:
            notes.setdefault(key, _unreadable(key, exc))
            return None

    def src(key: str, *cols: str) -> list[str]:
        return [f"{PATHS[key]}:{c}" for c in cols] if view.has(key) else []

    def add(name: str, label: str, value: int | str, estimate: str | None, sources: list[str],
            **extra: str) -> None:  # fmt: skip
        out.append({"name": name, "label": label, "value": str(value),
                    "estimate": estimate or UNESTIMATED, "sources": sources, **extra})  # fmt: skip

    def n_unmapped(kind: str) -> int:
        unmapped = table("unmapped")
        return 0 if unmapped is None else int((unmapped["kind"] == kind).sum())

    def genes(t: pd.DataFrame) -> None:
        tiers = sorted({str(x) for x in t["detection_tier"].dropna()})
        tier = tiers[0] if len(tiers) == 1 else None
        sources = src("genes_tsv", "detected", "detection_tier", "detected_2pep")
        n, total = _count(t, "detected"), f"{len(t)} genes"
        if tier == "two_peptide":  # §7.2: the rule has no FDR (lead decision 6)
            add("genes_detected", "genes detected", n, TWO_PEPTIDE, sources, detail=total)
            return
        fdr = None if tier is None else _tier_fdr(table("detection_summary"), a.identification)
        est = None
        if fdr is not None:
            sources += src("detection_summary", "value")
            what = _TIER_FDR.get(str(tier), "tier FDR")
            est = f"({what} {fdr} at q ≤ {format_value(a.identification)})"
        two = f"{_count(t, 'detected_2pep')} also pass the two-peptide rule; {total}"
        add("genes_detected", "genes detected", n, est, sources, detail=two)

    def peptides(hits: pd.DataFrame) -> None:
        placed = hits[pd.to_numeric(hits["placement_q"], errors="coerce") <= a.placement]
        n_in = hits["peptide"].nunique() + n_unmapped("peptide")
        strata = _strata(table("mapping_fdr"))
        q = format_value(a.placement)
        est = f"placement q ≤ {q} per family (0 / ≥ 1 difference, per decoy type)"
        if strata:
            est += f"; descriptive stratum FDR (not q-values): {', '.join(strata)}"
        sources = src("peptide_hits", "peptide", "placement_q") + src("unmapped", "kind")
        keys = ("decoy_type", "n_diff", "length_class", "fdr_stratum")
        sources += src("mapping_fdr", *keys) if strata else []
        add("peptides_placed", f"peptides placed of {n_in} input", placed["peptide"].nunique(),
            est, sources)  # fmt: skip

    def proteins(hits: pd.DataFrame) -> None:
        n = hits["protein_id"].nunique()
        add("proteins_placed", f"proteins placed of {n + n_unmapped('protein')} input", n, None,
            src("protein_hits", "protein_id") + src("unmapped", "kind"))  # fmt: skip

    def islands(t: pd.DataFrame) -> None:
        add("islands", "abundance islands", _count(t, "p_value", a.context),
            f"FWER ≤ {format_value(a.context)}", src("islands", "p_value"))  # fmt: skip

    def operons(t: pd.DataFrame) -> None:
        q = format_value(a.operon_flags)
        add(
            "operon_flags",
            f"operon pairs flagged at q ≤ {q} (discordant or co-abundant)",
            _count(t, "q_value", a.operon_flags),
            f"FDR ≤ {q}",
            src("operon_flags", "q_value"),
        )

    def variants(t: pd.DataFrame) -> None:
        alt = t[t["call"] == "ALT"]  # ``mixed`` is not an ALT call (lead decision 5)
        counts = alt["confidence"].astype(str).value_counts()
        known = [c for c in ("high", "candidate") if c in counts.index]
        order = known + sorted(c for c in counts.index if c not in known)
        by = ", ".join(f"{int(counts[c])} {c}" for c in order)
        add("alt_calls", f"ALT calls ({by or 'none'})", len(alt), _alt_estimate(alt),
            src("variants", "call", "confidence", "class_error_rate"))  # fmt: skip

    def novel(t: pd.DataFrame) -> None:
        add("novel_orfs", "novel ORFs called", _count(t, "called"), _novel_fdr(t),
            src("novel_orfs", "called", "class_fdr"))  # fmt: skip

    def model(t: pd.DataFrame) -> None:
        """ADR-0033: the elastic-net whole-model rows only; ADR-0034: "(uncalibrated)"."""
        en = t[(t["model"] == "elastic_net") & (t["feature_or_block"] == "all")]
        (r2, _), (sd, _) = _metric(en, "cv_r2_mean"), _metric(en, "cv_r2_sd")
        p, flags = _metric(en, "null_p")
        value = "NA" if r2 is None else f"{format_value(r2)} ± {format_value(sd)}"
        unc = " (uncalibrated)" if "null_uncalibrated" in flags else ""
        est = None if p is None else f"null p {format_value(p)}{unc}"
        add("model_r2", "(CV);", value, est,
            src("determinants_model", "model", "feature_or_block", "metric", "value", "flags"),
            prefix="Model R²")  # fmt: skip

    specs = [("genes_tsv", genes), ("peptide_hits", peptides), ("protein_hits", proteins),
             ("islands", islands), ("operon_flags", operons), ("variants", variants),
             ("novel_orfs", novel), ("determinants_model", model)]  # fmt: skip
    for key, card in specs:
        frame = table(key)
        if frame is None:
            continue
        try:
            card(frame)
        except Exception as exc:
            notes.setdefault(key, _unreadable(key, exc))
    return out, list(notes.values())


def build_report(
    view: RunView,
    plots: Mapping[str, PlotFiles],
    skipped: Mapping[str, str],
    cfg: GiapConfig,
    twins: Mapping[str, Twin] | None = None,
) -> str:
    """The report page (behaviour 2-6); deterministic, offline, no absolute path."""
    twins = twins or {}
    sections = [
        {"id": s.id, "heading": s.heading,
         "plots": [_plot_slot(p, plots.get(p), skipped, twins) for p in s.plots],
         "tables": [_table(view, k) for k in s.tables]}
        for s in SECTIONS
    ]  # fmt: skip
    summary, notes = cards(view, cfg)
    doc, note = content.read_manifest(view)
    shown = any(i.stem in twins for s in SECTIONS for p in s.plots if p in plots
                for i in plots[p].instances)  # fmt: skip
    plotlyjs = _plotlyjs() if shown else None
    version = None if plotlyjs is None else _BANNER.search(plotlyjs)
    return _ENV.get_template("report.html.j2").render(
        sections=sections, cards=summary, notes=notes, version=giap.__version__,
        provenance=content.provenance(view, doc, note), extras=content.extras(view, doc),
        snippets=content.snippets(), plotlyjs=plotlyjs and Markup(plotlyjs),
        plotly_version=version and version.group(1),
    )  # fmt: skip


def _plotlyjs() -> str:
    """plotly.js as bundled in the plotly package (no network)."""
    import plotly.offline

    text: str = plotly.offline.get_plotlyjs()
    return text
