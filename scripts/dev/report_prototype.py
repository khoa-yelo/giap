# ruff: noqa: RUF001, E501, E731, B023
"""Decision prototype of a simplified report.html (Khoa's request; targets 1.0.1). NOT shipped code.

    PYTHONPATH=src python scripts/dev/report_prototype.py OUT out.html [--titles]

Reads a finished run through the same readers and twin builders as ``gippa report``, then lays the
page out differently: key points first, one finding per section, grouped figures split into one
large plot per panel, everything else in collapsed detail. It computes no statistic: every number
is read from the run's tables, cards or finding titles. Long decimals are rounded for display only;
the full value stays in the cell's tooltip. Output links are absolute paths into OUT, so the page
is a local artefact: do not commit one built from pilot data.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.offline
import yaml

from gippa import __version__
from gippa.config import GippaConfig
from gippa.io.plot_files import read_plot_files, svg_title
from gippa.io.run_view import PATHS, RunView
from gippa.report import content
from gippa.report.build import SECTIONS, cards, embed_svg
from gippa.steps.plots import load_plot_modules, twin_meta, twins

TICK, AXIS, LEGEND, TITLE = 12, 13, 12, 14

#: how each grouped figure is split: panels matching ``show`` stay visible, the rest are collapsed.
#: "first" = the first panel (or the first replicon's panels); ``None`` = every panel.
PANELS: dict[str, dict[str, Any]] = {
    "V17": {"show": ["abundance", "peptide lengths"], "h": 400, "more_cols": 3, "more_h": 420},
    "V05": {"show": None, "cols": 1, "h": 500},
    "V06": {"show": None, "h": 430},
    "V07": {"show": None, "cols": 1, "h": 1040},
    "V18": {"show": ["reliability"], "h": 430, "more_cols": 3, "more_h": 270},
    "V02": {"show": "first", "cols": 1, "h": 430, "more_cols": 2},
    "V10": {"show": "first_replicon", "h": 420, "more_cols": 2},
    "V11": {"show": ["pair class", "partners"], "h": 440},
    "V08": {"show": None, "cols": 1, "h": 520, "top": 34, "legend": False},
    "V14": {"show": ["NES", "running"], "h": 520, "more_cols": 1, "more_h": 560},
    "V15": {"show": ["block", "observed"], "h": 450, "more_cols": 2, "more_h": 360},
}
#: table -> (sort column, ascending, key columns). Other columns stay behind "all columns".
TABLES: dict[str, tuple[str | None, bool, list[str]]] = {
    "genes_tsv": ("abundance_log2", False, ["locus_tag", "gene", "product", "detected", "detection_tier", "gene_q_value", "n_peptides_unique", "coverage_frac", "abundance_log2", "p_detect", "cog_category"]),
    "detection_summary": (None, True, ["section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n"]),
    "regions_stats": ("q_value", True, ["region_id", "region_type", "stat", "n_genes", "observed", "expected", "effect", "ci_low", "ci_high", "q_value"]),
    "islands": ("q_value", True, ["island_id", "direction", "contig", "start", "end", "n_genes", "mean_abundance_log2", "detected_frac", "q_value"]),
    "autocorrelation": ("q_value", True, ["contig", "statistic", "lag", "null_model", "value", "null_lo", "null_hi", "q_value"]),
    "pairs": ("q_value", True, ["gene_a", "gene_b", "class", "test", "intergenic_bp", "delta_log2", "q_value", "lfdr"]),
    "operon_flags": ("q_value", True, ["gene_a", "gene_b", "uniop_class", "abundance_pattern", "flag", "q_value", "lfdr"]),
    "novel_orfs": ("rank", True, ["orf_id", "class", "contig", "strand", "length_aa_best", "n_peptides", "best_peptide_q", "ostir_tir", "class_fdr", "called", "overlapping_gene"]),
    "variants": ("n_alt", False, ["gene", "aa_pos", "ref_aa", "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "class_error_rate"]),
}
LABELS = {
    "o_minus_e": "observed − expected", "q_value": "q", "abundance_log2": "log2 abundance",
    "p_detect": "expected P(detect)", "gene_q_value": "gene q", "coverage_frac": "coverage",
    "n_peptides_unique": "unique peptides", "cog_category": "COG", "delta_log2": "Δ log2",
    "n_genes": "genes", "ci_low": "CI low", "ci_high": "CI high", "lfdr": "local fdr",
    "n_ref": "ref peptides", "n_alt": "ALT peptides", "class_error_rate": "class error rate",
    "class_fdr": "class FDR", "best_peptide_q": "best peptide q", "ostir_tir": "OSTIR TIR",
    "length_aa_best": "length (aa)", "n_peptides": "peptides", "intergenic_bp": "gap (bp)",
    "null_lo": "null low", "null_hi": "null high", "mean_abundance_log2": "mean log2 abundance",
    "detected_frac": "detected fraction",
}
#: the sentence that opens each section: the finding title of this plot (generated from tables)
LEAD = {"inputs-qc": "V17", "mapping": "V05", "detection": "V07", "genome-context": "V01",
        "operons": "V11", "variants": "V12", "novel-orfs": "V13", "function": "V14",
        "determinants": "V15"}  # fmt: skip
NAV = {"summary": "Key points", "inputs-qc": "QC", "mapping": "Mapping", "detection": "Detection",
       "genome-context": "Genome", "operons": "Operons", "variants": "Variants",
       "novel-orfs": "Novel ORFs", "function": "Function", "determinants": "Determinants",
       "sensitivity": "Sensitivity", "methods": "Methods", "provenance": "Provenance"}  # fmt: skip
_LONG = re.compile(r"(?<![\w.])(\d+\.\d{4,}(?:[eE][+-]?\d+)?)(?![\w.])")
_LONGNUM = re.compile(r"^[+-]?\d*\.\d{4,}([eE][+-]?\d+)?$|^[+-]?\d+(\.\d+)?[eE][+-]?\d+$")
E = html.escape


def short(text: str) -> str:
    """Long decimals in a sentence rounded to 3 significant figures (display only)."""
    return _LONG.sub(lambda m: f"{float(m.group(1)):.3g}", text)


def cell(value: Any) -> str:
    raw = "NA" if value is None or (isinstance(value, float) and pd.isna(value)) else str(value)
    if re.fullmatch(r"-?\d+\.0", raw):
        raw = raw[:-2]
    if _LONGNUM.match(raw):
        return f'<td title="{E(raw)}">{float(raw):.3g}</td>'
    if len(raw) > 48:
        return f'<td title="{E(raw)}">{E(raw[:46])}…</td>'
    return f"<td>{E(raw)}</td>"


# --- figures: one Plotly figure per panel -----------------------------------------------------


def _suffix(ref: str) -> str:
    return ref[1:].split(" ")[0]


def explode(figure: dict[str, Any]) -> list[dict[str, Any]]:
    """One ``{title, data, layout}`` per subplot of a twin; blank panels are dropped."""
    lay, panels = figure["layout"], []
    keep = {k: lay[k] for k in KEEP_LAYOUT if k in lay}
    subtitles = [a for a in lay.get("annotations", []) if a.get("xref") == "paper" and a.get("yref") == "paper" and a.get("yanchor") == "bottom"]
    keys = sorted((k for k in lay if re.fullmatch(r"xaxis\d*", k)), key=lambda k: int(k[5:] or 1))
    for key in keys:
        xid = "x" + key[5:]
        x = dict(lay[key])
        yid = str(x.get("anchor", "y"))
        y = dict(lay.get("yaxis" + yid[1:], {}))
        traces = [dict(t) for t in figure["data"] if t.get("xaxis", "x") == xid and "domain" not in t]
        if not traces:
            continue
        centre, top = sum(x["domain"]) / 2 if "domain" in x else 0.5, (y.get("domain") or [0, 1])[1]
        title = ""
        for a in subtitles:
            if abs(a["x"] - centre) < 0.03 and abs(a["y"] - top) < 0.05:
                title = a["text"]
        for t in traces:
            t["xaxis"], t["yaxis"] = "x", "y"
        def local(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
            out = []
            for it in items:
                xr, yr = str(it.get("xref", "x")), str(it.get("yref", "y"))
                if xr.split(" ")[0] == xid or yr.split(" ")[0] == yid:
                    it = dict(it)
                    if xr.split(" ")[0] == xid:
                        it["xref"] = xr.replace(xid, "x", 1)
                    elif xr != "paper":
                        continue
                    if yr.split(" ")[0] == yid:
                        it["yref"] = yr.replace(yid, "y", 1)
                    elif yr != "paper":
                        continue
                    out.append(it)
            return out
        for ax in (x, y):
            for gone in ("domain", "anchor", "matches", "overlaying", "scaleanchor"):
                ax.pop(gone, None)
        panels.append({"title": title, "data": traces, "x": x, "y": y,
                       "annotations": local(lay.get("annotations", [])), "shapes": local(lay.get("shapes", [])), "keep": keep})  # fmt: skip
    for t in figure["data"]:
        if "domain" in t or t.get("type") in ("sankey", "table"):
            title = next((a["text"] for a in subtitles if t.get("type") == "table" and "warning" in a["text"]), "")
            panels.append({"title": title, "data": [{**t, "domain": {"x": [0, 1], "y": [0, 1]}}], "x": None, "y": None, "annotations": [], "shapes": [], "keep": {}})
    return panels


KEEP_LAYOUT = ("barmode", "barnorm", "bargap", "bargroupgap", "violinmode", "violingap", "boxmode", "hovermode")


def panel_figure(panel: dict[str, Any], template: Any, height: int, top: int = 12, legend: bool = True) -> str:
    """The panel as a figure at reading sizes (the twins use the print sizes, 6-7, as CSS px)."""
    named = any(t.get("showlegend", True) and t.get("name") for t in panel["data"])
    many = sum(1 for t in panel["data"] if t.get("showlegend", True) and t.get("name")) > 1
    for t in panel["data"]:
        if t.get("type") == "sankey":
            t["textfont"] = {"size": TICK}
    layout: dict[str, Any] = {
        "template": template, "height": height, "autosize": True,
        "font": {"size": TICK, "family": "DejaVu Sans, sans-serif"},
        "margin": {"l": 70, "r": 24, "t": top, "b": 60},
        "showlegend": bool(named and many and legend),
        **panel.get("keep", {}),
        "legend": {"orientation": "h", "y": -0.16, "yanchor": "top", "font": {"size": LEGEND}},
        "shapes": panel["shapes"],
        "annotations": [{**a, "font": {**a.get("font", {}), "size": max(a.get("font", {}).get("size", 0), 11)}} for a in panel["annotations"]],
    }
    if layout["showlegend"]:
        layout["margin"]["b"] = 110
    for name, ax in (("xaxis", panel["x"]), ("yaxis", panel["y"])):
        if ax is None:
            continue
        ax = {**ax, "automargin": True, "tickfont": {"size": TICK}}
        title = ax.get("title")
        if isinstance(title, dict):
            ax["title"] = {**title, "font": {"size": AXIS}}
        layout[name] = ax
    return json.dumps({"data": panel["data"], "layout": layout}, separators=(",", ":")).replace("<", "\\u003c")


class Page:
    def __init__(self) -> None:
        self.n = 0

    def plot(self, panel: dict[str, Any], template: Any, height: int, **opts: Any) -> str:
        self.n += 1
        spec = panel_figure(panel, template, height, **opts)
        title = re.sub(r"<br\s*/?>", " · ", panel["title"] or "")
        head = f'<div class="ptitle">{title}</div>' if title else ""
        return (f'<div class="panel">{head}<div class="pl" style="height:{height}px" data-spec="s{self.n}"></div>'
                f'<script type="application/json" id="s{self.n}">{spec}</script></div>')  # fmt: skip

    def grid(self, panels: list[dict[str, Any]], template: Any, cols: int, height: int, **opts: Any) -> str:
        body = "".join(self.plot(p, template, height, **opts) for p in panels)
        return f'<div class="grid c{cols}">{body}</div>' if body else ""


def links(inst: Any, out: Path) -> str:
    items = [(r.rsplit(".", 1)[1].upper(), r) for r in inst.images]
    if inst.data is not None:
        items.append(("data", inst.data))
    return "".join(f'<a href="file://{out}/plots/{r}">{label}</a>' for label, r in items)


def more(summary: str, body: str, open_: bool = False) -> str:
    if not body:
        return ""
    return f'<details class="more"{" open" if open_ else ""}><summary>{E(summary)}</summary>{body}</details>'


def split(pid: str, panels: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    show = PANELS[pid]["show"]
    if show is None:
        return panels, []
    if show == "first":
        return panels[:1], panels[1:]
    if show == "first_replicon":
        first = (panels[0]["title"] or "").split(" ")[-1]
        keep = [p for p in panels if first and (p["title"] or "").endswith(first)]
        return keep, [p for p in panels if p not in keep]
    keep = [p for p in panels if any(s in (p["title"] or "") for s in show)]
    return keep, [p for p in panels if p not in keep]


def twin_block(pid: str, inst: Any, spec: str, page: Page, out: Path) -> tuple[str, str]:
    conf = PANELS[pid]
    figure = json.loads(spec)
    template = figure["layout"].get("template")
    shown, rest = split(pid, explode(figure))
    cols, height = conf.get("cols", 2), conf["h"]
    opts = {"top": conf.get("top", 12), "legend": conf.get("legend", True)}
    visible = page.grid(shown, template, 1 if len(shown) == 1 else cols, height, **opts)
    hidden = ""
    if rest:
        tall = [p for p in rest if "SHAP" in (p["title"] or "")]
        hidden = page.grid(tall, template, 1, 560) + page.grid([p for p in rest if p not in tall], template, conf.get("more_cols", 2), conf.get("more_h", 340))
    label = f"{len(rest)} more panel{'s' if len(rest) != 1 else ''}"
    # figure-level notes (paper coordinates) are caveats, not panel furniture: keep them as text
    anns = figure["layout"].get("annotations", [])
    texts = [a["text"] for a in anns if a.get("xref") == "paper" and a.get("yref") == "paper" and a.get("yanchor") != "bottom"]
    dropped = [sh for sh in figure["layout"].get("shapes", []) if sh.get("xref") == "paper" and sh.get("yref") == "paper"]
    if dropped:
        print(f"{pid}: {len(dropped)} paper-level shapes not carried over", file=sys.stderr)
    notes = "".join(f'<p class="fnote">{t}</p>' for t in texts)
    tools = f'{notes}<div class="flinks">{pid} {links(inst, out)}</div>'
    return visible + (more(label, hidden) if hidden else "") + tools, inst.stem


def static_block(inst: Any, width: int, out: Path) -> str:
    svg = embed_svg(inst.svg, inst.stem)
    native = re.search(r'<svg[^>]*\swidth="([\d.]+)pt"', svg)
    if native:  # 1.8 x the print size: a single-column figure must not fill the page
        width = min(width, round(float(native.group(1)) * 4 / 3 * 1.8))
    return f'<figure class="static" style="max-width:{width}px">{svg}<div class="flinks">{E(inst.stem)} {links(inst, out)}</div></figure>'


def loci_by_rank(out: Path, instances: list[Any]) -> tuple[list[Any], list[Any]]:
    """Loci in the index's rank order (the current report shows them alphabetically); the index
    instance returned separately."""
    index = [i for i in instances if i.stem in ("V03_locus", "V13_novel_orfs_locus")]
    loci = [i for i in instances if i not in index]
    if index and index[0].data:
        tsv = pd.read_csv(out / "plots" / index[0].data, sep="\t", dtype=str, keep_default_na=False)
        if "locus_id" in tsv:
            rank = {f"{index[0].stem}_{v}": n for n, v in enumerate(tsv["locus_id"])}
            loci.sort(key=lambda i: rank.get(i.stem, len(rank)))
    return loci, index


# --- tables -----------------------------------------------------------------------------------


def table_block(view: RunView, key: str, top: int = 20) -> str:
    name = PATHS[key].rsplit("/", 1)[1]
    try:
        frame = view.table(key)
    except Exception as exc:
        return f'<p class="note">{E(name)}: could not be read: {E(type(exc).__name__)}</p>'
    if frame is None:
        return f'<p class="note">{E(name)}: not produced in this run.</p>'
    sort, asc, keep = TABLES[key]
    if not len(frame):
        return f'<p class="note">{E(name)}: no rows in this run.</p>'
    total = len(frame)
    if sort and sort in frame:
        frame = frame.assign(_s=pd.to_numeric(frame[sort], errors="coerce")).sort_values("_s", ascending=asc, na_position="last", kind="stable").drop(columns="_s")
    frame = frame.head(top)
    cols = [c for c in keep if c in frame] + [c for c in frame if c not in keep]
    head = "".join(f'<th class="{"" if c in keep else "x"}" title="{E(c)}">{E(LABELS.get(c, c.replace("_", " ")))}</th>' for c in cols)
    rows = []
    for rec in frame[cols].to_dict("records"):
        tds = [cell(rec[c]) for c in cols]
        rows.append("<tr>" + "".join(td if c in keep else td.replace("<td", '<td class="x"', 1) for td, c in zip(tds, cols, strict=True)) + "</tr>")
    order = f", sorted by {LABELS.get(sort, sort)}" if sort else ""
    n_x = len(cols) - len([c for c in keep if c in frame])
    toggle = f'<label class="xt"><input type="checkbox" data-allcols> all {len(cols)} columns</label>' if n_x else ""
    return (f'<div class="tbl"><div class="thead"><b>{E(name)}</b> <span>top {len(frame)} of {total} rows{order}</span>'
            f'<input type="search" data-filter="{key}" placeholder="filter rows" aria-label="filter {E(name)}">{toggle}'
            f'<a href="file://{view_root(view)}/{PATHS[key]}">full table</a></div>'
            f'<div class="table-wrap"><table data-table="{key}"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div></div>')  # fmt: skip


_OUT: list[Path] = []


def view_root(_: RunView) -> Path:
    return _OUT[0]


def extras_block(items: list[dict[str, Any]]) -> tuple[str, int]:
    body, n = "", 0
    for x in items:
        if "rows" in x:
            head = "".join(f"<th>{E(str(h))}</th>" for h in x["header"])
            rows = "".join("<tr>" + "".join(cell(v) for v in r) + "</tr>" for r in x["rows"])
            body += f'<p>{E(x["caption"])}:</p><div class="table-wrap"><table data-table="{E(x["key"])}"><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>'
        else:
            body += f'<p class="caveat">{E(x["text"])}</p>'
            n += 1
    return body, n


# --- key points -------------------------------------------------------------------------------


def tile(c: dict[str, Any]) -> str:
    label = (c.get("prefix", "") + " " + c["label"]).strip().rstrip(";")
    value, est = short(c["value"]), c["estimate"]
    empty = c["name"] == "proteins_placed" and c["value"] == "0"
    brief = est.split("; descriptive")[0] if c["name"] == "peptides_placed" else est
    detail = f'<div class="td">{E(short(c["detail"]))}</div>' if c.get("detail") else ""
    return (f'<div class="tile{" dim" if empty else ""}" data-card="{E(c["name"])}" title="{E(est)}"><div class="tv">{E(value)}</div>'
            f'<div class="tl">{E(label)}</div><div class="te">{E(short(brief))}</div>{detail}</div>')  # fmt: skip


def titles_of(files: dict[str, Any]) -> dict[str, str]:
    out = {}
    for pid, pf in files.items():
        for inst in pf.instances:
            if inst.svg and pid not in out:
                t = svg_title(inst.svg)
                if t:
                    out[pid] = t
    return out


def locus_titles(pf: Any) -> dict[str, str]:
    return {i.stem: svg_title(i.svg) or "" for i in pf.instances if i.svg}


def key_points(view: RunView, cfg: GippaConfig, files: dict[str, Any], titles: dict[str, str], summary: list[dict[str, Any]], notes: list[str], doc: Any) -> str:
    by = {c["name"]: c for c in summary}
    pre = view.preflight() or {}
    reps = pre.get("replicons", []) if isinstance(pre, dict) else []
    size = sum(int(r["length"]) for r in reps)
    who = f"{len(reps)} replicon{'s' if len(reps) != 1 else ''}, {size / 1e6:.2f} Mb" if reps else ""
    tier = pre.get("evidence_tier") if isinstance(pre, dict) else None
    meta = " · ".join(x for x in (who, f"evidence: {tier}" if isinstance(tier, str) else "", f"GIPPA {(doc or {}).get('gippa_version', __version__)} (run)") if x)

    finds: list[tuple[str, str, str]] = []
    def find(sec: str, label: str, pid: str) -> None:
        if pid in titles:
            finds.append((sec, label, titles[pid]))
    find("genome-context", "Overall", "V01")
    find("mapping", "Mapping", "V05")
    find("detection", "Detected vs expected", "V07")
    if "V03" in files:
        loci, _ = loci_by_rank(_OUT[0], files["V03"].instances)
        tops = [svg_title(i.svg) for i in loci[:2] if i.svg]
        if tops:
            finds.append(("genome-context", "Top-ranked loci", "; ".join(t for t in tops if t)))
    find("function", "Function", "V14")
    find("determinants", "Determinants", "V15")
    find("detection", "Calibration", "V18")

    care: list[str] = []
    checks = pre.get("checks", []) if isinstance(pre, dict) else []
    # the report must not carry absolute paths (M9-02 behaviour 6): a shipped version redacts properly
    care += [f'Preflight {c.get("level")}: ' + re.sub(r"/(?:[\w.-]+/)+([\w.-]+)", r"\1", str(c.get("message") or c.get("rule_id"))) for c in checks if isinstance(c, dict) and c.get("level") in ("warn", "error")]
    for c in summary:
        if c["name"] in ("alt_calls", "novel_orfs") and c["value"] != "0" and "unestimated" in c["estimate"]:
            care.append(f'{c["value"]} {c["label"]}: error rate unestimated.'.replace("1 novel ORFs", "1 novel ORF"))
    if "V06" in titles and ("No stratum" in titles["V06"] or "unestimated" in titles["V06"]):
        care.append(f'Tolerant mapping: {titles["V06"]}.')
    r2 = by.get("model_r2")
    if r2 and "uncalibrated" in r2["estimate"]:
        care.append("Model R²: the permutation null p is uncalibrated, not a discovery claim (ADR-0034).")
    care += [x["text"] for x in content.extras(view, doc).get("summary", []) if "text" in x]
    care += notes
    code = content.genetic_code_note(cfg.annotation.translation_table)
    if code:
        care.append(code["text"])
    li = "".join(f'<li><a href="#{s}">{E(lab)}</a> {E(short(t))}</li>' for s, lab, t in finds)
    warn = "".join(f"<li>{E(short(t))}</li>" for t in care) or "<li>No warnings from this run.</li>"
    order = ["genes_detected", "peptides_placed", "model_r2", "alt_calls", "novel_orfs", "islands", "operon_flags", "proteins_placed"]
    tiles = "".join(tile(by[n]) for n in order if n in by)
    return (f'<p class="meta">{E(meta)}</p><div class="tiles">{tiles}</div>'
            f'<div class="kp"><div><h3>What this run shows</h3><ul>{li}</ul></div>'
            f'<div class="care"><h3>Read with care</h3><ul>{warn}</ul></div></div>')  # fmt: skip


# --- sections ---------------------------------------------------------------------------------


def section_body(sid: str, plots: tuple[str, ...], tables: tuple[str, ...], view: RunView, files: dict[str, Any], live: dict[str, Any], page: Page, out: Path, titles: dict[str, str], extras: dict[str, Any]) -> tuple[str, str, str]:
    """(visible figures, collapsed figures, count summary)."""
    vis, hid = "", ""
    for pid in plots:
        pf = files.get(pid)
        if pf is None or not pf.instances:
            vis += f'<p class="note">{pid}: not produced in this run.</p>'
            continue
        if pid == "V03":
            loci, index = loci_by_rank(out, pf.instances)
            vis += "".join(static_block(i, 1240, out) for i in loci[:2] if i.svg)
            rest = "".join(static_block(i, 1240, out) for i in [*index, *loci[2:]] if i.svg)
            hid += more(f"{len(loci[2:])} more loci and the ranked index", rest)
            continue
        if pid == "V12":
            hid += more("Peptide-vs-genome alignments at the top ALT sites (V12)", "".join(static_block(i, 1100, out) for i in pf.instances if i.svg))
            continue
        if pid == "V13":
            first, *others = pf.instances
            if first.svg:
                vis += static_block(first, 760, out)
            loci, index = loci_by_rank(out, others)
            hid += more(f"Locus views of {len(loci)} novel ORFs", "".join(static_block(i, 1000, out) for i in [*index, *loci] if i.svg))
            continue
        for inst in pf.instances:
            if inst.stem in live and pid in PANELS:
                block, _ = twin_block(pid, inst, live[inst.stem][1], page, out)
                vis += block
            elif inst.svg:
                vis += static_block(inst, 1000 if pid == "V01" else 1100, out)
    return vis, hid, ""


CSS = """
body{font-family:"DejaVu Sans",Arial,sans-serif;margin:0;color:#222;font-size:15px;line-height:1.45;background:#fff}
main{max-width:1280px;margin:0 auto;padding:0 24px 60px}
nav{position:sticky;top:0;z-index:5;background:#123;color:#fff;padding:6px 24px;display:flex;flex-wrap:wrap;gap:4px 14px;font-size:13px}
nav a{color:#fff;text-decoration:none} nav a:hover{text-decoration:underline} nav b{margin-right:10px}
h1{font-size:1.5em;margin:18px 0 0} section{padding:10px 0 18px;border-top:1px solid #ddd;scroll-margin-top:40px}
#summary{border-top:0} h2{font-size:1.25em;margin:10px 0 4px} h3{font-size:1em;margin:0 0 6px}
.meta{color:#555;margin:2px 0 12px}
.tiles{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}
.tile{border:1px solid #ccd;border-radius:6px;padding:10px 14px;background:#fafbfd}
.tile.dim{opacity:.55}.tv{font-size:1.7em;font-weight:700;line-height:1.2}.tl{font-size:.92em}.te,.td{font-size:.8em;color:#556}
.kp{display:grid;grid-template-columns:3fr 2fr;gap:18px;margin-top:14px}
.kp ul{margin:0;padding-left:1.1em}.kp li{margin-bottom:5px}.kp li a{font-weight:600;margin-right:3px}
.care{border-left:4px solid #8a5a00;background:#fff8e6;padding:10px 14px;border-radius:0 6px 6px 0}
.lead{font-size:1.05em;margin:0 0 10px;font-weight:600;color:#123}
.grid{display:grid;gap:14px 18px;margin:8px 0}.grid.c1{grid-template-columns:1fr}.grid.c2{grid-template-columns:1fr 1fr}.grid.c3{grid-template-columns:1fr 1fr 1fr}
.panel{min-width:0}.ptitle{font-weight:600;font-size:14px;margin-bottom:2px}
.static{margin:10px 0}.static svg{width:100%;height:auto;display:block}
.flinks{font-size:.78em;color:#777;margin-top:2px}.flinks a{margin-left:8px;color:#369}
details.more{margin:10px 0;border:1px solid #dde;border-radius:6px;padding:0 14px}
details.more>summary{cursor:pointer;padding:9px 0;font-weight:600;color:#245}
details.more[open]{padding-bottom:10px}
.tbl{margin:8px 0 16px}.thead{display:flex;gap:12px;align-items:center;flex-wrap:wrap;font-size:.85em;margin-bottom:4px}.thead span{color:#667}
.thead input[type=search]{font-size:.95em;padding:2px 6px}.xt{color:#667}
.table-wrap{overflow-x:auto}table{border-collapse:collapse;font-size:.82em}
th,td{border-bottom:1px solid #e6e6ea;padding:3px 9px;text-align:left;white-space:nowrap}th{background:#f3f4f7;cursor:pointer}
.tbl:not(.allcols) .x{display:none}
.fnote{font-size:.85em;color:#555;margin:.2em 0}.note{color:#666;font-style:italic}.caveat{font-size:.9em;margin:.35em 0}
.prov td,.prov th{white-space:normal;overflow-wrap:anywhere}.prov{width:100%;table-layout:auto}
@media(max-width:900px){.tiles{grid-template-columns:1fr 1fr}.kp,.grid.c2,.grid.c3{grid-template-columns:1fr}}
"""

JS = """
(function(){
var done=new WeakSet();
var config={responsive:true,displayModeBar:false};
function draw(d){if(done.has(d))return;done.add(d);var s=JSON.parse(document.getElementById(d.dataset.spec).textContent);Plotly.newPlot(d,s.data,s.layout,config);}
var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){draw(e.target);io.unobserve(e.target);}});},{rootMargin:'900px'});
document.querySelectorAll('.pl').forEach(function(d){io.observe(d);});
window.__drawAll=function(){document.querySelectorAll('.pl').forEach(function(d){if(d.offsetParent!==null)draw(d);});};
document.querySelectorAll('[data-allcols]').forEach(function(c){c.addEventListener('change',function(){c.closest('.tbl').classList.toggle('allcols',c.checked);});});
})();
"""


def build(out: Path) -> str:
    _OUT[:] = [out]
    cfg = GippaConfig.model_validate(yaml.safe_load((out / "resolved_config.yaml").read_text()))
    load_plot_modules()
    view = RunView.from_outdir(out).with_modules(cfg.modules)
    files = read_plot_files(out / "plots")
    live = twins(out / "plots", files, twin_meta(view, cfg))
    summary, notes = cards(view, cfg)
    doc, note = content.read_manifest(view)
    titles = titles_of(files)
    extras = content.extras(view, doc)
    page = Page()
    tables_js = (Path(content.__file__).parent / "templates" / "tables.js").read_text()
    parts = []
    for s in SECTIONS:
        vis, hid, _ = section_body(s.id, s.plots, s.tables, view, files, live, page, out, titles, extras)
        lead_pid = LEAD.get(s.id)
        lead = f'<p class="lead">{E(short(titles[lead_pid]))}</p>' if lead_pid in titles else ""
        body = ""
        if s.id == "summary":
            demo = content.demo_banner(doc)
            body = (f'<p class="note">{E(demo)}</p>' if demo else "") + key_points(view, cfg, files, titles, summary, notes, doc)
        elif s.id == "variants":
            var = view.table("variants") if view.has("variants") else None
            if var is not None and len(var):
                counts = var.groupby(["call", "confidence"], dropna=False).size().reset_index(name="sites")
                rows = "".join("<tr>" + "".join(f"<td>{E(str(v))}</td>" for v in r) + "</tr>" for r in counts.itertuples(index=False))
                body = f'<table class="small"><thead><tr><th>call</th><th>confidence</th><th>sites</th></tr></thead><tbody>{rows}</tbody></table>'
        tkeys = [k for k in s.tables if k != "sensitivity"] + (["variants"] if s.id == "variants" else [])
        tbl = "".join(table_block(view, k) for k in tkeys)
        cav, n_cav = ("", 0)
        if s.id != "summary":
            cav, n_cav = extras_block(extras.get(s.id, []))
        foot = hid
        if tbl:
            foot += more(f"Tables ({len(tkeys)}): top rows, key columns", tbl, open_=s.id == "novel-orfs")
        if cav:
            foot += more(f"Caveats and assumptions ({n_cav})" if n_cav else "Diagnostics", cav)
        if s.id == "sensitivity":
            whole = view.table("sensitivity") if view.has("sensitivity") else None
            if whole is not None and len(whole):
                j = pd.to_numeric(whole["jaccard"], errors="coerce").dropna()
                if len(j):
                    lead = f'<p class="lead">Flagged sets under the missing-data variants: Jaccard {j.min():.2f} to {j.max():.2f} against the primary run.</p>'
                rows = [r for _, g in content.sensitivity_groups(whole) for r in g]
                head = "".join(f"<th>{E(c)}</th>" for c in whole.columns)
                tr = "".join("<tr>" + "".join(f"<td>{E(str(v))}</td>" for v in r) + "</tr>" for r in rows)
                body = f'<div class="table-wrap"><table data-table="sensitivity"><thead><tr>{head}</tr></thead><tbody>{tr}</tbody></table></div>'
        if s.id == "methods":
            snippets = "".join(f'<h3>{E(t)}</h3>{b}' for t, b in content.snippets())
            lines = f"<p>{E(content.genetic_code_methods(cfg.annotation.translation_table))}</p>"
            lic = content.licence_notes(doc)
            if lic:
                lines += f"<p>Licence notes: {E('; '.join(lic))}</p>"
            foot += more("Methods notes", lines + snippets + cav)
        if s.id == "provenance":
            prov = content.provenance(view, doc, note)
            rows = dict(prov.get("lines", []))
            lead = f'<p class="lead">GIPPA {E(rows.get("GIPPA version", __version__))}, commit {E(rows.get("commit", "not recorded")[:10])}, seed {E(rows.get("seed", "not recorded"))}.</p>'
            tab = lambda head, rs: f'<table class="prov"><thead><tr>{"".join(f"<th>{E(h)}</th>" for h in head)}</tr></thead><tbody>' + "".join("<tr>" + "".join(f"<td>{E(str(v))}</td>" for v in r) + "</tr>" for r in rs) + "</tbody></table>"
            det = tab(["item", "value"], prov.get("lines", [])) + f'<p>Resolved configuration SHA-256: {E(str(prov.get("config")))}</p>' + tab(["file", "role", "SHA-256"], prov["inputs"])
            if "tools" in prov:
                det += tab(["tool", "version", "digest"], prov["tools"]) + tab(["database", "version", "checksum", "retrieved"], prov.get("databases", []))
            foot += more("Full provenance: inputs, tools, databases, digests", det)
        heading = "Summary: key points" if s.id == "summary" else s.heading
        parts.append(f'<section id="{s.id}"><h2>{E(heading)}</h2>{lead}{body}{vis}{foot}</section>')
    nav = "".join(f'<a href="#{s.id}">{NAV[s.id]}</a>' for s in SECTIONS)
    return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>GIPPA report (prototype)</title><style>{CSS}</style>'
            f"<script>{plotly.offline.get_plotlyjs()}</script></head><body><nav><b>GIPPA report</b>{nav}</nav><main>"
            f'<h1>GIPPA report</h1>{"".join(parts)}</main><script>{tables_js}</script><script>{JS}</script></body></html>')  # fmt: skip


if __name__ == "__main__":
    out_dir, target = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(build(out_dir), encoding="utf-8")
    print(target, target.stat().st_size)
