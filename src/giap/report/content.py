# ruff: noqa: RUF001  (the report's wording is Unicode)
"""Report content beyond the M9-02 previews (M9-16; GIAP_SPEC.md §10.3 sections 11-13): the
mandatory caveats and small tables of each section, the sensitivity tables by family, the
provenance block and the methods snippets.

Pure apart from the reads a ``RunView`` does and the package-data snippets: a view in, template
data out. Every sentence below that an ADR or the spec mandates is quoted from it.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from importlib import resources
from typing import Any

import pandas as pd
from markupsafe import Markup

import giap
import giap.determinants
from giap.context.diagnostics import render_l_hat
from giap.io.run_view import RunView
from giap.io.tsv import format_value

#: GIAP_SPEC.md §7.4, verbatim
SIGMA70 = (
    "Promotech and Promoter Calculator see σ70 promoters only, and TerminatorNet and TransTermHP "
    "see intrinsic terminators only. A missing call therefore never proves a missing element."
)
DETERMINANTS = (
    "The univariate Spearman ρ and AUROC CIs (TU-block bootstrap) are uncalibrated (ADR-0032).",
    "The elastic-net ΔR² is the primary block attribution; LightGBM ΔR² and SHAP shares are "
    "secondary, nonlinear checks (ADR-0033).",
    "The permutation null p is uncalibrated: it is valid only for exchangeable TUs and is not a "
    "discovery claim (ADR-0034).",
)
GSEA = (
    "GSEA results are conditional on detection, with no p_detect adjustment (ADR-0026).",
    "NES is shown without an interval (ADR-0027).",
)
GSEA_CLUSTER = (
    "P-values are valid under TU exchangeability, and there is no level claim for sets flagged "
    "genomic_cluster (ADR-0025)."
)
#: docs/methods/islands.md, "Validity domain (ADR-0021 Amendment 1)"
ISLANDS = "The FWER is claimed for a dependence regime, not for a replicon size."
KEGG_TERMS = "KEGG data: academic use only (rest.kegg.jp terms)."
NOT_WRITTEN = "run_manifest.json not written by this GIAP version."
NOT_RECORDED = "not recorded"
_DIAG_HEADER = ["replicon", "ℓ̂ (l_hat)", "ℓ′ (block_len)", "Î₁ (lag-1 Moran's I)", "block rule"]
_BIN_COLUMNS = ("bin", "lower_bp", "upper_bp", "N_NO", "N_cv", "N_sw")


def _p(text: str) -> dict[str, Any]:
    return {"text": text}


def _bound(v: Any) -> str:
    if isinstance(v, int | float) and math.isinf(v):
        return "∞" if v > 0 else "−∞"
    return format_value(v)


def read_manifest(view: RunView) -> tuple[Mapping[str, Any] | None, str | None]:
    """The run manifest and, when it exists but cannot be read, a note naming it."""
    try:
        return view.manifest(), None
    except Exception as exc:
        return None, f"run_manifest.json could not be read: {type(exc).__name__}"


def _list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _databases(doc: Mapping[str, Any] | None) -> list[Mapping[str, Any]]:
    return [d for d in _list((doc or {}).get("databases")) if isinstance(d, Mapping)]


def kegg_date(doc: Mapping[str, Any] | None) -> str | None:
    """The ``retrieved`` dates of the manifest's ``kegg*`` databases (lead decisions 2, 8)."""
    kegg = [d for d in _databases(doc) if str(d.get("name", "")).startswith("kegg")]
    dates = sorted({str(d["retrieved"]) for d in kegg if d.get("retrieved")})
    return ", ".join(dates) or None


def _safe(view: RunView, key: str) -> pd.DataFrame | None:
    try:
        return view.table(key)
    except Exception:  # the table preview of the section already carries the note
        return None


def _function(view: RunView, doc: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    flags = view.run_flags()
    gsea = _safe(view, "enrichment_gsea") if view.has("enrichment_gsea") else None
    if view.has("enrichment_gsea"):
        out += [_p(t) for t in GSEA]
        if gsea is not None and gsea["flags"].astype(str).str.contains("genomic_cluster").any():
            out.append(_p(GSEA_CLUSTER))
    cats = _safe(view, "category_summary") if view.has("category_summary") else None
    if not view.has("category_summary"):
        return out
    cov = "" if cats is None else ", ".join(
        f"{c} {format_value(v)}"
        for c, v in cats.drop_duplicates("collection")[["collection", "coverage"]].itertuples(False)
    )  # fmt: skip
    out.append(_p(f"Annotation coverage leans toward well-studied genes (§7.7.1); coverage per "
                  f"collection: {cov or NOT_RECORDED}."))  # fmt: skip
    if "bakta_db_light" in flags:
        out.append(_p("Annotation used the Bakta light DB, which assigns fewer functions than "
                      "the full DB."))  # fmt: skip
    if "layer_missing:kegg" in flags:
        out.append(_p("KEGG not installed."))
    else:
        date = kegg_date(doc) or "not recorded in this run"
        out.append(_p(f"{KEGG_TERMS} KEGG retrieval date: {date}."))
    return out


def _diagnostics(view: RunView) -> list[dict[str, Any]]:
    if not view.has("island_diagnostics"):
        if not view.has("islands"):
            return []
        return [_p("island_diagnostics.parquet: not produced in this run, so ℓ̂, ℓ′ and Î₁ per "
                   "replicon are not shown.")]  # fmt: skip
    diag = _safe(view, "island_diagnostics")
    if diag is None:
        return [_p("island_diagnostics.parquet: could not be read.")]
    rows = [[str(r["contig"]), render_l_hat(r["l_hat"], str(r["block_rule"])),
             format_value(r["block_len"]), format_value(r["morans_i_lag1"]), str(r["block_rule"])]
            for r in diag.to_dict("records")]  # fmt: skip
    return [{"key": "island_diagnostics", "header": _DIAG_HEADER, "rows": rows,
             "caption": "Island diagnostics per replicon (ADR-0021 Amendment 1)"}]  # fmt: skip


def _missed(view: RunView) -> list[dict[str, Any]]:
    try:
        found = view.events("operons.missed_operons")
    except Exception:
        return [_p("operons.events.jsonl: could not be read.")]
    if not found:
        return []
    data = found[-1]
    bins = [b for b in data.get("bin_counts", []) if isinstance(b, Mapping)]
    rows = [[(_bound if c.endswith("_bp") else format_value)(b.get(c)) for c in _BIN_COLUMNS]
            for b in bins]  # fmt: skip
    caption = f"Missed-operon decoy counts per distance bin (π₀ = {format_value(data.get('pi0'))})"
    return [{"key": "missed_operons_bins", "header": list(_BIN_COLUMNS), "rows": rows,
             "caption": caption}]  # fmt: skip


def extras(view: RunView, doc: Mapping[str, Any] | None) -> dict[str, list[dict[str, Any]]]:
    """Section id -> its caveat paragraphs (``text``) and small tables (``key``, ``header``,
    ``rows``, ``caption``), each present when its condition holds (M9-16 behaviour 3, 4)."""
    out: dict[str, list[dict[str, Any]]] = {}
    if view.has("pairs") or view.has("operon_flags"):
        out["operons"] = [_p(SIGMA70), *_missed(view)]
    if view.has("determinants_model"):
        out["determinants"] = [_p(t) for t in (*giap.determinants.CAVEATS, *DETERMINANTS)]
    out["function"] = _function(view, doc)
    ctx = [_p(f"{ISLANDS} See the islands methods note.")] if view.has("islands") else []
    out["genome-context"] = ctx + _diagnostics(view)
    return out


def sensitivity_groups(frame: pd.DataFrame) -> list[tuple[str, list[list[str]]]]:
    """``sensitivity.tsv`` rows by ``family`` in first-seen order, cells formatted, ``jaccard``
    as ``%.2f``, ``NA`` or "no flagged set in either run" (M9-16 behaviour 5)."""
    groups: dict[str, list[list[str]]] = {}
    for rec in frame.to_dict("records"):
        cells = []
        for col in frame.columns:
            v = rec[col]
            if col == "jaccard":
                both0 = rec["n_primary"] == 0 and rec["n_variant"] == 0
                missing = v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))
                v = "no flagged set in either run" if both0 else "NA" if missing else f"{v:.2f}"
            cells.append(format_value(v))
        groups.setdefault(str(rec["family"]), []).append(cells)
    return list(groups.items())


def _tools(doc: Mapping[str, Any]) -> list[list[str]]:
    tools = [t for t in _list(doc.get("tools")) if isinstance(t, Mapping)]
    return [[str(t.get("name", "")), str(t.get("version") or NOT_RECORDED),
             str(t.get("container_digest") or t.get("conda_lock_hash") or NOT_RECORDED)]
            for t in tools]  # fmt: skip


def provenance(view: RunView, doc: Mapping[str, Any] | None, note: str | None) -> dict[str, Any]:
    """Provenance (M9-16 behaviour 6; lead decisions 1, 7, 8): from ``run_manifest.json`` when it
    is readable, else the pass-1 subset (preflight SHA-256s) with a note."""
    rows = doc.get("inputs") if doc is not None else (view.preflight() or {}).get("inputs")
    inputs = [[str(i.get("file", "")).replace("\\", "/").rsplit("/", 1)[-1],
               str(i.get("role", "")), str(i.get("sha256", ""))]
              for i in _list(rows) if isinstance(i, Mapping)]  # fmt: skip
    base = {"version": giap.__version__, "inputs": inputs, "config": view.digest("resolved_config")}
    if doc is None:
        return {**base, "note": note or NOT_WRITTEN}
    seeds = doc.get("seeds")
    dbs = [[str(d.get("name", "")), str(d.get("version", "")), str(d.get("checksum", "")),
            str(d.get("retrieved", ""))] for d in _databases(doc)]  # fmt: skip
    go = [str(d.get("version")) for d in _databases(doc) if d.get("name") == "go_basic"]
    dirty = " (uncommitted changes)" if doc.get("git_dirty") else ""
    lines = [
        ("GIAP version", str(doc.get("giap_version", NOT_RECORDED))),
        ("commit", f"{doc.get('git_commit') or NOT_RECORDED}{dirty}"),
        ("command line", " ".join(str(a) for a in doc.get("command_line") or []) or NOT_RECORDED),
        ("seed", str(seeds.get("master")) if isinstance(seeds, Mapping) else NOT_RECORDED),
        ("KEGG retrieval date", kegg_date(doc) or "not recorded in this run"),
        ("GO release", ", ".join(go) or NOT_RECORDED),
        ("manifest digest (SHA-256 of run_manifest.json)", view.digest("run_manifest") or ""),
    ]
    if doc.get("provenance_digest"):
        lines.append(("provenance digest", str(doc["provenance_digest"])))
    return {**base, "lines": lines, "tools": _tools(doc), "databases": dbs}


def snippets() -> list[tuple[str, Markup]]:
    """``(topic, html)`` of every methods snippet (package data, generated from
    ``docs/methods`` by scripts/dev/sync_report_snippets.py), by topic."""
    files = resources.files("giap.report").joinpath("snippets")
    found = [f for f in files.iterdir() if f.name.endswith(".html")]
    return [(f.name[: -len(".html")], Markup(f.read_text(encoding="utf-8")))
            for f in sorted(found, key=lambda f: f.name)]  # fmt: skip
