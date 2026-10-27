# ruff: noqa: RUF001  (the report's wording is Unicode)
"""Report content beyond the M9-02 previews (M9-16; GIPPA_SPEC.md §10.3 sections 11-13): the
mandatory caveats and small tables of each section, the sensitivity tables by family, the
provenance block and the methods snippets.

Pure apart from the reads a ``RunView`` does and the package-data snippets: a view in, template
data out. Every sentence below that an ADR or the spec mandates is quoted from it.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from importlib import resources
from pathlib import PurePosixPath
from typing import Any

import pandas as pd
from markupsafe import Markup

import gippa
import gippa.determinants
from gippa.context.diagnostics import render_l_hat
from gippa.data import calibration
from gippa.io.run_view import RunView
from gippa.io.tsv import format_value

#: GIPPA_SPEC.md §7.4, verbatim
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
#: ADR-0061 decision 1 (rows flagged ``tu_source_distance_model``), verbatim
MORAN_DISTANCE = (
    "Moran's I under the TU-block null is anti-conservative at short lags when TUs come from the "
    "distance model; rejection at α = 0.05 was 0.34 in simulation"
)
DISCORDANCE_DISTANCE = (
    "conservative under the distance-model fallback (merged TU boundaries widen the "
    "within-operon fit); fewer discordant pairs are found"
)
#: ADR-0061 decision 4 (rows flagged ``tu_block_tu_error``, every run), verbatim
MORAN_TU_ERROR = (
    "Moran's I under the TU-block null assumes exact TU boundaries and is anti-conservative at "
    "short lags when TUs are split; rejection at α = 0.05 was 0.67 in simulation with UniOP's "
    "K-12 error rates"
)
DISTANCE_MODEL = "distance_model"
NOT_WRITTEN = "run_manifest.json not written by this GIPPA version."
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


#: preflight rules whose ``warn`` opens the summary section (M11-14w; Q3: this one only in 1.0)
SUMMARY_WARNINGS = ("genome.genetic_code",)


def summary_warnings(view: RunView) -> list[dict[str, str]]:
    """The ``rule_id`` and ``message`` of each ``SUMMARY_WARNINGS`` preflight check at level
    ``warn`` in the run's ``preflight.json`` (M11-14w); empty when it is absent or unreadable."""
    try:
        doc = view.preflight()
    except Exception:
        return []
    checks = (doc or {}).get("checks") if isinstance(doc, Mapping) else None
    return [{"rule_id": str(c["rule_id"]), "message": str(c.get("message") or "")}
            for c in _list(checks)
            if isinstance(c, Mapping) and c.get("rule_id") in SUMMARY_WARNINGS
            and c.get("level") == "warn"]  # fmt: skip


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


def _num(value: Any) -> str:
    return f"{value:.3f}" if isinstance(value, int | float) else NOT_RECORDED


def operon_source_sentence(reason: str | None) -> str:
    """The M11-14d summary and Methods sentence of a distance-model run; the operating point
    (14s2) is read from ``calibration.yaml`` at render time."""
    op = calibration.load().get("operon_distance") or {}

    def ci(key: str) -> str:
        return f"{_num(op.get(key))} [{_num(op.get(key + '_ci_low'))}, " \
               f"{_num(op.get(key + '_ci_high'))}]"  # fmt: skip

    return (
        f"Operons/TUs: intergenic-distance model ({reason or NOT_RECORDED}); "
        f"precision {ci('precision_at_pair_threshold')}, recall {ci('recall_at_pair_threshold')}"
        " on RegulonDB K-12 at the pair threshold; TU-based tests use these units."
    )


def _tu_block_morans(view: RunView) -> bool:
    if not view.has("autocorrelation"):
        return False
    ac = _safe(view, "autocorrelation")
    if ac is None or not {"statistic", "null_model"} <= set(ac.columns):
        return False
    return bool(((ac["statistic"] == "morans_i") & (ac["null_model"] == "tu_block")).any())


def extras(view: RunView, doc: Mapping[str, Any] | None) -> dict[str, list[dict[str, Any]]]:
    """Section id -> its caveat paragraphs (``text``) and small tables (``key``, ``header``,
    ``rows``, ``caption``), each present when its condition holds (M9-16 behaviour 3, 4), plus
    the operon-source sentence and the ADR-0061 caveats (M11-14d)."""
    out: dict[str, list[dict[str, Any]]] = {}
    distance = (doc or {}).get("operon_source") == DISTANCE_MODEL
    if distance:
        sentence = _p(operon_source_sentence((doc or {}).get("operon_source_reason")))
        out["summary"] = [sentence]
        out["methods"] = [sentence]
    if view.has("pairs") or view.has("operon_flags"):
        out["operons"] = [_p(SIGMA70), *_missed(view)]
        if distance and view.has("pairs"):
            text = "Discordance tests (both_detected_t, flagged tu_source_distance_model) are"
            out["operons"].append(_p(f"{text} {DISCORDANCE_DISTANCE}."))
    if view.has("determinants_model"):
        out["determinants"] = [_p(t) for t in (*gippa.determinants.CAVEATS, *DETERMINANTS)]
    out["function"] = _function(view, doc)
    ctx = [_p(f"{ISLANDS} See the islands methods note.")] if view.has("islands") else []
    if _tu_block_morans(view):
        ctx.append(_p(f"{MORAN_TU_ERROR} (rows flagged tu_block_tu_error)."))
        if distance:
            ctx.append(_p(f"{MORAN_DISTANCE} (rows flagged tu_source_distance_model)."))
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


#: a path attached to a short option (``-o/x``)
_SHORT_PATH = re.compile(r"-[A-Za-z]/")


def _local_path(arg: str, outdir: str) -> str:
    """``arg`` without an absolute local path (M9-02 behaviour 6, §6.6; #323 B1): a path under
    the run's ``outdir`` as ``OUT/<relative>``, any other absolute path as its file name. The
    path may be the whole token, follow its last ``=``, be an item of a comma list, or be
    attached to a short option (``-o/x``); relative parts are unchanged."""
    if arg.startswith("/"):
        head, value = "", arg
    elif _SHORT_PATH.match(arg):
        head, value = arg[:2], arg[2:]
    else:
        cut = arg.rfind("=") + 1
        head, value = arg[:cut], arg[cut:]
    return head + ",".join(_shown_path(item, outdir) for item in value.split(","))


def _shown_path(value: str, outdir: str) -> str:
    """One path item: ``OUT/<relative>`` under ``outdir``, else the file name if absolute."""
    if not value.startswith("/"):
        return value
    path, out = PurePosixPath(value), PurePosixPath(outdir) if outdir.startswith("/") else None
    if out is not None and (path == out or out in path.parents):
        return str(PurePosixPath("OUT", path.relative_to(out)))
    return path.name or value


def _command_line(doc: Mapping[str, Any]) -> str:
    """The manifest's ``command_line`` joined by single spaces, absolute paths made local."""
    outdir = str(doc.get("outdir") or "")
    return " ".join(_local_path(str(a), outdir) for a in doc.get("command_line") or [])


def licence_notes(doc: Mapping[str, Any] | None) -> list[str]:
    """The run manifest's ``licence_notes`` (M11-15c2), for the Methods section."""
    return [str(n) for n in _list((doc or {}).get("licence_notes")) if isinstance(n, str)]


def provenance(view: RunView, doc: Mapping[str, Any] | None, note: str | None) -> dict[str, Any]:
    """Provenance (M9-16 behaviour 6; lead decisions 1, 7, 8): from ``run_manifest.json`` when it
    is readable, else the pass-1 subset (preflight SHA-256s) with a note."""
    rows = doc.get("inputs") if doc is not None else (view.preflight() or {}).get("inputs")
    inputs = [[str(i.get("file", "")).replace("\\", "/").rsplit("/", 1)[-1],
               str(i.get("role", "")), str(i.get("sha256", ""))]
              for i in _list(rows) if isinstance(i, Mapping)]  # fmt: skip
    base = {
        "version": gippa.__version__,
        "inputs": inputs,
        "config": view.digest("resolved_config"),
    }
    if doc is None:
        return {**base, "note": note or NOT_WRITTEN}
    seeds = doc.get("seeds")
    dbs = [[str(d.get("name", "")), str(d.get("version", "")), str(d.get("checksum", "")),
            str(d.get("retrieved", ""))] for d in _databases(doc)]  # fmt: skip
    go = [str(d.get("version")) for d in _databases(doc) if d.get("name") == "go_basic"]
    dirty = " (uncommitted changes)" if doc.get("git_dirty") else ""
    lines = [
        ("GIPPA version", str(doc.get("gippa_version", NOT_RECORDED))),
        ("commit", f"{doc.get('git_commit') or NOT_RECORDED}{dirty}"),
        ("command line", _command_line(doc) or NOT_RECORDED),
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
    files = resources.files("gippa.report").joinpath("snippets")
    found = [f for f in files.iterdir() if f.name.endswith(".html")]
    return [(f.name[: -len(".html")], Markup(f.read_text(encoding="utf-8")))
            for f in sorted(found, key=lambda f: f.name)]  # fmt: skip
