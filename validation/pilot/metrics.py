"""Pilot metrics for the arms A1 and A2 (GIAP_SPEC.md §15.5 metrics 1-5 and 7, §11.6 TE1 on real
data, §15.6 checklist, §15.7 hypotheses 1 and 3); spec docs/specs/M10-05-pilot-metrics.md.

Pure functions (frames and dicts in, dicts out) plus a thin CLI that reads and writes::

    python validation/pilot/metrics.py compute --a1 OUT_A1 --a2 OUT_A2 --cmp CMP \\
        --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --convert-a1 CONV_A1 --convert-a2 CONV_A2 \\
        [--summary-a1 SUM_A1 --summary-a2 SUM_A2] --bakta-json BAKTA_JSON --genome GENOME \\
        --k12-fasta K12_FASTA [--sacct SACCT] [--a3 OUT_A3] --out METRICS_DIR
    python validation/pilot/metrics.py render --metrics METRICS_DIR/metrics.json \\
        [--notes METRICS_DIR/notes.md] --out METRICS_DIR/REPORT.md

Exit codes: 0 written; 1 an internal inconsistency (metric 1 vs ``convert.json``, a render
placeholder missing, a rejected notes file); 2 a precondition. The outputs hold DIA-NN-derived
data and are never written inside the checkout (Khoa 2026-10-23, Q-K1 default).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import re
import subprocess
import sys
from collections.abc import Collection, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from validation.proxies.te_proxies import (
    INFRAME,
    MAX_EXCEPTIONS,
    non_atg_nterm,
    read_coding_genes,
    te1_result,
)

from giap.config import GiapConfig, resolve_config
from giap.mapping.digest import digest
from giap.model.genes import Gene
from giap.stats.cutoffs import wilson_ci
from giap.steps.mapping import select_queries

REPO = Path(__file__).resolve().parents[2]
SPEC = "docs/specs/M10-05-pilot-metrics.md"
SCHEMA = "giap-pilot-metrics/1"
Q_ACCEPT = 0.01  # accepted: q_value <= 0.01 (inclusive), the converter's written value
TE1_MIN = 0.995  # §11.6, point estimate (R3: not tuned)
H1_ALPHA, H1_MOSTLY = 0.05, 0.5
H3_RHO_MIN, H3_J_MIN, H3_UNION_MIN = 0.7, 0.5, 5
UNMAPPED_REASONS = ("no_hit", "rejected_placement_fdr", "below_identity", "too_short")
H3_FAMILIES = ("islands", "operon_flags", "enriched_categories")
OTHER_FAMILIES = ("novel_orfs", "alt_sites")
CODING = ("cds", "sorf")
SOURCE_SEP = " \u203a "


class Precondition(Exception):
    """Exit 2: a missing input, a failed run, a hash mismatch, ``--out`` inside the checkout."""


class Inconsistent(Exception):
    """Exit 1: two computations of the same number disagree."""


# ---------------------------------------------------------------- small pure helpers


def clean(v: Any) -> Any:
    """A JSON-ready scalar: None for NA, Python int/float/bool for NumPy ones."""
    if v is None or v is pd.NA or v is pd.NaT:
        return None
    if isinstance(v, bool | np.bool_):
        return bool(v)
    if isinstance(v, np.integer | int):
        return int(v)
    if isinstance(v, np.floating | float):
        return None if math.isnan(float(v)) else float(v)
    return str(v)


def share(k: int, n: int) -> dict[str, Any]:
    """k / n with its 95 % Wilson interval; ``value`` None and a reason when n = 0."""
    if n == 0:
        return {
            "value": None,
            "ci_low": None,
            "ci_high": None,
            "k": 0,
            "n": 0,
            "reason": "empty_denominator",
        }
    lo, hi = wilson_ci(k, n)
    return {"value": k / n, "ci_low": lo, "ci_high": hi, "k": k, "n": n}


def items(ids: Any) -> list[str]:
    return [] if clean(ids) is None else [x for x in str(ids).split(";") if x]


def contaminant_only(ids: Any, prefixes: Sequence[str]) -> bool:
    """Every ``;`` item is contaminant-prefixed (the preflight rule); no item: False."""
    its = items(ids)
    return bool(its) and all(i.startswith(tuple(prefixes)) for i in its)


def _target(P: pd.DataFrame) -> np.ndarray:
    if "is_decoy" not in P.columns:
        return np.ones(len(P), dtype=bool)
    out: np.ndarray = ~P["is_decoy"].fillna(False).astype(bool).to_numpy()
    return out


def _q(P: pd.DataFrame) -> np.ndarray:
    if "q_value" not in P.columns:
        return np.full(len(P), np.nan)
    out: np.ndarray = pd.to_numeric(P["q_value"], errors="coerce").to_numpy(dtype=float)
    return out


def _col(P: pd.DataFrame, name: str, default: Any) -> list[Any]:
    return list(P[name]) if name in P.columns else [default] * len(P)


# ---------------------------------------------------------------- metric 1


def _class_counts(keys: Mapping[Any, bool]) -> dict[str, int]:
    n_cont = sum(keys.values())
    return {"all": len(keys), "contaminant_only": n_cont, "non_contaminant": len(keys) - n_cont}


def _ident(rows: pd.DataFrame, prefixes: Sequence[str]) -> dict[str, dict[str, int]]:
    prec: dict[Any, bool] = {}
    peps: dict[str, bool] = {}
    groups: dict[tuple[str, ...], set[str]] = {}
    for pep, mod, z, ids in zip(
        rows["peptide"],
        _col(rows, "modified_peptide", None),
        _col(rows, "charge", None),
        rows["protein_ids"],
        strict=True,
    ):
        cont = contaminant_only(ids, prefixes)
        key = (str(pep) if clean(mod) is None else str(mod), clean(z))
        prec[key] = prec.get(key, True) and cont
        peps[str(pep)] = peps.get(str(pep), True) and cont
        groups.setdefault(tuple(sorted(set(items(ids)))), set()).add(str(pep))
    gcls = {g: bool(g) and all(i.startswith(tuple(prefixes)) for i in g) for g in groups}
    min2 = {g: c for g, c in gcls.items() if len(groups[g]) >= 2}
    return {
        "precursors": _class_counts(prec),
        "peptides": _class_counts(peps),
        "protein_groups": _class_counts(gcls),
        "protein_groups_min2": _class_counts(min2),
    }


def identifications(P: pd.DataFrame, prefixes: Sequence[str]) -> dict[str, Any]:
    """Metric 1: distinct precursors, peptides, protein groups (and groups with >= 2 accepted
    peptides) among accepted target rows, by contaminant class; union and per run."""
    acc = P.loc[_target(P) & (_q(P) <= Q_ACCEPT)]
    runs = sorted({str(r) for r in _col(acc, "run", None) if clean(r) is not None})
    per_run = {r: _ident(acc.loc[acc["run"].astype(str) == r], prefixes) for r in runs}
    return {"union": _ident(acc, prefixes), "per_run": per_run}


# ---------------------------------------------------------------- metric 2 and TE1 universes


def universe(P: pd.DataFrame, cfg: GiapConfig) -> tuple[list[str], dict[str, list[int]]]:
    """**U** = ``select_queries(P, cfg)`` and, per peptide of U, the indices of its target rows
    that ``select_queries`` keeps (q <= q_value_max or NA)."""
    U = [str(p) for p in select_queries(P, cfg)["peptide"]]
    keep = _target(P) & ~(_q(P) > cfg.input.q_value_max)
    wanted = set(U)
    rows: dict[str, list[int]] = {}
    for i, (pep, k) in enumerate(zip(P["peptide"].astype(str), keep, strict=True)):
        if k and pep in wanted:
            rows.setdefault(pep, []).append(i)
    return U, rows


def split_contaminants(
    P: pd.DataFrame, cfg: GiapConfig
) -> tuple[list[str], list[str], dict[str, list[int]]]:
    """(U, D, rows): D = U minus the peptides whose kept rows are all contaminant-only."""
    U, rows = universe(P, cfg)
    ids = list(P["protein_ids"]) if "protein_ids" in P.columns else [None] * len(P)
    pre = cfg.input.contaminant_prefixes
    D = [p for p in U if not all(contaminant_only(ids[i], pre) for i in rows.get(p, []))]
    return U, D, rows


def _routes(hits: pd.DataFrame) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for pep, route in zip(hits["peptide"].astype(str), hits["route"].astype(str), strict=True):
        out.setdefault(pep, set()).add(route)
    return out


def _reasons(unmapped: pd.DataFrame) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for seq, r in zip(unmapped["sequence"].astype(str), unmapped["reason"], strict=True):
        out.setdefault(seq, set()).add(str(r))
    return out


def unplaced_reason(reasons: Collection[str] | None) -> str:
    """The first of ``UNMAPPED_REASONS`` listed; ``not_listed`` when absent; else ``other``."""
    if not reasons:
        return "not_listed"
    return next((r for r in UNMAPPED_REASONS if r in reasons), "other")


def mapping_shares(
    P: pd.DataFrame, peptide_hits: pd.DataFrame, unmapped: pd.DataFrame, cfg: GiapConfig
) -> dict[str, Any]:
    """Metric 2: exact / tolerant / unmapped shares of **D** (U minus contaminant-only)."""
    U, D, _ = split_contaminants(P, cfg)
    routes, reasons = _routes(peptide_hits), _reasons(unmapped)
    counts = {"exact": 0, "tolerant": 0, "unmapped": 0}
    why = dict.fromkeys((*UNMAPPED_REASONS, "not_listed"), 0)
    for p in D:
        r = routes.get(p, set())
        kind = "exact" if "exact" in r else "tolerant" if "tolerant" in r else "unmapped"
        counts[kind] += 1
        if kind == "unmapped":
            key = unplaced_reason(reasons.get(p))
            why[key] = why.get(key, 0) + 1
    cont = unmapped[
        (unmapped["label"].astype(str) == "contaminant")
        & (unmapped["kind"].astype(str) == "peptide")
    ]
    n_c = len(U) - len(D)
    return {
        "n_universe": len(U),
        "n_contaminant_only": n_c,
        "n_denominator": len(D),
        "counts": counts,
        "shares": {k: share(v, len(D)) for k, v in counts.items()},
        "unmapped_reasons": why,
        "contaminant_difference": n_c - cont["sequence"].astype(str).nunique(),
    }


def genome_inconsistent(detection_summary: pd.DataFrame) -> dict[str, Any]:
    """The ``qc`` / ``genome_inconsistent_rate`` row: value, n and ``round(value * n)``."""
    ds = detection_summary
    row = ds[
        (ds["section"].astype(str) == "qc")
        & (ds["metric"].astype(str) == "genome_inconsistent_rate")
    ]
    if not len(row):
        return {"value": None, "n": None, "count": None, "reason": "missing_row"}
    value, n = clean(row["value"].iloc[0]), clean(row["n"].iloc[0])
    count = None if value is None or n is None else round(value * n)
    return {"value": value, "n": n, "count": count}


# ---------------------------------------------------------------- metric 3


def _digest_il(protein: str) -> set[str]:
    return {p.replace("I", "L") for _, _, p in digest(protein)}


def strain_specific_set(
    nissle_proteins: Mapping[str, str], k12_proteins: Mapping[str, str]
) -> dict[str, str]:
    """Locus tag -> ``strain_specific`` (no in-silico tryptic peptide in the K-12 digest),
    ``k12_shared`` or ``no_peptides``; I and L collapsed on both sides."""
    k12: set[str] = set()
    for seq in k12_proteins.values():
        k12 |= _digest_il(seq)
    out: dict[str, str] = {}
    for tag, seq in nissle_proteins.items():
        peps = _digest_il(seq)
        out[tag] = "no_peptides" if not peps else "k12_shared" if peps & k12 else "strain_specific"
    return out


def _region_types(v: Any) -> set[str]:
    return {i.split("_", 1)[0] for i in items(v)}


def _cov(sub: pd.DataFrame) -> dict[str, Any]:
    st = sub["detection_status"].astype(str).value_counts()
    det = int(st.get("detected", 0))
    return {
        "genes": len(sub),
        "detected": det,
        "share": share(det, len(sub)),
        "by_status": {str(k): int(v) for k, v in st.items()},
    }


def coverage(genes: pd.DataFrame, strain_specific: Mapping[str, str]) -> dict[str, Any]:
    """Metric 3: detected share of **G** (cds, sorf) per subset; pseudogenes on their own."""
    G = genes[genes["feature_type"].astype(str).isin(CODING)]
    rep = G["replicon_type"].astype(str)
    regions = [_region_types(v) for v in G["region_ids"]]
    cls = [strain_specific.get(str(t), "no_peptides") for t in G["locus_tag"]]
    subs: dict[str, pd.Series] = {"all": pd.Series(True, index=G.index)}
    subs["chromosome"] = rep == "chromosome"
    subs["plasmid"] = rep == "plasmid"
    for contig in sorted(set(G.loc[subs["plasmid"], "contig"].astype(str))):
        subs[f"plasmid:{contig}"] = subs["plasmid"] & (G["contig"].astype(str) == contig)
    subs["unknown"] = ~rep.isin(["chromosome", "plasmid"])
    for kind in ("prophage", "bgc"):
        subs[kind] = pd.Series([kind in r for r in regions], index=G.index)
    for c in ("strain_specific", "k12_shared", "no_peptides"):
        subs[c] = pd.Series([x == c for x in cls], index=G.index)
    out = {name: {"subset": name, **_cov(G[mask])} for name, mask in subs.items()}
    out["pseudogene"] = {
        "subset": "pseudogene",
        **_cov(genes[genes["feature_type"].astype(str) == "pseudogene"]),
    }
    return {"subsets": out, "order": list(out)}


# ---------------------------------------------------------------- metrics 4, 5


ABUNDANCE_FIELDS = (
    "n_shared",
    "n_units",
    "spearman_rho",
    "ci_low",
    "ci_high",
    "method",
    "n_boot",
    "ba_bias_log2",
    "ba_loa_low",
    "ba_loa_high",
    "flags",
)
STABILITY_FIELDS = ("n_a", "n_b", "n_shared", "jaccard", "flags")


def _pair(frame: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    a, b = frame["run_a"].astype(str), frame["run_b"].astype(str)
    return frame[[{x, y} == set(labels) for x, y in zip(a, b, strict=True)]]


def abundance_agreement(compare_abundance: pd.DataFrame, labels: Sequence[str]) -> dict[str, Any]:
    """Metric 4: the A1/A2 ``compare_abundance.tsv`` row, copied unchanged."""
    row = _pair(compare_abundance, labels)
    if not len(row):
        return {**dict.fromkeys(ABUNDANCE_FIELDS), "flags": "missing_row", "reason": "missing_row"}
    return {k: clean(row[k].iloc[0]) if k in row.columns else None for k in ABUNDANCE_FIELDS}


def _family_note(family: str, flags: Any) -> str | None:
    f = items(flags)
    if "both_empty" in f:
        return "none in either arm"
    miss = [x.split(":", 1)[1] for x in f if x.startswith("missing_in:")]
    if miss and family == "enriched_categories":
        return "function module off in " + ", ".join(miss)
    return "missing in " + ", ".join(miss) if miss else None


def stability(compare_stability: pd.DataFrame, labels: Sequence[str]) -> dict[str, Any]:
    """Metric 5: the A1/A2 ``compare_stability.tsv`` rows per family, copied unchanged (plus a
    plain-words ``note`` for the flags); ``novel_orfs`` and ``alt_sites`` under
    ``other_families``."""
    rows = _pair(compare_stability, labels)
    fam = rows["family"].astype(str)

    def one(family: str) -> dict[str, Any]:
        r = rows[fam == family]
        if not len(r):
            return {**dict.fromkeys(STABILITY_FIELDS), "flags": "missing_row", "note": "no row"}
        rec = {k: clean(r[k].iloc[0]) for k in STABILITY_FIELDS}
        return {**rec, "note": _family_note(family, rec["flags"])}

    out: dict[str, Any] = {f: one(f) for f in H3_FAMILIES}
    out["other_families"] = {f: one(f) for f in OTHER_FAMILIES}
    return out


# ---------------------------------------------------------------- metric 7


def _seconds(started: str, finished: str) -> float:
    return (
        dt.datetime.fromisoformat(finished) - dt.datetime.fromisoformat(started)
    ).total_seconds()


STEP_FIELDS = (
    "name",
    "status",
    "n_tasks",
    "n_retries",
    "wall_time_s",
    "cpu_time_s",
    "peak_rss_bytes",
)
SACCT_FIELDS = ("JobID", "JobName", "State", "Elapsed", "MaxRSS", "NodeList", "AllocCPUS")


def runtime(
    manifest: Mapping[str, Any],
    trace_rows: Sequence[Mapping[str, str]],
    sacct_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    """Metric 7: per-step timings of ``run_manifest.json``, the arm's wall time, sums, the max
    step peak RSS (not the arm's total memory), cached steps, trace hosts and the head job's
    ``sacct`` row (the first row whose JobID is not a job step)."""
    steps = [{k: s.get(k) for k in STEP_FIELDS} for s in manifest.get("steps", [])]
    walls = [s["wall_time_s"] for s in steps if s["wall_time_s"] is not None]
    rss = [s["peak_rss_bytes"] for s in steps if s["peak_rss_bytes"] is not None]
    cached = [str(s["name"]) for s in steps if s["status"] == "cached"]
    hosts = sorted({str(r.get("hostname")) for r in trace_rows} - {"", "-", "None"})
    head = [r for r in sacct_rows if "." not in str(r.get("JobID", ""))] or list(sacct_rows)
    started, finished = manifest.get("started_at"), manifest.get("finished_at")
    return {
        "steps": steps,
        "wall_time_s": _seconds(started, finished) if started and finished else None,
        "sum_task_wall_s": float(sum(walls)),
        "max_step_peak_rss_bytes": max(rss) if rss else None,
        "cached_steps": cached,
        "runtime_complete": not cached,
        "totals_lower_bound": bool(cached),
        "hosts": hosts,
        "sacct": {k: head[0].get(k) for k in SACCT_FIELDS} if head else None,
    }


# ---------------------------------------------------------------- TE1 on real data


def _te1_sources(
    P: pd.DataFrame, D: Sequence[str], rows: Mapping[str, list[int]], prefixes: Sequence[str]
) -> dict[str, set[str]]:
    ids = list(P["protein_ids"]) if "protein_ids" in P.columns else [None] * len(P)
    out: dict[str, set[str]] = {}
    for p in D:
        own = {
            i for r in rows.get(p, []) for i in items(ids[r]) if not i.startswith(tuple(prefixes))
        }
        out[p] = own or {"NA"}
    return out


def _miss_rows(
    sources: Mapping[str, set[str]],
    hits: pd.DataFrame,
    genes: Mapping[str, Any],
    circular: Mapping[str, int | None],
) -> list[dict[str, Any]]:
    """The proxy's miss rows, uncapped: one ``te1_result`` call per chunk of at most
    ``MAX_EXCEPTIONS`` peptides when the full call hit the cap."""
    peps = sorted(sources)
    hp = hits["peptide"].astype(str)
    out: list[dict[str, Any]] = []
    for i in range(0, len(peps), MAX_EXCEPTIONS):
        chunk = {p: sources[p] for p in peps[i : i + MAX_EXCEPTIONS]}
        res = te1_result(chunk, {}, hits[hp.isin(set(chunk))], genes, circular)
        out += [r for r in res.exceptions if r["kind"] == "miss"]
    return out


def _miss_cause(
    pep: str,
    own: Collection[str],
    hits: pd.DataFrame,
    why: Collection[str] | None,
    non_atg: Mapping[str, str],
) -> tuple[str, str | None]:
    """First match: non_atg_nterm, exact_other_gene, exact_other_class, tolerant_only,
    unplaced (sub-key: the class or the unmapped reason)."""
    if pep in non_atg:
        return "non_atg_nterm", None
    exact = hits[hits["route"].astype(str) == "exact"]
    if len(exact):
        cls = exact["locus_class"].astype(str)
        if any(
            c in INFRAME and g not in own
            for c, g in zip(cls, exact["gene"].astype(str), strict=True)
        ):
            return "exact_other_gene", None
        return "exact_other_class", min(cls)
    if len(hits):
        return "tolerant_only", None
    return "unplaced", unplaced_reason(why)


def te1_real(
    P: pd.DataFrame,
    peptide_hits: pd.DataFrame,
    unmapped: pd.DataFrame,
    genes: Iterable[Gene],
    circular: Mapping[str, int | None],
    non_atg: Mapping[str, str],
    cfg: GiapConfig,
    *,
    q_crossings: int = 0,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """TE1 on A2 (§11.6): share of **D₂** exact in frame at an own gene, through the proxy's
    ``te1_result`` with ``excluded = {}`` (L5); its 0.999 ``passed`` is ignored. Returns the TE1
    summary (no sequence) and the uncapped miss rows."""
    _, D, rows = split_contaminants(P, cfg)
    sources = _te1_sources(P, D, rows, cfg.input.contaminant_prefixes)
    intervals = {g.locus_tag: g.interval for g in genes}
    res = te1_result(sources, {}, peptide_hits, intervals, circular)
    exc = [r for r in res.exceptions if r["kind"] == "miss"]
    if len(exc) < res.counts["misses"]:
        exc = _miss_rows(sources, peptide_hits, intervals, circular)
    exc.sort(key=lambda r: (r["kind"], r["locus_tag"], r["peptide"]))
    missed = {r["peptide"] for r in exc}
    hp = peptide_hits["peptide"].astype(str)
    by_pep = {p: peptide_hits[hp == p] for p in missed}
    reasons = _reasons(unmapped)
    misses: dict[str, Any] = {
        "non_atg_nterm": 0,
        "exact_other_gene": 0,
        "exact_other_class": {},
        "tolerant_only": 0,
        "unplaced": {},
    }
    for p in sorted(missed):
        cause, sub = _miss_cause(p, sources[p], by_pep[p], reasons.get(p), non_atg)
        if sub is None:
            misses[cause] += 1
        else:
            misses[cause][sub] = misses[cause].get(sub, 0) + 1
    n, k = len(D), res.counts["exact_inframe"]
    s = share(k, n)
    na_peps = set(non_atg) & set(D)
    na_miss = len(na_peps & missed)
    qv = _q(P)
    at_cut = {
        p
        for p in D
        if (qs := [qv[i] for i in rows[p] if not math.isnan(qv[i])])
        and all(q == Q_ACCEPT for q in qs)
    }
    passed = s["value"] is not None and s["value"] >= TE1_MIN
    if q_crossings == 0:
        sq: dict[str, Any] = {
            "value": None,
            "ci_low": None,
            "ci_high": None,
            "reason": "no_crossings",
            "removed": None,
            "q_crossings": 0,
            "pass_becomes_fail": False,
        }
    else:
        sq = share(k - len(at_cut - missed), n - len(at_cut))
        sq |= {
            "removed": len(at_cut),
            "q_crossings": q_crossings,
            "pass_becomes_fail": passed and (sq["value"] or 0.0) < TE1_MIN,
        }
    te1 = {
        "n_peptides": n,
        "exact_inframe": k,
        "exact_inframe_masked": res.counts["exact_inframe_masked"],
        "share_exact_inframe": s["value"],
        "ci_low": s["ci_low"],
        "ci_high": s["ci_high"],
        "threshold": TE1_MIN,
        "passed": passed,
        "ci_low_below_threshold": s["ci_low"] is not None and s["ci_low"] < TE1_MIN,
        "target_q_na": sum(all(math.isnan(qv[i]) for i in rows[p]) for p in D),
        "non_atg_nterm": {"peptides": len(na_peps), "misses": na_miss},
        "misses": misses,
        "n_misses": len(missed),
        "share_excluding_non_atg": share(k - (len(na_peps) - na_miss), n - len(na_peps)),
        "share_without_q_at_cutoff": sq,
    }
    if n == 0:
        te1["reason"] = "empty_denominator"
    return te1, exc


def te1_diagnostic(P: pd.DataFrame, peptide_hits: pd.DataFrame, cfg: GiapConfig) -> dict[str, Any]:
    """A1 diagnostic (no pass flag): share of A1's **D** with an exact in-frame CDS/sORF row at
    any gene (A1's sources are UniProt accessions, not locus tags)."""
    _, D, _ = split_contaminants(P, cfg)
    h = peptide_hits
    ok = set(
        h.loc[
            (h["route"].astype(str) == "exact") & h["locus_class"].astype(str).isin(INFRAME),
            "peptide",
        ].astype(str)
    )
    k = sum(p in ok for p in D)
    return {"n_peptides": len(D), "exact_inframe_any_gene": k, "share": share(k, len(D))}


# ---------------------------------------------------------------- hypotheses


def _detected(genes: pd.DataFrame) -> set[str]:
    return set(
        genes.loc[genes["detection_status"].astype(str) == "detected", "locus_tag"].astype(str)
    )


def _p_greater(k: int, n: int) -> float | None:
    return float(binomtest(k, n, 0.5, alternative="greater").pvalue) if n else None


def hypothesis1(
    genes_a1: pd.DataFrame,
    genes_a2: pd.DataFrame,
    compare_detection: pd.DataFrame | None,
    subsets: Mapping[str, Collection[str]],
) -> dict[str, Any]:
    """H1 "A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes": gene
    effects, the exact one-sided TU-level sign test (verdict p) and the gene-level p
    ("assumes independent genes"), the "mostly" share, the verdict and failed clauses."""
    G2 = genes_a2[genes_a2["feature_type"].astype(str).isin(CODING)]
    G = [str(t) for t in G2["locus_tag"]]
    d1, d2 = _detected(genes_a1), _detected(genes_a2)
    b_set = {t for t in G if t in d2 and t not in d1}
    c_set = {t for t in G if t in d1 and t not in d2}
    b, c = len(b_set), len(c_set)
    units: dict[str, list[bool]] = {}
    for t, tu in zip(G, G2["tu_id"], strict=True):
        if t in b_set or t in c_set:
            u = units.setdefault(f"tu:{tu}" if clean(tu) else f"gene:{t}", [False, False])
            u[0 if t in b_set else 1] = True
    b_tu = sum(x and not y for x, y in units.values())
    c_tu = sum(y and not x for x, y in units.values())
    region = set(subsets.get("plasmid", ())) | set(subsets.get("prophage", ()))
    region |= set(subsets.get("bgc", ()))
    mostly = share(len(b_set & region), b)
    p_tu = _p_greater(b_tu, b_tu + c_tu)
    failed = [
        name
        for name, ok in (
            ("direction", b > c),
            ("test", p_tu is not None and p_tu <= H1_ALPHA),
            ("mostly", mostly["value"] is not None and mostly["value"] > H1_MOSTLY),
        )
        if not ok
    ]
    verdict = "not_evaluable" if b_tu + c_tu == 0 else "not_supported" if failed else "supported"
    warnings: list[dict[str, Any]] = []
    if compare_detection is not None:
        m = compare_detection[compare_detection["section"].astype(str) == "matrix"]
        for label, own in (("A1", d1), ("A2", d2)):
            r = m[m["run"].astype(str) == label]
            mat = set(
                r.loc[r["detection_status"].astype(str) == "detected", "locus_tag"].astype(str)
            )
            n_diff = len((mat ^ own) & set(G))
            if n_diff:
                warnings.append(
                    {"code": "compare_detection_mismatch", "run": label, "n_genes": n_diff}
                )
    return {
        "b": b,
        "c": c,
        "n_genes": len(G),
        "prop_a2_only": share(b, b + c),
        "coverage_diff": (b - c) / len(G) if G else None,
        "b_tu": b_tu,
        "c_tu": c_tu,
        "mixed_tu": sum(x and y for x, y in units.values()),
        "p_tu": p_tu,
        "p_gene": _p_greater(b, b + c),
        "unit_effect": share(b_tu, b_tu + c_tu),
        "mostly": mostly,
        "mostly_base_rate": share(len(set(G) & region), len(G)),
        "strain_specific_share": share(len(b_set & set(subsets.get("strain_specific", ()))), b),
        "verdict": verdict,
        "failed_clauses": [] if verdict == "not_evaluable" else failed,
        "warnings": warnings,
    }


def _status(evaluable: bool, holds: bool) -> str:
    return "not_evaluable" if not evaluable else "holds" if holds else "fails"


def hypothesis3(m4: Mapping[str, Any] | None, m5: Mapping[str, Any] | None) -> dict[str, Any]:
    """H3 rule (lead decision 5): rho ``ci_low >= 0.7``; per family with union >= 5,
    Jaccard >= 0.5; supported / not_supported / not_evaluable."""
    m4, m5 = m4 or {}, m5 or {}
    lo: float | None = m4.get("ci_low")
    rho = m4.get("spearman_rho")
    ev = rho is not None and lo is not None and "too_few_shared" not in items(m4.get("flags"))
    crit: dict[str, dict[str, Any]] = {
        "rho": {"status": _status(ev, ev and lo is not None and lo >= H3_RHO_MIN), "ci_low": lo}
    }
    for f in H3_FAMILIES:
        r = m5.get(f) or {}
        na, nb, ns = (r.get(k) for k in ("n_a", "n_b", "n_shared"))
        j: float | None = r.get("jaccard")
        union = None if na is None or nb is None or ns is None else int(na) + int(nb) - int(ns)
        flags = items(r.get("flags"))
        flagged = any(
            x in ("both_empty", "missing_row") or x.startswith("missing_in") for x in flags
        )
        e = union is not None and union >= H3_UNION_MIN and j is not None and not flagged
        crit[f] = {
            "status": _status(e, e and j is not None and j >= H3_J_MIN),
            "union": union,
            "jaccard": j,
        }
    st = [c["status"] for c in crit.values()]
    fam = st[1:]
    if "fails" in st:
        verdict = "not_supported"
    elif st[0] == "holds" and "holds" in fam:
        verdict = "supported"
    else:
        verdict = "not_evaluable"
    return {"verdict": verdict, "criteria": crit, "region_tests": "not compared in pass 1"}


def checklist(metrics: Mapping[str, Any], run_states: Mapping[str, Any]) -> list[dict[str, str]]:
    """§15.6 checklist in pass 1 (spec table order). ``run_states``: arm -> exit code (None when
    absent) and ``compare_ok`` (the three compare tables read under their schemas)."""
    ok = {a: run_states.get(a) == 0 for a in ("A1", "A2", "A3")}
    if all(ok.values()):
        runs = "pass"
    elif ok["A1"] and ok["A2"] and run_states.get("A3") is None:
        runs = "partial (A1, A2; A3 pending M10-10)"
    else:
        runs = "fail"
    te1 = metrics.get("te1") or {}
    te = "pass" if te1.get("passed") else "fail (no data)" if not te1.get("n_peptides") else "fail"
    m4 = metrics.get("m4") or {}
    cmp_ok = bool(run_states.get("compare_ok")) and "reason" not in m4
    return [
        {"bullet": "GIAP runs end to end on A1, A2 and A3", "status": runs},
        {"bullet": ">= 99.5 % of A2 peptides map exactly (TE1)", "status": te},
        {
            "bullet": "A1 and A2 detections and abundances are compared",
            "status": "pass" if cmp_ok else "fail",
        },
        {"bullet": "Measured runtimes replace the estimates", "status": "pending (M10-12)"},
        {"bullet": "TE2-TE5 are reported; failures explained", "status": "pending (M10-06)"},
    ]


# ---------------------------------------------------------------- render

PLACEHOLDER = re.compile(r"\{\{\s*([^}\s]+)\s*\}\}")
NOTES_MARK = "<!-- notes -->"
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
_MISSING = object()


def resolve(doc: Any, path: str) -> Any:
    cur = doc
    for part in path.split("."):
        if isinstance(cur, Mapping) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list) and part.isdigit() and int(part) < len(cur):
            cur = cur[int(part)]
        else:
            return _MISSING
    return cur


def fmt_value(key: str, v: Any) -> str:
    """Formats by key name: ``p_*`` 2 significant digits, ``*_s`` h:mm:ss, ``*_bytes`` GiB with
    1 decimal; other floats 3 significant digits; shares as ``v [lo, hi]``."""
    if isinstance(v, Mapping) and "value" in v:
        if v["value"] is None:
            return f"NA ({v.get('reason', 'no value')})"
        lo, hi = v.get("ci_low"), v.get("ci_high")
        ci = "" if lo is None else f" [{fmt_value('', lo)}, {fmt_value('', hi)}]"
        return fmt_value(key, v["value"]) + ci
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, list):
        return ", ".join(str(x) for x in v) if v else "none"
    if not isinstance(v, int | float):
        return str(v)
    if key.startswith("p_"):
        return f"{v:.2g}"
    if key.endswith("_s"):
        s = round(float(v))
        return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"
    if key.endswith("_bytes"):
        return f"{v / 1024**3:.1f} GiB"
    if isinstance(v, int) or float(v).is_integer():
        return str(int(v))
    return f"{v:.3g}" if abs(v) < 1000 else f"{v:.1f}"


def _numbers(doc: Any) -> set[str]:
    return set(re.findall(r"\d{4,}", json.dumps(doc, sort_keys=True)))


def check_notes(notes: str, metrics: Mapping[str, Any]) -> None:
    """Raise ``ValueError`` for a peptide-like string or a 4+ digit number not in metrics."""
    if m := PEPTIDE_LIKE.search(notes):
        raise ValueError(f"notes: peptide-like string at offset {m.start()}")
    known = _numbers(metrics)
    for num in re.findall(r"\d{4,}", notes):
        if num not in known:
            raise ValueError(f"notes: number {num} is not in metrics.json")


def missing_paths(template: str, metrics: Mapping[str, Any]) -> list[str]:
    return [p for p in PLACEHOLDER.findall(template) if resolve(metrics, p) is _MISSING]


def render(template: str, metrics: Mapping[str, Any], notes: str | None) -> str:
    """The report: each ``{{path}}`` becomes its formatted value and its source
    ``(metrics.json SOURCE_SEP parent)``; an unresolved path ``**missing: <path>**``."""
    if notes is not None:
        check_notes(notes, metrics)

    def sub(m: re.Match[str]) -> str:
        path = m.group(1)
        v = resolve(metrics, path)
        if v is _MISSING:
            return f"**missing: {path}**"
        parent = path.rsplit(".", 1)[0]
        return f"{fmt_value(path.rsplit('.', 1)[-1], v)} (metrics.json{SOURCE_SEP}{parent})"

    text = PLACEHOLDER.sub(sub, template)
    body = notes.strip() if notes and notes.strip() else "No notes."
    if NOTES_MARK in text:
        return text.replace(NOTES_MARK, body)
    return text + ("" if notes is None else f"\n## Notes\n\n{body}\n")


# ---------------------------------------------------------------- CLI (file I/O)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _need(path: Path) -> Path:
    if not path.exists():
        raise Precondition(f"missing input: {path}")
    return path


def _read_json(path: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads(_need(path).read_text(encoding="utf-8"))
    return doc


def _table(path: Path, stem: str) -> pd.DataFrame:
    from giap.io.run_view import _schema, read_schema_tsv

    schema = _schema(stem)
    assert schema is not None, stem
    return read_schema_tsv(path, stem, schema)


def _run_table(out: Path, name: str) -> Path | None:
    for rel in (f"tables/{name}", f"mapping/{name}"):
        if (out / rel).is_file():
            return out / rel
    return None


def _rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    text = path.read_text(encoding="utf-8")
    delim = "|" if "|" in text.split("\n", 1)[0] else "\t"
    return list(csv.DictReader(text.splitlines(), delimiter=delim))


def resolve_attempt(path: Path) -> tuple[Path, str | None]:
    """A run directory resolved through ``ATTEMPTS.tsv`` next to it (M10-04) to its final
    attempt (``.r<n>`` suffix included) and that attempt's head job ID."""
    rows = [
        r
        for r in _rows(path.parent / "ATTEMPTS.tsv")
        if re.fullmatch(re.escape(path.name) + r"(\.r\d+)?", Path(r.get("outdir", "")).name)
        and r.get("state") != "superseded"
    ]
    if not rows:
        return path, None
    return path.parent / Path(rows[-1]["outdir"]).name, rows[-1].get("job")


def _input_sha(manifest: Mapping[str, Any], role: str) -> str | None:
    return next(
        (i.get("sha256") for i in manifest.get("inputs", []) if i.get("role") == role), None
    )


def _cfg(manifest: Mapping[str, Any]) -> GiapConfig:
    raw = manifest.get("config")
    return GiapConfig.model_validate(raw) if raw else resolve_config()


def _spec_sha() -> str:
    res = subprocess.run(
        ["git", "-C", str(REPO), "log", "-1", "--format=%H", "--", SPEC],
        capture_output=True,
        text=True,
        check=False,
    )
    sha = res.stdout.strip()
    return sha if re.fullmatch(r"[0-9a-f]{40}", sha) else "0" * 40


def _arm(
    label: str,
    out: Path,
    pep: Path,
    conv: Path,
    summary: Path | None,
    genome_sha: str,
    sacct: list[dict[str, str]],
    job: str | None,
) -> dict[str, Any]:
    """Reads one arm and checks its preconditions."""
    from giap.io.peptides import read_peptide_table

    man = _read_json(out / "run_manifest.json")
    if man.get("exit_code") != 0:
        raise Precondition(f"{out / 'run_manifest.json'}: exit_code {man.get('exit_code')}")
    if _input_sha(man, "genome") != genome_sha:
        raise Precondition(f"{out / 'run_manifest.json'}: genome sha256 differs from --genome")
    if _input_sha(man, "peptides") != sha256(_need(pep)):
        raise Precondition(f"{pep}: sha256 differs from the run's peptides input ({out})")
    P, _ = read_peptide_table(pep)
    hits, unm, ds = (
        _run_table(out, n) for n in ("peptide_hits.tsv", "unmapped.tsv", "detection_summary.tsv")
    )
    trace = _rows(out / "logs" / "trace.tsv")
    head = [r for r in sacct if job is not None and r.get("JobID") == job]
    return {
        "label": label,
        "out": out,
        "manifest": man,
        "cfg": _cfg(man),
        "P": P,
        "genes": _table(_need(out / "tables" / "genes.tsv"), "genes"),
        "hits": None if hits is None else _table(hits, "peptide_hits"),
        "unm": None if unm is None else _table(unm, "unmapped"),
        "ds": None if ds is None else _table(ds, "detection_summary"),
        "conv": _read_json(conv),
        "sum": None if summary is None else _read_json(summary),
        "runtime": runtime(man, trace, head),
    }


def _missing(name: str) -> dict[str, str]:
    return {"reason": f"missing:{name}"}


def _cross_checks(arm: Mapping[str, Any], m1: Mapping[str, Any], warnings: list[Any]) -> None:
    n = m1["union"]["peptides"]["all"]
    counts = arm["conv"].get("counts", {})
    if n != counts.get("distinct_peptides_target_q01"):
        raise Inconsistent(
            f"{arm['label']}: metric 1 peptides {n} != convert.json "
            f"distinct_peptides_target_q01 {counts.get('distinct_peptides_target_q01')}"
        )
    if arm["sum"] is not None:
        diff = n - int(arm["sum"]["distinct_stripped_q01"])
        qc = int(counts.get("q_crossings", 0))
        if not 0 <= diff <= qc:
            warnings.append(
                {
                    "code": "raw_vs_written_q01_out_of_range",
                    "arm": arm["label"],
                    "raw_vs_written_q01": diff,
                    "q_crossings": qc,
                }
            )


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Every metric from the parsed arguments (reads only; raises on a precondition)."""
    out = a.out.resolve()
    if out == REPO or REPO in out.parents:
        raise Precondition(f"--out {a.out} is inside the repository")
    for p in (
        a.bakta_json,
        a.genome,
        a.k12_fasta,
        a.cmp,
        a.peptides_a1,
        a.peptides_a2,
        a.convert_a1,
        a.convert_a2,
        *(x for x in (a.summary_a1, a.summary_a2, a.sacct) if x),
    ):
        _need(p)
    genome_sha = sha256(a.genome)
    sacct = _rows(a.sacct)
    arms = {}
    for label, d, pep, conv, sm in (
        ("A1", a.a1, a.peptides_a1, a.convert_a1, a.summary_a1),
        ("A2", a.a2, a.peptides_a2, a.convert_a2, a.summary_a2),
    ):
        out_dir, job = resolve_attempt(d)
        arms[label] = _arm(label, _need(out_dir), pep, conv, sm, genome_sha, sacct, job)
    genes, circular = read_coding_genes(a.bakta_json, a.genome)
    a2 = arms["A2"]
    g2 = a2["genes"]
    tags = set(g2.loc[g2["feature_type"].astype(str).isin([*CODING, "pseudogene"]), "locus_tag"])
    if {g.locus_tag for g in genes} != {str(t) for t in tags}:
        raise Precondition(f"{a.bakta_json}: locus tags differ from {a2['out']}/tables/genes.tsv")
    from giap.io.protein_fasta import read_protein_fasta

    k12 = {r.id: r.sequence for r in read_protein_fasta(a.k12_fasta)}
    classes = strain_specific_set({g.locus_tag: g.aa for g in genes}, k12)
    warnings: list[Any] = []
    doc: dict[str, Any] = {"schema": SCHEMA, "spec": _spec_sha(), "warnings": warnings}
    doc["m1"], doc["m2"], doc["m3"], doc["m7"] = {}, {}, {}, {}
    for label, arm in arms.items():
        pre = arm["cfg"].input.contaminant_prefixes
        doc["m1"][label] = identifications(arm["P"], pre)
        _cross_checks(arm, doc["m1"][label], warnings)
        if arm["hits"] is None or arm["unm"] is None:
            doc["m2"][label] = _missing("peptide_hits.tsv/unmapped.tsv")
        else:
            m2 = mapping_shares(arm["P"], arm["hits"], arm["unm"], arm["cfg"])
            gi: dict[str, Any] = (
                genome_inconsistent(arm["ds"])
                if arm["ds"] is not None
                else _missing("detection_summary.tsv")
            )
            tu = m2["counts"]["tolerant"] + m2["counts"]["unmapped"]
            gi["consistent"] = None if gi.get("count") is None else gi["count"] == tu
            m2["genome_inconsistent"] = gi
            doc["m2"][label] = m2
        doc["m3"][label] = coverage(arm["genes"], classes)
        doc["m7"][label] = arm["runtime"]
    cmp_tables: dict[str, pd.DataFrame | None] = {}
    for stem in ("compare_detection", "compare_abundance", "compare_stability"):
        # giap compare writes CMP/tables/<stem>.tsv; a flat CMP/<stem>.tsv is read too
        f = next((x for x in (a.cmp / "tables", a.cmp) if (x / f"{stem}.tsv").is_file()), None)
        cmp_tables[stem] = None if f is None else _table(f / f"{stem}.tsv", stem)
    ab, st = cmp_tables["compare_abundance"], cmp_tables["compare_stability"]
    labels = ("A1", "A2")
    doc["m4"] = (
        abundance_agreement(ab, labels) if ab is not None else _missing("compare_abundance.tsv")
    )
    doc["m5"] = stability(st, labels) if st is not None else _missing("compare_stability.tsv")
    exceptions: list[dict[str, Any]] = []
    if a2["hits"] is None or a2["unm"] is None:
        doc["te1"] = {
            **_missing("peptide_hits.tsv"),
            "passed": False,
            "n_peptides": 0,
            "share_exact_inframe": None,
        }
    else:
        qc = int(a2["conv"].get("counts", {}).get("q_crossings", 0))
        doc["te1"], exceptions = te1_real(
            a2["P"],
            a2["hits"],
            a2["unm"],
            genes,
            circular,
            non_atg_nterm(genes),
            a2["cfg"],
            q_crossings=qc,
        )
    a1 = arms["A1"]
    doc["te1_a1_diagnostic"] = (
        _missing("peptide_hits.tsv")
        if a1["hits"] is None
        else te1_diagnostic(a1["P"], a1["hits"], a1["cfg"])
    )
    subsets = {
        "plasmid": set(g2.loc[g2["replicon_type"].astype(str) == "plasmid", "locus_tag"]),
        "strain_specific": {t for t, c in classes.items() if c == "strain_specific"},
    }
    for kind in ("prophage", "bgc"):
        subsets[kind] = {
            str(t)
            for t, r in zip(g2["locus_tag"], g2["region_ids"], strict=True)
            if kind in _region_types(r)
        }
    doc["h1"] = hypothesis1(a1["genes"], g2, cmp_tables["compare_detection"], subsets)
    warnings += doc["h1"]["warnings"]
    doc["h3"] = hypothesis3(
        doc["m4"] if "reason" not in doc["m4"] else None,
        doc["m5"] if "reason" not in doc["m5"] else None,
    )
    states: dict[str, Any] = {
        "A1": 0,
        "A2": 0,
        "A3": None,
        "compare_ok": all(t is not None for t in cmp_tables.values()),
    }
    if a.a3 is not None:
        states["A3"] = _read_json(resolve_attempt(a.a3)[0] / "run_manifest.json").get("exit_code")
    doc["checklist"] = checklist(doc, states)
    files = {
        "peptides_a1": a.peptides_a1,
        "peptides_a2": a.peptides_a2,
        "convert_a1": a.convert_a1,
        "convert_a2": a.convert_a2,
        "summary_a1": a.summary_a1,
        "summary_a2": a.summary_a2,
        "bakta_json": a.bakta_json,
        "genome": a.genome,
        "k12_fasta": a.k12_fasta,
        "sacct": a.sacct,
        "manifest_a1": arms["A1"]["out"] / "run_manifest.json",
        "manifest_a2": arms["A2"]["out"] / "run_manifest.json",
    }
    doc["inputs"] = {
        k: {"path": str(Path(p).resolve()), "sha256": sha256(Path(p))}
        for k, p in files.items()
        if p is not None
    }
    return doc, exceptions


def _dump(doc: Any) -> str:
    return json.dumps(doc, sort_keys=True, indent=2) + "\n"


def write_outputs(
    out: Path, doc: Mapping[str, Any], exceptions: Sequence[Mapping[str, Any]]
) -> None:
    from giap.io.tsv import write_tsv

    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(_dump(doc), encoding="utf-8")
    (out / "te1.json").write_text(_dump(doc["te1"]), encoding="utf-8")
    cols = ["kind", "locus_tag", "peptide", "hits"]
    write_tsv(out / "te1_exceptions.tsv", cols, [[r.get(c) for c in cols] for r in exceptions])


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for name in (
        "a1",
        "a2",
        "cmp",
        "peptides-a1",
        "peptides-a2",
        "convert-a1",
        "convert-a2",
        "bakta-json",
        "genome",
        "k12-fasta",
        "out",
    ):
        c.add_argument(f"--{name}", type=Path, required=True)
    for name in ("summary-a1", "summary-a2", "sacct", "a3"):
        c.add_argument(f"--{name}", type=Path, default=None)
    r = sub.add_parser("render")
    r.add_argument("--metrics", type=Path, required=True)
    r.add_argument("--notes", type=Path, default=None)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--template", type=Path, default=Path(__file__).with_name("report_template.md"))
    a = ap.parse_args(argv)
    try:
        if a.cmd == "compute":
            doc, exc = compute(a)
            write_outputs(a.out, doc, exc)
            print(
                f"wrote {a.out}/metrics.json; TE1 passed={doc['te1'].get('passed')}; "
                f"{len(doc['warnings'])} warning(s)"
            )
            return 0
        metrics = _read_json(a.metrics)
        template = _need(a.template).read_text(encoding="utf-8")
        notes = _need(a.notes).read_text(encoding="utf-8") if a.notes else None
        text = render(template, metrics, notes)
        if missing := missing_paths(template, metrics):
            print(f"render: missing paths in {a.metrics}: {', '.join(missing)}", file=sys.stderr)
            return 1
        a.out.write_text(text, encoding="utf-8")
        return 0
    except Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except (Inconsistent, ValueError) as err:
        print(f"inconsistent: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
