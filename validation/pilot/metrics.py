"""Pilot metrics for the arms A1 and A2 (GIAP_SPEC.md §15.5 metrics 1-5 and 7, §11.6 TE1 on real
data, §15.6 checklist, §15.7 hypotheses 1 and 3); spec docs/specs/M10-05-pilot-metrics.md.

Pure functions (frames and dicts in, dicts out) plus a thin CLI that reads and writes::

    python validation/pilot/metrics.py compute --a1 OUT_A1 --a2 OUT_A2 --cmp CMP \\
        --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --convert-a1 CONV_A1 --convert-a2 CONV_A2 \\
        [--summary-a1 SUM_A1 --summary-a2 SUM_A2] --bakta-json BAKTA_JSON --genome GENOME \\
        --k12-fasta K12_FASTA [--sacct SACCT] [--a3 OUT_A3 --peptides-a3 PEP_A3 \\
        --convert-a3 CONV_A3 --cmp3 CMP3 [--strain-db STRAIN_DB]] --out METRICS_DIR
    python validation/pilot/metrics.py render --metrics METRICS_DIR/metrics.json \\
        [--notes METRICS_DIR/notes.md] --out METRICS_DIR/REPORT.md

Exit codes: 0 written; 1 an internal inconsistency only (metric 1 vs ``convert.json``, a render
placeholder missing under a computed metric, a rejected notes file); 2 a precondition (a missing,
unreadable or schema-invalid input, a failed run, a hash mismatch, git unavailable, ``--out``
inside the checkout). Nothing is written on exit 1 or 2.
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime as dt
import hashlib
import json
import math
import re
import subprocess
import sys
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from importlib import resources
from itertools import pairwise
from pathlib import Path
from typing import Any, TypeVar

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from validation.proxies.te_proxies import INFRAME, non_atg_nterm, read_coding_genes, te1_result

from giap.config import GiapConfig
from giap.io.peptides import validate_peptides
from giap.mapping.digest import digest
from giap.model.genes import Gene
from giap.stats.cutoffs import wilson_ci
from giap.steps.mapping import select_queries

REPO = Path(__file__).resolve().parents[2]
SPEC = "docs/specs/M10-05-pilot-metrics.md"
SCHEMA = "giap-pilot-metrics/1"
ARMS = ("A1", "A2")
Q_ACCEPT = 0.01  # accepted: q_value <= 0.01 (inclusive), the converter's written value
TE1_MIN = 0.995  # §11.6, point estimate (R3: not tuned)
H1_ALPHA, H1_MOSTLY = 0.05, 0.5
H3_RHO_MIN, H3_J_MIN, H3_UNION_MIN = 0.7, 0.5, 5
UNMAPPED_REASONS = ("no_hit", "rejected_placement_fdr", "below_identity", "too_short")
H3_FAMILIES = ("islands", "operon_flags", "enriched_categories")
OTHER_FAMILIES = ("novel_orfs", "alt_sites")
CODING = ("cds", "sorf")
STRAIN_CLASSES = ("strain_specific", "k12_shared", "no_peptides")
SOURCE_SEP = " \u203a "
T = TypeVar("T")


class Precondition(Exception):
    """Exit 2: a missing or invalid input, a failed run, a hash mismatch, ``--out`` inside the
    checkout."""


class Inconsistent(Exception):
    """Exit 1: two computations of the same number disagree (or a rejected notes file)."""


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
        return {"value": None, "ci_low": None, "ci_high": None, "k": 0, "n": 0,
                "reason": "empty_denominator"}  # fmt: skip
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
        if group := tuple(sorted(set(items(ids)))):  # a row without protein_ids is no group
            groups.setdefault(group, set()).add(str(pep))
    gcls = {g: all(i.startswith(tuple(prefixes)) for i in g) for g in groups}
    min2 = {g: c for g, c in gcls.items() if len(groups[g]) >= 2}
    return {
        "precursors": _class_counts(prec),
        "peptides": _class_counts(peps),
        "protein_groups": _class_counts(gcls),
        "protein_groups_min2": _class_counts(min2),
    }


def identifications(P: pd.DataFrame, prefixes: Sequence[str]) -> dict[str, Any]:
    """Metric 1: distinct precursors, peptides, protein groups (and groups with >= 2 accepted
    peptides) among accepted target rows of the raw table, by contaminant class; union and per
    run."""
    acc = P.loc[_target(P) & (_q(P) <= Q_ACCEPT)]
    runs = sorted({str(r) for r in _col(acc, "run", None) if clean(r) is not None})
    per_run = {r: _ident(acc.loc[acc["run"].astype(str) == r], prefixes) for r in runs}
    return {"union": _ident(acc, prefixes), "per_run": per_run}


# ---------------------------------------------------------------- metric 2 and TE1 universes


def split_contaminants(
    P: pd.DataFrame, cfg: GiapConfig
) -> tuple[list[str], list[str], pd.DataFrame]:
    """(U, D, R) as GIAP's mapping step sees the table: **U** =
    ``select_queries(validate_peptides(P)[0], cfg)``; R its validated target rows that
    ``select_queries`` keeps (q <= q_value_max or NA); D = U minus the peptides whose rows of R
    are all contaminant-only (``validate_peptides``' ``is_contaminant``)."""
    V, _ = validate_peptides(P, cfg.input)
    U = [str(p) for p in select_queries(V, cfg)["peptide"]]
    keep = ~V["is_decoy"].to_numpy(dtype=bool) & ~(_q(V) > cfg.input.q_value_max)
    R = V.loc[keep & V["peptide"].astype(str).isin(U).to_numpy(dtype=bool)]
    cont = R.groupby(R["peptide"].astype(str))["is_contaminant"].all()
    return U, [p for p in U if not cont.get(p, False)], R


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


def gene_subsets(genes: pd.DataFrame, strain_classes: Mapping[str, str]) -> dict[str, set[str]]:
    """Locus tags of **G** (cds, sorf) per metric-3 subset, in the report's row order (a tag
    missing from ``strain_classes`` counts as ``no_peptides``)."""
    G = genes[genes["feature_type"].astype(str).isin(CODING)]
    tags, rep = list(G["locus_tag"].astype(str)), list(G["replicon_type"].astype(str))
    contig = list(G["contig"].astype(str))
    regions = [_region_types(v) for v in G["region_ids"]]
    cls = [strain_classes.get(t, "no_peptides") for t in tags]

    def pick(ok: Iterable[bool]) -> set[str]:
        return {t for t, o in zip(tags, ok, strict=True) if o}

    subs = {"all": set(tags), "chromosome": pick(r == "chromosome" for r in rep)}
    subs["plasmid"] = pick(r == "plasmid" for r in rep)
    for c in sorted({c for c, r in zip(contig, rep, strict=True) if r == "plasmid"}):
        subs[f"plasmid:{c}"] = subs["plasmid"] & pick(x == c for x in contig)
    subs["unknown"] = pick(r not in ("chromosome", "plasmid") for r in rep)
    for kind in ("prophage", "bgc"):
        subs[kind] = pick(kind in r for r in regions)
    for c in STRAIN_CLASSES:
        subs[c] = pick(x == c for x in cls)
    return subs


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
    tags = genes["locus_tag"].astype(str)
    coding = genes["feature_type"].astype(str).isin(CODING)
    out = {
        name: _cov(genes[coding & tags.isin(s)])
        for name, s in gene_subsets(genes, strain_specific).items()
    }
    out["pseudogene"] = _cov(genes[genes["feature_type"].astype(str) == "pseudogene"])
    return {"subsets": out, "order": list(out)}


# ---------------------------------------------------------------- metrics 4, 5


ABUNDANCE_FIELDS = ("n_shared", "n_units", "spearman_rho", "ci_low", "ci_high", "method", "n_boot",
                    "ba_bias_log2", "ba_loa_low", "ba_loa_high", "flags")  # fmt: skip
STABILITY_FIELDS = ("n_a", "n_b", "n_shared", "jaccard", "flags")


def _pair(frame: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    a, b = frame["run_a"].astype(str), frame["run_b"].astype(str)
    return frame[[{x, y} == set(labels) for x, y in zip(a, b, strict=True)]]


def abundance_agreement(compare_abundance: pd.DataFrame, labels: Sequence[str]) -> dict[str, Any]:
    """Metric 4: the A1/A2 ``compare_abundance.tsv`` row, copied unchanged (no row: every field
    None and ``reason = "missing_row"``)."""
    row = _pair(compare_abundance, labels)
    if not len(row):
        return {**dict.fromkeys(ABUNDANCE_FIELDS), "reason": "missing_row"}
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
            return {**dict.fromkeys(STABILITY_FIELDS), "reason": "missing_row", "note": "no row"}
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


STEP_FIELDS = ("name", "status", "n_tasks", "n_retries", "wall_time_s", "cpu_time_s",
               "peak_rss_bytes")  # fmt: skip
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
    head = [r for r in sacct_rows if "." not in str(r.get("JobID", ""))]
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
    _, D, R = split_contaminants(P, cfg)
    pre = tuple(cfg.input.contaminant_prefixes)
    own: dict[str, set[str]] = {p: set() for p in D}
    qs: dict[str, list[float]] = {p: [] for p in D}
    for pep, ids, q in zip(R["peptide"].astype(str), R["protein_ids"], _q(R), strict=True):
        if pep in own:
            own[pep] |= {i for i in items(ids) if not i.startswith(pre)}
            qs[pep].append(q)
    intervals = {g.locus_tag: g.interval for g in genes}
    # "NA": a peptide without a strain item still needs a source for the proxy's miss row
    sources = {p: s or {"NA"} for p, s in own.items()}
    res = te1_result(sources, {}, peptide_hits, intervals, circular, max_exceptions=None)
    exc = sorted(res.exceptions, key=lambda r: (r["kind"], r["locus_tag"], r["peptide"]))
    missed = {r["peptide"] for r in exc}
    hp = peptide_hits["peptide"].astype(str)
    reasons = _reasons(unmapped)
    misses: dict[str, Any] = {"non_atg_nterm": 0, "exact_other_gene": 0, "exact_other_class": {},
                              "tolerant_only": 0, "unplaced": {}}  # fmt: skip
    for p in sorted(missed):
        cause, sub = _miss_cause(p, sources[p], peptide_hits[hp == p], reasons.get(p), non_atg)
        if sub is None:
            misses[cause] += 1
        else:
            misses[cause][sub] = misses[cause].get(sub, 0) + 1
    n, k = len(D), res.counts["exact_inframe"]
    s = share(k, n)
    na_peps = set(non_atg) & set(D)
    na_miss = len(na_peps & missed)
    known = [[q for q in qs[p] if not math.isnan(q)] for p in D]
    at_cut = {p for p, v in zip(D, known, strict=True) if v and all(q == Q_ACCEPT for q in v)}
    passed = s["value"] is not None and s["value"] >= TE1_MIN
    if q_crossings == 0:
        sq = {**share(0, 0), "reason": "no_crossings", "removed": None, "q_crossings": 0,
              "pass_becomes_fail": False}  # fmt: skip
    else:
        sq = share(k - len(at_cut - missed), n - len(at_cut))
        # an empty superset denominator (value None) counts as a fail, with its reason
        sq |= {"removed": len(at_cut), "q_crossings": q_crossings,
               "pass_becomes_fail": passed and (sq["value"] or 0.0) < TE1_MIN}  # fmt: skip
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
        "target_q_na": sum(not v for v in known),
        "unknown_source": sum(not (own[p] & intervals.keys()) for p in D),
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
    """H3 rule (lead decision 5): rho ``ci_low >= 0.7`` (a rho without ``ci_low`` is not
    evaluable); per family with union >= 5, Jaccard >= 0.5; supported / not_supported /
    not_evaluable. Every criterion carries its numbers whatever its status (lead decision 6)."""
    m4, m5 = m4 or {}, m5 or {}
    lo: float | None = m4.get("ci_low")
    rho = m4.get("spearman_rho")
    ev = rho is not None and lo is not None and "too_few_shared" not in items(m4.get("flags"))
    rho_ok = ev and lo is not None and lo >= H3_RHO_MIN
    rec = {"status": _status(ev, rho_ok), "rho": rho, "ci_low": lo, "ci_high": m4.get("ci_high")}
    crit: dict[str, dict[str, Any]] = {"rho": rec}
    for f in H3_FAMILIES:
        r = m5.get(f) or {}
        na, nb, ns = (r.get(k) for k in ("n_a", "n_b", "n_shared"))
        j: float | None = r.get("jaccard")
        union = None if na is None or nb is None or ns is None else int(na) + int(nb) - int(ns)
        flags = items(r.get("flags"))
        flagged = any(x == "both_empty" or x.startswith("missing_in") for x in flags)
        e = union is not None and union >= H3_UNION_MIN and j is not None and not flagged
        crit[f] = {"status": _status(e, e and j is not None and j >= H3_J_MIN), "n_a": na,
                   "n_b": nb, "n_shared": ns, "union": union, "jaccard": j}  # fmt: skip
    st = [c["status"] for c in crit.values()]
    if "fails" in st:
        verdict = "not_supported"
    elif st[0] == "holds" and "holds" in st[1:]:
        verdict = "supported"
    else:
        verdict = "not_evaluable"
    return {"verdict": verdict, "criteria": crit, "region_tests": "not compared in pass 1"}


# ---------------------------------------------------------------- metric 6 (M10-10b)

#: the class enum of novel_orfs.schema.json
NOVEL_CLASSES = ("nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
                 "pseudogene", "sorf")  # fmt: skip
BUCKETS = ("called", "unestimated", "not_called")
ALT_TYPES = ("substitution", "swap", "block", "frameshift")  # variants.schema.json enum
CONFIDENCE = ("high", "candidate")
ORF_KEY = ("contig", "strand", "end")  # the stop codon (spec Q6): peptide-independent
ALT_KEY = ("contig", "codon_start", "strand", "alt_aa")  # compare_stability's alt_sites key
A3_NOT_RUN = {"reason": "a3_not_run"}  # lead decision after the M10-15a review


def _str_eq(df: pd.DataFrame, col: str, value: str) -> np.ndarray:
    out: np.ndarray = (df[col].astype(str) == value).to_numpy(dtype=bool)
    return out


def _called(df: pd.DataFrame) -> np.ndarray:
    out: np.ndarray = df["called"].fillna(False).astype(bool).to_numpy(dtype=bool)
    return out


def novel_summary(novel_orfs: pd.DataFrame) -> dict[str, Any]:
    """Metric 6 per arm: ``novel_orfs.tsv`` rows by class x bucket (``called``; ``unestimated``:
    not called and ``class_fdr = unestimated``; ``not_called``: the rest) and the totals."""
    called = _called(novel_orfs)
    unest = ~called & _str_eq(novel_orfs, "class_fdr", "unestimated")
    bucket = np.where(called, "called", np.where(unest, "unestimated", "not_called"))
    cls = novel_orfs["class"].astype(str).to_numpy()
    by = {c: {b: int(((cls == c) & (bucket == b)).sum()) for b in BUCKETS} for c in NOVEL_CLASSES}
    totals = {b: int((bucket == b).sum()) for b in BUCKETS}
    return {"n_rows": len(novel_orfs), "by_class": by, "totals": totals}


def variant_summary(variants: pd.DataFrame) -> dict[str, Any]:
    """Metric 6 per arm: ``call = ALT`` rows by ``alt_type`` x ``confidence``, the ``mixed``
    count and ``class_error`` (rows with a numeric ``class_error_rate`` vs ``unestimated``)."""
    alt = _str_eq(variants, "call", "ALT")
    t = variants["alt_type"].astype(str).to_numpy()
    c = variants["confidence"].astype(str).to_numpy()
    rate = variants["class_error_rate"]
    numeric = pd.to_numeric(rate, errors="coerce").notna().to_numpy(dtype=bool)
    return {
        "alt": {a: {k: int((alt & (t == a) & (c == k)).sum()) for k in CONFIDENCE}
                for a in ALT_TYPES},
        "mixed": int(_str_eq(variants, "call", "mixed").sum()),
        "class_error": {"numeric": int(numeric.sum()),
                        "unestimated": int(_str_eq(variants, "class_error_rate", "unestimated")
                                           .sum())},
    }  # fmt: skip


def _keys(df: pd.DataFrame, key: Sequence[str]) -> list[tuple[Any, ...]]:
    return [tuple(clean(v) for v in row) for row in zip(*(df[k] for k in key), strict=True)]


def overlap(a2: pd.DataFrame, a3: pd.DataFrame, key: Sequence[str]) -> dict[str, Any]:
    """A2/A3 overlap of one set (frames already filtered) on the ``key`` columns (keys
    de-duplicated): counts, union, Jaccard (None, ``reason = both_empty``, at union 0) and the
    two Wilson shares ``a3_confirmed_by_a2`` and ``a2_found_by_a3`` (descriptive)."""
    x, y = set(_keys(a2, key)), set(_keys(a3, key))
    shared, union = len(x & y), len(x | y)
    out: dict[str, Any] = {"n_a2": len(x), "n_a3": len(y), "n_shared": shared, "union": union,
                           "jaccard": shared / union if union else None,
                           "a3_confirmed_by_a2": share(shared, len(y)),
                           "a2_found_by_a3": share(shared, len(x))}  # fmt: skip
    if not union:
        out["reason"] = "both_empty"
    return out


def _g6(v: Any) -> str | None:
    return None if clean(v) is None else f"{float(v):.6g}"


def _a2_a3_row(compare_stability: pd.DataFrame | None, family: str) -> dict[str, Any] | None:
    """The (A2, A3) row of ``family`` oriented as A2 = ``n_a``, A3 = ``n_b`` (None if absent)."""
    if compare_stability is None:
        return None
    rows = _pair(compare_stability, ("A2", "A3"))
    rows = rows[rows["family"].astype(str) == family]
    if not len(rows):
        return None
    rec = {k: clean(rows[k].iloc[0]) for k in (*STABILITY_FIELDS, "run_a")}
    if rec["run_a"] == "A3":
        rec["n_a"], rec["n_b"] = rec["n_b"], rec["n_a"]
    return rec


def _novel_overlap(
    a2: pd.DataFrame, a3: pd.DataFrame, cmp_row: Mapping[str, Any] | None, warnings: list[Any]
) -> dict[str, Any]:
    c2, c3 = a2[_called(a2)], a3[_called(a3)]
    called = overlap(c2, c3, ORF_KEY)
    cls2: dict[tuple[Any, ...], str] = {}
    for k, c in zip(_keys(c2, ORF_KEY), c2["class"].astype(str), strict=True):
        cls2.setdefault(k, c)
    cls3: dict[tuple[Any, ...], str] = {}
    for k, c in zip(_keys(c3, ORF_KEY), c3["class"].astype(str), strict=True):
        cls3.setdefault(k, c)
    called["class_agree"] = sum(cls2[k] == c for k, c in cls3.items() if k in cls2)
    called["by_class_a3"] = {
        cl: {"n": sum(c == cl for c in cls3.values()),
             "shared": sum(c == cl and k in cls2 for k, c in cls3.items())}
        for cl in NOVEL_CLASSES
    }  # fmt: skip
    for arm, df in (("A2", a2), ("A3", a3)):
        if len(keys := _keys(df, ORF_KEY)) != len(set(keys)):
            warnings.append({"code": "duplicate_orf_key", "arm": arm,
                             "n_duplicates": len(keys) - len(set(keys))})  # fmt: skip
    orf_j = None if cmp_row is None else cmp_row["jaccard"]
    if cmp_row is None:
        warnings.append({"code": "compare_row_missing", "family": "novel_orfs"})
    else:
        if (cmp_row["n_a"], cmp_row["n_b"]) != (called["n_a2"], called["n_a3"]):
            warnings.append({"code": "duplicate_orf_key", "compare_n_a2": cmp_row["n_a"],
                             "compare_n_a3": cmp_row["n_b"], "n_a2": called["n_a2"],
                             "n_a3": called["n_a3"]})  # fmt: skip
        if _g6(orf_j) != _g6(called["jaccard"]):
            warnings.append({"code": "novel_orf_key_differs",
                             "stop_codon_jaccard": called["jaccard"],
                             "orf_id_jaccard": orf_j})  # fmt: skip
    return {"called": called, "any_row": overlap(a2, a3, ORF_KEY), "orf_id_jaccard": orf_j}


def _alt_overlap(
    a2: pd.DataFrame, a3: pd.DataFrame, cmp_row: Mapping[str, Any] | None
) -> dict[str, Any]:
    """ALT sites (``call = ALT``) and each arm's own ``high`` rows; the ALT set must equal
    compare's (A2, A3) ``alt_sites`` row (counts exactly, Jaccard after ``%.6g``)."""
    x2, x3 = a2[_str_eq(a2, "call", "ALT")], a3[_str_eq(a3, "call", "ALT")]
    alt = overlap(x2, x3, ALT_KEY)
    hi = overlap(x2[_str_eq(x2, "confidence", "high")], x3[_str_eq(x3, "confidence", "high")],
                 ALT_KEY)  # fmt: skip
    if cmp_row is None:
        raise Inconsistent("metric 6: compare_stability has no (A2, A3) alt_sites row")
    ours = (alt["n_a2"], alt["n_a3"], alt["n_shared"], _g6(alt["jaccard"]), alt["union"] == 0)
    theirs = (cmp_row["n_a"], cmp_row["n_b"], cmp_row["n_shared"], _g6(cmp_row["jaccard"]),
              "both_empty" in items(cmp_row["flags"]))  # fmt: skip
    if ours != theirs:
        raise Inconsistent(
            f"metric 6: ALT sites (n_a2, n_a3, n_shared, jaccard, both_empty) {ours} != "
            f"compare_stability alt_sites (A2, A3) {theirs}"
        )
    return {"alt": alt, "high": hi}


def metric6(
    tables: Mapping[str, Mapping[str, pd.DataFrame | None]],
    compare_stability: pd.DataFrame | None,
) -> dict[str, Any]:
    """Metric 6 (§15.5): per-arm summaries (an absent table ``missing:<table>.tsv``), the A2/A3
    overlap (``a3_not_run`` without an ``A3`` entry) with its cross-checks against compare's
    (A2, A3) rows, and ``warnings``. Raises ``Inconsistent`` when the ALT set differs from
    compare's."""
    out: dict[str, Any] = {}
    warnings: list[Any] = []
    for arm, t in tables.items():
        out[arm] = {
            "novel_orfs": _or_missing("novel_orfs.tsv", novel_summary, t.get("novel_orfs")),
            "variants": _or_missing("variants.tsv", variant_summary, t.get("variants")),
        }
    if "A3" not in tables:
        out["A3"] = dict(A3_NOT_RUN)
        out["overlap_a2_a3"] = dict(A3_NOT_RUN)
        out["warnings"] = warnings
        return out
    v3 = tables["A3"].get("variants")
    if v3 is not None and (n := int(pd.to_numeric(v3["class_error_rate"], errors="coerce")
                                    .notna().sum())):  # fmt: skip
        warnings.append({"code": "a3_class_error_estimated", "n_rows": n})
    a2, a3 = tables.get("A2") or {}, tables["A3"]
    ov: dict[str, Any] = {}
    n2, n3 = a2.get("novel_orfs"), a3.get("novel_orfs")
    ov["novel_orfs"] = (
        _missing("novel_orfs.tsv")
        if n2 is None or n3 is None
        else _novel_overlap(n2, n3, _a2_a3_row(compare_stability, "novel_orfs"), warnings)
    )
    w2, w3 = a2.get("variants"), a3.get("variants")
    ov["alt_sites"] = (
        _missing("variants.tsv")
        if w2 is None or w3 is None
        else _alt_overlap(w2, w3, _a2_a3_row(compare_stability, "alt_sites"))
    )
    out["overlap_a2_a3"] = ov
    out["warnings"] = warnings
    return out


# ---------------------------------------------------------------- hypothesis 4 (M10-10b)

H4_MIN_N = 20  # definition 6 (R3: not tuned)
H4_LEVEL = 0.95
H4_VERDICTS = ("not_supported", "consistent_with_lower_bound", "not_evaluable")
PLACEMENT_LOG = re.compile(r"placement: \d+ queries \((\d+) de novo\)")
SEARCHABLE_DIGEST = (
    "GIAP's digest (trypsin_p, at most 1 missed cleavage, 7-30 residues, no N-terminal Met "
    "excision; not DIA-NN's settings) of --strain-db, I folded to L"
)
H4_CAVEATS = (
    "The score bins span the de novo scores of (0.8, 1] only (Cascadia's reported range).",
    "c_value comes from A2's own exact matches and the confirmations from A2's identifications, "
    "so A2 enters both sides: a shared dependence, not an independent check.",
    "Class error rates of A3's ALT calls are unestimated (no spectrum pairing).",
    "Dependent rows: nested or overlapping de novo peptides of one locus share correctness, so "
    "the genome bin CI can be too narrow and a not_supported can come from that undercoverage; "
    "nested_share shows how clustered each bin is.",
    "Contaminant rows count in c_value (they never match the genome), which biases the genome "
    "precision up and lowers the power to refute; de novo contaminant peptides are not in A2's "
    "set and count as unconfirmed.",
    "A2's identifications are at 1 % FDR, so the lower bound can be off by up to that share "
    "(not corrected).",
    "A failure confined to a single score bin reads consistent_with_lower_bound: not_supported "
    "needs two consecutive reportable bins below the bound (the pre-registered fallback).",
)


def _fold(p: str) -> str:
    return p.replace("I", "L")


def _selected_denovo(P: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """The placement step's de novo queries (``validate_peptides``, ``aggregate_peptides``
    (score = max over a peptide's rows), ``select_queries``), before the length-class filter."""
    from giap.io.peptides import aggregate_peptides

    V, _ = validate_peptides(P, cfg.input)
    Q = select_queries(aggregate_peptides(V)[0], cfg)
    keep = (Q["source"].astype(str) == "denovo").to_numpy(dtype=bool)
    return Q.loc[keep, ["peptide", "score", "source"]].reset_index(drop=True)


def denovo_count(P: pd.DataFrame, cfg: GiapConfig) -> int:
    """The placement log's "(M de novo)": de novo queries before the length-class filter,
    unscored ones included (B1)."""
    return len(_selected_denovo(P, cfg))


def denovo_queries(P: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame:
    """Definition 1: the de novo queries as the mapping step bins them (``_selected_denovo``),
    keeping rows whose length falls in a ``mapping.length_classes`` class."""
    from giap.steps.mapping import length_label

    Q = _selected_denovo(P, cfg)
    classes = cfg.mapping.length_classes
    in_cls = [length_label(len(p), classes) is not None for p in Q["peptide"].astype(str)]
    return Q.loc[np.asarray(in_cls, dtype=bool)].reset_index(drop=True)


def placement_denovo_count(log_text: str) -> int | None:
    """The "(M de novo)" count of the last ``placement: N queries (M de novo)`` log line."""
    found = PLACEMENT_LOG.findall(log_text)
    return int(found[-1]) if found else None


def _nested(peps: Sequence[str]) -> set[str]:
    """Peptides that are a substring of another peptide of ``peps`` (equal strings included)."""
    lens = sorted({len(p) for p in peps})
    subs: set[str] = set()
    for p in set(peps):
        for n in (x for x in lens if x < len(p)):
            subs.update(p[i : i + n] for i in range(len(p) - n + 1))
    seen: dict[str, int] = {}
    for p in peps:
        seen[p] = seen.get(p, 0) + 1
    return {p for p in peps if p in subs or seen[p] > 1}


def h4_bin(n_b: int, k_b: int, genome: Mapping[str, Any], n_reportable: int) -> dict[str, Any]:
    """Definition 7 for one reportable bin: ``L_b`` (Wilson 95 %), ``L_b_adj`` (Wilson at
    1 - 0.05/B), ``U_b`` (GIAP's upper half-width rescaled to 1 - 0.025/B, capped at 1;
    descriptive) and the relations ``relation_95`` (descriptive) and ``relation_adj`` (the
    verdict's: the unrescaled genome upper limit below ``L_b_adj``'s lower limit, the adopted
    fallback)."""
    from scipy.stats import norm

    B = n_reportable
    p, hi = float(genome["precision"]), float(genome["ci_high"])
    lb = share(k_b, n_b)
    lo_adj, hi_adj = wilson_ci(k_b, n_b, level=1.0 - 0.05 / B)
    ratio = float(norm.ppf(1.0 - 0.025 / B) / norm.ppf(1.0 - (1.0 - H4_LEVEL) / 2.0))
    u = min(1.0, p + (hi - p) * ratio)
    return {
        "L_b": lb,
        "L_b_adj": lb | {"ci_low": lo_adj, "ci_high": hi_adj, "level": 1.0 - 0.05 / B},
        "U_b": u,
        "relation_95": "below_lower_bound" if hi < lb["ci_low"] else "consistent",
        "relation_adj": "below_lower_bound" if hi < lo_adj else "consistent",
    }


def _two_consecutive_below(rep: Sequence[Mapping[str, Any]]) -> bool:
    """Definition 8: two consecutive reportable bins (``rep`` in ``score_bin`` order, the
    non-reportable bins removed, so skipped labels do not break adjacency) both have
    ``relation_adj = below_lower_bound``. False when B < 2."""
    below = [x["relation_adj"] == "below_lower_bound" for x in rep]
    return any(a and b for a, b in pairwise(below))


def hypothesis4(
    queries: pd.DataFrame,
    d2: Collection[str],
    mapping_fdr: pd.DataFrame,
    *,
    searchable: Collection[str] | None = None,
    n_denovo_logged: int | None = None,
    n_denovo: int | None = None,
    a2_contaminant_rows: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """H4 (§15.7, definitions 2-10): per score bin of ``queries`` (``score_bins``, the step's
    own function) n_b, k_b (I→L folded on both sides), ``L_b`` and the genome row of
    ``mapping_fdr.tsv``; the relations of reportable bins (precision reported, n_b >= 20), the
    verdict, the descriptive extras (both over reportable bins), ``nested_share``,
    ``L_b_searchable`` and the caveats. Raises ``Inconsistent`` when the bin labels differ from
    ``mapping_fdr``'s or ``n_denovo_logged`` differs from ``n_denovo`` (the de novo queries
    before the length-class filter; default ``len(queries)``). Unscored (NaN) queries form no
    bin, as in the mapping step, but count in ``n_denovo``."""
    from giap.mapping.denovo_precision import score_bins

    peps = [_fold(p) for p in queries["peptide"].astype(str)]
    bins = score_bins(pd.to_numeric(queries["score"], errors="coerce").to_numpy(dtype=float))
    dn = mapping_fdr[mapping_fdr["section"].astype(str) == "denovo_precision"]
    recs = [r for r in dn.to_dict("records") if clean(r["score_bin"]) is not None]
    rows = {int(r["score_bin"]): r for r in recs}
    labels = sorted({int(b) for b in bins if b >= 0})
    if set(labels) != set(rows) or len(rows) != len(dn):
        raise Inconsistent(
            f"H4: score bins {labels} != mapping_fdr.tsv denovo_precision bins {sorted(rows)}"
        )
    warnings: list[Any] = []
    if n_denovo_logged is None:
        warnings.append({"code": "denovo_count_unchecked"})
    elif n_denovo_logged != (n_dn := len(peps) if n_denovo is None else n_denovo):
        raise Inconsistent(f"H4: {n_dn} de novo queries != placement log {n_denovo_logged}")
    D = {_fold(str(p)) for p in d2}
    if not D:
        warnings.append({"code": "a2_empty"})
    srch = None if searchable is None else {_fold(str(p)) for p in searchable}
    nested = _nested(peps)
    per: list[dict[str, Any]] = []
    for b in labels:
        S_b = [p for p, x in zip(peps, bins, strict=True) if x == b]
        r = rows[b]
        genome = {k: clean(r.get(k)) for k in ("precision", "ci_low", "ci_high", "c_value")}
        genome["flags"] = items(r.get("flags"))
        k_b = sum(p in D for p in S_b)
        if srch is None:
            searchable_b: dict[str, Any] = {**share(0, 0), "reason": "not_evaluable"}
        else:
            on = [p for p in S_b if p in srch]
            searchable_b = share(sum(p in D for p in on), len(on)) | {"label": "descriptive"}
        per.append({"score_bin": b, "n_b": len(S_b), "k_b": k_b, "genome": genome,
                    "reportable": genome["precision"] is not None and len(S_b) >= H4_MIN_N,
                    "nested_share": share(sum(p in nested for p in S_b), len(S_b)),
                    "L_b_searchable": searchable_b})  # fmt: skip
    rep = [x for x in per if x["reportable"]]
    for x in per:
        if x["reportable"]:
            x |= h4_bin(x["n_b"], x["k_b"], x["genome"], len(rep))
        else:
            x |= {"L_b": share(x["k_b"], x["n_b"]), "L_b_adj": None, "U_b": None,
                  "relation_95": None, "relation_adj": None}  # fmt: skip
    c_assumed = any("c_assumed" in items(r.get("flags")) for r in rows.values())
    if not rep or c_assumed:
        verdict = "not_evaluable"
    elif _two_consecutive_below(rep):
        verdict = "not_supported"
    else:
        verdict = "consistent_with_lower_bound"
    n_rep = sum(x["n_b"] for x in rep)
    mean = (
        {"value": sum(x["n_b"] * x["genome"]["precision"] for x in rep) / n_rep, "n": n_rep}
        if n_rep
        else {"value": None, "n": 0, "reason": "no_reportable_bin"}
    )
    caveats = list(H4_CAVEATS)
    if rep and all(x["genome"]["ci_high"] == 1.0 for x in rep):
        caveats.append("Every reportable bin has a genome CI clipped at 1: no refutation was "
                       "possible in those bins.")  # fmt: skip
    return {
        "verdict": verdict,
        "B": len(rep),
        "n_queries": len(peps),
        "c_assumed": c_assumed,
        "bins": per,
        "mean_genome_precision": mean | {"label": "descriptive"},
        "pooled_L": share(sum(x["k_b"] for x in rep), n_rep) | {"label": "descriptive"},
        "searchable_definition": SEARCHABLE_DIGEST if srch is not None else None,
        "a2_contaminant_rows": dict(a2_contaminant_rows) if a2_contaminant_rows else None,
        "caveats": caveats,
        "warnings": warnings,
    }


def contaminant_rows(P: pd.DataFrame, cfg: GiapConfig) -> dict[str, Any]:
    """S3: the share of A2's accepted target rows (q <= 0.01) that are contaminant rows, the
    rows ``c_by_class`` counts into c_value."""
    V, _ = validate_peptides(P, cfg.input)
    acc = V.loc[~V["is_decoy"].to_numpy(dtype=bool) & (_q(V) <= Q_ACCEPT)]
    return share(int(acc["is_contaminant"].astype(bool).sum()), len(acc))


def _bullet1_draft(codes: Mapping[str, Any]) -> str:
    """Bullet 1 as ``compute`` writes it (M10-05's strings; ``render --final`` relabels it)."""
    if all(codes.get(a) == 0 for a in ("A1", "A2", "A3")):
        return "pass"
    if codes.get("A1") == codes.get("A2") == 0 and codes.get("A3") is None:
        return "partial (A1, A2; A3 pending M10-10)"
    return "fail"


def bullet1(codes: Mapping[str, Any], a3_branch: str, feasibility: Mapping[str, Any]) -> str:
    """Bullet 1 at render time (S2): an arm's failure first, then the label from the runbook's
    ``a3_branch`` only (never guessed)."""
    if bad := [(a, codes.get(a)) for a in ("A1", "A2", "A3") if codes.get(a) not in (0, None)]:
        return f"fail ({bad[0][0]} exit {bad[0][1]})"
    if None in (codes.get("A1"), codes.get("A2")):
        return "fail (A1 or A2 not run)"
    if a3_branch not in A3_LABELS:
        return f"partial (A1, A2; A3 branch {a3_branch or 'absent'})"
    if codes.get("A3") is None:
        return f"fail (A3 not run: {feasibility.get('reason') or 'no reason given'})"
    if a3_branch == "cascadia" and feasibility.get("verdict") == "fail":
        return "fail (A3 feasibility fail)"
    return f"pass (A3 = {A3_LABELS[a3_branch]})"


def _bullet4(bench: Mapping[str, Any]) -> str:
    """Bullet 4: the 4 ``standard_run.wheat`` cells complete (non-null ``median_s``) and
    ADR-M10-D accepted; an incomplete cell wins over the ADR status."""
    if "reason" in bench:  # no --baseline
        return "pending (M10-12)"
    if bench["incomplete"]:
        return f"fail (incomplete: {', '.join(bench['incomplete'])})"
    status = bench.get("adr_status")
    if status in ("rejected", "superseded"):
        return f"fail (ADR-M10-D {status})"
    return "pass" if status in ADR_ACCEPTED else f"pending (ADR-M10-D {status or 'missing'})"


def _evidence_strings(metrics: Mapping[str, Any], bench: Mapping[str, Any]) -> list[str]:
    """One ``file SOURCE_SEP key`` string per checklist row; ``render --final`` appends the
    runbook's job IDs to rows 0-2 (``with_runbook``)."""
    if "reason" in bench:
        b4 = f"baseline.json{SOURCE_SEP}standard_run.wheat (not given, M10-12)"
    else:
        b4 = (
            f"baseline.json{SOURCE_SEP}standard_run.wheat.{{{','.join(CELLS)}}}.median_s; "
            f"{bench.get('adr') or 'ADR-M10-D'} {bench.get('adr_status') or 'missing'}; "
            "docs/user-guide/runtimes.md"
        )
    b5 = f"metrics.json{SOURCE_SEP}te_real.te2-te5.status" + "".join(
        f"; {t} not passed: report_notes.md{SOURCE_SEP}Diagnosis: {t}"
        for t in _not_passed(metrics)
        if t != "TE7"
    )
    return [
        f"metrics.json{SOURCE_SEP}evidence.exit_codes; reproduce.md{SOURCE_SEP}Evidence",
        f"metrics.json{SOURCE_SEP}te1.share_exact_inframe, te1.passed",
        f"metrics.json{SOURCE_SEP}m4, m5",
        b4,
        b5,
    ]


def checklist(metrics: Mapping[str, Any], run_states: Mapping[str, Any]) -> list[dict[str, str]]:
    """§15.6 checklist (spec table order) with one ``evidence`` string per row. ``run_states``:
    arm -> exit code (None when absent) and ``compare_ok`` (the three compare tables present and
    valid under their schemas). Bullet 4 reads ``metrics["evidence"]["benchmark"]`` (M10-15)."""
    bench = (metrics.get("evidence") or gate_evidence(None, None, None, {}))["benchmark"]
    te1 = metrics.get("te1") or {}
    te = "pass" if te1.get("passed") else "fail (no data)" if not te1.get("n_peptides") else "fail"
    m4 = metrics.get("m4") or {}
    cmp_ok = bool(run_states.get("compare_ok")) and "reason" not in m4
    statuses = [
        _bullet1_draft(run_states),
        te,
        "pass" if cmp_ok else "fail",
        _bullet4(bench),
        _bullet5(metrics),
    ]
    ev = _evidence_strings(metrics, bench)
    return [
        {"bullet": b, "status": s, "evidence": e}
        for b, s, e in zip(BULLETS, statuses, ev, strict=True)
    ]


# ---------------------------------------------------------------- M10-15 gate inputs

BULLETS = (
    "GIAP runs end to end on A1, A2 and A3",
    ">= 99.5 % of A2 peptides map exactly (TE1)",
    "A1 and A2 detections and abundances are compared",
    "Measured runtimes replace the estimates",
    "TE2-TE5 are reported; failures explained",
)
CELLS = ("k12_standard", "k12_core", "sim5m_standard", "sim5m_core")  # M10-12
ADR_ACCEPTED = ("accepted", "accepted-by-khoa")
A3_LABELS = {"cascadia": "Cascadia", "fallback": "Casanovo fallback on UTI89, different strain"}
EVIDENCE_STEPS = (*(f"S{n}" for n in range(8, 15)), "S-bench")
ROW_JOBS = {0: ("S8", "S9", "S11"), 1: ("S14",), 2: ("S12",)}  # checklist row -> runbook steps
ADR_STATUS = re.compile(r"^-?\s*Status:\s*([\w-]+)", re.M)
RUNBOOK_PENDING = "pending (render --final --runbook)"
RUNBOOK_STUB: dict[str, Any] = dict.fromkeys(("status", "freeze_commit", "a3_branch"))


def adr_status(text: str) -> str | None:
    """The first word of an ADR's ``Status:`` line (``accepted (lead)`` -> ``accepted``)."""
    m = ADR_STATUS.search(text)
    return m.group(1).lower() if m else None


def gate_evidence(
    baseline: Mapping[str, Any] | None,
    adr: tuple[str, str | None] | None,
    feasibility: Mapping[str, Any] | None,
    exit_codes: Mapping[str, Any],
) -> dict[str, Any]:
    """The ``evidence`` key of ``compute``: the arms' exit codes, the benchmark cells of
    ``baseline.json`` with ADR-M10-D (name, status) and M10-08's verdict; an absent input leaves
    a ``pending`` reason. The runbook is never recorded here (B1): ``runbook`` is a null stub
    that ``render --final --runbook`` fills."""
    wheat = ((baseline or {}).get("standard_run") or {}).get("wheat") or {}
    cells: dict[str, Mapping[str, Any]] = {
        c: e if isinstance(e := wheat.get(c), Mapping) else {} for c in CELLS
    }
    bench: dict[str, Any] = {
        "cells": {
            c: {k: e.get(k) for k in ("median_s", "commit", "cpu")} for c, e in cells.items()
        },
        "incomplete": [c for c, e in cells.items() if e.get("median_s") is None],
        "adr": None if adr is None else adr[0],
        "adr_status": None if adr is None else adr[1],
    }
    if baseline is None:
        bench["reason"] = "pending (M10-12)"
    feas = {k: (feasibility or {}).get(k) for k in ("verdict", "reason")}
    if feasibility is None:
        feas["reason"] = "pending (M10-08 verdict)"
    codes = {a: exit_codes.get(a) for a in ("A1", "A2", "A3")}
    return {"exit_codes": codes, "a3_feasibility": feas, "benchmark": bench,
            "runbook": RUNBOOK_STUB | {"reason": RUNBOOK_PENDING}}  # fmt: skip


def accepted_evidence(runbook: Any) -> dict[str, dict[str, str]]:
    """The runbook's Evidence rows S8-S14 and S-bench by step, superseded attempts skipped
    (B2); the latest accepted row of a step wins."""
    return {
        r["step"]: r
        for r in runbook.evidence
        if r.get("step") in EVIDENCE_STEPS and r.get("state") != "superseded"
    }


def with_runbook(metrics: Mapping[str, Any], runbook: Any | None) -> dict[str, Any]:
    """A copy of ``metrics`` with ``evidence.runbook``/``evidence.steps`` from the parsed
    runbook (a ``pending`` stub without one) and, with one, bullet 1 relabelled and the job IDs
    appended to the evidence strings of rows 0-2 (B1)."""
    doc = copy.deepcopy(dict(metrics))
    ev = doc.setdefault("evidence", {})
    if runbook is None:
        ev.setdefault("runbook", RUNBOOK_STUB | {"reason": RUNBOOK_PENDING})
        ev.setdefault("steps", {})
        return doc
    ev["runbook"] = {k: getattr(runbook, k) for k in ("status", "freeze_commit", "a3_branch")}
    steps = ev["steps"] = accepted_evidence(runbook)
    rows = doc.get("checklist")
    if isinstance(rows, list) and rows:
        codes, feas = ev.get("exit_codes") or {}, ev.get("a3_feasibility") or {}
        rows[0]["status"] = bullet1(codes, runbook.a3_branch, feas)
        for i, ids in ROW_JOBS.items():
            jobs = ", ".join(
                f"{s} {steps[s].get('job')} ({steps[s].get('state')}, exit "
                f"{steps[s].get('exit')}, {str(steps[s].get('commit'))[:7]})"
                for s in ids
                if s in steps
            )
            if jobs and i < len(rows) and str(rows[i].get("evidence") or "").strip():
                rows[i]["evidence"] += f"; reproduce.md{SOURCE_SEP}{jobs}"
    return doc


def pending_placeholders() -> dict[str, Any]:
    """``m6``, ``h4`` (M10-10) and ``m8`` (M10-07) before their hooks exist: every path of the
    final template resolves; the ``pending`` reason makes ``render --final`` refuse them."""
    level = {"fdp_combined": None, "fdp_lower": None, "assessment": None}
    return {
        "m6": {"overlap_a2_a3": {"novel_orfs": {"called": {"jaccard": None}},
                                 "alt_sites": {"high": {"jaccard": None}}},
               "reason": "pending (M10-10)"},
        "h4": {"verdict": None, "reason": "pending (M10-10)"},
        "m8": {"peptide": dict(level), "precursor": dict(level), "reason": "pending (M10-07)"},
    }  # fmt: skip


# ---------------------------------------------------------------- M10-06 hook (TE2-TE5, TE7)

TE_PENDING = "pending (M10-06)"
TE_KEYS = ("te2", "te3", "te4", "te5")
CAUSE_CLASSES = ("pipeline_bug", "annotation", "search", "biology", "definition")


def te_placeholders() -> dict[str, Any]:
    """``te_real`` and ``te7`` when ``compute`` has no ``--te-real`` / ``--te7`` (the template's
    paths resolve; bullet 5 stays pending)."""
    p = {"status": TE_PENDING}
    pairs = {k: {"log2_ratio": None} for k in ("rplL/rplJ", "atpA/atpD")}
    return {
        "te_real": {"te2": p | {"counts": {"violations": None}},
                    "te3": p | {"detection": {"share": None},
                                "abundance": {"percentile_rank": None}},
                    "te4": p | {"pairs": pairs},
                    "te5": p | {"A1": {"share_bucket": None}, "A2": {"share_bucket": None}}},
        "te7": {"giap": p | {"n_gate": None}, "diann": dict(p), "end_to_end": None},
    }  # fmt: skip


def _not_passed(metrics: Mapping[str, Any]) -> list[str]:
    """TEs (``TE3``) with a computed status other than ``pass``; pending ones are not listed."""
    te, te7 = metrics.get("te_real") or {}, metrics.get("te7") or {}
    out = [f"TE{k[2:]}" for k in TE_KEYS if "schema" in te and te[k].get("status") != "pass"]
    return out + (["TE7"] if "schema" in te7 and te7["giap"].get("status") != "pass" else [])


def _bullet5(metrics: Mapping[str, Any]) -> str:
    if "schema" not in (metrics.get("te_real") or {}):
        return TE_PENDING
    k = len([t for t in _not_passed(metrics) if t != "TE7"])
    return "pass" if k == 0 else f"reported ({k} not passed, explained)"


def diagnosis_problems(metrics: Mapping[str, Any], notes: str | None) -> list[str]:
    """Each not-passed TE needs a ``## Diagnosis: TE<n>`` notes section whose ``Cause class:``
    line names exactly one cause class."""
    out = []
    for te in _not_passed(metrics):
        m = re.search(rf"^## Diagnosis: {te}\b[^\n]*\n?(.*?)(?=^## |\Z)", notes or "", re.M | re.S)
        lines = re.findall(r"^Cause class:(.*)$", m.group(1), re.M) if m else []
        named = {c for ln in lines for c in re.findall(rf"\b({'|'.join(CAUSE_CLASSES)})\b", ln)}
        if m is None:
            out.append(f"{te}: no '## Diagnosis: {te}' section in the notes")
        elif len(named) != 1:
            out.append(f"{te}: the 'Cause class:' line names not one of {', '.join(CAUSE_CLASSES)}")
    return out


# ---------------------------------------------------------------- render

PLACEHOLDER = re.compile(r"\{\{\s*([^}\s]+)\s*\}\}")
NOTES_MARK = "<!-- notes -->"
M3_MARK = "<!-- m3 rows -->"
NOT_EVALUATED_MARK = "<!-- not evaluated -->"
DEFECT_REASONS = ("missing:", "pending")  # a defect in the final report, never a finding
M6_FINAL_ARMS = ("A2", "A3", "overlap_a2_a3")  # a ``missing:`` table here refuses --final (S3)
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
_MISSING = object()


class Gap(str):
    """An unresolved path under a metric that carries a ``reason``: rendered as a marker."""


def _walk(doc: Any, path: str) -> tuple[Any, Any]:
    """(value or ``_MISSING``, the nearest ``reason`` of the dicts walked through)."""
    cur, reason = doc, None
    for part in path.split("."):
        if isinstance(cur, Mapping):
            reason = cur.get("reason", reason)
        if isinstance(cur, Mapping) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list) and part.isdigit() and int(part) < len(cur):
            cur = cur[int(part)]
        else:
            return _MISSING, reason
    return cur, reason


def resolve(doc: Any, path: str) -> Any:
    """The value at a dotted path; ``Gap(reason)`` when it stops under a dict with a
    ``reason`` (the nearest one), else ``_MISSING``."""
    cur, reason = _walk(doc, path)
    if cur is _MISSING and reason is not None:
        return Gap(reason)
    return cur


def _fine(v: float) -> str:
    """4 significant figures, more when that would round a value below 1 up to "1"."""
    for d in range(4, 10):
        out = f"{v:.{d}g}"
        if v >= 1 or float(out) < 1:
            break
    return out


def fmt_value(key: str, v: Any) -> str:
    """Formats by key name: ``p_*`` 2 significant digits, ``*_s`` h:mm:ss, ``*_bytes`` GiB with
    1 decimal, shares (``value``, ``share*``, ``ci_*``) 4 significant figures; other floats 3;
    a share dict as ``v [lo, hi]``, another dict as ``k: v`` pairs (``none`` when empty)."""
    if isinstance(v, Mapping) and "value" in v:
        if v["value"] is None:
            return f"NA ({v.get('reason', 'no value')})"
        lo, hi = v.get("ci_low"), v.get("ci_high")
        ci = "" if lo is None else f" [{fmt_value('ci_low', lo)}, {fmt_value('ci_high', hi)}]"
        return fmt_value("value", v["value"]) + ci
    if isinstance(v, Mapping):
        return "; ".join(f"{k}: {fmt_value(str(k), x)}" for k, x in v.items()) or "none"
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
    if key == "value" or key.startswith(("share", "ci_")):
        return _fine(float(v))
    return f"{v:.3g}" if abs(v) < 1000 else f"{v:.1f}"


def _numbers(text: str) -> set[str]:
    """Digit runs with grouping commas removed ("9,999" -> "9999")."""
    return {n.replace(",", "") for n in re.findall(r"\d[\d,]*", text)}


def check_notes(notes: str, metrics: Mapping[str, Any]) -> None:
    """Raise ``Inconsistent`` for a peptide-like string or a 4+ digit number not in metrics."""
    if m := PEPTIDE_LIKE.search(notes):
        raise Inconsistent(f"notes: peptide-like string at offset {m.start()}")
    known = _numbers(json.dumps(metrics, sort_keys=True))
    for num in sorted(_numbers(notes)):
        if len(num) >= 4 and num not in known:
            raise Inconsistent(f"notes: number {num} is not in metrics.json")


def expand(template: str, metrics: Mapping[str, Any]) -> str:
    """The template with the metric-3 rows generated (one per subset of ``m3.<arm>.order``)."""
    m3 = metrics.get("m3") or {}
    names = dict.fromkeys(n for a in ARMS for n in (m3.get(a) or {}).get("order", []))
    rows = [f"| {n} | " + " | ".join(f"{{{{m3.{a}.subsets.{n}.share}}}}" for a in ARMS) + " |"
            for n in names]  # fmt: skip
    text = template.replace(M3_MARK, "\n".join(rows))
    for mark, lines in _m6_h4_rows(metrics).items():
        text = text.replace(mark, "\n".join(lines))
    return text


M6_ARMS = ("A1", "A2", "A3")
OVERLAP_SETS = (("novel_orfs", "called", "novel ORFs, called"),
                ("novel_orfs", "any_row", "novel ORFs, any row"),
                ("alt_sites", "alt", "ALT sites"),
                ("alt_sites", "high", "ALT sites, high confidence"))  # fmt: skip
OVERLAP_FIELDS = ("n_a2", "n_a3", "n_shared", "union", "jaccard", "a3_confirmed_by_a2",
                  "a2_found_by_a3")  # fmt: skip
H4_ROW = ("score_bin", "n_b", "L_b", "genome.precision", "genome.ci_low", "genome.ci_high",
          "genome.c_value", "U_b", "L_b_adj", "relation_95", "relation_adj", "nested_share",
          "L_b_searchable")  # fmt: skip


def _row(label: str, paths: Iterable[str]) -> str:
    return f"| {label} | " + " | ".join(f"{{{{{p}}}}}" for p in paths) + " |"


def _m6_computed(m6: Mapping[str, Any], table: str) -> tuple[list[str], list[str]]:
    """(arms whose ``table`` summary was computed, a reason placeholder per other arm): an arm
    or table with a reason is cited by its ``reason``, never by an unresolvable path."""
    done, why = [], []
    for a in M6_ARMS:
        arm = m6.get(a)
        if not isinstance(arm, Mapping) or "reason" in arm:
            why.append(f"{a} {{{{m6.{a}.reason}}}}" if isinstance(arm, Mapping) else "")
        elif isinstance(arm.get(table), Mapping) and "reason" in arm[table]:
            why.append(f"{a} {{{{m6.{a}.{table}.reason}}}}")
        else:
            done.append(a)
    return done, [w for w in why if w]


def _m6_h4_rows(metrics: Mapping[str, Any]) -> dict[str, list[str]]:
    """The generated lines of the metric-6 and H4 tables (per-arm columns for the arms whose
    tables were computed, one H4 row per bin, one line per caveat), keyed by template marker."""
    m6 = metrics.get("m6") if isinstance(metrics.get("m6"), Mapping) else {}
    arms, why = _m6_computed(m6, "novel_orfs")
    bk = [(a, b) for a in arms for b in BUCKETS]
    cls = [f"| class | {' | '.join(f'{a} {b}' for a, b in bk)} |", "|---" * (len(bk) + 1) + "|"]
    cls += [_row(c, (f"m6.{a}.novel_orfs.by_class.{c}.{b}" for a, b in bk)) for c in NOVEL_CLASSES]
    cls.append(_row("total", (f"m6.{a}.novel_orfs.totals.{b}" for a, b in bk)))
    if not arms:
        cls = []
    cls += [f"Not computed: {'; '.join(why)}."] if why else []
    arms, why = _m6_computed(m6, "variants")
    tc = [(a, c) for a in arms for c in CONFIDENCE]
    alt = [f"| ALT type | {' | '.join(f'{a} {c}' for a, c in tc)} |", "|---" * (len(tc) + 1) + "|"]
    alt += [_row(t, (f"m6.{a}.variants.alt.{t}.{c}" for a, c in tc)) for t in ALT_TYPES]
    alt.append("")
    alt += [f"{a}: mixed calls {{{{m6.{a}.variants.mixed}}}}; class error rate numeric "
            f"{{{{m6.{a}.variants.class_error.numeric}}}}, unestimated "
            f"{{{{m6.{a}.variants.class_error.unestimated}}}}." for a in arms]  # fmt: skip
    if not arms:
        alt = []
    alt += [f"Not computed: {'; '.join(why)}."] if why else []
    ov = [_row(label, (f"m6.overlap_a2_a3.{f}.{s}.{x}" for x in OVERLAP_FIELDS))
          for f, s, label in OVERLAP_SETS]  # fmt: skip
    h4 = metrics.get("h4") or {}
    bins = h4.get("bins") if isinstance(h4.get("bins"), list) else []
    caveats = h4.get("caveats") if isinstance(h4.get("caveats"), list) else []
    return {
        "<!-- m6 class rows -->": cls,
        "<!-- m6 alt rows -->": alt,
        "<!-- m6 overlap rows -->": ov,
        "<!-- h4 rows -->": [
            _row(f"{{{{h4.bins.{i}.score_bin}}}}", (f"h4.bins.{i}.{p}" for p in H4_ROW[1:]))
            for i in range(len(bins))
        ],
        "<!-- h4 caveats -->": [f"- {{{{h4.caveats.{i}}}}}" for i in range(len(caveats))],
    }


SECTION = re.compile(r"^<!-- (if|ifnot) ([\w.]+) -->\n(.*?)^<!-- end \2 -->\n", re.M | re.S)


def sections(template: str, metrics: Mapping[str, Any]) -> str:
    """Keeps ``<!-- if PATH -->`` blocks when PATH resolves, ``<!-- ifnot PATH -->`` blocks
    when it does not (M10-07: metric 8 only with ``--entrapment``); drops the markers."""

    def keep(m: re.Match[str]) -> str:
        v = resolve(metrics, m.group(2))
        has = v is not _MISSING and not isinstance(v, Gap)
        return m.group(3) if has == (m.group(1) == "if") else ""

    return SECTION.sub(keep, template)


def missing_paths(template: str, metrics: Mapping[str, Any]) -> list[str]:
    paths = PLACEHOLDER.findall(sections(expand(template, metrics), metrics))
    return [p for p in paths if resolve(metrics, p) is _MISSING]


def not_evaluated(metrics: Mapping[str, Any]) -> list[str]:
    """One line per metric or hypothesis (and its sub-blocks, to depth three) carrying a
    ``reason``."""
    out = []
    for key in sorted(k for k in metrics if re.fullmatch(r"m\d|h\d|te1", k)):
        blocks = [(key, metrics[key])]
        if isinstance(metrics[key], Mapping):
            for k, v in sorted(metrics[key].items()):
                blocks.append((f"{key}.{k}", v))
                if isinstance(v, Mapping):
                    blocks += [(f"{key}.{k}.{k2}", v2) for k2, v2 in sorted(v.items())]
        for path, v in blocks:
            if isinstance(v, Mapping) and isinstance(r := v.get("reason"), str):
                finding = not r.startswith(DEFECT_REASONS)
                out.append(f"- `{path}`: " + (f"not evaluable ({r})" if finding else r))
    return out or ["- none: every metric and hypothesis has a value."]


def final_problems(metrics: Mapping[str, Any]) -> list[str]:
    """What the final report refuses: a ``pending``/``partial`` checklist row or an empty
    evidence string, a ``pending`` or absent TE2-TE5 or TE7 status, a ``missing:`` table of
    metric 6 under A2, A3 or their overlap (S3; A1's stays a value) (each message names it)."""
    out = []
    for i, row in enumerate(metrics["checklist"]):
        if str(row.get("status")).startswith(("pending", "partial")):
            out.append(f"checklist.{i}.status is {row.get('status')!r}")
        if not str(row.get("evidence") or "").strip():
            out.append(f"checklist.{i}.evidence is empty")
    te, te7 = metrics.get("te_real") or {}, metrics.get("te7") or {}
    statuses = {f"te_real.{k}.status": (te.get(k) or {}).get("status") for k in TE_KEYS}
    statuses["te7.giap.status"] = (te7.get("giap") or {}).get("status")
    out += [f"{p} is {s!r}" for p, s in statuses.items() if s is None or s.startswith("pending")]
    m6 = metrics.get("m6")
    for arm in M6_FINAL_ARMS if isinstance(m6, Mapping) else ():
        block = m6.get(arm)
        for table, v in sorted(block.items()) if isinstance(block, Mapping) else ():
            if isinstance(v, Mapping) and str(v.get("reason", "")).startswith("missing:"):
                out.append(f"m6.{arm}.{table} is {v['reason']!r}")
    return out


def render(
    template: str, metrics: Mapping[str, Any], notes: str | None, final: bool = False
) -> str:
    """The report: each ``{{path}}`` becomes its formatted value and its source
    ``(metrics.json SOURCE_SEP parent)``; a path under a metric with a reason becomes
    ``**missing: <reason>**``, any other unresolved path ``**missing: <path>**``. ``final``
    (M10-15): a null or absent value under a finding reason (not ``missing:``/``pending``)
    renders as ``not evaluable (<reason>)`` (so does a mistyped path under such a block), and a
    null one under a defect reason as ``**missing: <reason>**``. In draft mode a null under a
    ``pending`` reason renders that reason (N4). Values of ``evidence.runbook`` and
    ``evidence.steps`` cite ``reproduce.md`` (``with_runbook`` fills them)."""
    if notes is not None:
        check_notes(notes, metrics)

    def sub(m: re.Match[str]) -> str:
        path = m.group(1)
        v, reason = _walk(metrics, path)
        if final and isinstance(reason, str) and (v is None or v is _MISSING):
            ok = not reason.startswith(DEFECT_REASONS)
            return f"not evaluable ({reason})" if ok else f"**missing: {reason}**"
        if v is None and isinstance(reason, str) and reason.startswith("pending"):
            return reason
        v = resolve(metrics, path)
        if v is _MISSING or isinstance(v, Gap):
            return f"**missing: {path if v is _MISSING else v}**"
        parent = path.rsplit(".", 1)[0]
        src = "reproduce.md" if path.startswith(("evidence.runbook", "evidence.steps")) else ""
        return (
            f"{fmt_value(path.rsplit('.', 1)[-1], v)} ({src or 'metrics.json'}{SOURCE_SEP}{parent})"
        )

    text = PLACEHOLDER.sub(sub, sections(expand(template, metrics), metrics))
    text = text.replace(NOT_EVALUATED_MARK, "\n".join(not_evaluated(metrics)))
    body = notes.strip() if notes and notes.strip() else "No notes."
    if NOTES_MARK in text:
        return text.replace(NOTES_MARK, body)
    return text + ("" if notes is None else f"\n## Notes\n\n{body}\n")


# ---------------------------------------------------------------- CLI (file I/O)

#: hashed input files (argparse dests); ``cmp``, ``a1``, ``a2``, ``a3`` are directories
#: M10-06 inputs: ``te_real.json`` and ``te7.json`` (optional; absent -> ``te_placeholders()``)
TE_INPUTS = ("te_real", "te7")
#: M10-15 inputs (committed, read-only): ``benchmarks/baseline.json``, ADR-M10-D and
#: ``a3_feasibility.json``; the runbook is an input of ``render --final`` only (B1)
GATE_INPUTS = ("baseline", "benchmark_adr", "a3_feasibility")
#: M10-10b inputs: ``--a3`` needs all of ``A3_COMPANIONS``; ``--strain-db`` is optional (S4)
A3_COMPANIONS = ("peptides_a3", "convert_a3", "cmp3")
A3_INPUTS = ("peptides_a3", "convert_a3", "strain_db")
INPUTS = ("peptides_a1", "peptides_a2", "convert_a1", "convert_a2", "summary_a1", "summary_a2",
          "bakta_json", "genome", "k12_fasta", "sacct", *TE_INPUTS, *GATE_INPUTS,
          "entrapment", *A3_INPUTS)  # fmt: skip
OPTIONAL = ("summary_a1", "summary_a2", "sacct", "a3", *TE_INPUTS, *GATE_INPUTS, "entrapment",
            *A3_INPUTS, "cmp3")  # fmt: skip
REQUIRED = ("a1", "a2", "cmp", *(k for k in INPUTS if k not in OPTIONAL), "out")
COMPARE = ("compare_detection", "compare_abundance", "compare_stability")
MAP_TABLES = "peptide_hits.tsv/unmapped.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _need(path: Path) -> Path:
    if not path.exists():
        raise Precondition(f"missing input: {path}")
    return path


def _load(path: Path, read: Callable[[Path], T]) -> T:
    """``read(path)``; a missing file or a read, parse, key or schema error is a precondition."""
    try:
        return read(_need(path))
    except (ValueError, KeyError, OSError) as err:
        raise Precondition(f"{path}: {err}") from err


def _json(path: Path) -> dict[str, Any]:
    doc: dict[str, Any] = _load(path, lambda p: json.loads(p.read_text(encoding="utf-8")))
    return doc


def table_schema(stem: str) -> dict[str, Any]:
    res = resources.files("giap.schemas").joinpath("tables", f"{stem}.schema.json")
    doc: dict[str, Any] = json.loads(res.read_text(encoding="utf-8"))
    return doc


def _table(path: Path, stem: str) -> pd.DataFrame:
    # read_schema_tsv keeps the written 1-based coordinates, which te1_result expects;
    # RunView.table would convert them to 0-based.
    from giap.io.run_view import read_schema_tsv

    return _load(path, lambda p: read_schema_tsv(p, stem, table_schema(stem)))


def _rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    text = path.read_text(encoding="utf-8")
    delim = "|" if "|" in text.split("\n", 1)[0] else "\t"
    return list(csv.DictReader(text.splitlines(), delimiter=delim))


def _run_table(out: Path, name: str) -> Path | None:
    """GIAP writes ``mapping/{peptide_hits,unmapped}.tsv`` and ``tables/detection_summary.tsv``."""
    return next((f for d in ("tables", "mapping") if (f := out / d / name).is_file()), None)


def _input_sha(manifest: Mapping[str, Any], role: str) -> str | None:
    return next((i.get("sha256") for i in manifest.get("inputs", []) if i.get("role") == role),
                None)  # fmt: skip


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


def _cfg(manifest: Mapping[str, Any], path: Path) -> GiapConfig:
    if not manifest.get("config"):
        raise Precondition(f"{path}: no resolved config")
    try:
        return GiapConfig.model_validate(manifest["config"])
    except ValueError as err:
        raise Precondition(f"{path}: invalid config: {err}") from err


def _spec_sha(spec: str = SPEC) -> str:
    """The commit of ``spec`` (git log); git unavailable or no commit is a ``Precondition``."""
    try:
        res = subprocess.run(
            ["git", "-C", str(REPO), "log", "-1", "--format=%H", "--", spec],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as err:
        raise Precondition(f"git log {spec}: {err}") from err
    if not re.fullmatch(r"[0-9a-f]{40}", sha := res.stdout.strip()):
        raise Precondition(f"git log {spec}: no commit found")
    return sha


def _git_ignored(path: Path) -> bool:
    """``git check-ignore``: ``path`` (existing or not) is ignored in the checkout."""
    try:
        res = subprocess.run(["git", "-C", str(REPO), "check-ignore", "-q", str(path)],
                             capture_output=True, check=False)  # fmt: skip
    except OSError:
        return False
    return res.returncode == 0


def _outside_repo(path: Path, flag: str, allow_ignored: bool = False) -> None:
    """``path`` lies outside the checkout, or (``allow_ignored``, M10-07 B2) in a git-ignored
    part of it such as ``local/``, where farmshare_env.sh puts ``GIAP_ROOT``."""
    p = path.resolve()
    if (p == REPO or REPO in p.parents) and not (allow_ignored and _git_ignored(p)):
        raise Precondition(f"{flag} {path} is inside the repository")


def _arm(label: str, out: Path, a: argparse.Namespace, genome_sha: str, job: str | None,
         sacct: list[dict[str, str]]) -> dict[str, Any]:  # fmt: skip
    """Reads one arm and checks its preconditions."""
    from giap.io.compare_io import genome_digest
    from giap.io.peptides import read_peptide_table

    x = label.lower()
    pep, summary = getattr(a, f"peptides_{x}"), getattr(a, f"summary_{x}")
    mp = out / "run_manifest.json"
    man = _json(mp)
    if man.get("exit_code") != 0:
        raise Precondition(f"{mp}: exit_code {man.get('exit_code')}")
    if genome_digest(out) != genome_sha:
        raise Precondition(f"{mp}: genome sha256 differs from --genome")
    if _input_sha(man, "peptides") != sha256(_need(pep)):
        raise Precondition(f"{pep}: sha256 differs from the run's peptides input ({out})")
    sm = None if summary is None else _json(summary)
    if sm is not None and not isinstance(sm.get("distinct_stripped_q01"), int):
        raise Precondition(f"{summary}: no integer distinct_stripped_q01")
    t: dict[str, pd.DataFrame | None] = {}
    for n in ("peptide_hits", "unmapped", "detection_summary"):
        f = _run_table(out, f"{n}.tsv")
        t[n] = None if f is None else _table(f, n)
    return {
        "label": label,
        "out": out,
        "manifest": man,
        "cfg": _cfg(man, mp),
        "P": _load(pep, lambda p: read_peptide_table(p)[0]),
        "genes": _table(out / "tables" / "genes.tsv", "genes"),
        "hits": t["peptide_hits"],
        "unm": t["unmapped"],
        "ds": t["detection_summary"],
        "conv": _json(getattr(a, f"convert_{x}")),
        "sum": sm,
        "runtime": runtime(man, _rows(out / "logs" / "trace.tsv"),
                           [r for r in sacct if job is not None and r.get("JobID") == job]),
    }  # fmt: skip


def _missing(name: str) -> dict[str, str]:
    return {"reason": f"missing:{name}"}


def _or_missing(name: str, fn: Callable[..., Any], *args: Any) -> Any:
    """``fn(*args)``, or the ``missing:<name>`` stub when a table argument is absent."""
    return _missing(name) if any(x is None for x in args) else fn(*args)


def _cross_checks(arm: Mapping[str, Any], m1: Mapping[str, Any], warnings: list[Any]) -> None:
    n = m1["union"]["peptides"]["all"]
    counts = arm["conv"].get("counts", {})
    if n != counts.get("distinct_peptides_target_q01"):
        raise Inconsistent(
            f"{arm['label']}: metric 1 peptides {n} != convert.json "
            f"distinct_peptides_target_q01 {counts.get('distinct_peptides_target_q01')}"
        )
    if arm["sum"] is not None:
        diff = n - arm["sum"]["distinct_stripped_q01"]
        qc = int(counts.get("q_crossings", 0))
        if not 0 <= diff <= qc:
            warnings.append({"code": "raw_vs_written_q01_out_of_range", "arm": arm["label"],
                             "raw_vs_written_q01": diff, "q_crossings": qc})  # fmt: skip


A3_NEEDS = "--a3 needs --peptides-a3, --convert-a3 and --cmp3"


def _a3_options(a: argparse.Namespace) -> None:
    """``--a3`` and its companions come together (example 13); ``--strain-db`` needs ``--a3``."""
    given = [k for k in (*A3_COMPANIONS, "strain_db") if getattr(a, k) is not None]
    if a.a3 is not None and any(getattr(a, k) is None for k in A3_COMPANIONS):
        raise Precondition(A3_NEEDS)
    if a.a3 is None and given:
        opts = ", ".join("--" + k.replace("_", "-") for k in given)
        raise Precondition(f"{opts} given without --a3 ({A3_NEEDS})")


def _arm_tables(out: Path) -> dict[str, pd.DataFrame | None]:
    """An arm's ``novel_orfs.tsv`` and ``variants.tsv`` (None when the module was off)."""
    t: dict[str, pd.DataFrame | None] = {}
    for stem in ("novel_orfs", "variants"):
        f = _run_table(out, f"{stem}.tsv")
        t[stem] = None if f is None else _table(f, stem)
    return t


def _placement_logged(out: Path) -> int | None:
    """The placement line's de novo count from the ``.command.log`` of the last COMPLETED
    placement task of ``OUT/logs/trace.tsv`` (B1); None without a trace, task dir or line."""
    from giap.runner.nextflow import _step_name, _task_dir

    trace = out / "logs" / "trace.tsv"
    if not trace.is_file():
        return None
    with trace.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    done = [
        r
        for r in rows
        if (r.get("status") or "").strip() == "COMPLETED"
        and _step_name(r.get("name") or "") == "placement"
    ]
    if not done:
        return None
    log = Path(_task_dir(done[-1], out / "work")) / ".command.log"
    if not log.is_file():
        return None
    return placement_denovo_count(log.read_text(encoding="utf-8", errors="replace"))


def _read_a3(a: argparse.Namespace, genome_sha: str) -> dict[str, Any]:
    """The A3 inputs and their preconditions (spec "metrics.py changes (10b)")."""
    from giap.io.compare_io import genome_digest
    from giap.io.peptides import read_peptide_table

    out = _need(resolve_attempt(a.a3)[0])
    mp = out / "run_manifest.json"
    man = _json(mp)
    if man.get("exit_code") != 0:
        raise Precondition(f"{mp}: exit_code {man.get('exit_code')}")
    if genome_digest(out) != genome_sha:
        raise Precondition(f"{mp}: genome sha256 differs from --genome")
    pep_sha = sha256(a.peptides_a3)
    if _input_sha(man, "peptides") != pep_sha:
        raise Precondition(f"{a.peptides_a3}: sha256 differs from the run's peptides input ({out})")
    if _input_sha(man, "reference_peptides") != sha256(a.peptides_a2):
        raise Precondition(f"{mp}: the reference_peptides input is not --peptides-a2")
    conv = _json(a.convert_a3)
    if conv.get("failures"):
        raise Precondition(f"{a.convert_a3}: failures {conv.get('failures')}")
    if ((conv.get("outputs") or {}).get("peptides") or {}).get("sha256") != pep_sha:
        raise Precondition(f"{a.convert_a3}: peptide table sha256 differs from --peptides-a3")
    cmp3: dict[str, pd.DataFrame] = {}
    for stem in COMPARE:
        cmp3[stem] = _table(a.cmp3 / "tables" / f"{stem}.tsv", stem)
    fdr = _run_table(out, "mapping_fdr.tsv")
    if fdr is None:
        raise Precondition(f"{out}: no mapping/mapping_fdr.tsv")
    return {"out": out, "manifest": man, "cfg": _cfg(man, mp),
            "P": _load(a.peptides_a3, lambda p: read_peptide_table(p)[0]),
            "tables": _arm_tables(out), "cmp3": cmp3,
            "mapping_fdr": _table(fdr, "mapping_fdr"),
            "logged": _placement_logged(out)}  # fmt: skip


def _searchable(path: Path | None) -> set[str] | None:
    """In-silico tryptic peptides of ``--strain-db`` (GIAP's digest, I→L), or None."""
    from giap.io.protein_fasta import read_protein_fasta

    if path is None:
        return None
    out: set[str] = set()
    for rec in _load(path, read_protein_fasta):
        out |= _digest_il(rec.sequence)
    return out


def _m6_h4(
    a: argparse.Namespace,
    arms: Mapping[str, Mapping[str, Any]],
    a3: Mapping[str, Any] | None,
    warnings: list[Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Metric 6 and H4 from the read inputs (``a3_not_run`` stubs without ``--a3``)."""
    tables = {x: _arm_tables(arms[x]["out"]) for x in ARMS}
    if a3 is None:
        m6 = metric6(tables, None)
        return m6, dict(A3_NOT_RUN)
    tables["A3"] = a3["tables"]
    m6 = metric6(tables, a3["cmp3"]["compare_stability"])
    a2 = arms["A2"]
    _, D2, _ = split_contaminants(a2["P"], a2["cfg"])
    h4 = hypothesis4(denovo_queries(a3["P"], a3["cfg"]), D2, a3["mapping_fdr"],
                     searchable=_searchable(a.strain_db), n_denovo_logged=a3["logged"],
                     n_denovo=denovo_count(a3["P"], a3["cfg"]),
                     a2_contaminant_rows=contaminant_rows(a2["P"], a2["cfg"]))  # fmt: skip
    if a3["manifest"].get("git_commit") != arms["A2"]["manifest"].get("git_commit"):
        warnings.append({"code": "giap_commit_differs", "A3": a3["manifest"].get("git_commit")})
    return m6, h4


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Every metric from the parsed arguments (reads only; raises on a precondition)."""
    from giap.io.output_checks import check_tsv
    from giap.io.protein_fasta import read_protein_fasta

    _a3_options(a)
    _outside_repo(a.out, "--out")
    for k in ("cmp", "cmp3", *INPUTS):
        if getattr(a, k) is not None:
            _need(getattr(a, k))
    genome_sha = sha256(a.genome)
    sacct = _rows(a.sacct)
    warnings: list[Any] = []
    arms = {}
    for label in ARMS:
        out_dir, job = resolve_attempt(getattr(a, label.lower()))
        arms[label] = _arm(label, _need(out_dir), a, genome_sha, job, sacct)
        if a.sacct is not None and arms[label]["runtime"]["sacct"] is None:
            code = "sacct_ignored_no_attempts" if job is None else "sacct_unmatched"
            warnings.append({"code": code, "arm": label})
    genes, circular = _load(a.bakta_json, lambda p: read_coding_genes(p, a.genome))
    a1, a2 = arms["A1"], arms["A2"]
    g2 = a2["genes"]
    tags = set(g2.loc[g2["feature_type"].astype(str).isin([*CODING, "pseudogene"]), "locus_tag"])
    if {g.locus_tag for g in genes} != {str(t) for t in tags}:
        raise Precondition(f"{a.bakta_json}: locus tags differ from {a2['out']}/tables/genes.tsv")
    k12 = _load(a.k12_fasta, lambda p: {r.id: r.sequence for r in read_protein_fasta(p)})
    classes = strain_specific_set({g.locus_tag: g.aa for g in genes}, k12)
    doc: dict[str, Any] = {"schema": SCHEMA, "spec": _spec_sha(), "warnings": warnings,
                           "giap_commit": a2["manifest"].get("git_commit")}  # fmt: skip
    if a1["manifest"].get("git_commit") != doc["giap_commit"]:
        warnings.append({"code": "giap_commit_differs", "A1": a1["manifest"].get("git_commit")})
    doc["m1"], doc["m2"], doc["m3"], doc["m7"] = {}, {}, {}, {}
    for label, arm in arms.items():
        doc["m1"][label] = identifications(arm["P"], arm["cfg"].input.contaminant_prefixes)
        _cross_checks(arm, doc["m1"][label], warnings)
        m2 = _or_missing(MAP_TABLES, mapping_shares, arm["P"], arm["hits"], arm["unm"], arm["cfg"])
        if "counts" in m2:
            gi = _or_missing("detection_summary.tsv", genome_inconsistent, arm["ds"])
            tu = m2["counts"]["tolerant"] + m2["counts"]["unmapped"]
            gi["consistent"] = None if gi.get("count") is None else gi["count"] == tu
            m2["genome_inconsistent"] = gi
        doc["m2"][label] = m2
        doc["m3"][label] = coverage(arm["genes"], classes)
        doc["m7"][label] = arm["runtime"]
    cmp: dict[str, pd.DataFrame | None] = {}
    compare_ok = True
    for stem in COMPARE:  # giap compare writes CMP/tables/<stem>.tsv
        f = a.cmp / "tables" / f"{stem}.tsv"
        cmp[stem] = _table(f, stem) if f.is_file() else None
        compare_ok &= cmp[stem] is not None
        try:  # enum, minimum, pattern, required: a fail and a warning, not an exit code
            if cmp[stem] is not None:
                check_tsv(stem, f, stem)
        except ValueError as err:
            compare_ok = False
            warnings.append({"code": "compare_table_invalid", "table": f"{stem}.tsv",
                             "error": str(err)})  # fmt: skip
    doc["m4"] = _or_missing("compare_abundance.tsv", abundance_agreement,
                            cmp["compare_abundance"], ARMS)  # fmt: skip
    doc["m5"] = _or_missing("compare_stability.tsv", stability, cmp["compare_stability"], ARMS)
    exceptions: list[dict[str, Any]] = []
    if a2["hits"] is None or a2["unm"] is None:
        doc["te1"] = _missing(MAP_TABLES)
    else:
        qc = int(a2["conv"].get("counts", {}).get("q_crossings", 0))
        doc["te1"], exceptions = te1_real(a2["P"], a2["hits"], a2["unm"], genes, circular,
                                          non_atg_nterm(genes), a2["cfg"],
                                          q_crossings=qc)  # fmt: skip
        for code in ("unknown_source", "target_q_na"):
            if doc["te1"][code]:
                warnings.append({"code": code, "arm": "A2", "n": doc["te1"][code]})
    doc["te1_a1_diagnostic"] = _or_missing("peptide_hits.tsv", te1_diagnostic, a1["P"],
                                           a1["hits"], a1["cfg"])  # fmt: skip
    doc["h1"] = hypothesis1(a1["genes"], g2, cmp["compare_detection"], gene_subsets(g2, classes))
    warnings += doc["h1"]["warnings"]
    doc["h3"] = hypothesis3(doc["m4"], doc["m5"])
    states: dict[str, Any] = {"A1": 0, "A2": 0, "A3": None, "compare_ok": compare_ok}
    a3 = None if a.a3 is None else _read_a3(a, genome_sha)
    if a3 is not None:
        states["A3"] = a3["manifest"].get("exit_code")
    doc |= te_placeholders() | pending_placeholders()
    doc["m6"], doc["h4"] = _m6_h4(a, arms, a3, warnings)
    warnings += doc["m6"]["warnings"] + doc["h4"].get("warnings", [])
    for key in TE_INPUTS:
        if getattr(a, key) is not None:
            doc[key] = _json(getattr(a, key))
    if a.entrapment is not None:  # M10-07 hook
        doc["m8"] = entrapment_block(_json(a.entrapment))
    doc["evidence"] = _gate_evidence(a, states)
    doc["checklist"] = checklist(doc, states)
    files = {k: getattr(a, k) for k in INPUTS} | {
        f"manifest_{x.lower()}": arms[x]["out"] / "run_manifest.json" for x in ARMS
    }
    if a3 is not None:
        files["manifest_a3"] = a3["out"] / "run_manifest.json"
    doc["inputs"] = {
        k: {"path": str(Path(p).resolve()), "sha256": sha256(Path(p))}
        for k, p in files.items()
        if p is not None
    }
    return doc, exceptions


M8_KEYS = ("peptide", "precursor", "universe", "caveats", "cross_checks", "warnings")


def entrapment_block(ent: Mapping[str, Any]) -> dict[str, Any]:
    """Metric 8: the blocks of ``entrapment.json`` (validation/pilot/entrapment.py, M10-07)."""
    if ent.get("schema") != "giap-pilot-entrapment/1" or any(k not in ent for k in M8_KEYS):
        raise Precondition(f"--entrapment: not a giap-pilot-entrapment/1 document with {M8_KEYS}")
    codes = sorted({str(w.get("code") if isinstance(w, Mapping) else w) for w in ent["warnings"]})
    return {k: ent[k] for k in M8_KEYS} | {"warning_codes": codes}


def _gate_evidence(a: argparse.Namespace, states: Mapping[str, Any]) -> dict[str, Any]:
    """Reads the M10-15 inputs for ``gate_evidence``."""
    adr = None
    if a.benchmark_adr is not None:
        text = _load(a.benchmark_adr, lambda q: q.read_text(encoding="utf-8"))
        adr = (a.benchmark_adr.stem, adr_status(text))
    baseline = None if a.baseline is None else _json(a.baseline)
    feas = None if a.a3_feasibility is None else _json(a.a3_feasibility)
    return gate_evidence(baseline, adr, feas, states)


def _dump(doc: Any) -> str:
    return json.dumps(doc, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write_outputs(
    out: Path, doc: Mapping[str, Any], exceptions: Sequence[Mapping[str, Any]]
) -> None:
    from giap.io.tsv import write_tsv

    text = _dump(doc)  # before any file is written (a NaN raises)
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(text, encoding="utf-8")
    (out / "te1.json").write_text(_dump(doc["te1"]), encoding="utf-8")
    cols = ["kind", "locus_tag", "peptide", "hits"]
    write_tsv(out / "te1_exceptions.tsv", cols, [[r.get(c) for c in cols] for r in exceptions])


def render_report(
    template: str,
    metrics: Mapping[str, Any],
    notes: str | None,
    final: bool = False,
    runbook: Any | None = None,
) -> str:
    """``render`` of ``with_runbook(metrics, runbook)`` plus its refusals (``Inconsistent``): an
    absent placeholder path, a not-passed TE without its diagnosis; ``final`` adds a runbook not
    ``status=final``, ``final_problems`` and any ``**missing:`` marker. A ``metrics.json``
    without ``checklist``, or no runbook, in final mode is a ``Precondition``."""
    if final and not isinstance(metrics.get("checklist"), list):
        raise Precondition("render --final: metrics.json has no checklist")
    if final and runbook is None:
        raise Precondition("render --final needs --runbook (bullet 1 and the job IDs)")
    doc = with_runbook(metrics, runbook)
    text = render(template, doc, notes, final)
    if missing := missing_paths(template, doc):
        raise Inconsistent(f"render: missing paths: {', '.join(missing)}")
    problems = diagnosis_problems(doc, notes)
    if final:
        if (status := getattr(runbook, "status", None)) != "final":
            problems.append(f"the runbook is status={status or 'absent'}, not final")
        problems += final_problems(doc)
        if marks := sorted(set(re.findall(r"\*\*missing: (.*?)\*\*", text))):
            problems.append("the report holds **missing: markers: " + ", ".join(marks))
    if problems:
        raise Inconsistent("render: " + "; ".join(problems))
    return text


def _render(a: argparse.Namespace) -> int:
    from validation.pilot.reproduce_check import parse as parse_runbook

    _outside_repo(a.out, "render --out")
    metrics = _json(a.metrics)
    template = _need(a.template).read_text(encoding="utf-8")
    notes = None if a.notes is None else _need(a.notes).read_text(encoding="utf-8")
    rb = None
    if a.runbook is not None:
        rb = parse_runbook(_load(a.runbook, lambda q: q.read_text(encoding="utf-8")))
    a.out.write_text(render_report(template, metrics, notes, a.final, rb), encoding="utf-8")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for dest in (*REQUIRED, *OPTIONAL):
        c.add_argument(f"--{dest.replace('_', '-')}", type=Path, required=dest in REQUIRED)
    r = sub.add_parser("render")
    r.add_argument("--metrics", type=Path, required=True)
    r.add_argument("--final", action="store_true", help="the gate report (M10-15)")
    r.add_argument("--runbook", type=Path, default=None, help="reproduce.md (required by --final)")
    r.add_argument("--notes", type=Path, default=None)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--template", type=Path, default=Path(__file__).with_name("report_template.md"))
    a = ap.parse_args(argv)
    try:
        if a.cmd == "render":
            return _render(a)
        doc, exc = compute(a)
        write_outputs(a.out, doc, exc)
        print(
            f"wrote {a.out}/metrics.json; TE1 passed={doc['te1'].get('passed')}; "
            f"{len(doc['warnings'])} warning(s)"
        )
        return 0
    except Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except Inconsistent as err:
        print(f"inconsistent: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
