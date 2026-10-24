"""Pilot metrics for the arms A1 and A2 (GIAP_SPEC.md §15.5 metrics 1-5 and 7, §11.6 TE1 on real
data, §15.6 checklist, §15.7 hypotheses 1 and 3); spec docs/specs/M10-05-pilot-metrics.md.

Pure functions (frames and dicts in, dicts out) plus a thin CLI that reads and writes::

    python validation/pilot/metrics.py compute --a1 OUT_A1 --a2 OUT_A2 --cmp CMP \\
        --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --convert-a1 CONV_A1 --convert-a2 CONV_A2 \\
        [--summary-a1 SUM_A1 --summary-a2 SUM_A2] --bakta-json BAKTA_JSON --genome GENOME \\
        --k12-fasta K12_FASTA [--sacct SACCT] [--a3 OUT_A3] --out METRICS_DIR
    python validation/pilot/metrics.py render --metrics METRICS_DIR/metrics.json \\
        [--notes METRICS_DIR/notes.md] --out METRICS_DIR/REPORT.md

Exit codes: 0 written; 1 an internal inconsistency only (metric 1 vs ``convert.json``, a render
placeholder missing under a computed metric, a rejected notes file); 2 a precondition (a missing,
unreadable or schema-invalid input, a failed run, a hash mismatch, git unavailable, ``--out``
inside the checkout). Nothing is written on exit 1 or 2.
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
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from importlib import resources
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


def checklist(metrics: Mapping[str, Any], run_states: Mapping[str, Any]) -> list[dict[str, str]]:
    """§15.6 checklist in pass 1 (spec table order). ``run_states``: arm -> exit code (None when
    absent) and ``compare_ok`` (the three compare tables present and valid under their
    schemas)."""
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
        {"bullet": "TE2-TE5 are reported; failures explained", "status": _bullet5(metrics)},
    ]


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
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
_MISSING = object()


class Gap(str):
    """An unresolved path under a metric that carries a ``reason``: rendered as a marker."""


def resolve(doc: Any, path: str) -> Any:
    """The value at a dotted path; ``Gap(reason)`` when it stops under a dict with a
    ``reason`` (the nearest one), else ``_MISSING``."""
    cur, reason = doc, None
    for part in path.split("."):
        if isinstance(cur, Mapping):
            reason = cur.get("reason", reason)
        if isinstance(cur, Mapping) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list) and part.isdigit() and int(part) < len(cur):
            cur = cur[int(part)]
        else:
            return _MISSING if reason is None else Gap(reason)
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
    return template.replace(M3_MARK, "\n".join(rows))


def missing_paths(template: str, metrics: Mapping[str, Any]) -> list[str]:
    paths = PLACEHOLDER.findall(expand(template, metrics))
    return [p for p in paths if resolve(metrics, p) is _MISSING]


def render(template: str, metrics: Mapping[str, Any], notes: str | None) -> str:
    """The report: each ``{{path}}`` becomes its formatted value and its source
    ``(metrics.json SOURCE_SEP parent)``; a path under a metric with a reason becomes
    ``**missing: <reason>**``, any other unresolved path ``**missing: <path>**``."""
    if notes is not None:
        check_notes(notes, metrics)

    def sub(m: re.Match[str]) -> str:
        path = m.group(1)
        v = resolve(metrics, path)
        if v is _MISSING or isinstance(v, Gap):
            return f"**missing: {path if v is _MISSING else v}**"
        parent = path.rsplit(".", 1)[0]
        return f"{fmt_value(path.rsplit('.', 1)[-1], v)} (metrics.json{SOURCE_SEP}{parent})"

    text = PLACEHOLDER.sub(sub, expand(template, metrics))
    body = notes.strip() if notes and notes.strip() else "No notes."
    if NOTES_MARK in text:
        return text.replace(NOTES_MARK, body)
    return text + ("" if notes is None else f"\n## Notes\n\n{body}\n")


# ---------------------------------------------------------------- CLI (file I/O)

#: hashed input files (argparse dests); ``cmp``, ``a1``, ``a2``, ``a3`` are directories
#: M10-06 inputs: ``te_real.json`` and ``te7.json`` (optional; absent -> ``te_placeholders()``)
TE_INPUTS = ("te_real", "te7")
INPUTS = ("peptides_a1", "peptides_a2", "convert_a1", "convert_a2", "summary_a1", "summary_a2",
          "bakta_json", "genome", "k12_fasta", "sacct", *TE_INPUTS)  # fmt: skip
OPTIONAL = ("summary_a1", "summary_a2", "sacct", "a3", *TE_INPUTS)
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


def _spec_sha() -> str:
    try:
        res = subprocess.run(
            ["git", "-C", str(REPO), "log", "-1", "--format=%H", "--", SPEC],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as err:
        raise Precondition(f"git log {SPEC}: {err}") from err
    if not re.fullmatch(r"[0-9a-f]{40}", sha := res.stdout.strip()):
        raise Precondition(f"git log {SPEC}: no commit found")
    return sha


def _outside_repo(path: Path, flag: str) -> None:
    p = path.resolve()
    if p == REPO or REPO in p.parents:
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


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Every metric from the parsed arguments (reads only; raises on a precondition)."""
    from giap.io.output_checks import check_tsv
    from giap.io.protein_fasta import read_protein_fasta

    _outside_repo(a.out, "--out")
    for k in ("cmp", *INPUTS):
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
    if a.a3 is not None:
        states["A3"] = _json(resolve_attempt(a.a3)[0] / "run_manifest.json").get("exit_code")
    doc |= te_placeholders()
    for key in TE_INPUTS:
        if getattr(a, key) is not None:
            doc[key] = _json(getattr(a, key))
    doc["checklist"] = checklist(doc, states)
    files = {k: getattr(a, k) for k in INPUTS} | {
        f"manifest_{x.lower()}": arms[x]["out"] / "run_manifest.json" for x in ARMS
    }
    doc["inputs"] = {
        k: {"path": str(Path(p).resolve()), "sha256": sha256(Path(p))}
        for k, p in files.items()
        if p is not None
    }
    return doc, exceptions


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


def _render(a: argparse.Namespace) -> int:
    _outside_repo(a.out, "render --out")
    metrics = _json(a.metrics)
    template = _need(a.template).read_text(encoding="utf-8")
    notes = None if a.notes is None else _need(a.notes).read_text(encoding="utf-8")
    text = render(template, metrics, notes)
    if missing := missing_paths(template, metrics):
        raise Inconsistent(f"render: missing paths in {a.metrics}: {', '.join(missing)}")
    if problems := diagnosis_problems(metrics, notes):
        raise Inconsistent("render: " + "; ".join(problems))
    a.out.write_text(text, encoding="utf-8")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for dest in (*REQUIRED, *OPTIONAL):
        c.add_argument(f"--{dest.replace('_', '-')}", type=Path, required=dest in REQUIRED)
    r = sub.add_parser("render")
    r.add_argument("--metrics", type=Path, required=True)
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
