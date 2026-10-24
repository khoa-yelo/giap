"""Real-data TE2-TE5 on the pilot arms and TE7 (a second A2 run, byte for byte); GIAP_SPEC.md
§11.6, §11.7, §15.6 bullet 5; spec docs/specs/M10-06-te-real-and-te7.md.

Pure functions (frames and dicts in, dicts out) plus a thin CLI::

    python validation/pilot/te_real.py compute --a2 OUT_A2 --a1 OUT_A1 \\
        --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --contaminants CONTAM --out METRICS_DIR
    python validation/pilot/te_real.py te7 --ref OUT_A --rep OUT_B --repo GIAP_REPO \\
        [--diann-a DIANN_A --diann-b DIANN_B --peptides-a PEP_A2 --peptides-b CONV_B] \\
        --out METRICS_DIR

Exit codes: 0 written, whatever the TE statuses; 1 an internal inconsistency; 2 a precondition.
``te_real_exceptions.tsv`` and ``te7_files.tsv`` hold row-level data and are never committed.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from itertools import zip_longest
from pathlib import Path, PurePosixPath
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from validation.pilot.metrics import (
    CODING,
    Q_ACCEPT,
    REPO,
    Inconsistent,
    Precondition,
    _cfg,
    _input_sha,
    _need,
    _q,
    _read_json,
    _routes,
    _rows,
    _run_table,
    _table,
    _target,
    clean,
    items,
    resolve_attempt,
    sha256,
    share,
)
from validation.proxies.te_proxies import TE2_EXEMPT, _flag_names, te2_result

from giap.annotate import bakta_layers as bl
from giap.config import GiapConfig
from giap.io.protein_fasta import resolve_protein_id
from giap.steps.mapping import select_queries

SPEC = "docs/specs/M10-06-te-real-and-te7.md"
TE3_MIN = 0.90  # §11.6, inclusive (R3: fixed before the real data)
TE4_PAIRS = (("rplL", "rplJ", 1, 3), ("atpA", "atpD", -1, 1))  # §11.6, closed intervals
RIBO_PRODUCT = re.compile(
    r"^(type [AB] )?(30S|50S) ribosomal protein [LS][0-9]+(/L[0-9]+)?( [A-Z][^ ]*)?$"
)
FAMILIES = {
    "trypsin": re.compile(r"^(Trypsin|Anionic trypsin|Cationic trypsin)$"),
    "keratin": re.compile(r"^Keratin, type (I|II) "),
    "keratin_associated": re.compile(r"^Keratin-associated protein"),
}
VERDICT_FAMILIES = ("trypsin", "keratin")
TE5_CLASSES = ("genome_exact", "bucket", "placed_tolerant", "labelled_unmapped", "absent")
TE5_VIOLATIONS = ("placed_tolerant", "labelled_unmapped", "absent")
TE2_COUNTS = ("alt_sites", "high", "high_start_codon", "high_stop_codon", "violations")
CODE_PATHS = ("src/", "workflow/", "pyproject.toml", "uv.lock", "scripts/slurm/giap_run.sbatch")
EXCLUDED_TOP = ("work", "logs", "report", "preflight")
Reader = Callable[[], bytes]


def _s(v: Any) -> str:
    c = clean(v)
    return "" if c is None else str(c)


def _true(v: Any) -> bool:
    c = clean(v)
    return c is True or (isinstance(c, str) and c.lower() in ("true", "1"))


def _tally(keys: Iterable[str]) -> dict[str, int]:
    c = Counter(keys)
    return {k: c[k] for k in sorted(c)}


def _status(parts: Sequence[str]) -> str:
    """Combined status: fail if any fails, pass if all pass, else not_evaluable."""
    if "fail" in parts:
        return "fail"
    return "pass" if parts and all(p == "pass" for p in parts) else "not_evaluable"


# ---------------------------------------------------------------- TE2


def _te2_bad(variants: pd.DataFrame) -> pd.DataFrame:
    """The ``high`` rows with neither exempt flag (te2_result's violations, all of them)."""
    high = (variants["confidence"].astype(str) == "high").to_numpy(dtype=bool)
    exempt = np.array([bool(TE2_EXEMPT & _flag_names(f)) for f in variants["flags"]], dtype=bool)
    return variants[high & ~exempt]


def _n_alt_bin(v: Any) -> str:
    c = clean(v)
    return "NA" if c is None else ">=3" if int(c) >= 3 else str(int(c))


def te2(
    variants: pd.DataFrame | None, genes: pd.DataFrame, detection_summary: pd.DataFrame | None
) -> dict[str, Any]:
    """TE2: no ``high`` ALT call without ``start_codon``/``stop_codon`` (te2_result), with
    breakdowns of every violation and the ``id_mixups`` cross-check."""
    warnings: list[dict[str, Any]] = []
    if variants is None:
        return {"status": "not_evaluable", "reason": ["missing:tables/variants.tsv"],
                "counts": dict.fromkeys(TE2_COUNTS), "violations": [], "vacuous": False,
                "breakdown": {}, "warnings": warnings}  # fmt: skip
    res = te2_result(variants)
    bad = _te2_bad(variants)
    tags = genes["locus_tag"].astype(str)
    rep = dict(zip(tags, genes["replicon_type"].map(_s), strict=True))
    reg = dict(zip(tags, genes["region_ids"], strict=True))
    region: list[str] = []
    for g in bad["gene"].map(_s):
        region += sorted({re.sub(r"_\d+$", "", r) for r in items(reg.get(g))}) or ["none"]
    breakdown = {
        "alt_type": _tally(bad["alt_type"].map(lambda v: _s(v) or "NA")),
        "call": _tally(bad["call"].map(lambda v: _s(v) or "NA")),
        "flags": _tally(
            ";".join(sorted(_flag_names(f) - TE2_EXEMPT)) or "none" for f in bad["flags"]
        ),
        "n_alt": _tally(_n_alt_bin(v) for v in bad["n_alt"]),
        "replicon": _tally(rep.get(g) or "not_in_genes" for g in bad["gene"].map(_s)),
        "region": _tally(region),
        "aa_pos_1_unflagged": int((pd.to_numeric(bad["aa_pos"], errors="coerce") == 1).sum()),
    }
    n_bad = res.counts["violations"]
    mix = None
    if detection_summary is not None:
        ds = detection_summary
        row = ds[(ds["section"].astype(str) == "qc") & (ds["metric"].astype(str) == "id_mixups")]
        mix = clean(row["value"].iloc[0]) if len(row) else None
    if mix is None:
        warnings.append({"code": "id_mixups_na"})
    elif n_bad > mix:
        warnings.append(
            {"code": "te2_id_mixups_inconsistent", "violations": n_bad, "id_mixups": mix}
        )
    out = {"status": "pass" if n_bad == 0 else "fail", "counts": dict(res.counts),
           "violations": list(res.exceptions), "vacuous": len(variants) == 0,
           "breakdown": breakdown, "id_mixups": mix, "warnings": warnings}  # fmt: skip
    if n_bad:
        out["reason"] = ["violations"]
    return out


# ---------------------------------------------------------------- TE3


def _coding(genes: pd.DataFrame) -> pd.DataFrame:
    return genes[genes["feature_type"].astype(str).isin(CODING)].reset_index(drop=True)


def _r_mask(G: pd.DataFrame) -> np.ndarray:
    """Set R (question L1): the strict product rule or §7.1.2's ``RIBOSOMAL_GENE``."""
    return np.array(
        [
            bool(RIBO_PRODUCT.match(_s(p))) or bool(bl.RIBOSOMAL_GENE.match(_s(g)))
            for p, g in zip(G["product"], G["gene"], strict=True)
        ],
        dtype=bool,
    )


def ribosomal_set(genes: pd.DataFrame) -> pd.DataFrame:
    """The genes of **G** (cds, sorf) in set R."""
    G = _coding(genes)
    return G[_r_mask(G)]


def _values(G: pd.DataFrame, col: str = "abundance_log2") -> pd.Series:
    return pd.to_numeric(G[col], errors="coerce").astype(float)


def _measured(G: pd.DataFrame) -> np.ndarray:
    imputed = np.array([_true(v) for v in G["abundance_imputed"]], dtype=bool)
    out: np.ndarray = _values(G).notna().to_numpy() & ~imputed
    return out


def _types(G: pd.DataFrame, mask: np.ndarray) -> str | None:
    return ";".join(sorted({_s(t) for t in G.loc[mask, "abundance_type"]} - {""})) or None


def _te3_clauses(G: pd.DataFrame, R: np.ndarray) -> dict[str, Any]:
    n = int(R.sum())
    k = int(((G["detection_status"].astype(str) == "detected").to_numpy() & R).sum())
    meas, vals = _measured(G), _values(G).to_numpy()
    qv, rv = vals[meas], vals[meas & R]
    atype = _types(G, meas)
    ab: dict[str, Any] = {"m": None, "q3": None, "percentile_rank": None, "holds": None,
                          "n_r_measured": len(rv), "n_measured": len(qv),
                          "abundance_type": atype}  # fmt: skip
    if len(rv) and atype is not None and "none" not in atype.split(";"):
        m, q3 = float(np.median(rv)), float(np.quantile(qv, 0.75))
        ab |= {"m": m, "q3": q3, "percentile_rank": float(np.mean(qv <= m)), "holds": m >= q3}
    return {"n_ribosomal": n, "abundance": ab,
            "detection": {"share": share(k, n), "holds": n > 0 and k / n >= TE3_MIN}}  # fmt: skip


def te3(genes: pd.DataFrame) -> dict[str, Any]:
    """TE3: R >= 90 % detected (clause 1) and median abundance of R >= the 75th percentile of the
    measured genes (clause 2); secondaries never change the status."""
    G = _coding(genes)
    R = _r_mask(G)
    out = _te3_clauses(G, R)
    type_b = G["product"].map(_s).str.startswith("type B").to_numpy(dtype=bool)
    s712 = [
        bl.is_ribosomal_protein(SimpleNamespace(gene=_s(g) or None, product=_s(p) or None))  # type: ignore[arg-type]
        for g, p in zip(G["gene"], G["product"], strict=True)
    ]
    out["secondary"] = {
        "without_type_b": _te3_clauses(G, R & ~type_b),
        "section_7_1_2": _te3_clauses(G, np.array(s712, dtype=bool)),
    }
    reason = [] if out["detection"]["holds"] else ["detection"]
    reason += ["abundance"] if out["abundance"]["holds"] is False else []
    if out["n_ribosomal"] == 0:
        out |= {"status": "not_evaluable", "reason": ["no_ribosomal_genes"]}
    elif reason:
        out |= {"status": "fail", "reason": reason}
    elif out["abundance"]["holds"] is None:
        out |= {"status": "not_evaluable", "reason": ["abundance_not_evaluable"]}
    else:
        out["status"] = "pass"
    return out


# ---------------------------------------------------------------- TE4


def te4(genes: pd.DataFrame) -> dict[str, Any]:
    """TE4: log2(RplL/RplJ) in [1, 3] and log2(AtpA/AtpD) in [-1, 1] on a molar scale."""
    G = _coding(genes)
    names = G["gene"].map(_s).to_numpy()
    meas, vals, top3 = _measured(G), _values(G), _values(G, "top3_log2")
    rows = {n: np.flatnonzero(names == n) for p in TE4_PAIRS for n in p[:2]}
    pairs: dict[str, Any] = {}
    for a, b, lo, hi in TE4_PAIRS:
        reason = [f"gene_{'missing' if not len(rows[x]) else 'ambiguous'}:{x}"
                  for x in (a, b) if len(rows[x]) != 1]  # fmt: skip
        ratio = t3 = None
        if not reason:
            i, j = int(rows[a][0]), int(rows[b][0])
            if any(_s(G["abundance_scale"].iloc[x]) != "molar" for x in (i, j)):
                reason.append("relative_scale")
            reason += [f"not_measured:{n}" for n, x in ((a, i), (b, j)) if not meas[x]]
            if meas[i] and meas[j]:
                ratio = float(vals.iloc[i]) - float(vals.iloc[j])
            if top3.iloc[[i, j]].notna().all():
                t3 = float(top3.iloc[i]) - float(top3.iloc[j])
        st = (
            "not_evaluable"
            if reason
            else "pass"
            if ratio is not None and lo <= ratio <= hi
            else "fail"
        )
        pairs[f"{a}/{b}"] = {"status": st, "reason": reason, "low": lo, "high": hi,
                             "log2_ratio": ratio, "top3_log2_ratio": t3}  # fmt: skip
    sec = {
        n: {k: clean(G[k].iloc[ix[0]]) if len(ix) == 1 else None
            for k in ("n_peptides_unique", "n_observable_peptides")}
        for n, ix in rows.items()
    }  # fmt: skip
    status = _status([p["status"] for p in pairs.values()])
    out = {"status": status, "pairs": pairs, "secondary": {"genes": sec},
           "abundance_type": _types(G, np.ones(len(G), dtype=bool))}  # fmt: skip
    if status != "pass":
        out["reason"] = [f"{k}:{p['status']}" for k, p in pairs.items() if p["status"] != "pass"]
    return out


# ---------------------------------------------------------------- TE5


def contaminant_families(headers: Iterable[str]) -> dict[str, set[str]]:
    """Family -> resolved IDs, by the rule on the description (header text up to `` OS=``)."""
    out: dict[str, set[str]] = {f: set() for f in FAMILIES}
    for h in headers:
        tok, _, rest = h.lstrip(">").strip().partition(" ")
        desc = rest.split(" OS=", 1)[0].strip()
        for fam, rule in FAMILIES.items():
            if rule.match(desc):
                out[fam].add(resolve_protein_id(tok)[0])
    return out


def _te5_judge(
    own: Mapping[str, set[str]],
    routes: Mapping[str, set[str]],
    labels: Mapping[str, Mapping[str, set[str]]],
    families: Mapping[str, Collection[str]],
    fams: Sequence[str],
) -> tuple[dict[str, Any], dict[str, str]]:
    ids = set().union(*(set(families.get(f, ())) for f in fams))
    S = sorted(p for p, o in own.items() if o & ids)
    counts = dict.fromkeys((*TE5_CLASSES, "match_other_entry"), 0)
    classes: dict[str, str] = {}
    for p in S:
        r, lab = routes.get(p, set()), labels.get(p, {})
        if "exact" in r:
            c = "genome_exact"
        elif "contaminant" in lab:
            c = "bucket"
            counts["match_other_entry"] += bool(lab["contaminant"] - {""} - own[p])
        elif "tolerant" in r:
            c = "placed_tolerant"
        elif "unmapped" in lab:
            c = "labelled_unmapped"
        else:
            c = "absent"
        counts[c] += 1
        classes[p] = c
    per = {f: sum(bool(own[p] & set(families.get(f, ()))) for p in S) for f in fams}
    if any(v > len(S) for v in per.values()) or sum(counts[c] for c in TE5_CLASSES) != len(S):
        raise Inconsistent("TE5: family or class counts exceed |S|")
    return {"n_s": len(S), "families": per, "counts": counts}, classes


def _te5(
    P: pd.DataFrame,
    peptide_hits: pd.DataFrame,
    unmapped: pd.DataFrame,
    families: Mapping[str, Collection[str]],
    cfg: GiapConfig,
) -> tuple[dict[str, Any], dict[str, str]]:
    U = set(select_queries(P, cfg)["peptide"].astype(str))
    acc = _target(P) & (_q(P) <= Q_ACCEPT)
    own: dict[str, set[str]] = {}
    for pep, ids, ok in zip(P["peptide"].astype(str), P["protein_ids"], acc, strict=True):
        if ok and pep in U:
            own.setdefault(pep, set()).update(items(ids))
    um = unmapped[unmapped["kind"].astype(str) == "peptide"]
    labels: dict[str, dict[str, set[str]]] = {}
    for sid, lab, m in zip(um["sequence_id"].astype(str), um["label"].astype(str),
                           um["contaminant_match"], strict=True):  # fmt: skip
        labels.setdefault(sid, {}).setdefault(lab, set()).add(_s(m))
    routes = _routes(peptide_hits)
    res, classes = _te5_judge(own, routes, labels, families, VERDICT_FAMILIES)
    ka, _ = _te5_judge(own, routes, labels, families, ("keratin_associated",))
    c = res["counts"]
    denom = res["n_s"] - c["genome_exact"]
    bad = sum(c[v] for v in TE5_VIOLATIONS)
    res |= {"share_bucket": share(c["bucket"], denom), "secondary": {"keratin_associated": ka},
            "status": "fail" if bad else "pass" if denom > 0 else "not_evaluable"}  # fmt: skip
    if res["status"] != "pass":
        res["reason"] = ["violations"] if bad else ["no_trypsin_keratin_peptides"]
    return res, {p: k for p, k in classes.items() if k in TE5_VIOLATIONS}


def te5_arm(
    P: pd.DataFrame,
    peptide_hits: pd.DataFrame,
    unmapped: pd.DataFrame,
    families: Mapping[str, Collection[str]],
    cfg: GiapConfig,
) -> dict[str, Any]:
    """TE5 on one arm: family peptides of **U** (accepted target rows) land in the bucket."""
    return _te5(P, peptide_hits, unmapped, families, cfg)[0]


# ---------------------------------------------------------------- TE7, GIAP layer


def path_class(path: str) -> str:
    top = path.split("/", 1)[0]
    if top in EXCLUDED_TOP or top.startswith(".nextflow") or path == "run_manifest.json":
        return "excluded"
    if path.startswith("annotation/tools/"):
        return "external"
    if path.startswith("plots/") and not path.startswith("plots/data/"):
        return "figures"
    return "gate"


def _numeric(la: list[str], lb: list[str]) -> dict[str, Any] | None:
    if len(la) != len(lb) or not la or la[0] != lb[0]:
        return None
    header = la[0].split("\t")
    cols: set[str] = set()
    n, max_abs, max_rel = 0, 0.0, 0.0
    for x, y in zip(la[1:], lb[1:], strict=True):
        if x == y:
            continue
        cx, cy = x.split("\t"), y.split("\t")
        if len(cx) != len(cy):
            return None
        for j, (u, v) in enumerate(zip(cx, cy, strict=True)):
            if u == v:
                continue
            try:
                fu, fv = float(u), float(v)
            except ValueError:
                return None
            n += 1
            cols.add(header[j] if j < len(header) else str(j))
            d, scale = abs(fu - fv), max(abs(fu), abs(fv))
            if math.isfinite(d):
                max_abs, max_rel = max(max_abs, d), max(max_rel, d / scale if scale else 0.0)
    order = [h for h in dict.fromkeys(header) if h in cols] + sorted(cols - set(header))
    return {"columns": order, "n_cells": n, "max_abs": max_abs, "max_rel": max_rel} if n else None


def _diagnose(
    path: str,
    ba: bytes | None,
    bb: bytes | None,
    labels: Sequence[str],
    known: Mapping[str, Collection[str]],
    producer: Mapping[str, Any],
) -> dict[str, Any]:
    d: dict[str, Any] = dict.fromkeys(
        ("sub_class", "columns", "n_cells", "max_abs", "max_rel", "first_differing_line")
    )
    d |= {"path": path, "producer": producer.get("producer", "unknown"),
          "cpu_classes": producer.get("cpu_classes")}  # fmt: skip
    if ba is None or bb is None:
        return d | {"class": f"only_in:{labels[0] if bb is None else labels[1]}"}
    if path.endswith(".gz"):
        try:
            ba, bb = gzip.decompress(ba), gzip.decompress(bb)
        except (OSError, EOFError):
            return d | {"class": "content"}
        if ba == bb:
            return d | {"class": "gzip_header_only"}
    try:
        la, lb = ba.decode("utf-8").splitlines(), bb.decode("utf-8").splitlines()
    except UnicodeDecodeError:
        return d | {"class": "content"}
    if la != lb and la[:1] == lb[:1] and sorted(la[1:]) == sorted(lb[1:]):
        return d | {"class": "row_order_only"}
    num = _numeric(la, lb)
    if num is not None:
        name = PurePosixPath(path).name
        listed = set(known.get(name, known.get(name.removesuffix(".gz"), ())))
        cc = producer.get("cpu_classes") or {}
        lapack = set(num["columns"]) <= listed and all(x in cc for x in labels[:2])
        lapack = lapack and set(cc[labels[0]]) != set(cc[labels[1]])
        return d | num | {"class": "numeric_only", "sub_class": "known_lapack" if lapack else None}
    line = next((i + 1 for i, (x, y) in enumerate(zip_longest(la, lb)) if x != y), None)
    return d | {"class": "content", "first_differing_line": line}


def _classify(
    files_a: Mapping[str, Reader],
    files_b: Mapping[str, Reader],
    known: Mapping[str, Collection[str]],
    producers: Mapping[str, Any],
    labels: Sequence[str],
) -> tuple[dict[str, Any], list[list[Any]]]:
    counts: Counter[str] = Counter()
    found: dict[str, list[Any]] = {"gate": [], "external": [], "figures": []}
    rows: list[list[Any]] = []
    for p in sorted(set(files_a) | set(files_b)):
        cls = path_class(p)
        counts[cls] += 1
        if cls == "excluded":
            continue
        ba = files_a[p]() if p in files_a else None
        bb = files_b[p]() if p in files_b else None
        ha, hb = (None if x is None else hashlib.sha256(x).hexdigest() for x in (ba, bb))
        rows.append([p, cls, ha, hb, ha is not None and ha == hb])
        if ha is not None and ha == hb:
            continue
        if cls == "gate":
            found[cls].append(_diagnose(p, ba, bb, labels, known, producers.get(p) or {}))
        else:
            found[cls].append({"path": p})
    giap = {
        "status": "fail" if found["gate"] else "pass",
        "differing": found["gate"],
        "external_differing": found["external"],
        "figures_differing": found["figures"],
        **{f"n_{c}": counts[c] for c in ("gate", "external", "figures", "excluded")},
    }
    return giap, rows


def classify_tree(
    files_a: Mapping[str, Reader],
    files_b: Mapping[str, Reader],
    known_columns: Mapping[str, Collection[str]],
    producers: Mapping[str, Any],
    *,
    labels: Sequence[str] = ("ref", "rep1"),
) -> dict[str, Any]:
    """te7.json's ``giap`` object: the byte comparison of the gate files decides TE7."""
    return _classify(files_a, files_b, known_columns, producers, labels)[0]


# ---------------------------------------------------------------- TE7, DIA-NN layer

L3_KEY = ["run", "modified_peptide", "charge", "is_decoy"]


def _keyed(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    for c in ("run", "modified_peptide", "charge"):
        d[c] = d[c].map(_s)
    d["is_decoy"] = [_true(v) for v in d["is_decoy"]]
    for c in ("q_value", "intensity"):
        d[c] = pd.to_numeric(d[c], errors="coerce").astype(float)
    if "peptide" not in d.columns:
        d["peptide"] = d["modified_peptide"].str.replace(r"\(.*?\)", "", regex=True)
    d["acc"] = ~d["is_decoy"] & (d["q_value"] <= Q_ACCEPT)
    return d.drop_duplicates(L3_KEY).set_index(L3_KEY)


def _jaccard(a: set[Any], b: set[Any]) -> float | None:
    return len(a & b) / len(a | b) if a | b else None


def _differs(x: pd.Series, y: pd.Series) -> int:
    u, v = x.to_numpy(dtype=float), y.to_numpy(dtype=float)
    return int((~((u == v) | (np.isnan(u) & np.isnan(v)))).sum())


def _l3(conv_a: pd.DataFrame, conv_b: pd.DataFrame) -> dict[str, Any]:
    ia, ib = _keyed(conv_a), _keyed(conv_b)
    both = ia.index.intersection(ib.index)
    A, B = ia.loc[both], ib.loc[both]
    ratio = np.log2(B["intensity"].to_numpy() / A["intensity"].to_numpy())
    ok = (A["intensity"].to_numpy() > 0) & (B["intensity"].to_numpy() > 0) & np.isfinite(ratio)
    return {
        "only_a": len(ia) - len(both),
        "only_b": len(ib) - len(both),
        "both": len(both),
        "jaccard_precursors": _jaccard(set(ia.index[ia["acc"]]), set(ib.index[ib["acc"]])),
        "jaccard_peptides": _jaccard(
            set(ia.loc[ia["acc"], "peptide"]), set(ib.loc[ib["acc"], "peptide"])
        ),
        "q_value_differs": _differs(A["q_value"], B["q_value"]),
        "acceptance_flips": int((A["acc"].to_numpy() != B["acc"].to_numpy()).sum()),
        "intensity_differs": _differs(A["intensity"], B["intensity"]),
        "max_abs_log2_intensity_ratio": float(np.abs(ratio[ok]).max()) if ok.any() else None,
    }  # fmt: skip


def diann_layer(
    sums_a: Mapping[str, str],
    sums_b: Mapping[str, str],
    conv_a: pd.DataFrame | None,
    conv_b: pd.DataFrame | None,
    *,
    conv_sha: tuple[str, str],
) -> dict[str, Any]:
    """L0 (file list), L1 (report sha256), L2 (converted table sha256), L3 (aggregates)."""
    names = sorted(set(sums_a) | set(sums_b))
    equal = [n for n in names if n in sums_a and sums_a[n] == sums_b.get(n)]
    l0 = {"files_equal": equal, "files_differ": [n for n in names if n not in equal]}
    report = [n for n in names if n.endswith("report.parquet")]
    if report and all(n in equal for n in report):
        return {"status": "identical_report", "l0": l0, "l3": None}
    if conv_sha[0] == conv_sha[1]:
        return {"status": "identical_converted", "l0": l0, "l3": None}
    if conv_a is None or conv_b is None:
        raise ValueError("diann_layer: converted tables are needed when L2 differs")
    return {"status": "nondeterministic", "l0": l0, "l3": _l3(conv_a, conv_b)}


# ---------------------------------------------------------------- TE7 pair checks


def _no_outdir(v: Any) -> Any:
    if isinstance(v, Mapping):
        return {k: _no_outdir(x) for k, x in v.items() if k != "outdir"}
    return [_no_outdir(x) for x in v] if isinstance(v, list) else v


def pair_problems(ref: Mapping[str, Any], rep: Mapping[str, Any]) -> list[str]:
    """The manifest fields that stop two runs from forming a TE7 pair (empty: a pair)."""
    bad = []
    for lab, m in (("ref", ref), ("rep", rep)):
        bad += [f"{lab}.exit_code"] if m.get("exit_code") != 0 else []
        bad += [f"{lab}.resume"] if m.get("resume") is not False else []
        bad += [f"{lab}.git_dirty"] if m.get("git_dirty") is not False else []
        bad += (
            [f"{lab}.steps.cached"]
            if any(s.get("status") == "cached" for s in m.get("steps", []))
            else []
        )
    pick: dict[str, Callable[[Mapping[str, Any]], Any]] = {
        "inputs.sha256": lambda m: sorted((i.get("role"), i.get("sha256"))
                                          for i in m.get("inputs", [])),
        "seeds": lambda m: m.get("seeds"),
        "databases": lambda m: sorted((d.get("name"), d.get("version"), d.get("checksum"))
                                      for d in m.get("databases") or []),
        "config": lambda m: _no_outdir(m.get("config")),
    }  # fmt: skip
    return bad + [k for k, f in pick.items() if f(ref) != f(rep)]


# ---------------------------------------------------------------- CLI (file I/O)


def _spec_sha() -> str:
    res = subprocess.run(["git", "-C", str(REPO), "log", "-1", "--format=%H", "--", SPEC],
                         capture_output=True, text=True, check=False)  # fmt: skip
    sha = res.stdout.strip()
    return sha if re.fullmatch(r"[0-9a-f]{40}", sha) else "0" * 40


def _out_dir(out: Path) -> Path:
    o = out.resolve()
    if o == REPO or REPO in o.parents:
        raise Precondition(f"--out {out} is inside the repository")
    return out


def _dump(doc: Any) -> str:
    return json.dumps(doc, sort_keys=True, indent=2) + "\n"


def _write(
    out: Path, name: str, doc: Any, tsv: str, cols: list[str], rows: list[list[Any]]
) -> None:
    from giap.io.tsv import write_tsv

    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(_dump(doc), encoding="utf-8")
    write_tsv(out / tsv, cols, sorted(rows, key=lambda r: [str(x) for x in r]))


def _manifest(out: Path) -> dict[str, Any]:
    man = _read_json(out / "run_manifest.json")
    if man.get("exit_code") != 0:
        raise Precondition(f"{out / 'run_manifest.json'}: exit_code {man.get('exit_code')}")
    return man


def _pinned(man: Mapping[str, Any], fasta: Path) -> bool:
    """False when the manifest records no contaminants checksum; raises when it differs."""
    db = next((d for d in man.get("databases") or [] if d.get("name") == "contaminants"), None)
    if not db or not db.get("checksum"):
        return False
    algo, _, hexv = str(db["checksum"]).rpartition(":")
    if hashlib.new(algo or "sha256", fasta.read_bytes()).hexdigest() != hexv:
        raise Precondition(
            f"{fasta}: checksum differs from the run's contaminants {db['checksum']}"
        )
    return True


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[list[Any]]]:
    """te_real.json and the exception rows (reads only; raises on a precondition)."""
    from giap.io.peptides import read_peptide_table

    _out_dir(a.out)
    contam = _need(a.contaminants)
    arms: dict[str, dict[str, Any]] = {}
    unpinned: list[str] = []
    for label, d, pep in (("A1", a.a1, a.peptides_a1), ("A2", a.a2, a.peptides_a2)):
        out = _need(resolve_attempt(d)[0])
        man = _manifest(out)
        if _input_sha(man, "peptides") != sha256(_need(pep)):
            raise Precondition(f"{pep}: sha256 differs from the run's peptides input ({out})")
        unpinned += [] if _pinned(man, contam) else [label]
        hits, unm = (_run_table(out, n) for n in ("peptide_hits.tsv", "unmapped.tsv"))
        arms[label] = {"out": out, "man": man, "pep": pep, "P": read_peptide_table(pep)[0],
                       "hits": hits and _table(hits, "peptide_hits"),
                       "unm": unm and _table(unm, "unmapped")}  # fmt: skip
    out2 = arms["A2"]["out"]
    genes = _table(_need(out2 / "tables" / "genes.tsv"), "genes")
    vpath, dpath = out2 / "tables" / "variants.tsv", _run_table(out2, "detection_summary.tsv")
    variants = _table(vpath, "variants") if vpath.is_file() else None
    doc: dict[str, Any] = {"schema": "giap-pilot-te-real/1", "spec": _spec_sha()}
    doc["te2"] = te2(variants, genes, dpath and _table(dpath, "detection_summary"))
    doc["te3"], doc["te4"] = te3(genes), te4(genes)
    heads = [ln[1:] for ln in contam.read_text(encoding="utf-8").splitlines() if ln.startswith(">")]
    fam = contaminant_families(heads)
    rows: list[list[Any]] = []
    doc["te5"] = {}
    for label, arm in arms.items():
        if arm["hits"] is None or arm["unm"] is None:
            doc["te5"][label] = {
                "status": "not_evaluable",
                "share_bucket": share(0, 0),
                "reason": ["missing:peptide_hits.tsv/unmapped.tsv"],
            }
            continue
        doc["te5"][label], bad = _te5(arm["P"], arm["hits"], arm["unm"], fam, _cfg(arm["man"]))
        rows += [["TE5", f"te5_{k}", p, label] for p, k in bad.items()]
    doc["te5"]["status"] = _status([doc["te5"]["A1"]["status"], doc["te5"]["A2"]["status"]])
    if variants is not None:
        for vid, g, pos in _te2_bad(variants)[["variant_id", "gene", "aa_pos"]].itertuples(
            index=False
        ):
            rows.append(["TE2", "te2_violation", _s(g), f"variant_id={_s(vid)};aa_pos={_s(pos)}"])
    G = _coding(genes)
    miss = G[_r_mask(G) & (G["detection_status"].astype(str) != "detected").to_numpy()]
    for r in miss.itertuples(index=False):
        detail = (f"gene={_s(r.gene)};product={_s(r.product)};"
                  f"detection_status={_s(r.detection_status)};"
                  f"n_observable={_s(r.n_observable_peptides)}")  # fmt: skip
        rows.append(["TE3", "te3_not_detected", r.locus_tag, detail])
    warnings = list(doc["te2"]["warnings"])
    warnings += [{"code": "contaminants_unpinned", "arms": unpinned}] if unpinned else []
    doc["warnings"] = warnings
    files = {"contaminants": contam, "genes_a2": out2 / "tables" / "genes.tsv",
             **({"variants_a2": vpath} if variants is not None else {})}  # fmt: skip
    for label, arm in arms.items():
        files[f"manifest_{label.lower()}"] = arm["out"] / "run_manifest.json"
        files[f"peptides_{label.lower()}"] = arm["pep"]
    doc["inputs"] = {k: {"path": str(Path(p).resolve()), "sha256": sha256(Path(p))}
                     for k, p in files.items()}  # fmt: skip
    return doc, rows


def _tree(root: Path) -> dict[str, Reader]:
    out: dict[str, Reader] = {}
    for d, _dirs, names in os.walk(root):
        for n in names:
            p = Path(d) / n
            if p.is_file():
                out[p.relative_to(root).as_posix()] = p.read_bytes
    return out


def _step_classes(out: Path) -> tuple[list[str], dict[str, list[str]]]:
    """Trace hosts and, per process (step), its sorted CPU classes."""
    from scripts.measure_runtimes import node_class, process_name

    hosts: set[str] = set()
    per: dict[str, set[str]] = {}
    for r in _rows(out / "logs" / "trace.tsv"):
        h = (r.get("hostname") or "").strip()
        if h and h != "-":
            hosts.add(h)
            per.setdefault(process_name(r.get("name", "")), set()).add(node_class(h))
    return sorted(hosts), {k: sorted(v) for k, v in sorted(per.items())}


def _producers(
    paths: Iterable[str], classes: Mapping[str, Mapping[str, list[str]]]
) -> dict[str, Any]:
    """OUT path -> its producing registered step (output key = file name with ``_`` for ``.``, or
    its stem) and that step's CPU classes per run."""
    import giap.steps
    from giap.steps.registry import all_steps

    giap.steps.load_step_modules()
    owner = {k: s.name for s in all_steps() for k in s.outputs}
    out: dict[str, Any] = {}
    for p in paths:
        name = PurePosixPath(p).name
        step = owner.get(name.replace(".", "_")) or owner.get(name.split(".", 1)[0])
        if step is not None:
            out[p] = {
                "producer": step,
                "cpu_classes": {k: list(v.get(step, [])) for k, v in classes.items()},
            }
    return out


def _known_columns() -> dict[str, set[str]]:
    rows = _rows_tsv(Path(__file__).with_name("te7_known_columns.tsv"))
    out: dict[str, set[str]] = {}
    for r in rows:
        out.setdefault(PurePosixPath(r["table"]).name, set()).add(r["column"])
    return out


def _rows_tsv(path: Path) -> list[dict[str, str]]:
    import csv

    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _git_equivalent(repo: Path, a: str, b: str) -> bool:
    res = subprocess.run(["git", "-C", str(repo), "diff", "--quiet", a, b, "--", *CODE_PATHS],
                         capture_output=True, text=True, check=False)  # fmt: skip
    if res.returncode not in (0, 1):
        raise Precondition(
            f"git_commit: cannot compare {a} and {b} in {repo}: {res.stderr.strip()}"
        )
    return res.returncode == 0


def _diann(a: argparse.Namespace) -> dict[str, Any]:
    given = [a.diann_a, a.diann_b, a.peptides_a, a.peptides_b]
    if not any(given):
        return {"status": "not_run", "reason": "not_given", "l0": None, "l3": None, "nodes": []}
    if not all(given):
        raise Precondition("--diann-a, --diann-b, --peptides-a and --peptides-b go together")
    sums, nodes = [], []
    for d in (a.diann_a, a.diann_b):
        text = _need(d / "outputs.sha256").read_text(encoding="utf-8")
        sums.append({ln.split(None, 1)[1].strip().lstrip("*"): ln.split(None, 1)[0]
                     for ln in text.splitlines() if ln.strip()})  # fmt: skip
        prov = d / "provenance.json"
        nodes.append(_read_json(prov).get("node") if prov.is_file() else None)
    shas = (sha256(_need(a.peptides_a)), sha256(_need(a.peptides_b)))

    def read(p: Path) -> pd.DataFrame:
        return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, compression="infer")

    lazy = shas[0] != shas[1] and not _same_report(sums[0], sums[1])
    conv = (read(a.peptides_a), read(a.peptides_b)) if lazy else (None, None)
    return diann_layer(sums[0], sums[1], conv[0], conv[1], conv_sha=shas) | {"nodes": nodes}


def _same_report(sa: Mapping[str, str], sb: Mapping[str, str]) -> bool:
    rep = [n for n in sa if n.endswith("report.parquet")]
    return bool(rep) and all(sa[n] == sb.get(n) for n in rep)


def te7(a: argparse.Namespace) -> tuple[dict[str, Any], list[list[Any]]]:
    """te7.json and the te7_files.tsv rows (reads only; raises on a precondition)."""
    _out_dir(a.out)
    (ref, ref_job), (rep, rep_job) = resolve_attempt(a.ref), resolve_attempt(a.rep)
    mans = [
        _read_json(_need(ref) / "run_manifest.json"),
        _read_json(_need(rep) / "run_manifest.json"),
    ]
    if problems := pair_problems(*mans):
        raise Precondition("not a TE7 pair: " + ", ".join(problems))
    commits = [str(m.get("git_commit") or "") for m in mans]
    if not all(commits) or not _git_equivalent(a.repo, *commits):
        raise Precondition(f"git_commit: {commits[0]} and {commits[1]} are not code-equivalent")
    m = re.search(r"rep\d+", rep.name)
    labels = ("ref", m.group(0) if m else "rep1")
    pair, classes = {}, {}
    for lab, out, man, job in (
        (labels[0], ref, mans[0], ref_job),
        (labels[1], rep, mans[1], rep_job),
    ):
        hosts, classes[lab] = _step_classes(out)
        pair[lab] = {"outdir": str(out.resolve()), "git_commit": man.get("git_commit"), "job": job,
                     "hosts": hosts, "cpu_classes_by_step": classes[lab]}  # fmt: skip
    fa, fb = _tree(ref), _tree(rep)
    gate = [p for p in sorted(set(fa) | set(fb)) if path_class(p) == "gate"]
    giap, rows = _classify(fa, fb, _known_columns(), _producers(gate, classes), labels)
    diann = _diann(a)
    e2e = giap["status"] == "pass" and diann["status"] in (
        "identical_report",
        "identical_converted",
    )
    doc = {"schema": "giap-pilot-te7/1", "spec": _spec_sha(), "pair": pair, "giap": giap,
           "diann": diann, "end_to_end": e2e}  # fmt: skip
    return doc, rows


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for name in ("a1", "a2", "peptides-a1", "peptides-a2", "contaminants", "out"):
        c.add_argument(f"--{name}", type=Path, required=True)
    t = sub.add_parser("te7")
    for name in ("ref", "rep", "repo", "out"):
        t.add_argument(f"--{name}", type=Path, required=True)
    for name in ("diann-a", "diann-b", "peptides-a", "peptides-b"):
        t.add_argument(f"--{name}", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "compute":
            doc, rows = compute(a)
            _write(
                a.out,
                "te_real.json",
                doc,
                "te_real_exceptions.tsv",
                ["te", "kind", "id", "detail"],
                rows,
            )
            st = {k: doc[k]["status"] for k in ("te2", "te3", "te4", "te5")}
        else:
            doc, rows = te7(a)
            cols = ["path", "class", "sha256_a", "sha256_b", "equal"]
            _write(a.out, "te7.json", doc, "te7_files.tsv", cols, rows)
            st = {"te7": doc["giap"]["status"], "diann": doc["diann"]["status"]}
        print(f"wrote {a.out}: " + ", ".join(f"{k} {v}" for k, v in st.items()))
        return 0
    except Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except (Inconsistent, ValueError) as err:
        print(f"inconsistent: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
