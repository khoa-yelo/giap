"""Real-data TE2-TE5 on the pilot arms; GIAP_SPEC.md §11.6, §15.6 bullet 5; spec
docs/specs/M10-06-te-real-and-te7.md. TE7 (a second A2 run, byte for byte) lives in te7.py; its
subcommand stays here.

Pure functions (frames and dicts in, dicts out) plus a thin CLI::

    python validation/pilot/te_real.py compute --a2 OUT_A2 --a1 OUT_A1 \\
        --peptides-a1 PEP_A1 --peptides-a2 PEP_A2 --contaminants CONTAM --out METRICS_DIR
    python validation/pilot/te_real.py te7 --ref OUT_A --rep OUT_B --repo GIAP_REPO \\
        [--diann-a DIANN_A --diann-b DIANN_B --peptides-a PEP_A2 --peptides-b CONV_B] \\
        --out METRICS_DIR

Exit codes: 0 written, whatever the TE statuses; 1 an internal inconsistency; 2 a precondition.
``te_real_exceptions.tsv`` holds row-level data and is never committed.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from collections.abc import Collection, Iterable, Mapping, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from validation.pilot import metrics as pm
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


def _s(v: Any) -> str:
    c = pm.clean(v)
    return "" if c is None else str(c)


def _true(v: Any) -> bool:
    c = pm.clean(v)
    return c is True or (isinstance(c, str) and c.lower() in ("true", "1"))


def _tally(keys: Iterable[str]) -> dict[str, int]:
    c = Counter(keys)
    return {k: c[k] for k in sorted(c)}


def _status(parts: Sequence[str]) -> str:
    """Combined status: fail if any fails, pass if all pass, else not_evaluable."""
    ok = bool(parts) and all(p == "pass" for p in parts)
    return "fail" if "fail" in parts else "pass" if ok else "not_evaluable"


# ---------------------------------------------------------------- TE2


def _te2_bad(variants: pd.DataFrame) -> pd.DataFrame:
    """The ``high`` rows with neither exempt flag (te2_result's violations, all of them)."""
    high = (variants["confidence"].astype(str) == "high").to_numpy(dtype=bool)
    exempt = np.array([bool(TE2_EXEMPT & _flag_names(f)) for f in variants["flags"]], dtype=bool)
    return variants[high & ~exempt]


def _n_alt_bin(v: Any) -> str:
    c = pm.clean(v)
    return "NA" if c is None else ">=3" if int(c) >= 3 else str(int(c))


def te2(
    variants: pd.DataFrame | None, genes: pd.DataFrame, detection_summary: pd.DataFrame | None
) -> dict[str, Any]:
    """TE2: no ``high`` ALT call without ``start_codon``/``stop_codon`` (te2_result), with
    breakdowns of every violation and the ``id_mixups`` cross-check; the exception rows ride
    along under ``_rows`` (popped before writing)."""
    warnings: list[dict[str, Any]] = []
    if variants is None:
        return {"status": "not_evaluable", "reason": ["missing:tables/variants.tsv"],
                "counts": dict.fromkeys(TE2_COUNTS), "violations": [], "vacuous": False,
                "breakdown": {}, "warnings": warnings}  # fmt: skip
    res, bad = te2_result(variants), _te2_bad(variants)
    n_bad = res.counts["violations"]
    if len(bad) != n_bad:  # S5: the breakdown rule must agree with te2_result's
        raise pm.Inconsistent(f"TE2: {len(bad)} breakdown rows, te2_result counts {n_bad}")
    tags = genes["locus_tag"].astype(str)
    rep = dict(zip(tags, genes["replicon_type"].map(_s), strict=True))
    reg = dict(zip(tags, genes["region_ids"], strict=True))
    region: list[str] = []
    for g in bad["gene"].map(_s):
        names = sorted({re.sub(r"_\d+$", "", r) for r in pm.items(reg.get(g))}) or ["none"]
        region += names if g in reg else ["not_in_genes"]
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
    mix = None
    if detection_summary is not None:
        ds = detection_summary
        row = ds[(ds["section"].astype(str) == "qc") & (ds["metric"].astype(str) == "id_mixups")]
        mix = pm.clean(row["value"].iloc[0]) if len(row) else None
    if mix is None:
        warnings.append({"code": "id_mixups_na"})
    elif n_bad > mix:
        warnings.append(
            {"code": "te2_id_mixups_inconsistent", "violations": n_bad, "id_mixups": mix}
        )
    rows = [
        ["TE2", "te2_violation", _s(g), f"variant_id={_s(v)};aa_pos={_s(p)}"]
        for v, g, p in bad[["variant_id", "gene", "aa_pos"]].itertuples(index=False)
    ]
    return {"status": "pass" if n_bad == 0 else "fail", "counts": dict(res.counts),
            "violations": list(res.exceptions), "vacuous": len(variants) == 0,
            "breakdown": breakdown, "id_mixups": mix, "warnings": warnings, "_rows": rows,
            **({"reason": ["violations"]} if n_bad else {})}  # fmt: skip


# ---------------------------------------------------------------- TE3


def _coding(genes: pd.DataFrame) -> pd.DataFrame:
    return genes[genes["feature_type"].astype(str).isin(pm.CODING)].reset_index(drop=True)


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
    det = {"share": pm.share(k, n), "holds": n > 0 and k / n >= TE3_MIN}
    return {"n_ribosomal": n, "abundance": ab, "detection": det}


def te3(genes: pd.DataFrame) -> dict[str, Any]:
    """TE3: R >= 90 % detected (clause 1) and median abundance of R >= the 75th percentile of the
    measured genes (clause 2); secondaries never change the status. ``_rows``: the genes of R not
    detected (popped before writing)."""
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
    miss = G[R & (G["detection_status"].astype(str) != "detected").to_numpy()]
    out["_rows"] = [["TE3", "te3_not_detected", r.locus_tag,
                     f"gene={_s(r.gene)};product={_s(r.product)};detection_status="
                     f"{_s(r.detection_status)};n_observable={_s(r.n_observable_peptides)}"]
                    for r in miss.itertuples(index=False)]  # fmt: skip
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
        ok = ratio is not None and lo <= ratio <= hi  # ratio is set whenever reason is empty
        st = "not_evaluable" if reason else "pass" if ok else "fail"
        pairs[f"{a}/{b}"] = {"status": st, "reason": reason, "low": lo, "high": hi,
                             "log2_ratio": ratio, "top3_log2_ratio": t3}  # fmt: skip
    sec = {
        n: {k: pm.clean(G[k].iloc[ix[0]]) if len(ix) == 1 else None
            for k in ("n_peptides_unique", "n_observable_peptides")}
        for n, ix in rows.items()
    }  # fmt: skip
    status = _status([p["status"] for p in pairs.values()])
    out: dict[str, Any] = {"status": status, "pairs": pairs, "secondary": {"genes": sec},
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
        raise pm.Inconsistent("TE5: family or class counts exceed |S|")
    return {"n_s": len(S), "families": per, "counts": counts}, classes


def te5_arm(
    P: pd.DataFrame,
    peptide_hits: pd.DataFrame,
    unmapped: pd.DataFrame,
    families: Mapping[str, Collection[str]],
    cfg: GiapConfig,
) -> dict[str, Any]:
    """TE5 on one arm: family peptides of **U** (accepted target rows) land in the bucket. The
    violating peptides ride along under ``_violations`` (popped before writing)."""
    U = set(select_queries(P, cfg)["peptide"].astype(str))
    acc = pm._target(P) & (pm._q(P) <= pm.Q_ACCEPT)
    own: dict[str, set[str]] = {}
    for pep, ids, ok in zip(P["peptide"].astype(str), P["protein_ids"], acc, strict=True):
        if ok and pep in U:
            own.setdefault(pep, set()).update(pm.items(ids))
    um = unmapped[unmapped["kind"].astype(str) == "peptide"]
    labels: dict[str, dict[str, set[str]]] = {}
    for sid, lab, m in zip(um["sequence_id"].astype(str), um["label"].astype(str),
                           um["contaminant_match"], strict=True):  # fmt: skip
        labels.setdefault(sid, {}).setdefault(lab, set()).add(_s(m))
    routes = pm._routes(peptide_hits)
    res, classes = _te5_judge(own, routes, labels, families, VERDICT_FAMILIES)
    ka, _ = _te5_judge(own, routes, labels, families, ("keratin_associated",))
    c = res["counts"]
    denom = res["n_s"] - c["genome_exact"]
    bad = sum(c[v] for v in TE5_VIOLATIONS)
    res |= {"share_bucket": pm.share(c["bucket"], denom), "secondary": {"keratin_associated": ka},
            "status": "fail" if bad else "pass" if denom > 0 else "not_evaluable",
            "_violations": {p: k for p, k in classes.items() if k in TE5_VIOLATIONS}}  # fmt: skip
    if res["status"] != "pass":
        res["reason"] = ["violations"] if bad else ["no_trypsin_keratin_peptides"]
    return res


# ---------------------------------------------------------------- CLI (file I/O)


def _spec_sha() -> str:
    res = subprocess.run(["git", "-C", str(pm.REPO), "log", "-1", "--format=%H", "--", SPEC],
                         capture_output=True, text=True, check=False)  # fmt: skip
    sha = res.stdout.strip()
    return sha if re.fullmatch(r"[0-9a-f]{40}", sha) else "0" * 40


def _out_dir(out: Path) -> Path:
    o = out.resolve()
    if o == pm.REPO or pm.REPO in o.parents:
        raise pm.Precondition(f"--out {out} is inside the repository")
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
    man = pm._read_json(out / "run_manifest.json")
    if man.get("exit_code") != 0:
        raise pm.Precondition(f"{out / 'run_manifest.json'}: exit_code {man.get('exit_code')}")
    return man


def _pinned(man: Mapping[str, Any], fasta: Path) -> bool:
    """False when the manifest records no contaminants checksum; raises when it differs or is
    not sha256."""
    db = next((d for d in man.get("databases") or [] if d.get("name") == "contaminants"), None)
    if not db or not db.get("checksum"):
        return False
    algo, _, hexv = str(db["checksum"]).rpartition(":")
    if algo not in ("", "sha256") or pm.sha256(fasta) != hexv:
        raise pm.Precondition(f"{fasta}: not the run's contaminants {db['checksum']}")
    return True


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[list[Any]]]:
    """te_real.json and the exception rows (reads only; raises on a precondition)."""
    from giap.io.peptides import read_peptide_table

    _out_dir(a.out)
    contam = pm._need(a.contaminants)
    arms: dict[str, dict[str, Any]] = {}
    unpinned: list[str] = []
    for label, d, pep in (("A1", a.a1, a.peptides_a1), ("A2", a.a2, a.peptides_a2)):
        out = pm._need(pm.resolve_attempt(d)[0])
        man = _manifest(out)
        if pm._input_sha(man, "peptides") != pm.sha256(pm._need(pep)):
            raise pm.Precondition(f"{pep}: sha256 differs from the run's peptides input ({out})")
        unpinned += [] if _pinned(man, contam) else [label]
        hits, unm = (pm._run_table(out, n) for n in ("peptide_hits.tsv", "unmapped.tsv"))
        arms[label] = {"out": out, "man": man, "pep": pep, "P": read_peptide_table(pep)[0],
                       "hits": hits and pm._table(hits, "peptide_hits"),
                       "unm": unm and pm._table(unm, "unmapped")}  # fmt: skip
    out2 = arms["A2"]["out"]
    gpath, vpath = out2 / "tables" / "genes.tsv", out2 / "tables" / "variants.tsv"
    genes, dpath = pm._table(pm._need(gpath), "genes"), pm._run_table(out2, "detection_summary.tsv")
    variants = pm._table(vpath, "variants") if vpath.is_file() else None
    doc: dict[str, Any] = {"schema": "giap-pilot-te-real/1", "spec": _spec_sha()}
    doc["te2"] = te2(variants, genes, dpath and pm._table(dpath, "detection_summary"))
    doc["te3"], doc["te4"] = te3(genes), te4(genes)
    heads = [ln[1:] for ln in contam.read_text(encoding="utf-8").splitlines() if ln.startswith(">")]
    fam = contaminant_families(heads)
    rows: list[list[Any]] = doc["te2"].pop("_rows", []) + doc["te3"].pop("_rows")
    doc["te5"] = {}
    for label, arm in arms.items():
        if arm["hits"] is None or arm["unm"] is None:
            doc["te5"][label] = {"status": "not_evaluable", "share_bucket": pm.share(0, 0),
                                 "reason": ["missing:peptide_hits.tsv/unmapped.tsv"]}  # fmt: skip
            continue
        doc["te5"][label] = te5_arm(arm["P"], arm["hits"], arm["unm"], fam, pm._cfg(arm["man"]))
        bad = doc["te5"][label].pop("_violations")
        rows += [["TE5", f"te5_{k}", p, label] for p, k in bad.items()]
    doc["te5"]["status"] = _status([doc["te5"]["A1"]["status"], doc["te5"]["A2"]["status"]])
    doc["warnings"] = list(doc["te2"]["warnings"]) + (
        [{"code": "contaminants_unpinned", "arms": unpinned}] if unpinned else []
    )
    files = {"contaminants": contam, "genes_a2": gpath,
             **({"variants_a2": vpath} if variants is not None else {})}  # fmt: skip
    for label, arm in arms.items():
        files |= {f"manifest_{label.lower()}": arm["out"] / "run_manifest.json",
                  f"peptides_{label.lower()}": arm["pep"]}  # fmt: skip
    doc["inputs"] = {k: {"path": str(Path(p).resolve()), "sha256": pm.sha256(Path(p))}
                     for k, p in files.items()}  # fmt: skip
    return doc, rows


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for name in ("a1", "a2", "peptides-a1", "peptides-a2", "contaminants", "out"):
        c.add_argument(f"--{name}", type=Path, required=True)
    t = sub.add_parser("te7")  # TE7 lives in te7.py; the subcommand stays here (spec interface)
    for name in ("ref", "rep", "repo", "out"):
        t.add_argument(f"--{name}", type=Path, required=True)
    for name in ("diann-a", "diann-b", "peptides-a", "peptides-b"):
        t.add_argument(f"--{name}", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "te7":
            from validation.pilot import te7

            st = te7.run(a)
        else:
            doc, rows = compute(a)
            cols = ["te", "kind", "id", "detail"]
            _write(a.out, "te_real.json", doc, "te_real_exceptions.tsv", cols, rows)
            st = {k: doc[k]["status"] for k in ("te2", "te3", "te4", "te5")}
        print(f"wrote {a.out}: " + ", ".join(f"{k} {v}" for k, v in st.items()))
        return 0
    except pm.Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except (pm.Inconsistent, ValueError) as err:
        print(f"inconsistent: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
