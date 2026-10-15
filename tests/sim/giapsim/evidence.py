"""Peptide evidence for sequence-free layouts (issue M4-09; GIAP_SPEC.md §11.3, §8.4, §8.7, §7.3.2).

The evidence half of the M4-06 contract (C9, C10, C11-C13) drawn from the M4-08 expression truth as
if preflight and mapping were perfect. True peptides: Trypsin/P, 7-30 aa, <= 1 missed cleavage,
observed with probability ``sigmoid(A + B (abundance_log2_true - 20) + N(0, 1))`` (>= 1 per
digestible detected gene), ``ln I = abundance_log2_true ln 2 + N(0, 0.5^2)``, score ~ N(mu1, 1).
False targets (``round(pi/(1 - pi) n_true)``) are unobserved tryptic peptides of genes drawn in
proportion to ``n_observable``, or unplaceable random sequences; decoys are reversed tryptic
peptides (C-terminal K/R kept) of genes drawn the same way; both score ~ N(mu0, 1). PEP is the
mixture posterior with ``pi0 = pi_wrong``; q(s) is the mean target PEP at scores >= s. ``emit``
only masks values (decoy rows are dropped without ``is_decoy``). Nothing is imported from ``giap``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

__all__ = ["EMIT_FIELDS", "Evidence", "EvidenceConfig", "simulate_evidence",
           "write_evidence_tables"]  # fmt: skip

EMIT_FIELDS = frozenset({"q_value", "pep", "score", "is_decoy", "intensity", "protein_ids"})
_CUT = re.compile(r"(?<=[KR])")
_LEN = (7, 30)
_NON_KR = np.array(list("ACDEFGHILMNPQSTVWY"))
_LOCUS_CLASS = {"cds": "cds_inframe", "sorf": "sorf", "pseudogene": "pseudogene"}
_A, _B, _AB_REF = -0.5, 0.6, 20.0  # peptide detectability logit (log2 abundance centred at 20)
_LN2 = math.log(2.0)

#: (column, logical dtype) of C9-C13 in contract order (literal copy of M4-06 ``TABLES``).
_TYPES: dict[str, tuple[tuple[str, str], ...]] = {
    "peptides": (
        ("peptide", "str"), ("q_value", "float"), ("pep", "float"), ("score", "float"),
        ("is_decoy", "bool"), ("is_contaminant", "bool"), ("protein_ids", "list"),
        ("modified_peptides", "list"), ("n_rows", "int"), ("n_runs", "int"),
        ("intensity", "float"), ("source", "str"), ("search_db", "str"),
    ),
    "abundance": (("protein_id", "str"), ("value", "float")),
    "peptide_hits": (
        ("peptide", "str"), ("peptide_len", "int"), ("source", "str"), ("best_q_input", "float"),
        ("pep_input", "float"), ("score_input", "float"), ("contig", "str"), ("start", "int"),
        ("end", "int"), ("strand", "str"), ("frame", "int"), ("crosses_origin", "bool"),
        ("route", "str"), ("n_diff", "int"), ("diff_types", "list"), ("placement_q", "float"),
        ("length_class", "str"), ("locus_class", "str"), ("uniqueness", "str"), ("gene", "str"),
        ("novel_candidate_id", "str"), ("flags", "list"),
    ),
    "protein_hits": (
        ("protein_id", "str"), ("protein_len_aa", "int"), ("route", "str"), ("contig", "str"),
        ("start", "int"), ("end", "int"), ("strand", "str"), ("rank", "int"),
        ("identity", "float"), ("query_coverage", "float"), ("reciprocal_overlap", "float"),
        ("class", "str"), ("gene", "str"), ("flags", "list"),
    ),
    "unmapped": (
        ("sequence_id", "str"), ("sequence", "str"), ("kind", "str"), ("label", "str"),
        ("reason", "str"), ("contaminant_match", "str"),
    ),
}  # fmt: skip
_ARROW = {"str": pa.string(), "int": pa.int64(), "float": pa.float64(), "bool": pa.bool_(),
          "list": pa.list_(pa.string())}  # fmt: skip


@dataclass(frozen=True)
class EvidenceConfig:
    pi_wrong: float = 0.01
    false_source: Literal["random_gene", "unplaceable"] = "random_gene"
    n_decoys: Literal["matched"] | int = "matched"
    mixture: tuple[float, float] = (0.0, 3.0)
    emit: frozenset[str] = frozenset({"q_value", "pep", "is_decoy", "intensity", "protein_ids"})
    protein_ids: Literal["locus_tag", "foreign"] = "locus_tag"
    multi_identical_pairs: int = 0
    shared_frac: float = 0.0
    tolerant_frac: float = 0.0
    proteins_only: bool = False
    abundance_table: Literal["none", "ibaq", "maxlfq"] = "none"

    def __post_init__(self) -> None:
        bad = []
        if not 0.0 <= self.pi_wrong < 1.0:
            bad.append(f"pi_wrong {self.pi_wrong} not in [0, 1)")
        for name in ("shared_frac", "tolerant_frac"):
            if not 0.0 <= getattr(self, name) <= 1.0:
                bad.append(f"{name} not in [0, 1]")
        if self.false_source not in ("random_gene", "unplaceable"):
            bad.append(f"false_source {self.false_source!r}")
        if self.protein_ids not in ("locus_tag", "foreign"):
            bad.append(f"protein_ids {self.protein_ids!r}")
        if self.abundance_table not in ("none", "ibaq", "maxlfq"):
            bad.append(f"abundance_table {self.abundance_table!r}")
        nd = self.n_decoys
        if nd != "matched" and (isinstance(nd, bool) or not isinstance(nd, int) or nd < 0):
            bad.append(f"n_decoys {nd!r}")
        if self.multi_identical_pairs < 0 or not set(self.emit) <= EMIT_FIELDS:
            bad.append("multi_identical_pairs < 0 or unknown emit field")
        if bad:
            raise ValueError("invalid EvidenceConfig: " + "; ".join(bad))


@dataclass(frozen=True, eq=False)
class Evidence:
    peptides: pd.DataFrame
    peptide_hits: pd.DataFrame
    protein_hits: pd.DataFrame | None
    abundance: pd.DataFrame | None
    truth: pd.DataFrame
    gene_truth: pd.DataFrame


def _digest(aa: object) -> tuple[list[tuple[str, int]], int]:
    """(distinct 7-30 aa Trypsin/P peptides with <= 1 missed cleavage and their first residue
    offset, in sequence order; n_observable = distinct 0-missed peptides)."""
    if not isinstance(aa, str) or not aa:
        return [], 0
    frags = [f for f in _CUT.split(aa) if f]
    pos = [0]
    for f in frags:
        pos.append(pos[-1] + len(f))
    out: dict[str, int] = {}
    zero: set[str] = set()
    for i in range(len(frags)):
        for j in (i + 1, i + 2):
            if j <= len(frags) and _LEN[0] <= pos[j] - pos[i] <= _LEN[1]:
                p = aa[pos[i] : pos[j]]
                out.setdefault(p, pos[i])
                if j == i + 1:
                    zero.add(p)
    return list(out.items()), len(zero)


def _reverse(p: str) -> str:
    return p[:-1][::-1] + p[-1] if p[-1] in "KR" else p[::-1]


def pep_q(scores: np.ndarray, is_target: np.ndarray, pi0: float,
          mixture: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """PEP = pi0 f0 / (pi0 f0 + (1 - pi0) f1) and q(s) = mean target PEP at scores >= s."""
    mu0, mu1 = mixture
    if pi0 <= 0:
        pep = np.zeros(len(scores))
    else:
        l0 = math.log(pi0) - 0.5 * (scores - mu0) ** 2
        l1 = math.log1p(-pi0) - 0.5 * (scores - mu1) ** 2
        pep = np.exp(l0 - np.logaddexp(l0, l1))
    ts = scores[is_target]
    if len(ts) == 0:
        return pep, pep.copy()
    cum = np.cumsum(pep[is_target][np.argsort(-ts, kind="stable")])
    k = np.maximum(len(ts) - np.searchsorted(np.sort(ts), scores, side="left"), 1)
    return pep, cum[k - 1] / k


def _span(g: dict[str, Any], a: int, b: int, length: int) -> tuple[int, int, bool]:
    """0-based half-open genomic span of residues [a, b) of gene row ``g``."""
    s, e = ((g["start"] + 3 * a, g["start"] + 3 * b) if g["strand"] == "+"
            else (g["end"] - 3 * b, g["end"] - 3 * a))  # fmt: skip
    if s >= length:
        s, e = s - length, e - length
    return int(s), int(e), bool(s < length < e)


def _frame(key: str, data: dict[str, Any]) -> pd.DataFrame:
    cols: dict[str, Any] = {}
    for name, dtype in _TYPES[key]:
        v = data[name]
        if dtype == "int" and any(x is None for x in v):
            cols[name] = pd.array(v, dtype="Int64")
        elif dtype in ("int", "float", "bool") and len(v):
            cols[name] = np.asarray(v, dtype={"int": np.int64, "float": float, "bool": bool}[dtype])
        else:
            cols[name] = pd.Series(list(v), dtype=object)
    return pd.DataFrame(cols, columns=[n for n, _ in _TYPES[key]])


def simulate_evidence(layout: Any, expr: Any, cfg: EvidenceConfig,
                      rng: np.random.Generator) -> Evidence:  # fmt: skip
    """Evidence tables for a ``giapsim.layout.Layout`` and its ``Expression`` (not mutated)."""
    genes = layout.genes.to_dict("records")
    tags = [str(g["locus_tag"]) for g in genes]
    n = len(tags)
    ex = expr.genes.set_index(expr.genes["locus_tag"].astype(str))
    ab = ex.loc[tags, "abundance_log2_true"].to_numpy(float)
    det = ex.loc[tags, "detected_true"].to_numpy(bool)
    length = dict(zip(layout.replicons["contig"], layout.replicons["length"], strict=True))
    dig = [_digest(g["aa"]) for g in genes]
    cands = [d[0] for d in dig]
    w = np.array([d[1] for d in dig], dtype=float)
    digestible = np.array([len(c) > 0 for c in cands])
    ids = tags if cfg.protein_ids == "locus_tag" else [f"UP_{100001 + i}" for i in range(n)]
    present = det & digestible
    aa_len = np.array([len(g["aa"]) if isinstance(g["aa"], str) else 0 for g in genes])
    # (peptide, gene index or -1, residue offset, kind) with kind true / false / decoy
    rows: list[tuple[str, int, int, str]] = []
    hits: list[tuple[int, int, int, str, str]] = []  # (row, gene, offset, route, uniqueness)
    partner: dict[int, int] = {}  # multi_identical pair: gene -> copy
    protein_hits = None
    if cfg.proteins_only:
        protein_hits = _frame("protein_hits", _protein_rows(genes, ids, np.flatnonzero(present)))
        rows_out = _rows_frames([], [], np.zeros(0), np.zeros(0), ids, genes, length, cfg)
    else:
        nc = np.array([len(c) for c in cands], dtype=np.int64)
        gi = np.repeat(np.arange(n), nc)
        logit = _A + _B * (ab[gi] - _AB_REF) + rng.normal(0.0, 1.0, len(gi))
        obs = rng.random(len(gi)) < 1.0 / (1.0 + np.exp(-logit))
        used: set[str] = set()
        off = np.concatenate([[0], np.cumsum(nc)])
        for i in np.flatnonzero(present).tolist():
            o = off[i]
            sel = [j for j in range(nc[i]) if obs[o + j]] or [int(np.argmax(logit[o : off[i + 1]]))]
            pick = [cands[i][j] for j in sel if cands[i][j][0] not in used]
            pick = pick or [c for c in cands[i] if c[0] not in used][:1]
            for p, a in pick:
                used.add(p)
                rows.append((p, i, a, "true"))
        n_true = len(rows)
        s_true, ln_true = rng.normal(cfg.mixture[1], 1.0, n_true), rng.normal(0.0, 0.5, n_true)
        ln_true += np.array([ab[r[1]] * _LN2 for r in rows])
        has_true = sorted({r[1] for r in rows})
        for _ in range(cfg.multi_identical_pairs):  # duplicated gene pairs (IS copies)
            free = [i for i in has_true if i not in partner and i not in partner.values()]
            if not free:
                break
            a = free[int(rng.integers(len(free)))]
            la = aa_len[a]
            taken = set(partner) | set(partner.values())
            bs = [j for j in np.flatnonzero(aa_len >= la).tolist() if j != a and j not in taken]
            if bs:
                partner[a] = bs[int(rng.integers(len(bs)))]
        u_tol, u_sh = rng.random(n_true), rng.random(n_true)
        for k, (p, i, a, _) in enumerate(list(rows)):
            if i in partner:
                hits.append((k, i, a, "exact", "multi_identical"))
                hits.append((k, partner[i], a, "exact", "multi_identical"))
                continue
            if u_tol[k] < cfg.tolerant_frac:
                x = int(rng.integers(len(p) - 1))
                mut = p[:x] + str(rng.choice(_NON_KR[p[x] != _NON_KR])) + p[x + 1 :]
                if mut not in used and mut not in genes[i]["aa"]:
                    used.add(mut)
                    rows[k] = (mut, i, a, "true")
                    hits.append((k, i, a, "tolerant", "unique"))
                    continue
            if u_sh[k] < cfg.shared_frac:
                js = np.flatnonzero(aa_len >= len(p))
                js = js[js != i]
                if len(js):
                    j = int(js[int(rng.integers(len(js)))])
                    b = int(rng.integers(len(genes[j]["aa"]) - len(p) + 1))
                    hits.append((k, i, a, "exact", "shared"))
                    hits.append((k, j, b, "exact", "shared"))
                    continue
            hits.append((k, i, a, "exact", "unique"))
        n_false = round(cfg.pi_wrong / (1.0 - cfg.pi_wrong) * n_true)
        _false_targets(rows, hits, used, cands, w, n_false, cfg, genes, rng)
        n_dec = n_false if cfg.n_decoys == "matched" else int(cfg.n_decoys)
        for _ in range(n_dec):
            for _try in range(1000):
                i = int(rng.choice(n, p=w / w.sum()))
                d = _reverse(cands[i][int(rng.integers(len(cands[i])))][0])
                if d not in used:
                    used.add(d)
                    rows.append((d, i, 0, "decoy"))
                    break
            else:
                raise ValueError("no free decoy sequence")
        n_null = len(rows) - n_true
        scores = np.concatenate([s_true, rng.normal(cfg.mixture[0], 1.0, n_null)])
        ln_int = np.concatenate([ln_true, (_AB_REF - 2) * _LN2 + rng.normal(0.0, 1.0, n_null)])
        rows_out = _rows_frames(rows, hits, scores, ln_int, ids, genes, length, cfg)
    peptides, peptide_hits, truth = rows_out
    n_tp = np.zeros(n, dtype=np.int64)
    for _, i, _, kind in rows:
        n_tp[i] += kind == "true"
    gene_truth = pd.DataFrame({"locus_tag": tags, "present": n_tp >= 1, "n_true_peptides": n_tp,
                               "detected_true": det & digestible})  # fmt: skip
    abundance = None
    if cfg.abundance_table != "none":
        idx = np.flatnonzero(present)
        vals = ab[idx] + rng.normal(0.0, 0.3, len(idx))
        abundance = _frame("abundance", {"protein_id": [ids[i] for i in idx], "value": vals})
    return Evidence(peptides, peptide_hits, protein_hits, abundance, truth, gene_truth)


def _false_targets(rows: list[tuple[str, int, int, str]], hits: list[Any], used: set[str],
                   cands: list[list[tuple[str, int]]], w: np.ndarray, n_false: int,
                   cfg: EvidenceConfig, genes: list[dict[str, Any]],
                   rng: np.random.Generator) -> None:  # fmt: skip
    """Append false target rows: unobserved tryptic peptides of genes drawn ∝ n_observable
    (placed, unique) or unplaceable random tryptic-like sequences."""
    allaa = "|".join(g["aa"] for g in genes if isinstance(g["aa"], str))
    for _ in range(n_false):
        for _try in range(10_000):
            if cfg.false_source == "random_gene":
                i = int(rng.choice(len(w), p=w / w.sum()))
                free = [c for c in cands[i] if c[0] not in used]
                if not free:
                    continue
                p, a = free[int(rng.integers(len(free)))]
                hits.append((len(rows), i, a, "exact", "unique"))
            else:
                body = "".join(rng.choice(_NON_KR, int(rng.integers(9, 20))))
                p, i, a = body + str(rng.choice(["K", "R"])), -1, 0
                if p in used or p in allaa:
                    continue
            used.add(p)
            rows.append((p, i, a, "false"))
            break
        else:
            raise ValueError("no free false-target sequence")


def _rows_frames(rows: list[tuple[str, int, int, str]], hits: list[tuple[int, int, int, str, str]],
                 scores: np.ndarray, ln_int: np.ndarray,
                 ids: list[str], genes: list[dict[str, Any]], length: dict[str, int],
                 cfg: EvidenceConfig) -> tuple[pd.DataFrame, ...]:  # fmt: skip
    """C9, C11 (internal 0-based) and truth frames from the drawn rows."""
    is_t = np.array([r[3] != "decoy" for r in rows], dtype=bool)
    pep, q = pep_q(scores, is_t, cfg.pi_wrong, cfg.mixture)
    keep = [k for k, r in enumerate(rows) if r[3] != "decoy" or "is_decoy" in cfg.emit]
    gid: dict[int, list[int]] = {}
    for k, i, *_ in hits:
        gid.setdefault(k, []).append(i)
    ids_of: list[list[str]] = []
    for k, (_, i, _, kind) in enumerate(rows):
        if "protein_ids" not in cfg.emit:
            ids_of.append([])
        elif kind == "decoy":
            ids_of.append(["rev_" + ids[i]])
        else:
            ids_of.append([ids[j] for j in gid[k]] if k in gid else [f"OTHER_{k + 1}"])
    nan = np.full(len(rows), np.nan)
    col = {"q_value": q, "pep": pep, "score": scores, "intensity": np.exp(ln_int)}
    col = {c: (v if c in cfg.emit else nan) for c, v in col.items()}
    m = len(keep)
    pdata: dict[str, Any] = {
        "peptide": [rows[k][0] for k in keep], "is_decoy": [rows[k][3] == "decoy" for k in keep],
        "is_contaminant": [False] * m, "protein_ids": [ids_of[k] for k in keep],
        "modified_peptides": [[] for _ in keep], "n_rows": [1] * m, "n_runs": [1] * m,
        "source": ["giapsim"] * m, "search_db": [None] * m,
        **{c: v[keep] for c, v in col.items()},
    }  # fmt: skip
    h: dict[str, list[Any]] = {c: [] for c, _ in _TYPES["peptide_hits"]}
    for k, i, a, route, uniq in hits:
        g, p = genes[i], rows[k][0]
        s, e, cross = _span(g, a, a + len(p), int(length[g["contig"]]))
        tol = route == "tolerant"
        for c, v in (("peptide", p), ("peptide_len", len(p)), ("source", "giapsim"),
                     ("best_q_input", col["q_value"][k]), ("pep_input", col["pep"][k]),
                     ("score_input", col["score"][k]), ("contig", g["contig"]), ("start", s),
                     ("end", e), ("strand", g["strand"]), ("frame", None),
                     ("crosses_origin", cross), ("route", route), ("n_diff", int(tol)),
                     ("diff_types", ["substitution"] if tol else []), ("placement_q", 0.0),
                     ("length_class", None), ("locus_class", _LOCUS_CLASS[g["feature_type"]]),
                     ("uniqueness", uniq), ("gene", g["locus_tag"]), ("novel_candidate_id", None),
                     ("flags", ["strain_difference"] if tol else [])):  # fmt: skip
            h[c].append(v)
    truth = pd.DataFrame({
        "peptide": [r[0] for r in rows], "is_false": [r[3] != "true" for r in rows],
        "gene_true": pd.Series([genes[r[1]]["locus_tag"] if r[3] == "true" else None
                                for r in rows], dtype=object),
        "decoy_gene": pd.Series([genes[r[1]]["locus_tag"] if r[3] == "decoy" else None
                                 for r in rows], dtype=object),
    })  # fmt: skip
    return _frame("peptides", pdata), _frame("peptide_hits", h), truth


def _protein_rows(genes: list[dict[str, Any]], ids: list[str], idx: np.ndarray) -> dict[str, Any]:
    """C12 rows (0-based): one identical rank-1 hit per detected digestible gene."""
    rs = [genes[i] for i in idx]
    k = len(rs)
    return {"protein_id": [ids[i] for i in idx], "protein_len_aa": [len(g["aa"]) for g in rs],
            "route": ["exact"] * k, "contig": [g["contig"] for g in rs],
            "start": [int(g["start"]) for g in rs], "end": [int(g["end"]) for g in rs],
            "strand": [g["strand"] for g in rs], "rank": [1] * k, "identity": [1.0] * k,
            "query_coverage": [1.0] * k, "reciprocal_overlap": [1.0] * k,
            "class": ["identical"] * k, "gene": [g["locus_tag"] for g in rs],
            "flags": [[] for _ in rs]}  # fmt: skip


def _cell(v: Any, dtype: str) -> str:
    if dtype == "list":
        return ";".join(str(x) for x in v) if len(v) else "NA"
    if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    if dtype == "bool":
        return "true" if v else "false"
    return repr(float(v)) if dtype == "float" else str(int(v)) if dtype == "int" else str(v)


def _write_tsv(path: Path, key: str, df: pd.DataFrame, one_based: bool) -> None:
    types = _TYPES[key]
    lines = ["\t".join(c for c, _ in types)]
    for r in df.to_dict("records"):
        if one_based:
            r["start"] = int(r["start"]) + 1
        lines.append("\t".join(_cell(r[c], t) for c, t in types))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_parquet(path: Path, key: str, df: pd.DataFrame) -> None:
    arrays = []
    for c, t in _TYPES[key]:
        v = df[c].tolist()
        if t == "list":
            v = [[str(x) for x in cell] for cell in v]
        else:
            v = [None if _cell(x, "str") == "NA" else x for x in v]
        arrays.append(pa.array(v, type=_ARROW[t]))
    schema = pa.schema([pa.field(c, _ARROW[t]) for c, t in _TYPES[key]])
    pq.write_table(pa.Table.from_arrays(arrays, schema=schema), path, compression="zstd")


def write_evidence_tables(ev: Evidence, work_dir: Path) -> None:
    """Write C9 (``peptides.norm.parquet``), C10 (``abundance.norm.parquet``, if any) under
    ``work_dir`` and C11-C13 (1-based TSVs, header-only when empty) under ``work_dir/mapping``."""
    mapping = work_dir / "mapping"
    mapping.mkdir(parents=True, exist_ok=True)
    _write_parquet(work_dir / "peptides.norm.parquet", "peptides", ev.peptides)
    if ev.abundance is not None:
        _write_parquet(work_dir / "abundance.norm.parquet", "abundance", ev.abundance)
    _write_tsv(mapping / "peptide_hits.tsv", "peptide_hits", ev.peptide_hits, True)
    ph = ev.protein_hits if ev.protein_hits is not None else _frame(
        "protein_hits", {c: [] for c, _ in _TYPES["protein_hits"]})  # fmt: skip
    _write_tsv(mapping / "protein_hits.tsv", "protein_hits", ph, True)
    p = ev.peptides
    placed = set(ev.peptide_hits["peptide"])
    seqs = [s for s, d in zip(p["peptide"], p["is_decoy"], strict=True)
            if not d and s not in placed]  # fmt: skip
    k = len(seqs)
    um = _frame("unmapped", {"sequence_id": [f"peptide_{i + 1}" for i in range(k)],
                             "sequence": seqs, "kind": ["peptide"] * k, "label": [None] * k,
                             "reason": ["no_match"] * k,
                             "contaminant_match": [None] * k})  # fmt: skip
    _write_tsv(mapping / "unmapped.tsv", "unmapped", um, False)
