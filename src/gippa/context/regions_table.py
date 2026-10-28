"""``regions_stats`` rows: plasmid replicon-level effects, prophage subsets, the table (M4-18).

GIPPA_SPEC.md §7.3.5(b), §8.5 attack 4, §9.2, §9.3; docs/specs/M4-18-regions-stats-table.md (lead
decisions 1-23); docs/methods/region-tests.md "Replicon-level effects" and "Descriptive rows".

- Plasmid rows: the plasmid's own O - E (E = sum of ``p_detect``) and the Hodges-Lehmann shift of
  the plasmid's raw log2 abundance against the reference replicons (``kind = chromosome``; else
  every replicon that is neither ``kind = plasmid`` nor a plasmid region's contig, flag
  ``no_reference_replicon``), imputed genes excluded. TU-block bootstrap CIs (HL: strata plasmid /
  reference); never a p-value (``null_model = none``). With depth, an extra ``hl_shift`` row
  minus log2(plasmid depth / length-weighted reference depth), depth treated as known.
- Subset rows: detected / total per prophage function subset, descriptive (no CI, no p).
- ``regions_stats``: M4-17 rows (``subset = all``, q-values copied, never recomputed), plasmid
  and subset rows, C7 coordinates, sorted.

Pure functions; randomness only from the generator passed in, from which each plasmid region gets
its own child stream ``for_step("context_stats.<region_id>", base)``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from gippa.config.model import GippaConfig
from gippa.context.regions import SMALL_REGION
from gippa.genetable import eligible
from gippa.rng import for_step
from gippa.stats.bootstrap import tu_block_bootstrap
from gippa.stats.robust import hodges_lehmann

__all__ = ["ADJUSTED_LABEL", "PLASMID_LABEL", "PLASMID_METHOD", "REGIONS_STATS_COLUMNS", "STATS",
           "SUBSETS", "depth_gaps", "plasmid_rows", "reference_contigs", "regions_stats",
           "subset_rows"]  # fmt: skip

#: §9.3 order (region, detection, abundance), the §9.2 test columns, then ``label``, ``flags``
REGIONS_STATS_COLUMNS: tuple[str, ...] = (
    "region_id", "region_type", "subset", "contig", "start", "end", "n_genes", "n_eligible",
    "observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model", "n_perm",
    "label", "flags",
)  # fmt: skip
STATS = ("o_minus_e", "hl_shift")  # schema enum order = row order within a region
SUBSETS = ("integrase", "repressor", "structural", "lysis", "cargo")  # §7.1.7
TYPE_ORDER = ("prophage", "plasmid", "bgc", "crispr", "amr")  # default context.region_types
PLASMID_METHOD = "replicon_contrast"
PLASMID_LABEL = "dosage_plus_regulation"
ADJUSTED_LABEL = "dosage_adjusted"
DESCRIPTIVE = "descriptive"
NO_REFERENCE = "no_reference_replicon"
DEPTH_FLAG = "depth_adjusted"
_FLOATS = ("observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2", "effect", "ci_low",
           "ci_high", "p_value", "q_value")  # fmt: skip
_INTS = ("n_genes", "n_eligible", "n_perm")
_NA = math.nan
B8 = npt.NDArray[np.bool_]


def _frame(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(REGIONS_STATS_COLUMNS))
    df = df.astype({**dict.fromkeys(_FLOATS, "float64"), **dict.fromkeys(_INTS, "int64"),
                    "start": "Int64", "end": "Int64"})  # fmt: skip
    for c in ("region_id", "region_type", "subset", "contig", "stat", "method", "null_model",
              "label", "flags"):  # fmt: skip
        df[c] = df[c].astype(object)
    return df


def _bool(genes: pd.DataFrame, col: str) -> B8:
    if col not in genes:
        return np.zeros(len(genes), dtype=bool)
    return np.asarray(genes[col].astype("boolean").fillna(False).to_numpy(bool), dtype=bool)


def _ids(v: Any) -> list[str]:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return []
    if isinstance(v, str):
        return [x for x in v.split(";") if x]
    return [str(x) for x in v]


def _plasmid_members(genes: pd.DataFrame) -> dict[str, B8]:
    """``plasmid_*`` region id (§6.5 naming, lead decision 10) -> member mask over ``genes``."""
    ids = [_ids(v) for v in genes["region_ids"]]
    out: dict[str, B8] = {}
    for rid in sorted({r for rs in ids for r in rs if r.startswith("plasmid_")}):
        out[rid] = np.array([rid in rs for rs in ids], dtype=bool)
    return out


def reference_contigs(
    replicons: pd.DataFrame, plasmid_contigs: Iterable[str]
) -> tuple[list[str], bool]:
    """(reference replicons, fallback used), lead decision 3: the ``kind = chromosome``
    replicons, else every replicon neither ``kind = plasmid`` nor in ``plasmid_contigs``."""
    contig = replicons["contig"].astype(str).tolist()
    kind = replicons["kind"].astype(str).tolist()
    chrom = [c for c, k in zip(contig, kind, strict=True) if k == "chromosome"]
    if chrom:
        return chrom, False
    pl = set(plasmid_contigs)
    return [c for c, k in zip(contig, kind, strict=True) if k != "plasmid" and c not in pl], True


def _plasmid_contigs(genes: pd.DataFrame, members: dict[str, B8]) -> dict[str, str]:
    contig = genes["contig"].astype(str).to_numpy()
    return {rid: str(pd.Series(contig[m]).mode().iloc[0]) for rid, m in members.items()}


def depth_gaps(
    genes: pd.DataFrame, replicons: pd.DataFrame, depth: pd.DataFrame, cfg: GippaConfig
) -> list[str]:
    """Plasmid-region and reference contigs absent from ``depth`` (their adjusted rows are
    omitted; the step logs a warning, lead decision 14); none when ``plasmid`` is not in
    ``context.region_types`` (decision 22)."""
    if "plasmid" not in cfg.context.region_types:
        return []
    pcontig = _plasmid_contigs(genes, _plasmid_members(genes))
    if not pcontig:
        return []
    ref, _ = reference_contigs(replicons, pcontig.values())
    have = set(depth["contig"].astype(str))
    return sorted({*pcontig.values(), *ref} - have)


def _base(rid: str, contig: str, n_genes: int, m: int) -> dict[str, Any]:
    return {"region_id": rid, "region_type": "plasmid", "subset": "all", "contig": contig,
            "start": None, "end": None, "n_genes": n_genes, "n_eligible": m, "p_value": _NA,
            "q_value": _NA, "method": PLASMID_METHOD, "null_model": "none", "n_perm": 0,
            "label": PLASMID_LABEL}  # fmt: skip


def _detection(o: float, e: float, hl: float) -> dict[str, float]:
    return {"observed": o, "expected": e, "o_minus_e": o - e,
            "o_over_e": o / e if e > 0 else _NA, "hl_shift_log2": hl}  # fmt: skip


def _hl_of(t: pd.DataFrame) -> float:
    a, w = t["a"].to_numpy(float), t["in"].to_numpy(float) > 0.5
    return hodges_lehmann(a[w], a[~w])


def _sum_r(t: pd.DataFrame) -> float:
    return float(t["r"].to_numpy().sum())


def _one_plasmid(
    genes: pd.DataFrame, mask: B8, rid: str, contig: str, ref: list[str], fallback: bool,
    depth: dict[str, float] | None, length: dict[str, float], B: int, rng: np.random.Generator,
) -> list[dict[str, Any]]:  # fmt: skip
    el = eligible(genes)
    pm = mask & el
    m = int(pm.sum())
    base = _base(rid, contig, int(mask.sum()), m)
    ref_flag = (NO_REFERENCE,) if fallback else ()
    if m < SMALL_REGION:  # lead decisions 5 and 13: NA statistics
        f = ("small_region", PLASMID_LABEL)
        na = {**_detection(_NA, _NA, _NA), "effect": _NA, "ci_low": _NA, "ci_high": _NA}
        return [{**base, **na, "stat": "o_minus_e", "flags": f},
                {**base, **na, "stat": "hl_shift", "flags": (*f, *ref_flag)}]  # fmt: skip
    d = _bool(genes, "detected")
    p = pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float)
    contig_col = genes["contig"].astype(str).to_numpy()
    tu = pd.factorize(pd.Series(list(zip(contig_col, genes["tu_id"].astype(str), strict=True))))[0]
    tu = tu.astype(np.float64)
    used = pm | np.isin(contig_col, ref)  # lead decision 22: plasmid and reference genes
    if genes.loc[used, "tu_id"].isna().any():
        raise ValueError(f"missing tu_id on plasmid region {rid!r} or its reference replicons "
                         "(C6: every gene has a tu_id)")  # fmt: skip
    r = d.astype(float) - p
    ci = tu_block_bootstrap(pd.DataFrame({"tu": tu[pm], "r": r[pm]}), "tu", _sum_r, B, rng)
    o, e = float(d[pm].sum()), float(p[pm].sum())

    a = pd.to_numeric(genes["abundance_log2"], errors="coerce").to_numpy(float)
    q = el & d & ~_bool(genes, "abundance_imputed") & np.isfinite(a)  # lead decision 2
    inside, outside = pm & q, np.isin(contig_col, ref) & ~mask & q
    hl = hodges_lehmann(a[inside], a[outside])
    lo = hi = _NA
    if not math.isnan(hl):
        both = inside | outside
        hd = pd.DataFrame({"tu": tu[both], "a": a[both], "in": inside[both].astype(np.float64)})
        hci = tu_block_bootstrap(hd, "tu", _hl_of, B, rng, strata_col="in")
        lo, hi = hci.ci_low, hci.ci_high
    det = _detection(o, e, hl)
    rows = [
        {**base, **det, "stat": "o_minus_e", "effect": o - e, "ci_low": ci.ci_low,
         "ci_high": ci.ci_high, "flags": (PLASMID_LABEL,)},
        {**base, **det, "stat": "hl_shift", "effect": hl, "ci_low": lo, "ci_high": hi,
         "flags": (PLASMID_LABEL, *ref_flag)},
    ]  # fmt: skip
    if depth is not None and not math.isnan(hl) and ref and {contig, *ref} <= set(depth):
        dc = sum(length[c] * depth[c] for c in ref) / sum(length[c] for c in ref)
        off = math.log2(depth[contig] / dc)  # depth treated as known: CI shifted by the constant
        rows.append({**base, **det, "stat": "hl_shift", "effect": hl - off, "ci_low": lo - off,
                     "ci_high": hi - off, "label": ADJUSTED_LABEL,
                     "flags": (DEPTH_FLAG, *ref_flag)})  # fmt: skip
    return rows


def plasmid_rows(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    depth: pd.DataFrame | None,
    cfg: GippaConfig,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Replicon-level rows of every ``plasmid_*`` region (found from the genes' ``region_ids``):
    ``o_minus_e`` and ``hl_shift`` (+ a ``dosage_adjusted`` ``hl_shift`` row when ``depth``
    covers the plasmid and every reference replicon and HL is defined). No row is tested.
    ``depth``: columns ``contig``, ``mean_depth`` (``gippa.io.depth.read_replicon_depth``)."""
    if "plasmid" not in cfg.context.region_types:
        return _frame([])
    members = _plasmid_members(genes)
    pcontig = _plasmid_contigs(genes, members)
    ref, fallback = reference_contigs(replicons, pcontig.values())
    length = dict(zip(replicons["contig"].astype(str), replicons["length"].astype(float),
                      strict=True))  # fmt: skip
    dmap = None
    if depth is not None:
        dmap = dict(zip(depth["contig"].astype(str), depth["mean_depth"].astype(float),
                        strict=True))  # fmt: skip
    base = int(rng.integers(0, 2**63 - 1))  # one draw: plasmid streams independent of their number
    rows: list[dict[str, Any]] = []
    for rid, mask in members.items():
        child = for_step(f"context_stats.{rid}", base)  # lead decision 7
        rows += _one_plasmid(genes, mask, rid, pcontig[rid], ref, fallback, dmap, length,
                             cfg.stats.bootstrap, child)  # fmt: skip
    return _frame(rows)


def subset_rows(
    genes: pd.DataFrame, regions: pd.DataFrame, region_genes: pd.DataFrame
) -> pd.DataFrame:
    """Descriptive detected / total rows per prophage and §7.1.7 subset present (C8 ``subset``):
    ``stat = o_minus_e``, ``effect = o_minus_e``, CI / p / q NA, ``null_model = none``; 0 eligible
    genes give NA statistics, fewer than ``SMALL_REGION`` add ``small_region``."""
    pro = regions[regions["region_type"].astype(str) == "prophage"]
    contig_of = dict(zip(pro["region_id"].astype(str), pro["contig"].astype(str), strict=True))
    rg = region_genes[region_genes["region_id"].astype(str).isin(contig_of)]
    rg = rg[rg["subset"].notna()]
    g = genes.assign(
        _el=eligible(genes),
        _d=_bool(genes, "detected"),
        _p=pd.to_numeric(genes["p_detect"], errors="coerce").astype(float),
    )
    g = g.drop_duplicates("locus_tag")
    g.index = pd.Index(g["locus_tag"].astype(str))
    rows: list[dict[str, Any]] = []
    for (rid, sub), grp in rg.groupby([rg["region_id"].astype(str), rg["subset"].astype(str)]):
        tags = [t for t in dict.fromkeys(grp["locus_tag"].astype(str)) if t in g.index]
        mem = g.loc[tags]
        el = mem["_el"].to_numpy(bool)
        m = int(el.sum())
        o, e = float(mem["_d"].to_numpy(bool)[el].sum()), float(mem["_p"].to_numpy(float)[el].sum())
        det = _detection(o, e, _NA) if m else _detection(_NA, _NA, _NA)
        flags = ("small_region",) * (m < SMALL_REGION) + (DESCRIPTIVE,)
        rows.append({"region_id": rid, "region_type": "prophage", "subset": sub,
                     "contig": contig_of[rid], "start": None, "end": None, "n_genes": len(tags),
                     "n_eligible": m, **det, "stat": "o_minus_e", "effect": det["o_minus_e"],
                     "ci_low": _NA, "ci_high": _NA, "p_value": _NA, "q_value": _NA,
                     "method": DESCRIPTIVE, "null_model": "none", "n_perm": 0,
                     "label": DESCRIPTIVE, "flags": flags})  # fmt: skip
    return _frame(rows)


def _test_rows(tests: pd.DataFrame) -> list[dict[str, Any]]:
    """M4-17 rows with ``subset = all`` and the detection columns of both stat rows."""
    key = list(zip(tests["region_id"].astype(str), tests["contig"].astype(str), strict=True))
    records = tests.to_dict("records")
    oe = {k: r for k, r in zip(key, records, strict=True) if r["stat"] == "o_minus_e"}
    hl = {k: float(r["effect"]) for k, r in zip(key, records, strict=True)
          if r["stat"] == "hl_shift"}  # fmt: skip
    rows = []
    for k, r in zip(key, records, strict=True):
        src = oe.get(k)
        o, e = (_NA, _NA) if src is None else (float(src["observed"]), float(src["expected"]))
        det = _detection(o, e, hl.get(k, _NA))
        det["o_minus_e"] = _NA if src is None else float(src["effect"])
        rows.append({**{c: r.get(c) for c in REGIONS_STATS_COLUMNS}, **det, "subset": "all",
                     "start": None, "end": None, "label": None,
                     "flags": tuple(_ids(r["flags"]))})  # fmt: skip
    return rows


def _check_rules(df: pd.DataFrame) -> None:
    """Lead decision 9: ``null_model = none`` => p, q NA and ``n_perm = 0``; a p-value =>
    ``null_model != none``; a q-value => a p-value."""
    none = df["null_model"].astype(str) == "none"
    p, q = df["p_value"].notna(), df["q_value"].notna()
    bad = (none & (p | q | (df["n_perm"] != 0))) | (q & ~p)
    if bad.any():
        bad_rows = df.loc[bad, ["region_id", "stat"]].to_dict("records")[:5]
        raise ValueError(f"regions_stats rows break the null_model rule: {bad_rows}")


def regions_stats(
    tests: pd.DataFrame, plasmids: pd.DataFrame, subsets: pd.DataFrame, regions: pd.DataFrame
) -> pd.DataFrame:
    """The ``regions_stats`` table (``REGIONS_STATS_COLUMNS``): M4-17 ``tests`` (q-values
    copied), ``plasmids`` and ``subsets`` rows with C7 coordinates (0-based half-open; AMR NA).
    Rows sort by region type (default ``region_types`` order), C7 contig (order of first
    appearance in ``regions``), C7 start, region_id; within a region ``subset = all`` first (then
    §7.1.7 order), ``stat`` enum order, test contig, unadjusted before adjusted."""
    rows = _test_rows(tests) if len(tests) else []
    rows += plasmids.to_dict("records") + subsets.to_dict("records")
    df = _frame([{c: r.get(c) for c in REGIONS_STATS_COLUMNS} for r in rows])
    if df.empty:
        return df
    reg = regions.drop_duplicates("region_id")
    reg.index = pd.Index(reg["region_id"].astype(str))
    rid = df["region_id"].astype(str)
    df["start"] = pd.array(rid.map(reg["start"]).to_numpy(object), dtype="Int64")
    df["end"] = pd.array(rid.map(reg["end"]).to_numpy(object), dtype="Int64")
    crank = {c: i for i, c in enumerate(dict.fromkeys(regions["contig"].astype(str)))}
    contig = df["contig"].astype(str)
    c7_contig = rid.map(reg["contig"].astype(str)).fillna(contig)
    sub_rank = {"all": 0} | {s: i + 1 for i, s in enumerate(SUBSETS)}
    keys = pd.DataFrame(
        {
            "type": df["region_type"].map({t: i for i, t in enumerate(TYPE_ORDER)}),
            "c7_contig": c7_contig.map(crank).fillna(len(crank)),
            "start": df["start"].astype("float64").fillna(-1.0),
            "region_id": rid,
            "subset": df["subset"].map(sub_rank),
            "stat": df["stat"].map({s: i for i, s in enumerate(STATS)}),
            "contig_rank": contig.map(crank).fillna(len(crank)),
            "contig": contig,
            "adjusted": (df["label"].astype(str) == ADJUSTED_LABEL).astype(int),
        }
    )
    order = keys.sort_values(list(keys.columns), kind="mergesort").index
    out = df.loc[order].reset_index(drop=True)
    _check_rules(out)
    return out
