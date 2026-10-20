"""Inputs and oracles for the M4-18 ``regions_stats`` tests.

docs/specs/M4-18-regions-stats-table.md.

Gene tables are the hand-built §9.3 ``genes.tsv`` frames of ``region_tests_data.gene_table`` (M4-17)
plus the M4-14 column ``abundance_imputed`` and the M4-16 summary columns ``hypothetical`` /
``cog_category`` (M4-15 is not merged; the column names follow §9.3 and the M4-06 contract). C7 /
C8 frames follow the contract columns; plasmid regions span the whole replicon (M2-19: "genes = all
genes on it"), prophage members carry a C8 ``subset`` from §7.1.7 (integrase, repressor,
structural, lysis, cargo). Region IDs are ``<type>_<n>`` (§6.5), so a region's type is readable
from the ``region_ids`` of a gene.

Oracles implement the lead decisions of the spec (2026-10-20) independently of ``giap.context``:

* ``reference_contigs`` (decision 3): ``kind = chromosome`` replicons; if none, every replicon that
  is neither ``kind = plasmid`` nor a contig of a plasmid region (flag ``no_reference_replicon``).
* ``plasmid_hl`` (decision 2): Hodges-Lehmann shift of the plasmid's detected, non-imputed genes
  with a finite raw ``abundance_log2`` against the same genes of the reference replicons (no
  residualisation, no imputed genes); NaN when a side is empty.
* ``detection_counts``: observed = detected eligible genes, expected = sum of ``p_detect`` over the
  eligible genes (not partial, ``p_detect`` present), O - E and O / E (NA when E = 0, decision 5).
* ``weighted_depth`` (spec "Plasmid rows"): length-weighted mean depth of a replicon set.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import region_tests_data as rd

#: §9.3 ``regions_stats.tsv`` order (region, detection, abundance), the §9.2 test columns in their
#: §9.2 order (``stat`` first), then ``label`` and ``flags`` (test-author reading of
#: "§9.3 order + the §9.2 test columns")
REGIONS_STATS_COLUMNS = (
    "region_id", "region_type", "subset", "contig", "start", "end", "n_genes", "n_eligible",
    "observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2",
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model", "n_perm",
    "label", "flags",
)  # fmt: skip
TEST_COLUMNS = ("stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method",
                "null_model", "n_perm")  # fmt: skip
DETECTION_COLUMNS = ("observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2")
#: §7.1.7 prophage function subsets (M4-06 C8 enum)
SUBSETS = ("integrase", "repressor", "structural", "lysis", "cargo")
#: default ``context.region_types`` order (sort key of the table)
TYPE_ORDER = ("prophage", "plasmid", "bgc", "crispr", "amr")
#: §9.4 "Region and context" + M4-17 / M4-index additions that rows of this table can carry
REGION_FLAGS = (
    "small_region", "dosage_plus_regulation", "descriptive", "no_oric", "uncalibrated",
    "few_shifts", "split_contigs", "hl_variance_nonpositive", "no_reference_replicon",
    "depth_adjusted",
)  # fmt: skip
#: lead decision 13: flag of the dosage-adjusted plasmid row; method / label of subset rows
DEPTH_FLAG = "depth_adjusted"
DESCRIPTIVE = "descriptive"
PLASMID_METHOD = "replicon_contrast"
PLASMID_LABEL = "dosage_plus_regulation"
ADJUSTED_LABEL = "dosage_adjusted"


@dataclass
class Genome:
    genes: pd.DataFrame
    replicons: pd.DataFrame
    regions: pd.DataFrame
    region_genes: pd.DataFrame


def genome(
    contigs: Sequence[rd.Contig],
    members: Mapping[str, tuple[str, Sequence[str]]] | None = None,
    *,
    plasmids: Iterable[str] = (),
    subsets: Mapping[str, str] | None = None,
    imputed: Iterable[str] = (),
) -> Genome:
    """Genes, C1, C7, C8 for ``contigs``.

    ``members``: ``{region_id: (region_type, locus_tags)}`` (M4-17 helper semantics);
    ``plasmids``: contig names that get a ``plasmid_<n>`` region over all their genes (numbered in
    the order given); ``subsets``: C8 ``subset`` per prophage member locus tag; ``imputed``: locus
    tags with ``abundance_imputed = true`` (they must be detected; their ``abundance_log2`` is
    whatever the contig gives them, i.e. the fill value).
    """
    genes, replicons = rd.gene_table(contigs)
    mem: dict[str, tuple[str, list[str]]] = {k: (t, list(v)) for k, (t, v) in
                                             (members or {}).items()}  # fmt: skip
    for i, c in enumerate(plasmids, start=1):
        tags = genes.loc[genes["contig"] == c, "locus_tag"].astype(str).tolist()
        mem[f"plasmid_{i}"] = ("plasmid", tags)
    genes, regions, rg = rd.region_frames(genes, mem)
    length = dict(zip(replicons["contig"], replicons["length"], strict=True))
    for j in regions.index[regions["region_type"] == "plasmid"]:
        regions.loc[j, "start"] = 0
        regions.loc[j, "end"] = int(length[str(regions.loc[j, "contig"])])
    sub = dict(subsets or {})
    ptype = dict(zip(regions["region_id"], regions["region_type"], strict=True))
    rg["subset"] = [
        sub.get(str(t)) if ptype[str(r)] == "prophage" else None
        for r, t in zip(rg["region_id"], rg["locus_tag"], strict=True)
    ]
    imp = set(map(str, imputed))
    genes["abundance_imputed"] = genes["locus_tag"].astype(str).isin(imp).to_numpy(bool)
    bad = genes["abundance_imputed"] & ~genes["detected"].astype(bool)
    if bad.any():
        raise ValueError(
            "imputed genes must be detected (M4-14: undetected genes are never imputed)"
        )
    genes["hypothetical"] = False
    genes["cog_category"] = None
    return Genome(genes, replicons, regions, rg)


# --------------------------------------------------------------------------- oracles


def eligible(genes: pd.DataFrame) -> np.ndarray:
    p = pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float)
    return ~genes["partial"].astype(bool).to_numpy() & np.isfinite(p)


def detection_counts(genes: pd.DataFrame, contigs: Iterable[str]) -> dict[str, float]:
    on = genes["contig"].astype(str).isin(set(contigs)).to_numpy()
    el = eligible(genes) & on
    d = genes["detected"].astype(bool).to_numpy() & el
    p = pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float)
    o = float(d.sum())
    e = float(p[el].sum())
    return {
        "n_genes": float(on.sum()),
        "n_eligible": float(el.sum()),
        "observed": o,
        "expected": e,
        "o_minus_e": o - e,
        "o_over_e": o / e if e > 0 else math.nan,
    }


def reference_contigs(
    replicons: pd.DataFrame, plasmid_contigs: Iterable[str]
) -> tuple[list[str], bool]:
    """(reference replicons, fallback used) per lead decision 3."""
    kind = dict(zip(replicons["contig"].astype(str), replicons["kind"].astype(str), strict=True))
    chrom = [c for c, k in kind.items() if k == "chromosome"]
    if chrom:
        return chrom, False
    pl = set(plasmid_contigs)
    return [c for c, k in kind.items() if k != "plasmid" and c not in pl], True


def hl_values(genes: pd.DataFrame, contigs: Iterable[str]) -> np.ndarray:
    """Raw log2 abundance of the detected, eligible, non-imputed genes of ``contigs`` (finite)."""
    on = genes["contig"].astype(str).isin(set(contigs)).to_numpy()
    a = pd.to_numeric(genes["abundance_log2"], errors="coerce").to_numpy(float)
    ok = (on & eligible(genes) & genes["detected"].astype(bool).to_numpy()
          & ~genes["abundance_imputed"].astype(bool).to_numpy() & np.isfinite(a))  # fmt: skip
    return a[ok]


def plasmid_hl(genes: pd.DataFrame, plasmid: str, reference: Iterable[str]) -> float:
    return rd.naive_hl(list(hl_values(genes, [plasmid])), list(hl_values(genes, reference)))


def weighted_depth(
    replicons: pd.DataFrame, depth: Mapping[str, float], contigs: Iterable[str]
) -> float:
    length = dict(zip(replicons["contig"].astype(str), replicons["length"].astype(float),
                      strict=True))  # fmt: skip
    cs = list(contigs)
    return sum(length[c] * depth[c] for c in cs) / sum(length[c] for c in cs)


def depth_frame(depth: Mapping[str, float]) -> pd.DataFrame:
    return pd.DataFrame({"contig": list(depth), "mean_depth": [float(v) for v in depth.values()]})


# --------------------------------------------------------------------------- row access


def rows(df: pd.DataFrame, region_id: str, stat: str | None = None, *, subset: str | None = None,
         label: str | None = None) -> pd.DataFrame:  # fmt: skip
    m = df["region_id"].astype(str) == region_id
    if stat is not None:
        m &= df["stat"].astype(str) == stat
    if subset is not None and "subset" in df:
        m &= df["subset"].astype(str) == subset
    if label is not None:
        m &= df["label"].astype(str) == label
    return df[m]


def one(df: pd.DataFrame, region_id: str, stat: str, **kw: Any) -> pd.Series:
    r = rows(df, region_id, stat, **kw)
    if len(r) != 1:
        raise AssertionError(f"{len(r)} rows for {region_id}/{stat}/{kw}: {r.to_dict('records')}")
    return r.iloc[0]


def unadjusted(df: pd.DataFrame, region_id: str, stat: str) -> pd.Series:
    """The plasmid row of ``stat`` that is not the depth-adjusted one."""
    r = rows(df, region_id, stat)
    r = r[r["label"].astype(str) != ADJUSTED_LABEL]
    if len(r) != 1:
        raise AssertionError(f"{len(r)} unadjusted rows for {region_id}/{stat}")
    return r.iloc[0]


def na(x: Any) -> bool:
    return rd.is_na(x) or (isinstance(x, float | np.floating) and math.isnan(float(x)))


def flags(cell: Any) -> set[str]:
    return rd.flag_set(cell)


# --------------------------------------------------------------------------- fixtures and runners

#: bootstrap B of the unit fixtures
B = 200
Contig = rd.Contig


def module() -> Any:
    """``giap.context.regions_table`` (imported lazily: the tests must fail, not error at
    collection, while it is absent)."""
    from giap.context import regions_table as T

    return T


def rng(seed: int = 1) -> np.random.Generator:
    from giap.rng import for_step

    return for_step("context_stats", seed)


def plasmid_rows(g: Genome, depth: pd.DataFrame | None = None, cfg: Any = None,
                 seed: int = 1) -> pd.DataFrame:  # fmt: skip
    cfg = cfg if cfg is not None else rd.config(bootstrap=B)
    return module().plasmid_rows(g.genes, g.replicons, depth, cfg, rng(seed))  # type: ignore[no-any-return]


def tests(g: Genome, seed: int = 1) -> pd.DataFrame:
    """M4-17 ``test_all_regions`` (merged) on the fixture."""
    from giap.context import regions as R

    return R.test_all_regions(g.genes, g.regions, g.region_genes, g.replicons,
                              rd.config(bootstrap=50), rng(seed))  # fmt: skip


def table(g: Genome, depth: pd.DataFrame | None = None, seed: int = 1) -> pd.DataFrame:
    """``regions_stats`` over M4-17 tests, plasmid rows and subset rows of the fixture."""
    T = module()
    pl = plasmid_rows(g, depth, seed=seed)
    sub = T.subset_rows(g.genes, g.regions, g.region_genes)
    return T.regions_stats(tests(g, seed), pl, sub, g.regions)  # type: ignore[no-any-return]


def abund(n: int, base: float, step: float = 0.5) -> list[float]:
    """Distinct log2 abundances base + step * (i mod 7) + 0.01 i."""
    return [base + step * (i % 7) + 0.01 * i for i in range(n)]


def pairs(name: str, n: int) -> list[str]:
    """Two-gene TUs."""
    return [f"TU_{name}_{i // 2}" for i in range(n)]


def chr_contig(name: str = "chr1", n: int = 120, base: float = 10.0, **kw: Any) -> rd.Contig:
    d = [0 if i % 3 == 0 else 1 for i in range(n)]
    kw.setdefault("tu", pairs(name, n))
    return rd.Contig(name, d, [0.6] * n, a=abund(n, base), **kw)


def plasmid_contig(name: str = "p1", n: int = 60, base: float = 14.3, **kw: Any) -> rd.Contig:
    d = [0 if i % 4 == 0 else 1 for i in range(n)]
    kw.setdefault("tu", pairs(name, n))
    kw.setdefault("kind", "plasmid")
    return rd.Contig(name, d, [0.7] * n, a=abund(n, base), **kw)


def example_1(**kw: Any) -> Genome:
    """Spec example 1: chromosome (120 genes) + plasmid p1 (60 genes) 4 log2 units above it (the
    0.01 i trend of the longer chromosome is offset by the plasmid base 14.3; HL = 4.0)."""
    return genome([chr_contig(), plasmid_contig()], plasmids=["p1"], **kw)


def prophage_genome() -> Genome:
    """Spec example 3: prophage_1 = 16 structural genes (8 detected) + 4 lysis genes (none
    detected) on a 120-gene chromosome, plus bgc_1."""
    n = 120
    d = [0] * n
    for i in [*range(20, 28), *range(50, 60)]:
        d[i] = 1
    p = [0.25 + 0.005 * i for i in range(n)]
    c = rd.Contig("chr1", d, p, tu=pairs("chr1", n))
    tags = genome([c]).genes["locus_tag"].tolist()
    structural, lysis = tags[20:36], tags[36:40]
    sub = {t: "structural" for t in structural} | {t: "lysis" for t in lysis}
    mem = {"prophage_1": ("prophage", structural + lysis), "bgc_1": ("bgc", tags[50:60])}
    return genome([c], mem, subsets=sub)


def mixed_genome() -> Genome:
    """Every region type: prophage_1 (with structural / cargo subsets), prophage_2 (lower start
    than prophage_1), bgc_1, crispr_1, scattered amr_1 on chr1 (200 genes) and plasmid_1 on p1."""
    n = 200
    d = [1 if (i * 7) % 5 < 3 else 0 for i in range(n)]
    p = [0.3 + 0.4 * ((i * 3) % 7) / 6 for i in range(n)]
    c = rd.Contig("chr1", d, p, tu=pairs("chr1", n))
    contigs = [c, plasmid_contig()]
    g0 = genome(contigs)
    tags = g0.genes.loc[g0.genes["contig"] == "chr1", "locus_tag"].tolist()
    mem = {
        "prophage_1": ("prophage", tags[100:130]),
        "prophage_2": ("prophage", tags[20:40]),
        "bgc_1": ("bgc", tags[60:75]),
        "crispr_1": ("crispr", tags[150:154]),
        "amr_1": ("amr", [tags[5], tags[90], tags[170]]),
    }
    sub = {t: ("structural" if k % 2 else "cargo") for k, t in enumerate(tags[100:130])}
    return genome(contigs, mem, plasmids=["p1"], subsets=sub)


SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"


def schema() -> dict[str, Any]:
    """The ``regions_stats`` table schema (``regions_stats.schema.json`` or ``.json``)."""
    found = [p for p in (SCHEMA_DIR / "regions_stats.json",
                         SCHEMA_DIR / "regions_stats.schema.json") if p.exists()]  # fmt: skip
    if len(found) != 1:
        raise AssertionError(f"expected one regions_stats schema, found {found}")
    return json.loads(found[0].read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def stat_order() -> list[str]:
    """Lead decision 11: rows of a region follow the schema's ``stat`` enum order."""
    f = next(f for f in schema()["fields"] if f["name"] == "stat")
    return [str(x) for x in f["constraints"]["enum"]]


def check_row_order(df: pd.DataFrame) -> None:
    """Lead decision 11 on a regions_stats frame: region types in the default ``region_types``
    order, a region's rows contiguous, ``subset = all`` (M4-17 / plasmid rows) before subset rows,
    and within each (region, subset) the ``stat`` enum order of the schema."""
    types = list(dict.fromkeys(df["region_type"].astype(str)))
    assert types == [t for t in TYPE_ORDER if t in types], types
    rid = df["region_id"].astype(str).to_numpy()
    seen: list[str] = []
    for r in rid:
        if not seen or seen[-1] != r:
            assert r not in seen, f"rows of {r} are not contiguous"
            seen.append(r)
    order = stat_order()
    for r in seen:
        grp = df[df["region_id"].astype(str) == r]
        subs = list(grp["subset"].astype(str))
        n_all = subs.count("all")
        assert subs[:n_all] == ["all"] * n_all, (r, subs)
        for _, g2 in grp.groupby(grp["subset"].astype(str), sort=False):
            rank = [order.index(str(s)) for s in g2["stat"]]
            assert rank == sorted(rank), (r, list(g2["stat"]))
