"""Inputs for the M4-19 sensitivity tests (docs/specs/M4-19-sensitivity-regions.md, lead decisions
1-4 of 2026-10-20; GIAP_SPEC.md §8.6, §2.4 G9, §7.3.2, §5.2).

* ``family_genes``: a hand-built §9.3 gene table (``region_tests_data.gene_table``) of a 200-gene
  circular chromosome ``chr1`` and a 100-gene circular plasmid ``p1`` with four regions:
  ``prophage_1`` (chr1 genes 0-9, every gene detected at p_detect 0.2: a strong O - E excess),
  ``bgc_1`` (chr1 genes 100-109, null-like), ``amr_1`` (split over chr1 150-154 and p1 0-4, both
  parts with an O - E excess; M4-17 tests it once per contig, flag ``split_contigs``) and
  ``crispr_1`` (chr1 60-61, 2 genes). The rest alternates detected / not detected at p 0.5.
  Region IDs follow §6.5 ``<type>_<n>``, so a region's type is its ID prefix; ``region_ids`` is
  filled (sorted) as C6 / genes.tsv carry it.
* ``fake_tests``: an M4-17 ``test_all_regions`` result frame (``TABLE_COLUMNS``) with chosen p / q
  / null_model, for the flagged-set rule (lead decision 4) without running the tests.
* ``Recorder``: a sensitivity family that records each call (genes, cfg, the generator's state
  at entry) and returns a fixed set, for the common-random-numbers checks (lead decision 2).
* ``sim_work_dir``: a giapsim work directory (M4-07/08/09 contract tables, flat C1-C10, M3 TSVs
  under ``mapping/`` and copied flat) with the M4-15 gene table at ``genetable/genes.parquet``
  (``run_gene_table`` then ``run_detectability``, PR #111).

``call_regions_family`` is the one place that encodes how ``regions_family`` receives its inputs:
``(genes, cfg, rng, *, regions, region_genes, replicons)`` (lead decisions 5 and 13: the step
binds C7 / C8 / C1 with ``functools.partial``; no fallback to ``region_ids``).
"""

from __future__ import annotations

import copy
import math
import shutil
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import region_tests_data as rd
from region_tests_data import Contig, gene_table, region_frames, tags_of

#: §8.6 sensitivity.tsv columns, in order (``agreement`` follows for the determinants family)
COLUMNS = ("family", "variant", "n_primary", "n_variant", "n_shared", "jaccard")
TABLE_COLUMNS = rd.TABLE_COLUMNS
SENS_B = 199  # §8.6 "Sensitivity reruns use B = 199"


def family_id(region_id: str, contig: str, stat: str) -> str:
    """Lead decision 3: ``<region_id>:<contig>:<stat>``."""
    return f"{region_id}:{contig}:{stat}"


def config(**over: Any) -> Any:
    """``GiapConfig`` with a small bootstrap (CIs are not compared) and nested overrides, e.g.
    ``config(stats={"permutations": 199}, output={"sensitivity": "full"})``."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {"seed": 20261020, "stats": {"bootstrap": 30}}
    for k, v in over.items():
        if isinstance(v, Mapping) and isinstance(raw.get(k), dict):
            raw[k] = {**raw[k], **v}
        else:
            raw[k] = v
    return GiapConfig.model_validate(raw)


# --------------------------------------------------------------------------- hand-built genes


def _chr1() -> Contig:
    n = 200
    d = [i % 2 for i in range(n)]
    p = [0.5] * n
    a = [10.0 + 0.1 * (i % 7) for i in range(n)]
    for i in range(10):  # prophage_1: detected, p 0.2
        d[i], p[i] = 1, 0.2
    for i in range(150, 155):  # amr_1 (chr1 part): detected, p 0.2
        d[i], p[i] = 1, 0.2
    return Contig("chr1", d, p, a)


def _p1() -> Contig:
    n = 100
    d = [i % 2 for i in range(n)]
    p = [0.5] * n
    a = [9.0 + 0.1 * (i % 5) for i in range(n)]
    for i in range(5):  # amr_1 (p1 part)
        d[i], p[i] = 1, 0.2
    return Contig("p1", d, p, a)


@dataclass
class FamilyTables:
    genes: pd.DataFrame
    regions: pd.DataFrame
    region_genes: pd.DataFrame
    replicons: pd.DataFrame


def family_tables() -> FamilyTables:
    genes, replicons = gene_table([_chr1(), _p1()])
    members = {
        "amr_1": ("amr", tags_of(genes, "chr1", range(150, 155)) + tags_of(genes, "p1", range(5))),
        "bgc_1": ("bgc", tags_of(genes, "chr1", range(100, 110))),
        "crispr_1": ("crispr", tags_of(genes, "chr1", [60, 61])),
        "prophage_1": ("prophage", tags_of(genes, "chr1", range(10))),
    }
    genes, regions, rg = region_frames(genes, members)
    genes["abundance_imputed"] = False
    return FamilyTables(genes, regions, rg, replicons)


def bound_tables(p1_topology: str = "linear") -> FamilyTables:
    """Lead decisions 5 / 13 (M27): C7 / C8 / C1 that differ from what ``region_ids`` implies.

    ``family_tables`` plus ``bgc_2`` = p1 genes 10-94 (85 genes, all detected at p_detect 0.2), a
    region present only in C7 / C8; every gene's ``region_ids`` is emptied; p1 is
    ``p1_topology``. On a circular p1 ``bgc_2`` O - E has p = 1/100 (flagged); on a linear p1 only
    shifts k = 1..5 keep the set inside the contig, so p >= 1/6 (never flagged)."""
    t = family_tables()
    g = t.genes.copy()
    on_p1 = (g["contig"] == "p1").to_numpy()
    pos = np.cumsum(on_p1) - 1
    hot = on_p1 & (pos >= 10) & (pos <= 94)
    g.loc[hot, "detected"] = True
    g.loc[hot, "detection_status"] = "detected"
    g.loc[hot, "p_detect"] = 0.2
    g.loc[hot & g["abundance_log2"].isna(), "abundance_log2"] = 9.5
    g["region_ids"] = [[] for _ in range(len(g))]
    tags = g.loc[hot, "locus_tag"].astype(str).tolist()
    reg = t.regions.copy()
    s, e = int(g.loc[hot, "start"].min()), int(g.loc[hot, "end"].max())
    extra = pd.DataFrame([("bgc_2", "bgc", "p1", s, e, "test", None)], columns=list(reg.columns))
    reg = pd.concat([reg, extra.astype(reg.dtypes.to_dict())], ignore_index=True)
    rg_extra = pd.DataFrame([("bgc_2", t_, 1.0, None, False) for t_ in tags],
                            columns=list(t.region_genes.columns))  # fmt: skip
    rg = pd.concat([t.region_genes, rg_extra.astype(t.region_genes.dtypes.to_dict())],
                   ignore_index=True)  # fmt: skip
    rep = t.replicons.copy()
    rep.loc[rep["contig"] == "p1", "topology"] = p1_topology
    return FamilyTables(g, reg, rg, rep)


def with_missing_abundance(genes: pd.DataFrame, tags: Iterable[str]) -> pd.DataFrame:
    """``genes`` with ``abundance_log2`` NaN on ``tags`` (what ``missing=exclude`` does to genes
    that the primary imputed); detection columns unchanged."""
    out = genes.copy()
    m = out["locus_tag"].isin(set(tags))
    out.loc[m, "abundance_log2"] = math.nan
    out.loc[m, "abundance_imputed"] = False
    return out


def call_regions_family(t: FamilyTables, cfg: Any, rng: np.random.Generator) -> frozenset[str]:
    """``giap.context.sensitivity.regions_family`` on the tables ``t`` (spec signature)."""
    from giap.context import sensitivity as S

    return S.regions_family(
        t.genes, cfg, rng, regions=t.regions, region_genes=t.region_genes, replicons=t.replicons
    )


def reported_flagged(t: FamilyTables, cfg: Any, rng: np.random.Generator) -> frozenset[str]:
    """Lead decision 4 applied to the reported M4-17 tests (what regions_stats.tsv holds):
    ``null_model = circular_shift``, p not NA, q <= ``alpha.context``."""
    from giap.context import regions as R

    res = R.test_all_regions(t.genes, t.regions, t.region_genes, t.replicons, cfg, rng)
    return flagged_from(res, cfg.stats.alpha.context)


def flagged_from(res: pd.DataFrame, alpha: float) -> frozenset[str]:
    keep = (
        (res["null_model"].astype(str) == "circular_shift")
        & res["p_value"].notna()
        & res["q_value"].notna()
        & (res["q_value"].astype(float) <= alpha)
    )
    return frozenset(
        family_id(str(r), str(c), str(s))
        for r, c, s in zip(
            res.loc[keep, "region_id"], res.loc[keep, "contig"], res.loc[keep, "stat"], strict=True
        )
    )


# --------------------------------------------------------------------------- fake M4-17 output


def fake_tests(rows: Iterable[tuple[str, str, str, str, float, float, str]]) -> pd.DataFrame:
    """M4-17 result frame from (region_id, region_type, contig, stat, p, q, null_model) rows."""
    out = []
    for rid, typ, contig, stat, p, q, null in rows:
        out.append(
            {
                "region_id": rid, "region_type": typ, "contig": contig, "stat": stat,
                "observed": 1.0, "expected": 0.5 if stat == "o_minus_e" else math.nan,
                "effect": 0.5, "ci_low": 0.1, "ci_high": 0.9, "p_value": p, "n_perm": 199,
                "n_genes": 5, "n_eligible": 5, "flags": (), "q_value": q,
                "method": rd.METHOD[stat], "null_model": null,
            }
        )  # fmt: skip
    return pd.DataFrame(out, columns=list(TABLE_COLUMNS)).astype(
        {"p_value": "float64", "q_value": "float64", "n_perm": "int64"}
    )


def patch_region_tests(monkeypatch: Any, frame: pd.DataFrame) -> list[dict[str, Any]]:
    """Replace M4-17 ``test_all_regions`` (module attribute and any name imported into
    ``giap.context.sensitivity``) by a stub returning ``frame``; returns the call log."""
    from giap.context import regions as R
    from giap.context import sensitivity as S

    calls: list[dict[str, Any]] = []

    def stub(*args: Any, **kwargs: Any) -> pd.DataFrame:
        calls.append({"args": args, "kwargs": kwargs})
        return frame.copy()

    monkeypatch.setattr(R, "test_all_regions", stub)
    monkeypatch.setattr(S, "test_all_regions", stub, raising=False)
    return calls


# --------------------------------------------------------------------------- recorder family


@dataclass
class Call:
    genes: pd.DataFrame
    cfg: Any
    state: dict[str, Any]
    first_draw: int


@dataclass
class Recorder:
    """Sensitivity family recording every call; returns ``result`` (fixed)."""

    result: frozenset[str] = frozenset()
    calls: list[Call] = field(default_factory=list)

    def __call__(self, genes: pd.DataFrame, cfg: Any, rng: np.random.Generator) -> frozenset[str]:
        state = copy.deepcopy(rng.bit_generator.state)
        probe = np.random.Generator(type(rng.bit_generator)())
        probe.bit_generator.state = copy.deepcopy(state)
        self.calls.append(Call(genes.copy(), cfg, state, int(probe.integers(0, 2**62))))
        return self.result


# --------------------------------------------------------------------------- simulated work dir

SIM_SEED = 20261020
SIM_GENES = 1000


def sim_work_dir(
    work: Path,
    *,
    evidence: Mapping[str, Any] | None = None,
    blank_intensity_genes: int = 60,
    cfg_over: Mapping[str, Any] | None = None,
    seed: int = SIM_SEED,
    layout: Mapping[str, Any] | None = None,
) -> Path:
    """giapsim work directory with the M4-15 gene table; returns ``work``.

    Expression: ``region_detect_shift = {prophage: -3}`` (prophage_1 O - E flagged at the default
    seed) and ``region_effect_log2 = {bgc: 3}``. ``blank_intensity_genes`` genes (drawn among the
    genes with peptide hits) lose every peptide intensity, so detected ones among them have a
    missing abundance (imputed under the primary ``median``; NA under ``exclude``).
    ``evidence``: ``giapsim.evidence.EvidenceConfig`` fields (e.g. ``emit = {"intensity"}`` gives
    the ``two_peptide`` tier); ``layout``: ``giapsim.layout.LayoutConfig`` fields (e.g.
    ``linear_contigs = (80,)`` adds a linear contig with partial end genes)."""
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la
    from giap.steps.detectability import run_detectability  # M4-15 (PR #111)
    from giap.steps.gene_table import run_gene_table  # M4-15 (PR #111)
    from giapsim.contract_tables import write_contract_tables

    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    lcfg = la.LayoutConfig(n_genes_chromosome=SIM_GENES, **dict(layout or {}))
    lay = la.simulate_layout(lcfg, rng)
    ecfg = ex.ExpressionConfig(
        region_detect_shift={"prophage": -3.0}, region_effect_log2={"bgc": 3.0}
    )
    expr = ex.simulate_expression(lay, ecfg, rng)
    e = ev.simulate_evidence(lay, expr, ev.EvidenceConfig(**dict(evidence or {})), rng)
    write_contract_tables(lay, work)
    ev.write_evidence_tables(e, work)
    if blank_intensity_genes:
        ph = pd.read_csv(work / "mapping" / "peptide_hits.tsv", sep="\t")
        pool = sorted(ph["gene"].dropna().astype(str).unique())
        pick = set(np.random.default_rng(3).choice(pool, blank_intensity_genes, replace=False))
        peps = set(ph.loc[ph["gene"].astype(str).isin(pick), "peptide"])
        pep_path = work / "peptides.norm.parquet"
        pep = pd.read_parquet(pep_path)
        if "intensity" in pep.columns:
            pep.loc[pep["peptide"].isin(peps), "intensity"] = np.nan
            pep.to_parquet(pep_path, index=False)
    for tsv in (work / "mapping").glob("*.tsv"):  # both layouts (M4-15 lead decision 8)
        shutil.copy(tsv, work / tsv.name)
    cfg = config(seed=seed, **dict(cfg_over or {}))
    run_gene_table(work, cfg)
    run_detectability(work, cfg)
    return work


def sim_genes(work: Path) -> pd.DataFrame:
    return pd.read_parquet(work / "genetable" / "genes.parquet")


def sim_m4(work: Path) -> Any:
    from giap.io.work_tables import default_paths, load_m4_inputs

    paths = {k: p for k, p in default_paths(work).items() if p.exists()}
    m4, _ = load_m4_inputs(paths)
    return m4


def sensitivity_tsv(work: Path) -> Path:
    return work / "tables" / "sensitivity.tsv"


def read_sensitivity(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
