"""Inputs and oracles for the M5-10 ``islands`` sensitivity family
(docs/specs/M5-10-sensitivity-islands.md, lead decision 38; GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.5(d)).

Interface assumed (spec "Public interface", decision 38 = M4-19 decision 5):
``giap.context.sensitivity.islands_family(genes, cfg, rng, *, tus, replicons) -> frozenset[str]``
registered as ``FAMILIES["islands"]``; the step binds C5 ``tus`` and C1 ``replicons`` by keyword
(``functools.partial``). The flagged set is ``{"<direction>:<locus_tag>"}`` over every gene of
every island row of M5-06 ``islands_frame`` (best-window gene lists, decision Q3 of M5-06), so a
gene in a high and a low island contributes both IDs.

* ``unit_tables``: a hand-built two-replicon gene table (``islands_data.sim_contig``): a
  400-gene circular ``chr1`` with a planted +3 SD high island of 12 TUs (found at p = 0.005 with
  ``unit_config``) and a 120-gene linear ``p1`` under the null.
* ``oracle_flagged``: the M5-06 pipeline as ``context_stats`` runs it (``tu_series`` per C1
  replicon, ``scan_replicons`` then ``islands_frame`` on one generator), turned into IDs.
* ``fake_islands`` / ``patch_islands_frame``: replace ``islands_frame`` by a stub returning chosen
  island rows, to pin the ID rule without depending on what a scan finds.
* ``sim_work_dir``: a giapsim work directory (M5-04 scenario ``island_power``: 1,000-gene
  chromosome, AR(1) rho = 0.5, one 15-TU island at +4 SD with abundance-dependent detection),
  peptide intensities blanked on ``blank_intensity_genes`` genes (imputed under the primary
  ``median``, NA under ``exclude``), then the M4-15 gene table. With ``work_config`` (LOESS span
  0.30, block length 1) the planted island is found at B = 199.
"""

from __future__ import annotations

import shutil
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import islands_data as isd

#: §8.6 "Sensitivity reruns use B = 199" (M5-10: both ``stats.bootstrap`` and
#: ``stats.permutations``, for the primary and for every variant)
SENS_B = 199
FAMILY = "islands"
#: M4-19 decision 2 (common random numbers): one stream per family
STREAM = "sensitivity/islands"
ISLAND_COLUMNS = isd.ISLAND_COLUMNS


def island_id(direction: str, locus_tag: str) -> str:
    """Spec "Public interface": ``<direction>:<locus_tag>``."""
    return f"{direction}:{locus_tag}"


# --------------------------------------------------------------------------- hand-built tables


@dataclass
class IslandTables:
    genes: pd.DataFrame
    tus: pd.DataFrame
    replicons: pd.DataFrame


def unit_tables(*, effect: float = 3.0, seed: int = 5) -> IslandTables:
    """``chr1`` (400 genes, circular, island at TUs 100-111 shifted by ``effect`` SD, 80 %
    detected) and ``p1`` (120 genes, linear, null)."""
    rng = np.random.default_rng(seed)
    c = isd.sim_contig(rng, "chr1", 400, rho=0.3, island=(100, 12, effect), detected_frac=0.8)
    p = isd.sim_contig(rng, "p1", 120, rho=0.3, circular=False)
    genes, tus, replicons = isd.frames([c, p])
    return IslandTables(genes, tus, replicons)


def tiny_tables() -> IslandTables:
    """One 40-gene circular replicon under the null (a cheap scan for the stubbed tests)."""
    rng = np.random.default_rng(11)
    genes, tus, replicons = isd.frames([isd.sim_contig(rng, "chr1", 40)])
    return IslandTables(genes, tus, replicons)


def unscanned_tables() -> IslandTables:
    """A 3-gene replicon: fewer than 2 ``min_tus`` TUs, never scanned (M5-06 decision 40)."""
    rng = np.random.default_rng(2)
    genes, tus, replicons = isd.frames([isd.sim_contig(rng, "chr1", 3)])
    return IslandTables(genes, tus, replicons)


def unit_config(**over: Any) -> Any:
    """``islands_data.config`` at B = 199 with the wide detrend span and fixed block of the
    M5-09 plumbing fixtures (the island of ``unit_tables`` is found at p = 0.005)."""
    kw: dict[str, Any] = {"bootstrap": SENS_B, "permutations": SENS_B, "span_frac": 0.30,
                          "block_len": 1}  # fmt: skip
    return isd.config(**(kw | over))


def call_islands_family(t: IslandTables, cfg: Any, rng: np.random.Generator) -> frozenset[str]:
    """``giap.context.sensitivity.islands_family`` on ``t`` (spec signature, keyword tables)."""
    from giap.context import sensitivity as S

    return S.islands_family(t.genes, cfg, rng, tus=t.tus, replicons=t.replicons)  # type: ignore[no-any-return]


def flagged_from_frame(frame: pd.DataFrame) -> frozenset[str]:
    """IDs of every gene of every island row (both directions kept)."""
    return frozenset(
        island_id(str(d), str(g))
        for d, genes in zip(frame["direction"], frame["genes"], strict=True)
        for g in genes
    )


def oracle_frame(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any,
    rng: np.random.Generator,
) -> pd.DataFrame:  # fmt: skip
    """M5-06 ``islands_frame`` as ``context_stats`` builds it, on ``rng``."""
    frame, _ = isd.step(genes, isd.build_series(genes, tus, replicons, cfg), cfg, rng)
    return frame


def oracle_flagged(
    genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame, cfg: Any,
    rng: np.random.Generator,
) -> frozenset[str]:  # fmt: skip
    return flagged_from_frame(oracle_frame(genes, tus, replicons, cfg, rng))


# --------------------------------------------------------------------------- stubbed islands


def fake_islands(rows: Iterable[tuple[str, Sequence[str]]], contig: str = "chr1") -> pd.DataFrame:
    """An ``islands_frame`` result from ``(direction, genes)`` rows (every row significant)."""
    out = []
    for n, (direction, genes) in enumerate(rows, start=1):
        out.append(
            {
                "island_id": f"ISL_{n}", "direction": direction, "contig": contig,
                "start": 100 * n, "end": 100 * n + 50, "tu_first": f"T{n}a",
                "tu_last": f"T{n}b", "n_tus": max(1, len(genes)), "n_genes": len(genes),
                "genes": [str(g) for g in genes], "mean_abundance_log2": 20.0,
                "detected_frac": 1.0, "scan_s": 5.0 if direction == "high" else -5.0,
                "block_len": 1, "stat": isd.STAT, "effect": 1.0, "ci_low": 0.5,
                "ci_high": 1.5, "p_value": 0.005, "q_value": float("nan"),
                "method": isd.METHOD, "null_model": isd.NULL_MODEL, "n_perm": SENS_B,
                "flags": (),
            }
        )  # fmt: skip
    df = pd.DataFrame(out, columns=list(ISLAND_COLUMNS))
    return df.astype({"p_value": "float64", "q_value": "float64", "scan_s": "float64"})


def patch_islands_frame(monkeypatch: Any, frames: Sequence[pd.DataFrame]) -> list[dict[str, Any]]:
    """Replace M5-06 ``islands_frame`` (module attribute and any name imported into
    ``giap.context.sensitivity``) by a stub returning ``frames`` in call order (the last one
    again once they run out); returns the call log."""
    from giap.context import islands as I
    from giap.context import sensitivity as S

    calls: list[dict[str, Any]] = []

    def stub(*args: Any, **kwargs: Any) -> pd.DataFrame:
        calls.append({"args": args, "kwargs": kwargs})
        return frames[min(len(calls), len(frames)) - 1].copy()

    monkeypatch.setattr(I, "islands_frame", stub)
    monkeypatch.setattr(S, "islands_frame", stub, raising=False)
    return calls


def record_scans(monkeypatch: Any) -> list[Any]:
    """Wrap M5-06 ``scan_replicons`` (module and any imported name) to record the config of
    every call; the real scan still runs."""
    from giap.context import islands as I
    from giap.context import sensitivity as S

    real = I.scan_replicons
    cfgs: list[Any] = []

    def wrapper(series: Any, cfg: Any, rng: Any) -> Any:
        cfgs.append(cfg)
        return real(series, cfg, rng)

    monkeypatch.setattr(I, "scan_replicons", wrapper)
    monkeypatch.setattr(S, "scan_replicons", wrapper, raising=False)
    return cfgs


# --------------------------------------------------------------------------- simulated work dir

SIM_SEED = 20261022
SIM_GENES = 1000


def work_config(**over: Any) -> Any:
    """``GiapConfig`` for the work directory: a small ``stats.bootstrap`` (50, so the B = 199
    override is visible), LOESS span 0.30 and block length 1, and nested overrides."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {
        "seed": SIM_SEED,
        "stats": {"bootstrap": 50},
        "context": {"islands": {"detrend_span_frac": 0.30, "block_len": 1}},
    }
    for k, v in over.items():
        raw[k] = {**raw[k], **v} if isinstance(v, Mapping) and isinstance(raw.get(k), dict) else v
    return GiapConfig.model_validate(raw)


def sim_work_dir(
    work: Path,
    *,
    blank_intensity_genes: int = 60,
    evidence: Mapping[str, Any] | None = None,
    seed: int = SIM_SEED,
) -> Path:
    """giapsim work directory (M5-04 ``island_power`` scenario) with the M4-15 gene table."""
    import giapsim.evidence as ev
    import giapsim.expression as ex
    import giapsim.layout as la
    from giap.steps.detectability import run_detectability
    from giap.steps.gene_table import run_gene_table
    from giapsim import scenarios as sc
    from giapsim.contract_tables import write_contract_tables

    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    lcfg, ecfg = sc.island_power(SIM_GENES, 0.5, 15, 4.0, detection="abundance")
    lay = la.simulate_layout(lcfg, rng)
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
    for tsv in (work / "mapping").glob("*.tsv"):
        shutil.copy(tsv, work / tsv.name)
    cfg = work_config()
    run_gene_table(work, cfg)
    run_detectability(work, cfg)
    return work


def sim_genes(work: Path) -> pd.DataFrame:
    from giap.io.genes_tsv import read_genes_parquet

    return read_genes_parquet(work / "genetable" / "genes.parquet")


def sim_m4(work: Path) -> Any:
    from giap.io.work_tables import default_paths, load_m4_inputs

    m4, _ = load_m4_inputs({k: p for k, p in default_paths(work).items() if p.exists()})
    return m4


def sensitivity_config(cfg: Any) -> Any:
    """``cfg`` with the two §8.6 rerun knobs at 199 (what the islands family runs under)."""
    from giap.steps.sensitivity import with_overrides

    return with_overrides(cfg, {"stats.permutations": SENS_B, "stats.bootstrap": SENS_B})


def oracle_sets(work: Path, cfg: Any) -> tuple[frozenset[str], dict[str, frozenset[str]]]:
    """(primary set, variant -> set) of the islands family on ``work`` for the ``basic``
    ``missing=*`` variants, each call on a fresh ``for_step(STREAM, seed)`` generator."""
    from giap.context import sensitivity as S
    from giap.rng import for_step
    from giap.steps.sensitivity import with_overrides

    genes = sim_genes(work)
    m4 = sim_m4(work)
    base = sensitivity_config(cfg)
    primary = oracle_flagged(genes, m4.tus, m4.replicons, base, for_step(STREAM, cfg.seed))
    out: dict[str, frozenset[str]] = {}
    for v in S.variants(cfg, str(genes["detection_tier"].iloc[0])):
        if v.single_peptide_detected or v.name == "top3":
            continue
        g = S.variant_genes(genes, v, cfg)
        vcfg = with_overrides(base, v.overrides)
        out[v.name] = oracle_flagged(g, m4.tus, m4.replicons, vcfg, for_step(STREAM, cfg.seed))
    return primary, out


def read_sensitivity(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
