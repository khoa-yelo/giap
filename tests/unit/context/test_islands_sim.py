"""M5-06 examples 4-5 on the simulator (docs/specs/M5-06-islands.md "Acceptance", ``slow``).

gippasim M4-07 layout + M4-08 expression with a planted ``IslandSpec`` (tests/sim/gippasim), the
M4-17 gene-table builder (tests/helpers/region_tests_data.py) and the real C5 ``tus`` frame.

Example 4 (a planted high island is found) needs ``detection = "abundance"`` and a large shift:
with gippasim's default ``covariates_only`` detection about 45 % of genes are censored
independently of abundance, which hides even an 8 SD island (ADR-0021 "Power": 0.18 for
20 TUs x 2 SD on 4,500 genes). Example 5 is the mirror case and is *not* recovered at all under
left-censoring - a low island is bounded below by the censored score, which the rest of the
genome shares - so the test asserts the part that must hold under any power: no island is
reported that misses the planted TUs (a false island). The found/not-found counts are printed;
the gate for power is M5-07.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as idat
import region_tests_data as rd

pytestmark = pytest.mark.slow

N_GENES = 1000
ISLAND_FIRST, ISLAND_TUS = 100, 20


def _abundance_beta() -> dict[str, float]:
    """Detection driven by abundance (slope 0.4 per log2), same mean detected share."""
    import gippasim.expression as ex

    slope = 0.4
    return {**dict(ex.DEFAULT_BETA), "log2_abundance": slope,
            "intercept": ex.DEFAULT_BETA["intercept"] - slope * 20.0}  # fmt: skip


def _simulate(
    seed: int, effect: float
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str]]:
    import gippasim.expression as ex
    import gippasim.layout as la

    rng = np.random.default_rng(seed)
    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=N_GENES), rng)
    cfg = ex.ExpressionConfig(
        ar1_rho=0.5,
        detection="abundance",
        beta=_abundance_beta(),
        islands=(ex.IslandSpec("chr1", ISLAND_FIRST, ISLAND_TUS, effect),),
    )
    expr = ex.simulate_expression(lay, cfg, rng)
    t = rd.tables_from(lay, expr, "true")
    chr_tus = lay.tus[lay.tus["contig"] == "chr1"]["tu_id"].astype(str).tolist()
    planted = set(chr_tus[ISLAND_FIRST : ISLAND_FIRST + ISLAND_TUS])
    return t.genes, lay.tus, t.replicons, planted


def _islands(genes: pd.DataFrame, tus: pd.DataFrame, reps: pd.DataFrame, seed: int) -> Any:
    from gippa.context import islands as M
    from gippa.context.regions import replicon_index
    from gippa.rng import for_step

    cfg = idat.config(bootstrap=199, permutations=99, alpha=0.05)
    rng = for_step(idat.STEP, seed)
    series = [
        M.tu_series(
            replicon_index(genes, str(r["contig"])),
            tus,
            contig=str(r["contig"]),
            circular=str(r["topology"]) == "circular",
            length=int(r["length"]),
            span_frac=cfg.context.islands.detrend_span_frac,
        )
        for _, r in reps.iterrows()
    ]
    scan = M.scan_replicons(series, cfg, rng)
    return M.islands_frame(scan, genes, cfg, rng), scan, cfg


def _covered(row: pd.Series, series: Any) -> set[str]:
    tu_ids = {s.contig: s.tu_ids for s in series}[str(row["contig"])]
    first = tu_ids.index(str(row["tu_first"]))
    return {tu_ids[(first + k) % len(tu_ids)] for k in range(int(row["n_tus"]))}


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_example4_planted_high_island_on_the_simulator(seed: int) -> None:
    """Example 4: one ``high`` island overlapping the planted TUs, ``p_value <= alpha``,
    ``method = scan_max_stat`` and ``q_value`` NA."""
    genes, tus, reps, planted = _simulate(seed, 6.0)
    out, scan, cfg = _islands(genes, tus, reps, seed)
    assert len(out) == 1, out[["direction", "tu_first", "tu_last", "p_value"]]
    row = out.iloc[0]
    assert row["direction"] == "high"
    assert float(row["p_value"]) <= cfg.stats.alpha.context
    assert row["method"] == idat.METHOD and idat.is_na(row["q_value"])
    assert row["null_model"] == idat.NULL_MODEL
    assert _covered(row, scan.series) & planted
    assert int(row["n_tus"]) <= cfg.context.islands.max_tus


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_example5_planted_low_island_reports_no_false_island(seed: int) -> None:
    """Example 5 under left-censoring: the planted low island is not recovered (recorded), and
    every island that is reported overlaps the planted TUs."""
    genes, tus, reps, planted = _simulate(seed, -6.0)
    out, scan, _ = _islands(genes, tus, reps, seed)
    found = [r for _, r in out.iterrows() if _covered(r, scan.series) & planted]
    false = [r for _, r in out.iterrows() if not _covered(r, scan.series) & planted]
    print(f"seed {seed}: {len(out)} islands, {len(found)} on the planted low island")
    assert not false, [(r["direction"], r["tu_first"], r["p_value"]) for r in false]
    for r in found:
        assert r["direction"] == "low" and float(r["scan_s"]) < 0
        assert float(r["detected_frac"]) < float(genes["detected"].mean())
