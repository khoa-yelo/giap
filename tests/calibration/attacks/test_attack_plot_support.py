"""Attack checks for the M9-17 plot-support nulls (GIAP_SPEC.md §8.5, §10.2 V14(a)/V18;
docs/specs/M9-17-plot-support-analysis.md). Statistician-owned.

The support files draw no new replicates and report no new p, so the attacks are reproduction
attacks: under the structure the region and GSEA tests face (TU blocks, AR(1) between genes,
covariate detectability, split and linear replicons, discrete scores with ties), the written file
must give back the statistic the table reports.

* **R1 region nulls reproduce p** (giapsim, AR(1) rho 0.5, 600-gene chromosome + 60-gene plasmid +
  80-gene linear contig; the chromosome AMR set merged with the plasmid AMR set, so the region is
  ``split_contigs``; a 10-gene BGC block on the linear contig, K < N - 1). Exhaustive and sampled
  shifts: ``perm_p(observed, null_value)`` read back from the written file equals each tested
  row's ``p_value`` in ``regions_stats.tsv`` and the non-NA count equals ``n_perm``.
* **R2 the percentile grid is bounded** (lead decision 2): for K >= 199 the grid's p is within
  0.02 of the table p (reviewer probe: <= 0.01 at K = 599, but up to 0.037 at K = 59-68).
* **R3 small nulls** (lead decision 8, statistician N1): a null of K <= 99 finite values
  (plasmids, linear contigs) keeps its exact values in the grid form, ``percentile`` the 1-based
  rank; it is not up-sampled to a 99-point interpolated grid. A hard check since decision 8 (the
  draft's strict xfail became this test; it fails until the small-null change).
* **G1 gsea_running near-tie** (§7.7.4 calcGseaStat tie rule): with discrete centred scores the
  engine's walk and the cumulative sum can break a ``maxP == -minP`` near-tie (|maxP + minP|
  ~ 1e-16) to opposite sides; the support check must accept the engine's ES there instead of
  raising (probe: 876 of 200,000 small discrete rankings raise at 1268c1a9). Run through
  ``giap.function.running.gsea_running`` itself (lead decision 8's tie clause).

Suites (tests/calibration/conftest.py): R1-R3 ``regions``, G1 ``gsea``; none draws replicates of
its own, so none reads ``tolerances.yaml`` (the bounds are the reproduction tolerances above).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_data as rd

pytestmark = pytest.mark.calibration
REGIONS = pytest.mark.calsuite("regions")
GSEA = pytest.mark.calsuite("gsea")


def _split_layout(seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    t = rd.sim_tables(rng, rd.calibration_layout_config(600, linear_contigs=(80,)),
                      {"ar1_rho": 0.5})  # fmt: skip
    rg, reg = t.region_genes.copy(), t.regions.copy()
    amr = reg[reg["region_type"] == "amr"].sort_values("region_id")
    chr_amr = amr[amr["contig"] == "chr1"]["region_id"].iloc[0]
    p_amr = amr[amr["contig"] != "chr1"]["region_id"].tolist()
    rg.loc[rg["region_id"].isin(p_amr), "region_id"] = chr_amr
    reg = reg[~reg["region_id"].isin(p_amr)]
    lin = t.genes[t.genes["contig"] == t.replicons["contig"].iloc[-1]]
    tags = lin["locus_tag"].iloc[30:40].tolist()
    extra = {"region_id": "bgc_lin", "region_type": "bgc", "contig": lin["contig"].iloc[0],
             "start": int(lin["start"].iloc[30]), "end": int(lin["end"].iloc[39]),
             "source_tool": "attack", "score": np.nan}  # fmt: skip
    reg = pd.concat([reg, pd.DataFrame([extra])], ignore_index=True)
    reg = reg.astype({"start": "Int64", "end": "Int64"})
    rg = pd.concat([rg, pd.DataFrame({"region_id": "bgc_lin", "locus_tag": tags,
                                      "overlap_frac": 1.0, "subset": None,
                                      "core_biosynthetic": False})], ignore_index=True)  # fmt: skip
    of: dict[str, list[str]] = {}
    for r, g in zip(rg["region_id"], rg["locus_tag"], strict=True):
        of.setdefault(g, []).append(r)
    genes = t.genes.copy()
    genes["region_ids"] = [sorted(set(of.get(g, []))) for g in genes["locus_tag"]]
    return genes, reg, rg, t.replicons


def _written(tmp: Path, seed: int, exhaustive: int, grid: bool, monkeypatch: pytest.MonkeyPatch
             ) -> tuple[pd.DataFrame, pd.DataFrame]:  # fmt: skip
    import giap.steps.context_stats as cs
    from giap.context.regions import region_tests
    from giap.context.regions_table import plasmid_rows, regions_stats, subset_rows
    from giap.io.regions_tsv import write_regions_stats

    genes, reg, rg, reps = _split_layout(seed)
    cfg = rd.config(bootstrap=20, permutations=199, exhaustive=exhaustive)
    tests, nulls = region_tests(genes, reg, rg, reps, cfg, np.random.default_rng(seed + 1))
    plas = plasmid_rows(genes, reps, None, cfg, np.random.default_rng(1))
    table = regions_stats(tests, plas, subset_rows(genes, reg, rg), reg)
    write_regions_stats(tmp / "regions_stats.tsv", table)
    if grid:
        monkeypatch.setattr(cs, "NULL_REGIONS_MAX_BYTES", 0)
    cs._write_null_regions(tmp, table, nulls)
    f = next((tmp / "support").glob("null_regions.tsv*"))
    kw = {"sep": "\t", "na_values": ["NA"], "keep_default_na": False}
    return pd.read_csv(f, **kw), pd.read_csv(tmp / "regions_stats.tsv", **kw)


def _tested(rs: pd.DataFrame) -> pd.DataFrame:
    return rs[(rs["subset"] == "all") & rs["p_value"].notna()]


def _null_of(nr: pd.DataFrame, row: pd.Series) -> pd.DataFrame:
    return nr[(nr["region_id"] == row["region_id"]) & (nr["contig"] == row["contig"])
              & (nr["stat"] == row["stat"])]  # fmt: skip


@REGIONS
@pytest.mark.parametrize("seed", [11, 12])
@pytest.mark.parametrize("exhaustive", [20000, 100], ids=["exhaustive", "sampled"])
def test_r1_region_nulls_reproduce_p(tmp_path: Path, seed: int, exhaustive: int,
                                     monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    from giap.stats.permutation import perm_p

    nr, rs = _written(tmp_path, seed, exhaustive, False, monkeypatch)
    tested = _tested(rs)
    assert tested["flags"].astype(str).str.contains("split_contigs").sum() == 4
    assert (tested["region_id"] == "bgc_lin").sum() == 2
    for _, row in tested.iterrows():
        sel = _null_of(nr, row)
        assert len(sel), (row["region_id"], row["contig"], row["stat"])
        z = sel["null_value"].to_numpy(float)
        assert int(np.count_nonzero(~np.isnan(z))) == int(row["n_perm"])
        p = perm_p(float(sel["observed"].iloc[0]), z)
        assert p == pytest.approx(float(row["p_value"]), rel=1e-6), (row["region_id"], row["stat"])


@REGIONS
@pytest.mark.parametrize("seed", [11, 12])
def test_r2_percentile_grid_p_is_close(tmp_path: Path, seed: int,
                                       monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    from giap.stats.permutation import perm_p

    nr, rs = _written(tmp_path, seed, 20000, True, monkeypatch)
    assert "percentile" in nr.columns and "shift_k" not in nr.columns
    small = []
    for _, row in _tested(rs).iterrows():
        sel = _null_of(nr, row)
        if int(row["n_perm"]) >= 199:
            assert len(sel) == 99
            p = perm_p(float(sel["observed"].iloc[0]), sel["null_value"].to_numpy(float))
            assert abs(p - float(row["p_value"])) <= 0.02
        else:
            small.append((row["region_id"], row["contig"], row["stat"], int(row["n_perm"])))
    assert small  # the K < 99 rows (plasmid, linear contig) are the ones R3 is about


@REGIONS
@pytest.mark.parametrize("seed", [11, 12])
def test_r3_small_nulls_keep_their_values_in_grid_form(
    tmp_path: Path, seed: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every tested (region, contig, stat) whose exhaustive null has K <= 99 finite values has,
    in the grid form, exactly those K values ascending with ``percentile`` 1..K; larger nulls
    keep the 99-point grid."""
    (tmp_path / "ex").mkdir()
    (tmp_path / "grid").mkdir()
    ex, rs = _written(tmp_path / "ex", seed, 20000, False, monkeypatch)
    grid, _ = _written(tmp_path / "grid", seed, 20000, True, monkeypatch)
    small = 0
    for _, row in _tested(rs).iterrows():
        e = _null_of(ex, row)["null_value"].to_numpy(float)
        fin = np.sort(e[np.isfinite(e)])
        g = _null_of(grid, row)
        key = (row["region_id"], row["contig"], row["stat"], fin.size)
        if fin.size <= 99:
            small += 1
            assert g["percentile"].astype(int).tolist() == list(range(1, fin.size + 1)), key
            np.testing.assert_allclose(g["null_value"].to_numpy(float), fin, rtol=0, atol=0,
                                       err_msg=str(key))  # fmt: skip
        else:
            assert g["percentile"].astype(int).tolist() == list(range(1, 100)), key
        assert set(g["observed"]) == set(_null_of(ex, row)["observed"]), key
    assert small, "vacuous: no null with K <= 99 (plasmid, linear contig)"


def _near_ties(n_trials: int) -> list[tuple[np.ndarray, np.ndarray]]:
    from giap.stats.fgsea import gsea_es, running_es

    rng = np.random.default_rng(1)
    out = []
    for trial in range(n_trials):
        n = int(rng.integers(4, 40))
        s = rng.choice([-1.0, 0.0, 1.0], size=n) if trial % 2 else rng.integers(-3, 4, n) * 0.1
        s = np.sort(s - np.median(s))[::-1].copy()
        m = int(rng.integers(1, n))
        S = np.sort(rng.choice(np.arange(1, n + 1), size=m, replace=False)).astype(np.int64)
        es, _ = gsea_es(s, S)
        w = running_es(s, S)
        hi, lo = max(0.0, float(w.max())), min(0.0, float(w.min()))
        if abs(hi + lo) <= 1e-9 and es != 0.0 and np.sign(es) != np.sign(hi if hi > -lo else lo):
            out.append((s, S))
    return out


def _step_inputs(s: np.ndarray, S: np.ndarray) -> tuple[object, list[object], pd.DataFrame]:
    """``gsea_running``'s inputs for one set on the ranking ``s`` (TU_i = rank i, gene g_i), the
    table's ``es`` the engine's (``gsea_es``) on that ranking."""
    from types import SimpleNamespace

    from giap.function.gsea import TuRanking
    from giap.stats.fgsea import gsea_es

    ids = tuple(f"TU_{i}" for i in range(1, len(s) + 1))
    rank = TuRanking(ids, np.asarray(s, np.float64), {t: (f"g_{t[3:]}",) for t in ids})
    colls = [SimpleNamespace(collection="go", sets={"T": frozenset(f"g_{int(m)}" for m in S)})]
    es, _ = gsea_es(s, S)
    gsea = pd.DataFrame({"collection": ["go"], "term": ["T"], "p_value": [0.01],
                         "q_value": [0.01], "nes": [1.0], "es": [float(es)]})  # fmt: skip
    gsea.attrs["rank_labels"] = list(ids)
    return rank, colls, gsea


@GSEA
def test_g1_gsea_running_accepts_engine_side_of_a_near_tie(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The step must not fail on a ranking whose ES the engine computed (the walk is right; only
    the side of a 1e-16 tie differs). Fails at 1268c1a9 (giap/function/running.py:73)."""
    from giap.function import running as R

    cases = _near_ties(20000)
    assert cases, "no near-tie found; the probe lost its teeth"
    raised = []
    for s, S in cases[:50]:
        rank, colls, gsea = _step_inputs(s, S)
        monkeypatch.setattr(R, "tu_ranking", lambda genes, rank=rank: rank)
        try:
            out = R.gsea_running(pd.DataFrame(), colls, gsea)  # type: ignore[arg-type]
        except ValueError as exc:
            raised.append(str(exc))
            continue
        assert out is not None and len(out) == len(s)
    assert not raised, f"{len(raised)} of {min(len(cases), 50)} near-ties raise: {raised[:3]}"
