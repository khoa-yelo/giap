"""M9-17: ``null_regions`` and ``null_islands`` from the ``context_stats`` step (V18;
docs/specs/M9-17-plot-support-analysis.md "Support files", "Null draws", examples 4 and 5, the
"region with K = 0 shifts" degenerate case, the tables-unchanged acceptance; lead decisions 1 and
2).

Fixture: sim1000 through the M4 chain (``helpers.determinants_step_data.build_work``: one
circular 1,000-gene chromosome, four regions, chr1 scanned for islands), then the step through
``run_step`` with the inputs ``run_context_stats`` gathers and the support key added.

Readings pinned here (test author; contract in tests/helpers/support_analysis_data.py):

* ``null_regions``: one row per (tested region, stat, shift) in ``regions_stats.tsv`` row order,
  ``shift_k`` ascending; ``null_value`` the statistic the test compares (O - E: ADR-0018's
  studentised sum(d - p) / sqrt(sum p (1 - p)); HL: the amendments' T = HL / sqrt(V)), NA where
  the step's null is NaN; ``observed`` its observed value, so ``perm_p(observed, null)`` is the
  table's ``p_value``. Example 4 on the fixture: a circular replicon with N eligible genes and
  every shift used gives N - 1 rows per stat, ``shift_k`` 1..N - 1.
* **Q2.** ``gippa.steps.context_stats.NULL_REGIONS_MAX_BYTES`` (5 x 2^20) bounds the written
  exhaustive file; at or above it the step writes the 1st-99th percentile grid instead, and the
  header says so: the third column is ``percentile`` (1..99) instead of ``shift_k`` (one §9.2
  header row; no comment line). Grid values are numpy's default (linear) percentiles of the
  non-NaN null.
* ``null_islands``: the genome-wide maximum statistic M_b the FWER uses (``IslandScan.null_max``)
  is the maximum over a replicate's rows of ``max_scan_s``; rows only for scanned replicons;
  ``replicate`` 1..B in each (contig, direction) group; ``observed_max`` constant in a group,
  its maximum over groups the largest observed |S| of any scanned replicon.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import support_analysis_data as sa

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_step_data as s
import m4_chain as mc

SEED = 20261022
KEY = sa.SUPPORT_KEYS["context_stats"]
NAMES = ("null_islands", "null_regions")
B = 100  # stats.bootstrap of the fixture config (island null replicates)


def _cfg(**over: Any) -> Any:
    return mc.config(SEED, stats={"permutations": 99, "bootstrap": B}, **over)


def _inputs(work: Path) -> dict[str, Path]:
    """What ``run_context_stats`` passes (contract tables, gene table, the summary parts, C15)."""
    from gippa.io.work_tables import OPTIONAL_KEYS, default_paths
    from gippa.steps.gene_table import contract_paths

    ins = contract_paths(work)
    ins |= {k: p for k, p in default_paths(work).items()
            if k in OPTIONAL_KEYS and k not in ins and p.exists()}  # fmt: skip
    ins["gene_table"] = work / "genetable" / "genes.parquet"
    parts = work / "genetable" / "detection_summary_parts"
    for key, name in (("detectability_summary", "detectability"), ("tier_summary", "tier")):
        if (parts / f"{name}.parquet").exists():
            ins[key] = parts / f"{name}.parquet"
    if (work / "features.parquet").exists():
        ins["features"] = work / "features.parquet"
    return ins


TABLE_OUTPUTS = {
    "regions_stats": "tables/regions_stats.tsv",
    "detection_summary": "tables/detection_summary.tsv",
    "autocorrelation": "tables/autocorrelation.tsv",
    "islands": "tables/islands.tsv",
    "genes_tsv": "tables/genes.tsv",
    "gene_context": "context/gene_context.parquet",
    "island_diagnostics": "context/island_diagnostics.parquet",
}


def _run(work: Path, out: Path, *, support: bool = True, cfg: Any = None) -> dict[str, Path]:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outs = {k: out / v for k, v in TABLE_OUTPUTS.items()}
    if support:
        outs[KEY] = out / KEY
    reg.run_step("context_stats", cfg or _cfg(), _inputs(work), outs)
    return outs


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return s.build_work(tmp_path_factory.mktemp("m917_ctx_work"), SEED)


@pytest.fixture(scope="module")
def full(work: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _run(work, tmp_path_factory.mktemp("m917_ctx_full"))


@pytest.fixture(scope="module")
def regions(full: dict[str, Path]) -> pd.DataFrame:
    t = mc.read(full["regions_stats"])
    return t[(t["subset"] == "all") & (t["method"] != "descriptive")].reset_index(drop=True)


@pytest.fixture(scope="module")
def nulls(full: dict[str, Path]) -> pd.DataFrame:
    return sa.read(sa.need(full[KEY], "null_regions"))


@pytest.fixture(scope="module")
def genes(work: Path) -> pd.DataFrame:
    return pd.read_parquet(work / "genetable" / "genes.parquet")


# ----------------------------------------------------------------------------- wiring


def test_registry_declares_the_optional_support_output() -> None:
    spec = sa.step_spec("context_stats")
    assert KEY in spec.outputs and KEY in spec.optional_outputs
    for key in ("regions_stats", "detection_summary"):
        assert key in spec.outputs and key not in spec.optional_outputs


def test_tables_bytes_unchanged_by_the_support_writer(
    work: Path, full: dict[str, Path], tmp_path: Path
) -> None:
    bare = _run(work, tmp_path / "bare", support=False)
    assert sa.tables_bytes({k: full[k] for k in TABLE_OUTPUTS}) == sa.tables_bytes(
        {k: bare[k] for k in TABLE_OUTPUTS}
    )


def test_published_under_plots_support(work: Path, tmp_path: Path) -> None:
    w = tmp_path / "work"
    sa.run_with_spec_outputs("context_stats", _cfg(), _inputs(work), w)
    view = sa.publish_and_view(w, tmp_path / "OUT")
    for name in NAMES:
        assert sa.find(tmp_path / "OUT" / "plots" / "support", name) is not None, name
        got = view.support(name)
        assert got is not None and list(got.columns) == list(sa.COLUMNS[name]), name


# ----------------------------------------------------------------------------- schema


def test_the_two_files_and_their_headers(full: dict[str, Path]) -> None:
    assert sa.names(full[KEY]) == sorted(NAMES)
    for name in NAMES:
        path = sa.need(full[KEY], name)
        assert sa.header(path) == list(sa.COLUMNS[name]), name
        sa.check_written_form(path)


@pytest.mark.parametrize("name", NAMES)
def test_validates_against_its_schema(full: dict[str, Path], name: str) -> None:
    import json

    frictionless = pytest.importorskip("frictionless")
    import gippa.schemas

    path = sa.need(full[KEY], name)
    d = json.loads((Path(gippa.schemas.__file__).parent / "plot_support" / f"{name}.json")
                   .read_text(encoding="utf-8"))  # fmt: skip
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ----------------------------------------------------------------------------- null_regions


def test_q2_size_bound_constant() -> None:
    from gippa.steps import context_stats

    assert context_stats.NULL_REGIONS_MAX_BYTES == sa.NULL_REGIONS_MAX_BYTES


def test_regions_are_the_tested_regions(nulls: pd.DataFrame, regions: pd.DataFrame) -> None:
    tested = regions[regions["n_perm"].astype(int) > 0]
    assert len(tested) >= 4, "vacuous without tested regions"
    got = list(dict.fromkeys(zip(nulls["region_id"], nulls["stat"], strict=True)))
    assert got == list(zip(tested["region_id"], tested["stat"], strict=True))


def test_example_4_n_minus_1_shifts_per_stat(
    nulls: pd.DataFrame, regions: pd.DataFrame, genes: pd.DataFrame
) -> None:
    """Circular chromosome, every shift used (N <= exhaustive_shift_max_genes): N - 1 rows per
    stat, ``shift_k`` = 1..N - 1; the O - E null has no NA, the HL null has ``n_perm`` values."""
    from gippa.context.regions import replicon_index

    for _, r in regions.iterrows():
        n = len(replicon_index(genes, r["contig"]))
        grp = nulls[(nulls["region_id"] == r["region_id"]) & (nulls["stat"] == r["stat"])]
        assert grp["shift_k"].astype(int).tolist() == list(range(1, n)), (r["region_id"], n)
        finite = int((grp["null_value"] != "NA").sum())
        assert finite == int(r["n_perm"]), (r["region_id"], r["stat"], finite, r["n_perm"])
        if r["stat"] == "o_minus_e":
            assert finite == n - 1


def test_observed_constant_and_perm_p_is_the_table_p(
    nulls: pd.DataFrame, regions: pd.DataFrame
) -> None:
    from gippa.stats.permutation import perm_p

    for _, r in regions.iterrows():
        grp = nulls[(nulls["region_id"] == r["region_id"]) & (nulls["stat"] == r["stat"])]
        assert grp["observed"].nunique() == 1
        obs = sa.num(grp["observed"].iloc[0])
        null = sa.nums(grp["null_value"])
        null = null[~np.isnan(null)]
        p = float(perm_p(obs, null))
        want = sa.num(r["p_value"])
        assert abs(p - want) <= 2 / (1 + len(null)) + 1e-12, (r["region_id"], r["stat"], p, want)


def test_o_minus_e_null_is_the_studentised_statistic(
    nulls: pd.DataFrame, regions: pd.DataFrame, genes: pd.DataFrame, work: Path
) -> None:
    """ADR-0018: the observed O - E statistic is (O - E) / sqrt(sum p (1 - p)), not the
    reported effect."""
    from gippa.context.regions import replicon_index

    rg = pd.read_parquet(work / "region_genes.parquet")
    for _, r in regions[regions["stat"] == "o_minus_e"].iterrows():
        ix = replicon_index(genes, r["contig"])
        members = set(rg.loc[rg["region_id"] == r["region_id"], "locus_tag"].astype(str))
        p = ix.loc[ix["locus_tag"].isin(members), "p"].to_numpy(float)
        t = sa.num(r["effect"]) / math.sqrt(float(np.sum(p * (1 - p))))
        grp = nulls[(nulls["region_id"] == r["region_id"]) & (nulls["stat"] == "o_minus_e")]
        assert sa.close(sa.num(grp["observed"].iloc[0]), t), (r["region_id"], t)


def test_k_zero_region_has_no_rows(work: Path, tmp_path: Path) -> None:
    """A region whose replicon gives K = 0 shifts (here: a region with no eligible gene, n_perm
    0) has no null rows; the others keep theirs."""
    from gippa.io.work_tables import default_paths

    w = tmp_path / "w"
    import shutil

    shutil.copytree(work, w)
    paths = default_paths(w)
    regs = pd.read_parquet(paths["regions"])
    rg = pd.read_parquet(paths["region_genes"])
    victim = str(regs["region_id"].iloc[0])
    rg = rg[rg["region_id"].astype(str) != victim]
    rg.to_parquet(paths["region_genes"], index=False)
    out = _run(w, tmp_path / "out")
    t = mc.read(out["regions_stats"])
    zero = t[(t["region_id"] == victim) & (t["subset"] == "all")]
    assert len(zero) and (zero["n_perm"].astype(int) == 0).all()
    nulls = sa.read(sa.need(out[KEY], "null_regions"))
    assert victim not in set(nulls["region_id"])
    assert len(set(nulls["region_id"])) >= 1


def test_q2_percentile_grid_above_the_bound(
    work: Path, nulls: pd.DataFrame, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Above the bound: 99 rows per (region, stat), ``percentile`` 1..99, the null's
    percentiles, the same observed value; the header names the grid."""
    from gippa.steps import context_stats

    monkeypatch.setattr(context_stats, "NULL_REGIONS_MAX_BYTES", 1)
    out = _run(work, tmp_path / "grid")
    path = sa.need(out[KEY], "null_regions")
    assert sa.header(path) == list(sa.GRID_COLUMNS)
    sa.check_written_form(path)
    grid = sa.read(path)
    pairs = list(dict.fromkeys(zip(nulls["region_id"], nulls["stat"], strict=True)))
    assert list(dict.fromkeys(zip(grid["region_id"], grid["stat"], strict=True))) == pairs
    for rid, stat in pairs:
        g = grid[(grid["region_id"] == rid) & (grid["stat"] == stat)]
        e = nulls[(nulls["region_id"] == rid) & (nulls["stat"] == stat)]
        assert g["percentile"].astype(int).tolist() == list(range(1, 100))
        v = sa.nums(e["null_value"])
        want = np.percentile(v[~np.isnan(v)], np.arange(1, 100))
        np.testing.assert_allclose(sa.nums(g["null_value"]), want, rtol=1e-4, atol=1e-6)
        assert set(g["observed"]) == set(e["observed"])


def test_q2_exhaustive_below_the_bound(full: dict[str, Path]) -> None:
    path = sa.need(full[KEY], "null_regions")
    assert path.stat().st_size < sa.NULL_REGIONS_MAX_BYTES
    assert sa.header(path)[3] == "shift_k"  # region_id, contig, stat, shift_k


# ----------------------------------------------------------------------------- null_islands


@pytest.fixture(scope="module")
def scan(work: Path, genes: pd.DataFrame) -> Any:
    """``scan_replicons`` as the step calls it (stream ``context_stats.islands``, drawn first)."""
    from gippa.context.islands import scan_replicons, tu_series
    from gippa.context.regions import replicon_index
    from gippa.io.work_tables import M4_INPUT_KEYS, load_m4_inputs
    from gippa.rng import for_step

    cfg = _cfg()
    ins = _inputs(work)
    m4, _ = load_m4_inputs({k: v for k, v in ins.items() if k in M4_INPUT_KEYS})
    span = float(cfg.context.islands.detrend_span_frac)
    series = [
        tu_series(replicon_index(genes, str(c)), m4.tus, contig=str(c),
                  circular=str(t) == "circular", length=int(n), span_frac=span)
        for c, t, n in zip(m4.replicons["contig"], m4.replicons["topology"],
                           m4.replicons["length"], strict=True)
    ]  # fmt: skip
    return scan_replicons(series, cfg, for_step("context_stats.islands", SEED))


def test_null_islands_replicates_and_genome_wide_max(full: dict[str, Path], scan: Any) -> None:
    t = sa.read(sa.need(full[KEY], "null_islands"))
    scanned = {ts.contig for ts, ok in zip(scan.series, scan.scanned, strict=True) if ok}
    assert scanned, "vacuous: no scanned replicon"
    assert set(t["contig"]) <= scanned
    assert (t["direction"] != "NA").all() and (t["direction"] != "").all()
    for _, grp in t.groupby(["contig", "direction"], sort=False):
        assert grp["replicate"].astype(int).tolist() == list(range(1, B + 1))
        assert grp["observed_max"].nunique() == 1
    t = t.assign(_r=t["replicate"].astype(int), _m=sa.nums(t["max_scan_s"]))
    got = t.groupby("_r")["_m"].max().reindex(range(1, B + 1)).to_numpy(float)
    np.testing.assert_allclose(got, np.asarray(scan.null_max, float), rtol=sa.REL)
    obs = max(float(np.nanmax(np.abs(w.s))) for w, ok in zip(scan.windows, scan.scanned,
                                                              strict=True) if ok)  # fmt: skip
    assert sa.close(float(np.nanmax(sa.nums(t["observed_max"]))), obs)


def test_no_scanned_replicon_no_null_islands(work: Path, tmp_path: Path) -> None:
    """``min_tus`` above the TU count: no replicon is scanned (``few_values``), no island null."""
    cfg = _cfg(context={"islands": {"min_tus": 900, "max_tus": 900}})
    out = _run(work, tmp_path / "noscan", cfg=cfg)
    assert sa.find(out[KEY], "null_islands") is None
    assert sa.find(out[KEY], "null_regions") is not None


# ----------------------------------------------------------------------------- determinism


def test_example_5_same_run_twice_byte_identical(
    work: Path, full: dict[str, Path], tmp_path: Path
) -> None:
    again = _run(work, tmp_path / "again")
    assert sa.tree_bytes(full[KEY]) == sa.tree_bytes(again[KEY])


# ----------------------------------------------------------------------------- review round
# Lead decision 8 (PR #272): split-contig regions (review B2, mutant M14), small nulls in the grid
# form (statistician N1, review N5), no HL null under hl_variance_nonpositive (N2 = N8), the
# island observed_max is |S| (review N2, mutant M7), and the grid schema.

GRID_SCHEMA = "null_regions_grid"
N_CHR, N_PL = 160, 130  # both K > 99 shifts, so the grid form takes percentiles on each


def _split_inputs() -> tuple[Any, ...]:
    """Two circular replicons (chr1, plasmid p1) and ``amr_1`` with members on both (tested once
    per contig, flagged ``split_contigs``), plus a one-contig ``bgc_1``."""
    import region_tests_data as rd

    contigs = [
        rd.Contig("chr1", [int(i % 3 == 0) for i in range(N_CHR)], [0.4] * N_CHR),
        rd.Contig("p1", [int(i % 2 == 0) for i in range(N_PL)], [0.6] * N_PL),
    ]
    genes, reps = rd.gene_table(contigs)
    tags = rd.tags_of(genes, "chr1", [0, 3, 9, 12, 40, 41]) + rd.tags_of(genes, "p1", [2, 4, 7, 30])
    genes, regions, rg = rd.region_frames(
        genes, {"amr_1": ("amr", tags), "bgc_1": ("bgc", rd.tags_of(genes, "chr1", range(60, 66)))}
    )
    return genes, regions, rg, reps, rd.config(bootstrap=20, permutations=199)


def _split_written(tmp: Path) -> tuple[pd.DataFrame, list[Any], Path]:
    """(regions_stats table, the step's RegionNulls, the support dir) through the step's writer."""
    import gippa.steps.context_stats as cs
    from gippa.context.regions import region_tests
    from gippa.context.regions_table import plasmid_rows, regions_stats, subset_rows

    genes, regions, rg, reps, cfg = _split_inputs()
    tests, nulls = region_tests(genes, regions, rg, reps, cfg, np.random.default_rng(3))
    plas = plasmid_rows(genes, reps, None, cfg, np.random.default_rng(4))
    table = regions_stats(tests, plas, subset_rows(genes, regions, rg), regions)
    cs._write_null_regions(tmp, table, nulls)
    return table, nulls, tmp


def _split_rows(table: pd.DataFrame) -> pd.DataFrame:
    t = table[(table["region_id"] == "amr_1") & (table["subset"] == "all")]
    assert len(t) == 4 and set(t["contig"]) == {"chr1", "p1"}, t[["contig", "stat"]]
    assert (t["n_perm"].astype(int) > 0).all(), "vacuous: a split row was not tested"
    return t


def _null_of(nulls: list[Any], rid: str, contig: str, stat: str) -> Any:
    hits = [n for n in nulls if (n.region_id, n.contig, n.stat) == (rid, contig, stat)]
    assert len(hits) == 1, (rid, contig, stat, len(hits))
    return hits[0]


def test_split_region_exhaustive_one_null_per_contig(tmp_path: Path) -> None:
    """Review B2 (kills M14): a region on two replicons gives one (region, contig, stat) group
    per contig, each with its own shifts 1..N_c - 1, its own null and its own observed value."""
    from gippa.stats.permutation import perm_p

    table, nulls, out = _split_written(tmp_path)
    path = sa.need(out, "null_regions")
    assert sa.header(path) == list(sa.COLUMNS["null_regions"])
    nr = sa.read(path)
    n_of = {"chr1": N_CHR, "p1": N_PL}
    observed = {}
    for _, r in _split_rows(table).iterrows():
        key = ("amr_1", str(r["contig"]), str(r["stat"]))
        g = nr[(nr["region_id"] == key[0]) & (nr["contig"] == key[1]) & (nr["stat"] == key[2])]
        assert len(g), f"no null rows for {key}"
        assert g["shift_k"].astype(int).tolist() == list(range(1, n_of[key[1]])), key
        mine = _null_of(nulls, *key)
        assert g["observed"].nunique() == 1
        obs = sa.num(g["observed"].iloc[0])
        assert sa.close(obs, mine.observed), (key, obs, mine.observed)
        v = sa.nums(g["null_value"])
        np.testing.assert_allclose(v, mine.null[np.argsort(mine.shifts)], rtol=sa.REL, atol=1e-9)
        assert int(np.count_nonzero(~np.isnan(v))) == int(r["n_perm"]), key
        p = float(perm_p(obs, v[~np.isnan(v)]))
        assert abs(p - float(r["p_value"])) <= 2 / (1 + int(r["n_perm"])) + 1e-12, key
        observed[key] = obs
    for stat in ("o_minus_e", "hl_shift"):
        assert observed[("amr_1", "chr1", stat)] != observed[("amr_1", "p1", stat)], stat


def test_split_region_grid_one_null_per_contig(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review B2, grid form: the same two groups per stat, each the 1st-99th percentiles of its
    own contig's finite null and its own observed value."""
    import gippa.steps.context_stats as cs

    monkeypatch.setattr(cs, "NULL_REGIONS_MAX_BYTES", 0)
    table, nulls, out = _split_written(tmp_path)
    path = sa.need(out, "null_regions")
    assert sa.header(path) == list(sa.GRID_COLUMNS)
    grid = sa.read(path)
    for _, r in _split_rows(table).iterrows():
        key = ("amr_1", str(r["contig"]), str(r["stat"]))
        g = grid[(grid["region_id"] == key[0]) & (grid["contig"] == key[1])
                 & (grid["stat"] == key[2])]  # fmt: skip
        mine = _null_of(nulls, *key)
        fin = mine.null[np.isfinite(mine.null)]
        assert fin.size > 99, "fixture: both contigs give more than 99 finite null values"
        assert g["percentile"].astype(int).tolist() == list(range(1, 100)), key
        np.testing.assert_allclose(sa.nums(g["null_value"]), np.percentile(fin, np.arange(1, 100)),
                                   rtol=sa.REL, atol=1e-9, err_msg=str(key))  # fmt: skip
        assert set(g["observed"]) == {g["observed"].iloc[0]}
        assert sa.close(sa.num(g["observed"].iloc[0]), mine.observed), key


def _exhaustive(groups: dict[tuple[str, str, str], tuple[list[float], float]]) -> pd.DataFrame:
    """An exhaustive ``null_regions`` frame: per key, the null values at shifts 1..K."""
    parts = [pd.DataFrame({"region_id": k[0], "contig": k[1], "stat": k[2],
                           "shift_k": np.arange(1, len(v) + 1, dtype=np.int64),
                           "null_value": np.asarray(v, dtype=float), "observed": obs})
             for k, (v, obs) in groups.items()]  # fmt: skip
    return pd.concat(parts, ignore_index=True)


def _grid_of(ex: pd.DataFrame, key: tuple[str, str, str]) -> pd.DataFrame:
    from gippa.context.support import null_regions_grid

    g = null_regions_grid(ex)
    assert list(g.columns) == list(sa.GRID_COLUMNS)
    return g[(g["region_id"] == key[0]) & (g["contig"] == key[1]) & (g["stat"] == key[2])]


@pytest.mark.parametrize("k", [1, 10, 98, 99])
def test_grid_small_null_keeps_its_exact_values(k: int) -> None:
    """Statistician N1: a group with K <= 99 finite values keeps them, ascending, with
    ``percentile`` the 1-based rank 1..K (not a 99-point interpolated grid)."""
    rng = np.random.default_rng(k)
    vals = rng.normal(size=k).round(3).tolist()
    key = ("r1", "chr1", "o_minus_e")
    ex = _exhaustive({key: (vals, 1.5), ("r2", "chr1", "o_minus_e"): (list(range(200)), 0.5)})
    g = _grid_of(ex, key)
    assert g["percentile"].astype(int).tolist() == list(range(1, k + 1))
    assert g["null_value"].to_numpy(float).tolist() == sorted(vals)
    assert (g["observed"].to_numpy(float) == 1.5).all()
    big = _grid_of(ex, ("r2", "chr1", "o_minus_e"))
    assert big["percentile"].astype(int).tolist() == list(range(1, 100))


def test_grid_small_null_counts_finite_values_only() -> None:
    """N1 with NaN and +-inf: K counts the finite values (7 here), and only they are kept."""
    vals = [3.0, np.nan, -np.inf, 1.0, 2.0, np.inf, 0.5, -1.0, np.nan, 4.0, -2.0]
    key = ("r1", "p1", "hl_shift")
    g = _grid_of(_exhaustive({key: (vals, 0.25)}), key)
    assert g["percentile"].astype(int).tolist() == list(range(1, 8))
    assert g["null_value"].to_numpy(float).tolist() == [-2.0, -1.0, 0.5, 1.0, 2.0, 3.0, 4.0]


def test_grid_drops_infinite_values_before_percentiles() -> None:
    """Lead decision 8: +-inf null values (an HL T with V <= 0) are dropped before the
    percentiles are taken; the grid is the percentiles of the finite values."""
    fin = np.linspace(-3.0, 3.0, 150)
    vals = [*fin.tolist(), np.inf, np.inf, -np.inf, np.nan]
    key = ("r1", "chr1", "hl_shift")
    g = _grid_of(_exhaustive({key: (vals, 0.0)}), key)
    assert g["percentile"].astype(int).tolist() == list(range(1, 100))
    got = g["null_value"].to_numpy(float)
    assert np.isfinite(got).all(), got[~np.isfinite(got)]
    np.testing.assert_allclose(got, np.percentile(fin, np.arange(1, 100)), rtol=0, atol=1e-12)


def test_grid_all_nan_null_writes_na_rows(tmp_path: Path) -> None:
    """Review N5: a group whose null is all NaN gives NA ``null_value`` rows (no crash); the
    other groups are unaffected and the file is written."""
    from gippa.context.support import null_regions_grid
    from gippa.io.plot_support import write_support

    nan_key = ("r1", "chr1", "hl_shift")
    ok_key = ("r1", "chr1", "o_minus_e")
    ex = _exhaustive({ok_key: (list(np.linspace(0, 1, 120)), 0.3), nan_key: ([np.nan] * 50, 0.7)})
    g = null_regions_grid(ex)
    bad = g[g["stat"] == "hl_shift"]
    assert len(bad) >= 1
    assert bad["null_value"].isna().all()
    assert (g[g["stat"] == "o_minus_e"]["percentile"].astype(int) == np.arange(1, 100)).all()
    path = write_support("null_regions", g, tmp_path)
    rows = sa.read(path)
    assert (rows.loc[rows["stat"] == "hl_shift", "null_value"] == "NA").all()


def test_grid_form_validates_against_the_grid_schema(
    work: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Lead decision 8: the grid form has its own schema, ``null_regions_grid.json``."""
    import json

    frictionless = pytest.importorskip("frictionless")
    import gippa.schemas
    from gippa.steps import context_stats

    monkeypatch.setattr(context_stats, "NULL_REGIONS_MAX_BYTES", 1)
    out = _run(work, tmp_path / "grid")
    path = sa.need(out[KEY], "null_regions")
    assert sa.header(path) == list(sa.GRID_COLUMNS)
    spath = Path(gippa.schemas.__file__).parent / "plot_support" / f"{GRID_SCHEMA}.json"
    assert spath.is_file(), f"missing grid schema {spath}"
    d = json.loads(spath.read_text(encoding="utf-8"))
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_no_hl_null_when_hl_variance_nonpositive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Statistician N2 = review N8: an HL row flagged ``hl_variance_nonpositive`` (p NA) has no
    null to plot, so ``region_tests`` returns none for it and ``null_regions`` writes none; the
    O - E null of the same region stays."""
    from gippa.context import regions as R
    from gippa.context.regions_table import plasmid_rows, regions_stats, subset_rows
    from gippa.context.support import null_regions

    orig = R._hl_t

    def nonpos(*args: Any, **kw: Any) -> Any:
        t, _ = orig(*args, **kw)
        return t, True

    monkeypatch.setattr(R, "_hl_t", nonpos)
    genes, regions, rg, reps, cfg = _split_inputs()
    tests, nulls = R.region_tests(genes, regions, rg, reps, cfg, np.random.default_rng(3))
    hl = tests[tests["stat"] == "hl_shift"]
    assert len(hl) and all("hl_variance_nonpositive" in str(f) for f in hl["flags"])
    assert hl["p_value"].isna().all()
    assert (hl["n_perm"].astype(int) > 0).any(), "vacuous: the HL null had values"
    assert [n.stat for n in nulls if n.stat == "hl_shift"] == []
    assert any(n.stat == "o_minus_e" for n in nulls)
    plas = plasmid_rows(genes, reps, None, cfg, np.random.default_rng(4))
    table = regions_stats(tests, plas, subset_rows(genes, regions, rg), regions)
    rows = null_regions(table, nulls)
    assert rows is not None and set(rows["stat"]) == {"o_minus_e"}


def test_null_islands_observed_max_is_the_largest_absolute_s() -> None:
    """Review N2 (kills M7): ``observed_max`` is max |S_w| of the replicon, also when its most
    extreme window is negative; ``direction`` is ``abs``."""
    from types import SimpleNamespace

    from gippa.context.support import null_islands

    scan = SimpleNamespace(
        series=(SimpleNamespace(contig="chr1"), SimpleNamespace(contig="p1")),
        scanned=(True, True),
        windows=(SimpleNamespace(s=np.array([1.0, -4.5, 2.0])),
                 SimpleNamespace(s=np.array([np.nan, 0.5, -0.25]))),
    )  # fmt: skip
    per = [np.array([3.0, 1.0]), np.array([0.5, 0.75])]
    t = null_islands(scan, per)  # type: ignore[arg-type]
    assert t is not None and list(t.columns) == list(sa.COLUMNS["null_islands"])
    obs = t.groupby("contig")["observed_max"].first().to_dict()
    assert obs == {"chr1": 4.5, "p1": 0.5}
    assert set(t["direction"]) == {"abs"}
