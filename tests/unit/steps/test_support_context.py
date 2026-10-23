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
* **Q2.** ``giap.steps.context_stats.NULL_REGIONS_MAX_BYTES`` (5 x 2^20) bounds the written
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
    from giap.io.work_tables import OPTIONAL_KEYS, default_paths
    from giap.steps.gene_table import contract_paths

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
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

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
    import giap.schemas

    path = sa.need(full[KEY], name)
    d = json.loads((Path(giap.schemas.__file__).parent / "plot_support" / f"{name}.json")
                   .read_text(encoding="utf-8"))  # fmt: skip
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ----------------------------------------------------------------------------- null_regions


def test_q2_size_bound_constant() -> None:
    from giap.steps import context_stats

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
    from giap.context.regions import replicon_index

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
    from giap.stats.permutation import perm_p

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
    from giap.context.regions import replicon_index

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
    from giap.io.work_tables import default_paths

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
    from giap.steps import context_stats

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
    assert sa.header(path)[2] == "shift_k"


# ----------------------------------------------------------------------------- null_islands


@pytest.fixture(scope="module")
def scan(work: Path, genes: pd.DataFrame) -> Any:
    """``scan_replicons`` as the step calls it (stream ``context_stats.islands``, drawn first)."""
    from giap.context.islands import scan_replicons, tu_series
    from giap.context.regions import replicon_index
    from giap.io.work_tables import M4_INPUT_KEYS, load_m4_inputs
    from giap.rng import for_step

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
