"""M4-18 example 8 and the step-level checks of ``giap.steps.context_stats.run_context_stats``.

docs/specs/M4-18-regions-stats-table.md: "fixture end to end through ``run_context_stats``:
``regions_stats.tsv`` validates against its schema; header = ``REGIONS_STATS_COLUMNS``", example 9
(byte-identical TSV for the same seed), acceptance "Attack null 4" (no plasmid row has a p-value)
and lead decision 7 (adding a plasmid leaves every M4-17 row unchanged).

Work directory (test-author reading; M4-15 is not merged): the M4-06 contract tables C1-C13 directly
under ``work_dir`` (``giap.io.work_tables.default_paths``, written by the M4-16 fixture
``detection_summary_fixture.sim_fixture``), the gene table at ``work_dir/genetable/genes.parquet``
(M4-15 ``run_detectability`` output; built here from the M4-07 layout and M4-08 truth with the §9.3
columns, ``abundance_imputed`` false) and outputs under ``work_dir/tables/``. The depth TSV
(§5.3) is passed as ``replicon_depth``.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_summary_fixture as fx
import region_tests_data as rd
import regions_stats_data as sd

TABLE_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"
SEED = 7
N_GENES = 300
COPY_20 = math.log2(20.0)


def _cfg() -> Any:
    from giap.config.model import GiapConfig

    return GiapConfig.model_validate({"stats": {"bootstrap": 50}, "seed": 11})


def _gene_table(inputs: Any, copy_log2: float) -> pd.DataFrame:
    from giapsim.expression import ExpressionConfig, simulate_expression
    from giapsim.layout import LayoutConfig, simulate_layout

    lay = simulate_layout(
        LayoutConfig(n_genes_chromosome=N_GENES, plasmids=(25,), linear_contigs=(15,)),
        np.random.default_rng(SEED),
    )  # the layout of fx.sim_fixture(seed=SEED, n_genes=N_GENES)
    ecfg = ExpressionConfig(detection="covariates_only", ar1_rho=0.5,
                            plasmid_copy_log2={"p1": copy_log2})  # fmt: skip
    expr = simulate_expression(lay, ecfg, np.random.default_rng(SEED + 1))
    t = rd.tables_from(lay, expr, "true", genes_c2=inputs.genes, layers=inputs.layers)
    genes = t.genes.copy()
    c2 = inputs.genes.set_index("locus_tag")
    tags = genes["locus_tag"].astype(str)
    genes["hypothetical"] = tags.map(c2["hypothetical"].astype(bool)).to_numpy(bool)
    genes["cog_category"] = fx.cog_categories(inputs.genes, inputs.xrefs).to_numpy()
    genes["abundance_imputed"] = False
    return genes


def _drop_plasmid_regions(work: Path, genes: pd.DataFrame) -> pd.DataFrame:
    """Remove every ``plasmid_*`` region from C7, C8, C6 ``region_ids`` and the gene table,
    keeping each Parquet file's Arrow schema."""
    from giap.io.work_tables import default_paths

    paths = default_paths(work)
    for key in ("regions", "region_genes"):
        tab = pq.read_table(paths[key])
        df = tab.to_pandas()
        df = df[~df["region_id"].astype(str).str.startswith("plasmid_")]
        pq.write_table(pa.Table.from_pandas(df, schema=tab.schema, preserve_index=False),
                       paths[key])  # fmt: skip
    tab = pq.read_table(paths["layers"])
    df = tab.to_pandas()
    df["region_ids"] = [[r for r in v if not str(r).startswith("plasmid_")] for v in
                        df["region_ids"]]  # fmt: skip
    pq.write_table(pa.Table.from_pandas(df, schema=tab.schema, preserve_index=False),
                   paths["layers"])  # fmt: skip
    out = genes.copy()
    out["region_ids"] = [[r for r in v if not str(r).startswith("plasmid_")] for v in
                         out["region_ids"]]  # fmt: skip
    return out


def _build(work: Path, *, copy_log2: float = 0.0, plasmid: bool = True,
           depth: dict[str, float] | None = None) -> Path:  # fmt: skip
    """Write the fixture to ``work``, run the step, return ``tables/regions_stats.tsv``."""
    from giap.steps.context_stats import run_context_stats

    sim = fx.sim_fixture(work, seed=SEED, n_genes=N_GENES)
    genes = _gene_table(sim.inputs, copy_log2)
    if not plasmid:
        genes = _drop_plasmid_regions(work, genes)
    (work / "genetable").mkdir(parents=True, exist_ok=True)
    genes.to_parquet(work / "genetable" / "genes.parquet", index=False)
    dpath = None
    if depth is not None:
        dpath = work / "replicon_depth.tsv"
        lines = ["contig\tmean_depth", *(f"{c}\t{v}" for c, v in depth.items())]
        dpath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    run_context_stats(work, _cfg(), dpath)
    out = work / "tables" / "regions_stats.tsv"
    assert out.exists(), sorted(p.name for p in (work / "tables").glob("*"))
    return out


def _read(path: Path) -> pd.DataFrame:
    from giap.io.tsv import read_tsv

    header, rows = read_tsv(path)
    return pd.DataFrame(rows, columns=header, dtype=object)


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    base = tmp_path_factory.mktemp("context_stats")
    depth = {"chr1": 50.0, "p1": 400.0, "ctg1": 50.0}
    return {
        "plain": _build(base / "plain"),
        "again": _build(base / "again"),
        "depth": _build(base / "depth", depth=depth),
        "copy20": _build(base / "copy20", copy_log2=COPY_20),
        "no_plasmid": _build(base / "no_plasmid", plasmid=False),
    }


def _descriptor() -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / "regions_stats.json",
                         TABLE_DIR / "regions_stats.schema.json") if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    d = json.loads(found[0].read_text(encoding="utf-8"))
    return {k: v for k, v in d.items() if k != "coordinates"}


@pytest.mark.parametrize("name", ["plain", "depth", "copy20"])
def test_example_8_validates_against_schema(runs: dict[str, Path], name: str) -> None:
    frictionless = pytest.importorskip("frictionless")
    out = runs[name]
    schema = frictionless.Schema.from_descriptor(_descriptor())
    resource = frictionless.Resource(path=out.name, basepath=str(out.parent), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_example_8_header_and_content(runs: dict[str, Path]) -> None:
    from giap.context.regions_table import REGIONS_STATS_COLUMNS

    df = _read(runs["plain"])
    assert tuple(df.columns) == tuple(REGIONS_STATS_COLUMNS) == sd.REGIONS_STATS_COLUMNS
    types = set(df["region_type"])
    assert {"prophage", "plasmid", "bgc", "crispr", "amr"} <= types
    pl = df[df["region_type"] == "plasmid"]
    assert set(pl["stat"]) == {"o_minus_e", "hl_shift"} and len(pl) == 2
    assert (pl["label"] == sd.PLASMID_LABEL).all()
    assert (df.loc[df["subset"] != "all", "region_type"] == "prophage").all()
    assert (df["subset"] != "all").any()  # giapsim draws §7.1.7 subsets for prophage members


def test_detection_summary_written(runs: dict[str, Path]) -> None:
    from giap.context.summary import SUMMARY_COLUMNS

    out = runs["plain"].parent / "detection_summary.tsv"
    assert out.exists()
    df = _read(out)
    assert tuple(df.columns) == tuple(SUMMARY_COLUMNS)
    assert "basic" in set(df["section"])


@pytest.mark.parametrize("name", ["plain", "depth", "copy20"])
def test_no_untested_row_has_a_p_value(runs: dict[str, Path], name: str) -> None:
    """Acceptance: no row has ``p_value`` without ``null_model != none``; no ``q_value`` without
    ``p_value``; plasmid rows are never tested (also not inside a chromosome + plasmid index)."""
    df = _read(runs[name])
    none = df["null_model"] == "none"
    assert df.loc[none, "p_value"].isna().all() and df.loc[none, "q_value"].isna().all()
    assert (df.loc[none, "n_perm"] == "0").all()
    assert df.loc[df["q_value"].notna(), "p_value"].notna().all()
    pl = df[df["region_type"] == "plasmid"]
    assert not pl.empty and (pl["null_model"] == "none").all() and pl["p_value"].isna().all()


def test_depth_adds_adjusted_plasmid_row(runs: dict[str, Path]) -> None:
    df = _read(runs["depth"])
    pl = df[df["region_type"] == "plasmid"]
    adj = pl[pl["label"] == sd.ADJUSTED_LABEL]
    assert len(adj) == 1 and adj.iloc[0]["stat"] == "hl_shift"
    assert sd.DEPTH_FLAG in str(adj.iloc[0]["flags"]).split(";")
    hl = pl[(pl["stat"] == "hl_shift") & (pl["label"] != sd.ADJUSTED_LABEL)].iloc[0]
    assert float(adj.iloc[0]["effect"]) == pytest.approx(float(hl["effect"]) - 3.0, abs=1e-9)


def test_attack_4_dosage_is_in_the_plasmid_hl(runs: dict[str, Path]) -> None:
    """With ``plasmid_copy_log2 = log2(20)`` (covariate-only detection, same seed) the plasmid HL
    and its CI move by exactly log2(20) and nothing else changes; no plasmid p-value."""
    a, b = _read(runs["plain"]), _read(runs["copy20"])
    assert list(a["region_id"]) == list(b["region_id"])
    is_hl = (a["region_type"] == "plasmid") & (a["stat"] == "hl_shift")
    for col in ("effect", "ci_low", "ci_high", "hl_shift_log2"):
        np.testing.assert_allclose(b.loc[is_hl, col].astype(float).to_numpy(),
                                   a.loc[is_hl, col].astype(float).to_numpy() + COPY_20,
                                   rtol=0.0, atol=1e-9)  # fmt: skip
    pl_oe = (a["region_type"] == "plasmid") & ~is_hl
    hl_cols = ["hl_shift_log2"]
    pd.testing.assert_frame_equal(a[pl_oe].drop(columns=hl_cols), b[pl_oe].drop(columns=hl_cols))
    chrom = a["region_type"] != "plasmid"
    pd.testing.assert_frame_equal(a[chrom], b[chrom])


@pytest.mark.parametrize("name", ["plain", "depth", "copy20", "no_plasmid"])
def test_row_order(runs: dict[str, Path], name: str) -> None:
    """Lead decision 11: type order, contiguous regions, subset ``all`` first, ``stat`` enum order
    of the schema within (region, subset)."""
    sd.check_row_order(_read(runs[name]))


def test_example_9_byte_identical(runs: dict[str, Path]) -> None:
    assert runs["plain"].read_bytes() == runs["again"].read_bytes()


def test_m4_17_rows_unchanged_by_the_plasmid(runs: dict[str, Path]) -> None:
    """Decision 7: plasmid rows draw from their own stream, so the M4-17 rows (CIs included) are
    identical with and without the plasmid region."""
    a, b = _read(runs["plain"]), _read(runs["no_plasmid"])
    assert (a["region_type"] == "plasmid").any()
    assert not (b["region_type"] == "plasmid").any()
    ka = a[a["region_type"] != "plasmid"].reset_index(drop=True)
    kb = b.reset_index(drop=True)
    pd.testing.assert_frame_equal(ka, kb)


def test_step_inputs_left_untouched(tmp_path: Path) -> None:
    """The step only reads its inputs: running it twice in the same work directory gives the
    same table."""
    out = _build(tmp_path / "w")
    first = out.read_bytes()
    from giap.steps.context_stats import run_context_stats

    run_context_stats(tmp_path / "w", _cfg(), None)
    assert out.read_bytes() == first
