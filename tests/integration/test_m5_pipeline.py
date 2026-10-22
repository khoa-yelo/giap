"""M5 end to end on giapsim fixtures (docs/specs/M5-11-integration-golden.md examples 1-5).

``build_m4_fixture`` (M4-20, extended by M5-11 with ``oric_frac`` and ``islands``) writes the
contract tables C1-C15, then ``run_gene_table`` -> ``run_detectability`` -> ``run_context_stats``
-> ``run_sensitivity`` (tests/helpers/m5_chain.py) write the six ``tables/`` files; each validates
with ``frictionless`` against its §9.5 resource.

* Example 1 -- ``sim50k`` at the M4 golden config: ``autocorrelation.tsv`` and ``islands.tsv``
  written and valid (``islands.tsv`` header-only, decision 39), and the M4 tables unchanged:
  ``detection_summary.tsv`` and ``sensitivity.tsv`` byte for byte against
  ``tests/data/golden/sim50k/m4/``, ``regions_stats.tsv`` with the M4-20 decision 9 float
  tolerance on its four lstsq columns, ``genes.tsv`` byte for byte outside the three columns M5
  fills (``island_id``, ``replichore``, ``replichore_position``: whether sim50k gets an oriC is a
  lead decision, see the report).
* Example 2 -- the island case ``m5_chain.ISLAND`` at the shipped defaults: exactly one ``high``
  island overlapping the planted TUs, the l' check of decision 39, three ori-ter rows labelled
  ``descriptive``, and the genes.tsv context columns consistent with islands.tsv. Byte equality
  with the golden files is tests/golden/test_m5_golden.py.
* Example 3 -- TE7: the same seed twice gives byte-identical ``tables/*`` and ``context/*.parquet``
  (sim50k here; the island case is ``slow``).
* Example 4 -- the island case with a different ``seed`` in the config (the fixture fixed): the
  planted island is still reported, sharing at least one TU (one seed here, five ``slow``).
* Example 5 -- the recorded K-12 Bakta JSON through the ``parse_bakta`` step and M5-08: the oriC is
  the circular mean of the two calls and every chromosome gene gets a replichore; no proteomics.
"""

from __future__ import annotations

import logging
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m5_chain as m5

REPO = Path(__file__).resolve().parents[2]
M4_GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "m4"
#: M4-20 lead decision 9 (tests/golden/test_m4_golden.py)
NUMERIC_COLUMNS = ("hl_shift_log2", "effect", "ci_low", "ci_high")
REL_TOL, ABS_TOL = 1e-9, 1e-12
#: extra config seeds of example 4 (``slow``)
MORE_SEEDS = (1, 2, 3, 4, 5)


@pytest.fixture(scope="module")
def db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return mc.test_db(tmp_path_factory.mktemp("m5_db") / "db")


@pytest.fixture(scope="module")
def sim50k(tmp_path_factory: pytest.TempPathFactory, db: Path, golden_seed: int) -> dict[str, Any]:
    work = m5.build_sim50k(tmp_path_factory.mktemp("m5_sim50k") / "work", golden_seed)
    return {"work": work, "tables": m5.run_m5(work, m5.sim50k_config(golden_seed), db)}


@pytest.fixture(scope="module")
def island(tmp_path_factory: pytest.TempPathFactory, db: Path, golden_seed: int) -> dict[str, Any]:
    work = m5.build_island(tmp_path_factory.mktemp("m5_island") / "work", golden_seed)
    return {"work": work, "tables": m5.run_m5(work, m5.island_config(golden_seed), db)}


def _close(want: str, got: str) -> bool:
    if want in ("NA", "") or got in ("NA", ""):
        return want in ("NA", "") and got in ("NA", "")
    a, b = float(want), float(got)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)


# --------------------------------------------------------------------------- example 1


def test_example_1_writes_and_validates_every_table(sim50k: dict[str, Any]) -> None:
    tables = sim50k["tables"]
    missing = [n for n, p in tables.items() if not p.is_file()]
    assert not missing, f"the chain wrote no {missing}"
    for name, stem in m5.ALL_TABLES.items():
        mc.validate(tables[name], stem)


def test_example_1_islands_is_header_only_on_50_kb(sim50k: dict[str, Any]) -> None:
    lines = sim50k["tables"]["islands.tsv"].read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1, lines[:3]
    assert lines[0].split("\t") == [f["name"] for f in mc.descriptor("islands")["fields"]]


def test_example_1_autocorrelation_has_both_statistics_and_nulls(sim50k: dict[str, Any]) -> None:
    ac = mc.read(sim50k["tables"]["autocorrelation.tsv"])
    moran = ac[ac["statistic"] != "oriter_slope"]
    assert len(moran) > 0, "no Moran's I / join-count rows on sim50k"
    assert set(moran["contig"]) == {"chr1"}
    assert set(moran["null_model"]) >= {"gene_label", "tu_block"}
    assert len(set(moran["statistic"])) == 2, set(moran["statistic"])


@pytest.mark.parametrize("name", ["detection_summary.tsv", "sensitivity.tsv"])
def test_example_1_m4_tables_are_byte_identical_to_the_m4_golden(
    sim50k: dict[str, Any], name: str
) -> None:
    want, got = M4_GOLDEN / name, sim50k["tables"][name]
    assert got.read_bytes() == want.read_bytes(), mc.first_diff(want, got)


def test_example_1_regions_stats_matches_the_m4_golden(sim50k: dict[str, Any]) -> None:
    """Byte for byte except decision 9's four lstsq columns (rel 1e-9 / abs 1e-12)."""
    want, got = (
        mc.read(M4_GOLDEN / "regions_stats.tsv"),
        mc.read(sim50k["tables"]["regions_stats.tsv"]),
    )
    assert list(got.columns) == list(want.columns) and len(got) == len(want)
    for col in want.columns:
        pairs = list(zip(want[col], got[col], strict=True))
        if col in NUMERIC_COLUMNS:
            bad = [(i, a, b) for i, (a, b) in enumerate(pairs) if not _close(a, b)]
        else:
            bad = [(i, a, b) for i, (a, b) in enumerate(pairs) if a != b]
        assert not bad, f"regions_stats {col} (row, want, got): {bad[:5]}"


def test_example_1_genes_tsv_is_unchanged_outside_the_m5_columns(sim50k: dict[str, Any]) -> None:
    want = m5.drop_columns(M4_GOLDEN / "genes.tsv", m5.CONTEXT_COLUMNS)
    got = m5.drop_columns(sim50k["tables"]["genes.tsv"], m5.CONTEXT_COLUMNS)
    assert got == want, m5.first_diff(want, got)


def test_example_1_no_island_ids_on_50_kb(sim50k: dict[str, Any]) -> None:
    g = mc.read(sim50k["tables"]["genes.tsv"])
    assert set(g["island_id"]) == {"NA"}


# --------------------------------------------------------------------------- example 2


def test_example_2_writes_and_validates_every_table(island: dict[str, Any]) -> None:
    for name, stem in m5.ALL_TABLES.items():
        assert island["tables"][name].is_file(), name
        mc.validate(island["tables"][name], stem)


def test_example_2_one_high_island_over_the_planted_tus(island: dict[str, Any]) -> None:
    islands = mc.read(island["tables"]["islands.tsv"])
    assert not m5.island_problems(island["work"], islands), islands.to_dict("records")
    row = islands.iloc[0]
    assert row["island_id"] == "ISL_1"
    assert row["contig"] == "chr1"
    assert float(row["p_value"]) <= 0.05, row["p_value"]
    assert int(row["n_tus"]) <= 20  # context.islands.max_tus


def test_example_2_l_prime_check(island: dict[str, Any]) -> None:
    """Decision 39 / ADR-0021: auto block length, l' = min(20, 2 l-hat + 1) < 20 TUs."""
    islands = mc.read(island["tables"]["islands.tsv"])
    assert not m5.l_prime_problems(island["work"], islands), m5.diagnostics(island["work"])


def test_example_2_oriter_rows_are_descriptive(island: dict[str, Any]) -> None:
    ac = mc.read(island["tables"]["autocorrelation.tsv"])
    ori = ac[ac["statistic"] == "oriter_slope"]
    assert sorted(ori["stat"]) == ["oriter_slope_left", "oriter_slope_pooled",
                                   "oriter_slope_right"]  # fmt: skip
    assert set(ori["contig"]) == {"chr1"}
    for _, r in ori.iterrows():
        flags = r["flags"].split(";")
        assert "descriptive" in flags and "no_oric" not in flags, r.to_dict()
        assert r["value"] not in ("NA", ""), r.to_dict()


def test_example_2_genes_tsv_context_columns(island: dict[str, Any]) -> None:
    """island_id is set exactly on the island's genes; every chr1 gene has a replichore and a
    position in [0, 1]; both arms are populated around the oriC at 0.8."""
    g = mc.read(island["tables"]["genes.tsv"])
    isl = mc.read(island["tables"]["islands.tsv"])
    members = m5.island_genes(isl.iloc[0])
    tagged = set(g.loc[g["island_id"] != "NA", "locus_tag"])
    assert tagged == members
    assert set(g.loc[g["island_id"] != "NA", "island_id"]) == {isl["island_id"].iloc[0]}
    chr1 = g[g["contig"] == "chr1"]
    assert set(chr1["replichore"]) == {"right", "left"}
    pos = chr1["replichore_position"].astype(float).to_numpy()
    assert np.all((pos >= 0.0) & (pos <= 1.0))


def test_example_2_features_carry_one_oric(island: dict[str, Any]) -> None:
    f = pd.read_parquet(island["work"] / "features.parquet")
    oric = f[f["type"].astype(str) == "oric"]
    assert len(oric) == 1 and str(oric["contig"].iloc[0]) == "chr1"
    length = int(pd.read_parquet(island["work"] / "replicons.parquet")
                 .set_index("contig").loc["chr1", "length"])  # fmt: skip
    mid = (int(oric["start"].iloc[0]) + int(oric["end"].iloc[0])) / 2 / length
    assert abs(mid - m5.ORIC_FRAC) < 0.02, mid


# --------------------------------------------------------------------------- example 3 (TE7)


def test_example_3_te7_sim50k(
    sim50k: dict[str, Any], tmp_path: Path, db: Path, golden_seed: int
) -> None:
    work = m5.build_sim50k(tmp_path / "again", golden_seed)
    m5.run_m5(work, m5.sim50k_config(golden_seed), db)
    first, second = m5.digests(sim50k["work"]), m5.digests(work)
    assert set(first) >= {f"tables/{n}" for n in m5.ALL_TABLES}, sorted(first)
    assert first == second, sorted(k for k in first if first[k] != second.get(k))


@pytest.mark.slow
def test_example_3_te7_island(
    island: dict[str, Any], tmp_path: Path, db: Path, golden_seed: int
) -> None:
    work = m5.build_island(tmp_path / "again", golden_seed)
    m5.run_m5(work, m5.island_config(golden_seed), db)
    first, second = m5.digests(island["work"]), m5.digests(work)
    assert first == second, sorted(k for k in first if first[k] != second.get(k))


# --------------------------------------------------------------------------- example 4


def _check_seed(tmp: Path, db: Path, fixture_seed: int, seed: int) -> None:
    """The ``island`` fixture rebuilt (same fixture seed), config ``seed`` changed."""
    work = m5.build_island(tmp / f"seed{seed}", fixture_seed)
    tables = m5.run_m5(work, m5.island_config(seed), db)
    islands = mc.read(tables["islands.tsv"])
    planted = set(m5.planted_tus(work))
    shared = [r for _, r in islands.iterrows()
              if str(r["direction"]) == "high" and set(m5.tu_span(r, work)) & planted]  # fmt: skip
    assert shared, f"config seed {seed}: planted island not reported: {islands.to_dict('records')}"


def test_example_4_another_seed_still_reports_the_island(
    tmp_path: Path, db: Path, golden_seed: int
) -> None:
    _check_seed(tmp_path, db, golden_seed, golden_seed + 1)


@pytest.mark.slow
@pytest.mark.parametrize("seed", MORE_SEEDS)
def test_example_4_more_seeds(tmp_path: Path, db: Path, golden_seed: int, seed: int) -> None:
    _check_seed(tmp_path, db, golden_seed, seed)


# --------------------------------------------------------------------------- oriC is M4-neutral


def test_the_oric_does_not_change_the_m4_tables(tmp_path: Path, db: Path) -> None:
    """Spec Interface: the oriC is a new file and the gradient stays 0, so a fixture with and one
    without the oriC give the same M4 tables (and the same genes.tsv outside the M5 columns)."""
    from giapsim.fixtures import build_m4_fixture

    cfg = mc.config(11, stats={"permutations": 199, "bootstrap": 200})
    out = {}
    for frac in (None, m5.ORIC_FRAC):
        work = tmp_path / f"oric_{frac}"
        work.mkdir()
        build_m4_fixture(work, seed=11, size="sim1000", oric_frac=frac)
        out[frac] = m5.run_m5(work, cfg, db)
    a, b = out[None], out[m5.ORIC_FRAC]
    for name in ("detection_summary.tsv", "regions_stats.tsv", "sensitivity.tsv"):
        assert a[name].read_bytes() == b[name].read_bytes(), name
    assert m5.drop_columns(a["genes.tsv"], m5.CONTEXT_COLUMNS) == m5.drop_columns(
        b["genes.tsv"], m5.CONTEXT_COLUMNS
    )
    assert set(mc.read(a["genes.tsv"])["replichore"]) == {"NA"}
    assert set(mc.read(b["genes.tsv"])["replichore"]) >= {"right", "left"}


# --------------------------------------------------------------------------- example 5


def test_example_5_k12_bakta_json_resolves_the_oric(tmp_path: Path) -> None:
    """Recorded K-12 (Bakta 1.12.2 full): ``parse_bakta`` writes C15 with two oriC calls; M5-08
    resolves one origin (their circular mean, 3,904,925.5 +- 1) and labels every chromosome gene
    ``right`` / ``left``; nothing needs proteomics."""
    from giap.config.model import GiapConfig
    from giap.context.oriter import gene_replichores, oric_position
    from giap.steps.annotation import parse_bakta
    from giap.steps.registry import StepContext

    rec = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full" / "k12"
    names = ("genes", "features", "replicons", "rrna16s", "trna_counts", "xrefs",
             "annotation_meta")  # fmt: skip
    ext = {"rrna16s": ".fna", "trna_counts": ".tsv", "annotation_meta": ".json"}
    outputs = {n: tmp_path / f"{n}{ext.get(n, '.parquet')}" for n in names}
    ctx = StepContext(
        cfg=GiapConfig.model_validate({"seed": 1}),
        inputs={"bakta_json": rec / "k12_mg1655.json.gz",
                "genome": REPO / "tests" / "data" / "genomes" / "k12_mg1655.fna.gz"},
        outputs=outputs, rng=np.random.default_rng(0), log=logging.getLogger("test"),
    )  # fmt: skip
    parse_bakta(ctx)
    feats = pd.read_parquet(outputs["features"])
    reps = pd.read_parquet(outputs["replicons"])
    genes = pd.read_parquet(outputs["genes"])
    assert int((feats["type"].astype(str) == "oric").sum()) == 2
    chrom = reps[reps["kind"].astype(str) == "chromosome"]
    assert list(chrom["contig"]) == ["U00096.3"]
    length = int(chrom["length"].iloc[0])
    got = oric_position(feats, "U00096.3", length)
    assert got is not None and got == pytest.approx((3_884_032.0 + 3_925_819.0) / 2, abs=1.0)
    rc = gene_replichores(genes, reps, feats)
    on_chr = rc[genes["contig"].astype(str).to_numpy() == "U00096.3"]
    assert len(on_chr) > 4000
    assert set(on_chr["replichore"]) == {"right", "left"}
    assert on_chr["replichore_position"].between(0.0, 1.0).all()
