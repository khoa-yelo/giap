"""M4-17 region tests: edge cases, flags, CIs, BH family, RNG order and real producer frames.

docs/specs/M4-17-region-tests.md "Edge cases", "Public interface" and the lead decisions of
docs/specs/M4-index.md (decision 5: shifts index eligible genes; BH within (region type, statistic);
``few_shifts`` flag when K < 19). Hand-built tables from tests/helpers/region_tests_data.py; the
round-trip test writes giapsim's C1-C8 and reads them with M4-06 ``load_m4_inputs``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import region_tests_data as rd
from region_tests_data import Contig, flag_set, gene_table, pick, region_frames, tags_of


def _regions() -> Any:
    from giap.context import regions as R

    return R


def _tables(
    contigs: list[Contig], members: dict[str, tuple[str, list[str]]]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes, replicons = gene_table(contigs)
    genes, regions, rg = region_frames(genes, members)
    return genes, regions, rg, replicons


def _run(
    contigs: list[Contig], members: dict[str, tuple[str, list[str]]], cfg: Any = None, seed: int = 1
) -> pd.DataFrame:
    genes, regions, rg, replicons = _tables(contigs, members)
    cfg = cfg if cfg is not None else rd.config()
    return _regions().test_all_regions(
        genes, regions, rg, replicons, cfg, np.random.default_rng(seed)
    )


def _alt(n: int, p: float = 0.5) -> tuple[list[int], list[float]]:
    return [i % 2 for i in range(n)], [p] * n


def _sorted(res: pd.DataFrame) -> pd.DataFrame:
    out = res.sort_values(["region_id", "contig", "stat"]).reset_index(drop=True)
    out["flags"] = [";".join(sorted(flag_set(f))) for f in out["flags"]]
    return out


# --------------------------------------------------------------------------- flags


@pytest.mark.parametrize(("n", "few"), [(19, True), (20, False)])
def test_few_shifts_flag(n: int, few: bool) -> None:
    """Lead decision 5 (Q3): flag ``few_shifts`` when K < 19 (the smallest p exceeds 0.05)."""
    d, p = _alt(n)
    c = Contig("p1", d, p)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "p1", [3, 4, 5]))})
    row = pick(res, "bgc_1", "o_minus_e")
    assert row["n_perm"] == n - 1
    assert ("few_shifts" in flag_set(row["flags"])) is few


def test_no_valid_shift_is_few_shifts() -> None:
    c = Contig("ctg1", [1, 0, 1, 1, 1], [0.5] * 5, circular=False)
    genes, _ = gene_table([c])
    res = _run([c], {"amr_1": ("amr", tags_of(genes, "ctg1", [0, 4]))})
    assert "few_shifts" in flag_set(pick(res, "amr_1", "o_minus_e")["flags"])


def test_flags_stay_in_region_vocabulary() -> None:
    """§9.4 "Region and context" plus the accepted ``few_shifts`` (lead decision 5)."""
    allowed = {"small_region", "dosage_plus_regulation", "descriptive", "no_oric", "uncalibrated",
               "few_shifts"}  # fmt: skip
    c = Contig("ctg1", [1, 0, 1, 1, 1], [0.5] * 5, circular=False)
    genes, _ = gene_table([c])
    res = _run([c], {"amr_1": ("amr", tags_of(genes, "ctg1", [0, 4]))})
    for f in res["flags"]:
        assert flag_set(f) <= allowed


# --------------------------------------------------------------------------- edge cases


def test_region_on_two_contigs_tested_per_contig() -> None:
    """Members on two contigs -> each contig tested separately, region_id kept (AMR sets spanning
    replicons are split the same way)."""
    d1, p1 = _alt(30)
    d2, p2 = _alt(25)
    contigs = [Contig("chr1", d1, p1), Contig("p1", d2, p2)]
    genes, _ = gene_table(contigs)
    tags = tags_of(genes, "chr1", [2, 9, 17, 23]) + tags_of(genes, "p1", [1, 5, 8])
    res = _run(contigs, {"amr_1": ("amr", tags)})
    assert sorted(res["contig"].unique()) == ["chr1", "p1"]
    assert (res["region_id"] == "amr_1").all() and len(res) == 4
    chr_row = pick(res, "amr_1", "o_minus_e", contig="chr1")
    p_row = pick(res, "amr_1", "o_minus_e", contig="p1")
    assert chr_row["n_genes"] == 4 and chr_row["n_perm"] == 29
    assert p_row["n_genes"] == 3 and p_row["n_perm"] == 24
    assert p_row["expected"] == pytest.approx(1.5)
    assert p_row["observed"] == pytest.approx(sum(d2[i] for i in [1, 5, 8]))


def test_all_partial_members() -> None:
    """Members all partial -> n_eligible 0 -> rows with NA statistics and small_region."""
    n = 10
    part = [True] + [False] * (n - 2) + [True]
    c = Contig("ctg1", [1] * n, [0.5] * n, circular=False, partial=part)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "ctg1", [0, n - 1]))})
    for stat in rd.STATS:
        row = pick(res, "bgc_1", stat)
        assert row["n_genes"] == 2 and row["n_eligible"] == 0
        assert rd.is_na(row["effect"]) and rd.is_na(row["p_value"])
        assert row["n_perm"] == 0
        assert "small_region" in flag_set(row["flags"])


def test_partial_member_counts_in_n_genes_only() -> None:
    n = 12
    part = [True] + [False] * (n - 2) + [True]
    d = [1] * n
    p = [0.25] * n
    c = Contig("ctg1", d, p, circular=False, partial=part)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "ctg1", [0, 1, 2, 3]))})
    row = pick(res, "bgc_1", "o_minus_e")
    assert row["n_genes"] == 4 and row["n_eligible"] == 3
    assert row["observed"] == pytest.approx(3.0) and row["expected"] == pytest.approx(0.75)
    # eligible index: N = 10, I = {0, 1, 2} -> rightward shifts 1..7 only
    assert row["n_perm"] == len(rd.expected_shifts(10, [0, 1, 2], circular=False)) == 7


def test_identical_statistics_ties_counted() -> None:
    """All p_detect equal and nobody detected: every shift ties the observed O - E -> p = 1."""
    n = 40
    c = Contig("chr1", [0] * n, [0.3] * n)
    genes, _ = gene_table([c])
    res = _run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", range(5, 12)))})
    assert pick(res, "prophage_1", "o_minus_e")["p_value"] == pytest.approx(1.0)


def test_detected_without_abundance_leaves_hl_both_sides() -> None:
    """Abundance policy ``exclude``: detected genes without abundance drop out of HL on both sides
    but still count in O - E."""
    n = 30
    d = [1] * n
    a = [float(10 + (i * 7) % 11) for i in range(n)]
    for i in (3, 4, 20, 21):
        a[i] = math.nan
    c = Contig("chr1", d, [0.5] * n, a=a)
    genes, _ = gene_table([c])
    members = {"bgc_1": ("bgc", tags_of(genes, "chr1", range(2, 8)))}
    genes, regions, rg = region_frames(genes, members)
    replicons = gene_table([c])[1]
    cfg = rd.config(abundance_missing="exclude")
    res = _regions().test_all_regions(genes, regions, rg, replicons, cfg, np.random.default_rng(1))
    inside = [a[i] for i in range(2, 8) if not math.isnan(a[i])]
    outside = [a[i] for i in range(n) if not 2 <= i < 8 and not math.isnan(a[i])]
    assert pick(res, "bgc_1", "hl_shift")["effect"] == pytest.approx(rd.naive_hl(inside, outside))
    oe = pick(res, "bgc_1", "o_minus_e")
    assert oe["observed"] == pytest.approx(6.0) and oe["n_eligible"] == 6


def test_hl_n_perm_counts_valid_shifts_only() -> None:
    """Shifts whose HL is NaN (no detected+observed gene inside) are dropped: K counts valid
    shifts. Circular N = 10, I = {0}; detected genes 0, 3, 6 -> HL valid at k = 3, 6 only."""
    d = [1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
    c = Contig("chr1", d, [0.3] * 10)
    genes, _ = gene_table([c])
    res = _run([c], {"crispr_1": ("crispr", tags_of(genes, "chr1", [0]))})
    assert pick(res, "crispr_1", "hl_shift")["n_perm"] == 2
    assert pick(res, "crispr_1", "o_minus_e")["n_perm"] == 9


def test_hl_compares_with_same_replicon_only() -> None:
    """HL "other detected" genes are the replicon's own: a 20x-copy plasmid region is compared with
    the plasmid's other detected genes, never with the chromosome (attack 4)."""
    rng = np.random.default_rng(4)
    ch = Contig("chr1", [1] * 50, [0.5] * 50, a=rng.normal(10, 1, 50).tolist())
    pa = rng.normal(14.3, 1, 30).tolist()
    pl = Contig("p1", [1] * 30, [0.5] * 30, a=pa)
    genes, _ = gene_table([ch, pl])
    res = _run([ch, pl], {"amr_1": ("amr", tags_of(genes, "p1", [4, 11, 19]))})
    want = rd.naive_hl(
        [pa[i] for i in (4, 11, 19)], [pa[i] for i in range(30) if i not in (4, 11, 19)]
    )
    row = pick(res, "amr_1", "hl_shift")
    assert row["effect"] == pytest.approx(want)
    assert row["n_perm"] == 29


def test_plasmid_regions_are_not_tested_here() -> None:
    """``plasmid`` regions belong to M4-18: no rows; only the four §7.3.5(b) types appear."""
    d, p = _alt(30)
    contigs = [Contig("chr1", d, p), Contig("p1", *_alt(20))]
    genes, _ = gene_table(contigs)
    res = _run(
        contigs,
        {
            "plasmid_1": ("plasmid", tags_of(genes, "p1", range(20))),
            "prophage_1": ("prophage", tags_of(genes, "chr1", range(3, 9))),
        },
    )
    assert set(res["region_type"]) == {"prophage"}
    assert "plasmid_1" not in set(res["region_id"])


def test_gene_in_two_regions_counts_in_both() -> None:
    """A gene in two regions (M2-20 example 4) is a member of each region's index set."""
    d, p = _alt(40)
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    pro = tags_of(genes, "chr1", range(10, 20))
    amr = tags_of(genes, "chr1", [12, 30, 35])
    res = _run([c], {"prophage_1": ("prophage", pro), "amr_1": ("amr", amr)})
    assert pick(res, "prophage_1", "o_minus_e")["observed"] == pytest.approx(5.0)
    assert pick(res, "amr_1", "o_minus_e")["n_eligible"] == 3


# --------------------------------------------------------------------------- CIs


def test_o_minus_e_ci_counts_only_region_members_of_a_tu() -> None:
    """O - E CI resamples the region's TUs, counting only region members in each TU: a region
    inside one TU (with outside members) has a degenerate CI at the effect."""
    n = 30
    d, p = _alt(n)
    d[0:5], p[0:5] = [0, 0, 0, 1, 1], [0.9, 0.9, 0.9, 0.5, 0.5]
    tu = [f"TU_{i}" for i in range(n)]
    tu[0:5] = ["TU_A"] * 5
    c = Contig("chr1", d, p, tu=tu)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", [3, 4]))})
    row = pick(res, "bgc_1", "o_minus_e")
    assert row["effect"] == pytest.approx(1.0)
    assert row["ci_low"] == pytest.approx(1.0) and row["ci_high"] == pytest.approx(1.0)


def test_o_minus_e_ci_inside_attainable_range() -> None:
    n = 40
    d, p = _alt(n, 0.3)
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(10, 16)))}, rd.config(200))
    row = pick(res, "bgc_1", "o_minus_e")
    # six single-gene TUs with r in {-0.3, +0.7}: any resample sum lies in [-1.8, 4.2]
    assert -1.8 - 1e-9 <= row["ci_low"] <= row["ci_high"] <= 4.2 + 1e-9
    assert row["ci_low"] < row["ci_high"]


def test_hl_ci_with_constant_sides_is_degenerate() -> None:
    """HL CI resamples inside and outside TUs as two strata: constant 20 inside and 10 outside
    gives the CI [10, 10]."""
    n = 30
    a = [10.0] * n
    a[5:10] = [20.0] * 5
    c = Contig("chr1", [1] * n, [0.5] * n, a=a)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(5, 10)))})
    row = pick(res, "bgc_1", "hl_shift")
    assert row["effect"] == pytest.approx(10.0)
    assert row["ci_low"] == pytest.approx(10.0) and row["ci_high"] == pytest.approx(10.0)


# --------------------------------------------------------------------------- BH, RNG


def _many_regions() -> tuple[list[Contig], dict[str, tuple[str, list[str]]]]:
    rng = np.random.default_rng(8)
    n = 200
    d = (rng.random(n) < 0.5).astype(int).tolist()
    p = rng.uniform(0.2, 0.8, n).tolist()
    a = rng.normal(15, 2, n).tolist()
    c = Contig("chr1", d, p, a=a, tu=[f"TU_{i // 2}" for i in range(n)])
    genes, _ = gene_table([c])
    members: dict[str, tuple[str, list[str]]] = {}
    for j, typ in enumerate(["prophage", "prophage", "prophage", "bgc", "bgc", "crispr"]):
        members[f"{typ}_{j}"] = (typ, tags_of(genes, "chr1", range(30 * j + 2, 30 * j + 14)))
    members["amr_1"] = ("amr", tags_of(genes, "chr1", [1, 50, 99, 150, 197]))
    return [c], members


def test_q_values_bh_within_type_and_stat() -> None:
    """q = BH within each (region_type, stat) family (lead decision 5)."""
    from giap.stats.multitest import bh

    contigs, members = _many_regions()
    res = _run(contigs, members)
    for (_, _), grp in res.groupby(["region_type", "stat"]):
        np.testing.assert_allclose(
            grp["q_value"].to_numpy(float), bh(grp["p_value"].to_numpy(float)), equal_nan=True
        )


def test_deterministic_and_independent_of_region_row_order() -> None:
    """Regions are processed in sorted region_id order: same seed -> identical results whatever the
    row order of C7/C8."""
    contigs, members = _many_regions()
    genes, regions, rg, replicons = _tables(contigs, members)
    cfg = rd.config(bootstrap=40)
    R = _regions()
    a = R.test_all_regions(genes, regions, rg, replicons, cfg, np.random.default_rng(5))
    b = R.test_all_regions(
        genes, regions.iloc[::-1].reset_index(drop=True),
        rg.sample(frac=1.0, random_state=2).reset_index(drop=True), replicons, cfg,
        np.random.default_rng(5),
    )  # fmt: skip
    pd.testing.assert_frame_equal(
        _sorted(a)[list(rd.TABLE_COLUMNS)], _sorted(b)[list(rd.TABLE_COLUMNS)]
    )


def test_gene_row_order_does_not_matter() -> None:
    contigs, members = _many_regions()
    genes, regions, rg, replicons = _tables(contigs, members)
    cfg = rd.config(bootstrap=40)
    R = _regions()
    a = R.test_all_regions(genes, regions, rg, replicons, cfg, np.random.default_rng(5))
    shuffled = genes.sample(frac=1.0, random_state=9).reset_index(drop=True)
    b = R.test_all_regions(shuffled, regions, rg, replicons, cfg, np.random.default_rng(5))
    cols = ["region_id", "stat", "observed", "expected", "effect", "p_value", "n_perm"]
    pd.testing.assert_frame_equal(_sorted(a)[cols], _sorted(b)[cols])


# --------------------------------------------------------------------------- real producer frames


def test_runs_on_frames_loaded_through_m4_06(tmp_path: Path) -> None:
    """giapsim M4-07 layout (plasmid with its own AMR set, origin crossing) + M4-08 truth:
    written as the C1-C8 contract and re-read with ``load_m4_inputs`` gives the same results as
    the in-memory frames (list cells, nullable dtypes)."""
    rng = np.random.default_rng(17)
    lay, expr = rd.simulate(rng, rd.calibration_layout_config(300), {"ar1_rho": 0.5})
    mem = rd.tables_from(lay, expr, "true")
    loaded, issues = rd.loaded_tables(lay, expr, tmp_path / "work")
    # only C1-C8 are written: the M3 mapping files are reported missing, nothing else may fail
    bad = [i for i in issues if i.severity == "error" and not i.code.endswith("missing_file")]
    assert not bad, bad
    cfg = rd.config(bootstrap=20)
    a = rd.run_region_tests(mem, cfg, np.random.default_rng(3))
    b = rd.run_region_tests(loaded, cfg, np.random.default_rng(3))
    cols = ["region_id", "contig", "stat", "effect", "p_value", "q_value", "n_perm", "n_eligible"]
    pd.testing.assert_frame_equal(_sorted(a)[cols], _sorted(b)[cols])
    assert {"prophage", "bgc", "crispr", "amr"} == set(a["region_type"])
    amr_contigs = set(a.loc[a["region_type"] == "amr", "contig"])
    assert amr_contigs == {"chr1", "p1"}
    assert len(a) == 2 * 5  # 4 chromosome regions + the plasmid AMR set, two statistics each
