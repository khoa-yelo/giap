"""M4-17 region tests: spec examples 1-10 and the public interface.

docs/specs/M4-17-region-tests.md.

GIAP_SPEC.md §7.3.5(b), §8.1, Appendix B.4. Inputs are hand-built §9.3 gene tables
(tests/helpers/region_tests_data.py); expected values come from the spec table, cross-checked
against the naive oracles there. ``giap.context.regions`` is imported as a module so pytest never
collects its ``test_region`` / ``test_all_regions`` functions.
"""

from __future__ import annotations

import dataclasses
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


def _run(
    contigs: list[Contig], members: dict[str, tuple[str, list[str]]], cfg: Any = None, seed: int = 1
) -> pd.DataFrame:
    genes, replicons = gene_table(contigs)
    genes, regions, rg = region_frames(genes, members)
    cfg = cfg if cfg is not None else rd.config()
    return _regions().test_all_regions(
        genes, regions, rg, replicons, cfg, np.random.default_rng(seed)
    )


def _alternating(n: int, p: float = 0.5) -> tuple[list[int], list[float]]:
    return [i % 2 for i in range(n)], [p] * n


# --------------------------------------------------------------------------- interface


def test_region_test_dataclass_fields() -> None:
    R = _regions()
    assert dataclasses.is_dataclass(R.RegionTest)
    assert tuple(f.name for f in dataclasses.fields(R.RegionTest)) == rd.RESULT_FIELDS
    params = R.RegionTest.__dataclass_params__
    assert params.frozen
    assert hasattr(R.RegionTest, "__slots__")


def test_table_columns_method_and_null_model() -> None:
    d, p = _alternating(30)
    genes, _ = gene_table([Contig("chr1", d, p)])
    res = _run(
        [Contig("chr1", d, p)], {"prophage_1": ("prophage", tags_of(genes, "chr1", range(5, 12)))}
    )
    assert set(rd.TABLE_COLUMNS) <= set(res.columns)
    assert sorted(res["stat"]) == ["hl_shift", "o_minus_e"]
    for _, row in res.iterrows():
        assert row["method"] == rd.METHOD[row["stat"]]
        assert row["null_model"] == "circular_shift"
        assert row["region_id"] == "prophage_1" and row["region_type"] == "prophage"
        assert row["contig"] == "chr1"
    oe = pick(res, "prophage_1", "o_minus_e")
    assert oe["n_genes"] == 7 and oe["n_eligible"] == 7
    assert pick(res, "prophage_1", "hl_shift")["expected"] is None or rd.is_na(
        pick(res, "prophage_1", "hl_shift")["expected"]
    )


def test_replicon_index_eligible_genes_in_genome_order() -> None:
    """Index = eligible genes of the contig in genome order, 0..N-1 (lead decision 5: partial
    genes are not indexed); columns idx, locus_tag, tu_id, d, p, a; row order of the input
    does not matter."""
    n = 8
    part = [True] + [False] * (n - 2) + [True]
    d = [1, 1, 0, 1, 0, 0, 1, 1]
    p = [0.1 * (i + 1) for i in range(n)]
    genes, _ = gene_table(
        [Contig("ctg1", d, p, circular=False, partial=part), Contig("p1", [1, 0], [0.5, 0.5])]
    )
    idx = _regions().replicon_index(genes.sample(frac=1.0, random_state=3), "ctg1")
    assert {"idx", "locus_tag", "tu_id", "d", "p", "a"} <= set(idx.columns)
    assert list(idx["idx"]) == list(range(n - 2))
    assert list(idx["locus_tag"]) == tags_of(genes, "ctg1", range(1, n - 1))
    assert [int(x) for x in idx["d"]] == d[1:-1]
    np.testing.assert_allclose(idx["p"].to_numpy(float), p[1:-1])
    a = idx["a"].to_numpy(float)
    assert np.isnan(a[np.array(d[1:-1]) == 0]).all()
    np.testing.assert_allclose(a[np.array(d[1:-1]) == 1], [11.0, 13.0, 16.0])  # 10 + gene index


def test_replicon_index_d_and_a_values() -> None:
    d = [1, 0, 1]
    genes, _ = gene_table([Contig("chr1", d, [0.2, 0.3, 0.4], a=[7.0, 99.0, 9.0])])
    idx = _regions().replicon_index(genes, "chr1")
    assert [int(x) for x in idx["d"]] == d
    a = idx["a"].to_numpy(float)
    assert a[0] == 7.0 and math.isnan(a[1]) and a[2] == 9.0


def test_o_minus_e_null_small() -> None:
    """S_k = sum of r over (I + k) mod N."""
    r = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    out = _regions().o_minus_e_null(r, np.array([0, 1]), np.array([1, 2, 3, 4]))
    np.testing.assert_allclose(out, [5.0, 7.0, 9.0, 6.0])


def test_hl_null_small() -> None:
    """HL(inside detected & observed, outside detected & observed) at every shift; NaN when a side
    is empty."""
    a = np.array([10.0, 12.0, 8.0, 9.0, np.nan, 11.0])
    det = np.array([True, True, True, True, False, True])
    shifts = np.array([1, 2, 3, 4, 5])
    out = _regions().hl_null(a, det, np.array([0, 1]), shifts)
    np.testing.assert_allclose(out, rd.naive_hl_null(a, det, [0, 1], shifts), equal_nan=True)
    empty_inside = _regions().hl_null(a, det, np.array([2]), np.array([2]))  # inside = {4}
    assert empty_inside.shape == (1,) and math.isnan(empty_inside[0])


# --------------------------------------------------------------------------- examples


def test_example_1_o_minus_e() -> None:
    """4 eligible genes, d = [1, 0, 1, 1], p = [0.9, 0.2, 0.5, 0.4] -> observed 3, expected 2.0,
    effect 1.0."""
    d, p = _alternating(20)
    d[5:9], p[5:9] = [1, 0, 1, 1], [0.9, 0.2, 0.5, 0.4]
    genes, _ = gene_table([Contig("chr1", d, p)])
    res = _run([Contig("chr1", d, p)], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(5, 9)))})
    row = pick(res, "bgc_1", "o_minus_e")
    assert row["observed"] == pytest.approx(3.0)
    assert row["expected"] == pytest.approx(2.0)
    assert row["effect"] == pytest.approx(1.0)
    assert row["n_eligible"] == 4 and row["n_genes"] == 4
    assert "small_region" not in flag_set(row["flags"])


def test_example_1_through_test_region() -> None:
    d, p = _alternating(20)
    d[5:9], p[5:9] = [1, 0, 1, 1], [0.9, 0.2, 0.5, 0.4]
    genes, _ = gene_table([Contig("chr1", d, p)])
    genes, _, _ = region_frames(genes, {"bgc_1": ("bgc", tags_of(genes, "chr1", range(5, 9)))})
    rows = _regions().test_region(
        genes, tags_of(genes, "chr1", range(5, 9)), "bgc", "chr1", rd.config(),
        np.random.default_rng(0),
    )  # fmt: skip
    by = {r.stat: r for r in rows}
    assert set(by) == {"o_minus_e", "hl_shift"}
    assert by["o_minus_e"].observed == pytest.approx(3.0)
    assert by["o_minus_e"].expected == pytest.approx(2.0)
    assert by["o_minus_e"].effect == pytest.approx(1.0)
    assert by["o_minus_e"].n_perm == 19
    assert by["hl_shift"].expected is None or math.isnan(by["hl_shift"].expected)


def test_example_2_hodges_lehmann() -> None:
    """Region detected log2 [10, 12], other detected [8, 9, 11] -> HL effect 1.5."""
    d = [1, 1, 0, 1, 1, 1, 0, 0]
    a = [10.0, 12.0, 0.0, 8.0, 9.0, 11.0, 0.0, 0.0]
    genes, _ = gene_table([Contig("chr1", d, [0.5] * 8, a=a)])
    res = _run(
        [Contig("chr1", d, [0.5] * 8, a=a)],
        {"crispr_1": ("crispr", tags_of(genes, "chr1", [0, 1, 2]))},
    )
    row = pick(res, "crispr_1", "hl_shift")
    assert row["effect"] == pytest.approx(1.5)
    assert row["effect"] == pytest.approx(rd.naive_hl([10, 12], [8, 9, 11]))


def test_example_3_appendix_b4() -> None:
    """B.4: circular replicon of 100 eligible genes, observed more extreme than all 99 shifts ->
    n_perm 99, p_value 0.01."""
    c = rd.hot_contig("chr1", 0)
    genes, _ = gene_table([c])
    res = _run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", [0, 1, 2]))})
    row = pick(res, "prophage_1", "o_minus_e")
    assert row["n_perm"] == 99
    assert row["p_value"] == pytest.approx(0.01)


@pytest.mark.parametrize("partial_ends", [False, True])
def test_example_4_linear_shifts(partial_ends: bool) -> None:
    """Linear contig, N = 10 eligible, I = {2, 3, 4} -> n_perm 7 (shifts 1-5, 8, 9). With partial
    genes at both ends (12 genes), the index still has N = 10 (lead decision 5)."""
    n = 12 if partial_ends else 10
    part = [True] + [False] * 10 + [True] if partial_ends else None
    off = 1 if partial_ends else 0
    a = [float(i) for i in range(n)]
    c = Contig("ctg1", [1] * n, [0.5] * n, a=a, circular=False, partial=part)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "ctg1", [2 + off, 3 + off, 4 + off]))})
    assert rd.expected_shifts(10, [2, 3, 4], circular=False) == [1, 2, 3, 4, 5, 8, 9]
    assert pick(res, "bgc_1", "o_minus_e")["n_perm"] == 7
    assert pick(res, "bgc_1", "hl_shift")["n_perm"] == 7


def test_example_5_linear_no_valid_shift() -> None:
    """Linear contig, N = 5, I = {0, 4} -> K = 0: p_value NA, n_perm 0 (both statistics)."""
    c = Contig("ctg1", [1, 0, 1, 1, 1], [0.5] * 5, circular=False)
    genes, _ = gene_table([c])
    res = _run([c], {"amr_1": ("amr", tags_of(genes, "ctg1", [0, 4]))})
    for stat in rd.STATS:
        row = pick(res, "amr_1", stat)
        assert row["n_perm"] == 0
        assert rd.is_na(row["p_value"])
        assert rd.is_na(row["q_value"])
    assert pick(res, "amr_1", "o_minus_e")["effect"] == pytest.approx(1.0)


@pytest.fixture
def shift_spy(monkeypatch: pytest.MonkeyPatch) -> list[tuple[np.ndarray, bool]]:
    """Records every ``choose_shifts`` result the region tests use (spec "Shifts")."""
    import giap.stats.permutation as perm

    R = _regions()
    orig = perm.choose_shifts
    calls: list[tuple[np.ndarray, bool]] = []

    def spy(*args: Any, **kw: Any) -> tuple[np.ndarray, bool]:
        out = orig(*args, **kw)
        calls.append((np.asarray(out[0]).copy(), bool(out[1])))
        return out

    monkeypatch.setattr(perm, "choose_shifts", spy)
    if hasattr(R, "choose_shifts"):
        monkeypatch.setattr(R, "choose_shifts", spy)
    return calls


def test_example_6_random_shifts_above_threshold(shift_spy: list[tuple[np.ndarray, bool]]) -> None:
    """N = 25,000 > exhaustive_shift_max_genes = 20,000, stats.permutations = 1000 -> n_perm
    1000, shifts drawn at random without replacement from 1..N-1."""
    n = 25_000
    rng = np.random.default_rng(6)
    d = (rng.random(n) < 0.002).astype(int).tolist()
    d[100:110] = [1] * 10
    p = rng.uniform(0.0, 0.01, n).tolist()
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    cfg = rd.config(bootstrap=10, permutations=1000)
    assert cfg.context.exhaustive_shift_max_genes == 20_000
    res = _run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", range(100, 130)))}, cfg)
    assert pick(res, "prophage_1", "o_minus_e")["n_perm"] == 1000
    assert pick(res, "prophage_1", "hl_shift")["n_perm"] <= 1000
    assert shift_spy, "choose_shifts was not used"
    for shifts, exhaustive in shift_spy:
        assert not exhaustive
        assert shifts.size == 1000 and np.unique(shifts).size == 1000
        assert shifts.min() >= 1 and shifts.max() <= n - 1


@pytest.mark.parametrize(("n", "K"), [(50, 49), (51, 20)])
def test_exhaustive_threshold_is_inclusive(n: int, K: int) -> None:
    """Exhaustive while N <= exhaustive_shift_max_genes, otherwise stats.permutations shifts."""
    d, p = _alternating(n, 0.4)
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    cfg = rd.config(permutations=20, exhaustive=50)
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(10, 15)))}, cfg)
    assert pick(res, "bgc_1", "o_minus_e")["n_perm"] == K


def test_example_7_small_region_still_tested() -> None:
    """2 eligible genes -> flag small_region; the test is still reported (K > 0)."""
    d, p = _alternating(40)
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    res = _run([c], {"crispr_1": ("crispr", tags_of(genes, "chr1", [3, 4]))})
    for stat in rd.STATS:
        row = pick(res, "crispr_1", stat)
        assert "small_region" in flag_set(row["flags"])
        assert row["n_eligible"] == 2
        assert row["n_perm"] > 0
    assert not rd.is_na(pick(res, "crispr_1", "o_minus_e")["p_value"])


def test_three_eligible_genes_not_small() -> None:
    d, p = _alternating(40)
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    res = _run([c], {"crispr_1": ("crispr", tags_of(genes, "chr1", [3, 4, 5]))})
    assert "small_region" not in flag_set(pick(res, "crispr_1", "o_minus_e")["flags"])


def test_example_8_bh_within_prophage_o_minus_e() -> None:
    """3 prophages with O - E p = [0.01, 0.04, 0.03] -> q = [0.03, 0.04, 0.04] (BH within
    prophage x o_minus_e). The p-values are exact: ``hot_contig`` plants 0, 3, 2 tying shifts."""
    contigs = [rd.hot_contig("chr1", 0), rd.hot_contig("p1", 3), rd.hot_contig("p2", 2)]
    genes, _ = gene_table(contigs)
    members = {
        f"prophage_{i + 1}": ("prophage", tags_of(genes, c.name, [0, 1, 2]))
        for i, c in enumerate(contigs)
    }
    res = _run(contigs, members)
    p = [pick(res, f"prophage_{i}", "o_minus_e")["p_value"] for i in (1, 2, 3)]
    q = [pick(res, f"prophage_{i}", "o_minus_e")["q_value"] for i in (1, 2, 3)]
    np.testing.assert_allclose(p, [0.01, 0.04, 0.03])
    np.testing.assert_allclose(q, [0.03, 0.04, 0.04])


def test_example_9_no_detected_gene() -> None:
    """Region with no detected gene -> HL row effect NA, p_value NA; O - E row normal."""
    d, p = _alternating(30)
    d[10:15] = [0] * 5
    c = Contig("chr1", d, p)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(10, 15)))})
    hl = pick(res, "bgc_1", "hl_shift")
    assert rd.is_na(hl["effect"]) and rd.is_na(hl["p_value"]) and rd.is_na(hl["q_value"])
    oe = pick(res, "bgc_1", "o_minus_e")
    assert oe["effect"] == pytest.approx(-2.5)
    assert oe["n_perm"] == 29
    assert 0.0 < oe["p_value"] <= 1.0


def test_example_10_planted_detection_shift() -> None:
    """d - p = +0.5 for 30 region genes on a 1,000-gene circular chromosome, iid elsewhere ->
    p_value <= 0.01."""
    rng = np.random.default_rng(10)
    n = 1000
    d = (rng.random(n) < 0.5).astype(int).tolist()
    d[400:430] = [1] * 30
    c = Contig("chr1", d, [0.5] * n)
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(400, 430)))})
    row = pick(res, "bgc_1", "o_minus_e")
    assert row["effect"] == pytest.approx(15.0)
    assert row["n_perm"] == 999
    assert row["p_value"] <= 0.01
