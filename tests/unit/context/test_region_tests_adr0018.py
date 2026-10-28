"""M4-17 region tests after the statistician's review: ADR-0018 and lead decisions 6 and 9.

docs/adr/ADR-0018-region-tests-confounding.md (accepted-by-khoa 2026-10-15),
docs/specs/M4-17-region-tests.md "Lead decisions after the statistician's review":

- decision 1 (ADR): HL is computed on log2 abundance residualised on the M4-13 detectability
  covariates (log length, n_obs, GRAVY) among the replicon's detected genes, before the
  circular-shift test. Test-author reading, pinned by the oracle ``rd.naive_residualise``: OLS with
  an intercept on log(length_aa), log1p(n_observable_peptides), gravy (the M4-13
  ``transform_covariates`` transforms; a non-finite transform -> the column median over the fitted
  genes), fitted once per replicon (residuals do not depend on the region, so the null is the same
  shift orbit). ``effect``, CI and null are all on the residual scale (one statistic).
- decision 2 (ADR): the O - E test statistic is the studentised sum
  T(I) = sum(d - p) / sqrt(sum p (1 - p)), recomputed at every shift; ``observed``, ``expected``
  and ``effect`` (= O - E) are unchanged.
- decision 3 (ADR): BH within (region_type, stat) is pinned in test_region_tests_edges.py.
- lead decision 9: the surviving mutants of the review (``split_contigs`` flag dropped, HL
  bootstrap without strata, index sorted by locus_tag before start) are killed here.
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
    from gippa.context import regions as R

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


# --------------------------------------------------------------------------- studentised O - E


def _hetero_contig(seed: int = 4, n: int = 40) -> tuple[Contig, list[int]]:
    """Replicon p = 0.9 except a mid-range region (p = 0.5): the attack-2 anti-conservative
    direction in miniature (statistician, limitation (i))."""
    rng = np.random.default_rng(seed)
    p = [0.9] * n
    idx = list(range(12, 20))
    for i in idx:
        p[i] = 0.5
    d = (rng.random(n) < np.array(p)).astype(int).tolist()
    return Contig("chr1", d, p), idx


def test_o_minus_e_p_uses_studentised_statistic() -> None:
    """p = two-sided centred perm p of T_k = sum(d - p) / sqrt(sum p (1 - p)) over the shifts;
    observed / expected / effect stay sum d / sum p / O - E."""
    c, idx = _hetero_contig()
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))})
    row = pick(res, "bgc_1", "o_minus_e")
    d, p = np.array(c.d, float), np.array(c.p, float)
    r, v = d - p, p * (1 - p)
    shifts = list(range(1, len(d)))
    want = rd.naive_perm_p(rd.studentised(r, v, idx), rd.naive_studentised_null(r, v, idx, shifts))
    raw = rd.naive_perm_p(float(r[idx].sum()), rd.naive_o_minus_e_null(r, idx, shifts))
    assert want != pytest.approx(raw), "fixture does not separate studentised from plain O - E"
    assert row["p_value"] == pytest.approx(want)
    assert row["n_perm"] == len(shifts)
    assert row["observed"] == pytest.approx(d[idx].sum())
    assert row["expected"] == pytest.approx(p[idx].sum())
    assert row["effect"] == pytest.approx(r[idx].sum())


def test_studentised_equals_plain_when_p_is_constant() -> None:
    """With one p on the whole replicon, sum p (1 - p) is the same at every shift: the
    studentised and plain statistics give the same p (spec example 10 keeps its value)."""
    rng = np.random.default_rng(11)
    n = 60
    d = (rng.random(n) < 0.4).astype(int).tolist()
    c = Contig("chr1", d, [0.4] * n)
    genes, _ = gene_table([c])
    idx = list(range(20, 30))
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))})
    r = np.array(d, float) - 0.4
    raw = rd.naive_perm_p(float(r[idx].sum()), rd.naive_o_minus_e_null(r, idx, range(1, n)))
    assert pick(res, "bgc_1", "o_minus_e")["p_value"] == pytest.approx(raw)


def test_o_minus_e_ci_stays_on_o_minus_e_scale() -> None:
    """The reported effect is O - E (ADR-0018 decision 2), so its TU-block CI is too: with every
    region gene in its own TU and r in {-0.5, +0.5}, any resample sum lies in [-4, 4]."""
    c, idx = _hetero_contig()
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}, rd.config(200))
    row = pick(res, "bgc_1", "o_minus_e")
    assert -4.0 - 1e-9 <= row["ci_low"] <= row["effect"] <= row["ci_high"] <= 4.0 + 1e-9
    assert row["ci_low"] < row["ci_high"]


# --------------------------------------------------------------------------- residualised HL


def _cov_contig(seed: int, n: int = 40) -> Contig:
    rng = np.random.default_rng(seed)
    d = (rng.random(n) < 0.7).astype(int).tolist()
    length = rng.integers(40, 800, n)
    n_obs = np.maximum(1, length // 25 + rng.integers(-3, 4, n))
    gravy = rng.normal(-0.3, 0.4, n)
    a = 18 + 1.5 * np.log(length) - 0.8 * gravy + rng.normal(0, 1.0, n)
    return Contig(
        "chr1", d, rng.uniform(0.3, 0.9, n).tolist(), a=a.tolist(),
        length_aa=length.tolist(), n_obs=n_obs.tolist(), gravy=gravy.tolist(),
    )  # fmt: skip


def _hl_oracle(c: Contig, idx: list[int]) -> tuple[float, float, int]:
    det = np.array(c.d, bool)
    a = np.where(det, np.array(c.a, float), np.nan)
    ar = rd.naive_residualise(a, det, c.length_aa, c.n_obs, c.gravy)  # type: ignore[arg-type]
    ok = det & np.isfinite(ar)
    x = [ar[i] for i in idx if ok[i]]
    y = [ar[i] for i in range(len(ar)) if ok[i] and i not in idx]
    shifts = list(range(1, len(ar)))
    # ADR-0018 amendment 2: p from T = HL / sqrt(V), V label-free (lag autocovariance)
    null = rd.naive_studentised_hl_null(ar, ok, idx, shifts)
    t_obs = float(rd.naive_studentised_hl_null(ar, ok, idx, [0])[0])
    eff = rd.naive_hl(x, y)
    return eff, rd.naive_perm_p(t_obs, null), int(np.count_nonzero(~np.isnan(null)))


def test_hl_on_residualised_abundance_matches_oracle() -> None:
    c = _cov_contig(21)
    genes, _ = gene_table([c])
    idx = list(range(8, 18))
    res = _run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", idx))})
    row = pick(res, "prophage_1", "hl_shift")
    eff, p, k = _hl_oracle(c, idx)
    det = np.array(c.d, bool)
    a = np.where(det, np.array(c.a, float), np.nan)
    raw = rd.naive_hl([a[i] for i in idx if det[i]],
                      [a[i] for i in range(len(a)) if det[i] and i not in idx])  # fmt: skip
    assert abs(eff - raw) > 0.05, "fixture does not separate residualised from raw HL"
    assert row["effect"] == pytest.approx(eff, abs=1e-9)
    assert row["p_value"] == pytest.approx(p)
    assert row["n_perm"] == k
    assert row["ci_low"] <= row["effect"] + 1e-9 and row["effect"] <= row["ci_high"] + 1e-9


def test_abundance_explained_by_covariates_gives_no_shift() -> None:
    """Attack 2 in miniature: the region's genes are short and log2 abundance is an exact linear
    function of the detectability covariates. Raw HL sees a large shift (the region is the most
    extreme position, p = 1/N); after residualising nothing is left: effect 0, p ~ 1."""
    n = 60
    idx = list(range(20, 30))
    length = [30 if i in idx else 300 + 7 * (i % 5) for i in range(n)]
    n_obs = [max(1, L // 20) for L in length]
    gravy = [-0.4 + 0.1 * (i % 7) for i in range(n)]
    a = [16 + 2.0 * math.log(L) + 0.7 * math.log1p(k) - 0.5 * g
         for L, k, g in zip(length, n_obs, gravy, strict=True)]  # fmt: skip
    c = Contig("chr1", [1] * n, [0.6] * n, a=a, length_aa=length, n_obs=n_obs, gravy=gravy)
    genes, _ = gene_table([c])
    raw_null = rd.naive_hl_null(np.array(a), np.ones(n, bool), idx, range(1, n))
    raw_eff = rd.naive_hl([a[i] for i in idx], [a[i] for i in range(n) if i not in idx])
    assert rd.naive_perm_p(raw_eff, raw_null) <= 0.05  # the confounded raw statistic rejects
    res = _run([c], {"prophage_1": ("prophage", tags_of(genes, "chr1", idx))})
    row = pick(res, "prophage_1", "hl_shift")
    assert abs(row["effect"]) < 1e-6
    assert row["p_value"] > 0.5


def test_constant_covariates_leave_hl_as_raw() -> None:
    """Constant covariates only centre the abundances: HL effect equals the raw HL (spec example 2
    keeps HL 1.5) and p is the studentised p of the raw values (ADR-0018 amendments 1, 2)."""
    d = [1, 1, 0, 1, 1, 1, 0, 0]
    a = [10.0, 12.0, 0.0, 8.0, 9.0, 11.0, 0.0, 0.0]
    c = Contig("chr1", d, [0.5] * 8, a=a)
    genes, _ = gene_table([c])
    res = _run([c], {"crispr_1": ("crispr", tags_of(genes, "chr1", [0, 1, 2]))})
    row = pick(res, "crispr_1", "hl_shift")
    assert row["effect"] == pytest.approx(1.5)
    det = np.array(d, bool)
    av = np.where(det, np.array(a), np.nan)
    null = rd.naive_studentised_hl_null(av, det, [0, 1, 2], list(range(1, 8)))
    t_obs = float(rd.naive_studentised_hl_null(av, det, [0, 1, 2], [0])[0])
    assert row["p_value"] == pytest.approx(rd.naive_perm_p(t_obs, null))


def test_residualisation_uses_detected_genes_only() -> None:
    """Covariates of undetected genes do not enter the fit: changing them leaves HL unchanged."""
    c = _cov_contig(23)
    idx = list(range(5, 14))
    genes, _ = gene_table([c])
    members = {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}
    und = [i for i, x in enumerate(c.d) if not x]
    assert und
    L, G = list(c.length_aa or []), list(c.gravy or [])
    for i in und:
        L[i], G[i] = 5000, 3.0
    c2 = Contig(c.name, c.d, c.p, a=c.a, length_aa=L, n_obs=c.n_obs, gravy=G)
    a = pick(_run([c], members), "bgc_1", "hl_shift")
    b = pick(_run([c2], members), "bgc_1", "hl_shift")
    for col in ("effect", "p_value", "ci_low", "ci_high"):
        assert a[col] == pytest.approx(b[col], abs=1e-9)


def test_residualisation_is_per_replicon() -> None:
    """The fit uses the tested replicon's genes only: a plasmid with other covariates and
    abundances does not change the chromosome's HL row."""
    c = _cov_contig(25)
    idx = list(range(10, 20))
    rng = np.random.default_rng(2)
    m = 30
    plasmid = Contig("p1", [1] * m, [0.5] * m, a=rng.normal(25, 3, m).tolist(),
                     length_aa=rng.integers(40, 900, m).tolist(),
                     n_obs=rng.integers(1, 40, m).tolist(),
                     gravy=rng.normal(0.5, 0.5, m).tolist())  # fmt: skip
    genes, _ = gene_table([c])
    members = {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}
    alone = pick(_run([c], members), "bgc_1", "hl_shift")
    both = pick(_run([c, plasmid], members), "bgc_1", "hl_shift")
    eff, p, _ = _hl_oracle(c, idx)
    assert both["effect"] == pytest.approx(eff, abs=1e-9)
    assert both["p_value"] == pytest.approx(p)
    assert both["effect"] == pytest.approx(alone["effect"], abs=1e-12)


def test_missing_covariate_imputed_by_median() -> None:
    """A detected gene whose gravy is NaN (M4-10: protein with nothing left after cleaning) is
    kept; its transform takes the column median of the fitted genes (M4-13 rule)."""
    c = _cov_contig(27)
    G = list(c.gravy or [])
    j = next(i for i, x in enumerate(c.d) if x and i not in range(6, 15))
    G[j] = math.nan
    c2 = Contig(c.name, c.d, c.p, a=c.a, length_aa=c.length_aa, n_obs=c.n_obs, gravy=G)
    genes, _ = gene_table([c2])
    idx = list(range(6, 15))
    row = pick(_run([c2], {"bgc_1": ("bgc", tags_of(genes, "chr1", idx))}), "bgc_1", "hl_shift")
    eff, p, k = _hl_oracle(c2, idx)
    assert not math.isnan(eff)
    assert row["effect"] == pytest.approx(eff, abs=1e-9)
    assert row["p_value"] == pytest.approx(p)
    assert row["n_perm"] == k


# --------------------------------------------------------------------------- surviving mutants


def test_split_contigs_flag_on_every_row_of_a_split_region_only() -> None:
    """Lead decision 2: a region whose members lie on two contigs is flagged ``split_contigs`` on
    each contig's rows; a region on one contig is not."""
    n1, n2 = 30, 25
    contigs = [
        Contig("chr1", [i % 2 for i in range(n1)], [0.5] * n1),
        Contig("p1", [i % 2 for i in range(n2)], [0.5] * n2),
    ]
    genes, _ = gene_table(contigs)
    tags = tags_of(genes, "chr1", [2, 9, 17, 23]) + tags_of(genes, "p1", [1, 5, 8])
    res = _run(contigs, {"amr_1": ("amr", tags),
                         "bgc_1": ("bgc", tags_of(genes, "chr1", range(10, 15)))})  # fmt: skip
    amr = res[res["region_id"] == "amr_1"]
    assert len(amr) == 4
    assert all("split_contigs" in flag_set(f) for f in amr["flags"])
    bgc = res[res["region_id"] == "bgc_1"]
    assert not any("split_contigs" in flag_set(f) for f in bgc["flags"])


def test_hl_bootstrap_resamples_inside_and_outside_tus_as_strata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Spec "CI": HL resamples inside and outside TUs as two strata. The HL bootstrap call gets a
    stratum column that splits the TUs into the region's 3 TUs and the 12 others."""
    R = _regions()
    calls: list[tuple[pd.DataFrame, str, str | None]] = []
    orig = R.tu_block_bootstrap

    def spy(df: pd.DataFrame, tu_col: str, *args: Any, **kw: Any) -> Any:
        strata = kw.get("strata_col", args[4] if len(args) > 4 else None)
        calls.append((df.copy(), tu_col, strata))
        return orig(df, tu_col, *args, **kw)

    monkeypatch.setattr(R, "tu_block_bootstrap", spy)
    n = 30
    c = Contig("chr1", [1] * n, [0.5] * n, a=[10.0 + (i % 4) for i in range(n)],
               tu=[f"TU_{i // 2}" for i in range(n)])  # fmt: skip
    genes, _ = gene_table([c])
    res = _run([c], {"bgc_1": ("bgc", tags_of(genes, "chr1", range(10, 16)))})
    assert not rd.is_na(pick(res, "bgc_1", "hl_shift")["ci_low"])
    hl_calls = [x for x in calls if len(x[0]) == n]
    assert len(hl_calls) == 1, [len(x[0]) for x in calls]
    df, tu_col, strata = hl_calls[0]
    assert strata is not None, "HL bootstrap must stratify inside / outside TUs"
    per_tu = df.groupby(tu_col)[strata].nunique()
    assert (per_tu == 1).all()
    sizes = sorted(df.groupby(strata)[tu_col].nunique().tolist())
    assert sizes == [3, 12]


def test_index_sorts_by_start_before_locus_tag() -> None:
    """Genome order is by ``start`` (locus tags need not sort in genome order, e.g. user
    annotations): tags reversed against position must not reorder the index."""
    n = 12
    d = [1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1]
    c = Contig("chr1", d, [0.1 * (1 + i % 9) for i in range(n)])
    genes, _ = gene_table([c])
    genes["locus_tag"] = [f"Z{n - i:03d}" for i in range(n)]  # descending in genome order
    ix = _regions().replicon_index(genes.sample(frac=1.0, random_state=4), "chr1")
    assert list(ix["locus_tag"]) == [f"Z{n - i:03d}" for i in range(n)]
    assert [int(x) for x in ix["d"]] == d


def test_linear_shift_count_follows_positions_not_tags() -> None:
    """Spec example 4 (linear, N = 10, I = {2, 3, 4} -> n_perm 7) with scrambled locus tags."""
    n = 10
    c = Contig("ctg1", [1] * n, [0.5] * n, a=[float(i) for i in range(n)], circular=False)
    genes, replicons = gene_table([c])
    perm = [7, 2, 9, 0, 5, 3, 8, 1, 6, 4]
    genes["locus_tag"] = [f"T{perm[i]}" for i in range(n)]
    genes, regions, rg = region_frames(genes, {"bgc_1": ("bgc", ["T9", "T0", "T5"])})
    res = _regions().test_all_regions(
        genes, regions, rg, replicons, rd.config(), np.random.default_rng(0)
    )
    assert pick(res, "bgc_1", "o_minus_e")["n_perm"] == 7
