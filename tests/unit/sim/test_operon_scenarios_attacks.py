"""giapsim operon scenarios, round 2: lead decisions 10, 11, 13 and 14 of
docs/specs/M6-03-sim-operon-scenarios.md (the statistician's block on PR #222; GIAP_SPEC.md §8.4
``missed_operons``, §8.5, ADR-0018). Decision 12 (the attack bites) and the R = 20 share of
decision 10 are in tests/calibration/test_sim_operon_attack_calibration.py.

Interface assumed (test-author choices, reported to the lead):

* **decision 10**: ``operon_attack_detection`` returns an ``OperonConfig`` with
  ``detection = "abundance"``, a positive ``beta["log2_abundance"]`` and ``operon_t`` set
  (``OperonConfig`` no longer rejects ``operon_t`` with abundance detection). The detection logit
  of every gene is the merged ``giapsim.expression`` logit (``beta`` intercept, ``log_length``,
  ``log1p_n_obs``, ``gravy``, the region shift) plus ``beta["log2_abundance"]`` times the gene's
  *final* (post-walk, post-planting) ``abundance_log2_true``; ``p_detect_true`` is its sigmoid and
  ``detected_true`` is drawn from it after the walk. Any centring offset lives in the intercept;
* **decision 11**: two new ``LayoutConfig`` fields next to ``short_gaps_in`` / ``uniop_pure``:
  ``short_convergent_gaps: bool = False`` (the opt-in) and
  ``convergent_gap_bp: tuple[int, int] = (-20, 300)`` (inclusive integer range, read only when the
  opt-in is on). When on, the gap of every convergent (``+`` then ``-``) TU boundary is drawn
  uniformly from that range instead of the 50-300 bp inter-TU range; other gaps are unchanged.
  Off, nothing changes: the golden digests below are those of fd3f3a3. ``missed_operon_power``
  and ``operon_attack_detection`` turn it on; ``operon_null``, ``discordance_power`` and
  ``operon_attack_ar1`` leave it off. The wrap pair of a circular replicon is not checked;
  convergent pairs with both genes inside a ``short_gaps_in`` block keep their intra-TU gap;
* **decision 13**: ``operon_attack_ar1(n_genes, rho, *, n_planted=0, effect_sd=0.0)`` and
  ``operon_attack_detection(n_genes, short_genes_in, *, n_planted=0, effect_sd=0.0)``; with
  planting, the expression config is the unplanted one with ``n_missed = n_planted`` and
  ``missed_effect_sd = effect_sd`` (the layout config unchanged), so planting is the
  ``missed_operon_power`` mechanism and ``write_operon_truth`` writes the rows to
  ``truth/missed_operons.tsv``;
* **decision 14**: the null's config itself carries rho = 0 and no planted pair.
"""

from __future__ import annotations

import csv
import dataclasses
import hashlib
import math
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

N = 1800


def _sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


@cache
def sim(name: str, args: tuple[Any, ...], seed: int, kw: tuple[tuple[str, Any], ...] = ()) -> Any:
    sc = _sc()
    return sc.simulate_operons(getattr(sc, name)(*args, **dict(kw)), seed)


def _b(s: pd.Series) -> np.ndarray:
    return np.asarray(s.astype(bool).to_numpy(), dtype=bool)


def _block_genes(s: Any, kind: str) -> set[str]:
    ids = set(s.layout.regions.loc[s.layout.regions["region_type"] == kind, "region_id"])
    rg = s.layout.region_genes
    return set(rg.loc[rg["region_id"].isin(ids), "locus_tag"].astype(str))


def _not_wrap(s: Any, p: pd.DataFrame) -> np.ndarray:
    """Pairs whose right gene starts after the left one (drops the wrap pair of a ring)."""
    start = dict(zip(s.layout.genes["locus_tag"].astype(str), s.layout.genes["start"].astype(int),
                     strict=True))  # fmt: skip
    return np.array([start[str(b)] > start[str(a)] for a, b in zip(p["gene_a"], p["gene_b"],
                                                                   strict=True)])  # fmt: skip


def _cv_gaps(s: Any, exclude: set[str] = frozenset()) -> np.ndarray:  # type: ignore[assignment]
    p = s.pairs[s.pairs["class"] == "CV"]
    p = p[_not_wrap(s, p)]
    out = ~(p["gene_a"].astype(str).isin(exclude) & p["gene_b"].astype(str).isin(exclude))
    return np.asarray(p.loc[out.to_numpy(), "intergenic_bp"].to_numpy(), dtype=np.int64)


def _with_layout(cfgs: tuple[Any, Any], **kw: Any) -> tuple[Any, Any]:
    return dataclasses.replace(cfgs[0], **kw), cfgs[1]


# --------------------------------------------------------------- decision 11: short CV gaps


def test_short_convergent_gaps_defaults_off_with_the_minus20_300_range() -> None:  # decision 23
    from giapsim.layout import LayoutConfig

    cfg = LayoutConfig()
    assert cfg.short_convergent_gaps is False
    assert tuple(cfg.convergent_gap_bp) == (-20, 300)


def test_short_convergent_gaps_off_explicitly_is_the_default_layout() -> None:
    from giapsim.layout import LayoutConfig, simulate_layout

    a = simulate_layout(LayoutConfig(), np.random.default_rng(7))
    b = simulate_layout(LayoutConfig(short_convergent_gaps=False), np.random.default_rng(7))
    pd.testing.assert_frame_equal(a.genes, b.genes)
    pd.testing.assert_frame_equal(a.tus, b.tus)


def _digest(*frames: pd.DataFrame) -> str:
    m = hashlib.sha256()
    for f in frames:
        # 9 significant digits: the last-ulp float noise differs across CPUs (lead decision 20)
        m.update(f.to_csv(index=False, float_format="%.9g").encode())
    return m.hexdigest()


def test_default_layout_is_byte_identical_to_fd3f3a3() -> None:
    from giapsim.layout import LayoutConfig, simulate_layout

    g = simulate_layout(LayoutConfig(), np.random.default_rng(7)).genes
    assert _digest(g) == "c889d561360d77be959edf0ad618d38bcffb8fdd635be93c8c8b97f4f895ea7e"


@pytest.mark.parametrize(
    "name,args,digest",
    [
        ("operon_null", (600,),
         "6af312e99ef72ed5920517156a7245c0a7bf1716e5a4feeb2d9894261158fad6"),
        ("discordance_power", (600, 10, 2.0),
         "84d1c4a4cbaeee098b6fb899d49c64bd767d0b8a684429d76c59c9221db70935"),
        ("operon_attack_ar1", (600, 0.9),
         "baf078c879e022522831b1cca029684acbe369931f2110fc9a448a14af7a8c22"),
    ],
)  # fmt: skip
def test_scenarios_without_short_gaps_are_byte_identical_to_fd3f3a3(
    name: str, args: tuple[Any, ...], digest: str
) -> None:
    """Decision 11 "default off so existing scenarios are byte-identical" (and decision 13's
    default ``n_planted = 0`` on ``operon_attack_ar1``): seed 7, layout genes, expression genes and
    the pair table as CSV (floats to 9 significant digits), digests equal at fd3f3a3 and on
    iron-01 and barley-04 (lead decision 20)."""
    s = sim(name, args, 7)
    assert _digest(s.layout.genes, s.expression.genes, s.pairs) == digest


@pytest.mark.parametrize(
    "name,args,on",
    [
        ("operon_null", (N,), False),
        ("discordance_power", (N, 40, 2.0), False),
        ("operon_attack_ar1", (N, 0.5), False),
        ("missed_operon_power", (N, 50, 0.1), True),
        ("operon_attack_detection", (N, "prophage"), True),
    ],
)
def test_which_scenarios_turn_short_convergent_gaps_on(
    name: str, args: tuple[Any, ...], on: bool
) -> None:
    lcfg, _ = getattr(_sc(), name)(*args)
    assert lcfg.short_convergent_gaps is on
    assert tuple(lcfg.convergent_gap_bp) == (-20, 300)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_short_convergent_gaps_on_draws_cv_gaps_uniformly_from_minus20_to_300(seed: int) -> None:
    """On the null layout with the opt-in on: every (non-wrap) CV gap lies in -20-300 bp, gaps below
    50 bp occur, and DV and same-strand inter-TU gaps keep the 50-300 bp range."""
    sc = _sc()
    s = sc.simulate_operons(_with_layout(sc.operon_null(N), short_convergent_gaps=True), seed)
    gaps = _cv_gaps(s)
    assert len(gaps) >= 100
    assert gaps.min() >= -20 and gaps.max() <= 300, (gaps.min(), gaps.max())
    assert (gaps < 0).any() and (gaps > 150).any()  # 3' overlaps and the 150-300 bin (decision 23)
    assert (gaps < 50).mean() >= 0.15, (gaps < 50).mean()
    p = s.pairs[_not_wrap(s, s.pairs)]
    dv = p.loc[(p["class"] == "DV").to_numpy(), "intergenic_bp"].to_numpy()
    assert dv.min() >= 50 and dv.max() <= 300
    across = p[(p["class"].isin(["OP", "NO"])).to_numpy() & ~_b(p["true_operon"])]
    assert across["intergenic_bp"].min() >= 50


def test_short_convergent_gaps_are_uniform_pooled() -> None:
    """Pooled over 5 seeds, the CV gaps follow the discrete uniform on -20..300 (KS p > 1e-3)."""
    sc = _sc()
    cfgs = _with_layout(sc.operon_null(N), short_convergent_gaps=True)
    gaps = np.concatenate([_cv_gaps(sc.simulate_operons(cfgs, s)) for s in range(5)])
    u = stats.randint(-20, 301)
    # continuity: compare x + U(0, 1) against the continuous uniform on [-20, 301)
    jit = gaps + np.random.default_rng(0).random(len(gaps))
    assert stats.kstest(jit, stats.uniform(-20, 321).cdf).pvalue > 1e-3
    assert abs(gaps.mean() - u.mean()) <= 4 * u.std() / math.sqrt(len(gaps))


def test_convergent_gap_range_is_honoured() -> None:
    sc = _sc()
    cfgs = _with_layout(sc.operon_null(600), short_convergent_gaps=True, convergent_gap_bp=(20, 30))
    gaps = _cv_gaps(sc.simulate_operons(cfgs, 3))
    assert len(gaps) and gaps.min() >= 20 and gaps.max() <= 30


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_detection_attack_has_short_cv_decoys_outside_the_block(seed: int) -> None:
    """B1(c): CV decoys outside the short-gap block populate the < 50 bp bins."""
    s = sim("operon_attack_detection", (N, "prophage"), seed)
    gaps = _cv_gaps(s, _block_genes(s, "prophage"))
    assert (gaps < 50).sum() >= 10, np.sort(gaps)[:20]
    assert gaps.min() >= -20 and gaps.max() <= 300  # decision 23


def test_missed_operon_power_has_short_cv_decoys() -> None:
    s = sim("missed_operon_power", (N, 50, 0.1), 0)
    gaps = _cv_gaps(s)
    assert (gaps < 50).sum() >= 10 and gaps.max() <= 300  # decision 23


# ----------------------------------------------------- decision 10: abundance-dependent detection


def test_operon_config_accepts_operon_t_with_abundance_detection() -> None:
    from giapsim.operons import OperonConfig

    sc = _sc()
    beta = {**dict(sc.OPERON_BETA), "log2_abundance": 0.3}
    cfg = OperonConfig(operon_t=sc.OPERON_T, detection="abundance", beta=beta)
    assert cfg.detection == "abundance" and cfg.operon_t == sc.OPERON_T


def test_detection_attack_config_uses_abundance_detection() -> None:
    _, ecfg = _sc().operon_attack_detection(N, "prophage")
    assert ecfg.detection == "abundance"
    assert ecfg.beta["log2_abundance"] > 0
    assert ecfg.operon_t is not None


def _sigmoid_inv(p: np.ndarray) -> np.ndarray:
    return np.log(p) - np.log1p(-p)


@pytest.mark.parametrize("seed", [0, 1])
def test_detection_attack_logit_uses_the_post_walk_abundance(seed: int) -> None:
    """``logit(p_detect_true) = covariate logit + region shift + b4 * final abundance`` gene by
    gene (``b4 > 0``), where the final abundance is ``abundance_log2_true`` after the within-TU
    walk; a detection computed before the walk leaves a residual that tracks the walk's t."""
    _, ecfg = _sc().operon_attack_detection(N, "prophage")
    s = sim("operon_attack_detection", (N, "prophage"), seed)
    b = ecfg.beta
    assert b["log2_abundance"] > 0
    lg = s.layout.genes.set_index(s.layout.genes["locus_tag"].astype(str))
    ex = s.expression.genes
    tags = ex["locus_tag"].astype(str).to_numpy()
    block = _block_genes(s, "prophage")
    shift = {k: float(v) for k, v in ecfg.region_detect_shift.items()}
    assert set(shift) == {"prophage"}
    want = np.full(len(tags), float(b["intercept"]))
    for i, t in enumerate(tags):
        aa = lg.at[t, "aa"]
        if isinstance(aa, str) and aa:
            want[i] += (b["log_length"] * math.log(len(aa))
                        + b["log1p_n_obs"] * math.log1p(float(ex["n_obs_true"].iloc[i]))
                        + b["gravy"] * float(ex["gravy_true"].iloc[i]))  # fmt: skip
        if t in block:
            want[i] += shift["prophage"]
    want += b["log2_abundance"] * ex["abundance_log2_true"].to_numpy(float)
    got = _sigmoid_inv(ex["p_detect_true"].to_numpy(float))
    ok = np.isfinite(got)
    assert ok.mean() > 0.99
    np.testing.assert_allclose(got[ok], want[ok], atol=1e-6)


@pytest.mark.parametrize("seed", [0, 1])
def test_detection_attack_pair_columns_follow_the_detection_truth(seed: int) -> None:
    s = sim("operon_attack_detection", (N, "prophage"), seed)
    ex = s.expression.genes
    det = dict(zip(ex["locus_tag"].astype(str), _b(ex["detected_true"]), strict=True))
    p = s.pairs
    assert [det[str(g)] for g in p["gene_a"]] == list(_b(p["detected_a"]))
    assert [det[str(g)] for g in p["gene_b"]] == list(_b(p["detected_b"]))


def _block_pairs(s: Any) -> pd.DataFrame:
    block = _block_genes(s, "prophage")
    p = s.pairs[s.pairs["class"] != "SW"]
    inb = p["gene_a"].astype(str).isin(block) & p["gene_b"].astype(str).isin(block)
    return p[inb.to_numpy()]


def test_detection_attack_block_pairs_survive_the_detection_filter() -> None:
    """Decision 10: >= 25 % of the block's adjacent pairs have both genes detected (seeds 0-4
    pooled); the block still loses detection relative to the rest (attack present)."""
    k = n = 0
    for seed in range(5):
        s = sim("operon_attack_detection", (N, "prophage"), seed)
        bp = _block_pairs(s)
        k += int((_b(bp["detected_a"]) & _b(bp["detected_b"])).sum())
        n += len(bp)
        ex = s.expression.genes
        inside = ex["locus_tag"].astype(str).isin(_block_genes(s, "prophage")).to_numpy()
        det = _b(ex["detected_true"])
        assert det[inside].mean() < det[~inside].mean(), seed
    assert n >= 50
    assert k / n >= 0.25, (k, n)


# ------------------------------------------------------- decision 13: planted pairs on attacks


ATTACKS = [("operon_attack_ar1", (N, 0.5)), ("operon_attack_detection", (N, "prophage"))]


@pytest.mark.parametrize("name,args", ATTACKS)
def test_attack_planting_defaults_to_nothing(name: str, args: tuple[Any, ...]) -> None:
    f = getattr(_sc(), name)
    assert f(*args, n_planted=0, effect_sd=0.0) == f(*args)
    _, ecfg = f(*args)
    assert ecfg.n_missed == 0 and ecfg.missed_effect_sd == 0.0
    s = sim(name, args, 0)
    assert len(s.truth["missed_operons"]) == 0


@pytest.mark.parametrize("name,args", ATTACKS)
def test_attack_planting_is_the_missed_operon_mechanism(name: str, args: tuple[Any, ...]) -> None:
    f = getattr(_sc(), name)
    base_l, base_e = f(*args)
    lcfg, ecfg = f(*args, n_planted=20, effect_sd=0.1)
    assert lcfg == base_l
    assert ecfg == dataclasses.replace(base_e, n_missed=20, missed_effect_sd=0.1)


@pytest.mark.parametrize("name,args", ATTACKS)
@pytest.mark.parametrize("seed", [0, 1])
def test_attack_planted_pairs_are_true_no_pairs_with_small_delta(
    name: str, args: tuple[Any, ...], seed: int
) -> None:
    s = sim(name, args, seed, (("n_planted", 20), ("effect_sd", 0.1)))
    m = s.truth["missed_operons"]
    assert list(m.columns[:3]) == ["pair_id", "gene_a", "gene_b"]
    assert len(m) == 20 and m["pair_id"].is_unique
    p = s.pairs.set_index("pair_id")
    rows = p.loc[m["pair_id"]]
    assert (rows["class"] == "NO").all()
    assert not _b(rows["true_operon"]).any()
    assert (rows["uniop_p"].astype(float) < 0.5).all()
    assert list(rows["gene_a"].astype(str)) == list(m["gene_a"].astype(str))
    g = s.expression.genes
    ab = dict(zip(g["locus_tag"].astype(str), g["abundance_log2_true"].astype(float), strict=True))

    def absd(fr: pd.DataFrame) -> np.ndarray:
        return np.array([abs(ab[str(b)] - ab[str(a)]) for a, b in zip(fr["gene_a"], fr["gene_b"],
                                                                      strict=True)])  # fmt: skip

    others = p[(p["class"] == "NO").to_numpy() & ~_b(p["true_operon"])]
    others = others.loc[~others.index.isin(m["pair_id"])]
    assert np.median(absd(rows)) <= 0.25 * np.median(absd(others))


@pytest.mark.parametrize("name,args", ATTACKS)
def test_attack_planting_writes_missed_operons_tsv(
    name: str, args: tuple[Any, ...], tmp_path: Path
) -> None:
    from giapsim.contract_tables import read_operon_truth, write_operon_truth

    s = sim(name, args, 0, (("n_planted", 20), ("effect_sd", 0.1)))
    paths = write_operon_truth(s, tmp_path)
    f = Path(paths["missed_operons"])
    assert f == tmp_path / "truth" / "missed_operons.tsv"
    with f.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0][:3] == ["pair_id", "gene_a", "gene_b"]
    assert len(rows) == 21
    back = read_operon_truth(tmp_path)["missed_operons"]
    assert list(back["pair_id"].astype(str)) == list(s.truth["missed_operons"]["pair_id"])


def test_attack_planting_too_many_pairs_raises() -> None:
    sc = _sc()
    with pytest.raises(ValueError, match="NO pairs"):
        sc.simulate_operons(sc.operon_attack_ar1(200, 0.5, n_planted=10_000, effect_sd=0.1), 0)


# ------------------------------------------------------------- decision 14: the null's config


def test_null_config_has_rho_zero_and_nothing_planted() -> None:
    lcfg, ecfg = _sc().operon_null(N)
    assert ecfg.ar1_rho == 0.0
    assert ecfg.n_missed == 0 and ecfg.missed_effect_sd == 0.0
    assert ecfg.n_discordant == 0 and ecfg.discordant_effect_log2 == 0.0
    assert ecfg.islands == () and len(ecfg.region_effect_log2) == 0
    assert lcfg.uniop_pure is True and lcfg.short_gaps_in is None
    s = sim("operon_null", (N,), 0)
    assert len(s.truth["missed_operons"]) == 0 and len(s.truth["discordant_pairs"]) == 0
    assert not _b(s.pairs["promoter_between_true"] & s.pairs["true_operon"]).any()
