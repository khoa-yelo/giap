"""gippasim null and power scenarios (issue M5-04; GIPPA_SPEC.md §8.4, §8.5 attack 1, §11.3).

Examples 8 and 9 and the scenario definitions of docs/specs/M5-04-sim-oric-gradient.md (the spec's
``tests/sim/test_scenarios.py``; lives in tests/unit like the M4-07/M4-08 tests), plus the
scenarios of lead decision 13 and the ADR-0020 attack list. The exact-level checks of the two
nulls are in tests/calibration/test_sim_scenarios_calibration.py.

Interface assumed for lead decision 13 (not in the spec's signatures; reported to the lead):
``label_null(n, *, oriter_log2=0.0)`` and ``tu_null(n, detection="abundance", *,
short_genes_in=None, oriter_log2=0.0)`` add the cosine gradient with ``oric_frac = 0.8`` when
``oriter_log2 != 0``; ``island_power(n, rho, n_tus, effect_sd, *, detection="covariates_only")``;
``draft_ar1(rho)`` (500-gene chr1 + 133 linear 30-gene contigs).

Example 8's within-TU pairs are consecutive genes of a TU (genome order); the between pairs are
the last gene of a TU and the first gene of the next TU on chr1 (genome order). Example 9's Moran's
I uses lag-1 neighbours along chr1 genome order (no wrap pair).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import oric_oracle as oo

SCENARIO_NAMES = ("label_null", "tu_null", "ar1_background", "island_power", "oriter_gradient",
                  "draft_ar1")  # fmt: skip


def _sc() -> Any:
    import gippasim.scenarios

    return gippasim.scenarios


def run(cfgs: tuple[Any, Any], seed: int) -> tuple[Any, Any]:
    from gippasim.expression import simulate_expression
    from gippasim.layout import simulate_layout

    lcfg, ecfg = cfgs
    lay = simulate_layout(lcfg, np.random.default_rng(seed))
    return lay, simulate_expression(lay, ecfg, np.random.default_rng(10_000 + seed))


def chr1_genome_order(lay: Any, ex: Any) -> pd.DataFrame:
    """chr1 genes in genome order with their TU and ``abundance_log2_true``."""
    g = lay.genes[lay.genes["contig"] == "chr1"][["locus_tag", "start"]]
    tu = dict(zip(lay.truth["locus_tag"], lay.truth["tu_id"], strict=True))
    ab = dict(zip(ex.genes["locus_tag"], ex.genes["abundance_log2_true"], strict=True))
    g = g.sort_values("start", kind="mergesort")
    return pd.DataFrame({"locus_tag": g["locus_tag"].to_numpy(),
                         "tu_id": [tu[t] for t in g["locus_tag"]],
                         "ab": [float(ab[t]) for t in g["locus_tag"]]})  # fmt: skip


def lag1_pairs(df: pd.DataFrame) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """(within-TU, between-TU) neighbour pairs along chr1 genome order."""
    within, between = [], []
    ab, tu = df["ab"].tolist(), df["tu_id"].tolist()
    for i in range(len(df) - 1):
        (within if tu[i] == tu[i + 1] else between).append((ab[i], ab[i + 1]))
    return within, between


def r_of(pairs: list[tuple[float, float]]) -> float:
    x, y = zip(*pairs, strict=True)
    return oo.pearson(x, y)


# ----------------------------------------------------------------------------- interface


def test_every_scenario_exists() -> None:
    m = _sc()
    for name in SCENARIO_NAMES:
        assert callable(getattr(m, name)), name


@pytest.mark.parametrize(
    "call",
    [
        lambda m: m.label_null(300),
        lambda m: m.tu_null(300),
        lambda m: m.tu_null(300, detection="covariates_only"),
        lambda m: m.ar1_background(300, 0.5),
        lambda m: m.island_power(300, 0.5, 5, 2.0),
        lambda m: m.oriter_gradient(300, 0.5, 2.0),
        lambda m: m.draft_ar1(0.5),
    ],
    ids=["label_null", "tu_null", "tu_null_cov", "ar1", "island", "oriter", "draft"],
)
def test_scenarios_return_config_pairs_that_run(call: Any) -> None:
    from gippasim.expression import ExpressionConfig
    from gippasim.layout import LayoutConfig

    cfgs = call(_sc())
    assert isinstance(cfgs, tuple) and len(cfgs) == 2
    lcfg, ecfg = cfgs
    assert isinstance(lcfg, LayoutConfig) and isinstance(ecfg, ExpressionConfig)
    assert call(_sc()) == cfgs  # a pure function of its arguments
    lay, ex = run(cfgs, 1)
    assert len(ex.genes) == len(lay.genes)


def _common(lcfg: Any, ecfg: Any, n: int) -> None:
    """Spec: all scenarios use the circular chr1 of ``n`` genes and no other replicon,
    ``covariates_only`` unless stated, no planted region effects."""
    assert lcfg.n_genes_chromosome == n
    assert tuple(lcfg.plasmids) == () and tuple(lcfg.linear_contigs) == ()
    assert dict(ecfg.region_effect_log2) == {}
    assert dict(ecfg.region_detect_shift) == {}
    assert dict(ecfg.plasmid_copy_log2) == {}


def test_label_null_definition() -> None:
    lcfg, ecfg = _sc().label_null(1000)
    _common(lcfg, ecfg, 1000)
    assert ecfg.tu_sd == 0.0 and ecfg.position_decay_log2 == 0.0
    assert ecfg.stoich_sd > 0 or ecfg.gene_sd > 0
    assert ecfg.detection == "covariates_only"
    assert tuple(ecfg.islands) == () and ecfg.oriter_log2 == 0.0
    assert lcfg.oric_frac is None


def test_abundance_beta_is_pinned() -> None:
    """Review B2: ``ABUNDANCE_BETA`` is ``DEFAULT_BETA`` with slope 0.5 on log2 abundance and the
    intercept offset ``-0.5 * mu_log2`` (so the detected share stays near the default)."""
    from gippasim.expression import DEFAULT_BETA, ExpressionConfig

    beta = dict(_sc().ABUNDANCE_BETA)
    mu = ExpressionConfig().mu_log2
    assert beta["log2_abundance"] == 0.5
    assert beta["intercept"] == pytest.approx(DEFAULT_BETA["intercept"] - 0.5 * mu)
    assert beta["intercept"] == pytest.approx(-18.65)
    rest = {k: v for k, v in beta.items() if k not in ("log2_abundance", "intercept")}
    assert rest == {k: v for k, v in dict(DEFAULT_BETA).items()
                    if k not in ("log2_abundance", "intercept")}  # fmt: skip


def test_abundance_detection_share_matches_the_default() -> None:
    """Review B2 behaviourally: the ``abundance`` detected share stays near the
    ``covariates_only`` one (0.534 vs 0.555 over 8 seeds); dropping the ``-b4 mu`` offset or
    doubling the slope moves it."""
    shares = {}
    for name, cfgs in (("abundance", _sc().tu_null(1000)),
                       ("covariates_only", _sc().tu_null(1000, "covariates_only"))):  # fmt: skip
        vals = [float(run(cfgs, seed)[1].genes["detected_true"].mean()) for seed in range(8)]
        shares[name] = float(np.mean(vals))
    assert shares["abundance"] == pytest.approx(0.534, abs=0.02)
    assert shares["covariates_only"] == pytest.approx(0.555, abs=0.02)
    assert abs(shares["abundance"] - shares["covariates_only"]) < 0.05


@pytest.mark.parametrize("detection", ["abundance", "covariates_only"])
def test_tu_null_definition(detection: str) -> None:
    m = _sc()
    lcfg, ecfg = m.tu_null(1000) if detection == "abundance" else m.tu_null(1000, detection)
    _common(lcfg, ecfg, 1000)
    assert ecfg.ar1_rho == 0.0 and ecfg.tu_sd == 2.0 and ecfg.position_decay_log2 == 0.2
    assert ecfg.detection == detection
    assert tuple(ecfg.islands) == () and ecfg.oriter_log2 == 0.0
    if detection == "abundance":
        assert dict(ecfg.beta) == dict(_sc().ABUNDANCE_BETA)  # review B2
    else:
        from gippasim.expression import DEFAULT_BETA

        assert dict(ecfg.beta) == dict(DEFAULT_BETA)


def test_tu_null_detection_inherits_tu_structure() -> None:
    """``detection = "abundance"``: the detection logit is correlated within TUs."""
    within, between = [], []
    for seed in range(5):
        lay, ex = run(_sc().tu_null(1000), seed)
        p = dict(zip(ex.genes["locus_tag"], ex.genes["p_detect_true"], strict=True))
        df = chr1_genome_order(lay, ex)
        q = np.clip([float(p[t]) for t in df["locus_tag"]], 1e-12, 1 - 1e-12)
        df["ab"] = np.log(q / (1 - q))
        w, b = lag1_pairs(df)
        within += w
        between += b
    assert r_of(within) > r_of(between) + 0.1


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
def test_ar1_background_definition(rho: float) -> None:
    from gippasim.expression import ExpressionConfig

    lcfg, ecfg = _sc().ar1_background(1000, rho)
    _common(lcfg, ecfg, 1000)
    d = ExpressionConfig()
    assert ecfg.ar1_rho == rho
    for f in ("mu_log2", "tu_sd", "gene_sd", "stoich_sd", "position_decay_log2"):
        assert getattr(ecfg, f) == getattr(d, f), f
    assert tuple(ecfg.islands) == () and ecfg.oriter_log2 == 0.0
    assert ecfg.detection == "covariates_only"
    assert lcfg.oric_frac is None


def test_ar1_background_09_has_strong_tu_autocorrelation() -> None:
    """Lead decision 13: ``ar1_background(1000, 0.9)``; the TU series keeps rho = 0.9."""
    import expression_oracle as eo

    lay, ex = run(_sc().ar1_background(1000, 0.9), 3)
    x = eo.tu_series(lay, ex.genes)["chr1"]
    assert eo.lag1_autocorr(x) > 0.75


@pytest.mark.parametrize("effect", [2.0, -1.5])
def test_island_power_definition(effect: float) -> None:
    lcfg, ecfg = _sc().island_power(1000, 0.5, 10, effect)
    _common(lcfg, ecfg, 1000)
    assert ecfg.ar1_rho == 0.5
    assert len(ecfg.islands) == 1
    isl = ecfg.islands[0]
    assert (isl.replicon, isl.n_tus, isl.effect) == ("chr1", 10, effect)  # sign = direction
    assert isl.first_tu >= 0
    assert ecfg.detection == "covariates_only" and ecfg.oriter_log2 == 0.0


def test_island_power_fits_every_layout_seed() -> None:
    """The seeded island start fits the chr1 TUs of any layout drawn from the scenario."""
    cfgs = _sc().island_power(1000, 0.5, 20, 2.0)
    for seed in range(10):
        _, ex = run(cfgs, seed)
        assert ex.genes["in_island"].any()
        assert ex.islands["n_tus"].tolist() == [20]


def test_island_power_abundance_detection() -> None:
    """Lead decision 13: ``island_power(..., detection = "abundance")`` (review B2: the pinned
    ``ABUNDANCE_BETA``, and the mode must reach the config)."""
    _, ecfg = _sc().island_power(1000, 0.5, 10, 2.0, detection="abundance")
    assert ecfg.detection == "abundance"
    assert dict(ecfg.beta) == dict(_sc().ABUNDANCE_BETA)


def test_island_first_tu_is_a_quarter_of_the_expected_tus() -> None:
    """S3: ``island_first_tu(1000) == 138`` (1,000 genes / (4 x mean TU size 1.8))."""
    m = _sc()
    assert m.island_first_tu(1000) == 138
    assert m.island_first_tu(4500) == 625
    assert m.island_first_tu(1000, mean_tu_size=2.5) == 100
    assert m.island_power(1000, 0.5, 10, 2.0)[1].islands[0].first_tu == 138


def test_oriter_gradient_definition() -> None:
    lcfg, ecfg = _sc().oriter_gradient(1000, 0.5, 2.0)
    _common(lcfg, ecfg, 1000)
    assert lcfg.oric_frac == 0.8
    assert ecfg.oriter_log2 == 2.0 and ecfg.oriter_shape == "cosine"
    assert ecfg.ar1_rho == 0.5 and tuple(ecfg.islands) == ()
    lay, _ = run((lcfg, ecfg), 1)
    assert (lay.features["type"] == "oric").sum() == 1


# ----------------------------------------------------------------------------- lead decision 13


def test_label_null_with_gradient() -> None:
    """``label_null`` + gradient: ``tu_sd = 0``, cosine 2 log2 (ADR-0020 attack b)."""
    lcfg, ecfg = _sc().label_null(1000, oriter_log2=2.0)
    assert lcfg.oric_frac == 0.8
    assert ecfg.tu_sd == 0.0 and ecfg.oriter_log2 == 2.0 and ecfg.oriter_shape == "cosine"
    _, ex = run((lcfg, ecfg), 2)
    assert ex.genes["replichore_true"].notna().all()


def test_tu_null_with_gradient() -> None:
    """ADR-0020 attack a: TU null + 2 log2 cosine."""
    lcfg, ecfg = _sc().tu_null(1000, oriter_log2=2.0)
    assert lcfg.oric_frac == 0.8
    assert ecfg.ar1_rho == 0.0 and ecfg.tu_sd == 2.0 and ecfg.oriter_log2 == 2.0


def test_tu_null_short_gene_prophage() -> None:
    """Lead decision 13 (ADR-0020 attack c): short prophage genes under the TU null with
    ``covariates_only`` detection."""
    lcfg, ecfg = _sc().tu_null(1000, detection="covariates_only", short_genes_in="prophage")
    assert lcfg.short_genes_in == "prophage"
    assert ecfg.detection == "covariates_only" and ecfg.ar1_rho == 0.0


def test_draft_ar1_layout() -> None:
    lcfg, ecfg = _sc().draft_ar1(0.5)
    assert lcfg.n_genes_chromosome == 500
    assert tuple(lcfg.linear_contigs) == (30,) * 133
    assert tuple(lcfg.plasmids) == ()
    assert ecfg.ar1_rho == 0.5 and tuple(ecfg.islands) == ()
    lay, _ = run((lcfg, ecfg), 1)
    assert (lay.replicons["topology"] == "linear").sum() == 133


# ----------------------------------------------------------------------------- examples 8, 9


def test_example8_tu_null_correlation_within_not_between() -> None:
    """``tu_null(1000)``: within-TU lag-1 r > 0.5; between TUs |r| < 0.1 (20 seeds pooled)."""
    within, between = [], []
    for seed in range(20):
        lay, ex = run(_sc().tu_null(1000), seed)
        w, b = lag1_pairs(chr1_genome_order(lay, ex))
        within += w
        between += b
    assert r_of(within) > 0.5
    assert abs(r_of(between)) < 0.1


def test_example9_label_null_morans_i_near_zero() -> None:
    """``label_null(1000)``: |Moran's I lag 1| < 0.08 in at least 95 of 100 seeds."""
    small = 0
    for seed in range(100):
        lay, ex = run(_sc().label_null(1000), seed)
        small += abs(oo.morans_i_lag1(chr1_genome_order(lay, ex)["ab"].to_numpy())) < 0.08
    assert small >= 95


def test_label_null_within_tu_uncorrelated() -> None:
    """Both nulls hold under ``label_null``: no within-TU correlation either (20 seeds)."""
    within = []
    for seed in range(20):
        lay, ex = run(_sc().label_null(1000), seed)
        within += lag1_pairs(chr1_genome_order(lay, ex))[0]
    assert abs(r_of(within)) < 0.1
