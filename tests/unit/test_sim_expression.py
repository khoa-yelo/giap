"""giapsim expression and detection truth (issue M4-08; GIAP_SPEC.md §11.3, §8.4, §8.5, §7.3.6).

Examples 1-8, edge cases and acceptance criteria of docs/specs/M4-08-sim-expression.md, checked
against tests/helpers/expression_oracle.py (written from the spec, not from the simulator):

* the detection probability is exactly the spec logistic in (log length, log1p n_obs, GRAVY
  [, log2 abundance]) + the region detection shift, with covariates from the M4-10 definitions;
* TU-level values are one per TU; plasmid copy number, region effects and islands are additive
  shifts (paired runs with the same seed differ by exactly the planted shift);
* the AR(1) background has the stated autocorrelation; the frames are a pure function of
  (layout, config, seed); the layout and the M4-07/M0 outputs are unchanged (pinned hashes);
* expression genes join 1:1 with the M4-06 contract genes loaded by ``load_m4_inputs``.

Example 3 and 4 tolerances (+/-0.5) are below the single-seed SD of a 40-/60-gene block under the
default AR(1) (about 0.8), so the literal examples are checked on the mean over 30 seeds and the
exact shift on paired runs. Lives in tests/unit (tests/sim is implementer tooling).

Ambiguities resolved here (reported to the lead): ``tu_sd`` is the AR(1) innovation SD sigma of the
spec formula; ``tu_log2`` is the TU-level value before gene-level terms (mu + replicon + region +
island + AR(1)); island ``first TU index`` is 0-based in replicon TU order; ``n_obs`` counts
distinct peptides (M4-10); ``log`` is natural; the expression frame has one row per layout gene.
"""

from __future__ import annotations

import hashlib
import math
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo
import layout_oracle as lo

REPO = Path(__file__).resolve().parents[2]
GIAPSIM = REPO / "tests" / "sim" / "giapsim"
GENE_COLUMNS = (
    "locus_tag", "tu_log2", "abundance_log2_true", "p_detect_true", "detected_true",
    "n_obs_true", "gravy_true", "in_island", "region_effect_log2",
)  # fmt: skip
BETA_KEYS = {"intercept", "log_length", "log1p_n_obs", "gravy", "log2_abundance"}

# Contract tables of ``simulate_layout(LayoutConfig(), default_rng(1))`` and of a mixed layout
# (main at b3846d7, before M4-08): M4-08 must not change the layout or its writer.
LAYOUT_SHA256 = {
    "default": {
        "annotation_meta.json": "0998816b81d3f81364ac87e43ef1ffc2d40c08bc620c535584e859619e215f72",
        "genes.parquet": "e0d4341901746a7444965ebe8070a3c382d5edc6f864bff8597bd8e697b65e65",
        "layers.parquet": "c7ae231f33b641c11a6b7b03ce95c36e5daf5c6b7fd5c1a46456d31e6d6c5ba5",
        "region_genes.parquet": "1197b721b0b595beaf182adb375dfbf0707c5cf8085b5d4598106988f1215fd7",
        "regions.parquet": "54ca39d8465b7720f4417364363d63ce4c6aa92cd893cad24c2f311c825f2ca5",
        "replicons.parquet": "baa9f94adf047f02a746031e33777abc727b0eea9e650454e0b7a868562bd85e",
        "tus.parquet": "c564de32b5fac4d9d1dc6a9ff5f885065ccdec0c2d7547b7ce48561f0e1ca825",
        "xrefs.parquet": "0fd3d2531b4d94b2dc6370de0378ecb870858dc2f6b091c37854771e474f592e",
    },
    "mixed": {
        "annotation_meta.json": "0998816b81d3f81364ac87e43ef1ffc2d40c08bc620c535584e859619e215f72",
        "genes.parquet": "36b633ae9875d5a51b30e2999070ff6dc35b03c6011ef0aa3a42df6fea54e2de",
        "layers.parquet": "1745a48214adc69e34c25c54c61906b07ac121d43a097bbfd10044afc0514171",
        "region_genes.parquet": "acea093cd168190d31138ba706167c04a25a5186727512189982206c4480369a",
        "regions.parquet": "c61cc5fa60051964a7435dc93e62a89687fb60c52ead2e53ff314b6055ecbe50",
        "replicons.parquet": "5e38eeaf9aa9a76a63410d12fdb7296f7d489d63ad5a41afe5c9adf211e6bab1",
        "tus.parquet": "7ad19f2cfb66e2c5f9bbe397f0d468fc7ce12ed1bbf47d0c3324305b965e1719",
        "xrefs.parquet": "ded2195445f5002254a0fbbdcd1abc402787ed8cf09621e98c66d9dac0cbed80",
    },
}
MIXED_KW: dict[str, Any] = {
    "n_genes_chromosome": 400, "plasmids": (25, 1), "linear_contigs": (12, 1),
    "origin_crossing": True,
}  # fmt: skip
# M0 + M3-19 proteomics (``SimConfig(proteomics=ProteomicsConfig())``, seed 1; main at b3846d7):
# the M3-19 hook swap (``expression_for_simulation``) must be opt-in.
M0_PROTEOMICS_SHA256 = {
    "abundance.tsv": "bdbe3b7329f9a369bb278ab48d82c925462027b7c719c6f2bfa4e922534fd34a",
    "annotation.gff3": "b9b0482bb56fb40a4069be0d9bce4629169779c41ebb2243506d8da57628a924",
    "bakta_like.json": "d312420eb3d404cb1abbb6d07844cccbc83c0ce5e9c912d3a02aa17da6c49064",
    "genome.fasta": "698373a45c787caa64c01d7050dd1d2a7b31e03771e4e389056ca510c6df1971",
    "peptides.tsv": "8489d448b3241131140f1ec59f28c0465aca8875e70f188550df4705328ee852",
    "proteins.faa": "138b2faa715a1f995fe510073c451253913d46589af8d2759a84df986286fb29",
    "truth/features.tsv": "14c5e4866fe3beb761ec7652455e4cc90e7496ac3370b6def7eb1c0fc03c6075",
    "truth/genes.tsv": "9c2c98152eed218d3855ccddf478bf49ea42bc69727fce466111f704e0387b0f",
    "truth/novel_orfs.tsv": "57b8bf862b97cc98740f605389ee2c51d5b7da03f2e03ce6bb42827c7a7da9f2",
    "truth/peptides.tsv": "4ddfb9b94088026d55e29e21550a030952d2f0fcc47bd2457816f62104cd0e4b",
}


# ----------------------------------------------------------------------------- helpers


def _expr_mod() -> ModuleType:
    import giapsim.expression

    return giapsim.expression


def layout(seed: int = 1, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def cfg(**kw: Any) -> Any:
    return _expr_mod().ExpressionConfig(**kw)


def express(lay: Any, seed: int = 1, **kw: Any) -> Any:
    m = _expr_mod()
    return m.simulate_expression(lay, m.ExpressionConfig(**kw), np.random.default_rng(seed))


def island(*args: Any) -> Any:
    return _expr_mod().IslandSpec(*args)


def by_tag(expr: Any) -> pd.DataFrame:
    return expr.genes.set_index("locus_tag")


@pytest.fixture(scope="module")
def lay1() -> Any:
    return layout(1)


@pytest.fixture(scope="module")
def lay4500() -> Any:
    return layout(1, n_genes_chromosome=4500)


@pytest.fixture(scope="module")
def ex1(lay1: Any) -> Any:
    return express(lay1, 1)


# ----------------------------------------------------------------------------- interface


def test_config_defaults_match_spec() -> None:
    c = cfg()
    assert c.mu_log2 == 20.0
    assert (c.tu_sd, c.gene_sd, c.stoich_sd) == (2.0, 0.5, 0.7)
    assert c.ar1_rho == 0.5
    assert c.position_decay_log2 == 0.2
    assert dict(c.plasmid_copy_log2) == {}
    assert dict(c.region_effect_log2) == {}
    assert dict(c.region_detect_shift) == {}
    assert tuple(c.islands) == ()
    assert c.detection == "covariates_only"
    assert dict(c.beta) == dict(_expr_mod().DEFAULT_BETA)
    assert set(_expr_mod().DEFAULT_BETA) == BETA_KEYS


def test_config_frozen_and_hashable() -> None:
    c = cfg(plasmid_copy_log2={"p1": 1.0}, islands=(island("chr1", 0, 5, 1.0),))
    with pytest.raises((AttributeError, TypeError)):
        c.mu_log2 = 1.0  # type: ignore[misc]
    hash(c)
    assert c == cfg(plasmid_copy_log2={"p1": 1.0}, islands=(island("chr1", 0, 5, 1.0),))


def test_default_beta_signs_follow_lead_proposal() -> None:
    """Q1 (proposed default): longer, more observable, less hydrophobic proteins detect better."""
    b = _expr_mod().DEFAULT_BETA
    assert b["log_length"] > 0 and b["log1p_n_obs"] > 0 and b["gravy"] < 0
    assert b["intercept"] < 0


def test_gene_frame_columns_and_one_row_per_gene(ex1: Any, lay1: Any) -> None:
    g = ex1.genes
    assert set(GENE_COLUMNS) <= set(g.columns)
    assert g["locus_tag"].is_unique
    assert sorted(g["locus_tag"]) == sorted(lay1.genes["locus_tag"])
    assert g["detected_true"].isin([True, False]).all()
    assert g["in_island"].isin([True, False]).all()
    for c in ("tu_log2", "abundance_log2_true", "p_detect_true", "gravy_true",
              "region_effect_log2"):  # fmt: skip
        assert pd.api.types.is_float_dtype(g[c]), c
    assert pd.api.types.is_integer_dtype(g["n_obs_true"]) or (
        (g["n_obs_true"] == np.round(g["n_obs_true"])).all()
    )


def test_tu_log2_is_one_value_per_tu(ex1: Any, lay1: Any) -> None:
    tl = by_tag(ex1)["tu_log2"]
    for tid, members in eo.tu_members(lay1).items():
        vals = {float(tl[g]) for g in members}
        assert len(vals) == 1, (tid, vals)


# ----------------------------------------------------------------------------- example 1


def test_example1_p_detect_in_open_unit_interval(ex1: Any) -> None:
    p = ex1.genes["p_detect_true"].to_numpy(dtype=float)
    assert np.isfinite(p).all()
    assert ((p > 0) & (p < 1)).all()


def test_example1_detected_share(ex1: Any) -> None:
    share = float(ex1.genes["detected_true"].astype(bool).mean())
    assert abs(share - 0.55) <= 0.15, share


def test_example1_detection_logit_exact(ex1: Any, lay1: Any) -> None:
    """``p_detect_true = sigmoid(b0 + b1 log len + b2 log1p n_obs + b3 GRAVY)``, covariates_only."""
    m = eo.merged(lay1, ex1.genes)
    m = m[~m["aa"].map(eo.is_na)]
    x = eo.expected_logit(m, _expr_mod().DEFAULT_BETA, "covariates_only", {}, {})
    np.testing.assert_allclose(m["p_detect_true"].to_numpy(float), eo.sigmoid(x), rtol=0,
                               atol=1e-9)  # fmt: skip


def test_covariate_columns_match_oracle(ex1: Any, lay1: Any) -> None:
    m = eo.merged(lay1, ex1.genes)
    m = m[~m["aa"].map(eo.is_na)]
    assert [int(v) for v in m["n_obs_true"]] == [eo.n_obs(a) for a in m["aa"]]
    np.testing.assert_allclose(m["gravy_true"].to_numpy(float), [eo.gravy(a) for a in m["aa"]],
                               atol=1e-9)  # fmt: skip


def test_covariates_agree_with_giap_m4_10(ex1: Any, lay1: Any) -> None:
    """M4-13 refits on giap's covariates: the simulator's truth must use the same definitions."""
    from giap.genetable.covariates import covariates

    m = eo.merged(lay1, ex1.genes)
    m = m[~m["aa"].map(eo.is_na)].head(300)
    for aa, n, gr in zip(m["aa"], m["n_obs_true"], m["gravy_true"], strict=True):
        c = covariates(aa, "trypsin_p", (7, 30))
        assert c.n_observable_peptides == int(n)
        assert c.gravy == pytest.approx(float(gr), abs=1e-9)


def test_partial_and_pseudogenes_keep_values() -> None:
    """Partial genes (linear-contig ends) and pseudogenes keep finite values (giap drops them)."""
    lay = layout(3, n_genes_chromosome=300, linear_contigs=(10, 10), frac_pseudogene=0.2)
    g = eo.merged(lay, express(lay, 3).genes)
    sub = g[g["partial"].astype(bool) | (g["feature_type"] == "pseudogene")]
    assert len(sub) >= 4
    for c in ("tu_log2", "abundance_log2_true", "p_detect_true"):
        assert np.isfinite(sub[c].to_numpy(float)).all(), c
    p = sub["p_detect_true"].to_numpy(float)
    assert ((p > 0) & (p < 1)).all()


def test_detected_true_is_bernoulli_of_p(lay4500: Any) -> None:
    """Total detections within 4 SD of sum p (Poisson-binomial); deeper check in calibration."""
    g = express(lay4500, 11).genes
    p = g["p_detect_true"].to_numpy(float)
    k = int(g["detected_true"].astype(bool).sum())
    assert abs(k - p.sum()) <= 4 * math.sqrt(float((p * (1 - p)).sum())) + 1


# ----------------------------------------------------------------------------- example 2, AR(1)


def test_example2_lag1_autocorrelation_rho_09(lay4500: Any) -> None:
    x = eo.tu_series(lay4500, express(lay4500, 1, ar1_rho=0.9).genes)["chr1"]
    r = eo.lag1_autocorr(x)
    assert abs(r - 0.9) <= 0.05, r


def test_rho_zero_gives_iid_tu_effects(lay4500: Any) -> None:
    x = eo.tu_series(lay4500, express(lay4500, 2, ar1_rho=0.0).genes)["chr1"]
    n = len(x)
    assert abs(eo.lag1_autocorr(x)) <= 4 / math.sqrt(n)


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
def test_tu_level_mean_and_innovation_sd(lay4500: Any, rho: float) -> None:
    """Whitened TU series ~ iid N(0, 1): mean mu, innovation SD ``tu_sd`` (no effects planted)."""
    x = eo.tu_series(lay4500, express(lay4500, 5, ar1_rho=rho).genes)["chr1"] - 20.0
    z = eo.ar1_whiten(x, rho, 2.0)
    n = len(z)
    assert abs(z.mean()) <= 4 / math.sqrt(n), z.mean()
    assert abs(z.var() - 1) <= 4 * math.sqrt(2 / n), z.var()
    assert abs(eo.lag1_autocorr(z)) <= 4 / math.sqrt(n)


@pytest.mark.parametrize("rho", [-0.1, 1.0, 1.5])
def test_rho_outside_unit_interval_rejected(lay1: Any, rho: float) -> None:
    with pytest.raises(ValueError):
        express(lay1, 1, ar1_rho=rho)


def test_ar1_restarts_per_replicon_and_single_tu_is_stationary() -> None:
    """A one-TU replicon gets the stationary draw: Var = sigma^2/(1 - rho^2)
    (rho = 0.9 -> 21.05 for sigma = 2),
    independent of the previous replicon's last TU (no chaining across replicons)."""
    rho, sigma = 0.9, 2.0
    vals, prev = [], []
    for s in range(15):
        lay = layout(100 + s, n_genes_chromosome=40, regions=(), plasmids=(1,) * 20)
        series = eo.tu_series(lay, express(lay, s, ar1_rho=rho).genes)
        last = float(series["chr1"][-1]) - 20.0
        for k in range(1, 21):
            x = series[f"p{k}"]
            assert len(x) == 1
            vals.append(float(x[0]) - 20.0)
            prev.append(last)
            last = vals[-1]
    v = np.array(vals)
    target = sigma**2 / (1 - rho**2)
    assert 0.75 * target <= float(np.mean(v**2)) <= 1.3 * target, np.mean(v**2)
    assert abs(float(np.corrcoef(v, prev)[0, 1])) <= 4 / math.sqrt(len(v))


# ----------------------------------------------------------------------------- gene level


def test_gene_level_position_decay_and_noise(lay4500: Any) -> None:
    """abundance - tu_log2 = -decay*position + stoichiometry + noise: slope -0.2 per TU position,
    residual SD sqrt(0.7^2 + 0.5^2)."""
    g = by_tag(express(lay4500, 4))
    pos = eo.gene_tu_position(lay4500)
    tags = list(g.index)
    k = np.array([pos[t][1] for t in tags], dtype=float)
    d = (g["abundance_log2_true"] - g["tu_log2"]).to_numpy(float)
    slope, icpt = np.polyfit(k, d, 1)
    assert abs(slope + 0.2) <= 0.06, slope
    resid = d - (slope * k + icpt)
    sd = math.sqrt(0.7**2 + 0.5**2)
    assert abs(resid.std() - sd) <= 0.06, resid.std()


def test_position_decay_knob_honoured(lay4500: Any) -> None:
    g = by_tag(express(lay4500, 4, position_decay_log2=1.0, stoich_sd=0.1, gene_sd=0.1))
    pos = eo.gene_tu_position(lay4500)
    k = np.array([pos[t][1] for t in g.index], dtype=float)
    d = (g["abundance_log2_true"] - g["tu_log2"]).to_numpy(float)
    slope = np.polyfit(k, d, 1)[0]
    assert abs(slope + 1.0) <= 0.03, slope
    assert np.std(d + k) <= 0.2  # sqrt(0.1^2 + 0.1^2) = 0.14


# ----------------------------------------------------------------------------- example 3, regions


def test_example3_prophage_effect_paired_exact(lay1: Any) -> None:
    """Same seed with and without the effect: prophage-only TUs shift by exactly -3."""
    base = by_tag(express(lay1, 1))
    eff = by_tag(express(lay1, 1, region_effect_log2={"prophage": -3.0}))
    status = eo.tu_region_status(lay1, "prophage")
    mem = eo.tu_members(lay1)
    assert any(s == "all" for s in status.values())
    for t, s in status.items():
        for g in mem[t]:
            da = eff.at[g, "abundance_log2_true"] - base.at[g, "abundance_log2_true"]
            dt = eff.at[g, "tu_log2"] - base.at[g, "tu_log2"]
            if s == "all":
                assert da == pytest.approx(-3.0, abs=1e-9) and dt == pytest.approx(-3.0, abs=1e-9)
                assert eff.at[g, "region_effect_log2"] == pytest.approx(-3.0)
            elif s == "none":
                assert da == pytest.approx(0.0, abs=1e-9) and dt == pytest.approx(0.0, abs=1e-9)
                assert eff.at[g, "region_effect_log2"] == pytest.approx(0.0)
            else:
                assert min(abs(da), abs(da + 3.0)) <= 1e-9
    assert (base["region_effect_log2"] == 0).all()


def test_example3_prophage_mean_difference_pooled() -> None:
    """mean(prophage) - mean(elsewhere) within -3 +/- 0.5 (mean of 30 seeds; one seed: SD ~ 0.8)."""
    diffs = []
    for s in range(30):
        lay = layout(500 + s)
        g = eo.merged(lay, express(lay, s, region_effect_log2={"prophage": -3.0}).genes)
        rt = eo.region_types_of_genes(lay)
        inp = g["locus_tag"].map(lambda t, rt=rt: "prophage" in rt.get(t, set()))
        diffs.append(g.loc[inp, "abundance_log2_true"].mean()
                     - g.loc[~inp, "abundance_log2_true"].mean())  # fmt: skip
    assert abs(float(np.mean(diffs)) + 3.0) <= 0.5, np.mean(diffs)


# ----------------------------------------------------------------------------- example 4, plasmids


def test_example4_plasmid_copy_paired_exact() -> None:
    lay = layout(2, plasmids=(60,))
    shift = math.log2(10)
    base = by_tag(express(lay, 2))
    eff = by_tag(express(lay, 2, plasmid_copy_log2={"p1": shift}))
    contig = dict(zip(lay.genes["locus_tag"], lay.genes["contig"], strict=True))
    for g in base.index:
        want = shift if contig[g] == "p1" else 0.0
        for c in ("tu_log2", "abundance_log2_true"):
            assert eff.at[g, c] - base.at[g, c] == pytest.approx(want, abs=1e-9), (g, c)


def test_example4_plasmid_shift_pooled() -> None:
    """Plasmid genes shifted by +3.32 +/- 0.5 against the chromosome (mean over 30 seeds)."""
    diffs = []
    for s in range(30):
        lay = layout(700 + s, plasmids=(60,))
        g = eo.merged(lay, express(lay, s, plasmid_copy_log2={"p1": math.log2(10)}).genes)
        on = g["contig"] == "p1"
        diffs.append(g.loc[on, "abundance_log2_true"].mean()
                     - g.loc[~on, "abundance_log2_true"].mean())  # fmt: skip
    assert abs(float(np.mean(diffs)) - math.log2(10)) <= 0.5, np.mean(diffs)


def test_unknown_replicon_in_plasmid_copy_has_no_effect_or_raises() -> None:
    lay = layout(2, plasmids=(20,))
    try:
        e = express(lay, 2, plasmid_copy_log2={"p9": 5.0})
    except (KeyError, ValueError):
        return
    pd.testing.assert_frame_equal(e.genes, express(lay, 2).genes)


# ----------------------------------------------------------------------------- example 5, BGC shift


def test_example5_bgc_detection_shift_exact(lay1: Any) -> None:
    shift = eo.logit(0.7) - eo.logit(0.5)
    e = express(lay1, 1, region_detect_shift={"bgc": shift})
    m = eo.merged(lay1, e.genes)
    m = m[~m["aa"].map(eo.is_na)]
    rt = eo.region_types_of_genes(lay1)
    x = eo.expected_logit(m, _expr_mod().DEFAULT_BETA, "covariates_only", {"bgc": shift}, rt)
    np.testing.assert_allclose(m["p_detect_true"].to_numpy(float), eo.sigmoid(x), atol=1e-9)
    # genes whose covariate-only p is 0.5 move to 0.7 inside the BGC
    base = eo.expected_logit(m, _expr_mod().DEFAULT_BETA, "covariates_only", {}, {})
    inb = m["locus_tag"].map(lambda t: "bgc" in rt.get(t, set())).to_numpy()
    assert inb.sum() >= 20
    p_in = m["p_detect_true"].to_numpy(float)[inb]
    lifted = eo.sigmoid(base[inb] + shift)
    np.testing.assert_allclose(p_in, lifted, atol=1e-9)
    assert float(eo.sigmoid(eo.logit(0.5) + shift)) == pytest.approx(0.7)


def test_region_detect_shift_leaves_abundance_unchanged(lay1: Any) -> None:
    a = express(lay1, 1).genes
    b = express(lay1, 1, region_detect_shift={"bgc": 0.85}).genes
    pd.testing.assert_series_equal(a["abundance_log2_true"], b["abundance_log2_true"])
    pd.testing.assert_series_equal(a["tu_log2"], b["tu_log2"])


# ----------------------------------------------------------------------------- example 6


def test_example6_covariates_only_independent_of_abundance(lay4500: Any) -> None:
    g = eo.merged(lay4500, express(lay4500, 1).genes)
    g = g[~g["aa"].map(eo.is_na)]
    X = eo.covariate_matrix(g, extra=("abundance_log2_true",))
    X[:, -1] = X[:, -1] - X[:, -1].mean()
    y = g["detected_true"].astype(float).to_numpy()
    b, se = eo.logistic_fit(X, y)
    z = b[-1] / se[-1]
    assert abs(z) < 3, (b, se)


def test_covariates_only_ignores_abundance_coefficient(lay1: Any) -> None:
    """b4 is used only when ``detection = "abundance"``."""
    beta = dict(_expr_mod().DEFAULT_BETA)
    beta["log2_abundance"] = 0.7
    beta["intercept"] = beta["intercept"] - 0.7 * 20.0
    g = eo.merged(lay1, express(lay1, 1, beta=beta, detection="covariates_only").genes)
    g = g[~g["aa"].map(eo.is_na)]
    x = eo.expected_logit(g, beta, "covariates_only", {}, {})
    np.testing.assert_allclose(g["p_detect_true"].to_numpy(float), eo.sigmoid(x), atol=1e-9)


def test_abundance_mode_logit_exact_and_dependent(lay4500: Any) -> None:
    beta = dict(_expr_mod().DEFAULT_BETA)
    beta["log2_abundance"] = 0.7
    beta["intercept"] = beta["intercept"] - 0.7 * 20.0
    g = eo.merged(lay4500, express(lay4500, 1, beta=beta, detection="abundance").genes)
    g = g[~g["aa"].map(eo.is_na)]
    x = eo.expected_logit(g, beta, "abundance", {}, {})
    np.testing.assert_allclose(g["p_detect_true"].to_numpy(float), eo.sigmoid(x), atol=1e-9)
    X = eo.covariate_matrix(g, extra=("abundance_log2_true",))
    X[:, -1] = X[:, -1] - 20.0
    b, se = eo.logistic_fit(X, g["detected_true"].astype(float).to_numpy())
    assert b[-1] / se[-1] > 5
    assert abs(b[-1] - 0.7) <= 4 * se[-1]


def test_unknown_detection_mode_rejected(lay1: Any) -> None:
    with pytest.raises(ValueError):
        express(lay1, 1, detection="both")


# ----------------------------------------------------------------------------- example 7


def test_example7_simple_covariates() -> None:
    aa = "MAGLLKAAGPEPTIDERPAAAGGSTLKAAAAAAAK"
    n, gr = _expr_mod().simple_covariates(aa)
    assert n == 3
    assert isinstance(n, int)
    assert gr == pytest.approx(eo.gravy(aa), abs=1e-12)


@pytest.mark.parametrize(
    ("aa", "n"),
    [
        ("MAGLLKAAGPEPTIDERPAAAGGSTLK" + "W" * 31 + "KAAAAAAAK", 3),  # 32-aa run excluded
        ("AAAAAAAKAAAAAAAKGGGGGGGR", 2),  # distinct peptides
        ("AAAAAAA", 1),  # C-terminal peptide without K/R counts; exactly 7 aa
        ("AAAAAAKAAAAAA", 1),  # 7 aa counts, 6 aa does not
        ("A" * 30, 1),
        ("A" * 31, 0),
        ("KKKKKKKKKK", 0),
    ],
)
def test_simple_covariates_n_obs_cases(aa: str, n: int) -> None:
    assert eo.n_obs(aa) == n  # oracle sanity
    assert _expr_mod().simple_covariates(aa)[0] == n


def test_simple_covariates_gravy_values() -> None:
    sc = _expr_mod().simple_covariates
    assert sc("MKTAYIAKQR")[1] == pytest.approx(-0.78, abs=1e-9)
    assert sc("IIII")[1] == pytest.approx(4.5)
    assert sc("RRRR")[1] == pytest.approx(-4.5)


# ----------------------------------------------------------------------------- example 8


def test_example8_same_seed_identical_frames(lay1: Any) -> None:
    kw = {"region_effect_log2": {"prophage": -3.0}, "islands": (island("chr1", 3, 6, 1.5),)}
    a, b = express(lay1, 8, **kw), express(lay1, 8, **kw)
    pd.testing.assert_frame_equal(a.genes, b.genes)
    pd.testing.assert_frame_equal(a.islands, b.islands)
    c = express(lay1, 9, **kw)
    assert not np.array_equal(a.genes["abundance_log2_true"], c.genes["abundance_log2_true"])


def test_rng_is_the_only_source_of_randomness(lay1: Any) -> None:
    m = _expr_mod()
    np.random.seed(0)  # noqa: NPY002 (legacy global state must not matter)
    a = m.simulate_expression(lay1, m.ExpressionConfig(), np.random.default_rng(4))
    np.random.seed(12345)  # noqa: NPY002
    b = m.simulate_expression(lay1, m.ExpressionConfig(), np.random.default_rng(4))
    pd.testing.assert_frame_equal(a.genes, b.genes)


def test_layout_not_mutated(lay1: Any) -> None:
    fresh = layout(1)
    express(lay1, 1, region_effect_log2={"prophage": -3.0}, islands=(island("chr1", 0, 5, 2.0),))
    for name in ("replicons", "genes", "tus", "operon_pairs", "regions", "region_genes", "truth",
                 "xrefs"):  # fmt: skip
        pd.testing.assert_frame_equal(getattr(lay1, name), getattr(fresh, name), obj=name)


# ----------------------------------------------------------------------------- islands


def _island_tus(lay: Any, contig: str, first: int, n: int) -> list[str]:
    return eo.tu_order(lay)[contig][first : first + n]


@pytest.mark.parametrize("effect", [2.0, -1.5])
def test_island_shift_paired_exact(lay1: Any, effect: float) -> None:
    """rho = 0 (stationary SD = sigma): island TUs shift by effect*tu_sd, nothing else moves."""
    spec = island("chr1", 10, 8, effect)
    base = by_tag(express(lay1, 3, ar1_rho=0.0))
    e = express(lay1, 3, ar1_rho=0.0, islands=(spec,))
    eff = by_tag(e)
    mem = eo.tu_members(lay1)
    in_tus = set(_island_tus(lay1, "chr1", 10, 8))
    in_genes = {g for t in in_tus for g in mem[t]}
    assert set(eff.index[eff["in_island"].astype(bool)]) == in_genes
    assert not base["in_island"].astype(bool).any()
    for g in eff.index:
        want = effect * 2.0 if g in in_genes else 0.0
        assert eff.at[g, "tu_log2"] - base.at[g, "tu_log2"] == pytest.approx(want, abs=1e-9)
        assert eff.at[g, "abundance_log2_true"] - base.at[g, "abundance_log2_true"] == (
            pytest.approx(want, abs=1e-9)
        )
    assert len(e.islands) == 1


def test_islands_truth_frame_rows_and_columns(lay1: Any) -> None:
    specs = (island("chr1", 0, 5, 1.0), island("chr1", 40, 20, -2.0))
    two = express(lay1, 1, islands=specs).islands
    none = express(lay1, 1).islands
    assert isinstance(none, pd.DataFrame) and len(none) == 0
    assert len(two) == 2
    assert list(none.columns) == list(two.columns) and len(two.columns) >= 3


def test_edge_empty_islands_frame_has_columns(ex1: Any) -> None:
    assert len(ex1.islands) == 0
    assert len(ex1.islands.columns) >= 3
    assert not ex1.genes["in_island"].astype(bool).any()


# ----------------------------------------------------------------------------- attack nulls


def test_attack2_short_prophage_genes_detected_less(lay1: Any) -> None:
    """§8.5 attack 2 input: short prophage genes have lower p_detect (confounded by covariates)."""
    lay = layout(1, short_genes_in="prophage")
    g = eo.merged(lay, express(lay, 1).genes)
    rt = eo.region_types_of_genes(lay)
    inp = g["locus_tag"].map(lambda t: "prophage" in rt.get(t, set()))
    assert g.loc[inp, "p_detect_true"].mean() < g.loc[~inp, "p_detect_true"].mean() - 0.1


def test_attack4_high_copy_plasmid_abundance() -> None:
    lay = layout(4, plasmids=(60,))
    g = eo.merged(lay, express(lay, 4, plasmid_copy_log2={"p1": math.log2(20)}).genes)
    on = g["contig"] == "p1"
    assert g.loc[on, "abundance_log2_true"].mean() > g.loc[~on, "abundance_log2_true"].mean() + 2


# ----------------------------------------------------------------------------- M4-06 contract


def test_contract_tables_still_valid_and_join_expression(tmp_path: Path, lay1: Any,
                                                         ex1: Any) -> None:  # fmt: skip
    """The written C1-C8 load with no error; every expression row is one contract gene."""
    from giap.genetable.contract import validate_inputs
    from giapsim.contract_tables import write_contract_tables

    write_contract_tables(lay1, tmp_path / "work")
    inputs, issues = lo.load_work(tmp_path / "work")
    assert lo.errors(issues) == []
    assert validate_inputs(inputs) == []
    tags = list(inputs.genes["locus_tag"])
    assert sorted(tags) == sorted(ex1.genes["locus_tag"])
    assert ex1.genes["locus_tag"].is_unique


# ----------------------------------------------------------------------------- backward compat


@pytest.mark.parametrize("name", ["default", "mixed"])
def test_layout_outputs_unchanged(tmp_path: Path, name: str) -> None:
    from giapsim.contract_tables import write_contract_tables

    lay = layout(1) if name == "default" else layout(1, **MIXED_KW)
    express(lay, 1, region_effect_log2={"prophage": -3.0})  # the feature must not touch it
    write_contract_tables(lay, tmp_path / "w")
    got = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
           sorted((tmp_path / "w").iterdir())}  # fmt: skip
    assert got == LAYOUT_SHA256[name]


def test_m0_proteomics_outputs_unchanged_with_hook_off(tmp_path: Path) -> None:
    """``expression_for_simulation`` is opt-in: the default M0 + M3-19 run is byte-identical."""
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    _expr_mod()  # the feature exists (so this test guards its default)
    giapsim.simulate(giapsim.SimConfig(proteomics=ProteomicsConfig()), 1, tmp_path)
    for rel, digest in M0_PROTEOMICS_SHA256.items():
        assert hashlib.sha256((tmp_path / rel).read_bytes()).hexdigest() == digest, rel


def test_expression_for_simulation_hook(tmp_path: Path) -> None:
    """The M3-19 hook: one finite log abundance per M0 gene, deterministic in the seed."""
    import giapsim
    from giapsim.layout import layout_from_simulation

    m = _expr_mod()
    res = giapsim.simulate(giapsim.SimConfig(), 1, tmp_path)
    lay = layout_from_simulation(res)
    a = np.asarray(m.expression_for_simulation(lay.genes, lay.tus, m.ExpressionConfig(),
                                               np.random.default_rng(2)), dtype=float)  # fmt: skip
    b = np.asarray(m.expression_for_simulation(lay.genes, lay.tus, m.ExpressionConfig(),
                                               np.random.default_rng(2)), dtype=float)  # fmt: skip
    assert a.shape == (len(lay.genes),)
    assert np.isfinite(a).all()
    np.testing.assert_array_equal(a, b)
    # genes of one TU share the TU term: within-TU spread is smaller than across TUs
    assert a.std() > 0


# ----------------------------------------------------------------------------- acceptance


def test_speed_4500_genes_under_one_second(lay4500: Any) -> None:
    m = _expr_mod()
    best = math.inf
    for s in range(3):
        t0 = time.perf_counter()
        m.simulate_expression(lay4500, m.ExpressionConfig(), np.random.default_rng(s))
        best = min(best, time.perf_counter() - t0)
    assert best < 1.0, best


def test_expression_module_imports_nothing_from_giap() -> None:
    """Shared grep test covers expression.py; also no ``giap`` module is loaded at import."""
    assert (GIAPSIM / "expression.py").is_file()
    code = (
        "import sys, giapsim.expression as e; "
        "bad = [m for m in sys.modules if m == 'giap' or m.startswith('giap.')]; "
        "print(bad); sys.exit(1 if bad else 0)"
    )
    env_path = f"{REPO / 'tests' / 'sim'}"
    env = {"PYTHONPATH": env_path, "PATH": "/usr/bin:/bin"}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env,
                       check=False)  # fmt: skip
    assert r.returncode == 0, r.stdout + r.stderr


def test_dataclass_replace_round_trip() -> None:
    c = cfg()
    assert replace(c, ar1_rho=0.9).ar1_rho == 0.9
