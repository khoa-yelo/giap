"""Small synthetic ``OUT/`` trees for the V15 ``determinants`` tests (M9-20).

Tables go through ``helpers.viz_out_tree.write_schema_table`` (schema header; absent cells
``NA``): ``tables/determinants_model.tsv``, ``tables/determinants_attribution.tsv`` and
``tables/determinants_features.tsv`` as the M8 ``determinants`` step writes them (every
``delta_r2`` row flagged ``ci_unestimated``, every ``cv_r2_mean`` row too, ``null_p`` rows
``null_uncalibrated``; M8-10). The M9-17 support files are written from the contract in
``helpers.support_analysis_data`` (not from M9-17 code) to ``plots/support/<name>.tsv``:
``determinants_shap`` (model rows x the top features, ``feature`` an encoded model column) and
``determinants_pd`` (``cai, tai, tir_ostir, promotech_max`` order, 21 ascending grid values on
the feature's own scale, ``pd_low <= pd_mean <= pd_high``; LightGBM only). Nothing here touches
``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from helpers import support_analysis_data as sa
from helpers.viz_out_tree import make_out, write_events, write_schema_table

BLOCKS = ("transcription", "translation", "codon", "protein", "context")
UNEXPLAINED = "abundance_not_explained_by_sequence"
#: 22 encoded model columns, by decreasing SHAP scale (the first 20 are the top 20)
SHAP_FEATURES: tuple[str, ...] = (
    "tir_ostir", "cai", "promotech_max", "mrna_dg", "tai", "gc3", "length_aa", "sd_dg",
    "tx_rate_max_tu", "enc", "pi", "gravy", "start_codon_ATG", "fop", "tu_position", "mw",
    "replichore_position", "gc_local_dev", "cpb", "spacer_bp", "at_2_10", "in_prophage",
)  # fmt: skip
N_MODEL = 120
RIBOSOMAL = ("RP_0001", "RP_0002")
UNEXPLAINED_TAGS = ("G_0007", "G_0042")
NA_PRED_TAG = "G_0099"  # an attribution row with no prediction (``unestimated``)


def mrow(
    model: str,
    what: str,
    metric: str,
    value: float | None,
    lo: float | None = None,
    hi: float | None = None,
    flags: str | None = None,
    n: int | None = N_MODEL,
    n_perm: int = 0,
) -> dict[str, object]:
    """One ``determinants_model.tsv`` row."""
    return {"model": model, "feature_or_block": what, "metric": metric, "value": value,
            "ci_low": lo, "ci_high": hi, "n": n, "n_perm": n_perm, "flags": flags}  # fmt: skip


def delta(model: str, block: str, value: float, lo: float | None, hi: float | None) -> dict:
    return mrow(model, block, "delta_r2", value, lo, hi, "ci_unestimated")


def cv_rows(model: str, mean: float, sd: float, null_p: float) -> list[dict[str, object]]:
    return [
        mrow(model, "all", "cv_r2_mean", mean, flags="ci_unestimated"),
        mrow(model, "all", "cv_r2_sd", sd),
        mrow(model, "all", "null_p", null_p, flags="null_uncalibrated", n_perm=49),
    ]


#: example 1: elastic net 0.21 ± 0.03, transcription 0.08 (0.05, 0.11) the largest; LightGBM's
#: largest block is codon (0.12), so a title naming codon read the wrong model
EN_DELTAS: Mapping[str, tuple[float, float | None, float | None]] = {
    "transcription": (0.08, 0.05, 0.11),
    "translation": (0.03, 0.01, 0.05),
    "codon": (0.02, -0.01, 0.05),
    "protein": (0.01, -0.02, 0.04),
    "context": (0.005, None, None),  # CI NA: point without a bar, "CI unestimated"
}
LGB_DELTAS: Mapping[str, tuple[float, float | None, float | None]] = {
    "transcription": (0.06, 0.02, 0.10),
    "translation": (0.04, 0.01, 0.07),
    "codon": (0.12, 0.08, 0.16),
    "protein": (0.01, -0.01, 0.03),
}  # no context row: a block with no feature left has no row (M8 rule)
#: example 2: every elastic-net CI includes 0 (LightGBM's still exclude it)
EN_DELTAS_NULL: Mapping[str, tuple[float, float | None, float | None]] = {
    "transcription": (0.02, -0.01, 0.05),
    "translation": (0.01, -0.02, 0.04),
    "codon": (0.0, -0.03, 0.03),
    "protein": (-0.01, -0.04, 0.02),
    "context": (0.005, -0.02, 0.03),
}


def univariate_rows() -> list[dict[str, object]]:
    return [
        mrow("univariate", "cai", "spearman_rho", 0.31, flags="ci_unestimated"),
        mrow("univariate", "tai", "spearman_rho", 0.22, flags="ci_unestimated"),
    ]


def model_table(
    en: Mapping[str, tuple[float, float | None, float | None]] = EN_DELTAS,
    lgb: Mapping[str, tuple[float, float | None, float | None]] = LGB_DELTAS,
    *,
    en_cv: tuple[float, float, float] = (0.21, 0.03, 0.02),
    lgb_cv: tuple[float, float, float] = (0.25, 0.04, 0.04),
    univariate: bool = True,
) -> list[dict[str, object]]:
    """``determinants_model.tsv`` rows in the M8 order (model, metric, block)."""
    rows = cv_rows("lightgbm", *lgb_cv)
    rows += [delta("lightgbm", b, *v) for b, v in lgb.items()]
    rows += cv_rows("elastic_net", *en_cv)
    rows += [delta("elastic_net", b, *v) for b, v in en.items()]
    return rows + (univariate_rows() if univariate else [])


def too_few_model_table() -> list[dict[str, object]]:
    return [
        mrow(m, "all", "cv_r2_mean", None, flags="too_few_rows;ci_unestimated", n=None)
        for m in ("lightgbm", "elastic_net")
    ] + univariate_rows()


# ---------------------------------------------------------------- genes


def model_tags(n: int = N_MODEL) -> list[str]:
    return [f"G_{i:04d}" for i in range(n)]


def tu_of(i: int) -> str:
    return f"TU{i // 3:04d}"  # 3 model genes per TU


def attribution_rows(n: int = N_MODEL, *, seed: int = 15) -> list[dict[str, object]]:
    """Model rows (``predicted_log2`` = observed + noise; two flagged unexplained; one with NA
    prediction), then two ``cai_reference`` ribosomal rows (not model rows)."""
    rng = np.random.default_rng(seed)
    obs = 20.0 + rng.normal(0.0, 1.5, n)
    pred = obs + rng.normal(0.0, 0.6, n)
    rows: list[dict[str, object]] = []
    for i, tag in enumerate(model_tags(n)):
        flags: list[str] = []
        p: float | None = round(float(pred[i]), 4)
        if tag in UNEXPLAINED_TAGS:
            p = round(float(obs[i]) - 4.0, 4)
            flags.append(UNEXPLAINED)
        if tag == NA_PRED_TAG:
            p = None
            flags.append("unestimated")
        a = round(float(obs[i]), 4)
        rows.append({"locus_tag": tag, "abundance_log2": a, "predicted_log2": p,
                     "residual": None if p is None else round(a - p, 4),
                     "top_decile": False, "flags": ";".join(flags) or None})  # fmt: skip
    for j, tag in enumerate(RIBOSOMAL):
        rows.append({"locus_tag": tag, "abundance_log2": 25.0 + j, "predicted_log2": 23.5 + j,
                     "residual": 1.5, "top_decile": False, "flags": "cai_reference"})  # fmt: skip
    return rows


def too_few_attribution(n: int = 60) -> list[dict[str, object]]:
    rng = np.random.default_rng(3)
    return [
        {"locus_tag": t, "abundance_log2": round(float(20 + rng.normal()), 4),
         "predicted_log2": None, "residual": None, "top_decile": False, "flags": "too_few_rows"}
        for t in model_tags(n)
    ]  # fmt: skip


def feature_rows(
    n: int = N_MODEL, *, tu: int = 3, extra_non_model: int = 0, seed: int = 16
) -> list[dict[str, object]]:
    """``determinants_features.tsv`` rows: ``n`` model rows (``tu`` genes per TU), the
    ribosomal genes, and ``extra_non_model`` undetected genes in TUs of their own."""
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for i, tag in enumerate(model_tags(n)):
        rows.append({
            "locus_tag": tag, "tu_id": f"TU{i // tu:04d}", "detected": True,
            "abundance_log2": 20.0, "abundance_imputed": False, "ribosomal_protein": False,
            "cai": round(float(rng.uniform(0.2, 0.9)), 4),
            "tai": round(float(rng.uniform(0.1, 0.6)), 4),
            "tir_ostir": round(float(rng.lognormal(6.0, 1.0)), 2),
            "promotech_max": round(float(rng.uniform(0.0, 1.0)), 4),
            "replicon_type": "chromosome",
        })  # fmt: skip
    for tag in RIBOSOMAL:
        rows.append({"locus_tag": tag, "tu_id": "TU_RP", "detected": True, "abundance_log2": 25.0,
                     "abundance_imputed": False, "ribosomal_protein": True, "cai": 0.95,
                     "tai": 0.5, "tir_ostir": 2000.0, "promotech_max": 0.9,
                     "replicon_type": "chromosome"})  # fmt: skip
    for k in range(extra_non_model):
        rows.append({"locus_tag": f"N_{k:04d}", "tu_id": f"TUN{k:04d}", "detected": False,
                     "abundance_log2": None, "abundance_imputed": False,
                     "ribosomal_protein": False, "replicon_type": "chromosome"})  # fmt: skip
    return rows


# ---------------------------------------------------------------- support files


def shap_frame(
    tags: Sequence[str] | None = None, features: Sequence[str] = SHAP_FEATURES[:20], seed: int = 17
) -> pd.DataFrame:
    """``determinants_shap``: every tag x every feature; feature ``j`` has SHAP scale
    ``1 / (j + 1)``, so mean |SHAP| strictly decreases along ``features``."""
    tags = [t for t in model_tags() if t != NA_PRED_TAG] if tags is None else list(tags)
    rng = np.random.default_rng(seed)
    rows = []
    for j, f in enumerate(features):
        vals = rng.uniform(0.0, 1.0, len(tags))
        shap = rng.normal(0.0, 1.0, len(tags))
        shap = shap / np.mean(np.abs(shap)) / (j + 1)  # mean |SHAP| exactly 1/(j+1)
        for t, v, s in zip(tags, vals, shap, strict=True):
            rows.append((t, f, round(float(v), 6), round(float(s), 6)))
    return pd.DataFrame(rows, columns=list(sa.COLUMNS["determinants_shap"]))


def pd_frame(features: Sequence[str] = sa.PD_FEATURES) -> pd.DataFrame:
    """``determinants_pd``: 21 ascending grid values per feature on its own scale."""
    ranges = {"cai": (0.2, 0.9), "tai": (0.1, 0.6), "tir_ostir": (50.0, 5000.0),
              "promotech_max": (0.0, 1.0)}  # fmt: skip
    rows = []
    for k, f in enumerate(features):
        lo, hi = ranges[f]
        grid = np.linspace(lo, hi, sa.PD_GRID_POINTS)
        for i, g in enumerate(grid):
            mean = 19.5 + 0.05 * (k + 1) * i
            rows.append((f, round(float(g), 6), round(mean, 6), round(mean - 0.2, 6),
                         round(mean + 0.15, 6)))  # fmt: skip
    return pd.DataFrame(rows, columns=list(sa.COLUMNS["determinants_pd"]))


def write_support(out: Path, name: str, frame: pd.DataFrame) -> Path:
    from giap.io.tsv import write_tsv

    sup = out / "plots" / "support"
    sup.mkdir(parents=True, exist_ok=True)
    path = sup / f"{name}.tsv"
    write_tsv(path, [str(c) for c in frame.columns], frame.itertuples(index=False, name=None))
    return path


# ---------------------------------------------------------------- the OUT tree


def build_out(
    tmp: Path,
    model: Sequence[Mapping[str, object]] | None = None,
    attribution: Sequence[Mapping[str, object]] | None = None,
    features: Sequence[Mapping[str, object]] | None = None,
    *,
    shap: pd.DataFrame | None = None,
    pd_: pd.DataFrame | None = None,
    support: bool = True,
    flags: Sequence[str] = (),
    no_model: bool = False,
) -> Path:
    """An ``OUT/`` with the three determinants tables (defaults: the main fixture) and, when
    ``support``, the two M9-17 support files; ``flags`` go to a step events file."""
    out = make_out(tmp)
    if not no_model:
        write_schema_table(out / "tables" / "determinants_model.tsv", "determinants_model",
                           model_table() if model is None else model)  # fmt: skip
    write_schema_table(
        out / "tables" / "determinants_attribution.tsv", "determinants_attribution",
        attribution_rows() if attribution is None else attribution,
    )  # fmt: skip
    write_schema_table(out / "tables" / "determinants_features.tsv", "determinants_features",
                       feature_rows() if features is None else features)  # fmt: skip
    if support:
        write_support(out, "determinants_shap", shap_frame() if shap is None else shap)
        write_support(out, "determinants_pd", pd_frame() if pd_ is None else pd_)
    if flags:
        write_events(
            out / "logs" / "steps" / "determinants.events.jsonl",
            [{"event": "determinants.layer_missing", "data": {"flags": list(flags)},
              "step": "determinants", "level": "WARNING",
              "ts": "2026-10-23T00:00:00+00:00"}],
        )  # fmt: skip
    return out


def many_attribution(n: int) -> list[dict[str, object]]:
    """``n`` model rows (for the > 20,000-gene binning rule)."""
    rng = np.random.default_rng(99)
    obs = 20.0 + rng.normal(0.0, 1.5, n)
    pred = obs + rng.normal(0.0, 0.6, n)
    return [
        {"locus_tag": f"G_{i:06d}", "abundance_log2": round(float(o), 4),
         "predicted_log2": round(float(p), 4), "residual": round(float(o - p), 4),
         "top_decile": False, "flags": None}
        for i, (o, p) in enumerate(zip(obs, pred, strict=True))
    ]  # fmt: skip


def many_features(n: int) -> list[dict[str, object]]:
    return [
        {"locus_tag": f"G_{i:06d}", "tu_id": f"TU{i // 3:06d}", "detected": True,
         "abundance_log2": 20.0, "abundance_imputed": False, "ribosomal_protein": False,
         "cai": 0.5, "tai": 0.3, "tir_ostir": 400.0, "promotech_max": 0.5,
         "replicon_type": "chromosome"}
        for i in range(n)
    ]  # fmt: skip
