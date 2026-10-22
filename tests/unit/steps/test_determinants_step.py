"""M8-10 step ``determinants``: registry entry, the three tables, flags, streams and threads
(docs/specs/M8-10-determinants-step-and-run.md examples 1 (below the minimum), 3, 4, 5, 6 and
steps 1-8; GIAP_SPEC.md §6.2, §7.8.2, §7.8.3, §9.3, §9.4, §6.6).

The fixture is the 1,000-gene M4-07 layout with M4-08 expression (``build_m4_fixture``
``sim1000``) through the M4 chain, so the step reads the real C1/C2/C5/C6/C7/C8 tables, the M4
gene table and M5-09's gene context. Example 3 runs with B = 49 (the spec's unit budget) and a
small CV (2 repeats x 3 folds, 30 boosting rounds): every row kind is produced, in seconds.

Example 1 (sim100k, ~90 genes, below the 100-row minimum) is checked here on the same fixture
with abundance kept on 60 detected genes (same code path, no 40-s ``giap run``); the real
``giap run --executor local`` on sim100k is tests/integration/test_determinants_run.py.
Contract choices: tests/helpers/determinants_step_data.py.
"""

from __future__ import annotations

import json
import logging
import math
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_step_data as s

SEED = 20261022
REPO = Path(__file__).resolve().parents[3]


@pytest.fixture(scope="module")
def step_mod() -> ModuleType:
    """First fixture, so a missing step fails before the fixture is built."""
    return s.mod()


@pytest.fixture(scope="module")
def work(step_mod: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> Path:
    return s.build_work(tmp_path_factory.mktemp("m810_work"), SEED)


def _events(records: list[logging.LogRecord], name: str) -> list[dict[str, Any]]:
    return [r.data for r in records if getattr(r, "event", None) == name]  # type: ignore[attr-defined]


class _Capture(logging.Handler):
    """Every record of ``giap.steps`` (module-scoped fixtures cannot use ``caplog``)."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _run_captured(cfg: Any, ins: dict[str, Path], out: Path) -> tuple[dict[str, Path], list]:
    cap = _Capture()
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(cap)
    log.setLevel(logging.DEBUG)
    try:
        paths = s.run(cfg, ins, out)
    finally:
        log.removeHandler(cap)
        log.setLevel(old)
    return paths, cap.records


@pytest.fixture(scope="module")
def full(work: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Example 3: enough model rows; threads 1."""
    ins = s.inputs(work)
    cfg = s.config(SEED, threads=1)
    paths, records = _run_captured(cfg, ins, tmp_path_factory.mktemp("m810_full"))
    return {"ins": ins, "cfg": cfg, "paths": paths, "records": records}


@pytest.fixture(scope="module")
def model(full: dict[str, Any]) -> pd.DataFrame:
    return s.read(full["paths"]["determinants_model"])


@pytest.fixture(scope="module")
def attribution(full: dict[str, Any]) -> pd.DataFrame:
    return s.read(full["paths"]["determinants_attribution"])


@pytest.fixture(scope="module")
def expected(full: dict[str, Any]) -> dict[str, Any]:
    """What the M8-07/M8-08 functions give on the same inputs (rows, encoded columns, blocks)."""
    from giap.determinants.features import model_matrix
    from giap.determinants.models import model_rows

    f = s.expected_features(full["ins"])
    X, blocks = model_matrix(f)
    rows = model_rows(f)
    return {"features": f, "X": X, "blocks": blocks, "rows": rows}


# ---------------------------------------------------------------------------- registry


def test_registry_entry(step_mod: ModuleType) -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("determinants")
    assert spec.kind == "py" and spec.func is not None
    assert spec.module == "core", "Module 8 is core (§6.2): never switched off by `modules`"
    assert spec.inputs == s.REQUIRED
    assert set(spec.optional_inputs) == set(s.OPTIONAL)
    assert set(spec.outputs) == set(s.OUTPUTS)
    assert spec.optional_outputs == ()
    assert spec.precomputed is None


def test_runs_after_gene_table_and_context_stats(step_mod: ModuleType) -> None:
    """Step 1: ``gene_context`` is published by CONTEXT, so ``determinants`` runs after
    ``context_stats`` locally too (topological order)."""
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    order = reg.topo_order(["determinants", "context_stats", "gene_table", "detectability"])
    assert order.index("determinants") > order.index("context_stats")
    assert order.index("determinants") > order.index("gene_table")


def test_local_executor_publishes_the_three_tables(step_mod: ModuleType) -> None:
    """Step 7: ``tables/determinants_*.tsv`` only, under the names Nextflow uses."""
    from giap.runner.local import PUBLISHED, PUBLISHED_FILES

    for key, rel in s.PUBLISHED.items():
        assert PUBLISHED_FILES.get(key) == rel, (key, PUBLISHED_FILES.get(key))
    det = {k: v for k, v in PUBLISHED_FILES.items() if k.startswith("determinants")}
    assert det == s.PUBLISHED
    assert not any(k.startswith("determinants") for k in PUBLISHED)


def test_caveats_are_the_spec_text_verbatim() -> None:
    """Step 8 / §7.8.2: three strings, verbatim, for the M9 report."""
    from giap.determinants import CAVEATS

    assert tuple(CAVEATS) == s.CAVEATS


# ---------------------------------------------------------------------------- schemas (ex. 6)


@pytest.mark.parametrize(
    ("stem", "columns"),
    [
        ("determinants_model", s.MODEL_COLUMNS),
        ("determinants_attribution", s.ATTRIBUTION_COLUMNS),
    ],
)
def test_schema_fields_are_the_table_columns(stem: str, columns: tuple[str, ...]) -> None:
    assert tuple(s.field_names(stem)) == columns


def test_features_schema_fields_are_the_feature_columns() -> None:
    from giap.determinants.features import FEATURE_COLUMNS

    assert tuple(s.field_names("determinants_features")) == tuple(FEATURE_COLUMNS)


def test_model_schema_enums() -> None:
    """Example 6: ``model`` enum {lightgbm, elastic_net, univariate}; ``metric`` as §9.3."""
    assert set(s.enum_of("determinants_model", "model") or ()) == set(s.MODEL_ORDER)
    assert set(s.enum_of("determinants_model", "metric") or ()) == set(s.METRIC_ORDER)


def test_attribution_schema_types() -> None:
    fields = {f["name"]: f for f in s.descriptor("determinants_attribution")["fields"]}
    assert fields["top_decile"]["type"] == "boolean"
    assert set(s.enum_of("determinants_attribution", "dominant_block") or ()) == set(s.BLOCKS)
    for c in ("abundance_log2", "predicted_log2", "residual", *(f"shap_{b}" for b in s.BLOCKS)):
        assert fields[c]["type"] == "number", c


@pytest.mark.parametrize("key", s.OUTPUTS)
def test_example_3_tables_validate(full: dict[str, Any], key: str) -> None:
    """Example 6: ``frictionless validate`` on what the step wrote."""
    path = full["paths"][key]
    assert path.is_file(), f"{key} not written"
    s.validate(path, s.SCHEMAS[key])


def test_writers_round_trip_na_and_empty_flags(tmp_path: Path) -> None:
    """The three writers (``giap.io.determinants_tables``): NA CI, empty flags and booleans
    written so the table validates; header-only attribution validates too."""
    from giap.io import determinants_tables as w

    tbl = pd.DataFrame(
        {"model": ["lightgbm", "univariate"], "feature_or_block": ["all", "cai"],
         "metric": ["cv_r2_mean", "auroc"], "value": [0.25, math.nan],
         "ci_low": [math.nan, math.nan], "ci_high": [math.nan, math.nan],
         "n": pd.array([120, 10], dtype="Int64"), "n_perm": pd.array([0, 0], dtype="Int64"),
         "flags": ["ci_unestimated", "unestimated;ci_unestimated"]}
    )  # fmt: skip
    p = tmp_path / "determinants_model.tsv"
    w.write_determinants_model(p, tbl)
    assert s.header(p) == list(s.MODEL_COLUMNS)
    got = s.read(p)
    assert got["value"].tolist() == ["0.25", "NA"]
    assert got["ci_low"].tolist() == ["NA", "NA"]
    s.validate(p, "determinants_model")

    a = tmp_path / "determinants_attribution.tsv"
    w.write_determinants_attribution(a, pd.DataFrame(columns=list(s.ATTRIBUTION_COLUMNS)))
    assert a.read_text(encoding="utf-8") == "\t".join(s.ATTRIBUTION_COLUMNS) + "\n"
    s.validate(a, "determinants_attribution")


# ---------------------------------------------------------------------------- example 3


def test_example_3_features_table(full: dict[str, Any], expected: dict[str, Any]) -> None:
    from giap.determinants.features import FEATURE_COLUMNS

    path = full["paths"]["determinants_features"]
    assert s.header(path) == list(FEATURE_COLUMNS)
    got = s.read(path)
    want = expected["features"]
    assert len(want) > 500, f"only {len(want)} eligible genes: a vacuous fixture"
    assert got["locus_tag"].tolist() == want["locus_tag"].astype(str).tolist()


def test_example_3_every_row_kind(model: pd.DataFrame, expected: dict[str, Any]) -> None:
    """CV, null (B from config), ΔR² and univariate rows; ``feature_or_block = all`` for the
    first three metrics; one ΔR² row per block that has a column (no empty block)."""
    n_model = int(np.count_nonzero(expected["rows"]))
    assert n_model >= 100, f"{n_model} model rows: the fixture must clear the minimum"
    present = sorted(set(expected["blocks"].values()))
    for m in s.MODELS:
        for metric in ("cv_r2_mean", "cv_r2_sd", "null_p"):
            r = s.rows(model, model=m, metric=metric)
            assert len(r) == 1, (m, metric, r.to_dict("records"))
            assert r["feature_or_block"].iloc[0] == "all"
            assert int(r["n"].iloc[0]) == n_model, (m, metric)
        d = s.rows(model, model=m, metric="delta_r2")
        assert sorted(d["feature_or_block"]) == present, (m, d["feature_or_block"].tolist())
    assert set(model["model"]) == set(s.MODEL_ORDER)


def test_example_3_rows_sorted(model: pd.DataFrame) -> None:
    keys = s.sort_key(model)
    assert keys == sorted(keys), "rows not sorted by model, metric, feature_or_block"


def test_example_3_null_rows(model: pd.DataFrame, full: dict[str, Any]) -> None:
    """Decision 4: ``null_uncalibrated`` on every null row; ``n_perm`` = the non-NaN null values
    (B = 49 here, all finite); p on the (1 + b)/(1 + B) grid; CI NA."""
    B = full["cfg"].determinants.permutations
    assert B == 49
    for m in s.MODELS:
        r = s.rows(model, model=m, metric="null_p").iloc[0]
        assert "null_uncalibrated" in s.flags(r["flags"]), r.to_dict()
        assert int(r["n_perm"]) == B, f"n_perm {r['n_perm']}: B must be determinants.permutations"
        p = s.num(r["value"])
        assert 1 / (1 + B) - 1e-9 <= p <= 1.0
        b = p * (1 + B) - 1
        assert abs(b - round(b)) < 1e-4, f"p = {p} is not (1 + b)/(1 + B)"
        assert r["ci_low"] == "NA" and r["ci_high"] == "NA"


def test_example_3_cv_r2_mean_ci_unestimated(model: pd.DataFrame) -> None:
    """Decision 4 (ADR-0031 pattern): the null p's effect carries ``ci_unestimated``."""
    for m in s.MODELS:
        r = s.rows(model, model=m, metric="cv_r2_mean").iloc[0]
        assert "ci_unestimated" in s.flags(r["flags"]), r.to_dict()
        assert r["ci_low"] == "NA" and r["ci_high"] == "NA"
        assert math.isfinite(s.num(r["value"]))
        sd = s.rows(model, model=m, metric="cv_r2_sd").iloc[0]
        assert math.isfinite(s.num(sd["value"])) and s.num(sd["value"]) >= 0


def test_example_3_delta_rows_have_t_intervals(model: pd.DataFrame) -> None:
    d = model[model["metric"] == "delta_r2"]
    assert len(d), "no delta_r2 rows"
    for _, r in d.iterrows():
        v, lo, hi = s.num(r["value"]), s.num(r["ci_low"]), s.num(r["ci_high"])
        assert math.isfinite(v) and lo <= v <= hi, r.to_dict()
        assert "ci_unestimated" in s.flags(r["flags"])


def test_example_3_univariate_rows(model: pd.DataFrame, expected: dict[str, Any]) -> None:
    """One ``spearman_rho`` and one ``auroc`` per encoded feature column, each with a CI."""
    cols = sorted(map(str, expected["X"].columns))
    uni = model[model["model"] == "univariate"]
    for metric in ("spearman_rho", "auroc"):
        assert sorted(s.rows(uni, metric=metric)["feature_or_block"]) == cols, metric
    assert set(uni["metric"]) == {"spearman_rho", "auroc"}
    assert all("ci_unestimated" in s.flags(f) for f in uni["flags"]), "ADR-0032"


def test_example_3_univariate_stream(model: pd.DataFrame, full: dict[str, Any],
                                     expected: dict[str, Any]) -> None:  # fmt: skip
    """§6.6: the univariate rows come from ``for_step("determinants.bootstrap", seed)`` with
    ``stats.bootstrap`` draws, on the M8-07 feature frame (no other stream consumed first)."""
    from giap.determinants.models import univariate
    from giap.rng import for_step

    cfg = full["cfg"]
    want = univariate(expected["features"], cfg.stats.bootstrap,
                      for_step("determinants.bootstrap", cfg.seed))  # fmt: skip
    got = model[model["model"] == "univariate"].set_index(["feature_or_block", "metric"])
    bad = []
    for r in want.itertuples(index=False):
        g = got.loc[(r.feature_or_block, r.metric)]
        for c in ("value", "ci_low", "ci_high", "n"):
            if g[c] != s.fmt(getattr(r, c)):
                bad.append((r.feature_or_block, r.metric, c, g[c], s.fmt(getattr(r, c))))
    assert not bad, bad[:5]


def test_example_3_model_flag_vocabulary(model: pd.DataFrame) -> None:
    seen = set().union(*(s.flags(f) for f in model["flags"]))
    assert seen <= s.MODEL_FLAGS, seen - s.MODEL_FLAGS
    assert "too_few_rows" not in seen


def test_example_3_attribution_rows(attribution: pd.DataFrame, expected: dict[str, Any],
                                    full: dict[str, Any]) -> None:  # fmt: skip
    """§9.3 columns in order; one row per model row plus the detected ribosomal genes
    (``cai_reference``); predicted values on every model row; residual = observed - predicted."""
    assert s.header(full["paths"]["determinants_attribution"]) == list(s.ATTRIBUTION_COLUMNS)
    f = expected["features"]
    tags = f["locus_tag"].astype(str).to_numpy()
    model_tags = set(tags[expected["rows"]])
    got = set(attribution["locus_tag"])
    assert model_tags <= got
    extra = got - model_tags
    for _, r in attribution[attribution["locus_tag"].isin(extra)].iterrows():
        assert "cai_reference" in s.flags(r["flags"]), r.to_dict()
    for _, r in attribution[attribution["locus_tag"].isin(model_tags)].iterrows():
        ab, pred, res = s.num(r["abundance_log2"]), s.num(r["predicted_log2"]), s.num(r["residual"])
        assert math.isfinite(pred), r.to_dict()
        assert res == pytest.approx(ab - pred, abs=1e-4 * max(1.0, abs(ab))), r.to_dict()


def test_example_3_attribution_flags_and_top_decile(attribution: pd.DataFrame,
                                                    expected: dict[str, Any]) -> None:  # fmt: skip
    """§7.8.2 step 6 / §9.4: dominant block only on the top decile; the residual flag on
    |residual| > 2 SD; flags from the controlled vocabulary; SHAP NA only for an empty block."""
    m = attribution[~attribution["flags"].map(lambda c: "cai_reference" in s.flags(c))]
    assert set(attribution["top_decile"]) <= {"true", "false"}
    top = m["top_decile"] == "true"
    assert 0.05 * len(m) <= top.sum() <= 0.15 * len(m), int(top.sum())
    assert (attribution.loc[attribution["top_decile"] == "false", "dominant_block"] == "NA").all()
    assert set(m.loc[top, "dominant_block"]) - {"NA"} <= set(s.BLOCKS)
    resid = m["residual"].map(s.num).to_numpy(float)
    sd = float(np.std(resid, ddof=1))
    flagged = m["flags"].map(lambda c: "abundance_not_explained_by_sequence" in s.flags(c))
    far = np.abs(resid) > 2 * sd
    clear = np.abs(np.abs(resid) - 2 * sd) > 1e-3 * sd  # %.6g rounding at the boundary
    assert (flagged.to_numpy()[clear] == far[clear]).all()
    assert flagged.any(), "no residual > 2 SD on 500+ genes: a vacuous check"
    seen = set().union(*(s.flags(c) for c in attribution["flags"]))
    assert seen <= s.ATTRIBUTION_FLAGS, seen - s.ATTRIBUTION_FLAGS
    present = set(expected["blocks"].values())
    for b in s.BLOCKS:
        col = attribution[f"shap_{b}"]
        if b in present:
            assert (col != "NA").any(), b
        else:
            assert (col == "NA").all(), f"shap_{b} filled for a block with no column"


def test_example_3_events(full: dict[str, Any]) -> None:
    """Step 8: the caveats and the dropped features go to the step's events for M9."""
    caveats = _events(full["records"], s.CAVEATS_EVENT)
    assert caveats == [{"caveats": list(s.CAVEATS)}] or (
        len(caveats) == 1 and list(caveats[0]["caveats"]) == list(s.CAVEATS)
    ), caveats
    dropped = _events(full["records"], s.DROPPED_EVENT)
    assert len(dropped) == 1, dropped
    names = dropped[0]["features"]
    assert isinstance(names, list) and all(isinstance(n, str) for n in names)
    tsv = s.read(full["paths"]["determinants_features"])
    from giap.determinants.features import FEATURE_BLOCKS

    all_na = [f for fs in FEATURE_BLOCKS.values() for f in fs if (tsv[f] == "NA").all()]
    assert set(all_na) <= set(names), (all_na, names)


# ---------------------------------------------------------------------------- example 4


def test_example_4_threads_do_not_change_the_tables(
    work: Path, full: dict[str, Any], tmp_path: Path
) -> None:
    """Same seed, threads 1 vs 8 → byte-identical tables (§6.6, fixed LightGBM ``num_threads``
    per fit; the null's permutations are collected in b order)."""
    cfg = s.config(SEED, threads=8)
    paths = s.run(cfg, s.inputs(work), tmp_path / "t8")
    for key in s.OUTPUTS:
        assert paths[key].read_bytes() == full["paths"][key].read_bytes(), key


# ---------------------------------------------------------------------------- below the minimum


@pytest.fixture(scope="module")
def few(work: Path, step_mod: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> dict:
    """Example 1's case (sim100k: < 100 model rows): abundance kept on 60 detected genes."""
    root = tmp_path_factory.mktemp("m810_few")
    gt = s.gene_table_variant(work, root / "genes.parquet", keep_abundance=60)
    paths = s.run(s.config(SEED), s.inputs(work, gene_table=gt), root / "OUT")
    return {"paths": paths, "model": s.read(paths["determinants_model"]),
            "attribution": s.read(paths["determinants_attribution"])}  # fmt: skip


def test_below_minimum_one_too_few_rows_row_per_model(few: dict[str, Any]) -> None:
    tbl = few["model"]
    for m in s.MODELS:
        r = s.rows(tbl, model=m)
        assert len(r) == 1, (m, r.to_dict("records"))
        row = r.iloc[0]
        assert row["metric"] == "cv_r2_mean" and row["feature_or_block"] == "all"
        assert row["value"] == "NA"
        assert {"too_few_rows", "ci_unestimated"} <= s.flags(row["flags"]), row.to_dict()
    uni = tbl[tbl["model"] == "univariate"]
    assert set(uni["metric"]) == {"spearman_rho", "auroc"}
    assert s.sort_key(tbl) == sorted(s.sort_key(tbl))


def test_below_minimum_attribution_has_na_predictions(few: dict[str, Any]) -> None:
    a = few["attribution"]
    assert len(a) == 60, len(a)
    assert (a["predicted_log2"] == "NA").all() and (a["residual"] == "NA").all()
    assert all("too_few_rows" in s.flags(f) for f in a["flags"])
    for b in s.BLOCKS:
        assert (a[f"shap_{b}"] == "NA").all()
    assert (a["dominant_block"] == "NA").all()


@pytest.mark.parametrize("key", s.OUTPUTS)
def test_below_minimum_tables_validate(few: dict[str, Any], key: str) -> None:
    s.validate(few["paths"][key], s.SCHEMAS[key])


# ---------------------------------------------------------------------------- example 5


@pytest.fixture(scope="module")
def no_abundance(work: Path, step_mod: ModuleType,
                 tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:  # fmt: skip
    root = tmp_path_factory.mktemp("m810_noab")
    gt = s.gene_table_variant(work, root / "genes.parquet", keep_abundance=None)
    ins = s.inputs(work, gene_table=gt)
    paths, records = _run_captured(s.config(SEED), ins, root / "OUT")
    return {"ins": ins, "paths": paths, "records": records}


def test_example_5_features_and_auroc_rows_only(no_abundance: dict[str, Any]) -> None:
    from giap.determinants.features import FEATURE_COLUMNS, model_matrix

    p = no_abundance["paths"]
    feats = s.read(p["determinants_features"])
    want = s.expected_features(no_abundance["ins"])
    assert s.header(p["determinants_features"]) == list(FEATURE_COLUMNS)
    assert feats["locus_tag"].tolist() == want["locus_tag"].astype(str).tolist()
    tbl = s.read(p["determinants_model"])
    assert len(tbl), "no AUROC rows"
    assert set(tbl["model"]) == {"univariate"} and set(tbl["metric"]) == {"auroc"}
    cols = sorted(map(str, model_matrix(want)[0].columns))
    assert sorted(tbl["feature_or_block"]) == cols
    assert s.sort_key(tbl) == sorted(s.sort_key(tbl))


def test_example_5_attribution_header_only(no_abundance: dict[str, Any]) -> None:
    text = no_abundance["paths"]["determinants_attribution"].read_text(encoding="utf-8")
    assert text == "\t".join(s.ATTRIBUTION_COLUMNS) + "\n"


def test_example_5_run_flag(no_abundance: dict[str, Any]) -> None:
    msgs = [r.getMessage() for r in no_abundance["records"]]
    assert any(s.RUN_FLAG in m for m in msgs), msgs


@pytest.mark.parametrize("key", s.OUTPUTS)
def test_example_5_tables_validate(no_abundance: dict[str, Any], key: str) -> None:
    s.validate(no_abundance["paths"][key], s.SCHEMAS[key])


# ---------------------------------------------------------------------------- optional inputs


def test_context_skipped_still_runs(work: Path, step_mod: ModuleType, tmp_path: Path) -> None:
    """No C7/C8 and no gene context (CONTEXT skipped): the step runs; ``replichore_position``,
    ``in_prophage`` and ``in_bgc`` are NA."""
    gt = s.gene_table_variant(work, tmp_path / "genes.parquet", keep_abundance=None)
    paths = s.run(s.config(SEED), s.inputs(work, gene_table=gt, optional=False), tmp_path / "O")
    feats = s.read(paths["determinants_features"])
    assert len(feats) > 500
    for c in ("replichore_position", "in_prophage", "in_bgc"):
        assert (feats[c] == "NA").all(), c


# ---------------------------------------------------------------------------- --publication


def test_publication_sets_the_null_to_10000_permutations() -> None:
    """Lead decision 3: ``--publication`` follows §8.2 (B = 10,000) for the determinants null,
    whose B is ``determinants.permutations``; ``--set`` still wins."""
    from giap.config import resolve_config

    assert resolve_config().determinants.permutations == 1000
    assert resolve_config(publication=True).determinants.permutations == 10_000
    cfg = resolve_config(publication=True, overrides=["determinants.permutations=49"])
    assert cfg.determinants.permutations == 49


def test_publication_help_states_the_runtime() -> None:
    """Degenerate cases: "the run takes ~10x longer, stated in the CLI help"."""
    from click.testing import CliRunner

    from giap.cli import cli

    res = CliRunner().invoke(cli, ["run", "--help"])
    assert res.exit_code == 0
    text = " ".join(res.output.split())
    m = re.search(r"--publication(.*?)(?=\s--[a-z]|$)", text)
    assert m, text
    assert "longer" in m.group(1).lower(), m.group(1)


# ---------------------------------------------------------------------------- giap _step


def test_giap_step_writes_events(
    work: Path, step_mod: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Nextflow's path (``giap _step determinants``): the caveats event reaches
    ``logs/events.jsonl`` (§6.7) for the M9 report."""
    from click.testing import CliRunner

    from giap.cli import cli

    gt = s.gene_table_variant(work, tmp_path / "in" / "genes.parquet", keep_abundance=None)
    monkeypatch.chdir(tmp_path)
    args = ["_step", "determinants", "--set", "stats.bootstrap=19"]
    for k, p in s.inputs(work, gene_table=gt).items():
        args += ["--in", f"{k}={p}"]
    for k in s.OUTPUTS:
        args += ["--out", f"{k}={s.PUBLISHED[k]}"]
    res = CliRunner().invoke(cli, args)
    assert res.exit_code == 0, res.output
    for k in s.OUTPUTS:
        assert (tmp_path / s.PUBLISHED[k]).is_file(), k
    lines = [json.loads(x) for x in (tmp_path / "logs" / "events.jsonl").read_text().splitlines()]
    ev = [e for e in lines if e["event"] == s.CAVEATS_EVENT]
    assert len(ev) == 1 and ev[0]["step"] == "determinants"
    assert ev[0]["data"]["caveats"] == list(s.CAVEATS)


# ---------------------------------------------------------------------------- Nextflow config


def test_base_config_gives_determinants_process_high_resources() -> None:
    """Lead decision 2: ``withName`` at ``process_high`` (16 CPU, 48 GB x attempt) in
    workflow/conf/base.config, so the null's permutations run in parallel.

    The spec's ``withName: 'DETERMINANTS'`` cannot be literal: a GIAP_STEP alias named
    DETERMINANTS inside the DETERMINANTS subworkflow is a Nextflow error ("`DETERMINANTS` is
    already included"). So the selector must match the process determinants.nf runs the step
    under -- its GIAP_STEP alias, or GIAP_STEP itself -- by simple or fully qualified name
    (``DETERMINANTS:<process>``), and must not match another py step's process."""
    nf = REPO / "workflow" / "subworkflows" / "local" / "determinants.nf"
    assert nf.is_file(), "workflow/subworkflows/local/determinants.nf missing"
    aliases = re.findall(r"GIAP_STEP\s+as\s+(\w+)", nf.read_text(encoding="utf-8"))
    assert len(aliases) <= 1, f"one GIAP_STEP process expected in determinants.nf: {aliases}"
    proc = aliases[0] if aliases else "GIAP_STEP"
    conf = (REPO / "workflow" / "conf" / "base.config").read_text(encoding="utf-8")
    sel = re.findall(r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{([^}]*(?:\}[^}]*)?)\}", conf)

    def matches(pat: str, simple: str, qualified: str) -> bool:
        return bool(re.fullmatch(pat, simple) or re.fullmatch(pat, qualified))

    hits = [(pat, body) for pat, body in sel if matches(pat, proc, f"DETERMINANTS:{proc}")]
    assert hits, (
        f"no withName selector for DETERMINANTS:{proc} in base.config: {[p for p, _ in sel]}"
    )
    for pat, _ in hits:
        for other in ("CONTEXT_STATS", "GIAP_STEP"):
            assert not matches(pat, other, f"CONTEXT:{other}"), f"{pat} also matches {other}"
    body = hits[0][1]
    assert re.search(r"cpus\s*=\s*16\b", body), body
    assert re.search(r"memory\s*=\s*\{\s*48\.GB\s*\*\s*task\.attempt\s*\}", body), body
