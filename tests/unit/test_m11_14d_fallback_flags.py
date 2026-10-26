# ruff: noqa: RUF001, RUF002  (the report's wording is Unicode)
"""M11-14d, ADR-0061 decision 1 (proposed): ``tu_source_distance_model`` on the TU-blocked rows.

Kept in its own file so that it can be dropped if Khoa rejects ADR-0061
(docs/adr/ADR-0061-distance-model-fallback-scope.md on branch m11/14s1-operon-fallback-attack;
14s1 outcome rule in docs/specs/M11-14-container-install.md: "the affected rows carry a flag
(``tu_source_distance_model``) with an R4 caveat"). When ``operon_source = distance_model``:

- every ``tables/autocorrelation.tsv`` row with ``statistic = morans_i`` and
  ``null_model = tu_block`` (all lags) carries the flag; the report states "Moran's I under the
  TU-block null is anti-conservative at short lags when TUs come from the distance model;
  rejection at α = 0.05 was 0.34 in simulation";
- every ``tables/pairs.tsv`` row with ``test = both_detected_t`` carries the flag; the report
  states "conservative under the distance-model fallback (merged TU boundaries widen the
  within-operon fit); fewer discordant pairs are found".

No other row changes, and nothing changes under ``uniop`` or ``precomputed``.

Interface pinned here (the ADR leaves it to the 14d implementer):

- ``giap.operons.flags.TU_SOURCE_FLAG == "tu_source_distance_model"``.
- ``giap.operons.flags.add_tu_source_flag(frame, table, operon_source) -> pandas.DataFrame``, pure:
  ``table`` is ``"autocorrelation"`` or ``"pairs"``; the frame's ``flags`` cells are token
  sequences (as ``autocorrelation_rows`` and ``merge_tests`` build them). It returns a new frame,
  appends the flag to the target rows only (existing tokens kept, no duplicate) and leaves the
  input frame untouched.
- The steps ``context_stats`` and ``operons`` take ``operon_meta`` (``operon_meta.json`` of
  ``operon_pairs``) as an optional input and apply the rule from its ``operon_source``.
- The table schemas' ``flags`` patterns accept the token.
- The report prints each caveat (Genome context and Operons sections) when the manifest's
  ``operon_source`` is ``distance_model`` and the table is present.
"""

from __future__ import annotations

import json
import logging
import math
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

FLAG = "tu_source_distance_model"
REPO = Path(__file__).resolve().parents[2]
TABLE_SCHEMAS = REPO / "src" / "giap" / "schemas" / "tables"
MORAN_CAVEAT = (
    "Moran's I under the TU-block null is anti-conservative at short lags when TUs come from the "
    "distance model; rejection at α = 0.05 was 0.34 in simulation"
)
DISCORDANCE_CAVEAT = (
    "conservative under the distance-model fallback (merged TU boundaries widen the "
    "within-operon fit); fewer discordant pairs are found"
)
SOURCES = ["uniop", "distance_model", "precomputed"]


def _add() -> Any:
    from giap.operons import flags

    fn = getattr(flags, "add_tu_source_flag", None)
    assert fn is not None, "giap.operons.flags.add_tu_source_flag is not implemented"
    return fn


def _tokens(cell: Any) -> list[str]:
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    if isinstance(cell, str):
        return [t for t in cell.split(";") if t]
    return [str(t) for t in cell]


def test_flag_constant() -> None:
    from giap.operons import flags

    assert getattr(flags, "TU_SOURCE_FLAG", None) == FLAG


# ---------------------------------------------------------------- pure rule (examples)


def _autocorr(rows: Sequence[tuple[str, str, int, Sequence[str]]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "contig": ["chr1"] * len(rows),
            "statistic": [r[0] for r in rows],
            "null_model": [r[1] for r in rows],
            "lag": pd.array([r[2] for r in rows], dtype="Int64"),
            "value": np.linspace(0.1, 0.2, len(rows)),
            "flags": pd.Series([tuple(r[3]) for r in rows], dtype=object),
        }
    )


AUTOCORR_ROWS = [
    ("morans_i", "tu_block", 1, ()),
    ("morans_i", "tu_block", 2, ("ci_unestimated",)),
    ("morans_i", "tu_block", 20, ()),
    ("morans_i", "gene_label", 1, ()),
    ("join_count_bb", "tu_block", 1, ()),
    ("oriter_slope", "bootstrap", None, ("descriptive",)),
]


def test_autocorrelation_example_distance_model() -> None:
    df = _autocorr(AUTOCORR_ROWS)
    got = _add()(df, "autocorrelation", "distance_model")
    want = [True, True, True, False, False, False]
    assert [FLAG in _tokens(c) for c in got["flags"]] == want
    assert _tokens(got["flags"].iloc[1]) == ["ci_unestimated", FLAG] or set(
        _tokens(got["flags"].iloc[1])
    ) == {"ci_unestimated", FLAG}
    assert _tokens(got["flags"].iloc[5]) == ["descriptive"]


def _pairs(tests: Sequence[str], flags: Sequence[Sequence[str]] | None = None) -> pd.DataFrame:
    n = len(tests)
    return pd.DataFrame(
        {
            "pair_id": [f"PAIR_{i + 1}" for i in range(n)],
            "class": ["OP"] * n,
            "test": list(tests),
            "p_value": np.linspace(0.01, 0.9, n),
            "flags": pd.Series([list(f) for f in (flags or [[]] * n)], dtype=object),
        }
    )


def test_pairs_example_distance_model() -> None:
    df = _pairs(["both_detected_t", "one_detected_pmiss", "missed_operon_decoy", "none",
                 "both_detected_t"], [[], [], [], [], ["rna_between"]])  # fmt: skip
    got = _add()(df, "pairs", "distance_model")
    assert [FLAG in _tokens(c) for c in got["flags"]] == [True, False, False, False, True]
    assert set(_tokens(got["flags"].iloc[4])) == {"rna_between", FLAG}


@pytest.mark.parametrize("source", ["uniop", "precomputed", None])
@pytest.mark.parametrize("table", ["autocorrelation", "pairs"])
def test_no_flag_unless_distance_model(source: str | None, table: str) -> None:
    df = _autocorr(AUTOCORR_ROWS) if table == "autocorrelation" else _pairs(["both_detected_t"])
    got = _add()(df, table, source)
    assert not any(FLAG in _tokens(c) for c in got["flags"])
    assert [_tokens(c) for c in got["flags"]] == [_tokens(c) for c in df["flags"]]


def test_empty_frames() -> None:
    for table, df in (("autocorrelation", _autocorr([])), ("pairs", _pairs([]))):
        got = _add()(df, table, "distance_model")
        assert len(got) == 0 and list(got.columns) == list(df.columns)


def test_unknown_table_raises() -> None:
    with pytest.raises((ValueError, KeyError)):
        _add()(_pairs(["both_detected_t"]), "islands", "distance_model")


# ---------------------------------------------------------------- pure rule (properties)

_STAT = st.sampled_from(["morans_i", "join_count_bb", "oriter_slope"])
_NULL = st.sampled_from(["gene_label", "tu_block", "bootstrap"])
_TOK = st.sampled_from(["ci_unestimated", "not_estimable", "descriptive", "few_values"])
_AC_ROW = st.tuples(_STAT, _NULL, st.integers(1, 30), st.lists(_TOK, unique=True, max_size=3))
_TEST = st.sampled_from(["both_detected_t", "one_detected_pmiss", "missed_operon_decoy", "none"])
_PTOK = st.sampled_from(["rna_between", "abundance_imputed", "ci_unestimated", FLAG])


def _target(table: str, df: pd.DataFrame) -> list[bool]:
    if table == "autocorrelation":
        return list((df["statistic"] == "morans_i") & (df["null_model"] == "tu_block"))
    return list(df["test"] == "both_detected_t")


def _check(table: str, df: pd.DataFrame, source: str) -> None:
    before = df.copy(deep=True)
    got = _add()(df, table, source)
    pd.testing.assert_frame_equal(df, before)  # input untouched
    assert list(got.columns) == list(df.columns) and len(got) == len(df)
    pd.testing.assert_frame_equal(got.drop(columns="flags"), df.drop(columns="flags"))
    hit = _target(table, df)
    for old, new, t in zip(df["flags"], got["flags"], hit, strict=True):
        o, n = _tokens(old), _tokens(new)
        assert len(n) == len(set(n)), n  # no duplicate token
        if source == "distance_model" and t:
            assert set(n) == set(o) | {FLAG}, (o, n)
        else:
            assert n == o, (o, n)
    again = _add()(got, table, source)  # idempotent
    assert [_tokens(c) for c in again["flags"]] == [_tokens(c) for c in got["flags"]]


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=st.lists(_AC_ROW, max_size=25), source=st.sampled_from(SOURCES))
def test_property_autocorrelation(rows: list[Any], source: str) -> None:
    _check("autocorrelation", _autocorr(rows), source)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    rows=st.lists(st.tuples(_TEST, st.lists(_PTOK, unique=True, max_size=3)), max_size=25),
    source=st.sampled_from(SOURCES),
)
def test_property_pairs(rows: list[Any], source: str) -> None:
    _check("pairs", _pairs([r[0] for r in rows], [r[1] for r in rows]), source)


# ---------------------------------------------------------------- schemas


@pytest.mark.parametrize(("stem", "other"), [("autocorrelation", "ci_unestimated"),
                                             ("pairs", "rna_between")])  # fmt: skip
def test_schema_flags_pattern_accepts_the_token(stem: str, other: str) -> None:
    doc = json.loads((TABLE_SCHEMAS / f"{stem}.schema.json").read_text(encoding="utf-8"))
    field = next(f for f in doc["fields"] if f["name"] == "flags")
    pattern = field["constraints"]["pattern"]
    for value in (FLAG, f"{other};{FLAG}", f"{FLAG};{other}"):
        assert re.fullmatch(pattern, value), (stem, value)
    assert not re.fullmatch(pattern, "tu_source_uniop")


# ---------------------------------------------------------------- steps declare operon_meta


@pytest.mark.parametrize("name", ["operons", "context_stats"])
def test_steps_take_operon_meta_as_optional_input(name: str) -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get(name)
    assert "operon_meta" in spec.optional_inputs, spec.optional_inputs
    assert "operon_meta" not in spec.inputs


# ---------------------------------------------------------------- operons step end to end


def _meta(path: Path, source: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"operon_source": source, "reason": "test"}), encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from helpers.layers_chain import build_chain

    return build_chain(Path(tmp_path_factory.mktemp("ops_chain_14d")))


def _gene_table(chain: Any, path: Path) -> Path:
    """A §9.3 gene table with abundance (as ``tests/unit/steps/test_operons_step.py``)."""
    from giap.genetable.assemble import GENES_COLUMNS, GENES_TYPES
    from giap.io.genes_tsv import write_genes_parquet

    c2 = chain.frame("genes")
    n = len(c2)
    default = {"list": lambda: [], "bool": lambda: False, "int": lambda: None,
               "float": lambda: math.nan, "str": lambda: None}  # fmt: skip
    df = pd.DataFrame({c: [default[GENES_TYPES[c]]() for _ in range(n)] for c in GENES_COLUMNS})
    for c in ("locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand",
              "partial", "crosses_origin", "hypothetical"):  # fmt: skip
        df[c] = c2[c].to_numpy()
    idx = np.arange(n)
    det = idx % 4 != 3
    df["detected"] = det
    df["detection_status"] = np.where(det, "detected", "not_detected")
    df["p_detect"] = np.linspace(0.2, 0.9, n)
    df["abundance_log2"] = np.where(det, 12.0 + np.sin(idx), np.nan)
    df["abundance_type"] = "ibaq"
    df["abundance_scale"] = "absolute"
    df["abundance_imputed"] = False
    path.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(path, df)
    return path


def _run_operons(chain: Any, root: Path, source: str | None) -> pd.DataFrame:
    from giap.config import GiapConfig
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    keys = ("genes", "replicons", "tus", "pairs_raw", "features", "promoter_windows",
            "terminators", "tx_rates")  # fmt: skip
    ins = {k: chain.inputs[k] for k in keys}
    ins["gene_table"] = _gene_table(chain, root / "genetable" / "genes.parquet")
    if source is not None:
        ins["operon_meta"] = _meta(root / "operon_meta.json", source)
    outs = {k: root / "tables" / f"{k}.tsv" for k in ("pairs", "operons", "operon_flags")}
    logging.getLogger("giap").setLevel(logging.WARNING)
    reg.run_step("operons", GiapConfig(seed=1), ins, outs)
    return pd.read_csv(outs["pairs"], sep="\t", dtype=str, keep_default_na=False)


@pytest.fixture(scope="module")
def pairs_by_source(chain: Any, tmp_path_factory: pytest.TempPathFactory) -> dict[Any, Any]:
    return {s: _run_operons(chain, tmp_path_factory.mktemp(f"ops_{s}"), s)
            for s in ("distance_model", "uniop", None)}  # fmt: skip


def test_operons_step_flags_every_both_detected_t_row(pairs_by_source: dict[Any, Any]) -> None:
    df = pairs_by_source["distance_model"]
    t = df["test"] == "both_detected_t"
    assert t.sum() >= 1, "fixture has no both_detected_t row"
    has = df["flags"].map(lambda v: FLAG in _tokens(v))
    assert has[t].all(), df.loc[t & ~has, ["pair_id", "flags"]]
    assert not has[~t].any(), df.loc[~t & has, ["pair_id", "test", "flags"]]


@pytest.mark.parametrize("source", ["uniop", None])
def test_operons_step_no_flag_otherwise(pairs_by_source: dict[Any, Any], source: Any) -> None:
    df = pairs_by_source[source]
    assert not df["flags"].map(lambda v: FLAG in _tokens(v)).any()


def test_operons_step_only_the_flags_column_differs(pairs_by_source: dict[Any, Any]) -> None:
    a, b = pairs_by_source["distance_model"], pairs_by_source["uniop"]
    pd.testing.assert_frame_equal(a.drop(columns="flags"), b.drop(columns="flags"))
    strip = b["flags"].map(_tokens)
    assert [sorted(set(x) - {FLAG}) for x in a["flags"].map(_tokens)] == [sorted(x) for x in strip]


# ---------------------------------------------------------------- context_stats step end to end

SEED = 20261022
CTX_OUT = {
    "regions_stats": "tables/regions_stats.tsv",
    "detection_summary": "tables/detection_summary.tsv",
    "autocorrelation": "tables/autocorrelation.tsv", "islands": "tables/islands.tsv",
    "genes_tsv": "tables/genes.tsv", "gene_context": "context/gene_context.parquet",
    "island_diagnostics": "context/island_diagnostics.parquet",
}  # fmt: skip


@pytest.fixture(scope="module")
def ctx_work(tmp_path_factory: pytest.TempPathFactory) -> Path:
    sys.path.insert(0, str(REPO / "tests" / "helpers"))
    import determinants_step_data as s

    return Path(s.build_work(tmp_path_factory.mktemp("ctx_14d"), SEED))


def _ctx_inputs(work: Path) -> dict[str, Path]:
    """What ``run_context_stats`` passes (as ``tests/unit/steps/test_support_context.py``)."""
    from giap.io.work_tables import OPTIONAL_KEYS, default_paths
    from giap.steps.gene_table import contract_paths

    ins = contract_paths(work)
    ins |= {k: p for k, p in default_paths(work).items()
            if k in OPTIONAL_KEYS and k not in ins and p.exists()}  # fmt: skip
    ins["gene_table"] = work / "genetable" / "genes.parquet"
    parts = work / "genetable" / "detection_summary_parts"
    for key, name in (("detectability_summary", "detectability"), ("tier_summary", "tier")):
        if (parts / f"{name}.parquet").exists():
            ins[key] = parts / f"{name}.parquet"
    if (work / "features.parquet").exists():
        ins["features"] = work / "features.parquet"
    return ins


@pytest.fixture(scope="module")
def autocorr_by_source(ctx_work: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[Any, Any]:
    sys.path.insert(0, str(REPO / "tests" / "helpers"))
    import m4_chain as mc

    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    cfg = mc.config(SEED, stats={"permutations": 99, "bootstrap": 100})
    out: dict[Any, Any] = {}
    for source in ("distance_model", "uniop"):
        root = tmp_path_factory.mktemp(f"ctx_{source}")
        ins = _ctx_inputs(ctx_work) | {"operon_meta": _meta(root / "operon_meta.json", source)}
        outs = {k: root / v for k, v in CTX_OUT.items()}
        reg.run_step("context_stats", cfg, ins, outs)
        out[source] = pd.read_csv(outs["autocorrelation"], sep="\t", dtype=str,
                                  keep_default_na=False)  # fmt: skip
    return out


def test_context_stats_flags_every_tu_block_morans_i_row(
    autocorr_by_source: dict[Any, Any],
) -> None:
    df = autocorr_by_source["distance_model"]
    t = (df["statistic"] == "morans_i") & (df["null_model"] == "tu_block")
    assert t.sum() >= 2, "fixture has fewer than two TU-block Moran's I lags"
    has = df["flags"].map(lambda v: FLAG in _tokens(v))
    assert has[t].all(), df.loc[t & ~has, ["lag", "flags"]]
    assert not has[~t].any(), df.loc[~t & has, ["statistic", "null_model", "lag"]]


def test_context_stats_no_flag_under_uniop(autocorr_by_source: dict[Any, Any]) -> None:
    df = autocorr_by_source["uniop"]
    assert not df["flags"].map(lambda v: FLAG in _tokens(v)).any()
    a = autocorr_by_source["distance_model"]
    pd.testing.assert_frame_equal(a.drop(columns="flags"), df.drop(columns="flags"))


# ---------------------------------------------------------------- report caveats (R4)


def _report(tmp_path: Path, source: str) -> Any:
    from helpers.report_content import manifest, page, write_manifest
    from helpers.report_fixtures import run_out
    from helpers.viz_out_tree import write_schema_table

    out = run_out(tmp_path)
    flags = FLAG if source == "distance_model" else None
    ac = [{"contig": "chr1", "statistic": "morans_i", "lag": lag, "null_model": "tu_block",
           "value": 0.1, "stat": "morans_i", "p_value": 0.01, "q_value": 0.02,
           "method": "tu_block_permutation", "n_perm": 999, "flags": flags}
          for lag in (1, 2)]  # fmt: skip
    write_schema_table(out / "tables" / "autocorrelation.tsv", "autocorrelation", ac)
    pairs = [{"pair_id": "PAIR_1", "class": "OP", "gene_a": "g1", "gene_b": "g2",
              "strand": "+", "intergenic_bp": 20, "test": "both_detected_t", "p_value": 0.01,
              "q_value": 0.02, "flags": flags}]  # fmt: skip
    write_schema_table(out / "tables" / "pairs.tsv", "pairs", pairs)
    write_manifest(out, manifest(operon_source=source, operon_source_reason="UniOP not installed"))
    return page(out)


def test_report_caveats_under_distance_model(tmp_path: Path) -> None:
    from helpers.report_content import section_text

    p = _report(tmp_path, "distance_model")
    assert MORAN_CAVEAT in section_text(p, "genome-context")
    assert DISCORDANCE_CAVEAT in section_text(p, "operons")


def test_report_no_caveats_under_uniop(tmp_path: Path) -> None:
    from helpers.report_content import words

    text = words(_report(tmp_path, "uniop").root.text())
    assert MORAN_CAVEAT not in text and DISCORDANCE_CAVEAT not in text
