"""M8-10 writers: any model table the step can build writes a schema-valid
``determinants_model.tsv`` with its rows in order, and any attribution table a schema-valid
``determinants_attribution.tsv`` (docs/specs/M8-10-determinants-step-and-run.md example 6;
GIPPA_SPEC.md §9.2, §9.3, §9.5). Contract choices: tests/helpers/determinants_step_data.py.
"""

from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import determinants_step_data as s

finite = st.floats(min_value=-1e6, max_value=1e6, allow_nan=False, allow_infinity=False)
maybe = st.one_of(st.just(math.nan), finite)
model_flags = st.lists(st.sampled_from(sorted(s.MODEL_FLAGS)), unique=True, max_size=3)
attr_flags = st.lists(st.sampled_from(sorted(s.ATTRIBUTION_FLAGS)), unique=True, max_size=3)


@st.composite
def model_rows(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(0, 12))
    rows = []
    for _ in range(n):
        rows.append({
            "model": draw(st.sampled_from(s.MODEL_ORDER)),
            "feature_or_block": draw(st.sampled_from(["all", *s.BLOCKS, "cai", "tir_ostir"])),
            "metric": draw(st.sampled_from(s.METRIC_ORDER)),
            "value": draw(maybe), "ci_low": draw(maybe), "ci_high": draw(maybe),
            "n": draw(st.one_of(st.none(), st.integers(0, 5000))),
            "n_perm": draw(st.integers(0, 10_000)),
            "flags": ";".join(draw(model_flags)),
        })  # fmt: skip
    df = pd.DataFrame(rows, columns=list(s.MODEL_COLUMNS))
    df[["value", "ci_low", "ci_high"]] = df[["value", "ci_low", "ci_high"]].astype(float)
    df[["n", "n_perm"]] = df[["n", "n_perm"]].astype("Int64")
    return df


@st.composite
def attribution_rows(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(0, 12))
    rows = []
    for i in range(n):
        r = {"locus_tag": f"G{i:04d}", "abundance_log2": draw(maybe),
             "predicted_log2": draw(maybe), "residual": draw(maybe)}  # fmt: skip
        for b in s.BLOCKS:
            r[f"shap_{b}"] = draw(maybe)
        r["dominant_block"] = draw(st.one_of(st.none(), st.sampled_from(s.BLOCKS)))
        r["top_decile"] = draw(st.booleans())
        r["flags"] = ";".join(draw(attr_flags))
        rows.append(r)
    df = pd.DataFrame(rows, columns=list(s.ATTRIBUTION_COLUMNS))
    df["dominant_block"] = df["dominant_block"].astype("string")
    df["top_decile"] = df["top_decile"].astype(bool)
    return df


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(model_rows())
def test_model_writer_valid_and_order_kept(df: pd.DataFrame) -> None:
    from gippa.io.determinants_tables import write_determinants_model

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "determinants_model.tsv"
        write_determinants_model(p, df)
        assert s.header(p) == list(s.MODEL_COLUMNS)
        got = s.read(p)
        assert len(got) == len(df)
        assert got["model"].tolist() == df["model"].tolist()
        assert got["metric"].tolist() == df["metric"].tolist()
        assert [s.flags(c) for c in got["flags"]] == [s.flags(c) for c in df["flags"]]
        s.validate(p, "determinants_model")


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(attribution_rows())
def test_attribution_writer_valid(df: pd.DataFrame) -> None:
    from gippa.io.determinants_tables import write_determinants_attribution

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "determinants_attribution.tsv"
        write_determinants_attribution(p, df)
        assert s.header(p) == list(s.ATTRIBUTION_COLUMNS)
        got = s.read(p)
        assert got["locus_tag"].tolist() == df["locus_tag"].tolist()
        assert got["top_decile"].tolist() == ["true" if v else "false" for v in df["top_decile"]]
        s.validate(p, "determinants_attribution")
