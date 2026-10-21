"""Hypothesis properties of the M5-09 table writers and readers (GIAP_SPEC.md §9.2, §6.5, §11.2;
docs/specs/M5-09-context-tables-step.md acceptance "round trips (dtype-aware) of both tables").

Invariants, over random tables:

* ``read(write(df))`` returns the same values (floats to the six significant digits ADR-0006's
  ``%.6g`` keeps) in the same row order, for both tables;
* ``islands.tsv`` coordinates are exactly ``start + 1`` and ``end``, for every span, including
  the origin crossers whose ``end`` exceeds the replicon length (§6.5);
* flags and the island gene list survive as sets / sequences, never as a pandas ``Flags`` object;
* writing is a pure function of the frame (R7): the same frame gives the same bytes.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import context_tables_data as cd

from giap.context.autocorrelation import AUTOCORR_COLUMNS
from giap.context.islands import ISLAND_COLUMNS

FAST = settings(max_examples=50, deadline=None,
                suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
NA = float("nan")


def _vocabulary() -> list[str]:
    from giap.context.vocab import CONTEXT_FLAGS

    return sorted(str(f) for f in CONTEXT_FLAGS)


SMALL_FLOATS = st.one_of(
    st.just(NA),
    st.floats(min_value=-1e6, max_value=1e6, allow_nan=False, allow_infinity=False),
)
PROBABILITIES = st.one_of(st.just(NA), st.floats(min_value=0.0, max_value=1.0))
TAGS = st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789", min_size=1, max_size=8)


#: deferred so importing this module does not import the (not yet written) vocabulary
VOCABULARY = st.deferred(lambda: st.sampled_from(_vocabulary()))


@st.composite
def flags(draw: st.DrawFn) -> tuple[str, ...]:
    chosen = draw(st.lists(VOCABULARY, max_size=3, unique=True))
    return tuple(sorted(chosen))


@st.composite
def autocorrelation_frame(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(min_value=0, max_value=6))
    rows = []
    for i in range(n):
        statistic = draw(st.sampled_from(["morans_i", "join_count_bb", "oriter_slope"]))
        oriter = statistic == "oriter_slope"
        p = NA if oriter else draw(PROBABILITIES)
        rows.append(
            {
                "contig": draw(st.sampled_from(["chr1", "p1", "ctg1"])),
                "statistic": statistic,
                "lag": NA if oriter else draw(st.integers(min_value=1, max_value=20)),
                "null_model": "bootstrap"
                if oriter
                else draw(st.sampled_from(["gene_label", "tu_block"])),
                "value": draw(SMALL_FLOATS),
                "null_mean": NA if oriter else draw(SMALL_FLOATS),
                "null_lo": NA if oriter else draw(SMALL_FLOATS),
                "null_hi": NA if oriter else draw(SMALL_FLOATS),
                "stat": f"{statistic}_{i}" if oriter else statistic,
                "effect": draw(SMALL_FLOATS),
                "ci_low": draw(SMALL_FLOATS),
                "ci_high": draw(SMALL_FLOATS),
                "p_value": p,
                "q_value": NA if cd.is_na(p) else draw(PROBABILITIES),
                "method": draw(st.sampled_from(["moran_i_permutation", "ols_replichore"])),
                "n_perm": draw(st.integers(min_value=0, max_value=10_000)),
                "flags": draw(flags()),
            }
        )
    return pd.DataFrame(rows, columns=list(AUTOCORR_COLUMNS))


@st.composite
def island_frame(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(min_value=0, max_value=5))
    rows = []
    for i in range(n):
        start = draw(st.integers(min_value=0, max_value=5_000_000))
        end = start + draw(st.integers(min_value=1, max_value=50_000))
        rows.append(
            {
                "island_id": f"ISL_{i + 1}",
                "direction": draw(st.sampled_from(["high", "low"])),
                "contig": draw(st.sampled_from(["chr1", "p1"])),
                "start": start,
                "end": end,
                "tu_first": draw(TAGS),
                "tu_last": draw(TAGS),
                "n_tus": draw(st.integers(min_value=2, max_value=20)),
                "n_genes": draw(st.integers(min_value=1, max_value=200)),
                "genes": tuple(draw(st.lists(TAGS, min_size=1, max_size=4, unique=True))),
                "mean_abundance_log2": draw(SMALL_FLOATS),
                "detected_frac": draw(st.floats(min_value=0.0, max_value=1.0)),
                "scan_s": draw(SMALL_FLOATS),
                "block_len": draw(st.integers(min_value=1, max_value=20)),
                "stat": "scan_max",
                "effect": draw(SMALL_FLOATS),
                "ci_low": draw(SMALL_FLOATS),
                "ci_high": draw(SMALL_FLOATS),
                "p_value": draw(st.floats(min_value=1e-6, max_value=1.0)),
                "q_value": NA,
                "method": "scan_max_stat",
                "null_model": "circular_block_bootstrap",
                "n_perm": draw(st.integers(min_value=0, max_value=10_000)),
                "flags": draw(flags()),
            }
        )
    return pd.DataFrame(rows, columns=list(ISLAND_COLUMNS))


def _canonical(df: pd.DataFrame, lists: tuple[str, ...]) -> list[tuple[Any, ...]]:
    out = []
    for row in df.itertuples(index=False, name=None):
        cells: list[Any] = []
        for column, value in zip(df.columns, row, strict=True):
            if column in lists:
                cells.append(tuple(str(v) for v in (value or ())))
            elif isinstance(value, str):
                cells.append(value)
            elif cd.is_na(value):
                cells.append(None)
            else:
                cells.append(float(f"{float(value):.6g}"))
        out.append(tuple(cells))
    return out


@given(df=autocorrelation_frame())
@FAST
def test_autocorrelation_round_trip(tmp_path_factory: Any, df: pd.DataFrame) -> None:
    from giap.io.context_tsv import read_autocorrelation, write_autocorrelation

    out = tmp_path_factory.mktemp("ac") / "autocorrelation.tsv"
    write_autocorrelation(out, df)
    assert _canonical(read_autocorrelation(out), ("flags",)) == _canonical(df, ("flags",))


@given(df=island_frame())
@FAST
def test_islands_round_trip(tmp_path_factory: Any, df: pd.DataFrame) -> None:
    from giap.io.context_tsv import read_islands, write_islands

    out = tmp_path_factory.mktemp("isl") / "islands.tsv"
    write_islands(out, df)
    lists = ("genes", "flags")
    assert _canonical(read_islands(out), lists) == _canonical(df, lists)


@given(
    start=st.integers(min_value=0, max_value=10_000_000),
    extra=st.integers(min_value=1, max_value=100_000),
)
@FAST
def test_island_start_is_always_one_based_and_end_unchanged(tmp_path_factory: Any, start: int,
                                                            extra: int) -> None:  # fmt: skip
    """§6.5 in both directions, origin crossers (``end`` beyond the length) included."""
    from giap.io.context_tsv import read_islands, write_islands

    base = _one_island(start, start + extra)
    out = tmp_path_factory.mktemp("coord") / "islands.tsv"
    write_islands(out, base)
    text = out.read_text(encoding="utf-8").splitlines()[1].split("\t")
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert text[header.index("start")] == str(start + 1)
    assert text[header.index("end")] == str(start + extra)
    back = read_islands(out)
    assert int(back.loc[0, "start"]) == start and int(back.loc[0, "end"]) == start + extra


def _one_island(start: int, end: int) -> pd.DataFrame:
    row = {
        "island_id": "ISL_1", "direction": "high", "contig": "chr1", "start": start, "end": end,
        "tu_first": "TU_1", "tu_last": "TU_2", "n_tus": 2, "n_genes": 3, "genes": ("G_1",),
        "mean_abundance_log2": 20.0, "detected_frac": 0.5, "scan_s": 4.0, "block_len": 3,
        "stat": "scan_max", "effect": 1.0, "ci_low": 0.5, "ci_high": 1.5, "p_value": 0.01,
        "q_value": NA, "method": "scan_max_stat", "null_model": "circular_block_bootstrap",
        "n_perm": 999, "flags": (),
    }  # fmt: skip
    return pd.DataFrame([row], columns=list(ISLAND_COLUMNS))


@given(df=island_frame())
@FAST
def test_writing_is_a_pure_function_of_the_frame(tmp_path_factory: Any,
                                                 df: pd.DataFrame) -> None:  # fmt: skip
    """R7: no clock, no path, no dict ordering in the bytes."""
    from giap.io.context_tsv import write_islands

    base = tmp_path_factory.mktemp("pure")
    write_islands(base / "a.tsv", df)
    write_islands(base / "b.tsv", df.copy())
    assert (base / "a.tsv").read_bytes() == (base / "b.tsv").read_bytes()


@given(df=autocorrelation_frame())
@FAST
def test_row_order_is_preserved(tmp_path_factory: Any, df: pd.DataFrame) -> None:
    from giap.io.context_tsv import write_autocorrelation

    out = tmp_path_factory.mktemp("order") / "autocorrelation.tsv"
    write_autocorrelation(out, df)
    lines = out.read_text(encoding="utf-8").rstrip("\n").split("\n")[1:]
    header = AUTOCORR_COLUMNS.index("stat")
    assert [line.split("\t")[header] for line in lines] == [str(v) for v in df["stat"]]


@given(chosen=st.lists(VOCABULARY, max_size=4, unique=True))
@FAST
def test_any_vocabulary_subset_round_trips(tmp_path_factory: Any, chosen: list[str]) -> None:
    from giap.io.context_tsv import read_islands, write_islands

    df = _one_island(10, 20)
    df.at[0, "flags"] = tuple(sorted(chosen))
    out = tmp_path_factory.mktemp("flags") / "islands.tsv"
    write_islands(out, df)
    assert cd.flag_set(read_islands(out).loc[0, "flags"]) == set(chosen)


@given(lag=st.one_of(st.none(), st.integers(min_value=1, max_value=20)))
@FAST
def test_lag_is_an_integer_or_missing(tmp_path_factory: Any, lag: int | None) -> None:
    from giap.io.context_tsv import read_autocorrelation, write_autocorrelation

    row = {
        "contig": "chr1", "statistic": "morans_i" if lag else "oriter_slope",
        "lag": NA if lag is None else lag,
        "null_model": "gene_label" if lag else "bootstrap", "value": 0.5, "null_mean": NA,
        "null_lo": NA, "null_hi": NA, "stat": "morans_i" if lag else "oriter_slope_pooled",
        "effect": 0.5, "ci_low": NA,
        "ci_high": NA, "p_value": NA, "q_value": NA, "method": "moran_i_permutation",
        "n_perm": 0, "flags": ("not_estimable",),
    }  # fmt: skip
    df = pd.DataFrame([row], columns=list(AUTOCORR_COLUMNS))
    out = tmp_path_factory.mktemp("lag") / "autocorrelation.tsv"
    write_autocorrelation(out, df)
    back = read_autocorrelation(out)
    if lag is None:
        assert pd.isna(back.loc[0, "lag"])
    else:
        assert int(back.loc[0, "lag"]) == lag
