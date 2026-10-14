"""Hypothesis properties for peptide-table aggregation (issue M1-06; GIAP_SPEC.md §5.3).

- aggregation is invariant to row order (shuffled input → identical frames);
- one aggregated row per distinct kept peptide;
- ``peptide_runs`` sums equal the per-run intensity sums of the input (NA when a run has none).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config.model import InputConfig
from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides

FAST = settings(max_examples=100, deadline=None)

HEADER = ["peptide", "modified_peptide", "protein_ids", "q_value", "score", "intensity", "run"]
VALID = ["PEPTIDEK", "PEPTLDEK", "AAAAKPEP", "LLLLIIIK", "MNQRSTVWYK"]
INVALID = ["PEPK", "PEPT(UniMod:21)IDEK", "peptidek", "A" * 41]

ROW = st.tuples(
    st.sampled_from(VALID + INVALID),
    st.sampled_from(["NA", "PEPTIDEK", "M(ox)PEP"]),
    st.sampled_from(["NA", "P1", "P2;P3", "rev_P1", "CON__X"]),
    st.one_of(st.just("NA"), st.integers(0, 100).map(lambda i: str(i / 100))),
    st.one_of(st.just("NA"), st.integers(-50, 50).map(str)),
    # integral intensities keep float sums exact regardless of order
    st.one_of(st.just("NA"), st.integers(0, 10**6).map(str)),
    st.sampled_from(["r1", "r2", "r3"]),
)


def _run(rows: list[tuple[str, ...]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    text = "\t".join(HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "peptides.tsv"
        p.write_text(text, encoding="utf-8")
        df_raw, _ = read_peptide_table(p)
    df, _ = validate_peptides(df_raw, InputConfig())
    agg, runs = aggregate_peptides(df)
    return df, agg, runs


@FAST
@given(rows=st.lists(ROW, min_size=1, max_size=25), data=st.data())
def test_aggregation_invariant_to_row_order(
    rows: list[tuple[str, ...]], data: st.DataObject
) -> None:
    shuffled = data.draw(st.permutations(rows))
    _, agg1, runs1 = _run(rows)
    _, agg2, runs2 = _run(list(shuffled))
    pd.testing.assert_frame_equal(agg1.reset_index(drop=True), agg2.reset_index(drop=True))
    pd.testing.assert_frame_equal(runs1.reset_index(drop=True), runs2.reset_index(drop=True))


@FAST
@given(rows=st.lists(ROW, min_size=1, max_size=25))
def test_one_row_per_distinct_kept_peptide(rows: list[tuple[str, ...]]) -> None:
    _, agg, _ = _run(rows)
    kept = {r[0] for r in rows if r[0] in VALID}
    assert len(agg) == len(kept)
    assert set(agg["peptide"]) == kept
    assert int(agg["n_rows"].sum()) == sum(r[0] in VALID for r in rows)


@FAST
@given(rows=st.lists(ROW, min_size=1, max_size=25))
def test_peptide_runs_sums_match_input(rows: list[tuple[str, ...]]) -> None:
    _, agg, runs = _run(rows)
    expected: dict[tuple[str, str], list[float]] = {}
    counts: dict[tuple[str, str], int] = {}
    for r in rows:
        if r[0] not in VALID:
            continue
        key = (r[0], r[6])
        counts[key] = counts.get(key, 0) + 1
        vals = expected.setdefault(key, [])
        if r[5] != "NA":
            vals.append(float(r[5]))
    got = {
        (p, run): (i, n)
        for p, run, i, n in runs[["peptide", "run", "intensity", "n_rows"]].itertuples(index=False)
    }
    assert set(got) == set(expected)
    for key, vals in expected.items():
        intensity, n = got[key]
        assert n == counts[key]
        if vals:
            assert intensity == sum(vals)
        else:
            assert pd.isna(intensity)
    # peptide intensity = median over quantified runs
    for pep, grp in runs.groupby("peptide"):
        q = grp["intensity"].dropna()
        val = agg.loc[agg["peptide"] == pep, "intensity"].iloc[0]
        if len(q):
            assert val == q.median()
        else:
            assert pd.isna(val)
