"""Placement-kind counts in ``detection_summary.tsv`` (M7-16b rule 17; spec examples 17-18; lead
decisions Q17-Q18; acceptance "exactly the 14 rows (16b, plus 16c sec_site and pyl_site) ...
Examples 17 and 18 hold.
Consistency checks ... ``qc_rows``' output is unchanged"; statistician review of the amendment,
item 4: counts, never rates).

Interface pinned here: ``giap.context.summary.placement_kind_rows(peptide_hits, unmapped)``, pure,
with ``peptide_hits`` the C11 frame (``peptide``, ``route``, ``diff_types``, ``flags``; the contract
reader gives lists, the TSV gives ``;``-joined strings, both are accepted) and ``unmapped`` the C13
frame or ``None``. It returns rows in ``SUMMARY_COLUMNS``: section ``qc``, metric ``n_peptides``
(count in ``value``), ``ci_low``/``ci_high`` NA, empty ``flags``, in this order:

* ``placement_kind``: ``exact``, ``tolerant``, ``frameshift`` (exclusive: exact if any row is
  exact; else tolerant if any non-frameshift row; else frameshift), ``n`` = placed peptides;
* ``tolerant_diff_type``: ``substitution``, ``swap``, ``block``, ``near_block`` (tolerant-kind
  peptides with a non-frameshift row carrying the type; not exclusive), ``n`` = tolerant-kind
  peptides;
* ``placement_flag``: ``readthrough``, ``n`` = placed peptides;
* ``unmapped_reason``: ``no_hit``, ``rejected_placement_fdr``, ``too_short``, ``contaminant_id``
  (C13 ``kind = peptide`` rows by ``reason``, both labels), ``n`` = peptide rows of C13.

``detection_summary`` sorts rows by section, ``group_type``, ``group``, ``metric``, so the rule's
order is the function's; the written table holds the same 14 rows.
"""

from __future__ import annotations

import gzip
import io
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[3]
ORDER: tuple[tuple[str, str], ...] = (
    ("placement_kind", "exact"),
    ("placement_kind", "tolerant"),
    ("placement_kind", "frameshift"),
    ("tolerant_diff_type", "substitution"),
    ("tolerant_diff_type", "swap"),
    ("tolerant_diff_type", "block"),
    ("tolerant_diff_type", "near_block"),
    ("placement_flag", "readthrough"),
    ("placement_flag", "pyl_site"),
    ("placement_flag", "sec_site"),
    ("unmapped_reason", "no_hit"),
    ("unmapped_reason", "rejected_placement_fdr"),
    ("unmapped_reason", "too_short"),
    ("unmapped_reason", "contaminant_id"),
)
GROUP_TYPES = {g for g, _ in ORDER}
SUMMARY_COLUMNS = (
    "section",
    "group_type",
    "group",
    "metric",
    "value",
    "ci_low",
    "ci_high",
    "n",
    "flags",
)


def _summary() -> Any:
    from giap.context import summary

    return summary


def kind_rows(hits: pd.DataFrame, unmapped: pd.DataFrame | None) -> pd.DataFrame:
    return _summary().placement_kind_rows(hits, unmapped)


def hit(
    peptide: str,
    route: str = "exact",
    diff_types: Sequence[str] = (),
    flags: Sequence[str] = (),
    start: int = 100,
    *,
    as_str: bool = False,
) -> dict[str, Any]:
    dt: Any = ";".join(diff_types) if as_str else list(diff_types)
    fl: Any = ";".join(flags) if as_str else list(flags)
    return {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "best_q_input": 0.001,
        "pep_input": 0.01,
        "score_input": 50.0,
        "contig": "chr1",
        "start": start,
        "end": start + 3 * len(peptide),
        "strand": "+",
        "frame": 1,
        "crosses_origin": False,
        "route": route,
        "n_diff": len(diff_types),
        "diff_types": dt if diff_types or as_str else [],
        "placement_q": 0.001,
        "length_class": "10-11",
        "locus_class": "cds_inframe",
        "uniqueness": "unique",
        "gene": "G1",
        "novel_candidate_id": None,
        "flags": fl,
    }


HIT_COLUMNS = tuple(hit("A"))
UNMAPPED_COLUMNS = ("sequence_id", "sequence", "kind", "label", "reason", "contaminant_match")


def hits_frame(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=list(HIT_COLUMNS))


def unmapped_frame(rows: Sequence[tuple[str, str, str, str]]) -> pd.DataFrame:
    """C13 rows from (sequence, kind, label, reason)."""
    recs = [
        (f"u{i}", s, k, lab, r, "CONT_1" if lab == "contaminant" else None)
        for i, (s, k, lab, r) in enumerate(rows)
    ]
    return pd.DataFrame(recs, columns=list(UNMAPPED_COLUMNS))


def example_17(as_str: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    s = as_str
    rows = [
        hit("EXACTAAAAK", as_str=s),
        hit("EXACTAAAAK", start=900, as_str=s),
        hit("EXACTCCCCK", as_str=s),
        hit("EXACTDDDDK", as_str=s),
        hit("APEPTIDEAAK", "tolerant", ["substitution"], ["readthrough"], as_str=s),
        hit("BPEPTIDEAAK", "tolerant", ["substitution", "swap"], as_str=s),
        hit("CPEPTIDEAAK", "tolerant", ["block"], as_str=s),
        hit("DPEPTIDEAAK", "tolerant", ["substitution"], as_str=s),
        hit("DPEPTIDEAAK", "tolerant", ["frameshift", "substitution"], start=700, as_str=s),
        hit("EPEPTIDEAAK", "tolerant", ["frameshift"], as_str=s),
    ]
    um = unmapped_frame(
        [
            ("NOHITAAAAK", "peptide", "unmapped", "no_hit"),
            ("NOHITCCCCK", "peptide", "unmapped", "no_hit"),
            ("SHORTK", "peptide", "unmapped", "too_short"),
            ("REJECTEDAAK", "peptide", "unmapped", "rejected_placement_fdr"),
            ("CONTAMAAAK", "peptide", "contaminant", "contaminant_id"),
            ("MPROTEINA", "protein", "unmapped", "below_identity"),
            ("MPROTEINC", "protein", "contaminant", "contaminant_id"),
        ]
    )
    return hits_frame(rows), um


def by_key(df: pd.DataFrame) -> dict[tuple[str, str], tuple[float, Any]]:
    return {
        (str(r.group_type), str(r.group)): (float(r.value), r.n) for r in df.itertuples(index=False)
    }


def _is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


# ---------------------------------------------------------------- shape


@pytest.mark.parametrize("as_str", [False, True], ids=["lists", "strings"])
def test_fixed_row_set_order_and_columns(as_str: bool) -> None:
    df = kind_rows(*example_17(as_str))
    assert tuple(df.columns) == SUMMARY_COLUMNS
    assert [(r.group_type, r.group) for r in df.itertuples(index=False)] == list(ORDER)
    assert set(df["section"]) == {"qc"} and set(df["metric"]) == {"n_peptides"}
    assert df["ci_low"].isna().all() and df["ci_high"].isna().all()
    assert all(f == "" for f in df["flags"])


def test_no_metric_is_a_rate_share_or_fdr() -> None:
    df = kind_rows(*example_17())
    assert not any(str(m).endswith(("_rate", "_share", "_fdr")) for m in df["metric"])


# ---------------------------------------------------------------- example 17


@pytest.mark.parametrize("as_str", [False, True], ids=["lists", "strings"])
def test_example_17_counts(as_str: bool) -> None:
    got = by_key(kind_rows(*example_17(as_str)))
    want = {
        ("placement_kind", "exact"): (3, 8),
        ("placement_kind", "tolerant"): (4, 8),
        ("placement_kind", "frameshift"): (1, 8),
        ("tolerant_diff_type", "substitution"): (3, 4),
        ("tolerant_diff_type", "swap"): (1, 4),
        ("tolerant_diff_type", "block"): (1, 4),
        ("tolerant_diff_type", "near_block"): (0, 4),
        ("placement_flag", "readthrough"): (1, 8),
        ("placement_flag", "pyl_site"): (0, 8),
        ("placement_flag", "sec_site"): (0, 8),
        ("unmapped_reason", "no_hit"): (2, 5),
        ("unmapped_reason", "rejected_placement_fdr"): (1, 5),
        ("unmapped_reason", "too_short"): (1, 5),
        ("unmapped_reason", "contaminant_id"): (1, 5),
    }
    assert {k: (v, int(n)) for k, (v, n) in got.items()} == {
        k: (float(v), n) for k, (v, n) in want.items()
    }


def test_example_17_qc_rows_unchanged_and_consistent() -> None:
    hits, um = example_17()
    qc = _summary().qc_rows(None, hits, um)
    q = {str(r.metric): float(r.value) for r in qc.itertuples(index=False) if r.group_type == "run"}
    assert (q["unmapped_peptides"], q["contaminant_peptides"]) == (4.0, 1.0)
    got = by_key(kind_rows(hits, um))
    reasons = sum(v for (g, _), (v, _n) in got.items() if g == "unmapped_reason")
    assert reasons == q["unmapped_peptides"] + q["contaminant_peptides"]
    # qc_rows does not carry the new groups (it is not changed)
    assert not set(qc["group_type"]) & GROUP_TYPES


def test_example_17_sums_and_v05_agreement() -> None:
    from giap.viz.mapping_sankey import flows

    hits, um = example_17()
    got = by_key(kind_rows(hits, um))
    kinds = {grp: v for (g, grp), (v, _n) in got.items() if g == "placement_kind"}
    assert sum(kinds.values()) == hits["peptide"].nunique()
    counts, _conflicts = flows(hits, None, um)
    assert kinds["exact"] == counts["peptide", "exact", "placed"]
    assert kinds["tolerant"] + kinds["frameshift"] == counts["peptide", "tolerant", "placed"]
    flagged = {p for p, f in zip(hits["peptide"], hits["flags"], strict=True) if "readthrough" in f}
    assert got["placement_flag", "readthrough"][0] == len(flagged)


# ---------------------------------------------------------------- rule details and edge cases


def test_exact_takes_precedence_over_tolerant_rows() -> None:
    hits = hits_frame([hit("MIXEDAAAAK"), hit("MIXEDAAAAK", "tolerant", ["swap"], start=500)])
    got = by_key(kind_rows(hits, unmapped_frame([])))
    assert got["placement_kind", "exact"][0] == 1 and got["placement_kind", "tolerant"][0] == 0
    assert got["tolerant_diff_type", "swap"][0] == 0  # not a tolerant-kind peptide


def test_diff_types_of_frameshift_rows_do_not_count() -> None:
    """``E``-like peptide: a frameshift row with a substitution is frameshift kind; ``D``-like
    peptide: its frameshift row adds nothing to the diff-type groups."""
    hits = hits_frame(
        [
            hit("FSONLYAAAK", "tolerant", ["frameshift", "swap"]),
            hit("MIXAAAAAAK", "tolerant", ["block"]),
            hit("MIXAAAAAAK", "tolerant", ["frameshift", "near_block"], start=800),
        ]
    )
    got = by_key(kind_rows(hits, None))
    assert got["placement_kind", "frameshift"][0] == 1
    assert got["placement_kind", "tolerant"][0] == 1
    assert got["tolerant_diff_type", "swap"][0] == 0
    assert got["tolerant_diff_type", "near_block"][0] == 0
    assert got["tolerant_diff_type", "block"][0] == 1


def test_near_block_and_multi_locus_flag_counted_once() -> None:
    hits = hits_frame(
        [
            hit("NEARBLOCKAK", "tolerant", ["near_block"]),
            hit("RTTWICEAAAK", "tolerant", ["substitution"], ["readthrough"]),
            hit(
                "RTTWICEAAAK",
                "tolerant",
                ["substitution"],
                ["readthrough", "strain_difference"],
                start=600,
            ),
            hit("RTTWICEAAAK", "tolerant", ["substitution"], ["strain_difference"], start=1200),
        ]
    )
    got = by_key(kind_rows(hits, None))
    assert got["tolerant_diff_type", "near_block"][0] == 1
    value, n = got["placement_flag", "readthrough"]
    assert (value, int(n)) == (1.0, 2)


def test_example_18_empty_run_writes_every_row_as_zero() -> None:
    df = kind_rows(hits_frame([]), unmapped_frame([]))
    assert [(r.group_type, r.group) for r in df.itertuples(index=False)] == list(ORDER)
    assert (df["value"] == 0).all()
    assert [int(n) for n in df["n"]] == [0] * len(ORDER)


def test_without_unmapped_the_reason_rows_are_na() -> None:
    hits, _um = example_17()
    df = kind_rows(hits, None)
    reason = df[df["group_type"] == "unmapped_reason"]
    assert len(reason) == 4
    assert reason["value"].isna().all()
    assert all(_is_na(n) for n in reason["n"])
    rest = df[df["group_type"] != "unmapped_reason"]
    assert not rest["value"].isna().any()


def test_rows_assemble_into_the_detection_summary() -> None:
    hits, um = example_17()
    s = _summary()
    table = s.detection_summary([s.qc_rows(None, hits, um), kind_rows(hits, um)])
    sub = table[table["group_type"].isin(GROUP_TYPES)]
    assert sorted(zip(sub["group_type"], sub["group"], strict=True)) == sorted(ORDER)
    assert set(sub["section"]) == {"qc"}


# ---------------------------------------------------------------- committed run outputs


def _read_tsv(path: Path) -> pd.DataFrame:
    data = path.read_bytes()
    text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode("utf-8")
    return pd.read_csv(
        io.StringIO(text), sep="\t", dtype=str, keep_default_na=False, na_values=["NA"]
    )


def _summary_rows(path: Path) -> pd.DataFrame:
    df = _read_tsv(path)
    return df[df["group_type"].isin(GROUP_TYPES)]


SNAPSHOT = REPO / "tests/data/snapshot/sim_plots/OUT"
GOLDEN_SUMMARIES = (
    REPO / "tests/data/golden/sim50k/m4/detection_summary.tsv",
    REPO / "tests/data/golden/sim100k_real/m4/detection_summary.tsv",
    SNAPSHOT / "tables/detection_summary.tsv",
)


@pytest.mark.parametrize("path", GOLDEN_SUMMARIES, ids=lambda p: str(p.relative_to(REPO)))
def test_committed_summaries_carry_the_12_rows_and_agree_with_the_qc_rows(path: Path) -> None:
    """After regeneration (lead): the 14 rows, metric ``n_peptides``, and within one table
    sum(unmapped_reason) = unmapped_peptides + contaminant_peptides and sum(placement_kind) =
    the distinct placed peptides (the ``n`` of ``engine_protein_mismatch_rate``)."""
    df = _read_tsv(path)
    sub = df[df["group_type"].isin(GROUP_TYPES)]
    assert sorted(zip(sub["group_type"], sub["group"], strict=True)) == sorted(ORDER), path
    assert set(sub["metric"]) == {"n_peptides"} and set(sub["section"]) == {"qc"}
    run = df[(df["group_type"] == "run") & (df["group"] == "all")].set_index("metric")
    reasons = sub[sub["group_type"] == "unmapped_reason"]["value"].astype(float).sum()
    unm = float(run.loc["unmapped_peptides", "value"]) + float(
        run.loc["contaminant_peptides", "value"]
    )
    assert reasons == unm
    kinds = sub[sub["group_type"] == "placement_kind"]["value"].astype(float).sum()
    n_placed = run.loc["engine_protein_mismatch_rate", "n"]
    assert kinds == (0.0 if _is_na(n_placed) else float(n_placed))


def test_snapshot_rows_equal_the_function_on_its_own_mapping_outputs() -> None:
    ph = _read_tsv(SNAPSHOT / "mapping/peptide_hits.tsv")
    um = _read_tsv(SNAPSHOT / "mapping/unmapped.tsv")
    want = kind_rows(ph, um)
    got = _summary_rows(SNAPSHOT / "tables/detection_summary.tsv")
    w = {(r.group_type, r.group): (float(r.value), int(r.n)) for r in want.itertuples(index=False)}
    g = {(r.group_type, r.group): (float(r.value), int(r.n)) for r in got.itertuples(index=False)}
    assert g == w


def test_snapshot_cross_checks_against_v05() -> None:
    from giap.viz.mapping_sankey import flows

    ph = _read_tsv(SNAPSHOT / "mapping/peptide_hits.tsv")
    um = _read_tsv(SNAPSHOT / "mapping/unmapped.tsv")
    got = by_key(kind_rows(ph, um))
    counts, _ = flows(ph, None, um)
    assert got["placement_kind", "exact"][0] == counts["peptide", "exact", "placed"]
    assert (
        got["placement_kind", "tolerant"][0] + got["placement_kind", "frameshift"][0]
        == counts["peptide", "tolerant", "placed"]
    )
    assert (
        sum(v for (g, _), (v, _n) in got.items() if g == "placement_kind")
        == ph["peptide"].nunique()
    )
    flagged = {
        p
        for p, f in zip(ph["peptide"], ph["flags"], strict=True)
        if "readthrough" in str(f).split(";")
    }
    assert got["placement_flag", "readthrough"][0] == len(flagged)
