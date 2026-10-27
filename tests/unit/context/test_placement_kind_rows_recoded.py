"""The 16c placement-flag rows and schema values (M7-16c; spec rule 17 with rule 19 "The
``placement_flag`` groups ``sec_site`` and ``pyl_site`` are added", rule 18 "16c adds ``pyl_site``,
``sec_site``", example 20; test plan C9).

``placement_kind_rows`` writes 14 rows: the 12 of 16b plus ``placement_flag`` ``sec_site`` and
``pyl_site`` (placed peptides with at least one row carrying the flag; ``n`` = placed peptides),
zeros included. The ``placement_flag`` rows stay together between ``tolerant_diff_type`` and
``unmapped_reason``; their order inside the group is not pinned here. The ``peptide_hits`` schema
``flags`` pattern accepts ``sec_site`` and ``pyl_site``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from test_placement_kind_rows import (
    ORDER,
    SUMMARY_COLUMNS,
    by_key,
    example_17,
    hit,
    hits_frame,
    kind_rows,
    unmapped_frame,
)

REPO = Path(__file__).resolve().parents[3]
SCHEMA = REPO / "src/gippa/schemas/tables/peptide_hits.schema.json"
NEW = (("placement_flag", "pyl_site"), ("placement_flag", "sec_site"))
ORDER_16C = ORDER  # d17a641f already added NEW to ORDER
assert set(NEW) <= set(ORDER_16C) and len(ORDER_16C) == 14
GROUP_TYPES = {g for g, _ in ORDER_16C}


def example_17_plus(as_str: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Example 17 plus ``F`` (``sec_site`` at one of two loci), ``G`` (``pyl_site``) and ``H``
    (``readthrough;sec_site``, a row consuming a terminal and a recoded cell)."""
    hits, um = example_17(as_str)
    s = as_str
    more = [
        hit("FSECPEPTIDK", "tolerant", ["substitution"], ["sec_site"], as_str=s),
        hit("FSECPEPTIDK", "tolerant", ["substitution"], [], start=400, as_str=s),
        hit("GPYLPEPTIDK", "tolerant", ["substitution"], ["pyl_site"], as_str=s),
        hit(
            "HBOTHPEPTIDEK",
            "tolerant",
            ["substitution", "substitution"],
            ["readthrough", "sec_site"],
            as_str=s,
        ),
    ]
    rows = [*hits.to_dict("records"), *more]
    return hits_frame(rows), um


def _rows(df: pd.DataFrame) -> list[tuple[str, str]]:
    return [(str(r.group_type), str(r.group)) for r in df.itertuples(index=False)]


def test_fourteen_rows_and_the_group_order() -> None:
    df = kind_rows(*example_17_plus())
    got = _rows(df)
    assert len(got) == 14 and sorted(got) == sorted(ORDER_16C)
    assert tuple(df.columns) == SUMMARY_COLUMNS
    types = [g for g, _ in got]
    want_types = ["placement_kind"] * 3 + ["tolerant_diff_type"] * 4 + ["placement_flag"] * 3
    assert types == want_types + ["unmapped_reason"] * 4
    assert set(df["metric"]) == {"n_peptides"} and set(df["section"]) == {"qc"}


def test_example_17_without_recoded_flags_has_zero_sec_and_pyl_rows() -> None:
    got = by_key(kind_rows(*example_17()))
    for key in NEW:
        value, n = got[key]
        assert (value, int(n)) == (0.0, 8)
    assert got["placement_flag", "readthrough"][0] == 1


def _expect(got: dict, key: tuple[str, str], value: int, n: int) -> None:  # type: ignore[type-arg]
    v, m = got[key]
    assert (v, int(m)) == (float(value), n), key


def test_sec_and_pyl_counts_and_n() -> None:
    for as_str in (False, True):
        got = by_key(kind_rows(*example_17_plus(as_str)))
        _expect(got, ("placement_kind", "tolerant"), 7, 11)
        _expect(got, ("placement_flag", "readthrough"), 2, 11)  # A and H
        _expect(got, ("placement_flag", "sec_site"), 2, 11)  # F (once) and H
        _expect(got, ("placement_flag", "pyl_site"), 1, 11)  # G
        _expect(got, ("tolerant_diff_type", "substitution"), 6, 7)


def test_example_18_empty_run_has_the_14_rows_at_zero() -> None:
    df = kind_rows(hits_frame([]), unmapped_frame([]))
    assert sorted(_rows(df)) == sorted(ORDER_16C)
    assert (df["value"] == 0).all() and [int(n) for n in df["n"]] == [0] * 14


def test_rows_assemble_into_the_detection_summary_with_the_new_groups() -> None:
    from gippa.context import summary as s

    hits, um = example_17_plus()
    table = s.detection_summary([s.qc_rows(None, hits, um), kind_rows(hits, um)])
    sub = table[table["group_type"].isin(GROUP_TYPES)]
    assert sorted(zip(sub["group_type"], sub["group"], strict=True)) == sorted(ORDER_16C)


# ---------------------------------------------------------------- the schema pattern


def _pattern() -> re.Pattern[str]:
    doc = json.loads(SCHEMA.read_text(encoding="utf-8"))
    (f,) = [f for f in doc["fields"] if f["name"] == "flags"]
    return re.compile(f["constraints"]["pattern"])


def test_schema_pattern_accepts_sec_site_and_pyl_site() -> None:
    pat = _pattern()
    for value in (
        "sec_site",
        "pyl_site",
        "readthrough;sec_site",
        "pyl_site;strain_difference",
        "multi;sec_site;strain_difference",
    ):
        assert pat.fullmatch(value), value


def test_schema_pattern_still_rejects_near_misses() -> None:
    pat = _pattern()
    for value in ("sec", "sec_sites", "pyl", "secsite", "sec_site;", "Sec_site"):
        assert not pat.fullmatch(value), value
