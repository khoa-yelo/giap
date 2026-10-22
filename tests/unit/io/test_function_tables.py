"""M8-06a item 4: writers of ``tables/category_summary.tsv``, ``enrichment_detection.tsv`` and
``enrichment_gsea.tsv`` (docs/specs/M8-06-function-step-and-run.md; GIAP_SPEC.md §9.2, ADR-0006).

Pinned: cells through ``giap.io.tsv`` (``NA``, ``%.6g``, ``true``/``false``, ``;``-joined
``leading_edge`` and ``flags``, an empty list is ``NA``); a None ``parent_group`` and a NA
``redundant`` are ``NA``; columns in each module's ``COLUMNS`` order whatever the frame's order;
GSEA rows sorted stably by (collection in ``COLLECTIONS`` order, term); a header-only frame writes
the header row only.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def _summary_row(**over: Any) -> dict[str, Any]:
    row = {"collection": "cog", "category": "J", "name": "Translation", "n_genes": 3,
           "n_detected": 2, "expected": 0.9 + 0.5 + 0.4, "o_minus_e": 2 - (0.9 + 0.5 + 0.4),
           "ci_low": math.nan, "ci_high": math.nan, "median_abundance_log2": math.nan,
           "coverage": 0.75, "flags": ["ci_few_units"]}  # fmt: skip
    return row | over


def _ora_row(**over: Any) -> dict[str, Any]:
    row = {"collection": "go_bp", "term": "GO:0000101", "name": "toy", "direction": "over",
           "n_universe": 40, "n_in_set": 12, "n_detected": 20, "k_detected_in_set": 9,
           "n_in_set_genes": 15, "k_detected_in_set_genes": 11, "odds_weight": 1.23456789,
           "p_classic": 0.0123456789, "p_noelim": 0.02, "parent_group": None, "coverage": 0.5,
           "stat": 9.0, "effect": 2.5, "ci_low": math.nan, "ci_high": math.nan,
           "p_value": 1e-300, "q_value": 0.04, "method": "cpb_tu_ora",
           "null_model": "conditional_poisson_binomial_tu", "n_perm": 0,
           "flags": ["ci_unestimated", "ora_uncalibrated"]}  # fmt: skip
    return row | over


def _gsea_row(**over: Any) -> dict[str, Any]:
    row = {"collection": "cog", "term": "J", "name": "Translation", "size_tus": 12, "es": 0.5,
           "nes": 1.5, "log2err": 0.25, "leading_edge": ["g1", "g2"], "redundant": False,
           "parent_group": None, "coverage": 0.8, "stat": 0.5, "effect": 1.5,
           "ci_low": math.nan, "ci_high": math.nan, "p_value": 0.01, "q_value": 0.03,
           "method": "fgsea_simple", "null_model": "tu_set_permutation", "n_perm": 1000,
           "flags": ["ci_unestimated"]}  # fmt: skip
    return row | over


def _frame(rows: list[dict[str, Any]], columns: tuple[str, ...], kind: str) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(columns))
    if kind == "gsea":
        df["redundant"] = df["redundant"].astype("boolean")
    if "parent_group" in df:
        df["parent_group"] = df["parent_group"].astype(object)
    return df


def _mods():  # type: ignore[no-untyped-def]
    from giap.function import gsea, ora, summary
    from giap.io import function_tables as ft

    return {
        "summary": (ft.write_category_summary_tsv, tuple(summary.COLUMNS), _summary_row),
        "ora": (ft.write_enrichment_detection_tsv, tuple(ora.COLUMNS), _ora_row),
        "gsea": (ft.write_enrichment_gsea_tsv, tuple(gsea.COLUMNS), _gsea_row),
    }


def _write(kind: str, rows: list[dict[str, Any]], path: Path) -> list[list[str]]:
    writer, cols, _ = _mods()[kind]
    writer(path, _frame(rows, cols, kind))
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    return [line.split("\t") for line in text.splitlines()]


def _cell(lines: list[list[str]], row: int, col: str) -> str:
    return lines[row + 1][lines[0].index(col)]


@pytest.mark.parametrize("kind", ["summary", "ora", "gsea"])
def test_header_only_frame_writes_the_header_row_only(kind: str, tmp_path: Path) -> None:
    writer, cols, _ = _mods()[kind]
    path = tmp_path / f"{kind}.tsv"
    writer(path, _frame([], cols, kind))
    assert path.read_text(encoding="utf-8") == "\t".join(cols) + "\n"


@pytest.mark.parametrize("kind", ["summary", "ora", "gsea"])
def test_columns_follow_the_module_columns_verbatim(kind: str, tmp_path: Path) -> None:
    writer, cols, row = _mods()[kind]
    df = _frame([row()], cols, kind)
    df = df[list(reversed(cols))]  # the writer, not the frame, fixes the order
    df["internal_extra"] = 1
    path = tmp_path / f"{kind}.tsv"
    writer(path, df)
    header = path.read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
    assert tuple(header) == cols


def test_summary_cells(tmp_path: Path) -> None:
    rows = [_summary_row(), _summary_row(category="K", flags=[], ci_low=-0.5, ci_high=1.25,
                                          median_abundance_log2=21.123456789)]  # fmt: skip
    lines = _write("summary", rows, tmp_path / "s.tsv")
    assert lines[1] == ["cog", "J", "Translation", "3", "2", "1.8", "0.2", "NA", "NA", "NA",
                        "0.75", "ci_few_units"]  # fmt: skip
    assert _cell(lines, 1, "flags") == "NA"  # an empty flag list
    assert _cell(lines, 1, "ci_low") == "-0.5"
    assert _cell(lines, 1, "median_abundance_log2") == "21.1235"  # %.6g


def test_detection_cells(tmp_path: Path) -> None:
    rows = [_ora_row(), _ora_row(direction="under", parent_group="GO:0000100",
                                 flags=["ci_unestimated", "unestimated", "null_degenerate",
                                        "tu_recal_skipped", "ora_uncalibrated"])]  # fmt: skip
    lines = _write("ora", rows, tmp_path / "d.tsv")
    assert _cell(lines, 0, "parent_group") == "NA"  # None
    assert _cell(lines, 1, "parent_group") == "GO:0000100"
    assert _cell(lines, 0, "flags") == "ci_unestimated;ora_uncalibrated"
    assert _cell(lines, 1, "flags") == ("ci_unestimated;unestimated;null_degenerate;"
                                        "tu_recal_skipped;ora_uncalibrated")  # fmt: skip
    assert _cell(lines, 0, "odds_weight") == "1.23457"
    assert _cell(lines, 0, "p_classic") == "0.0123457"
    assert _cell(lines, 0, "p_value") == "1e-300"
    assert _cell(lines, 0, "ci_low") == "NA"
    assert _cell(lines, 0, "n_perm") == "0"
    assert _cell(lines, 0, "n_universe") == "40"


def test_gsea_list_and_boolean_cells(tmp_path: Path) -> None:
    rows = [
        _gsea_row(term="J", leading_edge=["g1", "g2", "g3"], redundant=True),
        _gsea_row(term="K", leading_edge=[], redundant=False, parent_group="K"),
        _gsea_row(term="L", redundant=pd.NA,
                  flags=["ci_unestimated", "size_filtered", "genomic_cluster"]),
        _gsea_row(term="M", nes=math.nan, es=0.0,
                  flags=["ci_unestimated", "unestimated", "detection_confounded"]),
    ]  # fmt: skip
    lines = _write("gsea", rows, tmp_path / "g.tsv")
    assert [_cell(lines, i, "leading_edge") for i in range(4)] == ["g1;g2;g3", "NA", "g1;g2",
                                                                   "g1;g2"]  # fmt: skip
    assert [_cell(lines, i, "redundant") for i in range(4)] == ["true", "false", "NA", "false"]
    assert [_cell(lines, i, "parent_group") for i in range(4)] == ["NA", "K", "NA", "NA"]
    assert _cell(lines, 2, "flags") == "ci_unestimated;size_filtered;genomic_cluster"
    assert _cell(lines, 3, "nes") == "NA"
    assert _cell(lines, 3, "es") == "0"
    assert _cell(lines, 0, "n_perm") == "1000"


def test_gsea_rows_sorted_by_collection_order_then_term(tmp_path: Path) -> None:
    """``COLLECTIONS`` order (cog, kegg, go_bp, go_mf, go_cc, go_slim), not alphabetical."""
    rows = [
        _gsea_row(collection="go_slim", term="GO:0000101"),
        _gsea_row(collection="go_bp", term="GO:0000102"),
        _gsea_row(collection="kegg", term="map00020"),
        _gsea_row(collection="go_cc", term="GO:0000301"),
        _gsea_row(collection="cog", term="K"),
        _gsea_row(collection="go_bp", term="GO:0000101"),
        _gsea_row(collection="go_mf", term="GO:0000201"),
        _gsea_row(collection="kegg", term="map00010"),
        _gsea_row(collection="cog", term="J"),
    ]
    lines = _write("gsea", rows, tmp_path / "g.tsv")
    got = [(r[0], r[1]) for r in lines[1:]]
    assert got == [("cog", "J"), ("cog", "K"), ("kegg", "map00010"), ("kegg", "map00020"),
                   ("go_bp", "GO:0000101"), ("go_bp", "GO:0000102"), ("go_mf", "GO:0000201"),
                   ("go_cc", "GO:0000301"), ("go_slim", "GO:0000101")]  # fmt: skip


def test_gsea_sort_is_stable(tmp_path: Path) -> None:
    rows = [_gsea_row(term="K", name="second"), _gsea_row(term="J", name="a"),
            _gsea_row(term="K", name="first")]  # fmt: skip
    lines = _write("gsea", rows, tmp_path / "g.tsv")
    assert [(r[1], r[2]) for r in lines[1:]] == [("J", "a"), ("K", "second"), ("K", "first")]


def test_writers_do_not_modify_their_input(tmp_path: Path) -> None:
    for kind, (writer, cols, row) in _mods().items():
        rows = [row(), row()] if kind != "gsea" else [_gsea_row(term="K"), _gsea_row(term="J")]
        df = _frame(rows, cols, kind)
        before = df.copy(deep=True)
        writer(tmp_path / f"{kind}.tsv", df)
        pd.testing.assert_frame_equal(df, before)


_TERMS = st.sampled_from(["J", "K", "L", "map00010", "map00020", "GO:0000101", "GO:0000102"])
_COLLS = st.sampled_from(["cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim"])


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(keys=st.lists(st.tuples(_COLLS, _TERMS), min_size=1, max_size=12, unique=True),
       data=st.data())  # fmt: skip
def test_gsea_bytes_do_not_depend_on_input_row_order(
    keys: list[tuple[str, str]], data: st.DataObject, tmp_path: Path
) -> None:
    """Property: with unique (collection, term), any permutation of the rows writes the same
    file, and the rows come out in (COLLECTIONS order, term) order."""
    from giap.function.collections import COLLECTIONS

    rows = [_gsea_row(collection=c, term=t, es=0.1 * i) for i, (c, t) in enumerate(keys)]
    perm = data.draw(st.permutations(rows))
    a = _write("gsea", rows, tmp_path / "a.tsv")
    b = _write("gsea", list(perm), tmp_path / "b.tsv")
    assert a == b
    order = {c: i for i, c in enumerate(COLLECTIONS)}
    got = [(order[r[0]], r[1]) for r in a[1:]]
    assert got == sorted(got)
    assert len(a) == len(rows) + 1
