"""``read_detection_tier`` edge cases (#363 review S1, S2): the tier the manifest records.

genes.tsv ``detection_tier`` wins; ``detection_summary.tsv`` (``section = tier``,
``group_type = tier``) is the fallback; anything other than exactly one tier gives ``None``, so the
manifest falls back to preflight's provisional tier (MVP-03, row ``detection_tier``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gippa.io.manifest_sources import read_detection_tier

GENES_HEADER = "gene_id\tproduct\tdetection_tier\n"
SUMMARY_HEADER = "section\tgroup_type\tgroup\tvalue\n"


def _tables(tmp: Path, genes: str | None = None, summary: str | None = None) -> Path:
    tables = tmp / "tables"
    tables.mkdir()
    if genes is not None:
        (tables / "genes.tsv").write_text(genes, encoding="utf-8")
    if summary is not None:
        (tables / "detection_summary.tsv").write_text(summary, encoding="utf-8")
    return tables


def test_genes_tsv_tier(tmp_path: Path) -> None:
    genes = GENES_HEADER + "g1\tkinase\tmayu\ng2\tNA\tmayu\n"
    assert read_detection_tier(_tables(tmp_path, genes)) == "mayu"


def test_a_quote_leading_cell_does_not_hide_the_tier(tmp_path: Path) -> None:
    genes = GENES_HEADER + 'g1\t"quoted\tmayu\ng2\tplain\tmayu\n'
    assert read_detection_tier(_tables(tmp_path, genes)) == "mayu"


def test_two_tiers_give_none(tmp_path: Path) -> None:
    genes = GENES_HEADER + "g1\ta\tmayu\ng2\tb\tdecoy_picked\n"
    assert read_detection_tier(_tables(tmp_path, genes)) is None


def test_na_and_empty_cells_are_not_tiers(tmp_path: Path) -> None:
    genes = GENES_HEADER + "g1\ta\tNA\ng2\tb\t\ng3\tc\ttwo_peptide\n"
    assert read_detection_tier(_tables(tmp_path, genes)) == "two_peptide"


def test_genes_tsv_wins_over_the_summary(tmp_path: Path) -> None:
    genes = GENES_HEADER + "g1\ta\tmayu\n"
    summary = SUMMARY_HEADER + "tier\ttier\tdecoy_picked\t1\n"
    assert read_detection_tier(_tables(tmp_path, genes, summary)) == "mayu"


@pytest.mark.parametrize("genes", [None, GENES_HEADER, GENES_HEADER + "g1\ta\tNA\n"])
def test_summary_fallback(tmp_path: Path, genes: str | None) -> None:
    summary = (
        SUMMARY_HEADER
        + "tier\ttier\tmayu\t1\n"
        + "tier\tq_threshold\t0.01\t1\n"
        + "detected\ttier\tdecoy_picked\t3\n"
    )
    assert read_detection_tier(_tables(tmp_path, genes, summary)) == "mayu"


def test_no_tables_give_none(tmp_path: Path) -> None:
    assert read_detection_tier(_tables(tmp_path)) is None


def test_a_table_without_the_column_gives_none(tmp_path: Path) -> None:
    assert read_detection_tier(_tables(tmp_path, "gene_id\tproduct\ng1\ta\n")) is None


def test_crlf_and_bom_are_read(tmp_path: Path) -> None:
    genes = "﻿" + GENES_HEADER.replace("\n", "\r\n") + "g1\ta\tmayu\r\n"
    assert read_detection_tier(_tables(tmp_path, genes)) == "mayu"
