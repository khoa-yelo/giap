"""M6-06 review of #148: out-of-bin warning (B1), imputed label (N2), R1-R3 flags and the pi0
warning (code reviewer and statistician)."""

from __future__ import annotations

import logging

import pandas as pd
import pytest
from helpers.missed_data import block, cfg, example_frame, genes_for, stack

from gippa.operons.missed import TEST, missed_operons

EDGES_NO_OVERLAP = [0.0, 50.0, 150.0, 300.0, float("inf")]


def _scored(rows: pd.DataFrame) -> pd.DataFrame:
    return rows[rows["test"] == TEST]


def test_pairs_outside_every_bin_are_counted_in_a_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    pairs = example_frame()
    n_no = int(((pairs["class"] == "NO") & (pairs["intergenic_bp"] < 0)).sum())
    n_d = int(((pairs["class"] != "NO") & (pairs["intergenic_bp"] < 0)).sum())
    assert n_no and n_d
    with caplog.at_level(logging.WARNING, logger="gippa.operons.missed"):
        res = missed_operons(pairs, cfg(distance_bins_bp=EDGES_NO_OVERLAP))
    msg = [r.getMessage() for r in caplog.records if "outside every" in r.getMessage()]
    assert len(msg) == 1
    assert f"{n_no} eligible NO pair(s) and {n_d} decoy pair(s)" in msg[0]
    assert res.n_targets == 300 - n_no


def test_default_bins_never_warn_about_out_of_bin_pairs(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING, logger="gippa.operons.missed"):
        missed_operons(example_frame(), cfg())
    assert not any("outside every" in r.getMessage() for r in caplog.records)


@pytest.mark.parametrize("side", ["gene_a", "gene_b"])
def test_a_no_pair_with_an_imputed_gene_is_labelled(side: str) -> None:
    pairs = example_frame()
    no = pairs.index[pairs["class"] == "NO"]
    genes = genes_for(pairs, abundance_imputed=dict.fromkeys(pairs.loc[no[:5], side], True))
    rows = missed_operons(pairs, cfg(), genes=genes).rows.set_index("pair_id")
    for pid in pairs.loc[no, "pair_id"]:
        has = "abundance_imputed" in rows.at[pid, "flags"]
        assert has == (pid in set(pairs.loc[no[:5], "pair_id"]))


def test_imputed_label_is_kept_on_the_too_few_path() -> None:
    pairs = stack(block("NO", [0.1] * 5), block("CV", [0.2] * 30))
    genes = genes_for(pairs, abundance_imputed={pairs.at[0, "gene_a"]: True})
    rows = missed_operons(pairs, cfg(), genes=genes).rows
    assert rows.at[0, "flags"] == ["abundance_imputed", "too_few_no_pairs"]


def test_relative_scale_rows_are_unestimated() -> None:
    pairs = example_frame()
    tags = [*pairs["gene_a"], *pairs["gene_b"]]
    genes = genes_for(pairs, abundance_scale=dict.fromkeys(tags, "relative"))
    rows = missed_operons(pairs, cfg(), genes=genes).rows
    assert all({"relative_scale", "unestimated"} <= set(f) for f in rows["flags"])


def test_scored_rows_flag_the_missing_ci_and_unscored_rows_do_not() -> None:
    pairs = stack(example_frame(), block("NO", [0.1], detected=(True, False)))
    rows = missed_operons(pairs, cfg()).rows
    sc = _scored(rows)
    assert len(sc) == 300
    assert sc["ci_low"].isna().all() and sc["ci_high"].isna().all()
    assert all(f.count("ci_unestimated") == 1 and "unestimated" not in f for f in sc["flags"])
    assert all("ci_unestimated" not in f for f in rows.loc[rows["test"] == "none", "flags"])


def test_scored_rows_keep_their_input_flags() -> None:
    pairs = example_frame()
    first = pairs.index[pairs["class"] == "NO"][0]
    pairs.at[first, "flags"] = ["rna_between"]
    rows = missed_operons(pairs, cfg()).rows.set_index("pair_id")
    assert rows.at[pairs.at[first, "pair_id"], "flags"] == ["rna_between", "ci_unestimated"]


def test_pi0_below_one_warns(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING, logger="gippa.operons.missed"):
        missed_operons(example_frame(), cfg(pi0=0.8))
    assert any("operons.pi0 = 0.8" in r.getMessage() for r in caplog.records)
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="gippa.operons.missed"):
        missed_operons(example_frame(), cfg())
    assert not any("operons.pi0" in r.getMessage() for r in caplog.records)
