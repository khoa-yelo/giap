"""Per-gene promoter features and ``promoter_between`` (issue M2-13 examples 7-8; §7.1.4).

``window_scores`` is the ``promoter_layer`` frame (contig, strand, end3p, score). A gene's own
windows are those on its strand whose ``end3p`` lies at S + k, k in [-200, +20] (``+``:
``end3p`` in [S - 200, S + 20]; ``-``: ``end3p = E - k + 1`` in [E - 19, E + 201]).
``n_promoter_windows_tu``: windows with score >= yes on the TU strand whose 3' end lies between the
TU's first start codon - 200 and the gene's own start codon + 20 in transcription direction
(inclusive, the same convention as the k range).
"""

from __future__ import annotations

import math

import pandas as pd
import pytest
from helpers.promoter_oracle import api

from giap.annotate.operons import Tu
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

CHR1 = Replicon("chr1", 100_000, "circular", "chromosome")
COLUMNS = ["locus_tag", "promoter_score", "promoter_score_tu", "n_promoter_windows_tu"]


def _gene(tag: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval("chr1", start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


def _tu(tu_id: str, genes: list[Gene]) -> Tu:
    """``genes`` in transcription order."""
    lo = min(g.interval.start for g in genes)
    hi = max(g.interval.end for g in genes)
    strand = genes[0].interval.strand
    up = (
        Interval("chr1", max(lo - 300, 0), lo + 50, strand)
        if strand == "+"
        else (Interval("chr1", hi - 50, hi + 300, strand))
    )
    return Tu(tu_id, "chr1", Interval("chr1", lo, hi, strand), strand,
              tuple(g.locus_tag for g in genes), up)  # fmt: skip


def _scores(rows: list[tuple[str, int, float]]) -> pd.DataFrame:
    return pd.DataFrame(
        [("chr1", s, p, x) for s, p, x in rows], columns=["contig", "strand", "end3p", "score"]
    )


def _features(genes: list[Gene], tus: list[Tu], rows: list, yes: float = 0.5) -> pd.DataFrame:
    df = api("gene_promoter_features")(genes, tus, _scores(rows), yes=yes)
    assert list(df.columns) == COLUMNS
    return df.set_index("locus_tag")


# + TU: F (first) then G. F's windows end3p in [800, 1020]; G's in [1200, 1420].
F_PLUS, G_PLUS = _gene("F", 1000, 1300, "+"), _gene("G", 1400, 1700, "+")


def test_example7_gene_and_tu_scores() -> None:
    rows = [
        ("+", 900, 0.97), ("+", 1010, 0.30),  # F's own windows
        ("+", 1300, 0.83), ("+", 1250, 0.40),  # G's own windows
        ("-", 1300, 0.99),  # other strand: never counted
        ("+", 1100, 0.99),  # between F's and G's ranges: no gene owns it
    ]  # fmt: skip
    df = _features([F_PLUS, G_PLUS], [_tu("TU_1", [F_PLUS, G_PLUS])], rows)
    assert df.loc["G", "promoter_score"] == pytest.approx(0.83)
    assert df.loc["G", "promoter_score_tu"] == pytest.approx(0.97)
    assert df.loc["F", "promoter_score"] == pytest.approx(0.97)
    assert df.loc["F", "promoter_score_tu"] == pytest.approx(0.97)
    # G: windows >= 0.5 on + with end3p in [800, 1420] -> 900, 1100, 1300
    assert df.loc["G", "n_promoter_windows_tu"] == 3
    assert df.loc["F", "n_promoter_windows_tu"] == 1  # [800, 1020] -> 900


def test_n_promoter_windows_tu_boundaries_plus() -> None:
    rows = [
        ("+", 799, 0.9), ("+", 800, 0.9),  # TU first start - 200 = 800 (inclusive)
        ("+", 1420, 0.9), ("+", 1421, 0.9),  # G start + 20 = 1420 (inclusive)
        ("+", 1000, 0.5), ("+", 1001, 0.4999),  # threshold: >= yes
        ("-", 1200, 0.9),
    ]  # fmt: skip
    df = _features([F_PLUS, G_PLUS], [_tu("TU_1", [F_PLUS, G_PLUS])], rows)
    assert df.loc["G", "n_promoter_windows_tu"] == 3  # 800, 1420, 1000
    assert df.loc["F", "n_promoter_windows_tu"] == 2  # 800, 1000 (1020 is F's limit)


def test_n_promoter_windows_tu_respects_yes() -> None:
    rows = [("+", 900, 0.6), ("+", 1300, 0.7)]
    df = _features([F_PLUS, G_PLUS], [_tu("TU_1", [F_PLUS, G_PLUS])], rows, yes=0.65)
    assert df.loc["G", "n_promoter_windows_tu"] == 1


# - TU: F = (1400, 1700, -) first in transcription order, then G = (1000, 1300, -).
# F's windows end3p in [1681, 1901]; G's in [1281, 1501]; G's TU range [1281, 1901].
F_MINUS, G_MINUS = _gene("F", 1400, 1700, "-"), _gene("G", 1000, 1300, "-")


def test_minus_strand_tu_scores_and_boundaries() -> None:
    rows = [
        ("-", 1901, 0.95), ("-", 1902, 0.99),  # F's k = -200 limit, then outside
        ("-", 1281, 0.6), ("-", 1280, 0.99),  # G's k = +20 limit, then outside
        ("-", 1501, 0.7),  # G's k = -200
        ("+", 1500, 0.99),
    ]  # fmt: skip
    df = _features([G_MINUS, F_MINUS], [_tu("TU_1", [F_MINUS, G_MINUS])], rows)
    assert df.loc["F", "promoter_score"] == pytest.approx(0.95)
    assert df.loc["G", "promoter_score"] == pytest.approx(0.7)
    assert df.loc["G", "promoter_score_tu"] == pytest.approx(0.95)
    assert df.loc["G", "n_promoter_windows_tu"] == 3  # 1901, 1281, 1501
    assert df.loc["F", "n_promoter_windows_tu"] == 1  # 1901 (1681..1901)


def test_single_gene_tu_scores_equal() -> None:
    df = _features([F_PLUS], [_tu("TU_1", [F_PLUS])], [("+", 950, 0.61)])
    assert df.loc["F", "promoter_score"] == df.loc["F", "promoter_score_tu"] == 0.61
    assert df.loc["F", "n_promoter_windows_tu"] == 1


def test_features_are_deterministic() -> None:
    rows = [("+", 900, 0.97), ("+", 1300, 0.83)]
    args = ([F_PLUS, G_PLUS], [_tu("TU_1", [F_PLUS, G_PLUS])])
    a = api("gene_promoter_features")(*args, _scores(rows), yes=0.5)
    b = api("gene_promoter_features")(*args, _scores(rows), yes=0.5)
    pd.testing.assert_frame_equal(a, b)


# ------------------------------------------------------------------ promoter_between (example 8)


UP, DOWN = _gene("U", 100, 400, "+"), _gene("D", 420, 700, "+")  # D's windows end3p in [220, 440]


def _between(rows: list, up: Gene = UP, down: Gene = DOWN) -> float | None:
    return api("promoter_between")(up, down, _scores(rows), CHR1)


def test_example8_promoter_between_max_after_upstream_end() -> None:
    rows = [("+", 401, 0.2), ("+", 440, 0.71), ("+", 400, 0.99), ("+", 300, 0.95)]
    assert _between(rows) == pytest.approx(0.71)


def test_example8_none_when_no_window_after_upstream_end() -> None:
    assert _between([("+", 400, 0.99), ("+", 300, 0.95)]) is None


def test_promoter_between_ignores_other_strand_and_non_own_windows() -> None:
    rows = [("-", 420, 0.99), ("+", 441, 0.99), ("+", 430, 0.33)]  # 441 = D's S + 21
    assert _between(rows) == pytest.approx(0.33)


def test_promoter_between_minus_strand_mirror() -> None:
    # transcription runs right to left: up = (600, 900, -), down = (200, 580, -).
    # D's windows end3p in [561, 781]; "after up's end" = 3'-most base 0-based < 600 -> end3p <= 600
    up, down = _gene("U", 600, 900, "-"), _gene("D", 200, 580, "-")
    rows = [("-", 600, 0.4), ("-", 601, 0.99), ("-", 561, 0.3)]
    assert _between(rows, up, down) == pytest.approx(0.4)
    assert _between([("-", 601, 0.99), ("-", 700, 0.9)], up, down) is None


def test_promoter_between_returns_plain_float() -> None:
    v = _between([("+", 420, 0.5)])
    assert isinstance(v, float) and not math.isnan(v)
