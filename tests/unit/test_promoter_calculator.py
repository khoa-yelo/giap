"""Promoter Calculator: region FASTA, command, CSV parser, ``tx_layer``, TU features and
``tx_rate_between`` (M2-14; GIPPA_SPEC.md §7.1.1, §7.1.4, §9.3).

Examples 1-9 of docs/specs/M2-14-promoter-calculator.md and its edge cases, with hand-written CSVs
in the real 27-column Promoter Calculator 1.2.4 format (spec Q2: the recorded header).

Conventions pinned here where the spec is silent (test author; see the PR description):

* Q1 (proposed default adopted): a TU without a kept hit has ``tx_rate_max = tx_rate_sum = 0`` and
  ``tss_distance_bp`` NA.
* ``PcHit.tss`` is 1-based in forward region coordinates (spec). The tool's ``TSS`` column is the
  0-based index of the +1 base for ``+`` hits and the 1-based forward position of the +1 base for
  ``-`` hits (Promoter Calculator 1.2.4 ``wrapper.py``; verified on every recorded row by the
  contract test via the ``ITR`` column), so ``parse_pc_csv`` adds 1 to ``+`` rows only.
* ``parse_pc_csv`` reads CRLF (the tool) and LF (the nf stub) files, returns rows in file order,
  returns ``[]`` for a header-only or empty file (M2-07 writes the header when the tool wrote no
  CSV), and raises ``ValueError`` when ``TSS``, ``strand`` or ``Tx_rate`` is missing.
* ``tx_layer`` reduces genome TSSs modulo the replicon length on circular replicons (origin-crossing
  regions), so ``0 <= tss < length`` like the precomputed ``tx_rates`` layer; its ``strand`` is the
  TU's genome strand. TUs missing from ``hits`` (no CSV) contribute no rows.
* ``tu_tx_features`` takes ``replicons`` as well (like ``tx_layer`` and M2-13 decision 7): the
  distance to the start codon must be measured across the origin. One row per TU, in TU order. The
  first start codon is ``interval.start`` (+) / ``interval.end - 1`` (-) of the TU.
* ``tx_rate_between``: TSSs strictly between ``up``'s 3' end and ``down``'s start codon in
  transcription order; ``None`` when there is none. Circular replicons wrap. Boundary bases
  (lead decision 5): the first base after ``up``'s stop is kept, ``up``'s last base and ``down``'s
  start-codon base are dropped, on both strands; adjacent (gap 0) and overlapping genes give
  ``None`` without reading the layer.
* ``_signed`` reduces to ``(-length/2, length/2]``: a distance of exactly half an even contig is
  positive (upstream), and ``-length/2`` maps to ``+length/2`` (lead decision 5).
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import pandas as pd
import pytest
from helpers.promoter_calculator_fixtures import (
    HEADER,
    header_only,
    make_tu,
    revcomp,
    write_pc_csv,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

L = 10_000
_rng = random.Random(2014)
GENOME = {"chr1": "".join(_rng.choice("ACGT") for _ in range(L))}
REPS = {"chr1": Replicon("chr1", L, "circular")}
LIN = {"chr1": Replicon("chr1", L, "linear")}
# TU_1: + , first start codon S = 1000, region [850, 1050) (example 1)
TU1 = make_tu("TU_1", "chr1", "+", (1000, 1600), (850, 1050))
# TU_2: + , no hits (example 8)
TU2 = make_tu("TU_2", "chr1", "+", (3000, 3300), (2900, 3050))
# TU_3: - , gene interval end 5600 → start codon first base 5599, region [5550, 5800) (example 2)
TU3 = make_tu("TU_3", "chr1", "-", (5000, 5600), (5550, 5800))
TUS = [TU1, TU2, TU3]
TX_COLS = ["contig", "strand", "tss", "tx_rate", "tu_id"]
FEAT_COLS = ["tu_id", "tx_rate_max", "tx_rate_sum", "tss_distance_bp"]


def _mod():
    from gippa.annotate import promoter_calculator as pc

    return pc


def _hit(tss: int, strand: str, rate: float):
    return _mod().PcHit(tss=tss, strand=strand, tx_rate=rate)


def _rows(df: pd.DataFrame) -> set[tuple[object, ...]]:
    return {
        (str(c), str(s), int(t), float(x), str(u))
        for c, s, t, x, u in df[TX_COLS].itertuples(index=False, name=None)
    }


def _feat(df: pd.DataFrame, tu_id: str) -> dict[str, object]:
    sub = df[df["tu_id"] == tu_id]
    assert len(sub) == 1, f"{tu_id}: {len(sub)} feature rows"
    return sub.iloc[0].to_dict()


# ---------------------------------------------------------------- examples 1, 2, 9: region FASTA


def test_example1_plus_region_is_the_genome_slice_on_one_line() -> None:
    text = _mod().region_fasta_text(TU1, GENOME)
    assert text == ">TU_1\n" + GENOME["chr1"][850:1050] + "\n"
    assert len(text.splitlines()[1]) == 200


def test_example2_minus_region_is_the_reverse_complement() -> None:
    text = _mod().region_fasta_text(TU3, GENOME)
    assert text == ">TU_3\n" + revcomp(GENOME["chr1"][5550:5800]) + "\n"


@pytest.mark.parametrize("tu", TUS, ids=lambda t: t.tu_id)
def test_example9_exactly_one_record_and_two_lines(tu) -> None:
    text = _mod().region_fasta_text(tu, GENOME)
    assert text.count(">") == 1
    assert text.endswith("\n") and len(text.splitlines()) == 2
    assert text.splitlines()[0] == f">{tu.tu_id}"


def test_long_region_is_not_wrapped() -> None:
    """The tool joins all lines anyway, but the spec fixes one sequence line (no 60/80 wrap)."""
    tu = make_tu("TU_9", "chr1", "+", (2000, 2100), (1000, 2050))
    text = _mod().region_fasta_text(tu, GENOME)
    assert text == ">TU_9\n" + GENOME["chr1"][1000:2050] + "\n"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_crossing_region_is_one_contiguous_wrapped_slice(strand: str) -> None:
    tu = make_tu(
        "TU_7", "chr1", strand, (100, 700) if strand == "+" else (9300, 9900), (9800, 10150)
    )
    wrapped = GENOME["chr1"][9800:] + GENOME["chr1"][:150]
    want = wrapped if strand == "+" else revcomp(wrapped)
    assert _mod().region_fasta_text(tu, GENOME) == f">TU_7\n{want}\n"


def test_write_tu_regions_files_are_the_one_line_region_text(tmp_path: Path) -> None:
    """``tu_upstream/<tu_id>.fa`` (M2-12 writer) is exactly what the tool gets: one record, one
    sequence line (§7.1.1 gotcha), i.e. ``region_fasta_text``."""
    from gippa.annotate.operons import write_tu_regions

    long_tu = make_tu("TU_9", "chr1", "+", (2000, 2100), (1000, 2050))
    paths = write_tu_regions(tmp_path / "tu_upstream", [*TUS, long_tu], GENOME)
    for tu, path in zip([*TUS, long_tu], paths, strict=True):
        assert path.read_text() == _mod().region_fasta_text(tu, GENOME), tu.tu_id


# ---------------------------------------------------------------- example 3: command


def test_example3_command() -> None:
    from gippa.annotate.promoter_calculator_cmd import promoter_calculator_command

    got = promoter_calculator_command(Path("r.fa"), Path("r.csv"))
    assert got == ["promoter-calculator", "-i", "r.fa", "-o", "r.csv", "-j", "1"]
    assert all(isinstance(a, str) for a in got)


# ---------------------------------------------------------------- parser


def test_parser_converts_the_tool_tss_to_1_based_region_positions(tmp_path: Path) -> None:
    # tool TSS 100 on + = 0-based +1 base → 1-based 101; tool TSS 51 on - = 1-based already
    path = write_pc_csv(tmp_path / "TU_1.csv", [(100, "+", 1234.5), (51, "-", 80.0)], raw=True)
    hits = _mod().parse_pc_csv(path)
    assert hits == [_hit(101, "+", 1234.5), _hit(51, "-", 80.0)]


def test_parser_keeps_file_order_and_values(tmp_path: Path) -> None:
    rows = [(141, "+", 5.0), (101, "+", 10.0), (77, "-", 3.25), (121, "+", 30.0)]
    hits = _mod().parse_pc_csv(write_pc_csv(tmp_path / "x.csv", rows))
    assert [(h.tss, h.strand, h.tx_rate) for h in hits] == rows
    assert all(isinstance(h.tss, int) and isinstance(h.tx_rate, float) for h in hits)


def test_parser_reads_lf_files_like_the_nf_stub(tmp_path: Path) -> None:
    path = write_pc_csv(tmp_path / "x.csv", [(101, "+", 2.0)], crlf=False)
    assert _mod().parse_pc_csv(path) == [_hit(101, "+", 2.0)]


@pytest.mark.parametrize("crlf", [True, False])
def test_header_only_csv_has_no_hits(tmp_path: Path, crlf: bool) -> None:
    """No promoter: the tool writes no CSV and PROMOTER_CALCULATOR writes the header (CRLF)."""
    assert _mod().parse_pc_csv(header_only(tmp_path / "x.csv", crlf=crlf)) == []


def test_empty_csv_has_no_hits(tmp_path: Path) -> None:
    (tmp_path / "x.csv").write_bytes(b"")
    assert _mod().parse_pc_csv(tmp_path / "x.csv") == []


@pytest.mark.parametrize("drop", ["TSS", "strand", "Tx_rate"])
def test_parser_rejects_a_header_without_a_needed_column(tmp_path: Path, drop: str) -> None:
    cols = [c for c in HEADER.split(",") if c != drop]
    (tmp_path / "x.csv").write_text(",".join(cols) + "\r\n")
    with pytest.raises(ValueError, match=drop):
        _mod().parse_pc_csv(tmp_path / "x.csv")


def test_pchit_is_frozen() -> None:
    h = _hit(1, "+", 1.0)
    with pytest.raises((AttributeError, TypeError)):
        h.tss = 2  # type: ignore[misc]


# ---------------------------------------------------------------- examples 4-8: tx_layer, features


def test_example4_plus_tu_hit_maps_to_genome_and_distance() -> None:
    pc = _mod()
    tx = pc.tx_layer(TUS, {"TU_1": [_hit(101, "+", 1234.5)]}, REPS)
    assert list(tx.columns) == TX_COLS
    assert _rows(tx) == {("chr1", "+", 950, 1234.5, "TU_1")}
    f = _feat(pc.tu_tx_features(TUS, tx, replicons=REPS), "TU_1")
    assert f["tx_rate_max"] == 1234.5 and f["tx_rate_sum"] == 1234.5
    assert f["tss_distance_bp"] == 50


def test_example5_minus_tu_hit_maps_to_genome_and_distance() -> None:
    pc = _mod()
    tx = pc.tx_layer(TUS, {"TU_3": [_hit(51, "+", 80.0)]}, REPS)
    assert _rows(tx) == {("chr1", "-", 5749, 80.0, "TU_3")}
    f = _feat(pc.tu_tx_features(TUS, tx, replicons=REPS), "TU_3")
    assert f["tss_distance_bp"] == 150
    assert f["tx_rate_max"] == 80.0 and f["tx_rate_sum"] == 80.0


def test_example6_opposite_strand_hits_are_dropped() -> None:
    pc = _mod()
    tx = pc.tx_layer(TUS, {"TU_1": [_hit(60, "-", 999.0)], "TU_3": [_hit(51, "-", 7.0)]}, REPS)
    assert len(tx) == 0 and list(tx.columns) == TX_COLS
    feats = pc.tu_tx_features(TUS, tx, replicons=REPS)
    assert (feats["tx_rate_max"] == 0).all()


def test_example7_max_sum_and_distance_from_the_max_hit() -> None:
    pc = _mod()
    hits = {"TU_1": [_hit(101, "+", 10.0), _hit(121, "+", 30.0), _hit(141, "+", 5.0)]}
    tx = pc.tx_layer(TUS, hits, REPS)
    assert _rows(tx) == {
        ("chr1", "+", 950, 10.0, "TU_1"),
        ("chr1", "+", 970, 30.0, "TU_1"),
        ("chr1", "+", 990, 5.0, "TU_1"),
    }
    f = _feat(pc.tu_tx_features(TUS, tx, replicons=REPS), "TU_1")
    assert f["tx_rate_max"] == 30.0
    assert f["tx_rate_sum"] == pytest.approx(45.0)
    assert f["tss_distance_bp"] == 30  # 1000 - 970


def test_example8_tu_without_kept_hits_is_zero_zero_na() -> None:
    pc = _mod()
    hits = {"TU_1": [_hit(101, "+", 3.0)], "TU_3": [_hit(51, "-", 9.0)]}  # TU_2 absent: no CSV
    feats = pc.tu_tx_features(TUS, pc.tx_layer(TUS, hits, REPS), replicons=REPS)
    assert list(feats.columns) == FEAT_COLS
    assert list(feats["tu_id"]) == ["TU_1", "TU_2", "TU_3"]
    for tu_id in ("TU_2", "TU_3"):
        f = _feat(feats, tu_id)
        assert f["tx_rate_max"] == 0 and f["tx_rate_sum"] == 0
        assert pd.isna(f["tss_distance_bp"])
    assert not pd.isna(_feat(feats, "TU_1")["tss_distance_bp"])


def test_features_of_an_empty_layer_cover_every_tu() -> None:
    pc = _mod()
    feats = pc.tu_tx_features(TUS, pc.tx_layer(TUS, {}, REPS), replicons=REPS)
    assert list(feats["tu_id"]) == ["TU_1", "TU_2", "TU_3"]
    assert feats["tss_distance_bp"].isna().all()


def test_hits_outside_the_region_are_dropped() -> None:
    pc = _mod()
    hits = {
        "TU_1": [_hit(0, "+", 1.0), _hit(1, "+", 2.0), _hit(200, "+", 3.0), _hit(201, "+", 4.0)]
    }
    tx = pc.tx_layer(TUS, hits, REPS)
    assert _rows(tx) == {("chr1", "+", 850, 2.0, "TU_1"), ("chr1", "+", 1049, 3.0, "TU_1")}


def test_minus_region_ends_map_to_the_region_bounds() -> None:
    pc = _mod()
    hits = {"TU_3": [_hit(1, "+", 1.0), _hit(250, "+", 2.0), _hit(251, "+", 9.0)]}
    tx = pc.tx_layer(TUS, hits, REPS)
    assert _rows(tx) == {("chr1", "-", 5799, 1.0, "TU_3"), ("chr1", "-", 5550, 2.0, "TU_3")}


def test_hits_for_unknown_tus_are_ignored_or_rejected_but_never_mapped() -> None:
    pc = _mod()
    try:
        tx = pc.tx_layer(TUS, {"TU_99": [_hit(5, "+", 1.0)]}, REPS)
    except (KeyError, ValueError):
        return
    assert len(tx) == 0


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_crossing_tss_is_reduced_modulo_the_length(strand: str) -> None:
    """+ TU: S = 100, region [9800, 10150); - TU: start codon first base 9849, region same."""
    pc = _mod()
    if strand == "+":
        tu = make_tu("TU_7", "chr1", "+", (100, 700), (9800, 10150))
        hits = [_hit(101, "+", 1.0), _hit(251, "+", 2.0)]  # genome 9900 and 10050 → 50
        want = {("chr1", "+", 9900, 1.0, "TU_7"), ("chr1", "+", 50, 2.0, "TU_7")}
        dist = {1.0: 200, 2.0: 50}
    else:
        tu = make_tu("TU_7", "chr1", "-", (9250, 9850), (9800, 10150))
        hits = [_hit(101, "+", 1.0), _hit(301, "+", 2.0)]  # genome 10049 → 49 and 9849
        want = {("chr1", "-", 49, 1.0, "TU_7"), ("chr1", "-", 9849, 2.0, "TU_7")}
        dist = {1.0: 200, 2.0: 0}
    tx = pc.tx_layer([tu], {"TU_7": hits}, REPS)
    assert _rows(tx) == want
    for rate, d in dist.items():
        one = tx[tx["tx_rate"] == rate]
        f = _feat(pc.tu_tx_features([tu], one, replicons=REPS), "TU_7")
        assert f["tss_distance_bp"] == d, (strand, rate)


def test_tx_layer_on_a_linear_contig_keeps_positions() -> None:
    pc = _mod()
    tx = pc.tx_layer(TUS, {"TU_1": [_hit(101, "+", 1.0)]}, LIN)
    assert _rows(tx) == {("chr1", "+", 950, 1.0, "TU_1")}


def test_tx_layer_dtypes() -> None:
    pc = _mod()
    tx = pc.tx_layer(TUS, {"TU_1": [_hit(101, "+", 1.0)], "TU_3": [_hit(51, "+", 2.0)]}, REPS)
    assert pd.api.types.is_integer_dtype(tx["tss"])
    assert pd.api.types.is_float_dtype(tx["tx_rate"])


# ---------------------------------------------------------------- tx_rate_between (§9.3 pairs)


def _gene(tag: str, start: int, end: int, strand: str) -> Gene:
    return Gene(tag, "cds", Interval("chr1", start, end, strand), "M", None, "ATG")  # type: ignore[arg-type]


def _tx(rows: list[tuple[str, int, float]]) -> pd.DataFrame:
    return pd.DataFrame([("chr1", s, t, x, "TU_x") for s, t, x in rows], columns=TX_COLS).astype(
        {"tss": "int64", "tx_rate": "float64"}
    )


def test_tx_rate_between_plus_pair() -> None:
    up, down = _gene("a", 100, 400, "+"), _gene("b", 420, 700, "+")
    tx = _tx(
        [("+", 410, 7.0), ("+", 405, 2.0), ("+", 300, 99.0), ("-", 415, 50.0), ("+", 430, 80.0)]
    )
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) == 7.0


def test_tx_rate_between_minus_pair_uses_transcription_order() -> None:
    up, down = _gene("a", 800, 1000, "-"), _gene("b", 500, 700, "-")  # up transcribed first
    tx = _tx(
        [("-", 750, 4.0), ("-", 720, 6.0), ("-", 900, 99.0), ("-", 600, 88.0), ("+", 760, 77.0)]
    )
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) == 6.0


def test_tx_rate_between_none_without_tss() -> None:
    up, down = _gene("a", 100, 400, "+"), _gene("b", 420, 700, "+")
    tx = _tx([("+", 300, 99.0), ("-", 410, 5.0)])
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) is None
    assert _mod().tx_rate_between(up, down, _tx([]), REPS["chr1"]) is None


def test_tx_rate_between_wraps_the_origin() -> None:
    up, down = _gene("a", 9500, 9900, "+"), _gene("b", 50, 300, "+")
    tx = _tx([("+", 9950, 3.0), ("+", 20, 4.0), ("+", 9800, 100.0), ("+", 5000, 60.0)])
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) == 4.0


def test_tx_rate_between_ignores_other_contigs() -> None:
    up, down = _gene("a", 100, 400, "+"), _gene("b", 420, 700, "+")
    tx = _tx([("+", 410, 7.0)])
    tx.loc[len(tx)] = ["chr2", "+", 412, 70.0, "TU_y"]
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) == 7.0


# Window edges (lead decision 5). + pair: up [100, 400), down [420, 700) → window [400, 420).
# - pair: up [800, 1000), down [500, 700) → window [700, 800) in forward coordinates, i.e. from
# 799 (first base after up's stop) down to 700 (base after down's start codon 699).
_PLUS_PAIR = (("a", 100, 400, "+"), ("b", 420, 700, "+"))
_MINUS_PAIR = (("a", 800, 1000, "-"), ("b", 500, 700, "-"))


@pytest.mark.parametrize(
    ("pair", "tss", "kept"),
    [
        (_PLUS_PAIR, 399, False),  # up's last base
        (_PLUS_PAIR, 400, True),  # first base after up's stop
        (_PLUS_PAIR, 419, True),  # last base before down's start codon
        (_PLUS_PAIR, 420, False),  # down's start-codon base
        (_MINUS_PAIR, 800, False),  # up's last base (transcription order)
        (_MINUS_PAIR, 799, True),  # first base after up's stop
        (_MINUS_PAIR, 700, True),  # last base before down's start codon
        (_MINUS_PAIR, 699, False),  # down's start-codon base
    ],
    ids=[
        "+up_last",
        "+first_after_stop",
        "+before_start",
        "+start_codon",
        "-up_last",
        "-first_after_stop",
        "-before_start",
        "-start_codon",
    ],
)
def test_tx_rate_between_window_edges(pair, tss: int, kept: bool) -> None:
    up, down = (_gene(*g) for g in pair)  # type: ignore[arg-type]
    strand = up.interval.strand
    got = _mod().tx_rate_between(up, down, _tx([(strand, tss, 5.0)]), REPS["chr1"])
    assert got == (5.0 if kept else None)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("gap", [0, -1, -50])
def test_tx_rate_between_adjacent_or_overlapping_genes_is_none(strand: str, gap: int) -> None:
    if strand == "+":
        up, down = _gene("a", 100, 400, "+"), _gene("b", 400 + gap, 700, "+")
    else:
        up, down = _gene("a", 400 + gap, 700, "-"), _gene("b", 100, 400, "-")
    tx = _tx([(strand, t, 9.0) for t in range(90, 710)])
    assert _mod().tx_rate_between(up, down, tx, REPS["chr1"]) is None


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tx_rate_between_adjacent_genes_do_not_read_the_layer(strand: str) -> None:
    """Gap 0 returns ``None`` before touching ``tx`` (``span <= 0``, not ``span < 0``)."""
    if strand == "+":
        up, down = _gene("a", 100, 400, "+"), _gene("b", 400, 700, "+")
    else:
        up, down = _gene("a", 400, 700, "-"), _gene("b", 100, 400, "-")
    assert _mod().tx_rate_between(up, down, pd.DataFrame(), REPS["chr1"]) is None


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tx_rate_between_edges_across_the_origin(strand: str) -> None:
    # + : up [9900, 10000), down [30, 300) → window {0..29}; - : mirror image.
    if strand == "+":
        up, down = _gene("a", 9900, 10_000, "+"), _gene("b", 30, 300, "+")
        inside, outside = [0, 29], [9999, 30]
    else:
        up, down = _gene("a", 0, 100, "-"), _gene("b", 9700, 9970, "-")
        inside, outside = [9999, 9970], [0, 9969]
    for t in inside:
        assert _mod().tx_rate_between(up, down, _tx([(strand, t, 2.0)]), REPS["chr1"]) == 2.0, t
    for t in outside:
        assert _mod().tx_rate_between(up, down, _tx([(strand, t, 2.0)]), REPS["chr1"]) is None, t


@given(
    strand=st.sampled_from(["+", "-"]),
    a=st.integers(0, L - 1),
    gap=st.integers(1, 300),
    k=st.integers(-3, 303),
)
@settings(max_examples=200, deadline=None)
def test_tx_rate_between_keeps_exactly_the_gap_bases(strand: str, a: int, gap: int, k: int) -> None:
    """A TSS ``k`` bases past ``up``'s 3' end (transcription order) is kept iff ``0 <= k < gap``."""
    n = 50  # gene length; genes may not cross the origin (Interval), so skip such placements
    if strand == "+":
        us, ds = a, (a + n + gap) % L
        if us + n > L or ds + n > L:
            return
        up, down = _gene("a", us, us + n, "+"), _gene("b", ds, ds + n, "+")
        t = (us + n + k) % L
    else:
        us = a
        de = (us - gap) % L  # down's end (exclusive); its start codon is de - 1
        if us + n > L or de - n < 0:
            return
        up, down = _gene("a", us, us + n, "-"), _gene("b", de - n, de, "-")
        t = (us - 1 - k) % L
    got = _mod().tx_rate_between(up, down, _tx([(strand, t, 1.5)]), REPS["chr1"])
    assert got == (1.5 if 0 <= k < gap else None)


@pytest.mark.parametrize(
    ("d", "length", "expected"),
    [
        (5000, 10_000, 5000),  # exactly half an even contig: positive (upstream)
        (-5000, 10_000, 5000),  # -length/2 maps to +length/2
        (5001, 10_000, -4999),
        (4999, 10_000, 4999),
        (4999, 9999, 4999),  # odd contig: (-4999.5, 4999.5]
        (5000, 9999, -4999),
        (0, 10_000, 0),
        (10_000, 10_000, 0),
        (5000, None, 5000),
        (-7, None, -7),
    ],
)
def test_signed_at_half_the_contig_length(d: int, length: int | None, expected: int) -> None:
    """``_signed`` maps to ``(-length/2, length/2]``; the tie at exactly half is positive."""
    assert _mod()._signed(d, length) == expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tss_distance_at_exactly_half_the_contig_is_positive(strand: str) -> None:
    pc = _mod()
    if strand == "+":  # start codon 5000, TSS 0 → 5000 bp upstream either way round
        tu = make_tu("TU_h", "chr1", "+", (5000, 5300), (4900, 5050))
        tx = _tx([("+", 0, 3.0)])
    else:  # start codon 4999, TSS 9999 → 5000 bp either way round
        tu = make_tu("TU_h", "chr1", "-", (4700, 5000), (4950, 5100))
        tx = _tx([("-", 9999, 3.0)])
    tx["tu_id"] = "TU_h"
    f = _feat(pc.tu_tx_features([tu], tx, replicons=REPS), "TU_h")
    assert f["tss_distance_bp"] == 5000


# ---------------------------------------------------------------- properties


@st.composite
def _region_case(draw):
    size = draw(st.integers(1000, 3000))
    r = random.Random(draw(st.integers(0, 2**32 - 1)))
    seq = "".join(r.choice("ACGT") for _ in range(size))
    strand = draw(st.sampled_from(["+", "-"]))
    a = draw(st.integers(0, size - 1))
    tss = draw(st.integers(1, 350))
    return size, seq, strand, a, tss


def _geometry_tu(size: int, strand: str, a: int):
    """A TU as ``tu_upstream`` builds it with bp (-300, 50): region [a, a + 350)."""
    if strand == "+":
        s = (a + 300) % size  # start codon S
        return make_tu("TU_1", "c", "+", (s, s + 60), (a, a + 350))
    e = a + 50  # start codon first base e - 1
    lo = (e - 60) % size
    return make_tu("TU_1", "c", "-", (lo, lo + 60), (a, a + 350))


@settings(max_examples=150, deadline=None)
@given(_region_case())
def test_region_text_and_tss_mapping_agree(case) -> None:
    """The +1 base the tool reports in the region file is the base at the genome TSS (both
    strands, across the origin), and the distance to the start codon is 301 - tss."""
    size, seq, strand, a, tss = case
    pc = _mod()
    reps = {"c": Replicon("c", size, "circular")}
    tu = _geometry_tu(size, strand, a)
    region = pc.region_fasta_text(tu, {"c": seq}).splitlines()[1]
    assert len(region) == 350
    tx = pc.tx_layer([tu], {"TU_1": [pc.PcHit(tss=tss, strand="+", tx_rate=1.0)]}, reps)
    assert len(tx) == 1
    g = int(tx["tss"].iloc[0])
    assert 0 <= g < size
    base = seq[g] if strand == "+" else revcomp(seq[g])
    assert base == region[tss - 1]
    f = _feat(pc.tu_tx_features([tu], tx, replicons=reps), "TU_1")
    assert f["tss_distance_bp"] == 301 - tss


@settings(max_examples=100, deadline=None)
@given(
    st.lists(
        st.tuples(st.integers(-5, 260), st.sampled_from(["+", "-"]), st.floats(0.0, 1e5)),
        max_size=25,
    )
)
def test_features_are_max_and_sum_of_the_kept_hits(raw) -> None:
    pc = _mod()
    hits = [pc.PcHit(tss=t, strand=s, tx_rate=x) for t, s, x in raw]
    tx = pc.tx_layer(TUS, {"TU_3": hits}, REPS)
    kept = [h for h in hits if h.strand == "+" and 1 <= h.tss <= 250]
    assert len(tx) == len(kept)
    f = _feat(pc.tu_tx_features(TUS, tx, replicons=REPS), "TU_3")
    rates = [h.tx_rate for h in kept]
    assert f["tx_rate_max"] == (max(rates) if rates else 0)
    assert math.isclose(float(f["tx_rate_sum"]), sum(rates), rel_tol=1e-9, abs_tol=1e-9)
    assert f["tx_rate_max"] <= f["tx_rate_sum"] + 1e-9
    if rates:
        best = max(rates)
        dists = {(5800 - h.tss) - 5599 for h in kept if h.tx_rate == best}
        assert f["tss_distance_bp"] in dists
    else:
        assert pd.isna(f["tss_distance_bp"])


@settings(max_examples=100, deadline=None)
@given(st.lists(st.tuples(st.integers(1, 200), st.floats(0.1, 100.0)), min_size=1, max_size=12))
def test_parse_then_layer_round_trip(tmp_path_factory, raw) -> None:
    """CSV written in the tool's convention → parse → tx_layer gives a + 0-based genome TSS."""
    pc = _mod()
    path = write_pc_csv(tmp_path_factory.mktemp("pc") / "TU_1.csv", [(t, "+", x) for t, x in raw])
    tx = pc.tx_layer(TUS, {"TU_1": pc.parse_pc_csv(path)}, REPS)
    assert sorted(tx["tss"].tolist()) == sorted(850 + t - 1 for t, _ in raw)
