"""Plot support tables from run inputs (M9-03; docs/specs/M9-03-plot-support-run-inputs.md).

``gippa.viz.support.replicon_landmarks`` (oriC and terminus per replicon, V1 and V2) and
``gippa.viz.support.gc_skew`` (windowed GC skew, V1): spec examples 1-7, the degenerate-case
table, and Hypothesis properties -- the oriC agrees with ``gippa.context.oriter.oric_position`` on
gippasim layouts (acceptance criterion), the terminus is the antipode, windows tile each contig,
the cumulative skew is the running sum with NA windows carrying, case does not matter.

Interface assumed beyond the spec text (M9-03 test author):

- Both functions return the table **as written** (the "Tables" section): 1-based ``oric_position``
  / ``ter_position`` and 1-based closed ``start`` / ``end``, so ``write_support`` is a plain
  writer and ``RunView.support`` (which converts nothing) hands the plots the same frame.
- Column order is the spec's: landmarks ``contig, length_bp, topology, oric_position,
  ter_position, status``; GC skew ``contig, start, end, gc_skew, cumulative_gc_skew``.
- NA cells are anything ``pandas.isna`` accepts (NaN or ``pd.NA``).
- GC-skew rows follow the replicons' (C1 = FASTA) order, then ``start``; every replicon in the
  replicons frame gets windows, linear ones included (the spec excludes none).
- Where ``oric_position`` returns a non-integral midpoint (circular mean of several calls, or an
  odd-length call) the rounding to an integer position is not pinned: only |pos - (o + 1)| < 1.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

LANDMARK_COLUMNS = ["contig", "length_bp", "topology", "oric_position", "ter_position", "status"]
SKEW_COLUMNS = ["contig", "start", "end", "gc_skew", "cumulative_gc_skew"]
FEATURE_COLUMNS = [
    "feature_id", "type", "contig", "start", "end", "strand",
    "gene", "product", "amino_acid", "anticodon",
]  # fmt: skip

SETTINGS = settings(
    max_examples=30,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)


def _support() -> Any:
    from gippa.viz import support

    return support


def replicons(*rows: tuple[str, int, str, str]) -> pd.DataFrame:
    """C1 ``replicons`` from ``(contig, length, topology, kind)`` rows."""
    return pd.DataFrame(rows, columns=["contig", "length", "topology", "kind"]).astype(
        {"length": "int64"}
    )


def features(*rows: tuple[str, str, int, int]) -> pd.DataFrame:
    """C15 ``features`` from ``(type, contig, start, end)`` rows (0-based half-open)."""
    out = pd.DataFrame(
        [
            (f"F{i + 1:04d}", t, c, int(s), int(e), ".", None, None, None, None)
            for i, (t, c, s, e) in enumerate(rows)
        ],
        columns=FEATURE_COLUMNS,
    )
    return out.astype({"start": "int64", "end": "int64"})


def _row(df: pd.DataFrame, contig: str) -> pd.Series:
    sel = df[df["contig"].astype(str) == contig]
    assert len(sel) == 1, (contig, df)
    return sel.iloc[0]


def _na(v: object) -> bool:
    return bool(pd.isna(v))


# ================================================================ replicon landmarks


def test_example_1_origin_crossing_oric_gives_position_1_and_terminus_50001() -> None:
    reps = replicons(("c1", 100_000, "circular", "chromosome"))
    feats = features(("oric", "c1", 99_900, 100_100))
    out = _support().replicon_landmarks(reps, feats)
    assert list(out.columns) == LANDMARK_COLUMNS
    row = _row(out, "c1")
    assert int(row["length_bp"]) == 100_000
    assert row["topology"] == "circular"
    assert int(row["oric_position"]) == 1
    assert int(row["ter_position"]) == 50_001
    assert row["status"] == "called"


def test_example_2_circular_without_oric_is_no_oric() -> None:
    reps = replicons(("c2", 8_000, "circular", "plasmid"))
    out = _support().replicon_landmarks(reps, features(("crispr", "c2", 10, 50)))
    row = _row(out, "c2")
    assert _na(row["oric_position"]) and _na(row["ter_position"])
    assert row["status"] == "no_oric"
    assert int(row["length_bp"]) == 8_000


def test_example_3_linear_contig_is_linear_even_with_an_oric_call() -> None:
    reps = replicons(("lin", 5_000, "linear", "unknown"))
    out = _support().replicon_landmarks(reps, features(("oric", "lin", 1_000, 1_232)))
    row = _row(out, "lin")
    assert _na(row["oric_position"]) and _na(row["ter_position"])
    assert row["status"] == "linear"
    assert row["topology"] == "linear"


def test_no_features_table_behaves_like_no_oric_call() -> None:
    reps = replicons(("c1", 10_000, "circular", "chromosome"), ("lin", 900, "linear", "unknown"))
    out = _support().replicon_landmarks(reps, None)
    assert list(out["status"]) == ["no_oric", "linear"]
    assert out["oric_position"].isna().all() and out["ter_position"].isna().all()


def test_one_row_per_replicon_in_replicon_order() -> None:
    reps = replicons(
        ("zeta", 20_000, "circular", "chromosome"),
        ("alpha", 4_000, "circular", "plasmid"),
        ("mid", 3_000, "linear", "unknown"),
    )
    feats = features(("oric", "zeta", 500, 732), ("oric", "alpha", 100, 200))
    out = _support().replicon_landmarks(reps, feats)
    assert list(out["contig"]) == ["zeta", "alpha", "mid"]
    assert [int(x) for x in out["length_bp"]] == [20_000, 4_000, 3_000]
    assert list(out["status"]) == ["called", "called", "linear"]
    # zeta: midpoint 616 (0-based) -> 617; terminus (616 + 10000) mod 20000 + 1
    assert int(_row(out, "zeta")["oric_position"]) == 617
    assert int(_row(out, "zeta")["ter_position"]) == 10_617
    # alpha: midpoint 150 -> 151; terminus 150 + 2000 + 1
    assert int(_row(out, "alpha")["ter_position"]) == 2_151


def test_terminus_wraps_modulo_length() -> None:
    reps = replicons(("c1", 10_000, "circular", "chromosome"))
    out = _support().replicon_landmarks(reps, features(("oric", "c1", 7_900, 8_100)))
    row = _row(out, "c1")
    assert int(row["oric_position"]) == 8_001
    assert int(row["ter_position"]) == (8_000 + 5_000) % 10_000 + 1 == 3_001


def test_terminus_uses_integer_half_length_on_odd_lengths() -> None:
    reps = replicons(("c1", 10_001, "circular", "chromosome"))
    out = _support().replicon_landmarks(reps, features(("oric", "c1", 900, 1_100)))
    row = _row(out, "c1")
    assert int(row["oric_position"]) == 1_001
    assert int(row["ter_position"]) == (1_000 + 10_001 // 2) % 10_001 + 1 == 6_001


def test_several_close_oric_calls_give_their_circular_mean() -> None:
    # midpoints 99_000 and 1_000 on L = 100_000: 2_000 apart across the origin (< 0.05 L);
    # circular mean 0, not the arithmetic mean 50_000
    reps = replicons(("c1", 100_000, "circular", "chromosome"))
    feats = features(("oric", "c1", 98_900, 99_100), ("oric", "c1", 900, 1_100))
    row = _row(_support().replicon_landmarks(reps, feats), "c1")
    assert row["status"] == "called"
    assert int(row["oric_position"]) == 1  # 0 up to float error, nearest base (decision 4a)
    assert abs(int(row["ter_position"]) - 50_001) <= 1


def test_half_base_midpoint_rounds_up() -> None:
    # [100, 301) has midpoint 200.5: nearest base, ties up, gives 0-based 201 (decision 4a)
    reps = replicons(("c1", 1_000, "circular", "chromosome"))
    row = _row(_support().replicon_landmarks(reps, features(("oric", "c1", 100, 301))), "c1")
    assert (int(row["oric_position"]), int(row["ter_position"])) == (202, 702)


def test_quarter_base_midpoint_rounds_down() -> None:
    # midpoints 200 and 200.5 average to about 200.25, which rounds to 0-based 200, not up
    reps = replicons(("c1", 1_000, "circular", "chromosome"))
    feats = features(("oric", "c1", 100, 300), ("oric", "c1", 100, 301))
    row = _row(_support().replicon_landmarks(reps, feats), "c1")
    assert (int(row["oric_position"]), int(row["ter_position"])) == (201, 701)


def test_several_distant_oric_calls_give_no_oric() -> None:
    reps = replicons(("c1", 100_000, "circular", "chromosome"))
    feats = features(("oric", "c1", 9_900, 10_100), ("oric", "c1", 59_900, 60_100))
    row = _row(_support().replicon_landmarks(reps, feats), "c1")
    assert row["status"] == "no_oric"
    assert _na(row["oric_position"]) and _na(row["ter_position"])


def test_oric_on_another_contig_does_not_leak() -> None:
    reps = replicons(("c1", 50_000, "circular", "chromosome"), ("p1", 6_000, "circular", "plasmid"))
    out = _support().replicon_landmarks(reps, features(("oric", "c1", 100, 332)))
    assert _row(out, "c1")["status"] == "called"
    assert _row(out, "p1")["status"] == "no_oric"


def test_example_7_sim100k_replicons_two_rows_in_order() -> None:
    """The val-mvp-2026-10-22 replicons (chr1 100,000; plasmid_1 8,000, both circular), from the
    committed sim100k genome and recorded Bakta JSON, plus that run's chr1 oriC call."""
    from helpers import oriter_data as od

    recorded_feats, reps = od.recorded("sim100k")
    assert list(reps["contig"]) == ["chr1", "plasmid_1"]
    feats = pd.concat(
        [recorded_feats, features(("oric", "chr1", 12_445, 12_677))], ignore_index=True
    )
    out = _support().replicon_landmarks(reps, feats)
    assert list(out["contig"]) == ["chr1", "plasmid_1"]
    assert [int(x) for x in out["length_bp"]] == [100_000, 8_000]
    assert list(out["topology"]) == ["circular", "circular"]
    chr1 = _row(out, "chr1")
    assert (chr1["status"], int(chr1["oric_position"]), int(chr1["ter_position"])) == (
        "called",
        12_562,
        62_562,
    )
    assert _row(out, "plasmid_1")["status"] == "no_oric"


def test_inputs_are_not_mutated() -> None:
    reps = replicons(("c1", 100_000, "circular", "chromosome"))
    feats = features(("oric", "c1", 99_900, 100_100))
    r0, f0 = reps.copy(), feats.copy()
    _support().replicon_landmarks(reps, feats)
    pd.testing.assert_frame_equal(reps, r0)
    pd.testing.assert_frame_equal(feats, f0)


# ---------------------------------------------------------------- properties (gippasim)


@st.composite
def sim_layout(draw: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(replicons, features) of a small gippasim layout, with extra oriC calls sometimes added on
    chr1 (close to or far from the planted one) and on a plasmid."""
    from gippasim.layout import LayoutConfig, simulate_layout

    seed = draw(st.integers(0, 2**31 - 1))
    oric_frac = draw(st.one_of(st.none(), st.floats(0.0, 0.999)))
    cfg = LayoutConfig(
        n_genes_chromosome=draw(st.integers(20, 80)),
        plasmids=tuple(draw(st.lists(st.integers(1, 6), max_size=2))),
        linear_contigs=tuple(draw(st.lists(st.integers(1, 4), max_size=2))),
        regions=(),
        origin_crossing=draw(st.booleans()),
        oric_frac=oric_frac,
    )
    try:
        layout = simulate_layout(cfg, np.random.default_rng(seed))
    except ValueError:  # no intergenic gap fits an oriC (gippasim lead decision 14)
        assume(False)
    reps = layout.replicons.reset_index(drop=True)
    feats = layout.features.reset_index(drop=True)
    extra: list[tuple[str, str, int, int]] = []
    length = dict(zip(reps["contig"], reps["length"], strict=True))
    for contig, n in length.items():
        for _ in range(draw(st.integers(0, 2 if contig == "chr1" else 1))):
            s = draw(st.integers(0, int(n) - 1))
            w = draw(st.integers(2, 400))
            extra.append(("oric", str(contig), s, s + w))  # may be unwrapped past n (§6.5)
    if extra:
        feats = pd.concat([feats, features(*extra)], ignore_index=True)
    return reps, feats


@SETTINGS
@given(data=sim_layout())
def test_property_oric_agrees_with_oriter_on_gippasim_layouts(
    data: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    from gippa.context import oriter

    reps, feats = data
    out = _support().replicon_landmarks(reps, feats)
    assert list(out.columns) == LANDMARK_COLUMNS
    assert list(out["contig"]) == list(reps["contig"])
    assert [int(x) for x in out["length_bp"]] == [int(x) for x in reps["length"]]
    assert list(out["topology"]) == list(reps["topology"])
    for rep, row in zip(reps.itertuples(), out.itertuples(), strict=True):
        n = int(rep.length)
        if rep.topology != "circular":
            assert row.status == "linear"
            assert _na(row.oric_position) and _na(row.ter_position)
            continue
        o = oriter.oric_position(feats, str(rep.contig), n)
        if o is None:
            assert row.status == "no_oric"
            assert _na(row.oric_position) and _na(row.ter_position)
            continue
        assert row.status == "called"
        pos, ter = int(row.oric_position), int(row.ter_position)
        assert 1 <= pos <= n and 1 <= ter <= n
        assert pos - 1 == math.floor(o + 0.5) % n  # nearest base, ties up (decision 4a)
        assert ter == ((pos - 1) + n // 2) % n + 1


@SETTINGS
@given(
    length=st.integers(10, 500_000),
    mid=st.integers(0, 499_999),
    half=st.integers(1, 200),
)
def test_property_single_call_terminus_is_the_antipode(length: int, mid: int, half: int) -> None:
    mid %= length
    half = min(half, max(1, length // 50))  # keep the call short enough to be one call
    reps = replicons(("c1", length, "circular", "chromosome"))
    start = mid - half if mid >= half else mid - half + length  # unwrapped past L (§6.5)
    out = _support().replicon_landmarks(reps, features(("oric", "c1", start, start + 2 * half)))
    row = _row(out, "c1")
    pos, ter = int(row["oric_position"]), int(row["ter_position"])
    assert pos == mid % length + 1
    assert (ter - pos) % length == length // 2


# ================================================================ GC skew


def _skew(genome: Mapping[str, str], reps: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    out = _support().gc_skew(genome, reps, window_bp)
    assert list(out.columns) == SKEW_COLUMNS
    return out


def _circ(name: str, seq: str) -> tuple[str, int, str, str]:
    return (name, len(seq), "circular", "chromosome")


def test_example_4_one_window_skew_one_half() -> None:
    seq = "GGGC" + "A" * 9_996
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 10_000)
    assert len(out) == 1
    row = out.iloc[0]
    assert (row["contig"], int(row["start"]), int(row["end"])) == ("c1", 1, 10_000)
    assert row["gc_skew"] == pytest.approx(0.5)
    assert row["cumulative_gc_skew"] == pytest.approx(0.5)


def test_example_5_three_windows_last_one_shorter() -> None:
    seq = "ACGT" * 6_250  # 25,000 bp
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 10_000)
    assert [(int(s), int(e)) for s, e in zip(out["start"], out["end"], strict=True)] == [
        (1, 10_000),
        (10_001, 20_000),
        (20_001, 25_000),
    ]
    assert list(out["contig"]) == ["c1"] * 3
    assert np.allclose(out["gc_skew"].astype(float), 0.0)


def test_example_6_na_window_carries_the_cumulative_value() -> None:
    # window 8: G6 C2 -> 0.5; G3 C5 -> -0.25; no G/C -> NA
    seq = "GGGGGGCC" + "GGGCCCCC" + "AATTAATT"
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 8)
    skew = list(out["gc_skew"])
    assert skew[0] == pytest.approx(0.5) and skew[1] == pytest.approx(-0.25)
    assert _na(skew[2])
    assert [float(x) for x in out["cumulative_gc_skew"]] == pytest.approx([0.5, 0.25, 0.25])


def test_window_of_only_n_is_na_and_carries() -> None:
    seq = "GGGC" + "N" * 4 + "CCCG"
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 4)
    assert len(out) == 3
    assert out["gc_skew"].iloc[0] == pytest.approx(0.5)
    assert _na(out["gc_skew"].iloc[1])
    assert out["gc_skew"].iloc[2] == pytest.approx(-0.5)
    assert [float(x) for x in out["cumulative_gc_skew"]] == pytest.approx([0.5, 0.5, 0.0])


def test_contig_shorter_than_window_is_one_window() -> None:
    seq = "G" * 3 + "C" + "A" * 3_996  # 4,000 bp
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 10_000)
    assert [(int(s), int(e)) for s, e in zip(out["start"], out["end"], strict=True)] == [(1, 4_000)]
    assert out["gc_skew"].iloc[0] == pytest.approx(0.5)


def test_lowercase_is_counted_case_insensitively() -> None:
    upper = "GGGC" + "A" * 9_996
    lower = upper.lower()
    reps = replicons(_circ("c1", upper))
    a, b = _skew({"c1": upper}, reps, 10_000), _skew({"c1": lower}, reps, 10_000)
    pd.testing.assert_frame_equal(a, b)
    mixed = "GgGc" + "a" * 9_996
    assert _skew({"c1": mixed}, reps, 10_000)["gc_skew"].iloc[0] == pytest.approx(0.5)


def test_windows_do_not_wrap_on_circular_contigs() -> None:
    seq = "G" * 15  # 15 bp circular, window 10: [1, 10], [11, 15], never 11..5 across the origin
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 10)
    assert [(int(s), int(e)) for s, e in zip(out["start"], out["end"], strict=True)] == [
        (1, 10),
        (11, 15),
    ]


def test_cumulative_restarts_on_each_contig_and_rows_follow_replicon_order() -> None:
    a = "GGGG" * 2  # skew 1, 1
    b = "CCCC" * 3  # skew -1, -1, -1
    lin = "GGCC"  # skew 0
    reps = replicons(
        ("b", len(b), "circular", "chromosome"),
        ("a", len(a), "circular", "plasmid"),
        ("lin", len(lin), "linear", "unknown"),
    )
    out = _skew({"a": a, "b": b, "lin": lin}, reps, 4)
    assert list(out["contig"]) == ["b", "b", "b", "a", "a", "lin"]
    assert [int(s) for s in out["start"]] == [1, 5, 9, 1, 5, 1]
    assert [float(x) for x in out["cumulative_gc_skew"]] == pytest.approx(
        [-1.0, -2.0, -3.0, 1.0, 2.0, 0.0]
    )


def test_ambiguity_codes_are_not_counted_as_g_or_c() -> None:
    seq = "GGGCSSSS"  # S (G or C) is not counted
    out = _skew({"c1": seq}, replicons(_circ("c1", seq)), 8)
    assert out["gc_skew"].iloc[0] == pytest.approx(0.5)


# ---------------------------------------------------------------- properties


def _oracle(seq: str, w: int) -> list[tuple[int, int, float]]:
    """1-based closed windows and their skew, NaN when G + C = 0 (spec "Tables")."""
    out = []
    s = seq.upper()
    for i in range(0, len(s), w):
        chunk = s[i : i + w]
        g, c = chunk.count("G"), chunk.count("C")
        out.append((i + 1, i + len(chunk), (g - c) / (g + c) if g + c else math.nan))
    return out


_SEQ = st.text(alphabet="ACGTNacgtn", min_size=1, max_size=400)


@SETTINGS
@given(seqs=st.lists(_SEQ, min_size=1, max_size=3), window=st.integers(1, 120))
def test_property_windows_tile_and_skew_matches_the_oracle(
    seqs: Sequence[str], window: int
) -> None:
    names = [f"c{i}" for i in range(len(seqs))]
    reps = replicons(
        *((n, len(s), "circular", "chromosome") for n, s in zip(names, seqs, strict=True))
    )
    out = _skew(dict(zip(names, seqs, strict=True)), reps, window)
    for name, seq in zip(names, seqs, strict=True):
        sub = out[out["contig"] == name]
        want = _oracle(seq, window)
        assert len(sub) == len(want) == math.ceil(len(seq) / window)
        starts, ends = [int(x) for x in sub["start"]], [int(x) for x in sub["end"]]
        assert starts[0] == 1 and ends[-1] == len(seq)
        assert all(e - s + 1 == window for s, e in zip(starts[:-1], ends[:-1], strict=True))
        assert all(b == e + 1 for e, b in zip(ends[:-1], starts[1:], strict=True))
        running = 0.0
        assert list(zip(starts, ends, strict=True)) == [(ws, we) for ws, we, _ in want]
        for (_, _, wk), got, cum in zip(
            want, sub["gc_skew"], sub["cumulative_gc_skew"], strict=True
        ):
            if math.isnan(wk):
                assert _na(got)
            else:
                assert float(got) == pytest.approx(wk)
                assert -1.0 <= float(got) <= 1.0
                running += wk
            assert float(cum) == pytest.approx(running, abs=1e-9)


@SETTINGS
@given(seq=st.text(alphabet="ACGTN", min_size=1, max_size=300), window=st.integers(1, 60))
def test_property_swapping_g_and_c_negates_the_skew(seq: str, window: int) -> None:
    swapped = seq.translate(str.maketrans("GC", "CG"))
    reps = replicons(("c1", len(seq), "circular", "chromosome"))
    a, b = _skew({"c1": seq}, reps, window), _skew({"c1": swapped}, reps, window)
    sa, sb = a["gc_skew"].astype(float).to_numpy(), b["gc_skew"].astype(float).to_numpy()
    assert np.array_equal(np.isnan(sa), np.isnan(sb))
    assert np.allclose(sa[~np.isnan(sa)], -sb[~np.isnan(sb)])


@SETTINGS
@given(seq=st.text(alphabet="ACGTNacgtn", min_size=1, max_size=300), window=st.integers(1, 60))
def test_property_gc_skew_is_deterministic(seq: str, window: int) -> None:
    reps = replicons(("c1", len(seq), "linear", "unknown"))
    pd.testing.assert_frame_equal(
        _skew({"c1": seq}, reps, window), _skew({"c1": seq}, reps, window)
    )
