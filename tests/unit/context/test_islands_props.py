"""Hypothesis properties of the M5-06 island step (docs/specs/M5-06-islands.md "Acceptance").

- TU scores (and every number of the table) are invariant to the order of rows in the gene table.
- Relabelling TU IDs changes no number; only ``tu_first`` / ``tu_last`` follow the labels.
- An island never contains a TU outside its replicon, its TUs are consecutive in TU order and its
  span is the best window: ``min_tus <= n_tus <= max_tus``.
- On a linear replicon no island wraps (``start < end <= length``, TU indices increasing).
- ``effect`` is the mean detrended score of the reported window and ``p_value`` lies in (0, 1].
- Same seed, same islands frame and gene map.

The frames are small (12-40 single- or multi-gene TUs) with a 30 % LOESS span, so they are scanned
(decision 41: the span must give k >= 4) and fast; ``stats.alpha.context = 1`` makes every
non-zero window significant so that islands exist to check.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as idat
from islands_data import Contig

HEALTH = [HealthCheck.too_slow, HealthCheck.data_too_large]
SETTINGS = settings(max_examples=25, deadline=None, suppress_health_check=HEALTH)


def _cfg(**kw: Any) -> Any:
    return idat.config(bootstrap=19, permutations=9, span_frac=0.30, alpha=1.0, **kw)


def _run(
    genes: pd.DataFrame, tus: pd.DataFrame, reps: pd.DataFrame, cfg: Any, seed: int = 3
) -> tuple[pd.DataFrame, Any]:
    from giap.rng import for_step

    return idat.step(
        genes, idat.build_series(genes, tus, reps, cfg), cfg, for_step(idat.STEP, seed)
    )


@st.composite
def contigs(draw: Any, circular: bool | None = None) -> list[Contig]:
    """One replicon of 12-40 genes grouped into TUs, with detection, p_detect and abundance."""
    n = draw(st.integers(14, 40))
    d = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    a = draw(st.lists(st.floats(5.0, 25.0, allow_nan=False, width=32), min_size=n, max_size=n))
    sizes: list[int] = []
    while sum(sizes) < n:
        sizes.append(draw(st.integers(1, 3)))
    sizes[-1] -= sum(sizes) - n
    if sizes[-1] == 0:
        sizes.pop()
    circ = draw(st.booleans()) if circular is None else circular
    return [
        Contig(
            name="chr1",
            d=[int(v) for v in d],
            p=[0.4] * n,
            a=a,
            circular=circ,
            tu=idat.tu_labels(sizes, "chr1"),
        )
    ]


#: columns that may not change when the input is reordered or the TUs relabelled
NUMBERS = ("start", "end", "n_tus", "n_genes", "scan_s", "block_len", "effect", "ci_low",
           "ci_high", "p_value", "n_perm")  # fmt: skip
LABELS = ("direction", "contig", "island_id", "method", "null_model")


def assert_same_numbers(a: pd.DataFrame, b: pd.DataFrame) -> None:
    assert len(a) == len(b)
    for c in NUMBERS:
        assert np.allclose(a[c].to_numpy(float), b[c].to_numpy(float), equal_nan=True), c
    for c in LABELS:
        assert list(a[c].astype(str)) == list(b[c].astype(str)), c


@SETTINGS
@given(contigs(), st.integers(0, 2**31 - 1))
def test_gene_row_order_does_not_matter(cs: list[Contig], seed: int) -> None:
    genes, tus, reps = idat.frames(cs)
    shuffled = genes.sample(frac=1.0, random_state=seed % 2**31).reset_index(drop=True)
    cfg = _cfg()
    a, sa = _run(genes, tus, reps, cfg)
    b, sb = _run(shuffled, tus, reps, cfg)
    assert np.allclose(sa.series[0].score, sb.series[0].score, equal_nan=True)
    assert np.allclose(sa.series[0].detrended, sb.series[0].detrended, equal_nan=True)
    assert sa.series[0].tu_ids == sb.series[0].tu_ids
    assert_same_numbers(a, b)
    assert [list(g) for g in a["genes"]] == [list(g) for g in b["genes"]]


@SETTINGS
@given(contigs())
def test_relabelling_tus_changes_no_number(cs: list[Contig]) -> None:
    genes, tus, reps = idat.frames(cs)
    ren = {t: f"X{i:03d}" for i, t in enumerate(pd.unique(genes["tu_id"].astype(str)))}
    g2 = genes.copy()
    g2["tu_id"] = [ren[str(t)] for t in g2["tu_id"]]
    tus2 = idat.tu_frame(g2)
    cfg = _cfg()
    a, _ = _run(genes, tus, reps, cfg)
    b, _ = _run(g2, tus2, reps, cfg)
    assert_same_numbers(a, b)
    assert [ren[t] for t in a["tu_first"]] == list(b["tu_first"])
    assert [ren[t] for t in a["tu_last"]] == list(b["tu_last"])


@SETTINGS
@given(contigs())
def test_islands_stay_inside_their_replicon(cs: list[Contig]) -> None:
    genes, tus, reps = idat.frames(cs)
    cfg = _cfg()
    out, scan = _run(genes, tus, reps, cfg)
    ts = scan.series[0]
    known = set(ts.tu_ids)
    wmin, wmax = cfg.context.islands.min_tus, cfg.context.islands.max_tus
    for _, r in out.iterrows():
        assert r["contig"] == ts.contig
        assert r["tu_first"] in known and r["tu_last"] in known
        assert wmin <= int(r["n_tus"]) <= min(wmax, len(known) - (1 if ts.circular else 0))
        assert int(r["n_genes"]) == len(r["genes"]) >= int(r["n_tus"])
        assert set(r["genes"]) <= set(genes.loc[genes["contig"] == ts.contig, "locus_tag"])
        first = ts.tu_ids.index(r["tu_first"])
        pos = [(first + k) % len(ts.tu_ids) for k in range(int(r["n_tus"]))]
        assert ts.tu_ids[pos[-1]] == r["tu_last"]
        assert np.isclose(float(r["effect"]), float(np.mean([ts.detrended[j] for j in pos])))
        assert 0.0 < float(r["p_value"]) <= 1.0
        assert float(r["ci_low"]) <= float(r["ci_high"])


@SETTINGS
@given(contigs(circular=False))
def test_no_island_wraps_on_a_linear_replicon(cs: list[Contig]) -> None:
    genes, tus, reps = idat.frames(cs)
    length = int(reps["length"].iloc[0])
    out, scan = _run(genes, tus, reps, _cfg())
    ts = scan.series[0]
    assert not ts.circular
    for _, r in out.iterrows():
        assert ts.tu_ids.index(r["tu_first"]) <= ts.tu_ids.index(r["tu_last"])
        assert 0 <= int(r["start"]) < int(r["end"]) <= length
        assert int(r["n_tus"]) == ts.tu_ids.index(r["tu_last"]) - ts.tu_ids.index(r["tu_first"]) + 1


@settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(contigs(), st.integers(0, 2**31 - 1))
def test_same_seed_same_output(cs: list[Contig], seed: int) -> None:
    from giap.context import islands as M

    genes, tus, reps = idat.frames(cs)
    cfg = _cfg()
    a, _ = _run(genes, tus, reps, cfg, seed=seed)
    b, _ = _run(genes, tus, reps, cfg, seed=seed)
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(M.gene_island_map(a), M.gene_island_map(b))
