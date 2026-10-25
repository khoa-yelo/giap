"""``readthrough.consumed_sites`` must not scan every site for every row (M7-16b rule 15; code
review of PR #335, must-fix 1: 156 s for 100k rows against 2,300 sites per strand).

The check counts work instead of timing it: the sites are a ``StopSite`` subclass that counts
reads of ``codon_start``. A full scan reads it about ``rows x sites`` times; a bisect on the
sorted starts (``[start, end)`` plus ``[start - n, end - n)`` for the origin copy) reads it
``O(sites + rows log sites)`` times. The bound below is ``8 x sites + 100 x rows`` (about
208,000), far above any bisect (``log2(1000) ~ 10`` reads per bound, four bounds per row) and far
below the scan (about ``rows x sites`` = 2,000,000; 3.8 million measured). The answer is checked
against a brute-force oracle of rule 15 on the same rows, origin-crossing ones included.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from giap.mapping.readthrough import StopSite

N = 300_000
N_SITES, N_ROWS = 1000, 2000
_READS = [0]


class _CountingSite(StopSite):
    __slots__ = ()

    def __getattribute__(self, name: str) -> Any:
        if name == "codon_start":
            _READS[0] += 1
        return super().__getattribute__(name)


def _world(seed: int = 335) -> tuple[list[StopSite], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    starts = np.sort(rng.choice(N - 3, size=N_SITES, replace=False))
    sites: list[StopSite] = [
        _CountingSite("c1", "+", int(s), (f"G{k}",), "TGA", "terminal")
        for k, s in enumerate(starts)
    ]
    rows = []
    for k in range(N_ROWS):
        if k % 4 == 0:  # a window around a site, in frame or not
            s = int(starts[rng.integers(N_SITES)])
            a = s - int(rng.integers(0, 30))
        elif k % 4 == 1:  # across the origin
            a = N - int(rng.integers(1, 40))
        else:
            a = int(rng.integers(0, N))
        a %= N
        length = 3 * int(rng.integers(7, 30))
        rows.append(
            {
                "search": "decoy" if k % 50 == 0 else "target",
                "contig": "c1",
                "strand": "+",
                "start": a,
                "end": a + length,
                "blocks": "1-2;3-4" if k % 37 == 0 else None,
            }
        )
    # a site right after the origin, consumed through its +n copy
    sites.append(_CountingSite("c1", "+", 6, ("G_ori",), "TAG", "terminal"))
    rows.append({"search": "target", "contig": "c1", "strand": "+", "start": N - 9,
                 "end": N + 21, "blocks": None})  # fmt: skip
    return sites, pd.DataFrame(rows)


def _oracle(rows: pd.DataFrame, sites: list[StopSite]) -> list[tuple[StopSite, ...]]:
    """Rule 15 by brute force (vectorised over the sites; reads bypass the counter)."""
    order = sorted(sites, key=lambda s: object.__getattribute__(s, "codon_start"))
    cs = np.array([object.__getattribute__(s, "codon_start") for s in order])
    out: list[tuple[StopSite, ...]] = []
    for r in rows.to_dict("records"):
        if r["search"] != "target" or isinstance(r["blocks"], str):
            out.append(())
            continue
        a, b = int(r["start"]), int(r["end"])
        hit = np.zeros(len(cs), dtype=bool)
        for x in (cs, cs + N):
            hit |= (a <= x) & (x + 3 <= b) & ((x - a) % 3 == 0)
        out.append(tuple(order[k] for k in np.flatnonzero(hit)))
    return out


def test_consumed_sites_does_not_scan_every_site_per_row() -> None:
    from giap.mapping.readthrough import consumed_sites

    sites, rows = _world()
    want = _oracle(rows, sites)
    assert sum(map(len, want)) > 50  # the oracle is not vacuous
    assert want[-1] and want[-1][0].locus_tags == ("G_ori",)
    _READS[0] = 0
    got = consumed_sites(rows, sites, {"c1": N}, {"c1": True})
    reads = _READS[0]
    assert got == want
    n_rows, n_sites = len(rows), len(sites)
    assert reads <= 8 * n_sites + 100 * n_rows, (
        f"{reads} codon_start reads for {n_rows} rows x {n_sites} sites"
    )
