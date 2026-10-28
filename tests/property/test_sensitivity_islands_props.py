"""Properties of the M5-10 ``islands`` sensitivity family
(docs/specs/M5-10-sensitivity-islands.md "Public interface"; GIPPA_SPEC.md §8.6).

- I1 for any island rows M5-06 reports, the flagged set is exactly
  ``{<direction>:<locus_tag>}`` over every gene of every row (directions kept apart);
- I2 the set does not depend on the row order of the gene table (M5-06 orders genes itself);
- I3 the ``islands`` row of ``sensitivity_rows`` satisfies n_shared <= min(n_primary,
  n_variant), jaccard = n_shared / |union| (NA iff both empty), for family sets of I1.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import sensitivity_islands_data as sid

_TAGS = st.lists(st.sampled_from([f"G{i}" for i in range(15)]), min_size=1, max_size=6,
                 unique=True)  # fmt: skip
_ROWS = st.lists(st.tuples(st.sampled_from(["high", "low"]), _TAGS), max_size=5)
_FN = [HealthCheck.function_scoped_fixture]


def _rng() -> np.random.Generator:
    from gippa.rng import for_step

    return for_step(sid.STREAM, 1)


@settings(max_examples=60, deadline=None, suppress_health_check=_FN)
@given(_ROWS)
def test_i1_ids_of_any_island_rows(
    monkeypatch: pytest.MonkeyPatch, rows: list[tuple[str, list[str]]]
) -> None:
    frame = sid.fake_islands(rows)
    with monkeypatch.context() as m:
        sid.patch_islands_frame(m, [frame])
        got = sid.call_islands_family(sid.tiny_tables(), sid.unit_config(), _rng())
    assert got == frozenset(f"{d}:{g}" for d, gs in rows for g in gs)
    assert all(x.count(":") == 1 for x in got)


@settings(max_examples=8, deadline=None)
@given(st.integers(0, 2**31 - 1))
def test_i2_gene_row_order_does_not_matter(perm_seed: int) -> None:
    t = sid.unit_tables()
    cfg = sid.unit_config()
    want = sid.call_islands_family(t, cfg, _rng())
    shuffled = t.genes.sample(frac=1.0, random_state=perm_seed).reset_index(drop=True)
    got = sid.call_islands_family(sid.IslandTables(shuffled, t.tus, t.replicons), cfg, _rng())
    assert got == want


@settings(max_examples=60, deadline=None, suppress_health_check=_FN)
@given(_ROWS, _ROWS)
def test_i3_islands_row(
    monkeypatch: pytest.MonkeyPatch,
    a: list[tuple[str, list[str]]],
    b: list[tuple[str, list[str]]],
) -> None:
    from gippa.context import sensitivity as S

    with monkeypatch.context() as m:
        sid.patch_islands_frame(m, [sid.fake_islands(a), sid.fake_islands(b)])
        t = sid.tiny_tables()
        p = sid.call_islands_family(t, sid.unit_config(), _rng())
        v = sid.call_islands_family(t, sid.unit_config(), _rng())
    r: Any = S.sensitivity_rows({"islands": p}, {("islands", "missing=min"): v}).iloc[0]
    assert r["family"] == "islands"
    assert int(r["n_primary"]) == len(p) and int(r["n_variant"]) == len(v)
    assert int(r["n_shared"]) == len(p & v) <= min(len(p), len(v))
    j = r["jaccard"]
    if not p and not v:
        assert j is None or j is pd.NA or (isinstance(j, float) and math.isnan(j))
    else:
        assert float(j) == pytest.approx(len(p & v) / len(p | v))
