"""Properties of the M4-19 sensitivity framework (docs/specs/M4-19-sensitivity-regions.md acceptance
"Property"; lead decisions 3 and 4).

- J1 ``jaccard`` is symmetric, in [0, 1], 1 iff the sets are equal and non-empty, None iff both
  are empty, and equals |A n B| / |A u B|;
- R1 ``sensitivity_rows``: one row per rerun, n_primary = |P|, n_variant = |V|,
  n_shared = |P n V| <= min, jaccard as J1 (NA when both empty);
- F1 the regions family's flagged set is exactly the IDs ``<region_id>:<contig>:<stat>`` of the
  M4-17 rows with ``null_model = circular_shift``, p not NA and q <= ``alpha.context`` (any frame).
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
import sensitivity_data as sd

IDS = st.frozensets(st.sampled_from([f"r{i}:c{i % 3}:o_minus_e" for i in range(12)]), max_size=12)


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


@settings(max_examples=300, deadline=None)
@given(IDS, IDS)
def test_j1_jaccard(a: frozenset[str], b: frozenset[str]) -> None:
    S = _S()
    j, k = S.jaccard(a, b), S.jaccard(b, a)
    if not a and not b:
        assert j is None and k is None
        return
    assert j == k
    assert 0.0 <= j <= 1.0
    assert (j == 1.0) == (a == b)
    assert j == pytest.approx(len(a & b) / len(a | b), abs=1e-15)


@settings(max_examples=150, deadline=None)
@given(IDS, st.lists(IDS, min_size=1, max_size=4))
def test_r1_rows(primary: frozenset[str], variants: list[frozenset[str]]) -> None:
    S = _S()
    reruns = {("regions", f"v{i}"): v for i, v in enumerate(variants)}
    df = S.sensitivity_rows({"regions": primary}, reruns)
    assert len(df) == len(reruns)
    for _, r in df.iterrows():
        v = reruns[(r["family"], r["variant"])]
        assert int(r["n_primary"]) == len(primary)
        assert int(r["n_variant"]) == len(v)
        assert int(r["n_shared"]) == len(primary & v) <= min(len(primary), len(v))
        j = r["jaccard"]
        if not primary and not v:
            assert j is None or j is pd.NA or (isinstance(j, float) and math.isnan(j))
        else:
            assert float(j) == pytest.approx(len(primary & v) / len(primary | v), abs=1e-15)


_ROW = st.tuples(
    st.sampled_from(["prophage_1", "prophage_2", "bgc_1", "amr_1", "crispr_1"]),
    st.sampled_from(["chr1", "p1"]),
    st.sampled_from(["o_minus_e", "hl_shift"]),
    st.one_of(st.just(math.nan), st.floats(0.0, 1.0)),
    st.one_of(st.just(math.nan), st.floats(0.0, 1.0)),
    st.sampled_from(["circular_shift", "circular_shift", "none"]),
)


@settings(
    max_examples=120, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(st.lists(_ROW, max_size=12, unique_by=lambda r: (r[0], r[1], r[2])),
       st.sampled_from([0.01, 0.05, 0.10]))  # fmt: skip
def test_f1_flagged_set_rule(
    monkeypatch: pytest.MonkeyPatch, rows: list[tuple[Any, ...]], alpha: float
) -> None:
    S = _S()
    typ = {"prophage": "prophage", "bgc": "bgc", "amr": "amr", "crispr": "crispr"}
    full = [(r, typ[r.split("_")[0]], c, s, p, q, n) for r, c, s, p, q, n in rows]
    frame = sd.fake_tests(full)
    with monkeypatch.context() as m:
        sd.patch_region_tests(m, frame)
        got = S.regions_family(sd.family_tables().genes,
                               sd.config(stats={"alpha": {"context": alpha}}),
                               np.random.default_rng(0))  # fmt: skip
    assert got == sd.flagged_from(frame, alpha)
    assert all(len(x.split(":")) == 3 for x in got)
