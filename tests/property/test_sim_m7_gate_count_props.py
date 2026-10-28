"""Properties of ``novel_orf_power(..., n_per_class)`` (M7-12 lead decision 26; PR #225;
GIPPA_SPEC.md §7.6.6, §8.4). Interface in tests/unit/test_sim_m7_gate_count.py.

For any seed, size and count:

* ``n_per_class=None`` gives the same config as omitting it;
* an int ``n`` requests ``n`` of each of the four classes and changes nothing else.

Config-only (no simulation), so many examples are cheap.
"""

from __future__ import annotations

import dataclasses
from typing import Any

from hypothesis import given, settings
from hypothesis import strategies as st

seeds = st.integers(min_value=0, max_value=2**31 - 1)
sizes = st.sampled_from(["100k", "300k", "1M"])
CLASSES4 = ("intergenic", "cds_outframe", "antisense", "nterm_extension")


def sc() -> Any:
    import gippasim.scenarios

    return gippasim.scenarios


@settings(max_examples=50, deadline=None)
@given(seed=seeds, size=sizes, n=st.integers(min_value=1, max_value=1000))
def test_n_per_class_only_moves_the_counts(seed: int, size: str, n: int) -> None:
    default = sc().novel_orf_power(seed, size)
    assert sc().novel_orf_power(seed, size, n_per_class=None) == default
    cfg = sc().novel_orf_power(seed, size, n_per_class=n)
    got = {**dict(cfg.n_novel_orfs), "nterm_extension": cfg.n_nterm_extensions}
    assert got == dict.fromkeys(CLASSES4, n)
    back = dataclasses.replace(
        cfg, n_novel_orfs=default.n_novel_orfs, n_nterm_extensions=default.n_nterm_extensions
    )
    assert back == default
