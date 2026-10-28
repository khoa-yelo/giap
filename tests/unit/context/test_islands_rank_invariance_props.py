"""M5-07 invariant behind the island FWER calibration (GIPPA_SPEC.md §7.3.5(d) step 1; spec
docs/specs/M5-07-island-calibration-power.md): the island scan reads abundance only through Blom
scores of ranks, so the calibration on gippasim's log2 scale covers every strictly increasing
rescaling of abundance. Under one seed, a strictly increasing map of the detected abundances
leaves the TU scores, block lengths, null maximum and every window p unchanged.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as idat


def _scan(contigs: list[idat.Contig], seed: int) -> object:
    cfg = idat.config(bootstrap=49, permutations=49, alpha=0.05)
    genes, tus, replicons = idat.frames(contigs)
    _, scan = idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg,
                        np.random.default_rng(seed))  # fmt: skip
    return scan


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    seed=st.integers(0, 2**31 - 1),
    rho=st.sampled_from([0.0, 0.5, 0.9]),
    shift=st.floats(-10.0, 10.0),
    scale=st.floats(0.1, 10.0),
    cubic=st.booleans(),
)
def test_strictly_increasing_abundance_map_changes_nothing(
    seed: int, rho: float, shift: float, scale: float, cubic: bool
) -> None:
    rng = np.random.default_rng(seed)
    c = idat.sim_contig(rng, "chr1", 200, rho=rho, tu_p=0.6)

    def f(a: float) -> float:
        x = (a - 20.0) ** 3 if cubic else a
        return shift + scale * x

    mapped = idat.Contig(c.name, d=c.d, p=c.p, a=[f(a) for a in c.a], circular=c.circular,
                         tu=c.tu)  # fmt: skip
    a, b = _scan([c], seed), _scan([mapped], seed)
    sa, sb = a.series[0], b.series[0]  # type: ignore[attr-defined]
    assert sa.tu_ids == sb.tu_ids
    np.testing.assert_allclose(sa.score, sb.score, rtol=0, atol=1e-12)
    assert a.block_len == b.block_len  # type: ignore[attr-defined]
    assert a.l_hat == b.l_hat  # type: ignore[attr-defined]
    np.testing.assert_allclose(a.null_max, b.null_max, rtol=0, atol=1e-9)  # type: ignore[attr-defined]
    for pa, pb in zip(a.p_value, b.p_value, strict=True):  # type: ignore[attr-defined]
        np.testing.assert_array_equal(pa, pb)
