"""M3-08 review pin, large genome (docs/specs/M3-08-frameshift-splits.md, lead decision 1 after code
review): on a ~1 Mb genome the hashed seed prefilter triggers by itself (every query has far more
than 256 q-gram hits) and must not change the result compared with the unfiltered seeds
(``match_frameshift(..., prefilter_min_hits=NEVER)``, the keyword-only knob pinned in
tests/unit/test_frameshift_review.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from fs_genomes import SENSE_RESIDUES, background, codons_for, join_split, rc

from gippa.config.model import TolerantConfig
from gippa.mapping.decoys import DecoyType, build_decoys
from gippa.mapping.frameshift import match_frameshift
from gippa.mapping.sixframe import Frame, build_sixframe
from gippa.model import Replicon
from gippa.rng import for_step

NEVER = 10**12
SEGS = [("N", "GG"), ("GG", "N"), ("AG", "Q"), ("Q", "AG"), ("AS", "GT"), ("W", "AD"), ("GV", "R")]


def _genome(n: int, seed: int) -> tuple[str, list[str]]:
    """Background with 30 split peptides (12..40 aa; exact, + substitution, + block; δ = ±1;
    both strands; varied phases) on a circular replicon."""
    rng = np.random.default_rng(seed)
    g = background(n, seed)
    qs: list[str] = []
    for j in range(30):
        L = int(rng.integers(12, 41))
        p = int(rng.integers(4, L - 3))
        s = "".join(rng.choice(list(SENSE_RESIDUES), size=L))
        gp, gs = s[:p], s[p:]
        qp, qq = gp, gs
        r = j % 3
        if r == 1:
            i = int(rng.integers(len(qq)))
            qq = qq[:i] + ("W" if qq[i] != "W" else "Y") + qq[i + 1 :]
        elif r == 2:
            t, q = SEGS[j % len(SEGS)]
            if j % 6 < 3:
                gp, qp = gp + t, qp + q
            else:
                gs, qq = t + gs, q + qq
        if len(qp) + len(qq) > 40 or min(len(qp), len(qq)) < 4:
            continue
        nt = join_split(codons_for(gp, rng), codons_for(gs, rng), (1, -1)[j % 2])
        if nt is None:
            continue
        nt = nt if (j // 2) % 2 == 0 else rc(nt)
        pos = 5000 + (n // 31) * j + j % 8
        g = g[:pos] + nt + g[pos + len(nt) :]
        qs.append(qp + qq)
    return g[-30:] + g[:-30], qs


@pytest.mark.slow
@pytest.mark.parametrize("seed", [7])
def test_prefilter_on_a_one_megabase_genome_equals_unfiltered_seeds(seed: int) -> None:
    n = 1_000_000
    g, qs = _genome(n, seed)
    idx = build_sixframe({"c": g}, {"c": Replicon("c", n, "circular")})
    dec: dict[DecoyType, tuple[Frame, ...]] = {
        "reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))
    }
    decoy_of: dict[str, DecoyType] = dict.fromkeys(qs, "reverse")
    cfg = TolerantConfig()
    filtered = match_frameshift(qs, idx, dec, decoy_of, set(qs), cfg)
    unfiltered = match_frameshift(qs, idx, dec, decoy_of, set(qs), cfg, prefilter_min_hits=NEVER)
    pd.testing.assert_frame_equal(filtered, unfiltered)
    found = set(filtered.loc[filtered["search"] == "target", "peptide"])
    assert found == set(qs)
