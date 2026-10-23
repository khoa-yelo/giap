"""Attack nulls on the picked ranking key of the M7-11 novel-ORF class FDR (decision 11;
GIAP_SPEC.md §7.6.4, §8.1, §8.5; ADR-0009, ADR-0047; M4-24 shared-key rule; statistician,
PR #237 round 2).

Round 2 ranks a class on (q, then score) when every q is finite, on score alone when every score
is finite, else ``decoys_unscored``. These attacks probe the two branches:

- ``tied_q_unscored_decoys`` (block item B3): every row has a finite q, the decoys carry no score.
  Input q comes in plateaus (search engines report q = 0 for their whole top block and q is
  monotone, so ties are the rule). A NaN score is mapped to the worst secondary key, so inside a
  q plateau every scored target beats every unscored decoy: the tie no longer goes to the decoy and
  D collapses at the top. All-null world; claim: at most one ORF reaches q <= 0.05. M4-24 drops a
  key column that is incomplete on the paired decoys (``shared_key_columns``), which ranks on q
  alone here with ties to the decoy. ``xfail(strict=True)`` until the implementer applies the same
  rule to the secondary key.
- ``score_only_fallback``: decoys carry a score but no q, so the class ranks on score alone (the
  remedy chosen for B1). Mixed true/false world (intergenic, 40 % true); FDP within level, criterion
  ``fdr`` with ADR-0015 confirmation, as the other picked attacks.
- ``tied_q_scored``: the control for B3. Same plateaus, both sides scored; at level.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from helpers.novel_fdr_fixture import PickedWorld, q_of

pytestmark = pytest.mark.calsuite("novel_orfs")

SUITE = "novel_orfs"
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
Replicate = tuple[dict[str, float], list[str], set[str]]


def _fdr():  # type: ignore[no-untyped-def]
    import giap.novel.fdr as m

    return m


def _best(rng: np.random.Generator, mu: float) -> float:
    return float(rng.normal(mu, 1.0, rng.integers(1, 4)).max())


def _q_plateau(score: float) -> float:
    """Input q with plateaus: 0 for the engine's top block, else rounded to 1e-3."""
    return 0.0 if score > 0.0 else round(0.01 / (1.0 + math.exp(score)), 3)


def _world(
    rng: np.random.Generator,
    n: int,
    p_true: float,
    decoy_key: str,
) -> tuple[PickedWorld, list[str], set[str]]:
    """``decoy_key``: ``q_only`` (decoy score NaN), ``score_only`` (decoy q NaN) or ``both``."""
    w = PickedWorld()
    ids: list[str] = []
    false_ids: set[str] = set()
    for _ in range(n):
        true = rng.random() < p_true
        entry = w.entry("intergenic", "k0", "+")
        s_t = _best(rng, 3.0 if true else 0.0)
        orf, _ = w.target("intergenic", "k0", 0, entry=entry, q=_q_plateau(s_t), score=s_t)
        ids.append(orf.orf_id)
        if not true:
            false_ids.add(orf.orf_id)
        s_d = _best(rng, 0.0)
        q_d = math.nan if decoy_key == "score_only" else _q_plateau(s_d)
        sc_d = math.nan if decoy_key == "q_only" else s_d
        w.decoy(entry[0], 0, q=q_d, score=sc_d)
    return w, ids, false_ids


def _make(p_true: float, decoy_key: str) -> Callable[[np.random.Generator], Replicate]:
    def make(rng: np.random.Generator) -> Replicate:
        w, ids, false_ids = _world(rng, 1500, p_true, decoy_key)
        return dict(_fdr().picked_class_q(w.orfs, w.peptides, w.id_map)), ids, false_ids

    return make


def _gate(name: str, make: Callable[[np.random.Generator], Replicate], full: bool):  # type: ignore[no-untyped-def]
    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        qs, ids, false_ids = make(rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = {o for o in ids if q_of(qs, o) <= a}
            out[j] = len(called & false_ids) / len(called) if called else 0.0
        return out

    stream = f"{SUITE}.attack.picked.{name}" + (".full" if full else "")
    return calibrate_levels(SUITE, "fdr", levels, fdp, C.replicates(SUITE, full=full),
                            stream=stream)  # fmt: skip


@pytest.mark.xfail(strict=True, reason="PR #237 round-2 block B3: NaN decoy score loses q ties")
@pytest.mark.calibration
def test_picked_unscored_decoys_do_not_lose_q_ties() -> None:
    rng = C.replicate_rng(f"{SUITE}.attack.picked.tied_q_unscored_decoys", 0)
    w, ids, _ = _world(rng, 400, 0.0, "q_only")
    try:
        res = _fdr().picked_class_fdr(w.orfs, w.peptides, w.id_map)
    except ValueError:
        return
    n_called = sum(q_of(res.q, o) <= 0.05 for o in ids)
    assert n_called <= 1, f"{n_called} of {len(ids)} null ORFs at q <= 0.05 ({res.classes})"


@pytest.mark.parametrize("full", MODES)
def test_picked_score_only_fallback_holds_level(full: bool) -> None:
    _gate("score_only_fallback", _make(0.4, "score_only"), full).assert_ok()


@pytest.mark.parametrize("full", MODES)
def test_picked_tied_q_with_scored_decoys_holds_level(full: bool) -> None:
    _gate("tied_q_scored", _make(0.4, "both"), full).assert_ok()
