"""Genetic-code arms of the chance-hit suite (M11-16f; docs/specs/M11-16-translation-tables-4-25.md
item 9: "longer stop-free frames -> more chance hits; the chance model must stay calibrated"; D14).

Arm 11 is ``test_chance_hits.py`` (B.2 composition, s = 0.069), unchanged. At genetic code ``t``
the residue frequencies the model uses come from a six-frame index built at ``t``
(``residue_freqs``: TGA adds W at 4, G at 25, and the stop share drops). Two checks:

1. **Model vs Monte Carlo at the code-``t`` composition** (the B.2 claim of the arm-11 suite with
   ``f_t`` for the E. coli-like vector): a random target and random peptides drawn from ``f_t``,
   hits with <= 1 substitution counted exactly by ``test_chance_hits.hits_le1``, scaled to
   N = 10^7; |MC - E(L, 1)| <= 20 % of E(L, 1) and the MC SE <= 8 % (fast, L = 8) / 5 % (full,
   L = 8-10). Same streams and seed as arm 11 (``calibration.chance_hits.L{L}``,
   ``seeds.calibration_base``). The B.2 column itself is for s = 0.069 and is not compared.
2. **Six-frame target, one-sided** (statistician to review): random peptides from ``f_t`` against
   the real six-frame translation at ``t`` of a gippasim genome, windows holding a stop or ``X``
   excluded. E_i is a ranking score with N = all six-frame positions, so the model may
   over-count (stops cut windows) but must not under-count: MC <= 1.2 * E(L, 1). This is where
   the longer stop-free frames show (measured on 16f, 200 kb, L = 8: MC / E 0.73 at 11, 0.82 at
   4, 0.84 at 25). Arm 11 runs too, as the reference.

Genome: gippasim ``SimConfig(genome_length=100_000 | 1_000_000, genetic_code=t)`` at
``calibration_base``.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pytest
from test_chance_hits import hits_le1

from gippa.mapping.chance import ResidueFreqs, expected_hits_uniform, residue_freqs
from gippa.mapping.sixframe import STOP
from gippa.rng import for_step

pytestmark = pytest.mark.calsuite("placement")

FAST_TABLES = (4,)
FULL_TABLES = (4, 25)
N_PEPTIDES = {8: 200_000, 9: 1_000_000, 10: 4_000_000}


def _base(tolerances: dict[str, object]) -> int:
    return int(tolerances["seeds"]["calibration_base"])  # type: ignore[index]


@cache
def _index(table: int, length: int, seed: int):  # type: ignore[no-untyped-def]
    import tempfile

    from gippa.io.fasta import read_genome_fasta
    from gippa.mapping.sixframe import build_sixframe
    from gippa.model import Replicon
    from gippasim.simulate import SimConfig, simulate

    with tempfile.TemporaryDirectory() as d:
        res = simulate(SimConfig(genome_length=length, genetic_code=table), seed, Path(d))
        (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    return build_sixframe({rec.id: rec.sequence}, {rec.id: rep}, table=table)


def _iid(L: int, n_peptides: int, n_target: int, seed: int, f: np.ndarray) -> tuple[float, float]:
    """``test_chance_hits.simulate`` with frequency vector ``f`` (same stream and seed)."""
    rng = for_step(f"calibration.chance_hits.L{L}", seed)
    target = rng.choice(19, size=n_target, p=f).astype(np.uint8)
    peptides = rng.choice(19, size=(n_peptides, L), p=f).astype(np.uint8)
    h = hits_le1(target, peptides)
    scale = 1e7 / (n_target - L + 1)
    return float(h.mean() * scale), float(h.std(ddof=1) / np.sqrt(n_peptides) * scale)


def _model_vs_mc(table: int, L: int, n_pep: int, n_target: int, se_max: float, seed: int) -> None:
    fr = residue_freqs(_index(table, 100_000, seed))
    mean, se = _iid(L, n_pep, n_target, seed, fr.f)
    model = expected_hits_uniform(L, 1, fr.s, 1e7)
    assert se <= se_max * mean, f"code {table} L={L}: MC SE {se:.3g} > {se_max:.0%} of {mean:.3g}"
    assert abs(mean - model) <= 0.2 * model, f"code {table} L={L}: MC {mean:.4g} vs E {model:.4g}"


@pytest.mark.calibration
@pytest.mark.parametrize("table", FAST_TABLES)
def test_fast_model_matches_mc_at_code_composition(
    table: int, tolerances: dict[str, object]
) -> None:
    _model_vs_mc(table, 8, 100_000, 200_000, 0.08, _base(tolerances))


@pytest.mark.calibration_full
@pytest.mark.parametrize("table", FULL_TABLES)
@pytest.mark.parametrize("L", [8, 9, 10])
def test_model_matches_mc_at_code_composition_full(
    table: int, L: int, tolerances: dict[str, object]
) -> None:
    _model_vs_mc(table, L, N_PEPTIDES[L], 10**6, 0.05, _base(tolerances))


def _hits_valid(target: np.ndarray, peptides: np.ndarray) -> np.ndarray:
    """Windows within <= 1 substitution of each peptide, over windows free of codes >= STOP."""
    L = peptides.shape[1]
    n = len(target) - L + 1
    bad = np.concatenate([[0], np.cumsum(target >= STOP)])
    valid = (bad[L : L + n] - bad[:n]) == 0
    base = 21
    w = base ** np.arange(L, dtype=np.int64)
    tk = np.zeros(n, np.int64)
    for i in range(L):
        tk += target[i : i + n].astype(np.int64) * w[i]
    tk = tk[valid]
    pk = (peptides.astype(np.int64) * w).sum(axis=1)

    def count(keys: np.ndarray, q: np.ndarray) -> np.ndarray:
        return np.searchsorted(keys, q, "right") - np.searchsorted(keys, q, "left")

    total = -(L - 1) * count(np.sort(tk), pk)
    for j in range(L):
        masked = np.sort(tk - (tk // w[j] % base) * w[j])
        total += count(masked, pk - peptides[:, j].astype(np.int64) * w[j])
    return total


def _sixframe_target(index, L: int) -> np.ndarray:  # type: ignore[no-untyped-def]
    """Core residues of every frame (``chance._core_residues``), frames joined by L stops."""
    from gippa.mapping.chance import _core_residues

    sep = np.full(L, STOP, dtype=np.uint8)
    parts = []
    for fr in index.frames:
        parts += [_core_residues(fr, index.replicons[fr.contig].length), sep]
    return np.concatenate(parts)


def test_valid_window_kernel_matches_brute_force() -> None:
    rng = np.random.default_rng(0)
    target = rng.integers(0, 4, size=500).astype(np.uint8)
    target[rng.integers(0, 500, 40)] = STOP
    peptides = rng.integers(0, 3, size=(60, 5)).astype(np.uint8)
    win = np.lib.stride_tricks.sliding_window_view(target, 5)
    ok = (win < STOP).all(axis=1)
    brute = np.array([(((win != p).sum(axis=1) <= 1) & ok).sum() for p in peptides])
    np.testing.assert_array_equal(_hits_valid(target, peptides), brute)


def _sixframe_check(table: int, length: int, n_pep: int, seed: int) -> dict[str, float]:
    L = 8
    idx = _index(table, length, seed)
    fr: ResidueFreqs = residue_freqs(idx)
    target = _sixframe_target(idx, L)
    rng = for_step(f"calibration.chance_hits.sixframe.L{L}", seed)
    peptides = rng.choice(19, size=(n_pep, L), p=fr.f).astype(np.uint8)
    h = _hits_valid(target, peptides)
    mean, se = float(h.mean()), float(h.std(ddof=1) / np.sqrt(n_pep))
    model = expected_hits_uniform(L, 1, fr.s, fr.n_positions)
    stop_share = float((target >= STOP).mean())  # separators included (diagnostic only)
    out = {"mc": mean, "se": se, "model": model, "ratio": mean / model, "stop_share": stop_share}
    assert se <= 0.08 * mean, (table, out)
    assert mean <= 1.2 * model, f"code {table}: chance model under-counts six-frame hits {out}"
    return out


@pytest.mark.calibration
@pytest.mark.parametrize("table", (11, *FAST_TABLES))
def test_fast_sixframe_chance_hits_not_undercounted(
    table: int, tolerances: dict[str, object]
) -> None:
    _sixframe_check(table, 100_000, 100_000, _base(tolerances))


@pytest.mark.calibration_full
@pytest.mark.parametrize("table", (11, *FULL_TABLES))
def test_sixframe_chance_hits_not_undercounted_full(
    table: int, tolerances: dict[str, object]
) -> None:
    _sixframe_check(table, 1_000_000, 200_000, _base(tolerances))


@pytest.mark.calibration
@pytest.mark.parametrize("table", FAST_TABLES)
def test_code_index_has_fewer_stops_than_code_11(table: int, tolerances: dict[str, object]) -> None:
    """The premise of the arms: the code-``t`` six-frame index is cut by fewer stops."""

    def share(t: int) -> float:
        aa = np.concatenate([f.aa for f in _index(t, 100_000, _base(tolerances)).frames])
        return float((aa == STOP).mean())

    assert share(table) < 0.8 * share(11), (share(table), share(11))
