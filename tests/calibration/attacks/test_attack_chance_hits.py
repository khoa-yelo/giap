"""Attack nulls for the M3-09 chance-hit model (GIPPA_SPEC.md §7.2.5, Appendix B.2; statistician).

E(L, m) = N sum_{j<=m} C(L, j) (1 - s)^j s^(L - j) assumes that every one of the N six-frame
positions starts a searchable window of L residues drawn i.i.d. from the genome's composition f.
Real six-frame translations break this in three ways; these attacks measure the direction and size
of the error on genomes built through ``build_sixframe`` and ``residue_freqs``
(docs/methods/placement-fdr.md, "Chance hits and the ranking score"):

1. **Stops and X.** N counts every residue (spec: "N = number of six-frame positions"), but a
   window containing ``*`` or ``X`` can never match. Factor W_L / N < 1 (W_L = stop- and X-free
   windows), so E *over*-states the absolute count (conservative), more so for longer L.
2. **Heterogeneous composition.** Coding, antisense and intergenic frames have different residue
   laws g_k. For an i.i.d.-f peptide E[prod f(t_j)] = sum_k w_k (sum_a f_a g_k(a))^L, which by
   Jensen is >= s^L: E *under*-states the conditional count, more so for longer L.
3. **Low complexity.** Simple-sequence repeats put long runs of the same few residues into the
   target; windows of abundant residues have prod f(t_j) >> s^L, so E under-states again.

No peptides are sampled for the main quantities: the exact expected number of hit windows for a
peptide drawn i.i.d. from f is sum over searchable windows t of P(<= m mismatches), with mismatch
probability 1 - f(t_j) at position j (the same Poisson-binomial as E_i, by symmetry), computed in
closed form for m <= 2. A sampled exact-match count cross-checks that identity.

Reported ratios (methods note): ``stop`` = W_L / N, ``cond`` = observed / (W_L * tail),
``total`` = observed / E(L, m) = stop * cond. Assertions encode the *directions*, with margins:
i.i.d.-nucleotide genomes (any GC) must have cond within 3 % of 1 (the model is exact there once
stops are accounted for; replicate SE ~ 0.6 % at L = 8-10); gippasim and low-complexity genomes
must have cond > 1, with the low-complexity genome worse; total < 1 at L >= 8 on every genome
without low-complexity blocks and total > 1 at m = 0, L >= 8 with them. R = 8 genomes of 100 kb
per family (~10 s).
"""

from __future__ import annotations

from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest

from gippa.io.fasta import read_genome_fasta
from gippa.mapping.chance import expected_hits_uniform, residue_freqs
from gippa.mapping.sixframe import STOP, SixFrameIndex, build_sixframe
from gippa.model import Replicon
from gippa.rng import for_step

GENOME_NT = 100_000
LENGTHS = (7, 8, 10, 12)
MS = (0, 1, 2)
R_FAST = 8

_SSR_UNITS = (b"CAG", b"GCC", b"AGCGCA", b"AT", b"GA", b"A", b"CCGCCA", b"GGCTCA")


# --- genomes ---------------------------------------------------------------------------------


def _gippasim(seed: int, outdir: Path) -> str:
    import gippasim

    res = gippasim.simulate(gippasim.SimConfig(genome_length=GENOME_NT), seed=seed, outdir=outdir)
    (rec,) = read_genome_fasta(res.genome_fasta)
    return rec.sequence


def _iid_nt(rng: np.random.Generator, gc: float) -> str:
    p = np.array([1 - gc, gc, gc, 1 - gc]) / 2
    return bytes(np.frombuffer(b"ACGT", np.uint8)[rng.choice(4, size=GENOME_NT, p=p)]).decode()


def _low_complexity(seq: str, rng: np.random.Generator) -> str:
    """24 simple-sequence repeats of 360 nt (≈ 9 % of 100 kb), as in the decoy attack."""
    s = bytearray(seq.encode())
    for i in range(24):
        unit = _SSR_UNITS[i % len(_SSR_UNITS)]
        p = int(rng.integers(0, len(s) - 400))
        s[p : p + 360] = (unit * 400)[:360]
    return s.decode()


def _index(seq: str) -> SixFrameIndex:
    return build_sixframe({"c": seq}, {"c": Replicon("c", len(seq), "circular", "chromosome")})


# --- exact expected hits ---------------------------------------------------------------------


def _windows(index: SixFrameIndex, length: int) -> np.ndarray:
    """Searchable genomic windows, one per distinct genomic start (origin-crossing kept, the
    circular-extension duplicates dropped); rows of residue codes, ``*``/X windows included."""
    rows = []
    for fr in index.frames:
        n_rep = index.replicons[fr.contig].length
        w = np.lib.stride_tricks.sliding_window_view(fr.aa, length)
        k = np.arange(len(w))
        start = fr.offset + 3 * k if fr.strand == "+" else fr.ext_nt - fr.offset - 3 * (k + length)
        rows.append(w[start < n_rep])
    return np.concatenate(rows)


def _observed(rows: np.ndarray, f: np.ndarray) -> dict[int, float]:
    """sum over stop/X-free windows of P(i.i.d.-f peptide within <= m substitutions), m = 0..2."""
    ok = rows[(rows < STOP).all(axis=1)]
    ft = f[ok]
    prod = ft.prod(axis=1)
    r = (1.0 - ft) / ft
    e1 = r.sum(axis=1)
    e2 = (e1**2 - (r**2).sum(axis=1)) / 2
    return {
        0: float(prod.sum()),
        1: float((prod * (1 + e1)).sum()),
        2: float((prod * (1 + e1 + e2)).sum()),
    }


def _tail(length: int, m: int, s: float) -> float:
    return expected_hits_uniform(length, m, s, 1.0)


def _ratios(index: SixFrameIndex) -> dict[tuple[int, int], tuple[float, float, float]]:
    """(stop, cond, total) per (L, m)."""
    fr = residue_freqs(index)
    f = np.append(fr.f, [1.0] * 3)  # pad codes >= STOP (never used: those windows are dropped)
    out = {}
    for length in LENGTHS:
        rows = _windows(index, length)
        # one window per genomic codon start: 2 * length per circular replicon. M3-09 review: the
        # implementation drops the 2 origin-crossing minus-strand codons (N = 2 * length - 2);
        # accept either until fixed (relative effect 1e-5).
        assert 0 <= len(rows) - fr.n_positions <= 2
        w_l = int((rows < STOP).all(axis=1).sum())
        obs = _observed(rows, f)
        for m in MS:
            model = expected_hits_uniform(length, m, fr.s, fr.n_positions)
            out[(length, m)] = (
                w_l / fr.n_positions,
                obs[m] / (w_l * _tail(length, m, fr.s)),
                obs[m] / model,
            )
    return out


def _genomes(base: int, tmp: Path) -> dict[str, list[SixFrameIndex]]:
    out: dict[str, list[SixFrameIndex]] = {
        k: [] for k in ("iid_gc50", "iid_gc75", "gippasim", "gippasim_lowcx")
    }
    for i in range(R_FAST):
        rng = for_step(f"calibration.attack.chance_hits.{i}", base)
        out["iid_gc50"].append(_index(_iid_nt(rng, 0.50)))
        out["iid_gc75"].append(_index(_iid_nt(rng, 0.75)))
        seq = _gippasim(base + i, tmp / f"g{i}")
        out["gippasim"].append(_index(seq))
        out["gippasim_lowcx"].append(_index(_low_complexity(seq, rng)))
    return out


def _mean_ratios(indexes: list[SixFrameIndex]) -> dict[tuple[int, int], np.ndarray]:
    per = [_ratios(ix) for ix in indexes]
    return {k: np.mean([p[k] for p in per], axis=0) for k in per[0]}


# --- tests -----------------------------------------------------------------------------------


def _seed(tolerances: dict[str, object]) -> int:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    return int(seeds["calibration_base"])


def test_exact_observed_matches_sampled_exact_hits() -> None:
    """The closed form at m = 0 equals the mean sampled hit count (with multiplicity)."""
    rng = np.random.default_rng(7)
    ix = _index(_iid_nt(rng, 0.5)[:20_000])
    fr = residue_freqs(ix)
    rows = _windows(ix, 4)
    ok = rows[(rows < STOP).all(axis=1)].astype(np.int64)
    key = np.sort((ok * 19 ** np.arange(4)).sum(axis=1))
    peps = rng.choice(STOP, size=(200_000, 4), p=fr.f).astype(np.int64)
    pk = (peps * 19 ** np.arange(4)).sum(axis=1)
    hits = np.searchsorted(key, pk, "right") - np.searchsorted(key, pk, "left")
    exact = _observed(rows, np.append(fr.f, [1.0] * 3))[0]
    se = hits.std(ddof=1) / np.sqrt(len(hits))
    assert abs(hits.mean() - exact) <= 4 * se, (hits.mean(), exact, se)


@pytest.mark.calibration
def test_attack_chance_hit_model_directions(tolerances: dict[str, object], tmp_path: Path) -> None:
    res = {k: _mean_ratios(v) for k, v in _genomes(_seed(tolerances), tmp_path).items()}
    for name, table in res.items():  # printed for the methods note (pytest -s)
        for (length, m), (stop, cond, total) in sorted(table.items()):
            print(
                f"{name:14s} L={length:2d} m={m} stop={stop:.3f} cond={cond:.3f} total={total:.3f}"
            )
    for key in res["gippasim"]:
        length, m = key
        # 1. stops: W_L / N < 1 and decreasing in L (the conservative factor)
        for name in res:
            assert res[name][key][0] < 1.0
        # i.i.d. nucleotides (any composition): the model is exact given the stop factor
        if length <= 10:  # L = 12 conditional ratios at 100 kb carry ~ 5 % window-sampling noise
            for name in ("iid_gc50", "iid_gc75"):
                assert abs(res[name][key][1] - 1.0) <= 0.03, (name, key, res[name][key])
        # 2./3. heterogeneous frames and low complexity: E under-states the conditional count
        assert res["gippasim"][key][1] > 1.0, (key, res["gippasim"][key])
        assert res["gippasim_lowcx"][key][1] > res["gippasim"][key][1], key
        if length >= 8:  # net effect on realistic genomes without repeats: still conservative
            for name in ("iid_gc50", "iid_gc75", "gippasim"):
                assert res[name][key][2] < 1.0, (name, key)
            if m == 0:  # simple-sequence repeats: E is anti-conservative as an absolute count
                assert res["gippasim_lowcx"][key][2] > 1.0, key
    for name in res:
        stops = [res[name][(length, 0)][0] for length in LENGTHS]
        assert all(a > b for a, b in pairwise(stops)), (name, stops)
