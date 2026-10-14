"""Attack nulls for the M3-03 genome decoys (GIAP_SPEC.md §7.2.1, §8.2, §8.5; statistician).

Placement FDR (M3-10) estimates the chance target hits of wrong peptides by their decoy hits, so for
a wrong peptide that is unrelated to the genome, P(hit target) must equal P(hit decoy). These
attacks check that on the search spaces built by ``build_sixframe`` + ``build_decoys``
(docs/methods/placement-fdr.md, "Genome decoys"):

- ``composition null``: giapsim genomes (100 kb, one per replicate); wrong peptides drawn i.i.d.
  from the genome's six-frame residue frequencies (the null behind E_i, M3-09). The exact expected
  hit probability per peptide, sum over the *distinct* stop- and X-free L-mers w of prod f(w_j), is
  computed for targets and both decoys (no peptide sampling noise); reverse must equal target
  exactly (a theorem for any i.i.d. peptide law), shuffle must agree within replicate error. Sampled
  counts (i.i.d., tryptic with a C-terminal K/R, and <= 1 substitution at 8 aa) must agree within
  binomial error.
- ``low-complexity genome``: the same genomes with simple-sequence repeats (codon, dinucleotide and
  homopolymer runs) and ten copies of a 1.5-kb element. Repeats shrink the set of distinct target
  L-mers; the shuffled decoy destroys them, so it must not be *anti*-conservative (it is strongly
  conservative; the size of that power loss is recorded in the methods note).
- ``shuffle uniformity``: the window shuffle is uniform over permutations jointly across windows
  (around ``*`` and ``X``) and position by residue in a full 50-aa window.

Criteria: equality tests at p > 0.001 (``ks_min_p``); replicate means within ``se_multiplier`` SE
plus a 1 % equivalence margin. R = ``replicates_fast.placement`` (marker ``calibration``) or
``replicates_full.placement`` (``calibration_full``); seeds from ``seeds.calibration_base``.
"""

from __future__ import annotations

from collections.abc import Callable
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binomtest, chi2, chisquare

from giap.io.fasta import read_genome_fasta
from giap.mapping.decoys import build_decoys, shuffle_frame
from giap.mapping.sixframe import AMBIG, STOP, Frame, build_sixframe, encode_target
from giap.model import Replicon
from giap.rng import for_step

GENOME_NT = 100_000
EQUIV = 0.01  # equivalence margin on the decoy/target ratio of expected hits
KINDS = ("reverse", "shuffle")
K_R = (int(encode_target("K")[0]), int(encode_target("R")[0]))

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


# --- genomes ---------------------------------------------------------------------------------


def _giapsim_genome(seed: int, outdir: Path) -> str:
    import giapsim

    res = giapsim.simulate(giapsim.SimConfig(genome_length=GENOME_NT), seed=seed, outdir=outdir)
    (rec,) = read_genome_fasta(res.genome_fasta)
    return rec.sequence


_SSR_UNITS = (b"CAG", b"GCC", b"AGCGCA", b"AT", b"GA", b"A", b"CCGCCA", b"GGCTCA")


def _low_complexity(seq: str, rng: np.random.Generator) -> str:
    """24 simple-sequence repeats of 360 nt and 10 copies of a 1.5-kb element (≈ 24 % of 100 kb)."""
    s = bytearray(seq.encode())
    for i in range(24):
        unit = _SSR_UNITS[i % len(_SSR_UNITS)]
        p = int(rng.integers(0, len(s) - 400))
        s[p : p + 360] = (unit * 400)[:360]
    element = bytes(s[1000:2500])
    for _ in range(10):
        p = int(rng.integers(0, len(s) - 1600))
        s[p : p + 1500] = element
    return s.decode()


# --- search spaces ---------------------------------------------------------------------------


def _spaces(seq: str, seed: int) -> dict[str, tuple[Frame, ...]]:
    idx = build_sixframe({"c": seq}, {"c": Replicon("c", len(seq), "circular", "chromosome")})
    rng = for_step("sixframe_index", seed)
    out = {"target": idx.frames}
    for kind in KINDS:
        out[kind] = build_decoys(idx, kind, 50, rng)  # type: ignore[arg-type]
    return out


def _lmers(frames: tuple[Frame, ...], length: int) -> np.ndarray:
    """All stop- and X-free L-mers (rows of residue codes), with multiplicity."""
    rows = [
        w[(w < STOP).all(axis=1)]
        for f in frames
        if len(f.aa) >= length
        for w in [np.lib.stride_tricks.sliding_window_view(f.aa.astype(np.int64), length)]
    ]
    return np.concatenate(rows)


def _code(rows: np.ndarray) -> np.ndarray:
    c = np.zeros(len(rows), dtype=np.int64)
    for j in range(rows.shape[1]):
        c = c * 21 + rows[:, j]
    return c


def _freqs(frames: tuple[Frame, ...]) -> np.ndarray:
    aa = np.concatenate([f.aa for f in frames])
    return np.bincount(aa[aa < STOP], minlength=STOP) / int((aa < STOP).sum())


def _p_hit(rows: np.ndarray, f: np.ndarray) -> float:
    """P(an i.i.d.-f peptide equals some distinct L-mer) = sum over distinct w of prod f(w_j)."""
    u = np.unique(_code(rows))
    logp = np.zeros(len(u))
    for _ in range(rows.shape[1]):
        logp += np.log(f)[u % 21]
        u = u // 21
    return float(np.exp(logp).sum())


def _hits_exact(rows: np.ndarray, peps: np.ndarray) -> int:
    return int(np.isin(_code(peps), _code(rows)).sum())


def _hits_sub1(rows: np.ndarray, peps: np.ndarray) -> int:
    """Peptides with a hit at <= 1 substitution (mask each position in turn)."""
    hit = np.zeros(len(peps), dtype=bool)
    for j in range(peps.shape[1]):
        r, p = rows.copy(), peps.copy()
        r[:, j] = AMBIG
        p[:, j] = AMBIG
        hit |= np.isin(_code(p), _code(r))
    return int(hit.sum())


def _iid(rng: np.random.Generator, f: np.ndarray, n: int, length: int) -> np.ndarray:
    return rng.choice(STOP, size=(n, length), p=f)


def _tryptic(rng: np.random.Generator, f: np.ndarray, n: int, length: int) -> np.ndarray:
    peps = _iid(rng, f, n, length)
    peps[:, -1] = rng.choice(K_R, size=n)
    return peps


Sampler = Callable[[np.random.Generator, np.ndarray, int, int], np.ndarray]
#: name -> (peptide law, length, peptides per replicate, hit counter, run every k-th replicate)
SAMPLED: dict[str, tuple[Sampler, int, int, Callable[[np.ndarray, np.ndarray], int], int]] = {
    "iid_exact_7": (_iid, 7, 20_000, _hits_exact, 1),
    "tryptic_exact_7": (_tryptic, 7, 20_000, _hits_exact, 1),
    "iid_sub1_8": (_iid, 8, 10_000, _hits_sub1, 5),  # the masked search is the slow part
}


# --- attacks ---------------------------------------------------------------------------------


def _run(
    suite: str,
    tolerances: dict[str, object],
    tmp: Path,
    low_complexity: bool,
) -> tuple[dict[tuple[str, int], np.ndarray], dict[tuple[str, str], int], dict[str, int]]:
    reps, seeds = tolerances[suite], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    r = int(reps["placement"]) // (5 if low_complexity else 1)
    base = int(seeds["calibration_base"])
    ratio: dict[tuple[str, int], list[float]] = {}
    hits: dict[tuple[str, str], int] = {}
    hits_t: dict[str, int] = {}
    for i in range(r):
        rng = np.random.default_rng([base, i, int(low_complexity)])
        seq = _giapsim_genome(base + i, tmp / f"g{i}")
        if low_complexity:
            seq = _low_complexity(seq, rng)
        spaces = _spaces(seq, base + i)
        f = _freqs(spaces["target"])
        rows_by_len = {n: {k: _lmers(v, n) for k, v in spaces.items()} for n in (7, 8)}
        for length, rows in rows_by_len.items():
            assert all(len(rows[k]) == len(rows["target"]) for k in KINDS)  # same search space
            p_t = _p_hit(rows["target"], f)
            for kind in KINDS:
                ratio.setdefault((kind, length), []).append(_p_hit(rows[kind], f) / p_t)
        for name, (sampler, length, n, count, every) in SAMPLED.items():
            if i % every:
                continue
            peps = sampler(rng, f, n, length)
            rows = rows_by_len[length]
            hits_t[name] = hits_t.get(name, 0) + count(rows["target"], peps)
            for kind in KINDS:
                hits[(name, kind)] = hits.get((name, kind), 0) + count(rows[kind], peps)
    return {k: np.array(v) for k, v in ratio.items()}, hits, hits_t


@pytest.mark.parametrize("suite", SUITES)
def test_composition_null_decoy_hits_equal_target_hits(
    tolerances: dict[str, object], suite: str, tmp_path: Path
) -> None:
    k, min_p = float(tolerances["se_multiplier"]), float(tolerances["ks_min_p"])  # type: ignore[arg-type]
    ratio, hits, hits_t = _run(suite, tolerances, tmp_path, low_complexity=False)
    for (kind, length), x in ratio.items():
        if kind == "reverse":  # sum over rev(W) of an i.i.d. law equals the sum over W
            np.testing.assert_allclose(x, 1.0, rtol=1e-9, err_msg=f"reverse L={length}")
        else:
            se = x.std(ddof=1) / np.sqrt(len(x))
            assert abs(x.mean() - 1.0) <= k * se + EQUIV, (
                f"{kind} L={length}: mean D/T expected-hit ratio {x.mean():.4f} (SE {se:.4f})"
            )
    for (name, kind), d in hits.items():
        t = hits_t[name]
        p = binomtest(d, d + t, 0.5).pvalue
        assert p > min_p, f"{name} {kind}: target {t} vs decoy {d} hits, binomial p = {p:.2g}"


@pytest.mark.parametrize("suite", SUITES)
def test_low_complexity_genome_decoys_not_anticonservative(
    tolerances: dict[str, object], suite: str, tmp_path: Path
) -> None:
    k, min_p = float(tolerances["se_multiplier"]), float(tolerances["ks_min_p"])  # type: ignore[arg-type]
    ratio, hits, hits_t = _run(suite, tolerances, tmp_path, low_complexity=True)
    for (kind, length), x in ratio.items():
        if kind == "reverse":
            np.testing.assert_allclose(x, 1.0, rtol=1e-9, err_msg=f"reverse L={length}")
        else:
            se = x.std(ddof=1) / np.sqrt(len(x))
            assert x.mean() >= 1.0 - k * se - EQUIV, (
                f"{kind} L={length}: decoys under-count chance hits, ratio {x.mean():.4f}"
            )
    for (name, kind), d in hits.items():
        t = hits_t[name]
        p = binomtest(t, d + t, 0.5, alternative="greater").pvalue  # H1: target hit more often
        assert p > min_p, f"{name} {kind}: target {t} vs decoy {d} hits, one-sided p = {p:.2g}"


# --- exactness of the window shuffle ---------------------------------------------------------


def _codes(s: str) -> np.ndarray:
    return encode_target(s)


@pytest.mark.calibration
def test_shuffle_uniform_jointly_over_windows_around_stop_and_x() -> None:
    # window 4: "ACDE" | X | "FGHK" is one segment of 8 non-X residues -> windows ACDE and FGHK
    # (the X skipped, not a boundary); "*" ends it; "LMNP" is a new segment. 24^3 joint cells is
    # too many, so test (w1, w2) jointly (576 cells) and w3's marginal, plus w1-w3 independence.
    aa = _codes("ACDEXFGHK*LMNP")
    perms = {p: i for i, p in enumerate(permutations(range(4)))}
    n = 576 * 40
    rng = for_step("sixframe_index", 20261014)
    w1, w2, w3 = np.zeros(n, int), np.zeros(n, int), np.zeros(n, int)
    for i in range(n):
        out = shuffle_frame(aa, 4, rng)
        assert out[4] == AMBIG and out[9] == STOP
        for w, sl, ref in ((w1, slice(0, 4), aa[0:4]), (w2, slice(5, 9), aa[5:9]),
                           (w3, slice(10, 14), aa[10:14])):  # fmt: skip
            w[i] = perms[tuple(int(np.flatnonzero(ref == c)[0]) for c in out[sl])]
    joint = np.bincount(w1 * 24 + w2, minlength=576)
    assert chisquare(joint).pvalue > 0.001
    assert chisquare(np.bincount(w3, minlength=24)).pvalue > 0.001
    assert chisquare(np.bincount(w1 * 24 + w3, minlength=576)).pvalue > 0.001


@pytest.mark.calibration
def test_shuffle_full_window_position_by_residue_uniform() -> None:
    # 50 distinct labels are impossible with 19 residues; use residue identity of a 50-aa window
    # with known composition and test that every residue's count at every position matches
    # n * (its share of the window), pooled chi-square over 19 x 50 cells.
    rng0 = np.random.default_rng(5)
    aa = rng0.integers(0, STOP, size=50).astype(np.uint8)
    rng = for_step("sixframe_index", 20261014)
    n = 4000
    counts = np.zeros((STOP, 50))
    for _ in range(n):
        out = shuffle_frame(aa, 50, rng)
        counts[out, np.arange(50)] += 1
    share = np.bincount(aa, minlength=STOP) / 50
    keep = share > 0
    expected = n * np.repeat(share[keep, None], 50, axis=1)
    stat = ((counts[keep] - expected) ** 2 / expected).sum()
    dof = (int(keep.sum()) - 1) * 49

    assert chi2.sf(stat, dof) > 0.001, (stat, dof)
