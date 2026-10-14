# M3-09: Chance-hit model: E(L, m), Poisson-binomial ranking score

Milestone M3 · depends on: M3-02 (residue codes, frame lengths), M4-02 (`stats/poisson_binomial.py`;
whichever merges first creates it with the signatures below) · creates `src/giap/mapping/chance.py` · pure Python (+ numba for the Monte
Carlo helper in tests) · estimate ≈ 160 code lines · **statistical**: section "Chance hits and the
ranking score" of `docs/methods/placement-fdr.md`, signed off by the statistician (§14.1)

## Spec sections implemented

- §7.2.5 "Formula: E(L, m) = N · Σ_{j≤m} C(L, j) (1 − s)^j s^(L−j), with N = number of six-frame
  positions." Per-position match probability "s = Σ f_a² ≈ 0.069".
- Procedure step 2: "**Ranking score.** E_i = expected chance hits for peptide i at cost m: N
  times the Poisson-binomial probability of ≤ m mismatches, with per-position match probability
  f_a, the genome's six-frame frequency of the peptide's residue."
- §8.3 `poisson_binomial_sf(p, k)` — "exact dynamic programming".
- Appendix B.2 table and "A Monte Carlo simulation on random sequences with these frequencies MUST
  reproduce the ≤ 1 column for L = 8–10 within 20 %." §7.2.9 calibration: "the chance-hit table is
  reproduced by simulation within 20 %".

## Public interface

```python
# stats/poisson_binomial.py (pure, typed; owned jointly with M4-02, same names and contracts)
def poisson_binomial_pmf(p: ArrayLike) -> np.ndarray: ...             # length n + 1, exact DP
def poisson_binomial_cdf(p: ArrayLike, k: int) -> float: ...          # P(X ≤ k), summed directly (no 1 − sf)
def poisson_binomial_sf(p: ArrayLike, k: int) -> float: ...           # P(X ≥ k)  (§8.3 contract)

# mapping/chance.py
@dataclass(frozen=True)
class ResidueFreqs:
    f: np.ndarray            # length 19 over residue codes (I → L), sums to 1; excludes * and X
    n_positions: int         # N = Σ len(frame.aa) over the 6·contigs target frames, circular extension excluded
    @property
    def s(self) -> float: ...   # Σ f_a²
def residue_freqs(index: SixFrameIndex) -> ResidueFreqs: ...
def expected_hits_uniform(L: int, m: int, s: float, n_positions: float) -> float: ...   # E(L, m)
def expected_hits(query: str, m: int, freqs: ResidueFreqs) -> float: ...               # E_i
```

`expected_hits` uses per-position mismatch probability `1 − f_{a_i}`; a `NOMATCH` residue has
`f = 0`. It is a ranking score (monotone in the evidence, exchangeable between target and decoy);
calibration comes from the decoys (M3-10), which the methods note states.

## Examples

| # | Input | Expected (rel. tol 1e-9 unless rounded) |
|---|---|---|
| 1 | `expected_hits_uniform(L, m, 0.069, 1e7)` for L ∈ {7, 8, 9, 10, 12, 14}, m ∈ {0,1,2} | the B.2 table: e.g. (7,0) `0.0744635`, (8,1) `0.559742`, (9,2) `2.36692`, (10,1) `0.00332505`, (12,2) `0.00141835`, (14,1) `1.05295e-7` |
| 2 | same for L = 6 | `1.1`, `88`, `3.0e3` (rounded, B.2) |
| 3 | `poisson_binomial_pmf([0.5, 0.5])` | `[0.25, 0.5, 0.25]` |
| 4 | `poisson_binomial_sf([0.1, 0.2, 0.3], 2)` | `0.098` (0.1·0.2·0.7 + 0.1·0.8·0.3 + 0.9·0.2·0.3 + 0.006) |
| 5 | 8-aa query with f = 0.1 at 4 positions and 0.05 at 4, N = 10⁷ | E = `0.00625` (m = 0), `0.70625` (m = 1), `34.38125` (m = 2) |
| 6 | constant f = 0.069 at every position, L = 8, m = 1 | equals `expected_hits_uniform(8, 1, 0.069, 1e7)` = `0.559742` |
| 7 | `residue_freqs` of a frame set `AAAL` + `ILLK` | f(A) = 3/8, f(L) = 4/8 (I merged), f(K) = 1/8; `n_positions = 8` |
| 8 | a query containing `U` | E at m = 0 is 0 (`f_U = 0`) |

## Edge cases

- `m > L` is a `ValueError`; `p` outside [0, 1] is a `ValueError`; empty `p` → pmf `[1.0]`.
- Very small values (E ~ 1e-12 at L = 40, m = 0) stay positive and finite (cdf summed directly
  from the pmf; no subtraction from 1).
- Circular extension positions are not counted in N.

## Acceptance criteria

- [ ] `tests/unit/test_chance.py`: examples 1–2, 5–8; examples 3–4 join M4-02's
  `tests/unit/stats/test_poisson_binomial.py` if this issue creates the module.
- [ ] Property: pmf sums to 1 (≤ 1e-12), equals the binomial pmf for constant p, and is
  invariant to permuting `p`.
- [ ] Calibration (`tests/calibration/test_chance_hits.py`, marker `calibration`): random target
  sequence with a 19-letter frequency vector of Σf² = 0.069 (N = 10⁶ positions; seed from
  `tolerances.yaml`), random peptides of L = 8, 9, 10 drawn from the same frequencies; the mean
  number of ≤ 1-substitution hits per peptide, scaled to N = 10⁷, is within 20 % of the B.2 ≤ 1
  column. Enough peptides that the Monte Carlo SE is ≤ 5 % (counting kernel in numba).
- [ ] Methods note section exists (formula, why E is a ranking score, B.2 reproduction) and is
  signed off.

## Questions for the lead

1. Should E count the extra explanations available to the tolerant matcher (swaps, blocks), not
   just substitutions? **Proposed default:** no, as the spec defines it (substitutions only); it is
   a ranking score and the decoy search, which uses all operations, calibrates the FDR.

## Lead decisions after statistician review (PR #28, approved)
1. On a circular replicon N counts the minus-strand codons that cross the origin: the searchable six-frame
   positions of a circular replicon of length n with no stops/X total exactly 2·n (not 2·n − 2). Keep a
   minus-strand codon when `first_nt - 2 < length`. Tests pin N exactly for a stop-free circular replicon.
2. Forwarded to M3-10: E is a ranking score only (decoys calibrate FDR); never report E as an absolute
   expectation without a stop-free-window correction; simple-sequence repeats support the
   `low_complexity` flag (see docs/methods/placement-fdr.md, chance-hit section).
