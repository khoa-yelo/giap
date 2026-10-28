# M8-11: fgsea port II — the multilevel splitting stage (`sampleSize = 101`, `eps = 0`)

Status: superseded by ADR-0038 (GSEA runs R `fgsea::fgseaMultilevel` instead of a Python port; PR #193
closed)

Milestone M8 · track S (statistics engine) · **pass 2** · depends on: M8-02 · estimate ≈ 330 code
lines (+ ≈ 380 test lines) · **statistical** (label `stats`) · numba

## Spec sections implemented

- §7.7.4 "fgsea multilevel (Korotkevich 2021), ported to Python and numba from the fgsea source
  (MIT), with `sampleSize = 101`, `eps = 0`, `nPermSimple = 1000`. It estimates small p-values
  without a permutation floor."
- §8.2 "Where many tests share one correction, use analytical or multilevel p-values (Wallenius,
  fgsea multilevel) … which avoid the 1/(B + 1) floor."
- §8.3 `fgsea_multilevel(...)`; §6.6 determinism; §7.7.5 golden "|log₁₀ p − log₁₀ p_R| within
  3 × the combined `log2err` / log₂ 10" (M8-13).
- §17 "fgsea port misses the golden tolerances | golden tests | escalate; ADR to run R fgsea as an
  external tool".

## What is ported (licence as M8-02: MIT, attribution in every file header)

From **fgsea 1.36.2** (the same commit as M8-02): the multilevel half of `R/fgseaMultilevel.R`
(`multilevelImpl`, `multilevelError`, the `denomProb` correction, the `eps` handling, the
`isCpGeHalf` rule) and `src/fgseaMultilevel.cpp` + `src/fgseaMultilevelSupplement.{h,cpp}`
(`EsRuler`: the sample of `sampleSize` random sets of one size, levels at the sample median ES,
the MCMC perturbation that keeps samples above the current level, `extend` until the target ES or
`eps`, `getPvalue`, `calcLogCorrection`, `betaMeanLog`). The current GitHub master (1.39.x) adds a
`moveScale` argument and returns `log2err` from C++; **1.36.2's signatures are the ones to port**,
and the implementer records the diff in `docs/methods/fgsea-port.md`.

## Public interface

`src/gippa/stats/fgsea.py` (extends M8-02; kernels in `_fgsea_kernels.py`)
```python
@dataclass(frozen=True, slots=True)
class MultilevelResult:
    p: F64; log2err: F64; is_cp_ge_half: npt.NDArray[np.bool_]

def fgsea_multilevel(prep: Prepared, simple: SimpleResult, sample_size: int, eps: float,
                     rng: np.random.Generator) -> SimpleResult: ...      # simple, upgraded
```

1. **Which sets.** Exactly `simple.needs_multilevel` (M8-02 step 5), excluding the NaN
   (`modeFraction < 10`) sets. Grouped by set size; each size group runs one `EsRuler` per sign
   present (positive ES on the ranks, negative ES on the reversed ranks), as fgsea does.
2. **p.** `p = min(1, cppMPval / denomProb)` with `denomProb = (modeFraction + 1)/(nPerm + 1)`;
   per group `eps_group = eps · min(denomProb)`; with `eps = 0` there is no floor (the spec's
   setting). Sets with `p < eps` (only when eps > 0) get `p = eps`, `log2err` NaN.
3. **`log2err`** = `multilevelError(p, sampleSize)` =
   `sqrt(floor(−log2 p + 1) · (trigamma((sampleSize + 1)/2) − trigamma(sampleSize + 1))) / ln 2`;
   NaN where `isCpGeHalf` is false (fgsea's "P-values were likely overestimated" case).
4. **Unchanged from M8-02**: ES, NES, leading edge, and every set not needing multilevel (their
   simple p and simple `log2err`).
5. **Randomness.** One child seed per size group from `SeedSequence.spawn` in the group order
   (fgsea shuffles the group order with `sample()`; GIPPA iterates groups by size ascending, which
   only changes which random stream a group gets). Parallel over groups; results identical for any
   thread count. `sampleSize` even → +1 (fgsea's rule; GIPPA's 101 is odd).
6. **Call site.** M8-05's `gsea_table` calls `fgsea_multilevel` after `fgsea_simple` and sets
   `method = fgsea_multilevel` on every row of the run (one method per run, decision Q2). This
   switch is part of this PR (≤ 20 lines in `gippa.function.gsea`).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no set needs multilevel | result equals the simple result; `method` still `fgsea_multilevel` |
| `isCpGeHalf` false | p kept, `log2err` NaN, row flag `log2err_unestimated` (M8-05) |
| ES so extreme that `extend` cannot reach it before eps = 0 underflow | p = the smallest reachable value (≥ 0, never 0 by construction of the levels); `log2err` grows with depth — tested at ES of a top-30-of-3,000 set |
| one size group with both signs | two rulers, as fgsea |
| `sample_size < 3` | `ValueError` (config fixes 101) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | the top 50 of 2,000 N(0, 1) stats (seed 1) as a set; `nPerm` 1000 | simple p = 1/1001 → multilevel p < 10⁻¹⁰; `log2err` = `multilevelError(p, 101)` (0.64 at p = 10⁻⁶, 1.18 at p = 10⁻²⁰) |
| 2 | the same set, two seeds | log₂ p differs by < 3 · √2 · `log2err` (the M8-13 criterion applied to two Python runs) |
| 3 | a set with simple p = 0.3 | not multilevel; unchanged |
| 4 | a set of the bottom 50 | negative ES, the reversed-ranks ruler; p symmetric with example 1 within the tolerance |
| 5 | 8 threads vs 1 | identical p and `log2err` |
| 6 | eps = 10⁻¹⁰ (test only; config is 0) | a p below it is reported as 10⁻¹⁰ with `log2err` NaN |

## Acceptance criteria

- [ ] `tests/unit/stats/test_fgsea_multilevel.py` covers examples 1–6 and every degenerate row.
- [ ] Monte-Carlo check (fast calibration subset, `calibration` marker): over 50 seeds, for a set
  with a planted shift whose exact p is known by brute-force permutation at 10⁶ draws (p ≈ 10⁻⁴),
  |log₂ p̂ − log₂ p| ≤ 2 · `log2err` in ≥ 90 % of seeds.
- [ ] Performance: 2,500 TUs, 1,500 sets of which 100 need multilevel, in < 60 s on 4 CPUs.
- [ ] `docs/methods/fgsea-port.md` extended with the multilevel mapping and the 1.36.2 vs master
  differences.

## Calibration or attack work owed

- `gsea` suite at level (M8-15): multilevel p must not be anti-conservative in the tail; the
  statistician inspects the α = 0.01 rejection rate specifically.

## Hardening owed

- If M8-13 fails the p tolerance after 3 iterations (R9): stop, write the §7.7.4 fallback ADR (R
  fgsea as an external tool in `gippa-golden-r`-like runtime env) for Khoa, and leave
  `method = fgsea_simple` in place meanwhile — never widen `fgsea_log10p_k`.

## Questions for the lead

1. Group order: fgsea shuffles size groups; GIPPA iterates by size (proposed) — results are
   statistically equivalent and deterministic. Acceptable?
2. One `method` value per run (`fgsea_multilevel`, proposed), or per row (`fgsea_simple` for sets
   that did not need the multilevel stage)? fgsea itself reports no per-row method.
