# M8-02: fgsea port I — ES, leading edge, the simple permutation stage, NES and `log2err`

Milestone M8 · track S (statistics engine) · **pass 1** · depends on: nothing · estimate ≈ 320 code
lines (+ ≈ 400 test lines) · **statistical** (label `stats`) · numba

## Spec sections implemented

- §7.7.4 "Method (G15). fgsea multilevel (Korotkevich 2021), ported to Python and numba from the
  fgsea source (MIT), with `sampleSize = 101`, `eps = 0`, `nPermSimple = 1000`." This issue ports
  everything `fgseaMultilevel` does **before** the multilevel splitting: stat preparation
  (`prepareStats`, `preparePathways`), the observed ES and leading edge (`calcGseaStat`), the
  `nPermSimple` random-set stage (`fgseaSimpleImpl`), NES, the simple-stage p and `log2err`, the
  `modeFraction < 10` rule, and the decision which sets need the multilevel stage. M8-11 ports the
  splitting itself.
- §7.7.4 "Reported per set: ES, NES, p, `log2err` … and the leading edge".
- §2.4 G15; §8.3 row `fgsea_multilevel(...)`; §6.6 "Numba kernels that reduce floats in parallel
  MUST use a fixed chunking with a deterministic reduction order"; §6.6 seeds.

## What is ported, from where (licence)

fgsea is **MIT** (LICENCE: "YEAR: 2016-2019, COPYRIGHT HOLDER: Alexey Sergushichev"; licence field
"MIT + file LICENCE" on bioconductor.org, checked 2026-10-21). The port reads **fgsea 1.36.2** — the
version the golden env runs (M8-12) — from Bioconductor's git (`git.bioconductor.org/packages/fgsea`,
the `RELEASE_3_22` commit that carries 1.36.2): `R/fgsea.R` (`prepareStats`, `preparePathways`,
`calcGseaStat`, `fgseaSimpleImpl`), the simple-stage part of `R/fgseaMultilevel.R`, and
`src/fastGSEA.cpp` / `src/esCalculation.cpp` (`calcGseaStatCumulativeBatch`). Every ported file's
header carries fgsea's copyright line and the MIT permission notice; NOTICE.md gains a "Ported
code" row; `docs/methods/fgsea-port.md` records the commit hash and maps each Python function to its
origin. (The spec-writer read the current GitHub master, 1.39.4, to write this issue; the
implementer must diff it against 1.36.2 and port 1.36.2.)

## Public interface

`src/gippa/stats/fgsea.py` (pure; numba kernels in `src/gippa/stats/_fgsea_kernels.py`)
```python
@dataclass(frozen=True, slots=True)
class Prepared:
    stats: npt.NDArray[np.int64]                 # fgsea's integer-scaled |stats|, sorted
    order: npt.NDArray[np.int64]                 # input index of each sorted position
    names: tuple[str, ...]; members: tuple[npt.NDArray[np.int64], ...]   # 1-based ranks, sorted
    sizes: npt.NDArray[np.int64]; dropped: tuple[str, ...]

@dataclass(frozen=True, slots=True)
class SimpleResult:
    es: F64; nes: F64; p: F64; log2err: F64
    n_more_extreme: I64; mode_fraction: I64; needs_multilevel: npt.NDArray[np.bool_]
    leading_edge: tuple[npt.NDArray[np.int64], ...]      # input indices, fgsea's order

def prepare(stats: ArrayLike, labels: Sequence[str], sets: Mapping[str, Sequence[str]],
            min_size: int, max_size: int) -> Prepared: ...
def gsea_es(stats: npt.NDArray[np.int64], members: I64) -> tuple[float, I64]: ...
def fgsea_simple(prep: Prepared, n_perm: int, sample_size: int,
                 rng: np.random.Generator) -> SimpleResult: ...
```

1. **`prepare`** = `prepareStats` + `preparePathways` with `gseaParam = 1`, `scoreType = "std"`:
   non-finite stats raise; stats sorted **decreasing, stable** (R's radix `sort`: ties keep input
   order — verified against R by a tie case in the M8-12 inputs); `abs()`; then fgsea's **integer
   scaling**: `c = 2³⁰ / Σ|s|`, floored when ≥ 1, and `round(|s| · c)` stored as integers. ES is
   computed on these integers; without them the 10⁻⁹ ES tolerance cannot be met. Set members are
   mapped to ranks (names absent from `labels` dropped, duplicates removed), then sets outside
   `[min_size, max_size]` go to `dropped`. Duplicate labels raise.
2. **`gsea_es`** = `calcGseaStat` (`std`), on hit positions S (1-based, sorted), N stats, m = |S|:
   `NR = Σ r[S]`; `rCumSum = cumsum(r[S]) / NR` (or `(1…m)/m` when NR = 0);
   `tops = rCumSum − (S − (1…m)) / (N − m)`; `bottoms = tops − r[S]/NR` (or `− 1/m`). With
   `maxP = max(tops)`, `minP = min(bottoms)`: ES = maxP if maxP > −minP, minP if maxP < −minP, and
   **0 when they are equal**. Leading edge: `S[1 … argmax(tops)]` for ES > 0;
   `rev(S[argmin(bottoms) … m])` for ES < 0; **empty** on the tie. m = N raises (fgsea: "not
   defined when all genes are selected"). Accumulated in rank order, float64.
3. **`fgsea_simple`** = `fgseaSimpleImpl` via the cumulative batch: per permutation, one random
   sample of K = max set size ranks, whose prefixes of each set size give that size's random ES.
   Per set it counts `leEs`, `geEs` (random ES ≤ / ≥ the observed ES), `leZero`, `geZero` and sums
   `min(ES_rand, 0)`, `max(ES_rand, 0)`. Then:
   - `NES = ES / geZeroMean` if ES > 0, `ES / |leZeroMean|` if ES ≤ 0; NaN when that mean is 0;
   - `p = min((1 + nLeEs)/(1 + nLeZero), (1 + nGeEs)/(1 + nGeZero))` where NES exists;
   - `nMoreExtreme = nGeEs` if ES > 0 else `nLeEs`; `modeFraction = nGeZero` if ES ≥ 0 else
     `nLeZero`; sets with `modeFraction < 10` get p, NES and `log2err` NaN.
4. **Simple `log2err`** = `1/ln 2 · sqrt(trigamma(nMoreExtreme + 1) − trigamma(nPerm + 1))`.
5. **`needs_multilevel`** (ported literally): `crude = log2((nMoreExtreme + 1)/(nPerm + 1))`,
   Beta borders `left = log2(qbeta(0.025, nMoreExtreme, nPerm − nMoreExtreme + 1))`,
   `right = log2(qbeta(0.975, nMoreExtreme + 1, nPerm − nMoreExtreme))`,
   `simple_err = 0.5 · max(crude − left, right − crude)`, `multilevel_err = sqrt(floor(−log2 p₀ + 1)
   · (trigamma((sampleSize + 1)/2) − trigamma(sampleSize + 1))) / ln 2` with
   `p₀ = (nMoreExtreme + 1)/(nPerm + 1)`; `needs_multilevel = multilevel_err < simple_err`. Pass 1
   still reports the simple p for those sets (index decision 13); M8-11 replaces it.
6. **Randomness and parallelism.** One generator from the caller (stream `function.gsea`).
   Permutations run in fixed chunks of 100, each with a child seed from `SeedSequence.spawn`; the
   counts are integers and the two sums are added per chunk in chunk order, so results are
   identical for any thread count. R's stream cannot be reproduced (and need not be: M8-13's
   tolerances are Monte-Carlo-aware).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| maxP = −minP | ES = 0, empty leading edge (fgsea); NES uses the ES ≤ 0 branch |
| Σ_S \|s\| = 0 | equal hit weights 1/m (fgsea's NR = 0 branch) |
| a set covering all stats | `ValueError` from `gsea_es`; `prepare` never lets it through (`max_size` < N) |
| members absent from `labels` | dropped before the size filter |
| ties in stats | stable order (example 4); the methods note says ES depends on it, as fgsea warns |
| NES mean = 0 | NES and p NaN (`unestimated` in M8-05) |
| `n_perm < 1` or `sample_size < 3` | `ValueError` (fgsea silently raises `sampleSize` to 3; GIPPA's config is fixed at 101) |
| no set survives `prepare` | zero-length `SimpleResult` |

## Examples

ES examples use stats small enough that integer scaling is exact (c is an integer ≥ 1).

| # | Input | Expected |
|---|---|---|
| 1 | stats [5, 4, 3, 2, 1] (labels a–e), set {a, b} | tops [5/9, 1], bottoms [0, 5/9]: ES = **1.0**; leading edge [a, b] |
| 2 | same, set {d, e} | tops [−1/3, 0], bottoms [−1, −1/3]: ES = **−1.0**; leading edge **[e, d]** (fgsea reverses it) |
| 3 | stats [3, −1, 2, −2] (a–d), set {a, d} | sorted a, c, b, d; tops [0.6, 0], bottoms [0, −0.4]: ES = **0.6**, leading edge [a] |
| 4 | stats [2, 2, 1] labels [x, y, z], set {y} | stable order x, y, z: tops [0.5], bottoms [−0.5] → tie → ES = **0**, empty edge; with labels [y, x, z] ES = **1.0** |
| 5 | stats [1, 1, 1, 1] (a–d), set {a, c} | equal integer weights; tops [0.5, 0.5], bottoms [0, 0]: ES = **0.5**; leading edge **[a]** (the first maximum, R's `which.max`) |
| 6 | 1,000 standard-normal stats (seed 1), 20 random sets of size 15, `n_perm` 1000 | every p ≥ 1/1001; sign(NES) = sign(ES) where defined; same seed → identical arrays; 1 vs 8 numba threads → identical arrays |
| 7 | the top 50 of 2,000 normal stats as one set | `nMoreExtreme` = 0, p = 1/(1 + nGeZero), `needs_multilevel` **true** |
| 8 | a set whose random ES are almost all of one sign (`modeFraction` < 10) | p, NES, `log2err` NaN |
| 9 | `n_perm = 0` | `ValueError` |

## Acceptance criteria

- [ ] `tests/unit/stats/test_fgsea_es.py`: examples 1–5 to 10⁻¹², the leading edges and their
  order, the tie rule, the NR = 0 branch, the integer scaling (a hand-computed case where
  `round(|s|·c)` changes a weight).
- [ ] `tests/unit/stats/test_fgsea_simple.py`: examples 6–9; determinism across thread counts.
- [ ] Property: ES ∈ [−1, 1]; negating all stats negates ES (up to the stable tie order); the
  batch random ES of size k equals `gsea_es` on the same random prefix (exact).
- [ ] A slow pure-Python reference of `calcGseaStat` agrees on 200 random inputs to 10⁻¹².
- [ ] Performance (benchmark): 2,500 ranked TUs × 1,500 sets × 1,000 permutations in < 20 s on
  4 CPUs; entry in `benchmarks/baseline.json`.
- [ ] NOTICE.md and the module headers carry fgsea's MIT notice; `docs/methods/fgsea-port.md`
  exists with the commit and the function map.

## Calibration or attack work owed

- None alone; M8-15's `gsea` suite runs on this stage plus M8-11. The statistician checks that the
  simple p is a valid (+1, sign-conditional) permutation p.

## Hardening owed (pass 2)

- ES vs R fgsea to 10⁻⁹, NES within 0.02 (M8-13).
- M8-11 consumes `needs_multilevel`, `mode_fraction` and `denomProb = (modeFraction + 1)/(nPerm + 1)`.

## Questions for the lead

1. Port fgsea **1.36.2** (bioconda's version, so the golden env can run the same code) rather than
   Bioconductor's current 1.38.0? **Proposed default:** yes; the port and the golden reference must
   be one version, and bioconda is the only allowed R source (§3.2, §11.5).
2. `scoreType` `pos`/`neg` and `gseaParam ≠ 1`? **Proposed default:** not ported (R1); `prepare`
   accepts only the defaults.
3. Leading edge for the ES = 0 tie: fgsea returns NULL. **Proposed default:** keep it empty, and
   M8-05 writes an empty list (`NA` in the TSV).

## Lead decisions (2026-10-21)
1. Q1: port fgsea **1.36.2**, the version the golden env pins (M8-index lead decision 6). Read it from its official sources only: the Bioconductor 3.22 source tarball or the ctlab/fgsea tag. fgsea is MIT-licensed, so cite it in each file header and in NOTICE.md.
2. Q2: only `scoreType = "std"` and `gseaParam = 1` are ported. `prepare` rejects anything else.
3. Q3: the leading edge is empty when an ES = 0 tie occurs.
4. **Reproduce these 1.36.2 behaviours exactly:**
   - integer rounding of the statistics (scaled by 2³⁰/Σ) before ES is computed;
   - a sign-conditional simple-stage p;
   - an ES = 0 tie between the running-sum maximum and minimum;
   - the order of the leading edge for a negative ES.

   Check each against the 1.36.2 source, not master.
5. **Amends decision 4 (after the test step, 3d0d18b).** fgsea **1.36.2 has no integer scaling**: the 2³⁰/Σ rounding exists only on master (1.39.x), which is the version the designer read. GIPPA follows 1.36.2, the version the golden env runs, so statistics stay float |stats| and `Prepared.stats` is float64. This lead decision overrides spec item 1 for that one point.
6. The test-author's other resolutions are accepted:
   - the ES = 0 tie counts in both leZero and geZero;
   - p = (1 + nMoreExtreme)/(1 + modeFraction);
   - float ties are decided by IEEE rounding;
   - the `_fgsea_kernels.random_es_prefixes` hook is used for the prefix property.

   Owed to M8-13: R accumulates in long double, so near-ties may differ in rare cases, and the batch kernel's statEps weighting of zero stats is not pinned.
7. **After the statistician's review (395a99a):** M8-02 is approved. The port is valid and at level under random-set nulls (30/30 cells). Carried to M8-05:
   - **Tie order (lead-level; changes M8-05 rule 2 and M8-12's input):** TUs enter fgsea in a seeded random order, from `gippa.rng.for_step`, and the R golden gets the same order. The order is never by TU ID, which inflates the level to 0.165 under `abundance.missing = min` ties.
   - **Owed to Khoa (ADR):** the GSEA null under positional autocorrelation. With AR(1) φ = 0.5–0.8 and genome-adjacent sets, the rejection rate is 0.14–0.44 at α = 0.05. The options are a genome-rotation null, a TU-block permutation, or a documented scope limit. M8-15 gets this attack as a gate that must fail without the fix.
   - **Owed to Khoa (ADR):** the §8.2 95 % CI for NES. M8-05 leaves `ci_low`/`ci_high` empty, so it needs either an exemption or a CI method.
   - **Pass-1 power note:** rows that need multilevel keep a p ≥ ~1/501, so families with more than ~25 sets have little power until M8-11.
   - **The PR's calibration test has low power.** M8-15 switches it to the pooled criterion and adds a `calibration_full` arm at R = 300.

