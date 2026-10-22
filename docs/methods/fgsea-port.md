# fgsea port (G15; GIAP_SPEC.md §7.7.4)

GIAP's preranked GSEA is a Python/numba port of **fgsea 1.36.2** (MIT, © 2016-2019 Alexey
Sergushichev; licence text in NOTICE.md), the version the golden environment runs (M8-12).

## Source

- Bioconductor 3.22 source tarball `fgsea_1.36.2.tar.gz`, sha256
  `d6996ff50410d363b717b634baab22e91d1ee88c70355d71ffcb9f075a525f5a`; DESCRIPTION:
  `git_url https://git.bioconductor.org/packages/fgsea`, `git_branch RELEASE_3_22`,
  `git_last_commit b7862b7` (2026-01-05).
- Only scoreType `std` and gseaParam 1 are ported (M8-02 lead decision 2). No fgsea file is
  committed to GIAP.

## Function map (part I, M8-02)

| GIAP | fgsea 1.36.2 |
|---|---|
| `giap.stats.fgsea.prepare` | `R/fgsea.R` `preparePathwaysAndStats` (finite check, `sort(decreasing = TRUE)`, `abs`) and `R/util.R` `preparePathways` (`minSize >= 1`, `maxSize <= N - 1`, `unique(na.omit(fmatch(...)))`, duplicate-name check) |
| `giap.stats.fgsea.gsea_es`, `_fgsea_kernels.es_walk` | `R/fgsea.R` `calcGseaStat` (std, `returnLeadingEdge = TRUE`) |
| `_fgsea_kernels.random_es_prefixes` | `src/fastGSEA.cpp` `calcGseaStatCumulative` (std; statEps weighting of zero stats) |
| `_fgsea_kernels.random_es_batch` + the count loop of `fgsea_simple` | `src/fastGSEA.cpp` `calcGseaStatCumulativeBatch` |
| `giap.stats.fgsea.fgsea_simple` (NES, p, nMoreExtreme) | `R/fgsea.R` `fgseaSimpleImpl` |
| `fgsea_simple` (modeFraction, NaN rule, simple `log2err`, `needs_multilevel`) | `R/fgseaMultilevel.R` `fgseaMultilevel` up to the split into `dtSimpleFgsea` / `dtMultilevel`, and `multilevelError` |

## Behaviours reproduced and deviations

- **No integer scaling.** 1.36.2 computes ES on float |stats| (the 2³⁰/Σ rounding exists only on
  later masters; M8-02 lead decision 5).
- Stable decreasing order: ties keep input order, and ES depends on it (fgsea warns about ties).
- The ES = 0 tie (`maxP == -minP`) gives ES 0 and an empty leading edge; a negative ES's edge is
  `rev(S[which.min(bottoms) .. m])`; the positive edge stops at the first maximum.
- Random ES of exactly 0 count in both `leZero` and `geZero`.
- An even `sampleSize` is raised by one, as fgseaMultilevel does; `sample_size < 3` raises
  instead of being silently raised to 3.
- **Prefix walk instead of the convex hull.** fgsea's C++ finds each prefix's extremes with a
  block convex hull; GIAP walks every prefix directly (O(K²) per sample) with `calcGseaStat`'s
  arithmetic, so a random ES of size k equals `gsea_es` on the same ranks. R's random stream is
  not reproduced: each chunk of 100 permutations draws from a child of `rng.spawn`
  (`Generator.choice(N, K, replace = FALSE)`, random order like fgsea's shuffled `combination`).
- Owed to M8-13: R accumulates some sums in long double, so near-ties may be decided differently
  in rare cases; the statEps weighting of zero stats in the random kernel is not pinned by tests.

## Statistical review (M8-02, statistician, 2026-10-21)

Checked against the official 1.36.2 tarball (sha256 above): `R/fgsea.R` lines 666–695 (counts,
NES, p, nMoreExtreme) and `R/fgseaMultilevel.R` lines 100–200 (modeFraction, the `< 10` rule,
the Beta borders, `multilevelError`, the simple `log2err`). The port matches them line for line.

- **Validity of the simple p.** Under fgsea's null (the tested set is a uniform random set of its
  size), the observed ES and the `n_perm` random ES are exchangeable, and the observed and random
  ES both go through `es_walk`. Conditional on the sign of ES, (1 + nMoreExtreme)/(1 +
  modeFraction) is therefore a +1 permutation p among the same-signed draws, and the `min` in
  fgsea's formula always picks that side: for ES > 0 the le-side ratio is ≥ 1, and ES = 0 gives
  p = 1. Marginally, the p is a mixture of two conditionally valid p-values, so it is valid too.
  Dropping sets with `modeFraction < 10` depends only on the sign of ES and on the reference, so
  the sets that are kept stay conditionally valid, and BH over them (fgsea's `p.adjust` on non-NA)
  is not biased.
- **Evidence** (seeded simulations, not committed; scripts in
  `$GIAP_ROOT/stat/m802/`): N = 2,500 iid stats, 40 random sets per replicate, 500 replicates,
  `n_perm` 1,000, set sizes 10/15/50/200/500, under six nulls (normal, centred lognormal,
  uncentred lognormal, half-integer ties, 40 % zeros, t₂). All 30 cells are at level, with the
  rejection rate ≤ 0.0117 / 0.0521 / 0.1046 at α = 0.01 / 0.05 / 0.10 (20,000 p-values per cell).
  Exact ES = 0 ties: 0 in 36,000 random prefixes with k ≥ 10, except when every |stat| is equal.
  With equal |stats|, float rounding hides most exact ties, so the port gives ±x where exact
  arithmetic would give 0. That changes parity with R, not validity, because the observed and
  random ES share one rule.
- **The 1/(1 + modeFraction) floor.** With m tests per BH family, the floor (≈ 1/501 for balanced
  stats) is above the BH cut for one rejection, 0.05/m, once m ≳ 25. Pass 1 (simple stage only)
  can therefore reject only when several sets sit at the floor together. With 200 sets and 20
  shifted ones, power was 0 at shift 0.3 and 0.23 at shift 0.6 (mean FDP 0.019). Pass-1 GSEA is
  valid but has little power on large collections until M8-11 lands.
- **`log2err`** is fgsea's Beta(k + 1, n − k) posterior SD of log2 of the *unconditional* tail
  probability (it uses n = `n_perm`, not modeFraction). For the conditional p it understates the
  error slightly. It is ported as is, as a precision label.
- **Known structure that this null does not keep** (owed to M8-05/M8-15, not this port):
  (a) positional autocorrelation between TUs together with genomically clustered sets inflates the
  set-permutation p. With AR(1) φ = 0.5 TU scores and sets made of 5 clusters of 3 adjacent TUs,
  the rejection rate is 0.14 at α = 0.05; with one contiguous run of 15 TUs it is 0.22, and with
  φ = 0.8 it is 0.44. Scattered sets stay at level. (b) Stable tie order by TU ID (M8-05 rule 2)
  inflates level when many TUs tie at an *extreme* score, as with `abundance.missing = min`: with
  30 % of TUs tied at the 2nd percentile and contiguous 15-TU sets, the rate is 0.165 at
  α = 0.05. A seeded random input order restores 0.053. Ties at the centre (median imputation
  after centring) do not bite (0.048).
- **RNG.** `rng.spawn` children have distinct spawn keys and come from the SeedSequence, not from
  the parent's state. Draws made before the call therefore do not change the samples, and a
  second call on the same generator gets fresh children. Sampling without replacement is uniform
  (pair χ² p = 0.27 on N = 6, K = 3).

## Part II: the multilevel stage (M8-11)

`fgsea_multilevel(prep, simple, sample_size, eps, rng)` refines exactly the sets that
`simple.needs_multilevel` marks and that have a finite simple p. Every other field and every other set
comes back unchanged. Config: `sampleSize = 101`, `eps = 0`.

| GIAP | fgsea 1.36.2 |
|---|---|
| `fgsea_multilevel` (grouping by size, one ruler per sign, `denomProb`, `eps`, `p = min(1, cppMPval / denomProb)`) | `R/fgseaMultilevel.R` `fgseaMultilevel` after the split, `multilevelImpl` |
| `_multilevel_error` | `multilevelError` |
| `_beta_mean_log` | `betaMeanLog` |
| `_fgsea_kernels._ruler` | `EsRuler` (`extend`, `duplicateSamples`, `getPvalue`) |
| `_fgsea_kernels._perturb`, `_reaches` | `perturbate` |
| sign term of the p, `is_cp_ge_half` | `calcLogCorrection`, `isCpGeHalf` |

**Ruler.** A ruler draws `s` uniform random K-sets. Each set gets maxP, the positive deviation of
`calcGseaStat`'s walk. A round then does four things:

1. It records the median maxP as the level.
2. It keeps the top `h − 1 = (s − 1)/2` samples twice and the median sample once.
3. It runs MCMC swap moves until `s · K` moves have been accepted. A move swaps one member for one
   non-member, each chosen uniformly, and it is kept only when maxP stays ≥ the level. Each sample
   gets `max(1, K/10)` attempts per sweep.
4. It repeats until every target ES is below the median.

Negative ES runs the same ruler on the reversed |stats| with target |ES|. This is exact: reversing
the ranks maps the walk's minimum to minus the maximum.

**p of a target T.** T is resolved at the first round `r` whose median exceeds T. Let `c` be the
samples of that round with maxP > T, and `c₊` those of them whose ES is positive. Then:

- log cppMPval = r · betaMeanLog(h, s) + betaMeanLog(c + 1, s) + ψ(c₊ + 1) − ψ(c + 1);
- `isCpGeHalf` = 2c₊ ≥ c;
- p = min(1, cppMPval / denomProb);
- `log2err` = multilevelError(p, s), or NaN where `isCpGeHalf` is false;
- with eps > 0, a p < eps becomes eps with `log2err` NaN. Extension stops after
  ⌈−log₂(eps · min(denomProb) / 2)⌉ levels, using the minimum over the size group.

**Deviations, all GIAP's choices; part II was not line-checked against the 1.36.2 tarball.** The
implementer could not download the tarball in session. M8-13's golden runs and the statistician's
review own that check.

- `nPermSimple` is not a `SimpleResult` field. `denomProb = (modeFraction + 1)/(nPerm + 1)` takes
  nPerm from the simple `log2err` of the refined sets, because trigamma(nPerm + 1) =
  trigamma(nMoreExtreme + 1) − (ln 2 · log2err)². Newton steps invert it, then it is rounded. The
  recovery is exact to well under one permutation for nPerm ≤ 10⁷.
- **RNG.** Size groups get `rng.spawn` children in ascending size order (fgsea shuffles the
  groups). Each group spawns one child per sign. That child seeds a splitmix64 stream inside the
  numba kernel, which drives both the initial K-sets and the MCMC moves. Rulers run one per numba
  thread, so results do not depend on the thread count.
- **Stalled levels.** A level that accepts fewer than `s · K` moves within `10 · s · K` attempts
  (acceptance < 10 %) stops the extension, and so does reaching `MAX_LEVELS = 1000`. Targets still
  above the median are then resolved at that level, with `c` possibly 0, so p > 0 always. This only
  happens near the maximum ES of a size (the top-K set itself). There fgsea with eps = 0 keeps
  climbing at a cost of minutes per set. For example, the top 50 of 2,000 N(0, 1) stats stall
  near 1e-44 after about 140 levels.
- `_reaches` decides acceptance by walking from the last hit backwards and stopping once
  cum_j / nr < level. Its cumulative sums are nr minus the later weights, so they can differ from
  the forward walk by rounding. This matters only at exact ties with the level.
- **Not ported from master (1.39.x).** The `moveScale` argument and the C++-side `log2err` are
  absent, as the spec says. The full 1.36.2-vs-master diff is owed with the tarball check.
