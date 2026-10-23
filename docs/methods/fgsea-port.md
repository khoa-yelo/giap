# fgsea port (G15; GIAP_SPEC.md §7.7.4) — superseded by ADR-0038

> **Superseded.** Since M8-17b, `enrichment_gsea.tsv` comes from R fgsea 1.36.2
> (`fgsea::fgseaMultilevel`), run as a separate process (ADR-0038; docs/methods/function-gsea.md).
> The port below no longer produces any pipeline output. This note is kept as the history of
> M8-02 and as a reference: its function map names the fgsea 1.36.2 source behind each column of
> the recorded R calls that `ReplayEngine` replays (tests/helpers/fgsea_switch_data.py), and it
> documents `SIMPLE_ENGINE`, the in-process stand-in that unit tests inject in place of R. The hardening pass removes `giap.stats.fgsea` and its
> kernels unless Khoa keeps a no-R fallback (ADR-0038 item 6). Its 1/(1 + modeFraction) floor and
> the evidence below describe the port only, not the shipped engine.

GIAP's former preranked GSEA engine is a Python/numba port of **fgsea 1.36.2** (MIT, © 2016-2019
Alexey Sergushichev; licence text in NOTICE.md), the version the golden environment runs (M8-12).

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
