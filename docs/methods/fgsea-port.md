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
