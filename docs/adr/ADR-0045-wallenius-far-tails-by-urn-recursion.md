# ADR-0045: Wallenius far tails by the exact urn recursion
- Status: accepted-by-khoa (2026-10-23; lead recommendation)
- Date: 2026-10-22
- Spec sections: §7.7.3, §8.2–8.3; M8-01 lead decision 2 (hardening owed); docs/methods/function.md
  "Far-tail precision of scipy's pmf"
- Needs Khoa: no (no default, threshold, interface or dependency change; outputs become more
  accurate). The lead approves it.

## Context
`giap.stats.wallenius.wallenius_tails` sums scipy's `nchypergeom_wallenius.pmf`, which uses
Fog's precision constants. On the M8-13 golden data it is within 7.1e-12 of the exact
distribution on all 300 tails, including 6e-46 (ADR-0044). Elsewhere it is not.

The statistician checked 40 random urns (M 200–3,000, ω ∈ [e⁻¹, e]) against the exact urn
recursion in 80-bit precision (about 17,000 tails). The worst relative error per decade of the true
tail:

| true tail | worst relative error |
|---|---|
| 1e-5 – 1e-3 | ≤ 5.4e-10 |
| 1e-6 – 1e-5 | 6.7e-9 |
| 1e-7 – 1e-6 | 2.7e-8 |
| 1e-8 – 1e-7 | 6.2e-7 |
| 1e-9 – 1e-8 | 3.5e-6 |
| 1e-11 – 1e-10 | 1.5e-4 |
| 1e-14 – 1e-13 | 0.11 |
| 1e-19 – 1e-15 | 1 (**p = 0 returned**) |

The worst case is (M, n, N, ω) = (902, 151, 115, 1.98), lower tail. Its absolute error is a
near-constant ~8e-15, so scipy's pmf has an absolute floor rather than relative accuracy. This
matches the M8-01 measurement in docs/methods/function.md.

So GIAP reports p with relative error above 1e-6 below about 1e-8, and p = 0 for true values
near 1e-15. BH and elim decisions are unaffected (function.md "BH" and "elim"). But:
- the reported `p_classic` and `p_value` are wrong in the last decades;
- an exact-reference golden test (ADR-0044) passes only because the golden tails happen to be
  benign.

## Decision
1. `wallenius_tails` keeps scipy's pmf when the requested tail is ≥ `_EXACT_BELOW = 1e-6`
   (worst measured error there 6.7e-9). Below it, it recomputes the category's full pmf by the
   **exact urn recursion**. That recursion is the definition of Wallenius' distribution:
   - p₀(0) = 1;
   - each draw d moves mass from r red balls drawn to r + 1 with probability
     ω(m₁ − r) / (ω(m₁ − r) + m₂ − (d − r)).

   The tail is the `fsum` of the recursion's pmf.
2. Why the recursion: every term is a sum of products of positive factors, with no subtraction.
   So float64 is accurate to about N·ε relative per entry. It agreed with 80-bit precision to
   3.8e-15 on every golden tail and with the exact rational integral to < 1e-15 (ADR-0044 R2).
3. **Underflow.** If the float64 tail is < 1e-290, the recursion is rerun in log space
   (`np.logaddexp`) and exponentiated. A 0 can then only mean a true tail below float64's
   range (< 2e-308). No new column or flag is added.
4. **Odds range.** The recursion has no odds limit. `ODDS_RANGE` stays as it is for the scipy
   path. A later ADR may drop it once the recursion is the sole path.
5. **Cost.** The recursion is O(N · min(n, N)) per category, about 36 ms at N = 1,456 against
   1 ms for scipy. On the golden data it runs for 5 of 150 categories (+0.2 s). GIAP universes
   are a few thousand genes, so a GO run with ~5 % deep-tail categories adds seconds.

## Evidence required before acceptance (test-author, then implementer)
- A unit test against an independent oracle: the exact rational integral of ADR-0044, in
  `tests/helpers/`. It checks ≥ 10 urns, every k with a tail < 1e-6, and relative 1e-10.
- The xfail in `tests/unit/stats/test_wallenius.py` ("scipy's Wallenius pmf loses relative
  accuracy in some far tails") becomes a strict pass.
- `tests/calibration/test_wallenius_calibration.py` is unchanged and still passing.
- A benchmark of the function step on the largest M8 fixture, with < 10 % added wall time.

## Consequences
- The "Far-tail precision" section of docs/methods/function.md is replaced by the new error
  bound (≤ 1e-12 relative down to 1e-290).
- `enrichment_detection.tsv` no longer contains p = 0 from pmf truncation.
- It is independent of ADR-0044. ADR-0044 is needed for the gate whatever happens here. This ADR
  is needed for the reported numbers to be right where the golden data do not reach.
