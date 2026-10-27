# ADR-0044: The Wallenius golden reference is the exact distribution, not goseq
- Status: accepted-by-khoa (2026-10-23)
- Date: 2026-10-22
- Spec sections: §7.7.5 (criterion 1), §13 row M8, §11.4–11.5; M8-13 §1.1, question 3
- Needs Khoa: yes (changes the reference of a §7.7.5 criterion and adds a key to the protected
  `tests/calibration/tolerances.yaml`)

## Context
§7.7.5 criterion 1 compares GIPPA's Wallenius p-values with goseq's at relative 1e-6
(`golden_r.wallenius_rel_tol`). On the committed M8-12 golden data (goseq 1.62.0, BiasedUrn
2.0.12, default `precision = 1e-7`), 36 of the 300 tails differ by more than that (M8-13 F3). The
counts match for all 150 categories, and the odds match goseq's weight identity to 9e-16.

The statistician built a reference from the definition of Wallenius' distribution. No goseq or
BiasedUrn code was used or read. It was computed two independent ways:
- **R1, urn recursion.** The definition: N sequential draws without replacement, each red ball
  drawn with probability ω·r/(ω·r + w). Run in float64 and in 80-bit long double. Every term is
  a sum of products of positive factors, so nothing cancels. The two precisions agree to
  3.8e-15 relative on all 300 tails.
- **R2, exact rational integral.** f(x) = C(m₁,x) C(m₂,N−x) ∫₀¹ (1 − t^{ω/D})^x (1 − t^{1/D})^{N−x} dt,
  with D = ω(m₁ − x) + m₂ − N + x. The integral is expanded exactly as
  D Σⱼ (−1)ʲ C(x,j) B(ωj + D, N−x+1), where B(a, m+1) = m!/∏ᵢ₌₀ᵐ(a+i). This is exact in rationals
  because the float ω is a dyadic rational. The alternating sum is evaluated in fixed point with
  ≥ 200 surviving bits, so the cancellation is controlled.

  R2 equals R1 to every printed digit (rel < 1e-15) on all 9 tails checked:
  - CAT_002 under;
  - CAT_003 over;
  - CAT_004 over and under;
  - CAT_011 over;
  - CAT_012 over;
  - CAT_013 under;
  - CAT_105 over;
  - CAT_137 over.

## Evidence
On the 300 golden tails:

| | > 1e-6 against the exact reference | max relative error | max absolute error |
|---|---|---|---|
| GIPPA (`wallenius_tails`, scipy 1.18.1 pmf + `fsum`) | **0** | **7.1e-12** | 5.8e-12 |
| goseq / BiasedUrn | **36** | 0.49 (and 1.0 where it returns 0) | **2.6e-6** |

The failing tails, reference against goseq:

| tail | exact | GIPPA | goseq | goseq rel. error |
|---|---|---|---|---|
| CAT_004 over | 4.169558e-07 | 4.169558e-07 | 2.141283e-07 | 0.49 |
| CAT_011 over | 4.437464e-27 | 4.437464e-27 | 3.746333e-27 | 0.16 |
| CAT_012 over | 2.228390e-10 | 2.228390e-10 | 2.117715e-10 | 0.050 |
| CAT_013 under | 1.059794e-17 | 1.059794e-17 | 0 | 1 |
| CAT_137 over | 1.935414e-02 | 1.935414e-02 | 1.935152e-02 | 1.4e-4 |
| 31 others, p 0.1–0.99 | — | rel ≤ 5e-12 | — | 1.4e-6 – 1.8e-5 |

What goes wrong in goseq, from its outputs alone (goseq computes over = `dWNCHypergeo(k)` +
`pWNCHypergeo(k, lower.tail = FALSE)` and under = `pWNCHypergeo(k)`):
- **The point mass is accurate.** over + under − 1 equals the exact pmf(k) to 5.5e-8 relative
  on every category with pmf(k) > 1e-6. On CAT_004 and CAT_011, goseq's whole over-p equals
  pmf(k) to 1e-14 absolute.
- **The cdf has an absolute error of up to 2.6e-6.** The lower and the upper tail are off by the
  same amount with opposite signs, as they would be if upper = 1 − lower. This is about 26× the
  nominal `precision = 1e-7`.
- **Small tails inherit that absolute error and lose all relative accuracy:**
  - on CAT_004 and CAT_011 the P(X > k) part is exactly 0, so goseq's over-p is pmf(k) alone,
    about half of the true value;
  - on CAT_012 the P(X > k) part is only its first term;
  - on CAT_013 the lower tail is 0, although the true value is 1.06e-17.
- goseq is not uniformly bad. CAT_010's 6.07e-46 tail agrees to 5.5e-8. The error depends on
  BiasedUrn's internal method choice, which this ADR does not need to identify.

Conclusion: no GIPPA change can satisfy "relative 1e-6 against goseq" on these data. The
reference fails it, by up to a factor of 2, and it says 0 where the truth is 1e-17. **GIPPA is
right on every golden tail.**

## Decision
1. **Criterion 1's reference becomes the exact Wallenius distribution**, computed from the
   definition and never by the method GIPPA uses in the same regime:
   - **Tails < 1e-3:** by R2. These have few support points, so this is cheap. GIPPA's ADR-0045
     path for tails < 1e-6 is R1, so R2 keeps the comparison independent.
   - **Tails ≥ 1e-3:** by R1 in 80-bit precision. GIPPA uses scipy's pmf there.
   - **Larger tails:** taken as 1 − smaller tail + pmf(k), with pmf(k) by R2. This is exact
     because Σ f = 1, which is verified on ≥ 3 categories.
   - **R2 spot checks:** on ≥ 10 tails ≥ 1e-3, recorded.
   - **Cost:** R2 for one pmf is about x(N − x) big-integer products. Full R2 on CAT_001
     (n = 2,100) would take hours, which is why R2 is not used for every tail.

   The result is committed as `tests/data/golden/wallenius/exact_tails.tsv.gz` (category, tail,
   17-digit value) by `scripts/golden/make_wallenius_exact.py`. This is pure Python (`fractions`,
   integers, numpy long double), with no R and no scipy. It runs once on Slurm (≤ 4 CPUs, ≤ 1 h) and has a
   PROVENANCE entry. The comparison is
   `|p − p_exact| ≤ wallenius_rel_tol · p_exact` on all 300 tails, with the tolerance unchanged at
   1e-6. Today GIPPA passes this with a margin of 1.4e5.
2. **goseq stays as a second reference, at its demonstrated accuracy.** A new key
   `golden_r.wallenius_goseq_abs_tol: 1.0e-5` is added. The rule is
   `|p − p_goseq| ≤ max(wallenius_rel_tol · p_goseq, wallenius_goseq_abs_tol)`. The value is
   3.8× goseq's largest observed absolute error (2.6e-6), and it holds for p_goseq = 0 as well.
   It keeps the R-parity check, and it still catches gross faults. With an inverted ω (M8-13
   example 10), the goseq rule fails 288 of 300 tails (145 of 150 categories). The exact rule
   fails 293 tails (all 150 categories). The 5 categories the goseq rule misses are the
   deep-tail ones (CAT_010–013 and CAT_004), where only the exact rule has teeth. It replaces
   M8-13's interim `wallenius_abs_floor`; no separate floor is needed.
3. Counts and odds against goseq are unchanged (exact counts; ω to relative 1e-12).
4. **The vacuity floor** (M8-13 A5): ≥ 2 exact tails in (0, 1e-20) (CAT_010, CAT_011 today) and
   ≥ 1 in (1e-20, 1e-8) (CAT_012, CAT_013).

## Alternatives rejected
- **Loosen `wallenius_rel_tol`.** Any relative tolerance that goseq meets (> 0.49, or ∞ for
  the 0) makes criterion 1 vacuous in the tails it is meant to test.
- **Rerun goseq with a tighter BiasedUrn `precision`.** goseq's interface has no precision
  argument, so this needs a patched goseq call (a fork of the reference). Its accuracy would then have to be
  proved against the exact reference anyway.
- **Compare against R1 for every tail.** R1 is the algorithm proposed for GIPPA's far tails (ADR-0045).
  If ADR-0045 is adopted, R1 would test GIPPA against itself. R2 is a different method (closed
  form) and stays independent.

## Consequences
- M8-13 §1.1 item 4 is rewritten around the exact reference, plus the goseq rule above. A1 then
  depends only on this ADR.
- `tolerances.yaml` gains one key, applied in the test PR with `Gate-Change-Approved: ADR-0044`.
- §7.7.5 criterion 1 reads "against the exact distribution (and against goseq within its
  documented accuracy)". The spec text is not edited; this ADR overrides it.
- docs/methods/function.md "M8-13 golden" is updated: the 1e-6 comparison against goseq cannot
  hold *because of goseq*, while GIPPA's scipy pmf is within 1e-11 on every golden tail.
