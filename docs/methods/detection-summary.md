# Basic detection summary (O − E) and QC rows

Implements `giap.context.summary` (`basic_rows`, `qc_rows`, `detection_summary`) and
`giap.io.summary_tsv.write_detection_summary` (GIAP_SPEC.md §7.3.5(a), §8.2, §9.3; issue M4-16).
Statistician sign-off: see the end of this note (changes requested).

## O − E against Σ p_detect

For a group G of eligible genes (`detection_status != partial_excluded`; single-peptide genes are
eligible and undetected in the primary analysis), O = Σ d_g and E = Σ p_g with the cross-fitted
`p_detect` of M4-13, so E[O − E] = 0 when the detectability model is calibrated (§8.2: never
against raw counts). Reported: `n_genes`, `n_detected`, `pct_detected` (0–100), `expected`,
`o_minus_e`, `o_over_e` (NA when E = 0) and `n_tu`; `n` = n_genes on every row.

## TU-block CI

O − E = Σ_TU s_t with s_t = Σ_{g ∈ G ∩ t} (d_g − p_g): only the group's genes inside a drawn TU
count. The 95 % percentile CI comes from `tu_block_bootstrap` (docs/methods/bootstrap.md) over the
group's TU sums, B = `stats.bootstrap` (1000). Groups are processed in sorted
(group_type, group) order from one `for_step("context_stats", seed)` generator, so the draws for a
group depend only on the TU counts of the groups before it, not on gene labels outside it. A
group with fewer than 10 TUs carries the flag `ci_few_units`; a single-TU group has a degenerate
CI equal to the estimate. Calibration: tests/calibration/test_detection_summary_calibration.py
(coverage of 0 under a TU-comonotone null, where gene-level resampling under-covers).

What the CI does and does not cover: it is the sampling variability of Σ (d − p) over the group's
TUs **given** `p_detect`; the estimation error of the cross-fitted model is not propagated. For a
group that is (nearly) the whole eligible set (the main chromosome, `cog_category none` without
function annotation) Σ (d − p̂) ≈ 0 by construction (the logistic intercept fits calibration in the
large, up to cross-fitting noise), so its O − E is a sanity check, not a finding.

### Recorded coverage (review of PR M4-16; R = 2000 unless noted, B = 1000)

| Null | group | TUs | non-coverage (SE) |
|---|---|---|---|
| TU-comonotone, giapsim 400 genes, p = `p_detect_true` | replicon chr1 (gated) | 223 | 0.046 (0.005) |
| same | feature_class cds / hypothetical | 219 / 84 | 0.054 / 0.060 |
| same, **gene-level** resampling (teeth) | chr1 | 400 genes | 0.166 |
| independent detection within TU | chr1 | 223 | 0.054 |
| comonotone, p / 10 (rare detection) | chr1 | 223 | 0.080 |
| comonotone, k random TUs of chr1 | chr1 | 3 / 5 / 8 / 10 / 15 / 20 / 40 | 0.25 / 0.15 / 0.11 / 0.094 / 0.089 / 0.085 / 0.064 |
| same | feature_class sorf / pseudogene | ≈ 8 / ≈ 4 | **0.90 / 0.40** |
| single-gene TUs, independent (R = 1000) | 50 TUs, E = 1 | 50 | **0.37** |
| same | E = 5 (50 / 100 TUs) | 50 / 100 | 0.065 / 0.11 |
| same | E = 10 / 20 / 50 | 100–200 / 200–400 / 100 | 0.07–0.08 / 0.04–0.07 / 0.04 |
| 2-gene comonotone TUs | E = 2 / 10 / 20 | 50 / 50–100 / 100–200 | 0.36 / 0.05–0.12 / 0.07 |

Gate cost (ADR-0015): at the measured rate 0.046 the fast gate (R = 100, ≤ 9 misses) fails a block
with probability 0.017 and the full gate (R = 1000, ≤ 63) 0.006; two-block false alarms 3 × 10⁻⁴
and 3 × 10⁻⁵ per night, negligible against the ≈ 0.04 nightly budget. Power against gene-level
resampling (0.166): 0.96 (fast), 1.0 (full).

**The TU count alone does not govern accuracy.** The percentile CI of a sum of residuals d − p
fails when few TUs carry the minority outcome (few detections when p is small, few non-detections
when p is near 1): with no detection in the group every TU sum is negative and the CI excludes 0
whatever the TU count. sORF and pseudogene groups (low p) and small COG letters are exactly this
case. A CI should be read as 95 % only when ≥ 10 TUs **and** ≥ 10 expected minority events
(min(E, n_genes − E) ≥ 10); below, coverage is 0.63–0.95 (table). CIs from 10–19 TUs cover ≈ 0.91
(c(n) of bootstrap.md), so `ci_few_units` at < 10 is a lenient threshold.

## Groups

`replicon` (each contig), `feature_class` (`cds`, `sorf`, `pseudogene`, and `hypothetical` = any
feature type with `hypothetical = true`), `cog_category`: a gene is counted once in **every**
letter of its `cog_category` (COG convention; `KT` counts in K and in T), genes without a
category in `none`. Groups without eligible genes are omitted.

## QC rows (`section = qc`, `run`/`all`)

- `id_mixups`: distinct C11 peptides with flag `engine_protein_mismatch`.
- `genome_inconsistent_rate`: among distinct C9 target peptides with `search_db = same`, the share
  without an exact C11 placement (absent, or tolerant only); `n` = that denominator; NA without
  such peptides (lead decision 1).
- `contaminant_peptides`, `contaminant_proteins`, `unmapped_peptides`, `unmapped_proteins`: C13
  rows by `label` and `kind`; NA without C13.

## Table

Sections in order tier, basic, detectability, qc, then group_type, group, metric. Parts without
`n`/`flags` (M4-12, M4-13 rows) get NA / empty. Floats are written with round-trip precision
(shortest repr) rather than `%.6g`, so CIs read back exactly (writer tests require 1e-9).

## Statistician sign-off (M4-16)

Reviewed on branch m4/detection-summary at db7dae0. **CHANGES REQUESTED** (see the PR comment):

1. Flag rule: extend the CI flag to rare outcomes (min(E, n_genes − E) < 10, or fewer than 10
   TUs with a detection / a non-detection), not only < 10 TUs (table above); add a gated attack
   null for a low-E group.
2. `genome_inconsistent_rate` uses the wrong peptide subset: §7.5.1 defines it for peptides
   searched against **another** strain (`search_db = other`, `--genome-relation same`); for
   `search_db = same` it is ≈ 0 by construction. Exclude `is_contaminant` peptides from the
   denominator. `id_mixups` (§7.5.1: same-db ALT sites without `start_codon`/`stop_codon`) is
   approximated by `engine_protein_mismatch` peptides until M5 lands; this needs an ADR or a
   renamed metric.
3. Fail closed: an eligible gene with NaN `p_detect` gives `expected` = NaN but a finite CI
   (pandas groupby skips NaN); an NA `tu_id` is dropped from the bootstrap (pandas 3 keeps NaN
   through `astype(str)`). Both should raise.

Verified correct: O − E and O/E, eligibility, multi-letter COG counting, TU sums restricted to the
group's genes, sorted-group draws, writer round trip; 13 / 15 manual mutants killed.
