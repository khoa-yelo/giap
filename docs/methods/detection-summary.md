# Basic detection summary (O − E) and QC rows

Implements `gippa.context.summary` (`basic_rows`, `qc_rows`, `detection_summary`) and
`gippa.io.summary_tsv.write_detection_summary` (GIPPA_SPEC.md §7.3.5(a), §8.2, §9.3; issue M4-16).
Statistician sign-off: see the end of this note (approved after the review fixes, 2026-10-15).

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
group carries the flag `ci_few_units` (on the `o_minus_e` row only) when it has fewer than 10 TUs
or when min(E, n_genes − E) < 10 (lead decision 4); a single-TU group has a degenerate CI equal to
the estimate. Calibration: tests/calibration/test_detection_summary_calibration.py
(coverage of 0 under a TU-comonotone null, where gene-level resampling under-covers).

What the CI does and does not cover: it is the sampling variability of Σ (d − p) over the group's
TUs **given** `p_detect`; the estimation error of the cross-fitted model is not propagated. For a
group that is (nearly) the whole eligible set (the main chromosome, `cog_category none` without
function annotation) Σ (d − p̂) ≈ 0 by construction (the logistic intercept fits calibration in the
large, up to cross-fitting noise), so its O − E is a sanity check, not a finding.

### Recorded coverage (review of PR M4-16; R = 2000 unless noted, B = 1000)

| Null | group | TUs | non-coverage (SE) |
|---|---|---|---|
| TU-comonotone, gippasim 400 genes, p = `p_detect_true` | replicon chr1 (gated) | 223 | 0.046 (0.005) |
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

- `id_mixups` (§7.5.1: same-database ALT sites without `start_codon`/`stop_codon`): NA with flag
  `requires_variants` until the variant module (M7) exists (lead decision 6).
- `engine_protein_mismatch_rate`: distinct C11 peptides flagged `engine_protein_mismatch` over
  distinct placed (C11) peptides; `n` = that denominator; NA without C11 rows. A QC proxy reported
  under its own name, not `id_mixups`.
- `genome_inconsistent_rate` (§7.5.1): among distinct C9 target, non-contaminant peptides with
  `search_db = other`, the share without an exact C11 placement (absent, or tolerant only); `n` =
  that denominator; NA without such peptides (e.g. all `search_db = same`) (lead decision 5).
  §7.5.1 defines it under `--genome-relation same`; under `related` tolerant placements are
  expected and the rate is not a mismatch signal (it is still computed; read it accordingly).
- `contaminant_peptides`, `contaminant_proteins`, `unmapped_peptides`, `unmapped_proteins`: C13
  rows by `label` and `kind`; NA without C13.

## Table

Sections in order tier, basic, detectability, qc, then group_type, group, metric. Parts without
`n`/`flags` (M4-12, M4-13 rows) get NA / empty. Floats are written with round-trip precision
(shortest repr) rather than `%.6g`, so CIs read back exactly (writer tests require 1e-9).

### Re-review coverage after the fixes (Slurm job 1850360; R = 1000, B = 1000, TU-comonotone null, sORF group)
The E rule counts expected detections in genes; under within-TU dependence the effective count
is the number of TUs that carry a detection (≈ E / TU size when detections are comonotone).
Unflagged groups (≥ 10 TUs and min(E, n − E) ≥ 10) at the boundary:

| TUs × genes per TU, p | E | non-coverage (SE) | flagged |
|---|---|---|---|
| 300 × 1, 0.034 | 10.2 | 0.069 (0.008) | 0 % |
| 100 × 2, 0.05 / 200 × 2, 0.026 | 10 / 10.4 | 0.112 / 0.119 (0.010) | 0 % |
| 100 × 3, 0.04 | 12 | 0.089 (0.009) | 0 % |
| 50 × 4, 0.05 / 40 × 5, 0.05 | 10 | 0.089 / 0.113 (0.009–0.010) | 0 % |
| 20 × 10, 0.05 | 10 | 0.356 (0.015) | 0 % |
| 100 × 3, 0.10 | 30 | 0.059 (0.007) | 0 % |

For realistic operon sizes (2–5 genes) the unflagged boundary covers 0.88–0.93, the same order as
the accepted 10–19-TU regime (≈ 0.91), and coverage recovers once E is a few multiples of 10
(0.059 at E = 30). Only long, uniformly low-p TUs (20 × 10 genes) fall far below; that is not a
realistic group. Recommendation (non-blocking, backlog): flag in TU units as well, e.g. fewer than
10 TUs with an observed detection (or non-detection) in the group.

## Statistician sign-off (M4-16)

Re-review at 338f49f (2026-10-15): **APPROVED.** Lead decisions 4–9 are implemented and tested:
the E rule and flag only on `o_minus_e` (4, 9); `genome_inconsistent_rate` over `search_db = other`
non-contaminant targets (5); `id_mixups` NA + `requires_variants` and
`engine_protein_mismatch_rate` (6; the 3/5 example of 338f49f is the correct distinct-peptide
count); ValueError naming genes with NaN `p_detect` or no `tu_id` (7); COG `-`/blank → `none`,
non-boolean `detected` raises (8); group-order mutant killed (9). Unit + integration tests 80/80
pass locally; implementer's calibration with the low-E attack (job 1849314) passed: the sORF E = 1,
E = 2 (operons) and n − E = 1 attacks are flagged in every replicate and their unflagged miss rate
is a confirmed fail. No p- or q-values, so §8.1 families do not apply; every CI is a TU-block CI
(§8.2). Non-blocking: the TU-unit flag above; the `genome_relation` caveat of
`genome_inconsistent_rate`; `hypothetical` still uses `astype(bool)` (a string `"false"` would read
as True; the same fail-closed check as `detected` would cover it).

Earlier review at db7dae0 (history):

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
