# Basic detection summary (O − E) and QC rows

Implements `giap.context.summary` (`basic_rows`, `qc_rows`, `detection_summary`) and
`giap.io.summary_tsv.write_detection_summary` (GIAP_SPEC.md §7.3.5(a), §8.2, §9.3; issue M4-16).
Statistician sign-off: pending.

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
