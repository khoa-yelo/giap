# GIAP pilot report, pass 1 (arms A1 and A2)

Rendered by `validation/pilot/metrics.py render` from `metrics.json` (docs/specs/M10-05-pilot-metrics.md).
Every value below is followed by the file and key it comes from; none is typed by hand.
This template holds no numbers. Only aggregates appear here: no sequence, no gene list, no per-gene value.

## 1. Status

| §15.6 bullet | Status |
|---|---|
| GIAP runs end to end on A1, A2 and A3 | {{checklist.0.status}} |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | {{checklist.1.status}} |
| A1 and A2 detections and abundances are compared | {{checklist.2.status}} |
| Measured runtimes replace the estimates | {{checklist.3.status}} |
| TE2–TE5 are reported; failures explained | {{checklist.4.status}} |

GIAP commit: {{giap_commit}}. Spec SHA: {{spec}}.

## 2. Data and arms

PXD078459 wild-type runs (M10-01 manifest), DIA-NN 2.3.1 Academia (ADR-0051). A1 searched UniProt
K-12 (UP000000625) plus contaminants; A2 searched the Nissle strain database (`strain_db.faa`) plus
contaminants. Job IDs of M10-02/03/04 are in PROGRESS.md and the `SACCT` table.

## 3. Metrics

### Metric 1: identifications at 1 % FDR (union over runs)

Only the precursor count is at a controlled 1 % FDR (global precursor q-value). Peptide counts are
not FDR-controlled exactly. Protein groups are groups with ≥ 1 accepted precursor (no protein-level FDR control).
Per-run counts are "accepted at global q ≤ 0.01 among the rows DIA-NN reported for that run", not a per-run FDR.

| | A1 | A2 |
|---|---|---|
| precursors | {{m1.A1.union.precursors.all}} | {{m1.A2.union.precursors.all}} |
| peptides (all) | {{m1.A1.union.peptides.all}} | {{m1.A2.union.peptides.all}} |
| peptides (non-contaminant) | {{m1.A1.union.peptides.non_contaminant}} | {{m1.A2.union.peptides.non_contaminant}} |
| peptides (contaminant-only) | {{m1.A1.union.peptides.contaminant_only}} | {{m1.A2.union.peptides.contaminant_only}} |
| protein groups | {{m1.A1.union.protein_groups.all}} | {{m1.A2.union.protein_groups.all}} |
| protein groups with ≥ 2 accepted peptides | {{m1.A1.union.protein_groups_min2.all}} | {{m1.A2.union.protein_groups_min2.all}} |

### Metric 2: mapping shares (non-contaminant query peptides)

| | A1 | A2 |
|---|---|---|
| query peptides (U) | {{m2.A1.n_universe}} | {{m2.A2.n_universe}} |
| contaminant-only (C) | {{m2.A1.n_contaminant_only}} | {{m2.A2.n_contaminant_only}} |
| denominator (D) | {{m2.A1.n_denominator}} | {{m2.A2.n_denominator}} |
| exact | {{m2.A1.shares.exact}} | {{m2.A2.shares.exact}} |
| tolerant only | {{m2.A1.shares.tolerant}} | {{m2.A2.shares.tolerant}} |
| unmapped | {{m2.A1.shares.unmapped}} | {{m2.A2.shares.unmapped}} |

Genome-inconsistent peptides in A1 (`detection_summary.tsv`): rate {{m2.A1.genome_inconsistent.value}}, count {{m2.A1.genome_inconsistent.count}}; equal to tolerant + unmapped: {{m2.A1.genome_inconsistent.consistent}}.

### Metric 3: coverage of the strain CDS

Shares with nominal 95 % Wilson intervals, descriptive only. A1 and A2 are paired on the same genes,
so the arms are not compared by interval overlap; hypothesis 1 is the paired comparison.
"Strain-specific" is defined by peptide identity (no shared in-silico tryptic peptide with the K-12
digest, I/L collapsed), a proxy for "no K-12 counterpart"; it does not mean "undetectable by A1",
which can still detect such genes through tolerant placements.

| subset | A1 | A2 |
|---|---|---|
<!-- m3 rows -->

Subsets: all of G (cds, sorf), by replicon (`plasmid:<contig>` per plasmid; `unknown`: neither
chromosome nor plasmid), by region (prophage, bgc), by K-12 class (`no_peptides`: no in-silico
tryptic peptide); `pseudogene` is outside G.

### Metric 4: quantitative agreement (shared genes, from `giap compare`)

Spearman rho {{m4.spearman_rho}} [{{m4.ci_low}}, {{m4.ci_high}}] on {{m4.n_shared}} shared genes (TU-block bootstrap); Bland–Altman bias {{m4.ba_bias_log2}} log2; flags: {{m4.flags}}.

### Metric 5: stability of conclusions (Jaccard, from `giap compare`)

| family | Jaccard | note |
|---|---|---|
| islands | {{m5.islands.jaccard}} | {{m5.islands.note}} |
| operon flags | {{m5.operon_flags.jaccard}} | {{m5.operon_flags.note}} |
| enriched categories | {{m5.enriched_categories.jaccard}} | {{m5.enriched_categories.note}} |

### Metric 7: runtime and peak memory (pilot measurements on a shared cluster, not the benchmark)

| | A1 | A2 |
|---|---|---|
| wall time | {{m7.A1.wall_time_s}} | {{m7.A2.wall_time_s}} |
| sum of step wall times | {{m7.A1.sum_task_wall_s}} | {{m7.A2.sum_task_wall_s}} |
| max step peak RSS (one task, not the arm's total memory) | {{m7.A1.max_step_peak_rss_bytes}} | {{m7.A2.max_step_peak_rss_bytes}} |
| complete (no cached step; else totals are lower bounds) | {{m7.A1.runtime_complete}} | {{m7.A2.runtime_complete}} |

## 4. TE1 on real data (A2)

Every A2 target identification, true or false, is a substring of an annotated protein, so TE1 checks
that database, genome and mapping agree; it says nothing about identification error.
The expected misses are structural (N-terminal products of non-ATG starts, Sec sites); anything else is a mapping defect.

Share exact in frame at an own gene: {{te1.share_exact_inframe}} of {{te1.n_peptides}} peptides, nominal 95 % Wilson interval [{{te1.ci_low}}, {{te1.ci_high}}] (misses cluster by gene, so its coverage is optimistic; the interval is not a test).
Pass (point estimate at or above the threshold): {{te1.passed}}; lower end below the threshold: {{te1.ci_low_below_threshold}}.

Misses ({{te1.n_misses}}, first matching cause): non-ATG N-terminal {{te1.misses.non_atg_nterm}}, exact at another gene {{te1.misses.exact_other_gene}}, exact in another class {{te1.misses.exact_other_class}}, tolerant only {{te1.misses.tolerant_only}}, unplaced {{te1.misses.unplaced}}; masked exact in frame {{te1.exact_inframe_masked}}.
Non-ATG N-terminal products in the denominator: {{te1.non_atg_nterm.peptides}}. Peptides with no locus tag of `genes.tsv` among their sources: {{te1.unknown_source}}; with only NA q-values: {{te1.target_q_na}}.

Secondary, never the verdict: excluding non-ATG products {{te1.share_excluding_non_atg}}; without peptides whose accepted rows are all at the cut-off {{te1.share_without_q_at_cutoff}}.
Removed {{te1.share_without_q_at_cutoff.removed}} peptides (A2 q_crossings {{te1.share_without_q_at_cutoff.q_crossings}}); a pass that becomes a fail under this superset: {{te1.share_without_q_at_cutoff.pass_becomes_fail}}.
A1 diagnostic (exact in frame at any gene, no pass flag): {{te1_a1_diagnostic.share}}.

## 5. Hypotheses

Test family for this issue: {H1}, size 1, no BH.

**H1** (A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes): verdict {{h1.verdict}}; failed clauses: {{h1.failed_clauses}}.
`not_supported` means not supported by the pre-registered rule; `not_evaluable` means the rule could not be applied (no discordant TU).
A2-only genes b {{h1.b}}, A1-only genes c {{h1.c}}; A2-only share of discordant genes {{h1.prop_a2_only}}; coverage difference {{h1.coverage_diff}}.
These are nominal 95 % Wilson intervals; they describe the effects and are not the test (a
two-sided 95 % interval and the one-sided 0.05 test can disagree near the boundary).
TU-level test: {{h1.b_tu}} against {{h1.c_tu}} units ({{h1.mixed_tu}} mixed, left out), p_tu {{h1.p_tu}}, unit effect {{h1.unit_effect}}.
Gene-level p_gene {{h1.p_gene}} (assumes independent genes; decides nothing).
Share of A2-only genes in plasmid, prophage or BGC {{h1.mostly}} (base rate {{h1.mostly_base_rate}}); strain-specific share of A2-only genes {{h1.strain_specific_share}}.
The arms are two database searches of the same spectra, so discordance is largely fixed by sequence
divergence rather than sampled, and dependence at the scale of a region remains; p is not a replication probability.

**H3** (context conclusions agree between A1 and A2): verdict {{h3.verdict}}; rho criterion {{h3.criteria.rho.status}}; islands {{h3.criteria.islands.status}}; operon flags {{h3.criteria.operon_flags.status}}; enriched categories {{h3.criteria.enriched_categories.status}}.
The numbers are shown whatever the verdict: rho {{h3.criteria.rho.rho}} [{{h3.criteria.rho.ci_low}}, {{h3.criteria.rho.ci_high}}].

| family | n_a | n_b | n_shared | union | Jaccard | status |
|---|---|---|---|---|---|---|
| islands | {{h3.criteria.islands.n_a}} | {{h3.criteria.islands.n_b}} | {{h3.criteria.islands.n_shared}} | {{h3.criteria.islands.union}} | {{h3.criteria.islands.jaccard}} | {{h3.criteria.islands.status}} |
| operon flags | {{h3.criteria.operon_flags.n_a}} | {{h3.criteria.operon_flags.n_b}} | {{h3.criteria.operon_flags.n_shared}} | {{h3.criteria.operon_flags.union}} | {{h3.criteria.operon_flags.jaccard}} | {{h3.criteria.operon_flags.status}} |
| enriched categories | {{h3.criteria.enriched_categories.n_a}} | {{h3.criteria.enriched_categories.n_b}} | {{h3.criteria.enriched_categories.n_shared}} | {{h3.criteria.enriched_categories.union}} | {{h3.criteria.enriched_categories.jaccard}} | {{h3.criteria.enriched_categories.status}} |

**H2** and **H4**: deferred (they need A3, M10-07/M10-10 and the full validation).

## 6. Pending in pass 2

Metric 6 (M10-10), metric 8 (M10-07), TE2–TE5 and TE7 (M10-06), A3 (M10-08–11), benchmark (M10-12),
report-freeze re-runs (M10-14).

## 7. Deviations from the spec

- DIA-NN 2.3.1 Academia instead of 2.7.0 (ADR-0051).
- Inputs and metric files are not committed (Khoa's rule of 2026-10-23).
- Protein groups are groups with at least one accepted precursor, with no protein-level FDR control (metric 1).
- Strain-specific genes are defined by peptide identity, not orthology (metric 3).
- Region tests are not compared in pass 1, and the Jaccards are not restricted to genes detected in
  both arms (H3; a conservative stand-in, biased against H3).

## 8. Notes

<!-- notes -->
