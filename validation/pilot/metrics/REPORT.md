# GIAP pilot report (arms A1, A2 and A3)

Rendered by `validation/pilot/metrics.py render --final` from `metrics.json` (docs/specs/M10-05-pilot-metrics.md, docs/specs/M10-15-pilot-gate.md).
Every value below is followed by the file and key it comes from; none is typed by hand.
This template holds no numbers. Only aggregates appear here: no sequence, no gene list, no per-gene value.

## 1. Status

| §15.6 bullet | Status | Evidence |
|---|---|---|
| GIAP runs end to end on A1, A2 and A3 | pass (metrics.json › checklist.0) | metrics.json › evidence.exit_codes; reproduce.md › Evidence (metrics.json › checklist.0) |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | pass (metrics.json › checklist.1) | metrics.json › te1.share_exact_inframe, te1.passed (metrics.json › checklist.1) |
| A1 and A2 detections and abundances are compared | pass (metrics.json › checklist.2) | metrics.json › m4, m5 (metrics.json › checklist.2) |
| Measured runtimes replace the estimates | pending (M10-12) (metrics.json › checklist.3) | baseline.json › standard_run.wheat (not given, M10-12) (metrics.json › checklist.3) |
| TE2–TE5 are reported; failures explained | pending (M10-06) (metrics.json › checklist.4) | metrics.json › te_real.te2-te5.status (metrics.json › checklist.4) |

A failing bullet is explained in the notes below; no threshold was moved to make a bullet pass (R3).
GIAP commit: 6f969ee79d579d4bdb429be6f156679d644b2b51 (metrics.json › giap_commit). Spec SHA: e81c1559f9adba40cf8e53a43e4842a94c7492e5 (metrics.json › spec).
Submitted for Khoa's review (§14.8 item 4).

## 2. Data and arms

PXD078459 wild-type runs (M10-01 manifest), DIA-NN 2.3.1 Academia (ADR-0051). A1 searched UniProt
K-12 (UP000000625) plus contaminants; A2 searched the Nissle strain database (`strain_db.faa`) plus
contaminants. A3 is de novo: Cascadia 0.0.7 on the same runs or, when M10-08's verdict rules it
out, Casanovo 5.2.1 on PXD053727 / UTI89 (fallback, a different strain).
A3 branch of the runbook: pending (render --final --runbook); M10-08 verdict pass (metrics.json › evidence.a3_feasibility) (reason: NA (metrics.json › evidence.a3_feasibility)).
Job IDs are in `reproduce.md` (Evidence table) and in the status table above.

## 3. Metrics

### Metric 1: identifications at 1 % FDR (union over runs)

Only the precursor count is at a controlled 1 % FDR (global precursor q-value). Peptide counts are
not FDR-controlled exactly. Protein groups are groups with ≥ 1 accepted precursor (no protein-level FDR control).
Per-run counts are "accepted at global q ≤ 0.01 among the rows DIA-NN reported for that run", not a per-run FDR.

| | A1 | A2 |
|---|---|---|
| precursors | 28639 (metrics.json › m1.A1.union.precursors) | 33418 (metrics.json › m1.A2.union.precursors) |
| peptides (all) | 23792 (metrics.json › m1.A1.union.peptides) | 27952 (metrics.json › m1.A2.union.peptides) |
| peptides (non-contaminant) | 23704 (metrics.json › m1.A1.union.peptides) | 27827 (metrics.json › m1.A2.union.peptides) |
| peptides (contaminant-only) | 88 (metrics.json › m1.A1.union.peptides) | 125 (metrics.json › m1.A2.union.peptides) |
| protein groups | 2798 (metrics.json › m1.A1.union.protein_groups) | 3192 (metrics.json › m1.A2.union.protein_groups) |
| protein groups with ≥ 2 accepted peptides | 2433 (metrics.json › m1.A1.union.protein_groups_min2) | 2765 (metrics.json › m1.A2.union.protein_groups_min2) |

### Metric 2: mapping shares (non-contaminant query peptides)

| | A1 | A2 |
|---|---|---|
| query peptides (U) | 23792 (metrics.json › m2.A1) | 27952 (metrics.json › m2.A2) |
| contaminant-only (C) | 88 (metrics.json › m2.A1) | 125 (metrics.json › m2.A2) |
| denominator (D) | 23704 (metrics.json › m2.A1) | 27827 (metrics.json › m2.A2) |
| exact | 0.9898 [0.9884, 0.991] (metrics.json › m2.A1.shares) | 0.9985 [0.998, 0.9989] (metrics.json › m2.A2.shares) |
| tolerant only | 0.005948 [0.005046, 0.007011] (metrics.json › m2.A1.shares) | 0 [0, 0.000138] (metrics.json › m2.A2.shares) |
| unmapped | 0.004261 [0.003508, 0.005174] (metrics.json › m2.A1.shares) | 0.001509 [0.001117, 0.002039] (metrics.json › m2.A2.shares) |

Genome-inconsistent peptides in A1 (`detection_summary.tsv`): rate 0.0616 (metrics.json › m2.A1.genome_inconsistent), count 1540 (metrics.json › m2.A1.genome_inconsistent); equal to tolerant + unmapped: no (metrics.json › m2.A1.genome_inconsistent).

### Metric 3: coverage of the strain CDS

Shares with nominal 95 % Wilson intervals, descriptive only. A1 and A2 are paired on the same genes,
so the arms are not compared by interval overlap; hypothesis 1 is the paired comparison.
"Strain-specific" is defined by peptide identity (no shared in-silico tryptic peptide with the K-12
digest, I/L collapsed), a proxy for "no K-12 counterpart"; it does not mean "undetectable by A1",
which can still detect such genes through tolerant placements.

| subset | A1 | A2 |
|---|---|---|
| all | 0.5691 [0.5548, 0.5832] (metrics.json › m3.A1.subsets.all) | 0.6477 [0.6339, 0.6612] (metrics.json › m3.A2.subsets.all) |
| chromosome | 0.5708 [0.5566, 0.5849] (metrics.json › m3.A1.subsets.chromosome) | 0.6478 [0.634, 0.6614] (metrics.json › m3.A2.subsets.chromosome) |
| plasmid | 0.1111 [0.03102, 0.328] (metrics.json › m3.A1.subsets.plasmid) | 0.6111 [0.3862, 0.7969] (metrics.json › m3.A2.subsets.plasmid) |
| plasmid:CP171243.1 | 0.25 [0.04559, 0.6994] (metrics.json › m3.A1.subsets.plasmid:CP171243.1) | 0.25 [0.04559, 0.6994] (metrics.json › m3.A2.subsets.plasmid:CP171243.1) |
| plasmid:CP171244.1 | 0 [0, 0.3903] (metrics.json › m3.A1.subsets.plasmid:CP171244.1) | 0.6667 [0.3, 0.9032] (metrics.json › m3.A2.subsets.plasmid:CP171244.1) |
| plasmid:CP171245.1 | 0.125 [0.02242, 0.4709] (metrics.json › m3.A1.subsets.plasmid:CP171245.1) | 0.75 [0.4093, 0.9285] (metrics.json › m3.A2.subsets.plasmid:CP171245.1) |
| unknown | NA (empty_denominator) (metrics.json › m3.A1.subsets.unknown) | NA (empty_denominator) (metrics.json › m3.A2.subsets.unknown) |
| prophage | 0.02206 [0.00753, 0.06285] (metrics.json › m3.A1.subsets.prophage) | 0.2647 [0.1978, 0.3446] (metrics.json › m3.A2.subsets.prophage) |
| bgc | NA (empty_denominator) (metrics.json › m3.A1.subsets.bgc) | NA (empty_denominator) (metrics.json › m3.A2.subsets.bgc) |
| strain_specific | 0 [0, 0.004136] (metrics.json › m3.A1.subsets.strain_specific) | 0.3405 [0.3107, 0.3717] (metrics.json › m3.A2.subsets.strain_specific) |
| k12_shared | 0.7107 [0.696, 0.725] (metrics.json › m3.A1.subsets.k12_shared) | 0.725 [0.7105, 0.7391] (metrics.json › m3.A2.subsets.k12_shared) |
| no_peptides | 0 [0, 0.2775] (metrics.json › m3.A1.subsets.no_peptides) | 0 [0, 0.2775] (metrics.json › m3.A2.subsets.no_peptides) |
| pseudogene | 0 [0, 0.5615] (metrics.json › m3.A1.subsets.pseudogene) | 0 [0, 0.5615] (metrics.json › m3.A2.subsets.pseudogene) |

Subsets: all of G (cds, sorf), by replicon (`plasmid:<contig>` per plasmid; `unknown`: neither
chromosome nor plasmid), by region (prophage, bgc), by K-12 class (`no_peptides`: no in-silico
tryptic peptide); `pseudogene` is outside G.

### Metric 4: quantitative agreement (shared genes, from `giap compare`)

Spearman rho 0.978 (metrics.json › m4) [0.9729 (metrics.json › m4), 0.9825 (metrics.json › m4)] on 2625 (metrics.json › m4) shared genes (TU-block bootstrap); Bland–Altman bias -0.192 (metrics.json › m4) log2; flags: NA (metrics.json › m4).

### Metric 5: stability of conclusions (Jaccard, from `giap compare`)

| family | Jaccard | note |
|---|---|---|
| islands | NA (metrics.json › m5.islands) | none in either arm (metrics.json › m5.islands) |
| operon flags | NA (metrics.json › m5.operon_flags) | none in either arm (metrics.json › m5.operon_flags) |
| enriched categories | 0.551 (metrics.json › m5.enriched_categories) | NA (metrics.json › m5.enriched_categories) |

### Metric 6: novel ORFs and variants by class and confidence, A2/A3 overlap (M10-10)

Descriptive, unestimated: counts and Jaccards carry no error estimate (R4); the two shares carry
nominal 95 % Wilson intervals that describe the sample of calls, not an FDR. Per arm, novel ORFs
(`novel_orfs.tsv`) are counted as called, not called with an unestimated class FDR, and not called.

| class | A1 called | A1 unestimated | A1 not_called | A2 called | A2 unestimated | A2 not_called | A3 called | A3 unestimated | A3 not_called |
|---|---|---|---|---|---|---|---|---|---|
| nterm_extension | 0 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 15 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) | 1 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) |
| cterm_extension | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) |
| intergenic | 1 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 7 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) | 1 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) |
| cds_outframe | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) | 1 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) |
| antisense | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A3.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A3.novel_orfs.by_class.antisense) | 2 (metrics.json › m6.A3.novel_orfs.by_class.antisense) |
| pseudogene | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) |
| sorf | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) |
| total | 1 (metrics.json › m6.A1.novel_orfs.totals) | 22 (metrics.json › m6.A1.novel_orfs.totals) | 0 (metrics.json › m6.A1.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A3.novel_orfs.totals) | 0 (metrics.json › m6.A3.novel_orfs.totals) | 5 (metrics.json › m6.A3.novel_orfs.totals) |

ALT calls (`variants.tsv`, `call = ALT`) by type and confidence; for A3 every class error rate is
expected unestimated (no spectrum pairing).

| ALT type | A1 high | A1 candidate | A2 high | A2 candidate | A3 high | A3 candidate |
|---|---|---|---|---|---|---|
| substitution | 5 (metrics.json › m6.A1.variants.alt.substitution) | 131 (metrics.json › m6.A1.variants.alt.substitution) | 0 (metrics.json › m6.A2.variants.alt.substitution) | 0 (metrics.json › m6.A2.variants.alt.substitution) | 0 (metrics.json › m6.A3.variants.alt.substitution) | 42 (metrics.json › m6.A3.variants.alt.substitution) |
| swap | 0 (metrics.json › m6.A1.variants.alt.swap) | 0 (metrics.json › m6.A1.variants.alt.swap) | 0 (metrics.json › m6.A2.variants.alt.swap) | 0 (metrics.json › m6.A2.variants.alt.swap) | 1 (metrics.json › m6.A3.variants.alt.swap) | 124 (metrics.json › m6.A3.variants.alt.swap) |
| block | 0 (metrics.json › m6.A1.variants.alt.block) | 0 (metrics.json › m6.A1.variants.alt.block) | 0 (metrics.json › m6.A2.variants.alt.block) | 0 (metrics.json › m6.A2.variants.alt.block) | 0 (metrics.json › m6.A3.variants.alt.block) | 28 (metrics.json › m6.A3.variants.alt.block) |
| frameshift | 0 (metrics.json › m6.A1.variants.alt.frameshift) | 0 (metrics.json › m6.A1.variants.alt.frameshift) | 0 (metrics.json › m6.A2.variants.alt.frameshift) | 0 (metrics.json › m6.A2.variants.alt.frameshift) | 0 (metrics.json › m6.A3.variants.alt.frameshift) | 0 (metrics.json › m6.A3.variants.alt.frameshift) |

A1: mixed calls 1 (metrics.json › m6.A1.variants); class error rate numeric 0 (metrics.json › m6.A1.variants.class_error), unestimated 137 (metrics.json › m6.A1.variants.class_error).
A2: mixed calls 0 (metrics.json › m6.A2.variants); class error rate numeric 0 (metrics.json › m6.A2.variants.class_error), unestimated 0 (metrics.json › m6.A2.variants.class_error).
A3: mixed calls 53 (metrics.json › m6.A3.variants); class error rate numeric 0 (metrics.json › m6.A3.variants.class_error), unestimated 248 (metrics.json › m6.A3.variants.class_error).

A2/A3 overlap. Novel ORFs key on the stop codon (contig, strand, end), which does not depend on the
peptides; ALT sites key on (contig, codon start, strand, ALT residue) as `giap compare` does. The
`high` set is each arm's own high-confidence ALT rows.

| set | A2 | A3 | shared | union | Jaccard | A3 confirmed by A2 | A2 found by A3 |
|---|---|---|---|---|---|---|---|
| novel ORFs, called | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | NA (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.called) |
| novel ORFs, any row | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 5 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 5 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 [0, 0.4345] (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) |
| ALT sites | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 195 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 195 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 [0, 0.01932] (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.alt_sites.alt) |
| ALT sites, high confidence | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 1 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 1 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 [0, 0.7935] (metrics.json › m6.overlap_a2_a3.alt_sites.high) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.alt_sites.high) |

Novel ORFs called in both arms with the same class: 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called).
`giap compare`'s novel-ORF Jaccard (keyed on `orf_id`, which can split one ORF): NA (metrics.json › m6.overlap_a2_a3.novel_orfs).
High-confidence ALT Jaccard: 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high); called novel-ORF Jaccard: NA (metrics.json › m6.overlap_a2_a3.novel_orfs.called).

### Metric 8: free entrapment FDP on A1 (M10-07)

| level | combined FDP | lower bound | assessment |
|---|---|---|---|
| peptide (primary) | pending (M10-07) | pending (M10-07) | pending (M10-07) |
| precursor | pending (M10-07) | pending (M10-07) | pending (M10-07) |

## 4. TE1 on real data (A2)

Every A2 target identification, true or false, is a substring of an annotated protein, so TE1 checks
that database, genome and mapping agree; it says nothing about identification error.
The expected misses are structural (N-terminal products of non-ATG starts, Sec sites); anything else is a mapping defect.

Share exact in frame at an own gene: 0.9985 (metrics.json › te1) of 27827 (metrics.json › te1) peptides, nominal 95 % Wilson interval [0.998 (metrics.json › te1), 0.9989 (metrics.json › te1)] (misses cluster by gene, so its coverage is optimistic; the interval is not a test).
Pass (point estimate at or above the threshold): yes (metrics.json › te1); lower end below the threshold: no (metrics.json › te1).

Misses (42 (metrics.json › te1), first matching cause): non-ATG N-terminal 42 (metrics.json › te1.misses), exact at another gene 0 (metrics.json › te1.misses), exact in another class none (metrics.json › te1.misses), tolerant only 0 (metrics.json › te1.misses), unplaced none (metrics.json › te1.misses); masked exact in frame 18 (metrics.json › te1).
Non-ATG N-terminal products in the denominator: 42 (metrics.json › te1.non_atg_nterm). Peptides with no locus tag of `genes.tsv` among their sources: 0 (metrics.json › te1); with only NA q-values: 0 (metrics.json › te1).

Secondary, never the verdict: excluding non-ATG products 1 [0.9999, 1] (metrics.json › te1); without peptides whose accepted rows are all at the cut-off NA (no_crossings) (metrics.json › te1).
Removed NA (metrics.json › te1.share_without_q_at_cutoff) peptides (A2 q_crossings 0 (metrics.json › te1.share_without_q_at_cutoff)); a pass that becomes a fail under this superset: no (metrics.json › te1.share_without_q_at_cutoff).
A1 diagnostic (exact in frame at any gene, no pass flag): 0.9882 [0.9868, 0.9895] (metrics.json › te1_a1_diagnostic).

## 5. True examples (TE2–TE5) and determinism (TE7)

From `te_real.json` and `te7.json` (docs/specs/M10-06-te-real-and-te7.md); every definition was fixed
before the real data were read (R3). Shares carry nominal 95 % Wilson intervals, descriptive only.
Every TE that did not pass has a `## Diagnosis: TE<n>` section in the notes below, naming its cause class.

| TE | Status | Headline |
|---|---|---|
| TE2: no unexplained `high` ALT call in A2 | pending (M10-06) (metrics.json › te_real.te2) | violations NA (metrics.json › te_real.te2.counts) |
| TE3: ribosomal proteins detected (≥ 90 %) and in the top quartile | pending (M10-06) (metrics.json › te_real.te3) | detected NA (metrics.json › te_real.te3.detection); median percentile rank NA (metrics.json › te_real.te3.abundance) |
| TE4: operon stoichiometry (iBAQ, molar scale) | pending (M10-06) (metrics.json › te_real.te4) | log2 RplL/RplJ NA (metrics.json › te_real.te4.pairs.rplL/rplJ) in [1, 3]; log2 AtpA/AtpD NA (metrics.json › te_real.te4.pairs.atpA/atpD) in [−1, 1] |
| TE5: trypsin and keratin peptides in the contaminant bucket (A1 and A2) | pending (M10-06) (metrics.json › te_real.te5) | in the bucket: A1 NA (metrics.json › te_real.te5.A1), A2 NA (metrics.json › te_real.te5.A2) |
| TE7: a second A2 run gives byte-identical tables | pending (M10-06) (metrics.json › te7.giap) | gate files compared NA (metrics.json › te7.giap) |

TE7 is decided by GIAP's own tables only; the DIA-NN layer (a second search of the same spectra) is an
upstream finding, never the verdict: pending (M10-06) (metrics.json › te7.diann). Raw data to GIAP tables reproduce byte for byte: NA (metrics.json › te7).
A `known_lapack` difference is still a fail (M4-20 decision 9; writer rounding is the decision 10 backlog item).
TE7 is not a §15.6 bullet: a TE7 fail is shown here and in section 9 with its diagnosis and listed as an open item.

## 6. Hypotheses

Test family for this issue: {H1}, size 1, no BH.

**H1** (A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes): verdict not_supported (metrics.json › h1); failed clauses: mostly (metrics.json › h1).
`not_supported` means not supported by the pre-registered rule; `not_evaluable` means the rule could not be applied (no discordant TU).
A2-only genes b 409 (metrics.json › h1), A1-only genes c 40 (metrics.json › h1); A2-only share of discordant genes 0.9109 [0.881, 0.9339] (metrics.json › h1); coverage difference 0.0786 (metrics.json › h1).
These are nominal 95 % Wilson intervals; they describe the effects and are not the test (a
two-sided 95 % interval and the one-sided 0.05 test can disagree near the boundary).
TU-level test: 267 (metrics.json › h1) against 36 (metrics.json › h1) units (3 (metrics.json › h1) mixed, left out), p_tu 4.7e-45 (metrics.json › h1), unit effect 0.8812 [0.8399, 0.9129] (metrics.json › h1).
Gene-level p_gene 1.9e-78 (metrics.json › h1) (assumes independent genes; decides nothing).
Share of A2-only genes in plasmid, prophage or BGC 0.1027 [0.07687, 0.1359] (metrics.json › h1) (base rate 0.03282 [0.02809, 0.03831] (metrics.json › h1)); strain-specific share of A2-only genes 0.7702 [0.727, 0.8083] (metrics.json › h1).
The arms are two database searches of the same spectra, so discordance is largely fixed by sequence
divergence rather than sampled, and dependence at the scale of a region remains; p is not a replication probability.

**H3** (context conclusions agree between A1 and A2): verdict supported (metrics.json › h3); rho criterion holds (metrics.json › h3.criteria.rho); islands not_evaluable (metrics.json › h3.criteria.islands); operon flags not_evaluable (metrics.json › h3.criteria.operon_flags); enriched categories holds (metrics.json › h3.criteria.enriched_categories).
The numbers are shown whatever the verdict: rho 0.978 (metrics.json › h3.criteria.rho) [0.9729 (metrics.json › h3.criteria.rho), 0.9825 (metrics.json › h3.criteria.rho)].

| family | n_a | n_b | n_shared | union | Jaccard | status |
|---|---|---|---|---|---|---|
| islands | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | NA (metrics.json › h3.criteria.islands) | not_evaluable (metrics.json › h3.criteria.islands) |
| operon flags | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | NA (metrics.json › h3.criteria.operon_flags) | not_evaluable (metrics.json › h3.criteria.operon_flags) |
| enriched categories | 425 (metrics.json › h3.criteria.enriched_categories) | 287 (metrics.json › h3.criteria.enriched_categories) | 253 (metrics.json › h3.criteria.enriched_categories) | 459 (metrics.json › h3.criteria.enriched_categories) | 0.551 (metrics.json › h3.criteria.enriched_categories) | holds (metrics.json › h3.criteria.enriched_categories) |

**H2** (the strain-specific database lowers genome-inconsistent identifications): deferred (full validation, FDRBench); shown side by side only: A1 genome-inconsistent share 0.0616 (metrics.json › m2.A1.genome_inconsistent) and the free entrapment estimate (peptide level) pending (M10-07).

**H4** (genome-anchored de novo precision against the precision measured on A2): verdict consistent_with_lower_bound (metrics.json › h4).
DIA gives no spectrum pairing, so the only measure against A2 is sequence-level: the share of a
bin's de novo queries that A2 identified at 1 % FDR (I = L), a lower bound on the bin's precision.
The data can therefore count against H4 but never confirm it. A bin is reportable when GIAP reports
its precision and it holds at least 20 queries; the verdict uses the Bonferroni-adjusted relation
over the 10 (metrics.json › h4) reportable bins (`relation_adj`), the 95 % relation is descriptive.

| bin | n | A2-confirmed | genome precision | CI low | CI high | c_value | U (rescaled, descriptive) | L (adjusted) | relation 95 % | relation adjusted | nested share | A2 confirmation among searchable peptides (descriptive) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (metrics.json › h4.bins.0) | 137 (metrics.json › h4.bins.0) | 0.2555 [0.1898, 0.3345] (metrics.json › h4.bins.0) | 0.293 (metrics.json › h4.bins.0.genome) | 0.2201 (metrics.json › h4.bins.0.genome) | 0.3669 (metrics.json › h4.bins.0.genome) | 0.995 (metrics.json › h4.bins.0.genome) | 0.399 (metrics.json › h4.bins.0) | 0.2555 [0.1662, 0.3713] (metrics.json › h4.bins.0) | consistent (metrics.json › h4.bins.0) | consistent (metrics.json › h4.bins.0) | 0.0219 [0.007475, 0.0624] (metrics.json › h4.bins.0) | 0.9722 [0.8583, 0.9951] (metrics.json › h4.bins.0) |
| 1 (metrics.json › h4.bins.1) | 136 (metrics.json › h4.bins.1) | 0.3603 [0.2845, 0.4438] (metrics.json › h4.bins.1) | 0.346 (metrics.json › h4.bins.1.genome) | 0.2577 (metrics.json › h4.bins.1.genome) | 0.4348 (metrics.json › h4.bins.1.genome) | 0.995 (metrics.json › h4.bins.1.genome) | 0.473 (metrics.json › h4.bins.1) | 0.3603 [0.2553, 0.4806] (metrics.json › h4.bins.1) | consistent (metrics.json › h4.bins.1) | consistent (metrics.json › h4.bins.1) | 0 [0, 0.02747] (metrics.json › h4.bins.1) | 1 [0.9273, 1] (metrics.json › h4.bins.1) |
| 2 (metrics.json › h4.bins.2) | 136 (metrics.json › h4.bins.2) | 0.3603 [0.2845, 0.4438] (metrics.json › h4.bins.2) | 0.368 (metrics.json › h4.bins.2.genome) | 0.2839 (metrics.json › h4.bins.2.genome) | 0.4479 (metrics.json › h4.bins.2.genome) | 0.995 (metrics.json › h4.bins.2.genome) | 0.483 (metrics.json › h4.bins.2) | 0.3603 [0.2553, 0.4806] (metrics.json › h4.bins.2) | consistent (metrics.json › h4.bins.2) | consistent (metrics.json › h4.bins.2) | 0.007353 [0.001299, 0.04047] (metrics.json › h4.bins.2) | 1 [0.9273, 1] (metrics.json › h4.bins.2) |
| 3 (metrics.json › h4.bins.3) | 137 (metrics.json › h4.bins.3) | 0.4745 [0.3927, 0.5576] (metrics.json › h4.bins.3) | 0.492 (metrics.json › h4.bins.3.genome) | 0.4108 (metrics.json › h4.bins.3.genome) | 0.5726 (metrics.json › h4.bins.3.genome) | 0.995 (metrics.json › h4.bins.3.genome) | 0.608 (metrics.json › h4.bins.3) | 0.4745 [0.3594, 0.5923] (metrics.json › h4.bins.3) | consistent (metrics.json › h4.bins.3) | consistent (metrics.json › h4.bins.3) | 0 [0, 0.02728] (metrics.json › h4.bins.3) | 1 [0.9442, 1] (metrics.json › h4.bins.3) |
| 4 (metrics.json › h4.bins.4) | 136 (metrics.json › h4.bins.4) | 0.5221 [0.4387, 0.6042] (metrics.json › h4.bins.4) | 0.529 (metrics.json › h4.bins.4.genome) | 0.4436 (metrics.json › h4.bins.4.genome) | 0.6181 (metrics.json › h4.bins.4.genome) | 0.995 (metrics.json › h4.bins.4.genome) | 0.657 (metrics.json › h4.bins.4) | 0.5221 [0.4039, 0.6378] (metrics.json › h4.bins.4) | consistent (metrics.json › h4.bins.4) | consistent (metrics.json › h4.bins.4) | 0.02206 [0.00753, 0.06285] (metrics.json › h4.bins.4) | 1 [0.9487, 1] (metrics.json › h4.bins.4) |
| 5 (metrics.json › h4.bins.5) | 136 (metrics.json › h4.bins.5) | 0.5074 [0.4243, 0.59] (metrics.json › h4.bins.5) | 0.522 (metrics.json › h4.bins.5.genome) | 0.4361 (metrics.json › h4.bins.5.genome) | 0.604 (metrics.json › h4.bins.5.genome) | 0.995 (metrics.json › h4.bins.5.genome) | 0.64 (metrics.json › h4.bins.5) | 0.5074 [0.39, 0.6239] (metrics.json › h4.bins.5) | consistent (metrics.json › h4.bins.5) | consistent (metrics.json › h4.bins.5) | 0.01471 [0.004042, 0.05203] (metrics.json › h4.bins.5) | 1 [0.9473, 1] (metrics.json › h4.bins.5) |
| 6 (metrics.json › h4.bins.6) | 137 (metrics.json › h4.bins.6) | 0.6204 [0.5369, 0.6974] (metrics.json › h4.bins.6) | 0.638 (metrics.json › h4.bins.6.genome) | 0.5577 (metrics.json › h4.bins.6.genome) | 0.7191 (metrics.json › h4.bins.6.genome) | 0.995 (metrics.json › h4.bins.6.genome) | 0.754 (metrics.json › h4.bins.6) | 0.6204 [0.5005, 0.7272] (metrics.json › h4.bins.6) | consistent (metrics.json › h4.bins.6) | consistent (metrics.json › h4.bins.6) | 0 [0, 0.02728] (metrics.json › h4.bins.6) | 1 [0.9568, 1] (metrics.json › h4.bins.6) |
| 7 (metrics.json › h4.bins.7) | 136 (metrics.json › h4.bins.7) | 0.6985 [0.6168, 0.7693] (metrics.json › h4.bins.7) | 0.729 (metrics.json › h4.bins.7.genome) | 0.6504 (metrics.json › h4.bins.7.genome) | 0.8057 (metrics.json › h4.bins.7.genome) | 0.995 (metrics.json › h4.bins.7.genome) | 0.839 (metrics.json › h4.bins.7) | 0.6985 [0.5797, 0.7956] (metrics.json › h4.bins.7) | consistent (metrics.json › h4.bins.7) | consistent (metrics.json › h4.bins.7) | 0.01471 [0.004042, 0.05203] (metrics.json › h4.bins.7) | 1 [0.9611, 1] (metrics.json › h4.bins.7) |
| 8 (metrics.json › h4.bins.8) | 136 (metrics.json › h4.bins.8) | 0.7059 [0.6245, 0.776] (metrics.json › h4.bins.8) | 0.745 (metrics.json › h4.bins.8.genome) | 0.6726 (metrics.json › h4.bins.8.genome) | 0.8133 (metrics.json › h4.bins.8.genome) | 0.995 (metrics.json › h4.bins.8.genome) | 0.843 (metrics.json › h4.bins.8) | 0.7059 [0.5874, 0.8018] (metrics.json › h4.bins.8) | consistent (metrics.json › h4.bins.8) | consistent (metrics.json › h4.bins.8) | 0.02206 [0.00753, 0.06285] (metrics.json › h4.bins.8) | 1 [0.9615, 1] (metrics.json › h4.bins.8) |
| 9 (metrics.json › h4.bins.9) | 137 (metrics.json › h4.bins.9) | 0.8759 [0.8103, 0.9211] (metrics.json › h4.bins.9) | 0.887 (metrics.json › h4.bins.9.genome) | 0.8292 (metrics.json › h4.bins.9.genome) | 0.9446 (metrics.json › h4.bins.9.genome) | 0.995 (metrics.json › h4.bins.9.genome) | 0.97 (metrics.json › h4.bins.9) | 0.8759 [0.7759, 0.935] (metrics.json › h4.bins.9) | consistent (metrics.json › h4.bins.9) | consistent (metrics.json › h4.bins.9) | 0.0365 [0.01569, 0.08259] (metrics.json › h4.bins.9) | 1 [0.969, 1] (metrics.json › h4.bins.9) |

Descriptive, both over the reportable bins: n-weighted mean genome precision 0.555 (metrics.json › h4); pooled A2-confirmed share 0.5381 [0.5116, 0.5644] (metrics.json › h4).
Searchable peptides: GIAP's digest (trypsin_p, at most 1 missed cleavage, 7-30 residues, no N-terminal Met excision; not DIA-NN's settings) of --strain-db, I folded to L (metrics.json › h4). A2's contaminant share of its accepted target rows (counted in c_value): 0.003435 [0.003073, 0.003839] (metrics.json › h4).

Caveats:

- The score bins span the de novo scores of (0.8, 1] only (Cascadia's reported range). (metrics.json › h4.caveats)
- c_value comes from A2's own exact matches and the confirmations from A2's identifications, so A2 enters both sides: a shared dependence, not an independent check. (metrics.json › h4.caveats)
- Class error rates of A3's ALT calls are unestimated (no spectrum pairing). (metrics.json › h4.caveats)
- Dependent rows: nested or overlapping de novo peptides of one locus share correctness, so the genome bin CI can be too narrow and a not_supported can come from that undercoverage; nested_share shows how clustered each bin is. (metrics.json › h4.caveats)
- Contaminant rows count in c_value (they never match the genome), which biases the genome precision up and lowers the power to refute; de novo contaminant peptides are not in A2's set and count as unconfirmed. (metrics.json › h4.caveats)
- A2's identifications are at 1 % FDR, so the lower bound can be off by up to that share (not corrected). (metrics.json › h4.caveats)
- A failure confined to a single score bin reads consistent_with_lower_bound: not_supported needs two consecutive reportable bins below the bound (the pre-registered fallback). (metrics.json › h4.caveats)

## 7. Not evaluated

One line per metric or hypothesis whose value is null with a reason: the pilot data do not allow it (§15.7 "reported where the pilot data allow").

- `m8`: pending (M10-07)
- `te1.share_without_q_at_cutoff`: not evaluable (no_crossings)

## 8. Runtimes

### Metric 7: runtime and peak memory (pilot measurements on a shared cluster, not the benchmark)

| | A1 | A2 |
|---|---|---|
| wall time | 1:36:48 (metrics.json › m7.A1) | 1:43:28 (metrics.json › m7.A2) |
| sum of step wall times | 7:10:39 (metrics.json › m7.A1) | 6:33:46 (metrics.json › m7.A2) |
| max step peak RSS (one task, not the arm's total memory) | 7.5 GiB (metrics.json › m7.A1) | 7.5 GiB (metrics.json › m7.A2) |
| complete (no cached step; else totals are lower bounds) | yes (metrics.json › m7.A1) | yes (metrics.json › m7.A2) |

### Benchmark (M10-12, `benchmarks/baseline.json`, wheat nodes)

Medians of the standard run against the estimates of §13.2 (35 min, range 16–55, and 25 min); a median at an estimate is reported, not judged.

| cell | median wall time | commit | CPU |
|---|---|---|---|
| K-12, standard | pending (M10-12) | pending (M10-12) | pending (M10-12) |
| K-12, core | pending (M10-12) | pending (M10-12) | pending (M10-12) |
| simulated, standard | pending (M10-12) | pending (M10-12) | pending (M10-12) |
| simulated, core | pending (M10-12) | pending (M10-12) | pending (M10-12) |

Benchmark definition (ADR-M10-D): pending (M10-12), status pending (M10-12).

## 9. Reproducibility

Freeze commit of the runbook: pending (render --final --runbook) (runbook status pending (render --final --runbook)); every arm, the compare and the metrics ran at it (`validation/pilot/reproduce.md`, Evidence table).
The runbook's check writes `verify-<commit>.json` beside the freeze runs; it is part of the review packet.
TE7 at the freeze commit: pending (M10-06) (metrics.json › te7.giap).
The benchmark was not re-run at the freeze commit (M10-12 L12); its commit is in section 8.
Row-level outputs stay on FarmShare and are removed only by the runbook's guarded clean-up step.

## 10. Deviations from the spec

- DIA-NN 2.3.1 Academia instead of 2.7.0 (ADR-0051).
- Row-level inputs and outputs are not committed (Khoa's rule of 2026-10-23); aggregates and hashes are (2026-10-24).
- Protein groups are groups with at least one accepted precursor, with no protein-level FDR control (metric 1).
- Strain-specific genes are defined by peptide identity, not orthology (metric 3).
- Region tests are not compared, and the Jaccards are not restricted to genes detected in
  both arms (H3; a conservative stand-in, biased against H3).
- H2 is deferred to the full validation; metric 8 stands beside it without a verdict.
- The A3 branch when it is the fallback or absent, every deferral made at the freeze, and any
  accepted ADR that changes a pilot setting are stated in the notes.

## Notes

No notes.
