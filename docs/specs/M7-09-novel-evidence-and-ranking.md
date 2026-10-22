# M7-09: Class evidence rules, RBS support, promoter/terminator evidence and ranking

Milestone M7 · track N (novel ORFs) · **pass 1 (after the wiring)** · depends on: M7-06, M7-07,
M7-08 (fills columns of the already-wired `novel_orfs` step); merged M2-13 promoter tracks, M2-15
terminators · estimate ≈ 300 code lines (+ ≈ 450 test lines) · not statistical (the class FDR is
M7-11), but the quantile threshold is a data-derived cut-off — the statistician reviews decision 3

## Spec sections implemented

- §7.6.3 in full — the class table:

  | Class | Where the peptide sits | Extra evidence before calling |
  |---|---|---|
  | `nterm_extension` | in frame, upstream of the annotated start, after the previous in-frame stop | an upstream in-frame start codon and RBS support |
  | `intergenic` | no annotated CDS on either strand | start codon, ≥ 10 codons, RBS support; a homolog raises confidence later |
  | `cds_outframe` | same strand, other frame, inside a CDS | ≥ 2 distinct frame-specific peptides |
  | `antisense` | opposite strand, inside a CDS | ≥ 2 distinct peptides |
  | `pseudogene` | on a Bakta pseudogene | the frameshift or stop junction is checked: `junction_spanned` … and `homopolymer_junction` (a homopolymer ≥ 5 nt at the junction). These are labels only |
  | `sorf` | a Bakta sORF | detection evidence only |

  "**RBS support** means the OSTIR initiation score at the chosen start is ≥ the
  `novel_orfs.rbs_min_quantile` (default 0.10) quantile of scores at annotated starts of detected
  genes."
- ADR-0013 row for `cterm_extension`: "the in-frame continuation to the next stop and a peptide
  spanning the annotated stop; labelled `readthrough_candidate`".
- ADR-0012 decisions 1–5 (terminator 0–300 nt downstream on the same strand, best confidence;
  `promoter_supported` = `promoter_score ≥ promotech_yes`; the four columns and two flags;
  ranking key "terminator confidence, descending" after the promoter; "promoter/terminator evidence
  never gates a call").
- §7.6.5 (G13) ranking key; §7.6.4 "Neither: a call needs ≥ 2 distinct peptides and carries
  `class_fdr = unestimated`" (still the only branch until M7-11); config `novel_orfs.min_codons`,
  `novel_orfs.rbs_min_quantile`, `annotation.promoters.promotech_yes`,
  `annotation.terminators.terminatornet_yes`.

## Public interface

`src/giap/novel/evidence.py` (pure)
```python
def rbs_threshold(gene_table: pd.DataFrame, q: float, min_n: int = 20) -> float | None: ...
def terminator_downstream(orf: Orf, terminators: Sequence[Terminator],
                          rep: Replicon) -> tuple[float | None, int | None]: ...  # conf, distance
def pseudogene_junction(gene: Gene, genome: str, rep: Replicon) -> int | None: ...  # codon start
def class_evidence(row: Mapping[str, object], cfg: GiapConfig, rbs_min: float | None) -> tuple[bool, list[str]]: ...
def rank_key(row: Mapping[str, object]) -> tuple[object, ...]: ...
def apply_evidence(table: pd.DataFrame, ...) -> pd.DataFrame: ...   # fills rbs_supported, flags, called, rank
```

## Rules

1. **Chosen start** = `start_best_ostir`, else `start_longest` (M7-07 decision 2).
2. **RBS threshold.** The `rbs_min_quantile` quantile (`numpy.quantile`, linear interpolation) of
   `tir_ostir` over eligible (non-partial) gene-table rows with `detected = true` and `tir_ostir`
   not NA. Fewer than 20 such genes → threshold `None`: `rbs_supported` NA, flag
   `rbs_threshold_unestimated` on every ORF that needs RBS support, and those ORFs are **not
   called** (decision 1).
3. **`rbs_supported`** = `ostir_tir ≥ threshold` at the chosen start; NA when either is NA; flag
   `rbs_supported` when true (§9.4).
4. **Class rules** (each in addition to the error-control rule 8):
   - `nterm_extension`: `start_longest` exists and lies upstream of the annotated start (in
     transcription direction) **and** `rbs_supported`.
   - `intergenic`: a start codon exists, codons from the chosen start to the stop (stop excluded) ≥
     `novel_orfs.min_codons`, **and** `rbs_supported`.
   - `cds_outframe`: ≥ 2 distinct peptides (all its peptides are frame-specific: in-frame ones were
     excluded in M7-06).
   - `antisense`: ≥ 2 distinct peptides.
   - `pseudogene`: no extra rule; labels `junction_spanned` (some peptide's interval contains the
     junction codon — decision 2) and `homopolymer_junction` (a run of ≥ 5 identical nucleotides
     overlapping [junction − 5, junction + 8)).
   - `sorf`: no extra rule.
   - `cterm_extension`: an in-frame continuation to the next stop exists (M7-06 always gives one or
     `end_3p_open`) and a peptide spans the annotated stop codon; label `readthrough_candidate`.
     With the merged matcher no peptide can span a stop (M7-05 rule 6), so this class is reachable
     only through frameshift-split peptides and is not called until M7-16 (decision 4).
5. **Junction.** For a pseudogene, the first premature in-frame stop codon reading from the
   annotated start (both Bakta causes, internal stop and frameshift, end the start frame early);
   none found → no junction labels, flag `junction_unknown` (decision 2).
6. **Promoter (ADR-0012 d2).** `promoter_score` (M7-08) at the chosen start; `promoter_supported` =
   score ≥ `annotation.promoters.promotech_yes`; flag `promoter_supported`.
7. **Terminator (ADR-0012 d1).** Over the run's terminator layer (TerminatorNet, or TransTermHP
   under `--terminators transtermhp`), calls on the ORF's strand whose 5′ end (in transcription
   direction) lies in [stop codon end, stop codon end + 300] nt, wrapping on circular replicons:
   `terminator_conf` = the best confidence, `terminator_distance_nt` = its 5′ end − stop end,
   `terminator_downstream` = conf ≥ `annotation.terminators.terminatornet_yes`; flag
   `terminator_supported`. No call → conf and distance NA, `terminator_downstream` false. Layer
   missing → all three NA, run flag `layer_missing:terminators`.
8. **Error control (unchanged from M7-06 rule 10 until M7-11):** distinct peptides ≥
   `novel_orfs.min_peptides_without_decoys`; `class_fdr = unestimated`.
9. **`called`** = rule 8 **and** the class rule (rule 4). Promoter and terminator never gate.
10. **Rank** (G13 + ADR-0012 d4), over **called** ORFs only, 1 = best (decision 5): distinct
    peptides ↓; `best_peptide_q` ↑ (NA last); start codon of the chosen start ATG < GTG < TTG <
    other (none last); `ostir_tir` ↓; `promoter_score` ↓; `terminator_conf` ↓; `length_aa_best` ↓;
    then `orf_id` ↑ (deterministic tie-break). NA sorts last in every key. Uncalled rows: `rank`
    NA.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| < 20 detected genes with an OSTIR score | rule 2: `nterm_extension`/`intergenic` not called, flagged |
| `rbs_min_quantile` = 0 | threshold = the minimum score; every scored start passes |
| OSTIR missing at the chosen start | `rbs_supported` NA → RBS-requiring classes not called |
| intergenic ORF of 9 codons | not called (`min_codons` 10) |
| `nterm_extension` whose `start_longest` is the annotated start itself | not called (no upstream start) |
| terminator exactly 300 nt after the stop | found (closed interval) |
| terminator 301 nt after, or on the other strand | not found |
| terminator across the origin | found when within 300 nt along the wrapped strand |
| linear contig, window running off the end | clipped; no wrap |
| pseudogene with no premature stop | flag `junction_unknown`, no junction labels |
| all keys tie | `orf_id` decides |

## Examples

Gene-table `tir_ostir` of the 30 detected genes: 100, 200, …, 3000 → 10 % quantile (linear) =
100 + 0.1 · 29 · 100 = **390**.

| # | ORF | Expected |
|---|---|---|
| 1 | intergenic, 20 codons, chosen-start TIR 500, 2 peptides | `rbs_supported` true, **called** |
| 2 | as 1, TIR 300 | `rbs_supported` false, not called |
| 3 | as 1, 9 codons | not called |
| 4 | intergenic, 1 peptide, TIR 5000 | not called (rule 8) |
| 5 | `cds_outframe`, 2 distinct peptides, no OSTIR | called (no RBS rule) |
| 6 | `antisense`, 1 peptide | not called |
| 7 | `nterm_extension` (M7-06 ex. 6), TTG start upstream, TIR 800, 2 peptides | called |
| 8 | pseudogene with an internal TGA at [600, 603), a peptide over [590, 620), `AAAAA` at [595, 600) | labels `junction_spanned;homopolymer_junction` (hand-built interval: the label function only looks at coordinates) |
| 9 | `+` ORF stop ends at 1-based 123; terminators on `+` starting at 1-based 200 (conf 0.8) and 450 (0.95), one on `−` at 150 (0.99) | `terminator_conf` **0.8**, `terminator_distance_nt` **77**, `terminator_downstream` true (≥ 0.5) |
| 10 | only 15 detected genes with TIR | intergenic ORF of example 1: `rbs_supported` NA, flag `rbs_threshold_unestimated`, not called |
| 11 | ranking: A (3 pep, q 0.001, ATG), B (3 pep, q 0.001, GTG), C (2 pep, q 1e-5, ATG) | ranks A 1, B 2, C 3 (peptides first, then q, then start codon) |
| 12 | ranking tie on every key | lower `orf_id` first |

## Acceptance criteria

- [ ] `tests/unit/novel/test_evidence.py`: examples 1–12 and every degenerate row; §7.6.6 "each
  class", "the ≥ 2-peptide rule", "ranking order"; ADR-0012 consequences ("a terminator 0–300 nt
  downstream on the same strand is found, one on the other strand or beyond 300 nt is not, origin
  wrap, NA when no call").
- [ ] Property: `called` never depends on `promoter_score` or `terminator_conf` (permute them).
- [ ] Property: `rank` is a strict total order on called ORFs and invariant to input row order.
- [ ] The wired step now writes `rbs_supported`, the ADR-0012 columns, `called` and `rank`; the
  M7-08 executor-equality and TE7 tests still pass.

## Calibration and attacks owed

- The novel-ORF F1 gate (M7-13) runs on these rules.
- **Statistician:** the RBS threshold is estimated from the same run's detected genes; detected
  genes are enriched for strong RBSs, so the 10 % quantile is higher than the quantile over all
  genes and the rule is conservative. M7-13 records the share of planted true starts that pass.

## Questions for the lead

1. Minimum number of detected genes for the RBS quantile. **Proposed default:** 20 (the M6 floor
   for a fitted distribution); below it the RBS-requiring classes are not called and flagged, never
   called on an undefined threshold.
2. Pseudogene junction: Bakta's JSON carries no junction position. **Proposed default:** the first
   premature in-frame stop from the annotated start (rule 5); `junction_spanned` = a peptide
   interval contains that codon; homopolymer window [junction − 5, junction + 8). Labels only, so a
   simple, documented rule beats a reconstruction of Bakta's internals.
3. "≥ 2 distinct frame-specific peptides" for `cds_outframe` duplicates the §7.6.4 ≥ 2-peptide rule
   without decoys; with decoys (M7-11), is the class rule still ≥ 2? **Proposed default:** yes —
   the table's rule is class evidence, independent of the error-control branch.
4. `cterm_extension` calls before M7-16. **Proposed default:** report the rows with
   `readthrough_candidate` when a frameshift-split peptide crosses the annotated stop, but never
   `called = true` until M7-16 decides the stop-codon matching question.
5. Rank over called ORFs only (proposed) or all rows? **Proposed default:** called only; uncalled
   rows carry `rank` NA, so "top 20" (V13) means top 20 calls.

## Lead decisions (2026-10-22)

1. **Questions 1–5: the proposed defaults are accepted.**
   - **RBS quantile:** it needs at least 20 detected genes. Below that, the RBS-requiring classes are not called and the rows are flagged.
   - **Pseudogene junction:** the junction is the first premature in-frame stop from the annotated start. `junction_spanned` means that codon lies inside a peptide interval. The homopolymer window is [junction − 5, junction + 8). These fields are labels only.
   - **`cds_outframe`:** this class still needs at least 2 distinct frame-specific peptides when decoys are present (M7-11).
   - **`cterm_extension`:** these rows carry `readthrough_candidate` and are never `called = true` until M7-16.
   - **Rank:** only called ORFs are ranked; uncalled rows get rank NA.
2. **Branch:** m7/novel-evidence is stacked on m7/novel-steps (#168).

3. **Test contract accepted (TA 3ab4a3e, e295140, 774946d).** The test author resolved nine ambiguities, listed in the test docstrings, and they are accepted as written. In particular:
   - `class_evidence` excludes rule 8.
   - The fixed ≥ 2 peptides rule applies to `cds_outframe`/`antisense`.
   - Terminator distance = the terminator's 5′ base − the stop's last base, in a closed window 0..300.
   - Flags follow §9.4 vocabulary order, not alphabetical order.
4. **Flag order:** `apply_evidence` (or `orf_table`) emits flags in the novel_orfs schema vocabulary order. This replaces M7-07's "flags sorted"; the lead updates the stale `test_outputs.py` docstring.
5. **Peptide intervals for the junction and cterm labels** reach the step as additive columns of `novel_candidates.parquet`. The step's inputs are unchanged, so the M7-08 registration test still holds. If the frozen schema of the parquet needs it, the change is additive only.
