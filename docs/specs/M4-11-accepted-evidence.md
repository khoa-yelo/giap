# M4-11: Accepted evidence and per-gene evidence counts

Milestone M4 · track G · depends on: M4-06 · estimate ≈ 300 code lines (+ ≈ 350 test lines) ·
**can start now** (fixtures from M4-09, or hand-built C9/C11 frames)

## Spec sections implemented

- §7.3.1 quoted in full: "A peptide is **accepted evidence** for a gene when all of these hold: it
  passes the engine threshold (`q_value ≤ input.q_value_max` when a q-value exists); its placement is
  accepted (`placement_q ≤ 0.01`, or `same_db`); its locus class assigns it to the gene (§7.2.6); its
  uniqueness is `unique` or `multi_identical`; under `genome_relation = same`, it is an exact hit.
  **Distinct peptides** are distinct stripped sequences. Charge states, runs and modified forms of one
  sequence count once."
- §7.2.6: "`multi_identical` … counts toward detection of each copy, flagged, never toward abundance";
  "`shared`: different genes; counts toward neither detection nor abundance of a single gene";
  "`nterm_extension` counts toward the gene's detection (flag `nterm_extension_evidence`) but not its
  abundance".
- §7.2.7: "`related` (default): … Accepted tolerant hits count toward detection and abundance and are
  labelled `strain_difference`." "`same`: … Tolerant hits are reported but only exact hits count."
- §9.3 `genes.tsv` evidence columns: `n_peptides_unique`, `n_peptides_multi_identical`,
  `n_peptides_shared`, `n_peptides_tolerant`, `coverage_frac`, `best_peptide_q`, `detected_2pep`.
- §9.4 gene flags: `multi_identical`, `shared_evidence`, `tolerant_evidence`, `strain_difference`,
  `nterm_extension_evidence`.

## Public interface

`src/giap/genetable/evidence.py` (pure)
```python
ASSIGNING_CLASSES = ("cds_inframe", "sorf", "nterm_extension", "pseudogene")   # see M4-06 Q3

def accepted_peptides(peptide_hits: pd.DataFrame, peptides: pd.DataFrame | None,
                      cfg: GiapConfig) -> pd.DataFrame: ...
    # one row per (gene, peptide) accepted pair: gene, peptide, q_value, pep, score, intensity,
    # uniqueness, route, locus_class, abundance_eligible (bool), start, end, contig, crosses_origin
def gene_evidence(genes: pd.DataFrame, peptide_hits: pd.DataFrame, accepted: pd.DataFrame,
                  replicons: pd.DataFrame, min_peptides: int = 2) -> pd.DataFrame: ...
    # one row per C2 gene: the §9.3 evidence columns above + n_accepted, evidence_flags (list)
```
- Engine threshold: q from C9 joined on `peptide` (C11 `best_q_input` must equal it; C9 wins if both).
  NA q passes (the threshold applies "when a q-value exists").
- Placement accepted: `placement_q ≤ stats.alpha.placement` or `"same_db" in flags`.
- `abundance_eligible` = `uniqueness == "unique"` and `locus_class != "nterm_extension"` and
  (route exact, or route tolerant under `related`).
- Counts are of **distinct peptides** per gene: `n_peptides_unique` (accepted, unique),
  `n_peptides_multi_identical` (accepted, multi_identical), `n_peptides_shared` (C11 rows with
  `uniqueness = shared` on this gene; not accepted), `n_peptides_tolerant` (accepted, route tolerant).
- `n_accepted` = distinct accepted peptides; `detected_2pep` = `n_accepted ≥ min_peptides`
  (`detection.min_peptides_sequences_only`) for eligible genes, false otherwise.
- `best_peptide_q` = min C9 `q_value` over accepted peptides (NA if none has a q).
- `coverage_frac` = |union of accepted placement intervals ∩ gene interval| / gene length in bp
  (origin-crossing intervals split first, M1-01 `split_at_origin`); 0 when no accepted peptide.
- Flags: `multi_identical` (n_multi_identical > 0), `shared_evidence` (n_shared > 0),
  `tolerant_evidence` and `strain_difference` (an accepted tolerant hit under `related`),
  `nterm_extension_evidence` (an accepted `nterm_extension` peptide).

## Examples

Config: `q_value_max = 0.01`, `alpha.placement = 0.01`, `genome_relation = related`.

| # | C11 / C9 rows for gene G | Expected for G |
|---|---|---|
| 1 | `PEPTIDEK` (q 0.001, cds_inframe, unique, exact, placement_q 0.0) | n_unique 1, n_accepted 1, `detected_2pep` false, best_peptide_q 0.001 |
| 2 | #1 + `ANOTHERPEPK` (q 0.004, same) | n_unique 2, `detected_2pep` true, best_peptide_q 0.001 |
| 3 | #1 with C9 q = 0.02 | not accepted → n_accepted 0, best_peptide_q NA |
| 4 | #1 with `placement_q = 0.03`, flags `[same_db]` | accepted (same_db) |
| 5 | `PEPTIDEK` charge 2 and 3, two runs (one C9 row after aggregation) | counted once |
| 6 | peptide on G and its IS copy G′ (`multi_identical`) | G and G′ each: n_multi_identical 1, n_accepted 1, flag `multi_identical`; `abundance_eligible` false |
| 7 | peptide on G and unrelated H (`shared`) | G: n_shared 1, n_accepted 0, flag `shared_evidence` |
| 8 | tolerant hit (`route tolerant`, n_diff 1, placement_q 0.004) | `related`: accepted, n_tolerant 1, flags `tolerant_evidence;strain_difference`; `same`: not accepted, n_tolerant 0 |
| 9 | `nterm_extension` peptide | accepted (counts for detection), flag `nterm_extension_evidence`, `abundance_eligible` false |
| 10 | gene 300 bp, accepted peptides at [0, 30) and [15, 60) bp of the gene | coverage_frac = 60 / 300 = 0.2 |

## Edge cases

- No C11 / C9 (proteins only) → every gene gets zeros, NA best q, `detected_2pep` false.
- Partial genes get their counts but `detected_2pep` is false (not eligible).
- A peptide with C11 rows on two genes of which one is `cds_outframe`: only the assigning row counts
  (M3 already sets `gene` NA on non-assigning rows; the class check is a second guard).
- An origin-crossing gene with a peptide spanning the origin: coverage computed on split parts.
- C9 rows absent for a C11 peptide (contract violation) → `ValueError` naming the peptide.

## Acceptance criteria

- [ ] `tests/unit/genetable/test_evidence.py`: examples 1–10 as parametrised cases.
- [ ] Property: relabelling peptides (string permutation) leaves all counts unchanged; adding a
  duplicate C11 row for the same peptide and locus changes nothing; `n_accepted ≥ n_unique`.
- [ ] On the M4-09 default fixture: every `detected_true` gene has `n_accepted ≥ 1` (placement is
  perfect there); runs in < 1 s for 4,500 genes.

## Methods note (statistician sign-off)

`docs/methods/detection-tiers.md` §1 "Accepted evidence" (the five conditions, distinct peptides,
uniqueness handling).

## Questions for the lead

1. Peptides with no q-value but a PEP: apply a PEP threshold? **Proposed default:** no; §7.3.1 names
   only q. PEP enters through the tiers (M4-12).
2. Should `coverage_frac` use only accepted peptides (as above) or all placed peptides? **Proposed
   default:** accepted only, consistent with "evidence".

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted: pseudogene and `cterm_extension` placements count toward
   detection; NA `placement_q` without `same_db` is rejected; C9 decoy rows yield nothing; C9 wins over C11
   on q without raising; flags compared as sets; one row per gene in C2 order.
2. `cterm_extension` placements are not `abundance_eligible` (abundance from the canonical CDS only).
3. Contaminant-matching peptides that pass the five conditions are accepted; the gene row carries a separate
   bool column `contaminant_evidence` (accepted rows carry `is_contaminant`), not an `evidence_flags` entry
   (visible to M4-12/M4-14); M4-12 decides whether they count for tiers.
4. Several accepted placements of one peptide in one gene keep one row (exact first, then leftmost);
   coverage uses that row. Without C9, q/pep/score come from C11 (`best_q_input`, `pep_input`,
   `score_input`) and intensity is NA.
