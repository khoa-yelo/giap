# M7-17 G4-k3 flag gate, full tier: statistician's verdict (PR #291, ebfdd065)

**Verdict (b): a test-and-generator defect, not a method failure.** The flag follows its rule on
every reached row. The two unflagged "intergenic" rows are not false intergenic ORFs. They are
**planted true antisense ORFs** whose stop the injector reused. The pooled 300 records can be
rescored. No rerun is needed.

## What failed
`_check(True)` failed at its last line, `per_class["intergenic"]["flagged"] == ["reached"]`,
with 467 of 469. `summarise` bins rows by the **injected** class. Rule 4, and the test's own
violation check, use the row's **TSV** (observed) class.

## Class confusion in the 300 records (Slurm 1881039 for replicates 0–187, 1886178–87 for 188–299)
- Rows at injected keys, by injected → observed class:
  - intergenic → intergenic: 5,972;
  - intergenic → antisense: 20;
  - intergenic → cds_outframe: 6;
  - cds_outframe and antisense: 6,000 each, all matched.
- Of the 469 reached intergenic-injected rows, 443 are observed as `intergenic`, and all 443 are
  flagged. The other 26 are observed as antisense (20) or cds_outframe (6). 24 of them have no
  RBS support and are flagged. 2 have RBS support, so they are exempt (replicate 13,
  NOVEL_chr1_165770_m; replicate 253, NOVEL_chr1_791521_m).
- Violations of the rule, both ways: **0 of 18,448** reached rows.

## Cause (diagnostics: Slurm 1886794, 1886826, 1886919; census: 1886842)
- **Stop collision with a planted ORF: 25 of 26 rows, including both exempt rows.**
  - `inject_false_orfs` takes stop-closed segments from the non-coding pool. The pool masks the
    codons of planted ORFs, but a segment upstream of a planted ORF's start, in its frame, ends
    at that ORF's stop. Its key is then the planted ORF's key.
  - Example, replicate 13: the row is truth NOVEL_0578 (antisense, 165770–165982 −). It has 8
    peptides, all of them background reads of the true ORF.
  - Example, replicate 253: the row is truth NOVEL_0649 (antisense, stop 791521). The injected
    exact read LANMFWYK is placed `intergenic` at 791803–791826, upstream of the true start, and
    joins that ORF.
  - Census over all 300 replicates (simulation and injection only): 97 of 18,000 injected keys
    are planted truth keys. By injected class: intergenic 50 (25 at a true intergenic key, 19 at
    a true antisense key, 6 at a true cds_outframe key), cds_outframe 11, antisense 36. The
    per-replicate counts explain all 6 cds_outframe rows and 19 of the 20 antisense rows.
  - `null_segments` in the same harness already excludes truth keys and planted spans;
    `inject_false_orfs` does not.
- **Background wrong read joins the ORF: 1 of 26** (replicate 130, NOVEL_chr1_567837_p). A
  chance de novo read (PSSSLTK, score −4.05) is placed `antisense` over GSIM_02660 and shares the
  injected segment's stop. The row has no RBS support and is flagged.
- **Reclassification is correct behaviour.** An ORF is every placement to one stop (§7.6.2,
  "Grouping"). Its class is the highest-priority member class (§7.2.6; `CLASS_PRIORITY`), and
  antisense and cds_outframe rank above intergenic. GIAP cannot see what was injected.

## Q1: injected or observed class?
- The spec defines the observed class:
  - rule 4 reads "the class";
  - the edge-case row `class_mismatch` says the flag reads "its TSV class";
  - criterion 6 exempts "a `cds_outframe` or `antisense` **row**".
- "For `intergenic` this is every such row" follows from these only if the injected class equals
  the observed class. A stop collision breaks that assumption.
- The test's docstring also uses the TSV class ("unless its TSV class is …"). Only the last
  assertion bins by injected class.

## Fix that keeps the gate's intent
1. **Required.** The last assertion counts rows whose **observed** class is `intergenic`, and
   requires all of them to be flagged.
   - Why this is not a widening:
     - every reached row is still checked against rule 4 both ways, with zero tolerance;
     - R, q ≤ 0.05, the RBS threshold and the rule are unchanged;
     - the intergenic statement still has to hold exactly on every row GIAP calls `intergenic`.
   - What it gives up: exempting an injected-intergenic row that GIAP calls antisense with chance
     RBS support. That is the antisense chance-RBS escape the methods note already documents
     ("`intergenic` and `nterm_extension` are never exempt" refers to the TSV class), and it
     stays recorded.
   - The diagnostics should also record the injected → observed table above.
2. **Backlog; needs a rerun, not this PR.** `inject_false_orfs` should skip segments at a
   planted truth key, as `null_segments` does, so that every injected key is a false ORF.
   - This makes G4-k3 stricter, not looser: today 0.54 % of injections land on true ORFs and add
     no false ORF. It does not change the ADR-0050 conclusions.
   - It changes the G4-k3 stream (rng draws), so the lead should list it for Khoa's batch.
   - Do **not** exclude every segment whose span overlaps a CDS. That would narrow the
     intergenic attack territory, and the replicate-130 path is correct GIAP behaviour.

## Rescored, observed class, R = 300 (from the existing records; no rerun)

| Class | Reached | Flagged | Exempt | Escape share |
|---|---|---|---|---|
| intergenic | 443 | 443 | 0 | 0 |
| cds_outframe | 5,993 | 5,667 | 326 | 0.054 |
| antisense | 6,012 | 5,637 | 375 | 0.062 |

- With fix 1, the gate passes on the pooled records. The records store `class`, so rescoring is
  exact.
- The escape shares stay below the 0.2 line of lead decision 1.
- Under the injected class, the recorded escape shares are: intergenic 2/469 = 0.004,
  cds_outframe 326/5,987 = 0.054, antisense 373/5,992 = 0.062.
- Scripts and logs: `local/stat/m717-gate/` (`diag2.py`, `census.py`, `census_*.jsonl`, `*.out`).
