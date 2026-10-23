# ADR-0049: M7 variant-call precision under de novo errors (G1, A6a, A6b, A6c)
- Status: proposed
- Date: 2026-10-23
- Spec sections: §7.5 (call rule, flags), §8.4 (`f1_min.variants_high` 0.90), §8.5 attack 6, §13 gate M7
- Needs Khoa: yes. The fix changes the variant call rule, and the gate threshold must not be widened (R3).

## Context
The M7-13a test-author's reference run (Slurm 1872235) shows the M7 variants gate failing on the
current call rule:
- **G1:** F1 is 0.838 at 100k with R = 20 (P 0.776, Rc 0.910), and 0.852 at 1M over 10
  replicates. Both are below `f1_min` 0.90, and precision is the cause. The
  non-`mass_equivalent` subset has P 0.951 and Rc 0.675.
- **A6a:** 0–3 isobaric non-swap ALTs per replicate are not flagged `mass_equivalent`. An
  N-terminal block read such as AG→Q aligns as the substitution G>Q.
- **A6b:** at random adjacent swaps the share of replicates with a `high` call is 0.95. Two
  overlapping true peptides that carry the same swap give two distinct peptides.
- **A6c:** systematic swaps are expected to fail (M7-13 risk 1).
- **Rows outside G1 (statistician, #255).** Besides the 669 gated substitution predictions, GIAP
  emits 277 `high` rows that are not substitutions (block 218, swap 49, frameshift 10). §7.5.6
  speaks of "high ALT calls", so if these are false, precision over all `high` rows is 0.549
  (F1 0.685). Whether G1 covers them (M7-13 Q2 scopes G1 to substitutions) is part of the choice.

## Options (for Khoa)
1. **Call rule:** `mass_equivalent` substitutions and swap-type ALTs never reach `high` without
   support from a peptide that is not mass-equivalent (or a reference spectrum).
2. **Alignment:** treat N-terminal two-residue block reads as blocks, and flag them
   `mass_equivalent`.
3. **Distinct-peptide support:** require that the distinct supporting peptides disagree at more
   than the swapped pair (for example, different cleavage ends).
4. **Gate scope:** gate G1 on the non-`mass_equivalent` set. This changes what §8.4 measures.
   On current evidence it also fails: F1 0.790 (P 0.951, Rc 0.675). Option 1 is likely to cost
   recall the same way, since it removes the same calls from `high`.

## Decision
None yet. Until Khoa decides, M7-13a commits G1, A6a, A6b and A6c as strict xfails that name
this ADR. Each records its estimate in the nightly records, so a fix flips it to XPASS. Gate M7
is not passed.
