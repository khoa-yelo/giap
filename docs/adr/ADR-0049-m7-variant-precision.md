# ADR-0049: M7 variant-call precision under de novo errors (G1, A6a, A6b, A6c)
- Status: accepted-by-khoa (2026-10-23)
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

## Decision (Khoa, 2026-10-23)
1. **Mass-equivalent calls stay in `high`, flagged.** A substitution, block or swap whose ALT is
   mass-equivalent to the reference keeps its tier if it meets the current call threshold, and
   carries the `mass_equivalent` flag. The call rule's thresholds are unchanged.
2. **Flag coverage (option 2).** N-terminal two-residue block reads such as AG→Q are treated as
   blocks and flagged `mass_equivalent`. Swap-type ALTs are flagged too.
3. **Swap support (option 3 not taken).** Overlapping peptides that carry the same swap still count
   as distinct support. Swaps are flagged, not demoted.
4. **G1 scope.** G1 (`f1_min.variants_high` 0.90, unchanged) scores unflagged `high` calls against
   *identifiable* truth: truth variants that are mass-equivalent to the reference leave the recall
   denominator too. Flagged calls are reported, not scored. Before this applies, the statistician
   recomputes G1 under this scope (before and after the flag changes) and records it here; the
   tolerance is never widened (R3).
5. **Attacks.** A6a: every isobaric non-swap ALT is flagged `mass_equivalent`. A6b: every random
   adjacent swap call is flagged (replaces "no `high` call"). A6c: unchanged (expected fail).
6. **Rows outside substitutions** (block, swap, frameshift `high` rows): out of G1 (M7-13 Q2);
   their precision is recorded per type. Open: whether they get their own gate (backlog).

## Evidence (statistician, Slurm 1874972; local/stat/adr0049/)
| scope | 100k b1 | 100k b2 | 1M |
|---|---|---|---|
| (a) all high vs all truth | 0.848 | 0.838 | 0.852 |
| (b) unflagged vs all truth | 0.786 | 0.790 | 0.804 |
| (c) unflagged vs identifiable truth | 0.933 [0.916, 0.950] | 0.927 | 0.932 [0.924, 0.941] |
| (d) (c) + decided flags | 0.934 [0.917, 0.950] | 0.927 [0.905, 0.949] | 0.933 |
- Identifiable-truth rule: a truth row leaves the recall denominator when
  `giap.variants.flags.mass_equivalent(ref_aa, alt_aa, variants.mass_equivalent_tol_da)` is true
  (the same test, including the B.5 modification deltas, that sets the call flag). It removes 26 %
  of recall-eligible truth at 100k (800 → 593) and 24 % at 1M, mostly residue-plus-modification
  matches (V>L, N>K, D>E, S>T, G>A …), not plain isobaric residues.
- 37 % of `high` substitution calls are flagged and reported but not scored.
- Residual false positives are mostly unannotated start-codon M calls (ALT M at the peptide
  N-terminus over ref L) and reversed-stretch decompositions; residual false negatives are mostly
  `candidate`-tier rows (call threshold, not flags).
## Implementation rules (from the evidence)
- A6a: flag a substitution when an ALT peptide also aligns as the 2→1 block
  `ref[p-1:p+1] → alt` with its prefix shifted by one (q = 0, or `P[:q] == ref[p-1-q:p-1]`).
  This fixes 24 of 24 unflagged isobaric rows. No right-hand or looser form (they flag true
  positives by chance).
- A6b: `flags_for` (src/giap/variants/flags.py) flags every `swap` row `mass_equivalent`
  (75 of 75).

Open, pending Khoa (reversible default): "identifiable" uses the modification-inclusive
`mass_equivalent` test above; a strictly isobaric-residue rule is the alternative, to be recomputed
on request.

