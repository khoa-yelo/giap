# M7-16: Peptides that span a stop codon (Sec, read-through, `cterm_extension`) — decision-gated

Milestone M7 · track R (mapping) · **pass 2, gated on a lead/Khoa decision; not started before it**
· depends on: merged M3-06/M3-08 tolerant DP, M3-10 placement FDR, M3-14 loci, M3-23 makedb;
M7-05, M7-07, M7-09 · estimate ≈ 220 code lines (+ ≈ 400 test lines) · **stats** (changes what
the placement FDR searches)

## The gap

Three accepted parts of the design assume that a peptide can be placed across a genome stop codon:

- §7.5.3 G8 `stop_codon`: "The REF residue is a stop (e.g., Sec read through UGA)", and §7.5.1 "Sec
  sites … appear as ALT calls" under `search_db = same`;
- ADR-0013: `cterm_extension` evidence "a peptide spanning the annotated stop"; the
  `readthrough_candidate` label; the patch row with a new end;
- §11.6 TE2's exception list (`stop_codon` sites).

The merged matcher cannot do it. The six-frame index writes `*` at stops (§7.2.1), the exact
automaton never matches across one, and the M3-06 DP states "A target `*` (`STOP`) … is never
consumed". M3-22 found exactly this: the 4 TE1 misses are selenocysteine (fdnG) peptides. So
`stop_codon` and `cterm_extension` calls are unreachable in v1 as built (M7-05 rule 6, M7-09
rule 4).

## Options (for the decision; the spec-writer does not choose)

| Option | What changes | Cost | Risk |
|---|---|---|---|
| A. **Leave it** | nothing; document `stop_codon`/`cterm_extension` as reachable only through frameshift splits; TE2's `stop_codon` clause is vacuous | 0 | the G8 flag and ADR-0013 stay partly dormant |
| B. **Annotated-stop substitution** | the tolerant DP may consume a target `*` **only at the annotated stop codon of a CDS/sORF**, at cost 1, recorded as `substitution` with REF `*`; decoys get the same rule at the mirrored positions | ≈ 220 lines + re-calibration of M3-21 placement suites | small: one position per gene; decoys keep symmetry |
| C. **Any-stop substitution** | `*` consumable anywhere at cost 1 | ≈ 150 lines | large: every stop becomes a chance-hit site; the M3 chance-hit table (B.2) and placement FDR need re-derivation |

**Proposed default: B**, as an ADR (it changes Module 2's search space and therefore placement
FDR; Khoa accepts it). With B: `makedb --patch` writes the read-through residue of a
`cterm_extension` patch as `X` (decision for the ADR); the M7-12 simulator plants Sec/read-through
(the backlog line); TE1's fdnG misses become `stop_codon` ALT sites and TE2 exercises its second
exception.

## Acceptance criteria (if B is accepted)

- [ ] The M3 tolerant property tests (zero false negatives vs brute force) extended to the
  annotated-stop rule; the M3-21 placement calibration re-run green (fast + full).
- [ ] `stop_codon` ALT sites on the Nissle TE1/TE2 set at fdnG (and any other Sec genes).
- [ ] `cterm_extension` reachable and patchable; ADR-0013 consequences tested end to end.

## Questions for Khoa (through the lead)

1. Option A, B or C? **Proposed default:** B, by ADR.
2. With B, how does `makedb` translate the read-through codon of a patch? **Proposed default:** `X`
   (search engines treat it as any residue), since Sec (`U`) and Trp (`W`) read-through both occur.
