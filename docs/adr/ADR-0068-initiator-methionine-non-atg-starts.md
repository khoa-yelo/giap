# ADR-0068: Initiator methionine at annotated non-ATG starts
- Status: accepted-by-khoa (Accepted by Khoa, 2026-10-28, in session, Q-K56: "Yes")
- Date: 2026-10-28
- Spec sections: §2.2 improvement 6 and §15 ("non-ATG starts: No"), §7.2.1 (six-frame exact
  mapping); docs/specs/M11-20-initiator-methionine.md
- Needs Khoa: answered 2026-10-28 (a §2 decision changes)

## Context
Bacteria translate the first codon of a CDS as methionine whether it is ATG, GTG or TTG. The
six-frame index reads GTG and TTG as V and L, so the N-terminal peptide of a non-ATG-start gene
(M first) never maps exactly to its own gene: recall is lost, and a chance copy elsewhere can
become an unestimated false novel ORF or placement. On E. coli K-12 (Bakta) 402 of 4,373 CDS and
sORF starts are non-ATG. GIPPA_SPEC §15 accepted the miss ("non-ATG starts: No"), and the TE1 proxy
test pins it. The gap surfaced as code-25 G4 failures in the M11-16f calibration (M7-13
decision 14).

## Decision
Credit M at query position 0, in the exact route only, at the run's own annotated non-ATG
CDS/sORF starts whose codon is in `annotation.start_codons` (option A of M11-20). The start
sites become cells in the six-frame index (format 4), as ADR-0058 does for annotated stop
cells; decoys get mirrored cells (reverse: cell `n − 1 − p` at the last position of a mirrored
window; shuffle: cell `p`). Hits keep the table schema and gain the placement flag
`initiator_m`; the variants step treats the credited residue as an `initiator` op with no
pileup contribution. No catalogue of start codons is searched: only the run's annotation sets
cells, before any peptide is seen. Tolerant-route credit and credit at every GTG/TTG stay
follow-ups that need their own statistics review.

## Evidence (before → after: metrics, tests, simulations)
- Statistician review of M11-20 (2026-10-28): approved after the mirrored reverse window. With
  it, usable credited windows and i.i.d. chance mass are exactly equal for target and reverse
  decoy; E_i, strata and the §8.1 families are unchanged; the extra chance windows are 0.25 % of
  M-initial windows by count (0.30–0.47 % by chance mass) on K-12.
- Required before merge: `tests/calibration/attacks/test_attack_initiator_decoys.py` (exact mass,
  four mutant teeth) and arms I1 (placement FDR level) and I2 (code-25 G4 re-run, recurrence 0).

## Consequences
- Overrides GIPPA_SPEC §2.2 improvement 6 and §15 for annotated non-ATG starts; the TE1, TE2 and
  TE8 proxy expectations change (listed in M11-20).
- Six-frame index format 4; formats 2 and 3 load with no cells and a WARNING.
- M7-13 decision 14's exemption narrows to genes the site rule skips once 20i lands.
