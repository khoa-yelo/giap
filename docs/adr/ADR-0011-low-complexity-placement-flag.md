# ADR-0011: `low_complexity` placement flag (M3-10)
- Status: proposed
- Date: 2026-10-14
- Spec sections: §9.4 (placement flags), §7.2.5; M3-index lead decision (M3-10 `low_complexity`)
- Needs Khoa: yes (adds a value to an output vocabulary)

## Context
The M3-09 statistician attack showed that the chance-hit score E underestimates hits of simple-sequence
(repeat-like) peptides by 1.9× at L = 8 and 15× at L = 12. Placement FDR stays valid (decoys calibrate it),
but users should see which accepted placements rest on low-complexity peptides. §9.4's placement flag
vocabulary has no such value.

## Decision
Add `low_complexity` to the §9.4 placement flags: set when the peptide has ≤ 3 distinct residues
(I and L counted as one). Flags stay `;`-joined. No effect on q-values.

## Evidence
docs/methods/placement-fdr.md (chance-hit section, signed): E vs observed for simple-sequence repeats.

## Consequences
mapping outputs gain one flag value; report text explains it. No threshold or default changes.
