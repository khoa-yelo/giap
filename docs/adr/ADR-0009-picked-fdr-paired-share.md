# ADR-0009: Picked gene-level FDR rescales decoy wins by the paired share
- Status: accepted-by-khoa (2026-10-14)
- Date: 2026-10-14
- Spec sections: §7.3.2 (decision 4, picked target–decoy tier), §8.1
- Needs Khoa: answered

## Context
§7.3.2 uses picked gene-level target–decoy FDR when decoy rows can be paired with genes, and falls back
to MAYU only when more than 20 % of decoys are unpaired. Unpaired decoys (proteins absent from the
strain, IDs that do not map to the gene model) are dropped. The statistician's review of M4-01
(PR #19; docs/methods/target-decoy.md) showed that dropping 20 % of decoys at random inflates the true
picked FDR to about 1.2× nominal (0.059 at q = 0.05; 0.122 at q = 0.10).

## Decision
Picked gene-level FDR rescales the decoy-win count by the paired share: with a paired fraction
f = n_paired_decoys / n_decoys, each decoy win counts as 1/f. The 20 % fallback to MAYU stays as an
upper guard. The paired share and the rescaling factor are reported in `detection_summary.tsv`.

## Evidence (before → after: metrics, tests, simulations)
Before: 20 % random decoy loss → FDP ≈ 1.2 × q. After: M4-12 must show FDP ≤ q + 2 SE under
10 %, 20 % and 30 % random decoy loss (calibration suite, plus an attack null in
`tests/calibration/attacks/`).

## Consequences
Cross-strain searches (pilot arm A1) keep the picked tier with an honest estimate instead of falling
back to MAYU more often.
