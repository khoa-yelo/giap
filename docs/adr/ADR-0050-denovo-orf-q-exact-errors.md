# ADR-0050: ORF-level q in the de novo subgroup regime under exact systematic errors (G4-k3, attack GC)
- Status: proposed
- Date: 2026-10-23
- Spec sections: §7.6 (novel ORFs), §8.4 (class FDR), §8.5 (attacks), §13 gate M7; M7-11 decision 3; M7-13 rule 8
- Needs Khoa: yes. Options 4 and 5 change a gate generator or the call rule (R3); option 3 records a
  limitation of the class FDR.

## Context
M7-11 decision 3 sets a de novo ORF's class q to the minimum of its peptides' subgroup q. The
methods note flags this as anti-conservative ("k chances"). M7-13 rule 8 (G4-k3) measures it: 20
false ORFs per six-frame class, 3–5 wrong de novo reads each (k//2 reversed stretches, the rest
exact segment products).

**Full tier** (R = 300, Slurm 1875441 / 1876296 / 1877124), confirmed in both blocks:

| Check | Block 1 FDP | Block 2 FDP | Bound |
|---|---|---|---|
| intergenic, q 0.01 | 0.0168 | 0.0153 | 0.0114 |
| cds_outframe, q 0.01 and 0.05 | 0.0747 | 0.0749 | 0.0102 / 0.0502 |
| antisense, q 0.01 and 0.05 | 0.0719 | 0.0719 | 0.0102 / 0.0502 |

FDP is flat in q: the false ORFs sit at the top of the list (median minimum peptide q 5e-4).
Attack GC (exact products of GC-poor intergenic segments) fails at q 0.01 by the same mechanism
(0.0122 / 0.0126 against 0.0112).

**Decomposition** (100 paired replicates):
- exact products only: FDP 0.075 / 0.072, which is the whole failure;
- reversed stretches only: FDP 0.0098 / 0.0135 / 0.0140 against about 0.011–0.0125, a residual
  excess at q 0.01 only.

**Why no decoy can fix it:** a de novo read that exactly equals a non-coding translation cannot be
told apart from translation evidence by any decoy-based estimate, because a reversed or shuffled
genome decoy matches it exactly only by chance. The FDR estimand assumes wrong reads behave like
decoy hits. That holds for error-bearing reads, up to the k-chances effect, and fails for exact
systematic errors.

## Options
1. **Picked ORF-level competition on the genome decoy** (target ORF against its mirrored decoy ORF,
   aggregated score). Addresses the k-chances excess (the reversed-stretch arm). It cannot fix the
   exact arm, because an exact target ORF always beats its decoy.
2. **Per-ORF Bonferroni over peptides** (q_orf = min(1, k·min q)). Measured: no fix (cds_outframe
   0.076, antisense 0.074), and worse for intergenic (0.0243 against 0.0174 at q 0.01), because it
   drops true multi-peptide ORFs.
3. **Keep decision 3 and declare a documented limitation:** exact systematic de novo errors are
   outside what the class FDR controls.
   - G4-k3 (cds_outframe and antisense at 0.01 and 0.05, intergenic at 0.01) and GC at 0.01
     become strict xfails naming this ADR.
   - De novo-only calls with no database or RBS/TIR support are flagged.
4. **Re-specify rule 8 and attack GC to error-bearing wrong reads only** (a chance-hit null), and
   keep the exact arm as a recorded diagnostic. This changes a gate's generator, so it needs Khoa
   (R3 forbids relaxing a check). The reversed-only arm still shows an excess at q 0.01, so option
   1 or 3 would still be needed.
5. **Require independent support for de novo-only six-frame calls** (≥ 1 database PSM, or a
   planted RBS/TIR pass). This changes the call rule and power (G3) and needs an M7-11 amendment.

## Lead's interim default (reversible; pending Khoa)
Option 3's test side only: the failing G4-k3 and GC q 0.01 checks become strict xfails naming
ADR-0050, so a fix shows as XPASS. The flag for unsupported de novo-only calls is not added until
Khoa decides. The gate (M7-13c) lists these as open known failures, not as met.

Source: the M7-13b statistician's full-tier report (2026-10-23); its scripts and logs are in
`local/stat/m713b-full/`.
