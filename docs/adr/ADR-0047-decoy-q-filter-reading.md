# ADR-0047: Decoys above `q_value_max` leave the picked competition (reading of §5.2)
- Status: accepted-by-khoa (2026-10-22: the reading, and the defaults.yaml comment reworded)
- Date: 2026-10-22
- Spec sections: §5.2, §7.3.2; M4-12 decision 1; M4-24
- Needs Khoa: yes, for the wording change to the protected `src/giap/config/defaults.yaml`
  comment; Khoa approved it on 2026-10-22. The code already behaves this way.

## Context
§5.2 and the `input.q_value_max` comment in defaults.yaml say that rows above the cut-off are
"still used to estimate FDR". Since M4-12, `detection.pair_decoys` has excluded decoys with
q > `q_value_max` before pairing and picking. The M4-24 spec (unscored decoys) had to say which
reading holds.

## Decision
- Decoys above `q_value_max` are excluded before pairing and picking. They stay in C9 and in
  preflight's q-range summary. That is how "still used" is read: kept, not discarded at input.
- With q in the shared key (M4-24 rule 1), the exclusion changes no gene q of a target that has a
  q. It changes only f, the 80 % guard and the decoy counts.
- The filter applies only when some C9 target row carries q (M4-24 rule 6, a symmetry rule).

## Evidence
- M4-24 examples: the DIA-NN-shaped input and the sim50k golden are byte-identical under this
  reading.
- Slurm 1869554: on sim100k, scored decoys give the same detection counts with and without the
  filter at α 0.01 and 0.05.

## Consequences
- No code change for the filter itself.
- The `input.q_value_max` comment in defaults.yaml now reads "not evidence; kept in C9 and in
  preflight's q range", a protected-file change with `Gate-Change-Approved: ADR-0047`. The
  default value is unchanged.
