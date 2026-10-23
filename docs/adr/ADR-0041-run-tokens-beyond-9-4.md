# ADR-0041: Run-scope tokens beyond the §9.4 Run row
- Status: accepted-by-khoa
- Date: 2026-10-22
- Spec sections: §9.4 (Run row), §7.1, §6.7
- Needs Khoa: yes (vocabulary; decided by Khoa on 2026-10-22)

## Context
`run_manifest.json` (MVP-03, ADR-0040) collects run-scope flags. Four tokens that steps already emit, or that
the manifest adds, are not in the §9.4 Run row. ADR-0036 treated additions to that row as Khoa-level, so MVP-03a
kept them in `RUN_FLAGS_PENDING`. The schema accepts them, and the vocabulary test names them as pending.

## Decision
The §9.4 Run row gains these four tokens:

| Token | Emitted by | Meaning |
|---|---|---|
| `input_prefiltered` | `placement` | the peptide input was already FDR-filtered upstream, so the placement FDR is conditional on that filter |
| `cai_reference_missing` | `codon_usage` | no ribosomal-protein reference set was found, so CAI is NA |
| `no_detection_limit_model` | `operons` | the detection-limit model was not fitted, so missed-operon calls fall back (also a `pairs` row token) |
| `uncalibrated:<key>` | the manifest (MVP-03c), from `src/giap/data/calibration.yaml` | the claim `<key>` has no passing calibration on record (§7.1 "`uncalibrated` in the manifest", parameterised by claim) |

- `RUN_FLAGS_PENDING` becomes part of the accepted run vocabulary.
- Steps that only log `cai_reference_missing` or `no_detection_limit_model` now emit an event (MVP-03a).

## Consequences
- The manifest schema pattern is unchanged. The vocabulary test moves the four tokens from pending to accepted
  in a follow-up to #196.
- The report (M9-16) lists them with the other run flags.
