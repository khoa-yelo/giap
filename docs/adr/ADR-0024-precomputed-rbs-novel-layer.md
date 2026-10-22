# ADR-0024: `annotation.precomputed.rbs_novel`, a precomputed second-round OSTIR layer
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §6.4, §7.6 (second-round OSTIR), §13
- Needs Khoa: yes (a new config key, which changes the config interface). Answered: option (a).

## Context
The local executor reads every tool result from a precomputed layer, and a missing layer is an error (§6.4). Second-round OSTIR scores start codons that exist only after mapping, and these are not keyed by `locus_tag`.

## Decision
Add `annotation.precomputed.rbs_novel` to `defaults.yaml`, null by default.
- It holds a TSV keyed by position (contig, strand, 1-based first base of the start codon in transcription direction), with the OSTIR columns of the `rbs` layer.
- Nextflow honours the key the same way as the other precomputed layers, and runs OSTIR itself only when the key is null (M7-08 lead decision 6).
- Rejected: `START_<contig>_<pos1>_<p|m>` rows in the `rbs` layer. That needs no new key, but it overloads `locus_tag`.

## Evidence (before → after: metrics, tests, simulations)
- Before: the local executor could not run `NOVEL`.
- After: the M7-08 tests run it from a layer the simulator writes.

## Consequences
The input-path keys gain one entry. No threshold changes.
