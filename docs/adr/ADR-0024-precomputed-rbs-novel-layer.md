# ADR-0024: `annotation.precomputed.rbs_novel`, a precomputed second-round OSTIR layer
- Status: proposed
- Date: 2026-10-21
- Spec sections: §6.4, §7.6 (second-round OSTIR), §13
- Needs Khoa: yes (a new config key, which changes the config interface)

## Context
The local executor reads every tool result from a precomputed layer, and a missing layer is an error (§6.4). Second-round OSTIR scores start codons that exist only after mapping, and these are not keyed by `locus_tag`. Without a source for them, the local executor cannot run `NOVEL`.

## Decision
Add `annotation.precomputed.rbs_novel` to `defaults.yaml`, null by default.
- It holds a TSV keyed by position (contig, 1-based start-codon position, strand), with the OSTIR columns of the `rbs` layer.
- Nextflow honours the key the same way as the other precomputed layers, and runs OSTIR itself only when the key is null (M7-08 lead decision 6).
- The rejected alternative was to put `START_<contig>_<pos1>_<p|m>` rows into the `rbs` layer. It needs no new key, but it overloads `locus_tag`.

## Evidence (before → after: metrics, tests, simulations)
- Before: the local executor could not run `NOVEL`.
- After: the M7-08 tests run it from a layer the simulator writes.

## Consequences
- The input-path keys gain one entry, and no threshold changes.
- If Khoa rejects the key, M7-08 switches to the locus-keyed alternative. That needs no data change, because the simulator writes both forms of the layer.
