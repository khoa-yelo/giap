# UniOP format excerpts (unit-test fixtures, M2-12)

Source: a real UniOP run (pin 3ed14dfc, run as a process by the M2-12 reviewer, PR #49) on the GIAP
`uniop_input` GFF of giapsim SIM100K (seed `SIM100K_SEED`; full input sha256
`1d9c91e955d816df7b3e23890661af3e430f09317465dc6f479d233b7049053e`, 87 genes on `chr1` and `plasmid_1`).

- `multicontig.gff`: lines 1-7 (`chr1_1`..`chr1_6`) and the two `plasmid_1` lines of that input, verbatim.
- `multicontig_map.tsv`: the matching `uniop_map.tsv` rows, verbatim.
- `multicontig_uniop.pred`: the first 9 lines of the real `uniop.pred`, byte-for-byte (title line,
  header, rows 1-7 with UniOP's trailing space and blank-padded empty probabilities).

UniOP numbers genes after a stable sort by start across contigs, so rows 1-7 of the full run name
exactly the first 8 sorted genes, all of which are in the excerpt (plasmid_1_1 at 193 sorts first).
The probabilities were computed on the full genome; they are copied, not recomputed. UniOP itself is
not vendored (ADR-0008).
