# Fix: BacTermFinder parser reads the real output table

## Problem
`gippa.annotate.terminators.parse_bactermfinder` (M2-15, provisional) expects a header that names
contig, 1-based start/end, strand and probability columns. At the pinned commit `024870447f62`
(ADR-0008), `genome_scan.py` writes `<prefix><genome>_mean.csv` with the header
`SampleName,probability_binary,probability_ENAC,probability_PS2,probability_NCP,probability_mean`,
so every real run fails in the parser (backlog M2-15).

## Real format (recorded: `tests/data/recorded/bactermfinder/024870447f62/`)
- `SampleName` = `<u_id>_<seq_id>_<start>_<end>_<strand>` (genome_scan.py lines 81-83). `u_id` is an
  integer that restarts for each FASTA record. `seq_id` is the FASTA record id and may contain `_`.
- `start = i`, `end = i + 101` for the slice `seq[i:i + 101]` (lines 44-53). These are **0-based
  half-open** coordinates, as in GIPPA (§6.5). No shift is applied.
- Each window has a `+` row and a `-` row with the same coordinates. Every window is reported.

## Required behaviour
1. `parse_bactermfinder(path, replicons) -> list[Terminator]` keeps its signature. It returns one
   `Terminator` per data row, in file order, with
   `Interval(seq_id, int(start), int(end), strand)`, `confidence = float(probability_mean)` and
   `source = "BacTermFinder"`. The interval is validated against `replicons`.
2. Splitting: take `u_id` from the left at the first `_`, then split the rest from the right at
   its last three `_` (`rest.rsplit("_", 3)`). `u_id`, `start` and `end` must be integers and
   `strand` must be `+` or `-`.
3. Errors (`ValueError`): a missing `SampleName` or `probability_mean` column names the file. A
   malformed `SampleName`, a bad strand, or an unknown contig or out-of-range interval names `file:line`.
4. An empty file, or a file with only the header, gives `[]`.

## Out of scope / questions for the lead
- Q1: should the provisional contig/start/end/strand header still be accepted? It never matched a
  real output. **Proposed: drop it.** This means updating the `BTF_HEADER` fixtures in
  `tests/unit/test_terminators_review.py` and `test_bactermfinder_still_says_why` in
  `tests/unit/test_m2_21_contract_coverage.py` (lead).
- Q2: since every window is reported (about 2 rows per 3 bp), should the step keep only rows with
  `probability_mean` >= a cut-off, and merge them as the README suggests? The README gives 0.13,
  0.3 and 0.47 (F2/F1/F0.5). **Proposed: a follow-up config key, default 0.3.** The parser does not filter.
## Tests
`tests/contract/test_fix_bactermfinder_parser.py`: the recorded excerpt and smoke head; spot rows;
the 0-based convention; contig ids with `_` split from the right; a Hypothesis round trip; errors.
`tests/contract/test_bactermfinder_contract.py` loses `pending_sample` and reads the new recording.
