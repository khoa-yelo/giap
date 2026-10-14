# M3-08: Tolerant matching, ±1-nt frameshift splits

Milestone M3 · depends on: M3-02, M3-03, M3-06, M3-07 · creates `src/giap/mapping/frameshift.py`
· pure Python + numba · estimate ≈ 260 code lines

## Spec sections implemented

- §7.2.4 "a ±1-nt frameshift split: prefix in frame f, suffix in frame f ± 1 at the adjacent
  nucleotide offset, with ≥ 4 aa on each side" (costs one difference);
  "Frameshift splits: anchor the longer side (exact, or within the remaining budget), then verify
  the shorter side in the shifted frame."
- Config `mapping.tolerant.frameshift: {enabled: true, min_side_aa: 4}`.
- §9.4 placement flag `fs` is for proteins (miniprot); peptide frameshifts are reported through
  `diff_types = frameshift` and BED12 blocks (M3-16).

## Semantics

On one strand, in that strand's nucleotide coordinates `x`: a split at query position `p`
(`min_side_aa ≤ p ≤ L − min_side_aa`) places the prefix `q[0:p]` on codons `[x, x + 3p)` and the
suffix `q[p:L]` on codons `[x + 3p + δ, x + 3p + δ + 3(L − p))`, `δ ∈ {+1, −1}` (+1: one extra
genome nucleotide, −1: one nucleotide read twice). Cost = 1 + the M3-06 cost of each side against
its frame (the remaining budget `k − 1` may be spent on either side). The locus is
`[x, x + 3L + δ)` with two blocks; `frame` is the prefix frame.

**Decoys, searched identically.** Decoy frames have no nucleotides; a decoy split is two decoy
segments of the same strand whose **mirrored target** intervals (M3-03 `mirror_aa` + `to_interval`)
satisfy the same adjacency rule.

```python
def match_frameshift(queries: Sequence[str], index: SixFrameIndex,
                     decoys: Mapping[DecoyType, tuple[Frame, ...]], decoy_of: Mapping[str, DecoyType],
                     target_queries: set[str], cfg: TolerantConfig) -> pd.DataFrame: ...  # HIT_COLUMNS
```
Rows have `route = tolerant`, `diff_types` containing `frameshift`, and `blocks` =
`"s1-e1;s2-e2"` (0-based, genome coordinates, ascending).

## Examples

Sequence `GGCCGGAACCGACCATTGATAGAAAAACTGATGAACGCGTAA` (42 nt, linear `fs1`): frame `+3` reads
`PEPTIDRKTDERV`, frame `+1` reads `GRNRPLIEKLMNA*` (an `A` was inserted after `…GAT`).

| # | Input | Expected |
|---|---|---|
| 1 | query `PEPTIDEKLMNA` (12 aa, k = 2) | one target row: `Interval(fs1, 2, 39, +)`, `frame=+3`, `n_diff=1`, `diff_types=frameshift`, `blocks="2-20;21-39"` |
| 2 | same query with `min_side_aa=7` | no frameshift row (both sides are 6 aa) |
| 3 | query `PEPTIDEKLMNA` with the prefix carrying one substitution (`PEPTVDEKLMNA`) | `n_diff=2`, `diff_types="frameshift;substitution"` |
| 4 | query of 10 aa (k = 1) with a frameshift and a substitution | no row (cost 2) |
| 5 | reverse complement of the sequence, same query | one row on strand `-` covering `[3, 40)` with the mirrored blocks `3-21;22-40` |
| 6 | `δ = −1` case: a genome where one nucleotide is deleted at the junction | a row with `blocks` whose second block starts 1 nt before the first ends (`s2 = e1 − 1`) |
| 7 | frameshift disabled (`enabled: false`) | empty frame |

## Edge cases

- Overlapping blocks for `δ = −1` are legal in GIAP's model but BED12 forbids overlapping blocks:
  the writer (M3-16) shortens the first block by 1 nt and records the read-twice base in the GFF3
  attributes.
- Split points that land on a stop in either frame are rejected (no operation consumes `*`).
- Origin-crossing splits on circular contigs work within the 120-nt extension (M3-02 Q2).
- A query with both an ordinary tolerant hit and a frameshift hit at equal minimal cost: both
  loci are kept (they are different loci); at the same locus, the non-frameshift explanation wins
  (precedence, M3-05).

## Acceptance criteria

- [ ] `tests/unit/test_frameshift.py`: examples 1–7.
- [ ] Property (completeness, §7.2.4 "for every difference type"): random genomes ≤ 50 kb, query
  built from two adjacent frames with a random split and `δ = ±1` (plus 0..k−1 other planted
  differences) → found; equals a brute-force scan over all split points and positions.
- [ ] Property: strand symmetry under reverse complement (example 5 generalised).
- [ ] Benchmark: adds ≤ 50 % to the M3-07 benchmark runtime.

## Questions for the lead

1. Should frameshift loci count toward detection and abundance like other tolerant hits under
   `genome_relation = related`? **Proposed default:** yes, as the spec draws no distinction; they
   also carry `strain_difference` (M3-14).
