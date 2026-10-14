# M3-04: Route 1, exact peptide matching (Aho–Corasick)

Milestone M3 · depends on: M3-02, M3-03 · creates `src/giap/mapping/exact.py`,
`src/giap/mapping/hits.py` (shared hit schema) · pure Python (`pyahocorasick`, already a
dependency) · estimate ≈ 260 code lines

## Spec sections implemented

- §7.2.2 "Build an Aho–Corasick automaton (`pyahocorasick`) of all unique query peptides (I → L).
  Scan every target and decoy frame in one pass. Keep every locus. Peptides under 7 aa were dropped
  in preflight: a 6-mer has about one chance hit in a 5-Mb genome."
- §7.2 "Each sequence takes the cheapest route that can place it; only failures move on to the
  next step."
- Decision 5: "GIAP's own exact six-frame index for 7–40 aa." Config `mapping.exact: {min_len: 7,
  max_len: 40}`. §7.2.9 unit: "I = L; origin-spanning peptides on circular contigs; no wrap on
  linear contigs; Aho–Corasick multi-locus".

## Public interface

`src/giap/mapping/hits.py` — the hit table shared by M3-04, M3-07, M3-08 (work file
`hits_exact.parquet` / `hits_tolerant.parquet`, one row per peptide × locus × search):

| column | type | meaning |
|---|---|---|
| `peptide` | str | aggregated stripped peptide as in the input (I/L as written) |
| `query` | str | I → L normalised |
| `search` | enum{target, decoy} | |
| `decoy_type` | enum{reverse, shuffle} or NA | NA for target rows |
| `contig`, `frame_id`, `aa_start`, `aa_end` | | index coordinates (decoy rows: in the decoy frame) |
| `start`, `end`, `strand` | int, int, str | 0-based half-open, unwrapped; decoy rows: the **mirrored target** interval |
| `frame` | int | ±1..3 (M3-02 rule) |
| `route` | enum{exact, tolerant} | |
| `n_diff` | int | 0 for exact |
| `diff_types` | str | `;`-separated, sorted; empty for exact |
| `blocks` | str or NA | frameshift hits only: `start-end;start-end` (M3-08) |

```python
HIT_COLUMNS: tuple[str, ...]
def match_exact(queries: Sequence[str], index: SixFrameIndex,
                decoys: Mapping[DecoyType, tuple[Frame, ...]],
                decoy_of: Mapping[str, DecoyType]) -> pd.DataFrame: ...   # HIT_COLUMNS
```

Queries outside `[min_len, max_len]` are not searched (caller reports them). Duplicate queries
(after I → L, e.g. `PEPTIDEK`/`PEPTLDEK`) share one automaton key; each input peptide gets its own
rows. Target hits are de-duplicated with `dedupe_loci` (origin copies); decoy hits on their
mirrored interval. Each peptide is scanned against the decoy of its own source (`decoy_of`).
Rows are sorted by peptide, search, contig order, start, end, strand.

## Examples (genome of M3-02 unless stated)

| # | Input | Expected |
|---|---|---|
| 1 | query `MAKILW` on `c1` (test with `min_len=4`) | one target row: `Interval(c1,0,18,+)`, frame `+1`, `n_diff=0` |
| 2 | query `MAKLLW` | the same locus (I = L); `query` = `MAKLLW` |
| 3 | 50-kb sim genome, peptide cut from gene `GSIM_00005` | ≥ 1 target row whose interval lies inside the gene, frame of the gene |
| 4 | a peptide present twice in the genome (planted duplicate 30 nt) | two target rows ("keep every locus") |
| 5 | circular `r` of M3-02 (`wrap_nt=12`), query `MPKEL` (`min_len=4`) | one row `Interval(r,18,33,+)` (crosses origin), not two |
| 6 | the same genome declared linear | no target row |
| 7 | query containing `U` (`SEUKLLAAR`) | no rows (NOMATCH never matches) |
| 8 | query spanning a stop (`W*RX…`) | impossible (queries are A–Z); queries never match across `*` |
| 9 | database peptide whose reversed frame contains it | one `search=decoy, decoy_type=reverse` row with the mirrored interval |

## Edge cases

- Overlapping occurrences (`AAAAAAAA` in `AAAAAAAAAA`): every start offset is a locus (3 loci).
- Peptide on both strands: two loci with opposite strands.
- An empty query list returns an empty frame with `HIT_COLUMNS` and dtypes.

## Acceptance criteria

- [ ] `tests/unit/test_exact.py`: examples 1–9 and the edge cases.
- [ ] Property ("a peptide cut from a random genome maps back to its origin", §7.2.9): random
  genomes ≤ 50 kb (both topologies), random windows of 7–40 aa from any frame (origin-crossing
  ones on circular contigs) → the origin interval is among the target rows.
- [ ] Property: the result equals a brute-force `str.find` loop over every frame (target and
  decoy).
- [ ] Benchmark: 100,000 random tryptic peptides on a 5-Mb genome (target + reverse decoy)
  ≤ 30 s and ≤ 2 GB peak on one core.

## Questions for the lead

1. Should exact matching also run for peptides that will go to the tolerant route (≥ 10 aa)?
   **Proposed default:** yes; every peptide is first searched exactly (target and decoy). The
   tolerant route runs only for target leftovers, but the **decoy** tolerant search runs for every
   peptide ≥ 10 aa (see M3-10 Q1).

## Lead decisions after code review (PR #27)
1. `match_exact` validates up front: every value of `decoy_of` must be a key of `decoys`, and every searched
   peptide must be in `decoy_of`; otherwise `ValueError` naming the peptide or decoy type (a silent missing
   decoy search would make placement FDR optimistic).
2. Tests: two I/L-twin peptides sharing one automaton key with different decoy types; the replicon-length
   guard on a ~30-nt circular replicon; a `U` query is never added to the automaton.
3. `benchmarks/__init__.py` makes benchmark modules a package, so `pytest tests benchmarks` collects
   `benchmarks/test_exact.py` next to `tests/unit/test_exact.py` without an import-file mismatch.
4. `search`/`decoy_type` stay object columns for now; categorical dtypes are decided with M3-07/08 output.
