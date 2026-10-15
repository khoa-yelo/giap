# M3-07: Tolerant matching, q-gram seeds and the target/decoy search driver

Milestone M3 · depends on: M3-02, M3-03, M3-04 (hit schema), M3-06 · creates
`src/giap/mapping/qgram.py`, `src/giap/mapping/tolerant.py` · pure Python + numba ·
estimate ≈ 380 code lines

## Spec sections implemented

- §7.2.4 "**Difference budget k:** 1 at 10–11 aa, 2 at ≥ 12 aa." (config
  `mapping.tolerant.max_diff: {"10-11": 1, "12-40": 2}`, `min_len: 10`).
- "**Seed and verify (G4).** Index the six-frame (and decoy) translations by q-grams
  (q ∈ {3, 4, 5}). For each peptide, choose the largest q for which the q-gram lemma gives a
  lossless threshold t ≥ 1, where t = (L − q + 1) − k(q + 1); each difference destroys at most
  q + 1 q-grams. Count shared q-grams per diagonal band (diagonal ± 1, which absorbs 1↔2
  exchanges). Verify candidate diagonals with a banded dynamic program (band ± 2) … Implement the
  hot loops in numba."
- "**Completeness requirement.** Property tests MUST show zero false negatives against a
  brute-force matcher on random genomes ≤ 50 kb, for every difference type (§11.2)."
- "**Output.** For each peptide: all target loci at the minimal cost, the cost and the difference
  types; and the minimal decoy cost (search decoys identically)."
- Decision 6; §6.2 step `map_tolerant`: "leftover peptides ≥ 10 aa, index →
  `hits_tolerant.parquet` (target + decoy)".

## Public interface

```python
def choose_q(L: int, k: int) -> tuple[int, int]: ...          # (q, t); ValueError if none
@dataclass(frozen=True)
class QGramIndex:                                              # CSR over concatenated frames
    q: int; offsets: np.ndarray; positions: np.ndarray         # int32 global positions
    frame_starts: np.ndarray                                   # frame boundaries (no q-gram crosses one)
def build_qgram_index(frames: Sequence[Frame], q: int) -> QGramIndex: ...   # q-grams with * or X skipped
def candidate_diagonals(query: np.ndarray, idx: QGramIndex, t: int) -> np.ndarray: ...
def match_tolerant(queries: Sequence[str], index: SixFrameIndex,
                   decoys: Mapping[DecoyType, tuple[Frame, ...]], decoy_of: Mapping[str, DecoyType],
                   target_queries: set[str], cfg: TolerantConfig) -> pd.DataFrame: ...  # HIT_COLUMNS
```

- **Seeds.** For each query q-gram at offset `i` and each occurrence at global position `p`,
  diagonal `d = p − i`. A centre diagonal `c` is a candidate when the q-gram count summed over
  `c − 1, c, c + 1` is ≥ `t` (each query q-gram counted at most once per band). Candidates closer
  than 3 diagonals are merged; `verify` (M3-06) runs once per merged candidate with `d0 = c`.
- **Target search** only for `target_queries` (leftovers: no exact target hit, length ≥ `min_len`).
  **Decoy search** for every query ≥ `min_len` against its own decoy type (M3-10 Q1).
- **Per peptide and search:** keep all loci at the minimal cost (≤ k); report nothing above k.
  Target loci are converted with `to_interval`, de-duplicated with `dedupe_loci`; decoy loci keep
  their mirrored target interval (M3-03).
- Indexes are built lazily, only for the q values needed, and shared by target and decoy frames
  of one type. `numba.prange` over queries; each query's result is independent, then rows are
  sorted, so output does not depend on thread count (§6.6).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `choose_q(10, 1)` | `(4, 2)` (q = 5 gives t = 0) |
| 2 | `choose_q(11, 1)` | `(5, 1)` |
| 3 | `choose_q(12, 2)` | `(3, 2)` (q = 4 gives −1) |
| 4 | `choose_q(13, 2)` | `(3, 3)` |
| 5 | `choose_q(14, 2)`; `choose_q(20, 2)`; `choose_q(40, 2)` | `(4, 1)`; `(5, 4)`; `(5, 24)` |
| 6 | 50-kb sim genome; a 12-aa window of a gene with one N → GG change planted in the query | one target locus at the gene, `n_diff=1`, `diff_types="block"` |
| 7 | same window with a substitution and a swap planted | `n_diff=2`, `diff_types="substitution;swap"` |
| 8 | 10-aa query with 2 planted substitutions | no target row (k = 1) |
| 9 | 9-aa query | not searched (below `min_len`) |
| 10 | query whose window is duplicated in the genome, one copy exact and one with 1 substitution | in `match_exact` only (not a leftover): tolerant target search not run; decoy search still runs |

## Edge cases

- A diagonal range crossing a frame boundary: candidates are clipped to the frame; q-grams never
  span two frames (separator in the concatenation).
- Circular extension: windows found in the extension are reduced modulo L and de-duplicated
  (M3-02).
- Low-complexity queries (`AAAAAAAAAAAA`): many candidates; the result stays correct, and the
  benchmark below includes such queries.
- Queries with `NOMATCH` residues: their q-grams are skipped; `t` is computed from the full
  length, so a query with too many NOMATCH residues may find nothing (documented, not an error).

## Acceptance criteria

- [ ] `tests/unit/test_tolerant_seeds.py`: examples 1–10.
- [ ] Property, **completeness** (§7.2.4, §11.2): random genomes ≤ 50 kb (both topologies),
  queries cut from random frames with 0..k planted differences, each of `substitution`, `swap`,
  `block` (every pair of `BLOCKS`), `near_block` — the full result (set of loci at minimal cost,
  and the cost) equals a brute-force scan of every frame window with the M3-06 oracle. Zero
  false negatives; also zero false positives. Same for the decoy search.
- [ ] Property: thread count 1 vs 4 gives identical frames.
- [ ] Benchmark (`benchmarks/test_tolerant.py`): 10,000 leftover peptides (10–30 aa) on a 5-Mb
  genome, target + one decoy, ≤ 5 min on 8 CPUs and ≤ 2 GB peak; 1,000 peptides on 100 kb ≤ 5 s
  (the calibration suite depends on it).

## Questions for the lead

1. Peptides of 10–11 aa with k = 1 and a frameshift: handled by M3-08, not here. Confirm the
   split. **Proposed default:** yes; `match_tolerant` returns the non-frameshift result and M3-08
   adds frameshift candidates, the minimum of both is taken in M3-10.

## Lead decisions (test-author questions)
1. Batch entry point (M3-06 decision 6): `tolerant_dp.verify_batch(query, frame, d0s, k, ex)`, a numba
   dispatcher returning int64 rows (key, t_start, t_end), key −1 when cost > k; public
   `exchange_lookup(tables)` and `decode_key(key)`. Agrees with `verify` row by row, callable inside an
   njit `prange`, raises for L ≤ 4; `match_tolerant` compiles it only for uint8 C-contiguous inputs.
2. Queries of ≤ 4 aa yield no rows and no error, even if the config lowers `min_len`.
3. Recall: a planted peptide is recalled if a reported locus overlaps its origin in the same codon phase,
   or a strictly cheaper locus is reported (the planted one is then not the best explanation). Overlap
   recall must reach `mapping_recall_min`. Code review checks this cannot mask real misses.
