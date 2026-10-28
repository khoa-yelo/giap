# M3-06: Tolerant matching, banded verification DP (numba)

Milestone M3 · depends on: M3-02 (codes), M3-05 (operation tables) · creates
`src/gippa/mapping/tolerant_dp.py` · pure Python + numba · estimate ≈ 300 code lines

## Spec sections implemented

- §7.2.4 "**Each of these costs one difference:** a substitution; an adjacent swap (XY ↔ YX); a
  mass-equivalent block exchange …; a near-equal exchange …" "**I/L** is free."
- "Verify candidate diagonals with a banded dynamic program (band ± 2) whose only operations are
  those listed, each costing 1." "Implement the hot loops in numba."
- "**Output.** For each peptide: all target loci at the minimal cost, the cost and the difference
  types".
- G4: "banded verification that charges each allowed difference type as one difference".

(Seeding is M3-07; frameshift splits are M3-08.)

## Semantics (normative for implementer and test author)

Align the whole query `q[0:L]` (global) to a target substring `t[s:e]` (free ends) of one frame.
An alignment is a sequence of operations, each consuming `a` query and `b` target residues:

| operation | (a, b) | condition | cost |
|---|---|---|---|
| match | (1, 1) | `q_i == t_j` (codes already I → L) | 0 |
| substitution | (1, 1) | `q_i != t_j` | 1 |
| swap | (2, 2) | `q_i q_{i+1} == t_{j+1} t_j`, `q_i != q_{i+1}` | 1 |
| block | (1,2), (2,1), (2,2) | the two segments are the two sides of one `BLOCKS` pair, any order | 1 |
| near_block | (1,2), (2,1) | same with `NEAR_BLOCKS` | 1 |

No operation may consume a target `*` (`STOP`) or `X` (`AMBIG`). A query `NOMATCH` residue can
only be consumed by a substitution. The band: every DP cell `(i, j)` satisfies
`|(j − i) − d0| ≤ 2` for the candidate diagonal `d0`.

```python
@dataclass(frozen=True, slots=True)
class Alignment:
    cost: int; t_start: int; t_end: int; diff_types: tuple[DiffType, ...]   # sorted, with repeats
def verify(query: np.ndarray, frame: np.ndarray, d0: int, k: int,
           tables: BlockTables) -> list[Alignment]: ...
```

`verify` returns every alignment window of minimal cost ≤ `k` in the band, one per distinct
`(t_start, t_end)`, after **window de-duplication**: windows that overlap on this frame are one
locus; keep the one with the lowest cost, then the lowest diff-type precedence
(`substitution < swap < block < near_block`, M3-05 `DIFF_TYPES`), then the smallest `t_start`.
Within one window, ties between paths of equal cost are broken by the same precedence (the
reported `diff_types` is the lexicographically smallest multiset under that order). Returns `[]`
when the minimum exceeds `k`.

## Examples (strings shown before encoding; the target is a frame segment, `d0` = true diagonal)

| # | query (L, k) | target window | Expected |
|---|---|---|---|
| 1 | `AKLMNPQRSTVW` (12, 2) | `AKLMPNQRSTVW` | cost 1, `(swap,)` |
| 2 | `PEPTNDEKAR` (10, 1) | `PEPTGGDEKAR` | cost 1, `(block,)`, window length 11 |
| 3 | `VTEGTAPKLR` (10, 1) | `VTEASAPKLR` | cost 1, `(block,)` — not 2 substitutions |
| 4 | `SWAPKELLTR` (10, 1) | `SADAPKELLTR` | cost 1, `(near_block,)` |
| 5 | `MKRLLAEQTR` (10, 1) | `MKGVLLAEQTR` | cost 1, `(near_block,)` |
| 6 | `PEPTIDEKAR` (10, 1) | `PEPTLDEKAR` | cost 0, `()` (I = L) |
| 7 | `PEPTIDEKAR` (10, 1) | `PAPTIDEKGR` | `[]` (cost 2 > k) |
| 8 | `AKLMNPQRSTVW` (12, 2) | `AKLMPNQRSTVF` | cost 2, `(substitution, swap)` |
| 9 | `PEPTIDEKAR` (10, 1) | `PEPT*DEKAR` | `[]` (a stop is never consumed) |
| 10 | `AGLLKDEPTR` (10, 1) | `QLLKDEPTR` | cost 1, `(block,)` ({A,G} ↔ {Q}, window length 9) |

## Edge cases

- Two length-changing blocks in one direction move the end diagonal by 2: still inside the ±2 band.
- Block at the very first or last query residue (free target ends): handled (examples 2, 10).
- `AA` vs `AA` is a match, never a swap. `GG` vs `N` at the query end with the window at the frame
  end: the window may not run past the frame end.
- `k = 0` reduces to exact matching (used as a cross-check against M3-04).

## Acceptance criteria

- [ ] `tests/unit/test_tolerant_dp.py`: examples 1–10 and the edge cases.
- [ ] Property (Hypothesis): an **independent brute-force oracle** written by the test author
  (`tests/property/oracle_tolerant.py`: exhaustive recursion over operation sequences, no band,
  windows of length `L − 2k … L + 2k`) gives the same minimal cost on random (query, window) pairs
  built by planting 0–2 operations of each type into random sequences; every difference type is
  exercised (Hypothesis `@example`s per type).
- [ ] Property: cost is symmetric under reversing both query and target.
- [ ] Benchmark: 10⁶ `verify` calls on 12-aa queries ≤ 5 s on one core (numba `njit(cache=True)`;
  no Python objects in the kernel; traceback encodes diff types as small ints).

## Questions for the lead

1. Should a substitution be allowed against target `X` (ambiguous codon)? **Proposed default:** no;
   `X` never aligns (spec: "ambiguous codons → `X`, which never matches"), so gaps cannot create
   tolerant hits.
2. Tie-break precedence when two explanations of one window have equal cost. **Proposed
   default:** the precedence above (simplest explanation first).

## Lead decisions after code review (PR #39, changes requested)
1. `verify` raises `ValueError` for queries with L ≤ 2·BAND (= 4 aa): the ±2 band can then hold two
   disjoint best windows. Documented in the docstring. M3-08 (`min_side_aa: 4`) must not call `verify`
   with sides of ≤ 4 aa — recorded for M3-08 (use its own anchored check for such sides).
2. Tests for the stop/X/NOMATCH guards: a stop or X in each target slot and `#` in each query slot of every
   block shape (1↔1, 1↔2, 2↔1, 2↔2) is never consumed (e.g. `P` vs `*G`, `GS` vs `*S` are not blocks).
3. The kernel runs without numba: index values are cast to Python int (`int(q[i-1])`), and a test runs a
   representative subset with `NUMBA_DISABLE_JIT=1` in a subprocess, agreeing with the jitted results.
4. `verify` normalises inputs (C-contiguous uint8 via `np.ascontiguousarray(..., dtype=np.uint8)`) so other
   dtypes/strides never trigger a recompile. Test: int64 and strided inputs give the same result.
5. Test the `t_end` tie-break (equal cost, equal t_start → smallest t_end).
6. Deferred to the consuming issues (not in this PR): a public njit batch entry point for M3-07's `prange`,
   and anchored window ends for M3-08; each issue specifies its interface.
