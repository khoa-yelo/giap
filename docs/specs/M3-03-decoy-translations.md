# M3-03: Decoy translations and the `sixframe_index` work file

Milestone M3 · depends on: M3-02 · creates `src/giap/mapping/decoys.py`, `save_index`/`load_index`
in `src/giap/mapping/sixframe.py` · pure Python + numba · estimate ≈ 220 code lines ·
**statistical** (decoy construction): section "Genome decoys" of
`docs/methods/placement-fdr.md`, signed off by the statistician (§14.1) before merge

## Spec sections implemented

- §7.2.1 "**Decoy translations:** *reverse* (database input): each frame's amino-acid string
  reversed; *shuffle* (de novo input): residues shuffled within each stop-free segment, in windows
  of 50 aa, with a seeded RNG. Decoy hits are only ever counted, never placed."
- §8.2 "**Decoys mimic targets.** They keep composition, size and structure".
- §8.5 attack 5: "De novo errors that reverse a stretch of sequence → placement FDR with shuffled
  decoys (a reversed decoy would inflate D)" (the reason de novo uses shuffle).
- §7.6.4 (M7 consumer): "Each decoy hit is assigned the locus class of its mirrored target position
  (reversed or shuffled windows keep the frame)" → every decoy position must map to a mirrored
  target position.
- §6.2 step `sixframe_index`: "genome → six-frame and decoy translations (`.npz`)". §6.6 seeds.
- Config: `mapping.decoys: {database_input: reverse, denovo_input: shuffle, shuffle_window_aa: 50}`.

## Public interface

```python
DecoyType = Literal["reverse", "shuffle"]

def reverse_frame(aa: np.ndarray) -> np.ndarray: ...
def shuffle_frame(aa: np.ndarray, window: int, rng: np.random.Generator) -> np.ndarray: ...
def build_decoys(index: SixFrameIndex, kind: DecoyType, window: int,
                 rng: np.random.Generator) -> tuple[Frame, ...]: ...   # same order as index.frames
def mirror_aa(kind: DecoyType, n_aa: int, aa_start: int, aa_end: int) -> tuple[int, int]: ...
# reverse: (n_aa − aa_end, n_aa − aa_start); shuffle: identity

def save_index(path: Path, index: SixFrameIndex, decoys: Mapping[DecoyType, tuple[Frame, ...]]) -> None
def load_index(path: Path) -> tuple[SixFrameIndex, dict[DecoyType, tuple[Frame, ...]]]
```

- **Shuffle rule.** `*` (`STOP`) and `X` (`AMBIG`) stay in place. Each maximal run without `*` is
  a segment; within it, the non-`X` residues are cut into consecutive windows of `window` residues
  from the segment start (the last window may be shorter) and each window is permuted with
  `rng.permutation`. Frames are processed in index order, so the stream is deterministic.
- **RNG.** `giap.rng.for_step("sixframe_index", config.seed)`.
- **Which decoys are built.** `reverse` if any query has `source = database`; `shuffle` if any has
  `source = denovo`; both when the input is mixed (each peptide is searched against the decoy of
  its own source; M3-10 reports `decoy_type` per row).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `reverse_frame` of `MKV*AL` | `LA*VKM` |
| 2 | `shuffle_frame` of `MKVAL*GG*` (window 50) | positions 5 and 8 stay `*`; `[0,5)` is a permutation of `MKVAL`; `[6,8)` stays `GG` |
| 3 | `shuffle_frame` of 120 residues without `*`, window 50 | windows `[0,50)`, `[50,100)`, `[100,120)` each keep their multiset |
| 4 | `shuffle_frame` of `AXKL` | `X` stays at index 1; `{A,K,L}` permuted over indices 0, 2, 3 |
| 5 | same frame and seed twice | byte-identical output; another seed differs (for a 1,000-aa frame) |
| 6 | `mirror_aa("reverse", 100, 10, 18)` | `(82, 90)` |
| 7 | `mirror_aa("shuffle", 100, 10, 18)` | `(10, 18)` |
| 8 | `save_index` then `load_index` on the 50-kb simulated genome | frames, decoys and replicons equal |

## Edge cases

- Frame with no stop at all (short contig): one segment.
- Frame of length 0 (contig < 3 nt after the offset): empty arrays survive the round trip.
- The circular extension (`wrap_nt`) is decoyed like the rest of the frame; duplicated decoy
  windows near the origin do not matter because each peptide counts once at its best decoy cost
  (M3-10), and decoy hits are de-duplicated on their mirrored target interval.

## Acceptance criteria

- [ ] `tests/unit/test_decoys.py`: examples 1–8.
- [ ] Property ("decoy translations keep composition", §7.2.9): for random frames, reverse and
  shuffle keep the residue multiset of the frame, the positions of `*` and the per-window
  multisets.
- [ ] Property: `mirror_aa` is an involution for `reverse` and maps inside `[0, n_aa)`.
- [ ] Benchmark: both decoy sets for a 5-Mb genome ≤ 3 s (numba kernel for the shuffle, fixed
  window order; §6.6 determinism).
- [ ] Methods note section "Genome decoys" (what each decoy preserves, why shuffle for de novo,
  attack 5) exists and is signed off.

## Questions for the lead

1. `X` positions: fixed or shuffled? **Proposed default:** fixed (they never match either way;
   keeping them fixed preserves the "never match" geometry).
2. Mixed database + de novo input: one decoy per source (above) or one decoy for all?
   **Proposed default:** per source, with separate FDR rows per `decoy_type`.
