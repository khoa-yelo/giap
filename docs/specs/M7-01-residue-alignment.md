# M7-01: Residue-level alignment of every accepted placement

Milestone M7 · track V (variants) · **pass 1** · depends on: merged M3 mapping (M3-04 exact, M3-06
tolerant DP, M3-08 frameshift splits, M3-14 loci, M3-16 outputs) · estimate ≈ 240 code lines
(+ ≈ 400 test lines) · not statistical (exact recomputation), but every Module 5 number rests on it

## Why this issue exists

§7.5.2 step 1 assumes "the mapper's traceback gives a residue-level correspondence between each
peptide and the genome translation". The merged mapper keeps **no traceback**: the M3-06 DP
(`giap.mapping.tolerant_dp`) encodes the optimum as one additive integer key "so the DP keeps the
exact optimum without a traceback", and `mapping/peptide_hits.tsv` keeps only `n_diff` and
`diff_types`. The frameshift route keeps the two genomic blocks (`blocks = "a-b;c-d"`) in the work
table `hits_tolerant.parquet`, which `peptide_hits.tsv` drops. This issue recomputes the residue
correspondence for each accepted locus, outside the mapper, without changing M3 outputs.

## Spec sections implemented

- §7.5.2 step 1: "Residue alignment. The mapper's traceback gives a residue-level correspondence …
  Mass-equivalent block exchanges are recorded at the first codon of the block, with
  `alt_type = block` (for example `GG>N`). Frameshift splits are recorded as `alt_type = frameshift`
  events at the junction; they appear in `variants.tsv` but not in the VCF in version 1."
- §7.2.4 (the cost-1 operations: substitution, adjacent swap, the `BLOCKS` and `NEAR_BLOCKS`
  exchanges, ±1-nt frameshift split; "I/L is free"); §7.2.1 (table 11, `X` never matches, the
  circular 120-nt append and "map hits back modulo the replicon length").
- §6.5 (0-based half-open, unwrapped coordinates on circular replicons; minus strand).

## Public interface

`src/giap/variants/align.py` (pure; no file I/O)
```python
OpKind = Literal["match", "il", "substitution", "swap", "block", "near_block"]

@dataclass(frozen=True, slots=True)
class ResidueOp:
    kind: OpKind
    pep: tuple[int, int]             # peptide residues [i, j), 0-based
    codons: tuple[int, ...]          # genome codons consumed: lowest coordinate, 0-based, unwrapped
    ref: str                         # genome residues as translated (I and L kept distinct)
    alt: str                         # peptide residues as given (I and L kept distinct)

@dataclass(frozen=True, slots=True)
class FrameshiftEvent:
    pep_junction: int                # index of the first suffix residue
    last_prefix_codon: int           # lowest coordinate of the prefix's last codon
    first_suffix_codon: int
    shift: Literal[-1, 1]            # δ of M3-08

@dataclass(frozen=True, slots=True)
class HitAlignment:
    peptide: str                     # as in peptide_hits.tsv (I/L as observed)
    contig: str
    strand: Literal["+", "-"]
    frame: int
    ops: tuple[ResidueOp, ...]       # peptide order (N → C)
    frameshift: FrameshiftEvent | None

def translate_locus(seq: str, rep: Replicon, start: int, end: int, strand: str) -> str: ...
def align_hit(peptide: str, contig: str, start: int, end: int, strand: str, n_diff: int,
              blocks: str | None, seq: str, rep: Replicon) -> HitAlignment: ...
def align_hits(peptide_hits: pd.DataFrame, hits_tolerant: pd.DataFrame | None,
               genome: Mapping[str, str], replicons: Mapping[str, Replicon]) -> list[HitAlignment]: ...
```

## Rules

1. **Input rows.** One alignment per row of `mapping/peptide_hits.tsv` (1-based closed in the
   file; the step converts to 0-based before calling, §6.5). Frameshift rows are those whose
   `diff_types` contains `frameshift`; their blocks come from the matching `hits_tolerant` row
   (same peptide, search `target`, contig, start, end, strand). A frameshift row without a
   matching work row raises, naming the peptide.
2. **Genome residues.** `translate_locus` translates `[start, end)` on `strand` with table 11
   (M3-01 `translate11_text`), reverse-complementing on `−`, reading across the origin on circular
   replicons (`end > length`), raising on a linear contig when `end > length`. **No I → L
   folding** here: the REF residue is what the genome encodes.
3. **Exact route** (`route = exact`, `n_diff = 0`): residue-by-residue. Identical letters → `match`;
   letters that differ only as I/L → `il`. Any other difference raises (`ValueError`: the hit is not
   exact), because the exact route folded only I/L.
4. **Tolerant route, no frameshift.** A banded DP with traceback over the M3-06 alphabet (I folded
   into L **for the cost only**), the same operations at cost 1 each and the same band (±2). It
   must reproduce `n_diff`; a different optimum raises `ValueError` naming the peptide and both
   costs. Among co-optimal alignments the M3-06 key decides (cost, then the lexicographically
   smallest diff-type multiset under `DIFF_TYPES` precedence); remaining ties put each difference
   at the **leftmost** peptide position (decision 2). Matched positions are then labelled `match`
   or `il` by comparing the unfolded letters, as in rule 3.
5. **Block ops.** A `block` or `near_block` op consumes 1 or 2 genome codons and 1 or 2 peptide
   residues; `codons` lists all consumed codons in peptide order, so `codons[0]` is "the first
   codon of the block" (§7.5.2). `ref` = the genome residues, `alt` = the peptide residues
   (example 5 gives `ref = "N"`, `alt = "GG"`).
6. **Swap ops.** `XY ↔ YX` consumes two codons; `ref = "XY"` (genome), `alt = "YX"` (peptide).
7. **Frameshift route.** The prefix is aligned on block 1 and the suffix on block 2 (each with
   rule 4 on its own frame; `n_diff − 1` is shared between them). `FrameshiftEvent` records the
   junction. Both sides' ops are emitted. On `−` the prefix block is the higher-coordinate block
   (the 5′ end of a minus-strand hit is `end`, M3-14 note).
8. **Codon coordinates.** `codons` are the lowest coordinate of each codon, unwrapped (they may
   exceed the replicon length when the hit is stored unwrapped). On `+`, residue i of a hit
   starting at s sits at s + 3i; on `−`, residue i of a hit ending at e sits at e − 3(i + 1).
9. **Purity and speed.** No file I/O; numba is allowed for the DP but not required. The DP only
   runs on tolerant rows (exact rows are a zip).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| empty `peptide_hits` | `[]` |
| exact row whose letters differ beyond I/L | `ValueError` naming the peptide (mapper/aligner disagreement) |
| tolerant row whose recomputed cost ≠ `n_diff` | `ValueError` naming the peptide, `n_diff` and the recomputed cost |
| genome residue `X` (ambiguous codon) under a peptide residue | impossible for an accepted hit (`X` never matches, §7.2.1); raises if met |
| codon crossing the origin of a circular replicon | normal op; its `codons` entry is `length − 1` or `length − 2` (unwrapped coordinates continue past `length`) |
| linear contig, hit beyond the end | `ValueError` (cannot occur for an M3 hit) |
| frameshift row with no `blocks` | `ValueError` naming the peptide |
| peptide with several accepted loci | one `HitAlignment` per locus row |
| decoy rows | never passed (`peptide_hits.tsv` holds targets only) |

## Examples

Contig `c1`, circular, 1,000 nt unless stated.

| # | Input | Expected |
|---|---|---|
| 1 | exact `+` hit [100, 124) of `PEPTIDEK`, genome translates `PEPTLDEK` | 8 ops; op 4 `il` (ref `L`, alt `I`), codon 112; the others `match`; codons 100, 103, …, 121 |
| 2 | exact `−` hit [200, 224) of an 8-mer | op 0 codon **221**, op 7 codon **200** (rule 8) |
| 3 | tolerant `+` hit, `AGSDEFKLLR` on genome `AGSNEFKLLR`, `n_diff` 1 | op 3 `substitution`, ref `N`, alt `D`, codon start + 9; all others `match` |
| 4 | tolerant, `AGDSEFKLLR` on `AGSDEFKLLR`, `n_diff` 1 | one `swap` op on peptide [2, 4), ref `SD`, alt `DS`, two codons |
| 5 | tolerant, `AGGGEFKLLR` (10 aa) on the 9-codon locus `ANGEFKLLR`, `n_diff` 1 | one `block` op, peptide [1, 3), `codons` = (start + 3,), ref `N`, alt `GG` |
| 6 | tolerant, `ANGEFKLLRK` on `AGGGEFKLLRK`-encoding locus, `n_diff` 1 | `block` op, peptide [1, 2), two codons, ref `GG`, alt `N` (the spec's `GG>N`) |
| 7 | tolerant, `WEFKLLRAGK` on `ADEFKLLRAGK`, `n_diff` 1 | `near_block` op, ref `AD`, alt `W` |
| 8 | frameshift `+` hit, prefix `MKLLAE` on frame 0 at [300, 318), suffix `GSTRK` on frame +1 from 319, `blocks` "300-318;319-334", `n_diff` 1 | 11 `match` ops; `FrameshiftEvent(pep_junction=6, last_prefix_codon=315, first_suffix_codon=319, shift=+1)` |
| 9 | exact `+` hit [995, 1019) on the circular 1,000-nt contig | op 1's codon is **998** and spans the origin ([998, 1001)); op 2's codon is 1001 (unwrapped) |
| 10 | tolerant row with `n_diff` 1 whose locus actually aligns at cost 2 | `ValueError` naming the peptide, "n_diff 1, recomputed 2" |

## Acceptance criteria

- [ ] `tests/unit/variants/test_align.py`: examples 1–10 and every degenerate row.
- [ ] Property (Hypothesis): for random genomes ≤ 5 kb and peptides cut from them with ≤ 2 planted
  differences of every §7.2.4 type, on both strands and across the origin, the aligner's cost
  equals the M3-06/M3-08 cost of the same locus, and replaying the ops on the genome translation
  reproduces the peptide (I/L as given).
- [ ] Property: mirroring the genome (reverse complement) and the hit maps each op to the same
  peptide positions with mirrored codon coordinates.
- [ ] Conformance on the merged sim100k fixture: every row of the real `peptide_hits.tsv` aligns
  without error (the recomputation agrees with the mapper on every accepted hit).
- [ ] Benchmark in `benchmarks/`: 100,000 tolerant hits aligned in < 30 s on 1 CPU
  (`benchmarks/baseline.json` entry with node, CPU, commit).

## Calibration and attacks owed

None by itself: correctness is exact and tested by property. The F1 gate (M7-13) is where an
alignment bug would show.

## Hardening owed (pass 2)

- A recorded-hit contract test against the 50-kb golden fixture's `hits_tolerant.parquet`
  (M7-15).
- If the M3 mapper ever stores a traceback, replace the recomputation and keep this property test
  as the equivalence check.

## Questions for the lead

1. Recompute the traceback in Module 5 (proposed), or extend the M3-06/M3-08 kernels to emit one?
   **Proposed default:** recompute. The M3 kernels are calibrated and their key layout is
   load-bearing; a second, slower DP over accepted hits only (≤ 2 differences, ≤ 40 aa) is cheap
   and leaves M3 untouched. The cost equality is asserted on every hit.
2. Placement of co-optimal differences (same cost and same diff-type multiset, different
   positions, e.g. a substitution that could sit at either of two equal residues): §7.2.4 does not
   say. **Proposed default:** leftmost peptide position. It only affects which codon carries the
   ALT in rare repeats; the alternative codon still shows REF support from other peptides.
3. `hits_tolerant.parquet` becomes an input of the `variants` step (for frameshift blocks). Is a
   work-table dependency acceptable, or should `peptide_hits.tsv` gain a `blocks` column?
   **Proposed default:** read the work table; changing the published M3 schema is a breaking column
   change (§9.2 versioning) for a field only Module 5 needs.

## Lead decisions after the test step (2026-10-21)
1. The test-author's resolutions at cfc4fab are accepted.
   - **Co-optimal ties:** choose the lexicographically smallest tuple of difference start positions, in peptide order. This refines "leftmost" (M7-01 Q2).
   - **`frame`:** follows the `peptide_hits.tsv` convention.
   - **I/L folding:** applies inside swaps and blocks.
   - **Block codons:** listed in peptide order.
   - **Frameshift `shift`:** use the per-strand definitions in the tests.
   - **Tolerant route:** an end-to-end alignment within a ±2 band.
2. The benchmark (100k hits in under 30 s) and its `baseline.json` entry are owed by the implementer.
