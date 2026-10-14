# M3-02: Six-frame translation and coordinate maps

Milestone M3 · depends on: M1-01, M1-02 · creates `src/giap/mapping/sixframe.py` · pure Python +
numba · estimate ≈ 320 code lines

## Spec sections implemented

- §7.2.1 "**Translation.** Translate all six frames of every contig with table 11 (`*` at TAA, TAG,
  TGA; ambiguous codons → `X`, which never matches)."
- "**Alphabet.** Apply I → L in the index and in queries, so the alphabet has 19 letters plus `*`
  and `X`."
- "**Circular replicons.** Append the first 120 nt before translation, so peptides up to 40 aa can
  span the origin. Map hits back modulo the replicon length and deduplicate them."
- "**Storage.** Arrays of `uint8` per frame, with a function that maps (contig, frame, amino-acid
  index) to a genomic interval and strand."
- §6.5 origin-crossing intervals stored unwrapped; linear contigs never wrap. §7.2.9 unit tests:
  "table-11 translation; frame ↔ coordinate maps on both strands; I = L; origin-spanning peptides
  on circular contigs; no wrap on linear contigs".

Decoys are M3-03; this issue builds target frames only.

## Public interface (`src/giap/mapping/sixframe.py`)

```python
ALPHABET: str = "ACDEFGHKLMNPQRSTVWY*X"   # codes 0..18 residues (L includes I), 19 '*', 20 'X'
STOP: int = 19; AMBIG: int = 20; NOMATCH: int = 21  # NOMATCH: query letters B J O U X Z (never match)
WRAP_NT: int = 120

def encode_query(peptide: str) -> np.ndarray: ...           # uint8; I → L; others → NOMATCH
def encode_target(aa: str) -> np.ndarray: ...
def decode(codes: np.ndarray) -> str: ...
def translate11(nt: bytes) -> np.ndarray: ...               # numba; codons with non-ACGT → AMBIG

@dataclass(frozen=True, slots=True)
class Frame:
    contig: str; strand: Literal["+", "-"]; offset: int     # offset 0,1,2 from the strand's 5' end
    aa: np.ndarray                                          # uint8 codes
    ext_nt: int                                             # length of the translated (extended) sequence

@dataclass(frozen=True)
class SixFrameIndex:
    replicons: Mapping[str, Replicon]
    frames: tuple[Frame, ...]       # per contig in FASTA order: +0,+1,+2,-0,-1,-2
    def to_interval(self, frame_id: int, aa_start: int, aa_end: int) -> Interval: ...
    def frame_label(self, iv: Interval) -> int: ...          # ±1..3 (TSV `frame`)

def build_sixframe(genome: Mapping[str, str], replicons: Mapping[str, Replicon],
                   wrap_nt: int = WRAP_NT) -> SixFrameIndex: ...
def dedupe_loci(ivs: Iterable[Interval], replicons) -> list[Interval]: ...
```

**Coordinate map.** Let `E = S + S[:w]` (circular, `w = min(wrap_nt, L)`) or `E = S` (linear),
`n = len(E)`. Plus frame offset `f`, amino acids `[a, b)` → nt `[f + 3a, f + 3b)` of `E`. Minus
frame offset `f` translates `revcomp(E)`; amino acids `[a, b)` → `E` nt `[n − f − 3b, n − f − 3a)`.
Then reduce: `start' = start mod L`, `end' = start' + (end − start)`; on a linear contig
`end' ≤ L` always holds. The result is `Interval(contig, start', end', strand)` (unwrapped when
`end' > L`).

**Frame label** (TSV column `frame`, also used by Modules 5–6): plus `+(1 + start % 3)`; minus
`−(1 + ((L − end) % L) % 3)`. It equals `offset + 1` with the sign of the strand for every
interval fully inside `[0, L)`.

## Examples

Genome `g = ATGGCTAAAATTCTGTGGTAACGTNNA` (27 nt, linear contig `c1`):

| # | Call | Expected |
|---|---|---|
| 1 | `decode` of frame `+0` | `MAKLLW*RX` (I → L; `NNA` → `X`) |
| 2 | frames `+1`, `+2` | `WLKFCGNX`, `G*NSVVTX` (the last codon contains `N` → `X`) |
| 3 | frames `−0`, `−1`, `−2` | `XTLPQNFSH`, `XRYHRLLA`, `XVTTEF*P` |
| 4 | `to_interval(+0, 1, 4)` | `Interval("c1", 3, 12, "+")`, `frame_label = +1` |
| 5 | `to_interval(−0, 1, 2)` (`T`) | `Interval("c1", 21, 24, "-")`, `frame_label = −1` |
| 6 | `to_interval(−1, 2, 4)` (`YH`) | `Interval("c1", 14, 20, "-")`, `frame_label = −2` |
| 7 | circular `r` = `AAAGAACTGTAAGGCTCCATGCCG` (24 nt), test `wrap_nt=12`: frame `+0` | `KEL*GSMPKEL*`; `MPKEL` = aa `[6,11)` → `Interval("r",18,33,"+")`, crosses origin; TSV start 19, end 33 |
| 8 | same, `KEL` found at aa 0 and aa 8 | `dedupe_loci` → one locus `Interval("r",0,9,"+")` |
| 9 | same sequence declared linear | frame `+0` is `KEL*GSMP` (no wrap); no locus has `end > 24` |
| 10 | `encode_query("PEPTIDEK")` vs `encode_query("PEPTLDEK")` | identical arrays; `encode_query("SEUK")` has `NOMATCH` at index 2 |

## Edge cases

- Contig length not divisible by 3: frames have `⌊(n − f)/3⌋` amino acids; a trailing partial
  codon is dropped.
- Replicon shorter than `wrap_nt` (circular): append the whole replicon once (`w = L`).
- Lower-case and soft-masked input: upper-cased before translation (preflight normalises anyway).
- IUPAC codes other than ACGT in any codon position → `X` (`AMBIG`); `X` and `*` never match
  (enforced by matchers: no query code equals `STOP`/`AMBIG`).
- Minus-strand origin-crossing locus: example 7 mirrored by property test (reverse complement).

## Acceptance criteria

- [ ] `tests/unit/test_sixframe.py`: examples 1–10; all 64 codons of table 11 (against
  `giapsim.codons.CODE_11` and Biopython `Seq.translate(table=11)` as an independent oracle).
- [ ] Property (Hypothesis, random genomes ≤ 5 kb, both topologies): for every frame and every
  window `[a, b)` with `b − a ≤ 40`, the nucleotides of `to_interval(...)` (read with
  `split_at_origin`, reverse-complemented on minus) translate to the same codes.
- [ ] Property: reverse-complementing the genome swaps strands: a plus locus `[s, e)` becomes a
  minus locus `[L − e, L − s)` with the same amino acids ("reverse complement swaps the strand").
- [ ] Property: on circular replicons every window of ≤ 120 nt that crosses the origin is
  reachable exactly once after `dedupe_loci`.
- [ ] Benchmark (`benchmarks/test_sixframe.py`): `build_sixframe` on a 5-Mb random genome
  ≤ 3 s and ≤ 100 MB peak on one core (numba translation kernel; decision 20).

## Questions for the lead

1. The `frame` column convention is not fixed by the spec. **Proposed default:** the rule above
   (`±1..3`, phase relative to each strand's 5′ end of the replicon), documented in the schema and
   reused by Modules 5–6.
2. `wrap_nt` larger than `3 × max_len` would be needed for frameshift splits spanning the origin
   (40 aa + 1 nt = 121 nt). **Proposed default:** keep 120 (spec value); a frameshifted peptide
   whose split point is within 1 nt of the 120-nt edge is missed and listed as a known limitation.

## Lead decisions after code review (M3-02b, same PR)
- NOMATCH decodes to `#` (outside the protein alphabet), never `X`, so a query can never match a target X.
- `Frame.offset` is relative to the extended (wrapped) sequence; M3-04 uses `frame_label(iv)`, never
  `offset`; the docstring says so.
- `to_interval` rejects an interval longer than the replicon (`3*(b-a) > L`) with ValueError.
- The wall-clock/memory budget moves to `benchmarks/test_sixframe.py` (as the issue says), measuring peak
  memory portably (`ru_maxrss` is KiB on Linux, bytes on macOS — or use tracemalloc); it is not part of
  the per-commit suite.
- Missing tests: lowercase queries, the `i` → L fold, `errors="replace"` keeping coordinates for
  non-ASCII bytes, `dedupe_loci` with unwrapped starts ≥ L, the decode clamp.
