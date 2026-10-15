# M2-13: Promoter windows, Promotech wrapper and parser, per-gene promoter features, `promoters.bed`

Milestone M2 · depends on: M2-04, M2-12 (TUs) · estimate ≈ 380 code lines (tests not counted) ·
capabilities: none (recording the Promotech fixture needs `needs_conda`, `needs_promotech_models`)

## Spec sections implemented

- §7.1.4: "**Targeted Promotech.** For each start codon S of every CDS and sORF (and, in the second
  round, every novel-ORF candidate start), score the 40-nt windows on the gene's strand whose 3′ end
  lies at S + k for k ∈ [−200, +20]: 221 windows per start. Deduplicate identical windows, chunk them
  into files of 100,000 windows and run the chunks in parallel." Per-gene features: "`promoter_score`:
  the maximum score over the gene's own windows; `promoter_score_tu`: the same at the TU's first gene;
  the count of windows ≥ 0.5 between the TU start and the gene." "**Promoter between two adjacent
  genes:** the maximum score over windows of the downstream gene whose 3′ end lies after the upstream
  gene's end." "**`promoters.bed`** holds windows with score ≥ 0.5, merged when overlapping; the score is the maximum."
- §7.1.1 Promotech gotcha: "`-s` expects exactly 40-nt sequences. Window headers encode
  `contig|strand|end3p` (avoid characters the tool may mangle; test it). Genome-wide mode … is opt-in
  and preflight warns." Command: "targeted: `python promotech.py -s -m RF-HOT -f windows.fa -o out/`".
- §5.2 `promoters: {mode: targeted, window_nt: 40, window_3p_range: [-200, 20], promotech_yes: 0.5}`.
- Module 1 tests: "window generation property tests: count, strand, 3′-end range, origin wrap on
  circular replicons, no wrap on linear contigs".
- §6.5: BED is 0-based half-open; origin-crossing records split `_part1`/`_part2`; ADR-0006 BED scores
  are integers 0–1000.

## Public interface

`src/giap/annotate/promoters.py`, `src/giap/annotate/promotech_cmd.py`; steps `promoter_windows`,
`promoter_layer`.
```python
@dataclass(frozen=True, slots=True)
class Window:
    contig: str; strand: str; end3p: int     # 1-based genome coordinate of the window's 3'-most base
    seq: str                                  # 40 nt, 5'→3' on the gene strand
    owner: str; k: int                        # locus tag, offset

def gene_windows(gene: Gene, replicon: Replicon, seq: str, k_range=(-200, 20), width=40) -> list[Window]
def dedupe_and_chunk(windows: Iterable[Window], chunk_size=100_000) -> tuple[list[list[tuple[str, str]]], pd.DataFrame]
    # chunks of (header, seq); map frame: header, contig, strand, end3p, owner, k
def write_chunks(outdir: Path, chunks) -> list[Path]        # windows_000.fa, windows_001.fa, …
def promotech_command(chunk: Path, outdir: Path, promotech_home: Path) -> list[str]
def parse_promotech(path: Path) -> dict[str, float]          # header → score
def promoter_layer(window_map: pd.DataFrame, scores: Mapping[str, float]) -> pd.DataFrame  # contig, strand, end3p, score (unique)
def gene_promoter_features(genes, tus, window_scores, yes=0.5) -> pd.DataFrame
    # locus_tag, promoter_score, promoter_score_tu, n_promoter_windows_tu
def promoter_between(up: Gene, down: Gene, window_scores, replicon) -> float | None
def promoters_bed(window_scores, replicons, yes=0.5) -> list[BedRecord]
```
Window geometry (S = first base of the start codon, 0-based; E = gene interval end):
- `+`: window = genome `[S + k − 40, S + k)`; `end3p = S + k` (1-based of base `S+k−1`). k = 0 ends
  right before the start codon (Q1).
- `−`: window = revcomp of genome `[E − k, E − k + 40)`; `end3p = E − k + 1`.
- Circular replicons wrap; on linear contigs windows that leave `[0, length)` are dropped.
- Headers are `w<index>` over unique sequences in first-seen order (genes in §6.5 order, k ascending); the
  map frame keeps every (contig, strand, end3p, owner, k) occurrence (Q1).
- `n_promoter_windows_tu`: windows with score ≥ `yes` on the TU strand whose 3′ end lies between the
  TU's first start codon − 200 and the gene's own start codon + 20, in transcription direction.
- `promoters_bed`: windows ≥ `yes` merged per contig and strand when overlapping; BED6 name
  `promoter_<n>` (numbered §6.5), score `round(1000 × max)`, strand kept.

## Examples (chr1 circular, 100,000 bp)

| # | Input | Expected |
|---|---|---|
| 1 | gene `(1000, 1300, +)` | 221 windows; k = 0 → genome `[960, 1000)`, `end3p = 1000`; k = 20 → `[980, 1020)`, `end3p = 1020`; k = −200 → `[760, 800)` |
| 2 | gene `(2000, 2900, −)` | k = 0 → revcomp(genome `[2900, 2940)`), `end3p = 2901`; k = 20 → `[2880, 2920)`; k = −200 → `[3100, 3140)` |
| 3 | gene `(100, 400, +)` circular | k = −200 → `[−140, −100)` wrapped = genome `[99860, 99900)`, `end3p = 99900`; 221 windows |
| 4 | same gene on a linear contig | windows need `S + k − 40 ≥ 0` → k ∈ [−60, 20] → **81** windows |
| 5 | two genes whose windows share 30 identical 40-mers | 30 unique headers fewer; map frame keeps both owners |
| 6 | 250,001 unique windows | 3 chunk files with 100,000 / 100,000 / 50,001 records |
| 7 | scores: gene G windows max 0.83; TU first gene F max 0.97 | `promoter_score(G)=0.83`, `promoter_score_tu(G)=0.97` |
| 8 | up `(100,400,+)`, down `(420,700,+)`; down's windows with `end3p` > 400 scores {0.2, 0.71} | `promoter_between = 0.71`; none after 400 → `None` |
| 9 | windows ≥ 0.5 `[960,1000)` 0.6 and `[970,1010)` 0.9 on + | one BED record `chr1 960 1010 promoter_1 900 +` |
| 10 | Promotech header round trip on `w0 … w99999` | `parse_promotech` returns the same keys (no mangling) |

## Edge cases

- Window sequences containing `N` or other IUPAC codes are still scored (Promotech decides); a
  window shorter than 40 nt is never written (assert).
- Genes shorter than 20 nt downstream of S (sORFs): windows with k up to +20 still extend into or
  beyond the gene; allowed.
- Genome-wide mode (`promoters.mode: genome`) is out of scope here (opt-in, M2 only rejects it with
  a clear message until implemented, Q2).
- BED records crossing the origin → two records `promoter_<n>_part1/_part2`.

## Acceptance criteria

- [ ] `tests/unit/test_promoter_windows.py`: examples 1–6, 9; `tests/property/test_windows.py`:
  count = 221 on circular, every window 40 nt, strand = gene strand, `end3p − S` ∈ [−200, 20]
  (strand-adjusted), `revcomp` symmetry between a + gene and its mirrored − gene.
- [ ] `tests/unit/test_promoter_features.py`: examples 7–8, `n_promoter_windows_tu` boundaries.
- [ ] `tests/contract/test_promotech_contract.py` on the recorded sim100k chunk (example 10);
  `pending_sample` until recorded with the M2-10 harness (`needs_promotech_models`).

## Questions for the lead

1. §7.1.1 asks headers to encode `contig|strand|end3p`, but deduplication maps one sequence to
   several positions. **Proposed default:** `w<index>` headers + the Parquet map (no mangling risk);
   and "3′ end at S + k" means the window's last base is `S + k − 1` (k = 0 ends before the codon).
2. Genome-wide Promotech mode. **Proposed default:** not in M2 (opt-in, needs ≥ 24 GB and 5.8 GB per
   strand); config value `genome` exits 1 with "not implemented in this version" until an issue adds it.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (minus-strand `end3p = E−k+1`, BED windows `[p−1, p+39)`,
   inclusive bounds mirrored on the minus strand, unwrapped origin-crossing record split by the writer).
2. Errors are `ValueError` (never `AssertionError`): short windows, a Promotech header other than
   `CHROM SEQ PRED`.
3. (Corrected — GIAP_SPEC.md §7.1, `promoters.bed`): windows with score ≥ 0.5 that overlap are merged into
   one record with the maximum score, as in example 9; windows that only touch (end == start) are not merged.
4. Q2 (revised): `annotation.promoters.mode: genome` is a valid config — `giap preflight` accepts it and
   checks its resource needs (M1-11/M1-12 tests: exit 0, exit 3 with --strict). It is rejected at run time:
   `giap run` exits 1 before launching anything ("annotation.promoters.mode: genome is not implemented"),
   and the `promoter_windows` step raises the same error if reached. Not at config load.

## Lead decisions after code review (PR #76, changes requested)
5. ANNOTATE wiring (main, after #65): the `promoter_windows` call in annotate.nf passes `replicons` (the step
   needs it for circular wraps) and collects the chunks with
   `scatterSpec('windows', 'windows/windows_*.fa', …, true)` + `inDir`; the step writes chunks into the
   `windows/` directory. The M2-index row is updated. A test runs the step end to end through the
   registry, and the step-registry fixture (workflow/tests/data/annotate/step_registry.json) is regenerated.
6. Tests kill the surviving mutants: BED merge across the origin (and its max score); the `gap % n` wrap;
   circular normalisation in `promoter_between`; the `max(a, 0)` clamp; the `lo1` wrap (a window ending at
   end3p == n); the pseudogene filter in the step; the header check with a TAB-separated foreign header.
7. `gene_promoter_features(…, replicons)` requires `replicons` (no silent linear default).
8. Review finding 5 (`giap _step` never loads step modules) is fixed by the run-wiring PR #78.
