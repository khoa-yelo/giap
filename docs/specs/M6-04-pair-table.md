# M6-04: The pair table — classes, orientation, distances and between-gene evidence

Milestone M6 · track R (statistics) · **pass 1** · depends on: nothing new (merged M2-12 pairs/TUs,
M2-13/14/15 `*_between` helpers, M4-15 gene table) · estimate ≈ 340 code lines (+ ≈ 400 test lines)
· not statistical (geometry and joins), but it feeds every M6 test

## Spec sections implemented

- §7.4.1 in full: the four adjacent classes (OP, NO, CV, DV), the SW decoy ("genes i and i + 2 on
  the same strand, with gene i + 1 on the opposite strand … their distance is the sum of the two
  gaps", G16), `rna_between`, the orientation rule ("For OP, NO and SW pairs, gene *a* is upstream
  in the direction of transcription and Δ = log2 abundance(b) − log2 abundance(a) … For CV and DV
  pairs, *a* is the left gene in genome order") and the pair column list.
- §7.1.3 (adjacency and the wrap pair on circular replicons), §7.1.4 ("Promoter between two adjacent
  genes: the maximum score over windows of the downstream gene whose 3′ end lies after the upstream
  gene's end"), §7.1.5 ("Terminator between two same-strand adjacent genes …").
- §9.3 `tables/pairs.tsv` columns; §9.2 conventions; §6.5 coordinates and IDs.

## Public interface

`src/giap/operons/pairs.py` (pure)
```python
PairClass = Literal["OP", "NO", "CV", "DV", "SW"]

@dataclass(frozen=True, slots=True)
class Evidence:                      # one lookup object built once per run
    windows: pd.DataFrame | None     # promoter_windows.parquet
    terminators: tuple[Terminator, ...] | None
    tx: pd.DataFrame | None
    rna: pd.DataFrame                # C15 rows of RNA type

def adjacent_chain(genes: Sequence[Gene], replicons: Mapping[str, Replicon]
                   ) -> list[tuple[Gene, Gene]]: ...        # all adjacent pairs, wrap included
def sandwich_pairs(genes, replicons) -> list[tuple[Gene, Gene, Gene]]: ...
def classify(a: Gene, b: Gene, uniop_p: float | None, threshold: float) -> PairClass: ...
def pair_frame(genes_tbl: pd.DataFrame, genes: Sequence[Gene], replicons, pairs_raw: pd.DataFrame,
               ev: Evidence, cfg: GiapConfig) -> pd.DataFrame: ...
```

1. **Adjacency.** `giap.annotate.operons.adjacent_pairs` gains `same_strand: bool = True`
   (decision 1; the default keeps the M2-12 behaviour); `adjacent_chain` calls it with `same_strand=False`. One definition of genome order,
   of the wrap pair on circular replicons and of `intergenic_distance` (which already handles the
   origin) is therefore shared with M2-12.
2. **Members.** Every C2 gene of the replicon (cds, sorf, pseudogene), `partial` included
   (decision 2). Genes are ordered by `(start, end, locus_tag)` per contig, C1 contig order.
3. **Class.** Same strand and `uniop_p ≥ annotation.uniop.pair_threshold` → OP; same strand and
   below → NO; `a` on `+` and `b` on `−` → CV; `a` on `−` and `b` on `+` → DV. `uniop_p` comes from
   `pairs_raw` (which holds same-strand pairs only); a same-strand pair missing from `pairs_raw` is
   an error naming the pair, never a silent NO.
4. **SW.** For consecutive triples (i, i+1, i+2) on one contig (wrapping on circular replicons),
   emit a row when i and i+2 share a strand and i+1 is on the other; `intergenic_bp` = gap(i, i+1) +
   gap(i+1, i+2) (G16), `uniop_p` NA. SW rows are extra rows, never part of the chain, and a gene
   can appear in several of them.
5. **Orientation.** OP, NO, SW: on `+`, a is the left gene; on `−`, a is the **right** gene, so a is
   always upstream in transcription. CV and DV: a is the left gene. `strand` = the shared strand, or
   `.` for CV and DV. `delta_log2` = `abundance_log2(b) − abundance_log2(a)`, NA if either is NA.
6. **Evidence** (OP, NO, SW only; decision 4): `promoter_between` = `promoters.promoter_between(up,
   down, ev.windows, replicon)`; `tx_rate_between` = `promoter_calculator.tx_rate_between(...)`;
   `terminator_between` = `terminators.terminator_between(up, down, ev.terminators, replicon,
   yes=0.0)` — the raw maximum confidence (decision 5). A missing layer gives NA for the whole
   column and the existing `layer_missing:<layer>` run flag.
7. **`rna_between`** flag when any C15 feature of type `rrna`, `trna`, `ncrna` or `tmrna` overlaps
   the open interval between the two genes (for SW, the union of the two gaps).
8. **Per-gene columns** are joined from the gene table: `detected_a/b`, `p_detect_a/b`,
   `abundance_log2_a/b`. Genes absent from the gene table (impossible once M4-15 ran) raise.
9. **`pair_id`** = `PAIR_<n>`, n = 1, 2, … over the frame sorted by contig (C1 order), `start_a`,
   `start_b`, then class (decision 3); the sort key never uses `uniop_p`, so IDs are stable under
   the M6-10 threshold variants.
10. **Test columns** (`test`, `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`,
    `null_model`, `n_perm`, `lfdr`) are created here as NA with `test = "none"`; M6-05 and M6-06
    fill them.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| replicon with 1 gene | no pairs; no error |
| circular replicon with 2 genes | one adjacent pair (not two — the wrap pair equals the forward pair); no SW |
| linear replicon | no wrap pair, no SW triple across the ends |
| `partial` gene in a pair | row present, class and distance set, `detected`/`p_detect` as the gene table has them; the pair is excluded from every test by M6-05/06 (eligibility, §7.3) |
| overlapping genes | `intergenic_bp` negative; bin 0 |
| gene crossing the origin | `intergenic_distance` already unwraps (M2-12); coordinates in the TSV are 1-based closed and may exceed the replicon length (§6.5) |
| layer missing | evidence column all NA; no flag row later; run flag `layer_missing:<layer>` |
| CV / DV pair | evidence columns NA; `terminator_between` is **not** called (it raises on mixed strands by design, M2-15) |
| same-strand pair missing from `pairs_raw` | `ValueError` naming gene_a, gene_b |
| RNA gene between two CDS | `rna_between` in `flags`; the pair is otherwise ordinary |

## Examples

Genes on one circular 10 kb contig: G1 `+` [100, 400), G2 `+` [450, 900), G3 `−` [950, 1400),
G4 `−` [1500, 1900), G5 `+` [2000, 2400). `pair_threshold` 0.5.

| # | Input | Expected |
|---|---|---|
| 1 | G1–G2, `uniop_p` 0.8 | class **OP**, a = G1, b = G2, `strand +`, `intergenic_bp` 50 |
| 2 | G1–G2, `uniop_p` 0.2 | class **NO**, same geometry |
| 3 | G2–G3 | class **CV** (→ ←), a = G2 (left), `strand .`, `intergenic_bp` 50, `uniop_p` NA, evidence NA |
| 4 | G3–G4, `uniop_p` 0.9 | class **OP** on `−`: a = **G4** (upstream in transcription), b = G3, `intergenic_bp` 100, Δ = log2(G3) − log2(G4) |
| 5 | G4–G5 | class **DV** (← →), a = G4 (left) |
| 6 | G2, G3, G4 … G2 `+`, G3 `−`, G4 `−` | **no** SW row (i and i+2 are G2 `+` and G4 `−`: different strands) |
| 7 | G3 `−`, G4 `−`… replaced by G3 `+`, so G3 `+`, G4 `−`, G5 `+` | SW row a = G3, b = G5, `intergenic_bp` = 100 + 100 = **200** |
| 8 | wrap: G5 `+` [2000, 2400) and G1 `+` [100, 400) on a circular 10 kb contig, `uniop_p` 0.1 | one NO pair a = G5, b = G1, `intergenic_bp` = 10000 − 2400 + 100 = **7700** |
| 9 | a tRNA feature at [410, 440) | the G1–G2 row gains `rna_between` |
| 10 | abundance_log2 G1 = 18.0, G2 = 16.5 | `delta_log2` = **−1.5**; for the `−`-strand pair of example 4 with G4 = 16.5, G3 = 18.0, `delta_log2` = **+1.5** |
| 11 | `promoter_windows` absent | `promoter_between` NA in every row; run flag `layer_missing:promoters` |

## Acceptance criteria

- [ ] `tests/unit/operons/test_pairs.py` covers examples 1–11 and every degenerate row on
  hand-built frames.
- [ ] Property: every adjacent pair of the chain appears exactly once; the union of OP, NO, CV, DV
  rows equals `adjacent_chain`; SW rows are disjoint from it; classes are mutually exclusive.
- [ ] Property: reversing the strand of every gene and mirroring the coordinates maps OP↔OP, NO↔NO,
  CV↔CV, DV↔DV and SW↔SW. OP, NO and SW keep `gene_a`, `gene_b` and Δ (orientation is
  strand-aware); CV and DV swap `gene_a` and `gene_b` and negate Δ (a is the left gene). See lead
  decision 1 below.
- [ ] Property: `pair_id` assignment is unchanged when `annotation.uniop.pair_threshold` moves from
  0.3 to 0.7 (the M6-10 requirement), while `class` changes.
- [ ] Conformance: on the merged sim100k fixture, the OP pairs reproduce exactly the TU membership
  of `tus.parquet` — a pair is OP iff its two genes share a `tu_id` and are consecutive in it.
- [ ] Performance: 4,500 genes → pair frame in < 3 s on 1 CPU with all three layers present.

## Hardening owed (pass 2)

- A golden `pairs.tsv` for the 50-kb fixture (§11.2 "Golden … gene, pair and variant tables").
- The `partial`-gene decision re-examined against real Bakta output in the M10 pilot.
- Pair-level evidence for CV/DV, if M9's V11 panel turns out to need it (currently it does not).

## Questions for the lead

1. Extend the merged `giap.annotate.operons.adjacent_pairs` with `same_strand=False`, or write a
   second adjacency function in `giap.operons`? **Proposed default:** extend, so wrap handling and
   genome order have one definition (the M2 function is already imported by M6).
2. Do `partial` genes stay in the chain? **Proposed default:** yes for geometry (dropping them would
   make two genuine neighbours adjacent), no for tests (§7.3 eligibility).
3. `pair_id` scheme: §6.5 does not define one. **Proposed default:** `PAIR_<n>` with the
   threshold-independent ordering of step 9.
4. Evidence for CV/DV: NA (proposed) or computed with an arbitrary "downstream" choice?
   **Proposed default:** NA — no flag row in §7.4.5 uses them.
5. `terminator_between` called with `yes = 0.0` so the column is the raw maximum, with the cut-off
   applied in M6-07? **Proposed default:** yes; the column then survives a recalibration of
   `terminatornet_yes` (M6-12) without re-running the step.

## Lead decisions (2026-10-21, after the test-author)

1. **Mirror property.** The tests are right: OP, NO and SW keep `gene_a`, `gene_b` and Δ under
   mirroring; CV and DV swap `gene_a` and `gene_b` and negate Δ, because a is the left gene. The
   acceptance criterion's "negates no Δ" wording is corrected above.
2. **Circular replicon with 2 genes.** One adjacent pair (A, B); `pairs_raw` may hold both rows
   (M2-12 writes the wrap row too), and the extra row is ignored.
3. **`pair_frame` returns only the DataFrame.** The `layer_missing:<layer>` run flag (example 11,
   degenerate row "layer missing") is the `operons` step's job (M6-08).
4. **Test-author readings accepted:** `flags` is a list per row (`rna_between` the only value in
   pass 1); a gene absent from the gene table raises with a message naming the gene;
   `Evidence.rna` holds RNA rows only (the caller filters C15 by type); on a circular replicon
   with exactly 3 genes a sandwich pair (i, i + 2) is also the wrap-adjacent pair, so the
   property test does not check SW disjointness there (the row is still emitted by rule 4).
5. **Performance.** The 4,500-gene criterion is a benchmark in `benchmarks/` with a
   `benchmarks/baseline.json` entry (node, CPU, commit).
