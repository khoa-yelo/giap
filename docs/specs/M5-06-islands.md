# M5-06: Islands: TU scores, detrending, block length, circular block bootstrap, FWER, merged islands

Milestone M5 · track R (context statistics) · depends on: M5-01 (Blom, LOESS), M5-02
(`choose_block_len`), M5-03 (scan, FWER, merge); reuses merged `circular_block_bootstrap`,
`tu_block_bootstrap` (M4-04) · estimate ≈ 380 code lines (+ ≈ 400 test lines) · **statistical** (label
`stats`) · code can start now on hand-built frames; the gate calibration is M5-07

## Spec sections implemented

- §7.3.5(d) steps 1–7 quoted: "1. Gene scores. Within each replicon, rank eligible genes by log2
  abundance, with undetected genes tied at the lowest rank (left-censored). Convert to Blom normal scores
  … 2. TU score z_t = the mean of its genes' scores. One operon therefore cannot form an island alone;
  the minimum window is 2 TUs. 3. Detrend with LOESS along TU order (span 10 % of TUs, circular padding on
  circular replicons). 4. Scan all windows of 2–20 consecutive TUs (wrapping on circular replicons) … 5.
  Null. A circular block bootstrap of the detrended TU scores, with block length ℓ = the smallest lag at
  which TU-level Moran's I is not significant at 0.05 (gene-label null), bounded to [1, 20]. For each
  replicate b, M_b = max over windows of |S_w|. 6. FWER … 7. Report per island: ID, direction, TU and
  gene span, coordinates; mean log2 abundance, detected fraction; S, `p_value` (FWER), with `q_value = NA`
  and `method = scan_max_stat`; effect = mean detrended score, with a TU-bootstrap CI."
- §9.3 `tables/islands.tsv` columns; §9.3 genes.tsv `island_id`; §6.5 IDs `ISL_<n>`.
- §5.2 `context.islands: {min_tus: 2, max_tus: 20, detrend_span_frac: 0.10, block_len: auto}`,
  `stats.bootstrap`, `stats.permutations`, `stats.alpha.context`.

## Public interface

`src/giap/context/islands.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class TuSeries:
    contig: str; circular: bool
    tu_ids: tuple[str, ...]          # scored TUs in genome order
    genes: tuple[tuple[str, ...], ...]   # scored genes of each TU (transcription order)
    score: F64                       # z_t (mean Blom score)
    detrended: F64                   # z_t − LOESS trend

def gene_scores(ix: pd.DataFrame) -> pd.Series: ...              # Blom per scored gene (index locus_tag)
def tu_series(ix: pd.DataFrame, tus: pd.DataFrame, contig: str, circular: bool, span_frac: float) -> TuSeries: ...
def block_length(ts: TuSeries, cfg: GiapConfig, rng) -> int: ...
def scan_replicons(series: Sequence[TuSeries], cfg: GiapConfig, rng) -> IslandScan: ...
def islands_frame(scan: IslandScan, genes: pd.DataFrame, cfg: GiapConfig, rng) -> pd.DataFrame: ...
def gene_island_map(islands: pd.DataFrame) -> pd.DataFrame: ...   # locus_tag → island_id
```
1. **Scored genes.** From `replicon_index` (eligible: not partial, `p_detect` present). A gene is
   *censored* when not detected; a detected gene with finite `a` is uncensored; a detected gene with NaN
   abundance (possible only under `abundance.missing = exclude`) is **not scored** (dropped). n = scored
   genes of the replicon; `blom_scores`.
2. **TU score.** z_t = mean of the scored genes of TU t (`tu_id` from the gene table); TUs with no scored
   gene are dropped; TU order = genome order of each TU's first scored gene (C5 `start`, then `tu_id`).
3. **Detrend.** `loess_detrend(z, pos = 0 … T − 1, frac = detrend_span_frac, circular, period = T)`.
4. **Block length.** `context.islands.block_len` if an integer, else `choose_block_len(detrended,
   circular, B = stats.permutations, rng, alpha = 0.05, max_len = 20)`; then ℓ ≤ T.
5. **Replicates.** For b = 1 … B (`stats.bootstrap`): each scanned replicon is resampled by
   `circular_block_bootstrap(detrended, ℓ, B, rng)` (also on linear replicons: Politis–Romano wrap
   accepted, docs/methods/bootstrap.md), and M_b = **max over all scanned replicons** of
   `scan_null_max` (genome-wide family; Q1).
6. **Observed scan.** `scan_max(detrended, min_tus, max_tus, circular)` per replicon; p_w = `fwer_p(|S_w|,
   M)`; significant = p_w ≤ `stats.alpha.context` (0.05); `merge_windows` per replicon and direction.
7. **Islands.** One row per cluster, from its **best window**: `direction`, `contig`, `tu_first`,
   `tu_last` (TU IDs; on a wrapping window tu_first is the TU at the window start), `n_tus`, `genes` (all
   eligible genes of those TUs, genome order), `n_genes`, `start` = start of the first TU, `end` = end of
   the last TU (unwrapped, `end > length` when the window crosses the origin; §6.5), `mean_abundance_log2`
   = mean `a` of the detected genes with a value (NA if none), `detected_frac` = detected / eligible
   genes, `scan_s` = signed S of the best window, `block_len` = ℓ of the replicon, `stat = scan_s`,
   `effect` = mean detrended score of the best window's TUs, `ci_low`/`ci_high` = 95 %
   `tu_block_bootstrap` of that mean over the window's TUs (B = `stats.bootstrap`; flag `ci_few_units`
   when < 10 TUs, the existing M4-16 flag), `p_value` = p_w of the best window, `q_value` NA,
   `method = scan_max_stat`, `null_model = circular_block_bootstrap`, `n_perm` = B (valid M_b).
8. **IDs.** `ISL_<n>`, n = 1, 2, … by contig order, then `start`, then direction (high first).
9. **Gene map.** `island_id` per gene; a gene in two islands (only possible for a high and a low island
   overlapping) gets the one with the smaller p, then larger |S|, then smaller ISL number.
10. **RNG.** `giap.rng.for_step("context_stats.islands", seed)`: block-length permutations of every
    replicon (C1 order), then the bootstrap draws (replicons in C1 order), then the effect CIs (island
    order).

## Degenerate cases (the estimator does not exist)

| Situation | Behaviour |
|---|---|
| replicon with T < 2·`min_tus` (= 4) scored TUs | not scanned, contributes nothing to M_b; one `context_stats.islands.not_scanned` log event (no table row) |
| no replicon scanned | `islands.tsv` header only; every `island_id` NA |
| all genes undetected (all Blom scores 0) | z = 0 everywhere → S = 0, p_w = 1 → no island |
| B valid M_b = 0 (cannot happen once one replicon is scanned) | p NaN → no island |
| island with no detected gene (a low island) | `mean_abundance_log2` NA, `detected_frac` 0; still reported |
| best window of 2 TUs | CI from 2 units (bootstrap valid but coarse) → flag `ci_few_units` |
| `block_len` integer > T | bounded to T |
| `max_tus` ≥ T on a circular replicon | windows capped at T − 1 (M5-03) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | genes G1–G4 (one replicon), detected log2 [5, –, 3, –] (G2, G4 undetected) | gene scores [1.0491, −0.6289, 0.2993, −0.6289] |
| 2 | TU_1 = {G1, G2}, TU_2 = {G3}, TU_3 = {G4} with example 1's scores | z = [0.2101, 0.2993, −0.6289] |
| 3 | detected gene with NaN abundance under `missing = exclude` | not scored; its TU keeps its other genes; if none left the TU is dropped |
| 4 | 1,000-gene `ar1_background(0.5)` with a planted island of 20 TUs at +2 SD | one `high` island overlapping the planted TUs, `p_value ≤ 0.05`, `method scan_max_stat`, `q_value NA` |
| 5 | same with −2 SD | one `low` island; `scan_s < 0`; `detected_frac` below the replicon's |
| 6 | two significant overlapping high windows TU_10–TU_14 (S 5.1) and TU_13–TU_20 (S 5.6) | one island, best window TU_13–TU_20, `n_tus 8`, `scan_s 5.6` |
| 7 | circular replicon, best window wraps (last two TUs + first two) | `tu_first` = the TU at the window start (near the end), `end > length` |
| 8 | `block_len: 5` in config | `block_len` 5 in every row; no block-length permutations drawn |
| 9 | 30-TU plasmid + chromosome | plasmid scanned (T ≥ 4) and in the same M_b family |
| 10 | 3 scored TUs on a plasmid | plasmid not scanned; chromosome results unchanged by its presence (bitwise, same seed: its draws are skipped) |
| 11 | same seed twice | identical `islands` frame and gene map |

## Acceptance criteria

- [ ] `tests/unit/context/test_islands.py`: examples 1–3, 6–11 and every degenerate-table row (small
  hand-built frames; examples 4–5 as `slow` tests on the simulator).
- [ ] Property: TU scores of a replicon are invariant to the order of rows in the gene table; relabelling
  TU IDs does not change any number; islands never contain a TU outside the replicon; on linear replicons
  no island wraps.
- [ ] Hand-check of step 5: with `stats.bootstrap = 9` and a fixed generator, `p_value` equals
  (1 + #{M_b ≥ |S|})/10 recomputed in the test from `scan_null_max`.
- [ ] Performance: 4,500-gene chromosome (≈ 2,500 TUs), B = 1,000 bootstrap + 1,000 block-length
  permutations < 20 s on 4 CPUs (benchmark recorded).

## Methods note (statistician sign-off)

`docs/methods/islands.md`: left-censored scores and the effect of the missing-value policy, TU means,
detrending (and that it attenuates islands longer than about half the span), block length, the bootstrap
null, genome-wide maximum, merging, effect CI, degenerate cases; evidence from M5-07.

## Questions for the lead

1. FWER family: per replicon (M_b within each replicon, as the step text reads literally) or genome-wide
   (M_b = max over replicons)? §8.1 names one "Islands" claim with FWER 0.05. **Proposed default:**
   genome-wide maximum — FWER ≤ 0.05 for the whole run; per-replicon would give a run-level FWER up to
   0.05 × number of replicons.
2. Detected genes without abundance under `missing = exclude`: dropped from the scores (proposed), or
   ranked as censored? **Proposed default:** dropped (they are detected, so censoring them would be wrong,
   and their rank is unknown).
3. Island coordinates and gene list from the best window (spec: "represented by its best window") or the
   union of the merged windows? **Proposed default:** best window, as the spec says; the union is not
   reported (no new column, R1).
4. Significance threshold: `stats.alpha.context` (0.05) — §7.3.5 hard-codes 0.05. **Proposed default:**
   read it from `stats.alpha.context` ("the single source of every significance level", §5.2); the block
   length rule stays at the literal 0.05.

## Lead decisions (2026-10-20, after the statistician's attack check)

ADR-0021 accepted by Khoa 2026-10-20: the spec's block length gives FWER 0.205–0.27 at ρ ≥ 0.5. When accepted:
21. Q1: genome-wide maximum over replicons; each replicon keeps its own ℓ'.
22. Q2: detected genes without abundance are not scored.
23. Q3: report each island by its best window.
24. Q4: threshold `stats.alpha.context`; the block-length rule keeps the literal 0.05.
25. Block length ℓ' = min(20, 2ℓ̂ + 1, T); integer `context.islands.block_len` overrides.
26. Effect CI: the spec's TU-bootstrap CI; methods note says it is conditional on the selected window and
    ignores within-window autocorrelation; M5-07 records its coverage.
27. Keep the combined detection/abundance estimand; gate the short-gene prophage attack.
