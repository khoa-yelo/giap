# M5-05: Autocorrelation: Moran's I and join counts per replicon with both nulls; calibration

Milestone M5 · track R (context statistics) · depends on: M5-01 (LOESS), M5-02 (Moran, join count,
permutations) · calibration needs M5-04 scenarios (M4-05 harness, M4-07/08 simulator merged) ·
estimate ≈ 330 code lines (+ ≈ 400 test lines) · **statistical** (label `stats`) · code can start now on
hand-built frames; calibration full mode on **Slurm** (R = 1,000)

## Spec sections implemented

- §7.3.5(c) quoted: "Values. x = log2-abundance residuals of detected genes, after removing the replicon
  mean and a wide LOESS trend along gene index (span = 10 % of genes, circular padding on circular
  replicons). Moran's I at lag k = 1 … 20 … Join-count for detection. BB_k = the number of lag-k pairs of
  eligible genes where both are detected. Two nulls, B permutations each: gene label … TU block …
  Reporting. Two-sided p per lag and null (centred on the null mean); BH across lags within each null;
  2.5–97.5 % envelopes for V10."
- §8.1 "Autocorrelation | lag | gene-label and TU-block permutation | BH across lags | 0.05".
- §8.4 "Autocorrelation | gene-level noise (label null); TU-structured noise without inter-TU structure
  (TU null) | level | 100 | 1,000"; §8.5 attack 1; §7.3.6 "autocorrelation tests at level under the
  TU-block null".
- §9.3 `tables/autocorrelation.tsv` rows with `statistic ∈ {morans_i, join_count_bb}`; §9.2 test columns.
- §5.2 `context.moran_lags: [1, 20]`, `stats.permutations`.

## Public interface

`src/giap/context/autocorrelation.py` (pure)
```python
AUTOCORR_COLUMNS: tuple[str, ...]   # contig, statistic, lag, null_model, value, null_mean, null_lo,
                                    # null_hi, stat, effect, ci_low, ci_high, p_value, q_value, method,
                                    # n_perm, flags
def moran_values(ix: pd.DataFrame, circular: bool, span_frac: float = 0.10) -> F64: ...
def replicon_autocorrelation(ix: pd.DataFrame, contig: str, circular: bool, lags: Sequence[int],
                             B: int, rng: np.random.Generator) -> pd.DataFrame: ...
def autocorrelation_rows(genes: pd.DataFrame, replicons: pd.DataFrame, cfg: GiapConfig,
                         rng: np.random.Generator) -> pd.DataFrame: ...   # adds q_value
```
- **Index.** `giap.context.regions.replicon_index(genes, contig)` (merged M4-17): eligible genes = not
  partial and `p_detect` present, genome order, idx 0 … N − 1, columns `d`, `a` (log2 abundance of
  detected genes, NaN otherwise), `tu_id`. Same index as the region tests.
- **Values (`moran_values`).** Genes with finite `a`: subtract their mean, then `loess_detrend` with
  `pos` = idx, frac 0.10, `circular`, `period` = N; x = residuals at those genes, NaN elsewhere.
- **Lags.** `context.moran_lags = [lo, hi]` is an inclusive range (default 1 … 20); lags > `max_lag(N,
  circular)` are reported as rows with NA statistics and flag `not_estimable` (so every replicon has the
  same row grid).
- **Nulls** (B = `stats.permutations`, one permutation matrix per null and replicon, shared by all lags
  and both statistics):
  - `gene_label`: Moran — permute x among the genes with values (`movable` = has value, missingness
    fixed); join count — permute d among all eligible genes;
  - `tu_block`: `tu_block_permutations(tu_blocks(tu_id, circular))` applied to x (with its NaNs) and d
    jointly.
  The detrended x is permuted as is (no re-detrending per permutation).
- **Rows.** Per replicon (C1 order) × statistic (`morans_i`, `join_count_bb`) × `null_model`
  (`gene_label`, `tu_block`) × lag: `value` = I_k or BB_k; `lag_test` gives `null_mean`, `null_lo`,
  `null_hi`, `p_value`, `n_perm` = valid null draws; `stat` = statistic name; `effect` = value −
  null_mean; `ci_low`/`ci_high` = NA (see Q2); `method` = `moran_i_permutation` /
  `join_count_permutation`.
- **q_value.** BH (`giap.stats.multitest.bh`, NaN-aware) within each (contig, statistic, null_model) over
  its lags.
- **Flags.** `not_estimable` whenever the value or p is NA (invalid lag, n < 2 values, no pair, zero
  variance, degenerate null); additionally `few_values` when fewer than 10 genes have values (the test
  still runs). (Vocabulary additions: M5-09.)
- **RNG.** `giap.rng.for_step("context_stats.autocorrelation", seed)`; replicons in C1 order, within a
  replicon: label matrix for x, label matrix for d, then the TU-block matrix.

## Degenerate cases (the estimator does not exist)

| Situation | Rows |
|---|---|
| replicon with no eligible gene | all rows NA, `n_perm 0`, flag `not_estimable` |
| no detected gene | Moran rows NA (no values); join-count `value` 0, null constant → `p_value` NA, `not_estimable` |
| every eligible gene detected | Moran rows computed normally; join-count `value` = number of lag-k pairs, null constant → `p_value` NA, `not_estimable` |
| 1 gene with a value | Moran rows NA, `not_estimable` |
| all values equal after detrending | Moran NA (zero variance) |
| single TU on the replicon | `tu_block` null is the identity → `p_value` NA, `not_estimable`; `gene_label` rows normal |
| circular replicon with N = 30 | lags 1–14 computed; 15–20 NA, `not_estimable` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 10-gene circular replicon, every gene detected and its own TU, detrended x alternating ±1 | `morans_i` lag 1 value −1.0; lags 5–20 NA (`max_lag` 4) |
| 2 | 1,000-gene default layout, one replicon | 2 × 2 × 20 = 80 rows, sorted contig, statistic, null_model, lag |
| 3 | `tu_null` scenario (strong within-TU correlation, none across) | `morans_i` lag 1 `gene_label` p ≤ 0.01 and `tu_block` p > 0.01 in ≥ 80 % of 20 seeds |
| 4 | any replicon | `q_value` ≥ `p_value` row-wise; q NA exactly where p NA |
| 5 | `stats.permutations = 199` | `n_perm` 199 on every estimable row |
| 6 | 30-gene circular plasmid + 1,000-gene chromosome | plasmid lags 15–20 NA with `not_estimable`; BH family of the plasmid has 14 tests |
| 7 | same seed twice | identical frames (bitwise) |
| 8 | a replicon with no detected gene | Moran rows NA; join-count rows `value 0`, `p_value` NA, `not_estimable` |

## Acceptance criteria

- [ ] `tests/unit/context/test_autocorrelation.py`: examples 1, 2, 4–8 and every row of the degenerate
  table.
- [ ] **Calibration** `tests/calibration/test_autocorrelation.py` (suite `autocorrelation`, R = 100 fast /
  1,000 full, 1,000-gene layout fast, 4,500 full; `p_detect = p_detect_true`, abundance =
  `abundance_log2_true` for detected genes):
  - `label_null`: `check_level_pooled` over the 20 lags at α ∈ {0.01, 0.05, 0.10} for each
    (statistic, null_model) — both nulls must hold;
  - `tu_null`: the same for `null_model = tu_block` (both statistics);
  - teeth (recorded, not gated): `gene_label` rejection rate under `tu_null` at lag 1 ≥ 0.8.
- [ ] **Attack 1** (statistician, `tests/calibration/attacks/test_attack_autocorr.py`): `tu_null` plus an
  ori–ter cosine of 2 log2 (`oriter_gradient` with ρ = 0) → `tu_block` rows at level (the LOESS removes
  the gradient); the same without detrending is recorded to show the attack has teeth.
- [ ] Performance: 4,500-gene replicon, B = 1,000, both nulls, 20 lags < 30 s on 4 CPUs (benchmark).
- [ ] mypy, ruff; coverage ≥ 90 % of the module.

## Methods note (statistician sign-off)

`docs/methods/autocorrelation.md`: values and detrending, the two estimands, the missingness choice of the
label null, TU blocks, degenerate nulls, BH across lags within (contig, statistic, null), why no CI
(Q2), calibration and attack evidence, the known small negative bias of Moran's I after detrending
(≈ −1/k) and that it enters observed and null alike only for the TU-block null (it is checked by the
calibration).

## Questions for the lead

1. Gene-label null for Moran: permute values among genes **with** values (missingness fixed) or among all
   eligible genes (missingness moves)? **Proposed default:** among genes with values — the join count
   already tests where detected genes sit; Moran then asks only about abundance.
2. §8.2 asks for a 95 % CI with every p. Moran and join-count rows have the null envelope
   (`null_lo`/`null_hi`) but no CI of the statistic itself; a block-bootstrap CI for a lag-k
   autocorrelation needs blocks ≫ k and is uncalibrated. **Proposed default:** `effect` = value −
   null_mean, `ci_low`/`ci_high` = NA, documented in the methods note; record as an ADR (deviation from
   §8.2 for this table only) if the lead agrees.
3. BH family per (contig, statistic, null_model) — §8.1 says "across lags" without naming the contig.
   **Proposed default:** per contig; replicons answer separate questions.
4. Join count ignores detectability: a cluster of short genes (low `p_detect`) gives low BB under both
   nulls even with no abundance structure. The spec's nulls are kept; the statistician adds a recorded
   (not gated) scenario with detectability clustered beyond TUs. **Proposed default:** keep the spec;
   if that scenario shows the TU-block join count far off level, propose an ADR with the residual join
   count Σ (d_i − p_i)(d_j − p_j) (same idea as ADR-0018).

## Lead decisions (2026-10-20, after the statistician's attack check)

ADR-0020 accepted by Khoa 2026-10-20: the spec's "detrend then permute" is anti-conservative (pooled 0.096 /
0.155 at α 0.05) and the raw join count fails under clustered detectability. When accepted:
14. Q1: the gene-label null keeps the missingness pattern; it permutes the pre-detrended values among genes
    with values and detrends every replicate.
15. Q2: CI NA with flag `ci_unestimated` (ADR-0020 decision 4).
16. Q3: BH family per (contig, statistic, null_model); contigs with < 10 values get `few_values`.
17. Q4: residual join count now (ADR-0020 decision 3): `value` = raw BB_k, `stat`/`effect` = residual
    statistic − null mean.
18. Gates per (statistic, null), α ∈ {0.01, 0.05, 0.10}: lag-1 level, lag-20 level, BH-family
    any-rejection rate (FWER bound); pooled 20-lag rate recorded only.
19. Statistician attack file tests/calibration/attacks/test_attack_autocorr.py: (a) TU null + 2 log2
    cosine, (b) label null + gradient, (c) short-gene prophage under TU null (join count) — gated; (d)
    no-detrend teeth — recorded.
20. Methods note: large-lag negative I under strong short-range correlation is a span artefact.
21. (M5-01 decision 4) A replicon whose 10 % span gives k < 4 points has all-zero residuals: its rows are
    NA with flag `few_values` and do not enter a BH family.
22. (M5-01 statistician review) Supersedes the "Values" bullet: detrend on ordinals, `pos` = 0…n−1 over the
    n genes with values, `period` = n (circular), re-detrended per permutation as ADR-0020 states; one
    `loess_matrix` per (n, span, circular) serves the observed series and all B replicates (compute the
    observed fit as H @ y too). Consumers test `k < 4` (never `resid == 0`). Rows flagged `few_values`
    under decision 21 also carry `not_estimable`; the "< 10 values, test still runs" text in "Flags" and
    decision 16 is superseded for k < 4 replicons (≈ n ≤ 30).
23. (test-author report) `few_values` applies to the **morans_i** rows only: the join count never touches
    the detrend and stays estimable on a small plasmid. The trigger is the number of genes *with values*.
24. Row order is lexicographic on (contig in C1 order, statistic, null_model, lag), so `join_count_bb`
    precedes `morans_i` and `gene_label` precedes `tu_block`. The §9.3 enum order does not apply to rows.
25. The 0.30 wide pre-trend span and the 0.10 detrend span are **module constants**, not config keys:
    ADR-0020 fixes them and the calibration is conditional on them. No `context.autocorr` section.
26. `moran_values` on a degenerate replicon returns all-zero residuals, never NaN (consumers test k < 4).
27. Degenerate join-count rows ("no detected gene", "every gene detected") are written for a **constant**
    `p_detect`; with the Pearson statistic and varying p the null is not constant and the row is estimable.
28. **Join-count teeth (statistician question).** With giapsim's truth `p_detect` the Pearson residuals are
    exact Bernoulli residuals, so the join count has almost no power under any abundance-structured
    scenario (recorded 0.04 at lag 1 under the TU null, against 1.00 for Moran). Spec example 3's teeth
    claim therefore holds for `morans_i` only. Record both rates and gate neither for now; the
    statistician decides at review whether a miscalibrated-p scenario (constant `p_source`, or the
    short-gene prophage) must carry a join-count teeth gate.
29. The calibration file is `tests/calibration/test_autocorrelation_calibration.py` (basenames are unique
    across tests/); the spec's `test_autocorrelation.py` would collide with the unit file.
30. ADR-0020 decision-19 attack arms (a) and (b) need M5-04's `oriter_gradient` and wait for #119; arms
    (c) short-gene prophage and (d) the no-detrend teeth are written in this PR under
    tests/calibration/attacks/. The benchmark (4,500 genes, B = 1,000, < 30 s on 4 CPUs) is in this PR.
31. Arm (c) keeps **both** lag-1 and lag-20 gates in full mode, at 38 min, because ADR-0020 decision 5
    gates lag 20 explicitly; the cost is paid only by the nightly full array. This makes the existing
    backlog item (shard the full calibration tier across array tasks) a prerequisite for the next full
    arm, not a nice-to-have.
32. `benchmarks/baseline.json` gains `autocorrelation_rows_b1000_n4500_20_lags` only when the **real**
    implementation is measured (the implementer records it with node, CPU and job metadata). The
    stand-in's 2.57 s against a 30 s budget is a sanity check, not an accepted number.
33. #119 (M5-04) is merged, so `tests/helpers/autocorr_calib.py` and the attack file drop their local
    scenario definitions for `giapsim.scenarios`, and ADR-0020 attack arms (a) (TU null + a 2 log2
    ori–ter cosine) and (b) (label null + gradient) are written now.
34. (after the scenario switch) `tu_null_1000.morans_i.tu_block.lag1` at α 0.10 moved 0.100 → 0.190
    against a bound of 0.160 in the **fast** tier: block 1 fails, block 2 passes, so it records
    `unconfirmed`. At R = 100 one SE is 0.03, so it is about 1 SE over, and it is the only check in the
    suite needing confirmation — but it is the headline TU-null claim. The implementer runs the **full**
    tier (R = 1,000) for the main suite once the module exists and reports that number; if the full tier
    is also over, it goes to the statistician before anything is gated differently. No bound is touched
    meanwhile.
35. Arm (c) full lag 1 at α 0.01 now sits exactly on its bound (0.016 of 0.016, pass), up from 0.012;
    recorded so a later drift is visible.
