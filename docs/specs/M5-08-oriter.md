# M5-08: Origin-to-terminus gradient: oriC resolution, replichore columns, descriptive regressions

Milestone M5 · track R (context statistics) · depends on: M5-02 (`choose_block_len`), M5-04 (C15
contract, oriC in the simulator); reuses merged `circular_block_bootstrap`, `percentile_ci` (M4-04) ·
estimate ≈ 260 code lines (+ ≈ 300 test lines) · **statistical** (descriptive with a bootstrap CI; label
`stats`) · **real data, no tool run:** unit tests read the recorded Bakta 1.12.2 JSONs already in
`tests/data/recorded/bakta/1.12.2-full/` (K-12 has **two** oriC calls, Nissle one, sim100k none) · no
Slurm

## Spec sections implemented

- §7.3.5(e) quoted: "Position. oriC comes from Bakta; the terminus is taken as oriC + L/2. Replichore
  position = circular distance from oriC / (L/2), in [0, 1]. Regression. TU-level log2 abundance on
  replichore position, per replichore and pooled, with a circular block-bootstrap CI. Descriptive only:
  `p_value = NA` and the flag `descriptive`, because evolution places transcription and translation genes
  near oriC. Without an oriC call the analysis is skipped with the flag `no_oric`."
- §9.3 `autocorrelation.tsv` rows with `statistic = oriter_slope`, `null_model = bootstrap`; genes.tsv
  `replichore`, `replichore_position`.
- docs/methods/bootstrap.md "The replichore regression's CI uses the same rule" (block length).

## Public interface

`src/gippa/context/oriter.py` (pure)
```python
def oric_position(features: pd.DataFrame | None, contig: str, length: int) -> float | None: ...
def replichore(mid: ArrayLike, oric: float, length: int) -> tuple[npt.NDArray[np.str_], F64]: ...
def gene_replichores(genes, replicons, features) -> pd.DataFrame: ...
    # locus_tag, replichore {right, left} or NA, replichore_position or NA
def oriter_rows(genes, tus, replicons, features, cfg, rng) -> pd.DataFrame: ...   # AUTOCORR_COLUMNS
```
- **Replicons analysed:** `kind = chromosome` and `topology = circular`. Plasmids: no rows, NA gene
  columns. Linear or `unknown`-kind replicons: no rows (Q3).
- **oriC.** C15 rows with `type = oric` on the contig; midpoint of each = (start + end)/2 mod L (end may
  be unwrapped). One call → its midpoint. Several calls → the circular mean of the midpoints (angles
  2π·m/L) and a log event `context_stats.oriter.multiple_oric` listing them (Q1). No call, or no
  `features.parquet` → no oriC.
- **Replichore.** Gene midpoint m (mod L); d = (m − o) mod L; d ≤ L/2 → `right`, position d/(L/2);
  otherwise `left`, position (L − d)/(L/2). Position ∈ [0, 1]; 0 = oriC, 1 = terminus.
- **TU values.** y_t = mean `abundance_log2` of the TU's detected genes with a value (eligible genes of
  `replicon_index`); TUs without one are dropped. x_t = replichore position of the TU midpoint (first gene
  start to last gene end, circular midpoint), replichore of that midpoint.
- **Fits.** OLS y = a + b·x, separately for `right` TUs, `left` TUs and all TUs (`pooled`). Rows:
  `statistic = oriter_slope`, `stat ∈ {oriter_slope_right, oriter_slope_left, oriter_slope_pooled}`,
  `lag` NA, `null_model = bootstrap`, `value = effect = b` (log2 per replichore, i.e. ter minus ori of
  the fitted line), `null_mean`/`null_lo`/`null_hi` NA, `p_value` NA, `q_value` NA, `method =
  ols_replichore`, `n_perm` = valid bootstrap replicates, `flags` = `descriptive`.
- **CI.** ℓ = `context.islands.block_len` if an integer, else `choose_block_len` on the residuals of the
  pooled fit in genome order (gene-label null, B = `stats.permutations`); B = `stats.bootstrap` index
  replicates `circular_block_bootstrap(arange(T_s), ℓ, B, rng)` of the TU rows of each fit's sequence (the
  replichore's TUs in order from oriC to terminus; pooled: all TUs in genome order); refit; 95 %
  `percentile_ci` of the valid slopes; fewer than B/2 valid → CI NA.
- **RNG.** `gippa.rng.for_step("context_stats.oriter", seed)`; replicons in C1 order; pooled, right,
  left.

## Degenerate cases (the estimator does not exist)

| Situation | Behaviour |
|---|---|
| circular chromosome without an oriC call | one row `stat oriter_slope_pooled`, `value` NA, flags `descriptive;no_oric`; gene columns NA |
| a fit with < 3 TUs with values, or zero variance of x | that row `value`/CI NA, flags `descriptive;not_estimable` |
| bootstrap replicate with zero variance of x | slope NaN, dropped from the CI (counted out of `n_perm`) |
| fewer than B/2 valid replicates | CI NA, flag `not_estimable` added |
| no detected gene on the chromosome | all three rows NA with `not_estimable` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | L = 1,000, oriC midpoint 900; gene midpoints 100, 400, 500, 650 | (right, 0.4), (right, 1.0), (left, 0.8), (left, 0.5) |
| 2 | gene start 990, end 1,010 (unwrapped), oriC 0 | midpoint 0 → (right, 0.0) |
| 3 | recorded K-12 Bakta JSON (two oriC calls, 0-based midpoints 3,884,032 and 3,925,819 on L = 4,641,652) | one oriC at their circular mean 3,904,925.5 (± 1 bp); `multiple_oric` logged |
| 4 | recorded Nissle JSON | oriC midpoint 5,010,344 (0-based) on CP171242.1; the three plasmids get no rows (oriV/oriT ignored) |
| 5 | recorded sim100k JSON (no oriC) | one pooled row with `no_oric`; gene columns NA |
| 6 | right TUs x = [0, 0.25, 0.5, 0.75, 1], y = [10, 9.4, 9.1, 8.4, 8.1] | `oriter_slope_right` value −1.92 |
| 7 | `oriter_gradient` scenario, 2 log2 cosine, ρ = 0.5, 4,500 genes | pooled slope CI covers the true OLS slope of the planted curve in ≥ 90 % of 50 seeds (recorded, `slow`) |
| 8 | any run | no oriter row has a `p_value`; every one carries `descriptive` |

## Acceptance criteria

- [ ] `tests/unit/context/test_oriter.py`: examples 1–6, 8 and the degenerate table (examples 3–5 parse
  the recorded JSONs through `parse_bakta_json` → features frame; no tool runs).
- [ ] Property: replichore positions ∈ [0, 1]; rotating every coordinate and the oriC by the same offset
  leaves (replichore, position) unchanged; swapping strand never matters.
- [ ] Example 7 as a `slow` coverage check (recorded, not gated: the analysis is descriptive).
- [ ] mypy, ruff; coverage ≥ 90 %.

## Methods note (statistician sign-off)

`docs/methods/oriter.md`: why descriptive (gene placement is evolved, not a null), replichore geometry,
the oriC choice with several calls, TU values, the block-bootstrap CI and its block length.

## Questions for the lead

1. Several oriC calls on one chromosome (K-12: two calls 42 kb apart; the true oriC is the second). Options:
   circular mean (proposed), first call, longest call, or skip with `no_oric`. **Proposed default:**
   circular mean + log event; the difference is < 1 % of a replichore and the analysis is descriptive.
2. `replichore` values `right`/`left` (right = increasing coordinate from oriC) for genes.tsv, and the
   row label in the free `stat` column (`oriter_slope_right/left/pooled`) instead of a new column.
   **Proposed default:** as stated; it needs no change to the §9.3 column list.
3. Linear chromosomes (rare: e.g. Streptomyces, Borrelia) and `unknown`-kind contigs: skipped without a
   row. **Proposed default:** skip; the spec's geometry (terminus = oriC + L/2) assumes a circle.
4. Gradient on detected genes' abundance only (censoring ignored). **Proposed default:** yes, this is a
   descriptive regression of TU-level log2 abundance as §7.3.5(e) says.

## Lead decisions (2026-10-20, after the statistician's attack check)

31. Q1: circular mean of the oriC calls only if all lie within 0.05·L of each other (circular distance);
    otherwise `no_oric` plus a `multiple_oric` log event. (K-12: two calls 42 kb apart → mean 21 kb off
    the true oriC, 0.009 of a replichore; slope unaffected.)
32. Q2: `right`/`left` labels in `stat`.
33. Q3: linear and unknown-kind contigs are skipped.
34. Q4: detected genes only (descriptive).
35. CI block length ℓ' (ADR-0021, accepted); until accepted, tests use an integer `block_len`.

## Lead decisions after the test-author's report (2026-10-21)
1. **Per-arm block length:** `ℓ_s = block_len_prime(ℓ̂, T_s)` = min(20, 2ℓ̂+1, T_s) with `T_s` the *fit's
   own* TU count, and an integer `context.islands.block_len` used as `min(block_len, T_s)` — the same
   clamp the M5-06 oracle uses. `circular_block_bootstrap` raises when `block_len > n`, so a 20-TU arm
   cannot carry an ℓ' derived from 2,000 pooled TUs. `ℓ̂` is estimated **once per replicon** from the
   pooled fit's residuals. Re-estimating per arm is more faithful but costs more draws for no measured
   gain; revisit only if the statistician asks.
2. **Draw order:** a non-estimable fit consumes **no** bootstrap draws. Replicons in C1 order; per
   replicon the ℓ̂ permutations, then bootstrap replicates for pooled, right, left, skipping absent fits.
3. **Row order** within a replicon is pooled, right, left (M5-09 re-sorts for the table; the parity tests
   need a fixed order).
4. **`multiple_oric` is a logging record** from the pure module, using `gippa.io.events`' convention
   (`extra={"event": ..., "data": ...}`) rather than a diagnostics object threaded back through the step.
   Logging is not file I/O; this keeps the signature clean.
5. The event fires on **both** branches — when the calls are averaged and when they are > 0.05·L apart
   and refused (decision 31's `no_oric`).
6. The fast coverage arm stays in **NIGHTLY**: 47 s on a session node is under the 60 s rule but will
   exceed it on a 4-core CI runner.
7. **C15 routing confirmed** (M5-04 decision 19 stands): `features` does **not** join `load_m4_inputs` /
   `FILE_NAMES` / `M4Inputs`. `oriter_rows` takes `features: pd.DataFrame | None`, and `context_stats`
   gains it as an optional step input — the plumbing is M5-09's, as its spec already says. A missing
   features table gives `no_oric` rows and a successful run, not a contract issue.
8. **Example 7 stays recorded, not gated**, pending the statistician's review of the implementation PR:
   the test-author independently reproduced ADR-0021's 0.96 coverage at ρ 0.5 and measured the estimand
   at −2.42 against the closed form −2.432, so promotion is plausible — the statistician decides.
