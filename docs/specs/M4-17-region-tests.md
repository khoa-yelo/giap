# M4-17: Region tests: O − E and Hodges–Lehmann shift under circular gene-index shifts

Milestone M4 · track R · depends on: M4-01 (BH), M4-03 (shifts, perm_p), M4-04 (bootstrap, HL)
(calibration needs M4-05, M4-07, M4-08) · estimate ≈ 380 code lines (+ ≈ 400 test lines) ·
**can start now** · gate item: "circular-shift tests calibrated"

## Spec sections implemented

- §7.3.5(b) quoted: "Region tests (`regions_stats.tsv`), for the types `prophage`, `bgc`, `crispr` and
  `amr`: Detection statistic. O − E = Σ_{g∈R} (d_g − p_g) over the region's eligible genes. Abundance
  statistic. The Hodges–Lehmann shift: the median of pairwise differences in log2 abundance between
  the region's detected genes and the replicon's other detected genes. Null: circular gene-index shifts
  within the replicon. The region's gene-index set I becomes (I + k) mod N, for k = 1 … N − 1
  (exhaustive while N ≤ `exhaustive_shift_max_genes`, otherwise B random shifts). Recompute both
  statistics at every shift. On linear contigs, only shifts that keep the set inside [0, N) are used.
  Two-sided p = (1 + #{|S_k − S̄| ≥ |S_obs − S̄|}) / (1 + K) … CI. TU-block bootstrap: resample the
  region's TUs and, for the HL shift, the outside TUs too. Multiple testing. BH within each region
  type. Every row carries `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`,
  `null_model = circular_shift` and `n_perm = K`. Small regions. Regions with < 3 eligible genes are
  flagged `small_region`."
- §8.1 "Region effect | region | circular gene-index shifts within the replicon … | BH within type | 0.05".
- §8.4 "Region tests | AR(1) abundance, detection from p_detect only | level | 100 | 1,000".
- §8.5 attacks 2 (detectability confounded with region membership) and 4 (high-copy plasmids).
- §7.3.6: "exhaustive shift enumeration on circular and linear contigs"; property "O − E is unchanged by
  relabelling genes outside the region"; "circular shifts preserve region size"; calibration
  "region-test p-values uniform under autocorrelated nulls"; power "regions with detection shifts of
  10–50 %" (recorded, not gated).
- Appendix B.4.

## Public interface

`src/giap/context/regions.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class RegionTest:
    region_id: str; region_type: str; contig: str; stat: Literal["o_minus_e", "hl_shift"]
    observed: float; expected: float | None; effect: float; ci_low: float; ci_high: float
    p_value: float; n_perm: int; n_genes: int; n_eligible: int; flags: tuple[str, ...]

def replicon_index(genes: pd.DataFrame, contig: str) -> pd.DataFrame: ...
    # eligible genes of the contig in genome order with columns idx, locus_tag, tu_id, d, p, a (log2 or NaN)
def o_minus_e_null(r: npt.NDArray, idx: npt.NDArray, shifts: npt.NDArray) -> npt.NDArray: ...   # r = d − p
def hl_null(a: npt.NDArray, det: npt.NDArray, idx: npt.NDArray, shifts: npt.NDArray) -> npt.NDArray: ...
def test_region(genes, region_genes: Sequence[str], region_type: str, contig: str, cfg: GiapConfig,
                rng: np.random.Generator) -> list[RegionTest]: ...
def test_all_regions(genes, regions, region_genes, replicons, cfg, rng) -> pd.DataFrame: ...
    # adds q_value (BH per (region_type, stat)), method, null_model, test columns
```
- **Index:** the replicon's eligible genes in genome order, indices 0 … N − 1 (partial genes are not
  indexed, so shifts preserve the number of eligible region genes). Circular iff the replicon is
  circular. A region's index set I = indices of its eligible members on that replicon.
- **O − E:** S(I) = Σ_{i∈I} (d_i − p_i); `expected` = Σ p_i; `effect` = S_obs. Null over shifts by
  index arithmetic (vectorised: |I| × K gather).
- **HL:** `effect` = `hodges_lehmann(a[I ∩ detected ∩ observed], a[(¬I) ∩ detected ∩ observed])`,
  same at every shift; NaN when either side is empty. Shifts with NaN statistics are dropped
  (K counts valid shifts).
- **Shifts:** `choose_shifts(circular_shifts(N, I, circular), exhaustive_shift_max_genes, N,
  stats.permutations, rng)`. `n_perm = K`.
- **p:** `perm_p(S_obs, null, "two-sided")`; K = 0 → p NA.
- **CI:** `tu_block_bootstrap` with B = `stats.bootstrap`: O − E resamples the region's TUs (only region
  members in each TU count); HL resamples inside and outside TUs as two strata.
- **Flags:** `small_region` when n_eligible < 3 (tests still run if K > 0).
- RNG: `giap.rng.for_step("context_stats", seed)`, regions processed in sorted `region_id` order.
- `method` = `circular_shift_o_minus_e` / `circular_shift_hodges_lehmann`; `null_model = circular_shift`.
- Plasmid regions are **not** handled here (M4-18).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | region of 4 eligible genes, d = [1, 0, 1, 1], p = [0.9, 0.2, 0.5, 0.4] | `observed 3`, `expected 2.0`, `effect 1.0` |
| 2 | region detected log2 [10, 12], other detected [8, 9, 11] | HL `effect 1.5` |
| 3 | B.4: circular replicon of 100 eligible genes, observed more extreme than all 99 shifts | `n_perm 99`, `p_value 0.01` |
| 4 | linear contig, N = 10, region I = {2, 3, 4} | `n_perm 7` (shifts 1–5, 8, 9) |
| 5 | linear contig, N = 5, I = {0, 4} | K = 0 → `p_value NA`, `n_perm 0` |
| 6 | N = 25,000 > 20,000, `stats.permutations = 1000` | `n_perm 1000` (random, without replacement) |
| 7 | region with 2 eligible genes | flag `small_region`; test still reported |
| 8 | 3 prophages with p (O − E) = [0.01, 0.04, 0.03] | q = [0.03, 0.04, 0.04] (BH within prophage × o_minus_e) |
| 9 | region with no detected gene | HL row `effect NA`, `p_value NA`; O − E row normal |
| 10 | planted region d − p = +0.5 per gene for 30 genes on a 1,000-gene circular chromosome, iid elsewhere | `p_value ≤ 0.01` |

## Edge cases

- Region members on two contigs (should not happen for blocks) → each contig tested separately with
  `region_id` kept and flag recorded; AMR gene sets spanning replicons are split the same way.
- Region whose members are all partial → n_eligible 0 → rows with NA statistics and `small_region`.
- Identical statistics at many shifts (e.g. all p_detect equal): ties counted with tolerance (M4-03).
- Abundance policy `exclude`: detected genes without abundance drop out of HL on both sides.

## Acceptance criteria

- [ ] `tests/unit/context/test_regions.py`: examples 1–10.
- [ ] Property (`tests/property/test_region_props.py`): O − E unchanged by relabelling/permuting genes
  outside the region (values moved among outside indices); every shifted set has |I| members; on linear
  contigs no shifted set wraps; vectorised nulls equal a naive per-shift loop on random small inputs.
- [ ] **Calibration** `tests/calibration/test_region_tests.py` (suite `regions`, R = 100 fast / 1,000
  full): M4-07 layout (1,000 genes fast, 4,500 full; default regions + a 60-gene plasmid with its own
  AMR set) + M4-08 with `ar1_rho = 0.5` (full: ρ ∈ {0, 0.5, 0.9}), no region effects,
  `detection = "covariates_only"`; tests use `p = p_detect_true` (fast) and M4-13 cross-fitted
  `p_detect` (full). Criterion: `check_level` at α ∈ {0.01, 0.05, 0.10} for each (region_type, stat)
  separately (discrete p → level only).
- [ ] **Attack nulls** (statistician, `tests/calibration/attacks/`): (2) `short_genes_in="prophage"`
  with detection from covariates only and M4-13 `p_detect` → prophage O − E at level; (4)
  `plasmid_copy_log2 = {"p1": log2(20)}` → AMR-on-plasmid and chromosome region tests at level.
- [ ] Power (recorded, not gated): BGC detection logit shifts giving +10, 20, 30, 50 % detection,
  R = 100; rejection rates at α = 0.05 logged to `local/logs/calibration/regions_power.jsonl`.
- [ ] Performance: 4,500-gene chromosome, 30 regions, exhaustive shifts, B = 1,000 bootstrap: < 60 s
  on 4 CPUs (benchmark recorded).

## Methods note (statistician sign-off)

`docs/methods/region-tests.md`: index definition (eligible genes), both statistics, HL exact
computation, the circular-shift null and what it preserves (region size, spacing, autocorrelation),
linear restriction, two-sided centred p, TU-block CIs, BH family, small regions, calibration and
attack-null evidence, power table.

## Questions for the lead

1. Shift index over **eligible** genes only (not all rows)? **Proposed default:** yes: keeps
   n_eligible fixed under the null; partial genes only occur at linear contig ends.
2. BH family = (region_type, stat) rather than region_type alone (pooling O − E and HL)? **Proposed
   default:** (region_type, stat); detection and abundance answer different questions.
3. With few shifts (K < 19 the smallest p exceeds 0.05) add a flag? **Proposed default:** add
   `few_shifts` to the §9.4 "Region and context" vocabulary via the schema (needs the lead's OK, as
   §9.4 says "extend only through the schema").

## Note from M4-07 review (lead)
M4-17 adds overlapping regions and per-contig region placement to the simulator (`RegionSpec.contig`; a 60-gene plasmid with its own AMR set; a gene in two regions as in M2-20 example 4) — M4-07 does not provide them.

## Note from M4-08 review (lead)
Simulated TU expression follows an AR(1) that restarts per replicon and does not wrap at the origin of circular replicons (one uncorrelated join); this does not affect the level of the circular-shift test.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted ("+X % detection" = relative rise of the region's mean
   detection probability, capped at 0.99; `few_shifts` when K < 19; all-partial regions → NA effect/p,
   `n_perm` 0, `small_region`; HL `n_perm` counts valid shifts; bootstrap B = 20 in calibration).
2. A region whose members lie on two contigs is flagged `split_contigs` (and tested on each contig's
   own circular shifts); name added here.
3. The simulator gains `RegionSpec.contig` (M4-07 decision 6 deferred it to M4-17); overlap is covered by
   plasmid AMR genes also in `plasmid_1` — no extra overlap knob.
4. A2 (short prophage genes, gate on M4-13 p_detect) needs PR #83 merged; until then it may be marked
   xfail(strict=False) with reason "needs M4-13 (#83)".
5. The statistician reviews claims C1–C3, A2, A4 before merge.

## Lead decisions after the statistician's review (PR #89, changes requested)
6. ADR-0018 (accepted-by-khoa 2026-10-15): residualised HL, studentised O−E, BH within (type, stat).

7. Gate restructuring (within the M4-05 harness rules, no tolerance change): pool C1 over the four
   chromosome region types (`check_level_pooled`, m = 4); plasmid AMR gated at α 0.05/0.10 only; O−E at one
   ρ; A4 becomes a deterministic invariance test (bitwise-identical at copy number 1 and 20).
8. A2 xfail(strict=False) "needs M4-13 (#83)" until #83 merges; A2 gains the anti-conservative direction.
9. Tests kill the surviving mutants: `split_contigs`, HL bootstrap strata, start-before-locus_tag sort.
10. ADR-0018 residualisation (lead confirms the test-author's reading): per replicon, OLS with intercept
    of log2 abundance on log length_aa, log1p n_observable_peptides and GRAVY (§9.3 genes.tsv covariates),
    fitted on detected genes; the HL shift, its CI and p are reported on the residual scale (log2 units
    relative to the covariate expectation; report labels it "detectability-adjusted"); the O−E effect
    stays raw while its test uses the studentised sum.

## Open (lead, 12:55): A2-HL still fails after ADR-0018 (implementer 1186184)
Residualised HL rejects 0.144 at α 0.05 (full mode R=1000, job 1849427; fast 0.17/0.11), against the
statistician's 0.067 (R=300) with their residualisation. The implementer's side runs: the linear fit
roughly halves the raw bias (1000 genes 0.24 → 0.11; + (log length)² 0.09). Next (statistician): compare
their exp.py residualisation with decision 10 as implemented (covariates, fit set: detected-only vs all,
per-replicon vs pooled, imputation), then either amend decision 10 or the attack scenario; do not change
tolerances. Branch is not mergeable until A2-HL passes.
11. ADR-0018 amendment 1 (accepted-by-khoa 2026-10-16): the HL circular-shift p-value uses the studentised
    T = HL/√V at every shift (V: per-TU variance estimated once per replicon; formula in
    docs/methods/region-tests.md "Review round 2"); reported HL effect and CI unchanged. New attack null
    A2-HL-MAR (detection independent of abundance) gates it too. Resolves the "Open" item above.
12. (lead, after the test-author's b3f7d4c) Interface `tu_variance_components`, `hl_variance`,
    `studentised_hl_null` in giap.context.regions as pinned; a replicon with one TU has σ_u² = 0; V = 0
    gives T = 0; V counts as 0 below a data-scaled tolerance (relative to the residual variance, so
    fit rounding never makes T random).
13. (lead, statistician final review) Missing `tu_id` labels raise `ValueError` (C6 makes `tu_id`
    non-null); no per-gene fallback. Pinned by the statistician's
    test_atu_missing_region_tu_labels_do_not_change_the_hl_p (accepts the error).
14. Calibration records carry a run ID and commit (backlog: tests/calibration/_criteria.py writer).
