# M6-02: The detection-limit model — P(detect | abundance) and P(miss)

Milestone M6 · track S (statistics engine) · **pass 2 (hardening)** · depends on: nothing new
(reuses merged `giap.stats.poisson_binomial_sf`, `giap.genetable.covariates.digest`,
`giap.genetable.detectability` conventions) · estimate ≈ 300 code lines (+ ≈ 350 test lines) ·
**statistical** (label `stats`)

## Why it is pass 2

The one-detected discordance test is the only consumer, and §7.4.3 already defines the behaviour
without it: `unestimated (no_detection_limit_model)`. Pass 1 therefore ships the module with that
state; this issue turns the test on (together with M6-05H).

## Spec sections implemented

- §7.4.3 "One detected, one not" in full: the *peptide model* ("logistic regression of observed vs.
  not observed over the observable peptides of detected genes. Features: the gene's log2 abundance,
  peptide length and peptide GRAVY"), `P_det(a) = P(≥ k_min distinct peptides seen)` as "a
  Poisson-binomial over the gene's observable peptides", `k_min = 1, or 2 under two_peptide`; the
  *replicate model* ("when ≥ 2 runs exist … fitted on genes seen in some runs but not others");
  **P(miss)** = "the mean over the empirical within-operon log-ratios r of [1 − P_det,u(a_v + r)].
  The r values are oriented: Δ when u lies downstream, −Δ when upstream"; and the skip rule.
- §7.3.3 / M4-13 conventions (cross-fitting, separation handling, `detectability_constant`).
- §8.7 (inherited identification error — the model is fitted on observed peptides, which carry
  their own FDR).

## Public interface

`src/giap/operons/detection_limit.py`
```python
@dataclass(frozen=True, slots=True)
class PMissModel:
    kind: Literal["peptide", "replicate"]
    coef: Mapping[str, float]
    k_min: int
    n_fit: int
    auc: float
    usable: bool

def observable_peptides(aa: str, cfg: GiapConfig) -> list[tuple[int, float]]: ...  # (length, gravy)
def fit_peptide_model(genes, peptide_hits, cfg, rng) -> PMissModel: ...
def fit_replicate_model(genes, peptide_runs, cfg, rng) -> PMissModel: ...
def p_detect(model: PMissModel, a: float, peps: Sequence[tuple[int, float]]) -> float: ...
def p_miss(model: PMissModel, a_v: float, ratios: ArrayLike,
           peps: Sequence[tuple[int, float]]) -> float: ...
```

1. **Observable peptides** per gene: `digest(aa, cfg.abundance.enzyme, missed = 0)` filtered to
   `cfg.abundance.observable_len` (7–30), deduplicated within the gene; each keeps its length and
   GRAVY (the same GRAVY scale M4-10 uses).
2. **Peptide model.** Training rows = every observable peptide of every **detected** gene with a
   finite `abundance_log2`; the response is whether that peptide is assigned to that gene in C11
   `peptide_hits` (any route, any uniqueness — Q2). Features: gene `abundance_log2`, peptide length,
   peptide GRAVY. Logistic regression with the M4-13 separation and regularisation handling; the
   fit is **grouped by gene** when cross-fitting (`GroupKFold` by `tu_id`, as M4-13 does, so a
   partner's peptides never train the model that scores it).
3. **`p_detect(model, a, peps)`** = `poisson_binomial_sf([σ(β·x_j(a)) for j in peps], k_min)` with
   `k_min = 1`, or 2 under the `two_peptide` tier (`detection.min_peptides_sequences_only` and the
   tier in force; M4-12 already decides the tier).
4. **Replicate model.** Used when there are no peptides but ≥ 2 runs: per-run detection of a gene
   (from `peptide_runs`) regressed on the gene's mean `abundance_log2` plus the
   `detectability.covariates`, fitted on genes detected in some runs and not others. `p_detect` is
   then the model's probability directly (no Poisson-binomial).
5. **`p_miss`** = `mean_r [1 − p_detect(model, a_v + r, peps_u)]` over the oriented `ratios` of
   §7.4.3, computed by M6-05. Empty `ratios` → NaN.
6. **Usability.** `usable = False` when: no peptides and < 2 runs; < 50 training rows; the fitted
   model is constant (M4-13's `_constant` path); or the gene has no observable peptide. A
   non-usable model makes M6-05 emit `no_detection_limit_model`.
7. **RNG.** `giap.rng.for_step("operons.detection_limit", seed)` for the fold assignment only.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no peptide input and one run | `usable = False`; no fit attempted |
| peptides exist **and** ≥ 2 runs | peptide model wins (§7.4.3 "default when peptides exist"); the replicate model is not fitted (Q3) |
| gene with 0 observable peptides | `p_detect` NaN → `p_miss` NaN → M6-05 marks the pair |
| separated logistic | M4-13's handling; if it still fails, `usable = False` |
| `a_v + r` far outside the training range | the logistic extrapolates; `docs/methods/detection-limit.md` states the range and the report prints it |
| every peptide of a gene is also in another gene (`shared`) | still counted (Q2); the `shared_evidence` gene flag travels with it |
| `k_min = 2` with one observable peptide | `p_detect` = 0 → `p_miss` = 1 (the partner could never be called) |

## Examples

Toy model, `logit p_j = −12 + 0.6 · a − 0.02 · L_j`, gene u with peptides of length 10, 20, 30.

| # | Input | Expected |
|---|---|---|
| 1 | per-peptide p at a = 16 | **0.0691, 0.0573, 0.0474** |
| 2 | `p_detect(a = 17, k_min = 1)` | **0.2730** |
| 3 | `p_detect(a = 17, k_min = 2)` | **0.0281** |
| 4 | `p_miss(a_v = 17.0, ratios = [−1.0, −0.2, 0.3, 1.5])` | per-r 1 − P_det = 0.8359, 0.7524, 0.6853, 0.4820 → **0.6889** |
| 5 | the same with `ratios = [3.0]` (partner expected far higher) | small `p_miss` (< 0.25): the missing partner should have been seen |
| 6 | `poisson_binomial_sf([0.3, 0.4, 0.6], 1)` / `… , 2)` | **0.8320** / **0.3960** (the engine check) |
| 7 | fitted on 5,000 simulated peptides with known β | coefficients within 10 % and AUC ≥ 0.75 |
| 8 | a gene with no observable peptide | `p_detect` NaN, `p_miss` NaN |
| 9 | same seed twice | identical folds and coefficients |

## Acceptance criteria

- [ ] `tests/unit/operons/test_detection_limit.py`: examples 1–6, 8, 9 as hand-computed constants;
  example 7 as a `slow` test on M6-03's evidence simulator.
- [ ] Property: `p_detect` is non-decreasing in `a`, and non-increasing in `k_min`.
- [ ] Property: `p_miss` is non-increasing in `a_v`.
- [ ] Property: cross-fitting leaves no gene scored by a model trained on its own TU.
- [ ] Calibration (`detectability`-style, R = 20 fast / 100 full): with a logistic truth, the
  calibration slope of `p_detect` is 1 ± 0.2 — the same criterion M4-13 already meets.
- [ ] Integration: with M6-05H, `tables/pairs.tsv` on a simulated genome with planted undetected
  partners contains `one_detected_pmiss` rows and a single BH family.

## Hardening owed (its own pass, after the M10 pilot)

- Measure the one-detected branch's level on real data; P(miss) is a probability, not a p-value
  (the index risk register, M6-05 Q3).
- Mass-spectrometry-specific peptide features (missed cleavages, m/z range) if AUC on pilot data is
  poor; today the spec names exactly three features and this issue implements exactly those.

## Questions for the lead

1. `work/peptide_runs.parquet` is not in the M4-06 contract. **Proposed default:** add it as
   optional contract table **C16** with the M1-12 columns, so the replicate model has a validated
   input rather than an ad-hoc read.
2. "Observed" for the peptide model: any assignment of the peptide to the gene in C11 (proposed),
   or unique assignments only? **Proposed default:** any — restricting to unique would bias the
   model against genes with paralogues, which are exactly the undetected partners we score.
3. Peptide model preferred whenever peptides exist, even with many runs (proposed), or pick by
   cross-validated AUC? **Proposed default:** follow the spec literally; record both AUCs in the
   detection summary so a future ADR has evidence.

## Lead decisions (2026-10-21)
1. Q1: C16 `work/peptide_runs.parquet`, optional (M6-index decision 5).
2. Q2: "observed" = any C11 assignment of the peptide to the gene (proposed default).
3. Q3: follow the spec literally (peptide model whenever peptides exist); record both AUCs in the
   detection summary as evidence for a future ADR.
4. `p_miss` is a probability, never a p-value: it gets no q and never enters a BH family
   (M6-index decision 1). The M6-05H follow-up consumes it.
5. **Estimator (statistician review, e88e96e; docs/methods/detection-limit.md).** The model, its three
   features and its training rows are as specified, but the fit is a **conditional (truncated) likelihood**
   rather than plain logistic regression. Training only on detected genes selects on the response: the plain
   fit halves the abundance coefficient and over-predicts detection by +0.14 to +0.19, which makes P(miss)
   anti-conservative.
   - **Peptide model:** maximise Σ log-lik − Σ_g log P_det,g over training genes (finite abundance and
     ≥ k_min observed observable peptides). Fit unpenalised, with no ridge; separation or non-convergence
     means `usable = False`.
   - **Replicate model:** a binomial truncated at both ends (1 ≤ S ≤ R − 1) from 3 runs up. At exactly
     2 runs it is `usable = False` with reason `replicate_two_runs`, decided from the run count; spec
     item 6 already treats a constant model as unusable.
   - **Replicate `p_detect`:** stays the per-run q (spec item 4).
6. DL1 also **records** the calibration-in-the-large (CITL), meaning the intercept of d on logit(p̂) with
   the slope fixed at 1. Gating on it (±0.2) needs a `tolerances.yaml` key, which is Khoa's decision.
7. **For Khoa (pass 2, ADRs):**
   - widen the replicate training set to every gene seen in at least one run, which makes 2 runs usable;
   - switch the replicate `p_detect` to 1 − (1 − q)^R;
   - the DL1 CITL tolerance of ±0.2.
