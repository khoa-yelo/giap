# M6-13: Discordance noise scale as a function of mean abundance (ADR-0055 option (a))

Milestone M6 · track S + R · **pass 2** (hardening; closes the discordance half of the M6 gate) ·
depends on: M6-01, M6-05 (+05H), M6-09 (PR #310, the D-1 attack) · estimate ≈ 60–75 code lines
(≈ 40 in `giap.stats.tdist`, ≈ 20 in `giap.operons.discordance`, ≈ 10 in `giap.steps.operons`)
(+ ≈ 150 new test lines; ≈ 20 lines of existing-test migration owned by the lead, see "Test plan")
· **statistical** (label `stats`): the statistician reviews the methods-note draft below **before
code** (§14.1). Every open statistical choice is marked **[S-n]**.

## Spec sections implemented

- **ADR-0055, option (a)**, accepted by Khoa on 2026-10-24 ("1. Go with (a)"). It amends §7.4.2:
  the scale of the within-operon t becomes σ_i = σ₀ · exp(γ · z_i).
- §7.4.2 (the within-operon distribution), §7.4.3 "Both detected" (the two-sided p),
  §8.1 row "Operon discordance" (family and correction unchanged), §8.4 suite `discordance`
  (level), §8.5 (attacks D-1 and D-1b), §9.2 test columns, §9.3 `pairs.tsv` (no column change).
- Not changed: the one-detected P(miss) branch (§7.4.3, M6-05H), the BH family (ADR-0056), the
  missed-operon family (M6-06) and `MIN_OP_PAIRS` = 20.

## The method

Notation: i indexes the usable OP pairs, which are the M6-05 fit pairs: both genes eligible, both
detected, both `abundance_log2` finite, neither imputed. Δ_i = `delta_log2`.
m_i = (`abundance_log2_a` + `abundance_log2_b`) / 2 is the observed mean log2 abundance of the
pair.

1. **Covariate z.**
   - c = mean(m) and s = SD(m) with ddof = 1, both over the **fit pairs only**, in log2 units.
   - z_i = (m_i − c) / s.
   - A tested pair that is not a fit pair (today that means a pair with an imputed abundance,
     M6-05 example 5) gets z from its own observed (imputed) m. That z is **clamped** to
     [z_min, z_max] of the fit pairs **[S-2, S-3]**.
   - m is symmetric in a and b, so z does not depend on orientation, strand or which gene is
     upstream.
   - **Invariance.** The model family {σ₀ e^(γ z)} is the same for any affine standardisation of
     m. The fitted σ̂_i and every p are therefore invariant to c and s, up to optimiser
     tolerance. The standardisation only conditions the optimiser and gives γ its units: the
     change in log σ per SD of pair-mean log2 abundance **[S-1]**.
2. **Model.** Δ_i ~ μ + σ_i · t_ν, with σ_i = σ₀ · exp(γ · z_i).
   - μ, σ₀, γ and ν are fitted jointly by maximum likelihood over the fit pairs.
   - ν ∈ [2, 50], as in `fit_student_t`.
   - γ ∈ [−2, 2] **[S-4]**. At γ = ±2 the σ ratio between z = −2 and z = +2 is e⁸, far
     outside any plausible proteomics noise.
   - γ = 0 is exactly the current model.
3. **Fit procedure** (`giap.stats.tdist.fit_student_t_scaled`).
   1. Run the constant fit `fit_student_t(Δ)` exactly as today. If it is unestimated, or
      `converged = False`, every OP pair is `too_few_op_pairs`, exactly as today, and no
      scaled fit is tried.
   2. **Identification gate.** The scaled fit is attempted only when n_fit ≥ `MIN_GAMMA_PAIRS`
      = 50 and s ≥ `MIN_ABUNDANCE_SD_LOG2` = 0.1 **[S-5]**. Otherwise the null is the constant
      fit, and the reason is `few_pairs` or `flat_abundance` (step 4).
   3. **Start and simplex.** Nelder–Mead on (μ, log σ₀, γ, ν) starts from the constant optimum
      (μ̂_c, log σ̂_c, 0, ν̂_c), with an **explicit initial simplex**: steps of 0.1 · σ̂_c in μ,
      0.1 in log σ₀, 0.1 in γ, and 0.1 · ν̂_c in ν. The ν step is taken inward when ν̂_c is on a
      bound. scipy's default 0.00025 step for a zero start would leave γ almost frozen
      **[S-6]**. Use bounds ν ∈ [2, 50] and γ ∈ [−2, 2], `maxiter` = 4,000 and `maxfev` =
      40,000 (the current 2,000 per 3 parameters, scaled).
   4. **Snap to bounds.** As `_snap_to_bound` does today, score ν at both bounds, and γ at both
      bounds, at the fitted remaining parameters. Keep a bound value only when it does not fit
      worse.
   5. **Convergence check.** The scaled fit is **accepted** only if all of these hold:
      - `result.success`;
      - all four parameters are finite;
      - min σ̂_i over the fit pairs ≥ 10⁻⁸ · max(1, range(Δ)), the existing collapse rule;
      - NLL_scaled ≤ NLL_const + 10⁻⁹ · max(1, |NLL_const|). The start is the constant optimum,
        so the scaled fit can never be worse; if it is, the optimiser failed.

      If any check fails, the null is the constant fit, with the reason `not_converged`,
      `collapsed` or `worse_than_constant`, and a WARNING. This is **not** `too_few_op_pairs`,
      because a valid null (γ = 0) exists. ADR-0055's sentence "a non-converged fit → the
      current `too_few_op_pairs` rule" is read as applying to the constant fit (Q3).
   6. **No model selection.** When the gate in step 2 passes and step 5 accepts, the scaled fit
      is used, whatever γ̂ is (ADR-0055). Two values are **reported** but never used to choose:
      the likelihood ratio LR = 2 (NLL_const − NLL_scaled) ≥ 0, and its χ²₁ p **[S-7]**.
4. **Two-sided p.** For each both-detected OP pair:
   - σ̂_i = σ̂₀ · exp(γ̂ · clamp(z_i));
   - F = F_t((Δ_i − μ̂) / σ̂_i; ν̂);
   - p_i = 2 · min(F, 1 − F).

   This is the ADR's 2 · sf_t(|Δ − μ̂| / σ̂_i; ν̂). It is exactly 1 at Δ = μ̂, and NaN in gives
   NaN out.
5. **Effect and CI** (default unchanged, Q7):
   - `stat` = Δ and `effect` = Δ − μ̂, in log2 units, as today;
   - `ci_low` and `ci_high` stay NA, with `ci_unestimated` (M6-05 lead decision 2);
   - the standardised residual (Δ − μ̂) / σ̂_i is **not** a new column.
6. **BH family**: unchanged. The `both_detected_t` rows only (ADR-0056), and `q_value` is written
   back.
7. **The one-detected P(miss) branch**: unchanged. `ratios` is still `fit_rows(...)`, the raw
   sorted Δ of the fit pairs, and `p_miss(u, a_v, ratios, v)` is called exactly as today. This
   conflicts with an ADR-0055 consequence line; see Q1.
8. **Strings**: unchanged, so there is no schema change. `method = "student_t"` and
   `null_model = "within_operon_t"` (Q8). The run log records which σ model was in force.
9. **Determinism**: no RNG. The fit sorts (Δ, m) lexicographically before fitting, so it does
   not depend on row order or `pair_id`.

### Interpretation note (for the methods note, [S-10])
The fit pairs are both-detected pairs. At low m, a large |Δ| makes one gene more likely to be
undetected, so detection truncation shrinks the Δ spread there and biases γ̂ **upward**, against
the low-abundance widening. The test is also conditional on both genes being detected and on m, so
the conditional σ(z) is the right null for the tested pairs. γ̂ is a description of that
conditional spread, not of the instrument's noise curve. Under the gate null the ADR's oracle
found mean γ̂ = −0.004, so the truncation effect in the simulator is small.

## Public interface

`src/giap/stats/tdist.py` (pure; added to `__all__`):
```python
@dataclass(frozen=True, slots=True)
class ScaledStudentT:
    loc: float; scale0: float; gamma: float; df: float      # NaN when unestimated
    n: int; converged: bool
    model: Literal["loglinear", "constant", "none"]           # the null actually in force
    reason: str | None   # None | few_pairs | flat_abundance | not_converged | collapsed | worse_than_constant | disabled
    center: float; sd: float; z_min: float; z_max: float     # standardisation of m (log2 units)
    nll: float; nll_constant: float; lr_stat: float; lr_p: float

def fit_student_t_scaled(x: ArrayLike, covariate: ArrayLike, *,
                         df_bounds=(2.0, 50.0), gamma_bounds=(-2.0, 2.0),
                         min_pairs: int = 50, min_sd: float = 0.1) -> ScaledStudentT: ...
def scaled_t_two_sided(x: ArrayLike, covariate: ArrayLike, fit: ScaledStudentT) -> F64: ...
def scale_at(covariate: ArrayLike, fit: ScaledStudentT) -> F64: ...   # σ̂_i, clamped z
```
- `covariate` is raw m. Standardisation and clamping live in the engine, so they cannot drift
  between the fit and the test.
- Pairs where either x or the covariate is not finite are dropped together.
- `model = "constant"` returns gamma = 0.0 and the `fit_student_t` loc, scale and df **bit for
  bit**.

`src/giap/operons/discordance.py`:
- `discordance(..., *, sigma_model: Literal["loglinear", "constant"] = "loglinear")`.
  `"constant"` reproduces current main bit for bit (`reason = "disabled"`).
- `Discordance` gains `scaled: ScaledStudentT`, the null used for every p. `fit: StudentT`
  stays the constant fit, which existing consumers (`_operon_gate.DiscRun.fit`, attack C) read.
- `fit_rows` is unchanged. A new `fit_pairs(pairs) -> tuple[F64, F64]` returns (Δ, m) sorted
  lexicographically.

`src/giap/steps/operons.py`:
- The fit is threaded out of `run_operons`.
- One INFO event, `operons.discordance_fit`, is emitted. A WARNING is added when `reason` is
  `not_converged`, `collapsed` or `worse_than_constant`.
- Payload (`data`, strict JSON; NaN → null via `_plain`):
  - `sigma_model` (requested), `sigma_model_used` (`loglinear` | `constant` | `none`) and
    `reason`;
  - `n_fit`, `loc_log2`, `scale0_log2`, `gamma_per_sd`, `df`, `df_at_bound` and
    `gamma_at_bound`;
  - `abundance_center_log2`, `abundance_sd_log2`, `z_min` and `z_max`;
  - `lr_stat` and `lr_p`.
- Names follow §9.2/§9.3 style: snake_case, with units in the name. No table gains a column.

### Config
- **Proposed key:** `operons.discordance_sigma`, enum {`loglinear`, `constant`}, default
  `loglinear`.
  - It is **not** a protected tolerance: `tests/calibration/tolerances.yaml` and
    `src/giap/data/calibration.yaml` are untouched.
  - It is a method switch, also useful as an M6-10 sensitivity variant later.
  - **However**, `src/giap/config/defaults.yaml` is a protected file, and
    `test_model_defaults_equal_defaults_yaml` requires the key to be added there (Q2).
- `MIN_GAMMA_PAIRS` (50), `MIN_ABUNDANCE_SD_LOG2` (0.1) and `GAMMA_BOUNDS` (−2, 2) are **module
  constants**, like `MIN_OP_PAIRS`. They are not config keys.

## Examples (input → expected output)

E1–E5 use a **given** fit: μ = −0.3, σ₀ = 0.5, γ = −0.5, ν = 5, c = 20, s = 2, z range
[−2.5, 2.0]. The expected values are computed with `scipy.stats.t.sf`; assert to 1e-6 relative.

| # | Input | Expected |
|---|---|---|
| E1 | Δ = 0.7, m = 18 (z = −1) | σ̂_i = **0.824361**, p = **0.279284** |
| E2 | Δ = 0.7, m = 22 (z = +1) | σ̂_i = **0.303265**, p = **0.0215361**: the same Δ is 13× more significant at high abundance |
| E3 | Δ = 0.7, m = 20 (z = 0) | σ̂_i = 0.5, p = **0.101939**, equal to the constant-model p at σ = 0.5 |
| E4 | Δ = −1.3, m = 18 (\|Δ − μ\| = 1.0, z = −1) | p = **0.279284** (two-sided, symmetric about μ, not about 0) |
| E5 | an imputed tested pair, Δ = 0.7, m = 26 (z = 3) | z clamped to 2.0: σ̂_i = 0.183940, p = **0.00285671** (unclamped it would be 0.000288207) |
| E6 | 30 usable OP pairs (20 ≤ n < 50), varied m | `scaled.model = "constant"`, `reason = "few_pairs"`, γ = 0; every p **bit-identical** to `student_t_two_sided(Δ, fit_student_t(Δ_fit))`; no `too_few_op_pairs` |
| E7 | 200 usable pairs, every pair with m = 20.0 exactly (s = 0) | `model = "constant"`, `reason = "flat_abundance"`; p as on current main |
| E8 | recovery: n = 2,000, z ~ N(0, 1), Δ = −0.3 + 0.5 · e^(−0.5 z) · t₅; `rng = numpy.random.default_rng(11)`, z = `rng.standard_normal(2000)` then t₅ = `rng.standard_t(5, 2000)`, m = z | `model = "loglinear"`, converged; γ̂ ∈ [−0.62, −0.38]; μ̂ ∈ [−0.35, −0.25]; σ̂₀ within ±10 % of 0.5; ν̂ ∈ [3, 9]; LR > 20 (scratch oracle at this seed: γ̂ −0.517, μ̂ −0.282, σ̂₀ 0.507, ν̂ 5.28, LR 575) |
| E9 | homoscedastic: n = 2,000, `rng = default_rng(12)`: Δ = 0.05 + 0.4 · `rng.standard_t(6, 2000)`, then m = `rng.normal(20, 2, 2000)` | \|γ̂\| ≤ 0.08; max \|p_scaled − p_constant\| over the pairs ≤ 0.02 (oracle: γ̂ 0.0046, max difference 0.0066) |
| E10 | E8's data with m replaced by 3m + 7, and separately by m in reversed pair order | γ̂ and every p equal to E8's within 1e-6 absolute (affine and order invariance) |
| E11 | E8's frame with every pair's a and b swapped and Δ negated (minus-strand orientation) | z is unchanged; μ̂ → −μ̂ and γ̂, σ̂₀, ν̂ unchanged, and p unchanged, to optimiser tolerance (1e-3 absolute; the simplex is not reflection-equivariant) |
| E12 | E8's frame with `sigma_model = "constant"` | rows and `fit` **bit-identical** to current main on the same frame; `scaled.reason = "disabled"` |
| E13 | 19 usable pairs | unchanged: every OP row `none` + `too_few_op_pairs` + `unestimated`; `scaled.model = "none"` |
| E14 | E8's frame plus one-detected rows and a `p_miss` callable | the `ratios` passed to `p_miss` equal `fit_rows(...)` (raw Δ, not rescaled); the one-detected rows equal current main |

## Edge cases

| Case | Behaviour |
|---|---|
| empty input / no OP pair | empty rows; `scaled.model = "none"`, `n = 0`; event still emitted |
| n < 5 finite values or n < 20 | unchanged (`too_few_op_pairs`); scaled fit not tried |
| 20 ≤ n < 50 | constant null, `few_pairs` (E6) |
| flat m (s < 0.1 log2), including a single distinct m | constant null, `flat_abundance` (E7) |
| ties in m (quantised abundances) | no special case; fine if s ≥ 0.1. Ties in (Δ, m) are order-free by the lexicographic sort |
| dominant ties in Δ | the constant fit collapses → `too_few_op_pairs` (unchanged) |
| the scaled fit collapses (min σ̂_i < floor; a near-zero-noise high-abundance cluster) | constant null, `collapsed`, WARNING |
| γ̂ on ±2 | accepted, `gamma_at_bound: true` (as ν on a bound, M6-01); the statistician reads it as misspecification evidence |
| ν̂ on a bound | unchanged rule (`df_at_bound: true`) |
| optimiser stops on iterations | constant null, `not_converged`, WARNING |
| NLL_scaled > NLL_const | constant null, `worse_than_constant`, WARNING |
| tested pair outside the fit's z range (imputed pairs) | z clamped (E5) |
| a fit pair with extreme m (high leverage on γ) | no special case; it is inside [z_min, z_max] by definition (risk noted for the statistician, [S-3]) |
| NaN abundance or Δ on a tested row | untested, as today |
| partial genes | ineligible, as today; they never enter c, s or the fit |
| relative scale | unchanged: `relative_scale`; no fit is used |
| minus strand / orientation | m is symmetric; only Δ's sign flips (E11) |
| circular replicon, origin-crossing pair | no special case: m and Δ come from the gene table; adjacency is M6-04's |
| one-detected rows | unchanged (E14) |
| same input twice | identical output (no RNG) |

## Calibration and the gate

R comes from `tolerances.yaml` `discordance` (100 fast / 1,000 full); that file is unchanged (R3).
All level checks go through `calibrate_levels` / `_criteria` (ADR-0015 two-block confirmation),
with the bounds of §8.4.

| Check | Data | Criterion | Status |
|---|---|---|---|
| **gate null** (existing `test_discordance_pooled_level`, `_p_uniform`, `_bh_family_fdr`) | `operon_null(1800)` | unchanged; must pass with `loglinear` as the default | gated |
| **D-1** (`test_attack_d1_level`) | Δ − μ scaled by e^(−0.5 z) | pooled level + lowest quartile, as now; **strict xfail removed** | gated, must pass |
| D-1, every quartile | same genomes | rate per quartile ≤ `level_bound(α, n_q)` for each of the 4 quartiles **[S-9]** | gated (default) |
| **D-1b-step** (new) | as D-1, but log σ shifted by +log(2)/2 below the median of the planted z and −log(2)/2 above (σ ratio 2, centred) | pooled + every quartile | gated |
| **D-1b-floor** (new) | log σ = ½ log(0.15 + e^(−z)), centred (a variance floor at high abundance: the Rocke–Durbin shape) | pooled + every quartile | gated |
| teeth (new) | D-1 with `sigma_model = "constant"` | must be a confirmed **fail** (proves the attack still bites) | gated (fast tier) |
| **power non-inferiority** (new) | `discordance_power(1800, 40, 2.0)`, both models on the same genomes | over the R replicates, the mean paired difference (loglinear − constant) of the planted pairs' per-pair rejection rate at α = 0.05 is ≥ −2 SE **[S-8]** | gated |
| power, BH recall at q 0.1 | same | recorded; it is ~0 under both models (methods note, D-3) | recorded |
| γ̂ distribution | gate null and D-1 | mean and SD of γ̂, the share of fallbacks, the LR | recorded |

- D-1b uses the same `_noise_by_abundance` machinery with a shape argument. The quartiles are
  cut on the planted z, as in D-1.
- The fast-tier runtime of the whole `discordance` suite plus the attacks must stay under 10
  minutes on 4 vCPUs (§8.4). Record it in `benchmarks/baseline.json`.
- The full tier runs as the `calibration_full.sbatch` task `discordance` on Slurm, never in the
  session allocation.

### Predicted outcome: D-1b may fail (read before implementing)
A scratch oracle was run (spec-writer, 2026-10-24, outside the repo and not reproducible from it,
like the ADR's).
- Setup: Δ = −0.3 + 0.5 · e^(h(z)) · t₅ with h centred, z ~ N(0, 1), n = 500 pairs, R = 150,
  and the four-parameter fit as specified above.
- Rates are at α = 0.01 / 0.05 / 0.10; the quartile columns are Q1 (lowest) to Q4.

| shape | model | pooled | Q1 (lowest) | Q2 | Q3 | Q4 | γ̂ mean (SD) |
|---|---|---|---|---|---|---|---|
| none (null) | constant | .0103/.0502/.1005 | .0108/.0501/.0997 | — | — | — | — |
| none (null) | loglinear | .0102/.0504/.1005 | .0111/.0510/.1004 | — | — | — | 0.005 (0.039) |
| D-1 log-linear | loglinear | .0101/.0504/.1004 | .0110/.0513/.1006 | — | — | — | −0.495 (0.043) |
| step, ratio 2 | constant | .0084/.0531/.1077 | .0142/.0938/.1846 | .0157/.0925/.1804 | .0018/.0137/.0348 | .0009/.0102/.0280 | — |
| step, ratio 2 | loglinear | .0095/.0518/.1027 | .0079/.0418/.0929 | **.0201/.0885/.1646** | .0038/.0238/.0526 | .0085/.0458/.0956 | −0.276 (0.047) |
| floor | loglinear | .0101/.0504/.1001 | .0099/.0523/.1061 | .0097/.0465/.0946 | .0102/.0458/.0945 | .0118/.0515/.1049 | −0.416 (0.048) |

- The **step arm is predicted to fail**. Q2, just below the step, runs at about 1.7× nominal,
  and the pooled rate (0.0518 / 0.1027) is above the full-R bounds (0.0510 / 0.1013).
- The **floor arm is marginal** at full R (Q1 at α = 0.10, Q4 at α = 0.01).
- The correctly specified D-1 is at level, but Q1 at α = 0.01 sits about 1.8 SE above nominal.
  Estimation error in γ̂ widens the error of σ̂_i at extreme z, so a small anti-conservative
  excess in the outer quartiles is expected. Whether the full-R quartile bound (n ≈ 124k,
  treated as independent) tolerates it must be measured.
- **Decision path** (R3: no bound is widened and D-1b is not softened). If D-1b fails at full R:
  - D-1b becomes a strict xfail naming an **ADR-0055 amendment** (the ADR-0049 pattern);
  - the claim stays NOT MET for that shape;
  - the lead takes these options to Khoa: a more flexible σ(z) (a quadratic in z or a 3-df
    natural spline in log σ: 1–2 more parameters, the same machinery); the variance-floor model
    σ_i² = σ_f² + σ₁² e^(2γz); or narrowing the claim to "smooth, monotone, log-linear-like
    abundance noise".

## Acceptance criteria

- [ ] `fit_student_t_scaled`, `scaled_t_two_sided` and `scale_at` exist with the interface
  above. `mypy` and `ruff` are clean, with no new dependency.
- [ ] E1–E14 pass as unit tests at the stated tolerances.
- [ ] Property (hypothesis): for every input, NLL_scaled ≤ NLL_const + 1e-9 · max(1,
  |NLL_const|) when `model = "loglinear"`; the LR is ≥ 0; the fit is invariant to row order
  and to an affine map of m (1e-6).
- [ ] Property: `sigma_model = "constant"` leaves `Discordance.rows` and `.fit` bit-identical to
  current main on every existing unit fixture.
- [ ] Property: the fallback reasons are exhaustive and mutually exclusive. `model` is
  `loglinear` only when n ≥ 50, s ≥ 0.1 and every step-5 check holds.
- [ ] The step emits exactly one `operons.discordance_fit` event per run with the payload keys
  above. The JSON is strict (no NaN or inf).
- [ ] `pairs.tsv`, `operons.tsv`, `operon_flags.tsv` and their schemas are unchanged
  (byte-identical schema files).
- [ ] The calibration table above holds at **full R on Slurm**: the gate null passes; D-1
  passes with the strict xfail removed; D-1b-step and D-1b-floor pass, or are converted under
  the decision path above with numbers; the teeth arm is a confirmed fail; power is
  non-inferior. Logs go under `local/logs/calibration/m6_13_<sha>/`, with the job IDs in
  PROGRESS.md.
- [ ] `docs/methods/operons.md` carries the amendment below, signed off by the statistician
  **before** the implementation PR leaves draft. After the full run, the D-1 and D-1b rows of
  the attack table and the discordance verdict are updated.
- [ ] ADR-0055's Decision section records the outcome (lead edit).

## Test plan

**Unit tier** (every commit, `-m "not slow ..."`; ≈ 110 lines, all in new files written by the
test-author):
- `tests/unit/stats/test_tdist_scaled.py`:
  - E1–E5 with a hand-built `ScaledStudentT`;
  - E8–E11 (recovery, homoscedastic, invariance, orientation);
  - the degenerate rows (n < 5, flat, collapse, γ on a bound via a planted γ = −3,
    `worse_than_constant` via a monkeypatched optimiser result);
  - the hypothesis properties.
- `tests/unit/operons/test_discordance_sigma.py`: E6, E7, E12–E14, the clamp on an imputed pair
  through `discordance`, the `scaled` field, and the bit-identity of `sigma_model = "constant"`.
- `tests/unit/steps/test_operons_step_sigma.py`: the event payload on a sim100k-sized step run.
  It is expected to report `too_few_op_pairs` / `model = "none"` there, which also counts as a
  pass (M6-08).

**Calibration tier** (≈ 40 lines; fast on the allocation at `-n 2`, full on Slurm):
- `tests/calibration/attacks/test_attack_operons.py`:
  - D-1b-step, D-1b-floor and teeth, as new test functions;
  - the shape argument for `_noise_by_abundance` (new helper `_noise_by_shape`, so the
    existing helper and M-2 are untouched);
  - the every-quartile check, as a new function reading `_d1_cached`.
- `tests/calibration/test_operon_sigma_calibration.py`: power non-inferiority, plus the γ̂
  records.

**Existing tests that change** (role guard: implementer and test-author agents do not edit
existing tests; **the lead** makes these edits in its own commit, ≈ 20 lines):
- Remove `@D1_XFAIL` from `test_attack_d1_level` and update the module docstring. A strict
  XPASS fails the run, so this edit must land in the same PR as the implementation.
- Existing `tests/unit/operons/test_discordance.py` assertions that compare against
  `student_t_two_sided(·, fit_student_t(·))`: the `chain` fixtures are random-walk abundances, so
  m varies, and the `loglinear` default changes those p.
  - Fixtures with fewer than 50 pairs fall back to the constant fit and are unaffected.
  - Larger fixtures get `sigma_model="constant"` pinned. M6-05's examples 1–3 are statements
    about the constant fit.
- The same audit covers `test_discordance_one_detected.py`, `test_attack_operons_pmiss.py`,
  `test_operons_step.py` and `_operon_gate.py` (which needs `DiscRun.scaled` for the γ̂ records).

## Methods-note amendment draft (for `docs/methods/operons.md`; statistician review before code)

*Insert after "Two-sided p" in "The within-operon distribution", and replace the first paragraph
of "Discordance" › "Test".*

> **Abundance-dependent scale (M6-13, ADR-0055 option (a)).** Proteomics noise grows at low
> abundance. One σ for every OP pair is therefore too narrow for low-abundance pairs and too wide
> for high-abundance ones. Attack D-1 measured the consequence through the step: the pooled level
> failed at α = 0.05 and 0.10, and the lowest-abundance quartile ran at about 3× nominal
> (Slurm 1884493). The scale is therefore modelled as σ_i = σ₀ · exp(γ z_i), where z_i is the
> pair's mean log2 abundance, standardised over the fit pairs.
> - μ, σ₀, γ and ν are fitted jointly by maximum likelihood, by Nelder–Mead from the constant
>   fit, with ν ∈ [2, 50] and γ ∈ [−2, 2].
> - γ = 0 is the constant model of §7.4.2, and the p of a pair is 2 · sf_t(|Δ − μ̂| / σ̂_i; ν̂).
> - The fitted σ̂_i does not depend on how z is centred or scaled. γ is reported per SD of
>   pair-mean log2 abundance.
>
> Tested pairs outside the fit's abundance range (imputed pairs) are evaluated at the nearest end
> of that range, not extrapolated.
>
> The constant model is used instead, and the run log says why, when:
> - fewer than 50 pairs are usable (γ is then poorly determined, and ν is already unidentified
>   at n = 20);
> - the pair means vary by less than 0.1 log2 (SD);
> - the four-parameter fit fails to converge, collapses, or fits worse than the constant one.
>
> There is no model-selection step: when the four-parameter fit is available it is used. The
> likelihood ratio against γ = 0 is reported, not acted on.
>
> γ̂ describes the spread of Δ **among pairs whose two genes were both detected**. Detection
> truncation narrows that spread at low abundance, so γ̂ understates the instrument's noise slope.
> That is the right null for the test, which is conditional on detection.
>
> The one-detected P(miss) still averages over the raw within-operon ratios. Under
> abundance-dependent noise that sample is too narrow for a low-abundance partner (open, M6-13
> Q1).
>
> The log-linear form is an assumption. Attack D-1b checks a step and a variance floor. *[After
> the full run: the D-1/D-1b numbers, the γ̂ under the null and under D-1, and the power
> comparison.]* The claim covers abundance-dependent noise whose log scale is close to linear in
> mean abundance; shapes that D-1b fails are named here.

*In "The M6 gate, simulation half" › "Attacks", replace the two D-1 rows, add D-1b-step,
D-1b-floor, teeth and power rows, and update the verdict line and the sign-off when the full run
lands.*

## Questions for the lead (each with a reversible default)

1. **P(miss) ratio sample.** The brief says the one-detected branch is unchanged. ADR-0055's
   Consequences says that M6-05H's P(miss) "takes its ratio sample from the scaled fit:
   standardised residuals mapped back at the partner's z".
   - The two disagree.
   - Statistically, the raw r sample is too narrow for low-abundance partners. That makes P(miss)
     smaller and so anti-conservative there.
   - **Default:** unchanged in this issue (E14). A follow-up issue, "M6-05H: P(miss) over
     r_j = μ̂ + (Δ_j − μ̂) · σ̂(z_v) / σ̂(z_j)", is filed in the backlog and measured by attack C.
     Reversing it means one extra `ratios` transform in `_one_detected`.
2. **Config key in a protected file.** `operons.discordance_sigma` needs a line in
   `src/giap/config/defaults.yaml`, which is protected (and test_config enforces equality).
   **Default:** the implementer PR ships the `sigma_model` function parameter and a module
   constant `SIGMA_MODEL = "loglinear"`, with no config key. The key lands as a lead-made,
   Khoa-approved protected-file edit (one line in `OperonsConfig` and one in `defaults.yaml`).
   It is not a tolerance either way.
3. **"Non-converged" in ADR-0055.** Read as applying to the constant fit: a failed scaled fit
   falls back to γ = 0 with a WARNING, not to `too_few_op_pairs`. **Default:** as specified.
   The alternative, "unestimated", would discard a valid null.
4. **Row-level visibility of a fallback.** **Default:** run-level event only. No §9.4 flag is
   added, because that would be an ADR-0036-style interface change. Reversing it means adding a
   flag (e.g. `sigma_constant_fallback`) through the schema.
5. **D-1b predicted failure.** The step arm is predicted to fail (table above). **Default:**
   implement as specified, measure at full R and follow the decision path. Do not pre-emptively
   change the σ model in this issue; that needs an ADR amendment and Khoa.
6. **Existing-test migration** is a lead edit (role guard). **Default:** pin
   `sigma_model="constant"` on M6-05's constant-fit assertions, and remove `@D1_XFAIL` in the
   implementation PR.
7. **`effect` and `stat` units.** **Default:** unchanged (log2). ADR-0055 calls an effect in σ̂_i
   units "optional". Reversing it changes what `effect` means in `pairs.tsv` (a documentation
   change, plus M6-07's `abundance_pattern` audit).
8. **`null_model` string.** **Default:** `within_operon_t`, unchanged; the event says which σ
   model was used. The alternative is `within_operon_t_abundance` when the log-linear model is
   used (a schema description edit).

## Open statistical choices (statistician, before code)

- **[S-1]** Standardise by mean and SD (ddof = 1) over the fit pairs. Median and MAD are an
  alternative; the p does not depend on the choice.
- **[S-2]** An imputed tested pair takes z from its imputed abundances. The alternative is
  z = 0 for any imputed pair.
- **[S-3]** Clamp z to the fit range rather than extrapolate. Leverage of extreme fit pairs is
  not controlled.
- **[S-4]** γ ∈ [−2, 2].
- **[S-5]** `MIN_GAMMA_PAIRS` = 50 and `MIN_ABUNDANCE_SD_LOG2` = 0.1. The 50 is argued from the
  M6-05 identifiability table, not measured. A recorded check is suggested: the plug-in level of
  the four-parameter fit at n ∈ {20, 50}.
- **[S-6]** Start at the constant optimum with an explicit simplex (γ step 0.1).
- **[S-7]** No selection; the LR and its χ²₁ p are reported only.
- **[S-8]** Power non-inferiority: the mean paired difference is ≥ −2 SE, with no absolute
  margin, on homoscedastic `discordance_power(1800, 40, 2.0)`. Power under D-1 data with planted
  pairs is recorded only. There, the constant model's low-abundance "power" is inflated by the
  level failure.
- **[S-9]** The quartile criterion covers every quartile, for D-1 and D-1b. The current D-1 xfail
  checks only the lowest; the alternative is lowest-only for D-1 and every quartile for D-1b.
- **[S-10]** γ̂ is conditional on both genes being detected (interpretation note above).
- **[S-11]** D-1b shapes: a ratio-2 step at the planted median, and the floor
  ½ log(0.15 + e^(−z)). Both are centred so the average log σ matches D-1.
