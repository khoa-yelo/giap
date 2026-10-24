# M6-13: Discordance noise scale as a function of mean abundance (ADR-0055 option (a))

Milestone M6 · track S + R · **pass 2** (hardening; closes the discordance half of the M6 gate) ·
depends on: M6-01, M6-05 (+05H), M6-09 (PR #310, the D-1 attack) · estimate ≈ 60–75 code lines
(≈ 40 in `giap.stats.tdist`, ≈ 20 in `giap.operons.discordance`, ≈ 10 in `giap.steps.operons`)
(+ ≈ 150 new test lines; ≈ 20 lines of existing-test migration owned by the lead, see "Test plan")
· **statistical** (label `stats`): the statistician reviews the methods-note draft below **before
code** (§14.1). Every open statistical choice is marked **[S-n]**.

**Statistician design review (2026-10-24): approved with edits.** The decisions on S-1 to S-11
and the rulings on Q1–Q8 are in "Statistician decisions" at the end. The edits are already folded
into the sections below: the separate 50-pair gate is dropped (S-5), the power criterion gets an
absolute margin (S-8), and the every-quartile checks become records with a 2α escalation (S-9).
D-1c (per-gene measurement noise) is added, and the measured table replaces the prediction. The
measurements are scratch oracles through the M6-09 simulator and the merged constant step, at
R = 1,000 (Slurm 1885387, 1885388, 1885429, 1885466 and 1885467). Logs and scripts are in
`local/logs/calibration/m6_13_statreview{,_a}/`.

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
   2. **Flat-abundance guard.** The scaled fit is attempted whenever the constant fit is
      estimated (so n_fit ≥ `MIN_OP_PAIRS` = 20) and s ≥ `MIN_ABUNDANCE_SD_LOG2` = 0.1
      **[S-5]**. Otherwise the null is the constant fit, with the reason `flat_abundance`.
      There is **no separate pair-count gate**: the statistician measured the four-parameter
      fit at level from n = 20 (S-5).
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
    reason: str | None   # None | flat_abundance | not_converged | collapsed | worse_than_constant | disabled
    center: float; sd: float; z_min: float; z_max: float     # standardisation of m (log2 units)
    nll: float; nll_constant: float; lr_stat: float; lr_p: float

def fit_student_t_scaled(x: ArrayLike, covariate: ArrayLike, *,
                         df_bounds=(2.0, 50.0), gamma_bounds=(-2.0, 2.0),
                         min_sd: float = 0.1) -> ScaledStudentT: ...   # n floor: MIN_FIT_N / the caller's MIN_OP_PAIRS
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
  - `n_fit`, `loc_log2`, `scale0_log2`, `gamma_per_sd`, `gamma_per_log2` (= γ̂ / s, the
    physical slope, comparable across runs), `df`, `df_at_bound` and `gamma_at_bound`;
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
- `MIN_ABUNDANCE_SD_LOG2` (0.1) and `GAMMA_BOUNDS` (−2, 2) are **module constants**, like
  `MIN_OP_PAIRS`. They are not config keys. There is no `MIN_GAMMA_PAIRS` (S-5).

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
| E6 | 30 usable OP pairs (20 ≤ n < 50), varied m | `scaled.model = "loglinear"`, `reason = None` (no pair-count gate, S-5); every p finite; no `too_few_op_pairs`. With `sigma_model = "constant"` every p is **bit-identical** to `student_t_two_sided(Δ, fit_student_t(Δ_fit))` |
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
| 20 ≤ n < 50 | log-linear null, as at any n ≥ 20 (E6; S-5: at level, conservative at n = 20) |
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
| D-1, every quartile | same genomes, quartiles on the planted z | per-quartile rate recorded; **escalation: a rate > 2α in any quartile fails the test** (the attack-C pattern) **[S-9]** | recorded + escalation |
| quartile null check (new) | `operon_null(1800)`, both models | per-quartile rates on the planted pair mean m **and** on the upstream gene's abundance a, recorded. This documents that the m-strata are not a valid level check (below) | recorded |
| **D-1b-step** (new) | as D-1, but log σ shifted by +log(2)/2 below the median of the planted z and −log(2)/2 above (σ ratio 2, centred) | pooled level gated as D-1; every quartile recorded + 2α escalation. **Measured to fail pooled (below): a strict xfail naming the ADR-0055 amendment from the start** | gated, strict xfail |
| **D-1b-floor** (new) | log σ = ½ log(0.15 + e^(−z)), centred (a variance floor at high abundance: the Rocke–Durbin shape) | pooled level gated; every quartile recorded + 2α escalation | gated (pooled) |
| **D-1c** (new; statistician) | per-gene measurement noise: obs_g = true_g + 0.35 · e^(−0.5 z_g) · ε_g with ε ~ N(0, 1), z_g the standardised true abundance of the detected genes; the pair columns are rebuilt from the genes. The covariate carries the same noise as Δ (errors in variables) | pooled level gated; every quartile on the upstream gene's **true** abundance recorded + 2α escalation | gated (pooled) |
| teeth (new) | D-1 with `sigma_model = "constant"` | must be a confirmed **fail** (proves the attack still bites) | gated (fast tier) |
| **power non-inferiority** (new) | `discordance_power(1800, 40, 2.0)`, both models on the same genomes | per replicate, d = (loglinear − constant) planted-pair rejection rate at α = 0.05. Pass when mean(d) − 2 SE(d) ≥ **−0.02** (absolute margin) **[S-8]**. α = 0.01 and 0.10 and effect 4.0 are recorded | gated |
| power, BH recall at q 0.1 | same | recorded; it is ~0 under both models (methods note, D-3) | recorded |
| γ̂ distribution | gate null, D-1, D-1b, D-1c | mean and SD of γ̂ and of γ̂ / s, the share of fallbacks by reason, the LR | recorded |

- D-1b uses a new `_noise_by_shape(pairs, shape)` with the D-1 machinery. The quartiles are cut
  on the planted z, as in D-1. D-1c is a gene-level transform (the `_co_impute` pattern: perturb
  the genes, then rebuild `abundance_log2_{a,b}` and `delta_log2`), seeded through the replicate
  rng.
- The 2α escalation is not a level bound and not a tolerance. It is the threshold at which a
  record fails the test so a breach cannot pass silently, as for attack C. It is a module
  constant (`ESCALATE_QUARTILE = 2.0`, ratio to α).
- The fast-tier runtime of the whole `discordance` suite plus the attacks must stay under 10
  minutes on 4 vCPUs (§8.4). Record it in `benchmarks/baseline.json`.
- The full tier runs as the `calibration_full.sbatch` task `discordance` on Slurm, never in the
  session allocation.

### Statistician measurement (2026-10-24; replaces the spec-writer's prediction)
**Setup.**
- The oracle reproduces the specified procedure: Nelder–Mead from the constant optimum with the
  explicit simplex, both bounds and the snap.
- It runs through `G.simulate(operon_null(1800))` on the D-1 stream seeds, with the merged
  constant step as the reference, at R = 1,000 genomes and about 494 tested OP pairs each.
- The pooled rate here uses every tested pair. The harness subsamples 200 per genome, so its
  estimate is noisier but has the same expectation.
- Bounds: pooled 0.0104 / 0.0510 / 0.1013 (the harness's 200 · R); quartile 0.0106 / 0.0512 /
  0.1017 (n ≈ 124k, pairs treated as independent).
- Jobs: Slurm 1885387 (noise planted on the pair mean m, quartiles on m) and 1885429 (planted on
  the upstream abundance a, quartiles on a); D-1c in 1885466 and 1885467; power and small n in
  1885388 and 1885387.

Rates at α = 0.01 / 0.05 / 0.10, loglinear unless stated otherwise:

| arm | pooled | Q1 (lowest) | Q2 | Q3 | Q4 | γ̂ mean (SD) |
|---|---|---|---|---|---|---|
| null, constant (quartiles on m) | .0100/.0505/.1002 | **.0117/.0538/.1045** | .0084/.0474/.0960 | .0090/.0472/.0958 | **.0109/.0534/.1045** | — |
| null, loglinear (quartiles on m) | .0100/.0503/.1002 | **.0115/.0536/.1048** | .0085/.0476/.0961 | .0092/.0474/.0963 | **.0109/.0528/.1035** | 0.001 (0.041) |
| null, loglinear (quartiles on a) | .0100/.0503/.1002 | .0101/.0503/.1005 | .0103/.0511/.1003 | .0099/.0504/.1006 | .0097/.0496/.0993 | 0.001 (0.041) |
| D-1, constant (teeth) | .0075/**.0519/.1065** | **.0273/.1566/.2767** | .0021/.0368/.0997 | .0005/.0114/.0392 | .0001/.0022/.0098 | — |
| D-1 | .0101/.0505/.1002 | .0093/.0478/.0969 | .0099/.0506/.0994 | .0106/.0509/.1005 | .0106/**.0525/.1038** | −0.513 (0.041) |
| D-1b-step | .0095/**.0518/.1024** | .0079/.0472/.0976 | **.0181/.0910/.1662** | .0022/.0174/.0432 | .0097/**.0515/.1024** | −0.286 (0.042) |
| D-1b-floor | .0100/.0506/.1004 | .0105/**.0529/.1053** | .0086/.0466/.0932 | .0088/.0451/.0913 | **.0122/.0576/.1115** | −0.430 (0.042) |
| D-1c, constant (quartiles on true a) | .0102/.0506/.1004 | **.0289/.1206/.2083** | .0055/.0398/.0917 | .0035/.0240/.0585 | .0027/.0177/.0426 | — |
| D-1c (quartiles on true a) | .0098/.0480/.0973 | .0077/**.0514/.1086** | .0066/.0373/.0819 | .0083/.0409/.0846 | **.0168/.0624/.1139** | −0.229 |

Bold marks a value above its bound.

**What it shows.**
- **The spec-writer's prediction is confirmed.** D-1b-step fails the pooled level (0.0518 /
  0.1024 against 0.0510 / 0.1013), and Q2 runs at 1.8× nominal. It is a strict xfail from the
  start.
- **D-1 passes** the existing gate (pooled + lowest quartile). The teeth arm is a confirmed fail.
- **D-1b-floor and D-1c pass the pooled level.** Their worst quartile is the top one, at 1.1–1.7×
  nominal (largest at α = 0.01). That is the variance floor: at high abundance the spread of Δ
  stops shrinking, and a log-linear σ keeps shrinking it.
- **The quartile check on m is not a valid level check.** It fails even for the shipped
  constant method under the homoscedastic gate null: Q1 0.0117 and Q4 0.0109 at α = 0.01.
  - The reason is that m = a + Δ/2 contains the response. Pairs with a large |Δ − μ| are pushed
    into the outer m quartiles, so E[Δ | m] ≠ μ there.
  - No σ(m) model can remove this, because the excess comes from a location shift.
  - Cut on a, the upstream gene's abundance (independent of Δ under the null), every quartile
    is at level for both models.
  - Hence S-9: every-quartile checks are records with a 2α escalation, not level gates.
- **More flexible σ(z) was measured and does not pay its way now.**
  - Quadratic, cubic and 3-df natural-spline log σ, and the variance floor
    σ² = σ_f² + σ₁² e^(2γz), all reach the same pooled level as log-linear.
  - None fixes the step. Its Q2 is 0.0157–0.0220 at α = 0.01 for every smooth model: no smooth
    curve fits a jump.
  - Under D-1c the quadratic, the floor model and an additive two-gene variance model cut the top
    quartile from 1.68× to 1.36–1.45× at α = 0.01 (the spline was not run on D-1c). Under D-1b-floor
    the quadratic, the spline and the floor model cut it from 1.22× to 1.04–1.06×.
  - Cut on a, they **cost level under the null**: Q2 / Q3 at 1.03–1.10× (quadratic, spline and
    floor alike), against log-linear's clean quartiles.
- **Decision path** (R3: no bound is widened, and D-1b is not softened):
  - D-1b-step is a strict xfail naming an **ADR-0055 amendment** (the ADR-0049 pattern), and the
    claim is NOT MET for discontinuous σ(abundance);
  - the lead takes the amendment to Khoa. It narrows the claim to the pooled level under smooth
    monotone σ(abundance), with the per-quartile excess recorded (≤ 1.7× in the measured shapes,
    against about 3× for the constant model);
  - a backlog item covers the variance-floor model, measured on D-1b-floor and D-1c, for a later
    ADR.

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
  `loglinear` only when the constant fit is estimated (n ≥ 20), s ≥ 0.1 and every step-5 check
  holds.
- [ ] The step emits exactly one `operons.discordance_fit` event per run with the payload keys
  above. The JSON is strict (no NaN or inf).
- [ ] `pairs.tsv`, `operons.tsv`, `operon_flags.tsv` and their schemas are unchanged
  (byte-identical schema files).
- [ ] The calibration table above holds at **full R on Slurm**: the gate null passes; D-1
  passes with the strict xfail removed; D-1b-floor and D-1c pass the pooled level and stay
  under the 2α quartile escalation; D-1b-step is a strict xfail naming the ADR-0055 amendment
  (if it XPASSes at full R, the statistician re-reads the measurement before the xfail is
  removed); the teeth arm is a confirmed fail; power is non-inferior at the −0.02 margin. Logs go under `local/logs/calibration/m6_13_<sha>/`, with the job IDs in
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
  - D-1b-step (strict xfail, ADR-0055 amendment), D-1b-floor, D-1c and teeth, as new test
    functions;
  - the quartile null check on `operon_null` (cuts on m and on a, recorded);
  - the shape argument for `_noise_by_abundance` (new helper `_noise_by_shape`, so the
    existing helper and M-2 are untouched);
  - the every-quartile record with the 2α escalation, as a new function reading `_d1_cached`.
- `tests/calibration/test_operon_sigma_calibration.py`: power non-inferiority, plus the γ̂
  records.

**Existing tests that change** (role guard: implementer and test-author agents do not edit
existing tests; **the lead** makes these edits in its own commit, ≈ 20 lines):
- Remove `@D1_XFAIL` from `test_attack_d1_level` and update the module docstring. A strict
  XPASS fails the run, so this edit must land in the same PR as the implementation.
- Existing `tests/unit/operons/test_discordance.py` assertions that compare against
  `student_t_two_sided(·, fit_student_t(·))`: the `chain` fixtures are random-walk abundances, so
  m varies, and the `loglinear` default changes those p.
  - With no pair-count gate (S-5), **every** fixture with ≥ 20 usable pairs and varied m is
    affected. Those fixtures get `sigma_model="constant"` pinned. M6-05's examples 1–3 are
    statements about the constant fit.
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
>   It applies from the same 20 usable pairs as the constant model. At n = 20–50 it was measured
>   at level, and more conservative than the constant fit (statistician, M6-13 S-5).
> - The fitted σ̂_i does not depend on how z is centred or scaled. γ is reported per SD of
>   pair-mean log2 abundance, and per log2 unit.
>
> Tested pairs outside the fit's abundance range (imputed pairs) are evaluated at the nearest end
> of that range, not extrapolated.
>
> The constant model is used instead, and the run log says why, when:
> - the pair means vary by less than 0.1 log2 (SD). This is a relevance guard, not an
>   identifiability one: over so narrow a range no noise curve changes σ appreciably;
> - the four-parameter fit fails to converge, collapses, or fits worse than the constant one.
>
> There is no model-selection step: when the four-parameter fit is available it is used. The
> likelihood ratio against γ = 0 is reported, not acted on.
>
> γ̂ describes the spread of Δ **among pairs whose two genes were both detected, given their
> observed mean abundance**. It is not the instrument's noise curve, for two reasons.
> - Detection truncation can narrow the spread at low abundance (in the simulator the effect is
>   not measurable: γ̂ is 0.001 under the null).
> - m carries the same measurement noise as Δ, which attenuates γ̂.
>
> The conditional spread is the right null for the test, which is conditional on detection.
>
> **What the level claim covers.** The claim is the pooled per-pair level. m = a + Δ/2 contains
> the response, so even under homoscedastic noise the pairs in the outer quartiles of m are
> enriched for large |Δ − μ| (1.05–1.17× nominal at α = 0.01, for the constant and the
> log-linear model alike). Within quartiles of an abundance that is independent of Δ (the
> upstream gene's), the null is at level.
>
> The one-detected P(miss) still averages over the raw within-operon ratios. Under
> abundance-dependent noise that sample is too narrow for a low-abundance partner (open, M6-13
> Q1).
>
> The log-linear form is an assumption. Attack D-1b checks a step and a variance floor, and D-1c
> checks per-gene measurement noise. *[After the full run: the D-1/D-1b/D-1c numbers, the γ̂
> under the null and under D-1, and the power comparison.]*
>
> The claim covers the **pooled** level under smooth, monotone abundance-dependent noise.
> - Under the variance floor and under per-gene measurement noise, the top-abundance quartile
>   still runs at up to about 1.7× nominal at α = 0.01. The constant model runs at about 3×,
>   in the lowest quartile.
> - A **discontinuous** σ (D-1b-step) fails even the pooled level. No smooth σ(z) removes that,
>   and it is outside the claim (ADR-0055 amendment).

*In "The M6 gate, simulation half" › "Attacks", replace the two D-1 rows, add D-1b-step,
D-1b-floor, D-1c, the quartile null check, teeth and power rows, and update the verdict line and the sign-off when the full run
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

## Statistician decisions (design review, 2026-10-24)

**Verdict: approved with edits** (already folded in above). The implementer may start once the
lead accepts the edits. The methods note is signed off for the design; the numbers are filled in
after the full run. Measurements: "Statistician measurement" above, R = 1,000 (Slurm 1885387,
1885388, 1885429, 1885466, 1885467; `local/logs/calibration/m6_13_statreview{,_a}/`).

| | Decision | Reason |
|---|---|---|
| S-1 | **Mean and SD (ddof = 1) over the fit pairs.** Also report γ̂ / s (`gamma_per_log2`) | The p is invariant to the choice (affine invariance). Per-SD units make γ readable within a run; per-log2 units make it comparable across runs. MAD buys nothing, because no p depends on it |
| S-2 | **An imputed tested pair takes z from its own imputed m** | z = 0 would give a low-filled pair (MinProb-style fills) the average σ, which is anti-conservative exactly where the fill puts it. Its own m is the honest input. These pairs are tested but never fitted (M6-05 example 5) |
| S-3 | **Clamp to [z_min, z_max] of the fit pairs; no leverage control** | Extrapolating e^(γz) is unsupported by data. A fit pair's influence on γ̂ is linear in z_i but bounded in its residual (the t score is bounded by ν + 1), so one extreme pair moves γ̂ by roughly z_i² / n (an order-of-magnitude bound). That is ≈ 0.07 at z = 6, n = 500, inside γ̂'s own SD of 0.04–0.05 for z ≤ 4. No winsorising |
| S-4 | **γ ∈ [−2, 2]** | Not one of 10,000 fits at n ≈ 494 reached a bound, nor any of the 24,000 fits at n = 20–100 (D-1 planted −0.5). A γ̂ on the bound is reported as misspecification evidence |
| S-5 | **No pair-count gate: the scaled fit runs from `MIN_OP_PAIRS` = 20. `MIN_ABUNDANCE_SD_LOG2` = 0.1 kept** | Measured: in-sample rates, R = 2,000 per cell; z ~ N(0, 1); Δ = −0.3 + 0.5 e^(h) t_ν; table below. The four-parameter fit is at level or conservative at every n ≥ 20, while the constant fit under D-1 at n = 20 is anti-conservative pooled (0.064 at α = 0.05) and 3× in the lowest quartile. A 50-pair gate would hand small genomes the failing method. The SD guard is about relevance, not identifiability (the fit is affine invariant); it never binds on real proteomes (s ≈ 2.3 log2 in the simulator) |
| S-6 | **Approve the start (constant optimum) and the explicit simplex** | All 10,000 oracle fits converged, none was worse than the constant fit, and none fell back |
| S-7 | **No selection; report the LR and its χ²₁ p as descriptive** | A pre-test on γ makes a pre-test estimator, with its own level distortion near the detection threshold. The χ²₁ reference is nominal only for an interior γ̂ |
| S-8 | **Non-inferiority with an absolute margin: mean(d) − 2 SE(d) ≥ −0.02** at α = 0.05 | The spec's "mean ≥ −2 SE" is a test of zero difference, not non-inferiority. It fails on the measured cost: d = −0.0060 (SE 0.0011) at α = 0.05, with planted power 0.647 → 0.641. At α = 0.01 power *rises* (+0.0027, SE 0.0004); at effect 4, −0.0008 at α = 0.05 and +0.0138 at α = 0.01. 0.02 is about 3 % of the power, so it allows the one-df cost and catches a real loss |
| S-9 | **Every quartile is recorded, with a 2α escalation that fails the test; the level gates stay pooled (all arms) + lowest quartile (D-1, existing)** | The independent-pairs quartile bound on m fails for the *shipped constant method under the gate null* (Q1 0.0117 against 0.0106), because m contains Δ/2. So it is not a valid level check, and gating it would block on an artefact no σ model can remove. Cut on the Δ-independent upstream abundance, the null is at level in every quartile. The worst measured smooth-shape quartile is 1.7× at α = 0.01 (D-1c, top quartile); 2α catches anything near the constant model's 3× |
| S-10 | **Approve, with the reading made conditional on the observed m as well** (methods draft edited) | Detection truncation is not measurable in the simulator (null γ̂ 0.001 ± 0.041). Under D-1c, γ̂ = −0.23 for a planted per-gene slope of −0.5: errors in m attenuate it, and Δ_bio dilutes it. γ̂ is a description of the conditional spread, not of the instrument |
| S-11 | **Keep both D-1b shapes, and add D-1c** (per-gene measurement noise; the most realistic generator, and the only one where the covariate carries the noise) | The step is the stress test; the floor and D-1c are the realistic shapes. The step fails pooled (pre-registered strict xfail); the floor and D-1c pass pooled with a recorded top-quartile excess |

S-5 measurement. Rates are α = 0.01 / 0.05 / 0.10; "pooled" is over all pairs, "low Q" over the
lowest quartile of z.

| n | data | constant, pooled | log-linear, pooled | constant, low Q | log-linear, low Q | γ̂ SD |
|---|---|---|---|---|---|---|
| 20 | null, ν = 5 | .0016/.0535/.1060 | .0016/.0484/.1032 | .0014/.0541/.1074 | .0008/.0404/.0954 | 0.25 |
| 20 | D-1, ν = 5 | .0044/.0637/.1172 | .0018/.0480/.1047 | .0152/.1733/.2803 | .0013/.0364/.0916 | 0.27 |
| 30 | D-1, ν = 5 | .0070/.0608/.1144 | .0061/.0514/.1023 | .0219/.1647/.2707 | .0042/.0459/.0936 | 0.21 |
| 50 | null, ν = 5 | .0090/.0519/.1024 | .0089/.0513/.1017 | .0089/.0501/.1005 | .0077/.0467/.0958 | 0.14 |
| 50 | D-1, ν = 5 | .0078/.0575/.1122 | .0088/.0512/.1018 | .0252/.1577/.2696 | .0072/.0455/.0942 | 0.15 |
| 100 | D-1, ν = 5 | .0074/.0546/.1094 | .0103/.0512/.1011 | .0260/.1591/.2755 | .0096/.0508/.0997 | 0.10 |

(ν = 30 rows are in the log; they are conservative for both models.)

**D-1b recommendation.**
- **Proceed with the log-linear σ now**, with D-1b-step a strict xfail under an **ADR-0055
  amendment** (lead → Khoa). The amendment states the narrowed claim: the pooled level under
  smooth monotone σ(abundance), with the per-quartile excess recorded.
- **Do not adopt a more flexible σ(z) in this issue.** A quadratic or a 3-df spline in log σ
  would cost about 15 lines plus 1–2 parameters. The variance floor would cost about 25 lines:
  a 5th parameter, a boundary at σ_f → 0 needing its own snap, and a second start.
  - They do not fix the step (no smooth curve can), so the step xfail and the amendment are
    needed either way.
  - They reduce the floor / D-1c top-quartile excess (1.7× → about 1.4× at α = 0.01 for D-1c,
    1.22× → about 1.05× for D-1b-floor).
  - They cost level in Δ-independent quartiles under the null (1.03–1.10×).
- **Backlog:** "σ² = σ_f² + σ₁² e^(2γz) (the variance floor), gated on D-1b-floor and D-1c
  per-quartile at the 2α escalation and on the null quartile check cut on a". If it is adopted
  later, gate the same table.

**Rulings on Q1–Q8.**
- **Q1 agree.** P(miss) is unchanged here. The backlog item must name its measurement: attack C
  (oracle and through the step) on D-1 and D-1c noise, where the raw-ratio P(miss) is expected to
  be anti-conservative for low-abundance partners.
- **Q2 agree.**
- **Q3 agree.** The fallback never fired in 10,000 fits, so the WARNING path is reachable only
  through unit fixtures.
- **Q4 agree.** The run log is enough, given how rare a fallback is. Revisit only if a real run
  falls back.
- **Q5 agree,** with the pre-registered path above: D-1b-step is written as a strict xfail from
  the start, because the measurement already shows the pooled failure.
- **Q6 agree.** Note that with S-5, every fixture with 20 or more usable pairs and varied m needs
  `sigma_model="constant"`, not only those with 50 or more.
- **Q7 agree.** The p already carries the σ̂_i standardisation; `effect` in log2 keeps M6-07's
  `abundance_pattern` meaning. The standardised residual can be derived from the event (σ̂₀, γ̂,
  c, s).
- **Q8 agree.** One string per family; the event records the σ model.

**§8.1 / §8.4 check.**
- The family is unchanged: one BH family per run over the `both_detected_t` rows (ADR-0056).
- P(miss) stays out of the family.
- Every p keeps `effect` = Δ − μ̂ with `ci_unestimated` (M6-05 lead decision 2).
- Level is gated pooled through `calibrate_levels` with the ADR-0015 two-block confirmation, and
  no bound or tolerance changes (R3).
