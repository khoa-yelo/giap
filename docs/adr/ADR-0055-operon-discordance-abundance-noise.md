# ADR-0055: Operon discordance under abundance-dependent Δ noise (attack D-1)
- Status: accepted-by-khoa (2026-10-24, in session: "1. Go with (a)"): option (a), σ as a function of mean abundance, with attack D-1b gated before the claim
- Date: 2026-10-24
- Spec sections: §7.4.2 (within-operon distribution), §7.4.3 (discordance test), §8.1 row "Operon
  discordance", §8.4 suite `discordance`, §8.5, §13 gate M6
- Needs Khoa: yes. Either the method changes (§7.4.2) or the claim is narrowed. A failing attack
  blocks the gate, and no threshold is widened (R3).

## Context
§7.4.2 fits one Student-t (μ, σ, ν) to Δ over all usable OP pairs, and §7.4.3 tests every
both-detected OP pair against it. Real proteomics has more noise at low abundance, so the spread of
Δ shrinks as the mean abundance of the pair grows. The PR #147 statistician review owed attack D-1
to M6-09. In D-1, Δ − μ is scaled by exp(−0.5 z), with z the standardised mean log2 abundance of
the pair, on `operon_null(1800)`. It is implemented in `_noise_by_abundance` in
tests/calibration/attacks/test_attack_operons.py.

**M6-09 measured it through the merged step** (Slurm 1884412 at 90a0e709, rerun as Slurm 1884493 at
24267fd3 with identical pooled numbers; R = 1,000 genomes).
The bounds are the §8.4 level bounds; the quartile bound treats the pairs as independent.

| check | α = 0.01 | 0.05 | 0.10 |
|---|---|---|---|
| pooled level (200 random tested pairs per genome) | 0.0076 (≤ 0.0104) | **0.0521** (≤ 0.0510) | **0.1064** (≤ 0.1013) |
| lowest-abundance quartile, cut on the planted z (123,766 pairs) | **0.0273** (≤ 0.0106) | **0.157** (≤ 0.0512) | **0.277** (≤ 0.1017) |
| same data, no D-1 (the gate null) | 0.0098 | 0.0500 | 0.0998 |

- The pooled level fails at α = 0.05 and 0.10, 4.2 % and 6.4 % above nominal, and the failure is
  confirmed in both blocks (ADR-0015).
- The lowest quartile runs at about 3× nominal, as the PR #147 oracle found (2.4–3×).
- The fast tier fails the same way: R = 100, pooled 0.110 against 0.104 at α = 0.10.
- **Mechanism:** one σ for all pairs is too narrow for noisy low-abundance pairs and too wide for
  precise high-abundance pairs. Null calls therefore land on low-abundance pairs. BH over the
  family does not repair the per-pair error.

Both checks now form one strict xfail naming this ADR (`test_attack_d1_level`). This follows the ADR-0049 pattern for the M7 gate. **The M6
discordance claim is not met** until this ADR is decided.

## Options

**(a) Model σ as a function of mean abundance (method change, §7.4.2).**
- The model is σ_i = exp(log σ₀ + γ z_i), with z_i the standardised mean log2 abundance of pair i.
  Location and ν stay as now. The fit is maximum likelihood over (μ, log σ₀, γ, ν), using the same
  Nelder–Mead routine and ν bounds as `fit_student_t`.
- The p-value is 2 · sf_t(|Δ − μ̂| / σ̂_i; ν̂).
- γ = 0 recovers the current method. A likelihood-ratio test of γ = 0 can be reported, but it is
  not used to choose between the models (no selection step).
- **Statistician check** (a scratch oracle fit outside the code, not reproducible from the repo,
  R = 150 genomes; the implementer PR must re-measure it through D-1 and D-1b at the gate's R):

  | data | pooled level, α = 0.01 / 0.05 / 0.10 | lowest quartile | γ̂ |
  |---|---|---|---|
  | D-1 data | 0.0103 / 0.0510 / 0.1004 | 0.0085 / 0.047 / 0.095 | mean −0.508 (true −0.5) |
  | the null (no D-1) | 0.0100 / 0.0507 / 0.1008 | 0.011 / 0.054 / 0.105 | mean −0.004 |

  So it **restores the level in the lowest quartile** and costs nothing measurable under the null.
- Caveat: the log-linear σ is the shape D-1 plants, so this is a correctly specified check. A
  misspecified shape (a step or a curve) has to be added as attack D-1b before the gate is
  claimed.
- **Effort:** about 60 lines of implementation:
  - a `fit_student_t_scaled(delta, z)` in `giap.stats.tdist`;
  - wiring in `giap.operons.discordance`, plus `effect` in σ̂_i units, optional;
  - the methods note, about 40 lines;
  - about 150 test lines (unit tests, the γ = 0 reduction, D-1 and D-1b).

  That is one implementer PR plus a statistician review, roughly one feature-loop cycle. The z of
  a pair uses only the two observed abundances, so the step needs no new inputs.
- With fewer than 20 usable pairs, or with a non-converged fit, the current `too_few_op_pairs`
  rule still applies.

**(b) Narrow the claim.**
- Discordance calls are claimed only under Δ noise that does not depend on abundance.
- `operons.tsv` and `operon_flags.tsv` carry a flag such as `homoscedastic_assumption`, and the
  report says the calls are anti-conservative for low-abundance pairs when noise grows at low
  abundance.
- The D-1 tests stay strict xfails, and the gate is claimed as narrowed.
- Effort is small (a flag, a sentence, the methods note), but the per-pair claim is then false
  exactly where proteomics noise is largest, and nothing in a run can tell the user whether the
  assumption holds.

**(c) Stratify the fit by abundance quartile** (a variant of (a)).
- Fit one t per quartile of pair-mean abundance, each with ≥ 20 pairs, otherwise falling back to
  the pooled fit, and test each pair against the fit of its own stratum.
- No functional form is needed, but each fit uses about a quarter of the pairs: about 115 at
  1,800 genes, close to the n ≈ 50 regime in which ν̂ sits on a bound (docs/methods/operons.md,
  PR #147 table).
- Strata boundaries also create discontinuities: a pair on either side of a boundary gets a
  different null.

## Recommendation
**(a)**, with (b)'s label as the interim until (a) merges:
- (a) is the only option that restores the per-pair level where the failure is, at no
  measurable cost under the null.
- It reduces exactly to the current method when γ̂ = 0.
- It is also the remedy the PR #147 review listed: "σ modelled as a function of mean abundance".

Gate D-1b (a non-log-linear σ, for example a step at the median abundance) together with D-1 at
the same R before the M6 discordance claim is made.

## Decision
Pending Khoa.

## Consequences
- Until the decision: the M6 gate (simulation half) is **not met for discordance**. The
  missed-operon half and §7.4.7 criteria 1–3 on the gate null pass (docs/methods/operons.md,
  "The M6 gate, simulation half").
- With (a): an implementer issue on `giap.stats.tdist` and `giap.operons.discordance`. The D-1
  xfails become plain gated tests (an XPASS of a strict xfail fails the run, which forces the
  swap), and D-1b is added. M6-05H's P(miss) takes its ratio sample from the scaled fit:
  standardised residuals mapped back at the partner's z.
- With (b): the flag goes into the §9.3 enum (an interface change, ADR-0036 style), the report
  text changes, and the xfails stay with the narrowed claim recorded in the spec sign-off.
