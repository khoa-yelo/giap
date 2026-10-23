# M8-13: Gate, half 1 — GIAP's Wallenius against goseq, and GIAP's R fgsea call against a direct R fgsea run, on the golden data

Milestone M8 · track C (calibration and golden) · **pass 2** · depends on: M8-01 (Wallenius,
merged), M8-12 (#192, committed golden data, merged), M8-17 (#204 and #207 merged, #214 open),
**M8-12b** (median-exact golden GSEA inputs, proposed below) · estimate ≈ 0–20 code lines
(+ ≈ 380 test lines), see "Estimate" · **statistical** (items marked **[statistician]**)

> **Amended 2026-10-22 for ADR-0038.** GIAP no longer ports fgsea: `enrichment_gsea.tsv` comes
> from `fgsea::fgseaMultilevel` 1.36.2 run by `workflow/bin/giap_fgsea.R` (ADR-0038 items 1, 3).
> So the fgsea half of this gate cannot mean "the Python port agrees with R" any more. It now
> means **"GIAP's call into R fgsea reproduces a direct R fgsea run on the golden data"**: it
> checks the wiring (TU ranking, ties, TU sets, size filter, seed, parameters, parsing, row
> mapping), not the algorithm, which is fgsea's own by construction (ADR-0038 "Evidence"). The
> Wallenius half is unchanged in intent (ADR-0038 item 5). M8-11 (the multilevel port, #193) is
> cancelled. M8-17 lead decision 10 had marked this half "superseded by the M8-17a contract
> test". This amendment narrows that: A1 covers the R side, not `gsea_table`; see §2.3.

## What the gate says

§13 row M8: "GSEA and ORA match fgsea and goseq on golden data". §7.7.5 states the criteria and
`tests/calibration/tolerances.yaml` holds them. It is **protected and unchanged**; changing any
key needs a Khoa ADR (R3):
`golden_r: {wallenius_rel_tol: 1.0e-6, fgsea_es_abs_tol: 1.0e-9, fgsea_nes_abs_tol: 0.02,
fgsea_log10p_k: 3.0}`.

1. "Wallenius p-values against goseq run with a supplied PWF, relative difference ≤ 10⁻⁶";
2. "fgsea ES identical to 10⁻⁹";
3. "NES within 0.02";
4. "|log₁₀ p − log₁₀ p_R| within 3 × the combined `log2err` / log₂ 10".

Under ADR-0038, criteria 2–4 compare two runs of the **same** fgsea build. ES is deterministic
given the ranking, so criterion 2 still has teeth. NES and p differ only through the seed and the
tie order, so criteria 3–4 become Monte Carlo agreement checks (§3).

## Pre-implementation findings (spec-writer, 2026-10-22, on main e8ab5cd)

Each of these was a read-only check of the committed data: a few seconds of Python, no R.

- **F1. The golden GSEA scores are not exactly median-centred.** `inputs/gsea_stats.tsv.gz` has
  median −0.005705. `make_golden_inputs.py` centres the scores and then overwrites 106 of them
  with ties and zeros, accepting `|median| ≤ 0.05`. `giap.function.gsea.tu_ranking` subtracts
  the median over TUs (ADR-B). So `gsea_table` on the golden inputs ranks `score + 0.005705`, not
  the scores R saw, and every ES differs (`scoreType = "std"` weights by |score|). **An exact
  wiring check through `gsea_table` is impossible on the current data.** Remedy: M8-12b
  (question 1).
- **F2. Ties.** There are 106 tied TUs: two values × 50 TUs, plus 6 zeros. R saw them in TU-ID
  file order. GIAP sends them in its seeded tie order (`rng.permutation`, M8-02 lead decision 7),
  and fgsea's stable radix sort keeps that order (M8-17 §5, T5). So even after F1 is fixed, ES
  and the leading edges of sets whose peak falls inside a tie block legitimately differ between
  GIAP's own run and the golden run.
- **F3. The Wallenius half fails on current main [statistician].** I ran
  `giap.stats.wallenius.wallenius_ora` over the full 3,000-gene universe (`make_golden_r.R`
  passes `use_genes_without_cat = TRUE`), and the counts `numInCat` and `numDEInCat` match for
  all 150 categories. But **36 of the 300 tails exceed relative 10⁻⁶**:
  - CAT_137 over, p_R = 0.01935: rel 1.4 × 10⁻⁴;
  - about 30 tails with p_R in 0.1–0.2: rel 1–2 × 10⁻⁵;
  - CAT_012, p_R = 2.1 × 10⁻¹⁰: rel 0.05;
  - CAT_011, p_R = 3.7 × 10⁻²⁷: rel 0.18;
  - CAT_004, p_R = 2.1 × 10⁻⁷ against GIAP's 4.2 × 10⁻⁷: a factor of 2;
  - CAT_013, p_R = 0 against GIAP's 1.06 × 10⁻¹⁷, which the old rule "p_R = 0 requires
    p ≤ 10⁻³⁰⁰" fails.

  Errors of 10⁻⁵ at p ≈ 0.1 are not the deep-tail cancellation that docs/methods/function.md
  ("M8-13 golden") predicts. One side's pmf or tail sum is imprecise, and nothing here says which
  one (question 3). The test-author must reproduce these numbers before writing the test.
- **F4. The p-decade criterion cannot be met on the current data.** The old acceptance criterion
  "every decade 10⁻² … 10⁻³⁰ populated" fails. The fgsea golden p values fill decades 0–8, 10–12,
  16, 17, 26, 29 and 34. The goseq tails fill 0–3, 6, 9, 26, 45, plus one 0. The decades in
  between are empty. §6 replaces the criterion.

## 1. What half 1 checks now

### 1.1 Wallenius ORA against goseq (unchanged in intent)

`tests/golden/test_golden_wallenius.py` (layer "Golden", §11.2: every PR, CI, no R) reads
`inputs/ora_genes.tsv.gz`, `inputs/ora_categories.tsv.gz` and `outputs/goseq.tsv.gz`.

1. **Universe.** All 3,000 genes of `ora_genes` (`use_genes_without_cat = TRUE`): M = 3000,
   N = Σ detected, and `p_detect` used as the PWF verbatim.
2. **Counts.** For every category, `wallenius_ora`'s n and k equal goseq's `numInCat` and
   `numDEInCat`. Every one of the 150 categories is checked and none is skipped.
3. **Odds.** ω equals goseq's weight identity (M8-01 example 3) to relative 10⁻¹².
4. **Both tails** **[statistician]**: `|p − p_R| / p_R ≤ wallenius_rel_tol` where p_R > 0. Where
   p_R = 0 (BiasedUrn underflow, or cancellation in `1 − cdf`), the rule is
   `p ≤ wallenius_abs_floor`, defined in question 3. Until question 3 is decided it is **10⁻¹⁴**,
   GIAP's documented precision floor (docs/methods/function.md), **not** 10⁻³⁰⁰.
5. **Pipeline wiring of `p_classic`** (new, cheap). `giap.function.ora.detection_enrichment`
   runs on a gene table built from the golden ORA inputs: one gene per TU, all eligible, no
   folds. Its `p_classic` must equal `wallenius_ora` on the same category **exactly** (float `==`).
   The primary detection p is the conditional Poisson-binomial (ADR-0029) and has no goseq
   counterpart, so `p_classic` is the only ORA number the golden data can check.

### 1.2 GSEA: GIAP's R fgsea call against a direct R fgsea run

The reference is the direct R run in `make_golden_r.R`: one `set.seed(GOLDEN_SEED)` and one
`fgseaMultilevel` over the golden stats in file order, with sampleSize 101, eps 0,
nPermSimple 1000, minSize 10, maxSize 500, `scoreType = "std"` and `nproc = 1`. Its output is
`outputs/fgsea.tsv.gz`: 119 rows, because SET_027 drops below minSize. **The golden fixture**
turns the golden inputs into `gsea_table`'s inputs:
- one gene per TU (`GENE_<tu>` in TU `<tu>`), detected and eligible, with
  `abundance_log2 = score`, so the TU mean is the score exactly;
- one collection `golden` whose sets are the golden sets mapped to genes;
- the `TU_ABSENT_*` members as genes **absent** from the gene table, so that they drop exactly as
  R's name matching drops them;
- linear replicons and the positions the flags need.

It lives in `tests/helpers/` and is the test-author's to write.

The checks, from cheapest to strongest:

- **G1 call construction (CI, pure).** A capturing engine records the `EngineCall` that
  `gsea_table` builds from the fixture with `rng = for_step("function.gsea", GOLDEN_SEED)`:
  - `{tu_id: score}` of the call equals the golden stats **exactly** (float `==`, and `repr`
    round-trip through `write_call`). This needs M8-12b (F1);
  - the sorted-decreasing scores equal the golden sorted scores;
  - one job `main:golden` with min 10 and max 500. Its sets are the golden sets restricted to
    ranked TUs, in golden set order, with SET_027 absent and the other 119 present;
  - params 101 / 0 / 1000; `scoreType` is not a parameter, so it is fixed in R (A1 covers it);
  - the seed is `int` in `[1, 2³¹ − 2]` and equals the §5 derivation (the second draw after the
    permutation);
  - within each tie block, the TU order is the seeded permutation's order, not TU-ID order.
- **G2 result mapping (CI, replay).** A replay engine returns the golden `fgsea.tsv.gz` rows for
  `main:golden`. `gsea_table`'s table must then satisfy:
  - `es`, `nes`, `p_value` and `log2err` equal the golden values exactly, after parsing the
    17-digit strings;
  - `size_tus` equals `size`;
  - `leading_edge` equals the golden leading-edge TUs mapped to genes, **in fgsea's order**;
  - `q_value` equals `giap.stats.multitest.bh` over the golden p values;
  - `method = fgsea_multilevel`, `null_model = tu_set_permutation`, `n_perm = 1000`, and every
    row has the flag `ci_unestimated`.

  The written `enrichment_gsea.tsv` validates against its schema.
- **G3a exact wiring (Slurm, `slow` + `requires_r`).** `gsea_table` runs on the fixture with the
  real `RFgseaEngine` wrapped by a **normalising engine**. The wrapper does exactly two things to
  the `EngineCall` before delegating: it sets `seed = GOLDEN_SEED`, and it reorders the TUs inside
  each tie block to golden file order. Columns `pathway`, `pval`, `log2err`, `ES`, `NES`, `size`
  and `leadingEdge` of the engine's rows must equal `fgsea.tsv.gz` **as strings, row for row**,
  and G2's table rules must hold on the live table. This is the gate's wiring check. Every GIAP
  transformation runs except the two that G3a overrides, and those two are covered elsewhere:
  - seed derivation by M8-17 §5 unit tests;
  - tie order by M8-17 T5 on real R.
- **G3b GIAP's own run (Slurm, `slow` + `requires_r`) [statistician].** The same, **without**
  normalisation: GIAP's own seed and tie order. The §7.7.5 criteria are applied as Monte Carlo
  agreement, with the rules in §3:
  - **ES:** a test-side ES oracle, written from the GSEA definition with gseaParam 1 (about 20
    lines, no fgsea source), computes each set's ES under GIAP's order and under golden order.
    Then three things must hold to `fgsea_es_abs_tol`:
    - the oracle on golden order equals `calc_gsea_stat.tsv.gz`, which validates the oracle;
    - the live ES equals the oracle on GIAP's order;
    - the live ES equals the golden ES for every set where the two oracle values agree. Such a set
      is *tie-insensitive*. The others are counted and reported.
  - **NES** and **p**, per §3.
  - **Leading edges** are compared as exact sets on tie-insensitive sets (§4).

### 1.3 What the existing contract tests already cover

| Test | PR | Input | What it proves | What it does not |
|---|---|---|---|---|
| `tests/contract/test_fgsea_r_contract.py::test_golden_contract_strings` (M8-17 A1) | #204, merged | golden stats and sets written straight to `giap_fgsea.R` as one job, GOLDEN_SEED | `giap_fgsea.R` ≡ `make_golden_r.R`'s fgsea call, exact strings, 119 rows | `gsea_table`: TU collapse, centring, tie order, set construction, seed derivation |
| `…::test_golden_contract_through_wrapper` | #204 | the same, through `r_fgsea_engine` | `write_call`, the subprocess and `read_results` on golden data | the same |
| `…::test_seed_is_set_before_every_job` | #204 | the golden job behind another job | per-job `set.seed` | the same |
| `tests/contract/test_fgsea_r_switch_contract.py::test_live_engine_reproduces_the_replay_table` | #214, open | the M8-06 **simulated** step fixture (4 collections, seed 7) | live `gsea_table` + R = the replay table, byte for byte: live and recorded agree end to end | the golden data, and agreement with an *independent* direct R run (the replay was recorded through the same wrapper) |
| `…::test_function_step_on_real_r_is_byte_identical` | #214 | the same, through the `function` step | `method = fgsea_multilevel`; reruns byte-identical | the same |

So A1 already proves the R side on the golden data, and #214 proves `gsea_table` + R on
simulated data. **Nothing yet runs `gsea_table` on the golden data.** G1–G3 close that gap. G3a
is the composition of G1, G2 and A1 run live in one test. It stays because it catches what the
composition argument would miss, such as a set-order or member-order effect on fgsea's RNG
consumption.

## 2. Which tests run where

| Test | Where | Markers | Data |
|---|---|---|---|
| `tests/golden/test_golden_wallenius.py` (§1.1) | CI, every PR | none (golden layer) | committed goseq outputs |
| `tests/golden/test_golden_gsea_wiring.py::G1`, `::G2` | CI, every PR | none | committed fgsea outputs (replay) |
| `tests/golden/test_golden_gsea_wiring.py::G3a`, `::G3b` | Slurm only: `sbatch scripts/slurm/test_branch.sbatch <worktree> -m requires_r tests/golden` (≤ 2 CPUs of work, inside the 16-CPU giap-test job), and the weekly `gsea` task of the M8-15 full-calibration array (`calsuite gsea`) on main | `slow`, `requires_r` | live `giap-fgsea` env (`sbatch scripts/slurm/build_envs.sbatch fgsea`) |

- CI cannot run R (§11.4). G1 and G2 replay the committed outputs, so CI keeps a check of
  everything outside R.
- G3 skips with a reason when `real_prefix()` is `None`. That happens only off Slurm. In the
  weekly task a skip is a failure: the task asserts that the env is built first, as M8-15 does.
- The `requires_r` marker is registered where M8-17 lead decision 7 put it (`pyproject.toml` after
  #214).

## 3. Tolerances and NA rows [statistician]

All keys come from `tolerances.yaml` `golden_r` and are read, never copied. A change to any of
them needs a **Khoa ADR** (protected file, R3). This spec changes none.

| Criterion | G3a (seed and ties normalised) | G3b (GIAP's own seed and ties) |
|---|---|---|
| ES | exact strings (stronger than 1e-9) | `fgsea_es_abs_tol` = 1e-9, tie-insensitive sets; oracle equality for all |
| NES | exact strings | `fgsea_nes_abs_tol` = 0.02 where both are finite |
| p | exact strings | `|log10 p − log10 p_R| ≤ k · sqrt(e² + e_R²) / log2(10)`, k = `fgsea_log10p_k` = 3, with e the `log2err` of each side |
| NA | NA ↔ NA exactly | rules below |

**NA rules (G3b).** The golden output has **0** NA `pval`, **0** NA `NES` and **0** NA `log2err`
among its 119 rows.
- `pval` or `NES` NA on exactly one side fails the row; there is no cap.
- **`log2err` NA (`isCpGeHalf`) on exactly one side:** the row uses the finite side's `log2err`
  alone, `k · e / log2(10)`.
- **On both sides:** the row is not p-checked.
- **Recommendation for question 1 (the old Q1):** count the rows in either NA class and cap them
  at **5 % of compared sets**, i.e. ⌊0.05 × 119⌋ = **5** on the current data. They are reported
  in the JSON and never silently passed. Reasons:
  - the cap keeps a vacuous pass visible;
  - the golden side contributes none, so the cap only constrains GIAP's own run;
  - the cap is a test rule, not a `tolerances.yaml` key, so it needs no ADR.
- G3a needs no cap: NA must match NA exactly.

**Risks for the statistician (no tolerance is changed here):**
- **NES across seeds.** NES divides ES by the mean of the nPermSimple = 1000 null ES of the same
  size. The Monte Carlo sd of that mean is roughly sd(null ES)/mean · NES / √1000. For NES ≈ 2,
  that gives about 0.01 per run and about 0.015–0.02 for a difference of two runs. **The maximum
  over 119 sets can exceed 0.02.** The old port spec had the same exposure (it never matched
  R's RNG), but it was never measured. G3b reports the maximum and the 95th percentile. A miss
  is an evidence item for a Khoa ADR, not a wiring failure, because G3a has already shown the
  wiring exact.
- **Multiplicity of criterion 4.** At k = 3 per set, if `log2err` is calibrated, P(≥ 1 of 119
  exceeds) ≈ 1 − 0.9973¹¹⁹ ≈ 27 %. docs/methods/function-gsea.md notes that `log2err` "slightly
  understates" the error, which makes this worse. The run is deterministic (fixed seeds), so
  this means "about one fixed seed in four fails", not flakiness. G3b reports the number of
  exceedances.

**Default gating, proposed in question 4:** the gate is G1 + G2 + G3a, plus G3b's ES, oracle and
leading-edge rules, which are deterministic. G3b's NES and p rules are **reported** with their
§7.7.5 thresholds. A miss on them goes to the statistician and then to a Khoa ADR. It does not
block the gate.

## 4. Leading edges (old Q2: yes)

- **G2 and G3a:** the leading edge must equal the golden one as an **exact ordered list**: TUs
  in fgsea's order, mapped to genes in G2. Given the ranking and ES this is deterministic, so
  order is free to check.
- **G3b:** exact **set** equality on tie-insensitive sets, plus equality with the oracle's
  leading edge under GIAP's order on all sets. A tie-sensitive set may legitimately differ.
- It is not in §7.7.5, so it is a test rule and needs no ADR.

## 5. The report

`local/logs/golden/m8_golden.json` (G3 writes the `gsea_live` part on Slurm; CI writes the rest):

- **wallenius:** categories, max relative deviation per tail, the arg-max category, the count
  > 1e-6, the p_R = 0 rows, and the decade histogram of p_R;
- **gsea_wiring:** G1/G2 pass, number of sets 119;
- **gsea_live:**
  - G3a pass;
  - G3b max |ΔES| on tie-insensitive sets and the tie-sensitive count;
  - max and 95th percentile of |ΔNES|;
  - the number of p exceedances at k = 3 and the max of |Δlog10 p| / bound;
  - the NA `log2err` counts by class;
  - the p decade histogram;
  - the GIAP seed;
  - the R and fgsea versions from `EngineProvenance`.

The numbers and the Slurm job ID go into PROGRESS.md as gate evidence.

## Examples (what a passing run looks like)

| # | Check | Input | Expected |
|---|---|---|---|
| 1 | Wallenius counts | 150 golden categories, M = 3000 | n and k equal `numInCat` and `numDEInCat` for all 150 (they do today) |
| 2 | Wallenius tails | 300 tails | max relative deviation ≤ 1e-6 where p_R > 0. **Today 36 fail (F3); the gate is red until question 3 is resolved** |
| 3 | `p_classic` wiring | `detection_enrichment` on the golden ORA fixture | `p_classic` == `wallenius_ora` p for every category and direction |
| 4 | G1 | fixture → captured `EngineCall` | scores == golden (after M8-12b); 119 sets in golden order, SET_027 absent; 101/0/1000; [10, 500] |
| 5 | G2 | replayed `fgsea.tsv.gz` | table es/nes/p/log2err == golden; q == bh(golden p); leading edges are the golden TUs as genes |
| 6 | G3a | live R, seed → 50000, ties → file order | 119 rows equal `fgsea.tsv.gz` as strings |
| 7 | G3b ES | live R, GIAP's seed and ties | oracle(golden order) == `calc_gsea_stat` ≤ 1e-9; live == oracle(GIAP order) for 119/119; tie-insensitive sets == golden ≤ 1e-9 |
| 8 | G3b NES and p | the same | reported: max |ΔNES|, the number of p exceedances, NA counts ≤ 5 |
| 9 | mutation: scores not centred (F1 reintroduced, or `tu_ranking` skips the median) | G1 | G1 fails on the score equality; G3a fails on every ES |
| 10 | mutation: ω inverted (`category_odds` returns out/in) | §1.1 | counts pass, tails fail on all non-trivial categories |

## Degenerate cases

| Case | Rule |
|---|---|
| goseq p_R = 0 (CAT_013 under) | `p ≤ wallenius_abs_floor` (question 3; 1e-14 until decided), counted in the report |
| a category where `wallenius_ora` fails closed (`odds_out_of_range`) | the test fails. The golden ω all lie inside `ODDS_RANGE` today; a golden category must never be skipped |
| a set below minSize after name matching (SET_027, 9 `TU_ABSENT_*` members) | absent from G1's job and from both outputs, and asserted absent (M8-12 degenerate row) |
| a set covering every ranked TU | none in the golden data. Covered by M8-17 example 10 and not repeated here |
| tie blocks (2 × 50 TUs, 6 zeros) | G3a normalises them; G3b classifies sets as tie-sensitive or tie-insensitive via the oracle |
| `pval` NA on one side (unbalanced) | the row fails; the golden data have none |
| `log2err` NA | §3: one side uses the finite side; both sides are not p-checked; the class counts are capped at 5 |
| an empty leading edge | the golden data have none. The comparison is `()` == `()` (an empty cell) |
| giap-fgsea env missing | G3 skips with a reason off Slurm, and the weekly task fails |
| golden data regenerated (M8-12b) | ORA inputs and outputs must stay byte-identical. A1 (#204) and G1–G3 must pass on the new files in the same PR |

## Acceptance criteria

- [ ] **A1** `tests/golden/test_golden_wallenius.py` passes in CI. It covers examples 1, 2, 3
      and 10, with every category checked. This depends on question 3.
- [ ] **A2** `tests/golden/test_golden_gsea_wiring.py` G1 and G2 pass in CI without R (examples
      4, 5, 9).
- [ ] **A3** G3a passes on Slurm (`test_branch.sbatch -m requires_r`), and the job ID is in
      PROGRESS.md (example 6).
- [ ] **A4** G3b's deterministic rules (ES oracle, tie-insensitive ES and leading edges) pass on
      Slurm. Its NES and p numbers are recorded against the §7.7.5 thresholds (examples 7, 8).
- [ ] **A5** The JSON report (§5) exists with every field. A *vacuity floor* replaces the
      unsatisfiable decade criterion (F4):
      - ≥ 5 fgsea sets with p_R < 1e-10 (9 today);
      - ≥ 2 goseq tails with 0 < p_R < 1e-20 (2 today);
      - ≥ 1 tie-insensitive set with p_R < 1e-10.
- [ ] **A6** The mutation checks of examples 9 and 10 are documented in the PR: both injected
      faults are caught.
- [ ] **A7** PROGRESS.md records gate half 1 as met, with the max deviations and the test names.
- [ ] **A8** No change to `tolerances.yaml`, `src/` or `GIAP_SPEC.md` (a test-only issue, apart
      from M8-12b's generator change).

## Hardening owed / escalation

- **Wallenius (F3):** first decide which side is imprecise (question 3). If it is GIAP, an ADR
  on the computation (log-space or exact pmf) follows, never on `wallenius_rel_tol` (R3). If it is
  goseq or BiasedUrn, the reference itself misses 1e-6. That needs a Khoa ADR on the reference or
  on the tolerance, because no GIAP change can pass a test against an imprecise reference.
- **G3b NES or p miss after 3 iterations (R9):** evidence to the statistician, then a Khoa ADR on
  `fgsea_nes_abs_tol` / `fgsea_log10p_k` or on the multiplicity rule. Not a wiring bug once G3a
  passes.
- **The old fallback ADR** ("run R fgsea as an external tool", M8-index "Questions that need
  Khoa" 1) is moot: ADR-0038 adopted it.

## Dependencies (updated)

- M8-01 (Wallenius): merged.
- **M8-11: cancelled by ADR-0038**. PR #193 stays a reference draft (ADR-0038 item 6).
- M8-02 (`fgsea_simple`): no longer a dependency. `SIMPLE_ENGINE` is irrelevant here.
- **M8-12: merged (#192)**. Golden data in `tests/data/golden/r/`, generated by Slurm job
  1868639; the lock is `scripts/golden/golden_r.explicit.txt`.
- **M8-17: #204 (17a, env, `giap_fgsea.R`, wrapper, A1 contract) merged; #207 (17a-2,
  recording tooling) merged; #214 (17b, engine switch in `gsea_table`, `requires_r` marker in
  `pyproject.toml`) open.** G1–G3 need #214 (`engine=` in `gsea_table`). The Wallenius half needs
  nothing new.
- **M8-12b (new, proposed in question 1):** regenerate the golden GSEA inputs with an
  exactly-zero median, plus one `golden_r.sbatch` run (≤ 4 CPUs, short), a PROVENANCE update and a
  re-check of A1 (#204). G1 and G3 need it.
- M8-index row 13 should read "01, 12, 12b, 17 (#214)" instead of "01, 11, 12". That is the
  lead's edit, not part of this file.

## Estimate

- **Code:** about 0 lines in `src/`. Up to about 20 if the lead wants a public "normalising
  engine" helper; the default is test-only.
- **Tests:** about 380 lines:
  - the Wallenius golden test, about 90;
  - the golden GSEA fixture helper, about 80;
  - G1 and G2, about 110;
  - G3a and G3b with the ES oracle and the report, about 100.
- **M8-12b:** about 15 lines in `make_golden_inputs.py`, an updated `test_golden_r_inputs.py`
  expectation, and one Slurm job.
- **Wallenius investigation (question 3):** about 60 lines of an mpmath oracle in tests, ≤ 1 h of
  statistician time. Any ADR is extra.
- One PR for the tests (≤ 400 lines) and one small PR for M8-12b, which goes first.

## Questions for the lead

1. **M8-12b: regenerate the golden GSEA inputs with an exactly-zero median (F1)?**
   **Proposed default: yes.** Make the two central order statistics sum to exactly 0 in micro
   units (for example, force them to be equal to −x and +x before `_make_unique`). Keep the ties
   and the zeros. Leave the ORA inputs byte-identical (separate RNG stream or unchanged draw
   order). Rerun `golden_r.sbatch` and commit through the normal PROVENANCE path. *Alternative:*
   give up on the golden data for the wiring check and rely on A1 + #214. That leaves
   `gsea_table`'s centring and set construction checked only on simulated data.
2. **NA `log2err` rows (old Q1) [statistician].** **Proposed default:** a row NA on one side uses
   the finite side's `log2err`; a row NA on both sides is not p-checked. Both classes are counted
   and capped together at 5 % of compared sets (5 of 119). The cap is a test rule, so no ADR is
   needed. NA `pval` or `NES` on exactly one side fails with no cap.
3. **The Wallenius half fails today (F3) [statistician].** Who adjudicates, and what is the floor
   for p_R = 0? **Proposed default:** before writing the gate test, the statistician compares
   both sides with a high-precision oracle on the 36 failing tails. The oracle is the Wallenius
   pmf by mpmath quadrature from the definition (no goseq or BiasedUrn code). Then:
   - if GIAP is wrong, an ADR on the computation;
   - if goseq is wrong (likely: BiasedUrn's default `precision = 1e-7` and its own tail sum), a
     Khoa ADR: either compare against the oracle with `wallenius_rel_tol` and goseq at a looser
     documented tolerance, or change the reference.

   `wallenius_abs_floor = 1e-14` for p_R = 0 meanwhile. It is a test constant, not a
   `tolerances.yaml` key; it becomes one only by Khoa ADR.
4. **What gates in G3b [statistician].** **Proposed default:** the gate is G1 + G2 + G3a plus
   G3b's deterministic rules (ES, oracle, leading edges). G3b's NES ≤ 0.02 and p at k = 3 are
   reported and escalated, not gating, because across seeds they measure fgsea's Monte Carlo error
   rather than GIAP's wiring (§3 risks). *Alternative:* gate on them literally per §7.7.5 and
   accept a Khoa ADR if they miss.
5. **Leading edges (old Q2).** **Proposed default: yes.** An exact ordered list in G2 and G3a,
   and an exact set on tie-insensitive sets in G3b (§4).
6. **The normalising engine (G3a): test helper or `src/` helper?** **Proposed default:** a test
   helper in `tests/helpers/`. It rewrites a call only for comparison and has no pipeline use.
7. **Is a G3b failure a CI concern?** G3 never runs in CI. **Proposed default:** a weekly G3
   failure opens an issue labelled `stats` with the JSON attached, and the gate status in
   PROGRESS.md goes back to "not met" until a rerun passes or an ADR lands.

## Lead decisions (2026-10-22, on 4260941)

1. M8-12b goes first: regenerate the golden GSEA inputs with an exactly-zero median, ties and zeros
   kept, ORA files byte-identical, A1 rechecked.
2. Leading edges as proposed (ordered lists in G2/G3a, sets on tie-insensitive sets in G3b).
3. The G3a seed/tie-order wrapper lives in tests/helpers/, not src/.
4. A weekly G3 failure opens a `stats` issue with the JSON report and puts the gate back to "not
   met" in PROGRESS.md.
5. Questions 2–4 (NA `log2err`, the Wallenius disagreement, what gates) go to the statistician
   before any test is written; the Wallenius investigation comes first because it decides whether
   the gate can pass at all.
