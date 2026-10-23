# M7-10: Class-specific false-ALT rates from paired de novo and database results

Milestone M7 · track S (statistics) · **pass 2** · depends on: M7-02, M7-05 (classes need
`mass_equivalent`), M7-04 (the wired step), and a spectrum-level de novo table (decision 1) ·
estimate ≈ 180 code lines (+ ≈ 350 test lines) · **stats** (an error-rate estimate shipped on every
ALT row; statistician reviews the methods note and the simulation before the code, §14.1)

## Spec sections implemented

- §7.5.4: "Decoys control placement only (Module 2). Most de novo errors are near-correct and land
  at the true locus, where no decoy can see them. The false-ALT rate is therefore measured where
  truth is known: **With paired results.** When `--reference-peptides` contains `spectrum_id`s that
  pair with de novo rows, every de novo peptide whose spectrum has a confident (q ≤ 0.01) reference
  identification counts as an error if the two disagree (ignoring I/L). **Rates.** Class-specific
  error rates (by `alt_type` and by `mass_equivalent`) are written as `class_error_rate` and
  `class_error_n` on each ALT row. **Without paired results,** ALT confidence carries
  `class_error_rate = unestimated`."
- §8.1 row "ALT call is real | site | de novo errors on spectra with known database identity |
  class error rate | reported, else `unestimated`"; G2 (`spectrum_id`, `--reference-peptides`);
  §7.5.5 VCF `CLASS_ERR`.

## Public interface

`src/giap/variants/error_rates.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class ClassRate:
    alt_type: str
    mass_equivalent: bool
    n: int                   # paired ALT-supporting de novo observations
    errors: int
    rate: float | None       # errors / n; None when n == 0

def paired_observations(denovo: pd.DataFrame, reference: pd.DataFrame,
                        alt_support: pd.DataFrame, q_max: float = 0.01) -> pd.DataFrame: ...
def class_rates(obs: pd.DataFrame) -> list[ClassRate]: ...
def annotate_sites(sites: pd.DataFrame, rates: Sequence[ClassRate]) -> pd.DataFrame: ...
```

## Rules

1. **Pairs.** A de novo row pairs with a reference row when `(run, spectrum_id)` are equal. The
   reference row must be a target (not decoy) with `q_value ≤ 0.01` (§7.5.4; the same `REF_MAX_Q`
   as M3-11). One reference identification per spectrum (preflight already rejects a
   `(run, spectrum_id)` with two different peptides, ADR-0005).
2. **Observations.** One observation per de novo **row** (spectrum) whose stripped peptide supports
   at least one ALT site (M7-02) and which pairs (rule 1). Its classes are the `(alt_type,
   mass_equivalent)` of each ALT it supports (a peptide supporting two ALTs contributes to both
   classes).
3. **Error.** The observation is an error when the de novo and reference peptides differ after I/L
   folding.
4. **Rate.** Per class: `rate = errors / n`, `class_error_n = n`. On every `variants.tsv` ALT row
   (and VCF `CLASS_ERR`), the rate of its class; `n = 0` or no reference table → `unestimated`.
   Frameshift rows: their own class `(frameshift, false)`.
5. **What the rate means** (methods note): when the sample truly carries a variant, the reference
   search (against a database without it) disagrees with a *correct* de novo read, so the rate is an
   upper bound on the error rate in strain-difference runs. The note says so; the gate (M7-13)
   measures the bias on simulations with planted variants.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no `--reference-peptides` | every ALT row `unestimated` |
| reference rows without `spectrum_id` | `unestimated`, warning `no_spectrum_pairing` |
| database-only input (no de novo rows) | `unestimated`: §7.5.4 measures de novo error only |
| a class with paired observations but no error | rate 0, `n` reported |
| reference q above 0.01 | not a pair |
| reference peptide equals de novo up to I/L | not an error |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 10 paired substitution observations without `mass_equivalent`, 2 disagree | rate **0.2**, n 10 on every such ALT row |
| 2 | 4 paired block observations (all `mass_equivalent`), 4 disagree | rate **1**, n 4 |
| 3 | the same spectrum: de novo `PEPTIDEK`, reference `PEPTLDEK` | not an error |
| 4 | an ALT class with no paired observation | `unestimated`, `class_error_n` NA |
| 5 | reference q 0.02 | ignored → the class may become `unestimated` |
| 6 | VCF of example 1 | `CLASS_ERR=0.2` per allele |

## Acceptance criteria

- [ ] `tests/unit/variants/test_error_rates.py`: examples 1–6 and every degenerate row.
- [ ] `tests/calibration/test_variant_class_error.py` (`calibration`, R = 20 fast / 100 full,
  decision 3): on the M7-12 scenario **without** planted variants, the class rate is unbiased for
  the true per-class error rate (|bias| ≤ 0.05, the `denovo_precision_max_abs_bias` constant); with
  planted variants, the measured upward bias is recorded, not gated.
- [ ] `docs/methods/variants.md` §"Class error rates" (rule 5); statistician sign-off.

## Questions for the lead

1. **Spectrum pairing input.** C9 `peptides.norm.parquet` and `peptide_runs.parquet` drop
   `spectrum_id`; only the raw tables keep it. **Proposed default:** preflight writes
   `work/peptide_spectra.parquet` (`run`, `spectrum_id`, `peptide`, `source`) for de novo and
   reference rows, registered as an optional contract table (C17, after M6's C16), and the
   `variants` step takes it as an optional input. Alternative: the step re-reads the raw files,
   which puts input parsing in a Module 5 step.
2. A class rate is pooled over sites, so a site's own reliability is not measured. Is per-class
   pooling enough for v1? **Proposed default:** yes (it is what §7.5.4 asks); the methods note
   says a site-level estimate needs many more paired spectra than a pilot provides.
3. No `replicates_*` key exists for a variants suite in `tests/calibration/tolerances.yaml`
   (protected). **Proposed default:** ask Khoa to add `variants: 20` (fast) / `100` (full), shared by
   this check and the M7-13 F1 gate.

## Lead decisions (2026-10-22)

1. **Spectrum pairing input:** the proposed default. Preflight writes `work/peptide_spectra.parquet`
   (`run`, `spectrum_id`, `peptide`, `source`) for de novo and reference rows, registered as optional
   contract table C17; the `variants` step takes it as an optional input. Without it, the class
   rates are `unestimated` (flag), never computed from guessed pairs.
2. **Pooling:** per-class pooling over sites is enough for v1; the methods note says why a
   site-level rate is out of reach at pilot depth.
3. **Replicates:** `replicates_fast.variants = 20` and `replicates_full.variants = 100` already exist
   (ADR-0042); the check reads them from tolerances.yaml and does not hard-code its own counts.
4. **Test-author resolutions (e7d8cb50) accepted:**
   - one observation is one de novo row per class;
   - a site's `mass_equivalent` comes from its M7-05 `flags`;
   - `no_spectrum_pairing` may be reported by preflight or by the step.
5. **C17 stays at four columns.** Preflight writes only confident target reference rows
   (q ≤ 0.01, not decoy), with `source = "reference"`, plus every de novo row (`source = "denovo"`).
   The pure `paired_observations` still filters on q and decoy status, for callers that pass raw
   tables.
6. **Nextflow:** `VARIANTS` gains a `peptide_spectra` take after `reference_peptides`, staged with
   `kv` and left out of `needed`. The lead updated the existing key sets in
   test_variants_step.py `OPTIONAL`, test_contract.py `TABLE_KEYS` and
   test_nf_variants_registry.py `TAKE`.
7. **Calibration scenario:** build from giapsim pieces now and switch to M7-12 `variant_f1` after
   #225 merges (backlog). The planted arm is diagnostic only.
