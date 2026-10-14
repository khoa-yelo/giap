# M4-14: Gene-level abundance: iBAQ, Top3, protein-level route, transform, missing values

Milestone M4 · track A (abundance) · depends on: M4-06, M4-10 (`n_observable_peptides`), M4-11
(`accepted`, `abundance_eligible`), M4-12 (`detected` and the protein→gene map) · estimate ≈ 320 code
lines (+ ≈ 350 test lines) · **can start now** (fixtures; the M4-12 dependency is only on the
`Detection.table` and `protein_gene_map` interfaces, which tests can stub)

## Spec sections implemented

- §2.1 decision 2: "`median` by default, since these proteins were detected. `min` and `exclude` run as
  sensitivity checks; `zero` and `mean` are also accepted. Genes never identified are *not detected*
  and are never imputed."
- §2.1 decision 3: "With peptide intensities, GIAP computes gene-level iBAQ (primary) and Top3 itself,
  after mapping, from unique peptides. Protein-level values are used only when no peptide intensities
  are given: `ibaq` and `top3` count as molar; `maxlfq` and `intensity` count as relative".
- §7.3.4 quoted: "Use the accepted peptides of gene g with uniqueness `unique` (plus tolerant hits under
  `related`), with intensity aggregated as in §5.3. iBAQ_g = Σ intensity / max(1,
  `n_observable_peptides_g`) … Top3_g = the mean of the 3 most intense peptides (with fewer than 3, the
  mean of those present, flag `top3_lt3`). `abundance_raw` = iBAQ by default (`abundance.method`);
  `abundance_scale = molar`." "Protein-level input only. `abundance.tsv` joined through the protein →
  gene map. `abundance_type` = the declared type. Scale: `molar` for `ibaq` and `top3`; `relative`
  for `maxlfq`, `intensity` and `unknown`." "Transform. Linear values → log2; zeros → missing." The
  missing-value policy table (median, mean, min, zero = min − 1, exclude). "Filled values set
  `abundance_imputed = true`. Undetected genes are never imputed."
- §7.2.6: "Two input proteins on one gene → flag `two_proteins_one_gene`. With protein-level abundance
  only, the gene takes the larger value." "A cross-strain protein … carries its abundance to the strain
  gene it maps to."
- §9.3 abundance columns: `abundance_raw`, `abundance_log2`, `abundance_imputed`, `imputation_policy`,
  `abundance_type` enum{ibaq, top3, maxlfq, intensity, unknown, none}, `abundance_scale`
  enum{molar, relative}, `top3_log2`. §9.4 flags `imputed_abundance`, `top3_lt3`,
  `two_proteins_one_gene`, `cross_strain_protein`.
- §5.2 `abundance: {method: ibaq, type: unknown, scale: auto, missing: median, …}`.
- §7.3.6 unit: "the missing-value policies".

## Public interface

`src/giap/genetable/abundance.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class AbundanceResult:
    table: pd.DataFrame     # locus_tag + §9.3 abundance columns + abundance_flags (list)
    source: Literal["peptide_intensity", "protein_table", "none"]
    fill_value: float | None

def peptide_abundance(accepted: pd.DataFrame, n_observable: pd.Series) -> pd.DataFrame: ...
    # locus_tag, ibaq (linear), top3 (linear), n_quant_peptides
def protein_abundance(abundance: pd.DataFrame, pmap: pd.DataFrame, declared_type: str) -> pd.DataFrame: ...
    # locus_tag, value_log2, n_proteins, flags
def impute(log2: pd.Series, detected: pd.Series, policy: AbundanceMissing) -> tuple[pd.Series, pd.Series, float | None]: ...
    # (filled log2, imputed mask, fill value)
def gene_abundance(genes, accepted, covariates, detection, abundance, pmap,
                   cfg: GiapConfig, policy: AbundanceMissing | None = None) -> AbundanceResult: ...
```
- Source: `peptide_intensity` if any `abundance_eligible` accepted peptide has a non-NA intensity;
  else `protein_table` if C10 exists; else `none` (`abundance_type = none`, scale NA, all values NA).
- Peptide route: only `abundance_eligible` rows with intensity > 0; a peptide's intensity is the C9
  median-across-runs value (already aggregated). `abundance_raw` = iBAQ or Top3 per `abundance.method`
  (linear); `top3_log2` always when source is peptide_intensity. `abundance_type` = method;
  `abundance_scale = molar`. Top3 = arithmetic mean of the (up to) 3 largest linear intensities.
- Protein route: C10 values joined via `protein_gene_map` (M4-12); scale from C10 `scale` (or
  `abundance.scale` when not auto); linear values → log2 (`abundance_raw` keeps the input value);
  several proteins on one gene → the largest value + flag `two_proteins_one_gene`; mapping by
  `protein_hit` with `cross_strain_protein` in C12 flags → gene flag `cross_strain_protein`.
  `abundance_type = abundance.type`; scale `molar` iff type ∈ {ibaq, top3}.
- Transform: values ≤ 0 on a linear scale → missing (log2 undefined).
- Imputation (`policy` overrides `abundance.missing`, for sensitivity reruns): fill statistics over
  observed `abundance_log2` of **eligible detected** genes; filled only for `detected` genes with
  missing values; `exclude` fills nothing. `imputation_policy` = the policy on every row when the
  source is not `none`. Flag `imputed_abundance` on filled genes; `top3_lt3` when 1–2 quantified peptides.
- Undetected genes with an observed value (e.g. single-peptide genes) keep it; only imputation is
  restricted to detected genes.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | gene with unique accepted peptides of intensity 100, 300, 600, 1,000; `n_observable = 4` | iBAQ = 2,000/4 = 500 → `abundance_log2 = 8.965784`; Top3 = (1,000 + 600 + 300)/3 = 633.333 → `top3_log2 = 9.306821` |
| 2 | 2 quantified peptides 100, 300; `n_observable = 0` | iBAQ = 400/max(1, 0) = 400; Top3 = 200, flag `top3_lt3` |
| 3 | a `multi_identical` peptide (intensity 10⁶) plus a unique one (100) | only 100 used |
| 4 | tolerant accepted hit with intensity 50, `related` | included; under `same` it is not accepted at all (M4-11) |
| 5 | protein route: P1 → G (log2 20), P2 → G (log2 22), type `ibaq`, scale log2 | G = 22, flag `two_proteins_one_gene`; `abundance_scale = molar` |
| 6 | protein route, type `maxlfq` | `abundance_scale = relative` |
| 7 | observed detected log2 `[10, 12, 14, 20]`, a detected gene missing | median → 13; mean → 14; min → 10; zero → 9; exclude → stays NA |
| 8 | an undetected gene with NA abundance, policy `median` | stays NA, `abundance_imputed = false` |
| 9 | linear protein value 0 | `abundance_log2` NA (then imputed if detected) |
| 10 | no intensities and no C10 | `abundance_type = none`; all abundance columns NA; `imputation_policy` NA |

## Edge cases

- All detected genes missing → no fill statistic: fill value None, nothing imputed, warning
  `no_observed_abundance`.
- Intensities NA for some peptides of a gene: those peptides are skipped (not zeros).
- A peptide intensity of exactly 0 is skipped (it carries no quantity).
- Partial genes: values computed and reported, but they are not in the fill statistic (not eligible).

## Acceptance criteria

- [ ] `tests/unit/genetable/test_abundance_gene.py`: examples 1–10, with log2 values to 1e-6.
- [ ] Property: scaling every intensity by c > 0 shifts every `abundance_log2` (non-imputed) by
  log2(c); iBAQ ≤ Σ intensities; imputed genes ⊆ detected genes; `exclude` never sets
  `abundance_imputed`.
- [ ] On M4-09 fixtures (`abundance_table = "none"`): Spearman(`abundance_log2`, true abundance) ≥ 0.7
  over genes with ≥ 3 quantified peptides (sanity, not a gate).

## Methods note (statistician sign-off)

`docs/methods/abundance.md`: iBAQ and Top3 from unique peptides, protein-level route and molar vs
relative, transform, the five missing policies and why undetected genes are never imputed, the
sensitivity policies.

## Questions for the lead

1. Top3 averaged on the linear scale then log2 (standard) rather than the mean of log2 values?
   **Proposed default:** linear mean.
2. Fill statistic over eligible detected genes only (excluding partial genes)? **Proposed default:** yes.
3. Undetected genes keep an observed abundance in `genes.tsv` (e.g. single-peptide genes)? **Proposed
   default:** yes; statistics use detected genes only, so this is reporting, not inference.
