# M3-11: De novo precision estimated from the genome

Milestone M3 · depends on: M3-04 (exact target and shuffled-decoy hits; synthetic inputs suffice
for unit tests) · creates `src/giap/mapping/denovo_precision.py` · pure Python · estimate
≈ 220 code lines · **statistical**: `docs/methods/denovo-precision.md`, signed off by the
statistician (§14.1)

## Spec sections implemented

- §7.2.5 "**De novo precision, estimated from the genome.** For score bin b (deciles of `score`)
  and length class ℓ: m_bℓ = share of de novo peptides that match the target translation exactly;
  r_bℓ = share that match the shuffled decoy exactly; c_ℓ = share of reference peptides
  (`--reference-peptides`, q ≤ 0.01) that match the target exactly.
  prec_b = Σ_ℓ w_bℓ (m_bℓ − r_bℓ)/(c_ℓ − r_bℓ). Here w_bℓ is the share of bin-b peptides in class
  ℓ. Clip to [0, 1]. Take the CI from 1,000 bootstrap resamples within the bin. Do not report bins
  where c_ℓ − r_bℓ < 0.05. Without reference peptides, c_ℓ defaults to 0.99 and carries the flag
  `c_assumed`."
- §7.2.5 step 4: π̂_wrong (m = 0, de novo) = "1 − the lower 95 % bound of the de novo precision
  for the peptide's score bin; 1 if that cannot be estimated" (consumed by M3-10).
- §8.1 "De novo peptide correct | score bin | genome match vs. shuffled-decoy match, by length |
  precision with bootstrap CI". §8.4 "De novo precision | known precision per bin | |bias| ≤ 0.05
  | 50 | 300". §9.3 `mapping_fdr.tsv`: "De novo precision rows carry `section = denovo_precision`
  with `score_bin`, `precision`, `ci_low`, `ci_high`, `c_value`, `flags`." §9.4 run flag
  `c_assumed`. §5.3 reference peptides (G2).

## Public interface

```python
@dataclass(frozen=True, slots=True)
class BinPrecision:
    score_bin: int; n: int; precision: float | None; ci_low: float | None; ci_high: float | None
    c_value: float; flags: tuple[str, ...]          # c_value = Σ_ℓ w_bℓ c_ℓ

def score_bins(scores: np.ndarray, n_bins: int = 10) -> np.ndarray: ...          # 0..9
def c_by_class(ref_exact: pd.DataFrame | None, classes: Sequence[Sequence[int]]) -> tuple[dict[str, float], bool]
def denovo_precision(df: pd.DataFrame, c: Mapping[str, float], c_assumed: bool,
                     rng: np.random.Generator, n_boot: int = 1000) -> list[BinPrecision]: ...
#  df columns: peptide, score, length_class, target_exact (bool), decoy_exact (bool)
def pi_wrong_denovo(bins: Sequence[BinPrecision], bin_of_peptide: np.ndarray) -> np.ndarray: ...
```

- **Bins.** Deciles of `score` over the de novo peptides entering mapping; edges from
  `np.quantile(scores, [0.1, …, 0.9])`, bins right-closed, so tied scores share a bin; fewer than
  10 bins when edges coincide. Peptides without `score` form no bin (π̂ = 1).
- **Unreportable classes.** Within a bin, classes with `c_ℓ − r_bℓ < 0.05` are dropped and the
  remaining weights renormalised; if none remain, `precision = NA` with flag `unreported` and π̂ = 1
  (Q1).
- **Bootstrap.** 1,000 resamples of the bin's peptides with replacement (`c_ℓ` fixed), percentile
  95 % CI; RNG `giap.rng.for_step("placement", seed)`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | one class; m = 0.60, r = 0.05, c = 0.99 | prec = 0.55/0.94 = `0.585106` |
| 2 | two classes, w = 0.7/0.3; class 1 as in ex. 1; class 2 m = 0.40, r = 0.10, c = 0.90 | prec = 0.7·0.585106 + 0.3·0.375 = `0.522074` |
| 3 | m = 0.02, r = 0.05, c = 0.99 | raw −0.0319 → clipped `0` |
| 4 | m = 0.99, r = 0.0, c = 0.95 | raw 1.042 → clipped `1` |
| 5 | c = 0.08, r = 0.05 in the only class | `precision = NA`, flag `unreported`, π̂ = 1 |
| 6 | no reference peptides | every c_ℓ = 0.99, every row flag `c_assumed`, run flag `c_assumed` |
| 7 | reference peptides: 200 in class `8`, 198 match the target exactly | c_8 = `0.99` (no `c_assumed`) |
| 8 | bin with precision 0.9, CI [0.85, 0.94] | π̂_wrong for its peptides = 1 − 0.85 = `0.15` |
| 9 | 100 peptides with scores 1..100 | 10 bins of 10 peptides; scores 1–10 → bin 0 |

## Edge cases

- All scores equal → one bin. A bin with one peptide: CI from the bootstrap is degenerate; report
  it (n = 1 is visible in the row).
- Reference peptides with q > 0.01 or decoys are excluded from c.
- `m_bℓ`, `r_bℓ` use **exact** matches only (I = L), never tolerant ones.

## Acceptance criteria

- [ ] `tests/unit/test_denovo_precision.py`: examples 1–9.
- [ ] Property: precision ∈ [0, 1]; invariant to row order; CI contains the point estimate for
  bins with n ≥ 20.
- [ ] Calibration, statistics-only (`tests/calibration/test_denovo_precision.py`, `calibration`,
  R from `tolerances.yaml`): per bin a known precision p_b; correct peptides match the target with
  probability c_ℓ, wrong ones with probability r; mean |bias| per bin ≤ 0.05
  (`denovo_precision_max_abs_bias`). The genome-level version with giapsim de novo errors is M3-21.
- [ ] `docs/methods/denovo-precision.md` exists (estimator, bins, bootstrap, the reportability
  rule, how π̂ feeds M3-10) and is signed off.

## Questions for the lead

1. "Do not report bins where c_ℓ − r_bℓ < 0.05": per class or per bin? **Proposed default:** drop
   the class from that bin and renormalise; the bin is unreported only when no class remains.
2. Bootstrap: resample `c_ℓ` (reference peptides) too? **Proposed default:** no, "resamples
   within the bin" — c is held fixed; its uncertainty is small with hundreds of reference
   peptides.

## Lead decisions (test-author questions)
1. A length class with no reference peptides has no c_ℓ: its bins are not reportable and their π̂ = 1
   (conservative, ADR-0010). The class still appears in the output with c_ℓ = NA.
2. Accepted: M3-10's length-class labels; `c_by_class` reads `peptide`, `q_value` (≤ 0.01 counts),
   `is_decoy`, `target_exact`; peptides with no score get a negative bin and π̂ = 1; c − r = 0.05 exactly
   is reported; bin labels need not be consecutive; π̂ = 1 − ci_low (ADR-0010: conservative bound).
3. The calibration file is `tests/calibration/test_denovo_precision_calibration.py` (pytest cannot
   import two test modules with the same basename).
