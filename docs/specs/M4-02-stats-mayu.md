# M4-02: Stats engine II: MAYU-style gene FDR (closed form)

Milestone M4 · track S · depends on: M4-01 (`monotone_q`) · estimate ≈ 140 code lines (+ ≈ 220 test
lines) · **can start now**

`poisson_binomial_sf` (§8.3) is **not** in this issue: M3-09 creates it in `src/gippa/stats/poisson_binomial.py`
(it needs it first, for E_i in §7.2.5). M4 does not use it.

## Spec sections implemented

- §8.3: "`mayu_gene_fdr(best_q, pep, sizes, n_bins)` | closed form of §7.3.2".
- §7.3.2 MAYU-style closed form (G7), quoted in full:
  - "Each gene's weight a_g is `n_observable_peptides` (floor 1). Genes are grouped into
    `mayu_size_bins` bins (deciles of a_g). For a single false peptide, the per-gene probability of
    being hit is w_b = mean(a_g in bin b) / Σ_h a_h."
  - "P(t) = the accepted peptides with q ≤ t; F(t) = Σ PEP over P(t) when PEPs exist, otherwise
    t · |P(t)|; n_b(t) = genes in bin b with ≥ 1 peptide in P(t); N_b = eligible genes in bin b."
  - "p_b = 1 − (1 − w_b)^F(t)" and "FP_b(t) = (N_b − n_b(t)) p_b / (1 − p_b); FDR(t) = Σ_b FP_b(t) /
    Σ_b n_b(t)". "q_g = min over t ≥ best_q_g of FDR(t). Detected at q ≤ 0.01."
  - "≥ 2-peptide variant. Repeat the computation with p_b⁽²⁾ = 1 − e^(−λ)(1 + λ), λ = F · w_b. It is
    reported in `detection_summary.tsv` only."
- Appendix B.3.

## Public interface

`src/gippa/stats/mayu.py`
```python
@dataclass(frozen=True, slots=True)
class MayuResult:
    gene_q: npt.NDArray[np.float64]            # per gene, NaN for genes with no accepted peptide
    curve: pd.DataFrame                         # columns: t, n_peptides, F, n_genes, fp, fdr, fdr_2pep
    bin_of_gene: npt.NDArray[np.int64]

def mayu_gene_fdr(best_q: ArrayLike, peptide_q: ArrayLike, peptide_pep: ArrayLike | None,
                  sizes: ArrayLike, n_bins: int) -> MayuResult: ...
```
- `best_q[g]`: best accepted peptide q of eligible gene g (NaN = no accepted peptide).
- `peptide_q`, `peptide_pep`: one entry per distinct accepted peptide (defines P(t) and F(t)).
- `sizes[g]`: `n_observable_peptides` of eligible gene g (floored at 1 inside).
- Bins: equal-count bins of a_g by rank (`method="first"` order = input order), `n_bins` bins, fewer
  when there are fewer genes than bins. N_b counts every eligible gene in bin b.
- Thresholds t = the distinct finite `best_q` values. `fdr_2pep` uses p_b⁽²⁾ with the same n_b, N_b.
- Degenerate cases: p_b = 1 → FP_b = N_b − n_b (FDR then clipped to 1); Σ n_b = 0 → FDR NaN.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | B.3: one bin, N_b = 5,000, n_b = 3,000, F = 50, w = 1/5,000 (all sizes equal) | p_b = 0.0099512, FP_b = 20.1024, FDR = 0.0067008 (rel tol 1e-4) |
| 2 | B.3 inputs, ≥ 2-peptide variant: λ = 0.01 | p_b⁽²⁾ = 4.9668e-5, FP = 0.099341, `fdr_2pep` = 3.3114e-5 |
| 3 | two bins: 3,000 genes with a = 100 (1,000 hit), 2,000 genes with a = 400 (1,500 hit), F = 40 | w = (9.0909e-5, 3.6364e-4); FP = (7.2863, 7.3272); FDR = 14.6135 / 2,500 = 0.0058454 |
| 4 | F from PEPs: 3 accepted peptides with PEP `[0.01, 0.02, 0.5]`, threshold covers all | F = 0.53 |
| 5 | F without PEPs: 200 peptides with q ≤ t = 0.01 | F = 0.01 · 200 = 2.0 |
| 6 | gene q: FDR curve at t = (0.001, 0.005, 0.01) = (0.004, 0.002, 0.006) | gene with best_q 0.001 → 0.002; best_q 0.005 → 0.002; best_q 0.01 → 0.006 |

## Edge cases

- Genes with `best_q` NaN still count in N_b (they are the "genes without a true peptide").
- All sizes equal → bin edges collapse; equal-count rank binning still yields `n_bins` bins of equal
  w_b, which gives the one-bin answer (example 1 holds for any `n_bins`).
- n_b = N_b → FP_b = 0. F(t) = 0 → p_b = 0 → FDR 0.
- `peptide_pep` partly NaN → `ValueError` (the tier logic decides PEP vs t·|P| for the whole run).

## Acceptance criteria

- [ ] `tests/unit/stats/test_mayu.py`: examples 1–6.
- [ ] Property: `gene_q` non-decreasing in `best_q`; increasing F never lowers FDR; FDR ∈ [0, 1].
- [ ] Monte Carlo cross-check (unit, seeded, 2,000 draws): drop F false peptides on genes with
  probability ∝ a_g; the mean number of newly hit genes among the N_b − n_b unhit genes matches
  Σ_b FP_b-model expectation (N_b − n_b)·p_b within 10 % — confirms the closed form's model, not its
  conservativeness.
- [ ] Calibration of the gene FDR built on this function is in M4-12 (suite `gene_fdr`, R = 50 / 300).

## Methods note (statistician sign-off)

`docs/methods/mayu-gene-fdr.md`: the model (false peptides land ∝ observable-peptide count), the fixed
point FP_b = (N_b − n_b + FP_b)·p_b, binning by rank, the ≥ 2-peptide variant, PEP vs t·|P|.

## Questions for the lead

1. Binning "deciles of a_g" with many ties: equal-count rank bins (ties split by input order) or
   quantile edges with merged duplicates? **Proposed default:** equal-count rank bins (deterministic,
   always `n_bins` bins when ≥ `n_bins` genes; ties carry equal a_g so w_b barely changes).
2. `mayu_gene_fdr` takes the per-peptide q and PEP arrays rather than §8.3's `(best_q, pep, sizes,
   n_bins)` (F(t) needs every accepted peptide, not one per gene). **Proposed default:** the signature
   above; §8.3's argument list is read as shorthand.

## Statistician review (PR #34, approved) — forwarded to M4-12
1. Validate inputs at the M4-12 boundary: PEP in [0, 1], finite gene sizes, integer `n_bins` ≥ 1 (a negative
   PEP would silently lower F).
2. The estimate has no CI and the base case has no margin at q = 0.01: report `fp` and `n_genes` next to
   the q threshold; warn for semi-/unspecific searches (wrong hits spread uniformly → FDP ≈ 1.24 × estimate).
3. `fdr_2pep` is not an FDR for the set of genes with ≥ 2 peptides; never label it as one.
See docs/methods/mayu-gene-fdr.md and tests/calibration/attacks/test_attack_mayu.py.
