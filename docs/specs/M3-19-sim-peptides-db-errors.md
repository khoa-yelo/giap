# M3-19: giapsim, proteins, peptides and database-search errors

Milestone M3 · depends on: M3-18 · changes `tests/sim/giapsim/` (new `proteomics.py`) ·
simulator tooling · pure Python · estimate ≈ 360 lines (not gated by ADR-0007)

## Spec sections implemented (§11.3 "Expression and proteomics")

- "**Peptides:** tryptic digest (7–30 aa, ≤ 1 missed cleavage); peptide detectability logistic in
  peptide features and log abundance; log-normal intensities."
- "**Database-search errors:** false identifications at rate π_wrong, drawn from a shuffled
  proteome or from a related "other strain" proteome. Scores come from a two-component mixture,
  from which q-values and PEPs are computed consistently. Decoy rows (reversed sequences) carry
  null scores."
- Outputs: "`peptides.tsv`; `reference_peptides.tsv`; `proteins.faa`; `abundance.tsv`" and
  `truth/*.tsv`.
- Consumers: placement FDR calibration (M3-10/M3-21, "wrong peptides at known rate; database …
  error model"), recall gate (≥ 99.9 % of planted peptides), protein routes (M3-12/13).
- Gene-level abundance here is a plain log-normal per gene (the TU/region/AR(1) expression model
  of §11.3 belongs to M4's issue, which will replace `gene_log_abundance` behind the same
  function).

## New config (`ProteomicsConfig`, field of `SimConfig`, default `None` = M0 behaviour)

```python
digest_len: tuple[int, int] = (7, 30); missed_cleavages: int = 1     # trypsin: after K/R (Trypsin/P)
detect_logit: tuple[float, float, float] = (-2.0, 0.08, 0.9)          # intercept, per-aa length, log-abundance
pi_wrong: float = 0.01; wrong_source: Literal["shuffled", "other_strain"] = "shuffled"
score_mixture: tuple[float, float, float, float] = (3.0, 1.0, 0.0, 1.0)  # μ_correct, σ_correct, μ_wrong, σ_wrong
n_decoy_rows: Literal["match_wrong"] | int = "match_wrong"
search_db: Literal["same", "other"] = "same"
n_runs: int = 3; proteins: ProteinSimConfig | None = None             # proteins.faa: identity range for cross-strain copies
```

## Behaviour

- Peptides of detected (logistic draw) tryptic products of all non-pseudo CDS and sORFs; a
  fraction `pi_wrong / (1 − pi_wrong)` of extra **wrong** identifications is drawn from a shuffled
  proteome (per-protein residue shuffle) or the other-strain proteome (M3-18 `reference`) and
  digested the same way.
- **Scores:** correct ~ N(3, 1), wrong and decoy ~ N(0, 1). **PEP** = π f₀(s) / (π f₀(s) + (1 − π)
  f₁(s)); **q-value** = the mixture-implied FDR of the score threshold, π S₀(s) / (π S₀(s) +
  (1 − π) S₁(s)), made monotone. Decoy rows: reversed sequences, `protein_ids = rev_<locus>`,
  `is_decoy = true`, null scores.
- Per run: intensity = exp(log abundance + N(0, 0.5)), charge 2 or 3, `spectrum_id` unique per run.
- `proteins.faa`: Bakta proteins (`>locus_tag`) and, if configured, cross-strain copies at a
  chosen identity (80–99 %) under UniProt-like IDs, for the miniprot/tblastn routes.
- **Truth:** `truth/peptides.tsv`: `peptide`, `is_wrong`, `origin` (`gene` | `shuffled` |
  `other_strain` | `decoy`), `locus_tag`, `contig`, `start`, `end`, `strand` (1-based closed; NA
  for wrong), `n_loci_genome`.

## Examples

| # | Config | Expected |
|---|---|---|
| 1 | protein `MKAAKPRGGR…` digest | products split after K/R (also before P), lengths 7–30, ≤ 1 missed cleavage |
| 2 | `pi_wrong=0.0` | every peptide maps exactly to its truth locus |
| 3 | `pi_wrong=0.05`, 10⁴ peptides | wrong share within 0.05 ± 0.007 (3 SE) |
| 4 | a peptide with score s = 1.5, π = 0.05 | PEP = 0.05·φ(1.5) / (0.05·φ(1.5) + 0.95·φ(−1.5)) = `0.05` (equal densities at the midpoint) |
| 5 | q-values | non-decreasing in −score; the q of the lowest-scoring correct-dominated row ≤ 1 |
| 6 | decoy rows | `is_decoy=true`, `protein_ids` start with `rev_`, sequences are reversed tryptic peptides |
| 7 | `peptides.tsv` | passes M1-06's `read_peptide_table` + `validate_peptides` with no error-level issue |
| 8 | `wrong_source="other_strain"` with the strain config | wrong peptides differ from the genome's translation by ≥ 1 residue (or map elsewhere) |

## Acceptance criteria (test author)

- [ ] `tests/unit/test_sim_proteomics.py`: examples 1–8.
- [ ] Property: the empirical FDP at q ≤ 0.01 (vs truth `is_wrong`) over 20 seeds is ≤ 0.01 + 2 SE
  (the simulated q-values are themselves calibrated).
- [ ] Speed: 100-kb genome with peptides < 1 s; 5 Mb < 30 s (`slow`).

## Questions for the lead

1. Trypsin rule: cleave before P or not? **Proposed default:** yes (Trypsin/P, matching
   `abundance.enzyme: trypsin_p`).
2. Detectability coefficients are placeholders until M4 defines the detectability truth.
   **Proposed default:** the values above; M4 replaces them.
