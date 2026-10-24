# Free entrapment: the FDP of A1's identifications

Implements `validation/pilot/entrapment.py` (GIAP_SPEC.md §15.5 metric 8, §8.7; issue M10-07,
docs/specs/M10-07-free-entrapment.md with its binding statistician decisions). Validation only:
metric 8 gates nothing (§15.6) and is in no §8.1 testing family.

## Source

Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U. "Assessment of false discovery rate
control in tandem mass spectrometry analysis using entrapment." *Nature Methods* 22 (2025),
doi:10.1038/s41592-025-02719-x (open access, CC BY 4.0). N_T and N_E are the numbers of original
target and entrapment discoveries, and r is "the effective ratio of the entrapment to original
target database size".

- **Equation (1), the combined method.** FDP̂(T ∪ E) = N_E (1 + 1/r) / (N_T + N_E). Under an
  assumption analogous to the equal-chance assumption of target–decoy competition, the paper
  proves (Supplementary Note 2.2) that it is "an estimated upper bound, that is, on average it
  overestimates the true FDP". This is `fdp_combined`, §15.5's formula unchanged.
- **Equation (2), the lower bound.** N_E / (N_T + N_E). If every entrapment discovery is false, it
  is a lower bound on the FDP (Supplementary Note 2.1). This is `fdp_lower`.
- **Why the sample method (3) is not used.** N_E (1/r) / N_T estimates the FDP among the original
  targets only; the paper calls it "inherently flawed", so it is not computed.
- **The paired method** needs one shuffled partner per target peptide (r = 1). A free entrapment
  has no such pairs, so it does not apply.

The statement about (1) is about expectation: a single estimate can fall below the realised FDP.
The report therefore says "estimated upper bound", never "upper bound".

## Mapping onto the pilot

A1 searched UniProt K-12 plus contaminants. Against Nissle spectra its peptides split into the
original targets T, the K-12 peptides **present** in Nissle's translation (`shared`), and the
entrapment E, the K-12-only peptides (`k12_only`). The entrapment is free: it is already in A1's
database, and DIA-NN competed and filtered T and E together without knowing which is which. So the
combined method's quantity is the FDP of A1's accepted set after removing contaminant-sequence and
`not_in_universe` discoveries. r is computed in silico: r = |k12_only| / |shared|.

## The universe and "present"

- **Digest.** D(F) = the distinct `trypsin_p` products of 7–30 aa with ≤ 1 missed cleavage
  (initiator M kept, `giap.mapping.digest`), I → L, with every product holding a residue outside
  the 20 standard letters dropped (`non_standard`). The universe is D(K-12) \ D(contaminants):
  a peptide a contaminant also yields could be a true contaminant identification, so it is removed
  on both sides.
- **Present.** One function, applied unchanged to the in-silico universe and to the discoveries,
  with one index. Rule 1 (`exact`): an exact match in Nissle's six-frame translation
  (`match_exact`; I → L, both strands, the origin overlap of circular replicons), or (`strain_db`)
  a substring of a Bakta protein after I → L, which covers N-terminal products of non-ATG starts.
  Rule 2 (`mass_equivalent`): a tolerant target row (`match_tolerant`, GIAP's `mapping.tolerant`
  defaults, ≥ 10 aa) whose `diff_types` is non-empty and only `block`, the exact-composition
  exchanges of §7.2.4. `near_block`, `swap`, `substitution` and `frameshift` rows do not make a
  peptide present. Using one function on both sides keeps every asymmetry of the rules (for
  example that 7–9-aa peptides can only be present through rule 1) the same in r and in N_T, N_E.
- **Distinctness.** I/L-collapsed on both sides and at both levels. Two K-12 library entries that
  differ only by I/L give one collapsed element two chances of a false match; this is rare.

## The searchable space

Precursors carry carbamidomethyl on every C and 0 or 1 oxidised M (each oxidised position its own
form), at charge 1–4 with 300 ≤ m/z ≤ 1800, m/z = (M + z · 1.007276) / z and M the residue masses
plus 18.010565. At peptide level these bounds remove nothing: a 7–30-aa peptide weighs between
7 G + H₂O and 30 W + H₂O, carbamidomethyl-C and oxidised M lie inside the G–W residue range, and
charges 1 and 4 together already cover that whole mass range. The bounds are still applied (they
are parameters), and at precursor level they decide how many charge states a peptide has.
r_prec = Σ_{k12_only} n_prec / Σ_shared n_prec is taken from these configured bounds, not from
DIA-NN's library.

## Levels and the assessment rule

`q_value` is DIA-NN's experiment-wide precursor `Global.Q.Value`, so the set of distinct accepted
precursors (`norm(modified_peptide)`, charge), pooled over the wild-type runs, is the only unit
whose FDR the search claims to control at α = 0.01. The assessment is issued only there:

- `evidence_of_failure` if `fdp_lower` > α (checked first);
- `consistent_with_control` if `fdp_combined` ≤ α (inclusive);
- `inconclusive` otherwise; `not_evaluable` when there are no discoveries or no shared peptides.

Both comparisons are exact rationals (`fractions.Fraction`, α = 1/100): `fdp_combined` is
computed as N_E (k + s) / (k (N_T + N_E)) with k = |k12_only| and s = |shared|, which equals
equation (1) because 1 + 1/r = (k + s)/k. These are the paper's three readings (Fig. 1b).

The peptide level is §15.5's number but is **not assessed** (`not_assessed`): DIA-NN claims no FDR
for the distinct peptides derived from its precursor list, and false precursors spread over
distinct peptides while true ones cluster, so a peptide-level FDP above α is not evidence that the
search failed. Per-run blocks are descriptive and not assessed either, because the rows were
selected by the global q-value, not a run-level one. Estimates are reported unclipped and flagged
`above_one` when equation (1) exceeds 1, which small counts allow.

## Uncertainty

- **`fdp_lower` is a realised bound, not an estimate.** If every entrapment discovery is false,
  the realised FDP of this list is at least `fdp_lower`, with no sampling error. The condition is
  that the genome is complete and is the sample's strain, and that no K-12 material is in the
  sample. The report calls it a realised lower bound, conditional on these caveats.
- **`ucb95`**, a one-sided 95 % upper confidence bound for the FDP. Model: the equal chance
  assumption in its plainest form, with independent false matches. Given M false discoveries among
  N = N_T + N_E, N_E ~ Binomial(M, π) with π = k/(k + s) = r/(1 + r); then N_E/π estimates M
  without bias, and `fdp_combined` = N_E/(π N). The bound inverts the exact one-sided test:
  M_max = max{M ∈ {N_E, …, N} : BinomCDF(N_E; M, π) > 0.05} and `ucb95` = M_max / N. It bounds the
  realised FDP M/N, which lies on the grid {0, 1/N, …, 1}, so it is at most 1 and can sit below an
  unclipped `fdp_combined` (above one, or π close to 1); it is never below min(`fdp_combined`, 1)
  by a full grid step 1/N. It is not valid under clustered (protein-level) errors, and no
  two-sided interval is given. At precursor level π comes from the precursor sums.
- A protein-block bootstrap was rejected: the discoveries are one realised list, and the
  resampling unit for false matches is not the protein of the discovered peptides.

## Bias

| Source | Direction |
|---|---|
| Near-homolog entrapment: every K-12-only peptide near a Nissle variant has its neighbour's signal in the sample, so neighbour-driven errors concentrate on E (§15.5 "errs high") | high |
| 7–9-aa variant peptides are not tolerantly matched, so they stay in E | high |
| The tolerant route reports only minimal-cost loci; a block-only locus behind a cheaper non-block locus is not seen | high, negligible |
| An incomplete genome, a different strain, or K-12 material in the sample (carryover, lab contamination) makes some K-12-only discoveries true; the lower bound then stops being a bound | high |
| Six-frame chance hits move short K-12-only peptides into `shared` | lower N_E (a more conservative lower bound); negligible for equation (1) |
| E and T differ in length and charge composition, so equal chance per precursor mis-weights the peptide-level r (`sensitivity.r_prec_weighted`) | either |
| False matches attracted to T, e.g. interference from abundant co-eluting shared proteins (equal-chance violation in the low direction) | low, unquantified |
| A protein-level filter (`PG.Q.Value`) would remove K-12-only proteins first; the pilot applies none, and adding one is a lead finding | low |
| K-12-derived constructs in engineered strains | excluded: only wild-type runs are accepted |

Sensitivities, pre-registered and descriptive (no assessment): `present_broad` widens rule 2 to
rows of `block`, `swap` and `near_block` only; `r_prec_weighted` uses r_prec in equation (1);
`far` restricts E to K-12-only peptides of ≥ 10 aa with no tolerant target row of any type, and T
to shared peptides of ≥ 10 aa, and reports the split of `k12_only` and N_E into `lt10`,
`near_homolog` and `far`.

## Simulation evidence

`tests/calibration/attacks/test_attack_free_entrapment.py` (statistician; seed
`seeds.calibration_base` + 707) simulates a protein-grouped universe with r = 0.25 and checks: the
mean of `fdp_combined` is at least the mean true FDP under equal chance (edge and uniform), and
with whole-protein clustered errors; the excess grows with homolog attraction λ ∈ {1, 2, 4};
attraction to T biases low; `fdp_lower` never exceeds the true FDP; no list with true FDP ≤ α is
labelled `evidence_of_failure`; and `ucb95` covers at ≥ 95 % under equal chance. Two teeth checks
show the tests can fail: `fdp_lower` used as an estimate fails the mean check, and a median bound
fails the coverage check.
