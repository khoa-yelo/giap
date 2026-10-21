# M7-12: giapsim scenarios for variants and novel ORFs

Milestone M7 · track F (fixtures and simulator) · **pass 2** · depends on: merged giapsim (M3-18
features and strains, M3-19 proteomics, M3-20 de novo errors, M2-05 layers); M7-08's `rbs_novel`
layer (pass 1) · estimate ≈ 40 code lines (+ ≈ 320 simulator lines, + ≈ 350 test lines) · stats
(light): the scenarios define what the gate measures, so the statistician reviews them

## Spec sections implemented

- §11.3 "Planted features: … unannotated expressed ORFs (intergenic, out-of-frame, antisense);
  N-terminal extensions (annotated start set downstream of the true start)"; "Strain differences.
  A 'sample' genome differs from the 'reference' genome by SNVs (ALT residues) and indels"; "De novo
  errors: substitutions, mass-equivalent swaps and blocks (N ↔ GG, Q ↔ AG, …), I/L, adjacent swaps,
  reversed stretches, modifications, and random peptides … A `spectrum_id` pairs each de novo row
  with a reference identification"; "Signals with tunable strength: RBS motifs (AGGAGG at spacer
  5–10); σ70 −35/−10 motifs upstream of TU starts".
- §7.5.6 calibration: "classification F1 on simulated de novo errors (substitutions, swaps, blocks,
  modifications), with variants planted between 'sample' and 'reference' genomes"; §7.6.6
  calibration: "class FDR on simulations with planted novel ORFs and de novo errors"; §8.4 suite
  "Novel-ORF class FDR | planted ORFs + false peptides + de novo errors"; §8.5 attack 6.
- §11.6 TE2's in-silico same-strain set, simulated (the real-genome proxy is M7-14).

## What exists and what is missing (checked on c4d8841)

| Need | giapsim today | Added here |
|---|---|---|
| planted strain differences | `StrainConfig` SNVs + indels, `truth/strain_differences.tsv` (VCF orientation, `aa_change`) | `truth/alt_sites.tsv`: one row per non-synonymous SNV codon on the **reference** genome (the genome GIAP is given): contig, codon_start, strand, gene, ref_aa, alt_aa, `n_true_peptides` covering it |
| de novo errors with pairing | `simulate_denovo` (7 error types, `spectrum_id`, `reference_peptides.tsv`) | an `error_weights` preset per scenario; error truth per (spectrum, site) for the class-error check |
| novel ORFs | `n_novel_orfs` (intergenic planted, 10–30 aa; out-of-frame/antisense found inside CDSs), `n_nterm_extensions` | configurable length range (default stays 10–30; scenarios use 30–80 aa so ≥ 2 tryptic peptides are typical); an RBS motif at spacer 5–10 before planted intergenic starts and N-terminal true starts; their TIR rows in `rbs_novel` drawn from the "true" distribution (M7-08 rule 6 already plants the background); optional σ70 motif and a terminator after planted intergenic ORFs (ADR-0012 columns exercised) |
| truth for novel calls | `truth/novel_orfs.tsv` (id, class, coordinates, aa) | `stop` (0-based stop codon start), `true_start`, `n_true_peptides`, `rbs_planted` |
| same-strain peptides | M3-19 database peptides with `search_db = same` | scenario `same_strain_insilico`: every tryptic peptide (7–30 aa, ≤ 1 missed cleavage) of every annotated protein, Bakta-like M at non-ATG starts, no errors |
| isobaric attack | — | scenario `isobaric_attack`: de novo errors restricted to block, near-block, modification and swap types, **no** strain differences |

## Public interface

`tests/sim/giapsim/scenarios.py` gains
```python
def variant_f1(seed: int, size: Literal["100k", "1M"] = "100k") -> SimConfig: ...
def novel_orf_power(seed: int, size: Literal["100k", "1M"] = "100k") -> SimConfig: ...
def same_strain_insilico(seed: int) -> SimConfig: ...
def isobaric_attack(seed: int) -> SimConfig: ...
```
and writes the truth tables above. `src/` gets only the ≈ 40 lines of loaders the calibration
tests share (`tests/calibration/_m7.py` is test code; the 40 lines are a truth-table reader in
`giapsim`, counted as tooling).

## Rules

1. **Coordinates.** Truth is written on the genome GIAP receives. In strain scenarios that is
   `reference.fasta` (the related strain); peptides come from the sample. `truth/alt_sites.tsv`
   lifts each SNV to the reference codon and strand, 1-based closed as every truth table.
2. **Truth denominators** are written, not computed by the tests: `n_true_peptides` counts distinct
   (I/L-folded) simulated *correct* peptides covering the site (variants) or the ORF (novel), so
   M7-13 can restrict the gate to "planted with ≥ 2 peptides" exactly as §7.6.6 says.
3. **Determinism.** New draws use new `SeedSequence` keys (≥ 5), so every existing output stays
   byte-identical (the M3-20 pattern; a test pins it).
4. **Speed.** Each scenario at 100 kb in < 2 s, at 1 Mb in < 20 s (§11.3 budget).

## Examples

| # | Scenario | Expected |
|---|---|---|
| 1 | `variant_f1(seed)` | `reference.fasta`, sample-derived `peptides.tsv` (de novo) and `reference_peptides.tsv` (database, paired `spectrum_id`); `truth/alt_sites.tsv` with ≥ 20 sites covered by ≥ 2 true peptides at 1 Mb |
| 2 | `novel_orf_power(seed)` | ≥ 10 planted ORFs per class among intergenic, cds_outframe, antisense, nterm_extension at 1 Mb, with `n_true_peptides ≥ 2` for most; RBS motifs before the intergenic ones |
| 3 | `same_strain_insilico(seed)` | peptides with `search_db = same`, no errors; non-ATG starts give N-terminal peptides beginning with M |
| 4 | `isobaric_attack(seed)` | `alt_sites.tsv` empty; de novo errors only of the four isobaric kinds |
| 5 | any existing config (no new field set) | byte-identical outputs to c4d8841 |

## Acceptance criteria

- [ ] `tests/sim/test_m7_scenarios.py` (test-author): examples 1–5; truth tables parse; the counts
  of rule 2 match a brute-force recount on the 100-kb outputs.
- [ ] Calibration sanity in `tests/calibration/test_sim_m7_calibration.py`: planted RBS TIRs sit
  above the 10 % detected-gene quantile in ≥ 90 % of planted starts (so the gate tests GIAP, not
  the simulator).
- [ ] Runtime within rule 4 (benchmarks entry).

## Questions for the lead

1. Stop-codon readthrough / Sec (backlog line "the simulator plants stop-codon readthrough … for
   `cterm_extension` calibration", re-filed to M7 by the M6 lead): plant it here? **Proposed
   default:** not until M7-16 decides whether the matcher can place stop-spanning peptides; with the
   merged matcher such peptides are unmapped and the scenario would test nothing.
2. Indels in `strain_differences.tsv`: in-frame indels change peptide length and map as tolerant
   blocks or not at all. **Proposed default:** excluded from `alt_sites.tsv` (truth = non-synonymous
   SNVs only) and reported separately, so the F1 gate measures substitutions, the bulk of real ALT
   calls.
