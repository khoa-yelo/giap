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

## Lead decisions (2026-10-22)

1. No stop-codon read-through / Sec planting here; it waits for M7-16's decision on stop-spanning
   peptides (the backlog line stays under M7-16).
2. In-frame indels are excluded from `alt_sites.tsv` (truth = non-synonymous SNVs) and reported in a
   separate truth file, so the F1 gate measures substitutions.
3. Replicate counts for the variants suite come from ADR-0042 (#201, merged), not from a new
   tolerances key.
4. Feature loop: test-author → implementer (giapsim is tooling the implementer writes) → the
   statistician signs off the generators M7-13 gates on → code review.

## Lead decisions on the tests (2026-10-22, after 070a167)

5. The benchmark (`benchmarks/test_sim_m7_bench.py`: each scenario best of 3 under 2 s at 100 kb,
   under 20 s at 1 Mb, recorded in `benchmarks/baseline.json` from a Slurm run) is written by the
   implementer, as giapsim tooling.
6. No site-level error truth file in M7-12: `truth/denovo.tsv` (error type per spectrum) joined
   with `truth/alt_sites.tsv` by peptide is enough; if M7-10's tests need a `(spectrum, site)` file,
   M7-10 adds it.
7. `n_true_peptides` counts database-correct peptides (`truth/peptides.tsv`, `is_wrong = false`),
   as the tests read it; de novo reads are not counted.
8. **RBS bar.** The ≥ 0.90 point-estimate bar stays. Planted ORF RBSs may be drawn from a stronger
   distribution than first-round TU RBSs (documented in the scenario docstring) so the bar is met by
   design, not by false-positive `rbs.tsv` rows; the statistician reviews the generator at sign-off.
9. File names keep M3-20's `peptides_denovo.tsv` and `reference_peptides.tsv` (example 1's wording
   is loose).
10. The optional ADR-0012 σ70 promoter / terminator after planted ORFs is not in M7-12 (backlog).

## Lead decisions on the statistician's block (2026-10-22, PR #225 comment 6356933329)

11. **ALT truth in the sample frame.** After a frameshifting indel upstream in the same gene, the
    ALT codon is taken in the sample's frame and lifted to the reference; rows that become
    synonymous are dropped; `alt_sites.tsv` gains `frameshift_upstream` (bool). Frameshift indels
    are listed in `indel_sites.tsv` (with their junction), so M7-13 can mask them.
12. **False novel candidates.** A new opt-in wrong-peptide source draws peptides from non-coding
    six-frame translations (with matched decoys); `novel_orf_power` turns it on, so the novel-ORF
    FDR gate has false candidates to control.
13. **Genome size for novel-ORF gates.** `novel_orf_power` gates use a ≥ 300-kb genome (a new
    `sim300k` size, ≥ 200 genes) so gene detection is not at the 1/n q floor; 100 kb stays for fast
    unit tests. The GIAP small-genome q floor goes to the backlog.
14. **Planted TIR is a `SimConfig` field.** `novel_orf_power` defaults to the first-round
    distribution lognormal(8, 0.8); the calibration sanity arm keeps lognormal(9, 0.5) (decision 8).
    M7-13 adds the matched arm with an SE allowance.
15. **No crashes.** `features.novel_orfs` fixes the `min_gap` off-by-one; when a genome has too few
    intergenic gaps, the scenario plants as many as fit and the truth records the realised count
    (power denominators read the truth).
16. **Coordinates.** `novel_orfs.tsv` `stop` is written 1-based like the rest of that table.
17. Feature loop: test-author (tests for 11–16) → implementer → statistician re-check → code review.

## Lead decisions on the round-2 tests (2026-10-22, after c2cd9102)

18. The size key is `"300k"` (matching `"100k"` / `"1M"`); `junction = pos_reference + 1` (first
    reference base after the VCF anchor) is accepted.
19. Realised planted counts: the `novel_orfs.tsv` rows plus a WARNING are the record; no extra
    truth file. The per-class expected false-candidate count and the I/L `ms_observable` column go
    to M7-13 (backlog).
20. **Lead edits to existing tests** for decisions 11, 14 and 16: the alt_sites header test accepts
    the appended `frameshift_upstream`; the reference-frame alt_sites checks skip genes with a
    frameshifting indel (the round-2 lift test covers them); the indel_sites set check compares
    the in-frame rows only; `stop` is checked 1-based; the old 100k RBS calibration runs the
    decision 8 sanity arm (`planted_tir = (9, 0.5)`).

## Lead decisions on the statistician's round 2 (2026-10-22, PR #225 comment 6358601933)

21. **B1, scored makedb-style decoys in the simulator.** This applies to the six-frame mode only:
    - Each wrong target row gets a decoy with ID `rev_GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`.
    - Decoy scores come from the same null as the wrong targets, and q/pep from the same mixture.
    - Each null PSM goes to the target or to its entry's decoy with probability ½.
    - The `rev_<locus>` protein decoys stay, so gene-level `decoy_picked` is unchanged.
    - M7-13 keeps the simulator's decoys; `giap makedb` stays out of the calibration loop.
22. **B2, `nterm_extension` territory.**
    - A fourth wrong-peptide territory covers products in frame upstream of each CDS's true start,
      split evenly with the other three.
    - The 1M per-class test covers all four classes.
    - The excluded products of choice 2 are still excluded, except the in-frame upstream extension
      of a CDS, which is this class by definition.
23. **ORF-entry decoys leave gene-level pairing (GIAP change, test-author round 3, 70bf13bc).**
    - In a real makedb search, as in the six-frame scenarios, `rev_GIAPORF_*` decoys reach
      gene-level detection, where they are unpaired; at about 50 % passing q ≤ 0.01 they drop the
      paired share below 80 %.
    - They belong to the M7-11 class FDR. `detection.pair_decoys` therefore drops a decoy whose
      every stripped ID is a makedb ORF entry, as it drops all-contaminant decoys (M4-12
      decision 12); a mixed decoy pairs on its gene IDs.
    - Pinned by tests/unit/genetable/test_detection_orf_decoys.py (lead).
    - Protein decoys are `rev_<annotated locus_tag>` only (no `rev_NOVEL_*`), per test-author
      interface 5.
24. **Test-author interfaces 1–6 (70bf13bc) accepted:**
    - a decoy's class is the makedb entry's class;
    - N(mu0, sd0) scores with mixture q/pep at `pi_wrong_realised`;
    - a per-row ½ coin;
    - `origin = sixframe_noncoding` is kept;
    - the fourth territory is drawn in turn in `_sixframe_draw`.
25. **Round-2 decoy-count test retired (implementer round 3, 9d91aa93).**
    `test_sixframe_no_shortfall_and_matched_decoys` asserted that the decoy rows equal the wrong
    rows. Under decision 21 they cannot: six-frame decoys come from a ½ coin, and the
    `rev_<locus>` protein decoys are added on top (seed 7, 100k: T = 63, D6 = 54, Dp = 117).
    - The test becomes `test_sixframe_no_shortfall_and_null_scores`: no shortfall, both kinds of
      row present, and the wrong-ID score gap.
    - The pooled coin is checked by the round-3 calibration file.
    - Implementer deviations 1–4 are accepted:
      - protein decoys count `n_decoy_rows`;
      - giapsim restates the makedb rules, held honest by the round-3 oracle;
      - the speed is acceptable and the baseline is not re-recorded;
      - the `alt_sites` initiator-codon fix.
