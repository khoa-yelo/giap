# M7-13: The M7 gate — variant and novel-ORF F1, the class-FDR suite, attack 6

Milestone M7 · track C (calibration) · **pass 2 — the gate** · depends on: M7-05, M7-09, M7-10,
M7-11, M7-12 (and the wired steps M7-04, M7-08) · estimate ≈ 90 code lines (+ ≈ 600 test lines) ·
**stats**

## Spec sections implemented

- §13 row M7 gate: "F1 ≥ 0.9 on simulated variants and novel ORFs; TE2 proxy passes" (TE2 is
  M7-14).
- §7.5.6 "**Gate M7:** simulated F1 ≥ 0.9 for `high` ALT calls; no `high` ALT calls on in-silico
  same-strain peptides except at `start_codon` or `stop_codon` sites."
- §7.6.6 "**Gate M7:** F1 ≥ 0.9 for planted novel ORFs with ≥ 2 peptides"; "**Calibration:** class
  FDR on simulations with planted novel ORFs and de novo errors".
- §8.4 (FDR criterion "Mean FDP over R replicates ≤ q + 2 SE(FDP)", fixed seeds) and suite row
  "Novel-ORF class FDR | … | FDR | 50 | 300"; `tolerances.yaml` `f1_min: {variants_high: 0.90,
  novel_orfs: 0.90}`, `replicates_fast/full.novel_orfs` (50 / 300); ADR-0015 (two-block
  confirmation).
- §8.5 attack 6: "Isobaric swaps everywhere → `mass_equivalent` flags present; no `high` ALT calls
  from swaps alone at q ≤ 0.01."

## Checks

| Check | Scenario (M7-12) | Prediction set | Truth set | Criterion | R fast / full |
|---|---|---|---|---|---|
| **G1 variant F1** | `variant_f1` | `variants.tsv` rows with `confidence = high`, `alt_type = substitution` | `alt_sites.tsv` rows with `n_true_peptides ≥ 2` | pooled F1 (counts summed over replicates) ≥ `f1_min.variants_high` | 20 / 100 (Q1) |
| **G2 same-strain** | `same_strain_insilico` | `high` ALT rows | none | zero `high` rows without `start_codon` or `stop_codon`, in every replicate | as G1 |
| **G3 novel F1** | `novel_orf_power` | `called = true` ORFs | planted ORFs (intergenic, cds_outframe, antisense, nterm_extension) with `n_true_peptides ≥ 2` | pooled F1 ≥ `f1_min.novel_orfs`; match = same `(contig, strand, stop)` and class | 50 / 300 |
| **G4 class FDR** | `novel_orf_power` + false peptides + de novo errors, both regimes | ORFs with class q ≤ q | planted ORFs | mean FDP ≤ q + 2 SE at q ∈ {0.01, 0.05, 0.10}, per regime | 50 / 300 |
| **A6 isobaric attack** | `isobaric_attack` | `high` ALT rows | none | (a) every ALT row from an isobaric error carries `mass_equivalent` (swaps: see M7-05 Q1); (b) share of replicates with ≥ 1 `high` ALT ≤ 0.05 + 2 SE | 50 / 300 |
| recorded | `variant_f1` | — | — | class-error bias with planted variants (M7-10 rule 5); start accuracy of called ORFs (true start ∈ {longest, best}); share of planted RBS starts passing the threshold | recorded only |

## Rules

1. All checks go through `_criteria.calibrate()` (ADR-0015 confirmation; seeds from
   `SeedSequence([calibration_base, crc32(suite), rep])`).
2. F1 is pooled: TP, FP, FN summed over replicates, then F1 = 2TP / (2TP + FP + FN). A replicate
   with no truth sites still contributes its FPs.
3. Fast tier runs at 100 kb within the §8.4 10-minute budget; full tier at 1 Mb through
   `scripts/slurm/calibration_full.sbatch` tasks `variants` and `novel_orfs` (≤ 2 array tasks,
   ≤ 32 CPUs; sum the running `giap-*` CPUs first, R14).
4. **A failure is fixed by an ADR on the method, never by loosening `tolerances.yaml` (R3).** The
   likely candidates are written down now (risks).

## Risks (decided before the numbers)

- **A6 is expected to be hard.** Under `related`/de novo, a systematic misread (N read as GG by every
  spectrum of two overlapping peptides) is two distinct peptides with the same ALT → `high`. §7.5.3
  says flags never change the call, so (b) can fail by design. ADR options, for the statistician:
  (i) `high` requires ≥ 2 distinct peptides whose ALT is *not* `mass_equivalent`; (ii) `high` also
  requires `class_error_rate` ≤ a bound when it is estimated; (iii) restate the attack as
  "`mass_equivalent` present on every such call". Khoa decides an ADR that changes a call rule.
- **G3 depends on the RBS threshold**: if planted intergenic ORFs fail RBS support, the gate
  measures the simulator (M7-12 adds a sanity check for that).
- **G4 de novo regime**: the ORF-level minimum of peptide q (M7-11 decision 3) may be
  anti-conservative for ORFs with many peptides.

## Acceptance criteria

- [ ] `tests/calibration/test_variants_gate.py` (G1, G2, A6) and `test_novel_orfs_gate.py`
  (G3, G4), fast and full markers.
- [ ] Full-tier results (Slurm job IDs) and the recorded metrics in PROGRESS.md and in
  `docs/methods/variants.md` / `novel-orfs.md` (§8.4 calibration log).
- [ ] The M7 gate in PROGRESS.md lists G1–G4, A6 and TE2 (M7-14) separately; any `unconfirmed`
  check is visible.

## Questions for the lead

1. `tolerances.yaml` has `f1_min` but no `replicates_*.variants`. **Proposed default (needs Khoa:
   protected file):** add `variants: 20` (fast) / `100` (full).
2. G1 restricts the prediction set to substitutions, because the truth (non-synonymous SNVs,
   M7-12 Q2) has no other `alt_type`. **Proposed default:** yes; blocks and swaps planted by de novo
   errors are false positives *if* they are `high` — they are counted in FP when they hit a truth
   codon, and their rate is reported.
3. Match rule for G3: `(contig, strand, stop)` + class, start ignored (start accuracy recorded).
   **Proposed default:** yes — the stop identifies the ORF (M7-06 rule 3); the start is
   OSTIR-dependent and §7.6.6 asks for ORF recovery.
4. A6 criterion (b) as an FWER-style "share of replicates with ≥ 1 `high` call". **Proposed
   default:** yes; zero-tolerance per replicate is the stated attack, and the 0.05 + 2 SE form is
   the §8.4 FWER criterion.
