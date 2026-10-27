# M7 "Variants and novel ORFs": issue index

Spec: GIPPA_SPEC.md §13 row M7. Deliverables: "Modules 5 and 6, VCF, patch GFF3, second-round
Promotech and OSTIR". **Gate:** "F1 ≥ 0.9 on simulated variants and novel ORFs; TE2 proxy passes",
i.e. §7.5.6 ("simulated F1 ≥ 0.9 for `high` ALT calls; no `high` ALT calls on in-silico same-strain
peptides except at `start_codon` or `stop_codon` sites"), §7.6.6 ("F1 ≥ 0.9 for planted novel ORFs
with ≥ 2 peptides") and §11.6 TE2. Nothing in M7 is blocked on an external source.

Sources: §7.5 (whole module: §7.5.1 the source table, §7.5.2 pileup and calls, §7.5.3 flags,
§7.5.4 error control, §7.5.5 outputs and the VCF, §7.5.6 tests and gate), §7.6 (whole module:
§7.6.1 candidates, §7.6.2 reconstruction and the second round, §7.6.3 classes, §7.6.4 error
control, §7.6.5 ranking, §7.6.6 outputs and tests), §7.1.4 and §7.1.6 (targeted Promotech, OSTIR
input), §7.1.1 (tool pins and gotchas), §7.2.1–§7.2.6 (the index, the tolerant operations, placement
FDR, locus classes), §7.9.2 (`makedb`, class-tagged IDs, `--patch`), §2.1 decisions 13, 15, 16, §2.4
G2, G8, G13, G14, §5.2 (`variants.*`, `novel_orfs.*`, `sample_name`, `annotation.precomputed`), §6.2
rows `variants`, `novel_candidates`, `promotech_novel`/`ostir_novel`, `novel_orfs`, §6.3–§6.7,
§8.1 rows "ALT call is real" and "Novel ORF", §8.2, §8.4 suite "Novel-ORF class FDR", §8.5 attacks 5
and 6, §9.1–§9.5 (`pileup.tsv.gz`, `variants.tsv`, `variants.vcf`, `novel_orfs.tsv`, Variants and
Novel-ORF flag scopes), §11.2–§11.7 (golden "variant tables", the simulator's planted ORFs, strain
differences and de novo errors, TE2, TE7), Appendix B.5 (masses) and B.7 (the VCF record), the
`tolerances.yaml` template (`f1_min`, `replicates_*.novel_orfs`).

ADRs read: ADR-0005 (VCF writer strictness; "a codon crossing the origin raises in the VCF writer
(M7 decides how to report it)"; the generic FILTER line M7 should replace), ADR-0006 (value
encoding), ADR-0007 (size counts code), ADR-0008 (Promotech @ `56251ad` GPL-3.0, OSTIR 1.1.3
GPL-3.0, conda envs, no upstream images), ADR-0009 (picked FDR rescaled by the paired share),
ADR-0010 (π̂-weighted decoy FDR), ADR-0011 (vocabulary additions went to Khoa), **ADR-0012**
(promoter and terminator evidence for novel ORFs: four columns, two flags, a ranking key; M7 must
implement it), **ADR-0013** (`cterm_extension`: candidate class, `readthrough_candidate`, patch
rows), ADR-0015 (two-block confirmation), ADR-0019 (stratum-0 π̂), ADR-0021 Amendment 1 (no
`uncalibrated` flag anywhere).

One issue = one PR; size counts **code** only (ADR-0007; `tests/sim/` lines are tooling and listed
separately). "Stats" = statistical change: label `stats`, the statistician reviews the methods note
and the simulations before the code (§14.1).

## How M7 is sequenced (Khoa, 2026-10-20/21)

"Write the spec as usual, but prioritize end to end integration first, then work on the rigour
after." This index carries the **whole** of M7 — every table, file, test, calibration and attack
the spec owes — split into two passes:

- **Pass 1 (integration):** two tracks, each wired into `gippa run` **as soon as its minimal
  producers exist**, and completed inside the already-wired step afterwards.
  - Variants: 01 alignment → 02 pileup and calls → 03 tables and VCF → **04 wiring + real run** →
    05 flags (fills columns of the wired step).
  - Novel ORFs: 06 candidates and ORFs (with the minimal ≥ 2-peptide call) → 07 the four files →
    **08 steps, second-round Promotech and OSTIR, wiring + real run** → 09 class evidence, RBS
    support, ADR-0012 columns and ranking (fills columns of the wired step).
  - End state: `gippa run` on sim100k (both executors) writes `tables/pileup.tsv.gz`,
    `tables/variants.tsv`, `tables/variants.vcf`, `tables/novel_orfs.tsv`, `tables/novel_orfs.gff3`,
    `tables/novel_orfs.faa`, `tables/annotation_patch.gff3`, with ALT error rates and novel-ORF
    class FDR in their spec-sanctioned `unestimated` states.
- **Pass 2 (hardening and the gate):** 10 class error rates, 11 class FDR, 12 simulator scenarios,
  **13 the gate** (F1 × 2, class-FDR suite, attack 6), **14 TE2 proxy** and `id_mixups`, 15
  integration hardening (goldens, in-pipeline contracts, resume/failure, real-tool arm).
- **Decided (Khoa, 2026-10-24: option B, ADR-0058):** 16 stop-spanning peptides (read-through,
  `cterm_extension`; Sec at internal codons only if Khoa accepts ADR-0058 item 6). Three PRs:
  16a matcher + decoys, 16b M7 calls + `makedb` `X` + simulator, 16c recoded Sec/Pyl (gated).

Nothing is dropped; it is scheduled. Each producer issue states how its output reaches
`gippa run`: 04 for the three Module 5 files, 08 for the four Module 6 files.

## How this design lands (docs only, one branch `m7/design`, three commits)

| Commit | Files |
|---|---|
| 1/3 | this index, M7-01 … M7-05 (the variants track of pass 1) |
| 2/3 | M7-06 … M7-09 (the novel-ORF track of pass 1) |
| 3/3 | M7-10 … M7-16 (pass 2 and the gated issue) |

At ≈ 1,900 lines the design exceeds the ≈ 800-line review budget; the lead may open one PR per
commit (the M6 pattern) after deciding the questions below.

## Issues

| Order | Issue | Track | Pass | Scope | Depends on | Code lines | Stats |
|---|---|---|---|---|---|---|---|
| 1 | [M7-01](M7-01-residue-alignment.md) | V | 1 | residue-level traceback of every accepted hit (exact, tolerant, frameshift); the M3 mapper keeps none | merged M3 | 240 | no |
| 2 | [M7-02](M7-02-pileup-and-calls.md) | V | 1 | distinct-peptide pileup, REF/ALT/mixed/no_call, high/candidate, §7.5.1 exclusions | 01 | 260 | no |
| 3 | [M7-03](M7-03-variant-tables-and-vcf.md) | V | 1 | `pileup.tsv.gz`, `variants.tsv`, VCF 4.3 (G14 minimal codon, minus strand, multi-allelic), schemas | 02 | 330 | no |
| 4 | [M7-04](M7-04-variants-step-and-run.md) | T | 1 | **`variants` step, VARIANTS subworkflow, `sample_name` resolution, real Slurm run** | 01–03 | 300 | no |
| 5 | [M7-05](M7-05-variant-flags.md) | V | 1 | the seven §7.5.3 flags, `build_mass_equivalents.py` + TSV, exhaustive snv table | 02–04 | 260 | no |
| 6 | [M7-06](M7-06-novel-candidates-and-orfs.md) | N | 1 | candidates, exclusions, ORF reconstruction (both strands, wrap), grouping, IDs, second-round windows and OSTIR rows, minimal call | merged M3 | 360 | no |
| 7 | [M7-07](M7-07-novel-orf-outputs.md) | N | 1 | `novel_orfs.tsv` (+ ADR-0012 columns), `.gff3`, `.faa`, `annotation_patch.gff3`, makedb round trip | 06 | 300 | no |
| 8 | [M7-08](M7-08-novel-steps-and-second-round.md) | T | 1 | **five steps, PROMOTECH_NOVEL / OSTIR_NOVEL, `rbs_novel` layer + simulator, NOVEL subworkflow, real Slurm run** | 06, 07 | 380 (+90 sim) | no |
| 9 | [M7-09](M7-09-novel-evidence-and-ranking.md) | N | 1 | §7.6.3 class rules, RBS quantile, ADR-0012 promoter/terminator, ADR-0013 label, G13 rank | 06–08 | 300 | light (Q) |
| 10 | [M7-10](M7-10-class-error-rates.md) | S | 2 | §7.5.4 paired-spectrum class error rates | 02, 04, 05 | 180 | yes |
| 11 | [M7-11](M7-11-novel-class-fdr.md) | S | 2 | §7.6.4 picked class FDR (makedb IDs) and de novo subgroup FDR | 06, 09 | 300 | yes |
| 12 | [M7-12](M7-12-sim-variant-novel-scenarios.md) | F | 2 | gippasim: `variant_f1`, `novel_orf_power`, `same_strain_insilico`, `isobaric_attack`, truth tables | 08 (layer) | 40 (+320 sim) | light |
| 13 | [M7-13](M7-13-calibration-gate.md) | C | 2 | **the gate**: variant F1, same-strain check, novel F1, class-FDR suite, attack 6 (3 PRs: 13a variants, 13b novel, 13c evidence) | 05, 09–12, 11b | 60 harness (+≈ 750 tests) | yes |
| 14 | [M7-14](M7-14-te2-proxy-and-id-mixups.md) | V | 2 | **TE2 proxy** (Nissle, weekly Slurm) and the `id_mixups` QC row | 04, 05 | 170 | no |
| 15 | [M7-15](M7-15-integration-hardening.md) | T | 2 | goldens, in-pipeline contracts, resume/failure, real-tool arm, VCF FILTER, `novel_candidate_id` | 04, 08, 09 | 200 | no |
| 16 | [M7-16](M7-16-stop-codon-readthrough.md) | R | 2 | option B (ADR-0058): annotated terminal stop consumable at cost 1, mirrored in decoys; 16a matcher/decoys/index, 16b variants U/O, cterm call + patch, `makedb` `X`, simulator; 16c recoded Sec/Pyl (needs Khoa) | merged M3; 01, 05, 07, 09, 11b, 12, 18 | 370 (+130 sim) | yes |
| 18 | [M7-18](M7-18-variants-ambiguous-residues.md) | V | 2 | variants on B/J/O/U/X/Z peptide residues: skip ambiguous/rare, U@TGA and O@TAG as REF, one event (K-12 crash fix) | 02–05 | 75 | no |

Total ≈ 4,080 code lines (M7-16 re-estimated 220 → 370 after option B) (+ ≈ 540 simulator lines, ≈ 5,500 test lines); pass 1 ≈ 2,730. The
spec's 2 weeks fits if the two pass-1 tracks run in parallel (they share no code until 08/09) and
M7-12 is written while pass 1 finishes.

## Dependency DAG

```text
Pass 1  V: 01 ─► 02 ─► 03 ─► 04 (gippa run: 3 Module 5 files) ─► 05
        N: 06 ─► 07 ─► 08 (gippa run: 4 Module 6 files, second round) ─► 09
Pass 2  02,04,05 ─► 10 ┐
        06,09 ──► 11 ──┤
        08 ─────► 12 ──┴─► 13 (the gate)        04,05 ─► 14 (TE2)      04,08,09 ─► 15
B (ADR-0058)  16a ─► 16b ─► (re-run 13 G2 and 14 with stop_codon sites; M7-07 cterm patch rows)
              16a ─► 16c (only after Khoa accepts ADR-0058 item 6)
```

01 and 06 have no M7 dependency and start first; 06 is the longer chain (its wiring issue 08 is
the largest), so N is the critical path.

## Waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01, 06 | independent tracks |
| B | 02, 07 | |
| C | 03, 08 | 08 includes the OSTIR range recording (one small `sbatch`) |
| D | 04, 09, 05 | 04 is the Module 5 integration milestone (real `gippa_run.sbatch`); 08's real run ends wave C |
| E (hardening) | 10, 11, 12 in parallel; then 13, 14, 15 | 13 submits `calibration_full.sbatch` tasks `variants`, `novel_orfs` (≤ 2 jobs, ≤ 32 CPUs); 14 adds a weekly entry |
| F | 16a, then 16b (and 16c if Khoa accepts item 6) | option B accepted; statistician reviews M7-16 before 16a code; 16a re-runs the M3-21 placement suites fast + full on Slurm |

## What M7 consumes (checked against merged main, c4d8841)

| ID | Artefact | Producer (status) | What M7 relies on |
|---|---|---|---|
| C11 | `mapping/peptide_hits.tsv` | M3-14/M3-16 (merged) | accepted loci at best cost: `peptide` (I/L as observed), interval, strand, frame, `route`, `n_diff`, `diff_types`, `placement_q`, `locus_class` (incl. `cterm_extension`), `uniqueness`, `gene`; `novel_candidate_id` is NA ("novel-locus clustering is Module 6") |
| W1 | `work/hits_exact.parquet`, `work/hits_tolerant.parquet` | M3-04/06/08 (merged) | frameshift `blocks`; decoy rows with mirrored target intervals (M7-11). **No traceback** is stored anywhere (M7-01) |
| C9 | `peptides.norm.parquet` | M1-12 (merged) | `source`, `search_db`, `q_value`, `is_decoy`, `is_contaminant`, `protein_ids`; **no `spectrum_id`** (M7-10 Q1) |
| C1/C2/C15 | `replicons`, `genes`, `features` | M2-04 (merged) | topology, gene intervals, `start_codon`, `feature_type` (pseudogene, sorf); no pseudogene junction position (M7-09 Q2) |
| C4 | `annotation_meta.json` | M2-04 (merged) | `anti_sd` for the second-round OSTIR rows |
| G | `genetable/genes.parquet` | M4-15 + M4-13 (merged; wired by MVP-01) | `detected`, `partial`, `tir_ostir` for the RBS quantile |
| L | `promoter_windows.parquet`, `terminators.parquet` | M2-13, M2-15 (merged) | first-round scores by position; terminators for ADR-0012. `annotate.nf` must **emit** them — M6-08 does that (lead decision 3 of M6); if M6-08 has not merged, M7-08 adds the emission |
| M | `gippa.mapping.masses` | M3-05 (merged) | `mass_equivalents(tol)`, B.5 masses (no `build_mass_equivalents.py` or `mass_equivalents.tsv` yet: M7-05) |
| V | `gippa.io.vcf.write_vcf`, `gippa.io.gff3`, `gippa.io.protein_fasta` | M1-03/M1-05 (merged) | strict VCF 4.3 writer (equal-length alleles, origin-crossing REF raises), GFF3 with `Is_circular` |
| D | `gippa makedb --patch`, class-tagged six-frame ORFs, `db_id_map.tsv` | M3-23 (merged) | the round trip; **a patch CDS with an internal stop raises** (M7-07 Q4, M7-16) |
| T | `validation/proxies/te_proxies.py` | M3-22 (merged) | `run_chain`, the Nissle digest, the TE1 exception lists (431 non-ATG N-termini; 4 Sec misses) |
| S | `gippa.stats.tda.picked_fdr`, `monotone_q`; `gippa.mapping.placement_fdr.stratum_qvalues` | merged | M7-11 |
| H | calibration harness | M4-05 (merged) | `calibrate`, `check_fdr`; `tolerances.yaml` has `f1_min` and `replicates_*.novel_orfs`, **no `variants` replicate key** (Khoa Q1) |
| K | config | merged `VariantsConfig`, `NovelOrfsConfig`, `stats.alpha.{variants,novel_orfs}`, `annotation.start_codons`, `annotation.ostir.*`, `annotation.promoters.*`, `sample_name` | `sample_name: auto` is resolved nowhere (M7-04); `annotation.precomputed` has no second-round OSTIR key (M7-08 Q1) |
| Sim | gippasim | M3-18/19/20, M2-05 (merged) | strain differences with truth, planted novel ORFs and N-terminal extensions with truth, de novo errors with `spectrum_id` pairing; no RBS before planted ORFs, no read-through (M7-12) |
| W | `workflow/subworkflows/local/*.nf`, `PROMOTECH`, `OSTIR` modules | MVP-01, M2-13, M2-16 (merged) | the CONTEXT pattern; `PROMOTECH`/`OSTIR` reused by alias; `mapping.nf` must emit `hits_tolerant` (M7-04) |
| R | `gippa.runner.local.PUBLISHED` | merged | publishes only `annotation/` and `mapping/` on c4d8841; M7-04/08 add the `tables/` keys unless the table-publishing change has landed first |

## External tools (no new tool, env or image)

| Tool | Pin | Licence | Env / container | Used by | Slurm? |
|---|---|---|---|---|---|
| Promotech RF-HOT | commit `56251ad9b883ef831b4753fc623d5ec970fe65e0` (ADR-0008) | GPL-3.0 — separate process, never imported (`check_imports.py`) | `workflow/envs/promotech.yml` (Python 3.6.10, TF 2.1.0 wheel; lock in `workflow/envs/locks/`); RF-HOT model from `$GIPPA_DB_DIR/promotech` (installed, PROGRESS.md); no upstream image, conda under `-profile apptainer` | M7-08 `PROMOTECH_NOVEL` | under Nextflow `farmshare` profile each chunk is a Slurm task; sim100k second round ≤ 2 chunks |
| OSTIR | 1.1.3 + ViennaRNA 2.6.4 (ADR-0008) | GPL-3.0 — separate process | `workflow/envs/ostir.yml` (lock present); no upstream image | M7-08 `OSTIR_NOVEL`; the range recording | one task per run; the recording is a minutes-long `sbatch` |
| VCF validation | `pysam` (already a test dependency, MIT) | MIT | main env | M7-03 round trip | no |
| GFF3 validation | GIPPA's M1-05 reader | MIT | main env | M7-07 | no |

No `bcftools`, `vcf-validator` or other new binary is introduced. The TE2 proxy (M7-14) and the
full calibration tier (M7-13) are Slurm jobs; the per-PR fast tier fits the 8-CPU session and the
10-minute §8.4 budget.

## What happens when an estimator does not exist (summary; details per issue)

| Estimate | Does not exist when | Output |
|---|---|---|
| residue alignment | recomputed cost ≠ `n_diff` | the step fails (a mapper/aligner disagreement is a bug, not data) |
| ALT confidence | — (always defined) | `high` / `candidate` |
| class error rate | no `--reference-peptides`, no `spectrum_id`, or no paired observation in the class | `class_error_rate = unestimated`, `class_error_n` NA (§7.5.4) |
| VCF record | origin-crossing codon; 1↔2 block; frameshift | `variants.tsv` row only; `vcf_unrepresentable` for the first |
| `start_best_ostir` | no second-round OSTIR row at any candidate | NA; `*_best` columns fall back to the longest start, flag `ostir_missing` |
| RBS support | < 20 detected genes with an OSTIR score, or no OSTIR at the chosen start | `rbs_supported` NA; `nterm_extension`/`intergenic` not called; flag `rbs_threshold_unestimated` |
| class FDR | no class-tagged decoys and not de novo; or a class's paired share < 80 % | `class_fdr = unestimated`, flag `class_fdr_unestimated`, ≥ 2-peptide rule (§7.6.4) |
| terminator evidence | no terminator layer | NA, run flag `layer_missing:terminators` |
| `id_mixups` | no peptides | the M4-16 NA + `requires_variants` stays |
| `stop_codon` / `cterm_extension` calls | until M7-16a/16b merge; afterwards only at annotated terminal stops (internal Sec codons need ADR-0058 item 6) | implemented, unit-tested, unreachable until M7-16 |

No silent fallbacks: every NA carries a flag, and no `unestimated` value is ever turned into a
number. No `uncalibrated` flag anywhere (ADR-0021 Amendment 1).

## Calibration suites and attacks owed by M7

| Suite / check | Issue | Generator (M7-12) | Criterion | R fast / full |
|---|---|---|---|---|
| variant F1 (`high`) | M7-13 G1 | `variant_f1` (100k / 1M) | pooled F1 ≥ 0.90 (`f1_min.variants_high`) | 20 / 100 (ADR-0042) |
| same-strain `high` ALT | M7-13 G2 | `same_strain_insilico` | zero non-start/stop `high` calls per replicate | as G1 |
| novel-ORF F1 | M7-13 G3 | `novel_orf_power` 1M, `NOVEL_PER_CLASS_GATE` | pooled F1 ≥ 0.90 (`f1_min.novel_orfs`) | 10 / 300 (ADR-0048, pending Khoa) |
| `novel_orfs` class FDR | M7-13 G4 | `novel_orf_power` 1M, `NOVEL_PER_CLASS_GATE` + six-frame false peptides + de novo errors | FDR ≤ q + 2 SE, q ∈ {0.01, 0.05, 0.10}, per regime and class; `vacuous` below the class q floor | 10 / 300 (ADR-0048, pending Khoa) |
| class error rate bias | M7-10 | `variant_f1` without variants | \|bias\| ≤ 0.05 | as G1 |
| attack 6 (§8.5) | M7-13 A6a–c | `isobaric_attack` (+ systematic transform) | `mass_equivalent` on every isobaric non-swap call; ≤ 0.05 + 2 SE replicates with a `high` swap | 20 / 100 |
| attack (new) reversed stretches in intergenic space | M7-11/13 | de novo errors | class FDR at level (shuffled decoys) | as G4 |
| attack (new) GC-poor intergenic false peptides | M7-11/13 | false peptides | class FDR at level | as G4 |
| TE2 (real genome) | M7-14 | Nissle in-silico same-strain peptides | zero violations | weekly Slurm |

**Attack 6 is expected to be the hard one** (M7-13 risks): systematic isobaric misreads produce
two distinct peptides with the same ALT, and §7.5.3 forbids flags from changing calls. If it fails,
the fix is an ADR on the call rule (Khoa), never a tolerance change.

## Methods notes (statistician sign-off, §14.1, §14.9)

New: `docs/methods/variants.md` (the alignment recomputation and its equality check, distinct
peptides and I/L, calls and confidence, flags and why they never change a call, the class error
rate and its upward bias with real variants — M7-10 rule 5, the VCF representation choices),
`docs/methods/novel-orfs.md` (reconstruction, grouping by stop, classes and evidence, the RBS
quantile and its conservativeness, the three FDR regimes and the ORF ↔ makedb matching, ranking,
ADR-0012 evidence and what a missing promoter or terminator does *not* prove). Updates:
`target-decoy.md` (class-level picked FDR), `placement-fdr.md` (subgroup FDR by locus class),
`calibration.md` (the M7 suites and attacks).

## ADRs M7 may owe

- **Second-round OSTIR under the local executor** (M7-08 Q1): a new `annotation.precomputed.rbs_novel`
  key is a config-interface change → ADR (lead; Khoa if the lead judges it a new key in Khoa's
  sense, cf. ADR-0005 item 4).
- **VCF representation** (M7-03 Q1–Q6: tie-break, which `alt_type`s, origin codons, `J`, GT, FILTER
  PASS only): gap-filling, no default changes → lead ADR, like ADR-0005/0006.
- **New flag values** (`vcf_unrepresentable`; Novel-ORF `no_start_codon`, `partial`,
  `ostir_missing`, `rbs_threshold_unestimated`, `junction_unknown`, plus ADR-0012's
  `promoter_supported`/`terminator_supported` and ADR-0013's `readthrough_candidate` already
  accepted): through the schema; ADR-0011/0013 took vocabulary additions to Khoa, M6 lead decision 6
  accepted them as lead — the lead decides which precedent applies.
- **Stop-spanning peptides** (M7-16): ADR-0058, accepted by Khoa (B, `X`); item 6 pending.
- An ADR on the call rule **only if** attack 6 or a gate check fails (R3).

## Backlog items folded in (docs/backlog.md, not duplicated)

| Backlog line | Where it lands |
|---|---|
| "M7 variants and novel ORFs (declare GIPPA's VCF FILTER values with descriptions; M1-05 writes generic ones)" | M7-03 decision 6 (PASS only) + M7-15 item 5; the checkbox closes when 13 and 14 are green |
| "M7: compute the M4-16 QC metric id_mixups …" | M7-14 rule 3 (line deleted there) |
| "M6: the simulator plants stop-codon readthrough (C-terminal extensions, ADR-0013) for cterm_extension calibration" (re-filed to M7 by M6 lead decision 7) | M7-16 (and M7-12 Q1): useless until the matcher can place stop-spanning peptides |
| "M3-17/M4: resolve makedb protein IDs (GIPPAORF_, NOVEL_, *_patch) through db_id_map.tsv in annotate_protein_hits" | **related, not M7 scope**: it matters for second-pass searches with protein-level input; M7-11 needs the peptide-level decoy pairing only. Stays in the backlog |
| "M2: move the step helpers (_genes, _path, _replicons) from gippa.steps.operons to a public shared module" | done by M6-08 if merged first; otherwise M7-04 does the move (both import them) |
| "M2: one shared IUPAC-aware reverse complement" | M7-01/06 call the existing `translate11_text`/reverse-complement helpers and add no new copy; the line stays |
| "gippa run --annotation <bakta.json> does not stop BAKTA from running under Nextflow" | not touched |

## Risks

- **No traceback in the mapper.** Module 5 recomputes it (M7-01). A disagreement fails the step;
  the conformance criterion on the real sim100k `peptide_hits.tsv` is the early warning.
- **Stop codons are unmatchable** until M7-16a/16b merge (option B, ADR-0058). After them only
  annotated terminal stops are consumable; the fdnG Sec sites (internal TGA) stay unmapped unless
  Khoa accepts ADR-0058 item 6, so on Nissle TE2's `stop_codon` clause stays vacuous and
  PROGRESS.md must say so.
- **Attack 6** (above).
- **RBS threshold from the same run.** It needs ≥ 20 detected genes with OSTIR scores; sim100k may
  sit near that floor, so pass-1 runs may call few intergenic ORFs. Pass 1 accepts that (M7-08 Q4).
- **Second-round runtime.** 221 Promotech windows per candidate start; false peptides create many
  candidates on real data. M7-08 measures it; deduplication against the first round is the
  documented lever.
- **Local-executor strictness.** Registering `ostir_novel` makes every local run with peptides need
  `annotation.precomputed.rbs_novel` (§6.4). M7-08 updates the simulator and the recorded configs in
  the same PR; any config it misses fails loudly, naming the key.

## Feasibility on FarmShare

Pure Python except the second round, which reuses two built conda envs. Memory < 4 GB for a 5-Mb
genome. Slurm is needed for: the two real `gippa_run.sbatch` runs (04, 08), the OSTIR range recording
(08), the full calibration tier (13; one array, ≤ 2 tasks), the TE2 proxy (14, weekly) and the
real-tool arm (15). All within R14 (≤ 32 CPUs at once for M7's own jobs).

## Decisions taken in these specs (for the lead to confirm)

1. Module 5 recomputes the residue traceback; the M3 kernels are untouched (M7-01 Q1).
2. Distinct peptides are counted after I/L folding; an I or L ALT is written `J` (M7-02 Q1–Q2).
3. `shared` peptides support neither REF nor ALT; `multi_identical` ones count at every copy
   (M7-02 Q3).
4. VCF: lexicographic G14 tie-break; substitutions and equal-length two-codon events only;
   origin-crossing codons omitted and flagged `vcf_unrepresentable`; FILTER always `PASS`
   (M7-03 Q1–Q6).
5. `sample_name: auto` is resolved by `gippa run` (M7-04 Q1).
6. Swaps never carry `mass_equivalent`; the mass list is recomputed at run time and the committed
   TSV is pinned to it by a test (M7-05 Q1–Q2).
7. Novel-ORF identity = `(contig, strand, stop)`; class = highest-priority peptide class; ID from
   the longest ORF's lowest coordinate; multi-locus candidate peptides excluded; placed peptides
   checked against the contaminant library (M7-06 Q1–Q4).
8. Second-round OSTIR under the local executor through a new position-keyed precomputed layer
   `rbs_novel` (M7-08 Q1).
9. RBS quantile needs ≥ 20 detected genes; below that, RBS-requiring classes are not called
   (M7-09 Q1); ranks are over called ORFs only (M7-09 Q5).
10. Pass 1 applies §7.6.4's "Neither" branch to every input (`class_fdr = unestimated`, ≥ 2
    peptides) until M7-11 lands; pass 1 never claims error control it does not have.
11. The M7 gate is claimed check by check (G1–G4, A6, TE2); PROGRESS.md records any vacuous clause
    (the `stop_codon` exception while M7-16 is undecided).

## Questions that need Khoa (not the lead)

1. ~~Replicate count for a variants suite~~ resolved by ADR-0042 (20 / 100). New: replicate count
   for the 1-Mb novel-ORF gate, `novel_orfs_gate: 10 / 300` (ADR-0048 proposed; M7-13 Q1).
2. ~~Stop-spanning peptides (M7-16)~~ answered 2026-10-24: B, patch residue `X` (ADR-0058). New:
   ADR-0058 item 6, also consume annotated internal Sec/Pyl codons (fdnG), which B as worded does
   not reach (M7-16 Q1).
3. If attack 6 fails: an ADR changing the `high` rule (M7-13 risks) is a call-rule change.
4. Possibly: the new `annotation.precomputed.rbs_novel` key (M7-08 Q1) and the new flag values, if
   the lead treats them as interface changes in the ADR-0005/0011 sense.

## Lead decisions (2026-10-21)
1. The 11 decisions above, and every per-issue proposed default listed in the designer's report (M7-01 Q1 … M7-15 Q2), are **accepted**.
2. **New flag values** (`vcf_unrepresentable`, `no_start_codon`, `partial`, `ostir_missing`, `rbs_threshold_unestimated`, `junction_unknown`): add them through the schema, as M6 decision 6 did. They are vocabulary, not interface.
3. **`annotation.precomputed.rbs_novel`** (M7-08 Q1) is accepted. It is a config-interface addition, so M7-08 writes an ADR (lead-accepted, Khoa informed); the existing `rbs` layer keeps its keying.
4. **Publishing tables:** the local executor already publishes `tables/` through `PUBLISHED_FILES` (M4-21, #145), and a test pins it to the Nextflow `'key=path'` names. M7-04 and M7-08 add their keys there.
5. **Questions for Khoa** (K1 replicate counts in the protected `tolerances.yaml`, K2 the M7-16 option, K3 the attack-6 rule) stay open. They are pass-2 work, and nothing in pass 1 waits on them. PROGRESS.md records that the `stop_codon` clause is vacuous until K2 is decided.
6. **Integration-pass discipline:** one review round per PR; nits go to the backlog; the wiring issues M7-04 and M7-08 come right after their producers, as sequenced.

## Update (2026-10-24): M7-16 decided
Khoa chose option B ("B"), with the `makedb` patch residue `X`; ADR-0058 records it. M7-16 is now a
buildable stats issue in three PRs (16a, 16b, 16c), about 370 code lines (the row above). Lead decision 5's K2 is answered.
One open Khoa item remains: ADR-0058 item 6 (internal Sec/Pyl codons such as fdnG). Option B as worded does not reach them.
