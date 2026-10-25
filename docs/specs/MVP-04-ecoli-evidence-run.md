# MVP-04: the E. coli K-12 evidence run (real genome, recorded Bakta, simulated peptides)

Track: **MVP / end-to-end** · depends on: M4-24 (#238, unscored decoys → `mayu`), M7-12 (#225,
makedb-style six-frame decoys, `rev_GIAPORF_*` leave gene pairing), M7-11 (#237, novel class FDR),
MVP-01..03 (merged) · related: M7-15 (K-12 real-tool arm, branch `m7/integration-hardening`), M10
(pilot) · estimate ≈ 350 giapsim lines (tooling, `tests/sim/`) + ≈ 120 script lines + ≈ 300 test
lines · no change to `src/giap` is in scope (see Q6).

## Why (Khoa, 2026-10-22, PROGRESS.md › "Answered by Khoa")

The MVP and golden evidence runs on sim100k cannot detect a gene at `stats.alpha.identification` =
0.01. The picked gene FDR is (D + 1)/T, so its floor is 1/T. With T ≈ 74–91 eligible target-won
genes the floor is 0.011–0.014 > 0.01 (MVP-01 runs 1861957/1861972: q = 1/85; MVP real run
1864792: 1/82). Slurm 1869554 (`local/scratch/qfloor/`, 6 seeds × π_wrong ∈ {0.01, 0.05}) showed
that this is not a simulator artefact: with scored null decoys (the `patched` arm) the result is the
same. `n_det_001` = 0, `min_q` = 1/85, `n_decoy_wins` = 0. The opposite failure showed on the
813-gene validation run (1865149): 804/813 detected and q = 1/804 everywhere, which is saturated
and just as uninformative.

**Decision:** the evidence run uses the real E. coli K-12 MG1655 genome and its recorded full-DB
Bakta annotation. giapsim simulates the peptides on that genome until the M10 pilot supplies real
MS data.

## Spec sections implemented

- §7.3.2 detection tiers and the picked gene FDR (steps 1–4, the 80 % pairing rule); ADR-0009 (1/f);
  ADR-0047 (decoys above `q_value_max` leave the competition); M4-24 rules 1–6 (shared key columns,
  `decoy_unscored`).
- §11.3 giapsim: "Strain differences", "Expression and proteomics", "Database-search errors …
  Decoy rows (reversed sequences) carry null scores". This issue gives the decoys scores from the
  null (below), so the null now enters the competition rather than staying outside it.
- §11.4 real genomes (`tests/data/genomes/`, provenance) and the full-DB Bakta fixtures (G17).
- §11.6 TE7 (byte-identical tables for the same seed), as far as it applies here (acceptance 9).
- §7.4 operons, §7.5 variants (`related`, `database_derived`, `strain_difference`), §7.6 novel ORFs
  (classes and class FDR, M7-11), §7.3.5 context. Each module must give a non-vacuous output.
- §3.8 and R14 (Slurm, `giap-*` names, ≤ 128 CPUs, ≤ 20 jobs). R6, R7, R11.
- M7-12 decisions 21, 23, 24 and 26 (makedb-style decoys; the class q floor needs ≥ ~100 candidates).

## Inputs (official sources, hashes; nothing new committed)

| Input | Path | Source and hash |
|---|---|---|
| Genome | `tests/data/genomes/k12_mg1655.fna.gz` (already committed by M2-10 under §11.4) | NCBI Datasets, GCA_000005845.2 (U00096.3, 4,641,652 bp, `[topology=circular]` appended to the header). File sha256 `4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c`. Upper-case sequence sha256 `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`. Cross-checked against the NCBI FTP copy in `$GIAP_DATA_DIR/genomes/` (md5 verified) |
| Annotation | `tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz` (recorded, committed) | Bakta 1.12.2, full DB v6.0, Slurm **1849096** (barley-01, 592 s). JSON sha256 `c08c6466…f083d` (`PROVENANCE.md`), `sequences.nt` stripped. 4,295 CDS + 78 sORF = **4,373 protein genes** (2,140 +, 2,233 −), 6 truncated CDS, no CDS crossing the origin, 2 CRISPR arrays, 2 oriC |
| TU truth (simulation only) | `$GIAP_DATA_DIR/regulondb/14.5.0/TUSet.txt` (not committed, R11) | RegulonDB 14.5.0 (M6-11 fetch, sha256 `c96c8a3b…` in `provenance.json`); see Q3 |
| Databases | `$GIAP_DB_DIR` | geNOMAD, Promotech RF-HOT, antiSMASH, contaminants, COG/GO (`giap setup`), as recorded in PROGRESS.md › Environment notes. Bakta is **not** run (the recorded JSON is `--annotation`) |

Rules:
1. Before simulating, the input step checks both genome sha256 values and the annotation JSON
   sha256 against the table. If one differs, the step exits non-zero with `input_hash_mismatch:<file>`
   and writes nothing (R6: verified artefacts only).
2. **Nothing new is committed:** no genome copy, no genome slice, no simulated peptide table, no run
   output, no RegulonDB-derived table, and no KEGG-derived content (the function module is off, Q8).
   The simulated inputs are regenerated from the seed. Everything lives under
   `$GIAP_ROOT/mvp04/<seed>/` (inputs) and `$GIAP_ROOT/runs/mvp04/<seed>_<jobid>/` (run).
3. Tests may build small excerpts at test time, for example the first 200 kb of the committed FASTA
   and the recorded features inside it. They never write those excerpts into the repository.

## giapsim on a real genome (`tests/sim/giapsim/`, tooling)

The work reuses what exists: `proteomics.simulate_proteomics` (M3-19), `strain.derive_reference`
(M3-18), `layout.layout_from_simulation` + `expression.expression_for_simulation` (M4-08, opt-in,
currently unused), the M7-12 six-frame pool and `rev_GIAPORF_*` naming, `features.novel_orfs`
(out-of-frame / antisense "found" path). The new code goes in one module,
`giapsim/real_genome.py`, plus small opt-in hooks. **Every existing default stays byte-identical**,
so no golden and no pinned test changes. Nothing is imported from `giap`, as the rest of giapsim
already avoids (M7-12, 9d91aa93).

### Interface (proposed names; the implementer may rename within giapsim)
```python
@dataclass(frozen=True)
class RealGenomeConfig:
    proteomics: ProteomicsConfig            # decoys="scored_null", search_db="other", tuned detect_logit
    strain: StrainConfig                    # SNVs only (indel_rate_per_kb = 0), Q4
    expression: ExpressionConfig            # TU-level AR(1), planted islands
    tu_source: Literal["regulondb", "gap_rule"] = "regulondb"
    n_novel_per_class: int = 300            # = M7-12 NOVEL_PER_CLASS_GATE
    novel_len_aa: tuple[int, int] = (30, 150)

def simulate_on_genome(cfg: RealGenomeConfig, seed: int, outdir: Path, *,
                       genome_fasta: Path, bakta_json: Path, tu_set: Path | None) -> RealGenomeResult
# scenarios.k12_evidence() -> RealGenomeConfig ; K12_EVIDENCE_SEED (one fixed seed)
```
Hooks: `ProteomicsConfig.decoys: Literal["unscored", "scored_null"] = "unscored"`, and an optional
`simulate_proteomics(..., log_abundance=None)`. `None` keeps `gene_log_abundance`.

### Pipeline (all draws from `numpy.random.default_rng([seed, <stream>])`, new stream keys)
1. **Read.** Sequence from the FASTA. Features from the recorded JSON (CDS and sORF become `Feat`,
   contig IDs through `contig_map.tsv`, topology from the header). The real genome is never modified.
2. **Sample strain (variants).** giap is given the official K-12 as `--genome`. The sample proteome
   is derived from it: `derive_reference(K-12)` is applied and its output is used as the *sample*.
   The truth rows are swapped so that `ref` = the K-12 base and `alt` = the sample base. With SNVs
   only, sample and K-12 coordinates are identical, so features and makedb IDs lift trivially.
   Written to `truth/alt_sites.tsv` in the M7-12 format. Rate: 0.5 SNV/kb (≈ 2,300 SNVs).
3. **TUs and abundance (context, operons).** TUs are the connected components of same-strand genes
   that share a RegulonDB TU. A Bakta gene matches a RegulonDB gene on the same strand with the same
   stop coordinate. Unmatched genes are singletons. If `TUSet.txt` is absent, the TUs follow the
   giapsim gap rule (adjacent, same strand, gap ≤ 40 bp) and `provenance.json` records
   `tu_source: gap_rule` with a WARNING. A `truth/genes.tsv` (`locus_tag`, `tu_id`) feeds
   `layout_from_simulation`. `expression_for_simulation` then gives the TU-level AR(1) (ρ = 0.5) plus
   **6 planted islands** (3 at +3, 3 at −3 stationary SD, 8 TUs each, spread over the chromosome, not
   overlapping, truth in `truth/islands.tsv`). The natural-log abundance passed to
   `simulate_proteomics` is centred to mean `LOG_ABUNDANCE[0]`.
4. **Novel ORFs (found, not planted).** For each class (`intergenic`, `cds_outframe`, `antisense`),
   `n_novel_per_class` existing ORFs are drawn from the K-12 six-frame. Each starts with ATG, GTG or
   TTG, ends at the first in-frame stop, has a length in `novel_len_aa`, overlaps no annotated protein
   with the same peptide (I = L), and crosses no origin (Q7). These ORFs are declared expressed.
   `intergenic` uses a find-only path, because `features.novel_orfs` writes sequence and the real
   genome must stay unchanged. For `nterm_extension`, CDSs are drawn that have an in-frame upstream
   start codon ≥ 5 codons before the annotated start, with no in-frame stop between. They are
   expressed from that start (`Feat.true_start`), and the annotation is unchanged. On the minus strand
   "upstream" means higher coordinates. If a class runs short of supply, it takes as many as exist and
   logs a WARNING. The realised count is in `truth/novel_orfs.tsv` (M7-12 decision 15/19). Novel ORF
   abundance is drawn from the upper tercile of the gene distribution (Q5), so a class can reach the
   M7-12 decision-26 count of ≥ 100 ORFs with a true peptide.
5. **Peptides (the search that is simulated).** The simulated search mimics DIA-NN/DDA against
   `giap makedb --six-frame-orfs` on the **sample** (Bakta proteins + six-frame ORFs, §15.2 A2-like,
   but on the related sample): `search_db = "other"`, `--genome-relation related`. True peptides use
   the M3-19 logistic, with the intercept tuned (rule 7). Wrong identifications use π_wrong = 0.01.
   They are split between protein entries (shuffled source, attributed to a Bakta locus) and
   noncoding six-frame entries (`sixframe_noncoding`, attributed to a `GIAPORF_*` entry) in
   proportion to the tryptic products of each kind in the search space (Q2). There are 3 runs.
6. **Scored null decoys (`decoys = "scored_null"`).** This extends M7-12 decision 21 from six-frame
   entries to protein entries. For **every** wrong (null) target row, an independent ½ coin per row
   (decision 24) decides whether the PSM stays a target or becomes the decoy of its entry:
   - **Decoy row:** the peptide is reversed (`peptide[::-1]`, as M7-12). `protein_ids` becomes
     `rev_<locus_tag>` for a protein entry and `rev_GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`
     for a six-frame entry. `q_value`, `pep` and `score` are those of the same null draw,
     `is_decoy = true`, and intensity is NA.
   - **Collisions:** if the reversed peptide equals a target peptide (I = L), the row stays a target
     and the count goes to `provenance.json` (`decoy_collisions`).
   - **Same columns:** decoys and targets carry the same key columns (q, PEP, score). By M4-24 the
     shared key is then S = C_T = {q, pep, score} and there is no `decoy_unscored`.
   - **No extra decoys:** the M3-19 `match_wrong` reversed-random-protein decoys are **not** added in
     this mode. They would be decoys with no target competitor, which inflates D. Gene-level pairing
     uses only `rev_<locus_tag>` (f ≈ 1). The `rev_GIAPORF_*` decoys leave gene pairing (M7-12
     decision 23) and feed the M7-11 class FDR.
   - **Known bias:** target q/PEP are computed from the pre-coin wrong share (M7-12 backlog N1,
     ≈ 2× conservative at PSM level). This is accepted and documented, because the gene-level
     competition does not use the mixture q as a probability.
7. **Tuning, once, before the acceptance run.** A giapsim-only Slurm job (no giap) chooses the
   detection intercept b0 (other `detect_logit` terms unchanged). Its target: the share of eligible
   protein genes with ≥ 1 true peptide falls in **[0.55, 0.70]**, near the 2,000–2,600 proteins a
   good E. coli DDA run reports. b0 and the measured share are recorded in the scenario docstring and
   PROGRESS.md. The default (−2.0) is known to saturate (1865149).
8. **Outputs** (`RealGenomeResult`): `peptides.tsv` (C9), `abundance.tsv`, `config.yaml` (only
   `seed` and `sample_name`; no `annotation.precomputed` layers, since the real tools run),
   `truth/{peptides,genes,islands,alt_sites,novel_orfs}.tsv`, `provenance.json` (input hashes,
   giapsim git SHA, seed, b0, tu_source, realised counts, decoy_collisions).
   `sample_genome.fna` is written only under `outdir/truth/` for debugging and is never given to giap.

## The run

```bash
sbatch scripts/slurm/mvp04_inputs.sbatch            # giapsim (1 job, ≤ 4 CPU, ≤ 16 GB, ≤ 30 min)
sbatch scripts/slurm/giap_run.sbatch \
  --genome tests/data/genomes/k12_mg1655.fna.gz \
  --annotation tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz \
  --peptides $GIAP_ROOT/mvp04/$SEED/peptides.tsv --config $GIAP_ROOT/mvp04/$SEED/config.yaml \
  --seed $SEED --db-dir $GIAP_DB_DIR --outdir $GIAP_ROOT/runs/mvp04/${SEED}_<job>
python scripts/dev/mvp04_summary.py <outdir> <inputs dir>   # → <outdir>/evidence_summary.json
```
`giap_run.sbatch` runs with `--profile farmshare,apptainer`. The run uses the defaults:
`stats.alpha.identification` = 0.01, `input.q_value_max` = 0.01, `--genome-relation related`,
modules `core,mobile,bgc` (no `function`, Q8). If `--annotation` does not accept `.json.gz`, the
input step writes the decompressed JSON under `$GIAP_ROOT/mvp04/` (not committed).
`mvp04_summary.py` reads only published outputs plus the giapsim truth, and writes the numbers
below. It is a script, so no GIAP code changes (Q6).

## Where the run sits

- **Acceptance (once, this issue's DoD):** the `giap_run.sbatch` run above, under
  `farmshare,apptainer`. Job IDs go to PROGRESS.md and the measured runtimes to the table there.
- **Recurring: weekly, not nightly.** The nightly array (6 h; slot 3 = real run on sim100k) stays as
  it is. The K-12 run is hours long and its real tools are already exercised weekly. The current
  weekly task 0 (K-12) feeds a 60 % Bakta-protein subset as "evidence stand-in" in step (c). MVP-04
  replaces that stand-in with the MVP-04 inputs at `K12_EVIDENCE_SEED`, runs
  `tests/e2e/test_k12_evidence.py` on the output, and keeps `farmshare_local,apptainer` inside the
  16-CPU / 60 GB / 12 h allocation (weekly decisions 12 and 22: the `farmshare` profile is refused
  there). See Q9.
- **No byte golden.** The full tables of a 4.6-Mb real-tool run depend on container digests and would
  carry Bakta-DB-derived text, so byte-level regression stays with `sim50k` / `sim100k_real` /
  `sim50k_m7`. MVP-04 checks **ranges** instead (acceptance below), as constants in the e2e test (Q10).
- **Runtime budget (R14).** The head job takes 2 CPU / 8 GB. The `farmshare` profile caps the run
  at queueSize 8 × 16 CPU = 128 CPUs, so **before submitting**, the lead sums the running + pending
  `giap-*` CPUs (M7-15's K-12 arm included) and submits only if the total stays ≤ 128. Budget: wall
  ≤ 8 h, CPU time ≤ 150 CPU-h, peak task RSS ≤ 64 GB, giapsim step ≤ 30 min. Known costs:
  TerminatorNet on K-12 took 60 min (1847778). Promoter Calculator runs one task per TU (288 tasks
  for 0.9 Mb, so ≈ 1,400 on K-12), which at 8 concurrent tasks could alone exceed the budget: the
  backlog item "batch it" is a soft dependency (degenerate case D10). If the budget is exceeded, the
  run still counts as evidence but the overrun is recorded and becomes a backlog item. The budget is
  never a reason to drop a module.

## Acceptance criteria and expected numbers

Gates are what the e2e test asserts. "Expected" values are estimates to check against the first
measured run. **After that run, the lead re-sets every expected range from it (± the stated
margin) in a recorded decision and does not widen a gate** (R3).

| # | Quantity | How measured | Expected | Gate |
|---|---|---|---|---|
| 1 | run completes | `giap run` exit code; `run_manifest.json` status | 0, complete | exit 0; no step `failed` |
| 2 | detection tier | `genes.tsv` `detection_tier`; manifest warnings | `decoy_picked` | `decoy_picked`; warnings contain none of `decoy_unscored`, `decoy_pairing_low`, `q_floor_above_alpha`, `no_eligible_genes` |
| 3 | paired share f | `detection_summary.tsv` tier `decoy_paired_share` | ≥ 0.98 | ≥ 0.8 |
| 4 | eligible target-won genes T | 1 / tier `q_floor` | 2,300–3,100 | T ≥ 110 (q_floor ≤ 0.00909) |
| 5 | decoys compete | summary: eligible genes with a paired decoy (q ≤ `q_value_max`) and gene q NA (decoy won) | 20–300 | ≥ 1 decoy win (so it is not the D ≡ 0 picture of 1869554) |
| 6 | genes detected at α 0.01 | tier `n_detected` | 0.85–1.0 × truth-detected genes (≈ 2,100–3,000) | ≥ 1,000 and ≤ 0.90 × eligible (not saturated) |
| 7 | realised gene FDP | detected genes with no true peptide in `truth/peptides.tsv` / detected | ≤ 0.01 | ≤ 0.02 (single run; the calibration suites own the mean) |
| 8 | truth-detected share | genes with ≥ 1 true peptide / eligible (giapsim truth) | in the tuned band | [0.55, 0.70] |
| 9 | determinism (TE7, inputs) | two giapsim calls with the same seed | identical | byte-identical `peptides.tsv`, `truth/*.tsv` |
| 10 | mapping | `mapping_summary` exact / tolerant / unmapped shares | exact ≥ 0.95 | reported; tolerant > 0 (SNV peptides) |
| 11 | variants | `variants.tsv` `high` ALT rows vs `truth/alt_sites.tsv` (same codon, same ALT residue) | 60–300 `high`; precision ≥ 0.97 | ≥ 30 `high`; precision ≥ 0.95; recall ≥ 0.90 among truth sites with ≥ 2 distinct true covering peptides; every row carries `strain_difference` and `database_derived`; VCF validates (pysam) |
| 12 | novel ORFs | `novel_orfs.tsv` calls at class q ≤ 0.01 vs `truth/novel_orfs.tsv` (contig, strand, stop) | per class with ≥ 100 truth ORFs having a true peptide: ≥ 30 calls | ≥ 2 classes non-vacuous (≥ 1 call); called precision ≥ 0.95 per non-vacuous class; a class below 100 is reported as "below 1/α by construction" (M7-12 decision 26) |
| 13 | context | `regions_stats.tsv`, `autocorrelation`, `islands.tsv` vs `truth/islands.tsv` | Moran p = 1/(B+1), I ≈ 0.5; the planted high island recovered (11/11 seeds, Slurm 1872881); lows 0/3 with ceiling ≈ 2.2–2.5 against a critical value ≥ 5; spurious 0 | lag-1 Moran's I > 0 under `tu_block` with p ≤ 0.01; the planted high island (20 TUs, +3 SD) recovered by a `high` call covering ≥ 50 % of its genes; ≤ 2 island calls overlapping no planted island; `regions_stats` has a row for every region type the tools call, with non-NA statistics where the type has ≥ the minimum gene count. The three planted low islands are **exempt (structural cap)**: reported with the ceiling √w·\|Φ⁻¹(q/2)\| (w = scored TUs, q = undetected share of eligible genes from `genes.tsv`), not gated (decision 20) |
| 14 | operons | `pairs.tsv`, `operons.tsv`, `operon_flags.tsv` | ≈ 4,373 adjacent pairs (circular, one replicon); OP pairs ≥ 1,500 | tests not all `none`; no `too_few_op_pairs`; `operons.tsv` ≥ 1,000 TUs; RegulonDB TU recovery reported (not gated) |
| 15 | report | `report/report.html`, manifest | every module card has content | the detection, context, operons, variants, novel ORF cards are present and not the "no data" placeholder; the q-floor caveat is absent |
| 16 | runtime | `scripts/measure_runtimes.py` on the trace | 3–8 h wall | recorded under `mvp04_k12_evidence` in PROGRESS.md › Measured numbers; budget overrun → backlog, not failure |

Checklist:
- [ ] `giapsim/real_genome.py`, the `decoys` and `log_abundance` hooks, `scenarios.k12_evidence`;
  every existing giapsim pin unchanged (the existing tests pass untouched).
- [ ] Unit/property tests (test-author) on a test-time 200-kb excerpt: hash check (example 1), decoy
  coin (examples 2–3), strain swap (example 6), N-term extension on both strands (example 7), find-only
  intergenic ORFs leave the sequence unchanged, unscored mode reproduces main byte-for-byte.
- [ ] A slow test: on the excerpt, `scored_null` decoys under `giap.genetable.detection.detect` give
  `decoy_picked`, no `decoy_unscored`, f = 1 (example 4 shape).
- [ ] `scripts/slurm/mvp04_inputs.sbatch` (`giap-mvp04-inputs`), `scripts/dev/mvp04_summary.py`,
  `tests/e2e/test_k12_evidence.py` (marker `e2e`, reads `$GIAP_MVP04_OUT`, asserts the gates).
- [ ] The tuning job, then the acceptance run: job IDs in PROGRESS.md; `evidence_summary.json` numbers
  in the PR body; every gate passes.
- [ ] `docs/methods/` note (statistician): scored null decoys, why the ½ coin gives the picked FDR its
  null, why `match_wrong` decoys are off, and the 1/T floor and why K-12 clears it.
- [ ] Backlog entries closed: "giapsim real-run scenario needs decoys with null scores" (two entries)
  and "M4: warn when eligible genes < 1/α" (done by `q_floor_above_alpha`). The MVP-01 q-floor
  question is answered.

## Examples

Locus tags, peptides and coordinates in examples 2, 3, 6 and 7 are illustrative. The test-author
picks real ones from the excerpt; the rule each row shows is the one being tested.

| # | Input | Expected |
|---|---|---|
| 1 | `k12_mg1655.fna.gz` with one byte changed in the sequence (test-time copy) | input step exits non-zero with `input_hash_mismatch:k12_mg1655.fna.gz`; outdir has no `peptides.tsv` |
| 2 | one wrong target row: peptide `AGLVTEPK` (shuffled source, attributed to locus `K12_00420`), run2, score 0.513, q 0.0049, PEP 0.31; coin → decoy | that row is replaced by `KPETVLGA`, `rev_K12_00420`, q 0.0049, PEP 0.31, score 0.513, `is_decoy = true`, intensity NA; the target row is gone from run2 (other runs draw their own coin) |
| 3 | six-frame null PSM on entry `GIAPORF_antisense_U00096.3_120001_120150_m`, coin → decoy | decoy `rev_GIAPORF_antisense_U00096.3_120001_120150_m`; `pair_decoys` drops it from gene pairing (M7-12 decision 23); it counts in the M7-11 antisense class FDR |
| 4 | `decoys = "unscored"` (main behaviour) on the same seed | gene tier `mayu` with `decoy_unscored` (M4-24 rule 3); acceptance 2 fails: this is the guard that the scored mode is the one in use |
| 5 | T = 2,600 eligible target-won genes, D = 0 above the top 2,000 | q_floor = 1/2,600 = 0.000385 ≤ 0.01; no `q_floor_above_alpha`. The sim100k contrast, T = 85: q_floor 0.0118 > 0.01 and the warning |
| 6 | K-12 base C at U00096.3:337,021 (+ strand CDS, codon 2 = `CCG`, Pro); the sample SNV makes it `CTG` (Leu) | `truth/alt_sites.tsv`: ref C, alt T, pos 337021 on both sides, `aa_change` P?L with codon index. Two detected sample peptides cover it, so `variants.tsv` has an ALT `high` call `AA_REF=P`, `AA_ALT=L` with flags `snv;strain_difference;database_derived`, and the VCF record has `REF=C`, `ALT=T` at the codon's lowest position |
| 7 | a minus-strand CDS annotated on [s, e) whose in-frame start codon lies 12 codons above e (higher coordinate) with no in-frame stop between | `nterm_extension` truth ORF [s, e + 36); the protein is expressed from the upstream start; an N-terminal peptide maps across the annotated start; the recorded annotation is unchanged |
| 8 | b0 = −2.0 (the M3-19 default) | truth-detected share ≈ 0.94 (measured 0.944), `n_detected` ≥ 0.95 × eligible: gates 6 and 8 fail (saturated, as in 1865149). The tuned b0 is what the scenario pins |
| 9 | `TUSet.txt` absent | `tu_source = gap_rule`, WARNING, run proceeds; acceptance 14's RegulonDB recovery line is "not available" |
| 10 | island 4 planted at −3 SD over 8 TUs (≈ 15 genes), few detected peptides there | `islands.tsv` has a down island overlapping ≥ 50 % of those genes (M5-11: a down island over peptide-less genes is allowed) |

## Degenerate cases

- **D1, circular vs linear.** K-12 is one circular replicon (the header tag). If the tag is lost, giap
  reads the replicon as linear (PROGRESS.md environment note), the wrap pair disappears and the
  input-hash check fails first, because the file hash covers the header. No CDS crosses the origin.
  Six-frame candidates that cross it are excluded from the planted set (Q7) but remain in the
  simulated search space, so null hits may still land on them.
- **D2, minus strand.** About half the genes (2,233), the minus-strand extensions (example 7), the VCF
  REF/ALT reverse-complement (§7.5.5) and antisense ORFs are all covered by the run's truth checks.
- **D3, empty or near-empty classes.** If a novel class's supply or detection leaves < 100 truth
  ORFs with a true peptide, it is reported as vacuous by construction, not as a failure (gate 12
  needs ≥ 2 non-vacuous classes). Zero SNVs covered (impossible at 0.5/kb) would fail gate 11.
- **D4, no decoy wins.** D ≡ 0 is legitimate in principle, but here it signals that the decoys are not
  competing (the 1869554 picture), so gate 5 fails.
- **D5, pairing below 80 %.** This would come from `match_wrong` decoys left on, or from
  `rev_GIAPORF_*` decoys reaching gene pairing (M7-12 not merged): the tier is `mayu` with
  `decoy_pairing_low`, and gate 2 fails.
- **D6, ties.** Scores are rounded to 6 decimals, so target–decoy ties are rare. A tie goes to the
  decoy (§7.3.2 step 3). A decoy and its target can never come from the same draw, because the coin
  sends each null PSM to one side only. Palindromic reversals are handled by the collision rule
  (rule 6).
- **D7, shared peptides.** Paralogs (for example, the IS elements and the *tuf* genes) give peptides
  with several `protein_ids`. A decoy of a multi-locus wrong peptide uses the first locus only. The
  mapping module's ambiguity handling is unchanged.
- **D8, partial genes.** The 6 truncated CDS are `partial_excluded`. They do not count in T or in the
  eligible denominator.
- **D9, empty peptide input** (a giapsim bug writes a header only): `giap run` takes the no-peptide
  path, and gates 2 and 6 fail. The input step refuses to write a `peptides.tsv` with no target rows.
- **D10, runtime.** If Promoter Calculator is not batched, the run may exceed 8 h. That is recorded,
  not failed. If FarmShare holds the jobs ("user env retrieval failed"), the watchdog releases them.
- **D11, decoys above `q_value_max`.** Under ADR-0047 they leave pairing and picking. At π_wrong 0.01
  most nulls have mixture q ≤ 0.01, so most decoys survive. The surviving count is in the summary.

## Relation to M7-15 and M10

- **M7-15 item 4** (implementer running now) is a *real-tool* K-12 arm with the in-silico TE8 peptide
  set (UniProt K-12 tryptic peptides). Those peptides have no scores and no decoys, so that run
  exercises the tools and records `novel_second_round_k12` runtimes, but it says nothing at α. MVP-04
  is the *evidence* arm: the same genome, annotation, sbatch path and tool envs, with scored peptides
  and truth. The two do not merge.
  - MVP-04 reuses M7-15's run shape and its `$GIAP_M7_15_K12_OUT` lessons.
  - MVP-04 records its own runtime key (`mvp04_k12_evidence`) and does not re-record M7-15's.
  - The two runs are never submitted while their combined CPUs could pass 128.
- **M10.** When the pilot's A1/A2/A3 data exist (Nissle, DIA-NN, Cascadia), the real-MS runs become
  the evidence. MVP-04 then stays as the **truth-bearing** weekly regression (real genome, known
  answers), which the pilot cannot provide. The `scored_null` decoy mode is giapsim's model of DIA-NN
  `--report-decoys` output, and the methods note says so. If the pilot shows a different decoy shape,
  for example decoys with q but no PEP, M4-24 already handles it. giapsim would then gain a matching
  mode.

## Questions for the lead

1. **"Never commit the FASTA" vs §11.4.** The K-12 FASTA is already committed (M2-10,
   `tests/data/genomes/`, as §11.4 requires; NCBI places no restrictions on it). **Default:** read the
   committed file by hash and commit nothing new that is derived from the genome (no slices, no
   simulated tables). Reversible.
2. **How null PSMs are split between protein entries and six-frame entries.** **Default:** proportional
   to the tryptic products of each kind in the simulated makedb search space (the realistic choice;
   most nulls then land on six-frame entries). Alternative: a fixed 50/50. Either way, gate 5 needs
   enough gene-level nulls. If the proportional split gives fewer than ~50 gene-level null PSMs, the
   default switches to 50/50 and the switch is recorded.
3. **TU truth from RegulonDB (local, not committed) or the gap rule.** **Default:** RegulonDB 14.5.0
   TUs (independent truth, already fetched under ADR-0022's citation), with the gap rule as the
   fallback. The giapsim output then depends on a non-committed file, but only on the cluster, and its
   hash goes into `provenance.json`.
4. **SNVs only, or indels too?** **Default:** SNVs only at 0.5/kb, which keeps coordinates identical
   (as M7-15 decision 3 does). Indels and frameshifts stay with the `sim50k_m7` golden and the M7
   calibration.
5. **Novel ORF abundance.** Short ORFs at the tuned intercept rarely reach the ≥ 100-with-a-peptide
   count. **Default:** novel ORFs draw abundance from the upper tercile and lengths of 30–150 aa,
   documented as a simulation choice. Alternative: plant 1,000 per class (supply-limited).
6. **`n_decoy_wins` as a tier row?** Gate 5 is computed by the summary script from inputs and outputs.
   **Default:** no GIAP change. A `tier` row `n_decoy_wins` in `detection_summary.tsv` would be new
   §9.3 vocabulary, so it would go to an ADR and Khoa, and is **pending Khoa** if the lead wants it.
7. **Origin-crossing novel ORFs** are excluded from the planted set (the makedb naming of such
   entries, with `end1 < start1`, is untested). **Default:** exclude them and add a backlog item.
8. **The function module (GO/KEGG/COG, GSEA).** It is not in the list of modules the run must show,
   and its outputs carry KEGG-derived content that may never be committed. **Default:** off in the
   acceptance run and in the weekly slot. An optional second weekly arm may turn it on, with outputs
   kept under `$GIAP_ROOT` only.
9. **Weekly slot.** **Default:** replace weekly task 0's step (c) evidence stand-in (the 60 % Bakta
   protein subset) with the MVP-04 inputs and the e2e range test. If `weekly.sbatch` is protected or
   Khoa-owned, add a separate array task 2 instead. **Pending Khoa** only if the scrontab must change.
10. **Range gates as protected constants.** **Default:** after the first measured run, the lead sets
    the expected ranges and gates (table above) in the e2e test, and the file goes on
    `ci/protected.txt`, so that a later widening needs Khoa (R3). Adding the protection tightens; it
    does not weaken. Reversible until then.
11. **Seed.** **Default:** one fixed `K12_EVIDENCE_SEED` for acceptance and weekly. A 3-seed spread
    (to see the variance of gates 5–7 and 11–13) is one extra Slurm array, run once and recorded, not
    gated.
12. **TE7 on the full run** (two runs, same seed, byte-identical `tables/`) would double the cluster
    time. **Default:** TE7 on the inputs here (gate 9). The full-run TE7 stays with the integration
    runs (sim100k) and M7-15 item 7.

## Lead decisions (2026-10-23)
1. **Questions 1–5, 7, 8, 10–12:** the proposed defaults are accepted as written.
   - The committed `tests/data/genomes/k12_mg1655.fna.gz` (M2-10, §11.4) is used by hash, and
     nothing derived from it is committed.
   - Null hits split by search space, with a 50/50 fallback below about 50 gene-level nulls.
   - RegulonDB 14.5.0 truth TUs.
   - SNVs only, at 0.5/kb.
   - Upper-tercile novel ORFs of 30–150 aa.
   - Origin-crossing novel ORFs are excluded, with a backlog item.
   - The function module is off.
   - Range gates become constants after the first measured run.
   - One fixed seed; a 3-seed spread is recorded, not gated.
   - TE7 is checked on the inputs only.
2. **Question 6 (`n_decoy_wins` tier row): pending Khoa.**
   - Default for now: no GIAP change. The summary script computes the number.
   - Logged in PROGRESS.md › Questions for Khoa; not blocking.
3. **Question 9 (weekly slot):**
   - Default: add weekly array task 2 if `weekly.sbatch` is protected; otherwise replace the
     step (c) stand-in.
   - Changing the scrontab is pending Khoa; not blocking, and the first run is a one-off sbatch.
4. **Order:** the test-author starts after #225 (M7-12) merges, because `scored_null` extends its
   decision-21 coin. #238 and #237 must merge before the implementer's end-to-end run, not before
   the tests.
5. **Runtime:** before submitting, sum the giap-* CPUs together with M7-15's K-12 arm. Promoter
   Calculator batching (backlog) is needed if the first run exceeds 8 h; that is measured, not
   assumed.
6. **Test-author part 1 (54ac5d5e), spec gaps 1–8:** the test-author's readings are accepted.
   1. RegulonDB gene coordinates come from `GeneProductSet.txt` next to `TUSet.txt`.
   2. "Same stop" is read strand-aware: the right end on +, the left end on −.
   3. `truth/alt_sites.tsv` uses the M7-12 columns. The SNV-only checks compare
      `truth/sample_genome.fna` with K-12.
   4. `simulate_on_genome(..., check_hashes=True)`; test excerpts pass `False`. The CLI is
      `python -m giapsim.real_genome --genome --annotation --seed --outdir`.
   5. Island truth names layout `TU_<n>` ids; excerpts pass their own island specs.
   6. Upper tercile means at least the 2/3 quantile of gene abundances in `abundance.tsv`.
   7. The `nterm_extension` count is bounded by `n_novel_per_class`.
   8. The branch is stacked on #225 (merged in by the lead). When #225 squash-merges, merge its old
      tip and then main, following the stacked-merge rule.
7. **Implementer part 1 (4a16aa79):** the tuned b0 is −5.25, which gives a truth-detected share of
   0.631 (2,757 / 4,367), and every novel class has more than 200 ORFs with a true peptide. The
   implementer's unpinned choices 1–6 are accepted:
   1. Wrong identifications stay `shuffled` until part 2.
   2. `alt_sites` leaves out SNVs in the first codon (giap reads it as M).
   3. The replicon is passed to the layout as linear; no TU joins across the origin.
   4. Novel ORFs also exclude overlaps with an RNA gene or a same-strand pseudogene, and must stay
      intact on the sample (as in makedb).
   5. `proteins.faa` is dropped; `truth/novel_orfs.json` and `truth/proteomics.json` are added.
   6. Provenance records `decoys: "unscored"` until part 2.

   Example 8's expected share is corrected to 0.94; the test's ≥ 0.90 gate stands.
8. **Part 2** (`decoys="scored_null"`, question 2's split between protein and six-frame entries):
   test-author next, then the implementer, then the end-to-end run once #238, #225 and #237 have
   merged.
9. **Test-author part 2 (915a565e), spec gaps 1–7:**
   1. **Gap 1, accepted as a design change.** Protein-entry nulls must be able to give false gene
      evidence, otherwise the gene-level null is one-sided (conservative only) and gate 7 is vacuous.
      - In `scored_null` mode a protein-entry null is a real tryptic product (within the digest
        rules) of a recorded locus that is *not* detected in truth. Its target row carries that
        locus tag; on the coin, its decoy is `rev_<that locus>`.
      - The locus is drawn uniformly over loci with an undetected eligible product.
      - Shuffled peptides are not used in this mode.
      - The test-author adds tests for this (target maps to the genome and the locus; that locus has
        no true peptide; gate 7's false-discovery share becomes estimable). Then the implementer.
   2. "About 50" counts gene-level null PSMs before the coin. The two clear cases are tested.
   3. The search-space share uses distinct tryptic products, protein entries against six-frame-only
      products, with makedb `min_orf_aa` = 7 and ±0.04.
   4. Examples 4 and 3b activate once #238 and #237 merge (backlog).
   5. `decoy_collisions` counts PSMs.
   6. M7-12 decisions 28–29 exist on #225's current head; this branch merges it again before the
      implementer starts.
   7. Add `.hypothesis/` to `.gitignore` if missing (lead).
10. **Statistician review (09bb875d, review on #247, Slurm 1871592 and 1871626): BLOCKED on B1–B3.**
    Decisions 10–14 answer it.
    - **B1: the null split becomes a fixed 50/50 (question 2's alternative). The proportional
      default and its "~50 PSMs" fallback are withdrawn.**
      - Under the search-space split, K-12 gets about 21 protein nulls and 11 decoy wins.
      - The α cut never binds, and an engine that ignores the decoys still passes gate 7 (FDP 0.0076).
      - The fallback was also non-monotone in π_wrong.
      - With 50/50 (10 seeds): about 83 protein nulls and 41 decoy wins. The α 0.01 cut binds on
        10/10 seeds, with FDP 0.0104 ± 0.0013. Ignoring the decoys gives 0.030, which fails gate 7.
      - Six-frame entries get the other half. The class-level null was already dense and stays
        informative.
      - Implementer: change `_scored_nulls` / `k12_evidence`.
      - Test-author: update the split test and the gate-5 expected range (with 50/50, about 41 wins
        fits 20–300), and remove the statistician's two strict xfails ("gate 7 can fail", "α cut
        binds"). These become plain slow tests that must pass.
11. **B3: what gate 7 measures.** A "true peptide" is a `truth/peptides.tsv` row with `is_wrong =
    false`. A protein null's truth row carries its locus tag but has `is_wrong = true`, so it never
    makes a gene "true".
    - `mvp04_summary.py` uses this definition. It also reports false-detected genes by channel:
      protein null, six-frame null on a `GIAPORF_nterm_extension_*` entry, and tolerant-route
      (one-mismatch) six-frame null. That covers non-blocking item 9.
    - Test-author: the `gene_competition` oracle gains the nterm-extension channel and
      `multi_identical` acceptance.
12. **Test-author items from the review:**
    - (B2) The `gene_fdr` suite's replicate counts follow the suite convention: the full tier is
      R = 300, not `replicates_fast`.
    - (5) The two binomial tests of counts that are fixed by construction become exact-equality
      asserts.
    - (6) Teeth for the KS checks and `picked_null` (a score-dependent coin), and for `coin_pooled`
      (a p = 0.55 coin). Raise `picked_null`'s power against P = 0.6, or record its power.
    - (7) Use the calibration harness (`confirm` with records, `suite_seed`) rather than the file's
      own `confirm()`. Print each estimate with its 95 % CI.
13. **Recorded deviation (item 8).** Decision 9's protein nulls are drawn uniformly over undetected
    loci. The §8.4 gene-FDR generator draws size-proportionally. This is an evidence run, not a
    calibration gate, and the calibration suites keep the §8.4 generator. So the deviation is noted
    here and in the methods note, and not changed.
14. **Backlog (item 4, a GIAP issue, not MVP-04).**
    - nterm-extension decoys pair with a gene only after a target row is accepted (M7-12 decision
      28). A null gene whose rows all went to the decoy side then vanishes instead of counting as a
      decoy win. That is anti-conservative: about +0.8 false genes at K-12, more on real searches.
      The fix is to pair through makedb's entry-to-gene relation.
    - Also from #238's statistician: tighten ADR-0047's q-filter condition from "some target
      carries q" to "every target carries q".
    - Order: the test-author goes first, then the implementer, then the code review. The end-to-end
      run follows after #238 and #237 merge.
15. **Test-author round 3 (fc29cb8d): the readings of gaps 1–4 are accepted.**
    - Gap 1: the summary interface is `gate7(outdir, inputs) -> dict`, and its channels are not
      exclusive.
    - Gap 2: the protein side is capped by the pool and the six-frame side takes the rest.
    - Gap 3: the mode is `fifty_fifty`, and the search-space share is recorded for information only.
    - Gap 4: the pooled checks skip under sharding.
    - Gate 5's expected range stays at 20–300.
    - Next, the implementer: `_scored_nulls` is always 50/50, and a new
      `scripts/dev/mvp04_summary.py` provides `gate7()`, with the result under `gate7` in
      `evidence_summary.json`. After that, the code review.
16. **Implementer round 3 (225a2b1c): accepted.**
    - **Wider protein-null pool:** loci the sample strain truncates with a premature stop (giapsim
      `pseudo`) are now in the protein-null pool. Products spanning the stop are skipped.
      - GIAP treats these loci as eligible genes, and decision 9's "recorded locus not detected in
        truth" covers them. The test oracle already counted them.
      - This changes which loci each scored_null seed draws. No pins depend on it.
    - **K-12 measurement** (Slurm 1872190):
      - 174 null peptides, 87 protein and 87 six-frame;
      - 46 decoy wins out of 77 paired genes;
      - T = 2,744, of which 2,721 are detected at α (the cut binds);
      - picked FDP 0.0077; a decoy-blind engine gives 0.0284, which fails gate 7 as intended.
    - **Tests:**
      - 1872171: the fast tier of the K-12 integration, attack and calibration tests, 17 passed.
      - 1872199: the two full-tier attack tests, 2 passed.
      - 1872172: the fast suite, 14,035 passed.
    - The lead fixed a stale docstring in the round-3 tests (the mode is `fifty_fifty`).
    - Next: the statistician re-checks B1–B3, and the code review runs in parallel.
17. **Statistician re-check (6fcc3d95, Slurm 1872374): APPROVE.**
    - **Gate 7 catches only gross failures.** Over 10 seeds a correct engine's FDP is 0.0109 ±
      0.0013, and a decoy-blind engine reaches 0.023–0.029, which is 1.2–1.4 × the gate. A 1.5 × FDR
      inflation would usually pass. The calibration suites and the dense-null attack own the mean,
      so gate 7 is not widened.
    - **Q11's 3-seed spread** records both the FDP and the decoy-blind FDP.
    - The lead refreshed the numbers in the attack docstring.
18. **Code review (#247): CHANGES REQUESTED on B1–B3, all fixed by the lead.**
    - B1: the E501 in the lead's docstring refresh is fixed.
    - B2: `fifty_fifty_violations` asserts that the realised null count equals
      `n_wrong_requested`. This kills the six-frame shortfall mutant.
    - B3: `test_protein_nulls_skip_products_across_a_stop` and
      `test_every_peptide_uses_the_twenty_amino_acids` kill the removed-stop-filter mutant.
    - Non-blocking:
      - (2) Decision 9 is reworded. A protein null is a tryptic product of a recorded locus *as
        expressed in the sample*, which includes the planted N-terminal extensions. Two K-12 nulls
        lie in extensions, and gate 7 counts them as `protein_null`.
      - (3) Size: about 1,150 lines in tests/sim and scripts. This deviation is recorded rather
        than split, because splitting would discard the reviews already done.
      - This PR is **MVP-04 part 1**: the simulator, the scored nulls and `gate7()`. **Part 2**
        is a new PR: `scripts/slurm/mvp04_inputs.sbatch`, the full `mvp04_summary.py` (gates
        1–16), `tests/e2e/test_k12_evidence.py`, and the acceptance run once #238 and #237
        merge.
      - (1, 5) go to the backlog.
19. **Part 2 test-author (2a95422f, #258): the readings of gaps 1, 2 and 4–10 are accepted.**
    - Gap 2: detection warnings are logged as plain text and carry no flag or event. This goes to
      the backlog as a GIAP change: gene_table logs `decoy_pairing_low` / `decoy_unscored` /
      `q_floor_above_alpha` with an event, so the manifest carries them.
    - **Gap 3 (gate 13 islands, "≥ 4 of 6 recovered") may be infeasible.** Low islands are
      structurally capped once more than about a quarter of genes are undetected, and K-12 has about
      37 % undetected.
      - The statistician rules on it before the acceptance run. That covers the islands
        (direction-matched, as the test-author reads it) and the lag-1 Moran's I under the
        `tu_block` null (gap 4).
      - If the gate is infeasible, the lead re-states it before the first measured run (R3: the
        gate is changed by a recorded decision before any number is seen, not widened after).
    - Next steps, in parallel: the implementer (the full `mvp04_summary.py` and
      `scripts/slurm/mvp04_inputs.sbatch`) and the statistician's gate-13 feasibility check.
20. **Gate 13 re-stated before any acceptance number (R3), from the statistician's feasibility
    ruling (Slurm 1872831, 1872881; #258 comment 6365158721).**
    - "≥ 4 of 6 islands" was infeasible: 0/11 seeds (Wilson [0, 0.26]). The 8-TU lows are capped
      (ceiling 2.2–2.5 against a critical value of 6.3–7.4), and three 8-TU highs mask each other
      in the block-bootstrap maxT null (3/33 recovered).
    - `K12_ISLANDS` becomes one high island of 20 TUs at +3 SD from TU 970, plus three 8-TU lows
      at −3 SD from TUs 560, 1380 and 2200 (structural-zero exhibits). Measured: the high island
      recovered 11/11 in both detection arms; pooled with the no-lows variant 21/22 (Wilson
      [0.78, 0.99]), so a single-seed gate carries about a 5 % risk of a miss by masking. Recorded,
      not hidden. The detected share is unchanged, so gate 8 is unaffected.
    - Gate 13's criterion is the table row above. The summary gains `n_recovered_high`,
      `low_exempt` and the per-low ceilings; the e2e `MIN_ISLANDS` follows.
    - Rejected fallback: recovery recorded and not gated (a no-op islands step would pass).
    - Methods: the masking paragraph goes in `docs/methods/islands.md`; the method fix (ℓ̂ with
      top windows masked, or a step-down scan) is backlog, needing an ADR-0021 amendment.
    - Reversible (test scenario only); FYI for Khoa in PROGRESS › Questions.
21. **Part 2 implementer (cac56aa7): readings accepted, size deviation recorded.**
    - Gate 12 matches a call against truth ORFs of any class on (contig, strand, stop), per the
      spec wording; gate 13 fails `regions_ok` when `annotation/regions.tsv` is missing or empty;
      gate 13 keys planted TUs by the `TU_<n>` suffix without a replicon filter (K-12 has one);
      gate 16 labels node class from the trace hostname, "?" otherwise; gate 9 with no second
      call returns `None` and fails.
    - Size: 494 counted script lines against an estimate of about 120. Accepted as for part 1
      (decision 18): the 16 gate functions are one unit and splitting them would leave gates
      untested in the acceptance path.
    - Next: a test-author round for decision 20, then the implementer (`K12_ISLANDS`, the
      `k12_evidence` docstring, the island part of `gate13()` and its summary fields), then code
      review and the statistician's check of the summary, then the acceptance run.
22. **Gate-13 test-author (8b3895c9): gaps 1–6 accepted as pinned.**
    - `value` = `n_recovered_high`, expected [1, 1]; `in_expected` also needs `n_spurious == 0`
      (1–2 spurious calls are `outside_expected`, not a failure).
    - `low_ceilings` = {island_id: {w, q, ceiling}}, lows only.
    - w = the island's TUs with at least one eligible gene (`tables/genes.tsv`, `detection_status`
      not `partial_excluded`); q = the replicon's share of eligible genes not detected, from GIAP's
      `tables/genes.tsv`; the ceiling is the raw `low_direction_ceiling`, before LOESS attenuation.
    - The high-island part means every planted high island is recovered.
    - Next: the implementer (`K12_ISLANDS`, the `k12_evidence` docstring, the island part of
      `gate13()`).
23. **Statistician review of part 2 (d14f0468): APPROVE** (#258 comment 6365837248). Slurm 1873211:
    the scored-null calibration files pass (14), since they use their own 8-TU island excerpt.
    - The four methods-note corrections (M1–M4) are applied to `docs/methods/islands.md`.
    - N1: gate 7's expected band [0, 0.01] sits below a correct engine's 10-seed mean of
      0.0109 ± 0.0013. It stays as pre-registered; `outside_expected` is a warning, not a failure.
      When the ranges are re-set after the first run, the statistician's [0.008, 0.014] is the
      candidate.
    - N2–N6 (eligible vs scored genes in q and w, detrended S past the ceiling, gate 5's
      decoy-only wins, gate 12 precision risk at few calls, no CI on Moran's I and gate-7 FDP):
      recorded here, not changed.
24. **Code review of part 2 (d14f0468): APPROVE** (#258 comment 6366000222). 31 of 49 summary
    mutants were killed; the surviving threshold mutants are re-asserted by
    `tests/e2e/test_k12_evidence.py`'s own constants. One review round: the nits go to the backlog.
    - Backlog: boundary tests for gates 4, 6, 8, 11 and 13 (p) and a gate-13 case with no planted
      high island; the sbatch seed from `K12_EVIDENCE_SEED`, renamed `MVP04_SEED`, and tested as
      passed through; gate 15 counts the `empty_state` placeholder SVG as content (N2, a weaker
      check than "a figure or a table"); N3–N9 tidying.
    - **N10, binding for the acceptance run:** `mvp04_inputs.sbatch` runs the code in
      `$GIAP_REPO`. The acceptance inputs are made only after #258 merges and the main checkout is
      fast-forwarded, so they carry decision 20's `K12_ISLANDS`.

## Amendment (2026-10-24, Khoa: "bug fix"): how gate 11 is scored
Statistician audit stat/m7-k12-flags 2e9ecbf4 (docs/methods/variants-mvp04-k12-audit.md). Khoa
ruled these four scorer and truth fixes plain bug fixes, with no ADR and no sign-off caveat.
Thresholds and gate text are unchanged.
1. Residues compare with I/L folded on both sides (ALT `J` = truth `I` or `L`).
2. Initiator rows leave the precision denominator and are reported as `n_initiator_m_excluded`,
   because giapsim leaves start codons out of truth. An initiator row is an ALT M row that either
   carries `start_codon` or sits at a GTG/TTG/CTG truth novel-ORF start.
3. Truth is `truth/alt_sites.tsv` plus `truth/novel_alt_sites.tsv`: interior SNV codons of the
   planted intergenic, cds_outframe and antisense ORFs.
4. A recall site needs ≥ 2 distinct (folded) correct peptides of at least
   `mapping.tolerant.min_len` (from `OUTDIR/resolved_config.yaml`) whose interval covers the
   codon. With no readable config, recall is NA and the gate fails.

On run 20261023_a1 the literal text gives precision 55/70 and recall 55/65; the scored values are
55/55 and 55/58. Ten of the 15 excluded initiator rows carry no `start_codon` flag. The ADR
candidate "flag ALT M at unannotated start-capable codons" is in the backlog.

## Clarification (lead, 2026-10-24; statistician note local/stat/mvp04-gates-2-12.md): gate 2 tokens
Gate 2's forbidden warnings are the gene-level detection tokens. They match as whole tokens, so the
novel-ORF step's per-class `novel_orfs.class_q_floor_above_alpha` (the class FDR floor, §7.6;
expected for small classes such as nterm_extension) does not fail gate 2. This is a clarification,
not a gate change.
