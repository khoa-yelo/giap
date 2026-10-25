# GIAP preprint outline (M11-11)

An **outline**, not a draft (`docs/specs/M11-11-preprint-outline.md`): headings, one-line intents,
the figure and table plan, and a claim → evidence map. No result appears here that is not in a
committed evidence file; every other number is `TBD` and its row is `pending`, naming the issue
that will produce it. Kept in the repository and left out of the site navigation (Q-L6 style).
Aggregates only (Khoa's 2026-10-24 pilot data rule): no peptide sequences, no DIA-NN rows.

**Rule for the paper (R3).** The text may claim only what the gates and the readiness table
(`validation/release/1.0.0/readiness.json`, M11-03 and M11-13) say. A gate that is waived, vacuous,
`deferred` or `pending` is reported as such, never as a pass. A threshold is never restated more
loosely than `ci/gates.yaml` or `tests/calibration/tolerances.yaml` define it.

## 1. Front matter

| Item | Value |
|---|---|
| Working title | GIAP: genome-informed analysis of bacterial proteomes, with error control for every claim |
| Venue | **TBD (Khoa, Q-K11).** Default in use: bioRxiv preprint, then *Bioinformatics* (Application Note or Original Paper) |
| Authors and affiliations | **TBD (Khoa).** Default in use: Khoa Hoang, sole author; no ORCID or e-mail until Khoa gives them (R12) |
| AI-assistance statement | **Placeholder: Khoa to write.** The repository records every agent role in commit trailers (`Agent-Role:`); the statement's wording and scope are Khoa's |
| Figure permissions | **TBD (Khoa, Q-K11):** which pilot figures may appear (pilot data are PXD078459, public) |
| Competing interests, funding | **TBD (Khoa)** |

## 2. Sections (intent only)

1. **Abstract** (structure only): problem; what GIAP does (one gene table keyed by Bakta locus
   tags); how errors are controlled; what the Nissle 1917 pilot showed; availability. Numbers only
   from rows marked `supported` in §4 at submission time.
2. **Introduction.** The question: how much of what a bacterial proteome detects, and how much it
   holds, is explained by genome context (operons, islands, mobile elements, codon usage,
   promoters), and how to ask that without the database mismatch of a generic reference.
3. **Design.** Inputs (genome, Bakta annotation, peptide or protein tables from any search engine);
   one gene table keyed by Bakta locus tags; the modules (annotation layers, mapping, gene table and
   detectability, region tests, autocorrelation and islands, operons, variants, novel ORFs,
   function, determinants, report, compare). Tools run as separate programs (§5 below).
4. **Error control.** One paragraph per family with its estimand and its null, each pointing at
   its methods note: placement FDR (`docs/methods/placement-fdr.md`), gene-level FDR, picked and
   MAYU (`docs/methods/detection-tiers.md`, `docs/methods/mayu-gene-fdr.md`), circular-shift
   region tests (`docs/methods/region-tests.md`, `docs/methods/permutation-and-shifts.md`),
   island FWER (`docs/methods/islands.md`), operon decoys (`docs/methods/operons.md`), BH
   families (`docs/methods/multiple-testing.md`), the attack nulls and the calibration harness
   (`docs/methods/calibration.md`).
5. **Validation.** (a) Simulation calibration with the giapsim simulator (the four one-sided
   criteria of `docs/methods/calibration.md`); (b) the Nissle 1917 pilot with arms A1 (generic K-12
   reference), A2 (matched strain) and A3 (de novo); (c) the true-example assertions TE1–TE8
   (GIAP_SPEC.md §11.6).
6. **Performance.** Measured wall time and peak memory on FarmShare for the standard run on
   16 threads and 64 GB, with the node class; per-tool runtimes.
7. **Availability.** PyPI and Bioconda packages, containers pinned by digest, Zenodo DOI, the
   documentation site, licences (§5 below).
8. **Limitations.** What version 1 defers (§6 below), stated before the reader has to ask.

## 3. Figure and table plan

Run directories are on FarmShare under `$GIAP_ROOT` (not committed, R11); the figures use the
pilot freeze runs of M10-14 (`$GIAP_ROOT/runs/pilot/freeze/`).

| Figure | Plot ID | Comes from | Claim it supports (§4 row) | Status |
|---|---|---|---|---|
| Fig. 1 | new (workflow schematic) | drawn by Khoa or adapted from the docs site (M11-07, M11-08) | the design (no claim) | pending (Khoa; M11-08) |
| Fig. 2 | V01 `circos` | the A2 freeze run, `plots/` | where in the genome the A2 proteome lies (descriptive; C19) | pending (M10-04 #296, M10-14) |
| Fig. 3 | new panels: QQ of null p-values, observed FDP vs. nominal q | the `calibration_full` records on the release-candidate commit (`$GIAP_ROOT/logs/calibration/`) | C3, C4, C11–C14, C17 | pending (M11-13); no plotting code for these panels exists yet |
| Fig. 4 | V07 `obs_vs_exp`, V08 `cog_abundance` | the A1 and A2 freeze runs | detection against detectability, by region and category (C19, C20) | pending (M10-04 #296) |
| Fig. 5 | V11 `operon_pairs` | the A2 freeze run | operon partners and missed operons on the pilot (descriptive; C19); not evidence for the simulation claim C14 | pending (M10-04 #296) |
| Fig. 6 | V13 `novel_orfs` | the A2 and A3 freeze runs | novel ORFs by class and call status (C16, C21) | pending (M10-04 #296, M10-10) |
| Supp. Fig. S1 | V16 `compare` | `giap compare` of A1, A2, A3 (`$GIAP_ROOT/runs/pilot/pass1/cmp_a1_a2` and the freeze successor) | C20 | pending (M10-04 #296) |
| Table 1 | the claim map (§4), reduced to the rows `supported` at submission | this file | all | pending (M11-13) |
| Table 2 | runtimes per arm and per step | `benchmarks/baseline.json` and the pilot report | C8, C9, C22 | pending (M10-12) |

## 4. Claim → evidence map

Columns: the claim in one line; the number (`TBD` until evidence exists); its error estimate
(R4: a q-value, FDR, FWER, CI or FDP, or `unestimated`; `n/a` only for a claim that is not
statistical); the evidence, by repository path (`›` names a key or a section); the status.
`pending` rows name the issue that produces the evidence. Rows C1–C10 are the §13.2 bullets
(GIAP_SPEC.md), rows C19–C23 the §15.6 pilot bullets, rows C25–C28 the §15.7 hypotheses fixed
before running (reported where the pilot allows; none may be presented post hoc).

| ID | Source | Claim | Number | Error estimate | Evidence | Status |
|---|---|---|---|---|---|---|
| C1 | §13.2 b1 | ≥ 99.5 % of A2 peptides map exactly (TE1) | TBD | Wilson 95 % CI of the share | `validation/release/1.0.0/readiness.json` › R1; `validation/pilot/metrics/metrics.json` › te1 | pending (M10-04 #296, M10-15b, M11-13) |
| C2 | §13.2 b1 | the decoy-estimated mapping FDR is ≤ 1 % at defaults (A2) | TBD | placement FDR (decoy estimate) | `validation/release/1.0.0/readiness.json` › R2; `docs/methods/placement-fdr.md` | pending (M10-04 #296, M11-13) |
| C3 | §13.2 b2 | every statistical method passes its calibration suite (§8.4), **or ships flagged uncalibrated and is listed** (pre-stated: ORA `ora_uncalibrated`, determinants `null_uncalibrated`, the D-1b-step and M7 G1/A6 strict xfails, the island ρ = 0.9 breaches) | TBD | FDR and FWER at level, four one-sided criteria; flagged methods unestimated | `validation/release/1.0.0/readiness.json` › R3; `docs/methods/calibration.md` | pending (M11-13: `calibration_full` on the release candidate) |
| C4 | §13.2 b2 | every attack null (§8.5) stays at level, **or the method ships flagged uncalibrated and is listed** (as C3) | TBD | rejection rate at α (FWER/FDR per suite) | `validation/release/1.0.0/readiness.json` › R4; `docs/methods/calibration.md`; `tests/calibration/attacks/` | pending (M11-13) |
| C5 | §13.2 b3 | every true-example assertion the data allow (§11.6) passes | TBD | n/a (pass or fail per TE; shares carry Wilson CIs) | `validation/release/1.0.0/readiness.json` › R5; `validation/pilot/metrics/metrics.json` › te_real | pending (M10-06 real data, M10-15b; Q-K5 strictness) |
| C6 | §13.2 b3 | TE2 in-silico proxy: no `high` ALT call off `start_codon`/`stop_codon` sites. The proxy run predates M7-16a/b (stop-codon read-through), so `high_stop_codon` is 0 by construction; **re-run the proxy before submission** | passed (0 violations) | n/a (a count against zero) | `validation/proxies/2026-10-23_te2.json` | supported (pre-M7-16; re-run owed) |
| C7 | §13.2 b4 | the same seed gives byte-identical tables (TE7) | TBD | n/a (deterministic) | `validation/release/1.0.0/readiness.json` › R6; `validation/pilot/metrics/metrics.json` › te7; `validation/pilot/te7.py` | pending (M10-06, M11-13) |
| C8 | §13.2 b5 | a standard run on 16 threads and 64 GB takes about 35 min (16–55) | TBD | unestimated (median of 3, range reported) | `benchmarks/baseline.json` › standard_run.wheat.k12_standard; `validation/release/1.0.0/readiness.json` › R7 | pending (M10-12; Q-K6) |
| C9 | §13.2 b5 | core modules alone take about 25 min | TBD | unestimated (median of 3) | `benchmarks/baseline.json` › standard_run.wheat.k12_core; `validation/release/1.0.0/readiness.json` › R8 | pending (M10-12) |
| C10 | §13.2 b6, b7 | the pilot passes (§15.6); co-abundance and K-12 operon AUROC move to the full validation | TBD; the AUROC claim is not made in v1 | n/a | `validation/release/1.0.0/readiness.json` › R9, R10; `validation/pilot/REPORT.md` | pending (M10-15, Khoa's review; R10 is deferred by the spec) |
| C11 | gate M3 | placement FDR is calibrated on simulations at list level (mean FDP ≤ q + 2 SE); the 100-kb genome-level lists are smoke checks (no false peptides), stratum-0 evidence comes only from the stressed arm and the ADR-0019 estimate checks, and m ≥ 1 rests on the list-level suites | see the note | placement FDR ≤ q (simulation FDP, mean ± SE) | `docs/methods/placement-fdr.md` › Calibration (M3-10), › Genome-level calibration (M3-21); `docs/adr/ADR-0010-placement-fdr-weighted-decoys.md`; `docs/adr/ADR-0019-stratum0-pi-hat-exact-share.md` | supported |
| C12 | gate M4 | gene-level FDR (picked, MAYU) is calibrated; circular-shift region tests hold level. In giapsim the gene-FDR levels 0.05 and 0.10 are weak or vacuous (peptides capped at q ≤ 0.01; only PEP-only input informs 0.10); FDP is bimodal at a 50 % NA-q target share (backlog) | see the notes | gene-level FDR ≤ q; type I error of shift p-values | `docs/methods/detection-tiers.md` › Calibration evidence, › Sign-off (M4-24); `docs/methods/mayu-gene-fdr.md`; `docs/methods/region-tests.md` › Calibration evidence; `docs/adr/ADR-0009-picked-fdr-paired-share.md` | supported |
| C13 | gate M5 | island FWER ≤ 5 % under simulated nulls, within a stated dependence regime (outside it, ρ = 0.9 with short ℓ′, the FWER is breached) | see the note | FWER at α = 0.05, Wilson 95 % CI per cell, R = 500 | `docs/methods/islands.md` › Validity domain (ADR-0021 Amendment 1); `docs/adr/ADR-0021-island-oriter-block-length.md` | supported |
| C14a | gate M6 | missed-operon decoy FDR ≤ q in every arm that makes calls (convergent-only decoys, not shipped); the shipped configuration (convergent + sandwich decoys) makes **no call** in the gate or attacks A, A′, B and M-1 (recorded vacuous pass, near-powerless by construction); M-1 co-imputation is untested | no calls (shipped); see the note for the convergent-only arm | decoy FDR ≤ q (missed operons) | `docs/methods/operons.md` › The M6 gate, simulation half (M6-09); `PROGRESS.md` › Answered by Khoa (2026-10-22, M6-09) | supported (vacuous for the shipped decoys) |
| C14b | gate M6 | operon discordance is a BH family over per-pair t p-values that holds level only under the narrowed claim (smooth, monotone σ(abundance)); D-1b-step stays a confirmed strict-xfail failure outside the claim | see the note | BH q and pooled per-pair level (discordance) | `docs/adr/ADR-0055-operon-discordance-abundance-noise.md` › Amendment 1; `docs/methods/operons.md` › The M6 gate, simulation half (M6-09) | supported (narrowed claim only) |
| C15 | gate M6 | UniOP AUROC on K-12 against RegulonDB pairs: a recorded number for K-12 only, not a transfer claim; possibly in-sample (UniOP's training data unverified); RegulonDB label noise biases it toward the base rate | 0.911 [0.898, 0.923] | directon-bootstrap 95 % percentile CI | `docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md` | supported (K-12 only; in-sample status unverified; a recorded number, not a transfer claim) |
| C16 | gate M7 | variant and novel-ORF F1 ≥ 0.9 on simulations; TE2 proxy passes | TBD | F1 with replicate SE; class FDR | `docs/adr/ADR-0049-m7-variant-precision.md`; `docs/methods/variants.md`; `docs/methods/novel-orfs.md`; `docs/methods/denovo-flag-m717-gate-verdict.md` | pending (M7-13c, M7-16; gate M7 not recorded as passed) |
| C17 | gate M8 | ORA Wallenius tails match the exact reference (relative error ≤ 1e-6 on 300 golden tails; goseq is a second reference only, ADR-0044); GSEA runs R fgsea 1.36.2 itself (ADR-0038, so the fgsea match holds by construction and is not validation evidence) and its single-set null p are uniform under TU exchangeability (restricted claim, ADR-0025, ADR-0026); ORA p are shipped flagged `ora_uncalibrated` (ADR-0035) | see the notes | KS uniformity p (GSEA); golden relative tolerance (ORA); ORA level unestimated (flagged) | `docs/methods/function-gsea.md` › Evidence; `docs/methods/function.md` › Pass-1 limit (decision 13; ADR-0035 decision 3); `docs/adr/ADR-0038-gsea-runs-r-fgsea.md`; `docs/adr/ADR-0044-wallenius-golden-reference.md`; `docs/adr/ADR-0035-ora-tu-detection-recalibration.md`; `docs/adr/ADR-0025-gsea-clustered-sets.md`; `docs/adr/ADR-0026-gsea-conditional-on-detection.md` | supported (GSEA level; ORA tails only) |
| C18 | gate M8 | the determinants permutation p-value is shipped flagged `null_uncalibrated` | n/a | unestimated (fails the §8.5 attacks; flagged) | `docs/adr/ADR-0034-determinants-null-uncalibrated.md`; `docs/methods/determinants.md` | supported |
| C19 | §15.6 b1 | GIAP runs end to end on A1, A2 and A3 | TBD | n/a | `validation/pilot/REPORT.md`; `validation/pilot/metrics/metrics.json` | pending (M10-04 #296, M10-10) |
| C20 | §15.6 b3 | A1 and A2 detections and abundances are compared (`giap compare`) | TBD (Spearman ρ on shared genes; Jaccard of islands, operon flags, categories) | TU-block bootstrap 95 % percentile CI of ρ; Jaccard unestimated | `validation/pilot/metrics/metrics.json` › m4, m5; `docs/methods/compare.md` | pending (M10-04 #296, M10-15b) |
| C21 | §15.6 b1 (A3) | the A3 feasibility check (M10-08 rule, run 1 = WT_1) passes: Cascadia 0.0.7 peptides match the Nissle genome | share 0.636 (549 of 863), verdict pass against min share 0.20 | Wilson 95 % CI 0.604–0.668, optimistic because peptides of one protein co-vary (harmless at this margin); shuffled-decoy share 0.006. The contaminant share is also 5 of 863: checked in the detail file, the two sets are disjoint (a coincidence) | `validation/pilot/metrics/a3_feasibility.json`; `validation/pilot/a3_feasibility.py`; `docs/adr/ADR-0054-a3-external-artefacts.md` | supported (one run only; not a GIAP A3 run) |
| C22 | §15.6 b4 | measured runtimes replace the estimates | TBD | unestimated (median of 3) | `benchmarks/baseline.json` › standard_run; `PROGRESS.md` › Measured numbers | pending (M10-12) |
| C23 | §15.6 b2, b5 | TE1 passes; TE2–TE5 are reported on real data, failures explained | TBD | Wilson 95 % CI per share | `validation/pilot/metrics/metrics.json` › te1, te_real; `validation/pilot/REPORT.md` | pending (M10-06, M10-15b) |
| C24 | §15.5 m8 | A1's free entrapment: `fdp_lower` (a realised lower bound on the FDP) and `ucb95` (one-sided 95 % upper bound, invalid under clustered errors). The estimate errs high except for two unquantified low-direction sources: DIA-NN q refinement without `--no-refine-q`, and false matches attracted to targets | TBD | FDP bound (`fdp_lower`); one-sided 95 % UCB (`ucb95`) | `docs/methods/free-entrapment.md` › Uncertainty, › Bias; `validation/pilot/metrics/metrics.json` | pending (M10-04 #296, M10-07 on real data) |
| C25 | §15.7 b1 (pre-registered) | A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes | TBD | Wilson 95 % CI of the A2-only share; paired test | `validation/pilot/metrics/metrics.json` › h1; `validation/pilot/REPORT.md` | pending (M10-04 #296, M10-15b) |
| C26 | §15.7 b2 (pre-registered) | A1's genome-inconsistent FDP is at least its entrapment FDP | TBD; needs FDRBench | unestimated (full validation) | `docs/methods/free-entrapment.md` | pending (full validation, post-1.0; spec §16) |
| C27 | §15.7 b3 (pre-registered) | context conclusions (islands, region tests, operon concordance) agree between A1 and A2 on shared genes | TBD | Jaccard unestimated; ρ with TU-block bootstrap CI | `validation/pilot/metrics/metrics.json` › h3; `docs/methods/compare.md` | pending (M10-04 #296, M10-15b) |
| C28 | §15.7 b4 (pre-registered) | the genome-anchored de novo precision agrees with precision measured against A2 identifications | TBD | bootstrap CI of the de novo precision | `docs/methods/denovo-precision.md`; `validation/pilot/REPORT.md` | pending (M10-10) |

## 5. Data, code and licence availability

- [ ] **Proteomics data:** PXD078459 (iProX, Nissle 1917, three wild-type runs), cited by
  accession; never redistributed by GIAP. Only aggregates and file checksums are committed
  (`validation/pilot/inputs/pxd078459_files.tsv`).
- [ ] **Genomes and proteomes:** Nissle 1917 GCA_043228065.1; K-12 MG1655 GCA_000005845.2; UniProt
  K-12 proteome UP000000625 (A1's database). Later work: UTI89 GCA_000013265.1 / PXD053727.
- [ ] **Not redistributed, and why:** RegulonDB (its terms; fetched by
  `validation/fetch_regulondb.py`, never committed); DIA-NN output rows (Khoa's pilot data rule);
  DIA-NN itself (proprietary academic licence; the version used is 2.3.1 Academia, not the spec's
  2.7.0, `docs/adr/ADR-0051-diann-version.md`, proposed); the Cascadia checkpoint and vendor SDKs
  (`docs/adr/ADR-0054-a3-external-artefacts.md`).
- [ ] **Code:** GIAP is MIT-licensed; PyPI and Bioconda packages (TBD, M11-04, M11-10); Zenodo
  DOI (TBD, M11-13); documentation URL (TBD, M11-07, Q-K9); containers by digest
  (`workflow/conf/containers.config`, M11-06).
- [ ] **Third-party licences** (`NOTICE.md`, `docs/adr/ADR-0008-tool-pins.md`), to be stated in the paper:
  - GPL-3.0 tools (Bakta, Promotech, Promoter Calculator, OSTIR, TransTermHP, BacTermFinder) and
    AGPL-3.0 antiSMASH run as **separate programs**; GIAP never imports them
    (`scripts/ci/check_imports.py`). R and fgsea run as a separate process (ADR-0038).
  - **geNOMAD is licensed for non-commercial use only**; a commercial user must drop the module
    or obtain a licence.
  - **KEGG is for academic use**; non-academic users need a KEGG licence or the GO and COG
    collections only.
  - UniOP has no licence file: fetched at a pinned commit, never vendored or pushed as an image.
  - DIA-NN is proprietary (academic licence) and used only for users' own searches.

## 6. Limitations (what version 1 defers, GIAP_SPEC.md §16)

- **Single sample per run**: no co-expression or differential context tests across samples.
- **Full validation not done**: K-12 PXD075568, UTI89 beyond the A3 fallback, FDRBench entrapment,
  and the K-12 co-abundance AUROC (§13.2 b7). Release 1.0 is calibrated on simulations and one
  light pilot.
- **Sealed holdout not used**: the EDL933 holdout (PXD068911, §14.5) stays sealed for the full
  validation; no holdout claim is made.
- **A3 de novo arm** is exploratory; its feasibility check covers WT_1 only (C21).
- Deferred modules: assembly checks from frameshift-spanning peptides, CheckM2 in preflight,
  proteome mass fractions and the V9 treemap, RO-Crate export, H-NS GC test, prophage induction,
  RNA-seq and Ribo-seq inputs, the V4 locus browser, taxonomy of unmapped sequences, homology
  evidence for novel ORFs (`changelog.d/+known-limitations.misc.md`).
- **Out of scope by design (§16 "No")**: the recoding and start catalogue (Sec, prfB +1, dnaX −1,
  non-ATG starts, Met removal), start-codon correction from N-terminal peptides, signal-peptide or
  transmembrane prediction, PTM naming or classification, Plass, eggNOG-mapper and MMseqs2/UniRef
  function annotation, and spectrum searching (GIAP only writes the `makedb` FASTA). The paper
  must not imply any of them.

Statistical caveats a reader must apply (each recorded in a cited note):

- **ORA p-values are uncalibrated**: every ORA row carries `ora_uncalibrated`
  (`docs/methods/function.md` › Pass-1 limit (decision 13; ADR-0035 decision 3)).
- **Gate M7 is not passed**: variant F1 0.838 < 0.9 (`docs/adr/ADR-0049-m7-variant-precision.md`).
- **Operons**: the missed-operon family is near-powerless and vacuous with the shipped decoys;
  discordance level is claimed only under smooth, monotone σ(abundance)
  (`docs/adr/ADR-0055-operon-discordance-abundance-noise.md`).
- **GSEA restricted claim**: no level claim for genomically clustered sets; results are
  conditional on detection; BH spillover from flagged rows is not re-measured on the R engine
  (`docs/adr/ADR-0025-gsea-clustered-sets.md`, `docs/adr/ADR-0026-gsea-conditional-on-detection.md`).
- **Islands outside the validity domain**: the FWER is breached at ρ = 0.9 with short ℓ′, and the
  low direction is a structural zero, so a null result is not evidence of no silenced island
  (`docs/methods/islands.md`).
- **Region tests**: the HL shift is relative to random positions on this replicon; the ori → ter
  dosage gradient is not modelled in giapsim (`docs/methods/region-tests.md` › Limitations
  (statistician)).
- **Genome decoys under-count chance hits of protein-like wrong peptides** (D/T 0.96 at realistic
  dipeptide structure; `docs/methods/placement-fdr.md` › Limitations (documented, not CI gates));
  the real-protein TE8 check is still owed.
- **Determinants**: the permutation p-value ships flagged `null_uncalibrated`; it fails the §8.5
  dependence attacks (C18, `docs/adr/ADR-0034-determinants-null-uncalibrated.md`).
- **K-12 calibration**: the TerminatorNet precision against RegulonDB is a lower bound, and the
  UniOP AUROC is a K-12 number whose in-sample status is unverified (UniOP's training data are not
  verified), not a transfer claim (`docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md`).
- **P(miss)** is a probability, not a p-value, and sits outside the BH family
  (`docs/adr/ADR-0056-pmiss-outside-the-bh-family.md`).

## 7. Open questions for Khoa (batched; defaults in use)

- Q-K11: venue, author list, AI-use statement, figure permissions (defaults in §1).
- Whether this outline may stay in the public repository before submission (default: yes; it holds
  no number beyond committed evidence).
- Whether C21 (one-run A3 feasibility) belongs in the paper, or only GIAP's own A3 run (M10-10).
