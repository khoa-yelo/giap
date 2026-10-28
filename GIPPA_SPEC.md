# GIPPA v1: build specification for a Claude Code agent on FarmShare (Slurm HPC)

**GIPPA** is short for *genome integrated prokaryotic proteomics analysis*.

| | |
|---|---|
| Spec version | 1.1, 2026-10-14. Retargeted from the Claude Code cloud VM to Stanford FarmShare (Slurm, Apptainer); the cloud version is archived outside the repository and MUST NOT be followed |
| Owner | Khoa Hoang (the "user" in this file) |
| GIPPA licence | MIT |
| Design source | "GIPPA: design & development plan" (Claude Docs, revision 267). This file restates every decision needed to build version 1, so the agent never needs the design doc. |
| Status | All 20 design decisions and 18 proposed improvements are settled; the plots are chosen. |

---

## 0. How to use this file

### 0.1 For Khoa: one-time setup (about 15 minutes)

1. **Repository.** The GitHub repository is `git@github.com:khoa-yelo/gippa.git` (SSH; FarmShare's SSH key is already authorised for the `khoa-yelo` account). Its `main` branch holds only an initial `README.md` commit. The working copy is `/home/classes/bios/270/khoa/giap`, which already holds this file. In M0 the agent connects the folder to the remote and opens the first pull request (§13.1).
2. **Session.** On a FarmShare login node, start `tmux new -s gippa` (so a dropped SSH connection does not kill the session), then an interactive allocation and Claude Code:

   ```bash
   srun --partition=interactive --qos=interactive --cpus-per-task=8 --mem=32G --time=1-00:00:00 --pty bash
   cd /home/classes/bios/270/khoa/giap && module load claude-code && claude
   ```

   Paste the kickoff prompt from §0.3. The agent bootstraps its own user-space tools (§3.4); nothing needs root.
3. **GitHub login.** `git push` and `git fetch` already work over SSH. For pull requests and issues, type `! gh auth login` in the session when the agent asks (M0), choosing SSH as the git protocol, and set your git name and email for the repository (§3.10).
4. **Licences.** When the agent asks (M3 at the earliest), accept the DIA-NN 2.7.0 licence and unpack its Linux binary under `$GIPPA_TOOLS_DIR/diann-2.7.0` (§15.2). This is the only tool the agent may not install itself.
5. **Recurring jobs.** At M2 the agent proposes `scrontab` entries for the nightly and weekly FarmShare suites (§12). Review and install them yourself with `scrontab -e`.
6. **Check in.** Whenever you return, read `PROGRESS.md` and the open pull requests. Your checkpoints are listed in §14.8; the agent waits for you only there. When the 24-hour allocation expires, start a new one and use the resume prompt (§0.3); batch jobs keep running in the meantime.

**Everything stays in this folder.** Code is tracked in git; all data, databases, environments, containers, caches, Nextflow work directories and job logs live in the gitignored subfolder `local/` of the repository (`$GIPPA_ROOT`, §3.3). Nothing goes to `$HOME` or `/scratch`.

Everything in this spec runs on FarmShare: the full Bakta database, antiSMASH, the DIA-NN searches, de novo sequencing on the `gpu` partition (L40S) and the 16-thread benchmark. The agent submits these as Slurm batch jobs from `scripts/slurm/` (§3.7, §3.8, §15).

### 0.2 For the agent: rules that override everything else

These rules are numbered so that commits, ADRs and reviews can cite them (for example "R3").

- **R1 Scope.** Build exactly the version 1 scope in §2. Do not build anything listed in §16, even partially, and do not add features the spec does not ask for.
- **R2 Test first.** Every behaviour starts as a failing test. Tests and the code that satisfies them come from different subagents, following the feature loop in §14.3.
- **R3 Never weaken a gate.** Do not edit protected files (`ci/protected.txt` lists them: gates, coverage thresholds, calibration tolerances, fixed seeds), skip tests, mark them `xfail`, or widen a tolerance to get CI green. Such a change needs Khoa's written approval, recorded in an ADR.
- **R4 Error estimates.** Every statistical claim GIPPA writes carries an error estimate (q-value, FDR or FWER) or the literal label `unestimated` (§8).
- **R5 Licences.**
  - GPL, AGPL and unlicensed tools (Bakta, Promotech, Promoter Calculator, OSTIR, TransTermHP, antiSMASH, UniOP) run only as separate processes in their own environments.
  - Never `import` them from GIPPA code, and never copy their source into the repository.
  - geNOMAD is non-commercial: document that in the user guide.
- **R6 Network.**
  - Download pinned artefacts only from the sources in §3.2 and verify their checksums.
  - When a source fails (timeout, HTTP error or checksum mismatch) after one later retry, stop that sub-task. Record it under "Blocked" in `PROGRESS.md` and continue with other work.
  - Never substitute unofficial mirrors or re-hosted copies for a pinned artefact, and never disable TLS verification.
- **R7 Determinism.** All randomness flows from `config.seed`. The same inputs and seed must give byte-identical files under `tables/`, `mapping/` and `annotation/*.tsv`.
- **R8 Continuity.** Update `PROGRESS.md` at the end of every work unit (template in Appendix A.3). A new session must be able to resume from it alone.
- **R9 Escalation.** After 3 failed iterations on one problem or loop, stop, write the question under "Questions for Khoa" in `PROGRESS.md`, tell Khoa in chat, and move to other work.
- **R10 Decide or ask.** Ask Khoa only about the checkpoints in §14.8. Otherwise make the call, record it (an ADR when it changes a default, threshold, interface or dependency) and continue.
- **R11 Data.**
  - All data, databases, environments and caches live in the gitignored `local/` folder (§3.3); nothing under it is ever committed.
  - Never commit data you may not redistribute: RegulonDB files, raw MS data, DIA-NN binaries, Bakta or geNOMAD databases, or UniOP code.
  - No single file may exceed 5 MB, and `tests/data/` stays under 50 MB in total.
  - Compress text fixtures with gzip.
- **R12 Privacy.** Never put secrets, tokens or personal email addresses in code, commits, logs or reports.
- **R13 Spec changes.** Never edit `GIPPA_SPEC.md`. Deviations go into an ADR (`docs/adr/`), which wins over this file once Khoa accepts it.
- **R14 Shared cluster.** FarmShare is shared, and the same account runs Khoa's other projects.
  - Never run tests, tools or builds on a login node; work only inside the session allocation or a batch job.
  - Work that exceeds about 10 minutes or the session allocation goes to `sbatch` (§3.8).
  - Every GIPPA job name starts with `gippa-`. Cancel only jobs GIPPA submitted, by job ID; never `scancel -u`, and never touch, inspect in depth or reprioritise other jobs.
  - GIPPA uses at most 128 CPUs, 4 GPUs and 20 queued or running jobs at once (Nextflow `executor.queueSize = 20`).
  - All large state lives under `$GIPPA_ROOT` (`local/` in the repository folder, gitignored); stay within the §3.9 budget and clean Nextflow work directories after successful runs.

### 0.3 Kickoff prompt (Khoa pastes this into the first session)

```text
You are the lead agent building GIPPA v1. Read GIPPA_SPEC.md completely before acting.
Follow §0.2 at all times. You run on Stanford FarmShare inside a Slurm allocation (§3).
Start with §13, milestone M0: run the bootstrap (§3.4), create AGENTS.md, CLAUDE.md,
PROGRESS.md, .claude/ and ci/ from Appendix A, run the environment probe (§3.6), then
work through the milestones in order using the feature loop (§14.3). Keep pull requests
small. Update PROGRESS.md at the end of every work unit. When you need me, put the
question under "Questions for Khoa" in PROGRESS.md and tell me in chat; otherwise keep going.
```

For later sessions:

```text
Resume GIPPA v1. Read CLAUDE.md, then PROGRESS.md, then the GIPPA_SPEC.md sections it points to.
Check the GIPPA Slurm jobs listed under "In progress" (squeue/sacct). Continue from "Next actions".
```

### 0.4 Contents

| § | Section | Read when |
|---|---|---|
| 1 | What GIPPA is | first session |
| 2 | Version 1 scope and settled decisions | first session; before any design choice |
| 3 | The FarmShare HPC environment: Slurm, Apptainer, storage, GitHub | first session; when something cannot run; before submitting jobs |
| 4 | Repository layout | M0 |
| 5 | Interfaces: CLI, config, input formats, converters | M0–M1, then as needed |
| 6 | Architecture: steps, Nextflow, internal model | M0, M2 |
| 7 | Module specifications (Modules 0–8, utilities) | the milestone that builds each module |
| 8 | Statistical framework and calibration | every statistical change |
| 9 | Output standard and table schemas | every module |
| 10 | Visualisations and report | M9 |
| 11 | Testing strategy, simulator, fixtures | M0, then every feature |
| 12 | CI/CD | M0, M11 |
| 13 | Milestones and gates | every session |
| 14 | Agent workflow, roles, loops, guardrails | first session; every PR |
| 15 | Pilot validation kit (Nissle 1917) | M10 |
| 16 | Deferred and excluded: do not build | before adding anything |
| 17 | Risks and fallbacks | when blocked |
| 18 | References | as needed |
| A | File templates: AGENTS.md, CLAUDE.md, PROGRESS.md, subagents, hooks, ADR, PR | M0 |
| B | Worked numbers for tests | when writing tests |

Keywords: **MUST**, **MUST NOT**, **SHOULD** and **MAY** follow RFC 2119. "Khoa" is the user.

---
## 1. What GIPPA is

GIPPA takes one assembled bacterial genome plus proteomics evidence and answers one question: **what does genome context explain about which proteins are detected, and how abundant they are?**

It runs after the proteomics search and does not search spectra itself. It then:

- places peptides and proteins on the genome, with its own error control;
- builds a single gene table keyed by Bakta locus tags;
- tests genome context: regions, autocorrelation along the chromosome, abundance islands, operons, function, and the sequence determinants of abundance;
- calls REF/ALT residues at covered codons and finds novel ORFs;
- writes standardised tables, plots and a self-contained HTML report.

### 1.1 Inputs

| Input | Required | Format | Notes |
|---|---|---|---|
| `genome.fasta` | yes | nucleotide FASTA | 1 to a few contigs; circular replicons are flagged from headers |
| `peptides.tsv` | recommended | generic peptide table (§5.3) | from any search engine or de novo tool, via `gippa convert` |
| `proteins.faa` | yes, unless peptides are given | protein FASTA | protein hits from a database search, or de novo sequences; ID = first token of the header |
| `abundance.tsv` | optional | TSV: `protein_id` + one value column | protein-level values; used only when peptides carry no intensities |
| `--annotation` | optional | Bakta JSON (preferred), GFF3 or GBFF | reuses the gene model of the search database, so IDs match |
| `--reference-peptides` | optional | generic peptide table | same-strain database results for the same spectra; calibrates de novo precision and ALT error rates |
| `config.yaml` | optional | YAML (§5.2) | every parameter; resolved values are written to the run manifest |

### 1.2 Design rules (apply everywhere)

1. **One gene model.** Bakta locus tags are the key in every table. Every tool that accepts gene calls receives Bakta's.
2. **Every claim carries an error estimate**, or is labelled `unestimated`.
3. **Core modules always run.** The optional modules (`mobile` = geNOMAD, `bgc` = antiSMASH, `function` = Module 7) can be switched off without breaking later steps.
4. **Isolation.** Each external tool runs in its own pinned environment or container, as a separate process.
5. **Traceability.** `run_manifest.json` records tool versions, database checksums, parameters, input SHA-256 hashes, seeds and per-step time and peak memory.
6. **Say "detected", never "expressed".** All output text, column names and plot titles follow this.

### 1.3 Modules

| # | Module | Core? | Output tables |
|---|---|---|---|
| 0 | Preflight | core | `preflight/preflight.json`, `preflight.md` |
| 1 | Annotation and context layers | core; `mobile` and `bgc` optional | `annotation/*` |
| 2 | Mapping proteins and peptides | core | `mapping/*` |
| 3 | Gene table, detection, abundance, context statistics | core | `genes.tsv`, `detection_summary.tsv`, `regions_stats.tsv`, `autocorrelation.tsv`, `islands.tsv` |
| 4 | Operons | core | `pairs.tsv`, `operons.tsv`, `operon_flags.tsv` |
| 5 | Variants (REF/ALT) | core | `pileup.tsv.gz`, `variants.tsv`, `variants.vcf` |
| 6 | Novel ORFs and start sites | core | `novel_orfs.tsv`, `.gff3`, `.faa`, `annotation_patch.gff3` |
| 7 | Function: COG, GO, KEGG, GSEA | optional (`function`) | `category_summary.tsv`, `enrichment_detection.tsv`, `enrichment_gsea.tsv` |
| 8 | Determinants of abundance | core | `determinants_features.tsv`, `determinants_model.tsv`, `determinants_attribution.tsv` |

Utilities: `gippa setup`, `preflight`, `run`, `makedb`, `convert`, `compare`, `report` (§5.1).

---

## 2. Version 1 scope and settled decisions

### 2.1 The 20 decisions (all accepted by Khoa)

| # | Topic | Final decision |
|---|---|---|
| 1 | Name | GIPPA; package and command `gippa`. |
| 2 | Missing abundance (identified, not quantified) | `median` by default, since these proteins were detected. `min` and `exclude` run as sensitivity checks; `zero` and `mean` are also accepted. Genes never identified are *not detected* and are never imputed. |
| 3 | Abundance scale | With peptide intensities, GIPPA computes gene-level iBAQ (primary) and Top3 itself, after mapping, from unique peptides. Protein-level values are used only when no peptide intensities are given: `ibaq` and `top3` count as molar; `maxlfq` and `intensity` count as relative, which disables operon abundance ratios and mass fractions. |
| 4 | Detection rule | Three tiers, chosen from what the input carries:<br>• decoy rows → picked gene-level target–decoy FDR;<br>• q-values or PEPs only → MAYU-style gene FDR;<br>• sequences only → genes with ≥ 2 distinct peptides, single-peptide genes `unestimated`.<br>The ≥ 2-peptide tier is always reported too. |
| 5 | Peptide mapper | GIPPA's own exact six-frame index for 7–40 aa. miniprot handles proteins > 40 aa only, with tblastn as fallback. |
| 6 | Tolerant matching | For peptides ≥ 10 aa: 1 difference at 10–11 aa, 2 at ≥ 12 aa, accepted at placement q ≤ 0.01. Proteins are assigned at ≥ 80 % identity over ≥ 80 % of their length. `--genome-relation related` is the default. |
| 7 | Operon threshold | UniOP p ≥ 0.5 to start, recalibrated on RegulonDB to about 90 % precision. An intergenic-distance model is the fallback. |
| 8 | Promoters | Promotech RF-HOT on targeted 40-nt windows whose 3′ end lies −200 to +20 bp around every start codon, novel ORFs included. Promoter Calculator gives σ70 strength. Genome-wide Promotech is opt-in. |
| 9 | Terminators | TerminatorNet by default, with its cut-off set on RegulonDB. TransTermHP is the fallback and cross-check. BacTermFinder is opt-in for Rho-dependent terminators. |
| 10 | Function annotation | Bakta full-database cross-references only (COG, GO, EC, KEGG). No eggNOG-mapper. Every enrichment reports its annotation coverage. |
| 11 | Workflow engine | Nextflow DSL2 behind the `gippa` Click CLI, one container or conda environment per tool, nf-core conventions where they fit. |
| 12 | Preflight contiguity | Warn on N50 < 100 kb or largest contig < 1 Mb. The median contig length is recorded for information only. |
| 13 | Significance levels | 0.01 for identification, placement, novel ORFs and ALT calls; 0.05 for context and enrichment; 0.1 for operon flags. |
| 14 | Permutations | 1,000 by default; 10,000 with `--publication`. Besag–Clifford, Wallenius and fgsea multilevel avoid the permutation floor where many tests share one correction. |
| 15 | Variant calls | REF, ALT, `mixed` or `no_call` per covered codon, from distinct-peptide counts. ALT is `high` with ≥ 2 distinct peptides and `candidate` with 1. Mass look-alikes are flagged, never reclassified. No PTM naming. |
| 16 | Coordinates | 1-based closed in TSV, GFF3 and VCF; 0-based half-open in BED and proBed. Code works 0-based internally and converts only in the file writers. |
| 17 | Licence | MIT. GPL and AGPL tools are only ever called as separate programs. |
| 18 | UniOP | A standard dependency, fetched at a pinned commit when environments are built; never vendored (it has no licence file). The distance model is the fallback. |
| 19 | Validation | Pilot first: 3 Nissle 1917 wild-type runs (PXD078459), searched against K-12 UniProt (A1), against Nissle's Bakta proteins (A2), and de novo (A3). Cut-offs are calibrated on the K-12 MG1655 genome against RegulonDB. The full validation comes later. |
| 20 | Data stack | Python ≥ 3.11 with pandas + pyarrow, numba for hot loops, pydantic for config. Revisit polars only if peptide tables reach tens of millions of rows. |

### 2.2 Improvements

| # | Improvement | Version 1? |
|---|---|---|
| 1 | Generic peptide table from any engine or de novo tool | **Yes** |
| 2 | Reuse an existing annotation (`--annotation`, or Bakta `--regions`) | **Yes** |
| 3 | `gippa makedb`: optional strain search database for the user's own search; no FDR is built on it | **Yes** |
| 5 | Detectability model in every detection statistic | **Yes** |
| 9 | Contaminant check of the unmapped bucket (replaces Plass); taxonomy lookup later | **Yes** |
| 12 | `gippa compare` run-vs-run subcommand | **Yes** |
| 13 | Interactive self-contained HTML report (igv locus browser later) | **Yes** |
| 4 | Several samples in `abundance.tsv` (true co-expression) | Later |
| 7 | Assembly check from frameshift-spanning peptides | Later |
| 10 | CheckM2 genome QC in preflight | Later |
| 11 | Proteome mass fractions by category | Later |
| 14 | RO-Crate provenance export | Later |
| 15 | Low-abundance islands tested against GC (H-NS) | Later |
| 16 | Prophage induction signature | Later |
| 18 | RNA-seq and Ribo-seq inputs | Later |
| 6 | Recoding and start catalogue (Sec, prfB, dnaX, non-ATG starts) | No |
| 8 | Start-codon correction from N-terminal peptides | No |
| 17 | Signal-peptide and transmembrane predictor | No |

### 2.3 Plots in version 1

V1–V3, V5–V8 and V10–V18 (16 plots) are in version 1; V4 (igv browser) and V9 (mass-fraction treemap) come later (§10).

### 2.4 Implementation choices made by this spec

The design doc left the choices below open; this spec fixes them. The agent MAY revise any of them only through an ADR that Khoa accepts.

| # | Choice | Why |
|---|---|---|
| G1 | New subcommands `gippa setup` (download and verify databases and models) and `gippa convert` (engine converters) | The doc needs both but did not name them. |
| G2 | Optional peptide-table columns `spectrum_id` and `search_db`, plus `--reference-peptides` | De novo precision (c) and ALT error rates need same-spectrum database results. |
| G3 | Placement q-values are computed within each mismatch stratum (0, 1, 2), with hits ranked by expected chance hits, which pools length classes. The per-length × mismatch table is still reported. | With the +1 correction, a small stratum (T < 100) could never pass 0.01 on its own. Keeping mismatch strata separate preserves the subgroup control the doc asked for. |
| G4 | Tolerant matcher = q-gram seeds plus banded verification that charges each allowed difference type as one difference | Covers the doc's "canonicalised index", including 1↔2-residue exchanges that change length. |
| G5 | Detectability GAM = cubic B-spline basis per covariate plus L2-penalised logistic regression (scikit-learn) | Avoids `pygam` fragility; it is the same model class. |
| G6 | Transmembrane segments are counted from Kyte–Doolittle hydropathy (window 19, threshold 1.6) | Improvement 17 (a predictor) was declined, but the covariate is still needed. |
| G7 | MAYU-style gene FDR has a closed form per protein-size bin (§7.3.2) | Simple, testable, faithful to Reiter 2009. |
| G8 | ALT calls get the informational flags `start_codon` and `stop_codon` | No recoding catalogue (improvement 6 = No), but the true-example tests need to recognise non-ATG starts and Sec sites. |
| G9 | `tables/sensitivity.tsv`, with Jaccard stability under `min` and `exclude` | The doc requires sensitivity reruns but listed no file. |
| G10 | A developer-only local executor plus precomputed-layer injection (`annotation.precomputed`) | Fast end-to-end tests on simulated genomes without external tools. |
| G11 | Promoter Calculator runs once per transcription-unit (TU) upstream region, in parallel | Its CLI joins all FASTA lines into one sequence. |
| G12 | TerminatorNet always runs with `-s` (both strands); its scikit-learn version is pinned to the one that pickled `model.pkl` | The default searches the forward strand only. |
| G13 | Novel ORFs are ranked by a lexicographic key, with no weighted score | Interpretable, with no arbitrary weights. |
| G14 | VCF ALT = the minimal-change codon for the ALT residue; the residue itself goes in INFO | VCF needs nucleotides; peptides give amino acids. |
| G15 | fgsea multilevel = a Python port of fgsea (MIT), checked against R fgsea | The doc asks for a Python implementation cross-checked with fgsea. |
| G16 | For the → ← → decoy, the "intergenic distance" is the sum of its two gaps | Makes it comparable with NO-pair distances when reweighting. |
| G17 | Full-database Bakta outputs for test genomes are recorded by the agent on FarmShare (`scripts/slurm/record_bakta_full.sbatch`) and committed as fixtures | GitHub runners cannot hold the 84 GB database, so CI tests read the recorded outputs. |
| G18 | Development uses Python 3.12; CI covers 3.11 and 3.13 | numpy ≥ 2.5 and scipy ≥ 1.18 need Python ≥ 3.12; 3.11 stays supported through the lock file. |
| G19 | Development and heavy runs target FarmShare: Nextflow's Slurm executor with Apptainer (`-profile farmshare,apptainer`), conda as fallback; nightly and weekly suites run as Slurm jobs, GitHub Actions keeps the per-push checks | No root and no Docker on FarmShare; Apptainer runs the same OCI images by digest, and the cluster has the CPU, memory, disk and GPUs the cloud VM lacked. |

### 2.5 Out of scope for version 1

Everything marked Later or No in §2.2, V4 and V9, and the full validation (K-12 PXD075568, UTI89 PXD053727, the EDL933 holdout, FDRBench entrapment) are out of scope. So are multi-sample statistics, PTM classification, taxonomy lookup of unmapped sequences, `possible_assembly_error`, prophage induction flags and mass fractions. §16 is the binding list.

---
## 3. The FarmShare HPC environment

GIPPA is developed on **Stanford FarmShare**, a shared Slurm cluster. The agent is a Claude Code session running inside a Slurm interactive allocation; heavy work goes to Slurm batch jobs.

### 3.1 What FarmShare provides

Facts below were probed on 2026-10-14; the probe in §3.6 re-checks them at every session start.

| Resource | Value | Consequence for GIPPA |
|---|---|---|
| OS and rights | Ubuntu 24.04.5, x86_64, Slurm. A normal user account: no root, no `sudo`, no `apt` | Install everything in user space: micromamba, uv, Apptainer images |
| Session host | Claude Code runs inside an interactive Slurm allocation (current session: node `iron-01`, 8 CPU, 32 GB, 24 h limit). QOS `interactive` allows ≤ 16 CPU, 64 GB and 4 jobs | Lint, unit tests and small simulations run in the allocation. Anything longer than ~10 min or larger than the allocation goes to `sbatch` (§3.8). When the allocation ends, the session ends; batch jobs keep running |
| Partition `normal` (default, 2 days) | `wheat-01..08`: 80 CPU, 191 GB · `barley-01..04`: 128 CPU, 384 GB · `oat-01..06`: 64 CPU, 256 GB · `rye-01..02`: 512 CPU, 1.5 TB | All CPU work. `DefMemPerCPU` 2,400 MB and `MaxMemPerCPU` 4,000 MB, so always pass `--mem`; Slurm raises the CPU count when memory per CPU exceeds 4 GB |
| Partition `bigmem` | `rye-01..02`, up to 1.5 TB | not needed for bacterial genomes |
| Partition `gpu` | `oat-01..06`, 4 × NVIDIA L40S (48 GB) each; QOS `gpu` ≤ 4 GPUs | de novo sequencing (Cascadia, Casanovo; §15) |
| QOS limits per user | `normal`: ≤ 512 CPU, ≤ 128 jobs · `long`: 7 days, ≤ 128 CPU, 4 jobs · `dev`: 8 h, 8 CPU, 32 GB, 1 job · `interactive`: 16 CPU, 64 GB, 4 jobs | The same account also runs Khoa's other projects, so GIPPA stays far below these limits (R14) |
| Pre-installed | Python 3.12.3, OpenJDK 21, git, curl; Lmod modules include `apptainer/1.4.4`, `r/4.5.2`, `python/3.12–3.14`, `blast-plus/2.16.0`, `pandoc`, `texlive`, `claude-code`; micromamba 2.3.3 in `~/.local/bin` | Not installed: Docker, uv, gh, Nextflow, nf-test. The bootstrap (§3.4) installs them in user space |
| Containers | Apptainer 1.4.4; no Docker daemon | Images are referenced by digest and run as `docker://…@sha256:…` through Apptainer. The `docker` profile is kept for users elsewhere but is not tested here |
| Network | login, interactive and compute nodes reach the internet directly; there is no proxy (zenodo.org and github.com checked from a compute node) | Downloads may run inside batch jobs. Download only from the sources in §3.2 |
| GitHub | ordinary internet access; REST and GraphQL both work | `gh pr …`, `gh issue …` and `gh run …` work once `gh` is authenticated (§3.10) |
| Claude Code commands | foreground default 2 min, up to the limits set in §3.3 | Never block a foreground command on a long job; submit and poll (§3.8) |

**Storage.**

| Path | Kind | Size and notes | GIPPA use |
|---|---|---|---|
| `/home/classes/bios/270/khoa/giap` (also `/farmshare/home/classes/bios/270/khoa/giap`) | NFS class space | shared pool, 208 TB, about 199 TB free | the git repository **and**, in its gitignored subfolder `local/` (`$GIPPA_ROOT`), every database, dataset, environment, container image, cache, Nextflow work directory and job log. Everything under `local/` must be re-creatable by `gippa setup` and the scripts |
| `$HOME` (`/home/users/khoang99`) | NFS home | small quota, already about 32 GB used | nothing; every cache is redirected into `local/` (§3.3) |
| `/scratch/users/khoang99` | NFS scratch | shared pool | not used by GIPPA |
| `/tmp` on each node | node-local SSD (179 GB on `iron-01`, about 1.1 TB on compute nodes) | private to the node and the job | `TMPDIR` for tools inside a job; never for outputs or job logs |

### 3.2 Download sources (pinned)

FarmShare has no allowlist. The agent still downloads pinned artefacts only from these official sources, checks every checksum, and records each in `db_manifest.json` (R6):

| Source | Needed for | First needed |
|---|---|---|
| PyPI, `conda.anaconda.org` (conda-forge, bioconda; never the `defaults` channel), GitHub releases, `astral.sh` (uv installer) | Python and tool environments, Nextflow, gh | M0 |
| Docker Hub, `ghcr.io`, optionally `quay.io` (BioContainers) | container images, pulled by digest into Apptainer | M2 |
| `zenodo.org` | Bakta DB v6.0 (record 14916843), geNOMAD DB v1.9 (record 14886553) | M2 |
| `ftp.ncbi.nlm.nih.gov`, `api.ncbi.nlm.nih.gov` | AMRFinderPlus DB (pulled by `bakta_db`), reference genomes, COG category definitions | M2 |
| `www.cs.mun.ca` | Promotech RF models (`/~lourdes/public/PromoTech_models/`) | M2 |
| `dl.secondarymetabolites.org` | antiSMASH databases, when not already in the image | M2 |
| `rest.uniprot.org` | K-12 proteome UP000000625 (in-silico tests, pilot A1) | M3 |
| `purl.obolibrary.org`, `*.geneontology.org` | `go-basic.obo` | M8 |
| `rest.kegg.jp` | KO → pathway links (academic use) | M8 |
| `regulondb.ccg.unam.mx` | RegulonDB v12 transcription units and terminators (calibration, M6/M10) | M6 |
| iProX / ProteomeXchange (`proteomecentral.proteomexchange.org`, `www.iprox.cn` and the FTP host the dataset page names) | raw files of PXD078459 (pilot) | M10 |

When a source fails (timeout, HTTP 4xx/5xx, checksum mismatch), retry once later; if it still fails, the sub-task is blocked, not the project. The agent records the block (R6), uses recorded fixtures (§11.4), and continues.

### 3.3 Environment file (`scripts/farmshare_env.sh`, committed)

Every GIPPA shell sources this file: the SessionStart hook (§3.5), every `sbatch` script (§3.8) and the bootstrap (§3.4). It keeps all state inside the repository folder, under the gitignored `local/`, and separates GIPPA from Khoa's other projects (for example, the existing micromamba root under `/scratch/users/khoang99/bios270` is never touched).

```bash
# GIPPA environment on FarmShare. Source it; never execute it.
export GIPPA_REPO="${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}"
export GIPPA_ROOT="$GIPPA_REPO/local"                   # gitignored; all data and environments live here
export GIPPA_DB_DIR="$GIPPA_ROOT/db"
export GIPPA_DATA_DIR="$GIPPA_ROOT/data"
export GIPPA_TOOLS_DIR="$GIPPA_ROOT/tools"            # licence-gated binaries (DIA-NN), installed by Khoa
export GIPPA_LOG_DIR="$GIPPA_ROOT/logs"
export MAMBA_ROOT_PREFIX="$GIPPA_ROOT/micromamba"    # GIPPA's own root, separate from other projects
export CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"
export UV_CACHE_DIR="$GIPPA_ROOT/cache/uv"
export UV_PYTHON_INSTALL_DIR="$GIPPA_ROOT/uv-python"
export UV_PROJECT_ENVIRONMENT="$GIPPA_ROOT/venv"
export UV_LINK_MODE=copy
export PIP_CACHE_DIR="$GIPPA_ROOT/cache/pip"
export APPTAINER_CACHEDIR="$GIPPA_ROOT/cache/apptainer"
export NXF_HOME="$GIPPA_ROOT/nextflow"
export NXF_APPTAINER_CACHEDIR="$GIPPA_ROOT/containers"
export NXF_CONDA_CACHEDIR="$MAMBA_ROOT_PREFIX/nxf-envs"
export NXF_WORK="$GIPPA_ROOT/work"
export NXF_OPTS="-Xms512m -Xmx4g"
export PATH="$HOME/.local/bin:$MAMBA_ROOT_PREFIX/envs/gippa-wf/bin:$PATH"
module load apptainer/1.4.4 >/dev/null 2>&1 || true
```

Claude Code's own limits go in the `env` block of `.claude/settings.json` (Appendix A.4): `BASH_DEFAULT_TIMEOUT_MS=600000` and `BASH_MAX_TIMEOUT_MS=3600000`.

### 3.4 One-time bootstrap (`scripts/bootstrap_farmshare.sh`, run by the agent in M0)

User space only, idempotent, safe to re-run. It runs inside the session allocation (a few minutes); if it would take longer, submit it with `sbatch`.

```bash
#!/bin/bash
# GIPPA one-time bootstrap on FarmShare. User space only; idempotent.
set -euo pipefail
source "$(dirname "$0")/farmshare_env.sh"
mkdir -p "$GIPPA_ROOT"/{db,data,tools,logs,cache,containers,work,worktrees,ci} "$MAMBA_ROOT_PREFIX"
grep -qx '/local/' "$GIPPA_REPO/.gitignore" 2>/dev/null || echo '/local/' >> "$GIPPA_REPO/.gitignore"
# uv: the official standalone installer, into ~/.local/bin; record the version in PROGRESS.md
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
# Nextflow, nf-test, Java and gh from conda-forge and bioconda only, never the 'defaults' channel
micromamba env list | grep -q '^ *gippa-wf ' || \
  micromamba create -y -n gippa-wf -c conda-forge -c bioconda --override-channels \
    nextflow=26.10.0 nf-test=0.9.5 openjdk=21 gh
micromamba clean -a -y >/dev/null
nextflow -version && nf-test version && gh --version && uv --version
```

If `nextflow=26.10.0` is unavailable, pin the newest stable release at or above 25.04 and record it in an ADR. Heavy downloads (Bakta databases, container images, tool environments) do not belong here. The agent fetches them with `gippa setup` and `scripts/build_envs.sh`, as batch jobs (§3.8).

### 3.5 SessionStart hook (committed in the repository)

`.claude/settings.json` registers `scripts/session_start.sh` for `startup|resume` (Appendix A.4):

```bash
#!/bin/bash
# Runs at every Claude Code session start or resume on FarmShare.
cd "${CLAUDE_PROJECT_DIR:-.}"
source scripts/farmshare_env.sh
if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "WARNING: not inside a Slurm allocation. Do not run tests or tools on a login node (R14)."
fi
if [ -f pyproject.toml ] && command -v uv >/dev/null; then
  uv sync --frozen --all-extras >"$GIPPA_LOG_DIR/uv-sync.log" 2>&1 || echo "uv sync failed (see $GIPPA_LOG_DIR/uv-sync.log)"
fi
if [ -f scripts/probe_env.py ]; then uv run --no-sync python scripts/probe_env.py --quiet || true; fi
# Remind the lead of GIPPA jobs still queued or running
squeue --me --noheader -o '%i %j %T %M %l' 2>/dev/null | grep ' gippa-' || true
exit 0
```

### 3.6 Capability probe and test markers

`scripts/probe_env.py` writes `.gippa/capabilities.json` (gitignored) and prints a one-screen summary. It records:

- host, CPU count and RAM of the current allocation (`SLURM_CPUS_ON_NODE`, `SLURM_MEM_PER_CPU`) and the time left in it (`squeue -h -j $SLURM_JOB_ID -o %L`);
- free space under `$GIPPA_ROOT` and the size of each subdirectory (warn above the §3.9 budget), and free space in `/tmp`;
- Python, Java, Nextflow, nf-test, micromamba, uv, gh and Apptainer versions;
- whether `sbatch` works (`sbatch --test-only`), the GIPPA jobs in the queue, and the QOS limits (`sacctmgr show qos`);
- whether the `gpu` partition is visible;
- reachability of each §3.2 source (HEAD or small GET, 5 s timeout);
- presence and checksums of databases under `$GIPPA_DB_DIR`, and of DIA-NN under `$GIPPA_TOOLS_DIR`.

Tests that need a capability carry a pytest marker. Pytest skips them with an explicit reason when the capability is missing, and the skip list is printed at the end of every run.

| Marker | Requires |
|---|---|
| `needs_apptainer` | `apptainer` on PATH and able to run a pulled image |
| `needs_docker` | Docker daemon running (never satisfied on FarmShare; GitHub runners only) |
| `needs_conda` | micromamba and the relevant tool environment |
| `needs_nextflow` | `nextflow` and `nf-test` on PATH |
| `needs_slurm` | `sbatch --test-only` succeeds (tests that submit a real job; nightly only) |
| `needs_gpu` | running inside a job with a GPU (`CUDA_VISIBLE_DEVICES` set) |
| `needs_bakta_db` | Bakta light or full DB under `$GIPPA_DB_DIR/bakta` (`needs_bakta_full` for the full one) |
| `needs_genomad_db` | geNOMAD DB v1.9 |
| `needs_promotech_models` | `RF-HOT.model` present |
| `needs_diann` | DIA-NN 2.7.0 under `$GIPPA_TOOLS_DIR` |
| `needs_net` | the named host is reachable |
| `slow` | > 60 s; excluded from the per-commit run |
| `calibration`, `calibration_full` | statistical calibration (fast subset or full) |
| `e2e` | end-to-end runs on real genomes |

In CI and in the FarmShare nightly and weekly jobs, every marker a job declares as required MUST be satisfied; a missing capability fails the job instead of skipping it (`--require-capabilities` option of the `tests/conftest.py` plugin). Silent skips can therefore never hide failures.

### 3.7 What runs where

| Work | Session allocation (8 CPU, 32 GB) | Slurm batch jobs on FarmShare | GitHub Actions | Khoa |
|---|---|---|---|---|
| Unit, property, golden, contract, snapshot tests | yes | — | yes (every push) | — |
| Fast calibration subset | yes | — | yes (every PR) | — |
| Full calibration suite | no | nightly, as a job array (§12) | — | — |
| Tool integration on the 100-kb synthetic genome: Bakta (light and full DB), UniOP, TerminatorNet, OSTIR, Promoter Calculator, Promotech targeted, miniprot, tblastn, TransTermHP, codon-bias, geNOMAD, antiSMASH | short single-tool checks only | nightly, `-profile farmshare,apptainer` | stub mode only | — |
| Bakta with the full DB on K-12 MG1655, Nissle 1917 and the 100-kb genome (recorded fixtures, G17) | no | `scripts/slurm/record_bakta_full.sbatch` (16 CPU, 64 GB) | — | — |
| End-to-end runs on K-12 and Nissle (live full-DB Bakta) | no | weekly | — | — |
| 16-thread / 64-GB runtime benchmark | no | `scripts/slurm/benchmark_standard.sbatch`, pinned to one node class | — | — |
| DIA-NN searches (A1, A2) | no | `scripts/slurm/diann_pilot.sbatch` | — | accepts the DIA-NN licence and installs the binary once (§14.8) |
| De novo (A3: Cascadia, fallback Casanovo) | no | `gpu` partition, 1 × L40S | — | — |
| R golden references (fgsea, goseq) | yes, conda R env | or a short batch job | on demand | — |

### 3.8 Long-running work goes to Slurm

Anything that may exceed 10 minutes, or needs more than the session allocation, MUST run as a batch job. Every batch script lives in `scripts/slurm/`, starts from this header and sources the environment file:

```bash
#!/bin/bash
#SBATCH --job-name=gippa-<step>
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIPPA_REPO"
# ... the work ...
```

- **Submit and record.** `jobid=$(sbatch --parsable scripts/slurm/<x>.sbatch)`, then write the job ID, purpose and log path under "In progress" in `PROGRESS.md`, so a new session can pick it up.
- **Poll** with short commands every few minutes: `squeue -j $jobid -o '%i %T %M %l %R'`; when it leaves the queue, `sacct -j $jobid --format=JobID,State,Elapsed,MaxRSS,ExitCode` and `tail -n 50` of the log. Never wait on a long job in the foreground (`sbatch --wait`, `srun` of long work).
- **Quick checks** that need more than the session allocation use `--qos=dev` (≤ 8 h, 8 CPU, 32 GB, one job at a time).
- **Nextflow.** The Nextflow head process runs as its own batch job (`scripts/slurm/gippa_run.sbatch`: 2 CPU, 8 GB, up to 2 days) with `-profile farmshare,apptainer`; Nextflow then submits each task as a Slurm job. Small runs on the 100-kb genome may use `-profile farmshare_local` inside the session allocation.
- **Resume.** Prefer `-resume`; the work directory under `local/work` survives the session.
- **Clean up.** After a successful run, `nextflow clean -f -k` and delete `OUT/work/`.
- **Measured numbers.** `sacct` elapsed time and MaxRSS go into `PROGRESS.md` › Measured numbers, with the node class.

### 3.9 Disk budget in `local/` (`$GIPPA_ROOT`, target ≤ 400 GB)

| Item | Size (approx.) |
|---|---|
| uv environment and uv-managed Python | 2–3 GB |
| micromamba tool environments (Bakta, miniprot, BLAST+, TransTermHP, OSTIR, Promoter Calculator, UniOP, TerminatorNet, Promotech, geNOMAD, codon-bias, R golden); run `micromamba clean -a` after builds | 6–9 GB |
| Apptainer images, pinned by digest (antiSMASH alone is about 9 GB) | 15–25 GB |
| Bakta full DB v6.0 | about 84 GB unpacked; delete the archive after its checksum passes |
| Bakta light DB v6.0 (kept for fast integration tests) | 3.9 GB unpacked; delete the 1.3 GB archive after unpacking (MD5 `4a6e059ded39e9c5537ef4137d2f5648`) |
| geNOMAD DB v1.9 | about 0.8 GB |
| antiSMASH databases, Promotech RF-HOT model | measure on first download; record in `PROGRESS.md` |
| Raw MS data, 3 runs of PXD078459 | 5–10 GB (estimate) |
| DIA-NN reports, spectral libraries, de novo outputs | measure; keep only the converted generic tables after the pilot |
| Nextflow work directories | up to tens of GB per run (genome-wide Promotech writes about 5.8 GB per strand); clean after each run |
| Real genomes (K-12, Nissle; FASTA, gz) | < 10 MB |

The probe warns above 400 GB. `local/` is never committed (R11), so treat it as disposable: anything that cannot be regenerated by `gippa setup`, `scripts/build_envs.sh` or a recorded command belongs in git (within R11) or in a report.

### 3.10 GitHub from FarmShare

- **Remote.** `origin` = `git@github.com:khoa-yelo/gippa.git`. Git itself uses FarmShare's SSH key, which is already authorised (`ssh -T git@github.com` greets `khoa-yelo`).
- **Authentication for `gh`.** `gh` comes from the `gippa-wf` environment (§3.4). Khoa authenticates once by typing `! gh auth login` in the session (git protocol SSH, browser or token) and sets `git config user.name` / `user.email` for the repository. The token stays in `~/.config/gh`; never print, log or commit it (R12).
- **Branches.** Push to feature branches only. Never push to `main` directly, never push tags, never delete branches that are not yours.
- **Pull requests:** `gh pr create --base main --head <branch> --title … --body-file …`; merge with `gh pr merge <n> --squash` only when the merge rules in §14.3 are met.
- **Issues:** `gh issue create`; the backlog may also live in `docs/backlog.md`.
- **Actions:** `gh run list --branch <b>`, `gh run view <id> --log-failed`, `gh pr checks <n>`. Every workflow also posts a short result summary (failed tests, gate metrics) as a PR comment, so results are readable without logs.
- **FarmShare jobs report to GitHub.** Nightly and weekly Slurm jobs (§12) post their summary as a comment on the tracking issue "Nightly (FarmShare)" with `gh issue comment`, and on the PR that triggered them, if any.
- **Releases.** Tags and releases come from `release.yml` (workflow_dispatch) or from Khoa.
- **Actions minutes.** Heavy suites run on FarmShare, so GitHub Actions only carries the per-push and per-PR checks. Keep them under about 15 minutes per run on a private repository.

---
## 4. Repository layout

```text
gippa/
├── GIPPA_SPEC.md                this file (agents never edit it; R13)
├── AGENTS.md  CLAUDE.md        agent instructions (Appendix A.1, A.2)
├── PROGRESS.md                 session-to-session state (Appendix A.3)
├── README.md  LICENSE (MIT)  CHANGELOG.md  CITATION.cff  NOTICE.md (third-party tools and licences)
├── pyproject.toml  uv.lock     src layout, `gippa` console script, extras: dev, docs, viz
├── .claude/
│   ├── settings.json           permissions, PreToolUse guard, SessionStart hook (Appendix A.4)
│   ├── agents/                 subagent definitions (Appendix A.5)
│   └── hooks/protect_paths.py
├── ci/
│   ├── gates.yaml              coverage, budgets, benchmark tolerance (protected)
│   └── protected.txt           globs of protected files (protected)
├── .github/
│   ├── workflows/              ci.yml, nightly.yml, weekly.yml, containers.yml, golden-r.yml, release.yml
│   ├── pull_request_template.md
│   └── ISSUE_TEMPLATE/feature.md
├── src/gippa/
│   ├── cli.py                  Click entry point
│   ├── config/                 pydantic models, defaults.yaml, JSON-schema export
│   ├── model/                  Interval, Replicon, Gene, Feature, coordinates (0-based internally)
│   ├── io/                     fasta, gff3, bed, probed, vcf, tsv writer, peptide table, schemas
│   │   └── converters/         diann.py, casanovo.py, instanovo.py, cascadia.py, generic.py
│   ├── preflight/
│   ├── annotate/               one module per tool: build inputs + parse outputs (never import the tool)
│   ├── mapping/                sixframe, exact (Aho–Corasick), tolerant (numba), decoys, miniprot,
│   │                           tblastn, placement_fdr, loci, contaminants
│   ├── genetable/              gene table, detection tiers, detectability, abundance
│   ├── context/  operons/  variants/  novel/  function/  determinants/
│   ├── stats/                  multiple testing, target–decoy, MAYU, permutation, bootstrap,
│   │                           scan statistic, Moran/join-count, Wallenius ORA, fgsea multilevel
│   ├── compare/  makedb/  setupdb/
│   ├── viz/                    style.py + one module per plot ID
│   ├── report/                 Jinja2 templates, HTML assembly
│   ├── steps/                  step registry: one function per pipeline step, files in → files out
│   ├── runner/                 Nextflow launcher + developer-only local executor
│   ├── data/                   bundled small data: residue masses, mass equivalents, COG categories,
│   │                           calibration.yaml (cut-offs with provenance)
│   └── schemas/                Frictionless table schemas, manifest JSON Schema
├── workflow/
│   ├── main.nf  nextflow.config
│   ├── conf/                   base.config, test.config, farmshare.config, farmshare_local.config, standard.config
│   ├── modules/local/<tool>/   main.nf + meta.yml (one process per external tool or Python step group)
│   ├── subworkflows/local/     annotate.nf, mapping.nf, analysis.nf, report.nf
│   ├── envs/<tool>.yml         conda env per tool (conda-forge + bioconda, exact pins)
│   ├── containers/<tool>/      Dockerfile per tool (built in CI, pushed to GHCR, referenced by digest; run with Apptainer on FarmShare)
│   └── tests/                  nf-test suites
├── tests/
│   ├── unit/ property/ calibration/ golden/ contract/ snapshot/ integration/ e2e/
│   ├── sim/                    synthetic-truth generator (package `gippasim`, §11.3)
│   ├── data/                   fixtures ≤ 50 MB: tiny genomes, recorded/<tool>/<version>/…, golden/…
│   └── conftest.py             capability markers (§3.6), seeds, fixtures
├── benchmarks/                 pytest-benchmark suites, baseline.json
├── validation/                 pilot kit (§15), RegulonDB calibration scripts, reports; data/ gitignored
├── scripts/                    farmshare_env.sh, bootstrap_farmshare.sh, session_start.sh, probe_env.py, build_envs.sh, ci/,
│                               slurm/ (sbatch scripts: gippa_run, record_bakta_full, nightly, weekly, benchmark_standard, diann_pilot, denovo_pilot)
├── docs/                       MkDocs: user guide, methods/ (one note per statistical method), adr/, specs/
└── local/                      GITIGNORED ($GIPPA_ROOT, §3.3): db/, data/, tools/, micromamba/, venv/, containers/, cache/,
                                work/, logs/, worktrees/, ci/. Never committed (R11)
```

---

## 5. Interfaces

### 5.1 Command line

All commands share `--config FILE`, `--set key=value` (repeatable, dotted keys), `--log-level`, `--version`.

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | invalid input (preflight errors) |
| 2 | environment error: a tool, database or model is missing |
| 3 | warnings were raised and `--strict` was given |
| 4 | a pipeline step failed |
| 5 | internal error (a bug; the message asks to report it) |

```bash
# databases and models (wrappers around each tool's own downloader; writes DIR/db_manifest.json)
gippa setup --db-dir DIR [--bakta light|full|skip] [--genomad/--no-genomad] [--antismash/--no-antismash]
           [--promotech-models/--no-promotech-models] [--go] [--kegg] [--cog] [--contaminants] [--check]

# input validation only (seconds)
gippa preflight --genome G.fasta [--proteins P.faa] [--peptides T.tsv] [--abundance A.tsv]
               [--annotation X] [--outdir DIR] [--strict]

# the analysis
gippa run --genome G.fasta [--peptides T.tsv] [--proteins P.faa] [--abundance A.tsv]
         [--annotation bakta.json|ann.gff3|ann.gbff] [--reference-peptides R.tsv] [--replicon-depth D.tsv]
         --outdir OUT [--db-dir DIR] [--threads 16] [--memory 64.GB]
         [--modules core,mobile,bgc,function] [--genome-relation related|same]
         [--abundance-missing median|exclude|min|zero|mean]
         [--abundance-type ibaq|top3|maxlfq|intensity|unknown] [--abundance-scale auto|log2|linear]
         [--promoters targeted|genome] [--terminators terminatornet|transtermhp] [--add-bactermfinder]
         [--publication] [--seed INT] [--parquet] [--sensitivity basic|full|off] [--strict]
         [--profile farmshare,apptainer|farmshare,conda|farmshare_local|apptainer|conda|docker] [--resume] [--executor nextflow|local]

# optional strain search database for the user's own search (no FDR is built on it)
gippa makedb --genome G.fasta [--annotation bakta.json] --outdir DB [--db-dir DIR]
            [--six-frame-orfs] [--min-orf-aa 20] [--patch annotation_patch.gff3] [--no-contaminants]

# engine output → generic peptide table
gippa convert --from diann|casanovo|instanovo|cascadia|generic INPUT [INPUT ...] --out peptides.tsv
             [--abundance-out abundance.tsv] [--source database|denovo] [--search-db same|other]

# run-vs-run comparison (same genome; checked through the genome SHA-256 in each manifest)
gippa compare OUT_A OUT_B [OUT_C ...] --outdir CMP [--labels A1,A2,A3]

# rebuild plots and report from existing tables
gippa report OUT [--no-interactive] [--plots V01,V02,...]
```

Hidden commands: `gippa _step NAME --config resolved.yaml --in key=path ... --out key=path ...` (used by Nextflow, §6.3) and `gippa version --json`.

`--executor local` is a developer feature (§6.4): it runs only the Python steps and requires `annotation.precomputed` layers. It is documented under "Developer guide", not in the user guide.

`--publication` sets `stats.permutations: 10000` and `stats.bootstrap: 10000`, and writes the setting into the manifest.

### 5.2 Configuration (defaults)

Config is a pydantic v2 model. The YAML below is `src/gippa/config/defaults.yaml` and MUST stay in sync with the model (a unit test enforces it). CLI flags map to keys, and `--set` overrides any key. The resolved config is written to `run_manifest.json`.

```yaml
schema_version: "1.0"
seed: 20261014
sample_name: auto                                # auto = stem of the peptide (or protein) file
threads: 16
modules: [core, mobile, bgc, function]        # core is always on
genome_relation: related                       # related | same (§7.2.8)

input:
  peptide_source: auto            # auto | database | denovo; a per-row `source` column wins
  search_db: auto                 # auto | same | other; auto = same if ≥95 % of protein_ids match the gene model
  decoy_prefixes: ["rev_", "REV_", "DECOY_", "decoy_", "XXX_", "##"]
  contaminant_prefixes: ["Cont_", "CON__", "CONT_", "cRAP"]
  q_value_max: 0.01               # peptides with q above this are not evidence (still used to estimate FDR)
  peptide_len: [7, 40]            # shorter are dropped; longer protein records go to miniprot
  i_equals_l: true

preflight:
  max_contigs_warn: 10
  min_total_bp_warn: 100000
  min_n50_warn: 100000
  min_largest_contig_warn: 1000000
  small_replicon_bp: 500000
  gc_range_warn: [0.25, 0.75]
  n_run_min_bp: 10
  abundance_missing_warn_frac: 0.20

annotation:
  bakta: {db: auto, complete: auto, extra_args: []}   # complete=auto: --complete if all contigs circular
  translation_table: 11
  start_codons: [ATG, GTG, TTG]
  region_min_overlap_frac: 0.5                        # gene in region if ≥50 % of its CDS lies inside
  uniop: {pair_threshold: 0.5, fallback: distance}    # threshold replaced by calibration (§15.4)
  tu_upstream_bp: [-300, 50]                          # TU upstream region around the first start codon
  promoters:
    mode: targeted                                    # targeted | genome
    window_nt: 40
    window_3p_range: [-200, 20]                       # 3′ ends of 40-nt windows relative to each start codon
    promotech_yes: 0.5                                # "promoter between" threshold (Module 4)
  promoter_calculator: {enabled: true}
  terminators:
    tool: terminatornet                               # terminatornet | transtermhp
    cross_check_transtermhp: true
    bactermfinder: false
    terminatornet_yes: 0.5                            # placeholder; replaced by RegulonDB calibration (§15.4)
  ostir: {upstream_nt: 50, downstream_nt: 50, anti_sd: auto}   # auto = 3′ end of the strain's 16S rRNA
  codon: {cai_reference: ribosomal_proteins}
  precomputed:                                        # developer/test injection (§6.4); null = run the tool
    bakta_json: null
    operon_pairs: null
    promoters: null
    tx_rates: null
    terminators: null
    rbs: null
    regions: null

mapping:
  exact: {min_len: 7, max_len: 40}
  proteins:
    min_identity: 0.80
    min_coverage: 0.80
    near_identity: 0.95
    miniprot_args: "-S -N 5 --outs=0.95"              # -T 11 is applied at indexing time
    tblastn_args: "-db_gencode 11 -seg no -evalue 1e-5"
  tolerant:
    enabled: true
    min_len: 10
    max_diff: {"10-11": 1, "12-40": 2}
    swaps: true
    mass_equivalents: true
    frameshift: {enabled: true, min_side_aa: 4}
  decoys: {database_input: reverse, denovo_input: shuffle, shuffle_window_aa: 50}
  length_classes: [[7, 7], [8, 8], [9, 9], [10, 11], [12, 14], [15, 19], [20, 40]]

detection:
  min_peptides_sequences_only: 2
  mayu_size_bins: 10

abundance:
  method: ibaq                     # ibaq | top3, when peptide intensities exist
  type: unknown                    # protein-level input type: ibaq | top3 | maxlfq | intensity | unknown
  scale: auto                      # auto | log2 | linear
  missing: median                  # median | exclude | min | zero | mean
  sensitivity: [min, exclude]
  enzyme: trypsin_p                # trypsin_p (cleaves before P too; DIA-NN default) | trypsin
  observable_len: [7, 30]
  run_aggregation: median          # several runs → median per peptide across runs

detectability:
  covariates: [n_observable_peptides, length_aa, gravy, pi, tm_segments]
  spline_knots: 5
  folds: 5                         # GroupKFold by TU, cross-fitted

stats:
  alpha:                           # the single source of every significance level
    identification: 0.01           # gene-level detection q
    placement: 0.01                # placement q (Module 2)
    novel_orfs: 0.01
    variants: 0.01
    context: 0.05
    enrichment: 0.05
    operon_flags: 0.10
    determinants: 0.05
  permutations: 1000
  bootstrap: 1000
  besag_clifford_h: 10

context:
  region_types: [prophage, plasmid, bgc, crispr, amr]
  exhaustive_shift_max_genes: 20000
  moran_lags: [1, 20]
  islands: {min_tus: 2, max_tus: 20, detrend_span_frac: 0.10, block_len: auto}

operons:
  decoys: [convergent, sandwich]                      # sandwich = → ← →
  distance_bins_bp: [-.inf, 0, 50, 150, 300, .inf]
  pi0: 1.0

variants:
  high_min_peptides: 2
  mass_equivalent_tol_da: 0.05                       # includes K/Q (0.036 Da)

novel_orfs:
  min_codons: 10
  min_peptides_without_decoys: 2
  rbs_min_quantile: 0.10                              # RBS support threshold (§7.6.3)

function:
  go_slim: generic
  ora_size: [5, 500]
  gsea_size: [10, 500]
  gsea: {sample_size: 101, eps: 0.0, n_perm_simple: 1000}

determinants:
  outer_folds: 5
  repeats: 5
  permutations: 1000
  lightgbm: {n_estimators: 500, learning_rate: 0.03, num_leaves: 15, min_child_samples: 20,
             subsample: 0.8, subsample_freq: 1, colsample_bytree: 0.8, deterministic: true}
  elastic_net: {l1_ratio: [0.1, 0.5, 0.9], inner_folds: 3}

report:
  interactive: true
  max_locus_plots: 50
  formats: [svg, png, pdf]
  dpi: 300

output:
  parquet: false
  sensitivity: basic                # basic | full | off (§8.6)
```

### 5.3 Input formats

**Genome FASTA.**

- The first token of each header is the contig ID.
- `circular=true` or `[topology=circular]` in a header marks a circular replicon. Contigs without such a tag are linear unless `annotation.bakta.complete: true` declares every contig a complete, circular replicon.
- Unsafe ID characters are mapped to `[A-Za-z0-9._-]`; the mapping is saved.
- gzip input is accepted.

**Proteins FASTA.**

- The ID is the first header token. UniProt `sp|ACC|NAME` and `tr|ACC|NAME` resolve to `ACC`. A protein group written `A;B` is split, keeping the group ID.
- Records of ≤ 40 aa are routed to the peptide path, and records < 7 aa are dropped with a warning.
- Allowed letters: the 20 amino acids plus U. X, B, Z and J are counted. An internal `*` warns.

**Generic peptide table (`peptides.tsv`).** UTF-8 TSV, one header row. `gippa convert` writes it; users may also write it by hand.

| Column | Type | Required | Meaning and validation |
|---|---|---|---|
| `peptide` | str | **yes** | stripped sequence, uppercase A–Z only. Modifications are removed by the converter. |
| `modified_peptide` | str | no | engine notation or ProForma; kept for provenance |
| `protein_ids` | str | no | `;`-separated engine protein IDs. Decoy and contaminant prefixes are recognised from config. |
| `q_value` | float [0, 1] | no | precursor or peptide q-value; for DIA-NN the global q-value |
| `pep` | float [0, 1] | no | posterior error probability |
| `score` | float | no | engine or de novo score; higher is better |
| `is_decoy` | bool | no | decoy row; requires `protein_ids` for gene-level pairing |
| `intensity` | float ≥ 0 | no | linear precursor or peptide intensity, per run |
| `run` | str | no | run or sample name |
| `charge` | int | no | precursor charge |
| `spectrum_id` | str | no | scan or spectrum identifier, unique within `run` (G2) |
| `source` | enum | no | `database` or `denovo`; default from `input.peptide_source` |
| `search_db` | enum | no | `same` or `other`; default from `input.search_db` |

Aggregation, done once after reading:

- **Group by stripped `peptide`.** I/L are kept as written; equivalence is applied only in matching.
- **Combining rows of one peptide:**
  - `q_value` and `pep` take the minimum over rows;
  - `score` takes the maximum;
  - `is_decoy` is true only if all rows are decoys;
  - `protein_ids` is the union.
- **Intensity:** summed over charges and modified forms within a run, then the median across runs (`abundance.run_aggregation`).
- **Per-run detail:** per-run presence is kept in `work/peptide_runs.parquet` for the detection-limit model (§7.4.3).

**Abundance table (`abundance.tsv`).** The first column is a protein ID that matches `proteins.faa`, a locus tag, or a gene ID. Exactly one numeric value column is allowed in version 1; extra columns warn and are ignored. The scale is detected automatically (`log2` if any value is negative or all values are ≤ 50; otherwise `linear`) unless `--abundance-scale` is given. Negative values are an error only on a linear scale. `--abundance-type` declares iBAQ, Top3, MaxLFQ or raw intensity (default `unknown` = relative).

**Existing annotation (`--annotation`).**

- *Bakta JSON:* used as is, and Bakta is not run. The Bakta version and DB version are taken from the JSON. This route is how recorded full-database fixtures and Khoa's own annotation enter.
- *GFF3 or GBFF from another annotator:* Bakta runs with `--regions` on that file's CDS, so the gene coordinates match. `annotation/id_map.tsv` maps the source IDs to Bakta locus tags by identical CDS coordinates.

**Reference peptides (`--reference-peptides`).** Same schema as above, from a same-strain database search on the same spectra. Used only for de novo precision (c in §7.2.6) and ALT error rates (§7.5.5).

**Replicon depth (`--replicon-depth`).** TSV `contig`, `mean_depth`. Gives an optional copy-number offset for plasmid contrasts (§7.3.5).

### 5.4 Converters (`gippa convert`)

Each converter is a pure function from input files to the generic table, with a contract test against a recorded sample file (§11.4). Until a real sample exists, the converter is written from the tool's documentation and its contract test is marked `pending_sample`; the pilot (§15) supplies DIA-NN and Cascadia samples.

| Tool | Input | Mapping to the generic table |
|---|---|---|
| DIA-NN 2.x | `report.parquet` (or `.tsv`) | `peptide` ← `Stripped.Sequence`<br>`modified_peptide` ← `Modified.Sequence`<br>`protein_ids` ← `Protein.Ids`<br>`q_value` ← `Global.Q.Value`, else `Q.Value`<br>`pep` ← `PEP`<br>`is_decoy` ← `Decoy` (present with `--report-decoys`)<br>`intensity` ← `Precursor.Normalised`, else `Precursor.Quantity`<br>`run` ← `Run`; `charge` ← `Precursor.Charge`; `source` = database<br>`--abundance-out` writes `PG.MaxLFQ` per protein group, with type `maxlfq` |
| Casanovo 5.x | mzTab | `peptide` ← stripped `sequence`<br>`score` ← `search_engine_score[1]`<br>`spectrum_id` ← `spectra_ref`; `run` ← the `ms_run` location; `charge` ← `charge`<br>`source` = denovo |
| InstaNovo 1.x | CSV | `peptide` ← stripped prediction<br>`score` ← sequence log-probability (or confidence)<br>`spectrum_id` ← scan or spectrum column; `run` ← experiment name; `charge` ← precursor charge<br>`source` = denovo |
| Cascadia 0.0.x | its result file (mzTab-like) | as Casanovo; `source` = denovo |
| generic | TSV | validates, normalises and reports problems |

Stripping removes bracketed or parenthesised modifications (`(UniMod:35)`, `[+15.995]`, `[UNIMOD:35]`, `+15.995`, `n[42]`) and lower-case residues, and keeps only uppercase letters.

---
## 6. Architecture

### 6.1 Data flow

```text
inputs ─► preflight ─┬─► Bakta (or --annotation JSON) ─► gene model (locus tags) ─┬─► UniOP → pairs → TUs
                     │                                                             ├─► Promotech (targeted windows)
                     │                                                             ├─► Promoter Calculator (TU upstream)
                     │                                                             ├─► TerminatorNet (+ TransTermHP; BacTermFinder opt-in)
                     │                                                             ├─► OSTIR (anti-SD from 16S) ; codon-bias
                     │                                                             └─► antiSMASH [bgc]
                     ├─► geNOMAD [mobile] (parallel with Bakta)
                     └─► six-frame index ─► exact peptides ─► tolerant (+ decoys) ─► placement FDR ─┐
                         proteins > 40 aa ─► miniprot ─► tblastn (leftovers) ───────────────────────┤
                                                                                                    ▼
        context layers + loci ─► gene table ─► detection tiers ─► detectability (cross-fitted p_detect)
                                     │
        ┌──────────────┬─────────────┼──────────────┬───────────────┬───────────────┐
        ▼              ▼             ▼              ▼               ▼               ▼
   context stats    operons      variants      novel ORFs     function [fn]    determinants
   (M3)             (M4)         (M5)          (M6; 2nd-round Promotech/OSTIR)   (M8)
        └──────────────┴─────────────┴──────────────┴───────────────┴───────────────┘
                                     ▼
                  sensitivity reruns ─► tables + schemas ─► plots ─► report.html ─► run_manifest.json
```

The only join between the mapping track and the annotation track is the Bakta locus tag. Analysis modules read only the gene table and the layer tables, so each optional module can be switched off on its own.

### 6.2 Step registry

Every step is a Python function or an external process with declared file inputs and outputs. Python steps live in `src/gippa/steps/` and are registered with a decorator that declares names, inputs, outputs and the module flag. Intermediate files are Parquet under the work directory; final files follow §9.

| Step | Kind | Main inputs → outputs | Module |
|---|---|---|---|
| `preflight` | py | raw inputs → `preflight.json`, `genome.norm.fa`, `contig_map.tsv`, `replicons.tsv`, `peptides.norm.parquet`, `peptide_runs.parquet`, `proteins.long.faa`, `abundance.norm.parquet` | core |
| `bakta` | ext | `genome.norm.fa`, `replicons.tsv` [, `regions.gff3`] → Bakta JSON, GFF3, GBFF, FAA, TSV | core (skipped when a Bakta JSON is given) |
| `parse_bakta` | py | Bakta JSON → `genes.parquet`, `features.parquet` (rRNA, tRNA, CRISPR, oriC, gaps), `rrna16s.fa`, `trna_counts.tsv`, `xrefs.parquet` | core |
| `uniop` | ext | Prodigal-style GFF of Bakta CDS → `uniop.pred` | core |
| `operon_pairs` | py | `uniop.pred` (or the distance model) → `pairs_raw.parquet`, `tus.parquet` | core |
| `promoter_windows` | py | genes, genome → `windows_NNN.fa` (40-nt chunks) | core |
| `promotech` | ext | window chunk → score chunk | core |
| `promoter_calculator` | ext | one FASTA per TU upstream region → one CSV per region | core |
| `terminatornet` | ext | genome → terminator table (both strands, `-s`) | core |
| `transtermhp` | ext | genome + `.coords` → `.tt` | core (cross-check, or main tool if `--terminators transtermhp`) |
| `bactermfinder` | ext | genome → calls | opt-in |
| `ostir` | ext | `ostir_input.csv` → `ostir.csv` | core |
| `codon_usage` | py | genes, `trna_counts.tsv` → `codon.parquet` (uses codon-bias in-process; MIT) | core |
| `genomad` | ext | genome → summary TSVs | mobile |
| `antismash` | ext | Bakta GBFF → antiSMASH JSON | bgc |
| `context_layers` | py | parsed layers → `annotation/*` (§9.1), `layers.parquet` | core |
| `sixframe_index` | py | genome → six-frame and decoy translations (`.npz`) | core |
| `map_exact` | py | peptides, index → `hits_exact.parquet` | core |
| `miniprot` | ext | genome, `proteins.long.faa` → GFF | core |
| `tblastn` | ext | genome, leftover proteins → tabular hits | core |
| `map_tolerant` | py | leftover peptides ≥ 10 aa, index → `hits_tolerant.parquet` (target + decoy) | core |
| `placement` | py | all hits, decoy hits, input q/PEP → `placement.parquet`, `mapping_fdr.tsv` | core |
| `loci` | py | placements, genes → `peptide_hits.tsv`, `protein_hits.tsv`, `mappings.*`, `unmapped.tsv` | core |
| `novel_candidates` | py | loci → candidate ORFs, extra Promotech and OSTIR windows | core |
| `promotech_novel`, `ostir_novel` | ext | candidate windows → scores | core |
| `gene_table` | py | genes, layers, loci, abundance → `genes_pre.parquet` | core |
| `detectability` | py | gene table → `p_detect` (cross-fitted) | core |
| `context_stats` | py | gene table → `regions_stats`, `autocorrelation`, `islands`, `detection_summary`, ori–ter | core |
| `operons` | py | pairs, gene table → `pairs`, `operons`, `operon_flags` | core |
| `variants` | py | hits, genome, genes → `pileup`, `variants.tsv`, `variants.vcf` | core |
| `novel_orfs` | py | candidates + second-round scores → `novel_orfs.*`, `annotation_patch.gff3` | core |
| `function` | py | gene table, xrefs, GO/KEGG resources → `category_summary`, `enrichment_*` | function |
| `determinants` | py | gene table, features → `determinants_*` | core |
| `sensitivity` | py | reruns the `context_stats`, `operons`, `function` and `determinants` logic under `min` and `exclude` → `sensitivity.tsv` | core (`--sensitivity off` disables it) |
| `finalize` | py | everything → final TSV/GFF/BED/VCF, `genes.tsv`, `schemas/datapackage.json` | core |
| `plots` | py | tables → `plots/*` | core |
| `report` | py | tables, plots → `report/report.html` | core |
| `manifest` | py | versions, trace, config, hashes → `run_manifest.json` | core |

### 6.3 Nextflow workflow

- **Launch.** `gippa run` validates inputs and config in Python, writes `OUT/work/params.json` and `resolved_config.yaml`, then launches:
  `nextflow run <pkg>/workflow/main.nf -params-file params.json -profile <profile> -work-dir OUT/work -with-trace OUT/logs/trace.tsv -with-timeline OUT/logs/timeline.html -with-report OUT/logs/nf_report.html [-resume]`.
  Nextflow never parses user parameters on its own. Do not use Nextflow plugins in version 1, so runs never depend on the plugin registry.
- **Processes.** One process per external tool, and one per Python step or step group. Python processes call `gippa _step <name>`.
- **Process conventions.** Each process declares:
  - `tag`;
  - a `label` from `process_single` (1 CPU, 2 GB), `process_low` (2 CPU, 6 GB), `process_medium` (8 CPU, 16 GB) or `process_high` (16 CPU, 48 GB);
  - `conda` (`workflow/envs/<tool>.yml`) and `container` (`ghcr.io/<owner>/gippa-<tool>@sha256:…`);
  - an nf-core style `versions.yml` emission;
  - a `stub:` block that writes minimal valid outputs, used by nf-test in stub mode.
- **Resources.** `process.resourceLimits = [cpus: params.threads, memory: params.memory]`. Memory-heavy tools (Bakta, geNOMAD, antiSMASH, genome-wide Promotech) declare memory, so Nextflow never oversubscribes RAM.
- **Profiles:**
  - `apptainer` (`apptainer.enabled = true`, `apptainer.autoMounts = true`, images `docker://…@sha256:…` cached in `$NXF_APPTAINER_CACHEDIR`), `conda` (micromamba: `conda.useMicromamba = true`), `docker` (for users elsewhere; not tested on FarmShare);
  - `farmshare`: `executor = 'slurm'`, `queue = 'normal'`, `clusterOptions = '--qos=normal'`, `executor.queueSize = 20`, `executor.submitRateLimit = '10/1min'`, job names `gippa-<process>`, `process.scratch = '/tmp'` (node-local), `resourceLimits = [cpus: 64, memory: 240.GB, time: 48.h]`; work directory `$NXF_WORK` = `local/work`;
  - `farmshare_local`: local executor inside the session allocation, capped at `SLURM_CPUS_ON_NODE` CPUs and 30 GB;
  - `test` (the 100-kb synthetic genome with stubs, or with real tools when available);
  - `standard` (local executor, all cores).
- **Error strategy.** Retry once on exit codes 137 or 143 with memory × 2; otherwise terminate, and `gippa run` exits 4 naming the failed step and its log.
- **Optional modules.** `mobile`, `bgc` and `function` are switched with `params.modules`. Downstream steps accept a missing optional layer and record `layer_missing` in the manifest.
- **Publishing.** `publishDir` copies final files into `OUT/` (mode `copy`); `OUT/work/` can be deleted after success (`gippa run` prints how).

### 6.4 Developer-only local executor

`gippa run --executor local` runs the Python steps in topological order in one process. Every external layer must come from `annotation.precomputed` (Bakta-like JSON, operon pairs, promoters, Tx rates, terminators, RBS, regions), and a missing layer is an error naming the key. The simulator (§11.3) writes all precomputed layers, so end-to-end tests on simulated genomes run in seconds without Java, Docker or conda. Both executors call the same step functions, and an integration test checks that they produce byte-identical tables on the synthetic genome.

### 6.5 Internal model and coordinates

- **`Replicon`:** `id`, `length`, `topology` (`circular` | `linear`), `kind` (`chromosome` | `plasmid` | `unknown`; from Bakta's replicon table or size).
- **`Interval`:** `contig`, `start`, `end`, `strand`, 0-based half-open. On a circular replicon, a feature crossing the origin is stored with `end > length` (an unwrapped coordinate) and `crosses_origin = True`; helpers split it into `[start, length)` and `[0, end − length)`. Linear contigs never wrap.
- **Writers convert once:**
  - TSV and GFF3 use 1-based closed coordinates (`start + 1`, `end`). Origin-crossing features keep `end > length`, as the GFF3 convention for circular sequences allows, and the `##sequence-region` line and region feature carry `Is_circular=true`.
  - BED and proBed are 0-based half-open; origin-crossing features are written as two records with `name` suffixes `_part1` and `_part2`.
  - VCF positions are 1-based within `[1, length]`.
- **IDs:**
  - genes: Bakta locus tags;
  - novel ORFs: `NOVEL_<contig>_<start1>_<strand>` (strand `p` or `m`);
  - regions: `<type>_<n>`, numbered by contig order then start;
  - TUs: `TU_<n>`;
  - islands: `ISL_<n>`.
- **Partial genes.** Genes cut by a linear contig end are `partial` and excluded from detection statistics; they still appear in `genes.tsv`.
- **Sort order.** All tables sort by contig (as in the input FASTA), then start, end, strand and ID.

### 6.6 Determinism

- **Seeds.** Each step draws from `numpy.random.default_rng(SeedSequence([config.seed, crc32(step_name)]))`, so adding a step never changes another step's random stream.
- **Number formats.** Floats are written with 6 significant digits (`%.6g`); integers are written as integers; booleans as `true` and `false`.
- **Threading.** Numba kernels that reduce floats in parallel MUST use a fixed chunking with a deterministic reduction order, or integer arithmetic. LightGBM uses `deterministic=true`, `force_row_wise=true`, a fixed `num_threads` and a seed.
- **No timestamps or absolute paths** in `tables/`, `mapping/` or `annotation/*.tsv`; these belong only in `run_manifest.json` and `logs/`.

### 6.7 Logging and manifest

- **Logs.** `logs/<step>.log` (human-readable) and `logs/events.jsonl`, one JSON object per event: `ts`, `step`, `level`, `event`, `data`.
- **Manifest.** `run_manifest.json` follows a JSON Schema in `src/gippa/schemas/` and holds:
  - the GIPPA version and git commit; the command line; the resolved config;
  - inputs (path, SHA-256, size);
  - tools (name, version, container digest or conda lock hash);
  - databases (name, version, checksum);
  - seeds;
  - per-step wall time, CPU time and peak RSS, from the Nextflow trace;
  - `schema_version`, warnings and flags, layers missing, and the sensitivity summary.

---
## 7. Module specifications

Each module section lists inputs, parameters (config keys in §5.2), algorithm, outputs (schemas in §9), edge cases and required tests. Test IDs such as `T-MAP-03` are suggestions; test authors may add more.

### 7.0 Module 0: preflight

**Purpose.** Validate inputs and environment in seconds. Write `preflight/preflight.json` and `preflight.md`. Errors stop the run (exit 1). Warnings become flags carried into every table and the report.

| Input | Check | Rule (default) | Severity |
|---|---|---|---|
| Genome | parses as FASTA | ≥ 1 record, IUPAC nucleotides only | Error |
| Genome | contig IDs | unique first tokens; unsafe characters mapped to `[A-Za-z0-9._-]`, map saved | Error (duplicates) / Warn (renamed) |
| Genome | contig count | > 10 | Warn |
| Genome | total length | < 100 kb | Warn |
| Genome | contiguity | N50 < 100 kb or largest contig < 1 Mb; median length recorded for information | Warn |
| Genome | gaps and ambiguity | runs of ≥ 10 N and other IUPAC codes, with positions | Warn |
| Genome | GC content | outside 25–75 % | Warn |
| Genome | topology | `circular=true` or `[topology=circular]` → Bakta replicon table | Info |
| Genome | small replicons | contigs < 500 kb listed as putative plasmids or fragments | Info |
| Proteins | IDs | non-empty, unique; UniProt `sp\|ACC\|NAME` and `A;B` groups parsed | Error |
| Proteins | alphabet | 20 aa + U; X, B, Z, J counted; internal `*` | Warn |
| Proteins | length classes | ≤ 40 aa → peptide path, > 40 aa → proteins; < 7 aa dropped | Info / Warn |
| Proteins | duplicate sequences | same sequence under several IDs | Info |
| Peptides | schema | required `peptide`; types and ranges of optional columns (§5.3) | Error |
| Peptides | content | non-letters after stripping; length < 7 dropped; share of decoys; q-value range | Warn |
| Peptides | evidence tier | which detection tier the columns allow (§7.3.2) | Info |
| Abundance | shape | first column = ID, one numeric column; extra columns warn (multi-sample later) | Error / Warn |
| Abundance | join | IDs missing from proteins, and proteins missing from abundance | Warn |
| Abundance | values | non-numeric → error; negatives → error only on a linear scale (checked after scale detection); zeros and NA counted; > 20 % missing → warn | Error / Warn |
| Abundance | scale | log or linear auto-detected; `--abundance-scale` overrides | Info |
| Environment | tools and databases | binaries, environments and DB versions for enabled modules (Bakta DB v6 for Bakta 1.12) | Error (exit 2) |
| Environment | resources | free disk and RAM vs. module needs (genome-wide Promotech needs ≥ 24 GB) | Warn |

**Two kinds of missing abundance, kept apart:**

1. *Identified but not quantified*: detected, but with no abundance value. `--abundance-missing` applies: `median` (default), `exclude` (counted as detected but left out of abundance statistics), `min` (lowest observed value), `zero` or `mean`.
2. *Never identified*: an annotated gene with no evidence. This is not detected. It is never imputed, and the detection models handle it.

**Outputs.**

- `preflight.json`: checks with severity, value and rule; derived routing; the evidence tier; replicon table.
- `preflight.md`: a human-readable summary.
- Normalised inputs in the work directory (§6.2).

**Tests.**

- `T-PRE-*` covers every rule with one passing and one failing example.
- The exit codes are 0, 1 and 3 (`--strict`).
- The JSON report validates against its schema.
- A complete genome with two small plasmids raises no contiguity warning.

### 7.1 Module 1: annotation and context layers

Bakta runs first and defines the gene model. Every other layer either reads Bakta's genes or is joined to them by coordinates.

#### 7.1.1 Tool table (pinned)

| Layer | Tool, version | Environment | Command (essentials) | What GIPPA keeps |
|---|---|---|---|---|
| Genes, sORFs, tRNA, rRNA, ncRNA, CRISPR, oriC, pseudogenes, COG/GO/EC/KEGG cross-references | Bakta 1.12.2, DB v6.0 (full; light only for fast integration tests) | conda `bakta=1.12.2` or container `oschwengers/bakta` (pin by digest; Apptainer on FarmShare) | `bakta --db DB --output bakta --prefix P --threads N --keep-contig-headers (--complete \| --replicons replicons.tsv) [--regions regions.gff3] genome.norm.fa` | JSON (primary), GFF3, GBFF (for antiSMASH), FAA, TSV |
| Operons | UniOP, pinned commit | conda env with Python 3.11, pandas, numpy, scikit-learn; source cloned at a pinned commit when the env or image is built | `python <UniOP>/src/UniOP -f cds.gff --gff_source prodigal_gff -t outdir` (confirm the entry-point name on the pinned commit) | probability for each adjacent same-strand pair (`uniop.pred`) |
| Promoter presence | Promotech RF-HOT, pinned commit | frozen env: Python 3.6.10, TensorFlow 2.1.0, scikit-learn 0.23.1; models from `www.cs.mun.ca/~lourdes/public/PromoTech_models/` | targeted: `python promotech.py -s -m RF-HOT -f windows.fa -o out/`; genome-wide (opt-in): `-pg` then `-g` | score per 40-nt window |
| Promoter strength | Promoter Calculator 1.2.4 (PyPI `promotercalculator`, GPL) | own env | `promoter-calculator -i region.fa -o region.csv -j 1`, once per TU upstream region; regions run in parallel (G11) | σ70 Tx rate, TSS, strand |
| Terminators | TerminatorNet (Tjaden 2026), pinned commit | own env; scikit-learn pinned to the version that pickled `model.pkl` (G12) | `python TerminatorNet.py -f genome.fa -s -n N -o terms.tsv` | intrinsic terminators, confidence `pred`, hairpin ΔG |
| Terminators (fallback, cross-check) | TransTermHP 2.09 | conda `transtermhp=2.09` | `transterm -p expterm.dat genome.fa genes.coords` | terminators with confidence |
| Rho-dependent (opt-in) | BacTermFinder, pinned commit | own env | per its README | extra terminator layer, source-tagged |
| RBS / translation initiation | OSTIR 1.1.3 (GPL) | own env with ViennaRNA (conda-forge/bioconda) + `pip install ostir==1.1.3` | `ostir -i ostir_input.csv -o ostir.csv -j N` | `expression` (TIR), `dG_total`, `dG_rRNA:mRNA`, `dG_mRNA`, `dG_spacing`, `dG_standby`, `dG_start_codon`, `RBS_distance_bp` |
| Codon usage | codon-bias 0.5.0 (MIT, PyPI) | main GIPPA env, in-process | library calls | CAI, tAI, ENC, Fop, GC3, codon-pair bias |
| Prophages, plasmids | geNOMAD 1.12.0, DB v1.9 (non-commercial) | conda `genomad=1.12.0` | `genomad end-to-end --cleanup --splits 8 --threads N genome.fa out DB/genomad_db` | provirus coordinates, plasmid scores, virus gene annotations |
| BGCs | antiSMASH 8.0.4 (AGPL) | container `antismash/standalone` (pin by digest; Apptainer on FarmShare) or conda `antismash=8.0.4` | `antismash --genefinding-tool none --minimal --cpus N --output-dir out bakta.gbff` | regions, products, core biosynthetic genes |
| Peptide → genome, proteins | miniprot 0.18; BLAST+ (tblastn) | conda `miniprot=0.18`, `blast` (pinned) | §7.2 | alignments |

Tool gotchas the wrappers MUST handle, each with a test:

- **TerminatorNet.**
  - Its argument parser matches flags by substring. Pass only the documented flags, use staged file names without spaces, and always pass `-s`.
  - Output lines are `name, TERM_i, A_tract, hairpin, U_tract, structure, GH, pred, start, stop[, strand]`, 1-based, with minus-strand coordinates already converted to the forward strand.
  - Skip lines reading `NO TERMINATORS FOUND`.
- **Promoter Calculator.** It ignores FASTA headers and joins all lines into one sequence, so pass exactly one region per call.
- **Promotech.** `-s` expects exactly 40-nt sequences. Window headers encode `contig|strand|end3p` (avoid characters the tool may mangle; test it). Genome-wide mode writes about 5.8 GB of intermediate files per strand for E. coli, so it is opt-in and preflight warns.
- **miniprot.** `-T 11` is an indexing option. Build the index with `miniprot -t N -T 11 -d genome.mpi genome.fa`, then align with `miniprot -t N -S -N 5 --outs=0.95 --gff genome.mpi proteins.faa`. miniprot indexes ORFs ≥ 30 aa only, which is why short and divergent proteins fall to tblastn.
- **OSTIR.** `start_position` is 1-based in the input sequence. Restrict to the annotated start with CSV `start`/`end` columns.
- **Bakta.** `--regions` accepts CDS only, with at most 30 bp overlap with de novo CDS. Parse the JSON defensively: the contig field may be `sequence` or `contig`, and pseudogenes may be flagged CDS. A contract test on recorded 1.12.2 JSON fixes the field names.

#### 7.1.2 Gene model and replicons

- **`replicons.tsv` for Bakta** is built by preflight:
  - the largest contig is the chromosome if it is circular or ≥ 1 Mb;
  - other circular contigs < 500 kb are plasmids;
  - everything else is a contig.
- **Parsing.** `parse_bakta` produces `genes.parquet`, one row per CDS, sORF or pseudogene, with: locus tag, gene, product, type, interval, strand, `partial`, `crosses_origin`, `hypothetical`, amino-acid sequence, nucleotide sequence, start codon and `db_xrefs`. Non-coding features go to `features.parquet`.
- **Ribosomal-protein genes** (the CAI reference) are those whose product matches `^(30S|50S) ribosomal protein` or whose gene name matches `^(rps|rpl|rpm)[A-Z]`.
- **Anti-Shine–Dalgarno sequence.** Take the last 9 nt of each 16S rRNA (the most common 9-mer across copies). If no 16S exists, or the 9-mer lacks `CCTCC`, use the E. coli default `ACCTCCTTA` and flag `anti_sd_default`.
- **tRNA gene copy numbers** are counted per anticodon from Bakta tRNA features, for tAI.

#### 7.1.3 Operon pairs and transcription units (TUs)

- **UniOP adapter.** Write the Bakta CDS and sORFs as a Prodigal-style GFF in genome order, keep a map back to locus tags, run UniOP, and parse `uniop.pred`.
- **Fallback.** When UniOP fails, the distance model gives p(operon | d) = 1 / (1 + exp((d − d₀)/s)), where d is the intergenic distance in bp (negative for overlaps). Defaults are d₀ = 50, s = 20; both are refitted on RegulonDB K-12 pairs during calibration (§15.4) and stored in `src/gippa/data/calibration.yaml`. The manifest records `operon_source = distance_model`.
- **TUs** are maximal runs of same-strand adjacent genes linked by p ≥ `uniop.pair_threshold`.
- **TU upstream region.** From −300 to +50 bp around the first start codon, strand-aware, clipped at the end of the nearest upstream gene on either strand.

#### 7.1.4 Promoters

- **Targeted Promotech.** For each start codon S of every CDS and sORF (and, in the second round, every novel-ORF candidate start), score the 40-nt windows on the gene's strand whose 3′ end lies at S + k for k ∈ [−200, +20]: 221 windows per start. Deduplicate identical windows, chunk them into files of 100,000 windows and run the chunks in parallel.
- **Per-gene features:**
  - `promoter_score`: the maximum score over the gene's own windows;
  - `promoter_score_tu`: the same at the TU's first gene;
  - the count of windows ≥ 0.5 between the TU start and the gene.
- **Promoter between two adjacent genes:** the maximum score over windows of the downstream gene whose 3′ end lies after the upstream gene's end.
- **`promoters.bed`** holds windows with score ≥ 0.5, merged when overlapping; the score is the maximum.
- **Promoter Calculator.** Run once per TU upstream region. Keep promoters on the gene's strand whose TSS lies in the region. Features: `tx_rate_max`, `tx_rate_sum`, `tss_distance_bp`.

#### 7.1.5 Terminators

- **TerminatorNet** (default) runs on the whole genome with `-s`.
- **Terminator between two same-strand adjacent genes:** a terminator on that strand whose interval lies within [upstream gene stop − 50 bp, downstream gene start], in transcription direction, with confidence ≥ `terminatornet_yes`.
- **Calibration.** Until RegulonDB calibration (§15.4) the cut-off is a placeholder flagged `uncalibrated` in the manifest.
- **TransTermHP cross-check.** Report the share of TerminatorNet calls within ±20 bp of a TransTermHP call in `annotation/terminators_crosscheck.tsv`.
- **`--terminators transtermhp`** makes TransTermHP the main tool, with confidence = CONF/100.
- **BacTermFinder calls** (opt-in) are added to `terminators.bed` with `source=BacTermFinder` but are not used in Module 4 flags in version 1.

#### 7.1.6 RBS and codon usage

- **OSTIR input.** One row per gene: `id` = locus tag; `seq` = −50 to +50 nt around the start codon, strand-aware (clipped at contig ends); `start = end = 51`, or the clipped equivalent; `anti-Shine-Dalgarno` from §7.1.2.
- **Per-gene features:**
  - `tir_ostir` = `expression`;
  - `sd_dg` = `dG_rRNA:mRNA`;
  - `mrna_dg` = `dG_mRNA`;
  - `spacing_dg`, `standby_dg`;
  - `spacer_bp` = `RBS_distance_bp`;
  - `start_codon`.
- **codon-bias.**
  - CAI (reference = ribosomal-protein genes);
  - tAI (tRNA gene copy numbers);
  - ENC;
  - Fop (optimal codons from the reference set);
  - GC3;
  - codon-pair bias (reference = all CDS).
- **Computed directly:**
  - AT content of codons 2–10;
  - `rare_codons_2_10`: codons whose share within their synonymous family is < 10 % genome-wide.

#### 7.1.7 Regions

| Type | Source | Notes |
|---|---|---|
| `prophage` | geNOMAD proviruses | gene membership by ≥ 50 % CDS overlap. Prophage genes get a function subset by keyword rules on the Bakta product and geNOMAD markers: `integrase` (integrase, recombinase, excisionase), `repressor` (repressor, cI, Cro, antirepressor), `structural` (capsid, portal, tail, baseplate, head, terminase, scaffold, tape measure, fiber), `lysis` (holin, lysin, lysozyme, endolysin, spanin, Rz), `cargo` (everything else) |
| `plasmid` | Bakta replicon table, geNOMAD plasmid score | the whole replicon (replicon-level effect, §7.3.5) |
| `bgc` | antiSMASH regions | core biosynthetic genes flagged |
| `crispr` | Bakta CRISPR arrays | with flanking `cas` genes within 1 kb, if any |
| `amr` | genes Bakta annotates via AMRFinderPlus | a gene set rather than a contiguous region |

#### 7.1.8 Outputs of Module 1

- `annotation/genes.gff3`: the merged model with source-tagged layers.
- `annotation/regions.tsv` and `regions.bed`.
- `annotation/promoters.bed`, `annotation/terminators.bed`.
- `annotation/operon_pairs.tsv`, `annotation/tus.tsv`.
- `annotation/id_map.tsv`, when `--annotation` is a GFF3 or GBFF.
- `annotation/terminators_crosscheck.tsv`.
- `annotation/tools/<tool>/`: native outputs, untouched.

**Tests (Module 1):**

- a contract test per tool on recorded native output (§11.4);
- an adapter round-trip for UniOP (Bakta → Prodigal-style GFF → back to locus tags);
- window generation property tests: count, strand, 3′-end range, origin wrap on circular replicons, no wrap on linear contigs;
- a Promoter Calculator one-region-per-call test;
- a TerminatorNet `-s` and coordinate test;
- anti-SD extraction, tAI counts, ribosomal-protein regex;
- the distance-model fallback;
- region membership at the 50 % boundary.

### 7.2 Module 2: mapping proteins and peptides to the genome

Each sequence takes the cheapest route that can place it; only failures move on to the next step.

#### 7.2.1 Six-frame index

- **Translation.** Translate all six frames of every contig with table 11 (`*` at TAA, TAG, TGA; ambiguous codons → `X`, which never matches).
- **Alphabet.** Apply I → L in the index and in queries, so the alphabet has 19 letters plus `*` and `X`.
- **Circular replicons.** Append the first 120 nt before translation, so peptides up to 40 aa can span the origin. Map hits back modulo the replicon length and deduplicate them.
- **Storage.** Arrays of `uint8` per frame, with a function that maps (contig, frame, amino-acid index) to a genomic interval and strand.
- **Decoy translations:**
  - *reverse* (database input): each frame's amino-acid string reversed;
  - *shuffle* (de novo input): residues shuffled within each stop-free segment, in windows of 50 aa, with a seeded RNG.
  - Decoy hits are only ever counted, never placed.

#### 7.2.2 Route 1: exact peptides (7–40 aa)

Build an Aho–Corasick automaton (`pyahocorasick`) of all unique query peptides (I → L). Scan every target and decoy frame in one pass. Keep every locus. Peptides under 7 aa were dropped in preflight: a 6-mer has about one chance hit in a 5-Mb genome.

#### 7.2.3 Route 2: proteins > 40 aa (miniprot), then tblastn

- **Alignment.** miniprot as in §7.1.1. Parse GFF attributes `Identity`, `Positive`, `Rank`, `Target`, `Frameshift` and `StopCodon`. Keep the flags `fs` and `st`.
- **Assignment to a Bakta CDS** requires all of:
  - identity ≥ 0.80;
  - query coverage ≥ 0.80;
  - reciprocal overlap with the CDS ≥ 0.80 (overlap / CDS length and overlap / alignment span).
- **Similarity class:** `identical` (identity 1.0 and coverage 1.0), `near` (≥ 0.95) or `divergent` (0.80–0.95). Below the thresholds a protein goes to the next route.
- **Multiple loci.** Secondary alignments (Rank > 1) are kept as additional loci.
- **tblastn leftovers.** Proteins with no assignment go to `tblastn -db_gencode 11 -seg no -evalue 1e-5`, with the same thresholds after merging collinear HSPs per locus.

#### 7.2.4 Route 3: tolerant matching (leftover peptides ≥ 10 aa)

- **Difference budget k:** 1 at 10–11 aa, 2 at ≥ 12 aa.
- **Each of these costs one difference:**
  - a substitution;
  - an adjacent swap (XY ↔ YX);
  - a mass-equivalent block exchange between the residue compositions {N} ↔ {G,G}, {Q} ↔ {A,G}, {G,T} ↔ {A,S}, {A,V} ↔ {G,L} or {A,N} ↔ {G,Q} (both orders of each pair);
  - a near-equal exchange: {W} ↔ {A,D} or {E,G}, {R} ↔ {G,V};
  - a ±1-nt frameshift split: prefix in frame f, suffix in frame f ± 1 at the adjacent nucleotide offset, with ≥ 4 aa on each side.
- **I/L** is free.
- **Seed and verify (G4).**
  - Index the six-frame (and decoy) translations by q-grams (q ∈ {3, 4, 5}).
  - For each peptide, choose the largest q for which the q-gram lemma gives a lossless threshold t ≥ 1, where t = (L − q + 1) − k(q + 1); each difference destroys at most q + 1 q-grams.
  - Count shared q-grams per diagonal band (diagonal ± 1, which absorbs 1↔2 exchanges).
  - Verify candidate diagonals with a banded dynamic program (band ± 2) whose only operations are those listed, each costing 1.
  - Frameshift splits: anchor the longer side (exact, or within the remaining budget), then verify the shorter side in the shifted frame.
  - Implement the hot loops in numba.
- **Completeness requirement.** Property tests MUST show zero false negatives against a brute-force matcher on random genomes ≤ 50 kb, for every difference type (§11.2).
- **Output.** For each peptide: all target loci at the minimal cost, the cost and the difference types; and the minimal decoy cost (search decoys identically).

#### 7.2.5 Placement FDR (G3)

**Why decoys are needed.** Expected chance hits of one wrong peptide in about 10⁷ six-frame positions, at E. coli amino-acid frequencies with I = L (per-position match probability s = Σ f_a² ≈ 0.069):

| Peptide length | exact | ≤ 1 substitution | ≤ 2 substitutions |
|---|---|---|---|
| 6 aa (dropped) | 1.1 | 90 | 3,000 |
| 7 aa | 0.08 | 7 | 300 |
| 8 aa | 0.005 | 0.6 | 28 |
| 9 aa | 0.0004 | 0.05 | 2.4 |
| 10 aa | 0.00003 | 0.004 | 0.2 |
| 12 aa | < 10⁻⁶ | 0.00002 | 0.0015 |
| 14 aa | < 10⁻⁶ | < 10⁻⁶ | 0.00001 |

Formula: E(L, m) = N · Σ_{j≤m} C(L, j) (1 − s)^j s^(L−j), with N = number of six-frame positions. Appendix B reproduces the table.

**Procedure.**

1. Each peptide counts once, at its best cost, in the target search (m_T) and in the decoy search (m_D). The strata are m ∈ {0, 1, 2}.
2. **Ranking score.** E_i = expected chance hits for peptide i at cost m: N times the Poisson-binomial probability of ≤ m mismatches, with per-position match probability f_a, the genome's six-frame frequency of the peptide's residue.
3. **Lists.** Within stratum m, take the target list (peptides with m_T = m) and the decoy list (peptides with m_D = m), both ordered by E ascending. At each threshold e:

   ```latex
   \widehat{\mathrm{FDR}}_m(e) = \hat{\pi}_{\mathrm{wrong}} \, \frac{D_m(e) + 1}{T_m(e)}
   ```

   The placement q-value of peptide i is q_i = min over e ≥ E_i of FDR_m(e). Accept when q ≤ `stats.alpha.placement` (0.01).
4. **π̂_wrong:**
   - *m = 0, database input:* the mean PEP of the input peptides if PEP is present; otherwise the maximum input q-value (≤ `input.q_value_max`); otherwise 1.
   - *m = 0, de novo input:* 1 − the lower 95 % bound of the de novo precision for the peptide's score bin (below); 1 if that cannot be estimated.
   - *m ≥ 1:* always 1 (conservative).
5. **Same-strain input.** Peptides from a search against this genome's own proteins (`search_db = same`) skip step 3 for exact hits (`placement_q = 0`, flag `same_db`). Their only ambiguity is multi-locus placement. Tolerant hits still go through steps 1–4.
6. **`mapping/mapping_fdr.tsv`** reports one row per length class × m: T, D, π̂_wrong, the stratum FDR π̂(D + 1)/T and n_accepted. Per-m totals are reported too.

**De novo precision, estimated from the genome.** For score bin b (deciles of `score`) and length class ℓ:

- m_bℓ = share of de novo peptides that match the target translation exactly;
- r_bℓ = share that match the shuffled decoy exactly;
- c_ℓ = share of reference peptides (`--reference-peptides`, q ≤ 0.01) that match the target exactly.

```latex
\widehat{\mathrm{prec}}_b = \sum_{\ell} w_{b\ell} \, \frac{m_{b\ell} - r_{b\ell}}{c_\ell - r_{b\ell}}
```

Here w_bℓ is the share of bin-b peptides in class ℓ. Clip to [0, 1]. Take the CI from 1,000 bootstrap resamples within the bin. Do not report bins where c_ℓ − r_bℓ < 0.05. Without reference peptides, c_ℓ defaults to 0.99 and carries the flag `c_assumed`.

#### 7.2.6 Locus classes, uniqueness and gene assignment

- **Locus class** of each placement relative to Bakta features. When several apply, the first in this priority order wins:
  1. `rna_overlap` (suspicious);
  2. `pseudogene`;
  3. `sorf`;
  4. `cds_inframe`;
  5. `nterm_extension` (same frame, upstream of the annotated start, after the previous in-frame stop);
  6. `cds_outframe`;
  7. `antisense`;
  8. `intergenic`.

  A peptide with loci of different classes also gets `multi`.
- **Uniqueness:**
  - `unique`: one locus;
  - `multi_identical`: several loci in genes with identical protein sequences, such as IS copies; counts toward detection of each copy, flagged, never toward abundance;
  - `shared`: different genes; counts toward neither detection nor abundance of a single gene.
- **Peptides → genes.**
  - `cds_inframe` and `sorf` placements assign a peptide to that gene.
  - `nterm_extension` counts toward the gene's detection (flag `nterm_extension_evidence`) but not its abundance, and goes to Module 6.
- **Proteins → genes** through §7.2.3.
  - A cross-strain protein (for example a K-12 accession) carries its abundance to the strain gene it maps to.
  - Two input proteins on one gene → flag `two_proteins_one_gene`. With protein-level abundance only, the gene takes the larger value. With peptide intensities, GIPPA's own gene-level iBAQ and Top3 already avoid double counting.
- **Consistency check.** With same-strain input, a peptide whose engine `protein_ids` disagree with the mapped gene is flagged `engine_protein_mismatch`.

#### 7.2.7 Genome relation

- **`related`** (default): the genome may come from a different isolate than the sample. Accepted tolerant hits count toward detection and abundance and are labelled `strain_difference`.
- **`same`**: the genome is the sample's own assembly. Tolerant hits are reported but only exact hits count. Differences are treated as candidate variants or sequencing errors (Module 5).

#### 7.2.8 Unmapped bucket and contaminants

Sequences that land nowhere are matched against the universal contaminant library (Frankenfield 2022, 381 entries), downloaded by `gippa setup --contaminants` (vendor it only if its licence allows):

- peptides by exact substring (I = L);
- proteins by `blastp` at ≥ 90 % identity over ≥ 50 % of the query;
- engine protein IDs carrying a contaminant prefix also count as contaminants.

Matches are labelled `contaminant`. Everything else is `unmapped`, with a reason: `no_hit`, `rejected_placement_fdr`, `below_identity` or `too_short`. Plass is not used. A taxonomy lookup comes later.

#### 7.2.9 Outputs of Module 2

`mapping/peptide_hits.tsv`, `protein_hits.tsv`, `mappings.probed` (proBed 1.0, BED12+13), `mappings.bed`, `mappings.gff3`, `unmapped.tsv` and `mapping_fdr.tsv`. Schemas are in §9.3.

**Tests (Module 2):**

- **Unit:** table-11 translation; frame ↔ coordinate maps on both strands; I = L; origin-spanning peptides on circular contigs; no wrap on linear contigs; Aho–Corasick multi-locus.
- **Property:** a peptide cut from a random genome maps back to its origin; reverse complement swaps the strand; tolerant matcher = brute force; decoy translations keep composition.
- **Calibration:** placement FDR ≤ q + Monte Carlo error on simulated wrong peptides (database and de novo error models); the chance-hit table is reproduced by simulation within 20 %.
- **Golden:** the hit tables for the 50-kb fixture.
- **Gate M3:** ≥ 99.9 % recall of planted peptides; mapping FDR calibrated.

---
### 7.3 Module 3: gene table, detection, detectability, abundance, context statistics

All statistics run on one gene table (one row per annotated CDS, sORF or pseudogene) and compare what was detected against what each gene's detectability predicts. **Eligible genes** for detection statistics are the non-partial protein-coding rows. Outputs say "detected", never "expressed".

#### 7.3.1 Accepted evidence

A peptide is **accepted evidence** for a gene when all of these hold:

- it passes the engine threshold (`q_value ≤ input.q_value_max` when a q-value exists);
- its placement is accepted (`placement_q ≤ 0.01`, or `same_db`);
- its locus class assigns it to the gene (§7.2.6);
- its uniqueness is `unique` or `multi_identical`;
- under `genome_relation = same`, it is an exact hit.

**Distinct peptides** are distinct stripped sequences. Charge states, runs and modified forms of one sequence count once.

#### 7.3.2 Detection tiers (decision 4)

The tier is chosen once per run from the input and recorded in the manifest and in `genes.tsv` (`detection_tier`).

| Tier | When | Rule | Gene error |
|---|---|---|---|
| `decoy_picked` | decoy rows with `protein_ids` exist and ≥ 80 % of them pair with a gene | picked target–decoy over genes (below) | gene q-value |
| `mayu` | otherwise, q-values or PEPs exist | MAYU-style estimate (below) | gene q-value |
| `two_peptide` | otherwise, peptide sequences only | detected = ≥ 2 distinct accepted peptides; single-peptide genes get `detection_status = single_peptide_unestimated` | `unestimated` |
| `protein_list` | proteins only, no peptides | detected = a mapped input protein | `unestimated` (inherited from the engine) |

`detected_2pep` (≥ 2 distinct accepted peptides) is always reported, whatever the tier.

For downstream statistics, single-peptide genes under `two_peptide` count as not detected in the primary analysis and as detected in a sensitivity rerun.

**Picked gene-level FDR** (Savitski 2015, adapted from proteins to genes):

1. **Target score.** s_g is the best accepted target peptide of gene g, by −log₁₀ q (ties broken by lower PEP, then higher `score`).
2. **Decoy score.** Strip the decoy prefix from a decoy row's protein ID and map the resulting target protein to its gene. The decoy score s̃_g is the best decoy peptide paired with gene g. Unpaired decoys are counted and reported; if more than 20 % of decoys are unpaired, fall back to the `mayu` tier with a warning.
3. **Picking.** Per gene, keep max(s_g, s̃_g); a tie counts as a decoy.
4. **FDR and q.** FDR(t) = (#picked decoys ≥ t + 1) / #picked targets ≥ t, made monotone into gene q-values. A gene is detected at gene q ≤ `stats.alpha.identification` (0.01).

**MAYU-style closed form (G7)** (Reiter 2009):

- **Weights.** Each gene's weight a_g is `n_observable_peptides` (floor 1). Genes are grouped into `mayu_size_bins` bins (deciles of a_g). For a single false peptide, the per-gene probability of being hit is w_b = mean(a_g in bin b) / Σ_h a_h.
- **Thresholds.** For a threshold t on the best peptide q:
  - P(t) = the accepted peptides with q ≤ t;
  - F(t) = Σ PEP over P(t) when PEPs exist, otherwise t · |P(t)|;
  - n_b(t) = genes in bin b with ≥ 1 peptide in P(t); N_b = eligible genes in bin b.
- **Per bin.** p_b = 1 − (1 − w_b)^F(t) and

  ```latex
  \mathrm{FP}_b(t) = \frac{(N_b - n_b(t))\, p_b}{1 - p_b}
  \qquad
  \widehat{\mathrm{FDR}}(t) = \frac{\sum_b \mathrm{FP}_b(t)}{\sum_b n_b(t)}
  ```

  This is the solution of FP_b = (N_b − n_b + FP_b) · p_b: false genes are hit only among genes without a true peptide.
- **Gene q.** q_g = min over t ≥ best_q_g of FDR(t). Detected at q ≤ 0.01.
- **≥ 2-peptide variant.** Repeat the computation with p_b⁽²⁾ = 1 − e^(−λ)(1 + λ), λ = F · w_b. It is reported in `detection_summary.tsv` only.

#### 7.3.3 Detectability model (cross-fitted)

- **Response.** `detected` for eligible genes.
- **Covariates (physicochemical only):**
  - log1p(`n_observable_peptides`): unique tryptic 7–30 aa under the configured enzyme;
  - log(`length_aa`);
  - GRAVY;
  - pI (Biopython ProtParam);
  - `tm_segments`: Kyte–Doolittle, window 19, mean ≥ 1.6, runs merged, capped at 10 (G6).

  Expression determinants (CAI, promoter, RBS) are excluded so Module 8 is not circular. No tested category is ever a covariate: its observed-minus-expected would be zero by construction.
- **Pseudogenes.** Covariates come from Bakta's amino-acid sequence when present, otherwise from the in-frame translation up to the first stop.
- **Model (G5).**
  - Per covariate: `SplineTransformer(n_knots=5, degree=3, knots="quantile")`.
  - Then: `LogisticRegression(penalty="l2")`, with C chosen from {0.01, 0.1, 1, 10} by inner `GroupKFold(3)`.
- **Cross-fitting.** `GroupKFold(5)` with groups = TU (folds from seeded TU order). Each gene's `p_detect` comes from the model that did not see its TU.
- **Calibration check.** 10 equal-count bins, Brier score, and the calibration intercept and slope (logistic regression of d on logit p̂), written to `detection_summary.tsv` (section `detectability`) and plot V18.
- **Edge case.** With fewer than 50 detected or 50 undetected eligible genes, `p_detect` = the overall rate and the flag is `detectability_constant`.

#### 7.3.4 Abundance (decision 3)

- **With peptide intensities.** Use the accepted peptides of gene g with uniqueness `unique` (plus tolerant hits under `related`), with intensity aggregated as in §5.3.
  - iBAQ_g = Σ intensity / max(1, `n_observable_peptides_g`), where observable peptides are counted on the strain gene sequence with the configured enzyme, 0 missed cleavages and length 7–30.
  - Top3_g = the mean of the 3 most intense peptides (with fewer than 3, the mean of those present, flag `top3_lt3`).
  - `abundance_raw` = iBAQ by default (`abundance.method`); `abundance_scale = molar`.
- **Protein-level input only.** `abundance.tsv` joined through the protein → gene map. `abundance_type` = the declared type. Scale: `molar` for `ibaq` and `top3`; `relative` for `maxlfq`, `intensity` and `unknown`.
  - With a relative scale, Module 4's abundance tests are marked `unestimated (relative_scale)`, and stoichiometry checks are skipped.
- **Transform.** Linear values → log2; zeros → missing.
- **Missing values for detected genes** (`abundance.missing`):

  | Policy | Fill value |
  |---|---|
  | `median` (default) | median of observed log2 values |
  | `mean` | mean of observed log2 values |
  | `min` | minimum of observed log2 values |
  | `zero` | minimum of observed log2 values − 1 (zero cannot be represented on a log scale; documented) |
  | `exclude` | none: stays detected but leaves the abundance statistics |

  Filled values set `abundance_imputed = true`. Undetected genes are never imputed.
- **Rank-based where possible.** Single-sample analyses use ranks wherever they can: Spearman, Hodges–Lehmann, normal scores.

#### 7.3.5 Context statistics

**(a) Basic detection summary** (`detection_summary.tsv`, section `basic`). For each replicon, feature class (`cds`, `sorf`, `pseudogene`, `hypothetical`) and COG category, report:

- n_genes, n_detected, % detected;
- E = Σ p_detect;
- O − E, with a 95 % CI from a TU-block bootstrap (resample the group's TUs with replacement);
- O/E.

**(b) Region tests** (`regions_stats.tsv`), for the types `prophage`, `bgc`, `crispr` and `amr`:

- **Detection statistic.** O − E = Σ_{g∈R} (d_g − p_g) over the region's eligible genes.
- **Abundance statistic.** The Hodges–Lehmann shift: the median of pairwise differences in log2 abundance between the region's detected genes and the replicon's other detected genes.
- **Null: circular gene-index shifts** within the replicon. The region's gene-index set I becomes (I + k) mod N, for k = 1 … N − 1 (exhaustive while N ≤ `exhaustive_shift_max_genes`, otherwise B random shifts). Recompute both statistics at every shift.
  - On linear contigs, only shifts that keep the set inside [0, N) are used.
  - Two-sided p = (1 + #{|S_k − S̄| ≥ |S_obs − S̄|}) / (1 + K), where S̄ is the null mean and K the number of shifts.
- **CI.** TU-block bootstrap: resample the region's TUs and, for the HL shift, the outside TUs too.
- **Multiple testing.** BH within each region type. Every row carries `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`, `null_model = circular_shift` and `n_perm = K`.
- **Small regions.** Regions with < 3 eligible genes are flagged `small_region`.
- **Plasmids are replicon-level effects.** Report O − E and the HL shift of the plasmid against the chromosome, with TU-block bootstrap CIs. Set `p_value = NA` and `null_model = none`, and label the result `dosage_plus_regulation`. With `--replicon-depth`, also report the shift minus log2(plasmid depth / chromosome depth).
- **Prophage function subsets.** Per prophage, report detected / total for each subset of §7.1.7 (rows with `subset` ≠ `all`). These rows are descriptive and untested in version 1.

**(c) Autocorrelation** (`autocorrelation.tsv`):

- **Values.** x = log2-abundance residuals of detected genes, after removing the replicon mean and a wide LOESS trend along gene index (span = 10 % of genes, circular padding on circular replicons).
- **Moran's I at lag k = 1 … 20.** Pairs (i, i + k) along the gene index (circular where applicable) where both genes have values:

  ```latex
  I_k = \frac{\tfrac{1}{W_k}\sum_{(i,j):\,|i-j|=k} z_i z_j}{\tfrac{1}{n}\sum_i z_i^2}, \qquad z = x - \bar{x}
  ```

  W_k is the number of such pairs and n the number of genes with values.
- **Join-count for detection.** BB_k = the number of lag-k pairs of eligible genes where both are detected.
- **Two nulls**, B permutations each:
  - *gene label:* shuffle values among the replicon's genes; asks whether there is any structure;
  - *TU block:* shuffle whole TUs along the replicon, keeping internal order and orientation; asks whether there is structure beyond operons.
- **Reporting.** Two-sided p per lag and null (centred on the null mean); BH across lags within each null; 2.5–97.5 % envelopes for V10.

**(d) Islands** (`islands.tsv`):

1. **Gene scores.** Within each replicon, rank eligible genes by log2 abundance, with undetected genes tied at the lowest rank (left-censored). Convert to Blom normal scores Φ⁻¹((r − 3/8)/(n + 1/4)), using average ranks for ties.
2. **TU score** z_t = the mean of its genes' scores. One operon therefore cannot form an island alone; the minimum window is 2 TUs.
3. **Detrend** with LOESS along TU order (span 10 % of TUs, circular padding on circular replicons).
4. **Scan** all windows of 2–20 consecutive TUs (wrapping on circular replicons): S_w = Σ z_t / √n_w, for both high (S) and low (−S) islands.
5. **Null.** A circular block bootstrap of the detrended TU scores, with block length ℓ = the smallest lag at which TU-level Moran's I is not significant at 0.05 (gene-label null), bounded to [1, 20]. For each replicate b, M_b = max over windows of |S_w|.
6. **FWER.** p_w = (1 + #{M_b ≥ |S_w|}) / (1 + B). Islands are windows with p_w ≤ 0.05; overlapping significant windows with the same direction merge into one island, represented by its best window.
7. **Report per island:**
   - ID, direction, TU and gene span, coordinates;
   - mean log2 abundance, detected fraction;
   - S, `p_value` (FWER), with `q_value = NA` and `method = scan_max_stat`;
   - effect = mean detrended score, with a TU-bootstrap CI.

**(e) Origin-to-terminus gradient** (rows in `autocorrelation.tsv`, `statistic = oriter_slope`):

- **Position.** oriC comes from Bakta; the terminus is taken as oriC + L/2. Replichore position = circular distance from oriC / (L/2), in [0, 1].
- **Regression.** TU-level log2 abundance on replichore position, per replichore and pooled, with a circular block-bootstrap CI.
- **Descriptive only:** `p_value = NA` and the flag `descriptive`, because evolution places transcription and translation genes near oriC. Without an oriC call the analysis is skipped with the flag `no_oric`.

#### 7.3.6 Tests and gates (Module 3)

- **Unit:**
  - the picked-FDR toy example;
  - the MAYU closed form (Appendix B.3);
  - Blom scores with ties;
  - Moran's I on an alternating series (−1);
  - HL shift;
  - exhaustive shift enumeration on circular and linear contigs;
  - scan window enumeration with wrap-around;
  - the missing-value policies.
- **Property:**
  - O − E is unchanged by relabelling genes outside the region;
  - circular shifts preserve region size;
  - TU-block permutation preserves TU sizes and internal order.
- **Calibration:**
  - gene-level FDR ≤ q + Monte Carlo error under simulated identification errors, for both the picked and MAYU tiers;
  - the detectability calibration slope within 1 ± 0.2;
  - region-test p-values uniform under autocorrelated nulls;
  - island FWER ≤ 5 % under AR(1) backgrounds without islands;
  - autocorrelation tests at level under the TU-block null.
- **Power** (recorded, not gated): islands of 5, 10 and 20 TUs with effects of 1, 1.5 and 2 SD; regions with detection shifts of 10–50 %.
- **Gates:** M4 = circular-shift tests calibrated; M5 = island FWER ≤ 5 % and the power curve recorded.

---
### 7.4 Module 4: operons

A single sample measures **co-abundance, not co-expression**. Operon partners can legitimately differ 1–10-fold; the ATP synthase subunits, for example, come from one mRNA at very different rates. Calls therefore compare against empirical distributions, with built-in decoys.

#### 7.4.1 Pair table

Every adjacent pair of protein-coding genes on a replicon, in genome order, gets one of four classes. RNA genes between the two are noted with the flag `rna_between`.

| Class | Definition | Role |
|---|---|---|
| OP | same strand, UniOP p ≥ threshold | operon pairs; define the within-operon distribution |
| NO | same strand, UniOP p < threshold | candidates for missed operons |
| CV | convergent (→ ←) | decoys: no shared promoter or transcript |
| DV | divergent (← →) | reported separately (shared regulatory region) |

There is also a second decoy class, **SW** ("sandwich", → ← →): genes i and i + 2 on the same strand, with gene i + 1 on the opposite strand. They cannot share an operon. Their distance is the sum of the two gaps (G16).

**Orientation.** For OP, NO and SW pairs, gene *a* is upstream in the direction of transcription and Δ = log2 abundance(b) − log2 abundance(a), that is, downstream minus upstream. For CV and DV pairs, *a* is the left gene in genome order.

**Pair columns:**

- intergenic distance (bp; negative when genes overlap), UniOP p;
- `promoter_between` (maximum Promotech score, §7.1.4), `tx_rate_between`;
- `terminator_between` (maximum confidence, §7.1.5);
- both detection states and both `p_detect` values;
- both log2 abundances and Δ;
- test results (below).

#### 7.4.2 Within-operon distribution

Fit a Student-t (location, scale and df, with df in [2, 50]) by maximum likelihood to Δ over OP pairs where both genes are detected, starting from the median and MAD. Pairs with imputed abundance are excluded from the fit.

#### 7.4.3 Expected co-expressed, but not: discordant OP pairs

- **Both detected.** Two-sided p = 2 · min(F_t(Δ), 1 − F_t(Δ)) under the fitted t.
- **One detected, one not.** Uses a **detection-limit model** of P(detect | abundance, covariates):
  - *Peptide model* (default when peptides exist): logistic regression of observed vs. not observed over the observable peptides of detected genes. Features: the gene's log2 abundance, peptide length and peptide GRAVY. A gene at abundance a is then detected with probability P_det(a) = P(≥ k_min distinct peptides seen), where k_min = 1, or 2 under `two_peptide`; this is a Poisson-binomial over the gene's observable peptides.
  - *Replicate model* (when ≥ 2 runs exist): logistic regression of per-run detection on the gene's mean log2 abundance plus the detectability covariates, fitted on genes seen in some runs but not others.
  - **P(miss)** for the undetected partner u, given the detected partner's abundance a_v, is the mean over the empirical within-operon log-ratios r of [1 − P_det,u(a_v + r)]. The r values are oriented: Δ when u lies downstream, −Δ when upstream.
  - p = P(miss). A small value means the missing partner should have been seen.
  - With one run and no peptides, this test is skipped and marked `unestimated (no_detection_limit_model)`.
- **Multiple testing.** BH across all tested OP pairs (both kinds together). Pairs are discordant at q ≤ `stats.alpha.operon_flags` (0.1).
- **Relative-scale abundance** (§7.3.4): abundance tests are marked `unestimated (relative_scale)`.

#### 7.4.4 Not expected, but co-abundant: missed operons

- **Score.** s = −|Δ| for NO pairs where both genes are detected.
- **Decoys.** CV pairs and SW pairs, both genes detected, each reweighted to the NO pairs' distribution over intergenic-distance bins b (`operons.distance_bins_bp`). With π₀ = 1:

  ```latex
  \widehat{\mathrm{FDR}}(t) = \frac{\pi_0 \sum_b N_{\mathrm{NO},b} \, \dfrac{\#\{\mathrm{decoy}_b : s \ge t\} + 1}{N_{\mathrm{decoy},b} + 1}}{\#\{\mathrm{NO} : s \ge t\}}
  ```

- **Which decoy.** Compute the FDR with each decoy set and report the **larger** one, because convergent partners can look less alike than same-strand neighbours (transcriptional interference). Make the FDR monotone into q-values, and flag pairs at q ≤ 0.1.
- **Local fdr** for ranking: lfdr(s) = min(1, π₀ f₀(s) / f_NO(s)). Here f₀ is the distance-reweighted decoy density of |Δ| (Gaussian KDE, boundary reflection at 0, Silverman bandwidth) and f_NO the observed density of NO pairs.

#### 7.4.5 Promoter and terminator evidence decides the label

"Yes" means a Promotech score ≥ `annotation.promoters.promotech_yes` (0.5), or a TerminatorNet terminator at or above the calibrated confidence cut-off, between the two genes.

| UniOP | Promoter between | Terminator between | Abundance pattern | Flag |
|---|---|---|---|---|
| operon | yes | no | downstream higher | `internal_promoter_supported` |
| operon | no | yes | upstream higher | `internal_terminator_supported` |
| operon | no | no | discordant | `unexplained_discordance` |
| operon | any | yes | concordant | `conflict_terminator_in_operon` |
| non-operon | no | no | co-abundant | `novel_operon_candidate` |
| non-operon | any | yes | co-abundant | `coabundant_despite_terminator` |

How the abundance patterns are read:

- **"downstream higher" / "upstream higher":** discordant (q ≤ 0.1) with Δ > 0 or Δ < 0. For one-detected pairs, the detected gene is the higher one.
- **"discordant":** q ≤ 0.1 in either direction.
- **"concordant":** an OP pair with both genes detected and q > 0.1.
- **"co-abundant":** an NO pair flagged at q ≤ 0.1.

Any other combination for a discordant or co-abundant pair is labelled `other`.

Promotech and Promoter Calculator see σ70 promoters only, and TerminatorNet and TransTermHP see intrinsic terminators only. A missing call therefore never proves a missing element. The report states this next to the table.

#### 7.4.6 Truth check

On K-12 MG1655, UniOP pair probabilities are scored against RegulonDB v12 transcription units: AUROC, and the cut-off that reaches about 90 % precision becomes the default (§15.4). Comparing UniOP alone with UniOP plus co-abundance needs K-12 proteomics, so it belongs to the deferred full validation.

#### 7.4.7 Outputs and tests

- **Outputs:**
  - `tables/pairs.tsv`: all classes, SW included;
  - `operons.tsv`: per TU: genes, size, detected fraction, median log2 abundance, maximum |Δ| between adjacent members, flag counts;
  - `operon_flags.tsv`: flagged pairs only, with evidence columns.
- **Unit tests:**
  - class assignment, including SW and `rna_between`;
  - orientation of Δ on both strands;
  - the FDR formula on toy counts with empty bins;
  - the Student-t fit;
  - Poisson-binomial P_det;
  - P(miss) on a toy distribution;
  - every row of the flag table.
- **Calibration:**
  - missed-operon FDR ≤ 0.1 + Monte Carlo error on simulations with planted co-abundant NO pairs;
  - the same on the attack null (autocorrelated abundance with no operon effect);
  - discordance p-values uniform for OP pairs simulated without discordance.
- **Gate M6:** decoy FDR calibrated on simulations; UniOP and TerminatorNet cut-offs recalibrated on RegulonDB (§15.4); UniOP-alone AUROC recorded.

### 7.5 Module 5: variants, REF or ALT at each covered codon

At every genome codon covered by accepted peptides, GIPPA reports whether the peptides support the genome's residue (REF), a different residue (ALT), or both. It never names modifications.

#### 7.5.1 What a discrepancy means, by peptide source

| Peptide source | What a discrepancy usually means | Treatment |
|---|---|---|
| Search against this strain's Bakta proteins (`search_db = same`) | cannot occur, except at Sec sites and non-ATG starts, which appear as ALT calls | consistency check: ALT sites without the `start_codon` or `stop_codon` flag are counted as `id_mixups` in `detection_summary.tsv` |
| Search against another strain (e.g., K-12) | the database residue is absent from this genome: a false identification, or genome and sample do not match | with `--genome-relation same`: no ALT calls from these peptides; report the rate of genome-inconsistent peptides as QC in `detection_summary.tsv`. With `related`: ALT calls flagged `database_derived` |
| De novo sequencing | sequencing error, modification, or true variant | full classification below |
| Variant-aware or open search (future) | variant or modification | full classification below |

With `related`, any ALT may simply be a strain difference: such hits count toward detection and abundance and carry `strain_difference`.

#### 7.5.2 Pileup and calls

1. **Residue alignment.** The mapper's traceback gives a residue-level correspondence between each peptide and the genome translation. Mass-equivalent block exchanges are recorded at the first codon of the block, with `alt_type = block` (for example `GG>N`). Frameshift splits are recorded as `alt_type = frameshift` events at the junction; they appear in `variants.tsv` but not in the VCF in version 1.
2. **Pileup.** For each covered codon, count distinct peptides supporting REF (exact residue identity), each ALT residue, and `IL` (the residues differ only by I versus L).
3. **Call:**
   - `REF`: only REF support;
   - `ALT`: only ALT support;
   - `mixed`: both;
   - `no_call`: only `IL` support (mass spectrometry cannot tell I from L).

   `IL` support is ignored when REF or ALT support exists.
4. **Confidence.** An ALT (or the ALT part of a mixed call) is `high` with ≥ 2 distinct peptides supporting the same ALT residue, and `candidate` with 1. Every supporting peptide needs placement q ≤ 0.01.

#### 7.5.3 Flags (they never change the call)

- **`mass_equivalent`.** An equal- or near-equal-mass reading explains the ALT: \|m(ALT) − (m(REF) + Δmod)\| ≤ `variants.mass_equivalent_tol_da`, for Δmod in {0} ∪ a curated list of modifications. Blocks are included. The list is generated by `src/gippa/data/build_mass_equivalents.py` from the residue and modification masses in Appendix B.5 into `mass_equivalents.tsv`. Expected hits include:
  - N vs GG, Q vs AG, K vs Q (0.036 Da), W vs AD/EG, R vs GV;
  - N → D and Q → E (deamidation, +0.984016);
  - A → S and F → Y (+15.994915);
  - G → A, S → T, V → L, D → E and N → Q (methylation, +14.01565);
  - M(oxidised) vs F (0.033 Da).
- **`snv` / `multi_nt`.** Whether one nucleotide change of the genome codon (table 11) encodes the ALT residue.
- **`strain_difference`.** Any ALT under `related`.
- **`database_derived`.** An ALT from database peptides with `search_db = other`.
- **`start_codon`** (G8). The codon is the annotated start of a CDS and the ALT is M.
- **`stop_codon`** (G8). The REF residue is a stop (e.g., Sec read through UGA).

#### 7.5.4 Error control

Decoys control placement only (Module 2). Most de novo errors are near-correct and land at the true locus, where no decoy can see them. The false-ALT rate is therefore measured where truth is known:

- **With paired results.** When `--reference-peptides` contains `spectrum_id`s that pair with de novo rows, every de novo peptide whose spectrum has a confident (q ≤ 0.01) reference identification counts as an error if the two disagree (ignoring I/L).
- **Rates.** Class-specific error rates (by `alt_type` and by `mass_equivalent`) are written as `class_error_rate` and `class_error_n` on each ALT row.
- **Without paired results,** ALT confidence carries `class_error_rate = unestimated`.

#### 7.5.5 Outputs

- `tables/pileup.tsv.gz`: one row per covered codon, with REF, ALT and IL distinct-peptide counts.
- `tables/variants.tsv`: sites with an ALT or `mixed` call, plus frameshift events.
- `tables/variants.vcf`, VCF 4.3:
  - **Header:** `##fileformat=VCFv4.3`; `##source=gippa-<version>`; `##reference=` with the genome SHA-256; `##contig` lines with lengths; `##INFO` and `##FORMAT` definitions.
  - **Records:** CHROM = contig; POS = the lowest 1-based coordinate of the codon. REF = the codon bases on the forward strand. ALT = the minimal-change codon encoding each ALT residue, on the forward strand (reverse-complemented for minus-strand genes), multi-allelic when needed (G14).
  - **INFO:** `GENE`, `AA_POS`, `AA_REF`, `AA_ALT`, `CALL`, `CONF`, `NT_CHANGES`, `ALT_TYPE`, `FLAGS`, `CLASS_ERR`.
  - **FORMAT:** `GT:AD:DP`. GT is haploid: `1`, `2`, … for ALT, `.` for `mixed`. AD = REF then each ALT, in distinct peptides. DP = total distinct peptides.
  - **Sample column:** `config.sample_name` (default: the peptide file stem).

#### 7.5.6 Tests and gate

- **Unit:**
  - distinct-peptide counting (charges and runs collapse);
  - IL → `no_call`;
  - `mixed`;
  - the confidence thresholds;
  - an exhaustive `snv` / `multi_nt` table over 64 codons × 21 residues;
  - mass-equivalent flags for the listed pairs;
  - minus-strand codon conversion;
  - a VCF round trip with `pysam.VariantFile`;
  - the `start_codon` and `stop_codon` flags.
- **Calibration:** classification F1 on simulated de novo errors (substitutions, swaps, blocks, modifications), with variants planted between "sample" and "reference" genomes.
- **Gate M7:** simulated F1 ≥ 0.9 for `high` ALT calls; no `high` ALT calls on in-silico same-strain peptides except at `start_codon` or `stop_codon` sites. A2 itself is checked in the pilot.

### 7.6 Module 6: novel ORFs and start sites

#### 7.6.1 Candidates

A candidate is an accepted placement whose locus class is `nterm_extension`, `intergenic`, `cds_outframe`, `antisense`, `pseudogene` or `sorf`.

Peptides that also map in frame to an annotated gene (`cds_inframe` at any locus), or to a contaminant, never support a novel ORF.

#### 7.6.2 ORF reconstruction

- **Downstream.** Each candidate ORF extends in frame to the next stop codon.
- **Upstream.** Candidate start codons (`annotation.start_codons`) lie between the previous in-frame stop and the first peptide.
- **Two starts are reported:**
  - `start_longest`: the furthest upstream start;
  - `start_best_ostir`: the candidate with the highest OSTIR initiation score.
- **Second-round layers.** Every candidate start gets targeted Promotech windows (−200 to +20) and an OSTIR window (−50 to +50; for several candidates, one window spans them all with `start`/`end` ranges).
- **Grouping.** Peptides are grouped per ORF, and every ORF contains all its peptides.

#### 7.6.3 Classes and the evidence needed to call them

| Class | Where the peptide sits | Extra evidence before calling |
|---|---|---|
| `nterm_extension` | in frame, upstream of the annotated start, after the previous in-frame stop | an upstream in-frame start codon and RBS support |
| `intergenic` | no annotated CDS on either strand | start codon, ≥ 10 codons, RBS support; a homolog raises confidence later |
| `cds_outframe` | same strand, other frame, inside a CDS | ≥ 2 distinct frame-specific peptides |
| `antisense` | opposite strand, inside a CDS | ≥ 2 distinct peptides |
| `pseudogene` | on a Bakta pseudogene | the frameshift or stop junction is checked: `junction_spanned` (a peptide spans it) and `homopolymer_junction` (a homopolymer ≥ 5 nt at the junction). These are labels only; `possible_assembly_error` comes later |
| `sorf` | a Bakta sORF | detection evidence only |

**RBS support** means the OSTIR initiation score at the chosen start is ≥ the `novel_orfs.rbs_min_quantile` (default 0.10) quantile of scores at annotated starts of detected genes. The "downstream start" class from the draft is removed.

#### 7.6.4 Error control

Novel peptides have a lower prior than annotated ones, so the global search FDR understates their error.

- **With decoy rows in the input** (for example after a search against the `makedb` database): a picked target–decoy FDR at ORF level, within each class. Decoys pair through the class-tagged `makedb` IDs (§7.9.2). A call needs class q ≤ 0.01.
- **De novo input:** a subgroup placement FDR. Each decoy hit is assigned the locus class of its mirrored target position (reversed or shuffled windows keep the frame), then FDR_class = π̂_wrong (D_class + 1) / T_class. Calls need q ≤ 0.01.
- **Neither:** a call needs ≥ 2 distinct peptides and carries `class_fdr = unestimated`.

#### 7.6.5 Ranking (G13)

Novel ORFs are sorted by a lexicographic key:

1. distinct peptides, descending;
2. best peptide q, ascending;
3. start codon: ATG < GTG < TTG < other;
4. OSTIR score, descending;
5. promoter score, descending;
6. length, descending.

Conservation comes in a later version.

#### 7.6.6 Outputs and tests

- **Outputs:**
  - `tables/novel_orfs.tsv`;
  - `novel_orfs.gff3`;
  - `novel_orfs.faa` (both start variants);
  - `annotation_patch.gff3`: new ORFs as CDS (IDs `NOVEL_*`), and N-terminal extensions as the existing locus tag with the new start and `Note=nterm_extension`. `gippa makedb --patch` uses it to rebuild the search database for a second-pass search.
- **Tests:**
  - ORF reconstruction on both strands, including circular wrap;
  - each class;
  - the ≥ 2-peptide rule;
  - exclusion of annotated-gene and contaminant peptides;
  - ranking order;
  - GFF3 validity;
  - patch → `makedb` round trip.
- **Calibration:** class FDR on simulations with planted novel ORFs and de novo errors.
- **Gate M7:** F1 ≥ 0.9 for planted novel ORFs with ≥ 2 peptides.

---
### 7.7 Module 7: function (COG, GO, KEGG, GSEA), optional module `function`

Functional tests ask two separate questions:

- which categories are **detected** more often than detectability predicts;
- which are more **abundant** among detected proteins.

#### 7.7.1 Annotation and resources

- **Source.** Bakta full-database `db_xrefs` are the only annotation source (decision 10): `COG:` IDs (and category letters if present), `GO:`, `EC:` and `KEGG:` KOs. Genes without a close match stay unannotated. With the light DB, Module 7 runs on whatever cross-references exist, reports low coverage, and the manifest notes `bakta_db = light`.
- **Resources** (from `gippa setup`):
  - COG 2020 `cog-20.def.tab` (COG → category) and `fun-20.tab` (category names);
  - `go-basic.obo` and `goslim_generic.obo`;
  - KEGG KO → pathway links and pathway names (KEGG REST is free for academic use; the user guide says so).
- **GO propagation.** Terms are propagated to their ancestors with goatools.
- **Coverage.** The share of eligible genes with ≥ 1 term, per ontology, reported on every enrichment row and in `category_summary.tsv`. Coverage leans toward well-studied genes.

#### 7.7.2 Category table (`category_summary.tsv`)

One row per COG category, KEGG pathway and GO slim term, with:

- n_genes, n_detected;
- expected (Σ p_detect), O − E with a TU-bootstrap CI;
- median log2 abundance of detected genes;
- coverage.

Mass fraction is deferred (improvement 11).

#### 7.7.3 Detection enrichment: Wallenius ORA with p_detect weights

- **Setup.** The universe is the eligible genes. Per category C:
  - M = |universe|;
  - n = |C ∩ universe|;
  - N = number of detected genes;
  - k = number detected in C;
  - odds ω = mean(p_detect in C) / mean(p_detect outside C), the goseq weighting applied to detectability.
- **p-values:**
  - over-representation: `scipy.stats.nchypergeom_wallenius(M, n, N, ω).sf(k − 1)`;
  - under-representation: `.cdf(k)`.

  Wallenius is analytical, so there is no permutation floor.
- **Filters.** Category sizes are limited to `function.ora_size` ([5, 500]). BH runs within each ontology and direction.
- **GO redundancy (elim, Alexa 2006).** Test GO terms bottom-up; when a term is significant at the raw-p threshold 0.01, remove its genes from all its ancestors before testing those. Report `p_classic` and `p_value` (elim); BH applies to the elim p-values. Group surviving terms under their parents (`parent_group`).

#### 7.7.4 Abundance GSEA, preranked and TU-collapsed

- **Ranking unit.** Genes in one TU share a transcript, so permuting gene labels would be too liberal. The ranked units are TUs: each TU's score is the mean log2 abundance of its detected genes, and a TU belongs to a set when any of its detected genes does. Set sizes are counted in TUs and limited to `function.gsea_size` ([10, 500]).
- **Method (G15).** fgsea multilevel (Korotkevich 2021), ported to Python and numba from the fgsea source (MIT), with `sampleSize = 101`, `eps = 0`, `nPermSimple = 1000`. It estimates small p-values without a permutation floor.
- **Reported per set:** ES, NES, p, `log2err`, BH q within each collection (COG, KEGG, GO BP, GO MF, GO CC, GO slim), and the leading edge (genes).
- **Redundancy.** A significant GO parent is flagged `redundant` when re-testing it without the genes of its significant descendants gives q > α.
- **Fallback.** If the port cannot meet the golden tolerances (§11.5) within 3 iterations (R9), escalate. The proposed ADR would run R fgsea in its own environment as an external tool; fgsea is MIT, so this does not affect GIPPA's licence.

#### 7.7.5 Outputs and tests

- **Outputs:** `tables/category_summary.tsv`, `enrichment_detection.tsv`, `enrichment_gsea.tsv`.
- **Golden tests** against R on frozen data (§11.5):
  - Wallenius p-values against goseq run with a supplied PWF, relative difference ≤ 10⁻⁶;
  - fgsea ES identical to 10⁻⁹;
  - NES within 0.02;
  - |log₁₀ p − log₁₀ p_R| within 3 × the combined `log2err` / log₂ 10.
- **Unit:** elim on a toy DAG.
- **Calibration:**
  - ORA p-values uniform when detection depends only on p_detect;
  - the attack null "detectability confounded with category membership" (membrane-like genes with low p_detect) stays at level with weights and fails without them;
  - the GSEA attack null "gene sets packed into one operon" stays at level after TU collapse.

### 7.8 Module 8: what explains abundance (promoter, RBS, codon usage)

The question is attribution: for each abundant protein, is the predicted promoter, the RBS or the codon usage the best explanation? A model with grouped cross-validation answers it per gene, and a permutation null says whether any of it beats chance.

#### 7.8.1 Features (`determinants_features.tsv`, one row per eligible gene)

| Block | Features | Source |
|---|---|---|
| Transcription | `promotech_max`; `promotech_tu`; `tx_rate_max_tu`, `tx_rate_sum_tu`; `n_promoters_tu_to_gene`; `tu_position`; `dist_to_tu_start_bp` | Promotech, Promoter Calculator, UniOP |
| Translation initiation | `tir_ostir` (log); `sd_dg`; `spacer_bp`; `start_codon` (one-hot ATG, GTG, TTG, other); `mrna_dg` | OSTIR with the strain's 16S 3′ end |
| Codon usage | `cai` (ribosomal-protein reference); `tai`; `enc`; `fop`; `gc3`; `cpb`; `at_2_10`; `rare_codons_2_10` | codon-bias |
| Protein | `length_aa` (log); `mw`; `pi`; `gravy`; `nterm_residue` (Met removed when residue 2 ∈ {A, C, G, P, S, T, V}); `tm_segments` | sequence |
| Context | `replicon_type`; `replichore_position`; `in_prophage`; `in_bgc`; `gc_local_dev` (10-kb window GC − genome GC) | Modules 1 and 3 |

#### 7.8.2 Analysis

1. **Rows.**
   - Abundance models use detected genes with a non-imputed log2 abundance.
   - Detection AUROC uses all eligible genes.
   - Ribosomal-protein genes, the CAI reference set, are excluded from model fitting and from univariate CAI statistics, and are reported separately.
2. **Univariate.** Spearman ρ of each feature with log2 abundance, and AUROC for detection, each with a TU-block bootstrap CI.
3. **Multivariable.**
   - Models: LightGBM regression (§5.2 parameters) and an elastic-net baseline (standardised features, one-hot categoricals, `l1_ratio` ∈ {0.1, 0.5, 0.9} with inner `GroupKFold(3)`).
   - Cross-validation: `GroupKFold(5)` by TU, repeated 5 times with seeded TU shuffles. Report cross-validated R² as mean ± SD, and keep out-of-fold predictions averaged over repeats.
4. **Null.**
   - Permute the response at TU level, B = 1,000 times: response blocks move among TUs with the same number of included genes. TU sizes that occur only once stay in place, which is conservative.
   - Refit with one CV repeat each time.
   - p = (1 + b)/(1 + B), where b counts permutations with CV R² ≥ the observed mean CV R². Computed separately for each model.
5. **Block attribution.**
   - Drop-one-block ΔR² per repeat, with a t-interval over the repeats.
   - Per-gene SHAP values from LightGBM `pred_contrib=True` on the out-of-fold models, summed by block. No `shap` package is needed.
6. **Per gene.**
   - The dominant block for each top-decile protein = the block with the largest positive SHAP sum.
   - Out-of-fold residuals > 2 SD are flagged `abundance_not_explained_by_sequence`: candidates for regulation or stability.

**Caveats** (printed in the report):

- In a synthetic GFP library, mRNA folding near the RBS explained over half the variance; native genes add regulation and turnover, so expect a lower R².
- All features are predictions, and the data come from one condition.
- CAI, tAI and GC3 are collinear, so interpret blocks rather than single features.

#### 7.8.3 Outputs and tests

- **Outputs:**
  - `tables/determinants_features.tsv`;
  - `determinants_model.tsv`: per model and metric: CV R² mean and SD, null p, n_perm, block ΔR² with CI, univariate rows;
  - `determinants_attribution.tsv`: per-gene SHAP by block, dominant block, residual, flags.
- **Unit:**
  - the Met-processing rule;
  - TU-size-preserving permutation blocks;
  - the p formula;
  - SHAP block sums equal the prediction minus the expected value.
- **Calibration:**
  - null p-values uniform when the simulated response is independent of the features but keeps TU structure;
  - recovery of planted block effects with known R².
- **Gate M8:** null p-values uniform; GSEA and ORA match the golden data.

### 7.9 Utility commands

#### 7.9.1 `gippa setup`

Each item is downloaded inside its own tool environment or container where the tool has a downloader:

- `bakta_db download --type light|full`, which also runs `amrfinder_update`;
- `genomad download-database`;
- antiSMASH `download-databases` (skipped when the container image already contains them);
- the Promotech models zip;
- `go-basic.obo` and `goslim_generic.obo`;
- KEGG KO → pathway links and pathway names;
- COG 2020 tables;
- the contaminant library (at a pinned commit).

`DIR/db_manifest.json` records name, version, URL, checksum, size and retrieval date for each. `--check` only verifies. An unreachable source or a checksum mismatch gives exit code 2, names the host, and leaves the item absent; it never falls back to another source. On FarmShare, `gippa setup` runs as a batch job (`--bakta full` downloads about 84 GB).

#### 7.9.2 `gippa makedb` (convenience; no FDR is built on it)

- **Writes `DB/strain_db.faa`:**
  - Bakta CDS and sORF proteins (`>{locus_tag} {product}`);
  - with `--patch annotation_patch.gff3`, the patched proteins;
  - with `--six-frame-orfs`, class-tagged six-frame ORFs ≥ `--min-orf-aa` (20) that do not overlap annotated CDS in frame: `>GIPPAORF_{class}_{contig}_{start1}_{end1}_{p|m}`, classes `intergenic`, `cds_outframe`, `antisense`, `nterm_extension`;
  - contaminants as `>CONT_{id}` (unless `--no-contaminants`).
- **Also writes** `DB/id_map.tsv` and `DB/makedb_manifest.json`.
- **No decoys.** Search engines generate their own; with DIA-NN, use `--report-decoys` so decoy rows reach GIPPA.
- If no annotation is given, Bakta runs first.

#### 7.9.3 `gippa compare`

- **Inputs.** Two or more run directories on the same genome, checked through the genome SHA-256 in each manifest; a mismatch is exit code 1.
- **Writes to `CMP/`:**
  - `tables/compare_detection.tsv`: a gene × run detection matrix plus the UpSet intersections;
  - `compare_abundance.tsv`: pairwise Spearman ρ on shared detected genes with TU-bootstrap CIs, and Bland–Altman bias and limits;
  - `compare_stability.tsv`: Jaccard overlap of island gene sets, operon-flag pairs, enriched categories (q ≤ 0.05), novel ORFs and ALT sites;
  - `plots/V16_compare.*`;
  - `report/compare.html`.

#### 7.9.4 `gippa report`

Rebuilds `plots/` and `report/` from existing tables without recomputing statistics. Running it twice MUST produce identical `plots/data/*.tsv`.

---
## 8. Statistical framework: decoys, FDR, nulls

Every claim carries an error estimate from a decoy or from a null that keeps the structure of the data. Claims that cannot have one are labelled `unestimated`, never silently trusted.

### 8.1 Claims and their error control

| Claim | Unit | Decoy or null | Error control | Default |
|---|---|---|---|---|
| Peptide or protein identified | precursor | search engine's target–decoy (audited with entrapment in validation) | q-value (inherited) | 0.01 |
| Gene detected | gene | from the input: decoy rows (picked), q-values or PEPs (MAYU-style), or ≥ 2 peptides | gene q-value, or `unestimated` | 0.01 |
| Peptide placed at a locus (exact) | peptide | reversed frames; shuffled frames for de novo | placement q, π̂_wrong (D + 1)/T, ranked by E within the stratum (G3) | 0.01 |
| Peptide placed with ≥ 1 difference | peptide | same decoys and budget, π̂_wrong = 1 | placement q within the mismatch stratum | 0.01 |
| ALT call is real | site | de novo errors on spectra with known database identity | class error rate | reported, else `unestimated` |
| De novo peptide correct | score bin | genome match vs. shuffled-decoy match, by length | precision with bootstrap CI | reported |
| Novel ORF | ORF | decoy rows in the input (picked, per class); genome decoy for de novo | subgroup FDR | 0.01 |
| Region effect | region | circular gene-index shifts within the replicon; plasmids as replicon-level effects | BH within type | 0.05 |
| Autocorrelation | lag | gene-label and TU-block permutation | BH across lags | 0.05 |
| Islands | window | detrended TU scores, circular block bootstrap, maximum statistic | FWER | 0.05 |
| Operon discordance | pair | within-operon t distribution; detection-limit model | BH | 0.1 |
| Missed operon | pair | convergent and → ← → decoys, distance-reweighted | decoy FDR (+ local fdr for ranking) | 0.1 |
| Category detection | category | Wallenius with p_detect weights | BH per ontology and direction | 0.05 |
| GSEA | gene set | TU-collapsed scores, fgsea multilevel | BH per collection | 0.05 |
| Determinant model | model | abundance permuted across TUs of equal size | permutation p | 0.05 |

### 8.2 Rules applied everywhere

- **Decoys mimic targets.** They keep composition, size and structure: reversed or shuffled translations, convergent pairs, shifted regions.
- **Resample TUs, not genes.** Bootstraps and permutations move whole TUs.
- **Adjust for detectability.** Detection tests compare against Σ p_detect from cross-fitted predictions, never against raw counts.
- **Pre-declared test families.** BH applies within the families listed in §8.1, which are fixed in code and documented in `docs/methods/`. All p and q values are reported.
- **Effect sizes with 95 % CIs** accompany every p or q.
- **Permutations.**
  - 1,000 by default and 10,000 with `--publication`, with p = (1 + b)/(1 + B) and seeds recorded.
  - Where many tests share one correction, use analytical or multilevel p-values (Wallenius, fgsea multilevel) or sequential Besag–Clifford permutation (stop after h = 10 exceedances), which avoid the 1/(B + 1) floor.
- **Calibration tests ship with each method** (§8.4) and are CI gates.
- **Sensitivity.** Key results are rerun under other choices (§8.6), and stability is reported as the Jaccard overlap of flagged sets.

### 8.3 Shared statistics engine (`src/gippa/stats/`)

Each function is pure, typed and documented, with its methods note in `docs/methods/`:

| Function | Contract |
|---|---|
| `bh(p) -> q` | Benjamini–Hochberg, monotone, NaN-aware |
| `tda_qvalues(target_scores, decoy_scores, pi=1.0, plus_one=True)` | separate-search q-values, (D + 1)/T |
| `picked_fdr(target_best, decoy_best)` | picked competition, ties go to the decoy |
| `mayu_gene_fdr(best_q, pep, sizes, n_bins)` | closed form of §7.3.2 |
| `perm_p(obs, null, alternative)` | (1 + b)/(1 + B), two-sided centred on the null mean |
| `besag_clifford(stat_fn, obs, h, max_B, rng)` | sequential Monte Carlo p |
| `tu_block_bootstrap(df, tu_col, stat_fn, B, rng)` | percentile CI |
| `circular_shifts(n, idx, circular)` | all valid shifts (§7.3.5) |
| `circular_block_bootstrap(x, block_len, B, rng)` | for scan nulls |
| `scan_max(z, wmin, wmax, circular)` | numba; returns the window table and the maximum |
| `morans_i(x, lag, circular)`, `join_count(d, lag, circular)` | §7.3.5 |
| `wallenius_ora(...)` | §7.7.3 |
| `fgsea_multilevel(...)` | §7.7.4 |
| `poisson_binomial_sf(p, k)` | exact dynamic programming |
| `blom_scores(x, censored_mask)` | ties → average ranks |

### 8.4 Calibration requirements (CI gates)

Pass criteria (constants in `tests/calibration/tolerances.yaml`, a protected file):

- **Level.** For α ∈ {0.01, 0.05, 0.10}, the null rejection rate is ≤ α + 2√(α(1 − α)/R), where R counts independent null replicates.
- **Uniformity** (continuous p only). KS test against U(0, 1), passing at p > 0.001. Discrete permutation p-values use the level check only.
- **FDR.** Mean FDP over R replicates ≤ q + 2 SE(FDP).
- **FWER.** The share of replicates with ≥ 1 false discovery is ≤ 0.05 + 2 SE.
- **Fixed seeds.** Seeds are fixed in the tolerances file, and tests are deterministic.

| Suite | Null or truth generator (§11.3) | Criterion | Per PR (fast) | Nightly (full) |
|---|---|---|---|---|
| Placement FDR | wrong peptides at known rate; database and de novo error models | FDR | R = 50, 100 kb | R = 300, 1 Mb |
| Gene FDR (picked, MAYU) | false peptides on random genes, size-proportional | FDR | 50 | 300 |
| Detectability | logistic truth on physicochemical covariates | calibration slope 1 ± 0.2 | 20 | 100 |
| Region tests | AR(1) abundance, detection from p_detect only | level | 100 | 1,000 |
| Autocorrelation | gene-level noise (label null); TU-structured noise without inter-TU structure (TU null) | level | 100 | 1,000 |
| Islands | AR(1) backgrounds without islands | FWER | 50 | 500 |
| Operon discordance | OP pairs from the t distribution | level | 100 | 1,000 |
| Missed operons | NO pairs independent; planted co-abundant pairs for power | FDR | 50 | 300 |
| ORA | detection ∼ p_detect only | level | 100 | 1,000 |
| GSEA | abundance independent of sets; operon-packed sets | level | 50 | 300 |
| Determinants null | response independent of features, TU structure kept | level | 20 | 100 |
| Novel-ORF class FDR | planted ORFs + false peptides + de novo errors | FDR | 50 | 300 |
| De novo precision | known precision per bin | \|bias\| ≤ 0.05 | 50 | 300 |

The per-PR subset MUST finish in under 10 minutes on a 4-vCPU GitHub runner (and therefore in the 8-CPU session allocation); set R and the genome sizes to fit, and record the timings in `benchmarks/baseline.json`. The full suite runs nightly on FarmShare as a Slurm job array, one suite per task.

### 8.5 Attack nulls (owned by the statistician subagent)

Each attack is built to fool a method; the method passes when it stays at level.

1. Autocorrelated abundance with no islands → islands, autocorrelation.
2. Detectability confounded with region membership (for example, short genes concentrated in prophages) → region tests and ORA, with p_detect weighting.
3. Gene sets packed into one operon → GSEA after TU collapse.
4. High-copy plasmids → region tests, reported as replicon-level effects only.
5. De novo errors that reverse a stretch of sequence → placement FDR with shuffled decoys (a reversed decoy would inflate D).
6. Isobaric swaps everywhere → `mass_equivalent` flags present; no `high` ALT calls from swaps alone at q ≤ 0.01.

New attacks are added as the statistician finds weaknesses. A failing attack blocks the merge.

### 8.6 Sensitivity analyses (`tables/sensitivity.tsv`, G9)

- **`basic`** (default) reruns the context statistics, operons, GSEA and the determinants CV R² (without its permutation null) under `abundance.missing ∈ {min, exclude}`. Under the `two_peptide` tier it also counts single-peptide genes as detected.
- **`full`** adds UniOP thresholds 0.3 and 0.7, and Top3 instead of iBAQ.
- **Permutations.** Sensitivity reruns use B = 199, since they only compare flagged sets.
- **Columns:** `family`, `variant`, `n_primary`, `n_variant`, `n_shared`, `jaccard`; plus `agreement` for the dominant determinant block.
- **Families:** islands (gene sets), regions (q ≤ 0.05 IDs), operon flags (pair IDs), enriched categories (q ≤ 0.05), determinants.

### 8.7 Inherited identification error

In entrapment benchmarks, DIA tools did not consistently control precursor-level FDR, and protein-level control often failed (Wen 2025). Precursor q-values also do not control gene-level error, because false precursors pile up on otherwise undetected genes. Detection therefore uses the gene-level tiers, and validation measures the real false discovery proportion (A1's free entrapment, §15.5).

---

## 9. Output standard

Every run writes the same tree, every table has a published schema, and every plot comes with the data behind it.

### 9.1 Output tree

```text
OUT/
├── run_manifest.json          versions, DB versions + checksums, parameters, input SHA-256, seeds, timings
├── preflight/                 preflight.json, preflight.md
├── annotation/
│   ├── genes.gff3             merged gene model with source-tagged layers
│   ├── regions.tsv, regions.bed
│   ├── promoters.bed, terminators.bed, terminators_crosscheck.tsv
│   ├── operon_pairs.tsv, tus.tsv, id_map.tsv (with --annotation GFF3/GBFF)
│   └── tools/                 native outputs: bakta/, uniop/, promotech/, promoter_calculator/,
│                              terminatornet/, transtermhp/, ostir/, genomad/, antismash/
├── mapping/                   peptide_hits.tsv, protein_hits.tsv, mappings.probed, mappings.bed,
│                              mappings.gff3, unmapped.tsv, mapping_fdr.tsv
├── tables/                    genes.tsv (master), detection_summary.tsv, regions_stats.tsv,
│                              autocorrelation.tsv, islands.tsv, pairs.tsv, operons.tsv,
│                              operon_flags.tsv, pileup.tsv.gz, variants.tsv, variants.vcf,
│                              novel_orfs.tsv, novel_orfs.gff3, novel_orfs.faa, annotation_patch.gff3,
│                              category_summary.tsv, enrichment_detection.tsv, enrichment_gsea.tsv,
│                              determinants_features.tsv, determinants_model.tsv,
│                              determinants_attribution.tsv, sensitivity.tsv
├── plots/                     V01_circos.svg|png|pdf, …  +  plots/data/V01_circos.tsv, …
├── report/                    report.html (self-contained, interactive)
├── schemas/datapackage.json   Frictionless descriptor: every table, column, type, unit, coordinates
└── logs/                      <step>.log, events.jsonl, trace.tsv, timeline.html, nf_report.html
```

Optional modules that are switched off write no files; the manifest lists them under `modules_off`.

### 9.2 Table conventions

| Aspect | Rule |
|---|---|
| Format | TSV, UTF-8, one header row, `NA` for missing, one entity per row, sorted (§6.5) |
| Column names | snake_case, with units in the name: `length_aa`, `distance_bp`, `abundance_log2`, `mw_da` |
| Coordinates | 1-based closed in TSV, GFF3 and VCF; 0-based half-open in BED and proBed; declared per table in the schema |
| IDs | Bakta locus tags; `NOVEL_<contig>_<start>_<p\|m>`; `<type>_<n>` regions; `TU_<n>`; `ISL_<n>` |
| Booleans | `true` / `false` |
| Flags | `flags` column: `;`-separated controlled vocabulary (§9.4), defined in the schema |
| Test results | every test row has `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`, `null_model`, `n_perm` |
| Large tables | optional Parquet twin with `--parquet` |
| Versioning | `schema_version` in the manifest; semantic versioning; breaking column changes only in major versions |

### 9.3 Column dictionaries

Types: `str`, `int`, `float`, `bool`, `enum{…}`, `list` (`;`-separated). Units are in the names. Coordinates in TSV are 1-based closed.

**`tables/genes.tsv` (master gene table).** One row per CDS, sORF and pseudogene, in these column groups:

- **Identity:**
  - `locus_tag`, `source_id`;
  - `gene`, `product`, `feature_type` enum{cds, sorf, pseudogene}, `hypothetical`;
  - `contig`, `start`, `end`, `strand`;
  - `length_aa`, `partial`, `crosses_origin`, `replicon_type`.
- **Context:**
  - `tu_id`, `tu_position`, `tu_size`;
  - `uniop_p_prev`, `uniop_p_next`;
  - `region_ids` (list), `island_id`;
  - `replichore`, `replichore_position`, `gc_local_dev`.
- **Regulatory and sequence features:**
  - promoters and transcription: `promoter_score`, `promoter_score_tu`, `tx_rate_max_tu`, `tx_rate_sum_tu`, `terminator_conf_downstream`;
  - translation initiation: `tir_ostir`, `sd_dg`, `mrna_dg`, `spacer_bp`, `start_codon`;
  - codon usage: `cai`, `tai`, `enc`, `fop`, `gc3`, `cpb`, `at_2_10`, `rare_codons_2_10`;
  - protein: `mw_da`, `pi`, `gravy`, `tm_segments`, `n_observable_peptides`.
- **Evidence:**
  - `detected`, `detection_tier` enum{decoy_picked, mayu, two_peptide, protein_list};
  - `detection_status` enum{detected, not_detected, single_peptide_unestimated, partial_excluded};
  - `gene_q_value`, `detected_2pep`;
  - `n_peptides_unique`, `n_peptides_multi_identical`, `n_peptides_shared`, `n_peptides_tolerant`;
  - `coverage_frac`, `best_peptide_q`, `p_detect`.
- **Abundance:**
  - `abundance_raw`, `abundance_log2`, `abundance_imputed`, `imputation_policy`;
  - `abundance_type` enum{ibaq, top3, maxlfq, intensity, unknown, none}, `abundance_scale` enum{molar, relative};
  - `top3_log2`.
- **Function:** `cog_id`, `cog_category`, `go_terms` (list), `kegg_ko` (list), `ec` (list), `annotation_source` enum{bakta_full, bakta_light, user}.
- **Flags:** `flags`.

**`mapping/peptide_hits.tsv`.** One row per peptide × accepted locus. Rejected placements stay in `unmapped.tsv`.

- **Peptide:** `peptide`, `peptide_len`, `source`, `best_q_input`, `pep_input`, `score_input`.
- **Locus:** `contig`, `start`, `end`, `strand`, `frame`, `crosses_origin`.
- **Match:** `route` enum{exact, tolerant}, `n_diff`, `diff_types` (list).
- **Placement:** `placement_q`, `length_class`, `locus_class`, `uniqueness`.
- **Assignment:** `gene` (locus tag or NA), `novel_candidate_id`.
- **Flags:** `flags` (e.g., `same_db`, `strain_difference`, `engine_protein_mismatch`, `nterm_extension_evidence`).

**`mapping/protein_hits.tsv`.**

- **Protein:** `protein_id`, `protein_len_aa`, `route` enum{miniprot, tblastn}.
- **Locus:** `contig`, `start`, `end`, `strand`, `rank`.
- **Alignment:** `identity`, `query_coverage`, `reciprocal_overlap`, `class` enum{identical, near, divergent}.
- **Assignment:** `gene`.
- **Flags:** `flags` (`fs`, `st`, `two_proteins_one_gene`, `cross_strain_protein`).

**`mapping/unmapped.tsv`:** `sequence_id`, `sequence`, `kind` enum{peptide, protein}, `label` enum{contaminant, unmapped}, `reason` enum{no_hit, rejected_placement_fdr, below_identity, too_short, contaminant_id}, `contaminant_match`.

**`mapping/mapping_fdr.tsv`:** `length_class`, `n_diff`, `targets`, `decoys`, `pi_wrong`, `fdr_stratum`, `n_accepted`, `decoy_type` enum{reverse, shuffle}. Per-m total rows carry `length_class = all`. De novo precision rows carry `section = denovo_precision` with `score_bin`, `precision`, `ci_low`, `ci_high`, `c_value`, `flags`.

**`tables/detection_summary.tsv`.** Long format: `section` enum{tier, basic, detectability, qc}, `group_type`, `group`, `metric`, `value`, `ci_low`, `ci_high`, `n`, `flags`. Covers:

- the tier and its FDR curve points;
- the basic O − E rows;
- detectability calibration bins, Brier score, slope and intercept;
- QC rows: `id_mixups`, `genome_inconsistent_rate`, contaminant counts.

**`tables/regions_stats.tsv`:**

- region: `region_id`, `region_type`, `subset` (`all` or a prophage function subset), `contig`, `start`, `end`, `n_genes`, `n_eligible`;
- detection: `observed`, `expected`, `o_minus_e`, `o_over_e`;
- abundance: `hl_shift_log2`;
- test columns (§9.2): one row per statistic (`stat` enum{o_minus_e, hl_shift});
- `label` (e.g., `dosage_plus_regulation`), `flags`.

**`tables/autocorrelation.tsv`:** `contig`, `statistic` enum{morans_i, join_count_bb, oriter_slope}, `lag`, `null_model` enum{gene_label, tu_block, bootstrap}, `value`, `null_mean`, `null_lo`, `null_hi`, plus the test columns.

**`tables/islands.tsv`:** `island_id`, `direction` enum{high, low}, `contig`, `start`, `end`, `tu_first`, `tu_last`, `n_tus`, `n_genes`, `genes` (list), `mean_abundance_log2`, `detected_frac`, `scan_s`, `block_len`, plus the test columns (`p_value` = FWER, `q_value = NA`).

**`tables/pairs.tsv`:**

- `pair_id`, `class` enum{OP, NO, CV, DV, SW}, `gene_a`, `gene_b`, `strand`, `intergenic_bp`, `uniop_p`;
- evidence between the genes: `promoter_between`, `tx_rate_between`, `terminator_between`;
- per gene: `detected_a`, `detected_b`, `p_detect_a`, `p_detect_b`, `abundance_log2_a`, `abundance_log2_b`, `delta_log2`;
- tests: `test` enum{both_detected_t, one_detected_pmiss, missed_operon_decoy, none}, plus the test columns, `lfdr`;
- `flags`.

**`tables/operons.tsv`:** `tu_id`, `contig`, `start`, `end`, `strand`, `genes` (list), `n_genes`, `detected_frac`, `median_abundance_log2`, `max_abs_delta_log2`, `n_flags`.

**`tables/operon_flags.tsv`:** `pair_id`, `gene_a`, `gene_b`, `uniop_class`, `promoter_between_yes`, `terminator_between_yes`, `abundance_pattern`, `flag` enum{internal_promoter_supported, internal_terminator_supported, unexplained_discordance, conflict_terminator_in_operon, novel_operon_candidate, coabundant_despite_terminator, other}, `q_value`, `lfdr`, `evidence_note`.

**`tables/pileup.tsv.gz`:** `contig`, `codon_start`, `strand`, `gene`, `aa_pos`, `ref_aa`, `n_ref`, `n_il`, `alt_counts` (`X:n;Y:m`), `call` enum{REF, ALT, mixed, no_call}.

**`tables/variants.tsv`:**

- site: `variant_id`, `contig`, `codon_start`, `strand`, `gene`, `aa_pos`, `ref_codon`, `ref_aa`, `alt_aa`;
- `alt_type` enum{substitution, swap, block, frameshift};
- `call`, `confidence` enum{high, candidate}, `n_ref`, `n_alt`;
- `nt_changes`, `alt_codon_vcf`;
- `class_error_rate`, `class_error_n`;
- `source`, `flags` (§7.5.3).

**`tables/novel_orfs.tsv`:**

- `orf_id`, `class`, `contig`, `strand`, `frame`;
- `start_longest`, `start_best_ostir`, `end`, `length_aa_longest`, `length_aa_best`, `start_codon_best`;
- `n_peptides`, `peptides` (list), `best_peptide_q`, `ostir_tir`, `promoter_score`, `rbs_supported`;
- `class_fdr` (float or `unestimated`), `called`, `rank`;
- `overlapping_gene`, `flags`.

**`tables/category_summary.tsv`:** `collection` enum{cog, kegg, go_slim}, `category`, `name`, `n_genes`, `n_detected`, `expected`, `o_minus_e`, `ci_low`, `ci_high`, `median_abundance_log2`, `coverage`.

**`tables/enrichment_detection.tsv`:** `collection`, `term`, `name`, `direction` enum{over, under}, `n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set`, `odds_weight`, `p_classic`, `parent_group`, `coverage`, plus the test columns.

**`tables/enrichment_gsea.tsv`:** `collection`, `term`, `name`, `size_tus`, `es`, `nes`, `log2err`, `leading_edge` (list), `redundant`, `parent_group`, `coverage`, plus the test columns.

**`tables/determinants_model.tsv`:** `model` enum{lightgbm, elastic_net, univariate}, `feature_or_block`, `metric` enum{cv_r2_mean, cv_r2_sd, null_p, delta_r2, spearman_rho, auroc}, `value`, `ci_low`, `ci_high`, `n`, `n_perm`.

**`tables/determinants_attribution.tsv`:** `locus_tag`, `abundance_log2`, `predicted_log2`, `residual`, `shap_transcription`, `shap_translation`, `shap_codon`, `shap_protein`, `shap_context`, `dominant_block`, `top_decile`, `flags`.

**`tables/sensitivity.tsv`:** see §8.6.

**`annotation/regions.tsv`:** `region_id`, `region_type`, `contig`, `start`, `end`, `source_tool`, `score`, `genes` (list), `notes`.

**`annotation/operon_pairs.tsv`:** `gene_a`, `gene_b`, `uniop_p`, `source` enum{uniop, distance_model}.

**`annotation/tus.tsv`:** `tu_id`, `contig`, `start`, `end`, `strand`, `genes` (list), `upstream_start`, `upstream_end`.

### 9.4 Flag vocabulary (controlled; extend only through the schema)

| Scope | Flags |
|---|---|
| Gene | `partial`, `crosses_origin`, `pseudogene`, `multi_identical`, `shared_evidence`, `two_proteins_one_gene`, `cross_strain_protein`, `imputed_abundance`, `tolerant_evidence`, `strain_difference`, `nterm_extension_evidence`, `cai_reference`, `top3_lt3`, `detectability_constant`, `abundance_not_explained_by_sequence` |
| Placement | `same_db`, `fs`, `st`, `engine_protein_mismatch`, `origin_spanning` |
| Region and context | `small_region`, `dosage_plus_regulation`, `descriptive`, `no_oric`, `uncalibrated` |
| Operons | `rna_between`, `relative_scale`, plus the seven flag values of §7.4.5 |
| Variants | `mass_equivalent`, `snv`, `multi_nt`, `strain_difference`, `database_derived`, `start_codon`, `stop_codon` |
| Novel ORFs | `junction_spanned`, `homopolymer_junction`, `class_fdr_unestimated`, `rbs_supported` |
| Run | `anti_sd_default`, `c_assumed`, `layer_missing:<layer>`, `bakta_db_light`, `operon_source_distance_model` |

### 9.5 Schemas and validation

- **Descriptor.** `schemas/datapackage.json` is a Frictionless Data Package with one resource per table. Each field has name, type, unit, description and constraints (enums from §9.4). Each resource declares a custom `coordinates` property, `1-based-closed` or `0-based-half-open`.
- **Source of truth.** The descriptors live in `src/gippa/schemas/` and are copied into each run.
- **Tests.** Every integration and end-to-end test validates outputs with `frictionless validate`, and checks GFF3 and VCF with format validators (`pysam` for VCF; a GFF3 parser round trip).
- **Manifest.** `run_manifest.json` validates against `src/gippa/schemas/manifest.schema.json`.

---
## 10. Visualisations and report

Version 1 builds **16 plots**: V1–V3, V5–V8 and V10–V18. V4 (igv locus browser) and V9 (mass-fraction treemap) are deferred.

### 10.1 Conventions (all plots)

- **Files.** `plots/Vnn_<name>.svg`, `.png` (300 dpi) and `.pdf`, plus the exact data in `plots/data/Vnn_<name>.tsv`. Multi-instance plots add a suffix: `V03_locus_<region_or_flag_id>`.
- **One style module,** `gippa.viz.style`, sets fonts, sizes and colours for the CLI and the report alike. Figure widths are 3.5 in (single column) or 7.2 in (double).
- **Colourblind-safe palettes:** Okabe–Ito for categories; viridis or cividis for continuous values.
- **A category keeps one colour across all plots:**

  | Category | Value | Colour |
  |---|---|---|
  | Detection state | detected | `#0072B2` |
  | | not detected | `#BBBBBB` |
  | | unestimated | `#E69F00` |
  | Region type | prophage | `#D55E00` |
  | | plasmid | `#009E73` |
  | | BGC | `#CC79A7` |
  | | CRISPR | `#56B4E9` |
  | | AMR | `#F0E442` |
  | Island | high | `#0072B2` with hatching |
  | | low | `#D55E00` with hatching |

  COG letters are grouped into the four COG super-categories (Okabe–Ito hues), with letter labels.
- **Titles state the finding** once it is computed, for example "Prophage genes: 12 % detected vs 31 % expected". They are generated from the tables, never hand-written.
- **Large genomes.** Plots aggregate into bins when a track would draw more than 20,000 elements.
- **Snapshot tests.** Every plot is snapshot-tested with pytest-mpl on the simulated genome (tolerance RMS ≤ 2). Khoa approves the baselines at M9, and after that a style change is a reviewed change.

### 10.2 Plot catalogue

| ID | Name | Question | Content | Library | Interactive in report |
|---|---|---|---|---|---|
| V1 | `circos` | Where in the genome is the proteome? | contigs as sectors with ticks; genes by strand coloured by detection; abundance heat track; peptide density (10-kb bins); TU boundaries; promoter and terminator ticks; prophage, plasmid, BGC and CRISPR arcs; islands; GC skew; markers for ALT sites and novel ORFs | pyCirclize | static SVG |
| V2 | `abundance_landscape` | Is abundance clustered along the chromosome? | per replicon: log2 abundance vs. position (detected genes as points, undetected as a rug), LOESS trend, shaded islands, oriC and terminus lines | matplotlib | Plotly |
| V3 | `locus_<id>` | What happens at this prophage, BGC, island, flagged operon or novel ORF? | gene arrows coloured by abundance (cividis; undetected grey); peptide pileup track; promoter (bent arrow) and terminator (lollipop) marks; region shading. At most `report.max_locus_plots` (50), ordered by q | pyGenomeViz | static SVG |
| V5 | `mapping_sankey` | Where did my sequences go? | input class (peptide or protein) → route (exact, miniprot, tblastn, tolerant, unmapped) → outcome (placed, contaminant, unexplained, rejected by FDR) | custom matplotlib alluvial | Plotly Sankey |
| V6 | `mapping_fdr` | Can I trust tolerant hits and de novo peptides? | (a) targets vs. decoys by length class × differences; (b) stratum FDR with the 0.01 line; (c) de novo precision by score bin with CI (if de novo) | matplotlib | Plotly |
| V7 | `obs_vs_exp` | Which regions and categories are detected less than detectability predicts? | forest plot of O − E (or O/E) with 95 % CI for replicons, feature classes, COG categories, region types and individual regions; q labels | matplotlib | Plotly |
| V8 | `cog_abundance` | Which functions dominate abundance? | strips (or violins when n > 50) of log2 abundance per COG category, sorted by median, annotated with detected fraction | matplotlib | Plotly |
| V10 | `correlogram` | Is there structure beyond operons? | Moran's I vs. lag 1–20 with gene-label and TU-block 95 % envelopes; join-count panel for detection | matplotlib | Plotly |
| V11 | `operon_pairs` | Do operon partners look alike, and which other neighbours do too? | (a) Δ densities for OP, NO, CV and SW; (b) local fdr curve; (c) partner scatter (upstream vs. downstream log2) coloured by class and flag | matplotlib | Plotly |
| V12 | `variants` | Where do peptides support ALT residues? | (a) counts by call × confidence × flag; (b) peptide-vs-genome text alignments at up to 30 top ALT sites | matplotlib | static SVG + table |
| V13 | `novel_orfs` | Which new genes or upstream extensions are supported? | (a) class counts by call status; (b) locus views of the top 20 (V3 renderer) | matplotlib + pyGenomeViz | static SVG + table |
| V14 | `enrichment` | Which functions are high or low? | (a) GSEA running-score plots for the top 6 sets; (b) NES vs. −log₁₀ q dot plot; (c) detection ORA, O/E by category | matplotlib | Plotly |
| V15 | `determinants` | What explains abundance? | (a) block ΔR² with CIs; (b) partial dependence for CAI, tAI, OSTIR score and promoter score; (c) SHAP beeswarm (top 20 features); (d) observed vs. out-of-fold predicted | matplotlib | Plotly |
| V16 | `compare` | How do search strategies differ? | UpSet of detected genes across runs; pairwise abundance concordance matrix (Spearman ρ, scatter) | UpSetPlot (fallback: custom matplotlib) | Plotly (concordance) |
| V17 | `qc` | Are the inputs sane? | preflight summary panel; contig lengths; peptide length distribution; abundance distribution; missingness | matplotlib | Plotly |
| V18 | `stat_diagnostics` | Are the tests calibrated on this dataset? | detectability reliability curve; p-value histograms per test family; null distributions against the observed statistic (regions, islands, determinants) | matplotlib | Plotly |

V16 is written only by `gippa compare`. Every other plot is written by `gippa run` and `gippa report`.

### 10.3 Report (`report/report.html`)

- **Self-contained.** One HTML file built with Jinja2. Plotly.js is inlined once, and CSS and the small table-sorting script are inline; no network access is needed to view it. The size target is < 25 MB: large tables are linked by relative path, not embedded.
- **Sections:**
  1. summary cards;
  2. inputs and QC (V17);
  3. mapping (V5, V6);
  4. detection and detectability (V7, V18);
  5. genome context (V1, V2, V10, islands, regions, V3);
  6. operons (V11);
  7. variants (V12);
  8. novel ORFs (V13);
  9. function (V8, V14);
  10. determinants (V15);
  11. sensitivity;
  12. methods and caveats (generated from `docs/methods/` snippets);
  13. provenance (manifest digest).
- **Summary cards:** genes detected (by tier), mapping rates, islands, flagged operons, ALT calls, novel ORFs, model R².
- **Tables.** Every table section shows the top rows, sortable and filterable, with a link to the full TSV.
- **Wording.** "detected", not "expressed"; every number with a statistic carries its error estimate or `unestimated`.
- **Tests.** A snapshot of the HTML structure (headings and IDs). A link checker finds no broken relative links. A test asserts that no external URLs are loaded at view time.

---
## 11. Testing strategy

### 11.1 The test-driven loop

1. Each function starts as a failing test written from an example in this spec (input and expected output).
2. Then comes the smallest code that passes.
3. Then refactoring with the tests green.

Every bug gets a regression test before its fix. Statistics are proven on simulated data (§8.4) before they touch real data, and a real E. coli run gates every release.

### 11.2 Test layers

| Layer | Marker or location | What it checks | Examples | When | Budget |
|---|---|---|---|---|---|
| Unit | `tests/unit` | one behaviour per test, pure functions | table-11 translation; I = L; origin-crossing coordinates; BED ↔ GFF; BH q-values; MAYU closed form | every commit | < 60 s total |
| Property (Hypothesis) | `tests/property` | invariants over random inputs | a peptide cut from a random genome maps back to its origin; reverse complement swaps the strand; coordinate round trips; tolerant matcher = brute force; TU permutation keeps TU sizes | every commit | < 2 min |
| Calibration | `tests/calibration`, `calibration` / `calibration_full` | statistical validity on simulated nulls (§8.4) | uniform p; FDR ≤ q; island FWER ≤ 5 %; power curves recorded | fast subset per PR, full nightly | < 10 min per PR |
| Golden | `tests/golden` | exact outputs on frozen inputs | gene, pair and variant tables for the 50-kb fixture; GSEA vs. fgsea and ORA vs. goseq | every PR | < 2 min |
| Contract | `tests/contract` | parsers against recorded native outputs | Bakta 1.12.2 JSON; UniOP pairs; Promotech, Promoter Calculator, TerminatorNet, TransTermHP, OSTIR, geNOMAD, antiSMASH, miniprot, tblastn outputs; DIA-NN 2.x and Cascadia files when available | every PR and every tool upgrade | < 1 min |
| Snapshot | `tests/snapshot` | plots | pytest-mpl image comparison with tolerance; report structure | every PR | < 2 min |
| Integration | `workflow/tests` (nf-test), `tests/integration` | the workflow with real tools | the 100-kb synthetic genome through every wrapper (Bakta light DB plus recorded full-DB fixtures for Module 7); resume and failure paths; stub mode for fast workflow tests; local and Nextflow executors give identical tables | nightly, and PRs touching `workflow/` or `annotate/` | < 30 min |
| End-to-end | `tests/e2e`, `e2e` | biology on real E. coli genomes | §11.6 | weekly and before release | hours |
| Benchmarks | `benchmarks/` | speed and memory | per-step time and peak memory against `baseline.json` | nightly | — |

### 11.3 Synthetic-truth generator (`tests/sim/`, package `gippasim`)

Entry point: `gippasim.simulate(config: SimConfig, seed: int) -> SimResult`, which writes a directory of inputs plus truth. It is fast: a 100-kb genome in < 1 s, a 5-Mb genome in < 30 s.

**Genome**

- **Replicons:** a circular chromosome (50 kb – 5 Mb) and optional plasmids.
- **Coding sequence:** a 3-periodic Markov model with an E. coli-like codon usage table (bundled). Intergenic regions are AT-rich.
- **Genes:**
  - TU sizes geometric (mean 1.8 genes); intra-operon gaps −20 to +40 bp; inter-TU gaps 50–300 bp;
  - strand runs with set convergent and divergent frequencies;
  - log-normal lengths (median about 280 aa);
  - starts ATG 90 %, GTG 8 %, TTG 2 %.
- **Signals with tunable strength:**
  - RBS motifs (AGGAGG at spacer 5–10);
  - σ70 −35/−10 motifs upstream of TU starts;
  - GC hairpin + U-tract terminators after TU ends.
- **Non-coding features:** rRNA operons whose 16S ends in `…GATCACCTCCTTA`; tRNA genes with anticodons; an oriC-like region next to a marker gene.
- **Planted features:**
  - a prophage block with shifted GC;
  - a BGC-like block;
  - plasmids;
  - a CRISPR array;
  - pseudogenes (internal stop or frameshift);
  - sORFs (< 30 aa);
  - unannotated expressed ORFs (intergenic, out-of-frame, antisense);
  - N-terminal extensions (annotated start set downstream of the true start).
- **Strain differences.** A "sample" genome differs from the "reference" genome by SNVs (ALT residues) and indels.
- **Annotation.** Written as a Bakta-like JSON containing only the fields GIPPA parses (validated by GIPPA's own parser), plus GFF3. All truth tables are written too.

**Expression and proteomics**

- **TU-level log abundance** = μ + replicon effect (plasmid copy number) + region effects (for example, prophage −3) + planted islands (± effect over 5–20 TUs) + an AR(1) background (ρ tunable) + noise.
- **Gene level** = TU level + position-in-TU decay + a stoichiometry factor + noise.
- **Determinants.** Planted additive effects of promoter motif score, RBS strength and CAI with known coefficients, so the true R² is known.
- **Peptides:**
  - tryptic digest (7–30 aa, ≤ 1 missed cleavage);
  - peptide detectability logistic in peptide features and log abundance;
  - log-normal intensities.
- **Database-search errors:** false identifications at rate π_wrong, drawn from a shuffled proteome or from a related "other strain" proteome. Scores come from a two-component mixture, from which q-values and PEPs are computed consistently. Decoy rows (reversed sequences) carry null scores.
- **De novo errors:** substitutions, mass-equivalent swaps and blocks (N ↔ GG, Q ↔ AG, …), I/L, adjacent swaps, reversed stretches, modifications, and random peptides. Scores are correlated with correctness. A `spectrum_id` pairs each de novo row with a reference identification.

**Outputs**

`genome.fasta`; `reference.fasta` (for strain-difference runs); `bakta_like.json`; `peptides.tsv`; `reference_peptides.tsv`; `proteins.faa`; `abundance.tsv`; precomputed layers (`operon_pairs.tsv` with UniOP-like probabilities; promoters, Tx rates, terminators, RBS and regions, each with true elements plus false positives at known rates); `truth/*.tsv`; `config.yaml` with `annotation.precomputed` filled in.

### 11.4 Fixtures and recorded outputs

- **Synthetic.**
  - Generated at test time from fixed seeds.
  - The 50-kb golden fixture (seed `GOLDEN_SEED` in the tolerances file) has its expected tables committed under `tests/data/golden/sim50k/`.
- **Recorded native outputs** live in `tests/data/recorded/<tool>/<version>/<input>/`, with a `PROVENANCE.md` giving the command, environment lock or digest, date and host.
  - The agent records them by running the tools on the 100-kb synthetic genome on FarmShare (a batch job, or the `record-fixtures` mode of the nightly Slurm suite).
  - Re-recording happens only on a tool upgrade, which requires an ADR and a contract-test diff.
- **Full-DB Bakta fixtures (G17).** `scripts/slurm/record_bakta_full.sbatch` (16 CPU, 64 GB) runs Bakta 1.12.2 with the full DB v6.0 on:
  - K-12 MG1655 (GCA_000005845.2);
  - Nissle 1917 (GCA_043228065.1);
  - the 100-kb synthetic genome.

  The agent records them in M2, as soon as the full DB is downloaded, and commits the gzipped JSON, GFF3 and GBFF (a few MB each) under `tests/data/recorded/bakta/1.12.2-full/`. Until then, Module 7 tests use the light-DB fixture plus a synthetic cross-reference table, and the full-DB contract test is marked `pending_sample`.
- **Real genomes.** The K-12 and Nissle FASTA files (gzipped, < 2 MB each) are committed under `tests/data/genomes/`, with a provenance note. NCBI places no restrictions on these sequences.
- **Never committed:** RegulonDB files, UniProt bulk downloads (> 5 MB), databases, raw MS files. Scripts download them into `$GIPPA_DATA_DIR`.

### 11.5 Golden references from R

`scripts/golden/make_golden_r.R` runs in a conda environment (`r-base`, `bioconductor-fgsea`, `bioconductor-goseq`; conda-forge and bioconda). It writes:

- fgsea multilevel results (ES, NES, p, log2err) for frozen ranked TU scores and gene sets;
- goseq Wallenius p-values for frozen sets with a supplied PWF equal to p_detect.

Outputs go to `tests/data/golden/r/` with the R and package versions. `golden-r.yml` (manual) can regenerate them and open a pull request.

### 11.6 True-example assertions on E. coli

| ID | Assertion | Where and when |
|---|---|---|
| TE1 | ≥ 99.5 % of A2 peptides map exactly, in frame, to annotated CDS | pilot (M10). In-silico proxy: in-silico tryptic peptides of Nissle Bakta proteins map ≥ 99.9 % exactly in frame, the exceptions being N-terminal peptides at non-ATG starts |
| TE2 | no `high` ALT calls from A2 except at sites flagged `start_codon` or `stop_codon` | pilot; in-silico proxy on same-strain peptides |
| TE3 | ribosomal proteins ≥ 90 % detected in A2, with median abundance in the top quartile | pilot |
| TE4 | stoichiometry within operons: log2(RplL/RplJ) in [1, 3] (about 4×, from the rplKAJL operon); log2(AtpA/AtpD) in [−1, 1] (α3β3) | pilot (iBAQ) |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pilot; in-silico proxy using digests of porcine trypsin and human keratins added to a peptide table |
| TE6 | on K-12, Bakta's oriC call sits next to mnmG (gidA): mnmG or mioC is the nearest gene, and mnmG is within 2 kb | FarmShare (live full-DB Bakta run or recorded JSON) |
| TE7 | two runs with the same seed give byte-identical tables | every integration and end-to-end run |
| TE8 | cross-strain sanity: UniProt K-12 tryptic peptides on the K-12 genome map ≥ 99.5 % exactly; on the Nissle genome, the decoy-estimated placement FDR at defaults is ≤ 1 % and tolerant plus unmapped shares are reported | FarmShare (needs `rest.uniprot.org`) |

### 11.7 Coverage, mutation testing, determinism

- **Coverage** ≥ 90 % of `src/`, excluding the thin wrappers listed in `ci/gates.yaml`.
- **Mutation testing.** Monthly mutmut runs on `stats/` and `mapping/`. Scores are tracked in `benchmarks/mutation.json`, and a drop of more than 5 points opens an issue.
- **Determinism.** TE7 runs in integration and end-to-end tests. Any non-deterministic output fails CI.

---

## 12. CI/CD

CI is split between GitHub Actions (fast checks on every push and PR) and FarmShare Slurm jobs (everything heavy, G19). Every workflow and every FarmShare suite ends by posting a compact summary (failed tests, gate metrics, links) as a comment, so the agent can read results without downloading logs (§3.10).

| Workflow | Trigger | Jobs |
|---|---|---|
| `ci.yml` | push, pull_request | `ruff check`, `ruff format --check`, `mypy` (strict on `src/gippa`); pytest unit, property, golden, contract and snapshot tests on ubuntu (Python 3.11, 3.13) and macOS (3.12); the fast calibration subset (ubuntu, 3.12); coverage gate; nf-test in stub mode (setup-java 21, Nextflow from GitHub releases) |
| `scripts/slurm/nightly.sbatch` (FarmShare) | `scrontab` daily at 02:00 (installed by Khoa), or submitted by the agent on demand | full calibration as a job array; tool integration with `-profile farmshare,apptainer` on the 100-kb genome (Bakta light and full DB); `record-fixtures` mode (`--export=MODE=record`); benchmarks pinned to one node class with `--constraint` (median of 5; fail on > 15 % slowdown against `baseline.json` on the same node class). Tests `main` and the newest open PR branches, each in its own `git worktree` under `$GIPPA_ROOT/ci/` |
| `scripts/slurm/weekly.sbatch` (FarmShare) | `scrontab` weekly (Sunday 03:00), or on demand | end-to-end runs on K-12 and Nissle with live full-DB Bakta and with the recorded JSON; antiSMASH live run; in-silico peptide sets; TE6–TE8 proxies; monthly mutation testing (first Sunday) |
| `containers.yml` | changes to `workflow/containers/**` or `workflow/envs/**`, manual | build per-tool images, push to `ghcr.io/<owner>/gippa-<tool>`, open a PR updating the digests in `workflow/conf/containers.config` |
| `golden-r.yml` | manual | regenerate the R golden data and open a PR (may also be done on FarmShare with the conda R env) |
| `release.yml` | manual (version input) | verify gates, changelog and the version; build sdist and wheel; publish to PyPI (trusted publishing); create the tag and GitHub release (the Zenodo DOI comes from the GitHub integration). Bioconda recipe PR as a manual follow-up |

**`ci/gates.yaml`** (protected) holds:

- the coverage minimum (0.90);
- the per-layer time budgets (§11.2);
- the benchmark slowdown tolerance (0.15);
- the required capabilities per job, including the FarmShare suites (§3.6);
- the FarmShare resource caps of R14 (`farmshare: {max_cpus: 128, max_gpus: 4, max_jobs: 20}`);
- the list of thin wrappers excluded from coverage.

**Releases** use semantic versioning, changelog fragments (towncrier), PyPI and Bioconda packages, containers pinned by digest, and a Zenodo DOI each time.

**Actions minutes.** Heavy suites run on FarmShare, so Actions minutes are rarely a concern. If the private repository runs short, trim the macOS job first, then ask Khoa.

**FarmShare suite results** are posted to the "Nightly (FarmShare)" tracking issue with `gh issue comment`; a nightly failure becomes a top-priority issue. Job logs stay under `$GIPPA_LOG_DIR`.

---
## 13. Milestones and gates

Effort is for one developer working with agents, about 15 weeks in total. After M4, milestones M5–M8 can run in parallel Claude Code sessions, each in its own interactive allocation and its own `git worktree` (under `local/worktrees/<milestone>`), one milestone per session on separate branches. Parallel sessions share the R14 caps. A milestone is done when its gate passes and `PROGRESS.md` records the evidence: test names, metrics and PR links.

| # | Milestone | Effort | Deliverables | Gate |
|---|---|---|---|---|
| M0 | Scaffold | 0.5 wk | repository connected to `git@github.com:khoa-yelo/gippa.git`, `.gitignore` with `/local/`; `farmshare_env.sh`, `bootstrap_farmshare.sh` (run); Appendix A files; `pyproject.toml` + `uv.lock`; CLI skeleton with all subcommands as stubs; config model + `defaults.yaml` + sync test; schema skeleton; `probe_env.py`, `session_start.sh`; the `scripts/slurm/` job header and one trivial test job; `ci.yml`; simulator skeleton (genome, genes, truth JSON) | CI green; a toy genome simulated; a `gippa-` test job ran through Slurm; the probe report committed in `PROGRESS.md` |
| M1 | Inputs and preflight | 0.5 wk | `io/` readers and writers (FASTA, GFF3, BED, proBed, VCF, TSV); generic peptide table; `gippa convert` (generic, DIA-NN from its documentation); every preflight rule | every preflight rule unit-tested; `preflight.json` validated against its schema |
| M2 | Annotation wrappers and benchmarks | 1.5 wk | `workflow/envs/*`, Nextflow skeleton with stubs, wrappers and parsers for every tool in §7.1.1; `gippa setup`; recorded fixtures for the 100-kb genome; UniOP adapter; TUs; promoter windows; anti-SD; codon usage; regions | a contract test per tool; the full Bakta DB downloaded and the full-DB fixtures recorded (G17); `nightly.sbatch` and `weekly.sbatch` proposed to Khoa for `scrontab`; measured FarmShare runtimes per tool (with node class) replace the estimates |
| M3 | Mapping and decoy FDR | 2 wk | six-frame index, exact matching, miniprot and tblastn, tolerant matcher, decoys, placement FDR, loci, contaminants, Module 2 outputs; `gippa makedb`; **pilot kit (§15) prepared**, and Khoa asked to install DIA-NN | ≥ 99.9 % recall of planted peptides; mapping FDR calibrated; TE1 and TE8 proxies pass |
| M4 | Gene table, detectability, region tests | 1.5 wk | gene table, detection tiers, detectability, abundance, basic statistics, region tests | circular-shift tests calibrated; gene-level FDR calibrated (picked, MAYU) |
| M5 | Autocorrelation and islands | 1 wk | Moran's I and join counts with both nulls, islands, ori–ter | island FWER ≤ 5 % under simulated nulls; power curve recorded |
| M6 | Operons | 1.5 wk | pair table, discordance, missed operons, flags; RegulonDB calibration on the K-12 genome (§15.4) | decoy FDR calibrated on simulations; UniOP and TerminatorNet cut-offs recalibrated (ADR); UniOP AUROC recorded |
| M7 | Variants and novel ORFs | 2 wk | Modules 5 and 6, VCF, patch GFF3, second-round Promotech and OSTIR | F1 ≥ 0.9 on simulated variants and novel ORFs; TE2 proxy passes |
| M8 | Function and determinants | 1.5 wk | Modules 7 and 8; fgsea port; Wallenius ORA; R golden data | GSEA and ORA match fgsea and goseq on golden data; null p-values uniform |
| M9 | Visualisation and report | 1.5 wk | 16 plots, style module, report, `gippa report`, `gippa compare` + V16 | **Khoa approves** the plot snapshots and report layout |
| M10 | Pilot validation | 1 wk | §15: raw data download, DIA-NN (A1, A2) and de novo (A3) searches as Slurm jobs, conversion, runs A1/A2/A3, `gippa compare`, pilot report, the 16-thread benchmark | the pilot passes (§15.6) and **Khoa has reviewed** its report |
| M11 | Release 1.0 | 0.5 wk | packages, containers by digest, MkDocs site, changelog, citation file, preprint outline | the acceptance criteria below; **Khoa approves** the release |

### 13.1 The first session (M0), in order

1. Read this whole file.
2. Connect the folder to GitHub without touching the files already here:

   ```bash
   git init -b main && git remote add origin git@github.com:khoa-yelo/gippa.git
   git fetch origin && git reset origin/main          # adopt the remote's initial commit; keeps working files
   git checkout origin/main -- README.md
   printf '/local/\n.gippa/\n' >> .gitignore
   git switch -c m0/scaffold && git add GIPPA_SPEC.md .gitignore && git commit -m "docs: add GIPPA_SPEC.md"
   git push -u origin m0/scaffold
   ```

   Then run the bootstrap (§3.4) and ask Khoa for `! gh auth login` (§0.1 step 3), so the M0 PR can be opened with `gh`.
3. Create `AGENTS.md`, `CLAUDE.md`, `PROGRESS.md`, `.claude/` (settings, hooks, agents), `ci/gates.yaml`, `ci/protected.txt`, and the PR and ADR templates from Appendix A.
4. Write `scripts/farmshare_env.sh`, `scripts/session_start.sh` and the probe script; run the probe and record the capabilities and any unreachable sources in `PROGRESS.md`. Submit one trivial `gippa-probe` job with the §3.8 header to prove the batch path.
5. Create `pyproject.toml`:
   - Python ≥ 3.11;
   - runtime dependencies: numpy, pandas, pyarrow, numba, scipy, scikit-learn, statsmodels, lightgbm, pyahocorasick, pydantic, click, pyyaml, jinja2, matplotlib, pycirclize, pygenomeviz, upsetplot, plotly, goatools, codon-bias, biopython;
   - dev extras: pytest, hypothesis, pytest-benchmark, pytest-mpl, pytest-xdist, pysam, frictionless, mutmut, ruff, mypy, towncrier;
   - docs extras: mkdocs-material.

   Lock with `uv lock`.
6. Build the CLI skeleton, the config model with `defaults.yaml` exactly as in §5.2, and a test that the two agree.
7. Add `ci.yml`, then open the "M0 scaffold" PR.
8. Write the simulator skeleton test-first: genome + genes + a Bakta-like JSON for a 50-kb genome, deterministic by seed.
9. Update `PROGRESS.md` with the next actions for M1.

### 13.2 Acceptance criteria for release 1.0

- ≥ 99.5 % of A2 peptides map exactly, and the decoy-estimated mapping FDR is ≤ 1 % at defaults.
- Every statistical method passes its calibration suite (§8.4), and every attack null (§8.5) stays at level.
- Every true-example assertion that the available data allow (§11.6) passes.
- The same seed gives byte-identical tables.
- A standard run on 16 threads and 64 GB takes about 35 min (range 16–55; measured by `benchmark_standard.sbatch` on FarmShare, node class recorded). Core modules alone take about 25 min.
- The pilot passes (§15.6).
- "Adding co-abundance does not lower K-12 operon AUROC" moves to the full validation, because it needs K-12 proteomics.

---

## 14. Agent workflow

Agents build GIPPA in small, test-first pull requests. Separate subagents write the tests, the code and the reviews, and four improvement loops run against frozen benchmarks.

### 14.1 Roles

The lead is the main Claude Code session; every other role is a subagent defined in `.claude/agents/` (Appendix A.5).

| Role | Owns | May change | May not change |
|---|---|---|---|
| Lead (orchestrator) | backlog from this spec, issues, merge decisions, `PROGRESS.md` | issues, ADR drafts, `PROGRESS.md`, `docs/backlog.md` | merge without green CI and reviews |
| Spec writer | issue specs: examples, edge cases, acceptance criteria | `docs/specs/` | code, tests |
| Test author | failing unit, property and calibration tests | `tests/` (except protected tolerances) | `src/`, `workflow/` |
| Implementer | the smallest code that passes | `src/`, `workflow/`, `scripts/`, new fixtures under `tests/data/recorded/` with provenance | existing tests, tolerances, gates |
| Code reviewer | correctness, clarity, complexity, mutation score of changed code | review comments | the code it reviews |
| Statistician (independent) | every statistical change; reads the methods note and simulations before the code; owns the attack nulls; can block | `docs/methods/`, `tests/calibration/attacks/` | implementation |
| Performance | profiling (py-spy, scalene); optimisation with golden outputs unchanged | hot paths, `benchmarks/` | outputs, tolerances |
| Visualisation | renders plots and inspects the images against the style checklist | `src/gippa/viz/`, snapshot proposals | approved snapshots without Khoa's sign-off |
| Validation | end-to-end runs and the pilot | `validation/` | code or thresholds; never tunes on the holdout |
| Docs | user guide, tutorials, changelog | `docs/`, `CHANGELOG.md` | code |

### 14.2 Branches, commits, pull requests

- **Branch names:** `m<k>/<topic>`, for example `m3/tolerant-matcher`.
- **Commits** are conventional (`feat:`, `fix:`, `test:`, `docs:`, `perf:`, `refactor:`, `ci:`) with a trailer `Agent-Role: <role>`.
- **Pull request size:** about 400 changed lines or fewer, excluding lock files and fixtures.
- **Every PR links** an issue (or a `docs/backlog.md` entry) and the spec section it implements, and includes a changelog fragment.
- **Labels:** `needs-khoa` for anything in §14.8; `stats` for statistical changes, which need the statistician's approval.

### 14.3 Feature loop, one issue at a time

1. **Spec.** The spec writer writes the issue: the spec section, examples as input and expected output, edge cases, acceptance criteria.
2. **Failing tests.** The test author commits failing tests, so CI is red as expected.
3. **Code.** The implementer writes the smallest change that turns them green, then refactors.
4. **Review.** The code reviewer and, for statistical code, the statistician review independently and write their verdicts as PR comments.
5. **CI gates pass:** tests, coverage, calibration, golden files, snapshots, benchmarks, the role check and the protected-file check.
6. **Merge.** Squash-merge with `gh pr merge --squash`. The lead may merge when steps 4–5 pass and the PR has no `needs-khoa` label. Nightly suites run the slow checks afterwards, and a nightly failure becomes a top-priority issue.

### 14.4 Improvement loops

| Loop | Metric | Data | Stop rule | When |
|---|---|---|---|---|
| Speed | wall time and peak memory per step (`sacct` MaxRSS, same node class) | the 5-Mb simulated genome and K-12 | < 5 % gain over 3 attempts, or the step budget is met | after M3, M8 and M10 |
| Accuracy | mapping recall and precision; operon AUROC; novel-ORF precision; variant-class F1 | simulations and development genomes; never the holdout | no gain over 3 attempts without regressions elsewhere | after M3, M6 and M7 |
| Statistical rigour | uniform p-values, FDR ≤ q, FWER, power | null simulations including attack nulls | every test calibrated and every attack fails | continuous; owned by the statistician |
| Visualisation | style-checklist score, then Khoa's preference | plots rendered from the simulated and E. coli runs | Khoa approves and the snapshot is locked | M9 |

Each accepted improvement updates `benchmarks/baseline.json` through an ADR with before-and-after numbers.

**Style checklist (visualisation loop):**

- the title states the finding;
- axes carry units;
- palette as in §10.1;
- no overlapping labels;
- legible at single-column width;
- data TSV written;
- colourblind check (simulated deuteranopia);
- ≤ 20,000 drawn elements.

### 14.5 Sealed holdout (for the later full validation)

E. coli O157:H7 EDL933 (PXD068911, TMT DDA) plus one non-E. coli genome with public proteomics, to be chosen. They are never used in development; a second condition of a development strain would not be independent. In version 1 the holdout is not downloaded at all. When the full validation starts, it is checked only at milestone gates, with access logged and checksummed in `validation/holdout/ACCESS_LOG.md`.

### 14.6 Guardrails and how they are enforced

| Guardrail | Enforcement |
|---|---|
| Tests and implementation never come from the same agent in one PR | `Agent-Role` trailers. `scripts/ci/check_roles.py` (CI) fails a PR if the same role changes both `src/` and test files, if implementer commits modify existing test files, or if test-author commits touch `src/` or `workflow/` |
| Agents cannot edit CI gates, coverage thresholds, calibration tolerances or fixed seeds without Khoa's approval | `.claude/settings.json` deny rules and a PreToolUse hook (`.claude/hooks/protect_paths.py`) block Edit and Write on the globs in `ci/protected.txt`. `scripts/ci/check_protected.py` fails CI when a protected file changes unless the PR carries a `Gate-Change-Approved: ADR-NNNN` trailer naming an ADR whose status is `accepted-by-khoa` |
| Any change to a default or threshold needs an ADR with before-and-after metrics | review checklist; `check_protected.py` also guards `src/gippa/config/defaults.yaml` and `src/gippa/data/calibration.yaml` |
| PRs stay under about 400 changed lines | CI warning above 400 lines, failure above 800 (lock files and fixtures excluded) |
| External tool versions are pinned; upgrades go through contract tests | env files with exact pins; `containers.config` by digest; contract tests |
| GPL, AGPL and unlicensed tools are only called through their command lines | `scripts/ci/check_imports.py` fails on `import` of bakta, ostir, promoter_calculator, antismash, promotech or UniOP modules in `src/` |
| Seeds live in the config; non-deterministic output fails CI | TE7 in integration and end-to-end tests |
| After 3 failed iterations in one loop, stop and escalate | R9; `PROGRESS.md` "Questions for Khoa" |

### 14.7 What agents read and keep current

| File | Contents |
|---|---|
| `AGENTS.md` | setup, test commands per layer, conventions, definition of done, where specs live |
| `PROGRESS.md` | the current milestone, done, in progress, next actions, blocked, questions for Khoa, environment notes |
| `docs/adr/` | one record per decision, starting with ADR-0001, which imports §2.1 and §2.4 |
| `docs/methods/` | one note per statistical method with its simulation evidence; the statistician signs off |
| `benchmarks/baseline.json` | the accepted runtime and accuracy numbers |

### 14.8 Khoa's checkpoints: the only things the agent waits for

1. ADRs that change a default, threshold, interface, dependency or this spec's decisions, including G1–G19.
2. Any change to a protected file.
3. Plot snapshot baselines and the report layout (M9).
4. The pilot report (M10).
5. Each release (M11).
6. Licence and account steps only Khoa can take:
   - accepting the DIA-NN licence and installing its binary under `$GIPPA_TOOLS_DIR` (§15.2);
   - `gh auth login` and the git identity (§3.10).
7. Environment changes only Khoa can make: installing or changing `scrontab` entries, requests to FarmShare support (quota, QOS, software), raising the R14 caps, making the repository public, Actions budget.

Everything else, the agent decides, records and continues (R10).

### 14.9 Definition of done for a module

- Tests pass at every layer.
- Coverage is ≥ 90 %.
- A methods note exists for each statistical method, signed off by the statistician.
- The schema is updated and validated.
- The benchmark is recorded.
- The code and statistics reviews approve.
- `PROGRESS.md` is updated.

### 14.10 Working in long sessions on FarmShare

- Use subagents for exploration, reviews and test writing to keep the lead's context small. Subagents share the session allocation; never run more than one heavy test suite at a time in it.
- Commit and push at least at the end of every feature-loop step, so an expired allocation loses nothing.
- Check the time left in the allocation (the probe prints it). With under an hour left, finish the current step, update `PROGRESS.md` (including running job IDs), commit and push.
- Start a new session per milestone if the context grows large; the resume prompt (§0.3) plus `PROGRESS.md` is enough.
- Run long work as Slurm jobs (§3.8), and prefer `nextflow -resume`.

---
## 15. Pilot validation kit: E. coli Nissle 1917 (M10)

A light pilot (decision 19). The full validation follows only if the pilot looks good.

### 15.1 Data and arms

- **Sample data.** Three wild-type runs from **PXD078459** (iProX; Nissle 1917; timsTOF HT). The deposited Spectronaut file suggests diaPASEF; confirm the acquisition mode on the first run. The download is about 5–10 GB (estimate).
- **Genome.** **GCA_043228065.1**: chromosome, pMUT1 and pMUT2. Nissle 1917 is phylogroup B2 and K-12 is phylogroup A, so a K-12 database is a genuine mismatch.

| Arm | Search | Purpose |
|---|---|---|
| A1, generic reference | DIA-NN 2.7.0, library-free, against UniProt K-12 **UP000000625** (4,403 proteins) + contaminants | what most users do; measures cross-strain losses and artefacts |
| A2, matched strain | DIA-NN 2.7.0, same settings, against Nissle's Bakta proteins (sORFs included) + contaminants, from `gippa makedb` | best case for an annotated search |
| A3, de novo | Cascadia 0.0.7 on the same runs (exploratory), on the `gpu` partition (1 × L40S). If DIA de novo fails on timsTOF data: Casanovo 5.2.1 on 3 UTI89 DDA runs (PXD053727; genome GCA_000013265.1) | database-free evidence for novel ORFs and variants |

### 15.2 Who does what

| Step | Who | Action |
|---|---|---|
| 1 | agent | Run `scripts/slurm/record_bakta_full.sbatch` (Bakta 1.12.2, full DB v6.0) on Nissle and K-12 if M2 has not already done so. Commit the gzipped JSON/GFF3/GBFF under `tests/data/recorded/bakta/1.12.2-full/` |
| 2 | agent | Run `gippa makedb --annotation <nissle JSON> --outdir validation/pilot/db_a2`. Commit `strain_db.faa.gz` and a `README.md` with the exact DIA-NN settings for A1 and A2 (`validation/pilot/diann_A1.cfg`, `diann_A2.cfg`) |
| 3 | Khoa | Accept the DIA-NN 2.7.0 licence and unpack the Linux binary under `$GIPPA_TOOLS_DIR/diann-2.7.0`; check every flag of the template below against its documentation |
| 4 | agent | Download the 3 wild-type runs of PXD078459 to `$GIPPA_DATA_DIR/pxd078459/` (a batch job; record checksums). Confirm the acquisition mode on the first run |
| 5 | agent | Run DIA-NN for A1 and A2 with `scripts/slurm/diann_pilot.sbatch` (`normal`, 16 CPU, 64 GB), including decoy reporting, so GIPPA can use the picked tier |
| 6 | agent | Run Cascadia on the same runs with `scripts/slurm/denovo_pilot.sbatch` (`--partition=gpu --gres=gpu:1`); fall back as in §15.1 |
| 7 | agent | Run `gippa convert --from diann …` and `--from cascadia …` on FarmShare; the large reports stay in `local/data/`. Commit the gzipped generic tables (a few MB) to `validation/pilot/inputs/`, plus a 1,000-row sample of each raw report under `tests/data/recorded/diann/2.7.0/` and `…/cascadia/0.0.7/` for contract tests (DIA-NN output tables are not binaries, so R11 allows the sample) |
| 8 | agent | Run `gippa run` for A1, A2 and A3 as batch jobs; then `gippa compare out_a1 out_a2 out_a3 --labels A1,A2,A3`; then write `validation/pilot/REPORT.md` (§15.5–15.7) |
| 9 | agent | Run `scripts/slurm/benchmark_standard.sbatch` (16 CPU, 64 GB, one pinned node class, median of 3); the measured runtimes replace the estimates in the docs |

DIA-NN template. Khoa checks every flag against the DIA-NN 2.7.0 documentation before the agent's first run (step 3).

```bash
diann --f run1.d --f run2.d --f run3.d --lib "" --fasta-search --predictor --gen-spec-lib \
      --fasta <A1: UP000000625.fasta | A2: strain_db.faa> --fasta contaminants.fasta \
      --cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --met-excision \
      --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4 \
      --qvalue 0.01 --report-decoys --threads 16 --out report.parquet
```

GIPPA runs:

```bash
gippa run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a1.tsv.gz --outdir out_a1
gippa run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a2.tsv.gz --outdir out_a2
gippa run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a3.tsv.gz \
         --reference-peptides a2.tsv.gz --outdir out_a3
```

All three runs use the default `--genome-relation related`.

### 15.3 Data that is not mass spectrometry

| Item | Source | Use |
|---|---|---|
| K-12 MG1655 genome | GCA_000005845.2 | TE6, TE8, RegulonDB calibration |
| UniProt K-12 proteome | UP000000625 | A1 database (Khoa); in-silico TE8 |
| RegulonDB v12 | transcription units and terminators | §15.4. Not committed; downloaded to `$GIPPA_DATA_DIR` by `validation/fetch_regulondb.py` |

### 15.4 Cut-off calibration on K-12 (no MS data; runs at M6)

1. **Inputs.** The K-12 Bakta JSON (recorded), UniOP and TerminatorNet run on the K-12 genome.
2. **Pair truth.** A same-strand adjacent pair is positive when both genes lie in one RegulonDB transcription unit. It is negative when they lie in different TUs and no TU contains both. Everything else is ambiguous and excluded.
3. **UniOP.**
   - Record the AUROC.
   - Choose the smallest threshold with precision ≥ 0.90 (if none reaches 0.90, the threshold with maximum F1).
   - Refit the distance-model parameters d₀ and s.
4. **TerminatorNet.**
   - Precision = the share of calls in intergenic regions downstream of genes that lie within ±50 bp of a RegulonDB terminator on the same strand.
   - Choose the smallest confidence with precision ≥ 0.90 (else maximum F1).
   - RegulonDB is incomplete, so this precision is a lower bound; say so in the ADR.
5. **Record** the cut-offs in `src/gippa/data/calibration.yaml`, with RegulonDB version, date, AUROC, precision and recall. The ADR goes to Khoa (protected file).

### 15.5 Metrics in the pilot report

1. Identifications at 1 % FDR per arm: precursors, peptides, protein groups.
2. Mapping: exact, tolerant and unmapped shares; genome-inconsistent peptides in A1.
3. Coverage: share of strain CDS detected; strain-specific, plasmid, prophage and BGC genes detected per arm.
4. Quantitative agreement: Spearman ρ between A1 and A2 on shared genes.
5. Stability of conclusions: Jaccard overlap of islands, operon flags and enriched categories across arms (`gippa compare`).
6. Novel ORFs and variants by class and confidence, and their overlap between A2 and A3.
7. Runtime and peak memory per arm and per step.
8. **A1's free entrapment.** Peptides absent from Nissle's translation (no exact match after the I/L and mass-equivalent rules) cannot be true if the genome is complete. With N_E such peptides, N_T peptides present, and r the in-silico ratio of K-12-only to shared tryptic peptides (computed with the same exclusions, within the searchable length, charge and m/z range), the combined estimator of Wen 2025 gives:

   ```latex
   \widehat{\mathrm{FDP}} = \frac{N_E \,(1 + 1/r)}{N_T + N_E}
   ```

   These entrapment peptides are near-homologs, easier to mismatch than random sequences, so the estimate errs high. Only wild-type runs count, because engineered strains may carry K-12-derived constructs. Implemented in `validation/pilot/entrapment.py`.

### 15.6 Pilot pass criteria

- GIPPA runs end to end on A1, A2 and A3.
- ≥ 99.5 % of A2 peptides map exactly (TE1).
- A1 and A2 detections and abundances are compared (`gippa compare`).
- Measured runtimes replace the estimates.
- TE2–TE5 are reported. A failure is investigated and explained in the report before Khoa's review.

### 15.7 Hypotheses fixed before running

These are reported where the pilot data allow; the rest wait for the full validation.

- A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes.
- A1's genome-inconsistent FDP is at least its entrapment FDP (full validation, FDRBench), since it errs high.
- Context conclusions (islands, region tests, operon concordance) agree between A1 and A2 on shared genes.
- The genome-anchored de novo precision agrees with precision measured against A2 identifications.

---

## 16. Deferred and excluded: do not build in version 1

| Status | Item |
|---|---|
| Later | Several samples in `abundance.tsv`, true co-expression and differential context tests (improvement 4) |
| Later | Assembly check from frameshift-spanning peptides; the `possible_assembly_error` flag (7) |
| Later | CheckM2 completeness, contamination and species check in preflight (10) |
| Later | Proteome mass fractions and the V9 treemap (11) |
| Later | RO-Crate export (14) |
| Later | GC test of low-abundance islands for H-NS silencing (15) |
| Later | Prophage induction signature and flag (16) |
| Later | RNA-seq and Ribo-seq inputs (18) |
| Later | V4 igv-reports locus browser |
| Later | Taxonomy lookup of unmapped sequences for mixed cultures |
| Later | Homology and conservation evidence for novel ORFs |
| Later | Full validation: K-12 PXD075568, UTI89 PXD053727 (beyond the A3 fallback), EDL933 holdout, FDRBench entrapment, K-12 co-abundance AUROC |
| No | Recoding and start catalogue (Sec, prfB +1, dnaX −1, non-ATG starts, Met removal) (6) |
| No | Start-codon correction from N-terminal peptides (8) |
| No | Signal-peptide or transmembrane predictor (17); Kyte–Doolittle segments are computed internally (G6) |
| No | Plass; eggNOG-mapper; MMseqs2/UniRef function annotation |
| No | PTM naming or classification |
| No | GIPPA running spectrum searches itself; it only writes the `makedb` FASTA |

---

## 17. Risks and fallbacks

| Risk | Signal | Fallback |
|---|---|---|
| Promotech's frozen environment (Python 3.6, TF 2.1, scikit-learn 0.23.1) will not build | conda solve or import fails | conda-forge `python=3.6` + pip wheels; else a `python:3.6-slim` image with pip, built in `containers.yml` and run with Apptainer. After 3 failed iterations, escalate with an ADR proposing promoter evidence from Promoter Calculator only (`promoter_source = calculator_only`) |
| Promotech models unreachable (`www.cs.mun.ca`) | probe | `needs_promotech_models` tests skip; retry later; ask Khoa for a copy |
| TerminatorNet pickle incompatible with recent scikit-learn | load error | pin scikit-learn in its environment to the pickle's version |
| UniOP repository changes or disappears | clone fails at the pinned commit | distance-model fallback; ask Khoa for a copy |
| Bakta full DB download fails or exceeds the `local/` budget | `gippa setup` exit 2, or the probe's disk warning | retry as a batch job; light DB for live integration in the meantime; recorded fixtures for tests |
| antiSMASH image or databases fail to pull or run under Apptainer | pull or run error | conda `antismash=8.0.4` profile; recorded fixtures; escalate after 3 attempts (R9) |
| A tool's image misbehaves under Apptainer (writes into its own image, expects root, needs `$HOME`) | permission or read-only file-system errors | bind a writable directory, set `--containall` / `HOME`, or use the conda profile; record the fix in the module's `meta.yml` |
| fgsea port misses the golden tolerances | golden tests | escalate; ADR to run R fgsea as an external tool |
| Cascadia fails on timsTOF data | A3 empty or degenerate | Casanovo on 3 UTI89 DDA runs |
| PXD078459 acquisition mode differs | first-run check | adapt the DIA-NN settings; record in the pilot README |
| DIA-NN 2.7.0 flags differ from the template | Khoa's check | Khoa edits the configs; the converter's contract test covers column names |
| Actions minutes run out (private repository) | billing warnings | trim the macOS job; ask Khoa about making the repository public |
| Session allocation expires (24 h) mid-task | session ends | batch jobs, `-resume`, frequent pushes, job IDs in `PROGRESS.md`; Khoa restarts the session with the resume prompt |
| Long queue waits on a busy cluster | jobs `PENDING` for hours (`squeue` reason `Resources`/`Priority`) | smaller requests, the `dev` QOS for short checks, a different node class; never exceed the R14 caps to get ahead |
| Class share full or over quota | "no space left on device" or "Disk quota exceeded"; missing files | delete work directories, archives and caches; re-create databases with `gippa setup`; never conclude it is unrecoverable |
| Home quota full (caches written to `$HOME`) | write errors in `~/.cache`, `~/.conda`, `~/.apptainer` | check that `farmshare_env.sh` was sourced; move the cache to `$GIPPA_ROOT/cache` |
| Docker Hub pull rate limit | Apptainer pull fails with HTTP 429 | retry later; prefer `ghcr.io` images built by `containers.yml` |
| Nextflow on NFS (stale file handles, slow `.command.*` writes) | sporadic task failures | `process.scratch = '/tmp'`, retry once (§6.3), keep the work directory in `local/work` |
| DIA-NN not yet installed | `needs_diann` missing | continue with other work; the pilot waits on Khoa's checkpoint (§14.8) |
| UpSetPlot breaks with pandas 3 | V16 tests fail | small custom matplotlib UpSet |
| Download source down or changed | timeout, HTTP error, checksum mismatch | R6: retry once later, log it, skip with a reason, ask Khoa; never substitute an unofficial copy |

---

## 18. References

Versions and figures were checked on 2026-10-13/14.

**Tools and formats**

- Bakta: https://github.com/oschwengers/bakta; paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544; DB on Zenodo: https://doi.org/10.5281/zenodo.14916843
- UniOP: https://github.com/hongsua/UniOP; preprint: https://doi.org/10.1101/2024.11.11.623000
- Promotech: https://github.com/BioinformaticsLabAtMUN/PromoTech; paper: https://doi.org/10.1186/s13059-021-02514-9
- Promoter Calculator (Barrick Lab fork): https://github.com/barricklab/promoter-calculator; paper: https://www.nature.com/articles/s41467-022-32829-5
- OSTIR: https://github.com/barricklab/ostir
- TerminatorNet: https://github.com/btjaden/TerminatorNet; paper: https://doi.org/10.1093/bioinformatics/btag116. TransTermHP: https://transterm.cbcb.umd.edu/. BacTermFinder: https://pmc.ncbi.nlm.nih.gov/articles/PMC11890068
- geNOMAD: https://github.com/apcamargo/genomad; DB v1.9: https://zenodo.org/records/14886553
- antiSMASH: https://docs.antismash.secondarymetabolites.org/
- miniprot: https://lh3.github.io/miniprot/miniprot.html
- codon-bias: https://github.com/alondmnt/codon-bias; pyCirclize: https://pypi.org/project/pycirclize/; pyGenomeViz: https://pypi.org/project/pygenomeviz/
- proBed: https://www.psidev.info/probed; format paper: https://doi.org/10.1186/s13059-017-1377-x
- Contaminant library: https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics; paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10040255
- DIA-NN: https://github.com/vdemichev/DiaNN; FDRBench: https://github.com/Noble-Lab/FDRBench
- Casanovo: https://casanovo.readthedocs.io/; InstaNovo: https://instadeepai.github.io/InstaNovo/; Cascadia: https://cascadia.readthedocs.io/
- RegulonDB v12: https://academic.oup.com/nar/article/52/D1/D255/7424428; UniProt UP000000625: https://rest.uniprot.org/proteomes/UP000000625.json
- Nextflow: Di Tommaso 2017, https://doi.org/10.1038/nbt.3820; nf-core: Ewels 2020, https://doi.org/10.1038/s41587-020-0439-x
- FarmShare documentation: https://docs.farmshare.stanford.edu/
- Nextflow Slurm executor and Apptainer: https://www.nextflow.io/docs/latest/executor.html#slurm, https://www.nextflow.io/docs/latest/container.html#apptainer
- Apptainer: https://apptainer.org/docs/user/1.4/

**Datasets**

PXD078459 (Nissle 1917), PXD053727 (UTI89), PXD075568 (K-12) and PXD068911 (EDL933 holdout), all at https://proteomecentral.proteomexchange.org/.

**Statistics and FDR**

- Multiple testing and local fdr: Benjamini & Hochberg 1995, https://doi.org/10.1111/j.2517-6161.1995.tb02031.x; Storey & Tibshirani 2003, https://doi.org/10.1073/pnas.1530509100; Efron 2004, https://doi.org/10.1198/016214504000000089
- Spatial statistics: Kulldorff 1997, https://doi.org/10.1080/03610929708831995; Moran 1950, https://doi.org/10.1093/biomet/37.1-2.17
- Target–decoy and FDR in proteomics: Elias & Gygi 2007, https://doi.org/10.1038/nmeth1019; Fu & Qian 2014, https://doi.org/10.1074/mcp.O113.030189; Nesvizhskii 2014, https://doi.org/10.1038/nmeth.3144; Wen 2025, https://doi.org/10.1038/s41592-025-02719-x
- Picked protein FDR: Savitski 2015, https://doi.org/10.1074/mcp.M114.046995
- MAYU: Reiter 2009, https://doi.org/10.1074/mcp.m900317-mcp200
- Missing values: Lazar 2016, https://doi.org/10.1021/acs.jproteome.5b00981
- goseq: Young 2010, https://doi.org/10.1186/gb-2010-11-2-r14
- GSEA: Subramanian 2005, https://doi.org/10.1073/pnas.0506580102; fgsea: Korotkevich 2021, https://doi.org/10.1101/060012
- topGO elim: Alexa 2006, https://doi.org/10.1093/bioinformatics/btl140

**Proteomics and biology**

- Quantification: Schwanhäusser 2011 (iBAQ), https://doi.org/10.1038/nature10098; Silva 2006 (Top3), https://doi.org/10.1074/mcp.M500230-MCP200; Liebermeister 2014, https://doi.org/10.1073/pnas.1314810111
- Genome context:
  - Li 2014, https://doi.org/10.1016/j.cell.2014.02.033
  - Bryant 2014, https://doi.org/10.1093/nar/gku828
  - Couturier & Rocha 2006, https://doi.org/10.1111/j.1365-2958.2006.05046.x
  - Salgado 2000, https://doi.org/10.1073/pnas.110147297
  - Navarre 2006, https://doi.org/10.1126/science.1128794
  - Lucchini 2006, https://doi.org/10.1371/journal.ppat.0020081
- Sequence determinants:
  - Sharp & Li 1987 (CAI), https://doi.org/10.1093/nar/15.3.1281
  - dos Reis 2004 (tAI), https://doi.org/10.1093/nar/gkh834
  - Wright 1990 (ENC), https://doi.org/10.1016/0378-1119(90)90491-9
  - Goodman 2013, https://doi.org/10.1126/science.1241934
  - Kudla 2009, https://doi.org/10.1126/science.1170160
  - Kyte & Doolittle 1982, https://doi.org/10.1016/0022-2836(82)90515-0
- Small proteins and strains: Weaver 2019, https://doi.org/10.1128/mBio.02819-18; Grozdanov 2004, https://doi.org/10.1128/JB.186.16.5432-5441.2004; Touchon 2009, https://doi.org/10.1371/journal.pgen.1000344

---
## Appendix A. File templates (create in M0)

### A.1 `AGENTS.md`

````markdown
# AGENTS.md: how agents build, test and review GIPPA

## Source of truth
- GIPPA_SPEC.md (read-only). Accepted ADRs in docs/adr/ override it.
- Rules R1–R14 (GIPPA_SPEC.md §0.2) bind every agent and subagent.

## Setup (FarmShare; inside a Slurm allocation, never on a login node)
source scripts/farmshare_env.sh               # every shell and every sbatch script
bash scripts/bootstrap_farmshare.sh           # once: uv, Nextflow, nf-test, gh (user space)
uv sync --frozen --all-extras
uv run python scripts/probe_env.py            # writes .gippa/capabilities.json
sbatch scripts/slurm/build_envs.sbatch <tool> [...]   # micromamba envs (workflow/envs/*.yml) and Apptainer pulls
sbatch scripts/slurm/gippa_setup.sbatch        # gippa setup --bakta light,full --genomad --promotech-models --go --kegg --cog --contaminants

## Test commands
- every commit:   uv run pytest -m "not slow and not calibration_full and not e2e" -n auto
- calibration:    uv run pytest tests/calibration -m calibration -n auto      (full: -m calibration_full)
- snapshots:      uv run pytest tests/snapshot --mpl
- integration:    uv run pytest tests/integration ; nf-test test workflow/tests
- end-to-end:     sbatch scripts/slurm/weekly.sbatch    (never in the session allocation)
- pipeline run:   sbatch scripts/slurm/gippa_run.sbatch <args>   (-profile farmshare,apptainer)
- benchmarks:     uv run pytest benchmarks --benchmark-only
- lint and types: uv run ruff check . && uv run ruff format --check . && uv run mypy src/gippa

## Conventions
- src layout, type hints everywhere, pure functions for statistics; file I/O only in io/ and steps/.
- 0-based half-open coordinates internally; convert only in writers (spec §6.5).
- Randomness only through gippa.rng.for_step(step_name, seed) (spec §6.6).
- Never import GPL/AGPL/unlicensed tools (scripts/ci/check_imports.py).
- Commit trailer `Agent-Role: <role>`; conventional commits; PRs ≤ ~400 lines.
- Shared cluster (R14): jobs named `gippa-*`, ≤ 128 CPUs / 4 GPUs / 20 jobs, large files only under $GIPPA_ROOT.

## Definition of done (module): spec §14.9
## Where things live: specs docs/specs/, methods docs/methods/, ADRs docs/adr/, backlog docs/backlog.md, state PROGRESS.md
````

### A.2 `CLAUDE.md`

````markdown
@AGENTS.md

# Claude Code notes
- At session start: read PROGRESS.md, then .gippa/capabilities.json (run scripts/probe_env.py if missing).
- The main session is the lead; use the subagents in .claude/agents/ for the other roles (spec §14.1).
- You run on Stanford FarmShare inside a Slurm interactive allocation (spec §3). Source scripts/farmshare_env.sh.
- GitHub: `gh pr`, `gh issue`, `gh run` work (spec §3.10). Never push tags, never push to main.
- Commands longer than ~10 min or bigger than the allocation: sbatch from scripts/slurm/, record the job ID in PROGRESS.md, poll with squeue/sacct (spec §3.8).
- Never run work on a login node; cancel only GIPPA jobs by ID; respect the R14 caps.
- Unreachable source: log it under PROGRESS.md › Blocked and continue elsewhere; never use unofficial copies (R6).
- End every work unit: update PROGRESS.md, commit, push.
````

### A.3 `PROGRESS.md`

````markdown
# GIPPA progress
_Last updated: <YYYY-MM-DD> · session <link>_

## Current milestone
M<k>: <name>. Status: in progress | gate pending | done

## Done (with evidence)
- [M0] CI green at <sha>; toy genome simulated (tests/unit/test_sim_basic.py); PR #<n>

## In progress
- <branch>: <what>; next step: <…>
- Slurm jobs: <job ID> <job name> <purpose> <log path> <submitted at>

## Next actions (ordered)
1. …

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|

## Questions for Khoa
| # | Question | Options | Default if no answer | Asked |
|---|---|---|---|---|

## Environment notes
- capabilities: allocation cpus=…, mem=…, time left=…, $GIPPA_ROOT used=…, apptainer=…, micromamba=…, nextflow=…, gh=…, sbatch=…, gpu partition=…
- databases present: …

## Measured numbers
| Step | Genome | Wall | Peak RSS | Node class (Slurm job ID) |
|---|---|---|---|---|
````

### A.4 `.claude/settings.json`, hook and protected list

```json
{
  "env": {
    "BASH_DEFAULT_TIMEOUT_MS": "600000",
    "BASH_MAX_TIMEOUT_MS": "3600000"
  },
  "permissions": {
    "deny": [
      "Edit(/GIPPA_SPEC.md)",
      "Edit(/ci/gates.yaml)",
      "Edit(/ci/protected.txt)",
      "Edit(/tests/calibration/tolerances.yaml)"
    ]
  },
  "hooks": {
    "SessionStart": [
      { "matcher": "startup|resume",
        "hooks": [ { "type": "command", "command": "bash \"$CLAUDE_PROJECT_DIR\"/scripts/session_start.sh" } ] }
    ],
    "PreToolUse": [
      { "matcher": "Edit|Write|MultiEdit|NotebookEdit",
        "hooks": [ { "type": "command", "command": "python3 \"$CLAUDE_PROJECT_DIR\"/.claude/hooks/protect_paths.py" } ] }
    ]
  }
}
```

These `Edit(...)` rules cover every built-in file-editing tool. A leading `/` anchors the path at the project root. Do not write `Write(...)` rules: Claude Code ignores them. Shell commands can still bypass these rules, which is why CI also checks protected files (`scripts/ci/check_protected.py`).

`.claude/hooks/protect_paths.py`:

```python
#!/usr/bin/env python3
"""Block edits to protected files listed in ci/protected.txt (GIPPA_SPEC.md R3, §14.6)."""
import fnmatch, json, os, sys

event = json.load(sys.stdin)
tool_input = event.get("tool_input") or {}
path = tool_input.get("file_path") or tool_input.get("notebook_path") or ""
root = os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())
if path:
    rel = os.path.relpath(os.path.abspath(path), root)
    try:
        with open(os.path.join(root, "ci", "protected.txt")) as fh:
            globs = [l.split("#")[0].strip() for l in fh if l.split("#")[0].strip()]
    except FileNotFoundError:
        globs = []
    # New files may be created; existing protected files may not be changed.
    if os.path.exists(os.path.join(root, rel)) and any(fnmatch.fnmatch(rel, g) for g in globs):
        print(f"Blocked: {rel} is protected (GIPPA_SPEC.md R3). Draft an ADR and ask Khoa.", file=sys.stderr)
        sys.exit(2)
sys.exit(0)
```

`ci/protected.txt`:

```text
# Changing an existing protected file needs an ADR with status accepted-by-khoa and the
# commit trailer "Gate-Change-Approved: ADR-NNNN" (GIPPA_SPEC.md R3, §14.6). Creating them is allowed.
GIPPA_SPEC.md
ci/gates.yaml
ci/protected.txt
tests/calibration/tolerances.yaml
src/gippa/data/calibration.yaml
src/gippa/config/defaults.yaml        # adding new keys is allowed with an ADR; changing values needs Khoa
.claude/settings.json
.claude/hooks/*
```

`scripts/ci/check_protected.py` diffs the PR against its base:

- creating a protected file is allowed;
- modifying one requires the trailer, plus the referenced ADR with `Status: accepted-by-khoa`;
- for `defaults.yaml`, a diff that only adds keys passes with any ADR reference.

`ci/gates.yaml` (initial):

```yaml
coverage_min: 0.90
coverage_exclude: ["src/gippa/annotate/*_cmd.py"]       # thin command builders only
budgets_s: {unit: 60, property: 120, golden: 120, contract: 60, snapshot: 120,
            calibration_fast: 600, integration: 1800}
benchmark_slowdown_max: 0.15
pr_lines_warn: 400
pr_lines_fail: 800
required_capabilities:
  ci: []
  farmshare_nightly: [needs_apptainer, needs_conda, needs_nextflow, needs_slurm, needs_bakta_db, needs_bakta_full,
                      needs_genomad_db, needs_promotech_models]
  farmshare_weekly: [needs_apptainer, needs_conda, needs_nextflow, needs_slurm, needs_bakta_full]
farmshare: {max_cpus: 128, max_gpus: 4, max_jobs: 20}
```

`tests/calibration/tolerances.yaml` (initial):

```yaml
levels: [0.01, 0.05, 0.10]
se_multiplier: 2.0
ks_min_p: 0.001
detectability_slope: [0.8, 1.2]
denovo_precision_max_abs_bias: 0.05
mapping_recall_min: 0.999
f1_min: {variants_high: 0.90, novel_orfs: 0.90}
golden_r: {wallenius_rel_tol: 1.0e-6, fgsea_es_abs_tol: 1.0e-9, fgsea_nes_abs_tol: 0.02, fgsea_log10p_k: 3.0}
seeds: {GOLDEN_SEED: 50000, calibration_base: 20261014}
replicates_fast: {placement: 50, gene_fdr: 50, detectability: 20, regions: 100, autocorrelation: 100,
                  islands: 50, discordance: 100, missed_operons: 50, ora: 100, gsea: 50,
                  determinants: 20, novel_orfs: 50, denovo_precision: 50}
replicates_full: {placement: 300, gene_fdr: 300, detectability: 100, regions: 1000, autocorrelation: 1000,
                  islands: 500, discordance: 1000, missed_operons: 300, ora: 1000, gsea: 300,
                  determinants: 100, novel_orfs: 300, denovo_precision: 300}
```

### A.5 Subagents (`.claude/agents/*.md`)

Each file has YAML frontmatter followed by instructions. All of them inherit R1–R14 and read `AGENTS.md` first.

```markdown
---
name: spec-writer
description: Turns a GIPPA_SPEC.md section into an issue with input/expected-output examples, edge cases and acceptance criteria. Use at step 1 of the feature loop.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Write docs/specs/<issue>.md: the spec section, 3–10 concrete examples (input → expected output),
edge cases (circular/linear, minus strand, empty input, ties), and acceptance criteria that a test
can check. Do not write code or tests. Flag any gap in the spec as a question for the lead.
```

```markdown
---
name: test-author
description: Writes failing unit, property and calibration tests from an issue spec. Use at step 2 of the feature loop. Never edits src/ or workflow/.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Write tests only under tests/ (never tolerances.yaml). Tests must fail before the implementation
exists; run them and report the failure output. Prefer small examples from the spec and Appendix B;
add Hypothesis properties for invariants and calibration tests for every statistical claim (spec §8.4).
Commit with the trailer "Agent-Role: test-author".
```

```markdown
---
name: implementer
description: Writes the smallest code in src/ or workflow/ that makes the failing tests pass, then refactors. Use at step 3 of the feature loop.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Never modify existing tests, tolerances or gates. Keep functions pure and typed; follow spec §6.
Run the fast test suite before committing. Commit with the trailer "Agent-Role: implementer".
If a test looks wrong, report it to the lead instead of changing it.
```

```markdown
---
name: code-reviewer
description: Independent review of a PR for correctness, clarity, complexity and test quality. Use at step 4 of the feature loop.
tools: Read, Grep, Glob, Bash
---
Check the diff against the spec section and the issue. Look for coordinate off-by-one errors,
strand handling, circular wrap, determinism, error handling and missing edge-case tests. Run
mutation testing on changed stats/ or mapping/ code when it is feasible. Post the verdict
(approve / changes requested, with reasons) as a PR comment. Never edit the code under review.
```

```markdown
---
name: statistician
description: Independent statistical reviewer with blocking power. Reviews methods notes, simulations and calibration before code; owns attack nulls. Use for every PR labelled stats.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Read docs/methods/<method>.md and the calibration results first, then the code. Verify:
null models keep the data's structure (TU blocks, autocorrelation, detectability); multiple-testing
families match spec §8.1; every p/q has an effect size and CI; level/FDR/FWER pass with the
tolerances. Add attack nulls under tests/calibration/attacks/ (spec §8.5). Approve or block in a PR
comment. Never edit the implementation.
```

```markdown
---
name: performance
description: Profiles and speeds up hot paths without changing outputs. Use in the speed loop.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Profile with py-spy or scalene on the benchmark genomes; change only hot paths; golden outputs and
calibration results must be byte-identical or within tolerance. Record before/after in
benchmarks/ and propose a baseline update by ADR. Stop after 3 attempts with < 5 % gain.
```

```markdown
---
name: viz
description: Builds and inspects plots and the report against the style checklist (spec §10, §14.4). Use in M9 and the visualisation loop.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Render each plot from the simulated run, look at the PNG, and score it against the checklist.
Propose snapshot baselines; never approve them yourself (Khoa does at M9).
```

```markdown
---
name: validation
description: Runs end-to-end and pilot analyses and writes validation reports. Never changes code or thresholds.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Run the pilot (spec §15) and the weekly end-to-end runs as Slurm jobs (spec §3.8); write validation/pilot/REPORT.md with the
metrics of §15.5 and the pass criteria of §15.6. Never tune on the holdout; never edit src/.
```

```markdown
---
name: docs
description: Maintains the MkDocs user guide, methods index, tutorials and changelog.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Document every CLI option and output table from the code and schemas, not from memory. Include the
licence notes for third-party tools (geNOMAD non-commercial, KEGG academic use, GPL tools run as
separate programs).
```

### A.6 ADR template (`docs/adr/ADR-NNNN-<slug>.md`)

```markdown
# ADR-NNNN: <title>
- Status: proposed | accepted | accepted-by-khoa | superseded by ADR-MMMM
- Date: YYYY-MM-DD
- Spec sections: §…
- Needs Khoa: yes | no   (yes when a default, threshold, interface, dependency, protected file or §2 decision changes)

## Context
## Decision
## Evidence (before → after: metrics, tests, simulations)
## Consequences
```

ADR-0001 records §2.1 and §2.4 as accepted.

### A.7 Pull request template (`.github/pull_request_template.md`)

```markdown
## What and why
## Spec section(s)
## Issue / backlog item
## Roles (Agent-Role trailers) and which tests were red first
## Gates: coverage · calibration · golden · snapshot · benchmark
## ADR (if a default, threshold, interface or dependency changed)
- [ ] No protected file changed (or `Gate-Change-Approved: ADR-NNNN`)
- [ ] Changelog fragment added
- [ ] PROGRESS.md updated
```

---

## Appendix B. Worked numbers for tests

### B.1 Benjamini–Hochberg

p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]

q = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.105714, 0.216, 0.216, 0.216]

### B.2 Chance hits of a wrong peptide

E(L, m) = N Σ_{j≤m} C(L, j)(1 − s)^j s^(L−j), with N = 10⁷ and s = 0.069 (E. coli composition, I = L merged). The table in §7.2.5 shows these values rounded.

| L | exact | ≤ 1 | ≤ 2 |
|---|---|---|---|
| 6 | 1.1 | 88 | 3.0 × 10³ |
| 7 | 0.074 | 7.1 | 2.9 × 10² |
| 8 | 0.0051 | 0.56 | 27 |
| 9 | 3.5 × 10⁻⁴ | 0.043 | 2.4 |
| 10 | 2.4 × 10⁻⁵ | 0.0033 | 0.20 |
| 12 | 1.2 × 10⁻⁷ | 1.9 × 10⁻⁵ | 0.0014 |
| 14 | 5.5 × 10⁻¹⁰ | 1.1 × 10⁻⁷ | 9.3 × 10⁻⁶ |

A Monte Carlo simulation on random sequences with these frequencies MUST reproduce the ≤ 1 column for L = 8–10 within 20 %.

### B.3 MAYU closed form

One bin with N_b = 5,000 genes, n_b = 3,000 detected, F = 50 false peptides, uniform w = 1/5,000:

- p_b = 1 − (1 − 1/5000)⁵⁰ = 0.00995
- FP_b = 2000 × 0.00995 / 0.99005 = 20.10
- FDR = 20.10 / 3000 = 0.0067

### B.4 Autocorrelation and shifts

- Moran's I of the alternating series +1, −1, … (n = 10, circular, lag 1) is −1.
- Exhaustive circular shifts on a 100-gene circular replicon give K = 99 shifts. If the observed statistic is more extreme than all 99, p = (1 + 0)/(1 + 99) = 0.01.

### B.5 Masses for the mass-equivalence list (monoisotopic, Da)

**Residues.**

| Residue | Mass | Residue | Mass |
|---|---|---|---|
| G | 57.02146 | Q | 128.05858 |
| A | 71.03711 | K | 128.09496 |
| S | 87.03203 | E | 129.04259 |
| P | 97.05276 | M | 131.04049 |
| V | 99.06841 | H | 137.05891 |
| T | 101.04768 | F | 147.06841 |
| C | 103.00919 | R | 156.10111 |
| L, I | 113.08406 | Y | 163.06333 |
| N | 114.04293 | W | 186.07931 |
| D | 115.02694 | U | 150.95364 |
| | | O | 237.14773 |

**Modifications (Δ).**

| Modification | Δ | Modification | Δ |
|---|---|---|---|
| Oxidation | +15.994915 | Carbamyl | +43.005814 |
| Dioxidation | +31.989829 | Formyl | +27.994915 |
| Deamidated | +0.984016 | Phospho | +79.966331 |
| Methyl | +14.015650 | Amidated | −0.984016 |
| Dimethyl | +28.031300 | Dehydrated | −18.010565 |
| Trimethyl | +42.046950 | Ammonia-loss | −17.026549 |
| Acetyl | +42.010565 | Cation:Na | +21.981943 |
| Carbamidomethyl | +57.021464 | | |

**Checks the generated list MUST pass** at the default tolerance of 0.05 Da:

- N = GG (Δ 0);
- Q = AG (Δ 0);
- GT = AS (Δ 0);
- AV = GL (Δ 0);
- AN = GQ (Δ 0);
- W vs AD/EG (Δ 0.0153);
- R vs GV (Δ 0.0112);
- K vs Q (Δ 0.0364);
- A + Oxidation = S;
- N + Deamidated = D;
- G + Methyl = A;
- M + Oxidation vs F (Δ 0.0329).

### B.6 Placement q toy example

In stratum m = 1 with π̂ = 1, 150 targets precede the first decoy, so FDR = (0 + 1)/150 = 0.0067 and those 150 are accepted at 0.01. With only 50 targets and no decoys, FDR = 1/50 = 0.02, and none is accepted. The +1 makes small strata conservative by design.

### B.7 VCF record example

```text
#CHROM  POS   ID  REF  ALT  QUAL  FILTER  INFO                                                                                         FORMAT    pilot_A3
chr     1001  .   GCT  TCT  .     PASS    GENE=NISS_00012;AA_POS=45;AA_REF=A;AA_ALT=S;CALL=ALT;CONF=high;NT_CHANGES=1;ALT_TYPE=substitution;FLAGS=mass_equivalent|snv;CLASS_ERR=unestimated   GT:AD:DP  1:0,3:3
```

A→S (GCT → TCT) is one nucleotide change, and it is also explained by oxidation (+15.994915), hence `mass_equivalent`.
