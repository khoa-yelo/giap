# M2 "Annotation wrappers and benchmarks": issue index

Spec: GIAP_SPEC.md §13 row M2. Deliverables: "`workflow/envs/*`, Nextflow skeleton with stubs, wrappers
and parsers for every tool in §7.1.1; `giap setup`; recorded fixtures for the 100-kb genome; UniOP
adapter; TUs; promoter windows; anti-SD; codon usage; regions". **Gate:** "a contract test per tool;
the full Bakta DB downloaded and the full-DB fixtures recorded (G17); `nightly.sbatch` and
`weekly.sbatch` proposed to Khoa for `scrontab`; measured FarmShare runtimes per tool (with node
class) replace the estimates".

One issue = one PR. Size budget ≈ 400 **code** lines under `src/`, `workflow/`, `scripts/`
(ADR-0007: tests, fixtures, docs not counted). Estimates below are code lines only. Command builders
go in `src/giap/annotate/<tool>_cmd.py` (thin, coverage-excluded by `ci/gates.yaml`); input builders
and parsers go in `src/giap/annotate/<tool>.py` (full coverage). GPL/AGPL/unlicensed tools are only
ever run as processes (R5); `scripts/ci/check_imports.py` already forbids importing them.

## Issues

| # | Issue | Scope | Depends on | Code est. | Capabilities (§3.6) |
|---|---|---|---|---|---|
| 01 | [M2-01](M2-01-tool-envs-containers.md) | tool registry `annotate/tools.py`; `workflow/envs/*.yml` exact pins; `containers.config` digests; Dockerfiles for git-pinned tools; `scripts/build_envs.sh` + `slurm/build_envs.sbatch`; NOTICE rows | M0 | 380 | build: `needs_slurm`, `needs_conda`, `needs_apptainer`, `needs_net`; unit: none |
| 02 | [M2-02](M2-02-step-registry-precomputed.md) | step registry + `giap _step`; local executor (§6.4); `annotation.precomputed` formats + loaders | M0, M1-04 | 330 | none |
| 03 | [M2-03](M2-03-gene-model-bakta-json.md) | extends M3-01's `model/genes.py`; defensive Bakta 1.12.2 JSON parser; `bakta_cmd.py` | M1-01, M1-02; shares file with M3-01 | 330 | none (contract: recorded JSON, `pending_sample` until M2-10) |
| 04 | [M2-04](M2-04-parse-bakta-step.md) | `parse_bakta` step: `genes.parquet`, `features.parquet`, `rrna16s.fa`, `trna_counts.tsv`, `xrefs.parquet`; anti-SD; ribosomal set; TE6 | 02, 03 | 300 | none (TE6 on recorded K-12) |
| 05 | [M2-05](M2-05-simulator-v2.md) | `giapsim` additions on top of M3-18: signals (RBS, σ70, terminators), CRISPR, oriC, prophage/BGC blocks, real rRNA/tRNA option, precomputed layers, annotation truth, `SIM100K` | 02 (formats), 03 (field names), M3-18 | 400 (tooling) | none |
| 06 | [M2-06](M2-06-nextflow-skeleton.md) | `main.nf`, `nextflow.config`, `conf/*`, profiles, labels, `GIAP_STEP` process, `runner/nextflow.py` launcher, stub nf-test | 02 | 330 | `needs_nextflow` |
| 07 | [M2-07](M2-07-nextflow-tool-processes.md) | one process per external tool (10) with `stub:` + `versions.yml` + `meta.yml`; `subworkflows/local/annotate.nf`; precomputed bypass; module switches | 01, 06 | 400 | `needs_nextflow` (stub); real runs in 21 |
| 08 | [M2-08](M2-08-giap-setup-core.md) | `giap setup` core: items, `db_manifest.json`, `--check`, exit 2; Bakta light/full, geNOMAD, antiSMASH, Promotech models; `slurm/giap_setup.sbatch` | 01 | 360 | unit: none (fake runner); live: `needs_slurm`, `needs_net`, `needs_conda` |
| 09 | [M2-09](M2-09-giap-setup-resources.md) | `giap setup --go --kegg --cog --contaminants` | 08 | 220 | unit: none; live: `needs_net` |
| 10 | [M2-10](M2-10-recording-harness-wave-a.md) | recording harness (`scripts/record_fixture.py`, `slurm/record_fixtures.sbatch`, `slurm/record_bakta_full.sbatch`); genomes K-12/Nissle/sim100k committed; wave-A fixtures (Bakta light+full, TerminatorNet, geNOMAD, antiSMASH, BacTermFinder) | 01, 03 (cmd), 05, 08 | 300 | `needs_slurm`, `needs_conda`/`needs_apptainer`, `needs_bakta_db`, `needs_bakta_full`, `needs_genomad_db`, `needs_net` |
| 11 | [M2-11](M2-11-user-annotation-id-map.md) | `--annotation` GFF3/GBFF → Bakta `--regions` file and `annotation/id_map.tsv` | 03 | 240 | none |
| 12 | [M2-12](M2-12-uniop-pairs-tus.md) | UniOP adapter (Prodigal GFF ↔ locus tags), distance-model fallback, `operon_pairs.tsv`, TUs, TU upstream regions, `calibration.yaml` | 04 | 380 | none (UniOP fixture via harness) |
| 13 | [M2-13](M2-13-promoter-windows-promotech.md) | 40-nt windows (221/start, wrap via `replicons`, dedupe, 100k chunks in `windows/`), Promotech parser, per-gene promoter features, `promoter_between`, `promoters.bed` | 04, 12 | 380 | none (fixture via harness, `needs_promotech_models` to record) |
| 14 | [M2-14](M2-14-promoter-calculator.md) | one-region-per-call Promoter Calculator, parser, `tx_rates`, TU features | 12 | 200 | none |
| 15 | [M2-15](M2-15-terminators.md) | TerminatorNet (`-s`), TransTermHP (`.coords`, `.tt`), cross-check, BacTermFinder (opt-in), `terminator_between`, `terminators.bed` | 04 | 390 | none |
| 16 | [M2-16](M2-16-ostir-rbs.md) | OSTIR input (−50/+50, start = end = 51), parser → RBS features | 04 | 180 | none |
| 17 | [M2-17](M2-17-codon-usage.md) | codon-bias in-process: CAI, tAI, ENC, Fop, GC3, CPB; `at_2_10`, `rare_codons_2_10` | 04 | 220 | none |
| 18 | [M2-18](M2-18-genomad-antismash.md) | geNOMAD and antiSMASH command builders + parsers | 03 | 300 | none (fixtures from 10) |
| 19 | [M2-19](M2-19-regions.md) | `regions.tsv`/`.bed`: prophage (+ function subsets), plasmid, bgc, crispr (+cas), amr; 50 % membership | 04, 18 | 350 | none |
| 20 | [M2-20](M2-20-context-layers.md) | `context_layers` step: `layers.parquet`, `annotation/genes.gff3`, all `annotation/*` files, `tools/` natives, `layer_missing` | 02, 04, 11–19 | 380 | none |
| 21 | [M2-21](M2-21-integration-nightly-runtimes.md) | wave-B recordings, real-tool integration on sim100k, local ≡ Nextflow, `nightly.sbatch` / `weekly.sbatch`, measured runtimes | 06, 07, 10, 20 | 300 | all: `needs_slurm`, `needs_nextflow`, `needs_apptainer`, `needs_conda`, DB markers |

Total ≈ 6,700 code lines, 21 PRs.

## Dependency DAG

```text
M0/M1 ──┬─► 01 envs ──┬─► 08 setup core ──► 09 setup resources
        │             ├─► 07 NF tool processes ◄── 06 NF skeleton ◄── 02 registry
        │             └─► 10 recording (wave A) ◄── 05 sim v2 ◄── {02, 03, M3-18}
        ├─► 02 registry/_step/local/precomputed ──► 04
        └─► 03 gene model + Bakta JSON ──┬─► 04 parse_bakta ──┬─► 12 UniOP/pairs/TUs ──┬─► 13 promoters
                                         │                    │                        └─► 14 Promoter Calculator
                                         │                    ├─► 15 terminators
                                         │                    ├─► 16 OSTIR
                                         │                    ├─► 17 codon usage
                                         │                    └─► 19 regions ◄── 18 geNOMAD/antiSMASH ◄── 03
                                         └─► 11 user annotation / id_map
{02, 04, 11–19} ──► 20 context_layers ──► 21 integration + nightly/weekly + runtimes ◄── {06, 07, 10}

03 shares src/giap/model/genes.py with M3-01; 05 builds on M3-18 (see "Overlap with the M3 specs").
```

## Parallel tracks (separate agents, separate worktrees under `local/worktrees/`)

| Wave | Track | Issues (in order) | Notes |
|---|---|---|---|
| 1 | **T-ENV** (infra) | 01 → 08 → 09 | needs Slurm + network for the live parts; submit `build_envs.sbatch` early, it is long-latency |
| 1 | **T-WF** (workflow) | 02 → 06 → 07 | Groovy + registry; stub mode only, no tools needed |
| 1 | **T-MODEL** (gene model) | 03 → 04 | pure Python on the M0 sim JSON + hand-written JSON; contract test `pending_sample` until 10 |
| 1→2 | **T-SIM** | 05 (after 02's format table, 03's field list and M3-18 are merged) | tooling; unblocks recording |
| 2 | **T-REC** (lead or one agent with Slurm) | 10, then wave-B recordings as 12–16 land | the only track that needs the DBs; keep ≤ 20 jobs (R14) |
| 3 | **T-OPERON** | 12 → 13 → 14 | promoters need TUs |
| 3 | **T-TERM** | 15 | independent after 04 |
| 3 | **T-SEQ** | 16, 17 | two small PRs, either order |
| 3 | **T-REGION** | 18 → 19 | 18 can start right after 03 (wave 2) |
| 3 | **T-USER** | 11 | can start right after 03 (wave 2) |
| 4 | **T-INT** | 20 → 21 | 21 also proposes `nightly.sbatch` / `weekly.sbatch` to Khoa |

Maximum useful parallelism: 4 agents in wave 1, 6 in wave 3. Every issue in tracks T-MODEL,
T-OPERON, T-TERM, T-SEQ, T-REGION, T-USER and T-INT(20) is developed **against recorded fixtures,
hand-written doc-derived mini fixtures, or the simulator only**; none of them needs Slurm, a
database, a conda env or Apptainer to go green in CI. Their contract tests on real recorded output
are marked `pending_sample` until the fixture lands, then switched on in the recording PR.

## Issues that need cluster capabilities

| Issue | What needs it | Marker(s) | Where it runs |
|---|---|---|---|
| 01 | building envs, pulling images, resolving digests | `needs_net`, `needs_conda`, `needs_apptainer`, `needs_slurm` | `sbatch scripts/slurm/build_envs.sbatch` |
| 06, 07 | nf-test in stub mode | `needs_nextflow` | session allocation, and CI (`ci.yml` installs Nextflow) |
| 08, 09 | live downloads (84 GB full Bakta DB) | `needs_net`, `needs_conda`, `needs_slurm` | `sbatch scripts/slurm/giap_setup.sbatch` |
| 10 | recording wave A | `needs_bakta_db`, `needs_bakta_full`, `needs_genomad_db`, `needs_conda`/`needs_apptainer`, `needs_slurm` | `record_bakta_full.sbatch` (16 CPU, 64 GB), `record_fixtures.sbatch` |
| 12–16 | recording each tool's wave-B fixture (acceptance item; can be done by T-REC) | `needs_conda`, `needs_promotech_models` (13) | `record_fixtures.sbatch TOOLS=<tool>` |
| 21 | real-tool integration, nightly/weekly, runtimes | all of the above + `needs_nextflow` | Slurm only |

The Bakta DB prefetch (`scripts/slurm/fetch_bakta_db.sbatch`, job 1845787) already places
`local/db/bakta/db-light` and `local/db/bakta/db`; M2-08 adopts them (checksum + manifest) instead of
downloading again.

## Annotation step contract (single source for 02, 06, 07 and 20)

Work-directory file names. "py" steps run through `giap _step`; "ext" steps are Nextflow processes.
A `precomputed` key replaces the listed steps and supplies the listed output (formats in M2-02).

| Step | Kind | Inputs | Outputs | Replaced by `precomputed.` | Issue |
|---|---|---|---|---|---|
| `bakta` | ext | `genome.norm.fa`, `replicons.tsv` [, `regions.gff3`] | `bakta/<P>.{json,gff3,gbff,faa,tsv}` | `bakta_json` | 03, 07 |
| `parse_bakta` | py | Bakta JSON, `genome.norm.fa` | `genes.parquet`, `features.parquet`, `replicons.parquet`, `rrna16s.fa`, `trna_counts.tsv`, `xrefs.parquet`, `annotation_meta.json` | — | 04 |
| `uniop_input` | py | `genes.parquet` | `uniop_cds.gff`, `uniop_map.tsv` | `operon_pairs` | 12 |
| `uniop` | ext | `uniop_cds.gff` | `uniop.pred` | `operon_pairs` | 07, 12 |
| `operon_pairs` | py | `uniop.pred` (optional), `genes.parquet`, `uniop_map.tsv` | `pairs_raw.parquet`, `tus.parquet`, `tu_upstream/<tu_id>.fa` | `operon_pairs` (pairs only; TUs still built) | 12 |
| `promoter_windows` | py | `genes.parquet`, `genome.norm.fa` | `windows_NNN.fa`, `windows_map.parquet` | `promoters` | 13 |
| `promotech` | ext (per chunk) | `windows_NNN.fa` | `promotech_NNN/` (native) | `promoters` | 07, 13 |
| `promoter_layer` | py | `windows_map.parquet`, `promotech_*` | `promoter_windows.parquet` | `promoters` | 13 |
| `promoter_calculator` | ext (per TU) | `tu_upstream/<tu_id>.fa` | `pc/<tu_id>.csv` | `tx_rates` | 07, 14 |
| `tx_layer` | py | `pc/*.csv`, `tus.parquet` | `tx_rates.parquet` | `tx_rates` | 14 |
| `terminatornet` | ext | `genome.norm.fa` | `terminatornet.tsv` | `terminators` | 07, 15 |
| `transterm_input` | py | `genes.parquet` | `genes.coords` | `terminators` | 15 |
| `transtermhp` | ext | `genome.norm.fa`, `genes.coords` | `transtermhp.tt` | `terminators` | 07, 15 |
| `bactermfinder` | ext (opt-in) | `genome.norm.fa` | `bactermfinder/` | `terminators` | 07, 15 |
| `terminator_layer` | py | the three above | `terminators.parquet`, `terminators_crosscheck.tsv` | `terminators` | 15 |
| `ostir_input` | py | `genes.parquet`, `annotation_meta.json`, genome, replicons | `ostir_input.csv` | `rbs` | 16 |
| `ostir` | ext | `ostir_input.csv` | `ostir.csv` | `rbs` | 07, 16 |
| `rbs_layer` | py | `ostir.csv`, `ostir_input.csv` | `rbs.parquet` | `rbs` | 16 |
| `codon_usage` | py | `genes.parquet`, `trna_counts.tsv` | `codon.parquet` | — (always in-process) | 17 |
| `genomad` | ext [mobile] | `genome.norm.fa` | `genomad/` | `regions` | 07, 18 |
| `antismash` | ext [bgc] | Bakta GBFF | `antismash/` | `regions` | 07, 18 |
| `region_layer` | py | `genes.parquet`, `features.parquet`, `replicons.parquet`, `genomad/`, `antismash/` (optional) | `regions.parquet`, `region_genes.parquet` | `regions` | 19 |
| `context_layers` | py | all `*.parquet` layers above | `annotation/*` (§9.1), `layers.parquet` | — | 20 |

New py step names (`uniop_input`, `promoter_layer`, `tx_layer`, `transterm_input`,
`terminator_layer`, `ostir_input`, `rbs_layer`, `region_layer`) refine §6.2's rows without changing
any output; see Questions.

## Overlap with the M3 specs already written

- **M3-01** defines `src/giap/model/genes.py` (`Gene`, `NcFeature`, `GeneModel`,
  `gene_model_from_bakta_json`) and requires M2's `parse_bakta` to return the same `GeneModel`.
  M2-03 extends that module (extra optional fields, full parser) instead of creating a second gene
  model; whichever PR merges second rebases. Assign M2-03 and M3-01 to the same agent, or merge M3-01 first.
- **M3-18** adds plasmids, rRNA/tRNA, pseudogenes, sORFs and origin-crossing genes to `giapsim`
  (`features.py`). M2-05 adds only the annotation-specific parts (`signals.py`, `blocks.py`,
  `layers.py`) and reuses M3-18's `SimConfig` fields. Merge M3-18 before M2-05.
- **M4-21** checks M2's `parse_bakta`, TU and region outputs against the simulator's truth; the
  TU rule (singletons are TUs) and region membership in M2-12/M2-19 match what M4-06 assumes.

## Decisions taken in these specs (for the lead to confirm)

1. A precomputed layer = the **output of the layer's parse step** (one TSV format per key, M2-02),
   so tools and injection share one downstream code path.
2. Nextflow tool processes call the tool directly (GIAP is not installed in tool images); each
   `*_cmd.py` builder is kept in sync with its `main.nf` by a test that checks every flag.
3. Recording runs tools **directly in their micromamba envs** through `scripts/record_fixture.py`,
   not through Nextflow, so recordings do not wait for the workflow (wave A in 10, wave B with 12–16).
4. miniprot and BLAST+ get env files and container pins in M2-01; their wrappers, parsers and
   recordings move to M3 (they belong to §7.2).
5. `calibration.yaml` is created in M2-12 with spec placeholders (`d0 = 50`, `s = 20`,
   `terminatornet_yes = 0.5`, all `uncalibrated`); creating a protected file is allowed, changing it
   later needs an ADR (M6).

## Key questions for the lead (details in each issue)

- M2-01 Q1/Q3: git-pinned tools (UniOP, Promotech, TerminatorNet, BacTermFinder) pinned to the
  default-branch HEAD on the build date (ADR); `apptainer` profile uses conda for tools without a
  GHCR image until `containers.yml` publishes them.
- M2-02 Q1: precomputed TSV formats (table in M2-02); local executor missing layer → exit 1.
- M2-10 Q1/Q2: K-12/Nissle full-DB Bakta JSON may exceed 5 MB gzipped → drop `sequences[].nt`
  (genome committed separately); commit JSON + GFF3 for real genomes, GBFF only for sim100k.
- M2-13 Q1: Promotech headers become `w<index>` + sidecar map (dedupe makes `contig|strand|end3p`
  ambiguous); window `k = 0` ends right before the start codon.
- M2-19 Q1–Q3: AMR rows one per contig with NA coordinates; geNOMAD plasmid score cut-off 0.7;
  CRISPR region extended to include flanking `cas` genes.
- M2-04 Q2: `annotation_source` from the Bakta DB type; `user` only when unknown.

## Lead decisions on the questions (2026-10-14)
1. Source-only tools (UniOP, Promotech, TerminatorNet, BacTermFinder) are pinned to the default-branch
   commit on the build date and recorded in an ADR; tools without a published image run from conda
   until `containers.yml` pushes GHCR images. Accepted.
2. Precomputed layers use the M2-02 TSV formats; a missing layer under the local executor exits 1. Accepted.
3. The eight extra step names are internal; no ADR. Accepted.
4. Fixture trimming (drop `sequences[].nt` from recorded Bakta JSON, trim antiSMASH JSON, GBFF only
   for sim100k) to respect R11's 5 MB limit. Accepted.
5. Promotech window headers `w<index>` plus a map file (deviation from §7.1.1's `contig|strand|end3p`,
   which §7.1.1 itself warns the tool may mangle); recorded in an ADR. Accepted.
6. **Changed:** no new geNOMAD plasmid cut-off key (Khoa declined new thresholds on 2026-10-14). Use
   geNOMAD's own default-filtered `*_plasmid_summary.tsv`. AMR regions: one row per contig with NA
   coordinates, left out of BED. CRISPR regions widened to nearby `cas` genes. Accepted.
7. `annotation_source` from the Bakta DB type; `user` only when unknown. Accepted.
8. New table columns defined in M2-11/M2-15; `check_imports.py` gains terminatornet/bactermfinder;
   miniprot/tblastn parsers live in M3 (envs only in M2). Accepted.
Shared foundations with M3: M3-01 (`giap.model.genes`) and M3-18 (simulator features) are built first;
M2-03 and M2-05 extend them.
- From the M4 design: M2-20's layers table also produces `gc_local_dev` (local GC deviation per gene),
  which M4 consumes; `cog_category` is mapped from Bakta COG cross-references via the COG 2020 table.
