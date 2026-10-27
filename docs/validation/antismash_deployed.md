# antiSMASH in a deployed run (2026-10-26/27)

Khoa asked us to check that antiSMASH runs in the deployed pipeline (`-profile farmshare,apptainer`). In every WT pilot pass it was skipped because the precomputed Bakta JSON had no GBFF next to it (`layer_missing:antismash`). This report covers two real runs on E. coli Nissle 1917 WT (`NZ_CP058217.1` plus pMUT1 and pMUT2, GCF_021559835.1). This genome is the positive control because it carries the pks (colibactin) island.

**Result: both cases PASS.** antiSMASH 8.0.4 ran from the pinned image and used the `giap setup` databases. It found 6 BGC regions, one of which holds `clbA`–`clbS`. The region layer parsed them, and `layers_missing` is `[]` in both run manifests.

## Runs

Both runs used code from main at f1ff5b2f (worktree `local/worktrees/val-antismash`, `GIAP_REPO` exported) and `scripts/slurm/giap_run.sbatch` (`--profile farmshare,apptainer`). Common arguments were `--genome $GIAP_ROOT/pilot/pass3/genome/nissle1917_wt.fna.gz --peptides $GIAP_ROOT/pilot/inputs/a1.tsv.gz --threads 12 --memory 48.GB`, with the default modules (core, mobile, bgc, function). A launch-directory `nextflow.config` in each outdir set `executor.queueSize = 2` (R14) and `process.clusterOptions = '--qos=normal --exclude=rye-02'`. Outputs are under `$GIAP_ROOT/validation/antismash_deployed/`.

| case | input | head job(s) | result | antiSMASH task |
|---|---|---|---|---|
| (a) Bakta in the pipeline | genome only; no `--db-dir` (uses `$GIAP_DB_DIR`); Bakta full db 6.0 | 1897203 (8 h 16 min, wheat-07) | COMPLETED 0:0, 79 tasks, 0 failed | 1897622, barley-04 |
| (b) precomputed JSON + GBFF | `--annotation .../recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.json.gz` (the GBFF `nissle1917_wt.gbff.gz` sits next to it); `--db-dir $GIAP_DB_DIR` | 1897099, then a resume as 1898463 | 1897099 FAILED at PROMOTER_CALCULATOR (a stalled node, see D3), after ANTISMASH and REGION_LAYER had completed; the resume 1898463 COMPLETED 0:0 (38 run, 40 cached) | 1897161, barley-02 |

## Checks

- **Databases.** `$GIAP_DB_DIR/antismash` (4.1 GB on disk, `.giap-complete`) is in `db_manifest.json` as `antismash 8.0.4`, `ok`, from dl.secondarymetabolites.org. Run inside the image, `antismash --check-prereqs --databases $GIAP_DB_DIR/antismash` printed "All prerequisites satisfied".
- **Container.** The digest-pinned `antismash/standalone@sha256:d3f390e6…` resolved from `$GIAP_ROOT/containers` (trace `container` column). Conda is off for this process (`container_common.config`).
- **Databases path.** `.command.sh` contains `antismash --genefinding-tool none --minimal --databases /home/classes/bios/270/khoa/giap/local/db/antismash --cpus 12 …`. The directory is bound read-only by `apptainer.runOptions`. In case (a), `db_dir` came from the `$GIAP_DB_DIR` fallback.
- **GBFF wiring.** In case (a), the GBFF came from `BAKTA.out.gbff`. In case (b), it was found next to `x.json.gz` as `x.gbff.gz` and decompressed by the module.
- **Parsing.** antiSMASH wrote exactly one top-level `genome.json`. All record ids matched the genome's replicons.

## Regions found (identical in both cases)

`regions.parquet` is byte-identical between the two runs (sha256 `4d3331d8…`). Its 6 `bgc` rows, all from source `antiSMASH`:

| region | contig:start-end | products | genes |
|---|---|---|---|
| bgc_1 | NZ_CP058217.1:880584-912999 | NI-siderophore | 33 |
| **bgc_2** | NZ_CP058217.1:1948603-2082184 | **NRPS, T1PKS, NRP-metallophore** | 84 |
| bgc_3 | NZ_CP058217.1:3082441-3126029 | arylpolyene | 44 |
| bgc_4 | NZ_CP058217.1:3232161-3258454 | azole-containing-RiPP | 20 |
| bgc_5 | NZ_CP058217.1:3524250-3579236 | NRP-metallophore, NRPS (enterobactin, `entA`–`entF`) | 48 |
| bgc_6 | NZ_CP058217.1:3689877-3710777 | terpene-precursor | 23 |

bgc_2 is the pks island. It contains `clbA`–`clbS` (all 19 genes), and its `core_biosynthetic` genes are `clbB clbH clbI clbJ clbK clbN irp1 irp2 ybtS`. antiSMASH merged it with the adjacent yersiniabactin HPI (`irp1/2`, `ybt*`, `fyuA`) into one 134 kb region. The BGC regions reach context stats and the report (region O − E rows for `bgc_2` and the others).

## Cost per antiSMASH task

| task | Slurm request (allocated) | realtime | %CPU | peak RSS |
|---|---|---|---|---|
| (a) 1897622 | cpu=12, mem=48G (14 CPUs) | 25.7 s | 151 % | 3.3 GB |
| (b) 1897161 | cpu=12, mem=48G (14 CPUs) | 25.9 s | 155 % | 3.2 GB |

`--minimal` on a 5 Mb genome needs about 1.5 CPUs and 3.3 GB for under 30 s, so the request is about 9× too many CPUs and 15× too much memory (D2). In case (a), Bakta (full db) took 21 min on 12 CPUs (14 allocated) at 148 % CPU, with 7.9 GB peak RSS.

## Defects (none blocks antiSMASH; routed to the lead)

**D1: the antiSMASH version is lost from the run manifest (provenance).**
- **Failure:** the head log of both runs shows `WARNING giap.io.manifest_sources: unreadable versions.yml skipped: …/case_b/work/78/62ed30…/versions.yml (ScannerError)`. `run_manifest.json` `tools` therefore has no `antismash` entry.
- **Cause:** `workflow/modules/local/antismash/main.nf:25` runs `antismash --version 2>&1`. On FarmShare, `$HOME` (`/home/users/…`) is read-only inside the container, so matplotlib prints a two-line `MPLCONFIGDIR` warning on stderr. That warning goes into the YAML value, and `src/giap/io/manifest_sources.py:174-180` drops the whole file.
- **Minimal fix:** capture only stdout and the version token, e.g. `antismash --version 2>/dev/null | sed 's/^antiSMASH //'`. Optionally also `export MPLCONFIGDIR=$TMPDIR` for this process, which removes the warning from `.command.err` too.
- **Tests:** a module-script test that the ANTISMASH `versions.yml` is valid YAML when the tool writes noise to stderr (a stub `antismash` that writes a warning to stderr). A `tool_reports` test that a valid antiSMASH `versions.yml` yields `{"name": "antismash", "version": "8.0.4"}`. The other modules that use `--version 2>&1` (bakta, genomad, miniprot, tblastn) are worth the same check, though bakta's was clean here.

**D2: ANTISMASH is over-provisioned, and its own memory directive is dead.**
- **Cause:** `workflow/modules/local/antismash/main.nf:4-5` sets `label 'process_high'` and `memory { 32.GB * task.attempt }`. `workflow/conf/base.config:19` (`withLabel: 'process_high'`, 48 GB) wins over the process directive, so the task asks for 48 GB and 12–16 CPUs to do a job that needs about 1.5 CPUs and 3.3 GB.
- **Minimal fix:** a `withName: 'ANTISMASH'` entry in `base.config` (e.g. `cpus = 4; memory = { 8.GB * task.attempt }`), or `process_medium` plus a withName memory, and drop the dead directive.
- **Test:** a config test (as for the other withName entries) that ANTISMASH resolves to about 8 GB on attempt 1 under `farmshare`.

**D3: Slurm time-limit exit 140 is not retried, and barley-01 stalls like rye-02.**
- **Failure:** in run 1897099, `ANNOTATE:PROMOTER_CALCULATOR (pc_TU_406_TU_496)` (Slurm 1897804, barley-01) hit its 2 h limit after 2 min 50 s of CPU, with exit 140. Log: `$GIAP_ROOT/validation/antismash_deployed/case_b/work/67/db2227fe53c9e9e4eeef6050f2c79a/.command.log`. The run terminated. On resume, away from barley-01, the same task took 19 min 38 s.
- **Also on barley-01:** case (a) `PARSE_BAKTA` (1897625) took 20 min at 1 % CPU, against 7 s in case (b).
- **Cause:** `workflow/conf/base.config:12` retries only on 137 and 143, but Slurm's time-limit kill shows up as 140 here.
- **Minimal fix:** treat 140 like 143 in `errorStrategy` (nf-core retries 130–145). Add barley-01 to the slow-node stopgap list until the slow-node guard lands.
- **Test:** an errorStrategy config test that exit 140 with `task.time < resourceLimits.time` retries.

**D4 (minor, R14 accounting): Slurm adds CPUs to some requests.**
- **What happens:** `farmshare.config` computes CPUs as `ceil(mem / 4 GB)`, with 4 GB = 4096 MB, but `MaxMemPerCPU` is 4000 MB. A request of `cpu=12, mem=48G` therefore became 13 CPUs per task and 14 allocated (BAKTA 1897521, ANTISMASH 1897161/1897622).
- **Effect:** the queueSize × cpus budget in `farmshare.config:2` can undercount. At `--threads 16`, a 64 GB retry would also take 17–18 CPUs.
- **Fix:** use 4000 MB (`3.9.GB`) as the divisor, or clamp memory to `cpus × 4000 MB`.

## For the H1 sensitivity run: official GBFFs exist

- **pass-3 WT JSON** (`$GIAP_ROOT/pilot/pass3/inputs/nissle1917_wt.bakta.json`, Bakta 1.12.2, db 6.0 full, run 2026-10-25 20:55:44). It is the same Bakta run as `$GIAP_ROOT/pilot/pass3/bakta/nissle1917_wt/bakta/nissle1917_wt.json`, with only `sequences[].nt` stripped (`runs/pilot/pass3/jobs/prep_db.sh`). The GBFF from that run is `…/pass3/bakta/nissle1917_wt/bakta/nissle1917_wt.gbff`. Its gzip copy, `…/pass3/recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.gbff.gz`, is byte-identical after decompression and already sits next to `nissle1917_wt.json.gz`.
- **Cheapest faithful input:** pass `--annotation $GIAP_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.json.gz`, which is exactly case (b) above. No Bakta rerun is needed.
- **pass-1/2 JSON** (`$GIAP_ROOT/pilot/inputs/nissle1917.bakta.json`, CP171242–5, run 2026-10-15 09:47:23) matches `$GIAP_ROOT/recordings/bakta_full_circ/nissle1917/bakta/nissle1917.json` apart from `nt`. Its GBFF is `…/bakta_full_circ/nissle1917/bakta/nissle1917.gbff`. To use it, put a link or copy named `nissle1917.bakta.gbff` next to the pilot JSON, or pass the recording's JSON.
- **If a rerun were ever needed:** case (a) shows that Bakta 1.12.2 / db 6.0 full, run in the pipeline, gives byte-identical regions. That run cost 21 min on 12 CPUs.
