<!-- reproduce: status=draft freeze_commit=pending a3_branch=pending -->
# Pilot reproduction runbook (GIPPA_SPEC.md §15.2 steps 4-9; M10-14)

This is the recorded command (§3.9) for every step of the pilot, in execution order: the inputs
it reads, the job it runs as and the job that ran it. The format is fixed by
docs/specs/M10-14-pilot-reproducibility.md ("The runbook format"), and
`validation/pilot/reproduce_check.py static` checks it on every commit
(`tests/unit/test_m10_14_reproduce.py`). **Status: draft (14a).** Steps S1-S7 carry their pass-1
evidence. The freeze steps stay `pending` until the lead declares the freeze window and
writes the freeze commit F into the header (14b).

Every command runs on FarmShare with the environment sourced, from the main checkout:

```bash
source "$GIPPA_REPO/scripts/farmshare_env.sh"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
P=$GIPPA_ROOT/pilot/inputs
F7=$(git -C "$GIPPA_REPO" rev-parse --short=7 HEAD)
R=$GIPPA_ROOT/runs/pilot/freeze/$F7
```

The job classes are the scripts' `#SBATCH` headers unless the command overrides them. No step
runs on a login node, and none runs `uv sync` or `pip`.

## Inputs

A `.d` run has no single hash, so each run is listed twice: through its archive (the hash in the
committed manifest; the archives were deleted after unpacking, so `verify` reports them `absent`)
and through its `analysis.tdf` (M10-02 `provenance.json`). `analysis.tdf_bin` is not hashed (L11).

| id | path | sha256 | recorded_in |
|---|---|---|---|
| genome | tests/data/genomes/nissle1917.fna.gz | fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee | tests/data/genomes/PROVENANCE.md;tests/data/recorded/bakta/1.12.2-full/nissle/PROVENANCE.md |
| bakta_json | tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz | a012798400eaa71587dc7954473a201ef046b4907bdfa31936a2927a1cb6c97c | validation/pilot/db_a2/makedb_manifest.json;tests/data/recorded/bakta/1.12.2-full/nissle/PROVENANCE.md |
| strain_db | validation/pilot/db_a2/strain_db.faa.gz | 80dd01f3eeacf209131c7e14ebc65c12174aad89d9a9fd5a285704843faa887e | $GIPPA_ROOT/pilot/A2/provenance.json |
| cfg_a1 | validation/pilot/diann_A1.cfg | 0546a8905bcf092675991a793a45ba46f8e93a875834546587c07cc7d97ae608 | $GIPPA_ROOT/pilot/A1/provenance.json |
| cfg_a2 | validation/pilot/diann_A2.cfg | c0d145e87e582e450295517ac393d2c86093587605dd3da66c6f6a923c1f98e4 | $GIPPA_ROOT/pilot/A2/provenance.json |
| acquisition | validation/pilot/inputs/acquisition.json | a7af58ae2a31392f199dbca2259b899ab75f4a613ad020005ccf322e286bacb3 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| proteome_a1 | $GIPPA_DATA_DIR/proteomes/UP000000625.fasta | a8052a607a98d2444c9d5dd685bd1370150df8c5570c3116811544bbf3d093a2 | $GIPPA_ROOT/pilot/A1/provenance.json |
| contaminants | $GIPPA_DB_DIR/contaminants/contaminants.fasta | c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| diann_image | $GIPPA_ROOT/images/diann_2.3.1.sif | 1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381 | validation/pilot/diann_image.sha256;validation/pilot/README.md |
| px_record | $GIPPA_DATA_DIR/pxd078459/record/PXD078459.xml | 5d3f16c3137fb11ecae7f79ad02063d8e3008b21380784fcfc8aa77b7fb533f9 | validation/pilot/README.md;$GIPPA_DATA_DIR/pxd078459/record/PXD078459.xml.info |
| wt1_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914689_WT_P1_1_Y-D7_1_13119.d.zip | 58d09a41d6642f35f0cb25e88e4a715fbfe0c30f7b4f839c9fd30601bd99f013 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt2_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914690_WT_P2_2_Y-D8_1_13120.d.zip | 05d6d6499ac7ec3cf1af10d7eeabc28d741e826d81f4768667386828554e3102 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt3_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914691_WT_P3_3_Y-D9_1_13121.d.zip | c45618952c709de617934ca63eb29b9329286f41397c00341671056cdda1cd11 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt1_tdf | $GIPPA_DATA_DIR/pxd078459/WT_1.d/analysis.tdf | f4d5fd823d21d4d6c979bed74c74024d196e2393696bade787ae42b1bde66c5f | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| wt2_tdf | $GIPPA_DATA_DIR/pxd078459/WT_2.d/analysis.tdf | f77a283e669422688c6e7fe88d8d0d40c3921fa67ee0d02e4645d175abd71892 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| wt3_tdf | $GIPPA_DATA_DIR/pxd078459/WT_3.d/analysis.tdf | 512187712c5a523b60155c0f746ac098826dfbd6c282a78567291eeb23babf9a | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |

## Steps

### S1 resolve (M10-01)

Job: script=scripts/slurm/fetch_pxd078459.sbatch name=gippa-fetch-pxd078459 partition=normal qos=dev cpus=1 gpus=0 mem=2G time=00:30:00 constraint=- array=-

```bash
sbatch --qos=dev --cpus-per-task=1 --mem=2G --time=00:30:00 \
  "$GIPPA_REPO/scripts/slurm/fetch_pxd078459.sbatch" --resolve
```

Reads: -
Writes: $GIPPA_DATA_DIR/pxd078459/record/

The draft manifest was reviewed by hand and committed as
`validation/pilot/inputs/pxd078459_files.tsv` (README › Raw data). Job 1879851 was submitted with
the normal QOS (4 CPUs, 8G, 1 h). It ended with exit 2 after writing the record and the draft.

### S2 download (M10-01)

Job: script=scripts/slurm/fetch_pxd078459.sbatch name=gippa-fetch-pxd078459 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=24:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/fetch_pxd078459.sbatch"
```

Reads: px_record, wt1_zip, wt2_zip, wt3_zip
Writes: $GIPPA_DATA_DIR/pxd078459/

### S3 acquisition mode (M10-01)

Job: script=- name=gippa-acqmode partition=normal qos=dev cpus=2 gpus=0 mem=2G time=00:15:00 constraint=- array=-

```bash
srun -J gippa-acqmode -p normal --qos=dev -c 2 --mem=2G -t 00:15:00 \
  "$PY" -I "$GIPPA_REPO/validation/pilot/acquisition_mode.py" \
  "$GIPPA_DATA_DIR"/pxd078459/WT_1.d "$GIPPA_DATA_DIR"/pxd078459/WT_2.d "$GIPPA_DATA_DIR"/pxd078459/WT_3.d \
  --out "$TMPDIR/acquisition.json"
cmp "$TMPDIR/acquisition.json" "$GIPPA_REPO/validation/pilot/inputs/acquisition.json"
```

Reads: wt1_tdf, wt2_tdf, wt3_tdf
Writes: validation/pilot/inputs/acquisition.json (committed)

Job 1879882 ran it as a dev-QOS batch job. The `srun` above is the equivalent command. It writes
to `$TMPDIR` and compares, so a repeat never rewrites the committed file.

### S4 image install (M10-02)

Job: script=scripts/slurm/diann_image.sbatch name=gippa-diann-image partition=normal qos=normal cpus=1 gpus=0 mem=2G time=00:30:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/diann_image.sbatch"
```

Reads: diann_image
Writes: $GIPPA_ROOT/images/diann_2.3.1.sif

Job 1879883 ended FAILED (exit 2) on the banner check fixed in e7edc457. The copy was verified by
hand afterwards (mode 444, the pin). S-check re-verifies it at F. The image is never re-pinned.

### S5 searches (M10-02)

Job: script=scripts/slurm/diann_pilot.sbatch name=gippa-diann-A1 partition=normal qos=normal cpus=16 gpus=0 mem=64000M time=24:00:00 constraint=- array=-

```bash
sbatch -J gippa-diann-A1 "$GIPPA_REPO/scripts/slurm/diann_pilot.sbatch" A1
# after the A1 job has ended (two jobs of one arm must never run at once):
sbatch -J gippa-diann-A2 "$GIPPA_REPO/scripts/slurm/diann_pilot.sbatch" A2
```

Reads: proteome_a1, contaminants, strain_db, cfg_a1, cfg_a2, acquisition, diann_image, wt1_tdf, wt2_tdf, wt3_tdf
Writes: $GIPPA_ROOT/pilot/A1/, $GIPPA_ROOT/pilot/A2/

**Not re-run at the freeze** (index decision 4). DIA-NN may not repeat byte for byte (M10-06 L3).
A later reader who repeats S5 should take the converted-table hashes below as the target. The
pass-1 jobs asked for `--mem=64G` and got 17 CPUs; the script now asks for 64000M (16 CPUs).

### S6 report check (M10-02)

Job: script=scripts/slurm/diann_report.sbatch name=gippa-diann-report-A1 partition=normal qos=normal cpus=4 gpus=0 mem=16000M time=01:00:00 constraint=- array=-

```bash
sbatch -J gippa-diann-report-A1 "$GIPPA_REPO/scripts/slurm/diann_report.sbatch" A1
sbatch -J gippa-diann-report-A2 "$GIPPA_REPO/scripts/slurm/diann_report.sbatch" A2
```

Reads: a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/A1/check.json, $GIPPA_ROOT/pilot/A2/check.json, $GIPPA_ROOT/pilot/recorded/diann/2.3.1/

### S7 conversion (M10-03)

Job: script=scripts/slurm/diann_convert.sbatch name=gippa-diann-convert-A1 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=00:30:00 constraint=- array=-

```bash
sbatch -J gippa-diann-convert-A1 "$GIPPA_REPO/scripts/slurm/diann_convert.sbatch" A1
sbatch -J gippa-diann-convert-A2 "$GIPPA_REPO/scripts/slurm/diann_convert.sbatch" A2
```

Reads: a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/inputs/

S-check re-converts both reports at F and compares the hashes (example 7). A changed table is
re-converted with these commands at F, and the old Outputs row becomes `superseded` (example 8).

### S-check reproducibility check (M10-14)

Job: script=scripts/slurm/pilot_reproduce_check.sbatch name=gippa-pilot-repro-check partition=normal qos=normal cpus=2 gpus=0 mem=8G time=01:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/pilot_reproduce_check.sbatch"
```

Reads: genome, bakta_json, strain_db, cfg_a1, cfg_a2, acquisition, proteome_a1, contaminants, diann_image, px_record, wt1_tdf, wt2_tdf, wt3_tdf, a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/reproduce/

It runs at F, before S8, and must exit 0 (freeze gate item 4). It runs `static`, then `verify`
into `verify-<F7>.json`, then the re-conversion into `$TMPDIR`.

### S8 freeze A2 (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a2 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
G=$GIPPA_REPO/tests/data/genomes/nissle1917.fna.gz
J=$GIPPA_REPO/tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz
sbatch -J gippa-run-freeze-a2 "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P/a2.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2
```

Reads: genome, bakta_json, a2_table
Writes: $R/out_a2

The main checkout stays at F until the last freeze head job has ended (freeze gate item 3). Each
attempt goes into `$R/ATTEMPTS.tsv`. An accepted attempt records `git_commit = F` and
`git_dirty: false`.

### S9 freeze A1 (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a1 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
# after S8's head job has ended (any state):
sbatch -J gippa-run-freeze-a1 "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P/a1.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a1
```

Reads: genome, bakta_json, a1_table
Writes: $R/out_a1

### S10 freeze TE7 replicate (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a2-rep partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
sbatch -J gippa-run-freeze-a2-rep "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P/a2.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2_rep
```

Reads: genome, bakta_json, a2_table
Writes: $R/out_a2_rep

### S15 clean-up (M10-14)

Job: script=scripts/slurm/pilot_cleanup.sbatch name=gippa-pilot-cleanup partition=normal qos=normal cpus=2 gpus=0 mem=8G time=02:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/pilot_cleanup.sbatch"            # dry run: plan.tsv only
CLEAN=1 sbatch "$GIPPA_REPO/scripts/slurm/pilot_cleanup.sbatch"    # after reading the plan
```

Reads: -
Writes: $GIPPA_ROOT/pilot/cleanup/

It runs once, last. It deletes the Nextflow work directories and the mzML files. Until Khoa
answers K2, the DIA-NN reports and libraries are listed as "kept pending K2" and stay. After the
clean-up, M10-03's `slow` full-report test still finds the reports, and the 1,000-row sample test
keeps running.

## Steps that enter with their PRs

These steps are not yet headings, because their scripts are not on main. Each one gets a `### S`
heading, a `Job:` line and its command when its PR merges:

| Step | Issue | What | Re-run at freeze? |
|---|---|---|---|
| S11 A3 | M10-08 to M10-11 | the GIPPA run on the active A3 branch (below) | the GIPPA run yes; the de novo search no |
| S12 compare | M10-04/10 | `pilot_compare.sbatch` (`gippa-pilot-compare`) | yes |
| S13 evidence | M10-04 | `run_evidence.py arm --expect-commit F` per arm | yes |
| S14 TE and metrics | M10-05/06/07 | `te_real.py`, `entrapment.sbatch`, `metrics.py compute` and `render` into `$GIPPA_ROOT/pilot/metrics/freeze-<F7>/` | yes |
| S-bench | M10-12 | `benchmark_standard.sbatch` (`gippa-bench`) | no; `BENCH_COMMIT` recorded |

## A3

### A3-cascadia (status: inactive)

Active when `validation/pilot/metrics/a3_feasibility.json` has verdict `pass`. The Cascadia
searches and their conversion are not re-run. S11 runs `gippa run` on `$P/a3.tsv.gz` with
`--reference-peptides $P/a2.tsv.gz` into `$R/out_a3`, and S12 compares A1, A2 and A3. The
commands come from M10-09 and M10-10.

### A3-fallback (status: inactive)

Active when the verdict is `fail`, or when an E1-E3 outcome sent A3 to M10-11. Casanovo is not
re-run. S11 runs `gippa run` on the UTI89 genome without `--reference-peptides` into
`$R/out_a3_uti89`, and S12 compares only A1 and A2. The inputs join the Inputs table from
M10-11a's manifests.

## Verify a reproduction

Run the S8 command at F into a new outdir. Then run M10-06's
`te_real.py te7 --ref $R/out_a2 --rep <new> --repo $GIPPA_REPO --out <dir>`. "Reproduced" means
TE7 `pass`: every gate file is byte-identical, so the gate digest matches the Outputs row. A
`known_lapack` diagnosis on another CPU class is the documented non-portability (M4-20 decision
9). It does not count as reproduced.

## Evidence

Pass-1 jobs (`sacct`). The `commit` column is the commit the job ran at, where it was recorded.

| step | label | job | state | exit | elapsed | maxrss | node | commit |
|---|---|---|---|---|---|---|---|---|
| S1 | resolve (first try) | 1879840 | superseded | 2:0 | 00:00:07 | 14M | wheat-02 | - |
| S1 | resolve | 1879851 | FAILED | 2:0 | 00:00:07 | 15M | wheat-08 | - |
| S2 | download | 1879859 | COMPLETED | 0:0 | 00:23:21 | 5.7G | wheat-07 | - |
| S3 | acquisition mode | 1879882 | COMPLETED | 0:0 | 00:00:07 | 12M | wheat-05 | - |
| S4 | image install | 1879883 | FAILED | 2:0 | 00:00:16 | 654M | wheat-05 | - |
| S5 | A1 step 1 | 1879913 | COMPLETED | 0:0 | 00:08:52 | 8.0G | barley-03 | 26adc18c40e6eceb7f30b42f84e27a491d606b82 |
| S5 | A2 step 1 | 1879914 | superseded | 1:0 | 00:08:06 | 8.1G | barley-04 | 26adc18c40e6eceb7f30b42f84e27a491d606b82 |
| S5 | A1 | 1879955 | COMPLETED | 0:0 | 00:14:21 | 6.6G | rye-02 | 81272fc0c68dd45fc4f82309d9eb3e6c770ee0df |
| S5 | A2 | 1879956 | COMPLETED | 0:0 | 00:14:21 | 6.4G | rye-02 | 81272fc0c68dd45fc4f82309d9eb3e6c770ee0df |
| S6 | A1 | 1880002 | COMPLETED | 0:0 | 00:00:11 | 203M | wheat-01 | - |
| S6 | A2 | 1880003 | COMPLETED | 0:0 | 00:00:07 | 203M | wheat-08 | - |
| S7 | A1 (before acquisition.json) | 1880208 | superseded | 0:0 | 00:00:18 | - | wheat-01 | b051c2b75f1daf341dcf13dc491b4cbbee402903 |
| S7 | A2 (before acquisition.json) | 1880209 | superseded | 0:0 | 00:00:20 | - | wheat-05 | b051c2b75f1daf341dcf13dc491b4cbbee402903 |
| S7 | A1 | 1880219 | COMPLETED | 0:0 | 00:00:17 | 266M | wheat-01 | fdf586be0b1310a484ac51c8483f2e3b9b6d663c |
| S7 | A2 | 1880220 | COMPLETED | 0:0 | 00:00:17 | 290M | wheat-05 | fdf586be0b1310a484ac51c8483f2e3b9b6d663c |
| S-check | at F | pending | pending | pending | pending | pending | pending | pending |
| S8 | A2 freeze | pending | pending | pending | pending | pending | pending | pending |
| S9 | A1 freeze | pending | pending | pending | pending | pending | pending | pending |
| S10 | A2 replicate | pending | pending | pending | pending | pending | pending | pending |
| S15 | clean-up | pending | pending | pending | pending | pending | pending | pending |

## Outputs

The DIA-NN reports (kept pending K2) and the converted tables. The freeze `report.html` files,
`compare.html`, `metrics.json` and the per-arm gate digests join in 14b.

| id | path | sha256 | recorded_in | kept |
|---|---|---|---|---|
| a1_report | $GIPPA_ROOT/pilot/A1/A1_report.parquet | f217a40fc0777e99783f3d410743bc80d301290dcdb219a6c69574ee6a4926d7 | validation/pilot/README.md;$GIPPA_ROOT/pilot/A1/outputs.sha256 | yes |
| a2_report | $GIPPA_ROOT/pilot/A2/A2_report.parquet | 733a3884a15664f13385b241e48186e538308993ff2d5a8cc3ec9ff212866d36 | validation/pilot/README.md;$GIPPA_ROOT/pilot/A2/outputs.sha256 | yes |
| a1_table | $GIPPA_ROOT/pilot/inputs/a1.tsv.gz | 81818184f0cd222f091ae6b835ff72cc5743b6ce7e6403e7a2cb1740bbf24f3a | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a1.convert.json | yes |
| a1_abundance | $GIPPA_ROOT/pilot/inputs/a1_abundance.tsv.gz | 20284d4234055b26cf779dc405548087e4efc9bf1007dc3bebc031c5b6c74c54 | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a1.convert.json | yes |
| a2_table | $GIPPA_ROOT/pilot/inputs/a2.tsv.gz | 1eb5e1e9a8855fc8ac62bed4d60f2715f9270f2b96f40206dd35a3a754adf14f | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a2.convert.json | yes |
| a2_abundance | $GIPPA_ROOT/pilot/inputs/a2_abundance.tsv.gz | c61122f3afc9aa31986de8bbfef1c542719627d3d632ec9113e7719582b41c19 | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a2.convert.json | yes |
