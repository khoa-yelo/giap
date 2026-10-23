# M10-02: DIA-NN A1 and A2 as Slurm jobs: image check and read-only copy, provenance, report samples

Milestone M10 · **pass 1** (MVP) · depends on: M10-01 (the three `.d` runs, the committed
`validation/pilot/inputs/acquisition.json` with `dataset_mode = diaPASEF`) · estimate ≈ 330 code
lines (the index said 120; see "Estimate") + ≈ 450 test lines · not `stats` · roles: validation
(`validation/pilot/`), implementer (`scripts/slurm/`), test-author (`tests/`) · Khoa: ADR-0051
(DIA-NN 2.3.1 as the version of record, drafted here), and a re-check of the flags if the
acquisition mode forces a settings change.

## Spec sections implemented

- §15.1 arms A1 ("DIA-NN 2.7.0, library-free, against UniProt K-12 UP000000625 + contaminants")
  and A2 ("same settings, against Nissle's Bakta proteins (sORFs included) + contaminants, from
  `giap makedb`"). The version is 2.3.1 Academia (index decision 2, ADR-0051).
- §15.2 step 5: "Run DIA-NN for A1 and A2 with `scripts/slurm/diann_pilot.sbatch` (`normal`, 16
  CPU, 64 GB), including decoy reporting, so GIAP can use the picked tier".
- §15.2 step 7, the DIA-NN half only: "a 1,000-row sample of each raw report under
  `tests/data/recorded/diann/2.7.0/` … for contract tests (DIA-NN output tables are not binaries,
  so R11 allows the sample)". The directory is `2.3.1/` (index decision 6). The conversion itself
  is M10-03.
- §15.2 DIA-NN template, as amended by the M3-24 README flag table (Khoa accepted it on
  2026-10-20) and the two committed cfgs.
- §3.7 row "DIA-NN searches (A1, A2)"; §3.8 (header, submit and record, poll, measured numbers);
  §3.9 rows "DIA-NN reports, spectral libraries" (measure; keep only the converted tables after the
  pilot); R5 (DIA-NN runs as a separate process in its container); R6 (pinned artefact, checksum,
  never an unofficial copy); R11 (never DIA-NN binaries, ≤ 5 MB per file, `tests/data/` < 50 MB);
  R14 (`giap-*`, ≤ 128 CPUs, ≤ 20 jobs, large files only under `$GIAP_ROOT`).
- §17 rows "PXD078459 acquisition mode differs", "DIA-NN 2.7.0 flags differ from the template".
- Index lead decision 2 (verify the shared `.sif` against the pin, copy it to
  `$GIAP_ROOT/images/diann_2.3.1.sif` with mode 0444, and run only from that copy) and lead
  decision 3 (draft ADR-M10-A, numbered **ADR-0051** here).

## What exists (main d29a2ea5 + `m10/index`)

| Piece | State |
|---|---|
| `scripts/slurm/diann_pilot.sbatch` | M3-24: `giap-diann`, normal, 16 CPU, 64 GB, 24 h. It reads `validation/pilot/diann_<arm>.cfg` and runs the two calls (step 1 is skipped when `<arm>.predicted.speclib` exists). It runs **`${DIANN_SIF:-/farmshare/…/jordan/diann_2.3.1.sif}` directly**, checks no checksum, takes every `*.d` in `RAW_DIR`, and writes `diann_version.txt` but no provenance |
| `validation/pilot/diann_A{1,2}.cfg` | the settings; they differ only in `--fasta` (tested). Unchanged by this issue |
| `validation/pilot/README.md` | version of record, the image's path and sha256 `1946452d…6381`, the flag table |
| `tests/helpers/pilot_kit.py` | the sandbox: fake `apptainer`, a 12-byte stand-in image (`stub sif\n`), mode 0644, three runs `wt_rep{1,2,3}.d`, `DIANN_SIF`, `RAW_DIR` and `OUT_DIR` set |
| `tests/unit/test_m3_24_*.py` | static and stub tests of the kit. **Three of them assert that the container runs `DIANN_SIF` or its default directly**, which lead decision 2 now forbids (see "Tests that must change") |
| `tests/contract/test_diann_contract.py` | `pending_sample`; it looks for `tests/data/recorded/diann/report.{parquet,tsv}`, the top level, not `2.3.1/` (M10-03 adapts it) |
| `src/giap/io/converters/diann.py` | reads parquet or TSV. Required: `Stripped.Sequence`, `Protein.Ids`, `Global.Q.Value` or `Q.Value`, `Run`, `Precursor.Charge`. Optional: `Modified.Sequence`, `PEP`, `Decoy`, `Precursor.Normalised`/`Precursor.Quantity`, `Protein.Group`, `PG.MaxLFQ` |

Checked on disk on 2026-10-23, in the session allocation:

| Fact | Value |
|---|---|
| Shared image | `/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif`: 532,586,496 bytes, mode `-rwxrwxrwx`, owner `khoang99`, group `operator`, mtime 2026-07-09 18:20 −0700. Its directory is `drwxrwxrwx` too. (The index calls the image another user's; `stat` gives `khoang99`. It is world-writable either way.) |
| sha256 now | `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381`: **equals the pin** (2.5 s to hash) |
| `apptainer inspect` labels | build date `Tuesday_6_January_2026_19:18:54_PST`, bootstrap `docker`, from `khoayelo189/diann:2.3.1`, Apptainer 1.4.4 |
| `/diann-2.3.1/` in the image | `diann-linux` (25,747,216 bytes), `libtimsdata.so`, `ThermoFisher.CommonCore.*.dll`, `libtorch_cpu.so`, `models/` (`base-*.pt`, `light-*.pt`, `tiny-*.pt`), `diann-stats.py`. No licence file: `find / -maxdepth 3 -iname '*licen*'` finds nothing |
| `$GIAP_ROOT/images/`, `$GIAP_ROOT/pilot/` | do not exist yet |
| `$GIAP_TOOLS_DIR` | empty: no `diann-2.7.0` |

## Deliverables

| # | File | Role | What |
|---|---|---|---|
| 1 | `validation/pilot/diann_image.sha256` | validation | the pin, one line in `sha256sum` format: `1946452d…6381  diann_2.3.1.sif`. It is the **only** place the scripts read the pin from, and a test keeps it equal to the README and to `DIANN_SIF_SHA256` in the test helper |
| 2 | `scripts/slurm/diann_image.sbatch` | implementer | new script: verifies the shared image and installs the read-only copy |
| 3 | `scripts/slurm/diann_pilot.sbatch` | implementer | changed: it runs only the verified copy, checks the acquisition mode first, takes the runs from `acquisition.json`, uses a step-1 done marker, and writes provenance and `outputs.sha256` |
| 4 | `validation/pilot/diann_report.py` | validation | `check`, `sample` and `provenance` subcommands (below). `provenance` uses the stdlib only; `check` and `sample` use pyarrow |
| 5 | `scripts/slurm/diann_report.sbatch` | implementer | new script: `check`, then `sample`, for one arm |
| 6 | `tests/data/recorded/diann/2.3.1/` | validation | the committed samples and their `PROVENANCE.md` (layout below) |
| 7 | `validation/pilot/README.md` | validation | a "Searches (M10-02)" section: job IDs, node, elapsed time and MaxRSS, report row counts, the location of the copy; the pin line names deliverable 1 |
| 8 | `docs/adr/ADR-0051-diann-version.md` | spec-writer (this commit) | proposed; needs Khoa |

Nothing else is committed: not the image, not any file from inside it, not a speclib, not an
`*_empirical.parquet`, not a full report or a matrix.

### Jobs and resources (R14)

| Job (script) | `-J` | Partition / QOS | CPUs | Mem | Time limit | Node class | When |
|---|---|---|---|---|---|---|---|
| install the copy (`diann_image.sbatch`) | `giap-diann-image` | normal / normal | 1 | 2 GB | 00:30:00 | any | once, before the searches |
| search A1 (`diann_pilot.sbatch A1`) | `giap-diann-A1` | normal / normal | 16 | 64 GB | 24:00:00 | any; recorded | after the install |
| search A2 (`diann_pilot.sbatch A2`) | `giap-diann-A2` | normal / normal | 16 | 64 GB | 24:00:00 | any; recorded | concurrent with A1 |
| check and sample (`diann_report.sbatch A1`, `… A2`) | `giap-diann-report-A1`/`-A2` | normal / normal | 4 | 16 GB | 01:00:00 | any | after each search (`--dependency=afterok:<id>`) |

- The peak is 32 CPUs (both searches), or 40 if the report jobs overlap. Before each `sbatch`,
  the lead sums the CPUs of every running and pending `giap-*` job from every agent; if the sum
  plus 32 exceeds 128, A2 waits for A1. Memory: 64 GB / 16 CPUs = 4 GB per CPU, which equals
  `MaxMemPerCPU`, so Slurm adds no CPUs.
- The install and report jobs use the `normal` QOS, not `dev`. `dev` allows one job at a time and
  M10-01 or M10-08 may hold it.
- The searches do not pin a node class: runtimes are first measured here, and the benchmark node
  class belongs to M10-12. `provenance.json` records the node and its CPU model (`lscpu` model
  name) so the runtimes can be compared later.
- Submission, all from the worktree with `GIAP_REPO` pointing at it, so the cfgs and scripts of
  this branch are the ones that run:

  ```bash
  export GIAP_REPO=/home/classes/bios/270/khoa/giap/local/worktrees/m10-02
  img=$(sbatch --parsable -J giap-diann-image "$GIAP_REPO/scripts/slurm/diann_image.sbatch")
  a1=$(sbatch --parsable -J giap-diann-A1 --dependency=afterok:$img "$GIAP_REPO/scripts/slurm/diann_pilot.sbatch" A1)
  a2=$(sbatch --parsable -J giap-diann-A2 --dependency=afterok:$img "$GIAP_REPO/scripts/slurm/diann_pilot.sbatch" A2)
  sbatch -J giap-diann-report-A1 --dependency=afterok:$a1 "$GIAP_REPO/scripts/slurm/diann_report.sbatch" A1
  sbatch -J giap-diann-report-A2 --dependency=afterok:$a2 "$GIAP_REPO/scripts/slurm/diann_report.sbatch" A2
  ```

  Each job ID, its purpose and its log path go under PROGRESS › In progress at submission. After
  each job ends, `sacct -j <id> --format=JobID,State,Elapsed,MaxRSS,ExitCode,NodeList` goes to
  PROGRESS › Measured numbers (§3.8).

### Image verification and the read-only copy (`diann_image.sbatch`)

Variables: `DIANN_SIF_SHARED` (the source; default
`/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif`) and `DIANN_SIF` (the copy; default
`$GIAP_ROOT/images/diann_2.3.1.sif`). `PIN` is the hex field of deliverable 1.

1. The pin file must hold exactly one line, a 64-hex sha256 and the name `diann_2.3.1.sif`; else
   exit 2.
2. If `DIANN_SIF` already exists:
   - if its sha256 equals `PIN` and it has no write bit (`stat -c %a` gives `444`, `440` or
     `400`), print "already installed: <path>" and exit 0. The source is not read;
   - otherwise exit 2, "the installed copy <path> does not match the pin (or is writable);
     inspect it and delete it by hand". The job never replaces an existing copy itself, so a
     tampered copy is kept as evidence.
3. If `DIANN_SIF_SHARED` is not a readable regular file, exit 2, naming the path and
   "ask Khoa for the authentic 2.3.1 Academia image; never fetch DIA-NN elsewhere (R6)".
4. Free space under `$GIAP_ROOT/images` must be ≥ 2 × the source size; else exit 2.
5. `mkdir -p -m 0755 $GIAP_ROOT/images`. Copy with `cp --preserve=timestamps` to
   `$GIAP_ROOT/images/.diann_2.3.1.sif.tmp.$SLURM_JOB_ID`. Then hash **the copy**, not the
   source: the source is world-writable and could change during the copy. If the copy's sha256
   is not `PIN`, delete the temporary file and exit 2, naming both hashes and "the shared image
   no longer matches the pin". The lead then adds a Blocked row and asks Khoa.
6. `chmod 0444` the temporary file, then `mv -n` it onto `DIANN_SIF`, an atomic rename in the same
   directory. If another job won the race, `mv -n` leaves its file in place; the job re-checks
   rule 2 on the result.
7. Write `$GIAP_ROOT/images/diann_2.3.1.sif.json`:
   - the source: path, size, mode, owner, group, mtime and its sha256 (from step 5);
   - the copy: path, size, mode and sha256;
   - the pin and the pin file's path;
   - `installed_utc`, `job_id`, `hostname`;
   - the `apptainer inspect` labels and the Apptainer version.

   The banner check runs as well:
   `apptainer exec $DIANN_SIF /diann-2.3.1/diann-linux --version`. The first line must contain
   `DIA-NN 2.3.1 Academia`, else exit 2 (the README documents that 2.3.1 prints the banner, warns
   about `--version` and exits 0). The first two lines go into the JSON as `banner`.

Exit codes: 0 installed or already installed; 2 any precondition or mismatch. There is no exit 1:
nothing here is retried automatically.

### The searches (`diann_pilot.sbatch A1|A2`)

The script's argument handling, the cfg parsing (`load_step`, `expand`, `resolve_fastas`), the
binds and the two calls stay as M3-24 wrote them. The changes, in order:

1. **Image.** `SIF="${DIANN_SIF:-$GIAP_ROOT/images/diann_2.3.1.sif}"`. If it is missing, exit 2,
   naming the path, `DIANN_SIF`, and "submit scripts/slurm/diann_image.sbatch first". If its
   sha256 is not `PIN` (deliverable 1), exit 2 with both hashes. If any write bit is set, exit 2
   ("the copy must be read-only, mode 0444"). All three checks run **before any container
   call**, including the `--version` smoke check. The banner check of the install job is
   repeated, and its output goes into `diann_version.txt`, as now.
   The script never reads `DIANN_SIF_SHARED`.
2. **Acquisition mode.** It reads `ACQUISITION_JSON` (default
   `$GIAP_REPO/validation/pilot/inputs/acquisition.json`, from M10-01) with
   `$GIAP_PYTHON -I validation/pilot/diann_report.py runs ACQUISITION_JSON`. That prints one run
   label per line, in the file's order, or exits:
   - with 2 if the file is missing or unparsable;
   - with **3** if `dataset_mode` is not `diaPASEF`. The message is "§17: the acquisition mode
     is <mode>; the DIA-NN settings must change (ADR-0051 amendment, Khoa re-checks the flags)";
   - with **3** too if any run's `notes` names a diaPASEF variant (synchro-, slice- or
     midia-PASEF), until the lead or Khoa confirms the settings (question Q5).

   `GIAP_PYTHON` defaults to `$UV_PROJECT_ENVIRONMENT/bin/python`.
3. **Runs.** `RAW_FILES` becomes `--f $RAW_DIR/<run>.d` for each label in that order, not the
   `*.d` glob. A label with no `analysis.tdf` or no `analysis.tdf_bin` under its `.d` is exit 2,
   naming the run. A `*.d` in `RAW_DIR` that is not in the list is ignored, with a warning that
   names it. An empty list is exit 2 (as now).
4. **Step 1 resume.** Step 1 is skipped only if `$OUT/step1.done` exists and the sha256 it
   records equals the current `<arm>.predicted.speclib`. Otherwise any leftover
   `<arm>.parquet`/`<arm>.predicted.speclib` is deleted and step 1 runs again. After a step 1 that
   exits 0 and leaves a non-empty speclib, the script writes `step1.done`
   (`<sha256>  <arm>.predicted.speclib`). This stops a library truncated by a timeout from being
   reused.
5. **Step 2 is always fresh.** Before step 2 the script deletes `$OUT/<arm>_report*`,
   `$OUT/<arm>_empirical*`, the matrix files a previous attempt left (DIA-NN writes them next to
   `--out`; the names 2.3.1 uses are recorded on the first real run, see "Output layout"),
   `outputs.sha256` and `FAILED`. The logs of earlier attempts stay in `$GIAP_LOG_DIR`.
6. **After step 2.** A non-zero DIA-NN exit is exit 1. So is a missing or zero-byte
   `<arm>_report.parquet`. In both cases the script writes `$OUT/FAILED` (step, exit code, UTC
   time) and no `outputs.sha256`. On success it writes `outputs.sha256`: one line per regular
   file under `$OUT` that DIA-NN wrote, excluding `fasta/` and the provenance files, sorted by
   name.
7. **Provenance.** It writes `$OUT/provenance.json` with
   `$GIAP_PYTHON -I validation/pilot/diann_report.py provenance …` on every exit after the
   argument check, through an `EXIT` trap, so failures are recorded too (schema below). It also
   writes `command_step1.txt` and `command_step2.txt`: the exact argv of each call, one token per
   line, after expansion.

Exit codes: 0 report written; 1 DIA-NN failed or wrote no report; 2 a precondition (usage, image,
pin, mode, cfg, FASTA, runs, acquisition file); 3 the acquisition mode forbids these settings.

### The exact DIA-NN command lines

The command lines are built from the committed cfgs, with `GIAP_DATA_DIR`, `GIAP_DB_DIR` and
`GIAP_ROOT` at their FarmShare values, `SLURM_CPUS_PER_TASK=16`, and the runs `WT_1`–`WT_3`
(M10-01). `C` stands for
`apptainer exec -B /home/classes -B /scratch -B /farmshare /home/classes/bios/270/khoa/giap/local/images/diann_2.3.1.sif /diann-2.3.1/diann-linux`.

**A1, step 1** (`OUT=/home/classes/bios/270/khoa/giap/local/pilot/A1`):

```bash
C --fasta-search --gen-spec-lib --predictor \
  --fasta /home/classes/bios/270/khoa/giap/local/data/proteomes/UP000000625.fasta \
  --fasta /home/classes/bios/270/khoa/giap/local/db/contaminants/contaminants.fasta \
  --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --unimod4 \
  --var-mods 1 --var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4 \
  --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --peptidoforms --verbose 1 \
  --threads 16 --out-lib $OUT/A1.parquet
```

**A1, step 2**:

```bash
C --fasta /home/classes/bios/270/khoa/giap/local/data/proteomes/UP000000625.fasta \
  --fasta /home/classes/bios/270/khoa/giap/local/db/contaminants/contaminants.fasta \
  --gen-spec-lib --qvalue 0.01 --report-decoys --matrices --peptidoforms --verbose 1 \
  --f /home/classes/bios/270/khoa/giap/local/data/pxd078459/WT_1.d \
  --f /home/classes/bios/270/khoa/giap/local/data/pxd078459/WT_2.d \
  --f /home/classes/bios/270/khoa/giap/local/data/pxd078459/WT_3.d \
  --threads 16 --lib $OUT/A1.predicted.speclib \
  --out $OUT/A1_report.parquet --out-lib $OUT/A1_empirical.parquet
```

**A2**: the same two calls with `OUT=…/local/pilot/A2`, `A2` in every output name, and the first
`--fasta` of both steps replaced by
`/home/classes/bios/270/khoa/giap/local/pilot/A2/fasta/strain_db.faa`. That is the unpacked
`${GIAP_REPO}/validation/pilot/db_a2/strain_db.faa.gz`. The second `--fasta` is the same
contaminants file.

The quoting of `K*,R*` above is for the reader. The script passes the token unglobbed under
`set -f`, so the container receives `K*,R*` literally. The tests compare the argv token by token
against the cfg (`test_the_cfg_flags_reach_dia_nn`), so the lines above are a consequence of the
cfgs, not a second source of truth.

### Output layout under `$GIAP_ROOT`

```text
$GIAP_ROOT/images/
  diann_2.3.1.sif                 0444, the only image any DIA-NN job runs
  diann_2.3.1.sif.json            install provenance (source stat and hash, copy hash, labels, banner)
$GIAP_ROOT/pilot/<arm>/           <arm> = A1 | A2  (OUT_DIR overrides $GIAP_ROOT/pilot)
  fasta/strain_db.faa             A2 only: unpacked from the committed .gz
  diann_version.txt               banner of the image (smoke check)
  <arm>.parquet                   step 1 --out-lib target (if DIA-NN writes it)
  <arm>.predicted.speclib         step 1 library (never committed)
  step1.done                      sha256 of the speclib that step 1 finished
  <arm>_report.parquet            the report (never committed)
  <arm>_empirical.parquet         step 2 empirical library (never committed)
  <the other files DIA-NN writes next to --out: matrices, stats, log>   recorded as observed
  command_step1.txt, command_step2.txt
  outputs.sha256                  written only after a successful step 2
  provenance.json                 written on every exit
  FAILED                          present only after a failed attempt
  check.json                      diann_report.py check
  sample/report.sample.parquet, sample/schema.tsv, sample/summary.json   diann_report.py sample
```

The spec's §15.2 step 7 says "the large reports stay in `local/data/`". `$GIAP_ROOT/pilot` is
under `local/` as well and is the M3-24 default, which the stub tests use. The layout keeps it.
The sizes of the speclibs, reports and matrices are measured by `du` after the run and go into
PROGRESS › Measured numbers (§3.9: "measure"). M10-14 deletes the libraries and raw reports.

### `validation/pilot/diann_report.py`

CLI `python -I validation/pilot/diann_report.py <command> …`; every command is a pure function
plus a thin CLI. Nothing writes into the repository.

| Command | Inputs → outputs | Exit |
|---|---|---|
| `runs ACQUISITION_JSON` | the run labels, one per line, in file order | 0; 2 missing or unparsable; 3 not diaPASEF or a variant (rule 2 above) |
| `provenance --arm A --out OUT [--status …]` | `OUT/provenance.json` (schema below); stdlib only | 0; 2 usage |
| `check REPORT --acquisition ACQUISITION_JSON --sums OUT/outputs.sha256 --out OUT/check.json` | `check.json`: `ok`, `failures[]`, `rows`, `columns`, `runs` (distinct `Run` values with their row counts), `decoy_rows`, `max_q` (of `Q.Value` and `Global.Q.Value`), the pyarrow version | 0 ok; 1 partial or invalid; 2 usage |
| `sample REPORT --out DIR [--n 1000] [--decoys 100]` | `DIR/report.sample.parquet`, `DIR/schema.tsv`, `DIR/summary.json` | 0; 1 the report fails `check`'s structural rules; 2 usage |

**`check` rules** (each failure is one string in `failures`):

1. The file exists and is non-empty, and `pyarrow.parquet.read_metadata` succeeds. A file cut
   short mid-write has no footer, so it fails here. `num_rows > 0`.
2. Every column the converter requires is present, plus the columns this pilot depends on:
   `Stripped.Sequence`, `Modified.Sequence`, `Protein.Ids`, `Protein.Group`, `Run`,
   `Precursor.Charge`, `Global.Q.Value` or `Q.Value`, `Decoy` and `PG.MaxLFQ`.
3. Each run label of the acquisition file matches some `Run` value, and no `Run` value is left
   unmatched. A label matches when the value equals it, or its basename without a `.d` suffix
   equals it. The exact `Run` spelling 2.3.1 writes is recorded in `PROVENANCE.md`.
4. At least one row has `Decoy = 1`. This shows `--report-decoys` reached the output.
5. The report's sha256 equals its line in `outputs.sha256`: the file is the one the search job
   finished.

### The committed samples (`tests/data/recorded/diann/2.3.1/`)

```text
tests/data/recorded/diann/2.3.1/
  PROVENANCE.md
  A1/report.sample.parquet   ≤ 1,000 rows, every column of the report, the report's Arrow schema
  A1/schema.tsv              column, arrow_type, nullable: every column of the FULL report, in order
  A1/summary.json            aggregates of the full report (below); no row-level data
  A2/…                       the same three files
```

**Row selection** (deterministic; no random number generator, so `giap.rng` is not involved):

1. Read only `Run`, `Decoy` and `Precursor.Id` to select (if `Precursor.Id` is absent, use
   `Modified.Sequence` + `Precursor.Charge`). Give each row a file-order index `i` and a key
   `k = sha256(utf8(Run + "\t" + Precursor.Id)).hexdigest()`.
2. Strata are (`Run`, decoy). Runs are sorted by their `Run` string. There are `R` runs. The
   targets get `n − d` rows and the decoys `d` rows (defaults `n = 1000`, `d = 100`), each split
   evenly over the runs. The remainder goes one row each to the first runs in sorted order. With
   R = 3 that gives 300/300/300 targets and 34/33/33 decoys.
3. Within each stratum, take the rows with the smallest `(k, i)`, up to the stratum's quota.
4. If a stratum has fewer rows than its quota, fill the shortfall from the rows not yet taken, by
   smallest `(k, i)` over all strata, until `n` rows are taken or the report is exhausted.
5. Write the taken rows **in file order** (ascending `i`), with every column, using
   `pyarrow.parquet.write_table` on the source schema (`compression="zstd"`). The same report
   gives the same rows; the bytes depend on the pyarrow version, which `summary.json` records.

This replaces index decision 6 ("the first 1,000 rows in file order"), question Q3. Those rows
would most likely all come from the first run and might hold no decoy, and both the converter's
`run` handling and its `is_decoy` path need both. Hashing `(Run, Precursor.Id)` spreads the
sample over the q-value and retention-time range without a seed.

**`summary.json`** (full report): `rows`, `columns`, `rows_per_run`, `decoy_rows`,
`distinct_precursors`, `distinct_stripped_q01` (distinct `Stripped.Sequence` of target rows with
`Global.Q.Value ≤ 0.01`), `max_q`, `report_sha256`, `report_bytes`, `sample_rows`,
`sample_rule` (`"strata(run,decoy) sha256(Run\tPrecursor.Id) n=1000 d=100"`),
`duplicate_keys`, `pyarrow`.

**`PROVENANCE.md`** (written by hand from the JSON files, in the house style of
`tests/data/recorded/genomad/PROVENANCE.md`):

- the tool, `DIA-NN 2.3.1 Academia` (banner), the image copy's sha256 and the pin, and ADR-0051;
- the job IDs, node and elapsed time of both searches and both report jobs;
- the cfg sha256s and the git commit of `GIAP_REPO`;
- the input sha256s: both FASTAs per arm, `acquisition.json`, and the `analysis.tdf` of each run
  with the manifest sha256 of each archive;
- the exact command lines (a pointer to `command_step*.txt` and their content);
- the report's sha256, its row counts and the observed `Run` spelling;
- the selection rule;
- what is **not** here and why (the licence rules below).

**Licence and size rules** (R5, R11, index licence table). These are committed only:
`PROVENANCE.md`, `*/report.sample.parquet`, `*/schema.tsv` and `*/summary.json`. **Never**
committed:
- the image, or any file copied out of it (`diann-linux`, `*.so`, `*.dll`, `models/*.pt`,
  `diann-stats.py`);
- a spectral library (`*.speclib`, `*.predicted.speclib`, `*_empirical.parquet`, the step-1
  `<arm>.parquet`);
- a full report (`*_report.parquet`, any file over 1,000 data rows), a matrix (`*matrix*.tsv`) or
  a DIA-NN log containing paths of the image.

Each committed file is ≤ 5 MB (R11); the samples are expected to be well under 1 MB.
`tests/data/` stays < 50 MB (26 MB today). The image holds no licence text, so the permission to
commit output rows rests on §15.2 step 7 ("DIA-NN output tables are not binaries, so R11 allows
the sample"). ADR-0051 asks Khoa to confirm that this matches the Academia licence he accepted.

### `provenance.json` (per arm; also written on failure)

`arm`, `status` (`ok` | `failed` | `precondition` | `acquisition_mode`), `exit_code`, `step1`
(`ran` | `skipped` | `not_reached`) with its exit code, `step2` likewise, and `job`:

- `id`, `name`, `partition`, `cpus` (`SLURM_CPUS_PER_TASK`), `mem` (`SLURM_MEM_PER_NODE`);
- `node` (hostname) and `cpu_model`;
- `start_utc`, `end_utc`, `elapsed_s`.

The other fields:

- `image`: path, sha256, pin, mode, banner;
- `apptainer`: version;
- `cfg`: path, sha256;
- `repo`: `GIAP_REPO`, `git rev-parse HEAD`, dirty flag;
- `inputs`: each `--fasta` (path, sha256; for A2 both the `.gz` and the unpacked file),
  `acquisition_json` (path, sha256, `dataset_mode`), and each run (label, `.d` path, sha256 of
  `analysis.tdf`, size of `analysis.tdf_bin`);
- `commands`: the two argv arrays;
- `outputs`: the `outputs.sha256` path, or `null`.

`analysis.tdf_bin` is not hashed in the job, because it runs to gigabytes and the archive's sha256
is in the M10-01 manifest. `sacct` values (MaxRSS) do not exist until the job ends; the agent adds
them to PROGRESS and `PROVENANCE.md`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `diann_image.sbatch`; shared image hashes to the pin; no copy yet | `$GIAP_ROOT/images/diann_2.3.1.sif` exists, mode `444`, sha256 = pin; `diann_2.3.1.sif.json` records source mode `777`, owner, both hashes, labels, banner `DIA-NN 2.3.1 Academia…`; exit 0 |
| 2 | Example 1 run again | "already installed", exit 0; the shared image is not read (no new `cp` and no read of it; the stub logs it) |
| 3 | Shared image replaced by a different file (sha256 `ab12…`) | exit 2, message names `ab12…` and the pin and "never fetch DIA-NN elsewhere"; no new file and no temporary file left under `$GIAP_ROOT/images/`; no `apptainer exec` |
| 4 | `diann_pilot.sbatch A1` with the copy installed, acquisition `diaPASEF`, runs `WT_1`–`WT_3` | exactly two container calls on the copy, with the argv of "The exact DIA-NN command lines"; `step1.done`, `outputs.sha256`, `provenance.json` (`status ok`, `step1 ran`); exit 0 |
| 5 | Example 4, then re-submitted | step 1 skipped (`step1.done` matches), one container call (step 2), old `A1_report.parquet` deleted before it; `provenance.json` `step1 skipped`; exit 0 |
| 6 | `DIANN_SIF` points at a file whose bytes differ from the pin (the copy was altered) | exit 2 before any container call (not even `--version`); message has both hashes; `provenance.json` `status precondition` |
| 7 | The copy has mode `644` but the right hash | exit 2, "the copy must be read-only, mode 0444"; no container call |
| 8 | `acquisition.json` `dataset_mode = ddaPASEF` | exit 3, message cites §17 and ADR-0051; no container call; `provenance.json` `status acquisition_mode` |
| 9 | Step 2 killed at the 24 h limit (Slurm `TIMEOUT`), a 40 %-written `A1_report.parquet` left behind; re-submitted | the re-run skips step 1 (valid marker), deletes the partial report, runs step 2 again. If the first job was killed inside step 1, the marker is absent: the speclib is deleted and step 1 re-runs |
| 10 | `diann_report.sbatch A1` on a report whose parquet footer is missing (truncated) | `check.json` `ok = false`, `failures` names "parquet footer unreadable"; exit 1; no `sample/` |
| 11 | `sample` on a report with 3 runs, 600,000 target rows and 12,000 decoy rows | 1,000 rows: 300 targets and 34/33/33 decoys per run (runs sorted), in file order; `summary.json` `sample_rows 1000`, `duplicate_keys 0` |
| 12 | `sample` where run `WT_3` has only 120 target rows | `WT_3` contributes 120 targets; the 180 missing rows come from the remaining rows of any stratum by smallest key; total 1,000 |

## Edge cases

- **sha mismatch.** For the shared source, see example 3: the lead adds a Blocked row
  "DIA-NN image no longer matches the pin" and asks Khoa; work moves to M10-08. For the copy,
  see example 6: the agent inspects and deletes the copy, then re-runs the install, which
  succeeds only if the source still matches. A copy is never re-pinned to fit: re-pinning
  needs a new image from Khoa and an ADR-0051 amendment.
- **Image missing.** A missing copy is exit 2 with "submit diann_image.sbatch first". A missing
  shared source is exit 2 and a Blocked row. A missing pin file, or a malformed one, is exit 2
  in both scripts.
- **Not diaPASEF.** `ddaPASEF`, `mixed`, `ms1_only`, `unknown` or `inconsistent` is exit 3 before
  any container call, as is a diaPASEF variant named in `notes` (Q5). Example 8 is the `ddaPASEF`
  case. The settings change in an ADR-0051 amendment with Khoa's re-check; this issue does not
  guess other flags.
- **Timeout.** Example 9. One re-submission with the same limit. If it times out again, the second
  attempt asks for `--time=2-00:00:00`, the `normal` partition maximum. A third timeout is R9: a
  question to the lead, with the `long` QOS (7 days, ≤ 4 jobs) as the option.
- **Partial report.** DIA-NN exits non-zero, the report is missing or empty, or its footer is
  unreadable. The first two are exit 1 of the search job; the third is exit 1 of the report job
  (example 10). The rules of `check` cover runs absent from `Run` (a run DIA-NN skipped) and a
  report with no decoy rows. `FAILED` stays until the next attempt starts. M10-03 never converts
  a report whose `check.json` is not `ok`.
- **Runs.** The acquisition file names `WT_1`–`WT_3`, but `RAW_DIR` has an extra `old.d`: it is
  ignored with a warning. A listed run that lacks `analysis.tdf_bin` is exit 2.
- **Concurrency.** A1 and A2 share only the read-only image. Two install jobs racing: each copies
  to its own temporary name, and `mv -n` keeps the first; both then verify the result.
- **Empty input.** A report with 0 rows fails `check` rule 1. One with fewer than 1,000 rows is
  sampled whole (`sample_rows < 1000`, noted in `PROVENANCE.md`).
- **Ties.** Equal keys need duplicate (`Run`, `Precursor.Id`) pairs, which a well-formed report
  does not have. If they occur, `i` breaks the tie, `duplicate_keys` counts them and
  `PROVENANCE.md` reports them. Runs that sort equal cannot occur: `Run` values are distinct
  strings.
- **Stratum quota when R does not divide n − d.** The remainder rule (step 2) makes the split
  unique: it is a function of the sorted runs only.
- **Circular or linear genomes, strands:** not applicable. This issue handles search outputs, not
  coordinates; the A2 database's coordinates are `giap makedb`'s (M3-23).
- **Disk.** A copy needs 532 MB plus margin. Speclibs and reports are unmeasured; the search job
  does not pre-check their space (unknown sizes), but `$GIAP_ROOT` has 197 TB free in the class
  pool, and the 400 GB `local/` budget check in the probe warns.

## Acceptance criteria

Stub and static tests (CI, every commit; no DIA-NN, no container, no real data, no network):

1. `diann_image.sbatch` in the sandbox (fake `apptainer`, a stand-in source image, the sandbox pin
   file holding the stand-in's sha256), shown by examples 1–3:
   - it installs a mode-444 copy whose bytes equal the source;
   - it writes the JSON with both hashes;
   - it is idempotent;
   - a source mismatch is exit 2 with no copy left and no `apptainer exec`;
   - a missing source is exit 2;
   - an existing mismatched copy is exit 2 and is not replaced;
   - a banner without `DIA-NN 2.3.1 Academia` is exit 2.
2. `diann_pilot.sbatch`:
   - **every** container call (smoke check included) uses `DIANN_SIF`, which defaults to
     `$GIAP_ROOT/images/diann_2.3.1.sif`, never `DIANN_SIF_SHARED`;
   - a hash mismatch, a writable copy or a missing copy is exit 2 with zero container calls;
   - `ddaPASEF`/`mixed`/`unknown`/`inconsistent` and a variant in `notes` are exit 3 with zero
     container calls;
   - a missing acquisition file is exit 2;
   - the `--f` values are exactly the acquisition runs in their order;
   - a listed run without `analysis.tdf_bin` is exit 2.
3. Step-1 resume: with a valid `step1.done` step 1 is skipped. With a speclib but no marker, or
   a marker whose hash disagrees, step 1 re-runs and the old library is gone. Step 2 always
   starts with no `<arm>_report*` left from an earlier attempt.
4. Failure: `FAKE_DIANN_RC=1` in step 2 gives exit 1, a `FAILED` file, no `outputs.sha256`, and
   `provenance.json` with `status failed`, `step2` with exit code 1. A stub that exits 0 but
   writes no report gives the same.
5. `provenance.json` validates against the field list above, in a test-local schema. Its
   `commands` equal the argv the fake container logged, and its input hashes equal the sandbox
   files' sha256s.
6. `diann_report.py check` on parquet fixtures built at test time with pyarrow gives `ok` for a
   well-formed one. It names the rule for each of these: a truncated file (last 100 bytes cut), a
   missing `Decoy` column, no decoy rows, a run missing from `Run`, an unexpected `Run`, and a
   sha256 that disagrees with `outputs.sha256`.
7. `diann_report.py sample` gives exactly the selections of examples 11–12. Its tests:
   - small analogues (n = 10, d = 2, 3 runs) whose expected row indices are written out by hand
     from the rule;
   - the output is in file order, with every column and the source schema;
   - a run twice gives identical row sets;
   - a report under n rows is returned whole;
   - duplicate keys are counted.
8. `runs` gives the labels and the exit codes of rule 2, using JSON fixtures in M10-01's
   `acquisition.json` shape.
9. Pin consistency: `validation/pilot/diann_image.sha256`, the README and `DIANN_SIF_SHA256` in
   `tests/helpers/pilot_kit.py` hold the same 64-hex value.
10. Licence and size guard on the committed tree. It applies to whatever exists and is vacuous
    before the samples land:
    - `git ls-files tests/data/recorded/diann/` lists only the four allowed patterns;
    - no tracked file anywhere has a suffix `.sif`, `.speclib` or `.pt`, or starts with the ELF
      magic, under `tests/` or `validation/`;
    - every `report.sample.parquet` has ≤ 1,000 rows (by `read_metadata`), and its columns
      equal the rows of its `schema.tsv`;
    - every file is ≤ 5 MB;
    - `tests/data/` is < 50 MB.
11. House rules for the two new scripts (as `test_sbatch_house_rules` checks): a `giap-*` name,
    `--partition`, `--time`, `--cpus-per-task` and `--mem` as in the jobs table, the
    environment sourced, `set -euo pipefail`, `--output` under `local/logs/`. The kit's
    concurrent total stays ≤ 128 CPUs. No script writes under `$GIAP_REPO` (the sandbox repo
    is clean after every run).
12. The M3-24 tests pass, with only the changes listed under "Tests that must change".

Evidence on real data (PR description, PROGRESS, `PROVENANCE.md`; not CI):

13. The job IDs of the install, both searches and both report jobs, with `sacct` State, Elapsed,
    MaxRSS and NodeList, are in PROGRESS › Measured numbers. Each job's exit code is 0.
14. `$GIAP_ROOT/images/diann_2.3.1.sif` is mode 444 with sha256 = pin, and its JSON is in place.
15. `check.json` is `ok` for A1 and A2. The committed samples and summaries exist for both arms;
    `PROVENANCE.md` states every item listed above.
16. The README "Searches (M10-02)" section exists. ADR-0051 is in the PR, status `proposed`, and
    has a row under PROGRESS › Questions for Khoa.

Link to §15.6: these reports are the input of M10-03's conversion, which feeds "GIAP runs end to
end on A1, A2" (M10-04), TE1 (M10-05) and "A1 and A2 … are compared" (M10-04).

## Test plan

Owner: the test-author, writing tests first (R2). No test runs DIA-NN, starts a container, reads
the real image (except the existing `test_the_pinned_sha256_is_still_the_image_on_disk`, which
only hashes it and skips when it is absent) or touches `$GIAP_DATA_DIR`. Every script runs in
`helpers.pilot_kit.pilot_sandbox`.

New files:

- `tests/unit/test_m10_02_diann_image.py`: criteria 1, 11.
- `tests/unit/test_m10_02_diann_searches.py`: criteria 2–5, using the extended fake `apptainer`.
  It logs calls already; add `FAKE_DIANN_NO_REPORT=1`, which skips writing `--out`, and a banner
  override.
- `tests/unit/test_m10_02_diann_report.py`: criteria 6–8. Fixtures are built in `tmp_path` with
  pyarrow, a few dozen rows, column names from the converter's map.
- `tests/unit/test_m10_02_recorded_diann.py`: criteria 9–10.

Helper changes (`tests/helpers/pilot_kit.py`):

- the sandbox writes `validation/pilot/diann_image.sha256` with the stand-in's sha256 and makes
  the stand-in mode 0444;
- it writes `validation/pilot/inputs/acquisition.json` (`diaPASEF`, runs `wt_rep1`–`wt_rep3`,
  with `analysis.tdf_bin` stand-ins);
- it sets `GIAP_PYTHON` to `sys.executable`;
- it adds a stand-in `DIANN_SIF_SHARED`.

Each change is committed with a re-commit of the sandbox repository, as the helper does now.

### Tests that must change (lead decision 2 contradicts them; lead approval, question Q1)

| Test | Today | After |
|---|---|---|
| `test_m3_24_pilot_kit.py::test_diann_pilot_uses_the_image_from_decision_1` | `${DIANN_SIF:-…}` default must be the jordan path | the default in `diann_pilot.sbatch` is `$GIAP_ROOT/images/diann_2.3.1.sif`; the jordan path is the `DIANN_SIF_SHARED` default of `diann_image.sbatch`; binds and binary unchanged |
| `test_m3_24_diann_pilot_run.py::test_diann_sif_override_is_honoured` | an override stand-in (mode 644) is run | the override must be mode 444 and match the sandbox pin; then it is run |
| `test_m3_24_diann_pilot_run.py::test_the_default_image_is_the_path_of_decision_1` | with `DIANN_SIF` unset, runs the real jordan image path if present | with `DIANN_SIF` unset and no copy under the sandbox `$GIAP_ROOT/images`, exit 2 naming that path and `diann_image.sbatch`. It must never copy the real 532 MB image inside a unit test |

All other M3-24 tests pass unchanged once the helper changes above are in. The fake DIA-NN's
text "reports" stay acceptable, because the search job checks only existence and size; parquet
validity is `diann_report.py`'s job.

Targeted local runs only. `pytest tests/unit/test_m3_24_* tests/unit/test_m10_02_*` with the venv
Python (`--no-sync`; never `uv run` in a worktree) takes well under a minute. The full suite runs
on Slurm (`test_branch.sbatch`), not in the session.

## Estimate

| Piece | Code lines |
|---|---|
| `diann_image.sbatch` | 70 |
| `diann_pilot.sbatch` changes (image checks, acquisition gate, runs, marker, cleanup, trap) | 80 |
| `diann_report.sbatch` | 25 |
| `diann_report.py` (`runs`, `provenance`, `check`, `sample`, CLI) | 150 |
| pin file, README section | 5 (+ prose) |
| **Total** | **≈ 330** (index: 120; the provenance, report check and sampler were not counted there) |

Tests ≈ 450 lines. Under the ≈ 400-line PR guideline (ADR-0007 counts code only), this is one PR.
If the reviewer prefers, it splits into (a) image + searches and (b) `diann_report.py` + samples.

Wall time: the install takes about a minute (2.5 s to hash 532 MB plus the copy). The searches
are unmeasured: the index estimates 2–8 h each against a 24 h limit. The report jobs should take
minutes.

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| Q1 | Three M3-24 tests assert that the container runs `DIANN_SIF` or the jordan path directly, which lead decision 2 forbids. May the test-author amend them (table above) and extend the helper in M10-02's test step? | yes, in the test-author's commit, each change citing lead decision 2 |
| Q2 | Where does the pin live for the scripts? A committed `validation/pilot/diann_image.sha256` (the sandbox can rewrite it; no environment backdoor) or a hard-coded constant in each script (tests then need an override variable)? | the committed file, kept equal to the README and the test constant by a test |
| Q3 | The sample selection replaces index decision 6 ("first 1,000 rows in file order") with the stratified hash rule, because the first rows are probably one run and decoy-free. Accept? | stratified rule (n = 1000, d = 100); it can be re-run in minutes from the kept reports until M10-14 deletes them |
| Q4 | Sample format: parquet (the format 2.3.1 writes and the converter reads on the real data) instead of the index's "gzipped" TSV? | parquet only; the converter's TSV path keeps its hand-written fixtures |
| Q5 | M10-01 Q3 counts synchro-, slice- and midia-PASEF as `diaPASEF` with `notes`. Should the searches run on such a variant without a settings check? | no: exit 3 until the lead or Khoa confirms the cfgs suit the window scheme |
| Q6 | Timeouts: re-submit once at 24 h, then once at 2 days, then R9? | yes |
| Q7 | Lead decision 2 puts the copy in `$GIAP_ROOT/images/`. §3.3 says `$GIAP_TOOLS_DIR` holds "licence-gated binaries (DIA-NN), installed by Khoa", and §14.8 item 6 reserves installing DIA-NN to Khoa. Is an agent copy of Khoa's image under `images/` acceptable? (ADR-0051 asks Khoa.) | `images/` as decided; ADR-0051 records that the agent copies, it does not install, an image Khoa supplied |
| Q8 | `needs_diann` (tests/conftest.py) and the probe look for the jordan path or `$GIAP_TOOLS_DIR/diann-2.7.0*`, not the copy. Teach them the copy here or later? | later: a backlog chore; this issue's tests do not use the marker |
| Q9 | ADR number: 0050 is taken on `origin/m7/calibration-gate-novel` and 0022 is reserved for M6-12 part C, so this draft is ADR-0051. Renumber at merge if another branch lands 0051 first? | yes; the file name and the references here move together |
| Q10 | The search job runs from the worktree (`GIAP_REPO=…/worktrees/m10-02`) before the PR merges, so the evidence lands in the same PR. Is that acceptable, given `provenance.json` records the commit? | yes; the commit and dirty flag are in the provenance, and the PR must be clean (dirty = false) at run time |
| Q11 | The image's labels say it was bootstrapped from Docker Hub `khoayelo189/diann:2.3.1`, and it contains no licence file. The search proceeds on Khoa's 2026-10-20 acceptance, and ADR-0051 asks him to confirm the image's origin and the licence permission for the committed sample rows. Agreed? | proceed (Khoa supplied and accepted the image); the question goes in the next Khoa batch |

## Lead decisions (2026-10-23)
1. **Q1, the three M3-24 tests that pin the direct image path:**
   - `test_diann_pilot_uses_the_image_from_decision_1`;
   - `test_diann_sif_override_is_honoured`;
   - `test_the_default_image_is_the_path_of_decision_1`.

   The **lead** amends them and the `tests/helpers/pilot_kit.py` sandbox, citing index decision 2, in the
   implementation step (role guard: the test-author writes only new tests). The test-author's new tests
   may assume the sandbox additions listed in the spec. Spell them out in a helper under a new name if
   needed.
2. **Q2–Q11:** the proposed defaults are taken, all reversible:
   - the pin lives in a committed file;
   - the sample is stratified;
   - the sample is parquet only;
   - diaPASEF variants stop with exit 3;
   - on timeout, resubmit at 24 h, then at 2 days, then escalate;
   - the copy goes in `$GIAP_ROOT/images/`;
   - the `needs_diann` probe change goes to the backlog;
   - the ADR is renumbered at merge if needed;
   - the real searches run from a clean worktree;
   - work proceeds on the 2026-10-20 acceptance.

   ADR-0051's questions go to Khoa in the next batch.
3. **Ownership:** the shared image is owned by khoang99 (Khoa), not another user, and its directory is world-writable. The verified read-only copy stands. Khoa is asked whether to chmod the source.
4. **Recorded samples (Khoa, 2026-10-23): FarmShare only.**
   - `diann_report.py sample` writes to `$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/` and is never committed.
   - Under `tests/data/recorded/diann/2.3.1/` the repo holds only `PROVENANCE.md` and `<arm>/schema.tsv`, which carries column metadata only, with no rows or summaries.
   - Tests that need real rows read `$GIAP_ROOT` and skip on CI. This supersedes the committed-sample parts of the spec above.
5. **After the tests (d109b51e):** the interfaces the test-author fixed are accepted as the contract:
   - `select_rows`;
   - the provenance and install-record shapes;
   - the check failure strings;
   - `check.json`;
   - the exit codes.

   The free-space check and the `mv -n` race are untested and backlogged.
