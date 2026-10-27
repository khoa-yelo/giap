# Pilot validation kit: E. coli Nissle 1917 (GIAP_SPEC.md §15)

The kit for the light pilot of §15: three arms (A1 generic reference, A2 matched strain, A3 de
novo) over three wild-type runs of **PXD078459** (iProX; Nissle 1917; timsTOF HT, diaPASEF
expected). **M3 prepares the kit; M10 runs it.** Nothing here has been run on real data: the only
thing M3 executed is a `--version` smoke check of the DIA-NN container (M3-24 decision 6).

| File | What it is |
|---|---|
| `diann_A1.cfg` | DIA-NN settings for arm A1 (UniProt K-12 UP000000625 + contaminants) |
| `diann_A2.cfg` | the same settings for arm A2; only `--fasta` differs |
| `db_a2/strain_db.faa.gz` | the A2 database built by `giap makedb` from the recorded Nissle Bakta JSON |
| `db_a2/db_id_map.tsv` | one row per database record (1-based closed coordinates) |
| `db_a2/makedb_manifest.json` | inputs with their SHA-256, options and counts |
| `../../scripts/slurm/fetch_pxd078459.sbatch` | step 4: download the runs to `$GIAP_DATA_DIR/pxd078459/` |
| `../../scripts/slurm/diann_pilot.sbatch` | step 5: DIA-NN for A1 and A2 |
| `../../scripts/slurm/denovo_pilot.sbatch` | step 6: Cascadia (A3), with the Casanovo fallback |

## DIA-NN: the container, not an unpacked binary

Khoa supplied DIA-NN as an Apptainer image instead of the unpacked `$GIAP_TOOLS_DIR/diann-2.7.0`
build of §15.2 step 3 (reference invocation:
`/farmshare/home/classes/bios/270/jordan/run_2step.sh`, which this kit follows).

| | |
|---|---|
| Version string | `DIA-NN 2.3.1 Academia`, compiled 2025-12-05 |
| Image | `/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif` (532 MB, world-readable) |
| sha256 of the image | `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381` (2026-10-20); the scripts read it from `validation/pilot/diann_image.sha256` |
| Binary inside | `/diann-2.3.1/diann-linux` |
| Module | `apptainer/1.4.4`; binds `/home/classes`, `/scratch`, `/farmshare` |
| Override | `DIANN_SIF=/path/to/image.sif sbatch scripts/slurm/diann_pilot.sbatch A1`; the override must be mode 0444 and hash to the pin, or the job exits 2 |
| Smoke check | the job logs `… --version` into `$OUT_DIR/<arm>/diann_version.txt`; 2.3.1 prints the banner, warns `unrecognised option [--version]` and exits 0 |

The image is **never copied into the repository** (532 MB, R11; decision 3): the repository keeps
only the two cfg files, this README, `db_a2/` and the sbatch scripts.
`scripts/slurm/diann_image.sbatch` verifies the image above against the pin and installs a
read-only copy at `$GIAP_ROOT/images/diann_2.3.1.sif` (M10-02); `scripts/slurm/diann_pilot.sbatch`
reads `DIANN_SIF`, defaults to that copy and exits 2 when it is missing, writable or off the pin.
The pin above is checked by `tests/unit/test_m3_24_pilot_kit.py`; the scripts read it from
`validation/pilot/diann_image.sha256`, and the test fakes use `DIANN_SIF_SHA256`
(`tests/helpers/pilot_kit.py`). If the image is replaced, re-pin it in all of these and here.

## Version actually used

**Khoa reviewed and accepted the flag table and this section on 2026-10-20** ("looks fine"), including
the two deviations from the §15.2 template: `--reanalyse` / `--rt-profiling` left out, and
`--met-excision` absent. 2.3.1 Academia is the version of record unless he says otherwise.


The spec names DIA-NN **2.7.0** (§15.1, §15.2 step 3). The version actually available on FarmShare,
and the one this kit is written and recorded against, is **2.3.1** Academia. The spec text is not
rewritten (decision 2). If a 2.7.0 build arrives later, only `diann_A1.cfg`, `diann_A2.cfg` and this
section change; the recorded reports then carry the 2.7.0 version string, and
`tests/data/recorded/diann/` gets the sample of whichever version produced them.

## The two-step DIA-NN workflow

In DIA-NN 2.x, combining `--fasta-search --predictor` with raw files in one call disables the FASTA
digest. The pilot therefore runs two calls in order (decision 4), both in `diann_<arm>.cfg`:

1. `[step1]` predicts the in-silico spectral library from the FASTA files alone
   (`--fasta-search --gen-spec-lib --predictor`, `--out-lib $OUT/<arm>.parquet`, from which DIA-NN
   writes `<arm>.predicted.speclib`).
2. `[step2]` searches the runs against that library (`--lib <arm>.predicted.speclib`, `--f <run>.d`,
   `--out <arm>_report.parquet`).

Step 1 is skipped when its library already exists, so a re-queued job resumes. `--f`, `--lib`,
`--fasta`, `--out`, `--out-lib` and `--threads` (= `$SLURM_CPUS_PER_TASK`) are supplied by the
sbatch script; everything else comes from the cfg file, one flag per line, with `$GIAP_*`
variables expanded at submission time. A2's database is committed gzipped and DIA-NN reads plain
FASTA, so the script unpacks `strain_db.faa.gz` into `$OUT_DIR/<arm>/fasta/`, outside the
repository.

### Digest flags = GIAP's digest (decision 5)

`--cut K*,R*`, `--missed-cleavages 1`, `--min-pep-len 7`, `--max-pep-len 30` and `--unimod4` are
exactly the defaults of `giap.mapping.digest` (M3-22); the tests derive them from that module, so
the two cannot drift apart.

### §15.2 template flags that the cfg files do not use

| Flag | Why |
|---|---|
| `--met-excision` | **left out** (decision 7). GIAP's digest keeps the initiator Met (`src/giap/mapping/digest.py` decision 1), so letting DIA-NN excise it would produce peptides GIAP never proposes and break the A1/A2 comparison. |
| `--lib ""` | not needed: step 1 builds no library of its own, and step 2 gets `--lib <predicted>` from the script. |
| `--f`, `--fasta`, `--out`, `--out-lib`, `--threads` | supplied by `scripts/slurm/diann_pilot.sbatch` (run paths, output names and the Slurm CPU count). |

Everything else of the template is used verbatim: `--fasta-search --gen-spec-lib --predictor
--cut K*,R* --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --unimod4 --var-mods 1
--var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4 --qvalue 0.01
--report-decoys`.

Added beyond the template, from Khoa's reference invocation: the mass ranges `--min-pr-mz 300
--max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800`, `--peptidoforms`, `--verbose 1`, and (step 2
only) `--matrices` and `--gen-spec-lib --out-lib <arm>_empirical.parquet`. His `--reanalyse`
(MBR) and `--rt-profiling` are deliberately **not** used: match-between-runs weakens run-level
FDR control, so the identifications are no longer the single-pass q ≤ 0.01 set that §15.5 metric
8 estimates a false-discovery proportion for from entrapment peptides.

## How `db_a2` was built

```bash
giap makedb --annotation tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz \
            --genome local/data/genomes/nissle1917/GCA_043228065.1_ASM4322806v1_genomic.fna.gz \
            --outdir validation/pilot/db_a2 --no-contaminants
gzip -9n validation/pilot/db_a2/strain_db.faa
```

Counts: 4,620 CDS + 72 sORFs = 4,692 proteins (no six-frame ORFs, no decoys, no contaminants;
DIA-NN adds `contaminants.fasta` as a second `--fasta`, M3-23 Q2). `strain_db.faa.gz` is 989 KB.

The recorded Bakta JSON carries no nucleotides (R11), so `--genome` is required; `build_strain_db`
passes it through to the Bakta parser (decision 8). The genome is **GCA_043228065.1**
(`ASM4322806v1`: chromosome CP171242.1, plasmids CP171243.1, CP171244.1, CP171245.1) and is not
committed. Fetch it with `sbatch scripts/slurm/fetch_resources.sbatch` (it lands in
`$GIAP_DATA_DIR/genomes/nissle1917/`, md5 `ff062688cf999ade769f4a27c2cb0949`), or directly from
NCBI:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/228/065/GCA_043228065.1_ASM4322806v1/GCA_043228065.1_ASM4322806v1_genomic.fna.gz
```

## Raw data

PXD078459 (ProteomeXchange; hosted by iProX as IPX0017169000) is the pilot's raw data.

- **Record:** https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD078459&outputMode=XML,
  retrieved 2026-10-24 (02:43:27 UTC), sha256
  `5d3f16c3137fb11ecae7f79ad02063d8e3008b21380784fcfc8aa77b7fb533f9` (Slurm job 1879851, `--resolve`).
- **Licence:** the record states no data licence (no licence cvParam in the XML); ProteomeXchange
  data are public. To confirm with Khoa before any redistribution (none is planned: nothing raw is
  committed).
- **Instrument:** timsTOF HT (MS:1003404, from the record).
- **Selected runs:** the first three wild-type runs by natural sort (M10-01 Q2), WT_1–WT_3:
  `MLA1914689_WT_P1_1_Y-D7_1_13119.d.zip`, `MLA1914690_WT_P2_2_Y-D8_1_13120.d.zip`,
  `MLA1914691_WT_P3_3_Y-D9_1_13121.d.zip` (about 1.0 GB each).
- **Not selected** (resolve job 1879851, `file`, reason):

  | file | reason |
  |---|---|
  | `20260203_201309_sp_MJ20260114273_PP20260127494.sne` | not_raw (a Spectronaut file) |
  | `MLA1914692_DGC_P1_4_Y-D10_1_13122.d.zip` … `MLA1914694_DGC_P3_6_Y-D12_1_13124.d.zip` (3) | not_wild_type |
  | `MLA1914698_DGC_ROS_P1_10_Y-E4_1_13128.d.zip` … `MLA1914700_DGC_ROS_P3_12_Y-E6_1_13130.d.zip` (3) | not_wild_type |
  | `MLA1914695_WT_ROS_P1_7_Y-E1_1_13125.d.zip` … `MLA1914697_WT_ROS_P3_9_Y-E3_1_13127.d.zip` (3) | beyond_three |

  The `WT_ROS` runs are wild type under ROS treatment and also carry the `WT` token. The selector
  does not tell them apart by name: P1–P3 are the untreated wild-type runs only because the
  untreated set comes first in the MLA numbering (MLA1914689–91 before MLA1914695–97). This was
  checked by hand against the file names.
- **Transport (R6):** the record lists plain-HTTP links on `download.iprox.org`. That host's
  certificate does not name it, so HTTPS fails verification there. iProX's own domain
  `download.iprox.cn` (www.iprox.org redirects to www.iprox.cn) serves the same files over verified
  HTTPS: same sizes and same first MiB (checked 2026-10-24). The manifest's `url` is that HTTPS
  route; `url_record` keeps the record's link. TLS verification is never disabled.
- **Manifest:** `validation/pilot/inputs/pxd078459_files.tsv`. The record gives no checksums, so integrity rests on
  TLS only; `sha256` is filled by the download job (`SHA256SUMS` next to the data).
- **Download:** Slurm job 1879859 (2026-10-24, 23 min, COMPLETED) fetched the three archives over
  HTTPS (about 2.9 GB), unpacked them to `$GIAP_DATA_DIR/pxd078459/` as WT_1.d, WT_2.d and WT_3.d and wrote
  `SHA256SUMS`; the archive sha256s are in the manifest.
- **Acquisition mode:** diaPASEF for all three runs (`acquisition_mode.py`, Slurm job 1879882):
  734 MS1 frames and 8804–8805 diaPASEF frames (MsMsType 9) per run.
  `validation/pilot/inputs/acquisition.json` records it.

## Searches (M10-02)

DIA-NN 2.3.1 Academia ran from the read-only copy `$GIAP_ROOT/images/diann_2.3.1.sif`.
- **The install:** `scripts/slurm/diann_image.sbatch`, Slurm 1879883, copied and hashed the image
  but ended FAILED with exit 2: it ran before the banner fix (e7edc457), and its banner check read
  the blank line the binary prints first. The copy was verified afterwards: mode 444, sha256 equal
  to the pin in `validation/pilot/diann_image.sha256`. Its `diann_2.3.1.sif.json` keeps the
  pre-fix banner shape. The banner itself was checked by the searches (`diann_version.txt` of
  1879955/1879956). The install job now checks the banner on the temporary copy, before it is
  renamed into place, so a failed banner installs nothing.
- **Step 1 (the predicted libraries):** built by 1879913 (A1) and 1879914 (A2), submitted from
  commit `26adc18c`. 1879955/1879956 found a matching `step1.done` and reused them (`step1:
  skipped`), so their `command_step1.txt` is that of 1879913/1879914. Step 1 passes no raw files
  and its code is the same at both commits. `step1.done` now records the job and commit that built
  the library, and a skipped step 1 names them in `provenance.json`.
- **Step 2 (the searches):** both arms ran at commit `81272fc0` of `m10/diann-searches` from a
  clean worktree (spec Q10), on the three diaPASEF runs of PXD078459 (`acquisition.json`).
- **CPUs:** these jobs asked for `--mem=64G` (65,536 MB), more than 16 x `MaxMemPerCPU` (4,000 MB),
  so Slurm gave each 17 CPUs per task (18 allocated) and DIA-NN ran with `--threads 17`. The
  script now asks for `--mem=64000M`, which gives exactly 16.

| arm | Slurm job | node | state | elapsed | MaxRSS |
|---|---|---|---|---|---|
| A1 | 1879955 | rye-02 | COMPLETED, exit 0 | 00:14:21 | 6.6 GiB |
| A2 | 1879956 | rye-02 | COMPLETED, exit 0 | 00:14:21 | 6.4 GiB |

- **Outputs:** they stay under `$GIAP_ROOT/pilot/<arm>/`: report, library, `provenance.json`,
  `outputs.sha256`, and both command lines. Nothing derived from the reports is committed (Khoa,
  2026-10-23). Report row counts wait for Q-K1; the default is that they are not committed.
- **Superseded attempts:**
  - 1879913 (A1) and 1879914 (A2) ran concurrently without `--temp`;
  - A2 failed with "library size mismatch", because the two arms shared `<run>.d.quant` files next to the raw runs;
  - both arms were rerun with a per-arm `--temp`.
- **One submission per arm at a time:** two concurrent jobs of the *same* arm delete each other's
  `quant/` and outputs.

## Conversion (M10-03)

`scripts/slurm/diann_convert.sbatch <arm>` ran `giap convert --from diann` on each arm's DIA-NN
2.3.1 report, at commit `fdf586be` of `m10/diann-convert-tests` from a clean worktree. The job
converts twice and compares the hashes (determinism), then runs the five contract checks
(`tests/contract/test_m10_03_diann_contract.py`, 0 skipped).

```
giap convert --from diann $GIAP_ROOT/pilot/A1/A1_report.parquet --out $GIAP_ROOT/pilot/inputs/a1.tsv.gz \
  --abundance-out $GIAP_ROOT/pilot/inputs/a1_abundance.tsv.gz --search-db other
giap convert --from diann $GIAP_ROOT/pilot/A2/A2_report.parquet --out $GIAP_ROOT/pilot/inputs/a2.tsv.gz \
  --abundance-out $GIAP_ROOT/pilot/inputs/a2_abundance.tsv.gz --search-db same
```

| arm | Slurm job | node | state | elapsed | MaxRSS |
|---|---|---|---|---|---|
| A1 | 1880219 | wheat-01 | COMPLETED, exit 0 | 00:00:17 | 0.3 GiB |
| A2 | 1880220 | wheat-05 | COMPLETED, exit 0 | 00:00:17 | 0.3 GiB |

sha256 of each input report and of the two tables written from it:

| arm | file | sha256 |
|---|---|---|
| A1 | `A1_report.parquet` (input) | `f217a40fc0777e99783f3d410743bc80d301290dcdb219a6c69574ee6a4926d7` |
| A1 | `a1.tsv.gz` | `81818184f0cd222f091ae6b835ff72cc5743b6ce7e6403e7a2cb1740bbf24f3a` |
| A1 | `a1_abundance.tsv.gz` | `20284d4234055b26cf779dc405548087e4efc9bf1007dc3bebc031c5b6c74c54` |
| A2 | `A2_report.parquet` (input) | `733a3884a15664f13385b241e48186e538308993ff2d5a8cc3ec9ff212866d36` |
| A2 | `a2.tsv.gz` | `1eb5e1e9a8855fc8ac62bed4d60f2715f9270f2b96f40206dd35a3a754adf14f` |
| A2 | `a2_abundance.tsv.gz` | `c61122f3afc9aa31986de8bbfef1c542719627d3d632ec9113e7719582b41c19` |

- **Where the tables live:** they are derived from the DIA-NN reports, so they stay under
  `$GIAP_ROOT/pilot/inputs/` and are not committed (Khoa, 2026-10-23). This departs from §15.2
  step 7; the amendment goes to ADR-0051.
- **Row counts** (data rows, header excluded; aggregates are allowed, Khoa 2026-10-24):

  | arm | peptide rows (`aN.tsv.gz`) | protein groups (`aN_abundance.tsv.gz`) |
  |---|---|---|
  | A1 | 81,278 | 2,878 |
  | A2 | 93,781 | 3,257 |
- **Reproducing:** conversion is deterministic, so rerunning the job and matching these hashes is a
  full check (M10-14).
- **Superseded runs:** 1880208 and 1880209 ran at `b051c2b7`, before main's
  `acquisition.json`, so one contract check was skipped. Their outputs were overwritten by the
  runs above.

## A3 de novo (M10-08)

Arm A3 runs Cascadia 0.0.7 (Apache-2.0, PyPI wheel pinned by sha256 in `workflow/envs/cascadia.yml`)
on per-window mzML converted with TIMSCONVERT 2.0.0 and pyTDFSDK (Apache-2.0; the bundled Bruker
TDF-SDK is used under Bruker's licence, accepted by Khoa, and never redistributed: ADR-0054). The
weights (`cascadia.ckpt`, Apache-2.0, from the authors' Drive folder) and the Unimod table (from
psims 1.4.0) are pinned by sha256 and bytes in `a3_assets.tsv`; nothing of either is committed.
The env locks are `workflow/envs/locks/{cascadia,timsconvert}.{explicit,pip}.txt`.

```bash
# the env: the PyTorch 2.0 CUDA 11.8 build only solves on a CPU node with the virtual package set
CONDA_OVERRIDE_CUDA=11.8 sbatch scripts/slurm/build_envs.sbatch cascadia timsconvert
A3_WEIGHTS_HOST_OK=1 sbatch scripts/slurm/fetch_a3_assets.sbatch     # verified copies (R6)
sbatch scripts/slurm/mzml_pilot.sbatch WT_1                          # .d -> mzML, then mzml-check
RUN=WT_1 sbatch --time=04:00:00 scripts/slurm/denovo_pilot.sbatch cascadia   # 1 GPU, offline Unimod
python validation/pilot/a3_feasibility.py evaluate $GIAP_ROOT/pilot/A3_cascadia/WT_1/WT_1.ssl \
  --run WT_1 --out validation/pilot/metrics/a3_feasibility.json \
  --route B2 --route-version <commit> --weights-sha256 <sha> --mzml-sha256 <sha> \
  --mzml-check $GIAP_DATA_DIR/pxd078459/mzml/WT_1.mzml_check.json --job-id <id> \
  --date <YYYY-MM-DD> \
  --detail $GIAP_ROOT/pilot/A3_cascadia/WT_1/feasibility_top1000.tsv
```

The verdict is rule L6 on WT_1 (spec M10-08): fail when fewer than 500 distinct 7-30-mers pass
Cascadia's own cut-off (0.8), or fewer than 20 % of the top 1,000 match the Nissle six-frame
translation exactly; exit 0 pass, 3 fail. `early-fail --reason … --evidence …` records E1-E3.
`metrics/a3_feasibility.json` holds aggregates only; the job IDs, wall times, sizes and the
verdict itself come from the smoke run (pending). `--date` defaults to today (UTC); pass it so a
rerun is byte-identical.

## Pass-1 runs (M10-04)

§15.2 step 8 for arms A1 and A2, then `giap compare` (the first and third §15.6 criteria; spec
`docs/specs/M10-04-pilot-runs.md`). This section holds provenance only: job IDs, commits, states,
exit codes, runtimes, memory, sha256s and byte sizes. No biological number is committed here (Q-K1);
those belong to M10-05's `REPORT.md`. The run directories, the comparison, the evidence JSON and
`ATTEMPTS.tsv` stay under `$GIAP_ROOT/runs/pilot/pass1/` and are not committed (Khoa,
2026-10-23). Pass-1 outputs may carry `uncalibrated` flags; that is expected, because the
M6/M7/M8 hardening is still open.

```bash
P=$GIAP_ROOT/pilot/inputs; R=$GIAP_ROOT/runs/pilot/pass1
G=$GIAP_REPO/tests/data/genomes/nissle1917.fna.gz
sbatch -J giap-run-pilot-a2 scripts/slurm/giap_run.sbatch \
  --genome $G --annotation $P/nissle1917.bakta.json --peptides $P/a2.tsv.gz \
  --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2
sbatch -J giap-run-pilot-a1 scripts/slurm/giap_run.sbatch \
  --genome $G --annotation $P/nissle1917.bakta.json --peptides $P/a1.tsv.gz \
  --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a1
sbatch -J giap-pilot-compare scripts/slurm/pilot_compare.sbatch $R/out_a1 $R/out_a2 $R/cmp_a1_a2
python validation/pilot/run_evidence.py arm $R/out_a1 --arm A1 --job JOBID \
  --expect-commit SHA --expect-input genome=SHA --expect-input annotation=SHA \
  --expect-input peptides=SHA --out $R/evidence/a1.json
```

Both arms ran with the profile `farmshare,apptainer`, the default modules (`modules_off` empty),
the default seed 20261014, the default `--genome-relation related`, no `--abundance` and no
`--strict`. Both used the recorded Bakta 1.12.2 annotation (DB 6.0, full), decompressed once to
`$GIAP_ROOT/pilot/inputs/nissle1917.bakta.json` (read-only), with no live `BAKTA` task. Both
used detection tier `decoy_picked`. The accepted attempts ran at commit
`6f969ee79d579d4bdb429be6f156679d644b2b51` from the clean pinned worktree `pilot-6f969ee7`
(`git_dirty` false).

### Inputs (sha256 from `run_manifest.json`)

| role | file | sha256 |
|---|---|---|
| genome | `tests/data/genomes/nissle1917.fna.gz` | `fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee` |
| annotation | `nissle1917.bakta.json` (13,668,551 bytes) | `cf7de472a179ae122c085d33995e1b8399048a24d250b15e4c3d82e5e41f1ee3` |
| peptides A1 | `a1.tsv.gz` | `81818184f0cd222f091ae6b835ff72cc5743b6ce7e6403e7a2cb1740bbf24f3a` |
| peptides A2 | `a2.tsv.gz` | `1eb5e1e9a8855fc8ac62bed4d60f2715f9270f2b96f40206dd35a3a754adf14f` |

The peptide hashes equal the M10-03 conversion table above. The assembly has a chromosome and
three circular plasmid records, not the "pMUT1 and pMUT2" named in §15.1. The runs use the assembly
as committed, and M10-05 reports plasmid genes per replicon, not per plasmid name.

### Attempts (every submission)

| arm | attempt | head job | class | commit | state | exit | elapsed | MaxRSS (head) | node | rule set | outcome |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A2 | 1 | 1892236 | fresh | `0667e366` | COMPLETED | 0:0 | 01:52:31 | 1200512K | barley-03 | A: `commit_mismatch` | superseded; kept as `out_a2_0667e366` |
| A1 | 1 | 1892719 | fresh | `6f969ee7` | COMPLETED | 0:0 | 01:37:08 | 1384424K | wheat-02 | A: pass | **accepted**, `out_a1` |
| A2 | 2 | 1892882 | fresh | `1549a930` | COMPLETED | 0:0 | 03:08:09 | 1229680K | wheat-02 | A: `commit_mismatch` | superseded; kept as `out_a2_1549a930` |
| A2 | 3 | 1893676 | fresh | `6f969ee7` | COMPLETED | 0:0 | 01:43:45 | 1233216K | wheat-08 | A: pass | **accepted**, `out_a2` |
| compare | 1 | 1894360 | fresh | `5d1cf8d3` | COMPLETED | 0:0 | 00:00:17 | 268924K | wheat-08 | C: pass | **accepted**, `cmp_a1_a2` |

- **A2 attempt 1 (job 1892236):** a main merge landed between the two arms, so this attempt ran at
  `0667e366` while A1 ran at `6f969ee7`. Rule A3 fails, so the attempt is superseded.
- **A2 attempt 2 (job 1892882):** `giap_run.sbatch` set `GIAP_REPO` to the pinned worktree, but it
  imported the shared venv's editable `giap` from the main checkout, so the manifest records
  `1549a930`. The diff `6f969ee7..1549a930` touches only the cascadia and timsconvert lock files,
  but the gate stands (R3) and rule A3 fails. Since this run, `giap_run.sbatch` on main pins
  `PYTHONPATH` to `GIAP_REPO`.
- **A2 attempt 3 (job 1893676):** `PYTHONPATH` pinned to the worktree; rule set A passes.
- A2 used all three of its R9 attempts (two commit mismatches, then a pass); A1 used one. No
  attempt was a `--resume`.
- Both superseded directories were renamed from `out_a2` after their runs (their manifests'
  `command_line` says `out_a2`). Their rule set A evidence is `evidence/a2_0667e366.json` and
  `evidence/a2_1549a930.json`.
- The compare job was submitted from the `spec-m10-04` worktree with `GIAP_REPO` exported. Its
  head was `5d1cf8d3` at submission; the job itself does not record a commit.

### Resources (§15.5, the runtime and memory metric)

| arm | head job | task max AllocCPUS | task peak MaxRSS | top steps by wall time (`measure_runtimes.py`) |
|---|---|---|---|---|
| A1 | 1892719 | 14 CPUs | 10,597,531,648 bytes | terminatornet 1:08:26, promoter_calculator 0:17:17, determinants_step 0:10:24, sensitivity 0:10:07, function_step 0:02:10 |
| A2 | 1893676 | 14 CPUs | 10,301,911,040 bytes | terminatornet 1:10:54, promoter_calculator 0:18:47, determinants_step 0:15:06, sensitivity 0:09:22, promotech 0:03:39 |

No task failed or was retried in either accepted arm. Both head runs took under 2 h, well below
MVP-04's expected range of 3 h to 8 h; TerminatorNet dominates the wall time.

### Reports

| subject | file | sha256 | bytes |
|---|---|---|---|
| A1 | `out_a1/report/report.html` | `3005b2c2d2efab68eded2afb15d68c483c79461e8fcd05b86795a2e71d0f4c26` | 20,686,264 bytes |
| A2 | `out_a2/report/report.html` | `8b2ab25820ef2c551ea9d92b55d18b4bc82a9d64db390db9e46b779951b546bd` | 20,693,932 bytes |
| compare | `cmp_a1_a2/report/compare.html` | `b477574b65e4c6dd9b10971ce0923a57ae6908e904773dd1e10d882ff0f25f7d` | 5,199,817 bytes |
| compare | `cmp_a1_a2/plots/V16_compare.svg` | `91583ab11241ac0e160841da65b88fa8962915055a98fff5a41ccc9e77a6f27d` | 253,809 bytes |
| compare | `cmp_a1_a2/tables/compare_detection.tsv` | `ed2998b09c790f4053dcb6b4b6038dc3abcb82a8b69042312403ebaba2c2abb6` | 422,450 bytes |
| compare | `cmp_a1_a2/tables/compare_abundance.tsv` | `3aa654ebe3228517b27e8988e6ec551752cc0d034923f6504031e2467992e418` | 214 bytes |
| compare | `cmp_a1_a2/tables/compare_stability.tsv` | `ba83f391498ee8206c64b40ed0f1b3351731208090500a7a1b0d594b3a3bcc81` | 230 bytes |

The second `giap compare` (into `$TMPDIR/cmp2`) was byte-identical (rule C6).

### Evidence verdicts

`evidence/a1.json`, `evidence/a2.json` and `evidence/compare.json` all have `failures: []`.
Warnings are recorded, not failed:

- Both arms warn `novel_orfs.decoys_not_makedb` and `region_layer.layer_missing`.
- A2 also warns `novel_orfs.layer_missing`.
- `layers_missing` lists `antismash` in both arms. The head log says that the precomputed
  `bakta_json` has no GBFF next to it, so antiSMASH was skipped. The `bgc` module was on, but the
  antiSMASH layer did not run in these runs (a backlog finding, not a rule failure).
- A2 also lacks `rbs_novel`, and its novel-ORF promoter and RBS layers were skipped (`not_run`).

## A3 run (M10-10)

§15.2 steps 6-8 for arm A3 (Cascadia 0.0.7 de novo), then the three-way `giap compare` (the A3
part of the first §15.6 criterion; spec `docs/specs/M10-10-a3-run.md`). As in the pass-1 section,
this section holds provenance and A3-only aggregates: job IDs, commits, states, runtimes, memory,
sha256s, byte sizes, `.ssl` row counts and the per-run L6 numbers. Metric 6 (the A2/A3 overlap)
and hypothesis 4 involve A2 data and live in `metrics.json`/`REPORT.md` (M10-05's Q-K1, spec Q10).
The `.ssl` files, `a3.tsv.gz`, `a3.convert.*`, the per-run `feasibility.json`, `out_a3/`,
`cmp_a1_a2_a3/`, the evidence JSON and the `pass2/` attempts table stay under `$GIAP_ROOT` and are not
committed.

### Deviation: CPU-route `.ssl`, no GPU jobs

The spec plans three fresh GPU Cascadia jobs (steps 1-2, `gpu` partition, written to
`$GIAP_ROOT/pilot/A3/`). They were **not run**. A3 instead uses the `.ssl` files that M10-08's
CPU route (B2, Slurm work queue, no GPU) had already produced for all three runs, under
`$GIAP_ROOT/pilot/A3_cascadia_cpu/` (`A3_OUT_DIR` of the convert job). The command, weights
(sha256 `0f5e9141…`, `a3_assets.tsv`), Unimod table, batch size and the default cut-off 0.8 are
M10-08's, and no other checkpoint or threshold was tried (R3). The deviation is recorded in the
`convert` row of the `pass2/` attempts table and again in the A3 row. The consequences:

- The GPU sizing rule (step 2) and the GPU model were not exercised: no `giap-denovo-WT_*` job and
  no `nvidia-smi` line.
- WT_1 was not re-run, so there is no smoke-vs-rerun comparison. The WT_1 input is the M10-08
  smoke `.ssl` itself: its L6 aggregates below equal `metrics/a3_feasibility.json`, including
  the row count.
- Provenance is not uniform across runs (spec Q1). WT_1 came from the M10-08 work-queue arrays
  (1889426_0, 1890599, 1890876, 1890952, 1891287, 1891944; merge 1890902). WT_2 and WT_3 came from
  later work-queue jobs, merged by 1891845 and 1892227.

```bash
A3=$GIAP_ROOT/pilot/A3_cascadia_cpu; P=$GIAP_ROOT/pilot/inputs; R2=$GIAP_ROOT/runs/pilot/pass2
A3_OUT_DIR=$A3 sbatch -J giap-cascadia-convert scripts/slurm/cascadia_convert.sbatch
sbatch -J giap-run-pilot-a3 scripts/slurm/giap_run.sbatch \
  --genome $GIAP_REPO/tests/data/genomes/nissle1917.fna.gz --annotation $P/nissle1917.bakta.json \
  --peptides $P/a3.tsv.gz --reference-peptides $P/a2.tsv.gz \
  --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R2/out_a3
python validation/pilot/run_evidence.py arm $R2/out_a3 --arm A3 --job 1894682 \
  --expect-commit 242ce506396b2f03f54d70efd90309efe8cd8bc7 --expect-input genome=SHA \
  --expect-input annotation=SHA --expect-input peptides=SHA --expect-input reference_peptides=SHA \
  --a2 $GIAP_ROOT/runs/pilot/pass1/out_a2 --out $R2/evidence/a3.json
sbatch -J giap-pilot-compare3 --cpus-per-task=2 --mem=8000M --time=01:00:00 \
  scripts/slurm/pilot_compare.sbatch $GIAP_ROOT/runs/pilot/pass1/out_a1 \
  $GIAP_ROOT/runs/pilot/pass1/out_a2 $R2/out_a3 $R2/cmp_a1_a2_a3
```

The A3 run used the same settings as A1 and A2: profile `farmshare,apptainer`, default modules
(`modules_off` empty), seed 20261014, `--genome-relation related`, no `--abundance`, no
`--strict`, the same decompressed Bakta JSON, and no `--resume`. It ran at commit
`242ce506396b2f03f54d70efd90309efe8cd8bc7` from the clean pinned worktree `pilot-242ce506`
(`git_dirty` false; `PYTHONPATH` pinned by `giap_run.sbatch`). Its detection tier is
`two_peptide`, because de novo rows carry no decoys. A1 and A2 used `decoy_picked`.

### mzML and `.ssl` (per run)

| run | mzML job | mzML bytes | mzML sha256 | mzml-check | `.ssl` bytes | `.ssl` data rows | `.ssl` sha256 |
|---|---|---|---|---|---|---|---|
| WT_1 | 1884087 | 1,472,931,770 | `691556c139e9ebfc31f64413541ea39eaa965f78c2b1b0bde1510fb145bcc75f` | M1-M8 pass | 221,559 | 1,445 | `601344b23f89c3caf8142b4a2569187f1a40a8af2e3017377375267734a3712c` |
| WT_2 | 1889330 | 1,442,505,322 | `da70d21e28747d5af458cbeab518c68726304dba95254f9e472a42f248ec82f5` | M1-M8 pass | 212,204 | 1,383 | `59706a604ece87b6f4b6938eb5dca4093b56349aedcd6c9f97f2eb3876b8a433` |
| WT_3 | 1889332 | 1,472,217,635 | `8398a104311bfb54577ada2b4139d285b46e1997be41436d52861486b6b60c9b` | M1-M8 pass | 215,722 | 1,406 | `8f62a86bea3c2285cab26181eeddc3c34d9d5fbcd5aa8bf6f268b6a7d216e45a` |

### Per-run L6 (`a3_feasibility.py evaluate`, rule: ≥ 500 distinct 7-30-mers above 0.8 and ≥ 0.20 of the top 1,000 matching the six-frame translation)

| run | rows | 7-30-mers | n_pass (= n_top) | n_match | share (95 % CI) | decoy share | contaminant share | verdict |
|---|---|---|---|---|---|---|---|---|
| WT_1 | 1,445 | 1,424 | 863 | 549 | 0.636 (0.604-0.668) | 0.006 | 0.006 | pass |
| WT_2 | 1,383 | 1,365 | 848 | 546 | 0.644 (0.611-0.675) | 0.008 | 0.005 | pass |
| WT_3 | 1,406 | 1,383 | 853 | 531 | 0.623 (0.589-0.654) | 0.006 | 0.006 | pass |

All three runs pass L6, so no run is excluded and all three were converted.

### Conversion (`cascadia_convert.sbatch`, job 1894680)

The job completed (exit 0:0) in 00:00:25 on wheat-07 (2 CPUs, MaxRSS 237376K) at commit `242ce506`
(clean). `a3.convert.json` has `failures: []`, and the second conversion was byte-identical.

| item | value |
|---|---|
| `a3.tsv.gz` | 63,221 bytes, sha256 `743396d090db4f7b0e1075d63d68c4a81b903b8b6fd3a8b61d80fedcac1bf0b8` |
| rows (= sum of `.ssl` data rows) | 4,234 (WT_1 1,445, WT_2 1,383, WT_3 1,406) |
| distinct peptides (7-30-mers) | 1,386 (1,364) |
| charge 2 / 3 | 4,153 / 81 |
| score min / median / max | 0.800034 / 0.930794 / 0.999955; 0 rows at or below 0.8 |
| log issues | `short_peptide` warn (62 rows), `stripped` info (1 row) |

### Attempts (every submission, the `pass2/` attempts table)

| arm | attempt | job | class | commit | state | exit | elapsed | MaxRSS (head) | node | rule set | outcome |
|---|---|---|---|---|---|---|---|---|---|---|---|
| convert | 1 | 1894680 | fresh | `242ce506` | completed | 0:0 | 00:00:25 | 237376K | wheat-07 | L6 | **accepted**; CPU-route `.ssl` (deviation above) |
| A3 | 1 | 1894682 | fresh | `242ce506` | completed | 0:0 | 03:01:03 | 1303576K | wheat-05 | A: pass | **accepted**, `out_a3` |
| compare3 | 1 | 1895328 | fresh | `242ce506` | completed | 0:0 | 00:00:19 | 317544K | wheat-06 | C: pass | **accepted**, `cmp_a1_a2_a3` |

A3 used one of its three R9 attempts. Before each submission, the R14 check summed the CPUs of
running and pending `giap-*` jobs: at the compare submission they held 18 CPUs, and 18 + 2 ≤ 128.
The compare job was submitted from the main checkout at `242ce506`.

### Inputs (sha256 from `run_manifest.json`)

| role | file | sha256 |
|---|---|---|
| genome | `tests/data/genomes/nissle1917.fna.gz` | `fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee` |
| annotation | `nissle1917.bakta.json` (13,668,551 bytes) | `cf7de472a179ae122c085d33995e1b8399048a24d250b15e4c3d82e5e41f1ee3` |
| peptides A3 | `a3.tsv.gz` | `743396d090db4f7b0e1075d63d68c4a81b903b8b6fd3a8b61d80fedcac1bf0b8` |
| reference peptides | `a2.tsv.gz` (A2's `peptides` input, rule A11) | `1eb5e1e9a8855fc8ac62bed4d60f2715f9270f2b96f40206dd35a3a754adf14f` |

### Resources (§15.5)

| arm | head job | tasks | task max AllocCPUS | task peak MaxRSS | top steps by wall time (`measure_runtimes.py`) |
|---|---|---|---|---|---|
| A3 | 1894682 | 77 (0 failed, 0 retried) | 14 CPUs | 11,809,894,400 bytes | terminatornet 2:54:44, promotech (11 tasks) 0:29:46, promoter_calculator (27 tasks) 0:16:10, genomad 0:02:34, promotech_novel 0:01:10 |

The head run took 3 h 01 min, inside MVP-04's 3-8 h range. TerminatorNet dominates again. Its wall
time was about 2.5 times A1's (1:08:26) on the same genome and settings, so it depends on the node.
This is a runtime finding, not a rule.

### Reports

| subject | file | sha256 | bytes |
|---|---|---|---|
| A3 | `out_a3/report/report.html` | `289a13e262f510944e2f5bb38863ff01d73df9661222f1caecb12e43e8649737` | 11,121,121 bytes |
| compare | `cmp_a1_a2_a3/report/compare.html` | `9c793c495504cd4b19ef9f833a0c88902a2dd841a27b253f9b5ad728a7b56477` | 5,206,317 bytes |
| compare | `cmp_a1_a2_a3/plots/V16_compare.svg` | `6e7dc23b6e19ee1ae6083e385dc3308532670c97b60f61aebd88a60ac505cc4b` | 257,607 bytes |
| compare | `cmp_a1_a2_a3/tables/compare_detection.tsv` | `d99450795e5018bdb29413008d2a4377abcae3c2908475a9fa5a2fa95533bba3` | 650,081 bytes |
| compare | `cmp_a1_a2_a3/tables/compare_abundance.tsv` | `64732f1a1baa67b923dd2a50a390b77506c8629e5e6b01a7ead069c6854a7e70` | 354 bytes |
| compare | `cmp_a1_a2_a3/tables/compare_stability.tsv` | `8cbadf2f151fb1a07eeb823cbfc21ab8bc9474050c022774e64ea7858419c8da` | 599 bytes |

The second `giap compare` (into `$TMPDIR/cmp2`) was byte-identical (rule C6).

### Evidence verdicts

`pass2/evidence/a3.json` and `pass2/evidence/compare3.json` both have `failures: []`.

- **A3 (rule set A, A1-A11):**
  - A10 passes: `mapping_fdr.tsv` has 10 `denovo_precision` rows and none is flagged `c_assumed`.
  - A11 passes: the reference sha256 equals `out_a2`'s `peptides` input.
  - `expected_warnings` holds `no_spectrum_pairing`, as expected (DIA-NN writes no
    `spectrum_id`), so ALT class error rates are `unestimated` (§7.5.4).
  - The manifest warns `novel_orfs.class_q_floor_above_alpha`, `peptides.content`,
    `peptides.spectrum_pairing` and `region_layer.layer_missing`.
  - `layers_missing` lists `abundance` (A3 has no intensities) and `antismash` (no GBFF next to the
    precomputed `bakta_json`, as in A1 and A2).
- **Compare (rule set C, three labels):**
  - The two `compare_abundance` rows with A3 carry `missing_in:A3`, as expected.
  - C7 raised no warning: the (A1, A2) rows of `compare_abundance` and `compare_stability` equal
    those of `pass1/cmp_a1_a2`.
  - A1/A2 ran at `6f969ee7` and A3 at `242ce506`. This is recorded, not failed (spec Q4).

## Pass-3 runs: wild-type reference (M10 re-run)

Passes 1 and 2 used **GCA_043228065.1** as the Nissle reference. That assembly is the engineered
chassis XTC004-ILA (BioSample SAMN44054533), not wild-type Nissle 1917. Khoa decided on
2026-10-25 ("go with your recommendation") to re-run the pilot against an official wild-type
assembly. Passes 1 and 2 stay as earlier attempts. Pass 3 re-did every genome-dependent step
(Bakta, the A2 database, the A2 search and conversion), then the GIAP runs of A1 and A2 and their
`giap compare`. A3 was not re-run in pass 3.

Like the pass-1 section, this section holds provenance only. The run directories, the evidence
JSON, the compare, the attempts table and the interim metrics stay under `$GIAP_ROOT` and are not
committed:

- runs: `$GIAP_ROOT/runs/pilot/pass3/`;
- inputs: `$GIAP_ROOT/pilot/pass3/`;
- the interim `metrics.json` and `REPORT.md` (rendered without `--final`):
  `$GIAP_ROOT/pilot/metrics/pass3/`.

M10-15 owns the committed `metrics/metrics.json` and `REPORT.md`.

### Reference genome

| item | value |
|---|---|
| assembly | **GCF_021559835.1** (ASM2155983v1, RefSeq; GenBank twin GCA_021559835.1) |
| records | chromosome `NZ_CP058217.1`, plasmid pMUT1 `NZ_CP058218.1`, plasmid pMUT2 `NZ_CP058219.1` |
| source | NCBI (`ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/559/835/`), R6 |
| GCA download | md5 `3c1d874a2172b3bfd63a05fb651ab7d2`, equal to NCBI's `md5checksums.txt` |
| GCF copy | re-gzipped; decompressed, it is byte-identical to NCBI's GCF file (NCBI md5 `eaa414c72e26fc9dd49dacd4ca256283`, re-checked 2026-10-26). Its sequences equal the GCA download's |
| GIAP input | `nissle1917_wt.fna.gz`: the GCF records with ` [topology=circular]` appended to each header, sha256 `1a654e6b83292353f15be4e7a7973b896d728761a4c9dfde533953f76e7961ad` |

The two plasmids are the pMUT1 and pMUT2 that §15.1 names, so the pass-1 note about three unnamed
plasmid records does not apply to pass 3. The genome is not committed. A circular-tagged copy sits
untracked in the staging worktree, and the reference copy is `$GIAP_ROOT/pilot/pass3/genome/`.

### Annotation, strain database, A2 search and conversion

| step | Slurm job | node | state | elapsed | CPUs |
|---|---|---|---|---|---|
| Bakta 1.12.2 full DB 6.0 (`record_bakta_full.sbatch`) | 1895437 | rye-02 | COMPLETED, exit 0 | 01:42:40 | 18 CPUs |
| strip, decompress, `giap makedb`, L6 on WT_1 | 1895583 | wheat-05 | COMPLETED, exit 0 | 00:00:23 | 2 CPUs |
| DIA-NN 2.3.1 A2, steps 1 and 2 | 1895586 | rye-02 | COMPLETED, exit 0 | 00:42:00 | 16 CPUs |
| `diann_report.sbatch A2` checks | 1895667 | wheat-08 | COMPLETED, exit 0 | 00:00:08 | 4 CPUs |
| `diann_convert.sbatch A2` | 1895670 | wheat-01 | COMPLETED, exit 0 | 00:00:28 | 4 CPUs |

- **Bakta:**
  - Bakta 1.12.2 ran from image `sha256:1e24ae2f`, with `--keep-contig-headers --replicons` on the
    preflight-normalised genome and 16 threads.
  - The JSON was stripped of `sequences.nt` (`import_recording.py --strip-json-key`, as in M2-10):
    `nissle1917_wt.json.gz`, sha256 `b680d455d8760a23476edf4d624dd5e816bba0df387098812354ea13adc7f223`.
  - It was decompressed once to `nissle1917_wt.bakta.json` (13,680,327 bytes, sha256
    `dbe5e126ede2d027e63f1f57b7abc6475238718dcbfd4468ee0bbecef2b93a00`, read-only).
- **Strain database:**
  - Built as in "How `db_a2` was built" above, with `--genome` set to the GCF copy.
  - Contents: 4,613 CDS + 72 sORFs = 4,685 proteins, with no contaminants and no decoys.
    `strain_db.faa.gz` is 989,209 bytes, sha256
    `20e36adfa049110906bc688d20ea83f8ca251f9afcc7add18d78ca31029417b8`.
  - It is **not** committed: the committed `db_a2/` stays the pass-1 database.
  - DIA-NN and the conversion contract tests read the database from the staging worktree
    `pilot3-stage-4ccce8a532c3`, whose `db_a2/strain_db.faa.gz` was replaced and left
    uncommitted. That is why `a2.convert.json` records `dirty: true`.
- **A2 search:**
  - DIA-NN 2.3.1 used the same image pin and an unchanged `diann_A2.cfg` (sha256 `c0d145e8`).
  - Step 1 built a new predicted library from the pass-3 database, and step 2 searched the three
    runs. Report checks: `ok`, `failures: []`.
  - `A2_report.parquet` sha256 `c274bbb028b386c67d72410fd5ae48a6eaf0edeb7ac2a8030ab8bbf73ea6bc06`.
- **Conversion:** deterministic (the two conversions were equal under `cmp`), and the five contract
  tests passed with 0 skipped. `a2.tsv.gz` is 2,054,199 bytes, sha256
  `d7f1ddf5d2badb0038b50f55a6d99b43be61ccc9e2c8fadff9f5a7ab0c405fbf`, with 94,363 peptide rows.
  `a2_abundance.tsv.gz` has 3,250 protein groups, sha256
  `3970344b9c0edc0fcf756b267063a1d1c30c05f6b6bdb6654d973eae1e6f888e`.
- **A1** reuses the pass-1 `a1.tsv.gz`. The K-12 search does not depend on the Nissle genome.
- **A3 L6 on the wild-type genome:** WT_1's `.ssl` was re-evaluated against the GCF genome and the
  pass-3 Bakta JSON, and the verdict is **pass** (`$GIAP_ROOT/pilot/pass3/a3_feasibility.json`).
  Its aggregates equal the pass-2 WT_1 row above.

### GIAP runs

```bash
P=$GIAP_ROOT/pilot/pass3; R=$GIAP_ROOT/runs/pilot/pass3
G=$P/genome/nissle1917_wt.fna.gz; J=$P/inputs/nissle1917_wt.bakta.json
sbatch -J giap-run-pilot3-a2 --time=12:00:00 scripts/slurm/giap_run.sbatch --genome $G --annotation $J \
  --peptides $P/inputs/a2.tsv.gz --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2
sbatch -J giap-run-pilot3-a1 --time=12:00:00 scripts/slurm/giap_run.sbatch --genome $G --annotation $J \
  --peptides $GIAP_ROOT/pilot/inputs/a1.tsv.gz --db-dir $GIAP_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a1
sbatch -J giap-pilot3-compare --cpus-per-task=2 --mem=8000M scripts/slurm/pilot_compare.sbatch \
  $R/out_a1 $R/out_a2 $R/cmp_a1_a2
```

Both arms ran with the pass-1 settings:

- profile `farmshare,apptainer`;
- default modules (`modules_off` empty), seed 20261014 and `--genome-relation related`;
- no `--abundance`, no `--strict` and no `--resume`;
- the pass-3 Bakta JSON, with no live `BAKTA` task;
- detection tier `decoy_picked`.

Both ran at commit `4ccce8a532c365c9e657c981a3f61b3d0d5c3aa6` from the clean pinned worktree
`pilot-4ccce8a532c3` (`git_dirty` false), with `GIAP_REPO` exported and `PYTHONPATH` pinned by
`giap_run.sbatch`. The compare job was submitted from the same worktree.

| role | file | sha256 |
|---|---|---|
| genome | `nissle1917_wt.fna.gz` | `1a654e6b83292353f15be4e7a7973b896d728761a4c9dfde533953f76e7961ad` |
| annotation | `nissle1917_wt.bakta.json` (13,680,327 bytes) | `dbe5e126ede2d027e63f1f57b7abc6475238718dcbfd4468ee0bbecef2b93a00` |
| peptides A1 | `a1.tsv.gz` (pass 1) | `81818184f0cd222f091ae6b835ff72cc5743b6ce7e6403e7a2cb1740bbf24f3a` |
| peptides A2 | `a2.tsv.gz` (pass 3) | `d7f1ddf5d2badb0038b50f55a6d99b43be61ccc9e2c8fadff9f5a7ab0c405fbf` |

#### Attempts (the `pass3/` attempts table)

| arm | attempt | head job | class | commit | state | exit | elapsed | MaxRSS (head) | node | rule set | outcome |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A2 | 1 | 1895673 | fresh | `4ccce8a5` | COMPLETED | 0:0 | 04:00:04 | 1285072K | barley-03 | A: pass | **accepted**; OSTIR pre-#390 (below); kept as `out_a2_pre390` |
| A1 | 1 | 1896178 | fresh | `4ccce8a5` | COMPLETED | 0:0 | 04:29:54 | 1300748K | rye-02 | A: pass | **accepted**, `out_a1` |
| compare | 1 | 1896978 | fresh | `4ccce8a5` | COMPLETED | 0:0 | 00:00:28 | 270028K | wheat-08 | C: pass | **accepted**, `cmp_a1_a2` |
| A2 | 2 | 1897025 | resume (#390) | `8b6e93d4` | COMPLETED | 0:0 | 00:55:08 | 1243668K | barley-01 | A: pass | **accepted**, `out_a2`; supersedes attempt 1 |
| A1 | 2 | 1897026 | resume (#390) | `8b6e93d4` | COMPLETED | 0:0 | 00:23:48 | 1367036K | wheat-06 | A: pass | **accepted**, `out_a1` (commit parity) |
| compare | 2 | 1897034 | fresh | `8b6e93d4` | COMPLETED | 0:0 | 00:00:45 | 452932K | wheat-04 | C: pass | **accepted**, `o390/cmp_a1_a2` |

Attempt 1 of each arm passed rule set A, and attempt 2 is the #390 resume below. Neither attempt is a failed iteration, which is what R9 counts. While A2 1895673 was pending, its head time limit was
lowered from 2 days to 12 h so that it could backfill. R14 sums of the running and pending `giap-*`
CPUs were 1 + 72 at the A2 submission, 50 + 72 at A1 and 1 + 2 at the compare, all within 128.

#### Resources (§15.5)

| arm | head job | tasks | task max AllocCPUS | task peak RSS | top steps by wall time (`measure_runtimes.py`) |
|---|---|---|---|---|---|
| A2 | 1895673 | 73 (0 failed, 0 retried) | 14 CPUs | 7.5 GB | terminatornet 1:47:10, promotech (11 tasks) 0:28:06, variants_step 0:14:47, promoter_calculator (27 tasks) 0:13:22, determinants_step 0:10:07 |
| A1 | 1896178 | 77 (0 failed, 0 retried) | 14 CPUs | 7.5 GB | terminatornet 1:58:32, sensitivity 0:43:39, promotech_novel 0:42:08, promotech (11 tasks) 0:39:44, promoter_calculator (27 tasks) 0:10:53 |

Both head runs took about 4 h, inside MVP-04's range of 3 h to 8 h. TerminatorNet dominates again,
as it did in passes 1 and 2.

#### Reports

| subject | file | sha256 | bytes |
|---|---|---|---|
| A1 | `out_a1/report/report.html` | `aa002250be3819553a6b6c8a7e7c8f7af8ea9a79eaf72fc127ac785df2909914` | 20,649,053 bytes |
| A2 (pre-#390) | `out_a2_pre390/report/report.html` | `ccb0c1ebfdcf452ab498f8b5ae8e70052036eee5304807e7152bace01f33f66d` | 20,700,021 bytes |
| compare | `cmp_a1_a2/report/compare.html` | `efb58e2886ebd243e3b05b99705bd21ab3dc3117c1f46e4e97649e7d3458cfed` | 5,201,198 bytes |
| compare | `cmp_a1_a2/plots/V16_compare.svg` | `9d5de6e87eccccd8ca1d8cd2e4560be0eb399a12e7acd1dc5e729043ce019141` | 254,719 bytes |
| compare | `cmp_a1_a2/tables/compare_detection.tsv` | `180dca4afed2b511d5a4945fa11c20ec9c62fe51b6ef04c6abf59e8a834a4e8c` | 421,641 bytes |
| compare | `cmp_a1_a2/tables/compare_abundance.tsv` | `6a0aac4a5b77f00af8f1954ec9a4dd005db205c13ee13b6bf43cc06049300822` | 214 bytes |
| compare | `cmp_a1_a2/tables/compare_stability.tsv` | `1a60b9a2f747c55225b9f66f982bc953d432d0b8de7f8b4f4f71964c8125bf03` | 230 bytes |

The second `giap compare` (into `$TMPDIR/cmp2`) was byte-identical (rule C6).

#### Evidence verdicts

`pass3/evidence/a1.json`, `a2.json` and `compare.json` all have `failures: []`.

- **Warnings:**
  - Both arms warn `novel_orfs.decoys_not_makedb` and `region_layer.layer_missing`.
  - A2 also warns `novel_orfs.layer_missing`.
  - The compare has no warnings.
- **Missing layers:**
  - `layers_missing` lists `antismash` in both arms: the precomputed `bakta_json` has no GBFF
    next to it, as in passes 1 and 2.
  - A2 also lacks `rbs_novel`, and its novel-ORF promoter and RBS layers were `not_run`, as in
    pass 1.
  - **Cause of the missing antiSMASH layer:** the run input, not the tool or the profile.
    `annotate.nf` looks for the GBFF next to the precomputed JSON, as `<stem>.gbff[.gz]`. Every pass
    passed the decompressed `inputs/<name>.bakta.json` alone. The recorded GBFFs
    (`nissle1917.gbff.gz`, `nissle1917_wt.gbff.gz`) sit next to the recorded `<name>.json.gz`, under
    another stem and in another directory. The `bgc` module is on and the antiSMASH DB is
    installed, but the process was never called.
  - **What that means for H1:** with no antiSMASH there are no BGC genes, so H1's "mostly plasmid,
    prophage and BGC genes" clause was tested without BGC genes in every pass. Its verdict rule is
    unchanged.
  - **Islands and operon flags:** islands are empty in both arms of every pass, and no operon
    flag is significant in any pass. H3's "supported" therefore rests on the rho and
    enriched-categories criteria alone. The island and operon families are `not_evaluable`; they
    do not count as support.
- **§15.6 status in the interim `metrics.json`:**
  - "runs end to end" is partial: A1 and A2 ran, and A3 was not re-run.
  - TE1 passes.
  - "A1 and A2 compared" passes.
  - The runtime and TE2-TE5 bullets are pending (M10-12, M10-06), as in pass 2.
  - The metrics have no warnings. A2's metrics were computed from the `out_a2_pre390` snapshot.

### OSTIR and the Python user site (#390)

Before #390 (merged 2026-10-26 01:15 PDT), the py3.10 tool environments imported `~/.local`
packages. OSTIR therefore folded RNA with ViennaRNA 2.7.1 instead of the locked 2.6.4. See
`validation/reports/ostir-usersite-rerecord-2026-10-26.md`: three Nissle genes change, and the
gene-level TIR changes by up to about 65x.

- **A2 (1895673): pre-fix.** Its head job was submitted before #390 and ran code at `4ccce8a5`,
  which does not have the fix. Its OSTIR `versions.yml` has no `viennarna:` line, and its `ostir.csv` holds the user-site values for
  the three genes. `rbs`, the `tir_ostir` determinant and everything computed from them
  (`determinants`, `sensitivity`, the report) are pre-fix in this attempt.
- **A1 (1896178): post-fix in effect.** It ran at the same commit, but it was submitted at 05:02,
  after #390. Most likely `PYTHONNOUSERSITE=1` came in from the submitting shell; this was not
  recorded. Its `ostir_input.csv` is byte-identical to A2's, and its `ostir.csv` differs from A2's
  in exactly the three affected rows. A1's values in those rows match the #397 re-record (the
  isolated ViennaRNA 2.6.4). That match is the only evidence, because neither log nor version file
  records which ViennaRNA OSTIR imported. Every OSTIR `.command.log` prints `Vienna version: 2.6.4`
  both before and after the fix: the line comes from `RNAfold --version`, not from the imported
  module. From the `8b6e93d4` reruns on, the `viennarna:` line in OSTIR's `versions.yml` is the
  provenance.
- **A1 is re-run anyway, for commit parity.** Pass 1 superseded an A2 attempt for running at a
  different commit from A1 (rule A3). The lead decided on 2026-10-26 to resume A1 at `8b6e93d4` as
  well, so that both pass-3 arms carry one commit. The rerun counts as A1's second R9 attempt.
- **What the pre-fix A2 affects:** the compare and the interim metrics above read no OSTIR output.
  The stability families are islands, operon flags, enriched categories, novel ORFs and ALT sites.
  These numbers are expected to stand, and a re-compare after the rerun checks that.
- **The post-rerun job:** job 1897034 runs after both reruns. It re-applies rule set A at
  `8b6e93d4` (`evidence/a1_o390.json`, `a2_o390.json`) and re-runs the compare into
  `o390/cmp_a1_a2` (rule set C), checking the tables against the pre-fix compare. It then
  regenerates the interim `metrics.json` and `REPORT.md` into `$GIAP_ROOT/pilot/metrics/pass3/`;
  the pre-fix render moves to `pass3_pre390/`.
- **The reruns:** the A2 and A1 directories were copied to `out_a2_pre390` and `out_a1_pre390`,
  then resumed in place at the post-fix commit `8b6e93d4` (#390), so the unchanged tasks are reused
  from the cache. See "OSTIR
  reruns" below.

### Rerun results (pass 3)

The post-processing job 1897034 recorded the following.

- **OSTIR provenance:** both reruns' OSTIR `versions.yml` records `viennarna: 2.6.4`, and the
  three affected genes hold the #397 re-record values.
- **A1 (job 1897026):** 75 tasks were cached and 2 re-ran (`OSTIR`, `OSTIR_NOVEL`). Every
  published table and `rbs.parquet` is byte-identical to `out_a1_pre390`, which confirms that
  attempt 1 was already post-fix.
- **A2 (job 1897025):** 60 tasks were cached and 13 re-ran: OSTIR, the RBS and context layers, the
  gene table and every step after it. The changed outputs are `rbs.parquet`, `genes.tsv`,
  `sensitivity.tsv` and the three `determinants_*` tables. The A2 attempt-1 values of those outputs
  are pre-fix.
- **Rule sets:** A passes on both arms at `8b6e93d4` (`evidence/a1_o390.json`, `a2_o390.json`).
  The re-compare passes rule set C (`o390/evidence/compare.json`). Its detection, abundance and
  stability tables are byte-identical to the pre-fix compare.
- **Interim metrics:** `$GIAP_ROOT/pilot/metrics/pass3/` was regenerated, and the pre-fix render
  moved to `pass3_pre390/`. No metric key changed apart from paths, hashes and runtimes.

| subject | file | sha256 |
|---|---|---|
| A1 (attempt 2) | `out_a1/report/report.html` | `e173bdecd32a9dfa2e38dd53c499bec5c690aa4a017b102218d6f395f80f322e` |
| A2 (attempt 2) | `out_a2/report/report.html` | `38bc0161c22d0d7bf6e056fb58ed8f135cb7ee5cc494064dbcb380d1b20b9f96` |
| compare (attempt 2) | `o390/cmp_a1_a2/report/compare.html` | `bcb9179b1bf13a972b45b5635d2b3fc8f7c7160936521e61fb97e5a6671cbd02` |

## OSTIR reruns after #390

Every pipeline run whose OSTIR task predates #390 is re-run for `rbs` and the `tir_ostir`
determinant (#397). This is the minimal rerun:

- Each run directory is copied to `<dir>_pre390`, which keeps the accepted attempt and its evidence.
- Each run is then resumed **in place**. The out directory stays the same, so `params.config` and
  every task script are unchanged, and Nextflow reuses each task whose hash did not change.
- The runs use the pinned worktree `pilot-8b6e93d4`. Commit `8b6e93d4` is the #390 merge. Its
  `workflow/` differs from the pass-3 commit only in the OSTIR module, `nextflow.config` (the env
  scope) and `images.json`.
- OSTIR re-runs because its module changed. Its dependants re-run only where their inputs' content
  changed, because `GIAP_STEP` uses `cache 'deep'`.
- For the older commits, the changed container strings (`genomad`, `transtermhp`) and, for MVP-04,
  the `sixframe_index`/`loci` inputs also re-run.
- Each submission sources `farmshare_env.sh` from that worktree, which exports
  `PYTHONNOUSERSITE=1`, and keeps the recorded arguments (`--resume`, `--time=12:00:00`).
- At most two chains run at a time, to stay within this agent's R14 share.

| run | out directory | head job | after | original commit | result (rule set A at `8b6e93d4`) |
|---|---|---|---|---|---|
| pilot pass 3 A2 | `pass3/out_a2` | 1897025 | - | `4ccce8a5` | COMPLETED 00:55:08; pass |
| pilot pass 3 A1 (commit parity) | `pass3/out_a1` | 1897026 | job 1897025 | `4ccce8a5` | COMPLETED 00:23:48; pass; outputs unchanged |
| pilot pass 1 A1 | `pass1/out_a1` | 1897027 | job 1897026 | `6f969ee7` | COMPLETED 04:20:00 (about 3 h of it waiting for the OSTIR task's slot); pass; 63 cached, 14 re-ran |
| pilot pass 1 A2 | `pass1/out_a2` | 1897028 | job 1897027 | `6f969ee7` | COMPLETED 01:00:37; pass; 60 cached, 13 re-ran |
| pilot pass 2 A3 | `pass2/out_a3` | 1897029 | - | `242ce506` | COMPLETED 00:27:31; pass (A1-A11, `no_spectrum_pairing` expected); 63 cached, 14 re-ran |
| MVP-04 A1 | `mvp04/20261023_a1` | 1897030 | job 1897029 | `76591515` | not re-run: superseded by A2 (#326 inputs). Attempt 1897030 was cancelled at 1:06 because it resumed on A2's peptides (`3a65981d` instead of `0b6b2fd4`); the partial output is kept as `20261023_a1_o390_wronginput` |
| MVP-04 A2 | `mvp04/20261023_a2` | 1897031 | job 1897030 | `a6bf2e96` | COMPLETED 01:05:23; pass apart from A6 `cpu_cap:16` and `module_off:function`, both not applicable (below); 52 cached, 24 re-ran |

**Pass 2 A3 (job 1897029):** `evidence/a3_o390.json` has `failures: []`. Rule A11 was checked
against `pass1/out_a2_pre390`, which has the same `peptides` input as `out_a2`. OSTIR's
`versions.yml` records `viennarna: 2.6.4`. The changed outputs are `rbs.parquet`, `genes.tsv`,
`determinants_features.tsv` and `determinants_model.tsv`.

**Pass 1 A1 (job 1897027) and A2 (job 1897028):** `pass1/evidence/a1_o390.json` and
`a2_o390.json` have `failures: []`. The same holds for both arms below. OSTIR's
`versions.yml` records `viennarna: 2.6.4`, and the three affected genes hold the #397 re-record
values. The changed outputs are `rbs.parquet`, `genes.tsv`, `sensitivity.tsv` and the three
`determinants_*` tables. `detection_summary.tsv` also changed, but only by gaining the QC rows
that later code writes (`placement_flag` `pyl_site`/`sec_site`, all 0); no existing row changed.

**MVP-04 A2 (job 1897031):** the evidence is `$GIAP_ROOT/runs/mvp04/evidence/a2_o390.json`. Two
of the pilot rules do not apply to MVP-04, whose recorded command is `--threads 16 --modules
core,mobile,bgc`:

- A6 `cpu_cap:16` reflects `--threads 16`, against the pilot cap of 14.
- `module_off:function` reflects the module list, and the original run's manifest also lists
  `function` as off.

Every other rule passes. OSTIR records `viennarna: 2.6.4`, and the three K-12 genes hold the #397
values. The changed outputs are `rbs.parquet`, `genes.tsv`, `sensitivity.tsv`, the three
`determinants_*` tables and `novel_orfs.tsv`; in `novel_orfs.tsv` only the `ostir_tir` column
changed. `detection_summary.tsv` gained the newer QC rows (`placement_flag`, `placement_kind`,
`tolerant_diff_type`), and no existing row changed.

**MVP-04 A1 (job 1897030) was cancelled, and its orphaned tasks (jobs 1897137 and 1897139) too.**
The file `$GIAP_ROOT/mvp04/20261023/peptides.tsv` was regenerated by #326 (sha256 `3a65981d…`,
the MVP-04 A2 input). The original A1 input is `20261023.pre326/peptides.tsv` (sha256
`0b6b2fd4…`). The resume therefore re-mapped A1 on A2's peptides, and its rule A5 input check would
fail. The partial directory was kept as `20261023_a1_o390_wronginput`, and `20261023_a1` was
restored from `20261023_a1_pre390`. Lead decision (2026-10-26): it is not re-run, because A2
superseded A1 when #326 regenerated the inputs. The decision can be reversed: a later rerun would
pass `--peptides …/20261023.pre326/peptides.tsv` (and the matching abundance and config).

The first submission of these runs (jobs 1896986, 1897005, 1897008, 1897010, 1897011 and 1897014)
failed in 4 s with `uv: command not found`, before any run directory was touched. `GIAP_ROOT` was
not exported, so the worktree's `farmshare_env.sh` defaulted it to `<worktree>/local`. The wrapper
now exports it.

The superseded pass-1 A2 attempts (`out_a2_0667e366`, `out_a2_1549a930`) are not rerun: they
failed rule A3, and no metric reads them.

The table of jobs is `$GIAP_ROOT/runs/ostir390/JOBS.tsv`. Rule set A was re-applied to each rerun
with `--expect-commit 8b6e93d4`, and each verdict is in that pass's attempts table. Every rerun
completed with exit 0.

## Running the pilot (§15.2 steps 4-9)

Every step is a batch job; nothing runs on a login node or in the session allocation.

```bash
# step 4 - the three wild-type runs (~5-10 GB) into $GIAP_DATA_DIR/pxd078459/, with SHA256SUMS
sbatch scripts/slurm/fetch_pxd078459.sbatch
# step 5 - DIA-NN 2.3.1 for A1 and A2 (normal, 16 CPU, 64000 MB, two calls each)
sbatch -J giap-diann-A1 scripts/slurm/diann_pilot.sbatch A1
sbatch -J giap-diann-A2 scripts/slurm/diann_pilot.sbatch A2
# step 6 - A3 de novo on one GPU (Cascadia 0.0.7; `casanovo` selects the fallback)
sbatch scripts/slurm/denovo_pilot.sbatch
# step 7 - the generic peptide tables (the large reports stay under $GIAP_DATA_DIR)
giap convert --from diann    $GIAP_ROOT/pilot/A1/A1_report.parquet --out a1.tsv.gz
giap convert --from diann    $GIAP_ROOT/pilot/A2/A2_report.parquet --out a2.tsv.gz
giap convert --from cascadia $GIAP_ROOT/pilot/A3_cascadia/a3_cascadia --out a3.tsv.gz
```

```bash
# step 8 - the three GIAP runs (as batch jobs: sbatch scripts/slurm/giap_run.sbatch <args>)
giap run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a1.tsv.gz --outdir out_a1
giap run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a2.tsv.gz --outdir out_a2
giap run --genome nissle.fasta.gz --annotation nissle_bakta_full.json.gz --peptides a3.tsv.gz --reference-peptides a2.tsv.gz --outdir out_a3
giap compare out_a1 out_a2 out_a3 --labels A1,A2,A3
# step 9 - the measured runtimes that replace the estimates in the docs
sbatch scripts/slurm/benchmark_standard.sbatch
```

All three runs use the default `--genome-relation related`. Step 8 ends with
`validation/pilot/REPORT.md` (§15.5-15.7); the free-entrapment estimator of §15.5 metric 8 lives in
`validation/pilot/entrapment.py`.

The gate report is `metrics.py render --final --runbook validation/pilot/reproduce.md` (M10-15);
`review_packet.py` builds Khoa's review packet under `local/review/m10-pilot/`. One caveat of
`--final`: an absent path under a block that carries a finding reason (for example a mistyped
`m8.peptide.fdp_combinedd` under `reason: empty_discoveries`) renders `not evaluable (<reason>)`
rather than failing. Reason-only blocks need this, so check template paths by eye (review N3).

## Who does what

| Step | Who | What |
|---|---|---|
| 1-2 | agent | done in M2/M3: the recorded Bakta JSON and `db_a2/` above |
| 3 | **Khoa** | the DIA-NN licence. Khoa accepted it and supplied the 2.3.1 Academia image above, and checks the flag table of this README against the DIA-NN documentation before the first real run |
| 3b | **Khoa** | DIA-NN 2.7.0 under `$GIAP_TOOLS_DIR/diann-2.7.0`, if the pilot should use the version the spec names |
| 4-9 | agent | the batch jobs above, then `validation/pilot/REPORT.md` |

Also needed before step 5: `$GIAP_DATA_DIR/proteomes/UP000000625.fasta` (the A1 database, from
`sbatch scripts/slurm/fetch_resources.sbatch`) and `$GIAP_DB_DIR/contaminants/contaminants.fasta`
(from `giap setup --contaminants`) — the same contaminant library GIAP itself uses, so GIAP and
DIA-NN share one set.
