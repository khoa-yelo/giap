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
