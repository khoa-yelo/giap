# M10-01: PXD078459 raw data: file list from the official record, download, checksums, acquisition mode

Milestone M10 · **pass 1** (MVP) · depends on: none (main d29a2ea5) · estimate ≈ 250 code lines
(validation scripts + the sbatch change; + ≈ 300 test lines) · not `stats` · role: validation
(scripts under `validation/pilot/`, `scripts/slurm/`) · Khoa: only if the data are login-gated,
private, or the run selection is ambiguous.

## Spec sections implemented

- §15.1 "Three wild-type runs from **PXD078459** (iProX; Nissle 1917; timsTOF HT). The deposited
  Spectronaut file suggests diaPASEF; confirm the acquisition mode on the first run. The download
  is about 5–10 GB (estimate)."
- §15.2 step 4: "Download the 3 wild-type runs of PXD078459 to `$GIAP_DATA_DIR/pxd078459/` (a
  batch job; record checksums). Confirm the acquisition mode on the first run."
- §3.2 source row "iProX / ProteomeXchange (`proteomecentral.proteomexchange.org`, `www.iprox.cn`
  and the FTP host the dataset page names)"; R6 (pinned official sources, checksums, one retry,
  Blocked row, no mirrors, never TLS off); R11 (raw MS data never committed; ≤ 5 MB per file);
  R14 and §3.8 (batch job, `giap-*` name, job ID in PROGRESS); §3.9 (raw data 5–10 GB budget line).
- §15.5 metric 8 "Only wild-type runs count, because engineered strains may carry K-12-derived
  constructs": the selection rule must be explicit and recorded.
- §17 rows "PXD078459 acquisition mode differs" and "Download source down or changed".
- M3-24 decisions 6 and 11 (https/ftps only; URLs confirmed at M10), 16 (`set -f`).

This spec checks the acquisition mode of **all three** runs, not only the first: one mixed run
would change the DIA-NN settings for all.

## What exists and what is wrong with it

`scripts/slurm/fetch_pxd078459.sbatch` downloads `$PXD078459_BASE_URL/<name>` for each name in
`$PXD078459_FILES`, defaulting to `https://download.iprox.cn/dataset/PXD078459` and
`WT_1.d.zip WT_2.d.zip WT_3.d.zip`. Both are **unverified guesses** (decision 11). A HEAD request
on 2026-10-23 returned **404** for the base URL. The ProteomeXchange record
`https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD078459` answers 301 (to
`/ui?pxid=PXD078459`), and `https://www.iprox.cn/` answers 200. No file name has been seen.

The M3-24 tests (`tests/unit/test_m3_24_pilot_kit.py`, `tests/unit/test_m3_24_diann_pilot_run.py`)
run the script with no arguments against stubbed `curl`/`wget`. They expect at least three https
URLs, `SHA256SUMS` listing every created file, writes only under `$GIAP_DATA_DIR/pxd078459/`, and
no TLS-off flag. **They must keep passing unchanged.** The sandbox (`tests/helpers/pilot_kit.py`)
copies `validation/`, so a committed manifest can provide the defaults.

## Deliverables

| # | File | What |
|---|---|---|
| 1 | `validation/pilot/px_record.py` | stdlib only (`xml.etree`, `csv`, `re`). `parse_record(xml_bytes) -> list[RecordFile]` and `select_wild_type(files) -> Selection`; CLI `python -I validation/pilot/px_record.py RECORD.xml --out MANIFEST.tsv` |
| 2 | `scripts/slurm/fetch_pxd078459.sbatch` | gains `--resolve` (fetch the record and write a **draft** manifest) and reads the committed manifest by default; the env overrides stay |
| 3 | `validation/pilot/inputs/pxd078459_files.tsv` | committed manifest: one row per selected file (columns below) |
| 4 | `validation/pilot/acquisition_mode.py` | stdlib `sqlite3`, read-only: one JSON object per run from `<run>.d/analysis.tdf` |
| 5 | `validation/pilot/inputs/acquisition.json` | committed result for the three runs |
| 6 | `validation/pilot/README.md` | a "Raw data" section: accession, record URL, retrieval date, record sha256, the data licence the record states, selection rule and result, total size, acquisition mode, job IDs |

Nothing else is committed: no archive, no `.d` content, no XML record (index question L4).

### Two jobs, in order

1. **Resolve** (`sbatch --qos=dev scripts/slurm/fetch_pxd078459.sbatch --resolve`; 1 CPU, 2 GB,
   30 min). It fetches `https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD078459&outputMode=XML`
   (a small file from a §3.2 host) into `$GIAP_DATA_DIR/pxd078459/record/PXD078459.xml`, records
   its sha256 and the UTC date, runs `px_record.py`, and writes the draft manifest next to it. It
   downloads no raw data. The agent reads the draft, checks the selection, and commits it as
   deliverable 3.
2. **Download** (`sbatch scripts/slurm/fetch_pxd078459.sbatch`; 2 CPU, 8 GB, 24 h). It reads
   deliverable 3, downloads every row into `$GIAP_DATA_DIR/pxd078459/`, verifies, unpacks, and
   writes `SHA256SUMS`. Then `acquisition_mode.py` runs on the three `.d` directories (dev QOS,
   1 CPU, 4 GB, 30 min; reading SQLite metadata takes seconds).

Precedence: `PXD078459_FILES` / `PXD078459_BASE_URL` (explicit overrides, as in M3) **over** the
committed manifest. With neither, the script exits 2 with "no file list: run --resolve first".

### Manifest columns (`pxd078459_files.tsv`, tab-separated, header line, `NA` for missing)

| Column | Meaning |
|---|---|
| `run` | `WT_1`, `WT_2`, `WT_3`: GIAP's run label, in selection order |
| `file_name` | the name in the record (basename of the URI) |
| `url` | the URL the job downloads (https or ftps) |
| `url_record` | the URI exactly as the record gives it (may be `ftp://`, question Q1) |
| `bytes` | size the record states, else `NA` |
| `checksum_record` | `<algorithm>:<hex>` the record states (for example `md5:…`), else `NA` |
| `sha256` | filled by the download job's verification step (the committed draft has `NA`; the agent commits the filled file after the download) |
| `sample` | the record's sample or file-level description that the selection used |
| `selected_because` | the rule that matched, e.g. `token:WT`, `sdrf:genotype=wild type` |

Candidates that were **not** selected are listed in the README with their reason
(`not_wild_type`, `beyond_three`, `not_raw`), so the choice can be reviewed.

### Parsing the record (`parse_record`)

From the ProteomeXchange XML: every `DatasetFile` whose cvParams mark it as a raw file
(`Associated raw file URI` or the PRIDE raw-file type), with its URI, the size and checksum
cvParams when present, and any sample or description text attached to it. Search files
(the Spectronaut output), results and peak lists are type `not_raw` and never downloaded. The
instrument cvParam goes into the README; if it names no timsTOF instrument, the job stops (exit 2,
"the record does not describe timsTOF data").

### Selection (`select_wild_type`)

1. A raw file is **wild-type** when its sample annotation (SDRF genotype / characteristics) says
   wild type, or, failing an annotation, its file name contains a token `WT`, `wildtype`,
   `wild_type` or `wild-type` (case-insensitive, token-bounded: `WT` in `WT_1` matches, in
   `SWTX` it does not). A file whose name also carries a modification token (`Δ`, `delta`, `KO`,
   `mut`, `pMUT`-cured, `cured`, `::`) is not wild type.
2. Order candidates by natural sort of `file_name` (`WT_2` before `WT_10`).
3. Exactly three candidates: select them. More than three: select the first three, list the
   rest as `beyond_three` (question Q2). Fewer than three: exit 2 with the candidate list,
   "fewer than three wild-type runs"; nothing is downloaded; a PROGRESS question goes to Khoa.
4. A `.d` run stored as several files (`analysis.tdf`, `analysis.tdf_bin`, under a folder path in
   the URI) is one run: its files share a `run` label and keep their relative paths.

### Download rules

- Scheme: `https://` or `ftps://` only. An `ftp://` or `http://` URL in the manifest is exit 2
  ("R6 / M3-24 decision 11: plain FTP or HTTP is refused"), naming the file. Never
  `--insecure`/`-k`/`--no-check-certificate`.
- Disk: before the first byte, free space under `$GIAP_DATA_DIR` must be ≥ 3 × the sum of
  `bytes` (archive + unpacked + margin); else exit 2. With `bytes` `NA`, assume 5 GB per run.
- Resume: `curl --continue-at -` on a partial file; a file whose size equals `bytes` and whose
  record checksum matches is skipped.
- Verify: when `checksum_record` is present it must match; a mismatch deletes the file and
  retries once; a second mismatch is exit 1, no `SHA256SUMS` line for that file, and a Blocked
  row. `sha256` is computed for every file into `SHA256SUMS` (as M3 does) and into the manifest.
- Unpack: `.zip` (`unzip -t` first), `.tar`, `.tar.gz`. An archive member with an absolute path
  or a `..` component is refused (exit 1, "unsafe archive member"): the archives are untrusted.
  The run ends up at `$GIAP_DATA_DIR/pxd078459/<run>.d/`, whatever folder nesting the archive
  had; exactly one `analysis.tdf` per run. Any other format (`.7z`, `.rar`) is exit 2 with a
  question for the lead. The archive is deleted after a successful unpack (§3.9).
- Exit codes: 0 done; 1 download or verification failed after the retry; 2 a precondition
  (scheme, disk, manifest, selection, format).

### Acquisition mode (`acquisition_mode.py`)

Opens `analysis.tdf` with `sqlite3.connect("file:…?mode=ro", uri=True)` and reads:

- `Frames.MsMsType` counts. In TDF, 0 = MS1, 8 = ddaPASEF MS/MS, 9 = diaPASEF MS/MS (as
  documented by the open-source TDF readers, for example AlphaTims; the test-author takes the
  values from that documentation, not by importing the reader);
- whether `DiaFrameMsMsWindows` and `PasefFrameMsMsInfo` exist and have rows;
- the window scheme (number of window groups, m/z and 1/K0 ranges) from `DiaFrameMsMsWindows`;
- `GlobalMetadata` keys `InstrumentName`, `InstrumentFamily`, `AcquisitionSoftwareVersion`,
  `AcquisitionDateTime` when present; the gradient length = the maximum `Frames.Time`.

`mode` is `diaPASEF` (type 9 present, no type 8, windows table non-empty), `ddaPASEF` (type 8
present, no type 9, PASEF table non-empty), `ms1_only` (only type 0), `mixed` (8 and 9), or
`unknown` (anything else). The file gets a top-level `dataset_mode`: the common mode when all
runs agree, else `inconsistent`. The script exits 0 only when `dataset_mode = diaPASEF`; for
`ddaPASEF` it exits 3 ("settings must change: §17; stop before M10-02"); for anything else, 1.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | Record with six raw files `WT_1.d.zip`, `WT_2.d.zip`, `WT_3.d.zip`, `dLon_1.d.zip`, `dLon_2.d.zip`, `dLon_3.d.zip`, https URIs, sizes and MD5s | manifest rows `WT_1`–`WT_3`, `selected_because = token:WT`, `checksum_record = md5:…`; three `dLon_*` listed as `not_wild_type`; the Spectronaut `.sne` as `not_raw` |
| 2 | Record with an SDRF giving `genotype = wild type` for `S01.d.zip`–`S03.d.zip` and `ΔfliC` for `S04`–`S06` | `S01`–`S03` selected, `selected_because = sdrf:genotype=wild type` (annotation wins over names) |
| 3 | Four wild-type files `WT_1`, `WT_10`, `WT_2`, `WT_3` | natural sort → `WT_1`, `WT_2`, `WT_3` selected; `WT_10` listed `beyond_three` |
| 4 | Two wild-type files only | resolve exits 2, "fewer than three wild-type runs: WT_1.d.zip, WT_2.d.zip"; no manifest written |
| 5 | Record URIs `ftp://download.iprox.org/IPX…/WT_1.d.zip` only, no override | download exits 2 naming the file and decision 11; nothing downloaded (question Q1 decides a same-host https route) |
| 6 | `WT_2.d.zip` MD5 differs twice | exit 1; `SHA256SUMS` has `WT_1`, `WT_3` only; log says "checksum mismatch after one retry: WT_2.d.zip"; the agent adds a Blocked row |
| 7 | Job killed at the 24 h limit with `WT_3.d.zip` half written; re-submitted | `WT_1`, `WT_2` skipped ("already there"), `WT_3` resumed with `--continue-at -`, then verified |
| 8 | Archive `WT_1.d.zip` holding `../../evil.sh` | exit 1, "unsafe archive member: ../../evil.sh"; nothing extracted |
| 9 | `analysis.tdf` with `MsMsType` {0: 1,200, 9: 10,800} and 24 rows in `DiaFrameMsMsWindows` for all three runs | `mode = diaPASEF` per run, `dataset_mode = diaPASEF`, exit 0 |
| 10 | Runs 1–2 diaPASEF, run 3 {0, 8} with `PasefFrameMsMsInfo` rows | run 3 `ddaPASEF`, `dataset_mode = inconsistent`, exit 1 |

## Edge cases

- **Login-gated or private data:** the record lists no file URIs, or they answer 401/403 or a
  login page (HTML where an archive was expected: the content type or the zip test fails). The
  job exits 1; the agent adds a Blocked row "PXD078459 raw files not publicly downloadable" and
  asks Khoa. Never another copy (R6).
- **Record unreachable:** one retry later (R6), then Blocked.
- **`bytes` and `checksum_record` both missing:** download, compute sha256, and say in the README
  that the record published no checksum ("integrity: TLS only").
- **A `.d` already unpacked from an earlier run but its archive deleted:** skipped when its
  `analysis.tdf` and `analysis.tdf_bin` exist and the manifest `sha256` of the archive is recorded;
  otherwise downloaded again.
- **Nested or renamed folders in the archive** (`WT_1/20260101_EcN_WT_1.d/`): normalised to
  `<run>.d`; the original folder name goes into `acquisition.json` (`source_folder`).
- **No `analysis.tdf` in a run** (the archive held `.raw` or mzML): `mode = unknown`,
  "not a Bruker TDF run", exit 1.
- **`analysis.tdf_bin` missing:** acquisition mode still reads (metadata only) but the run is
  flagged `incomplete_run` and the script exits 1: DIA-NN needs the binary.
- **diaPASEF variants** (synchro-PASEF, slice-PASEF, midia-PASEF) still write type 9 frames:
  `mode = diaPASEF` plus `notes` naming the window scheme (question Q3).
- **Empty `Frames` table:** `unknown`, exit 1.
- Circular/linear genomes, strands and ties of coordinates do not apply to this issue; the
  ordering tie (equal names after natural sort cannot happen: names are unique in a record) is
  covered by example 3.

## Acceptance criteria

Unit and stub tests (CI, every commit):

1. `parse_record` on synthetic XML fixtures (written by the test-author under
   `tests/data/pilot/px_record_*.xml`, a few kB each) returns the files, types, sizes and checksums
   of examples 1–3; the Spectronaut file is `not_raw`.
2. `select_wild_type` gives exactly the selections of examples 1–4, including natural sort and
   the modification-token exclusion; fewer than three raises with the candidate list.
3. `fetch_pxd078459.sbatch` against the stubs: with the committed manifest and no overrides it
   requests exactly the manifest URLs; an `ftp://` or `http://` row is exit 2 before any request;
   a record checksum mismatch retries once then exits 1 without a `SHA256SUMS` line for that file;
   a partial file is resumed (`--continue-at -` in the stubbed call); an unsafe archive member is
   exit 1; too little disk (stubbed `df`) is exit 2.
4. The M3-24 tests in `tests/unit/test_m3_24_*.py` pass **unchanged**.
5. `acquisition_mode.py` on SQLite fixtures built at test time (tables `Frames`,
   `DiaFrameMsMsWindows`, `PasefFrameMsMsInfo`, `GlobalMetadata`) gives the modes and exit codes
   of examples 9–10 and of `ms1_only`, `mixed`, `unknown`, missing `analysis.tdf`, missing
   `analysis.tdf_bin`; it never writes to the database (the fixture's mtime and sha256 are unchanged).
6. No file under `validation/` exceeds 5 MB; `git ls-files` lists no `.d`, `.zip`, `.tdf`,
   `.tdf_bin`, `.mzML` or `.xml` record under `validation/pilot/`.

Evidence on real data (the M10-01 PR description and PROGRESS, not CI):

7. The resolve and download job IDs, their `sacct` elapsed time and MaxRSS, and the total size are
   in PROGRESS (In progress, then Measured numbers).
8. `validation/pilot/inputs/pxd078459_files.tsv` has three runs with `sha256` filled, and each
   matches `SHA256SUMS` on disk.
9. `validation/pilot/inputs/acquisition.json` has three runs and `dataset_mode = diaPASEF`, or
   the issue stops with the §17 path (exit 3, PROGRESS question, M10-02 waits).
10. The README "Raw data" section states the record URL, retrieval date, record sha256, data
    licence, selection, sizes and acquisition mode.

Link to §15.6: every pass criterion needs these runs; §15.1's "confirm the acquisition mode"
and §15.5 metric 8's "only wild-type runs" are satisfied here.

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| Q1 | The record may give only `ftp://` URIs. May the job use `https://` on the **same host and path** when a HEAD returns 200 and its `Content-Length` equals the record's `bytes`? It is the official host, not a mirror. | yes; `url` and `url_record` both recorded; without a matching HEAD, exit 2 and Blocked |
| Q2 | More than three wild-type runs: the first three by natural sort, or by the record's replicate numbering if it has one? | natural sort of file names; the rest listed; the choice is reviewable and reversible (only the download changes) |
| Q3 | Count diaPASEF variants (synchro-, slice-, midia-PASEF) as `diaPASEF`? DIA-NN 2.3.1 supports some of them. | yes, with `notes`; M10-02 checks the window scheme against the DIA-NN documentation |
| Q4 | Check all three runs or only the first (§15.2 says "the first run")? | all three; it costs seconds and catches a mixed dataset |
| Q5 | Should the resolve job's draft manifest be committed before the download, then again with `sha256` filled, or once at the end? | twice, two commits in the same PR: the reviewed selection first, the hashes after the download |

## Lead decisions (2026-10-23, after the tests at 0041b05e)
1. **The M3-24 fetch tests vs the committed manifest:** at implementation, the lead amends the M3-24 tests (`test_fetch_downloads_three_runs…` and the others that run the script bare against stub payloads) and their sandbox. They then pass an explicit sandbox manifest instead of relying on the defaults. The script gets no stub-detection path.
2. **The committed-data tests (acceptance 8–10):** the implementation PR is opened only after the real `--resolve` job has run and its manifest and README section are committed, so those tests are green on CI. The download stage's hashes are a second commit in the same PR, before it is marked ready.
3. **The test-author's readings are accepted:**
   - escaping symlinks are refused as an "unsafe archive member";
   - the zip-bomb limit gives exit 1 or 2;
   - the filled manifest is written in the data directory;
   - the HEAD upgrade may happen at resolve or at download.

   The 1/K0 range and the Q3 variant notes are backlogged.
