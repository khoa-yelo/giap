# M10-03: `giap convert --from diann` on the pilot A1/A2 reports (DIA-NN 2.3.1)

Milestone M10 · **pass 1** (MVP) · depends on: M10-02 (both reports under
`$GIAP_ROOT/pilot/<arm>/`, `check.json` `ok`, the stratified samples under
`$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/`) · estimate ≈ 170 code lines (index: 150) +
≈ 380 test lines · not `stats` · roles: test-author (new contract test, CI fixtures), implementer
(`src/giap/io/converters/diann.py` fixes, `scripts/slurm/diann_convert.sbatch`), validation
(`validation/pilot/diann_convert.py`, README section), lead (retires the old contract test) ·
Khoa: one question (Q-K1: are the converted tables "derived from the DIA-NN reports"?); the
default below does not wait for him.

## Spec sections implemented

- §5.4 converters: "Each converter is a pure function from input files to the generic table, with
  a contract test against a recorded sample file (§11.4). Until a real sample exists, the
  converter is written from the tool's documentation and its contract test is marked
  `pending_sample`; the pilot (§15) supplies DIA-NN and Cascadia samples." The DIA-NN 2.x row of
  the mapping table is unchanged; this issue confirms it against 2.3.1 output.
- §5.3 generic peptide table (columns, types, validation) and the abundance table
  (`--abundance-type maxlfq`).
- §5.1 `giap convert --from diann INPUT --out peptides.tsv [--abundance-out abundance.tsv]
  [--search-db same|other]`.
- §15.2 step 7, the DIA-NN half: "Run `giap convert --from diann …` on FarmShare; the large
  reports stay in `local/data/`. Commit the gzipped generic tables (a few MB) to
  `validation/pilot/inputs/` …". **The commit half is overridden** by Khoa's 2026-10-23 rule
  (M10-02 lead decisions 4 and 6: nothing derived from the DIA-NN reports is committed), read
  here as covering the converted tables (default; Q-K1).
- §11.4 contract tests; R2 (tests first), R5 (no DIA-NN code imported: the converter reads its
  output files only), R11 (≤ 5 MB per file; nothing derived from the reports committed), R14
  (`giap-*` jobs, CPU cap, large files under `$GIAP_ROOT`).
- Index row M10-03, index decision 5 (superseded, see "What is committed"), M10-02 edge case
  "M10-03 never converts a report whose `check.json` is not `ok`".

## What exists (main dfe3fe19, `m10/index`, `origin/m10/diann-searches` 6f1f1d6e)

| Piece | State |
|---|---|
| `src/giap/io/converters/diann.py` (M1-07) | written from the DIA-NN 2.x documentation. Loads only the mapped columns. Required: `Stripped.Sequence`, `Modified.Sequence`, `Protein.Ids`, `Global.Q.Value` or `Q.Value`, `Run`, `Precursor.Charge`. Optional: `PEP`, `Decoy`, `Precursor.Normalised`/`Precursor.Quantity`, `Protein.Group` + `PG.MaxLFQ` (abundance). Parquet numbers are widened to float64; `Decoy` must be 0/1; intensity 0 → NA; the first present column of each fallback pair wins **for the whole file** (no per-row mixing) |
| `CONVERTERS` (`src/giap/io/converters/__init__.py`) | `diann`, `generic`. Unchanged here |
| `giap convert` (`src/giap/cli.py`) | exit 0, or 1 (`INVALID_INPUT`) on any error issue, with no output written; prints every issue to stderr; `.gz` outputs have a fixed gzip header (`mtime=0`), so the bytes are deterministic. Floats are written `%.6g` (§6.6, ADR-0006) |
| `tests/unit/test_convert_diann.py`, `test_cli_convert.py` | hand-written TSV and parquet fixtures in the 2.x shape; all float64/int64 columns, no extra 2.x columns, no `Cont_` IDs |
| `tests/contract/test_diann_contract.py` | `pending_sample`; looks for `tests/data/recorded/diann/report.{parquet,tsv}` (the top level) and optional byte-golden `expected_*.tsv` next to it. **Both are now forbidden** in the repository (Khoa's rule), so as written it can never run |
| M10-02 outputs (not yet produced; `$GIAP_ROOT/pilot/` does not exist on 2026-10-23) | `$GIAP_ROOT/pilot/<arm>/<arm>_report.parquet`, `outputs.sha256`, `check.json`, `provenance.json`; the sample `$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/report.sample.parquet` (≤ 1,000 rows, stratified by run and decoy: 900 targets + 100 decoys), `schema.tsv`, `summary.json` |
| Search databases | A1: `$GIAP_DATA_DIR/proteomes/UP000000625.fasta` (`>sp|A5A616|MGTS_ECOLI …`). A2: `validation/pilot/db_a2/strain_db.faa.gz` (`>JMJJOG_00001 chromosomal replication initiator protein DnaA`). Both arms: `$GIAP_DB_DIR/contaminants/contaminants.fasta`, 381 records, every one `>sp|Cont_<ACC>|…` (358) or `>tr|Cont_<ACC>|…` (23) |
| Contaminant tagging | not the converter's job: preflight sets `is_contaminant` when **every** `protein_ids` item starts with an `input.contaminant_prefixes` entry (`Cont_`, `CON__`, `CONT_`, `cRAP`). `is_decoy` comes from the explicit column first and only falls back to `input.decoy_prefixes` where it is NA |
| Abundance use downstream | `genetable.abundance.gene_abundance` uses peptide intensities whenever accepted peptides carry them, and the protein-level table (C10) **only** when they do not. DIA-NN rows carry `Precursor.Normalised`, so a MaxLFQ table passed to `giap run` would be read and validated but not used for gene abundance (Q4) |

## What changes for DIA-NN 2.3.1

The mapping of §5.4 stays. What changes is (a) the contract test, which now reads the real sample
on FarmShare, (b) a set of converter adjustments that are **conditional**: each is made only if
the real sample shows the condition, and each has a CI test on a hand-made fixture, and (c) a
Slurm job that converts the full reports.

### Facts to confirm on the sample (the contract test asserts each)

| # | Area | Expected for 2.3.1 (DIA-NN 2.x documentation) | What the contract test asserts | Converter change if the sample differs |
|---|---|---|---|---|
| F1 | Columns | every column of the converter's `_MAP` plus `Protein.Group`, `PG.MaxLFQ`, present under these exact names; many extra 2.x columns (`Precursor.Id`, `Lib.Q.Value`, `Peptidoform.Q.Value`, `PG.Q.Value`, `RT`, `IM`, …) that the converter ignores | the first-choice column of every fallback pair is present: `Global.Q.Value` (not `Q.Value`), `Precursor.Normalised` (not `Precursor.Quantity`); `PEP`, `Decoy`, `Protein.Group`, `PG.MaxLFQ` present | a renamed column gets one more name at the **end** of its fallback tuple (the §5.4 name stays first), with a CI fixture using only the new name. A missing first-choice column is a finding for the lead, not a silent fallback |
| F2 | Types | parquet; scores and intensities may be `float` (float32); `Precursor.Charge` an integer type; `Decoy` an integer (0/1) or `bool`; text columns `string`/`large_string`, possibly dictionary-encoded (pandas `category`) | conversion has no error and no warning issue; `charge` is integral; `is_decoy` is boolean | float32: none (widening is exact; output is `%.6g`). `bool` `Decoy`: the `_decoy` path must accept `True/False` as 1/0 (CI fixture). Dictionary text: cast to string before `StringDtype` if pandas raises (CI fixture) |
| F3 | q-values | `Global.Q.Value` populated on every target row, in [0, 1]. On decoy rows it may be populated or NA | `q_value` equals `float64(Global.Q.Value)` row by row; targets never NA; the decoy NA count is recorded in the test output | none. **No per-row fallback** to `Q.Value`: one column of one definition (Q2). Decoys without q are M4-24's "partly scored decoys" (they carry `PEP`) |
| F4 | Decoys | `--report-decoys` gives rows with `Decoy = 1`; their `Stripped.Sequence` is the decoy sequence; their `Protein.Ids` are the IDs of the target they were derived from | `is_decoy` count equals the `Decoy = 1` count (100 in the sample); every decoy row has `protein_ids` (else `decoy_without_protein_ids` would be warned, and gene-level pairing, ADR-0009, fails) | none if the IDs are there. Empty decoy IDs are a finding for the lead (the picked tier needs them; ADR-0047), not something the converter can invent |
| F5 | Modifications | `Modified.Sequence` in `(UniMod:n)` notation; only `UniMod:4` (fixed C, `--unimod4`) and `UniMod:35` (M, `--var-mod`), no N-terminal acetyl (no `--met-excision`, no `UniMod:1` variable mod); `Stripped.Sequence` uppercase A–Z | `strip_modifications(modified_peptide) == peptide` on every row; the set of `UniMod:n` tokens ⊆ {4, 35}; `peptide` is 7–30 letters | a token outside {4, 35} is a finding on the cfgs (M10-02 / ADR-0051), not the converter |
| F6 | Run names | `Run` is the run stem (`WT_1`), not a path; one value per raw run, matching `acquisition.json` by M10-02 check rule 3 | `run` values = `Run` values verbatim; the set equals the three acquisition labels (or their `.d` stems) | **none: `run` stays verbatim** even if it is a path (Q3). Only if Q3 is answered otherwise: basename, then strip `.d` |
| F7 | Protein IDs | `Protein.Ids` `;`-separated; UniProt headers reduced to the accession (`P0A7V0`, `Cont_P00761`); Bakta headers reduced to the first token (`JMJJOG_01234`), not the description | every item of every `protein_ids` value, and of every abundance `protein_id`, is an ID of the arm's search FASTAs, as `giap.io.protein_fasta.resolve_protein_id` reads their first header tokens (A1: UP000000625 + contaminants; A2: `strain_db.faa.gz` + contaminants) | if DIA-NN keeps a description or the full `sp|ACC|NAME`, each `;` item is reduced to its first whitespace token and then through `resolve_protein_id` (CI fixture). Any other mismatch is a finding |
| F8 | Contaminant tag | the `Cont_` prefix survives in `Protein.Ids` (it sits in the accession field of `sp|Cont_ACC|NAME`) | at least one contaminant ID is checked through F7. Whether any contaminant row exists is recorded, not asserted (the sample may hold none); the full-report count is in `convert.json` | none. If DIA-NN drops the prefix (F7 fails with `P00761` instead of `Cont_P00761`), the converter cannot recover it: the contaminant FASTA headers must change and DIA-NN must re-run (finding for the lead; M10-02 re-run) |
| F9 | Abundance | one `PG.MaxLFQ` value per (`Protein.Group`, `Run`) repeated on each precursor row; 0 or NA for groups not quantified in a run | the abundance table has one row per non-decoy `Protein.Group` with a positive value in some run, its value the median over runs of the per-run value (independent oracle in the test), no duplicate `protein_id` | none expected |
| F10 | Row count | one output row per report row; nothing dropped or added by the converter (filters are preflight's) | `len(peptides) == sample rows`; the multiset of `(run, peptide, modified_peptide, charge, is_decoy)` equals the sample's `(Run, Stripped.Sequence, Modified.Sequence, Precursor.Charge, Decoy)` | none |

Each converter change ships with a CI test on a hand-made parquet fixture written from the DIA-NN
documentation and the condition the sample showed. **The fixture never copies rows, values or the
column list of the real report** (Khoa's rule; Q-K2): it names only the columns the change is
about, plus the converter's required ones.

### The contract test (`tests/contract/test_m10_03_diann_contract.py`, test-author)

- **Data.** For `arm` in `A1`, `A2`:
  `SAMPLE = $GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/report.sample.parquet`, and the arm's
  FASTAs (A1: `$GIAP_DATA_DIR/proteomes/UP000000625.fasta`; A2: the committed
  `validation/pilot/db_a2/strain_db.faa.gz`; both: `$GIAP_DB_DIR/contaminants/contaminants.fasta`)
  and `validation/pilot/inputs/acquisition.json` (M10-01).
- **Skip, never fail, when the data are absent.** `GIAP_ROOT` unset (as on CI) → skip with
  "`GIAP_ROOT` is not set: the DIA-NN samples live on FarmShare only (Khoa 2026-10-23)". Sample
  or FASTA missing → skip naming the path and "M10-02 report job not run". This is the pattern of
  `tests/contract/test_genomad_contract.py` and M10-02's `_sample_dir`. The test asserts that the
  sample path is **outside** the checkout.
- **No `pending_sample` marker** (the sample exists on FarmShare); no `slow` marker (1,000 rows
  convert in well under a second).
- **Assertions:** F1–F10 above, each its own test function, parametrized by arm, plus:
  - the CLI on the sample (`giap convert --from diann SAMPLE --out t/a.tsv.gz --abundance-out
    t/ab.tsv.gz --search-db other|same`) exits 0, and a second call gives byte-identical files;
  - the written table round-trips: `read_peptide_table` on it gives no error issue and the same
    row count; `read_abundance` on the abundance file gives no error issue;
  - no byte-golden file is read or written anywhere (golden outputs of the converter would be
    derived data and would only test the converter against itself).
- **One optional FarmShare test, marked `slow`:** the same F1, F2, F6, F7 and F10 checks on the
  full `$GIAP_ROOT/pilot/<arm>/<arm>_report.parquet`, skipped when absent. It runs in the
  conversion job (below), not in the fast suite.

The old `tests/contract/test_diann_contract.py` is **retired by the lead** in the implementation
step (role guard: the test-author does not edit existing tests; Q1). Kept, it would skip forever
on a path that may never exist.

### The conversion job (`scripts/slurm/diann_convert.sbatch <arm>`, implementer)

| `-J` | Partition / QOS | CPUs | Mem | Time limit (expected) | Count |
|---|---|---|---|---|---|
| `giap-diann-convert-A1`, `…-A2` | normal / normal | 2 | 8 GB | 00:30:00 (< 2 min) | 2, may run together: 4 CPUs |

Sizing. The converter is single-threaded pandas; pyarrow's parquet read uses the second CPU.
`MaxMemPerCPU` is 4 GB, so 2 CPUs give 8 GB. The reports are unmeasured. A bacterial
library-free search on 3 runs is of the order of 10⁵ precursor rows (M10-02 records the true
number in `summary.json`); M1-07 measured 500,000 rows in < 15 s. About 12 of the columns are
loaded, roughly 0.5 KB per row, so even 2 × 10⁶ rows stay under 2 GB. On `OUT_OF_MEMORY`, resubmit
once at 4 CPUs / 16 GB. The `normal` QOS is used, not `dev` (`dev` allows one job and M10-01/08
may hold it). Before `sbatch` the lead sums the running and pending `giap-*` CPUs (R14). The job
adds 4 at most.

Submission, from the worktree, after M10-02's report jobs (`--dependency=afterok:<id>` when they
are still queued):

```bash
export GIAP_REPO=/home/classes/bios/270/khoa/giap/local/worktrees/m10-03
sbatch -J giap-diann-convert-A1 "$GIAP_REPO/scripts/slurm/diann_convert.sbatch" A1
sbatch -J giap-diann-convert-A2 "$GIAP_REPO/scripts/slurm/diann_convert.sbatch" A2
```

The script, in order (header as `scripts/slurm/_header.sbatch`, `set -euo pipefail`, `set -f`):

1. **Arguments.** Exactly one of `A1`, `A2`, else exit 2. `SEARCH_DB` is `other` for A1 (K-12
   database, Nissle genome) and `same` for A2 (Q5). `PILOT="${OUT_DIR:-$GIAP_ROOT/pilot}"`.
2. **Preconditions** (each failure exits 2, names the path, and writes nothing):
   - `REPORT=$PILOT/<arm>/<arm>_report.parquet` exists and is non-empty;
   - `$PILOT/<arm>/check.json` exists and has `"ok": true`;
   - the sha256 of `REPORT` equals its line in `$PILOT/<arm>/outputs.sha256` (it is the file the
     search finished and the check passed);
   - `$PILOT/<arm>/FAILED` does not exist.
3. **Code under test.** `GIAP_PYTHON="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"`;
   giap runs as `PYTHONPATH="$GIAP_REPO/src" "$GIAP_PYTHON" -c 'from giap.cli import cli; cli()'
   …`, so the branch's converter runs, not the editable install of the main checkout. No
   `uv run`, no install.
4. **Convert.** `OUT=$PILOT/inputs`, created `0755`. Outputs go to temporary names
   (`.<name>.tmp.$SLURM_JOB_ID`) and are renamed only after success:
   - `<a>.tsv.gz` (`a1`/`a2`), the generic peptide table;
   - `<a>_abundance.tsv.gz`, the MaxLFQ table (`protein_id`, `maxlfq`);
   - `<a>.convert.log`, the converter's stderr (every issue line).

   Command: `giap convert --from diann "$REPORT" --out … --abundance-out … --search-db
   "$SEARCH_DB"`. A non-zero exit is exit 1: the temporary files are deleted and
   `$OUT/<a>.FAILED` is written (exit code, UTC time).
5. **Determinism.** The same command is run a second time into `$TMPDIR`; both outputs must be
   byte-identical to the first (`cmp`), else exit 1.
6. **Summary.** `"$GIAP_PYTHON" -I validation/pilot/diann_convert.py summary --arm <arm>
   --report "$REPORT" --peptides … --abundance … --log … --out $OUT/<a>.convert.json`
   (schema below). Exit 1 if the summary finds an inconsistency (its rules below).
7. **Full-report contract checks.** `"$GIAP_PYTHON" -m pytest -m slow -p no:cacheprovider
   tests/contract/test_m10_03_diann_contract.py -k full_report`, run with `PYTHONPATH` set as in
   step 3. A failure is exit 1. The pytest output goes to the job log.
8. Outputs are moved into place only after steps 5–7 pass, so the paths M10-04 reads hold only
   checked tables. A re-run overwrites them: same report, same bytes.

Exit codes: 0 converted and checked; 1 conversion, determinism, summary or contract failure; 2 a
precondition. The script never writes under `$GIAP_REPO`.

### `validation/pilot/diann_convert.py summary` (validation)

A pure function `summarize(report_cols, peptides, abundance, issues) -> dict` plus a thin CLI.
pyarrow and the giap readers; no network. `<a>.convert.json` holds:

- `arm`, `search_db`;
- `report`: path, sha256, bytes, rows;
- `outputs`: for the peptide table, the abundance table and the log: path, sha256, bytes, data
  rows;
- `giap`: `giap version --json`, `git -C $GIAP_REPO rev-parse HEAD`, dirty flag;
- `command`: the argv;
- `job`: id, node, start and end UTC;
- `issues`: code → {severity, rows}, from the converter's stderr;
- `counts`, all over the **converted table**:
  - `rows`, `target_rows`, `decoy_rows`, `decoy_rows_q_na`;
  - `contaminant_only_rows` (every `protein_ids` item starts with a configured contaminant
    prefix), `mixed_contaminant_rows` (some items but not all);
  - `rows_per_run`;
  - `distinct_peptides_target_q01` (distinct `peptide` of target rows with `q_value ≤ 0.01`, read
    back from the TSV);
  - `q_crossings`: rows whose report `Global.Q.Value` is > 0.01 but whose written `q_value` is
    ≤ 0.01 after `%.6g` (see "Edge cases");
  - `intensity_na_rows`;
  - `abundance_groups`.

Summary rules (each failure is a string in `failures`; any failure is exit 1):

1. The data rows of the peptide table equal the report rows (F10 on the full report).
2. `decoy_rows` equals the report's `Decoy = 1` count, and is > 0.
3. No issue of severity `error` (the CLI would have exited 1); any `warn` code is listed, not
   failed.
4. `rows_per_run` has the same keys as the report's distinct `Run` values.

`<a>.convert.json` holds paths, hashes and counts only, no row (Q-K1).

### What is committed (and what is not)

| Item | Where | Committed |
|---|---|---|
| `a1.tsv.gz`, `a2.tsv.gz`, `a{1,2}_abundance.tsv.gz`, `*.convert.log`, `*.convert.json` | `$GIAP_ROOT/pilot/inputs/` | **no** (default; Q-K1) |
| contract test, CI fixtures (hand-made), converter fixes, the sbatch, `diann_convert.py` | repository | yes |
| `validation/pilot/README.md` › "Conversion (M10-03)": the two commands, job IDs, `sacct` State/Elapsed/MaxRSS, and per arm the report sha256, the two output sha256s and their row counts | repository | yes, if Q-K1 allows counts and hashes (as M10-02's `PROVENANCE.md` does for the report); else only commands and job IDs |

**Default reading of Khoa's rule: the converted tables are derived from the DIA-NN reports**
(every row is a report row, renamed and typed), so they live under `$GIAP_ROOT` and are not
committed. Consequences, each a change to an earlier plan:

- index decision 5 and §15.2 step 7 ("commit the gzipped generic tables to
  `validation/pilot/inputs/`") no longer hold for DIA-NN. `validation/pilot/inputs/` keeps only
  M10-01's `pxd078459_files.tsv` and `acquisition.json`. The deviation from §15.2 step 7 is
  recorded as an amendment to ADR-0051 (the lead adds one paragraph; Q6);
- M10-04 reads `--peptides $GIAP_ROOT/pilot/inputs/a1.tsv.gz` (and `a2`). M10-14's
  `reproduce.md` regenerates the tables with this job and checks them against the sha256s in the
  README. Conversion is deterministic, so the hash is a full check;
- M10-14 deletes the raw reports (§3.9); the converted tables are then the only copy. They stay
  under `$GIAP_ROOT/pilot/inputs/` until Khoa says otherwise. They are a few MB.

A guard test (criterion 9) keeps this true.

## Examples

Rows are hand-made in the 2.3.1 shape (not taken from the real report). Output columns in §5.3
order; `score` and `spectrum_id` are absent because DIA-NN gives neither.

| # | Input | Expected |
|---|---|---|
| 1 | A1 report row: `Run=WT_1`, `Stripped.Sequence=PEPTMIDEK`, `Modified.Sequence=PEPTM(UniMod:35)IDEK`, `Protein.Ids=P0A7V0`, `Global.Q.Value=0.0021` (float32), `Q.Value=0.0009`, `PEP=0.003`, `Decoy=0`, `Precursor.Charge=2`, `Precursor.Normalised=15234.5`, `Precursor.Quantity=14000`; `--search-db other` | `PEPTMIDEK  PEPTM(UniMod:35)IDEK  P0A7V0  0.0021  0.003  false  15234.5  WT_1  2  database  other` (q from `Global.Q.Value`, not `Q.Value`; intensity from `Precursor.Normalised`) |
| 2 | Example 1 plus its decoy: `Stripped.Sequence=PLPTMIDDK`, `Decoy=1`, `Protein.Ids=P0A7V0`, `Global.Q.Value=0.4`, `PEP=0.97` | a second row with `is_decoy=true`, `protein_ids=P0A7V0`, `q_value=0.4`. Both rows are kept: q filtering is preflight's (ADR-0047 then drops the decoy before picking, q > 0.01) |
| 3 | Decoy row with `Global.Q.Value` NA (`Q.Value=0.3` present) | `q_value=NA`, `pep` kept; no fallback to `Q.Value` for that row; no issue from the converter. `convert.json` `decoy_rows_q_na=1` |
| 4 | Contaminant row: `Protein.Ids=Cont_P00761` (porcine trypsin) | `protein_ids=Cont_P00761` unchanged; at `giap run` preflight sets `is_contaminant=true`; `convert.json` `contaminant_only_rows` counts it |
| 5 | A1 row shared by the lacZ contaminant and its K-12 entry: `Protein.Ids=Cont_P00722;P00722` | `protein_ids=Cont_P00722;P00722`; not contaminant-only (not every item is prefixed), so it is target evidence; `mixed_contaminant_rows` counts it |
| 6 | A2 row: `Protein.Ids=JMJJOG_01234`, `Modified.Sequence=AC(UniMod:4)DLM(UniMod:35)TGK`, `Stripped.Sequence=ACDLMTGK`, `Precursor.Normalised=0`; `--search-db same` | `peptide=ACDLMTGK`, `intensity=NA` (0 → missing), `search_db=same` |
| 7 | Abundance: rows of `Protein.Group=JMJJOG_01234` with `PG.MaxLFQ` 100 (`WT_1`, three precursor rows), 300 (`WT_2`), 0 (`WT_3`); a decoy row of the same group with `PG.MaxLFQ=900` | one abundance row `JMJJOG_01234  200` (per run one value, 0 excluded, decoys excluded, median of 100 and 300) |
| 8 | `diann_convert.sbatch A1` with `check.json` `"ok": false` | exit 2, message names `check.json`; no file under `$GIAP_ROOT/pilot/inputs/`; no `giap convert` call |
| 9 | `diann_convert.sbatch A2` where the report was rewritten after the search (sha256 ≠ `outputs.sha256`) | exit 2 naming both hashes; nothing written |
| 10 | `diann_convert.sbatch A1` twice on the same report | both runs exit 0; `a1.tsv.gz` and `a1_abundance.tsv.gz` have the same sha256 both times (fixed gzip header, stable sort); `a1.convert.json` differs only in `job` |

## Edge cases

- **Empty input.** A report with 0 rows cannot reach this job: M10-02 `check` rule 1 fails it,
  and step 2 here exits 2. Called directly, the converter gives `empty_table` (error) and the CLI
  exits 1 with no output. A report with targets but **no decoy** fails M10-02 check rule 4 and so
  stops at step 2. A report of decoys only would convert, but summary rule 2 passes and M10-04's
  run would find no target evidence. This cannot happen at `--qvalue 0.01` with real targets, so
  it is not special-cased.
- **Ties.** Output rows are sorted stably (mergesort) by `run`, `peptide`, `modified_peptide`,
  `charge`. Rows equal on all four (the same precursor twice in one run, which a well-formed
  report does not have; or a target and a decoy with the same stripped sequence, modification
  and charge) keep their report order. The output therefore depends only on the report's bytes.
  A decoy whose stripped sequence equals a target's is aggregated with it at preflight
  (`is_decoy` only if all rows are decoys), so it counts as a target. The contract test records
  how many such sequences the sample holds (expected 0).
- **q-values at the cut-off.** The writer's `%.6g` can round a value just above 0.01 down to
  `0.01` (for example 0.01000004 → `0.01`), which then passes `q ≤ 0.01` when read back. This is
  the house TSV convention (ADR-0006) and is not changed here. `q_crossings` counts such rows per
  arm. A count above 0 is reported in the README and to M10-05, which computes TE1 at that cut-off
  (Q7). float32 widening itself is exact, and moves no value across the cut-off.
- **Fallback columns.** If `Global.Q.Value` were absent, the converter would use `Q.Value` for the
  whole file. That is legal for §5.4 but wrong for this pilot (run-level, not global). The
  contract test (F1) fails instead of letting it pass.
- **Multiple inputs.** One report per arm. The converter would concatenate several reports
  without de-duplicating (M1-07); the job passes exactly one.
- **Protein groups.** `Protein.Group` values like `A;B` are written as-is to the abundance table;
  `read_abundance` resolves the members. Duplicate group strings cannot occur (one row per group).
- **Large report.** Only mapped columns are loaded (`pq.read_schema` then `columns=`). A report
  of more rows than sized for gives `OUT_OF_MEMORY`, then one resubmission at 16 GB.
- **Partial outputs.** Outputs are written under temporary names and renamed only after every
  check passes, so a killed job never leaves a table M10-04 could read. A stale `<a>.FAILED`
  from an earlier attempt is removed when the next attempt starts.
- **Concurrency.** A1 and A2 write different names in one directory; no lock is needed.
- **Circular or linear genomes, strands, coordinates:** not applicable. This issue reads search
  output tables; no coordinate is read or written.

## Acceptance criteria

On CI (every commit; no DIA-NN, no real data, no `$GIAP_ROOT`):

1. The new contract test collects, and **skips** with the reasons above when `GIAP_ROOT` is unset
   or the sample is absent; it never errors at collection (no module-level read of
   `$GIAP_ROOT`).
2. For each converter change the sample forced (F1, F2, F7), a unit test on a hand-made parquet
   fixture shows the condition and the fixed output. All existing tests in
   `tests/unit/test_convert_diann.py` and `test_cli_convert.py` pass unchanged.
3. New CI fixtures in the 2.3.1 shape (built in `tmp_path` with pyarrow):
   - float32 `Global.Q.Value`, `PEP`, `Precursor.Normalised` and `PG.MaxLFQ`;
   - int8 and bool `Decoy`;
   - dictionary-encoded `Run` and `Protein.Ids`;
   - a handful of unmapped extra columns;
   - `Cont_` IDs;
   - A2 locus-tag IDs.

   They give examples 1–7 exactly. The written TSV has the §5.3 column order and `%.6g` floats
   (`0.0021`, not `0.0020999999`).
4. Example 3: a decoy with NA `Global.Q.Value` and a present `Q.Value` gives `q_value` NA for
   that row; target rows still take `Global.Q.Value`.
5. `diann_convert.sbatch` in a sandbox (a hand-made report, its `check.json` and
   `outputs.sha256`; `GIAP_ROOT`/`OUT_DIR` in `tmp_path`):
   - success gives the three outputs plus `convert.json` with exit 0;
   - examples 8 and 9 and a present `FAILED` give exit 2 with nothing written;
   - a forced converter failure (a report missing `Stripped.Sequence`, its sha256 listed) gives
     exit 1, no table, and an `<a>.FAILED`;
   - a stubbed non-deterministic second output gives exit 1;
   - nothing is written under the sandbox repo.
6. `summarize` on hand-made frames: every `counts` field computed by hand for a 12-row fixture.
   It includes one `q_crossings` row (`Global.Q.Value=0.01000004`) and one mixed contaminant row.
   Each summary rule fails on a fixture built to break it.
7. House rules for the new sbatch (`test_sbatch_house_rules`):
   - a `giap-*` name;
   - `--partition=normal`, `--qos=normal`, `--cpus-per-task=2`, `--mem=8G`,
     `--time=00:30:00`;
   - `--output` under `local/logs/`;
   - the environment sourced, `set -euo pipefail`;
   - no `uv run`, `uv sync` or `pip`.
8. The test helper's fixtures name only the columns they use; no test file contains the real
   report's column list (checked by review, Q-K2).
9. Guard (extends M10-02's recorded-diann guard):
   - `git ls-files validation/pilot/inputs/` ⊆ {`pxd078459_files.tsv`, `acquisition.json`,
     `README.md`};
   - no tracked file anywhere matches `a[1-3]*.tsv*`, `*_abundance.tsv*`, `*.convert.json` or
     `*.convert.log`;
   - nothing under `tests/data/recorded/diann/` other than `2.3.1/PROVENANCE.md`.

On FarmShare (evidence; PR description, PROGRESS, README; not CI):

10. The contract test passes for A1 and A2 on the real samples (`pytest
    tests/contract/test_m10_03_diann_contract.py` with the venv Python and `GIAP_ROOT` set). Its
    output (passed, skipped and recorded values: decoy q-NA count, contaminant rows in the sample)
    is pasted in the PR.
11. Both conversion jobs exit 0. Their job IDs and `sacct` State, Elapsed, MaxRSS and NodeList are
    in PROGRESS › Measured numbers. `a{1,2}.tsv.gz`, `a{1,2}_abundance.tsv.gz` and
    `a{1,2}.convert.json` exist under `$GIAP_ROOT/pilot/inputs/`, and `convert.json` has
    `failures: []`.
12. The README "Conversion (M10-03)" section exists (content per Q-K1), and PROGRESS points
    M10-04 at the two table paths.

Link to §15.6: these tables are the `--peptides` inputs of "GIAP runs end to end on A1, A2"
(M10-04), of TE1 (M10-05, "≥ 99.5 % of A2 peptides map exactly") and of the A1/A2 comparison
(M10-04).

## Test plan

Owner: the test-author, tests first (R2), in new files only:

- `tests/contract/test_m10_03_diann_contract.py`: criteria 1 and 10; F1–F10 and the round-trip,
  parametrized by arm; the `slow` full-report variant (`-k full_report`).
- `tests/unit/test_m10_03_diann_2_3_1_shape.py`: criteria 2–4. Fixtures are written with
  pyarrow in `tmp_path`, with types set explicitly (`pa.float32()`, `pa.int8()`, `pa.bool_()`,
  `pa.dictionary(pa.int32(), pa.string())`). The facts of F2/F7 come from the DIA-NN
  documentation. Once the sample is seen, they also come from the conditions it showed, stated in
  the test docstring as a condition ("`Decoy` stored as bool"), never as copied values.
- `tests/unit/test_m10_03_diann_convert_job.py`: criteria 5 and 7, in a sandbox like
  `helpers.pilot_kit.pilot_sandbox`. The real converter runs on a 10-row report; no fake `giap`.
  The determinism failure uses a `GIAP_CONVERT_CMD` test hook only if the implementer adds one;
  otherwise it is checked by review (Q8).
- `tests/unit/test_m10_03_convert_summary.py`: criterion 6.
- `tests/unit/test_m10_03_guard.py`: criterion 9.

The **sequence** follows from the conditional changes:

1. the test-author writes the contract test and the CI tests for examples 1–7 and the job, and
   commits and pushes;
2. the lead (or the implementer, under the lead's brief) runs the contract test on FarmShare
   against the M10-02 samples. Each failing F becomes a converter change with its own CI test,
   which the test-author adds;
3. the implementer makes the changes, the sbatch and `diann_convert.py`, and the lead retires the
   old contract test (Q1);
4. the jobs run.

If step 2 passes everything, the converter is unchanged and the issue is the test, the job and
the summary.

Local runs: targeted only, with the venv Python and `PYTHONPATH=$PWD/src:$PWD/tests/sim`, never
`uv run` in a worktree. Contract test and unit tests: seconds. The full suite goes to Slurm
(`test_branch.sbatch`), not the session.

## Estimate

| Piece | Code lines |
|---|---|
| converter changes (conditional F1, F2, F7; 0 if the sample matches) | 0–40 |
| `scripts/slurm/diann_convert.sbatch` | 60 |
| `validation/pilot/diann_convert.py` (`summarize`, CLI) | 90 |
| README section | prose |
| **Total** | **≈ 150–190** (index: 150) |

Tests ≈ 380 lines. One PR. Wall time: minutes for the jobs; the wait is on M10-02's searches.

## Questions

### For Khoa (in the next batch; the defaults do not block)

| # | Question | Default (reversible) |
|---|---|---|
| Q-K1 | Do the converted generic tables (`a1.tsv.gz`, `a2.tsv.gz`, the MaxLFQ tables) count as "derived from the DIA-NN reports" under your 2026-10-23 rule? §15.2 step 7 planned to commit them. And may the README carry their sha256s and row counts, as M10-02's `PROVENANCE.md` carries the report's? | derived: kept under `$GIAP_ROOT/pilot/inputs/`, never committed; the README carries commands, job IDs, sha256s and row counts but no rows. If you allow committing, one commit adds the two tables (a few MB, each ≤ 5 MB) and the guard drops that pattern |
| Q-K2 | CI tests need fixtures in the 2.3.1 shape. May a test encode a *fact* the real sample showed (for example "`Decoy` is stored as bool", "IDs are the accession only") without its values or column list? | yes, facts about types and naming only; no values, no column list |

### For the lead

| # | Question | Default (reversible) |
|---|---|---|
| Q1 | The existing `tests/contract/test_diann_contract.py` looks for a top-level `report.*` and byte-goldens that Khoa's rule now forbids. Retire it, and have the test-author write a new file instead of editing it (role guard)? | yes: new file by the test-author; the lead deletes the old one in the implementation step, citing this spec |
| Q2 | Decoy rows with NA `Global.Q.Value`: fill per row from `Q.Value`? | no; one q definition per column. Decoys keep `PEP` and fall under M4-24 (partly scored decoys) |
| Q3 | If 2.3.1 writes `Run` as a path, normalise it to the stem in the converter? | keep verbatim (§5.4 `run ← Run`); M10-02 check rule 3 already matches path or stem. Normalising would be a §5.4 change for all users |
| Q4 | `giap run` uses peptide intensities whenever they exist, so the MaxLFQ tables would be validated but not used for gene abundance. Write them anyway, and should M10-04 pass `--abundance`? | write them (index scope, cheap). M10-04 does **not** pass `--abundance`; M10-05 may use them as an external comparator for metric 4 |
| Q5 | Pass `--search-db other` (A1) and `same` (A2) explicitly instead of the `auto` rule (≥ 95 % of IDs in the gene model)? | explicit. A1's IDs are K-12 accessions and A2's are Nissle locus tags, so `auto` should agree; M10-04 can confirm from the manifest |
| Q6 | Not committing the tables departs from §15.2 step 7 and index decision 5. Record that as a paragraph in ADR-0051 (proposed, M10-02) or as its own ADR? | a paragraph in ADR-0051, added by the lead when Q-K1 is answered |
| Q7 | `%.6g` can round a q just above 0.01 to `0.01`, which then passes the cut-off. Act on it in the pilot? | count it (`q_crossings`) and report it; no writer change (ADR-0006 governs every table). If the count is above 0, M10-05 states it next to TE1 |
| Q8 | The determinism step is hard to fail in a test without a hook. Add a `GIAP_CONVERT_CMD` override to the sbatch for the test, or check the step by review only? | review only, unless the implementer adds the hook at no extra cost |
| Q9 | The job runs the converter from the worktree through `PYTHONPATH=$GIAP_REPO/src` before the PR merges, so the evidence lands in the same PR (as M10-02 Q10). Acceptable? | yes; `convert.json` records the commit and the dirty flag, and the worktree must be clean at run time |
