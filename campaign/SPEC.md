# giap-campaign: new genes and ORFs from public prokaryote proteomics

Status: draft v2 (2026-10-26). Owner: Khoa. This directory is **outside GIAP v1** (`src/giap` is not
touched): it catalogues public datasets, re-searches their raw spectra against a six-frame strain
database and de novo, runs `giap run` per dataset and search arm, and pools the per-run findings.

## 0. Decisions (Khoa, 2026-10-26)

| # | Decision |
|---|---|
| D1 | Goal: new genes and ORFs (GIAP novel ORFs and start sites). The other GIAP findings (ALT residues, operon flags, islands, detected hypothetical proteins) are pooled too. |
| D2 | All acquisition types (DIA and DDA). Only organisms with a **complete** genome of the **exact** taxon. |
| D3 | **Re-search the raw spectra** (v1 of this spec reused submitted results; superseded). Two arms per dataset: **db** = database search against `giap makedb --six-frame-orfs` with decoys; **denovo** = de novo sequencing, with the db results as GIAP `--reference-peptides`. |
| D4 | Separate directory: `campaign/`, package `giap_campaign` under `campaign/src/`. |

Engines (each runs as a separate process in its own environment, R5):

| Acquisition | db arm | denovo arm | Spectra input |
|---|---|---|---|
| DIA | DIA-NN 2.3.1 (the pilot image and two-step flags, `validation/pilot/diann_A2.cfg`), `--report-decoys` | Cascadia 0.0.7 (`giap-cascadia` env, pilot weights) | Bruker `.d` directly for DIA-NN; mzML for DIA-NN (Thermo) and Cascadia |
| DDA | Sage (MIT; new env `campaign/envs/sage.yml`), decoys generated and reported | Casanovo (`giap-casanovo` env, pilot weights) | mzML |

Raw → mzML: Thermo `.raw` via ThermoRawFileParser (Apache-2.0; new env
`campaign/envs/thermorawfileparser.yml`); Bruker `.d` via TIMSCONVERT (`giap-timsconvert`, pilot
route B2). SCIEX `.wiff` is not supported.

## 1. Layout and conventions

```
campaign/
  SPEC.md  README.md
  diann.cfg  pins.tsv
  envs/sage.yml  envs/thermorawfileparser.yml  envs/casanovo.yml   # pinned versions: Sage 0.15.0, ThermoRawFileParser 1.4.5, Casanovo 5.2.1
  src/giap_campaign/
    __init__.py  __main__.py  cli.py
    http.py      # §2
    pride.py     # §3
    ncbi.py      # §4
    catalog.py   # §5
    genome.py    # §6
    spectra.py   # §7 raw-file selection, acquisition mode
    engines.py   # §8 command builders (pure: argv lists and config files)
    sage.py      # §9 Sage results -> generic peptide table
    casanovo.py  # §9 Casanovo mzTab -> generic peptide table
    driver.py    # §10
    aggregate.py # §11
  slurm/campaign_array.sbatch
  tests/  conftest.py (puts campaign/src on sys.path), test_*.py, data/
```

- Python ≥ 3.11, type hints, stdlib + pandas (both in the GIAP venv). May import `giap` for
  `giap.io.converters.base.{ConvertResult, finish_peptides}`, `giap.io.converters.strip.strip_modifications`,
  `giap.io.issues.Issue`, `giap.io.peptides.PEPTIDE_COLUMNS`, `giap.annotate.bakta_cmd.bakta_db_path`.
- Pure functions for parsing and command building; network and subprocess calls only in `http.py`
  and `driver.py`, always injectable (`opener`, `run`, `download` parameters) — tests never touch the
  network or real tools.
- Network (R6, R12): https only, TLS verified, one retry, then the dataset fails (the campaign goes
  on). `NCBI_API_KEY`, if set, is sent as header `api-key` **only** to `api.ncbi.nlm.nih.gov` and is
  never logged or written.
- Output TSVs: UTF-8, tab, one header row, missing = empty cell (never `NA`), deterministic row order.
  Lists in a cell are `;`-joined.
- Commands: `uv run pytest campaign/tests`, `uv run ruff check campaign && uv run ruff format --check campaign`,
  `uv run mypy campaign/src`.

## 2. `http.py`

- `get_json(url, *, cache_dir=None, timeout=60, retries=1, opener=None) -> Any`. Non-`https://` →
  `ValueError` before any fetch. Builds a `urllib.request.Request` (Accept: application/json; the
  NCBI key header per §1) and calls `opener(request, timeout=timeout)` (default
  `urllib.request.urlopen`). With `cache_dir`, the body is stored as `<sha256(url)>.json` and reused;
  errors are never cached. Retries once on `URLError`, `TimeoutError` or HTTP ≥ 500; HTTP 404 raises
  `NotFound(LookupError)` at once; other 4xx raise `HttpError` at once; after the retries the last
  error is raised as `HttpError` (subclass of `OSError`).
- `download(url, dest, *, expected_bytes=None, checksum=None, run=subprocess.run) -> str`.
  Non-https → `ValueError`. Runs `curl -fsSL --proto =https --proto-redir =https --retry 1 -o <dest>.part <url>`
  (`<dest>.part` = `str(dest) + ".part"`); a non-zero exit → `DownloadError`. Then size check, then
  checksum check (32 hex = MD5, 40 = SHA-1, 64 = SHA-256, case-insensitive; anything else ignored),
  then atomic rename; returns the file's SHA-256. Any mismatch removes the `.part` file and raises
  `DownloadError`. An existing `dest` that passes the checks is not downloaded again; one that fails
  them is replaced.
- `safe_extract(archive, out_dir) -> list[Path]` for `.zip`, `.tar`, `.tar.gz`/`.tgz`: refuses absolute
  paths, `..` components, symlinks and hard links (`DownloadError`); returns the extracted **files**.

## 3. `pride.py` (PRIDE Archive REST v3, `https://www.ebi.ac.uk/pride/ws/archive/v3`)

Endpoints: `GET /projects?pageSize=N&page=P` (a JSON list; an empty list ends paging; a dict holding
the list under `_embedded.projects`, `list` or `content` is accepted), `GET /projects/{acc}`,
`GET /projects/{acc}/files/all` (JSON list of file records).

`parse_project(record) -> Project` (frozen dataclass): `accession`; `title` ("" if missing);
`submission_type`; `taxids` (from `organisms[]` items: a dict's `accession` such as `"562"`,
`"NEWT:562"`, `"NCBITaxon:562"`, or a string; the trailing integer; others skipped; order kept,
duplicates removed); `organism_names`; `instruments` (`instruments[].name`); `experiment_types`;
`software` (`softwares[].name`); `keywords`; `pubmed_ids` (`references[].pubmedID`, as str; 0, "0",
null, "" dropped); `dois` (non-empty `references[].doi`).

`acquisition(project) -> "DIA" | "DDA" | "unknown"`: case-insensitive whole-word search of experiment
types, keywords and title: DIA if any of `DIA`, `SWATH`, `SWATH-MS`, `data-independent`, `diaPASEF`;
else DDA if any of `DDA`, `data-dependent`, `ddaPASEF`, `shotgun`; else unknown.

`parse_file(record) -> PrideFile` (frozen): `name` ← `fileName`; `category` ← `fileCategory.value`
upper-cased ("" if missing); `size_bytes` ← `fileSizeBytes` (int or None); `url` ← the
`publicFileLocations[]` item named `FTP Protocol`, prefix `ftp://ftp.pride.ebi.ac.uk/` replaced by
`https://ftp.pride.ebi.ac.uk/`; an `https://` value is kept; anything else → None; `checksum` ←
`checksum` if a non-empty string.

## 4. `ncbi.py` (NCBI Datasets v2, `https://api.ncbi.nlm.nih.gov/datasets/v2`)

`GET /taxonomy/taxon/{taxid}` → `{"taxonomy_nodes": [{"taxonomy": {tax_id, organism_name, rank?,
lineage: [int], classification: {rank: {name, id}}}}]}`. `parse_taxon(report) -> Taxon` (frozen:
`tax_id`, `name`, `rank` upper-case or "", `lineage` tuple, `species_id` int | None, `domain`
str | None); a missing node or an `errors` entry → `LookupError`.
`is_prokaryote(t)`: 2 or 2157 is `t.tax_id` or in `t.lineage`.
`is_below_species(t)`: `t.species_id is not None and t.species_id != t.tax_id`.

`GET /genome/taxon/{taxid}/dataset_report?filters.assembly_level=complete_genome&tax_exact_match=true&page_size=1000`
(+ `page_token` while `next_page_token` is present) → `{"reports": [...]}` with `accession`,
`organism.{tax_id, organism_name, infraspecific_names.strain}`,
`assembly_info.{assembly_level, refseq_category, release_date}`.
`parse_assembly(report) -> Assembly` (frozen: `accession`, `tax_id`, `organism`, `strain`, `level`,
`refseq_category`, `release_date`). `pick_assembly(assemblies, tax_id) -> (Assembly | None, reason)`:
keep `level == "Complete Genome"` with that `tax_id`; none → `(None, "no_complete_genome")`; rank by
`refseq_category` in {`reference genome`, `representative genome`} first, then `GCF_` before `GCA_`,
then latest `release_date`, then accession ascending; reason `unique` or `chosen_of_<n>`.

## 5. `catalog.py` → sample sheet (`sheet.tsv`)

`build_sheet(projects, *, taxon, assemblies, files) -> pd.DataFrame` with injected `taxon(taxid) ->
Taxon`, `assemblies(taxid) -> list[Assembly]`, `files(accession) -> list[PrideFile]`; a
`LookupError`/`OSError` from them sets a status (below), never aborts. `taxon` is memoised per taxid
within one call. One row per (project, prokaryotic taxid); projects with none give no row and
`files` is not called for them; otherwise `files` is called once per project.

| column | value |
|---|---|
| `dataset` | accession |
| `tax_id` (int), `organism` | from the taxon |
| `below_species` | `true` / `false` |
| `title`, `submission_type` | §3 |
| `acquisition` | §3 metadata guess (`DIA`, `DDA`, `unknown`); the data decide at run time (§7) |
| `instruments`, `pubmed_ids`, `dois` | `;`-joined |
| `genome_accession` | picked assembly, else empty |
| `genome_choice` | `pick_assembly` reason, or `species_level` |
| `raw_format` | chosen raw format (§7), else empty |
| `n_raw`, `raw_gb` | count and total size (GB = 1e9 bytes, 2 decimals) of the chosen format's files |
| `raw_files`, `raw_urls`, `raw_bytes`, `raw_checksums` | `;`-joined, same order (empty items kept) |
| `formats_seen` | `;`-joined, sorted |
| `flags` | `;`-joined, sorted: `multi_organism` (> 1 taxid of any kind), `partial_submission` |
| `status` | first that applies: `taxon_error`, `no_raw_files` (no supported raw format; `sciex_wiff` only counts as unsupported), `needs_genome` (species level, no complete genome, or assembly lookup failed), else `ready` |

Rows sorted by `dataset`, then numeric `tax_id`. Curators may edit `genome_accession` and `status`;
the driver runs only `ready` rows.

`iter_projects(get_json, *, page_size=100, max_projects=None, accessions=None) -> Iterator[dict]`:
pages `/projects` until an empty page or `max_projects` records; with `accessions`, fetches
`/projects/{acc}` for each, in the given order.

## 6. `genome.py`

- NCBI zip URL: `/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GBFF`
  → `ncbi_dataset/data/{acc}/*_genomic.fna` and `ncbi_dataset/data/{acc}/genomic.gbff`.
- `find_genome_files(paths) -> (fna, gbff)`: exactly one of each among the paths, else `ValueError`.
- `gbff_topologies(gbff) -> dict[str, str]`: `VERSION` accession and `LOCUS` name → `circular` /
  `linear` (from the LOCUS line tokens; `linear` when neither occurs). gzip accepted.
- `tag_fasta(src, dst, topologies) -> {"circular", "linear", "unknown": counts}`: byte-for-byte copy,
  appending ` [topology=circular]` to headers whose first token is circular; unknown first tokens are
  written as is and counted. gzip input; plain output.

## 7. `spectra.py`

`classify_raw_file(name) -> str | None` on the lower-cased basename:

| format | rule |
|---|---|
| `thermo_raw` | ends with `.raw` |
| `bruker_d` | ends with `.d.zip`, `.d.tar`, `.d.tar.gz`, `.d.tgz` |
| `mzml` | ends with `.mzml` or `.mzml.gz` |
| `sciex_wiff` | ends with `.wiff` or `.wiff2` (unsupported) |
| None | anything else (including `.wiff.scan`, `.mgf`, result files) |

`choose_raw(files) -> (format | None, list[PrideFile], formats_seen)`: among files with a URL and a
supported format (`thermo_raw`, `bruker_d`, `mzml`), the format with the most files (ties:
`thermo_raw` > `bruker_d` > `mzml`); its files sorted by name. `formats_seen` includes `sciex_wiff`.

`select_runs(files, *, max_runs=12, max_gb=100.0) -> list[PrideFile]`: from the name-sorted list
take `k = min(max_runs, n)` evenly spaced files (indices `round(i * (n - 1) / (k - 1))`, i = 0…k−1;
`k = 1` → the first), then lower `k` by one until their total size ≤ `max_gb` (files of unknown
size count 0). Always at least one file. Deterministic.

Acquisition from the data (the metadata guess is only a fallback):
- `tdf_acquisition(analysis_tdf) -> "DIA" | "DDA" | "unknown"`: read-only sqlite,
  `SELECT DISTINCT MsMsType FROM Frames`; 9 → DIA, 8 → DDA, both or neither → unknown.
- `mzml_acquisition(mzml, max_spectra=200) -> "DIA" | "DDA" | "unknown"`: stream the first
  `max_spectra` MS2 spectra (`ms level` cvParam `MS:1000511` = 2); isolation width = lower offset
  (`MS:1000828`) + upper offset (`MS:1000829`); median width ≥ 5 Th → DIA, < 5 → DDA; no MS2 or no
  offsets → unknown. gzip accepted.
- `resolve_acquisition(data_modes, metadata) -> (mode, note)`: all data modes equal and known → that
  mode (note `metadata_disagrees` when the metadata says the other mode); mixed known modes →
  `("mixed", "mixed_acquisition")`; all unknown → metadata if known (note `from_metadata`), else
  `("unknown", "undetermined")`. `mixed` and `unknown` fail the dataset.

## 8. `engines.py` (pure command and config builders)

Each returns an argv list (`list[str]`) or file content; the driver runs them.

- `bakta_argv(genome, outdir, db, threads, complete) -> argv`: `micromamba run -n giap-bakta bakta
  --db <db> --output <outdir> --prefix genome --threads <t> --keep-contig-headers [--complete] <genome>`
  (`--complete` when every contig is circular; db from `bakta_db_path("auto", db_dir)`).
- `makedb_argv(annotation_json, genome, outdir, db_dir) -> argv`: `giap makedb --annotation … --genome …
  --six-frame-orfs --min-orf-aa 20 --db-dir … --outdir …`.
- `thermo_to_mzml_argv(raw, outdir)`: `micromamba run -n giap-thermorawfileparser ThermoRawFileParser
  -i <raw> -o <outdir> -f 2` (indexed mzML, vendor peak picking on).
- `bruker_to_mzml_argv(d_dir, outdir)`: `micromamba run -n giap-timsconvert timsconvert --input <d>
  --outdir <outdir>`.
- `diann_argvs(spectra, fasta, outdir, threads, cfg_path, sif)` → two argvs (library prediction, then
  search) built like `scripts/slurm/diann_pilot.sbatch`: flags read from `cfg_path` (one per line,
  `[step1]`/`[step2]` sections as in `validation/pilot/diann_A2.cfg`), plus `--fasta`, `--out-lib`,
  `--lib`, `--f` per spectra file, `--out <outdir>/report.parquet`, `--threads`, and always
  `--report-decoys`; run as `apptainer exec <sif> /diann-2.3.1/diann-linux …`. The campaign's own
  `campaign/diann.cfg` holds the A2 flags without the `--fasta` line (the strain database differs per
  dataset); `$VAR`/`${VAR}` in it are expanded from the environment.
- `sage_config(fasta, mzmls, outdir) -> dict` (written as JSON): `database.fasta`,
  `database.enzyme = {cleave_at: "KR", restrict: null, missed_cleavages: 1, min_len: 7, max_len: 30}`
  (GIAP's digest, `giap.mapping.digest`), `database.static_mods = {"C": 57.0215}`,
  `database.variable_mods = {"M": [15.9949]}`, `database.decoy_tag = "rev_"`,
  `database.generate_decoys = true`, `precursor_tol = {ppm: [-20, 20]}`,
  `fragment_tol = {ppm: [-20, 20]}`, `output_directory`, `mzml_paths`;
  `sage_argv(config_json) -> micromamba run -n giap-sage sage <config_json>`.
- `cascadia_argv(mzml, out_prefix, model)`: as `scripts/slurm/denovo_pilot.sbatch` (`cascadia sequence
  <mzml> <model> -o <prefix> -b 32` in `giap-cascadia`).
- `casanovo_argv(mzmls, out_mztab, model)`: `micromamba run -n giap-casanovo casanovo sequence
  --model <model> --output <out_mztab> <mzml …>`.
- `giap_run_argv(genome, annotation, peptides, outdir, db_dir, profile, threads, reference=None)`:
  `giap run --genome … --annotation … --peptides … [--reference-peptides …] --outdir … --db-dir …
  --profile … --threads … --resume`.
- `giap_compare_argv(db_out, denovo_out, outdir)`: `giap compare <db_out> <denovo_out> --outdir …
  --labels db,denovo`.

`campaign/pins.tsv` (`name`, `path`, `sha256`; `path` may use `$GIAP_DB_DIR`/`$GIAP_ROOT`) lists the
pinned external files: `diann_sif`, `cascadia_model`, `cascadia_unimod`, `casanovo_model`. The SHA-256
values are copied from the pilot pins (`validation/pilot/diann_image.sha256`, `validation/pilot/a3_assets.tsv`);
`casanovo_model` has no pin yet (Khoa records it when the weights are fetched): a job that needs a
file whose pin is empty or does not match fails with an environment error (exit 2) before running.
Model weights and the DIA-NN image are the pilot's pinned files (`$GIAP_DB_DIR/cascadia/…`,
`$GIAP_DB_DIR/casanovo/…`, `$GIAP_ROOT/images/diann_2.3.1.sif`); the driver verifies them against the
SHA-256 pins the pilot scripts use before the first use.

## 9. Converters (pure; return `ConvertResult` with `abundance=None`)

Common: gzip inputs accepted; several inputs concatenated; `peptide` =
`strip_modifications(raw).upper()`, empty results dropped with an `empty_peptide` warn issue
(count in the message); `charge` integer; missing required columns → one `missing_column` error and
an empty table; output through `finish_peptides`.

**Sage `results.sage.tsv`** (`convert_sage(paths)`): required `peptide`, `proteins`, `filename`,
`scannr`, `label`, `spectrum_q`. `modified_peptide` ← `peptide`; `protein_ids` ← `proteins`
(`;`-separated already); `is_decoy` ← `label == -1`; `q_value` ← `peptide_q` if present else
`spectrum_q`; `pep` ← `10 ** posterior_error`, clipped to [0, 1] (Sage writes log10(PEP); when every
value is ≥ 1, the placeholder Sage writes without its LDA rescoring, `pep` is left empty); `score` ← `sage_discriminant_score`, else
`hyperscore`; `run` ← `filename` without directory and without `.mzML`/`.mzML.gz`; `spectrum_id` ←
`scannr`; `charge` ← `charge`; `source` = `database`; `search_db` = `same`. Only `rank == 1` rows are
kept when `rank` is present. The contract test uses a hand-written sample and is marked
`pending_sample` until a recorded Sage output exists (GIAP §5.4 practice).

**Casanovo mzTab** (`convert_casanovo(paths)`): `MTD ms_run[i]-location` → run name (basename,
`.gz` then the extension removed); PSH/PSM rows; `peptide` ← `sequence` (stripped; Casanovo writes
mass-shift notation like `+42.011` and `M[Oxidation]`); `modified_peptide` ← `sequence`;
`score` ← `search_engine_score[1]`; `charge` ← `charge`; `run`, `spectrum_id` ← `spectra_ref`
`ms_run[i]:<id>`; `null` = missing; `source` = `denovo`; `search_db` = `same`. Required `sequence`,
`spectra_ref`.

DIA-NN and Cascadia outputs use GIAP's own converters (`giap convert --from diann|cascadia`).

## 10. `driver.py`

**plan** `plan(sheet, root, *, max_runs=12, max_gb=100.0, denovo_runs=3, disk_budget_gb=300.0) ->
(jobs DataFrame, throttle int)`: `ready` rows with a `genome_accession` and a `raw_format`; per row
`select_runs` on its raw files; columns `job` (0-based, sheet order), `dataset`, `tax_id`,
`organism`, `genome_accession`, `raw_format`, `acquisition_metadata`, `run_files`, `run_urls`,
`run_bytes`, `run_checksums` (`;`-joined, selected runs), `run_gb`, `denovo_runs` (the first
`denovo_runs` selected files by the same even spacing), `workdir` = `root/<dataset>_<tax_id>`.
`throttle = max(1, min(8, floor(disk_budget_gb / (3 × max run_gb))))` (raw + mzML + search work);
the CLI prints `sbatch --array=0-<N-1>%<throttle> …`.

**run_job** `run_job(job, *, db_dir, profile, threads, keep_inputs, run, download, env) -> int`.
Steps run in `workdir`, each recorded in `workdir/status.json`
(`{"steps": {name: {"state": "done"|"failed", "message": str}}, "acquisition": …}`); a `done` step
is skipped; the first failure records `failed` with its message and returns 1; all done → 0.

1. `fetch_genome` → `genome/genome.fna` (download zip, `safe_extract`, `find_genome_files`,
   `gbff_topologies`, `tag_fasta`).
2. `annotate` → `annotation/genome.json` (`bakta_argv`).
3. `makedb` → `searchdb/strain_db.faa` (`makedb_argv`).
4. `fetch_raw` → `raw/` (download each run with its bytes and checksum; Bruker archives
   `safe_extract`ed into `raw/<name>.d/`).
5. `acquisition` → per run `tdf_acquisition` / (after step 6 for Thermo) `mzml_acquisition`;
   `resolve_acquisition` with the metadata guess; stored in `status.json`.
6. `to_mzml` → `mzml/<run>.mzML` for every run that needs it: Thermo (all), Bruker (only the
   `denovo_runs` and, for DDA, all runs). mzML inputs are used as they are.
7. `search_db` → DIA: `diann_argvs` on `.d` (Bruker) or mzML (others) → `search/report.parquet`;
   DDA: `sage_config` + `sage_argv` → `search/results.sage.tsv`.
8. `denovo` → DIA: `cascadia_argv` per denovo run → `denovo/<run>.ssl`; DDA: `casanovo_argv` on the
   denovo runs → `denovo/denovo.mztab`.
9. `convert` → `peptides_db.tsv` (DIA: `giap convert --from diann … --search-db same`; DDA:
   `convert_sage` → `peptides_db.raw.tsv` → `giap convert --from generic … --search-db same`) and
   `peptides_denovo.tsv` (Cascadia: `giap convert --from cascadia …`; Casanovo: `convert_casanovo`
   → `giap convert --from generic …`). Error issues or a non-zero exit → failed.
10. `giap_db` → `giap_run_argv(…, peptides_db.tsv)` → `giap_db/`.
11. `giap_denovo` → `giap_run_argv(…, peptides_denovo.tsv, reference=peptides_db.tsv)` → `giap_denovo/`.
12. `compare` → `giap_compare_argv` → `compare/`.
13. `cleanup` unless `keep_inputs`: delete `raw/`, `mzml/`, `search/` work files except the report
    and results file, and `giap_*/work/`.

`run(argv, cwd) -> returncode` and `download(url, dest, expected_bytes, checksum)` are injected;
`env` maps the pinned paths (`db_dir`, `diann_sif`, `cascadia_model`, `casanovo_model`, `diann_cfg`).

**Slurm** `slurm/campaign_array.sbatch`: `--job-name=giap-campaign`, `--cpus-per-task=16`,
`--mem=64000M`, `--time=3-00:00:00`, `--output=…/local/logs/%x.%A_%a.out`; sources
`scripts/farmshare_env.sh`; runs `uv run --no-sync python -m giap_campaign run-job --jobs "$1"
--index "$SLURM_ARRAY_TASK_ID" --db-dir "$2" --profile farmshare_local,apptainer --threads 16`.
Submit as printed by `plan` (`%throttle` ≤ 8: 8 × 16 = 128 CPUs, R14).

## 11. `aggregate.py`

`aggregate(jobs, out_dir) -> dict[str, int]` (rows written per table). For each job and each arm
(`db` → `giap_db/`, `denovo` → `giap_denovo/`) whose `run_manifest.json` exists, read the tables below
when present; prefix rows with `dataset`, `tax_id`, `organism`, `genome_accession`, `arm`. Missing
tables are skipped; unreadable ones are noted in `campaign_runs.tsv`.

| output | source | rows | extra columns |
|---|---|---|---|
| `campaign_novel_orfs.tsv` | `tables/novel_orfs.tsv` | all | `arms_called` = `;`-joined sorted arms of this dataset that call the same `orf_id`; `n_datasets_called` = distinct datasets with the same `genome_accession` and `orf_id` called in any arm |
| `campaign_novel_orfs.faa` | `tables/novel_orfs.faa` | all records | header ID prefixed `<dataset>|<arm>|<genome_accession>|` |
| `campaign_variants.tsv` | `tables/variants.tsv` | all | `n_datasets` = distinct datasets with the same `genome_accession`, `contig`, `codon_start`, `strand`, `alt_aa` |
| `campaign_operon_flags.tsv`, `campaign_islands.tsv` | `tables/operon_flags.tsv`, `tables/islands.tsv` | db arm only | — |
| `campaign_detected_hypothetical.tsv` | `tables/genes.tsv` | db arm; `hypothetical` and `detected` true | `locus_tag`, `product`, `contig`, `start`, `end`, `strand`, `length_aa`, `detection_tier`, `gene_q_value`, `n_peptides_unique`, `abundance_log2` |
| `campaign_runs.tsv` | every job × arm | one row each | `state` (`done`, `failed:<step>` from `status.json`, `not_run`), `acquisition`, `giap_version` (manifest), `n_genes_detected`, `n_novel_called`, `n_alt_high`, `notes` |

Rows sorted by `dataset`, numeric `tax_id`, `arm`, then the table's own order. Per-run q-values and
FDRs are kept unchanged; no pooled error rate is claimed (R4).

## 12. CLI (`python -m giap_campaign`, argparse; `cli.main(argv) -> int`)

```
catalog   --out sheet.tsv [--cache-dir DIR] [--max-projects N] [--accessions PXD1,PXD2]
plan      --sheet sheet.tsv --root DIR --out jobs.tsv [--max-runs 12] [--max-gb 100] [--denovo-runs 3] [--disk-budget-gb 300]
run-job   --jobs jobs.tsv --index I --db-dir DIR [--profile farmshare_local,apptainer] [--threads 16] [--keep-inputs]
convert   --format sage|casanovo INPUT [INPUT ...] --out peptides.raw.tsv
aggregate --jobs jobs.tsv --out DIR
```

Exit codes: 0 success; 1 invalid input or a failed job; 2 network or environment error.

## 13. Later

Cross-species clustering of novel proteins; homology and literature checks of candidates (UniProt,
RefSeq, sORF catalogues, Europe PMC); GPU partition for de novo; Sage LFQ intensities.
