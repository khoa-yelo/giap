# giap-campaign: GIAP over many public prokaryote proteomics datasets

Status: draft v1 (2026-10-26). Owner: Khoa. This directory is **outside GIAP v1** (`src/giap` is not
touched): it is tooling that catalogues public datasets, converts the submitters' own search
results, runs `giap run` once per dataset on FarmShare, and pools the per-run tables.

## 0. Decisions (Khoa, 2026-10-26)

| # | Decision |
|---|---|
| D1 | "Novel" covers every GIAP finding: novel ORFs and start sites, protein variants (ALT), operon and regulatory flags, abundance islands, and detected hypothetical proteins. |
| D2 | All acquisition types (DIA and DDA). Only organisms with a **complete** genome of the **exact** taxon. |
| D3 | Reuse each submission's own search results; no re-search of raw spectra in v1 of the campaign. |
| D4 | Separate directory: `campaign/`, package `giap_campaign` under `campaign/src/`. |

Consequence of D3 (told to Khoa): submitters searched annotated protein databases, so their peptides
rarely fall outside annotated genes. Novel-ORF yield comes mainly from differences between the
submitter's database and Bakta's annotation; operon, island, detection and ALT findings are fully
supported. A re-search arm (`giap makedb --six-frame-orfs` + search) can be added later per dataset.

## 1. Layout and conventions

```
campaign/
  SPEC.md  README.md
  src/giap_campaign/
    __init__.py  __main__.py  cli.py
    http.py        # https-only GET with JSON cache, file download with checksum (§2)
    pride.py       # PRIDE Archive v3 records (§3)
    ncbi.py        # NCBI Datasets v2 taxonomy and genome records (§4)
    catalog.py     # sample sheet (§5)
    genome.py      # NCBI genome zip -> topology-tagged FASTA (§6)
    converters/    # __init__.py (CONVERTERS), common.py, maxquant.py, fragpipe.py, mztab.py, mzid.py (§7)
    driver.py      # plan + run-job (§8)
    aggregate.py   # pooled tables (§9)
  slurm/campaign_array.sbatch
  tests/  conftest.py (puts campaign/src on sys.path), test_*.py, data/ (small text fixtures, gzip if > 20 kB)
```

- Python ≥ 3.11, type hints everywhere, stdlib + pandas only (both in the GIAP venv). The package may
  import `giap` (installed in the same venv) for `giap.io.converters.base.{ConvertResult, finish_peptides}`,
  `giap.io.converters.strip.strip_modifications`, `giap.io.issues.Issue`, `giap.io.peptides.PEPTIDE_COLUMNS`.
- Pure functions for parsing; network and subprocess calls only in `http.py` and `driver.py`, always
  injectable (a `fetch`/`run` callable parameter) so tests never touch the network.
- Network rules (R6, R12): https only, TLS verified, one retry, then fail that dataset (never abort
  the whole campaign). An optional `NCBI_API_KEY` env var is sent as header `api-key` and never logged.
- Output TSVs: UTF-8, tab, one header row, `NA` never written (empty cell = missing), rows sorted
  deterministically as stated per table. Lists inside a cell are `;`-joined.
- Commands: `uv run pytest campaign/tests`, `uv run ruff check campaign && uv run ruff format --check campaign`,
  `uv run mypy campaign/src`.

## 2. `http.py`

- `get_json(url, *, cache_dir: Path | None = None, timeout=60, retries=1, opener=None) -> Any`.
  Rejects non-`https://` URLs (`ValueError`). With `cache_dir`, the response body is stored as
  `<sha256(url)>.json` and reused on later calls. Retries once on `URLError`/`HTTPError` ≥ 500/timeout;
  HTTP 404 raises `NotFound` (subclass of `LookupError`) without retry. `opener` defaults to
  `urllib.request.urlopen`; tests pass a fake.
- `download(url, dest: Path, *, expected_bytes: int | None, checksum: str | None, run=subprocess.run) -> str`:
  `curl -fsSL --proto =https --proto-redir =https --retry 1 -o dest.part url`, then size check
  (`expected_bytes` when given), checksum check (32 hex = MD5, 40 hex = SHA-1, 64 hex = SHA-256,
  case-insensitive; anything else ignored), then atomic rename to `dest`; returns the SHA-256 of the file.
  A mismatch deletes `dest.part` and raises `DownloadError`. An existing `dest` that passes the checks is
  not downloaded again.
- `safe_extract(zip_path, out_dir) -> list[Path]`: extracts a zip, refusing absolute paths, `..`
  components and symlinks (`DownloadError`).

## 3. `pride.py` (PRIDE Archive REST v3, `https://www.ebi.ac.uk/pride/ws/archive/v3`)

Endpoints: `GET /projects?pageSize=N&page=P` (JSON list of project records; an empty list ends
paging; also accept a dict holding the list under `_embedded.projects`, `list` or `content`),
`GET /projects/{acc}`, `GET /projects/{acc}/files/all` (JSON list of file records).

`parse_project(record) -> Project` (frozen dataclass):

| field | from the record |
|---|---|
| `accession` | `accession` |
| `title` | `title` (missing → "") |
| `submission_type` | `submissionType` (`COMPLETE`, `PARTIAL`, …; missing → "") |
| `taxids` | `organisms[]`: each item a dict with `accession` (e.g. `"562"`, `"NEWT:562"`, `"NCBITaxon:562"`) or a string; the trailing integer; items without one are skipped; order kept, duplicates removed |
| `organism_names` | `organisms[].name` |
| `instruments` | `instruments[].name` |
| `experiment_types` | `experimentTypes[].name` |
| `software` | `softwares[].name` |
| `keywords` | `keywords` (list of str) |
| `pubmed_ids` | `references[].pubmedID` (int or str); 0, "0", null and "" dropped; as str |
| `dois` | `references[].doi`, non-empty |

`acquisition(project) -> "DIA" | "DDA" | "unknown"`: case-insensitive search of experiment types,
keywords and title; DIA if any of `DIA`, `SWATH`, `data-independent`, `diaPASEF` occurs as a word;
else DDA if any of `DDA`, `data-dependent`, `shotgun` occurs; else unknown.

`parse_file(record) -> PrideFile` (frozen): `name` ← `fileName`; `category` ← `fileCategory.value`
upper-cased ("" if missing); `size_bytes` ← `fileSizeBytes` (int or None); `url` ← the
`publicFileLocations[]` item named `FTP Protocol`, with prefix `ftp://ftp.pride.ebi.ac.uk/` replaced by
`https://ftp.pride.ebi.ac.uk/`; any other URL scheme → `None`; `checksum` ← `checksum` if a non-empty string.

`classify_result_file(name) -> str | None` on the lower-cased basename, after removing one trailing
`.gz`:

| format | rule |
|---|---|
| `diann` | `report.parquet`, `report.tsv`, or ends with `.report.parquet` / `.report.tsv` / `_report.parquet` / `_report.tsv` and the name contains `diann` or `dia-nn`; never names containing `matrix`, `stats`, `.log` |
| `maxquant` | `evidence.txt` |
| `fragpipe` | `psm.tsv` |
| `mztab` | ends with `.mztab` |
| `mzid` | ends with `.mzid` |
| `archive` | ends with `.zip` and contains one of `txt`, `maxquant`, `combined`, `fragpipe`, `diann`, `result`, `search` |
| None | everything else |

`choose_results(files) -> tuple[str | None, list[PrideFile], list[str]]`: the chosen format, its
files (all files of that format, sorted by name) and all formats seen (sorted). Preference:
`diann` > `maxquant` > `fragpipe` > `mztab` > `mzid`; `archive` is chosen only when nothing else is
present (its contents are checked by the driver, §8). Files with `url is None` are ignored.

## 4. `ncbi.py` (NCBI Datasets v2, `https://api.ncbi.nlm.nih.gov/datasets/v2`)

`GET /taxonomy/taxon/{taxid}` → `{"taxonomy_nodes": [{"taxonomy": {...}}]}` with keys `tax_id`,
`organism_name`, `rank` (e.g. `SPECIES`, `STRAIN`; may be absent), `lineage` (list of int ancestors),
`classification` (dict of rank → `{"name", "id"}`, e.g. `domain`, `species`).

`parse_taxon(report) -> Taxon` (frozen: `tax_id`, `name`, `rank` upper-case or "", `lineage`
tuple, `species_id` int | None, `domain` str | None). Missing node or an `errors` entry →
`LookupError`.

- `is_prokaryote(t)`: 2 (Bacteria) or 2157 (Archaea) is `t.tax_id` or in `t.lineage`.
- `is_below_species(t)`: `t.species_id is not None and t.species_id != t.tax_id`.

`GET /genome/taxon/{taxid}/dataset_report?filters.assembly_level=complete_genome&tax_exact_match=true&page_size=1000`
(+ `page_token` while `next_page_token` is present) → `{"reports": [...], "next_page_token"?}`.
Report keys: `accession`, `organism.tax_id`, `organism.organism_name`,
`organism.infraspecific_names.strain`, `assembly_info.assembly_level`,
`assembly_info.refseq_category`, `assembly_info.release_date`.

`parse_assembly(report) -> Assembly` (frozen: `accession`, `tax_id`, `organism`, `strain`,
`level`, `refseq_category`, `release_date`).

`pick_assembly(assemblies, tax_id) -> tuple[Assembly | None, str]`: keep `level == "Complete Genome"`
and `tax_id` equal; none → `(None, "no_complete_genome")`. Rank by: `refseq_category` in
{`reference genome`, `representative genome`} first; then `GCF_` before `GCA_`; then the latest
`release_date`; then accession ascending. Return `(best, "unique")` when one candidate, else
`(best, f"chosen_of_{n}")`.

## 5. `catalog.py` → sample sheet

`build_sheet(projects, *, taxon, assemblies, files) -> pd.DataFrame`, where `taxon(taxid) -> Taxon`,
`assemblies(taxid) -> list[Assembly]` and `files(accession) -> list[PrideFile]` are injected
(production versions wrap `http.get_json`; a `LookupError` from them sets the row's status as below,
it never aborts). One row per (project, prokaryotic taxid); projects with no prokaryotic taxid give no
row. `files` is called once per project, only when it has a prokaryotic taxid.

| column | value |
|---|---|
| `dataset` | accession |
| `tax_id`, `organism` | from the taxon |
| `below_species` | `true`/`false` |
| `title`, `submission_type`, `acquisition` | §3 |
| `instruments`, `pubmed_ids`, `dois` | `;`-joined |
| `genome_accession` | picked assembly accession, else empty |
| `genome_choice` | the `pick_assembly` reason, or `species_level` when not below species |
| `result_format` | chosen format or empty |
| `result_files`, `result_urls`, `result_bytes`, `result_checksums` | `;`-joined, same order (empty cells kept as empty items) |
| `formats_seen` | `;`-joined |
| `flags` | `;`-joined, sorted: `multi_organism` (project has > 1 taxid of any kind), `partial_submission` (`submission_type == "PARTIAL"`) |
| `status` | first that applies: `taxon_error` (taxon lookup failed), `no_supported_results` (no format), `needs_inspection` (format `archive`), `needs_genome` (species level, or no complete genome, or assembly lookup failed), else `ready` |

Rows sorted by `dataset`, then `tax_id`. Taxon lookups are memoised per taxid within one build.
Curators edit `genome_accession` and `status` by hand; the driver only runs `ready` rows (§8).

`iter_projects(get_json, page_size=100, max_projects=None, accessions=None)`: pages `/projects`
until an empty page (or `max_projects`); with `accessions`, fetches `/projects/{acc}` for each instead.

## 6. `genome.py`

- NCBI download URL: `/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GBFF`
  → zip holding `ncbi_dataset/data/{acc}/*_genomic.fna` and `ncbi_dataset/data/{acc}/genomic.gbff`.
- `find_genome_files(paths) -> tuple[Path, Path]`: the single `.fna` and the single `.gbff` among
  extracted paths; anything else → `ValueError`.
- `gbff_topologies(gbff) -> dict[str, str]`: per record, `VERSION` accession (e.g. `NC_000913.3`)
  → `circular` or `linear` from the `LOCUS` line tokens (`linear` if neither word occurs). Also maps
  the `LOCUS` name. gzip accepted.
- `tag_fasta(src, dst, topologies) -> dict[str, int]`: copies the FASTA, appending
  ` [topology=circular]` to headers whose first token is circular (header kept otherwise); returns
  counts `{"circular": n, "linear": n, "unknown": n}` (unknown = not in the map; written as linear).
  gzip input accepted; output plain text.

## 7. Converters (`converters/`)

Each is `convert_<fmt>(paths: Sequence[Path]) -> giap.io.converters.base.ConvertResult` with
`abundance=None`. Common rules (`common.py`):

- inputs may be gzip; several inputs are concatenated in argument order;
- `peptide` = `strip_modifications(raw).upper()`; rows with an empty result are dropped with an
  `empty_peptide` warn issue (count in the message);
- `source = "database"` on every row; intensity 0 or negative → missing; `charge` integer;
- missing required columns → one `missing_column` error issue and an empty table;
- the result goes through `finish_peptides(df, columns)` (generic column order, stable sort).
- `CONVERTERS = {"maxquant": ..., "fragpipe": ..., "mztab": ..., "mzid": ...}`.

**MaxQuant `evidence.txt`.** Required: `Sequence`, `Raw file`, and `Proteins` or `Leading razor protein`.
`modified_peptide` ← `Modified sequence` with surrounding `_` removed; `protein_ids` ← `Proteins`
(already `;`-separated; empty → `Leading razor protein`); `pep` ← `PEP`; `score` ← `Score`;
`is_decoy` ← `Reverse == "+"` (absent → false); `intensity` ← `Intensity`; `run` ← `Raw file`;
`charge` ← `Charge`; `spectrum_id` ← `MS/MS scan number` as an integer string (empty if missing).
Rows whose `Type` is `MULTI-MATCH` (match-between-runs transfers without MS/MS) are dropped with an
info issue `mbr_dropped`.

**FragPipe `psm.tsv`.** Required: `Peptide`, `Spectrum`. `modified_peptide` ← `Modified Peptide`;
`protein_ids` ← `Protein` plus the comma-separated `Mapped Proteins`, de-duplicated, `;`-joined;
`pep` ← `1 − PeptideProphet Probability` (or `1 − Probability`), clipped to [0, 1]; `score` ←
`Hyperscore`; `intensity` ← `Intensity`; `charge` ← `Charge`; `spectrum_id` ← `Spectrum`;
`run` ← `Spectrum File` basename with a leading `interact-` and a trailing `.pep.xml` removed,
else the `Spectrum` value without its last three `.`-fields. No decoys (psm.tsv is pre-filtered).

**mzTab 1.0 (PSM section).** `MTD ms_run[i]-location` gives run names (basename without extension,
`.gz` removed first); `MTD psm_search_engine_score[i] [CV, ACC, name, ]` names score columns.
`PSH` is the header, `PSM` rows the data; `null` = missing. Required: `sequence`, `spectra_ref`.
`peptide` ← `sequence`; `protein_ids` ← union of `accession` over rows sharing (`PSM_ID`,
`sequence`, `spectra_ref`), which collapse to one row; `charge` ← `charge`;
`run`/`spectrum_id` ← `spectra_ref` `ms_run[i]:<id>`; `is_decoy` ← the first `opt_` column whose name
contains `decoy` (value `1`/`true`); `q_value` ← a score column whose CV accession is in Q (below),
else the first `opt_` column whose lower-cased name contains `q-value`, `q_value` or `qvalue`;
`pep` ← a score column with an accession in P, else an `opt_` column containing `posterior_error`
or ending in `_pep`; `score` ← the first score column whose accession is in S. If the MTD
`software` lines name Casanovo, emit a `use_giap_converter` error (use `giap convert --from casanovo`).

**mzIdentML 1.1/1.2.** Streamed with `xml.etree.ElementTree.iterparse` (namespaces stripped,
elements cleared after use). `DBSequence@id → accession`; `Peptide@id → PeptideSequence`;
`PeptideEvidence@id → (peptide_ref, dBSequence_ref, isDecoy)`; `SpectraData@id → location`
(run = basename without extension, `.gz` removed first). For each
`SpectrumIdentificationResult` (`spectraData_ref`, `spectrumID`) keep the items with `rank="1"` and
`passThreshold="true"`: `peptide` ← the peptide's sequence; `protein_ids` ← the accessions of its
`PeptideEvidenceRef`s (sorted unique); `is_decoy` ← every evidence `isDecoy="true"`; `charge` ←
`chargeState`; `spectrum_id` ← `spectrumID`; `run` ← its SpectraData; `q_value`, `pep`, `score`
from the item's `cvParam`s by accession lists Q, P, S (first match in list order).

CV accession lists (in priority order):
- Q (q-value): `MS:1002354` PSM-level q-value, `MS:1001491` percolator:Q value, `MS:1002054` MS-GF:QValue, `MS:1001868` distinct peptide-level q-value, `MS:1002055` MS-GF:PepQValue.
- P (PEP): `MS:1001493` percolator:PEP.
- S (higher is better): `MS:1002338` Andromeda:score, `MS:1001171` Mascot:score, `MS:1002252` Comet:xcorr, `MS:1002049` MS-GF:RawScore, `MS:1001950` PEAKS:peptideScore.

DIA-NN reports use GIAP's own converter (`giap convert --from diann`).

## 8. `driver.py`

**plan** `plan(sheet, root) -> pd.DataFrame` (jobs table, written as `jobs.tsv`): rows of the sheet with
`status == "ready"`, a non-empty `genome_accession` and `result_format` in
{`diann`, `maxquant`, `fragpipe`, `mztab`, `mzid`, `archive`}; columns `job` (0-based, in sheet
order), `dataset`, `tax_id`, `organism`, `genome_accession`, `result_format`, `result_files`,
`result_urls`, `result_bytes`, `result_checksums`, `workdir` = `root/<dataset>_<tax_id>`.

**run-job** `run_job(job_row, *, db_dir, profile, threads, keep_inputs, run, download)` performs the
steps below in `workdir`, each recorded in `workdir/status.json` as
`{"steps": {name: {"state": "done"|"failed", "message": str}}}`; a step already `done` is skipped
(resumable). The first failure stops the job, records `failed` with the message, and returns exit
code 1 (exit 0 when all steps are done).

1. `fetch_results`: `download` each result URL into `inputs/` with its bytes and checksum.
2. `unpack`: for `archive`, `safe_extract` into `inputs/unpacked/` and re-classify the members with
   `classify_result_file`; pick by the §3 preference; none → failed `no_supported_results_in_archive`.
3. `fetch_genome`: download the NCBI zip (§6) into `genome/`, extract, `tag_fasta` → `genome/genome.fna`.
4. `convert`: `diann` → run `giap convert --from diann <files> --out peptides.tsv`; otherwise
   the campaign converter writes `peptides.raw.tsv`, then run
   `giap convert --from generic peptides.raw.tsv --out peptides.tsv` (validation and normalisation).
   Error issues → failed.
5. `giap_run`: run `giap run --genome genome/genome.fna --peptides peptides.tsv --outdir giap_out
   --db-dir <db_dir> --profile <profile> --threads <threads> --resume`. Exit code ≠ 0 → failed with
   the code in the message.
6. `cleanup`: unless `keep_inputs`, delete `inputs/` and `giap_out/work/`.

`run` (subprocess runner) and `download` are injected; tests use fakes that record the commands.

**Slurm** (`slurm/campaign_array.sbatch`): job name `giap-campaign`, `--cpus-per-task=16`,
`--mem=64000M`, `--time=2-00:00:00`, output `local/logs/%x.%A_%a.out`; sources
`scripts/farmshare_env.sh`; runs `uv run --no-sync python -m giap_campaign run-job --jobs "$1"
--index "$SLURM_ARRAY_TASK_ID" --db-dir "$2"` with profile `farmshare_local,apptainer`.
Submit with `sbatch --array=0-<N-1>%8 campaign/slurm/campaign_array.sbatch jobs.tsv DB` (8 × 16 CPUs =
128 CPUs, the R14 cap).

## 9. `aggregate.py`

`aggregate(jobs, out_dir) -> dict[str, int]` (row counts per written table). For each job whose
`workdir/giap_out/run_manifest.json` exists, read the tables below when present; prefix every row
with `dataset`, `tax_id`, `organism`, `genome_accession`. Missing tables are skipped silently;
an unreadable table is listed in `campaign_runs.tsv` `notes`.

| output | source | rows kept | extra columns |
|---|---|---|---|
| `campaign_novel_orfs.tsv` | `tables/novel_orfs.tsv` | all | `n_datasets_called` = datasets with the same `genome_accession` and `orf_id` called (`called` true) |
| `campaign_variants.tsv` | `tables/variants.tsv` | all | `n_datasets` = datasets with the same `genome_accession`, `contig`, `codon_start`, `strand`, `alt_aa` |
| `campaign_operon_flags.tsv` | `tables/operon_flags.tsv` | all | — |
| `campaign_islands.tsv` | `tables/islands.tsv` | all | — |
| `campaign_detected_hypothetical.tsv` | `tables/genes.tsv` | `hypothetical` and `detected` both true | keeps `locus_tag`, `product`, `contig`, `start`, `end`, `strand`, `length_aa`, `detection_tier`, `gene_q_value`, `n_peptides_unique`, `abundance_log2` |
| `campaign_novel_orfs.faa` | `tables/novel_orfs.faa` | all records | header ID prefixed `<dataset>|<genome_accession>|` |
| `campaign_runs.tsv` | every job | one row per job | `state` (`done`, `failed:<step>`, `not_run`), `giap_version` (manifest `giap_version` if present), `n_genes_detected`, `n_novel_called`, `n_alt_high`, `n_operon_flags`, `n_islands`, `notes` |

Rows sorted by `dataset`, `tax_id`, then the table's own row order. GIAP's per-run q-values and
FDRs are kept unchanged; the campaign tables claim no pooled error rate (each row keeps its run's
error estimate, R4).

## 10. CLI (`python -m giap_campaign`, argparse)

```
catalog   --out sheet.tsv [--cache-dir DIR] [--max-projects N] [--accessions PXD1,PXD2]
plan      --sheet sheet.tsv --root DIR --out jobs.tsv
run-job   --jobs jobs.tsv --index I --db-dir DIR [--profile farmshare_local,apptainer] [--threads 16] [--keep-inputs]
convert   --format maxquant|fragpipe|mztab|mzid INPUT [INPUT ...] --out peptides.raw.tsv
aggregate --jobs jobs.tsv --out DIR
```

Exit codes: 0 success, 1 invalid input or a failed job, 2 network or environment error.

## 11. Later (not in this draft)

Re-search arm for novel ORFs (`giap makedb --six-frame-orfs`); cross-species clustering of novel
proteins; homology and literature checks of candidates (UniProt, RefSeq, sORF catalogues, Europe PMC).
