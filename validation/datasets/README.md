# Additional real-data validation datasets (beyond the E. coli K-12 / Nissle pilot)

Khoa (2026-10-26): "gather more data and run full pipeline". Goal: stress GIAP's assumptions on
bacteria that differ from E. coli in GC content, genome size, taxonomy and genetic code, with the
**full** pipeline (`giap run`, all four modules `core,mobile,bgc,function`, live Bakta,
`-profile farmshare,apptainer`, via `scripts/slurm/giap_run.sbatch`).

Nothing derived from these datasets at row level is committed (R11, Khoa 2026-10-23 rule for the
pilot applied here too): genomes, DIA-NN reports, converted peptide tables and run directories stay
under `$GIAP_ROOT/more_datasets/`. This directory holds only the plan, per-dataset summary reports
and small provenance summaries (accessions, sha256, exit codes, runtimes, counts).

## Selection criteria

1. Complete genome in NCBI RefSeq/GenBank (R6: official sources only).
2. Shotgun proteomics in PRIDE/ProteomeXchange **with a deposited search result that an existing
   GIAP converter reads** (`giap convert --from diann`: a DIA-NN `report.tsv` / `report.parquet`).
   No re-search of raw spectra: the deposited search is the evidence, as a user would bring it.
   Deposits with only Spectronaut, MaxQuant or FragPipe outputs were skipped (no converter).
3. Whole-cell proteome (not a pull-down), searched against a proteome of the same species.
4. Licence that allows use: all three are PRIDE CC0. NCBI genome records are public domain
   (NCBI policy; GenBank submitters' records carry no use restrictions).
5. Small enough to run in a few hours under the 32-CPU budget of this work unit.
6. Different from E. coli (50.8 % GC, 4.6-5.1 Mb, Gammaproteobacteria, table 11).

## Datasets

| Key | Organism | Genome (RefSeq / GenBank) | Size, GC, replicons | Proteomics | Search | Why |
|---|---|---|---|---|---|---|
| `msmeg` | *Mycolicibacterium smegmatis* mc²155 (Actinomycetota, Gram-positive) | GCF_000015005.1 / GCA_000015005.1 (ASM1500v1; CP000480.1 / NC_008596.1) | 6,988,209 bp, 67.5 % GC, 1 chromosome | PXD075873 (CC0), timsTOF Pro 2 diaPASEF; `Msm_report.tsv` 168 MB | DIA-NN 1.9, library-free, UniProt UP000000757 (built from GCA_000015005.1, so `search_db = same`), MBR on, no decoys reported | high GC (stresses the GC warning band, six-frame ORF density, promoter/terminator models trained on E. coli), 1.4x the K-12 genome, a different phylum |
| `mhyo` | *Mesomycoplasma hyorhinis* SK76 (Mycoplasmatota, Mollicutes; no cell wall) | GCF_000313635.1 / GCA_000313635.1 (ASM31363v1) | 836,897 bp, 26 % GC, 1 chromosome | PXD080153 (CC0), Orbitrap Fusion Lumos DIA; `2025NK081_ELFR_001to004_Table_results_report.parquet` 11 MB | DIA-NN 2.2, predicted library, "M. hyorhinis" proteome (`search_db = other`), MBR on | minimal genome, GC near the 0.25 warning edge, **translation table 4** (UGA = Trp): GIAP 1.0 assumes table 11 in Bakta, six-frame index, codon usage; the run shows what a user gets |
| `gsulf` | *Geobacter sulfurreducens* PCA (Thermodesulfobacteriota, Gram-negative) | GCF_000007985.2 / GCA_000007985.2 (ASM798v2; NCBI reference genome) | 3,814,128 bp, 61 % GC, 1 chromosome | PXD069850 (CC0), Q Exactive HF-X DIA; `report.tsv` 60 MB | DIA-NN 1.8.1, library-free, all *G. sulfurreducens* UniProt sequences (`search_db = other`) | Orbitrap (not timsTOF) DIA, older DIA-NN 1.8 report format, a non-Enterobacterales Gram-negative with >100 multi-haem cytochromes |

Candidates looked at and not taken: PXD065462 (*B. subtilis* 168; only `report.pg_matrix.tsv`
deposited), PXD067993 / PXD055092 (*S. aureus*; Spectronaut reports), PXD068569 (*M. agalactiae*;
co-IP, two small reports), PXD067466 (*Bacteroides*; searched against a non-reference proteome),
PXD053195 and other *B. subtilis* timsTOF sets (MaxQuant `combined.zip` only). No Firmicutes
(Bacillota) whole-proteome DIA-NN deposit was found in the PRIDE search; *M. hyorhinis*
(Mollicutes, a Bacillota-derived lineage) is the closest.

## Inputs and settings

- Download: `local/more_datasets/jobs/fetch.sbatch` (job `giap-fetch-more`): NCBI FTP
  `…_genomic.fna.gz` checked against NCBI's `md5checksums.txt`; PRIDE FTP files (PRIDE publishes
  no checksum; sha256 recorded on arrival in `$GIAP_ROOT/more_datasets/SHA256SUMS`).
- Conversion: `giap convert --from diann <report> --out peptides.tsv.gz --abundance-out
  abundance.tsv.gz --search-db same|other` (existing converter only).
- Run: `sbatch -J giap-run-more-<key> scripts/slurm/giap_run.sbatch --genome G --peptides P
  --db-dir $GIAP_DB_DIR --threads 3 --memory 12.GB --set annotation.bakta.complete=true --outdir
  $GIAP_ROOT/more_datasets/<key>/out`.
  - No `--annotation`: Bakta runs live (full database from `giap setup`), as for a new user.
  - `annotation.bakta.complete=true`: NCBI headers carry no topology tag; each genome is one
    complete circular chromosome.
  - R14 share of this work unit: 32 CPUs. `mhyo` ran with `--threads 3 --memory 12.GB` and the
    default `queueSize = 8` (at most 8 x 3 + head). After the coordinator's R14 alert (GIAP at
    20 jobs), the later runs cap Nextflow's queue with a launch-directory `OUT/nextflow.config`
    holding `executor { queueSize = 4 }`. `giap run` starts Nextflow with `cwd = OUT`, and
    `nextflow config` confirms that the override wins. They use `--threads 6 --memory 24.GB`, at
    most 4 x 6 + head = 28 CPUs and 5 jobs. Bakta full peaked at 11.8 GB on the 0.84 Mb genome,
    at the 12 GB clamp. Runs go one at a time.
  - Defaults otherwise (seed 20261014, `genome_relation = related`, all modules).

## Reports

`<key>/REPORT.md` per dataset: completion, exit codes, wall time per process (Nextflow trace),
failures, sanity of `genes.tsv`, mapping and FDR tier, operons and `operon_source`, context,
function, report rendering, and anything that looks wrong.
