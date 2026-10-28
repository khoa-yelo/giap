# From raw spectra to a GIPPA run

GIPPA does not search spectra. It reads the peptide table that a search engine or a de novo
sequencer produced. This page is for users who have mass spectra but no peptide or protein table
yet. It describes four routes from spectra to a `gippa run`, says when each one fits, and gives
the commands. In every route the steps are the same:

1. build a search database with `gippa makedb` (not needed for de novo);
2. run the search engine yourself;
3. turn its output into GIPPA's peptide table with `gippa convert`;
4. run `gippa run`.

The worked example throughout is GIPPA's pilot: three wild-type E. coli Nissle 1917 runs from
**PXD078459** (timsTOF HT, diaPASEF). Every time and size below was recorded in the pilot kit
(`validation/pilot/README.md` and `validation/pilot/reproduce.md` in the repository). Where the
pilot did not run a route, this page says so.

## Choose a route

| Route | Input spectra | Typical use | `gippa convert --from` | What it feeds in GIPPA | Pilot run time |
|---|---|---|---|---|---|
| [A. DIA-NN against the called genes](#a-dia-nn-against-the-called-genes) | DIA | the default for an annotated genome | `diann` | every analysis: gene-level FDR (`decoy_picked`), abundance, regions, islands, operons, function, variants | search 14–42 min; GIPPA 1.6–4.5 h |
| [B. DIA-NN against six-frame ORFs](#b-dia-nn-against-six-frame-orfs) | DIA | finding unannotated genes (novel ORFs) | `diann` | as A, plus novel-ORF class FDR | not run in the pilot |
| [C. De novo with Cascadia](#c-de-novo-with-cascadia) | DIA | no usable database: a poorly annotated or divergent strain | `cascadia` | detection (`two_peptide` tier), novel ORFs, variants; no abundance | GIPPA 3.0 h (de novo step not timed) |
| [D. Any other engine, including DDA](#d-other-engines-including-dda) | DDA or DIA | you already use another engine | `generic` (you write the table) | as A or C, depending on the columns you provide | not run in the pilot |

All four routes need the genome and, for routes A and B, its Bakta annotation. They feed the
`core` module, so they work with any `--modules` choice. The other modules (`mobile`, `bgc` and
`function`) work from the genome and do not depend on the route.

## Before you start: the annotation

`gippa makedb` builds the database from a **Bakta JSON**. `gippa makedb` does not run Bakta
itself (without `--annotation` it exits with code 2), so annotate the genome first. Run Bakta
1.12.2, the version GIPPA pins, with the database that `gippa setup` installed:

```bash
bakta --db ~/.local/share/gippa/db/bakta/db-light --output bakta --prefix strain \
  --keep-contig-headers genome.fna
```

Bakta is GPL-3.0 software. Install it from Bioconda, or use the pinned image
`oschwengers/bakta:v1.12.2`. Pass the same JSON to `gippa run --annotation`, so that the protein IDs
in your search results are the locus tags of the run. Leave `strain.gbff` next to `strain.json`:
without that GBFF file, the `bgc` module skips antiSMASH. On the pilot genome, Bakta with the
full database took 1 h 43 min on 18 CPUs.

## A. DIA-NN against the called genes

**When to use it.** You have DIA data and an assembled genome of the strain you measured. This is
the default route: the smallest search space, the best-calibrated FDR, and abundances from
DIA-NN's intensities. It cannot find genes that Bakta did not call: a peptide from an unannotated
ORF has no protein in the database to match.

**Database.** The Bakta CDS and sORF proteins, plus the contaminant library:

```bash
gippa makedb --annotation bakta/strain.json --outdir strain_db
```

The contaminants come from the database directory of `gippa setup`, and are added as
`>CONT_<id>` records. `--no-contaminants` leaves them out. The output is
`strain_db/strain_db.faa`, `db_id_map.tsv` and `makedb_manifest.json`
([Building a search database](../user-guide/makedb.md)).

**Search.** Run DIA-NN in two calls. DIA-NN 2.x turns the FASTA digest off when
`--fasta-search --predictor` and raw files are given in one call. So the first call predicts a
spectral library from the FASTA alone, and the second searches the runs against it. These are the
pilot's settings (`validation/pilot/diann_A2.cfg`). The digest flags equal GIPPA's own digest:
trypsin, one missed cleavage, peptides of 7 to 30 residues, fixed carbamidomethyl C.

```bash
# 1. predict the library from the FASTA (no raw files)
diann --fasta-search --gen-spec-lib --predictor --fasta strain_db/strain_db.faa \
  --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --unimod4 \
  --var-mods 1 --var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4 \
  --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --peptidoforms \
  --threads 16 --out-lib lib/strain.parquet
# 2. search the runs against the predicted library
diann --lib lib/strain.predicted.speclib --f run1.d --f run2.d --f run3.d \
  --fasta strain_db/strain_db.faa --gen-spec-lib --out-lib lib/strain_empirical.parquet \
  --qvalue 0.01 --report-decoys --matrices --peptidoforms --threads 16 \
  --out strain_report.parquet
```

- Keep `--report-decoys`. GIPPA needs the decoy rows for gene-level FDR (the `decoy_picked`
  tier).
- Do not add `--met-excision`. GIPPA's digest keeps the initiator Met, so DIA-NN would report
  peptides that GIPPA never proposes.
- The pilot also left out match-between-runs (`--reanalyse`) and `--rt-profiling`. MBR weakens
  run-level FDR control.

**Convert and run.**

```bash
gippa convert --from diann strain_report.parquet --out peptides.tsv.gz \
  --abundance-out abundance.tsv.gz --search-db same
gippa run --genome genome.fna --annotation bakta/strain.json --peptides peptides.tsv.gz \
  --outdir out
```

`--search-db same` says that the database was built from this genome. GIPPA then accepts exact
hits without a placement decoy search (the `same_db` flag). DIA-NN peptides carry intensities, so
GIPPA computes gene abundance from them; the protein table from `--abundance-out` is needed only
for peptides without intensities.

**A related reference proteome instead.** If you search a reference proteome of a related strain
(pilot arm A1: UniProt UP000000625, E. coli K-12, against the Nissle genome), pass
`--search-db other`. GIPPA then places every peptide on the genome with its own tolerant mapper
and placement FDR. Peptides that differ between the strains are scored as tolerant hits, or not
placed at all.

**Pilot (arm A2).**

| Step | Recorded |
|---|---|
| `gippa makedb` | 4,685 proteins (4,613 CDS and 72 sORFs), 989 KB gzipped; the job took 23 s |
| DIA-NN 2.3.1, both calls, 3 runs | 42 min on 16 CPUs (pass 3). In pass 1, 14 min with 17 threads, peak memory 6.4 GiB |
| `gippa convert` | 17–28 s, 0.3 GiB; 94,363 peptide rows, 2.1 MB gzipped; 3,250 protein groups |
| `gippa run` (all modules, 14 CPUs, 54 GB) | 1 h 44 min (pass 1) to 4 h 00 min (pass 3); TerminatorNet took most of the time; report 20.7 MB |

Arm A1 (the K-12 proteome) gave 81,278 peptide rows. Its runs took 1 h 37 min (pass 1) and
4 h 30 min (pass 3).

**Caveats.**

- **FDR.** DIA-NN's q-values control precursor FDR, not gene FDR. GIPPA estimates gene-level FDR
  itself, from the decoys. See
  [Placement FDR is not search-engine FDR](../user-guide/interpreting-results.md#placement-fdr-is-not-search-engine-fdr).
- **Search space.** It is only the called genes. A gene missing from the annotation can never be
  detected.
- **Sensitivity.** This is the highest of the four routes: a small database and spectral
  library-based scoring.

## B. DIA-NN against six-frame ORFs

**When to use it.** You want to find protein-coding genes that the annotation missed: short ORFs,
genes on the opposite strand, and N-terminal extensions. This route feeds GIPPA's novel-ORF
analysis (`tables/novel_orfs.tsv`).

**Database.** Add the six-frame ORFs to the database of route A:

```bash
gippa makedb --annotation bakta/strain.json --six-frame-orfs --min-orf-aa 20 --outdir strain_db_orfs
```

Each ORF of at least `--min-orf-aa` residues that is not in frame with an annotated CDS gets a
record `>GIPPAORF_<class>_<contig>_<start>_<end>_<p|m>`. The classes are `intergenic`,
`cds_outframe`, `antisense` and `nterm_extension`. On the SIM100K demo genome, the database grows
from 85 proteins (29 KB) to 85 proteins and 2,941 ORFs (272 KB).

**Search, convert and run.** As in route A, with `strain_db_orfs/strain_db.faa` as the FASTA, and
with the same `--search-db same`.

**Caveats.**

- **FDR.** GIPPA estimates FDR for novel ORFs per ORF class, from decoy rows whose protein IDs are
  the `GIPPAORF_…` IDs of this database. Keep `--report-decoys`. Without such decoys, the run warns
  `novel_orfs.decoys_not_makedb`, the class FDR is `unestimated`, and only the at-least-2-peptides
  rule applies. A class with fewer than about 1/alpha candidates (100 at alpha 0.01) cannot be
  called at all (warning `novel_orfs.class_q_floor_above_alpha`).
- **Search space.** The database is many times larger, as the SIM100K example shows. A larger
  search space makes DIA-NN's FDR threshold stricter, so fewer peptides of annotated genes pass.
  Library prediction and the search also take longer.
- **Sensitivity.** Expect fewer identifications of annotated genes than with route A on the same
  spectra. A useful pattern is to run route A for the gene table and route B for novel-ORF
  discovery, and to compare the two runs with `gippa compare`.
- **Not run in the pilot.** The pilot databases had no six-frame ORFs, so there is no recorded run
  time or size for this route on real data.

## C. De novo with Cascadia

**When to use it.** There is no good protein database for the strain: the genome is poorly
annotated, or divergent from every reference. A de novo sequencer reads peptide sequences
directly from the spectra. GIPPA then places them on the genome with its own tolerant mapper,
which allows mismatches. This makes the route useful for variant evidence too. GIPPA reads the
output of **Cascadia** 0.0.7, a de novo sequencer for DIA data.

**Commands.** No `makedb` step is needed. Cascadia reads mzML. For timsTOF `.d` data, the pilot
converted each run with TIMSCONVERT 2.0.0 first:

```bash
timsconvert --input run1.d --outdir mzml --mode centroid --compression zlib --exclude_mobility
cascadia sequence mzml/run1.mzML cascadia.ckpt -o denovo/run1 -b 32
```

Cascadia writes `denovo/run1.ssl`, and keeps predictions above its default score cut-off of 0.8.
It needs a GPU (the pilot also had a CPU route). Then convert all runs together and run GIPPA,
ideally with a database search of the same spectra as reference peptides:

```bash
gippa convert --from cascadia denovo/run1.ssl denovo/run2.ssl denovo/run3.ssl --out denovo.tsv.gz
gippa run --genome genome.fna --annotation bakta/strain.json --peptides denovo.tsv.gz \
  --reference-peptides peptides.tsv.gz --outdir out_denovo
```

`--reference-peptides` is a peptide table from a database search (route A) of the same spectra.
GIPPA uses it to estimate the precision of the de novo peptides, and the error rates of variant
calls. Without it, GIPPA assumes a precision and flags the run `c_assumed`.

**Pilot (arm A3).**

| Step | Recorded |
|---|---|
| timsTOF `.d` to mzML | 1.44–1.47 GB of mzML per run |
| Cascadia | about 1,400 rows per run (`.ssl` of about 215 KB); the pilot used `.ssl` files from a CPU route, whose wall time is not recorded |
| `gippa convert --from cascadia` | 25 s on 2 CPUs; 4,234 rows, 1,386 distinct peptides, 63 KB gzipped |
| `gippa run` (all modules, 14 CPUs, 54 GB) | 3 h 01 min, 77 tasks; report 11.1 MB |

**Caveats.**

- **FDR.** De novo rows carry no decoys and no q-values. Gene detection therefore uses the
  `two_peptide` tier: at least 2 distinct accepted peptides, with the gene error `unestimated`.
  GIPPA's placement FDR uses shuffled decoy frames for de novo input.
- **Search space.** There is none: any sequence can be read. The genome alone constrains which
  peptides are placed.
- **Sensitivity.** Much lower than a database search: the pilot got about 1,400 de novo rows per
  run, against about 94,000 peptide rows from DIA-NN on the same three runs. Cascadia gives no
  intensities, so the abundance analyses (GSEA, islands, determinants of abundance) do not run.
  The manifest lists `abundance` under `layers_missing`.
- **Other de novo tools.** `gippa convert --from casanovo` and `--from instanovo` are accepted
  choices, but not implemented: both exit with code 5. Until a converter exists, write the
  generic table yourself (route D) with `--source denovo`.

## D. Other engines, including DDA

**When to use it.** You use a search engine other than DIA-NN, for example on DDA data.
`gippa convert` has no converter for a DDA engine yet. Write the generic peptide table yourself,
then let `gippa convert --from generic` check it.

The table is a TSV with one header row. Only `peptide` is required, but each column you add
unlocks more of GIPPA. The columns are listed in
[Inputs › Peptide table](../user-guide/inputs.md#peptide-table-peptidestsv).

| Columns | What they give you |
|---|---|
| `peptide` | detection by the `two_peptide` tier |
| `q_value` or `pep` | the `mayu` tier: gene-level FDR without decoys |
| `is_decoy` and `protein_ids` on the decoy rows, plus a ranking column (`q_value`, `pep` or `score`) shared with the targets | the `decoy_picked` tier, the strongest |
| `intensity` and `run` | gene abundance and the abundance analyses |
| `spectrum_id` | spectrum pairing for variant error rates |

```bash
gippa convert --from generic my_engine_peptides.tsv --out peptides.tsv --source database --search-db same
gippa run --genome genome.fna --annotation bakta/strain.json --peptides peptides.tsv --outdir out
```

The generic converter strips modification syntax such as `(UniMod:35)` or `[+15.995]`, validates
every column and drops unknown ones with a warning. `--source` and `--search-db` fill those columns
where they are empty. Build the database with `gippa makedb` as in route A or B, so that the
protein IDs are the run's locus tags and `GIPPAORF_…` IDs.

**Caveats.** As for route A or B, depending on the database. Search-engine q-values from DDA
engines differ in level (PSM, precursor or peptide). GIPPA reads `q_value` as given, and
estimates gene-level FDR itself. Not run in the pilot.

## Run times and sizes from the pilot

| Step | Arm | Recorded | Where |
|---|---|---|---|
| download, 3 runs (about 1.0 GB each, zipped) | all | 23 min, 2.9 GB | Slurm 1879859 |
| Bakta 1.12.2, full database | all | 1 h 43 min, 18 CPUs | Slurm 1895437 |
| `gippa makedb` | A2 | 23 s (job total) | Slurm 1895583 |
| DIA-NN 2.3.1, 3 runs | A2 | 42 min, 16 CPUs (pass 3); 14 min, 17 threads, 6.4 GiB (pass 1) | Slurm 1895586, 1879956 |
| `gippa convert --from diann` | A1, A2 | 17–28 s, 0.3 GiB | Slurm 1880219, 1880220, 1895670 |
| `gippa convert --from cascadia` | A3 | 25 s, 2 CPUs | Slurm 1894680 |
| `gippa run` | A1 | 1 h 37 min (pass 1), 4 h 30 min (pass 3) | Slurm 1892719, 1896178 |
| `gippa run` | A2 | 1 h 44 min (pass 1), 4 h 00 min (pass 3) | Slurm 1893676, 1895673 |
| `gippa run` | A3 | 3 h 01 min | Slurm 1894682 |

Every `gippa run` used 14 CPUs and 54 GB, all modules, and the default seed. On FarmShare the
pilot ran these steps with `scripts/slurm/diann_pilot.sbatch` (DIA-NN),
`scripts/slurm/mzml_pilot.sbatch` and `scripts/slurm/denovo_pilot.sbatch` (de novo),
`scripts/slurm/diann_convert.sbatch` and `scripts/slurm/cascadia_convert.sbatch` (conversion),
and `scripts/slurm/gippa_run.sbatch`. The step-by-step runbook is `validation/pilot/reproduce.md`.
These scripts carry FarmShare paths, so copy them as a template, not as is.

## Licences of the search tools

GIPPA bundles none of these tools. It never runs them for you, and it reads only their output
files. You install each one and accept its terms yourself.

| Tool | Licence | Notes |
|---|---|---|
| DIA-NN | proprietary, its own terms (academic licence) | the pilot used DIA-NN 2.3.1 Academia; GIPPA never redistributes DIA-NN binaries |
| Cascadia | Apache-2.0 (code); the weights are under the terms given in the paper | the pilot pinned the 0.0.7 wheel and the weights by sha256 |
| TIMSCONVERT, pyTDFSDK | Apache-2.0; the Bruker TDF-SDK they bundle is under Bruker's licence | needed only for timsTOF `.d` data |
| Bakta | GPL-3.0 | runs as a separate program, also inside `gippa run` |

GIPPA itself is MIT-licensed. The licences of the tools and databases that `gippa run` uses are on
[Licences](../user-guide/licences.md).
