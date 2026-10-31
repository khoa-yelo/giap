# Quick start

This page takes you from nothing to a finished GIPPA run, one command at a time. For each step it
gives the command, what it prints, how long it takes and how much disk it uses. If you have mass
spectra but no peptide or protein table yet, read
[From raw spectra to a GIPPA run](from-spectra.md) after step 5.

## The short version

```bash
pip install gippa                         # or: conda install -c bioconda gippa
gippa demo                                # about 2 minutes; prints the report path
gippa setup                               # databases and images for the detected runtime
gippa run --genome g.fna --peptides p.tsv --outdir out
```

`gippa demo` needs nothing else: no Nextflow, no container runtime, no databases and no network.
`gippa setup` and `gippa run` choose the Nextflow profile themselves, from what is installed on
the machine. They print what they found (see [What was detected](#what-was-detected)), and
`--profile` overrides it. The run writes `out/report/report.html`.

## 1. Install

```bash
pip install gippa                                     # Python >= 3.11
conda install -c conda-forge -c bioconda gippa        # or Bioconda: also brings Nextflow and Java
```

`pip` installs the `gippa` command and the bundled Nextflow workflow. For `gippa run` you also
need these (the conda package brings the first two):

- **Nextflow**, version 25.04.8 or newer (26.10.0 is the tested version);
- **Java 17** or newer, which Nextflow needs;
- **a container runtime**: Apptainer, SingularityCE or Docker.

[Install › Prerequisites](install.md#prerequisites) has the commands. Installing takes a minute or
two. The package itself is small; the demo data inside it are under 1 MiB.

Check the installation:

```bash
gippa --version
```

## 2. Set up the databases and images

```bash
gippa setup
```

With no flags, `gippa setup` installs everything the four default modules need. That is the
Bakta light database, the geNOMAD database, the antiSMASH databases, the Promotech models, GO,
KEGG, COG and the contaminant library. It also pulls the pinned container images for the runtime
it detects. It writes them to one database directory:

- `$GIPPA_DB_DIR`, if you set it;
- else `$XDG_DATA_HOME/gippa/db`, which is usually `~/.local/share/gippa/db`.

Setup records the directory in `~/.config/gippa/db_dir`. `gippa run`, `gippa preflight` and
`gippa makedb` then find it without `--db-dir`, and they say so:
`gippa run: DB directory /x (from gippa setup)`.

**What it prints.** Before it downloads anything, setup prints a plan. The plan lists every item
with its size, the total, and the free space on each filesystem it writes to. For a new
Apptainer setup it starts like this (shortened):

```text
gippa setup: profile apptainer (auto: apptainer found; pass --profile to choose)
gippa setup: plan: bakta_light 5.0 GiB
gippa setup: plan: genomad 4.0 GiB
gippa setup: plan: antismash 10.0 GiB
gippa setup: plan: promotech_models 5.0 GiB
…
gippa setup: plan: image antismash 4.6 GiB
…
gippa setup: plan: temporary space in /tmp (APPTAINER_TMPDIR) for the largest image build 4.6 GiB
gippa setup: plan: total … GiB
gippa setup: free: … GiB on /home/u/.local/share/gippa/db (DB directory)
gippa setup: Proceed? [Y/n]
```

The plan also prints two licence lines, one for each restricted item it will install:

```text
geNOMAD (module mobile) is for non-commercial use only; for commercial use: gippa setup --no-genomad, gippa run --modules core,bgc,function
KEGG is free for academic use only; for other use: gippa setup --no-kegg
```

Every item has an opt-out flag, such as `--no-genomad`, and `--no-images` skips the images. See
[Set up the databases](setup-databases.md) for the flags, and [Licences](../user-guide/licences.md)
for the terms of each tool and database.

**UniOP.** UniOP, the operon predictor, has no licence file. GIPPA never distributes it: setup
fetches it from its authors and builds it on your machine. It does so only after you agree, either
with the prompt `Build UniOP for your own use? [y/N]` or with `--accept-uniop-terms`. Without your
consent, runs use the distance-based operon model instead.

**Disk.** The plan's sizes are the free space setup asks for before it downloads. On FarmShare,
the installed default databases take about 19.3 GB, and the images another 6.0 GB:

| Item | On disk |
|---|---|
| Bakta light | 4.0 GB |
| geNOMAD | 1.4 GB |
| antiSMASH | 9.4 GB |
| Promotech | 4.5 GB |
| GO, KEGG, COG and contaminants | 34 MB |
| the six pinned Apptainer images | 6.0 GB |

The Bakta full database (`--bakta full`) adds about 84 GB.

**Time.** Setup time depends on your network. Most of it is spent downloading about 20 GB.

**Downloading by hand.** On a host where setup cannot reach the sources, see
[Download the databases by hand](manual-databases.md).

### What was detected

`gippa setup`, `gippa preflight` and `gippa run` print the profile they chose, and why:

```text
gippa run: profile slurm,apptainer (auto: sbatch on PATH; pass --profile to choose)
```

| Your machine | Profile chosen |
|---|---|
| Stanford FarmShare | `farmshare,apptainer` |
| a Slurm cluster with Apptainer | `slurm,apptainer` |
| a workstation with only Apptainer | `apptainer` |
| a laptop with a working Docker | `docker` |
| `sbatch` and conda (micromamba) only | `slurm,conda` |

If the images were pulled for one runtime, a later run prefers that runtime. To choose yourself,
pass `--profile`, for example `gippa run --profile docker …`. On a cluster login node, setup also
prints this hint: `the Nextflow head process runs on this host for the whole run; submit gippa run
as a job`.

## 3. Check that it runs: `gippa demo`

```bash
gippa demo
```

The demo runs GIPPA on SIM100K, a small synthetic genome bundled in the package. SIM100K is a
100 kb circular chromosome and an 8 kb circular plasmid, with a simulated peptide table. The
demo runs the modules `core`, `mobile` and `bgc`. It does not use the real annotation tools.
Instead it reads precomputed Bakta, promoter, terminator, RBS, operon and region layers, which
GIPPA's test simulator (gippasim) wrote. Every step runs inside this one Python process, with no
Nextflow, no containers, no databases and no network.

**What it prints.** The last line names the report, for example:

```text
gippa demo: done in 48 s; report: /home/u/gippa-demo/report/report.html
```

**Time.** About 2 minutes. The target is 2 minutes on 2 CPUs; one test run took 121 s on a busy
8-CPU host. **Disk.** The output goes to `./gippa-demo` (use `--outdir` to put it elsewhere).
`gippa demo --open` also opens the report in a browser.

**What a good result looks like.** The command exits with code 0, and the report opens with the
banner "Synthetic demo data (SIM100K, gippasim seed 20261014)". `run_manifest.json` has
`"status": "completed"` and `"demo": true`. The numbers describe a simulated genome, so do not
interpret them as biology.

Once setup has finished, check the whole installation (preflight, images, databases, Nextflow)
with the same inputs and the real tools:

```bash
gippa demo --full --outdir gippa-demo-full
```

This is a real run, so it takes as long as step 4a.

## 4. Example runs

### 4a. A genome and a peptide table

The demo inputs are in the package. Copy them out:

```bash
DEMO=$(python -c "from importlib.resources import files; print(files('gippa.data') / 'demo')")
cp "$DEMO/sim100k.fna.gz" "$DEMO/sim100k_peptides.tsv.gz" .
```

From a checkout, the same files are `tests/data/genomes/sim100k.fna.gz` and
`docs/tutorials/data/sim100k_peptides.tsv.gz`.

**The genome** (`sim100k.fna.gz`, 34 KB) has two records. The `[topology=circular]` tag marks
each one circular:

```text
>chr1 [topology=circular]
TACACCGGGTTACTCCACAGGATAATATCTTACTGATAGCCATTCTATTCTATAATATAA…
>plasmid_1 [topology=circular]
…
```

**The peptide table** (`sim100k_peptides.tsv.gz`, 68 KB) holds 1,685 rows of one run: 1,668
target rows and 17 decoy rows. It is in the generic format that `gippa convert` writes
([Inputs › Peptide table](../user-guide/inputs.md#peptide-table-peptidestsv)). The first rows,
and one decoy:

```text
peptide                      protein_ids     q_value          pep              score     is_decoy  intensity  run   charge  spectrum_id  source    search_db
DAAHIFVNEARGAR               GSIM_00480      8.222418025e-05  0.0004403747136  2.550507  false     7.49312    run1  3       run1.1       database  same
GGPVGHVLWLAVVMHLAGFGTEPARLR  GSIM_00455      5.90062343e-06   1.864430642e-05  3.604676  false     2.1869     run1  3       run1.2       database  same
DLLFGDCFNIK                  GSIM_00200      1.016223144e-06  2.61801616e-06   4.259056  false     68.8833    run1  2       run1.3       database  same
RSIDQKNDANLTDTTHFIGVWML      rev_GSIM_00120  NA               NA               NA        true      NA         run1  3       run1.41      database  same
```

Run it:

```bash
gippa preflight --genome sim100k.fna.gz --peptides sim100k_peptides.tsv.gz   # seconds
gippa run --genome sim100k.fna.gz --peptides sim100k_peptides.tsv.gz --outdir sim100k_run
```

Preflight checks the inputs and the installed tools, and exits with code 0 when the run can
start. On the demo it gives one warning, `genome.contiguity`: a 100 kb chromosome is much smaller
than a bacterial genome.

**Time and disk.** On FarmShare, this run took 37 minutes (Slurm job 1894327, profile
`farmshare,apptainer`, 43 tasks). It wrote 20 MB, including 12 MB in Nextflow's `work/`
directory. The report alone is 5.8 MB. If a run stops, rerun the same command with `--resume`.

### 4b. A genome and a protein list

If your search engine gives you proteins rather than peptides, pass a protein FASTA with
`--proteins`. To try this, make a protein FASTA from the demo's Bakta layer with
`gippa makedb`, then run on it:

```bash
gippa makedb --annotation "$DEMO/layers/bakta.json.gz" --no-contaminants --outdir sim100k_db
gippa run --genome sim100k.fna.gz --proteins sim100k_db/strain_db.faa --outdir sim100k_proteins
```

`makedb` takes about a second and prints
`gippa makedb: wrote sim100k_db/strain_db.faa {'cds': 82, 'sorf': 3, 'patch': 0, 'orf': 0, 'contaminant': 0}`.
The FASTA is 29 KB:

```text
>GSIM_00005 hypothetical protein
MNATRTGRAMASLEAGDAPEKVLAADTEGANVHLDLAALLMQPQIAEKAIFRSQKVPNQI
…
```

GIPPA aligns proteins longer than 40 aa to the genome with miniprot, then tblastn. A protein list
gives the `protein_list` detection tier: a gene is detected when an input protein maps to it, and
its error is `unestimated`. This example lists every protein of the genome, so every mapped gene
is detected. The detection and region contrasts then have nothing to compare against. Use it to
learn the mechanics, not the statistics.

### A tour of the outputs

```text
sim100k_run/
├── run_manifest.json      tool and database versions, parameters, input SHA-256, seeds, timings
├── resolved_config.yaml   the configuration the run used
├── preflight/             preflight.json, preflight.md
├── annotation/            genes.gff3, regions.tsv, regions.bed, promoters.bed, terminators.bed,
│                          operon_pairs.tsv, tus.tsv
├── mapping/               peptide_hits.tsv, protein_hits.tsv, unmapped.tsv, mapping_fdr.tsv, …
├── tables/                genes.tsv (the master table), detection_summary.tsv,
│                          regions_stats.tsv, operons.tsv, variants.tsv, novel_orfs.tsv, …
├── plots/                 Vnn_<name>.svg|png|pdf, data/Vnn_<name>.tsv
├── report/report.html
└── logs/                  events.jsonl, trace.tsv, timeline.html, nf_report.html
```

The rows below come from the committed SIM100K golden tables
(`tests/data/golden/sim100k_real/m4/`). They were computed from the same genome and layers, but
with a test configuration (gene-level alpha 0.05) and gippasim's own peptides, so your numbers
will differ.

**`tables/genes.tsv`** has one row per gene, keyed by locus tag, with 73 columns. A few of them:

```text
locus_tag   contig     start  end    strand  region_ids  detected  detection_tier  detection_status  gene_q_value  n_peptides_unique  abundance_log2  flags
GSIM_00005  chr1       288    2231   +       NA          true      mayu            detected          1.66005e-08   21                 -1.05081        detectability_constant
GSIM_00065  chr1       13194  14465  -       NA          false     mayu            not_detected      NA            0                  NA              pseudogene;detectability_constant
GSIM_00120  chr1       25534  26205  +       bgc_1       true      mayu            detected          7.08563e-08   17                 2.31231         detectability_constant
GSIM_00505  plasmid_1  193    1854   -       plasmid_1   true      mayu            detected          2.57978e-08   54                 4.20376         detectability_constant
```

Coordinates in the tables are 1-based and closed. `detection_tier` says how "detected" was
decided ([Interpreting results › Detection tiers](../user-guide/interpreting-results.md#detection-tiers)).

**`annotation/regions.tsv`** lists the genomic regions: prophages (geNOMAD), biosynthetic gene
clusters (antiSMASH), CRISPR arrays and plasmids. SIM100K plants one of each. The demo's
bundled regions layer:

```text
region_id   region_type  contig     start  end    source_tool    score   genes
bgc_1       bgc          chr1       25380  35379  antiSMASH      NA      GSIM_00120;GSIM_00125;…;GSIM_00160
prophage_1  prophage     chr1       36616  48615  geNOMAD        0.7761  GSIM_00185;GSIM_00190;…;GSIM_00235
crispr_1    crispr       chr1       62388  62799  Bakta          NA      NA
plasmid_1   plasmid      plasmid_1  1      8000   Bakta+geNOMAD  0.9074  GSIM_00505;GSIM_00510
```

**`tables/regions_stats.tsv`** tests each region: are its genes detected more or less often than
their detectability predicts (`o_minus_e`), and are they more or less abundant (`hl_shift`)?
A few columns, rounded:

```text
region_id   subset  n_genes  observed  expected  o_minus_e  stat       p_value  q_value  method
prophage_1  all     11       11        10.747    0.253      o_minus_e  1        1        circular_shift_o_minus_e
bgc_1       all     9        9         8.793     0.207      hl_shift   0.2235   0.2235   circular_shift_hodges_lehmann
plasmid_1   all     2        NA        NA        NA         o_minus_e  NA       NA       replicon_contrast
```

**`report/report.html`** is one self-contained file that loads nothing from the network. Its
sections are: Summary, Inputs and quality control, Peptide mapping, Gene detection, Genome
context, Operons, Variants, Novel ORFs, Function, Determinants of abundance, Sensitivity, Methods
and Provenance. On a cluster, copy the whole output directory to your computer, so that the
report's links to the full tables keep working.

[Outputs](../user-guide/outputs.md) describes every file, and the
[table reference](../reference/tables.md) describes every column.

## 5. Your own data: which route?

Which route you take depends on what you have:

| You have | Route |
|---|---|
| a peptide table from a search engine | `gippa run --genome genome.fna --peptides peptides.tsv --outdir out` |
| DIA-NN output (`report.parquet`) | `gippa convert --from diann`, then `--peptides` ([Converting](../user-guide/convert.md)) |
| Cascadia de novo output (`.ssl`) | `gippa convert --from cascadia`, then `--peptides`, ideally with `--reference-peptides` |
| a protein list (FASTA) | `gippa run --genome genome.fna --proteins proteins.faa --outdir out` |
| an existing Bakta annotation | add `--annotation genome.bakta.json`; Bakta is then not run |
| raw spectra and no table yet | [From raw spectra to a GIPPA run](from-spectra.md) |

The genome must be one assembled bacterial genome in FASTA format. Add ` [topology=circular]` to
the header of each circular record, or pass `--set annotation.bakta.complete=true` for a complete
assembly ([Inputs › Genome FASTA](../user-guide/inputs.md#genome-fasta)). Keep the search
engine's decoy rows: with decoys, GIPPA can estimate gene-level FDR (the `decoy_picked` tier).

For example, the E. coli K-12 MG1655 assembly with the UniProt reference proteome as a protein
list:

```bash
curl -fLO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz
curl -fL -o UP000000625.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29"
gippa run --genome GCA_000005845.2_ASM584v2_genomic.fna.gz --proteins UP000000625.fasta \
  --set annotation.bakta.complete=true --outdir out_k12
```

On a cluster, submit `gippa run` as a batch job: the Nextflow head process runs for the whole
run. The run takes hours, not minutes. On FarmShare, an E. coli Nissle genome with a DIA-NN
peptide table took between 1.6 and 4.5 hours on 14 CPUs and 54 GB
([From raw spectra › Run times](from-spectra.md#run-times-and-sizes-from-the-pilot)).

## If this fails

!!! failure "No container runtime found"
    `gippa run` exits with code 2 and prints `no container runtime or conda found: install
    Apptainer (no root needed on most Linux systems) or Docker`. Install Apptainer (or
    SingularityCE, or Docker with a running daemon you can reach), then run `gippa setup` again,
    so that it pulls the images for that runtime. `gippa demo` still works without one. If setup
    ran on a host with no runtime, such as a transfer node, it downloaded the databases and
    skipped the images. Rerun `gippa setup --images --profile apptainer` on a host that has
    one, with the same database directory.

!!! failure "Java older than 17"
    Nextflow needs Java 17 or newer. Check with `java -version`. To install a newer Java without
    root, for example:
    `micromamba create -n gippa-java -c conda-forge openjdk=21`. The Bioconda package brings
    Java with it.

!!! failure "Not enough disk"
    Setup stops before it downloads anything, with exit code 2 and a message such as
    `needs 31.4 GiB, 12.0 GiB free on /home/u: pass --db-dir (or set GIPPA_DB_DIR) on a larger
    disk`. Point `--db-dir` at a project or scratch filesystem; clusters often give home
    directories a small quota, and setup warns when the default directory is there. Or install
    less: `gippa setup --minimal` installs only what the `core` and `function` modules need
    (no geNOMAD, no antiSMASH). Later runs on that directory then use `--modules core,function`.

## Next

- [From raw spectra to a GIPPA run](from-spectra.md): DIA-NN, six-frame ORFs and de novo.
- [Inputs](../user-guide/inputs.md): peptide tables, protein lists, abundance and annotations.
- [Running GIPPA](../user-guide/run.md): modules, seeds, `--publication` and exit codes.
- [Interpreting results](../user-guide/interpreting-results.md): what the q-values mean.
