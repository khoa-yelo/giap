# Changelog

All notable changes are recorded here by towncrier from `changelog.d/` fragments.

<!-- towncrier release notes start -->

## GIPPA 1.0.0 (2026-10-29)

### Features

- The genome-wide TerminatorNet scan is split into chunks of 250,000 windows per contig strand
  (`TERMINATORNET_PLAN` → `TERMINATORNET_CHUNK` → `TERMINATORNET_MERGE`), each on
  `--terminatornet_chunk_cpus` CPUs (default 4), with output byte-identical to the unsharded
  `-n ≥ 2` run. A timeout now re-runs one chunk instead of the whole scan; `--terminatornet_sharded
  false` runs the old single task. On FarmShare a chunk asks for 4 GB and 2 h per attempt, and the step
  uses at most 32 CPUs at a time. ([#455](https://github.com/khoa-yelo/gippa/pull/455))
- Preflight warns (`inputs.bakta_gbff`) when a precomputed Bakta JSON (`--annotation x.json[.gz]` or `annotation.precomputed.bakta_json`) has no `x.gbff` or `x.gbff.gz` next to it while the bgc module is on: antiSMASH would be skipped and the bgc layer would be empty. The warning names the expected paths, appears in `preflight.json`, `preflight.md` and the `run_manifest.json` warnings, and `--strict` turns it into exit code 3.
  `preflight.md` now shows check messages in full (the Value column keeps its 200-character cap). ([#467](https://github.com/khoa-yelo/gippa/pull/467))
- `gippa run` now publishes the Bakta JSON and GBFF it produced as `OUT/annotation/tools/bakta/genome.json` and `genome.gbff` (copied, untouched), so `gippa makedb --annotation OUT/annotation/tools/bakta/genome.json` works after a run without running Bakta by hand. A precomputed Bakta JSON is not copied into `OUT/`. Runs started with an earlier version re-run Bakta and the steps after it on `--resume`. ([#468](https://github.com/khoa-yelo/gippa/pull/468))
- Under `--profile docker`, Nextflow's tool containers no longer pull images at run time (`docker run --pull=never`), and preflight stops the run when a pinned image the selected modules need is not present locally, naming `gippa setup --images --profile docker`. It also reports a Docker client older than 20.10 (podman older than 3.0) and a Docker daemon that is not reachable. Docker 20.10 or later is required. ([#477](https://github.com/khoa-yelo/gippa/pull/477))
- <!-- m00 -->
  GIPPA installs as a Python package with a `gippa` command line, a validated YAML configuration
  and a Nextflow workflow. It ships with a built-in genome and proteomics simulator for testing, and
  continuous integration covers both.
- <!-- m01 -->
  Inputs and preflight: `gippa preflight` checks the genome, proteins, peptides, abundance table and
  environment before a run and writes `preflight.json` (validated against its JSON Schema) and
  `preflight.md`, with exit codes 0, 1, 2 and 3. `gippa convert --from generic|diann` turns a search
  engine report into GIPPA's peptide table. Readers and writers cover FASTA, GFF3, BED, proBed, VCF 4.3
  and TSV.
- <!-- m02 -->
  Genome annotation: `gippa setup` downloads and checks the Bakta (light or full), geNOMAD, antiSMASH,
  Promotech, GO, KEGG, COG and contaminant databases and records them in `db_manifest.json`. A run
  annotates the genome with Bakta and adds operon pairs and transcription units (UniOP), promoters,
  terminators, ribosome binding sites, codon usage and regions (prophage, plasmid, biosynthetic gene
  cluster, CRISPR, AMR). A Bakta JSON or your own GFF3/GBFF annotation can be supplied instead.
  The pipeline runs through Nextflow (`-profile farmshare,apptainer`) or a local executor.
- <!-- m03 -->
  Peptide and protein mapping: peptides are placed on the genome through a six-frame index (exact,
  tolerant and frameshift matches) and proteins through miniprot and tblastn, with target-decoy
  control of the placement FDR. Unmapped and contaminant evidence is reported with reasons.
  `gippa makedb` writes a strain search database (annotated proteins plus six-frame ORFs) for your
  spectrum search.
- <!-- m04 -->
  Gene table and region tests: `genes.tsv` has one row per Bakta locus tag with detection tiers
  (picked FDR, MAYU, two-peptide, protein list), a cross-fitted detectability model, iBAQ or Top3
  abundance and the genome context columns. `detection_summary.tsv` compares observed with expected
  detection, and `regions_stats.tsv` tests prophages, plasmids and other regions with circular-shift
  nulls.
- <!-- m05 -->
  Autocorrelation and islands: `autocorrelation.tsv` reports Moran's I and join counts along each
  replicon under label and transcription-unit block nulls, and `islands.tsv` lists runs of
  neighbouring genes with unusually high or low abundance, with family-wise error control. Abundance
  is also tested for a gradient from the origin to the terminus of replication.
- <!-- m06a -->
  Operons: `pairs.tsv` describes every pair of adjacent genes (predicted operon pair or not, strand,
  distance, promoter and terminator evidence), and `operons.tsv` summarises each transcription unit.
  GIPPA tests whether genes in one operon disagree in abundance (discordance), finds likely operons that
  the prediction missed with a decoy FDR, and lists flags in `operon_flags.tsv`.
- <!-- m06b -->
  Operon and terminator calls can be checked against RegulonDB: `validation/fetch_regulondb.py` (in the source repository)
  downloads the RegulonDB transcription units, terminators and gene coordinates with SHA-256
  provenance (TLS verification stays on), and `validation/regulondb_truth.py` builds the E. coli
  K-12 pair and terminator truth sets used for validation.
- <!-- m07 -->
  Variants and novel ORFs: peptides that disagree with the genome give REF/ALT residue calls in
  `variants.tsv` and `variants.vcf` (VCF 4.3), with flags and class error rates. Peptides outside
  annotated genes are assembled into candidate novel ORFs in `novel_orfs.tsv`, with GFF3 and protein
  FASTA files, promoter and ribosome-binding-site evidence from a second Promotech and OSTIR round,
  and `annotation_patch.gff3` to update the annotation.
- <!-- m08 -->
  Function and determinants: `category_summary.tsv` and `enrichment_detection.tsv` test GO, KEGG and
  COG categories for detection with a detectability-aware (Wallenius) over-representation test, and
  `enrichment_gsea.tsv` runs GSEA on abundance with genes collapsed to transcription units.
  `determinants_features.tsv`, `determinants_model.tsv` and `determinants_attribution.tsv` show which
  sequence and context features explain detection and abundance.
- <!-- m09 -->
  Plots and report: every run writes 16 plots (SVG, PNG and PDF, each with its data as TSV) and a
  self-contained, interactive `report.html` that works offline. `gippa report` rebuilds the report from
  a finished run directory, and `gippa compare` compares two or more runs of one genome (detection
  overlap, abundance agreement, stability) with its own plot and `compare.html`.
- <!-- m10 -->
  Real-data support for the pilot: DIA-NN reports from public E. coli Nissle 1917 proteomics data
  (iProX PXD078459) can be searched against a generic K-12 reference proteome or against the Nissle
  database from `gippa makedb`, and the two runs can be compared with `gippa compare`.
  `gippa convert --from cascadia` reads de novo peptides. A reproducibility runbook in the source
  repository lists every command and input hash. Pilot results are reported with the release.
- A Bioconda recipe for GIPPA is rendered from `pyproject.toml` by
  `scripts/release/bioconda_recipe.py` (PyPI sdist source, conda-forge and bioconda dependencies,
  Nextflow and OpenJDK included), ready for a `bioconda-recipes` pull request after each release.
- At genetic code 4 or 25 the variants VCF header records the table as `##gippa_translation_table=<t>` (right after `##reference=`); code-11 VCFs are unchanged.
- Bakta, miniprot and tblastn now run with the genetic code set by `annotation.translation_table`: Bakta gets `--translation-table` at tables 4 and 25 (its command line at table 11 is unchanged), miniprot indexes the genome with `-T <table>`, and tblastn gets `-db_gencode <table>` appended to `mapping.proteins.tblastn_args`. The default `mapping.proteins.tblastn_args` is now `-seg no -evalue 1e-5` and no longer contains `-db_gencode`; one unquoted `-db_gencode` equal to the table still works, and any other value (or a quoted one such as `-db_gencode "4"`) stops the run. A precomputed Bakta annotation JSON that records a different translation table is refused: its step stops with an input error. `annotation_meta.json` and the run manifest's `annotation` block record `translation_table` (output schema version 1.1). The report's Methods section names the genetic code of every run, and at tables 4 and 25 the summary opens with a note on how TGA is read.
- Codon usage, novel ORFs, variants and `gippa makedb` now follow `annotation.translation_table`. Under genetic codes 4 and 25, codon usage accepts the run's table, so TGA joins the Trp family (code 4) or the Gly family (code 25) and can be a rare codon. A CDS with an internal TGA is scored. `novel_orfs.faa` translates TGA as W or G, and a pseudogene's junction is only looked for at the table's stops. In variant calling, a peptide with W or G at a TGA codon matches the reference instead of being a stop-codon site, and a substitution there is an ordinary SNV. A `U` over TGA is skipped as a rare residue. `gippa makedb` translates six-frame ORFs, N-terminal extensions, Bakta proteins without a sequence and `--patch` CDSs with the table taken from the resolved configuration (`--config` or `--set annotation.translation_table=N`), so a TGA in a patch CDS is no longer rejected under 4 or 25. `makedb_manifest.json` records `translation_table`, and `gippa makedb` refuses a Bakta JSON recorded at a different table. At table 11 every output is unchanged apart from the new manifest key.
- Content-addressed image keys and tags (`<version>-<key12>`), `workflow/conf/images.json` written by `update_digests.py`, and the shared image-cache helper.
- N-terminal peptides of genes annotated with a GTG or TTG start codon (Bakta writes their first residue as `M`) now map exactly to their own gene, flagged `initiator_m` in `peptide_hits.tsv` (output schema 1.2), instead of being lost, mapped as a strain difference with an ALT `M` at the start codon, or placed only at a chance copy elsewhere. The six-frame index is now format 4 and is rebuilt on the next fresh run (without `-resume`); an older index still loads but gives no initiator credit and logs a warning to rebuild it.
- New `gippa demo` command: runs a bundled synthetic genome (SIM100K) with precomputed layers in about a minute, with no Nextflow, containers, databases or network, and prints the report path (`--open` opens it). `gippa demo --full` runs the same inputs through `gippa run` with the real tools after `gippa setup`, as a check of the whole install. The demo data (under 1 MiB) ships in the wheel; the manifest records `demo: true` and the report shows a synthetic-data banner.
- Preflight now treats a pinned container image as the install of its tool under the `docker` and
  `apptainer` profiles: under docker no local environment is needed, and under apptainer the image
  must be in the image cache (`NXF_APPTAINER_CACHEDIR`, else `$GIPPA_ROOT/containers`, else
  `<db_dir>/containers`), so a pip install without `$GIPPA_ROOT` works. A tool that has only a conda
  environment is reported missing under `docker` or plain `apptainer`, with the reason, instead of
  running on the host without its environment; a missing image names `gippa setup --images`.
- Preflight now warns (`genome.genetic_code`) when a genome left at the default translation table 11 probably uses genetic code 4 (TGA = Trp, as in Mycoplasma and Spiroplasma): long open reading frames cover far less of the genome with TGA as a stop than with TGA as a codon (ratio < 0.55, at least 20 such ORFs). The warning appears in `preflight.json`, the manifest `warnings` and at the top of the report summary; it never changes any result (`--strict` exits 3, as for every warning).
- R fgsea runs from the new `gippa-fgsea` image (R 4.5.3, fgsea 1.36.2, `gippa_fgsea.R` baked in) when no `GIPPA_FGSEA_PREFIX` and no prebuilt `gippa-fgsea` env is found: `gippa run` exports `GIPPA_CONTAINER_RUNTIME` and `GIPPA_DB_DIR`, preflight accepts the pinned image under `docker` and plain `apptainer`, and `containers.yml` builds the image.
- Releases are built and published by a manual release workflow: it checks the version, the release
  readiness record and the changelog, builds and smoke-tests the wheel and sdist, uploads them to PyPI
  by trusted publishing and creates the tagged GitHub release with these notes.
- Selenocysteine and pyrrolysine codons: annotated Sec (`U` over TGA) and Pyl (`O` over TAG) codons
  are now consumable sites for tolerant peptide placement, with decoys mirrored as for annotated
  stops. Rows over them are flagged `sec_site` / `pyl_site` in `peptide_hits.tsv` (GFF3
  `recoded_codon`), `variants.tsv` reports the `U`/`O` residue as an ALT site, and
  `detection_summary.tsv` counts the two new `placement_flag` groups.
- Step images `gippa-step` and `gippa-step-r` (`workflow/containers/step*/Dockerfile`) carry GIPPA
  and its hashed, permissive-only Python dependencies (`requirements.lock`, checked by
  `scripts/release/check_image_licences.py`), plus R fgsea in `gippa-step-r`; they are tagged with
  the GIPPA version (ADR-0065). `scripts/dev/build_step_image.sh` builds dev copies.
  The image licence check fails closed (missing site, empty run, `--require NAME`, every licence
  signal checked), and the sdist-only dependencies build on a hash-checked setuptools.
- The `docker` profile now matches `apptainer`: `params.db_dir` is bound read-only, tasks run as the calling user (`-u $(id -u):$(id -g)`), the per-task `TMPDIR` is passed into the container, the BAKTA and antiSMASH images get the `ps` shim, and pinned-container processes never use conda (also under `-profile docker,conda`). The TMPDIR `beforeScript` and the conda-off selector live once in `workflow/conf/container_common.config`, included by both container profiles; `update_digests.py --nextflow-config` follows `includeConfig` to edit that selector. docker-smoke gains an advisory BAKTA DB-visibility job.
- The full calibration launcher now submits one job per suite, sized from the checkout's
  `calibration_plan.tsv` (lanes x CPUs, memory and time per suite), instead of a throttled 4-CPU
  array. It validates the plan first, refuses a launch whose CPUs would pass the shared-cluster cap
  or whose longest suite would run into the nightly, and cancels every job already submitted when a
  submission fails. It prints `jobs=<ids> rows=<id> log=<dir> cpus=<total>`.
- The release check now reuses calibration evidence from an earlier commit only when every change
  since then lies in documentation, CI, report or plot code (or is the version line alone), and says
  so in the row. It also checks that each container image key is current in `images.json`, and the
  release workflow verifies that every GHCR image pin is publicly pullable.
- Tool images on GHCR now carry one immutable tag per build, `<version>-<key>`, where the key hashes
  the image's Dockerfile, env file and locks; an existing tag is never rebuilt or overwritten. The
  images install from the committed conda and pip locks, and every downloaded source archive is
  checked against its sha256. BacTermFinder now has locks too.
- `annotation.translation_table` now accepts genetic codes 4 (Mycoplasma and Spiroplasma) and 25 (Candidate Division SR1 and Gracilibacteria) as well as 11 in preflight; any other value stops the run with an input error. A new preflight check, `config.tool_genetic_code`, refuses tool options that set a genetic code themselves (Bakta `--translation-table`, miniprot `-T`, or a tblastn `-db_gencode` that differs from the table or is repeated). The `genome.genetic_code` warning is now evaluated at every supported table: at table 11 it names codes 4 and 25 when TGA looks like a sense codon, and at tables 4 and 25 it warns when the genome looks like it uses code 11.
- `annotation.uniop.enabled: true | false | auto` (default `true`). Preflight resolves it before the
  run: it runs UniOP's `--help` once to check that UniOP works, and `auto` turns into `true` or
  `false` in `resolved_config.yaml`. `gippa run` then passes `params.run_uniop`, and UNIOP is skipped
  when it is false. The distance-model fallback is always reported, with its reason:
  - preflight warns for `false`, for `auto` resolving to false and for a failed probe, and quotes
    the last failed `image:uniop` setup;
  - under `--publication`, `auto` without UniOP is an error;
  - `gippa run` prints a closing warning line;
  - the manifest gains `operon_source` and `operon_source_reason`;
  - the report states the distance model's precision and recall.

  ADR-0061 row flags:
  - TU-block Moran's I rows carry `tu_block_tu_error` on every run;
  - under the distance model, those rows and every `both_detected_t` pair row also carry
    `tu_source_distance_model`, with their caveats in the report.
- `annotation.uniop.enabled` now defaults to `auto` (M11-14d2, ADR-0061). When UniOP is installed for the profile the run uses it as before. When it is not, for example because a build failed, offline, or without build rights, operons come from the ADR-0022 distance model and the affected rows are flagged `tu_source_distance_model`, instead of preflight stopping with exit 2. Set `enabled: true` to require UniOP. Under `--publication`, `auto` must still resolve to `true`.
- `gippa compare` warns when the runs differ in GIPPA minor version, tool version or image digest,
  database version or operon source (on stderr, in `logs/events.jsonl` and in `compare.html`).
  Runs whose output schema majors differ stop with exit 1 unless `--allow-cross-major` is given.
  `run_manifest.json` and `datapackage.json` now record `output_schema_version` (`1.0`), the version
  of the output tables and manifest, apart from the config's `schema_version`.
- `gippa run`, `gippa preflight` and `gippa setup` choose a Nextflow profile when `--profile` is not given (new default `auto`): FarmShare is detected by its Slurm cluster name, `sbatch` on `PATH` selects the new generic `slurm` profile, and the container runtime follows what `gippa setup --images` recorded (`setup.runtime` in `db_manifest.json`, schema 1.1), else Apptainer, Singularity, Docker or conda. The choice is printed and recorded in `run_manifest.json` (`launch`). The accepted profile names are read from the bundled `nextflow.config`; two executor profiles or two container runtimes in `--profile` are a usage error.
- `gippa setup --images --profile apptainer|docker` pulls the pinned container images a config
  needs (`--images all`: every pinned image), checks each pulled digest against its pin and records
  the images in `db_manifest.json`; `--cache-dir` chooses the Apptainer image cache. `gippa setup`
  itself now runs the Bakta and antiSMASH downloaders through Docker when no environment or cached
  image is available, so it needs neither `$GIPPA_ROOT` nor a source checkout. `gippa setup --check`
  reports a database installed from another release as `installed <old>, this release pins <new>`.
- `gippa setup --uniop --profile docker --accept-uniop-terms` fetches UniOP from its authors (pinned commit, sha256 checked), builds the image locally with docker and records it in `local_tools.config` and `db_manifest.json`; `setup --images` includes the step (`--no-uniop` skips it), and `gippa run` adds the local build under docker, apptainer and conda (M11-14u1).
- `gippa setup --uniop` builds UniOP for every runtime and asks first:
  - the terms notice is always printed; without `--accept-uniop-terms` setup asks
    `Build UniOP for your own use? [y/N]` on a terminal and builds nothing otherwise (no is the
    default; `--uniop` then exits 2, `--images` warns and exits 0);
  - apptainer builds `uniop.def` into a `0600` `.sif` (`--fakeroot`, then
    `--fakeroot --ignore-fakeroot-command`); plain `conda` creates `<cache-dir>/envs/gippa-uniop`
    from the lock; `farmshare*` profiles check and record the prebuilt env;
  - an explicit `setup --uniop` installs UniOP only, and a rerun with a current build does nothing;
  - `db_manifest.json` records the build as `local_build`; preflight accepts it only when its
    runtime, key and image ID, `.sif` sha256 or env commit still match, and otherwise names the
    setup command;
  - the run manifest's `uniop` entry carries `local_build`, `runtime`, `source_sha256` and
    `container_digest`.

  Stub runs (`-stub-run`) write `operon_meta.json` with `operon_source: uniop` when UNIOP ran, so
  the stub run's manifest records the operon source.
- `gippa setup` installs the official, checksum-pinned Nextflow 26.10.0 `-dist` launcher as `DIR/bin/nextflow` when no `nextflow` is on `PATH` (`--no-nextflow` skips it), with a per-item redirect allowlist; `gippa run` uses `nextflow` on `PATH`, else that launcher, and exits 2 when Java is missing or older than 17 (`gippa setup` only warns). The Java check picks Java exactly as the Nextflow launcher does (`NXF_JAVA_HOME`, an executable `JAVA_CMD`, a `JAVA_HOME` directory, then `PATH`), and `run_manifest.json` records the launcher under `launch.nextflow` (`path`, `source`: `path` or `db_dir`; volatile, outside the provenance digest).
- `gippa setup` no longer needs `--db-dir`: it uses `$GIPPA_DB_DIR`, else
  `$XDG_DATA_HOME/gippa/db` (`~/.local/share/gippa/db`), and records the directory in
  `~/.config/gippa/db_dir`, which `gippa run`, `gippa preflight --profile` and `gippa makedb` read
  when neither `--db-dir` nor `$GIPPA_DB_DIR` is given. Before downloading, setup prints the plan
  (each item and image with its size, the apptainer build space, the total, the free space and
  user quota of each filesystem), refuses with exit 2 when less than 1.1 x the need is free, warns
  about a default DB directory in a shared home, and asks `Proceed? [Y/n]` on a terminal (`--yes`
  skips it). A full disk or quota during a download or pull marks the item `failed` and leaves no
  partial file; `apptainer pull` runs with `--disable-cache` unless `APPTAINER_CACHEDIR` is set.
- `gippa setup` now installs everything by default: antiSMASH, GO, KEGG, COG and the contaminant library join Bakta light, geNOMAD and the Promotech models, and the container images are pulled when `--profile` names a runtime. Each item has a `--no-<item>` opt-out (`--no-antismash`, `--no-go`, `--no-kegg`, `--no-cog`, `--no-contaminants`, `--no-images`); `--minimal` leaves out geNOMAD and antiSMASH for small disks, and a later `gippa run` on such a directory runs `core,function` with a notice. Setup prints the geNOMAD (non-commercial) and KEGG (academic use) terms before downloading, and `run_manifest.json` and the report's Methods section record the `licence_notes` of what each run used.
- gippasim gains `SimConfig.n_readthrough`: planted stop-codon read-through genes, with their truth in `truth/readthrough.tsv`.

### Fixes

- SVG plots no longer carry matplotlib's float-hashed clip-path, marker and hatch ids: `write_instance` renumbers them (`clip0`, `marker0`, `hatch0`, …) in order of appearance, so 1-ulp layout noise (seen on macOS for V18) cannot change the SVG bytes. ([#452](https://github.com/khoa-yelo/gippa/pull/452))
- The UniOP preflight probe and the R fgsea engine no longer pull container images: their `docker run` carries `--pull=never`. A missing fgsea image fails at once with a message to run `gippa setup` (detected from docker's "No such image" and podman's "image not known" errors), and under docker preflight's fgsea row checks the pinned image is present locally. Nextflow's own step images are not covered yet (N6). ([#465](https://github.com/khoa-yelo/gippa/pull/465))
- `gippa setup`: an OSError while installing an item (for example EXDEV, ESTALE or EIO renaming a
  verified download into place on NFS) now fails that item with a message naming the errno and the
  path, exit 2, and the remaining items still install, instead of aborting with exit 5 and no
  `db_manifest.json`; a full disk or quota gives the usual disk-full advice. A directory database is
  swapped by renaming the old copy aside to `<path>.gippa-old` first, so a failed swap keeps the old
  copy, or the next `gippa setup` recovers it. ([#466](https://github.com/khoa-yelo/gippa/pull/466))
- - `gippa report OUT` on a local run (no `resolved_config.yaml`) now takes the run's config, modules included, from `run_manifest.json`, so a `--modules core` run's page is byte-identical to the run's own report.
- BacTermFinder results are now read from the table the tool actually writes (`<genome>_mean.csv`, with `SampleName` and `probability_mean` columns), so runs with `annotation.terminators.bactermfinder` enabled no longer fail when their output is parsed. Contig names that contain underscores are handled, and window coordinates are kept as written. Every scanned window is still reported; no probability cut-off is applied yet.
- FarmShare profile: barley-03 joins barley-01 and rye-02 on the slow-node exclusion list (quick Python steps stalled there for about 20 minutes at 1 % CPU).
- The TransTermHP terminator cross-check no longer fails on a complete circular genome with a gene across the origin: such a gene is now passed to TransTermHP as two in-range pieces instead of one with an end past the replicon length.
- The antiSMASH version is now recorded in the run manifest: the antiSMASH, Bakta, geNomad, miniprot and tblastn version captures discard stderr, so a tool warning no longer corrupts `versions.yml`. antiSMASH is right-sized to 4 CPUs and 8 GB per attempt instead of the `process_high` 16 CPUs and 48 GB. Tasks killed at the Slurm time limit (exit 140) are retried like exit 143. On FarmShare, memory-derived CPU counts divide by 4000 MB to match Slurm's MaxMemPerCPU, and the slow nodes barley-01 and rye-02 are excluded from task submission.
- The opt-in BacTermFinder terminator step now calls `genome_scan.py` with the arguments its README gives (step 3, batch 10000) and runs it from a task-local directory, so it no longer fails at start-up or writes into the shared tool directory.
- Tool runs no longer import Python packages from the user site (`~/.local`): every workflow task,
  including Docker and Apptainer tasks, and every `gippa setup` tool call run with
  `PYTHONNOUSERSITE=1`, so a user-installed ViennaRNA or biopython can no longer change OSTIR and
  other tool results. OSTIR's `versions.yml` now records the ViennaRNA version it actually used.
- V03 locus plots are ordered within each test family instead of on one mixed scale. Loci
  significant at their own family's level come first, then non-significant loci, then loci with
  no estimate, and the families (prophage, BGC, island, operon discordance, missed operon) take
  turns within each group, so under the `report.max_locus_plots` cap one family's significant loci
  can no longer crowd out another's.
  Operon flags from the two flag families are no longer merged into one locus. The index figure
  names each estimate (`q`, `FWER p`, `decoy-FDR q`) and labels the groups.
- `-profile docker` runs the pinned images: `containers.config` drops the `docker://` prefix that
  `docker run` rejects.
- `gippa setup` from an installed package reads the pinned images from the bundled workflow and
  reports a missing `GIPPA_ROOT` clearly, instead of looking under site-packages. The report section
  for V15 is titled "Determinants of abundance".
- `gippa setup`, preflight and `gippa run` now use one Apptainer image cache: setup's Bakta and antiSMASH downloaders look for their image in `--cache-dir` (else the default cache), and `gippa run` exports `NXF_APPTAINER_CACHEDIR` to the cache preflight checked when it is not set, so `DIR/containers` works without any variable. `gippa setup --images` pulls the images before it downloads the databases, so an Apptainer-only host installs Bakta (and antiSMASH) in one run. `gippa preflight --profile` checks the installed tools as `gippa run --profile` does, and `gippa run` without Nextflow points at the install page.
- `run_manifest.json` now records as `detection_tier` the tier that detection actually used, the same as `genes.tsv`. It no longer records preflight's provisional tier. Detection warnings such as `decoy_unscored` now appear in the manifest's `warnings`.

### Documentation

- The user guide has a new "Genetic code" section under Inputs: when to set `annotation.translation_table` to 4 (Mycoplasma, Spiroplasma and other Mollicutes) or 25 (candidate division SR1, Gracilibacteria), what the table changes in Bakta, miniprot, tblastn, mapping, novel ORFs, variants, codon usage, `gippa makedb` and the report, how `-db_gencode` is added to `mapping.proteins.tblastn_args` and which values stop a run, why a precomputed Bakta JSON made at another table is refused, and the preflight checks `config.translation_table`, `config.tool_genetic_code` and `genome.genetic_code`. A new methods page explains the genetic-code check's evidence and limits: codes 4 and 25 cannot be told apart from the genome, on inputs with only 20–40 long open reading frames either warning is wrong about 1 % of the time, and the code-25 evidence rests on one assembly.

### Other

- Add an advisory `docker-smoke` workflow (real MINIPROT and TBLASTN under `-profile test,docker` on
  path-filtered pull requests and weekly) and the platforms support matrix in the user guide.
- Known limitations of version 1:
    - licences: geNOMAD is for non-commercial use only. `gippa setup` downloads its database
      unless you pass `--no-genomad` (or `--minimal`). KEGG is for academic use;
      `gippa setup --kegg` downloads it on your machine, and GIPPA does not redistribute it. The GPL and AGPL
      tools (Bakta, Promotech, Promoter Calculator, OSTIR, TransTermHP, BacTermFinder, antiSMASH)
      run as separate programs; GIPPA does not import them, and the wheel and sdist contain none
      of their code. UniOP has no licence file, so it is never published; `gippa setup --uniop`
      builds it on your machine (see `NOTICE.md`);
    - one sample per run: several samples in the abundance table, co-expression and differential
      context tests are not supported yet;
    - the V4 locus browser (igv-reports) and the V9 mass-fraction treemap are not built, so 16
      plots are made instead of 18;
    - no RO-Crate export;
    - full validation (K-12 PXD075568, UTI89 beyond the pilot, the EDL933 holdout, FDRBench
      entrapment) is planned after 1.0; 1.0 is calibrated on simulations;
    - stop codons: a peptide is placed across a stop codon only at a gene's annotated terminal
      stop, or at a selenocysteine (`U` over TGA) or pyrrolysine (`O` over TAG) codon that Bakta
      annotates, and only if it has at least 10 residues; other stop codons are never read
      through. Such a peptide has one difference from the genome, so it is judged in a small
      placement-FDR stratum and is rejected when that stratum holds few peptides. On E. coli
      Nissle 1917 the selenocysteine peptides of the two fdnG copies are rejected this way
      (q = 0.25), so their Sec sites are not reported. A C-terminal read-through extension
      (`cterm_extension`, labelled `readthrough_candidate`) is called only when a peptide crosses
      the annotated stop into the 3′ flank and the ORF has at least two distinct peptides; it is
      written to `annotation_patch.gff3`, where `gippa makedb --patch` reads the old stop as `X`.
      Under genetic codes 4 and 25 a `U` over TGA is not recoded;
    - genetic codes: the code comes from `annotation.translation_table` (11, 4 or 25) in the
      configuration and is never inferred. Preflight warns (`genome.genetic_code`) when TGA looks
      like a sense codon at table 11, or like a stop at tables 4 and 25, but the genome cannot
      tell code 4 (TGA = Trp) from code 25 (TGA = Gly), so you choose. With only 20–40 long open
      reading frames (a plasmid alone, a partial assembly) the warning is wrong about 1 % of the
      time. When Bakta calls no tRNA-Trp(UCA) gene, as on Mycoplasma hyorhinis, the tRNA
      adaptation index (tAI) at code 4 is approximate. `gippa compare` does not check that its
      runs used the same genetic code;
    - spectrum inputs: `gippa convert` reads DIA-NN, Cascadia and the generic table only;
      `--from casanovo` and `--from instanovo` are accepted but not implemented, and DDA search
      results must be written as the generic table;
    - variant calls: precision is checked on unflagged substitution calls; `mass_equivalent`
      calls and block, swap and frameshift calls are reported but not covered by that check;
    - low-power tests: missed operons (with the shipped decoys), low islands at normal detection
      rates and operon discordance at realistic effect sizes have little or no power by
      construction, so an empty result there is not evidence of absence;
    - older Linux (glibc < 2.28, such as CentOS 7): `pip install gippa` fails to build numpy,
      numba or pyarrow from source; use `pip install --prefer-binary gippa`. The older versions
      it picks are not in the tested matrix
      ([#502](https://github.com/khoa-yelo/gippa/issues/502));
    - `gippa setup` retries a failed Apptainer image pull once automatically; on some clusters
      `mksquashfs` exit status 139 still needs `APPTAINER_CACHEDIR` set
      ([#504](https://github.com/khoa-yelo/gippa/issues/504)).
- Refreshed the bactermfinder, cascadia and timsconvert env locks from fresh builds; the pip locks no longer list packages from the user's own site-packages.
- TerminatorNet asks for the process_high CPUs (clamped by `--threads`) and Promotech for 2 CPUs and 10 GB per chunk, matching their measured use; outputs are unchanged (ADR-0062).
- The README describes the released package: install from PyPI or Bioconda, a short quickstart, how
  to cite GIPPA, and a link to the documentation site.
- The determinants step refits the permutation null of both models on one pool of `threads` worker processes instead of a LightGBM thread pool and a serial elastic net; the null values, p and b are bit-identical (ADR-0063).
- The fixed SVG hash salt is now `gippa`; written SVGs are unchanged, because their hashed ids are renumbered on write.
- The miniprot conda environment now has a committed lock file, so `build_envs.sh` can rebuild it exactly.
- The project, Python package, command and environment variables are renamed to `gippa` / `GIPPA_*` (ADR-0064); `scripts/farmshare_env.sh` still exports the old variable names for a transition window.
- The six-frame index, exact and tolerant peptide mapping, the annotated read-through stop sites and the N-terminal extension walk now use `annotation.translation_table`: under genetic codes 4 and 25, TGA reads as W and G and is never a stop, a Sec (`U`) residue cannot be recoded (it is counted as a codon mismatch and never gives `sec_site`), and Pyl (`O` over TAG) works as under 11. The six-frame index work file records the table it was built with (index format 3), and a run with another table refuses it with a message to rebuild the six-frame index; an older format-2 index is still accepted at table 11. The `mapping.stop_sites` log event reports the table. Output at table 11 is unchanged.
- `calibration.yaml: operon_distance` records the distance-model fallback's precision and recall at
  `pair_threshold` 0.5 on RegulonDB 14.5.0 K-12 (0.828 [0.805, 0.847] and 0.916 [0.901, 0.931],
  directon-bootstrap CIs; ADR-0022 addendum).
