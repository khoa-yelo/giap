# M. hyorhinis SK76 at genetic code 4 (M11-16g, spec item 10)

*Mesomycoplasma hyorhinis* SK76 (GCF_000313635.1, 836,897 bp, one circular chromosome) uses NCBI
genetic code 4 (TGA = Trp). PR #411 ran it at the default table 11. This file compares that run
(`$GIPPA_ROOT/more_datasets/mhyo/out`, code 11) with a re-run at code 4
(`$GIPPA_ROOT/more_datasets/mhyo/out_code4`), as `docs/specs/M11-16-translation-tables-4-25.md`
item 10 asks. The thresholds are validation expectations, not CI gates (decision D12). A miss is
reported with its cause. No code, threshold or tolerance was changed.

## Runs

| Job | What | Result |
|---|---|---|
| 1903202 | first submission of the code-4 run | cancelled by the validation agent within seconds, before it ran anything: `GIPPA_REPO` pointed at the pinned worktree but `GIPPA_ROOT` was not exported, so `farmshare_env.sh` derived `GIPPA_ROOT` from the worktree and every input path was wrong |
| 1903204 | code-4 run, the spec's command | **FAILED** after 36 min: `ANNOTATE:TRANSTERMHP` exit 3, `Bad gene coordinates: 836688 .. 838112` (finding F1) |
| 1903516 | the same command plus `--set annotation.terminators.cross_check_transtermhp=false --resume` | COMPLETED, exit 0, 36 min (39 tasks run, 13 cached) |
| 1903218 | K-12 code-11 control at `6c8d90b9` (main just before M11-16a) | COMPLETED, 3 h 25 min |
| 1903219 | K-12 code-11 control at `77df81f5` (main with 16a-16e) | COMPLETED, 3 h 30 min |

- Code: a detached worktree at `origin/main` `77df81f5` (`local/worktrees/val-code4`). The
  manifest records `git_commit 77df81f5…`, `git_dirty false`.
- Command (spec item 10): `sbatch -J gippa-run-mhyo-code4 scripts/slurm/gippa_run.sbatch --genome
  …/GCF_000313635.1_ASM31363v1_genomic.fna.gz --peptides …/inputs/peptides.tsv.gz --db-dir
  $GIPPA_DB_DIR --threads 6 --memory 24.GB --set annotation.bakta.complete=true --set
  annotation.translation_table=4 --outdir $GIPPA_ROOT/more_datasets/mhyo/out_code4`, with
  `OUT/nextflow.config` = `executor { queueSize = 4 }` (as #411's later runs). Bakta 1.12.2, DB
  6.0 full, run live.
- The code-11 baseline `out` is #411's run (`--threads 3 --memory 12.GB`, otherwise the same
  inputs and settings). It predates M11-14w, so its preflight has no `genome.genetic_code` row;
  the row below was recomputed with `gippa.preflight.genetic_code.genetic_code_check` at
  `77df81f5` on the same FASTA.
- Deviation in the completed run: TransTermHP cross-check off (F1). It only feeds
  `terminators_crosscheck.tsv`; TerminatorNet (the terminator tool) ran as usual, and no measure
  below reads terminators.
- R14: at most 4 Nextflow task jobs per run (`queueSize = 4`). The GIPPA total before each
  submission was 13 to 48 CPUs, so it stayed under the 128-CPU cap.

## RefSeq annotation (R6)

`GCF_000313635.1_ASM31363v1_genomic.gff.gz` from
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/313/635/GCF_000313635.1_ASM31363v1/`
(2026-10-28 21:31 UTC), MD5 `e6700474126cbe028a241476b52101c9`, equal to the line in that
directory's `md5checksums.txt`. Stored in `$GIPPA_ROOT/more_datasets/mhyo/raw/refseq/`. The
fresh `md5checksums.txt` lists the same `_genomic.fna.gz` MD5 as #411's copy. RefSeq has 720
genes: **684 protein-coding** (684 non-pseudo CDSs), 27 pseudogenes, 30 tRNAs, 3 rRNAs and 4 other
ncRNAs. The spec's "about 750" is the gene count with RNAs and pseudogenes. The ±10 % band below
uses the 684 protein-coding genes: 616 to 752.

## Comparison (spec item 10 table)

| Measure | Code 11 (`out`) | Code 4 (`out_code4`) | Expected at code 4 | Verdict |
|---|---|---|---|---|
| preflight `genome.genetic_code` | `warn` (W11), ratio 0.227, `translation_table: 11` | `ok`, ratio 0.226994, `translation_table: 4` (`config.translation_table` ok 4, `config.tool_genetic_code` ok `[]`) | `ok`, ratio 0.227, `translation_table: 4` | pass |
| Bakta CDS / pseudogene | 1,200 / 274 | **754** / 6 | CDS 616-752 (±10 % of 684); pseudogenes ≤ 68 | **CDS miss by 2** (+10.2 %, F2); pseudogenes pass |
| CDSs ending in TGA | 801 (TAA 340, TAG 59) | 0 (TAA 632, TAG 122) | 0 | pass |
| median CDS protein | 80 aa (90 over all features) | 256.5 aa (255 over all features) | ≥ 250 aa | pass |
| RefSeq CDSs with a Bakta CDS on the same strand and stop codon | 319 of 684 (46.6 %) | 673 of 684 (98.4 %) | ≥ 90 % | pass |
| the 390 TGA-split adjacent pairs | 390 | 373 of 390 (95.6 %) inside one code-4 CDS, same strand and frame | ≥ 90 % | pass |
| organism peptides `no_hit` / those with W | 8,163 / 1,524 | 7,737 / 1,098 | W-containing `no_hit` ≤ 304 (20 % of 1,524) | **miss** (1,098 = 72 % of 1,524, F3) |
| distinct mapped peptides (target rows) | 10,214 (10,446 rows) | 10,839 (11,079 rows) | ≥ 10,214 | pass, gain +625 |
| flags `readthrough` / `strain_difference` / `nterm_extension_evidence` (rows) | 316 / 321 / 975 | 0 / 0 / 9 | ≤ 79 / 80 / 243 | pass |
| `sec_site` | 0 | 0 (`detection_summary.tsv` `placement_flag sec_site` 0; `pyl_site` 0) | 0 | pass |
| `trna_counts.tsv` anticodon `TCA` (tRNA-Trp UCA) | absent | **absent**: Bakta calls 27 tRNA rows, Trp only as `CCA` (1 copy) | report | reported (F4) |
| codon usage: Trp family synonymous shares | TGG 1.0 (TGA is a stop) | TGA 0.939, TGG 0.061 (2,378 and 154 internal codons over 754 CDSs) | TGA and TGG both present | pass |

Definitions used (script logic, reproducible from the files named):

- CDS and pseudogene counts: `genes.parquet` `feature_type`. Stop codon: the last 3 nt of `nt`.
  Protein length: `len(aa)` without a trailing `*`.
- RefSeq match: a RefSeq non-pseudo CDS (GFF 1-based, converted to 0-based half-open) whose stop
  end (plus: `end`; minus: `start`) and strand equal a Bakta CDS's.
- TGA-split pairs: code-11 CDSs adjacent by position on one strand, same frame (plus: `end mod 3`;
  minus: `start mod 3`), gap ≤ 300 bp, the upstream one (in transcription direction) ending in
  TGA. This reproduces the spec's 390. A pair is inside one code-4 CDS when a code-4 CDS on that
  strand covers both and has the same frame.
- Peptides: `mapping/unmapped.tsv` rows with `kind = peptide`, `label = unmapped`,
  `reason = no_hit`; distinct mapped peptides are the distinct `peptide` values of
  `mapping/peptide_hits.tsv` (it holds target rows only). Flags are counted per row of
  `peptide_hits.tsv` (the spec's 316 / 321 / 975 are row counts; distinct peptides are 314 / 319 /
  970).
- Trp shares: `gippa.annotate.codon_usage.synonymous_shares(cds_nt, table)` at the run's table
  (`codon.parquet` stores per-gene indices, not the shares). `codon.parquet` has CAI, tAI, ENC
  and Fop for all 754 CDSs.
- Raw numbers: `measures_code11.json` and `measures_code4.json` next to this file.

The report has the code-4 note (`data-genetic-code="4"`: "Genetic code 4 (Mycoplasma and
Spiroplasma): TGA is read as W.") and the Methods sentence; `run_manifest.json` has
`annotation.translation_table: 4` and `output_schema_version 1.1`; `sixframe.npz` meta is
`format 3`, `translation_table 4`.

## Findings

**F1 (bug, not code-4 specific): TransTermHP fails on an origin-spanning gene.** At code 4 with
`annotation.bakta.complete=true`, Bakta calls one CDS across the origin (`CAPMEI_0796`, 0-based
836,687-838,112 on the 836,897 bp chromosome; at code 11 the same locus was two partial CDSs cut
at TGA, `CAPMEI_01510` and `CAPMEI_00001`). `write_transterm_coords`
(`src/gippa/annotate/terminators.py:95`) writes the internal interval as is, `836688 838112`, past
the sequence end, and TransTermHP 2.09 exits 3 ("Bad gene coordinates"). Any complete circular
genome with a gene across the origin hits this at any table; K-12 and the #411 runs had none.
The completed run turned the cross-check off; a fix (wrap or skip origin-crossing genes in
`genes.coords`) is for an implementer. AC 16 ("the run completes, exit 0") is met only with this
deviation.

**F2 (expectation, not GIPPA): 754 CDSs is 2 above the +10 % band.** 673 Bakta CDSs share a stop
with a RefSeq CDS. **73 more lie on 27 RefSeq pseudogenes**, mostly IS3-family transposase and
"Mobile element protein" fragments plus one N-acetylneuraminate operon remnant: RefSeq records
each frameshifted element as one pseudogene, Bakta calls the intact pieces as CDSs. The other 8
are hypothetical proteins with no RefSeq feature (35-155 aa). Without the 73 pseudogene pieces
Bakta has 681 CDSs against 684. Bakta CDSs + pseudogenes (760) against RefSeq protein-coding +
pseudogenes (711) is +6.9 %. The code-4 gene model is right; the band does not count pseudogene
fragments.

**F3 (expectation, not GIPPA): W-containing `no_hit` stays at 1,098 (expected ≤ 304).** Of the
1,524 W-containing `no_hit` peptides at code 11, 425 map at code 4 and 1 is now
`rejected_placement_fdr`. The other 1,098 are still `no_hit`, and **every one of them is assigned
by DIA-NN only to accessions outside the M. hyorhinis ORF set** (UniProt entries such as P02766,
P05386, P26641: the host cell proteins of the deposit, which `contaminant` does not cover). No
code-4 `no_hit` peptide is attributed to an M. hyorhinis ORF. W occurs in 14.2 % of the code-4
`no_hit` peptides, the background rate for peptides of this length. The expectation assumed that
the W-containing `no_hit` peptides were TGA artefacts; only 425 were. All 425, plus 204 that were
`too_short` and 1 `rejected_placement_fdr` at code 11, map at code 4 (630 in all); 5 peptides
mapped at code 11 are `rejected_placement_fdr` at code 4, hence the net gain of 625.

**F4 (Q-L5, D10): Bakta does not call tRNA-Trp(UCA).** `trna_counts.tsv` has Trp only as `CCA`,
the same 27 rows as at code 11, although TGA is 94 % of Trp codons. tAI at code 4 therefore has no
anticodon for TGA from Bakta's tRNAs. Per D10 this is recorded, not patched; the tAI column is
approximate at code 4.

## Q-L8: the backlog's "628 unmapped organism peptides, all with W"

The figure was written in the M11-14w spec and backlog (commit `109e043b`) from #411's notes,
which do not define the set. Re-derived from `out`:

- Every unmapped organism peptide (`label = unmapped`, any reason) that occurs exactly in the
  code-4 six-frame translation of the genome (I = L, circular wrap): **630** = 425 `no_hit` + 204
  `too_short` + 1 `rejected_placement_fdr`. All 630 contain W, by construction: none occurs in the
  code-11 six-frame translation, so each spans a TGA. This is within 2 of 628, and is the only
  definition found that gives "all with W"; the exact 628 was not reproduced (the old note's
  filter is unknown).
- Over the whole `no_hit` set, the statement that holds is the spec's recomputed **1,524 of
  8,163** (18.7 %) contain W. "All unmapped organism peptides contain W" is false as written: there
  are 8,163 `no_hit`, 2,858 `too_short` and 8 `rejected_placement_fdr` unmapped organism peptides.
- At code 4 all 630 map (above), so the backlog number and the code-4 gain agree.

The backlog sentence should read "1,524 of the 8,163 `no_hit` organism peptides contain W; 630
unmapped organism peptides match the code-4 six-frame translation exactly, all through a TGA".

## Code-11 control: E. coli K-12 (AC 1)

No weekly run has ever written `$GIPPA_ROOT/ci/runs/weekly_*` (there is no `logs/weekly.tsv`);
the nightly runs only sim100k. The newest K-12 run on main is MVP-04 A2
(`$GIPPA_ROOT/runs/mvp04/20261023_a2`, `8b6e93d4`, 2026-10-26), which predates many other
merged changes, so a byte comparison against it would not isolate M11-16. The control instead ran
the same inputs at two commits: `6c8d90b9` (main just before 16a) and `77df81f5` (main with
16a-16e).

- Inputs: `$GIPPA_ROOT/data/genomes/k12_mg1655/GCA_000005845.2_ASM584v2_genomic.fna.gz` with
  `[topology=circular]` (as `weekly.sbatch`), MVP-04's Bakta JSON
  (`$GIPPA_ROOT/mvp04/k12_mg1655.bakta.json`), MVP-04's peptides and config, seed 20261023,
  the weekly 60 % protein stand-in from the recorded full-DB FAA (so miniprot and tblastn run),
  `--modules core,mobile,bgc --threads 6 --memory 24.GB`, `queueSize = 4`. Inputs and SHA-256 in
  `$GIPPA_ROOT/validation/code4_k12/inputs/`.
- Result: **2,875 files byte-identical, 12 different. Every file under `tables/` (19) and
  `mapping/` (7) is identical**, as are all parquet layers, `genes.coords`, `codon.parquet`,
  `trna_counts.tsv` and the novel-ORF outputs. The sixframe.npz arrays are equal; only `meta`
  differs.

The 12 differences, all documented:

| File | Difference |
|---|---|
| `preflight/preflight.json`, `preflight/preflight.md` | new `config.tool_genetic_code` row (ok, `[]`); `genome.genetic_code` rule text and `value.translation_table: 11`; `n_ok` 18 → 19; also the `config.translation_table` rule text ("in {4, 11, 25}", item 2.1) |
| `plots/V17_qc.{pdf,png,svg}`, `plots/data/V17_qc.tsv` | ok bar 18 → 19 (the brief's 14 → 15 is the same +1 on other inputs) |
| `schemas/datapackage.json` | `output_schema_version` 1.0 → 1.1 |
| `resolved_config.yaml` | `tblastn_args` `-db_gencode 11 -seg no -evalue 1e-5` → `-seg no -evalue 1e-5` (D1) |
| `annotation_meta.json` | gains `translation_table: 11` (the source of the manifest field, item 7) |
| `run_manifest.json` | `annotation.translation_table: 11`, `output_schema_version` 1.1, config snapshot `tblastn_args`, `provenance_digest`; plus the run-specific `git_commit`, `outdir`, times and step order/timings |
| `sixframe.npz` | meta `format` 2 → 3 with `translation_table: 11` |
| `report/report.html` | the Methods sentence "Genes were called and all six frames translated with NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)."; the V17 plot; the commit, manifest, provenance and config hashes; no `data-genetic-code` note |

No other difference.
