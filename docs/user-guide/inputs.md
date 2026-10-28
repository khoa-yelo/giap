# Inputs

GIPPA reads one assembled genome and proteomics evidence from the search you already ran. It never
searches spectra itself. Run `gippa preflight` to check every input in seconds before a run (see
[Running GIPPA](run.md)).

| Input | Option | Required | Format |
|---|---|---|---|
| Genome | `--genome` | yes | nucleotide FASTA, optionally gzipped |
| Peptides | `--peptides` | recommended | the generic peptide table (below) |
| Proteins | `--proteins` | yes, unless peptides are given | protein FASTA |
| Abundance | `--abundance` | optional | TSV: a protein ID column plus one value column |
| Annotation | `--annotation` | optional | Bakta JSON (preferred), GFF3 or GenBank (GBFF) |
| Reference peptides | `--reference-peptides` | optional | the generic peptide table |
| Replicon depth | `--replicon-depth` | optional | TSV with columns `contig` and `mean_depth` |
| Configuration | `--config` | optional | YAML; see the [configuration reference](../reference/configuration.md) |

`gippa preflight` takes `--genome`, `--proteins`, `--peptides`, `--abundance` and `--annotation`.
`gippa run` takes them all. Without `--peptides` and `--proteins`, preflight stops with the error
`inputs.evidence_present` (exit code 1). It does the same when every peptide row was dropped and
no proteins were given.

## Genome FASTA

- The first token of each header is the contig ID. Characters outside `[A-Za-z0-9._-]` are
  replaced by `_`. The renaming is saved and reported as the warning `genome.ids_renamed`.
- **Circular replicons.** A header that contains `circular=true` or `[topology=circular]` marks the
  replicon as circular. Contigs without either tag are linear. The exception is
  `--set annotation.bakta.complete=true`, which declares every contig a complete, circular
  replicon.
- Gzipped FASTA is accepted.
- **Genetic code.** GIPPA translates with NCBI genetic code 11 by default. Codes 4 and 25 are
  also supported. See [Genetic code](#genetic-code) below.

Topology matters. On a circular replicon, a gene or peptide can cross the origin, and the
statistics wrap around it. On a linear contig, nothing wraps. A gene cut by a contig end is
`partial` and is excluded from detection statistics.

**Plain NCBI FASTA.** NCBI headers have no topology tag. Preflight therefore reports every contig
as linear in the info row `genome.topology`, which lists each replicon's length, topology and kind.
For a complete assembly, add the tag or set `annotation.bakta.complete`. Preflight also lists
contigs shorter than `preflight.small_replicon_bp` (500,000 bp by default) as putative plasmids or
fragments, in the info row `genome.small_replicons`.

Other genome checks warn on the following:

- more than `preflight.max_contigs_warn` contigs;
- a total length or N50 below the configured minimum;
- runs of N;
- a GC content outside `preflight.gc_range_warn`.

## Genetic code

One key, `annotation.translation_table`, sets the NCBI genetic code for every step of a run. It
accepts three values:

| Value | NCBI name | TGA | Use it for |
|---|---|---|---|
| `11` (default) | Bacterial, Archaeal and Plant Plastid | stop (Sec when recoded) | most bacteria and archaea |
| `4` | Mycoplasma and Spiroplasma | W (Trp) | Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma and other Mollicutes |
| `25` | Candidate Division SR1 and Gracilibacteria | G (Gly) | candidate division SR1 and Gracilibacteria |

The three tables differ only at TGA. The other 63 codons translate the same way in all three.

```bash
gippa run --genome mhyo.fna.gz --peptides peptides.tsv --set annotation.translation_table=4 ...
```

You can also put it in a YAML file passed with `--config`:

```yaml
annotation:
  translation_table: 4
```

**Choose the table from the organism's taxonomy.** The genome sequence alone cannot tell code 4
from code 25, because both read TGA as a sense codon. Preflight can tell whether TGA acts as a stop
(see "Preflight checks" below), but not which of the two sense codes is right. A code-25 genome run
at table 4, or a code-4 genome run at table 25, gives no warning. Peptides with G (or W) at TGA
codons then map only through the tolerant route, with one difference, as `strain_difference`
under `related`.

### What the table changes

Bakta (GPL-3.0), miniprot (MIT) and BLAST+ (public domain) run as separate programs; GIPPA passes
the table on their command lines and imports none of them ([Licences](licences.md)).

- **Bakta** runs as `bakta --translation-table 4` (or `25`). At table 11 no flag is passed,
  because 11 is Bakta's own default, so the Bakta command line and its cached results do not
  change. A GFF3 or GBFF given with `--annotation` is re-called by `bakta --regions` with the
  configured table. Its `/transl_table` qualifier is not read.
- **miniprot** indexes the genome with `-T <table>`.
- **tblastn** gets `-db_gencode <table>`. See "tblastn arguments" below.
- **Six-frame translation and peptide mapping.** TGA reads as W under 4 and as G under 25, and is
  never a stop. A peptide with W or G at a TGA codon maps exactly. The six-frame index records the
  table it was built with. A run with another table refuses it and tells you to rebuild the
  six-frame index.
- **Read-through and recoded stops.** Under 4 and 25 only TAA and TAG are stops. Selenocysteine
  cannot be recoded over TGA. A `U` in a precomputed protein is counted as a codon mismatch and
  never gives `sec_site`. Pyrrolysine (`O` over TAG) works as under 11 and still gives `pyl_site`.
- **Novel ORFs** end only at the table's stops, and `novel_orfs.faa` translates TGA as W or G.
- **Variants.** A peptide with W or G at a TGA codon matches the reference, and a substitution
  there is an ordinary SNV. Under 4, a strain SNV TGA → TGG is synonymous (W → W), not a stop
  loss. A `U` over TGA is skipped as a rare residue.
- **Codon usage.** Under 4, TGA joins the Trp family (TGA, TGG). Under 25, it joins the Gly family
  (GGN and TGA). A CDS with an internal TGA is scored.
- **`gippa makedb`** reads the table from `--config` or `--set annotation.translation_table=N`
  (there is no separate flag). Six-frame ORFs, N-terminal extensions, Bakta proteins without a
  sequence and `--patch` CDSs are translated with it, and `makedb_manifest.json` records
  `translation_table`.
- **Outputs.** `annotation_meta.json` and the `annotation` block of `run_manifest.json` record
  `translation_table`. The report's Methods section names the code of every run: "Genes were
  called and all six frames translated with NCBI genetic code 4 (Mycoplasma and Spiroplasma)."
  At tables 4 and 25, the summary also opens with a note such as "Genetic code 4 (Mycoplasma and
  Spiroplasma): TGA is read as W."

These do not change:

- **Start codons.** `annotation.start_codons` (default ATG, GTG, TTG) does not follow the table.
  NCBI lists more starts for code 4, but Bakta's gene caller uses ATG, GTG and TTG at every table.
  You can add starts to `annotation.start_codons` yourself.
- **geNOMAD** (module `mobile`, non-commercial use only; see [Licences](licences.md)) chooses its
  own genetic code and reports it in its `genetic_code` column. GIPPA does not pass the table to
  it or compare the two.
- **antiSMASH** reads Bakta's GBFF, which carries `/transl_table`. Promoter, terminator, OSTIR and
  UniOP steps work on nucleotides only.
- **One table per run.** A plasmid or phage that uses a different code from its chromosome is not
  supported.

### tblastn arguments

The default `mapping.proteins.tblastn_args` is `-seg no -evalue 1e-5`. GIPPA appends
`-db_gencode <table>` to it, so the default runs with `-seg no -evalue 1e-5 -db_gencode 4` at
table 4. If you set your own `tblastn_args`, leave `-db_gencode` out and GIPPA adds it. A
`-db_gencode` equal to the table, for example `-db_gencode 11` at table 11, is passed through
unchanged. These stop the run:

- `-db_gencode` with a value other than the table, or a trailing `-db_gencode` with no value:
  preflight error `config.tool_genetic_code`;
- `-db_gencode` given more than once, even with the same value: preflight error
  `config.tool_genetic_code`;
- a quoted value such as `-db_gencode '4'`, even when it equals the table: preflight error
  `config.tool_genetic_code`. Write the number without quotes, or better, leave `-db_gencode` out.

### Precomputed Bakta JSON

Bakta records its table in the JSON (`genome.translation_table`). When you pass a Bakta JSON with
`--annotation`, or to `gippa makedb --annotation`, and it records a different table from the run,
the annotation is refused with this message:

```text
annotation JSON was made with translation table 11; this run uses 4: re-run Bakta with
--translation-table 4 or set annotation.translation_table: 11
```

`gippa makedb` exits with code 1. In `gippa run` the annotation step fails with an input error,
and the run stops with exit code 4; the manifest names the failed step and its log. Preflight does
not open the JSON, so this happens when the annotation is read, not in preflight.
A JSON without `genome.translation_table` (older or hand-made) is assumed to use the run's table,
and the log records the info event `annotation.table_unrecorded`.

### Preflight checks

| Check | Level | When |
|---|---|---|
| `config.translation_table` | error (exit 1) | the table is not 4, 11 or 25. Message: `annotation.translation_table=12: GIPPA supports genetic codes 4, 11 and 25 (the tables Bakta accepts)` |
| `config.tool_genetic_code` | error (exit 1) | a tool option sets a genetic code itself: Bakta's translation-table option (`bakta --translation-table`, or an abbreviation of it of at least 4 characters, such as `bakta --transl`) in `annotation.bakta.extra_args`; `-T` (also inside a short-option cluster such as `-uT4`) in `mapping.proteins.miniprot_args`; a `-db_gencode` in `mapping.proteins.tblastn_args` that differs from the table, is repeated, is quoted or has no value. An option string the shell cannot split (for example an unbalanced quote) is also an error. `value` lists the offending keys. Remove the option and set `annotation.translation_table` instead |
| `genome.genetic_code` | warn | the genome does not look like it uses the configured table (below) |
| `genome.genetic_code` | info | fewer than 20 open reading frames of at least 300 codons, so the check is not evaluated |

`genome.genetic_code` compares how many bases long open reading frames (at least 300 codons, six
frames) cover when TGA is a stop with how many they cover when TGA is a sense codon. The ratio is
low (below 0.55) for genomes that use code 4 or 25 and high for code-11 genomes:

- at table 11, a ratio below 0.55 warns that the genome probably uses code 4 or code 25. Set
  `annotation.translation_table` to 4 or 25 from the organism's taxonomy;
- at tables 4 and 25, a ratio of 0.55 or more warns that the genome probably uses code 11. Genes
  that end in TGA would be read through into the next gene.

The warning does not change any result, and the report's summary repeats it. `--strict` turns it
into exit code 3. On whole genomes the two classes are well separated. On partial inputs with only
20–40 long open reading frames, such as a plasmid alone or a small part of an assembly, either
warning is wrong about 1 % of the time. See
[Genetic codes 4, 11 and 25](../methods/genetic-code.md) for the evidence and its limits.

## Peptide table (`peptides.tsv`)

The generic peptide table is a UTF-8 TSV with one header row. `gippa convert` writes it from
DIA-NN or Cascadia output (see [Converting search-engine output](convert.md)). You can also write
it yourself. Only `peptide` is required.

| Column | Type | Meaning |
|---|---|---|
| `peptide` | text | stripped sequence, uppercase A–Z only |
| `modified_peptide` | text | engine notation or ProForma; kept for provenance |
| `protein_ids` | text | `;`-separated engine protein IDs; decoy and contaminant prefixes are recognised (`input.decoy_prefixes`, `input.contaminant_prefixes`) |
| `q_value` | number in [0, 1] | precursor or peptide q-value (for DIA-NN, the global q-value) |
| `pep` | number in [0, 1] | posterior error probability |
| `score` | number | engine or de novo score, higher is better |
| `is_decoy` | `true`/`false` | decoy row; needs `protein_ids` to pair the decoy with a gene |
| `intensity` | number ≥ 0 | linear precursor or peptide intensity, per run |
| `run` | text | run or sample name |
| `charge` | integer | precursor charge |
| `spectrum_id` | text | scan or spectrum ID, unique within `run` |
| `source` | `database` or `denovo` | default from `input.peptide_source` |
| `search_db` | `same` or `other` | was the search database this genome's own proteins? Default from `input.search_db` |

Columns not in this list are kept, with an info message. Rows are grouped by stripped `peptide`:

- `q_value` and `pep` take the minimum over a peptide's rows, and `score` the maximum;
- `is_decoy` is true only when every row is a decoy;
- `protein_ids` is the union.

Intensities are summed over charges and modified forms within a run. The median across runs is
then taken (`abundance.run_aggregation`). Rows whose stripped length is outside `input.peptide_len`
(7–40 aa by default) are dropped, and preflight reports them under `peptides.content`.

**Decoys help.** Keep the search engine's decoy rows, with their `protein_ids`. GIPPA can then use
picked target–decoy FDR at the gene level, the `decoy_picked` tier. With DIA-NN, search with
`diann --report-decoys`. See
[Interpreting results](interpreting-results.md#detection-tiers).

## Protein FASTA (`--proteins`)

- The ID is the first token of the header. UniProt `sp|ACC|NAME` and `tr|ACC|NAME` resolve to
  `ACC`. A protein group written `A;B` is split.
- Records of 40 aa or shorter go to the peptide path, and records shorter than 7 aa are dropped
  with a warning.
- Allowed letters are the 20 amino acids plus U. X, B, Z and J are counted, and an internal `*`
  warns.

Proteins longer than 40 aa are aligned to the genome with miniprot, then tblastn
(`mapping/protein_hits.tsv`).

## Abundance table (`--abundance`)

The first column is a protein ID: one that matches `--proteins`, a locus tag or a gene ID. Exactly
one numeric value column is used. Extra columns warn and are ignored. `NA` and empty cells are
missing values.

**Scale detection.** With `--abundance-scale auto` (the default), the scale is `log2` when any
value is negative, or when every value is ≤ 50. Otherwise it is `linear`. Give
`--abundance-scale log2` or `--abundance-scale linear` to override the detection. A negative value
is an error only on a linear scale.

**Abundance type.** `--abundance-type` declares what the values are: `ibaq`, `top3`, `maxlfq`,
`intensity` or `unknown`. The default `unknown` means relative values. Missing values are filled
by the policy `--abundance-missing` (`median` by default; or `exclude`, `min`, `zero`, `mean`).
Imputed genes carry `abundance_imputed = true` in `genes.tsv`.

The abundance table is used only when the peptides carry no intensities.

## Existing annotation (`--annotation`)

- **Bakta JSON** (`.json` or `.json.gz`). GIPPA uses it as is and does not run Bakta. The Bakta and
  database versions are read from the JSON. Use this when your search database came from the same
  Bakta annotation, so that the protein IDs match the locus tags.
- **GFF3 or GBFF** from another annotator. Bakta runs with `bakta --regions` on that file's CDS, so the
  gene coordinates match. `annotation/id_map.tsv` maps the source IDs to Bakta locus tags by
  identical CDS coordinates, and `genes.tsv` keeps them in `source_id`.

Preflight recognises the format by its extension: `.json`, `.gff3`, `.gff`, `.gbff` or `.gbk`,
optionally gzipped.

With a Bakta JSON, antiSMASH (the `bgc` module) reads Bakta's GBFF output from next to the JSON:
`<stem>.gbff`, else `<stem>.gbff.gz`, where `<stem>` is the JSON path without `.json` or
`.json.gz`. When neither file exists, antiSMASH is skipped and the bgc layer stays empty, so
preflight warns with `inputs.bakta_gbff` and names the expected paths. The same holds for
`--set annotation.precomputed.bakta_json=<path>`. The rule is not evaluated (info) when `bgc` is
not in `--modules` or `annotation.precomputed.regions` is set. It is a warning only, so `--strict`
turns it into exit code 3.

## Reference peptides (`--reference-peptides`)

A generic peptide table from a database search of the same spectra against this strain's
proteins. It is used for two things:

- the de novo precision of de novo peptides (without it, GIPPA assumes a value and flags the run
  `c_assumed`);
- the error rates of ALT calls in the variant module.

## Replicon depth (`--replicon-depth`)

A TSV with the columns `contig` and `mean_depth`, for example sequencing depth per replicon. It
gives a copy-number offset for plasmid contrasts. Rows flagged `depth_adjusted` used it.

## Coordinates in inputs

Input GFF3 and GBFF files use their formats' own 1-based closed coordinates. GIPPA converts them
once on reading. Every output coordinate convention is listed in
[Outputs › Coordinates](outputs.md#coordinates).
