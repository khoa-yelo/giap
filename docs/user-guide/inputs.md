# Inputs

GIAP reads one assembled genome and proteomics evidence from the search you already ran. It never
searches spectra itself. Run `giap preflight` to check every input in seconds before a run (see
[Running GIAP](run.md)).

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

`giap preflight` takes `--genome`, `--proteins`, `--peptides`, `--abundance` and `--annotation`.
`giap run` takes them all. Without `--peptides` and `--proteins`, preflight stops with the error
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
- **Genetic code.** This release supports genetic code 11 (bacteria and archaea) only. Any other
  `annotation.translation_table` stops the run in preflight with the input error
  `config.translation_table`. Support for genetic code 4 (Mycoplasma) is planned.
  A genome left at table 11 that probably uses code 4 (TGA read as tryptophan) gets the preflight
  warning `genome.genetic_code`: long open reading frames (at least 300 codons, six frames) cover
  less than 0.55 times as many bases with TGA as a stop as with TGA as a sense codon (evaluated
  with at least 20 such frames). It is a warning only, so `--strict` turns it into exit code 3.

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

## Peptide table (`peptides.tsv`)

The generic peptide table is a UTF-8 TSV with one header row. `giap convert` writes it from
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

**Decoys help.** Keep the search engine's decoy rows, with their `protein_ids`. GIAP can then use
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

- **Bakta JSON** (`.json` or `.json.gz`). GIAP uses it as is and does not run Bakta. The Bakta and
  database versions are read from the JSON. Use this when your search database came from the same
  Bakta annotation, so that the protein IDs match the locus tags.
- **GFF3 or GBFF** from another annotator. Bakta runs with `bakta --regions` on that file's CDS, so the
  gene coordinates match. `annotation/id_map.tsv` maps the source IDs to Bakta locus tags by
  identical CDS coordinates, and `genes.tsv` keeps them in `source_id`.

Preflight recognises the format by its extension: `.json`, `.gff3`, `.gff`, `.gbff` or `.gbk`,
optionally gzipped.

## Reference peptides (`--reference-peptides`)

A generic peptide table from a database search of the same spectra against this strain's
proteins. It is used for two things:

- the de novo precision of de novo peptides (without it, GIAP assumes a value and flags the run
  `c_assumed`);
- the error rates of ALT calls in the variant module.

## Replicon depth (`--replicon-depth`)

A TSV with the columns `contig` and `mean_depth`, for example sequencing depth per replicon. It
gives a copy-number offset for plasmid contrasts. Rows flagged `depth_adjusted` used it.

## Coordinates in inputs

Input GFF3 and GBFF files use their formats' own 1-based closed coordinates. GIAP converts them
once on reading. Every output coordinate convention is listed in
[Outputs › Coordinates](outputs.md#coordinates).
