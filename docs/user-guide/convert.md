# Converting search-engine output

`giap convert` turns a search engine's output into the generic peptide table that `giap run
--peptides` reads (see [Inputs](inputs.md#peptide-table-peptidestsv)). Every option is listed in
the [CLI reference](../reference/cli.md#giap-convert).

```bash
giap convert --from diann report.parquet --out peptides.tsv --abundance-out abundance.tsv
giap run --genome genome.fna --peptides peptides.tsv --outdir out
```

DIA-NN peptides carry `intensity`, so `giap run` computes gene abundance from the peptides. The
abundance table written by `--abundance-out` is then not used, even when you pass it with
`--abundance`. It is used only for peptides without intensities (see
[Inputs › Abundance table](inputs.md#abundance-table-abundance)). Passing it with
`--abundance-type maxlfq` makes sense only in that case.

## Options

| Option | Meaning |
|---|---|
| `INPUTS` | one or more input files; several files are read in the order given |
| `--from` | the engine: `diann`, `casanovo`, `instanovo`, `cascadia` or `generic` (required) |
| `--out` | the peptide table to write (required) |
| `--abundance-out` | also write a protein abundance table (DIA-NN only) |
| `--source` | `database` or `denovo`: the default `source` for `--from generic` |
| `--search-db` | `same` or `other`: was the search database this genome's own proteins? |

## Engines

| `--from` | Status | Input | Notes |
|---|---|---|---|
| `diann` | implemented | DIA-NN 2.x `report.parquet` or `report.tsv` | `source` = database |
| `cascadia` | implemented | Cascadia 0.0.7 result file (`.ssl`) | `source` = denovo; no q-value or PEP is derived from `score` |
| `generic` | implemented | a generic peptide TSV | validates and normalises it, and reports problems |
| `casanovo` | not implemented: exit code 5 | — | — |
| `instanovo` | not implemented: exit code 5 | — | — |

`--from casanovo` and `--from instanovo` are accepted choices, but GIAP 1.0 has no converter for
them. Both print "not implemented" and exit with code 5. Until a converter exists, write a generic
peptide table yourself (see below).

With `--from diann`, `--source` may only be `database`; with `--from cascadia`, only `denovo`.
Any other value exits with code 1.

### DIA-NN

| Generic column | DIA-NN column |
|---|---|
| `peptide` | `Stripped.Sequence` (modification syntax stripped) |
| `modified_peptide` | `Modified.Sequence` |
| `protein_ids` | `Protein.Ids` |
| `q_value` | `Global.Q.Value`, else `Q.Value` |
| `pep` | `PEP` |
| `is_decoy` | `Decoy` (present only when DIA-NN ran with `diann --report-decoys`) |
| `intensity` | `Precursor.Normalised`, else `Precursor.Quantity` (0 becomes missing) |
| `run`, `charge` | `Run`, `Precursor.Charge` |

`--abundance-out` writes one row per `Protein.Group`, with a column `maxlfq`. The value is the
median over runs of `PG.MaxLFQ`. Missing values, zeros and decoy rows are excluded, and groups with
no value are dropped. The command then prints a reminder to use `--abundance-type maxlfq` with
the table. That reminder applies only when the table is used (see above). With any
other engine, `--abundance-out` is an error (`no_abundance`, exit code 1).

DIA-NN is proprietary software under its own licence. GIAP never redistributes it, and you run it
yourself (see [Licences](licences.md)). Keep the decoy rows, so that GIAP can use the
`decoy_picked` detection tier.

### Cascadia

`peptide` is the stripped `sequence` and `modified_peptide` the original `sequence`. `score` is
copied from `score`, `run` is the stem of `file`, and `charge` is copied from `charge`.
`spectrum_id` is `cascadia:<scan>`. Note that Cascadia's `scan` is a prediction index, not a scan
number.

### Generic

The generic converter checks a table you wrote. It strips modifications, such as `(UniMod:35)`,
`[+15.995]`, `[UNIMOD:35]`, `+15.995` or `n[42]`, and lower-case residues. It then validates every
column. `--source` and `--search-db` fill the `source` and `search_db` columns where they are
empty. Booleans are canonicalised, and columns that are not in the generic table are dropped with
the warning `unknown_column`. (`giap run` keeps unknown columns instead.)

## Messages and exit codes

Problems are printed as `giap convert: <severity>: <code>: <message>`, with example row numbers.
Any `error` exits with code 1, and no table is written. Warnings and info messages do not stop the
conversion. For the full list of exit codes, see [Running GIAP](run.md#exit-codes).
