# GIPPA

GIPPA (genome integrated prokaryotic proteomics analysis) takes one assembled bacterial genome plus
proteomics evidence and asks: what does genome context explain about which proteins are detected,
and how abundant they are?

GIPPA runs after your proteomics search; it does not search spectra itself. It:

- places peptides and proteins on the genome, with its own error control (placement FDR);
- builds one gene table keyed by Bakta locus tags, with detection tiers and gene q-values;
- tests genome context: regions, autocorrelation along the chromosome, abundance islands,
  operons, function and the sequence determinants of abundance;
- calls REF/ALT residues at covered codons and finds novel ORFs and start sites;
- writes standard tables with schemas, plots with their data, and a self-contained HTML report.

Every claim carries a q-value, an FDR, an FWER or the label `unestimated`, and GIPPA says
"detected", never "expressed".

## Start here

- [Install](get-started/install.md) the `gippa` command, Nextflow and a container profile.
- [Set up the databases](get-started/setup-databases.md) once with `gippa setup`.
- Run the [Quickstart](get-started/quickstart.md) on public E. coli K-12 inputs.

## User guide

- [Inputs](user-guide/inputs.md), [Running GIPPA](user-guide/run.md) and
  [Outputs](user-guide/outputs.md).
- [Interpreting results](user-guide/interpreting-results.md): q-values, detection tiers and kinds
  of placement.
- [Plots and report](user-guide/plots-and-report.md), [Comparing runs](user-guide/compare.md),
  [Building a search database](user-guide/makedb.md) and
  [Converting search-engine output](user-guide/convert.md).
- [Licences](user-guide/licences.md): geNOMAD is non-commercial, KEGG is academic-use; read this
  before installing.
- [Platforms](user-guide/platforms.md) and [Runtimes](user-guide/runtimes.md).

## Reference

- [CLI](reference/cli.md), [configuration](reference/configuration.md),
  [output tables](reference/tables.md) and [flags](reference/flags.md), generated from the code and
  the schemas.
- [Methods](methods/index.md): one note per statistical method.
- [Decisions](adr/index.md): the architecture decision records.
