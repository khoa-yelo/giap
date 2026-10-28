# Licences

GIPPA itself is MIT-licensed. The tools and databases it uses keep their own licences. Some of
these restrict who may use them. Read this page before you install the databases or publish
results.

## What you need to know

- **geNOMAD is for non-commercial use only.** geNOMAD (module `mobile`) is licensed for
  non-commercial use. For commercial use, switch the module off: run with
  `--modules core,bgc,function`, and install with `gippa setup --no-genomad`.
- **KEGG is free for academic use only.** `gippa setup` downloads the KEGG KO-to-pathway
  links and pathway names from `rest.kegg.jp` by default (`--no-kegg` leaves them out). Academic
  users may use them freely. Other users need a licence from KEGG. Runs with the `function` module
  contain KEGG pathway names, and the manifest records the KEGG retrieval date.
- **Where the terms are recorded.** `gippa setup` prints the geNOMAD and KEGG terms before it
  downloads anything. Each run lists the terms of what it used in `run_manifest.json`
  (`licence_notes`, e.g. `geNOMAD: non-commercial`, `KEGG: academic use`), and the report's
  Methods section repeats them.
- **GPL and AGPL tools run as separate programs.** Bakta, Promotech, Promoter Calculator, OSTIR,
  TransTermHP and BacTermFinder (GPL-3.0), and antiSMASH (AGPL-3.0), each run in their own pinned
  container or Conda environment. GIPPA calls them as separate processes and never imports or
  vendors their code. A CI check (`scripts/ci/check_imports.py`) guards this. R and the fgsea runner run
  the same way ([ADR-0038](../adr/ADR-0038-gsea-runs-r-fgsea.md)).
- **UniOP has no licence.** UniOP's repository has no licence file. GIPPA fetches it from its
  upstream repository at a pinned commit (`3ed14dfc`) when the tool environments are built. GIPPA
  never redistributes it, and no UniOP container image is published.
- **DIA-NN is proprietary.** DIA-NN is used only for your own searches, under your own DIA-NN
  licence (NOTICE.md: "proprietary (academic licence)"). GIPPA never redistributes DIA-NN binaries, and
  `gippa convert --from diann` reads only its report files.
- **RegulonDB is not redistributed.** GIPPA's developers use RegulonDB only to calibrate cut-offs.
  `gippa setup` does not download it.
- **The contaminant library has no licence file.** The universal protein contaminant FASTA
  (Frankenfield et al. 2022) is therefore downloaded by `gippa setup` (on by default) at a pinned
  commit and never vendored.

## Tools

| Tool | Licence | How GIPPA uses it |
|---|---|---|
| Bakta | GPL-3.0 | separate process |
| UniOP | no licence | separate process; fetched at a pinned commit, never redistributed |
| Promotech | GPL-3.0 | separate process |
| Promoter Calculator | GPL-3.0 (or later) | separate process |
| OSTIR | GPL-3.0 | separate process |
| TerminatorNet | MIT | separate process |
| TransTermHP | GPL-3.0 (or later) | separate process |
| BacTermFinder | GPL-3.0 | separate process; opt-in (`--add-bactermfinder`) |
| geNOMAD | **non-commercial** | separate process; module `mobile` |
| antiSMASH | AGPL-3.0 | separate process; module `bgc` |
| miniprot | MIT | separate process |
| BLAST+ | public domain (NCBI) | separate process |
| fgsea (R) | MIT; R is GPL-2 | separate process (`workflow/bin/gippa_fgsea.R`); parts of fgsea are also ported to Python under its MIT licence |
| codon-bias | MIT | imported (the only in-process tool) |
| DIA-NN | proprietary | your own search only; never redistributed |

Versions, pins and image digests are in `NOTICE.md` and
[ADR-0008](../adr/ADR-0008-tool-pins.md).

## Data resources

| Resource | Terms |
|---|---|
| KEGG | **academic use** free; others need a KEGG licence |
| Gene Ontology | CC BY 4.0 |
| NCBI COG 2020 | NCBI FTP (public) |
| Bakta DB, geNOMAD DB | see each database's Zenodo record (14916843, 14886553) |
| RegulonDB | not redistributed; not downloaded by users |
| Contaminant library | no licence file; downloaded, never vendored |

## NOTICE.md

The notice file of the repository, reproduced verbatim:

--8<-- "NOTICE.md"
