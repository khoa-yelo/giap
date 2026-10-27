# Third-party tools and licences

GIPPA is MIT-licensed. It calls the tools below as separate programs, each in its own pinned
environment or container; it never imports or vendors their code (GIPPA_SPEC.md R5), except
codon-bias (MIT), which runs in-process. Pins, digests and licence checks: docs/adr/ADR-0008-tool-pins.md.

| Tool | Licence | Notes |
|---|---|---|
| Bakta | GPL-3.0 | 1.12.2; separate process |
| UniOP | unlicensed (no licence file) | git 3ed14dfc; fetched at a pinned commit when environments are built; never vendored |
| Promotech | GPL-3.0 | git 56251ad9; separate process |
| Promoter Calculator | GPL-3.0 (or later) | 1.2.4 (PyPI `promotercalculator`); separate process |
| OSTIR | GPL-3.0 | 1.1.3; separate process |
| TerminatorNet | MIT | git 72180309; separate process |
| TransTermHP | GPL-3.0 (or later) | 2.09; separate process |
| BacTermFinder | GPL-3.0 | git 02487044; opt-in; separate process |
| geNOMAD | non-commercial | 1.12.0; separate process; **non-commercial use only** |
| antiSMASH | AGPL-3.0 | 8.0.4; separate process |
| miniprot | MIT | 0.18; separate process |
| BLAST+ | public-domain (NCBI) | 2.16.0; separate process |
| fgsea | MIT | 1.36.2 (bioconda `bioconductor-fgsea`); run by `workflow/bin/gippa_fgsea.R` in R 4.5 (GPL-2) as a separate process (ADR-0038); R and fgsea are never imported |
| codon-bias | MIT | 0.5.0; imported in-process (PyPI) |
| DIA-NN | proprietary (academic licence) | used only by users' own searches; binaries never redistributed |

Container images: `.github/workflows/containers.yml` builds six images for GHCR
(`ghcr.io/khoa-yelo/gippa-<tool>`: promotech, promoter_calculator, terminatornet, bactermfinder, ostir, and
fgsea, described below).
Their digests enter `workflow/conf/containers.config` only through the digest-pin PR that the workflow opens;
until that PR merges, these tools run from their Conda environments. Four of them contain GPL-3.0 tools fetched from
upstream at build time, unmodified: Promotech (https://github.com/BioinformaticsLabAtMUN/PromoTech),
Promoter Calculator (GPL-3.0-or-later, https://pypi.org/project/promotercalculator/1.2.4/), OSTIR
(https://pypi.org/project/ostir/1.1.3/) and BacTermFinder (https://github.com/BioinformaticsLabAtMUN/BacTermFinder);
the corresponding source is available at those links at the pinned versions and commits above, and each image's
Dockerfile is in `workflow/containers/<tool>/`. TerminatorNet (MIT, https://github.com/btjaden/TerminatorNet) is
the fifth. UniOP (unlicensed) is never published.

The `gippa-fgsea` image (`ghcr.io/khoa-yelo/gippa-fgsea`, built by the same workflow from
`workflow/containers/fgsea/Dockerfile`) contains R 4.5.3 (GPL-2, https://www.r-project.org/), fgsea 1.36.2
(MIT, Bioconductor) and data.table 1.18.6.1 (MPL-2.0), installed unmodified from the conda packages of the
committed lock `workflow/envs/locks/fgsea.explicit.txt`, plus GIPPA's `workflow/bin/gippa_fgsea.R`. GIPPA starts it
only as a separate process (`docker run` / `apptainer exec`) and never imports R or fgsea; the corresponding
source is available from CRAN, Bioconductor and conda-forge/bioconda at those versions. Its ref is pinned in
`workflow/conf/images.json` only.

The step images (`ghcr.io/khoa-yelo/gippa-step` and `ghcr.io/khoa-yelo/gippa-step-r`, tagged with the GIPPA
version; ADR-0065) run GIPPA's own steps. `gippa-step` contains GIPPA (MIT) and its permissive Python
dependencies from the hashed lock `workflow/containers/step/requirements.lock` (every licence checked against
`workflow/containers/step/licences.allow` by `scripts/release/check_image_licences.py`), on a Debian
`python:3.12-slim` base with `procps` (GPL-2) and `ca-certificates` as separate system programs. The
`gippa-step-r` image (`workflow/containers/step-r/Dockerfile`) adds R 4.5.3 (GPL-2, https://www.r-project.org/),
fgsea 1.36.2 (MIT, Bioconductor) and data.table 1.18.6.1 (MPL-2.0), installed unmodified from the conda packages
of the committed lock `workflow/envs/locks/fgsea.explicit.txt`. GIPPA starts R only as a separate process
(`subprocess`) and never imports R or fgsea; the corresponding source is available from CRAN, Bioconductor and
conda-forge/bioconda at those versions, and Debian's from https://sources.debian.org.

Data resources: KEGG (academic use), RegulonDB (not redistributed), Gene Ontology (CC BY 4.0), NCBI COG 2020
(NCBI FTP).

Contaminant library: the universal protein contaminant FASTA of Frankenfield et al. 2022
(HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics, commit bddf6da9); no licence
file found in the repository (checked 2026-10-15), so it is downloaded by `gippa setup --contaminants` and
never vendored (R11).

## Inlined in the report

| Code | Licence | Notes |
|---|---|---|
| plotly.js | MIT | 4.1.1 (from the plotly 7.1.0 Python package); inlined once in `report.html` when interactive figures are shown; the page footer names it |

## Ported code

| Source | Licence | Notes |
|---|---|---|
| fgsea 1.36.2 | MIT | Bioconductor 3.22 tarball (git b7862b7); `calcGseaStat`, `fgseaSimpleImpl`, `calcGseaStatCumulative(Batch)`, `preparePathways` and the simple stage of `fgseaMultilevel` ported to Python/numba in `src/gippa/stats/fgsea.py` and `src/gippa/stats/_fgsea_kernels.py`; no fgsea file is vendored (docs/methods/fgsea-port.md) |

fgsea licence (MIT):

    Copyright (c) 2016-2019 Alexey Sergushichev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
