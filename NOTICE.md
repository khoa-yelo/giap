# Third-party tools and licences

GIAP is MIT-licensed. It calls the tools below as separate programs, each in its own pinned
environment or container; it never imports or vendors their code (GIAP_SPEC.md R5), except
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
| fgsea | MIT | 1.36.2 (bioconda `bioconductor-fgsea`); run by `workflow/bin/giap_fgsea.R` in R 4.5 (GPL-2) as a separate process (ADR-0038); R and fgsea are never imported |
| codon-bias | MIT | 0.5.0; imported in-process (PyPI) |
| DIA-NN | proprietary (academic licence) | used only by users' own searches; binaries never redistributed |

Container images: `.github/workflows/containers.yml` publishes five images to GHCR
(`ghcr.io/khoa-yelo/giap-<tool>`: promotech, promoter_calculator, terminatornet, bactermfinder, ostir),
referenced by digest in `workflow/conf/containers.config`. Four of them contain GPL-3.0 tools fetched from
upstream at build time, unmodified: Promotech (https://github.com/BioinformaticsLabAtMUN/PromoTech),
Promoter Calculator (GPL-3.0-or-later, https://pypi.org/project/promotercalculator/1.2.4/), OSTIR
(https://pypi.org/project/ostir/1.1.3/) and BacTermFinder (https://github.com/BioinformaticsLabAtMUN/BacTermFinder);
the corresponding source is available at those links at the pinned versions and commits above, and each image's
Dockerfile is in `workflow/containers/<tool>/`. TerminatorNet (MIT, https://github.com/btjaden/TerminatorNet) is
the fifth. UniOP (unlicensed) is never published.

Data resources: KEGG (academic use), RegulonDB (not redistributed), Gene Ontology (CC BY 4.0), NCBI COG 2020
(NCBI FTP).

Contaminant library: the universal protein contaminant FASTA of Frankenfield et al. 2022
(HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics, commit bddf6da9); no licence
file found in the repository (checked 2026-10-15), so it is downloaded by `giap setup --contaminants` and
never vendored (R11).

## Inlined in the report

| Code | Licence | Notes |
|---|---|---|
| plotly.js | MIT | 4.1.1 (from the plotly 7.1.0 Python package); inlined once in `report.html` when interactive figures are shown; the page footer names it |

## Ported code

| Source | Licence | Notes |
|---|---|---|
| fgsea 1.36.2 | MIT | Bioconductor 3.22 tarball (git b7862b7); `calcGseaStat`, `fgseaSimpleImpl`, `calcGseaStatCumulative(Batch)`, `preparePathways` and the simple stage of `fgseaMultilevel` ported to Python/numba in `src/giap/stats/fgsea.py` and `src/giap/stats/_fgsea_kernels.py`; no fgsea file is vendored (docs/methods/fgsea-port.md) |

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
