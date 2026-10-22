# ADR-0039: bioconductor.org as a download source for the R golden environment
- Status: accepted-by-khoa
- Date: 2026-10-22
- Spec sections: §3.2, R6, §11.5
- Needs Khoa: yes (download source; decided by Khoa on 2026-10-22)

## Context
The §11.5 golden environment installs goseq from bioconda (`scripts/golden/golden_r.explicit.txt`, PR #192).
goseq depends on two Bioconductor data packages, `bioconductor-genelendatabase` and `bioconductor-go.db`. Their
bioconda post-link scripts download the package tarballs from bioconductor.org and check them against md5
sums. bioconductor.org is not listed in §3.2, so #192 was held under R6.

## Decision
- `bioconductor.org` is added to the §3.2 sources, for one purpose only: the tarballs that bioconda
  post-link scripts fetch while building the R golden environment (§11.5, `golden-r.yml`,
  `scripts/slurm/golden_r.sbatch`).
- Everything else in R6 still applies:
  - The artefacts are pinned. The conda packages are fixed by the explicit lock, with md5 sums.
  - The post-link checksum must pass.
  - No mirrors and no disabled TLS.
  - A failure is logged under PROGRESS.md › Blocked.
- The pipeline's own environments never use this source. ADR-0038's fgsea environment needs no Bioconductor
  data packages, and the pipeline does not call goseq.
- Nothing fetched from bioconductor.org is committed; only the golden outputs are.

## Evidence
- The goseq golden data in #192 (Slurm 1868639) were produced this way. The md5 checks of both post-link
  downloads passed.

## Consequences
- #192 is unblocked. The goseq Wallenius golden comparison (§7.7.5) validates the scipy ORA (ADR-0038
  item 5).
- Any other use of bioconductor.org needs a new ADR.
