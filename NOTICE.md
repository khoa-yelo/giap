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
| codon-bias | MIT | 0.5.0; imported in-process (PyPI) |
| DIA-NN | proprietary (academic licence) | used only by users' own searches; binaries never redistributed |

Data resources: KEGG (academic use), RegulonDB (not redistributed), Gene Ontology (CC BY 4.0).
