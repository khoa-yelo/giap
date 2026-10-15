# Real Bakta 1.12.2 (DB v6.0, full) output excerpts (unit-test fixtures, M2-11)

Source: the M2-10 full-database recordings under `$GIAP_ROOT/recordings/bakta_full/` (Bakta run as a
process on FarmShare, 2026-10-14). Used as realistic `--annotation` GFF3/GBFF inputs from "another
annotator" whose CDS coordinates and locus tags are known.

| file | origin | sha256 of the full source file |
|---|---|---|
| `sim100k_bakta.gff3` | `sim100k/bakta/sim100k.gff3`, every line before `##FASTA`, verbatim | `f31f08227809f931cc08017c829f1bf1dfbe01bbd2ddbe569ad64b52dce7dec1` |
| `sim100k_bakta.tsv` | `sim100k/bakta/sim100k.tsv`, verbatim | `7f2371e81837e1992715d185a4d6c9b092daf3ec1e3996a8c7e9f5176de29f4c` |
| `sim100k_bakta.gbff.gz` | `sim100k/bakta/sim100k.gbff`, verbatim, `gzip -9 -n` | `c028b1c43b22260713a77796c6ce79d855773b682a68c2f4a2fe36124481c6fd` |
| `k12_excerpt.gff3` | `k12_mg1655/bakta/k12_mg1655.gff3`, lines 1-9 plus the CDS lines of `KOIAIB_00744`..`KOIAIB_00746` | `f6e22eb8f6fdc60cc3e4f0dfeab4e08311107cbdb9a63337b9647f427f6a5585` |
| `k12_excerpt.gbff` | `k12_mg1655/bakta/k12_mg1655.gbff`, lines 1-30 (header, `source`) and 21654-21727 (gene/CDS of `KOIAIB_00744`..`00746`), then an empty `ORIGIN` and `//` | `bc0ad5c0ece48585e913f96ec230f596c9b75238dbf2417451348768b71871d7` |

sim100k is giapsim `SIM100K` / `SIM100K_SEED` (`chr1` 100,000 bp circular, `plasmid_1` 8,000 bp
circular); 86 CDS. K-12 is `U00096.3` (4,641,652 bp, written `linear` by Bakta). `KOIAIB_00745` is a
Bakta pseudogene (`pseudogene=unitary`, location `complement(716388..>716597)`).
