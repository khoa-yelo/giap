# gippa demo data: provenance

Everything here is **synthetic**: gippasim (tests/sim/gippasim in the GIPPA repository)
simulated it from the SIM100K configuration with seed 20261014 (`SIM100K_SEED`). No file comes
from a real organism, a database or a third party. `gippa demo` reads these files; nothing else does.

| File | What | Written by |
|---|---|---|
| `sim100k.fna.gz` | SIM100K genome: a 100 kb circular chromosome and an 8 kb circular plasmid (byte-identical to `tests/data/genomes/sim100k.fna.gz`) | gippasim; copied by `scripts/docs/make_demo_layers.py` |
| `sim100k_peptides.tsv.gz` | simulated peptide table, one proteomics run (byte-identical to `docs/tutorials/data/sim100k_peptides.tsv.gz`) | `scripts/docs/make_demo_peptides.py`; copied by `make_demo_layers.py` |
| `layers/bakta.json.gz` | Bakta-like annotation JSON (`annotation.precomputed.bakta_json`) | `scripts/docs/make_demo_layers.py` |
| `layers/operon_pairs.tsv.gz`, `promoters`, `tx_rates`, `terminators`, `rbs`, `regions` | precomputed tool layers (`annotation.precomputed.<key>`) | `scripts/docs/make_demo_layers.py` |
| `layers/rbs_novel.tsv.gz` | novel-ORF RBS layer `rbs_novel` (ADR-0024) | `scripts/docs/make_demo_layers.py` |

Layers are gzip (level 9) with no file name, `mtime=0` and the header OS byte forced to `0xff`
(RFC 1952 "unknown"; Python 3.12 would write `0x03`), so they are reproducible byte for byte on
any Python version with the same zlib.
Regenerate with `python scripts/docs/make_demo_layers.py`; CI runs it with `--check`.
Licence: MIT, as GIPPA.
