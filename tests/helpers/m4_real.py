"""The real M2/M3 path for the M4-21 tests (docs/specs/M4-21-m2-m3-conformance.md).

``giap run --executor local`` (the MVP-01 harness of tests/integration/test_mapping_local.py and
test_gene_table_local.py) on a giapsim genome with every precomputed layer and its proteomics:
``parse_bakta``, ``operon_pairs`` (TUs), the regions layer, ``context_layers`` and the M3 mapping
steps (``sixframe_index`` ... ``loci``) run for real. Nothing here builds a contract table.

Layout of ``OUT/work`` (the local executor names every file by its step-output key): the M2
tables are ``work/<key>``, M3's final tables are ``work/mapping/*.tsv`` (C11-C13, as the M4-06
contract names them), the normalised peptides are the staged root ``work/peptides`` and C14 is
``work/annotation/id_map.tsv``. ``real_paths`` is that key mapping (M4-06 lead decision 4).

The simulator's ``config.yaml`` is used unchanged, so the regions arrive through
``annotation.precomputed.regions`` (the planted prophage, BGC, CRISPR and plasmid) -- the run
the spec names ("precomputed layers, local executor") and the one MVP-01 ran.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m4_chain as mc

M3_KEYS = ("peptide_hits", "protein_hits", "unmapped")


def sim100k(out: Path, seed: int | None = None, **proteomics: Any) -> Any:
    """giapsim ``SIM100K`` (the 100-kb genome of every sim100k recording) with M3-19 proteomics."""
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(giapsim.SIM100K, proteomics=ProteomicsConfig(**proteomics))
    return giapsim.simulate(cfg, giapsim.SIM100K_SEED if seed is None else seed, out)


def giap_run(sim: Any, out: Path, db: Path, seed: int, *sets: str) -> Path:
    """``giap run --executor local`` into ``out``; returns ``out/work``. ``sets`` are
    ``--set key=value`` overrides."""
    from click.testing import CliRunner

    from giap.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml),
            "--outdir", str(out), "--seed", str(seed), "--executor", "local"]  # fmt: skip
    for s in sets:
        args += ["--set", s]
    with mc.db_env(db):
        r = CliRunner().invoke(cli, args)
    assert r.exit_code == 0, r.output
    return out / "work"


def real_paths(work: Path) -> dict[str, Path]:
    """M4-06 key -> file of a ``giap run --executor local`` work directory."""
    from giap.io.work_tables import FILE_NAMES

    paths = {k: work / k for k in FILE_NAMES}
    paths.update({k: work / "mapping" / FILE_NAMES[k] for k in M3_KEYS})
    paths["id_map"] = work / "annotation" / "id_map.tsv"
    return paths
