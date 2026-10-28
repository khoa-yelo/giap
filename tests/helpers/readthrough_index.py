"""The small two-contig world of the M7-16 index-format tests (review change 7).

``world()`` is deterministic: contig ``a`` (linear, 3,000 nt) and ``b`` (circular, 1,200 nt) from
``numpy.random.default_rng(16)``, a handful of annotated CDS/sORF/pseudogene features and a
peptide table with one ``database`` and one ``denovo`` query (so ``sixframe_index`` builds both
decoy types). ``run_sixframe_index`` runs the registered step on it.

``tests/data/mapping/readthrough/sixframe_format1.npz`` is the **format-1** work file that
``sixframe_index`` wrote for this world on main at 508e86e0 (before M7-16a), with ``seed=7``:
``python tests/data/mapping/readthrough/make_format1_fixture.py`` (only meaningful on a checkout
whose ``save_index`` still writes format 1). Index format 2 must keep its decoy frames
byte-identical (spec rule 3). ``sixframe_format2.npz`` next to it is the format-2 file of the same
world *with* genes (main 1f7203d3, before M11-16b; ``make_format2_fixture.py``): the table-11
format-3 index must equal it apart from ``meta.format`` and ``meta.translation_table``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gippa.model import Replicon
from gippa.model.genes import GeneModel
from helpers.readthrough_toys import CODON, gene, model, write_annotation, write_genome

FIXTURE = Path(__file__).resolve().parents[1] / "data/mapping/readthrough/sixframe_format1.npz"
SEED = 7
_STOPS = ("TAA", "TAG", "TGA")


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def world() -> tuple[dict[str, str], dict[str, Replicon], GeneModel]:
    rng = np.random.default_rng(16)
    a = list("".join(rng.choice(list("ACGT"), size=3000)))
    b = list("".join(rng.choice(list("ACGT"), size=1200)))

    def put(seq: list[str], at: int, codon: str) -> None:
        n = len(seq)
        for i, c in enumerate(codon):
            seq[(at + i) % n] = c

    put(a, 237, "".join(CODON[x] for x in "MKTAYLAKQRQLSFVKSHFS"))  # stop-free before A1's stop
    put(a, 297, "TAA")  # A1 + [0, 300)
    put(a, 300, "".join(CODON[x] for x in "GSHMWEPCNL"))  # and after it
    put(a, 897, "TGA")  # A2 + [600, 900) and A3 + [750, 900): shared stop
    put(a, 1500, _rc("TAG"))  # A4 - [1500, 1800)
    put(a, 2097, "TAA")  # A5 sORF + [2001, 2100)
    put(a, 2497, "TGA")  # A6 pseudogene + [2200, 2500): never a site
    put(a, 2897, "GCT")  # A7 + [2600, 2900): last codon not a stop
    put(b, 1199, "TAG")  # B1 + [1000, 1202): stop crosses the origin
    put(b, 3, _rc("TAA"))  # B2 - [3, 153): stop at [3, 6), two copies
    genome = {"a": "".join(a), "b": "".join(b)}
    reps = {"a": Replicon("a", 3000, "linear"), "b": Replicon("b", 1200, "circular")}
    genes = [
        gene("A1", "a", 0, 300, "+"), gene("A2", "a", 600, 900, "+"),
        gene("A3", "a", 750, 900, "+"), gene("A4", "a", 1500, 1800, "-"),
        gene("A5", "a", 2001, 2100, "+", "sorf"), gene("A6", "a", 2200, 2500, "+", "pseudogene"),
        gene("A7", "a", 2600, 2900, "+"), gene("B1", "b", 1000, 1202, "+"),
        gene("B2", "b", 3, 153, "-"),
    ]  # fmt: skip
    return genome, reps, model(reps, genes)


def peptides_table(path: Path) -> Path:
    df = pd.DataFrame(
        {
            "peptide": ["PEPTLDEKAR", "SAMPLEWLDK"],
            "source": ["database", "denovo"],
            "q_value": [0.001, np.nan],
            "pep": [0.01, np.nan],
            "score": [np.nan, 0.9],
            "is_decoy": [False, False],
            "search_db": ["other", "other"],
        }
    )
    df.to_parquet(path, index=False)
    return path


def inputs(tmp: Path) -> dict[str, Path]:
    genome, _reps, gm = world()
    paths = write_annotation(gm, tmp / "ann")
    return {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": paths["replicons"],
        "genes": paths["genes"],
        "peptides": peptides_table(tmp / "peptides.parquet"),
    }


def run_sixframe_index(tmp: Path, with_genes: bool = True) -> Path:
    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    ins = inputs(tmp)
    spec = reg.get("sixframe_index")
    given = {k: v for k, v in ins.items() if k in spec.inputs and (with_genes or k != "genes")}
    out = tmp / "sixframe.npz"
    reg.run_step("sixframe_index", resolve_config(overrides=[f"seed={SEED}"]), given,
                 {"sixframe": out})  # fmt: skip
    return out
