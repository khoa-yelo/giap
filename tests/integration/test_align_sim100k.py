"""M7-01 conformance on the merged sim100k producer chain (docs/specs/M7-01-residue-alignment.md,
acceptance: "every row of the real peptide_hits.tsv aligns without error").

The genome is gippasim ``SIM100K`` (seed ``SIM100K_SEED``) with a strain layer (SNVs and short
indels, M3-18): the peptides are digested from the *sample* strain and mapped to the *reference*
genome, so the tolerant route (substitutions, frameshift splits) has rows to align. Preflight
and the merged mapping steps run through ``run_step`` (``helpers.mapping_steps.drive``); the
test reads the published ``peptide_hits.tsv`` (1-based closed -> 0-based half-open, as the
``variants`` step converts it, rule 1), the ``hits_tolerant`` work table, the normalised genome
and ``replicons.parquet``, and aligns every row. ``$GIPPA_DB_DIR`` is a test database (never the
real one).

Checks: one alignment per row; no ``ValueError`` (the recomputation agrees with the mapper on
every accepted hit); each alignment's cost equals the row's ``n_diff`` (``+1`` with a junction
event exactly on frameshift rows); the ops replay the peptide.
"""

from __future__ import annotations

import dataclasses
import importlib
from pathlib import Path

import pandas as pd
import pytest
from helpers.mapping_steps import drive, fresh_steps, key_like, module2_files, preflight, roots_of

from gippa.config import resolve_config
from gippa.io.fasta import read_genome_fasta
from gippa.model import Replicon

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SIM100K, SIM100K_SEED, simulate
    from gippasim.strain import StrainConfig

    tmp = Path(tmp_path_factory.mktemp("align_sim100k"))
    cfg_sim = dataclasses.replace(
        SIM100K,
        proteomics=ProteomicsConfig(),
        strain=StrainConfig(snv_rate_per_kb=5.0, indel_rate_per_kb=0.5, indel_len=(1, 1)),
        write_precomputed=False,
    )
    res = simulate(cfg_sim, SIM100K_SEED, tmp / "sim")
    assert res.reference_fasta is not None and res.reference_json is not None
    ref = dataclasses.replace(res, genome_fasta=res.reference_fasta, bakta_json=res.reference_json)
    run = preflight(ref, tmp / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={SIM100K_SEED}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots_of(run))
    files = module2_files(tmp / "work", avail)
    ph = pd.read_csv(files["peptide_hits.tsv"], sep="\t", keep_default_na=False, na_values=["NA"])
    ph["start"] = ph["start"].astype("int64") - 1  # §6.5: 1-based closed -> 0-based half-open
    ph["end"] = ph["end"].astype("int64")
    ht = pd.read_parquet(avail[key_like(list(avail), "hits_tolerant")])
    genome = {r.id: r.sequence for r in read_genome_fasta(run.genome)}
    rows = pd.read_parquet(avail[key_like(list(avail), "replicons")]).to_dict("records")
    reps = {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}
    return {"ph": ph, "ht": ht, "genome": genome, "reps": reps}


def test_every_accepted_hit_aligns(chain: dict[str, object]) -> None:
    mod = importlib.import_module("gippa.variants.align")
    ph = chain["ph"]
    assert isinstance(ph, pd.DataFrame)
    tolerant = ph[ph["route"] == "tolerant"]
    assert len(tolerant) > 0, "the fixture must exercise the tolerant route"
    got = mod.align_hits(ph, chain["ht"], chain["genome"], chain["reps"])
    assert len(got) == len(ph)
    for a, row in zip(got, ph.to_dict("records"), strict=True):
        assert (a.peptide, a.contig, a.strand) == (row["peptide"], row["contig"], row["strand"])
        assert "".join(op.alt for op in a.ops) == row["peptide"]
        cost = sum(op.kind not in ("match", "il") for op in a.ops)
        fs = isinstance(row["diff_types"], str) and "frameshift" in row["diff_types"]
        assert (a.frameshift is not None) == fs
        assert cost + int(fs) == int(row["n_diff"]), row["peptide"]
        codons = [c for op in a.ops for c in op.codons]
        assert min(codons) == row["start"] and max(codons) + 3 == row["end"], row["peptide"]
