"""Bakta light with ``--regions`` keeps the user's CDS coordinates (M2-11 acceptance).

docs/specs/M2-11-user-annotation-id-map.md; GIAP_SPEC.md §5.3 ("Bakta runs with ``--regions`` on
that file's CDS, so the gene coordinates match"), §7.1.1 (``--regions`` accepts CDS only).

Runs Bakta 1.12.2 (light DB) as a process (R5) from its micromamba env on the giapsim SIM100K
genome, with the simulator's GFF3 as ``--annotation``; at least 99 % of the id map rows must be
``identical``.
"""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
_GIAP_REPO = Path(os.environ.get("GIAP_REPO", REPO))
_ROOT = Path(os.environ.get("GIAP_ROOT", _GIAP_REPO / "local"))
PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", _ROOT / "micromamba"))
DB = Path(os.environ.get("GIAP_DB_DIR", _ROOT / "db"))

pytestmark = [pytest.mark.needs_conda, pytest.mark.needs_bakta_db, pytest.mark.slow]


def test_bakta_regions_keeps_given_cds(tmp_path: Path) -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.annotate.bakta_cmd import bakta_command
    from giap.annotate.tools import TOOLS
    from giap.annotate.user_annotation import (
        build_id_map,
        detect_format,
        read_user_cds,
        write_bakta_regions,
        write_id_map,
    )
    from giapsim import SIM100K, SIM100K_SEED, simulate

    env = TOOLS["bakta"].env_name
    if not (PREFIX / "envs" / env).is_dir():
        pytest.skip(f"{env} not built under {PREFIX / 'envs'}")
    db = DB / "bakta" / "db-light"
    if not db.is_dir():
        pytest.skip(f"no Bakta light DB at {db}")

    sim = simulate(SIM100K, SIM100K_SEED, tmp_path / "sim")
    reps = parse_bakta_json(sim.bakta_json).model.replicons
    assert sim.annotation_gff3 is not None
    assert detect_format(sim.annotation_gff3) == "gff3"
    cds = read_user_cds(sim.annotation_gff3, "gff3", reps, {k: k for k in reps})
    assert len(cds) > 50
    regions = tmp_path / "regions.gff3"
    write_bakta_regions(regions, cds, list(reps.values()))

    out = tmp_path / "bakta"
    argv = bakta_command(
        genome=sim.genome_fasta,
        outdir=out,
        prefix="sim100k",
        db=db,
        threads=4,
        all_circular=all(r.topology == "circular" for r in reps.values()),
        complete="auto",
        replicons_tsv=tmp_path / "replicons.tsv",
        regions_gff3=regions,
        extra_args=["--skip-plot"],
    )
    proc = subprocess.run(
        f"micromamba run -n {env} {shlex.join(argv)}",
        shell=True,
        capture_output=True,
        text=True,
        timeout=3600,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-4000:]

    genes = parse_bakta_json(out / "sim100k.json", reps).model.genes
    rows = build_id_map(cds, genes)
    write_id_map(tmp_path / "id_map.tsv", rows)
    identical = sum(r.match == "identical" for r in rows)
    assert len(rows) == len(cds)
    assert identical / len(rows) >= 0.99, f"{identical}/{len(rows)} identical"
