"""Live miniprot 0.18 on the 100-kb sim genome equals the recorded fixture (issue M3-12).

Needs the built ``miniprot`` micromamba env (``needs_conda``; built by
``sbatch scripts/slurm/build_envs.sbatch miniprot``) and the recorded run
``tests/data/recorded/miniprot/0.18/sim100k/`` (``pending_sample`` until then). Runs the two
§7.1.1 commands from the Python builders, then compares the parsed tables.
"""

from __future__ import annotations

import gzip
import os
import shlex
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from gippa.config.model import ProteinsConfig

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "miniprot" / "0.18" / "sim100k"
_GIPPA_REPO = Path(os.environ.get("GIPPA_REPO", REPO))
_ROOT = Path(os.environ.get("GIPPA_ROOT", _GIPPA_REPO / "local"))
PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", _ROOT / "micromamba"))

pytestmark = [
    pytest.mark.needs_conda,
    pytest.mark.slow,
    pytest.mark.pending_sample,
    pytest.mark.skipif(not RECORDED.is_dir(), reason=f"pending_sample: no recording in {RECORDED}"),
]


def _plain(path: Path, tmp: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = tmp / path.name.removesuffix(".gz")
    out.write_bytes(gzip.decompress(path.read_bytes()))
    return out


def _run(env: str, argv: list[str], stdout: Path | None = None) -> None:
    cmd = f"micromamba run -n {env} {shlex.join(argv)}"
    with open(stdout, "w") if stdout else open(os.devnull, "w") as out:
        proc = subprocess.run(
            cmd, shell=True, stdout=out, stderr=subprocess.PIPE, text=True, timeout=1800
        )
    assert proc.returncode == 0, proc.stderr


def test_live_equals_recorded(tmp_path: Path) -> None:
    from gippa.annotate.bakta import parse_bakta_json
    from gippa.annotate.tools import TOOLS
    from gippa.io.protein_fasta import read_protein_fasta
    from gippa.mapping.miniprot import miniprot_align_cmd, miniprot_index_cmd, parse_miniprot_gff
    from gippa.mapping.protein_assign import assign_alignments
    from gippasim import SIM100K, SIM100K_SEED, simulate

    env = TOOLS["miniprot"].env_name
    if not (PREFIX / "envs" / env).is_dir():
        pytest.skip(f"{env} not built under {PREFIX / 'envs'}")
    sim = simulate(SIM100K, SIM100K_SEED, tmp_path / "sim")
    faa = _plain(sorted(RECORDED.glob("proteins.faa*"))[0], tmp_path)
    recorded_gff = _plain(sorted(RECORDED.glob("*.gff*"))[0], tmp_path)
    mpi = tmp_path / "genome.mpi"
    _run(env, miniprot_index_cmd(sim.genome_fasta, mpi, 2, table=11))
    live_gff = tmp_path / "live.gff"
    _run(env, miniprot_align_cmd(mpi, faa, 2, ProteinsConfig().miniprot_args), stdout=live_gff)

    model = parse_bakta_json(sim.bakta_json).model
    lens = {r.id: len(r.sequence) for r in read_protein_fasta(faa)}
    cfg = ProteinsConfig()
    live = assign_alignments(parse_miniprot_gff(live_gff, model.replicons, lens), model, cfg)
    rec = assign_alignments(parse_miniprot_gff(recorded_gff, model.replicons, lens), model, cfg)
    pd.testing.assert_frame_equal(live[0].reset_index(drop=True), rec[0].reset_index(drop=True))
    assert live[1] == rec[1]
