"""Live makeblastdb + tblastn (BLAST+ 2.16.0) equal the recorded fixtures (issue M3-13).

Runs the commands built by ``giap.mapping.tblastn`` inside the pinned BLAST image
(``workflow/conf/containers.config``, pulled to ``$GIAP_ROOT/containers/``; ``needs_apptainer``):
- on the tiny fixture ``tests/data/tblastn/mini/`` (seconds): the raw tabular output equals the
  recorded ``tblastn.tsv`` line for line;
- on the 100-kb sim genome (``pending_sample`` until ``tests/data/recorded/blast/2.16.0/sim100k/``
  exists): the merged, assigned tables equal those of the recording.
"""

from __future__ import annotations

import gzip
import os
import re
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from giap.config.model import ProteinsConfig

REPO = Path(__file__).resolve().parents[2]
MINI = REPO / "tests" / "data" / "tblastn" / "mini"
RECORDED = REPO / "tests" / "data" / "recorded" / "blast" / "2.16.0" / "sim100k"
_GIAP_REPO = Path(os.environ.get("GIAP_REPO", REPO))
_ROOT = Path(os.environ.get("GIAP_ROOT", _GIAP_REPO / "local"))

pytestmark = [pytest.mark.needs_apptainer, pytest.mark.slow]


def _image() -> Path:
    conf = (REPO / "workflow" / "conf" / "containers.config").read_text()
    m = re.search(
        r"container = '(?:docker://)?quay\.io/biocontainers/blast@sha256:([0-9a-f]{64})'", conf
    )
    assert m, "BLAST image digest not pinned in workflow/conf/containers.config"
    img = _ROOT / "containers" / f"quay.io-biocontainers-blast@sha256-{m.group(1)}.img"
    if not img.is_file():
        pytest.skip(f"BLAST image not pulled: {img}")
    return img


def _run(img: Path, argv: list[str], cwd: Path, stdout: Path | None = None) -> None:
    cmd = ["apptainer", "exec", "-B", str(cwd), str(img), *argv]
    with open(stdout, "w") if stdout else open(os.devnull, "w") as out:
        proc = subprocess.run(
            cmd, cwd=cwd, stdout=out, stderr=subprocess.PIPE, text=True, timeout=1800
        )
    assert proc.returncode == 0, proc.stderr


def _search(img: Path, genome: Path, faa: Path, tmp: Path) -> Path:
    from giap.mapping.tblastn import makeblastdb_cmd, tblastn_cmd

    db = tmp / "db" / "genome"
    _run(img, makeblastdb_cmd(genome, db), tmp)
    out = tmp / "live.tsv"
    _run(img, tblastn_cmd(db, faa, 2, ProteinsConfig().tblastn_args), tmp, stdout=out)
    return out


def test_live_mini_equals_recorded(tmp_path: Path) -> None:
    img = _image()
    genome = tmp_path / "genome.fa"
    faa = tmp_path / "proteins.faa"
    genome.write_bytes((MINI / "genome.fa").read_bytes())
    faa.write_bytes((MINI / "proteins.faa").read_bytes())
    live = _search(img, genome, faa, tmp_path)
    assert live.read_text() == (MINI / "tblastn.tsv").read_text()


def _plain(path: Path, tmp: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = tmp / path.name.removesuffix(".gz")
    out.write_bytes(gzip.decompress(path.read_bytes()))
    return out


@pytest.mark.pending_sample
@pytest.mark.skipif(not RECORDED.is_dir(), reason=f"pending_sample: no recording in {RECORDED}")
def test_live_sim100k_equals_recorded(tmp_path: Path) -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.mapping.protein_assign import assign_alignments
    from giap.mapping.tblastn import merge_hsps, parse_tblastn
    from giapsim import SIM100K, SIM100K_SEED, simulate

    img = _image()
    sim = simulate(SIM100K, SIM100K_SEED, tmp_path / "sim")
    work = tmp_path / "work"
    work.mkdir()
    genome = work / "genome.fa"
    genome.write_bytes(sim.genome_fasta.read_bytes())
    faa = work / "proteins.faa"
    faa.write_bytes(_plain(sorted(RECORDED.glob("proteins.faa*"))[0], tmp_path).read_bytes())
    live_tsv = _search(img, genome, faa, work)
    rec_tsv = _plain(sorted(RECORDED.glob("*.tsv*"))[0], tmp_path)

    model = parse_bakta_json(sim.bakta_json).model
    cfg = ProteinsConfig()
    live = assign_alignments(merge_hsps(parse_tblastn(live_tsv), model.replicons), model, cfg)
    rec = assign_alignments(merge_hsps(parse_tblastn(rec_tsv), model.replicons), model, cfg)
    pd.testing.assert_frame_equal(live[0].reset_index(drop=True), rec[0].reset_index(drop=True))
    assert live[1] == rec[1]
