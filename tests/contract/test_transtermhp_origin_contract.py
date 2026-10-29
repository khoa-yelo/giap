"""Contract test: real TransTermHP 2.09 on ``write_transterm_coords`` output for a tiny circular
genome with genes across the origin (Slurm job 1903204, M. hyorhinis code 4: TransTermHP exited 3
with ``Bad gene coordinates`` on an origin-crossing gene; GIPPA_SPEC.md §7.1.1, §7.1.5).

Runs the pinned image (``workflow/conf/containers.config``, pulled under
``$GIPPA_ROOT/containers``) with the command of ``transterm_command`` and the ``expterm.dat``
lookup of ``workflow/modules/local/transtermhp/main.nf``; Slurm only (``slow``), skipped when the
image is not pulled. The fixture is a seeded random genome: contig ``c1`` (6,000 bp, circular)
with a ``+`` gene across the origin, contig ``c2`` (5,000 bp, circular) with a ``-`` gene across
it, plus ordinary genes on both strands.

Expected: exit 0; every parsed terminator lies in range on its contig.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

from gippa.annotate.terminators import parse_transterm_tt, write_transterm_coords
from gippa.annotate.transtermhp_cmd import transterm_command
from gippa.model.coords import Interval, Replicon, check_interval
from gippa.model.genes import Gene

REPO = Path(__file__).resolve().parents[2]
_GIPPA_REPO = Path(os.environ.get("GIPPA_REPO", REPO))
_ROOT = Path(os.environ.get("GIPPA_ROOT", _GIPPA_REPO / "local"))

pytestmark = [pytest.mark.slow]

REPS = {
    "c1": Replicon("c1", 6_000, "circular", "chromosome"),
    "c2": Replicon("c2", 5_000, "circular", "plasmid"),
}
HAIRPIN = "AAAAAAAAAAGCCCGCCTAATGAGCGGGCTTTTTTTTT"  # a canonical intrinsic terminator


def _image() -> Path:
    conf = (REPO / "workflow" / "conf" / "containers.config").read_text()
    m = re.search(
        r"container = '(?:docker://)?quay\.io/biocontainers/transtermhp@sha256:([0-9a-f]{64})'",
        conf,
    )
    assert m, "TransTermHP image digest not pinned in workflow/conf/containers.config"
    img = _ROOT / "containers" / f"quay.io-biocontainers-transtermhp@sha256-{m.group(1)}.img"
    if not img.is_file():
        pytest.skip(f"TransTermHP image not pulled: {img}")
    return img


def _gene(tag: str, contig: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * ((end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


GENES = [
    _gene("C1_0001", "c1", 200, 1_100, "+"),
    _gene("C1_0002", "c1", 2_000, 3_200, "-"),
    _gene("C1_0003", "c1", 5_400, 6_600, "+"),  # across the origin: 1-based 5401..6600
    _gene("C2_0001", "c2", 300, 1_200, "+"),
    _gene("C2_0002", "c2", 4_700, 5_600, "-"),  # across the origin: 1-based 4701..5600
]


def _genome(tmp: Path) -> Path:
    rng = np.random.default_rng(1903204)
    fa = tmp / "genome.fa"
    recs = []
    for r in REPS.values():
        seq = list("".join(rng.choice(list("ACGT"), r.length)))
        for pos in (1_120, 1_900):  # terminators after C?_0001 and before the minus gene
            seq[pos : pos + len(HAIRPIN)] = HAIRPIN
        recs.append(f">{r.id}\n{''.join(seq)}\n")
    fa.write_text("".join(recs))
    return fa


@pytest.mark.needs_apptainer
def test_transtermhp_runs_on_genes_across_the_origin(tmp_path: Path) -> None:
    img = _image()
    _genome(tmp_path)
    coords = tmp_path / "genes.coords"
    write_transterm_coords(coords, GENES, replicons=REPS)
    for ln in coords.read_text().splitlines():
        _tag, a, b, contig = ln.split(" ")
        assert 1 <= int(a) <= REPS[contig].length and 1 <= int(b) <= REPS[contig].length, ln
    argv = transterm_command(Path("genome.fa"), Path("genes.coords"), Path("$EXPTERM"))
    script = (
        'EXPTERM=$(dirname "$(command -v transterm)")/../data/expterm.dat; '
        + " ".join(argv)
        + " > transtermhp.tt"
    )
    proc = subprocess.run(
        ["apptainer", "exec", "-B", str(tmp_path), "--pwd", str(tmp_path), str(img),
         "sh", "-c", script],
        cwd=tmp_path, capture_output=True, text=True, timeout=600,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr[-2000:]
    tt = tmp_path / "transtermhp.tt"
    text = tt.read_text()
    assert "SEQUENCE c1" in text and "SEQUENCE c2" in text
    got = parse_transterm_tt(tt, REPS)
    for t in got:
        check_interval(t.interval, REPS[t.interval.contig])
        assert t.interval.end <= REPS[t.interval.contig].length
