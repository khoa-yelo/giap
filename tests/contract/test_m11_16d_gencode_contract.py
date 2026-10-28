"""Contract test (M11-16d; docs/specs/M11-16-translation-tables-4-25.md item 3.2, 3.3, Q-L4):
miniprot 0.18 ``-T 4`` / ``-T 25`` and tblastn (BLAST+ 2.16.0) ``-db_gencode 4`` on a 2-gene
fixture whose first protein has W (code 4) or G (code 25) at its TGA codons.

Runs the commands of the Python builders (``miniprot_index_cmd(..., table=t)``,
``tblastn_cmd(..., table=t)``) in the pinned envs, on Slurm only (``slow``; miniprot needs the
micromamba env, BLAST the pulled image). The fixture is built here, deterministically, from a
seeded generator (no network, no database):

- contig ``c1`` (linear, 4,000 nt); gene A on ``+`` at [500, 500 + 3 * 152): ATG, 150 sense
  codons of which 8 are TGA, then TAA; gene B on ``-`` at [2,500, 2,500 + 3 * 122): ATG, 120 sense
  codons without TGA, TAA;
- the query of gene A is its translation under the run's table (W or G at TGA).

Expected: under the matching table the whole protein A aligns on ``+`` from 500, without stop or
frameshift, at identity 1; under table 11 it does not (the TGA codons are stops). If miniprot
refuses ``-T 4`` or ``-T 25``, that is Q-L4's finding (the protein route at 4/25 needs a decision).
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pytest

from gippa.model.coords import Replicon

REPO = Path(__file__).resolve().parents[2]
_GIPPA_REPO = Path(os.environ.get("GIPPA_REPO", REPO))
_ROOT = Path(os.environ.get("GIPPA_ROOT", _GIPPA_REPO / "local"))
PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", _ROOT / "micromamba"))

pytestmark = [pytest.mark.slow]

BASES = "TCAG"
CODE11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
TGA_AS = {4: "W", 25: "G"}
A_START, B_START, LENGTH = 500, 2_500, 4_000
N_A, N_B, N_TGA = 150, 120, 8


def _translate(nt: str, table: int) -> str:
    code = dict(zip(CODONS, CODE11, strict=True))
    if table != 11:
        code["TGA"] = TGA_AS[table]
    return "".join(code[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fixture(tmp: Path, table: int) -> tuple[Path, Path, int]:
    rng = np.random.default_rng(1604)
    sense = [c for c, aa in zip(CODONS, CODE11, strict=True) if aa != "*"]
    genome = list("".join(rng.choice(list("ACGT"), LENGTH)))
    body_a = [str(c) for c in rng.choice(sense, N_A)]
    for k in sorted(rng.choice(np.arange(5, N_A - 5), N_TGA, replace=False)):
        body_a[int(k)] = "TGA"
    gene_a = "ATG" + "".join(body_a) + "TAA"
    gene_b = "ATG" + "".join(str(c) for c in rng.choice(sense, N_B)) + "TAA"
    genome[A_START : A_START + len(gene_a)] = gene_a
    genome[B_START : B_START + len(gene_b)] = _rc(gene_b)
    fa = tmp / "genome.fa"
    fa.write_text(">c1\n" + "".join(genome) + "\n", encoding="utf-8")
    faa = tmp / "proteins.faa"
    prot_a = _translate(gene_a[:-3], table)
    prot_b = _translate(gene_b[:-3], 11)
    assert "*" not in prot_a and prot_a.count(TGA_AS[table]) >= N_TGA
    faa.write_text(f">A\n{prot_a}\n>B\n{prot_b}\n", encoding="utf-8")
    return fa, faa, len(prot_a)


def _conda(
    env: str, argv: list[str], stdout: Path | None = None
) -> subprocess.CompletedProcess[str]:
    cmd = f"micromamba run -n {env} {shlex.join(argv)}"
    with open(stdout, "w") if stdout else open(os.devnull, "w") as out:
        return subprocess.run(
            cmd, shell=True, stdout=out, stderr=subprocess.PIPE, text=True, timeout=900
        )


def _miniprot_env() -> str:
    from gippa.annotate.tools import TOOLS

    env = TOOLS["miniprot"].env_name
    if not (PREFIX / "envs" / env).is_dir():
        pytest.skip(f"{env} not built under {PREFIX / 'envs'}")
    return env


def _miniprot(tmp: Path, table: int, index_table: int) -> list:  # type: ignore[type-arg]
    from gippa.config.model import ProteinsConfig
    from gippa.mapping.miniprot import miniprot_align_cmd, miniprot_index_cmd, parse_miniprot_gff

    env = _miniprot_env()
    fa, faa, n_a = _fixture(tmp, table)
    mpi = tmp / "genome.mpi"
    proc = _conda(env, miniprot_index_cmd(fa, mpi, 2, table=index_table))
    assert proc.returncode == 0, f"Q-L4: miniprot -T {index_table} refused: {proc.stderr[-2000:]}"
    gff = tmp / "miniprot.gff"
    proc = _conda(env, miniprot_align_cmd(mpi, faa, 2, ProteinsConfig().miniprot_args), gff)
    assert proc.returncode == 0, proc.stderr[-2000:]
    reps = {"c1": Replicon("c1", LENGTH, "linear")}
    lens = {"A": n_a, "B": N_B + 1}
    return [a for a in parse_miniprot_gff(gff, reps, lens) if a.protein_id == "A"]


def _full(a, n_a: int) -> bool:  # type: ignore[no-untyped-def]
    return (
        a.interval.start == A_START
        and a.interval.strand == "+"
        and a.q_start == 0
        and a.q_end == n_a
        and a.identity == 1.0
        and a.stop_codon == 0
        and a.frameshift == 0
    )


@pytest.mark.needs_conda
@pytest.mark.parametrize("table", [4, 25])
def test_miniprot_index_table_aligns_the_tga_protein(tmp_path: Path, table: int) -> None:
    hits = _miniprot(tmp_path, table, table)
    n_a = N_A + 1
    assert any(_full(a, n_a) for a in hits), hits


@pytest.mark.needs_conda
def test_miniprot_index_table_11_does_not(tmp_path: Path) -> None:
    hits = _miniprot(tmp_path, 4, 11)
    assert not any(_full(a, N_A + 1) for a in hits), hits


def _blast_image() -> Path:
    conf = (REPO / "workflow" / "conf" / "containers.config").read_text()
    m = re.search(
        r"container = '(?:docker://)?quay\.io/biocontainers/blast@sha256:([0-9a-f]{64})'", conf
    )
    assert m, "BLAST image digest not pinned in workflow/conf/containers.config"
    img = _ROOT / "containers" / f"quay.io-biocontainers-blast@sha256-{m.group(1)}.img"
    if not img.is_file():
        pytest.skip(f"BLAST image not pulled: {img}")
    return img


def _tblastn(tmp: Path, table: int, search_table: int) -> list:  # type: ignore[type-arg]
    from gippa.config.model import ProteinsConfig
    from gippa.mapping.tblastn import makeblastdb_cmd, parse_tblastn, tblastn_cmd

    img = _blast_image()
    fa, faa, _ = _fixture(tmp, table)
    db = tmp / "db" / "genome"
    out = tmp / "tblastn.tsv"
    for argv, stdout in (
        (makeblastdb_cmd(fa, db), None),
        (tblastn_cmd(db, faa, 2, ProteinsConfig().tblastn_args, table=search_table), out),
    ):
        cmd = ["apptainer", "exec", "-B", str(tmp), str(img), *argv]
        with open(stdout, "w") if stdout else open(os.devnull, "w") as fh:
            proc = subprocess.run(
                cmd, cwd=tmp, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=900
            )
        assert proc.returncode == 0, proc.stderr[-2000:]
    return [h for h in parse_tblastn(out) if h.qseqid == "A"]


@pytest.mark.needs_apptainer
def test_tblastn_db_gencode_4_aligns_the_w_protein(tmp_path: Path) -> None:
    hits = _tblastn(tmp_path, 4, 4)
    n_a = N_A + 1
    assert any(
        h.qstart == 1 and h.qend == n_a and h.nident == h.length == n_a and h.sstart == A_START + 1
        for h in hits
    ), hits


@pytest.mark.needs_apptainer
def test_tblastn_db_gencode_11_splits_the_w_protein(tmp_path: Path) -> None:
    hits = _tblastn(tmp_path, 4, 11)
    n_a = N_A + 1
    assert not any(h.qstart == 1 and h.qend == n_a and h.nident == n_a for h in hits), hits
