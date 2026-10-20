"""TE8 on the full UniProt K-12 proteome (issue M3-22 example 6; GIAP_SPEC.md §11.6 TE8, §12).

Weekly only (lead decision 3): ``e2e`` + ``slow`` + ``needs_net("rest.uniprot.org")``; without the
network the test is skipped by the capability marker (example 6). It never runs in CI or in the
session allocation (``scripts/slurm/weekly.sbatch``).

Assumed: ``python validation/fetch_uniprot_k12.py`` downloads UP000000625 into
``$GIAP_DATA_DIR/uniprot/`` (a FASTA named ``UP000000625*.fa*``) and records its SHA-256 and the
retrieval date in a file next to it (§11.4, R11: never committed). The test points
``$GIAP_DATA_DIR`` at a temporary directory, so it never touches the shared data directory, and
uses a test ``$GIAP_DB_DIR``.
"""

from __future__ import annotations

import hashlib
import importlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.mapping_steps import write_test_db

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "tests" / "data"
BAKTA = DATA / "recorded" / "bakta" / "1.12.2-full"
DATE = re.compile(r"\b20\d\d-\d\d-\d\d")

pytestmark = [pytest.mark.e2e, pytest.mark.slow, pytest.mark.needs_net("rest.uniprot.org")]


def test_te8_full_uniprot_k12_passes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    data_dir = tmp_path / "data"
    monkeypatch.setenv("GIAP_DATA_DIR", str(data_dir))
    monkeypatch.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path / "db")))
    script = REPO / "validation" / "fetch_uniprot_k12.py"
    run = subprocess.run(
        [sys.executable, str(script)], env=dict(os.environ), capture_output=True, text=True
    )
    assert run.returncode == 0, run.stderr
    fastas = sorted((data_dir / "uniprot").glob("UP000000625*.fa*"))
    assert len(fastas) == 1, fastas
    fasta = fastas[0]
    sha = hashlib.sha256(fasta.read_bytes()).hexdigest()
    sidecars = [p for p in (data_dir / "uniprot").iterdir() if p != fasta and p.is_file()]
    assert any(sha in p.read_text() and DATE.search(p.read_text()) for p in sidecars), sidecars

    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    te = importlib.import_module("validation.proxies.te_proxies")
    r = te.te8_proxy(
        DATA / "genomes" / "k12_mg1655.fna.gz",
        BAKTA / "k12" / "k12_mg1655.json.gz",
        DATA / "genomes" / "nissle1917.fna.gz",
        BAKTA / "nissle" / "nissle1917.json.gz",
        fasta,
    )
    # release 2026_03 (4,403 proteins): 151,424 distinct peptides, 354 non_atg_nterm (decision 6),
    # share 150,918 / 151,070 = 0.99899; not pinned because the live release changes
    c = r.counts
    assert c["peptides"] > 100_000
    assert 0 < c["non_atg_nterm"] < 0.01 * c["peptides"]
    assert r.shares["k12_exact"] == pytest.approx(
        c["k12_exact"] / (c["peptides"] - c["non_atg_nterm"])
    )
    assert r.shares["k12_exact"] >= 0.995
    assert r.shares["nissle_placement_fdr"] <= 0.01
    for key in ("nissle_tolerant", "nissle_unmapped"):
        assert 0.0 <= r.shares[key] <= 1.0
    assert r.passed is True
