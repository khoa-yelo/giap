"""`giap setup` adopts the databases the prefetch job left in $GIAP_DB_DIR (issue M2-08; §3.9).

Read-only: only ``verify`` runs (no download, no manifest written, $GIAP_ROOT untouched). The
layouts are the real ones (``bakta/db-light``, ``bakta/db``, ``genomad_db``,
``promotech/RF-HOT.model``), so ``SOURCES[...].path`` must point at them.
The live Promotech download test of the acceptance list is left out: tests never download (lead).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("GIAP_ROOT", REPO / "local"))
DB = Path(os.environ.get("GIAP_DB_DIR", ROOT / "db"))


def _verify(name: str) -> Any:
    from giap.setupdb import SOURCES, verify

    return verify(SOURCES[name], DB)


@pytest.mark.needs_bakta_db
def test_prefetched_bakta_light_verifies() -> None:
    if not (DB / "bakta" / "db-light").is_dir():
        pytest.skip("no bakta light DB")
    rec = _verify("bakta_light")
    assert (rec.status, rec.version, rec.path) == ("ok", "6.0", "bakta/db-light")


@pytest.mark.needs_bakta_full
def test_prefetched_bakta_full_verifies() -> None:
    rec = _verify("bakta_full")
    assert (rec.status, rec.version, rec.path) == ("ok", "6.0", "bakta/db")


@pytest.mark.needs_genomad_db
def test_prefetched_genomad_verifies() -> None:
    from giap.setupdb import SOURCES

    assert SOURCES["genomad"].path == "genomad_db"
    rec = _verify("genomad")
    assert (rec.status, rec.version, rec.path) == ("ok", "1.9", "genomad_db"), rec


@pytest.mark.needs_bakta_full
def test_real_amrfinderplus_version_file_is_where_setup_reads_it() -> None:
    v = DB / "bakta" / "db" / "amrfinderplus-db" / "latest" / "version.txt"
    assert v.is_file()
    assert v.read_text().strip() == "2026-08-07.1"


@pytest.mark.slow
@pytest.mark.needs_promotech_models
def test_prefetched_promotech_model_verifies() -> None:
    """Hashes the kept 1 GB RF-HOT.zip against the ADR-0014 SHA-256 pin."""
    rec = _verify("promotech_models")
    assert rec.status == "ok", rec
