"""TE2 in-silico proxy (docs/specs/M7-14-te2-proxy-and-id-mixups.md rules 1-2 and 5, examples 1-2;
GIAP_SPEC.md §11.6 TE2 "no ``high`` ALT calls from A2 except at sites flagged ``start_codon`` or
``stop_codon`` | pilot; in-silico proxy on same-strain peptides", §13 gate M7 "TE2 proxy passes").

Interface (see tests/unit/test_te2_proxy_result.py for the pinned ``ProxyResult`` fields)::

    validation.proxies.te_proxies.te2_proxy(genome: Path, bakta_json: Path) -> ProxyResult

Peptides = TE1's digest of every Bakta CDS/sORF protein, as ``database`` with ``search_db =
same``, through ``run_chain`` and M7-04's ``run_variants``. Marked as TE1's test is (``slow``);
the Nissle run (≈ 164,000 peptides on 5.4 Mb) is also ``e2e`` and runs only in the weekly Slurm
job (``scripts/slurm/weekly.sbatch``, R14), never in the session allocation. Every run gets a
test ``$GIAP_DB_DIR`` (never the real one). The committed Nissle genome is identified by its
SHA-256 only (tests/data/genomes/PROVENANCE.md).
"""

from __future__ import annotations

import hashlib
import importlib
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.mapping_steps import write_test_db

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full" / "nissle"
NISSLE_JSON = RECORDED / "nissle1917.json.gz"
NISSLE_GENOME = REPO / "tests" / "data" / "genomes" / "nissle1917.fna.gz"
NISSLE_SHA256 = "fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee"
COUNT_KEYS = {"alt_sites", "high", "high_start_codon", "high_stop_codon", "violations"}

pytestmark = pytest.mark.slow


def _te_proxies() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.proxies.te_proxies")


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db"))))


def _check_shape(r: Any) -> None:
    c = r.counts
    assert r.name == "TE2"
    assert set(c) == COUNT_KEYS
    assert c["violations"] <= c["high"] <= c["alt_sites"]
    assert r.passed is (c["violations"] == 0)
    assert len(r.exceptions) == min(50, c["violations"])
    for ex in r.exceptions:
        assert {"locus_tag", "aa_pos", "ref", "alt"} <= set(ex)


def _sim(tmp: Path, seed: int, **kw: Any) -> Any:
    from giapsim.simulate import SimConfig, simulate

    return simulate(SimConfig(genome_length=100_000, write_precomputed=True, **kw), seed, tmp)


# ------------------------------------------------------------------ the proxy on giapsim


def test_te2_proxy_on_an_all_atg_sim_has_no_alt_site(tmp_path: Path) -> None:
    """Same-strain peptides with ATG starts only: nothing contradicts the genome."""
    res = _sim(tmp_path / "sim", 2220, start_codons=(("ATG", 1.0),))
    r = _te_proxies().te2_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    assert r.counts["alt_sites"] == 0
    assert r.passed is True


def test_te2_proxy_on_gtg_starts_flags_every_high_site_start_codon(tmp_path: Path) -> None:
    """Rule 2: Bakta's M at GTG/TTG starts appears as ALT ``M`` flagged ``start_codon``; with
    50 % GTG starts some reach ``high`` and none is a violation."""
    res = _sim(tmp_path / "sim", 2221, start_codons=(("ATG", 0.5), ("GTG", 0.5)))
    r = _te_proxies().te2_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    assert r.counts["alt_sites"] > 0, "the GTG starts are meant to give ALT M sites"
    assert r.counts["high_start_codon"] > 0
    assert r.counts["high_stop_codon"] == 0
    assert r.counts["violations"] == 0
    assert r.passed is True


def test_example_1_same_strain_insilico_scenario(tmp_path: Path) -> None:
    """Example 1: M7-12's ``same_strain_insilico`` through ``run_variants``: zero ``high`` rows
    without ``start_codon``."""
    scenarios = pytest.importorskip("giapsim.scenarios")
    if not hasattr(scenarios, "same_strain_insilico"):
        pytest.skip("giapsim.scenarios.same_strain_insilico (M7-12) is not merged")
    from giapsim.simulate import simulate

    res = simulate(scenarios.same_strain_insilico(2222), 2222, tmp_path / "sim")
    r = _te_proxies().te2_proxy(res.genome_fasta, res.bakta_json)
    _check_shape(r)
    assert r.counts["violations"] == 0 and r.passed is True


# ------------------------------------------------------------------ example 2: Nissle (weekly)


@pytest.mark.e2e
@pytest.mark.pending_sample
@pytest.mark.skipif(
    not (NISSLE_JSON.is_file() and NISSLE_GENOME.is_file()),
    reason=f"pending_sample: no recorded Nissle full-DB Bakta JSON at {NISSLE_JSON}",
)
def test_example_2_te2_nissle_passes() -> None:
    assert hashlib.sha256(NISSLE_GENOME.read_bytes()).hexdigest() == NISSLE_SHA256
    r = _te_proxies().te2_proxy(NISSLE_GENOME, NISSLE_JSON)
    _check_shape(r)
    assert r.counts["violations"] == 0, r.exceptions
    assert r.passed is True
    # rule 2: the 4 TE1 misses (Sec, fdnG) stay unmapped until M7-16 -> no stop_codon site yet
    assert r.counts["high_stop_codon"] == 0
    # high_start_codon is recorded (expected > 0), not asserted (example 2)
    assert r.counts["high_start_codon"] >= 0
