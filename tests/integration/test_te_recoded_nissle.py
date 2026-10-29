"""TE1 and TE2 on Nissle with recoded Sec codons (M7-16c; spec "TE2 and TE1 expectations on
Nissle (16c)", rules 24-25, Q26; test plan C14; ``slow`` and ``e2e``: Slurm only, never in the
session allocation).

TE1: ``misses`` stays 4 and the share 0.99998; each of the 4 fdnG Sec peptides records
``placement_q``, ``m_target`` (1) and ``accepted`` in the result JSON.

M11-20 Nissle finding (docs/specs/M11-20-initiator-methionine.md, end of § What changes;
statistician 2026-10-29, accepted): the non-ATG N-terminal peptides are credited exact rows
(stratum 0), so the 4 fdnG Sec peptides are alone in placement stratum 1 (T1 = 4, D1 = 0) and
knockoff+ gives q = 1/4: all 4 are ``rejected_placement_fdr`` (``kind = miss``, ``hits =
unplaced``, ``accepted = False``, ``placement_q = 0.25``), TE1 ``sec_site`` 0, TE2
``alt_stop_codon`` 0 and ``alt_sites`` 0 (no ALT ``M`` at credited starts, no ``U`` site), and the
note is rule 25's text naming that rejection. These are pinned: a change to stratum 1 (a Sec
peptide accepted, or a different q) fails here instead of warning (R3: alpha and every criterion
stay as they are).
"""

from __future__ import annotations

import hashlib
import importlib
import json
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
FDNG = {"JMJJOG_02511", "JMJJOG_04733"}

pytestmark = [
    pytest.mark.slow,
    pytest.mark.e2e,
    pytest.mark.skipif(
        not (NISSLE_JSON.is_file() and NISSLE_GENOME.is_file()),
        reason=f"no recorded Nissle full-DB Bakta JSON at {NISSLE_JSON}",
    ),
]


def _te() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.proxies.te_proxies")


@pytest.fixture(scope="module")
def _db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return write_test_db(tmp_path_factory.mktemp("test_db"))


@pytest.fixture(scope="module")
def te1(_db: Path, tmp_path_factory: pytest.TempPathFactory) -> Any:
    assert hashlib.sha256(NISSLE_GENOME.read_bytes()).hexdigest() == NISSLE_SHA256
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(_db))
        r = _te().te1_proxy(NISSLE_GENOME, NISSLE_JSON)
    out = tmp_path_factory.mktemp("te1") / "te1.json"
    out.write_text(json.dumps(r.to_json(), indent=2, sort_keys=True))
    print(f"TE1 result JSON: {out}")
    return r


@pytest.fixture(scope="module")
def te2(_db: Path) -> Any:
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(_db))
        return _te().te2_proxy(NISSLE_GENOME, NISSLE_JSON)


def _sec_rows(r: Any) -> list[dict[str, Any]]:
    return [x for x in r.exceptions if x["kind"] in ("miss", "sec_site")]


def test_te1_misses_and_share_unchanged(te1: Any) -> None:
    assert te1.counts["misses"] == 4
    assert te1.shares["exact_inframe"] == pytest.approx(0.99998, abs=5e-6)
    assert te1.passed is True


def test_te1_records_the_q_values_of_the_four_sec_peptides(te1: Any) -> None:
    rows = _sec_rows(te1)
    assert len(rows) == 4
    doc = te1.to_json()
    for x in rows:
        assert x["locus_tag"] in FDNG, x
        assert {"placement_q", "m_target", "accepted"} <= set(x), x
        assert int(x["m_target"]) == 1, x
    assert all("placement_q" in x for x in doc["exceptions"] if x["kind"] != "non_atg_nterm")
    print("Sec peptides:", [(x["peptide"], x["placement_q"], x["accepted"]) for x in rows])


def test_te1_sec_peptides_rejected_by_the_stratum_1_placement_fdr(te1: Any) -> None:
    """The Nissle finding: T1 = 4, D1 = 0, q = (0 + 1) / 4 for each; none accepted."""
    rows = _sec_rows(te1)
    assert len(rows) == 4
    accepted = [x["peptide"] for x in rows if bool(x.get("accepted"))]
    assert not accepted, (
        f"placement stratum 1 changed on Nissle: Sec peptides {accepted} are now accepted; the "
        "M11-20 Nissle finding pins all 4 rejected at q = 0.25 (update the spec first)"
    )
    for x in rows:
        assert (x["kind"], x["hits"]) == ("miss", "unplaced"), x
        assert float(x["placement_q"]) == pytest.approx(0.25), x
    assert te1.counts["sec_site"] == 0


def test_te2_counts_and_note(te2: Any) -> None:
    c = te2.counts
    assert (c["high"], c["high_start_codon"], c["high_stop_codon"], c["violations"]) == (
        0,  # M11-20: no ALT M at credited starts (was 80, all start_codon)
        0,
        0,
        0,
    )
    # the Nissle finding: the Sec peptides are rejected, so no ALT U site and no ALT M site
    assert (c["alt_sites"], c["alt_stop_codon"]) == (0, 0), (
        f"TE2 alt_sites {c['alt_sites']}, alt_stop_codon {c['alt_stop_codon']}: placement "
        "stratum 1 changed on Nissle (the M11-20 Nissle finding pins both at 0)"
    )
    assert te2.passed is True
    assert te2.notes == [_te().TE2_NO_STOP_NOTE]
    note = te2.notes[0]
    assert "rejected by the stratum-1 placement FDR (4 targets, no decoys, q = 0.25)" in note
