"""TE1 and TE2 on Nissle with recoded Sec codons (M7-16c; spec "TE2 and TE1 expectations on
Nissle (16c)", rules 24-25, Q26; test plan C14; ``slow`` and ``e2e``: Slurm only, never in the
session allocation).

TE1: ``misses`` stays 4 and the share 0.99998; the 4 fdnG exceptions become ``kind = sec_site``
with ``hits`` ``tolerant:cds_inframe:<copy>`` and ``counts["sec_site"] == 4``; each records
``placement_q``, ``m_target`` (1) and ``accepted`` in the result JSON.

TE2: ``alt_sites`` 235, ``high`` 80, ``high_start_codon`` 80, ``high_stop_codon`` 0,
``violations`` 0, ``alt_stop_codon`` 2, and the note is rule 25's text.

"Rejection is a finding": the counts hold only if all 4 Sec peptides are accepted. If one is
rejected it stays ``kind = miss`` / ``unplaced`` with ``accepted = False``, TE2 has one site fewer,
and this test records the values (a warning naming the peptides, the JSON in the test's tmp
directory) instead of asserting those counts; alpha and every criterion stay as they are (R3).
"""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
import warnings
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
        mp.setenv("GIAP_DB_DIR", str(_db))
        r = _te().te1_proxy(NISSLE_GENOME, NISSLE_JSON)
    out = tmp_path_factory.mktemp("te1") / "te1.json"
    out.write_text(json.dumps(r.to_json(), indent=2, sort_keys=True))
    print(f"TE1 result JSON: {out}")
    return r


@pytest.fixture(scope="module")
def te2(_db: Path) -> Any:
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(_db))
        return _te().te2_proxy(NISSLE_GENOME, NISSLE_JSON)


def _sec_rows(r: Any) -> list[dict[str, Any]]:
    return [x for x in r.exceptions if x["kind"] in ("miss", "sec_site")]


def _all_accepted(r: Any) -> bool:
    return all(bool(x.get("accepted")) for x in _sec_rows(r))


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


def test_te1_sec_site_labels(te1: Any) -> None:
    rows = _sec_rows(te1)
    if not _all_accepted(te1):
        rejected = [x["peptide"] for x in rows if not bool(x.get("accepted"))]
        warnings.warn(
            f"finding: Sec peptides rejected by the placement FDR: {rejected}", stacklevel=1
        )
        for x in rows:
            if not bool(x.get("accepted")):
                assert (x["kind"], x["hits"]) == ("miss", "unplaced"), x
        assert te1.counts["sec_site"] == sum(bool(x.get("accepted")) for x in rows)
        return
    assert te1.counts["sec_site"] == 4
    for x in rows:
        assert x["kind"] == "sec_site", x
        parts = x["hits"].split(";")
        assert parts and all(p.startswith("tolerant:cds_inframe:JMJJOG_") for p in parts), x


def test_te2_counts_and_note(te2: Any, te1: Any) -> None:
    c = te2.counts
    assert (c["high"], c["high_start_codon"], c["high_stop_codon"], c["violations"]) == (
        80,
        80,
        0,
        0,
    )
    assert te2.passed is True
    assert te2.notes == [_te().TE2_NO_STOP_NOTE]
    if _all_accepted(te1):
        assert (c["alt_sites"], c["alt_stop_codon"]) == (235, 2)
    else:
        warnings.warn(
            f"finding: TE2 alt_sites {c['alt_sites']}, alt_stop_codon "
            f"{c['alt_stop_codon']} with a rejected Sec peptide",
            stacklevel=1,
        )
        assert c["alt_stop_codon"] <= 2 and 233 <= c["alt_sites"] <= 235
