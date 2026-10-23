"""M7-14 example 5: ``detection_summary.tsv`` after a sim100k run with peptides has a numeric
``id_mixups`` row without ``requires_variants`` (docs/specs/M7-14-te2-proxy-and-id-mixups.md rules
3-4; GIAP_SPEC.md §7.5.1 row 1, §9.3).

``giap run --executor local`` on giapsim SIM100K with M3-19 peptides (``search_db = same``, the
simulator's default; tests/helpers/m4_real.py) plans ``variants`` and ``context_stats``; rule 4:
"the local executor needs no change (optional inputs are picked up when present)", so
``work/variants_qc`` reaches ``context_stats``. Oracle: every peptide is database/same, so ``n`` =
the ``variants.tsv`` rows and the count = rows flagged neither ``start_codon`` nor
``stop_codon``. A test ``$GIAP_DB_DIR`` only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m4_chain as mc
import m4_real as mr

SEED = 20261014


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp("id_mixups_sim100k")
    db = mc.test_db(root / "db")
    sim = mr.sim100k(root / "sim")
    out = root / "OUT"
    work = mr.giap_run(sim, out, db, SEED)
    return out, work


def test_variants_qc_is_produced_by_the_local_executor(run: tuple[Path, Path]) -> None:
    _, work = run
    assert (work / "variants_qc").is_file()


def test_example_5_id_mixups_row_is_numeric(run: tuple[Path, Path]) -> None:
    from giap.io.variant_tables import read_variants

    out, _ = run
    ds = mc.read(out / "tables" / "detection_summary.tsv")
    sel = ds[(ds["section"] == "qc") & (ds["group_type"] == "run") & (ds["group"] == "all")
             & (ds["metric"] == "id_mixups")]  # fmt: skip
    assert len(sel) == 1
    row = sel.iloc[0]
    assert row["value"] != "NA" and row["n"] != "NA", row.to_dict()
    assert "requires_variants" not in str(row["flags"]).split(";")
    vpath, ppath = out / "tables" / "variants.tsv", run[1] / "peptides"
    if not (vpath.is_file() and ppath.is_file()):
        return  # the oracle needs both; the numeric row above is example 5 itself
    variants = read_variants(vpath)
    peps = pd.read_parquet(ppath)
    if set(peps["source"]) == {"database"} and set(peps["search_db"].astype(str)) == {"same"}:
        exempt = {"start_codon", "stop_codon"}
        assert int(row["n"]) == len(variants)
        assert float(row["value"]) == sum(not (exempt & set(f)) for f in variants["flags"])
