"""M7-03 review (#160): ``datapackage.json`` is a valid tabular Data Package. Laid out with the
written tables and their schemas, ``frictionless`` validates the package and both resources
(§9.5); parsing the JSON alone missed the missing resource profiles."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from helpers.variant_tables_fixture import GP, peptides_at, records, run, sub, vt
from unit.io.test_variant_tables import SCHEMAS, TABLES


def test_datapackage_validates_with_written_tables(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    pile, _ = run(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    _, table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    (tmp_path / "tables").mkdir()
    vt().write_pileup(tmp_path / "tables" / "pileup.tsv.gz", pile)
    vt().write_variants(tmp_path / "tables" / "variants.tsv", table)
    for name in ("pileup", "variants"):
        shutil.copy(TABLES / f"{name}.schema.json", tmp_path / "tables")
    # only the variant resources: other modules (M7-07 novel_orfs) share the descriptor
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    pkg["resources"] = [r for r in pkg["resources"] if r["name"] in ("pileup", "variants")]
    (tmp_path / "datapackage.json").write_text(json.dumps(pkg), encoding="utf-8")
    report = frictionless.validate(str(tmp_path / "datapackage.json"))
    assert report.valid, report.flatten(["type", "note"])
    assert len(report.tasks) == 2
