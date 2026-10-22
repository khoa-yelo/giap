"""M7-04 integration (#160 follow-up): a two-residue block reading Asn-Ala is ``alt_aa`` "NA",
the same text as the ADR-0006 missing token. The residue columns are required and never
missing, so the schema gives them no missing values and ``read_variants`` keeps the text."""

from __future__ import annotations

from pathlib import Path

from helpers.variant_tables_fixture import GP, peptides_at, records, sub, vt
from unit.io.test_variant_tables import _assert_valid


def test_asn_ala_residues_survive_write_validate_read(tmp_path: Path) -> None:
    _, table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    table = table.copy()
    table.loc[table.index[0], ["ref_aa", "alt_aa"]] = ["NA", "NA"]
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, table)
    _assert_valid(path, "variants")
    back = vt().read_variants(path)
    assert back.loc[0, "ref_aa"] == "NA" and back.loc[0, "alt_aa"] == "NA"
    out = tmp_path / "again.tsv"
    vt().write_variants(out, back)
    assert out.read_bytes() == path.read_bytes()
