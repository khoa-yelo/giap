"""PR #228 review B1 (lead decision 15): the peptide card states the placement claim per family
(q ≤ alpha, 0 / ≥ 1 difference, per decoy type) and lists the ``mapping_fdr`` stratum rows only as
descriptive values, each with its full key (decoy type, differences, length class)."""

from __future__ import annotations

from pathlib import Path

from helpers.report_content import page, words, write_mapping
from helpers.report_fixtures import run_out
from helpers.viz_out_tree import write_schema_table


def _rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for d, base in (("reverse", 0.001), ("shuffle", 0.002)):
        for n in (0, 1, 2):
            for k, lc in enumerate(("7-8", "9-12")):
                rows.append(
                    {
                        "section": "placement",
                        "length_class": lc,
                        "n_diff": n,
                        "targets": 100,
                        "decoys": 1,
                        "n_accepted": 90,
                        "decoy_type": d,
                        "fdr_stratum": round(base + 0.0001 * (10 * n + k), 6),
                    }
                )
    return rows  # fmt: skip


def _card(tmp: Path) -> str:
    out = run_out(tmp)
    write_mapping(out, placed_q=[0.001], unmapped_peptides=0, proteins=["P1"],
                  unmapped_proteins=0, strata={"7-8": 0.0072})  # fmt: skip
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", _rows())
    return words(page(out).cards()["peptides_placed"].text())


def test_peptide_card_claims_q_per_family_not_per_stratum(tmp_path: Path) -> None:
    text = _card(tmp_path)
    assert "q ≤ 0.01 per family" in text, text
    assert "not q-values" in text, text
    assert "at q ≤ 0.01" not in text.split("per family")[1], text


def test_every_stratum_row_carries_its_full_key(tmp_path: Path) -> None:
    text = _card(tmp_path)
    for r in _rows():
        key = f"{r['decoy_type']} {r['n_diff']} diff {r['length_class']}"
        assert key in text, (key, text)
    assert text.count("diff ") == len(_rows()), text
