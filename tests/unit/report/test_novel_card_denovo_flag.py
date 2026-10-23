"""The novel-ORF summary card names the de novo-only calls (docs/specs/M7-17-denovo-support-flag.md
"Report", acceptance criterion 3, example 13; ADR-0050 option 3).

Readings pinned by the test author: the suffix counts CALLED rows whose ``flags`` carry
``denovo_only_unsupported`` (uncalled flagged rows do not count); with no flagged call the card's
estimate is exactly the M7-11b text. Card texts are compared whitespace-normalised
(tests/helpers/report_fixtures.py hooks).
"""

from __future__ import annotations

from pathlib import Path

from helpers.report_fixtures import Page, cfg, parse, run_out
from helpers.viz_out_tree import write_schema_table

from giap.io.run_view import RunView

FLAG = "denovo_only_unsupported"
SUFFIX = ("; {k} of {m} calls de novo-only without database or RBS support (exact de novo errors"
          " are not covered by the class FDR)")  # fmt: skip


def _page(out: Path) -> Page:
    from giap.report.build import build_report

    return parse(build_report(RunView.from_outdir(out), {}, {}, cfg()))


def _estimate(page: Page) -> str:
    found = page.cards()["novel_orfs"].by_class("estimate")
    assert found, "no novel_orfs estimate"
    return " ".join(found[0].text().split())


def _novel(out: Path, rows: list[tuple[str, str, bool, str, str | None]]) -> None:
    body = [{"orf_id": i, "class": c, "called": k, "class_fdr": f, "flags": g}
            for i, c, k, f, g in rows]  # fmt: skip
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", body)


# Example 13: examples 1, 2, 3, 6, 7 (1 and 3 flagged and called, 7 flagged and not called).
EX13 = [
    ("N1", "intergenic", True, "0.004", f"rbs_supported;{FLAG}"),
    ("N2", "cds_outframe", True, "0.002", "rbs_supported"),
    ("N3", "cds_outframe", True, "0.002", FLAG),
    ("N6", "intergenic", True, "0.006", "rbs_supported"),
    ("N7", "cds_outframe", False, "0.03", FLAG),
]


def test_example_13_card_suffix(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, EX13)
    assert _estimate(_page(out)) == "class FDR ≤ 0.006" + SUFFIX.format(k=2, m=4)


def test_suffix_follows_unestimated(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", "sorf", True, "unestimated", f"class_fdr_unestimated;{FLAG}"),
                 ("N2", "sorf", True, "unestimated", "class_fdr_unestimated")])  # fmt: skip
    est = _estimate(_page(out))
    assert est.startswith("unestimated"), est
    assert est.endswith(SUFFIX.format(k=1, m=2)), est


def test_no_flagged_call_card_unchanged(tmp_path: Path) -> None:
    """Only an uncalled row carries the flag: the estimate is the M7-11b text exactly."""
    out = run_out(tmp_path)
    rows = [r if r[0] != "N1" and r[0] != "N3" else (*r[:4], "rbs_supported") for r in EX13]
    _novel(out, rows)
    assert _estimate(_page(out)) == "class FDR ≤ 0.006"
