"""The local executor publishes the function tables and the TU folds (M8-06 item 8, 06b).

docs/specs/M8-06-function-step-and-run.md, "Publishing (M6-08 decision 7)": only the three
``tables/*.tsv`` files and ``genetable/tu_folds.parquet`` are named, and
``giap.runner.local.PUBLISHED_FILES`` gains the same four keys (``test_local_publish_names.py``
pins them against the Nextflow names).
"""

from __future__ import annotations

from pathlib import Path

import pytest

WANT = {
    "category_summary": "tables/category_summary.tsv",
    "enrichment_detection": "tables/enrichment_detection.tsv",
    "enrichment_gsea": "tables/enrichment_gsea.tsv",
    "tu_folds": "genetable/tu_folds.parquet",
}


@pytest.mark.parametrize("key", sorted(WANT))
def test_published_files_name_the_function_outputs(key: str) -> None:
    from giap.runner.local import PUBLISHED_FILES

    assert PUBLISHED_FILES.get(key) == WANT[key], (key, PUBLISHED_FILES.get(key))


def test_publish_copies_function_outputs_to_their_names(tmp_path: Path) -> None:
    from giap.runner.local import publish

    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    for key in WANT:
        (work / key).write_bytes(f"{key}\n".encode())
    done = publish(work, out)
    for key, rel in WANT.items():
        assert (out / rel).read_bytes() == f"{key}\n".encode(), rel
        assert out / rel in done


def test_function_off_publishes_no_function_table(tmp_path: Path) -> None:
    """Spec §9.1 "Optional modules that are switched off write no files": nothing in work/,
    nothing published."""
    from giap.runner.local import publish

    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    (work / "operons").write_text("tu_id\n", encoding="utf-8")
    publish(work, out)
    for rel in WANT.values():
        assert not (out / rel).exists(), rel
