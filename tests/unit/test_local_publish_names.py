"""The local executor publishes the same files under the same names as Nextflow (§6.4, §9.1).

The ``key=path`` names live in the Nextflow subworkflows; ``giap.runner.local.PUBLISHED_FILES``
repeats them, so this pins the two together (lead fix after the code review of #145: the local
executor wrote no ``tables/`` at all).
"""

from __future__ import annotations

import re
from pathlib import Path

from giap.runner.local import PUBLISHED_FILES, publish

ROOT = Path(__file__).resolve().parents[2]


def _nextflow_names() -> dict[str, str]:
    names: dict[str, str] = {}
    for nf in sorted((ROOT / "workflow").rglob("*.nf")):
        if "tests" in nf.relative_to(ROOT / "workflow").parts:
            continue
        for key, rel in re.findall(r"'([a-z_0-9]+)=([a-z_]+/[^']+)'", nf.read_text()):
            names[key] = rel
    return names


def test_published_files_equal_the_nextflow_names() -> None:
    assert _nextflow_names() == PUBLISHED_FILES


def test_publish_copies_files_to_their_names_and_skips_absent_ones(tmp_path: Path) -> None:
    work, out = tmp_path / "work", tmp_path / "OUT"
    work.mkdir()
    (work / "regions_stats").write_text("region_id\n")
    done = publish(work, out)
    assert (out / "tables" / "regions_stats.tsv").read_text() == "region_id\n"
    assert out / "tables" / "regions_stats.tsv" in done
    assert not (out / "tables" / "sensitivity.tsv").exists()
