"""``workflow/main.nf`` resolves installed database items at the paths Python resolves (MVP-01).

``main.nf`` cannot import ``giap.setupdb``: it re-reads ``<db_dir>/db_manifest.json`` and builds
the path of the two items ``giap run`` passes to the workflow (the M3-15 contaminant library and
the M4-15 COG 2020 table) from a literal. These checks keep the literals equal to
``giap/data/db_sources.yaml``, so a moved file cannot silently drop ``cog_category`` from the gene
table. No Nextflow needed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from giap.setupdb.items import SOURCES

MAIN_NF = Path(__file__).resolve().parents[2] / "workflow" / "main.nf"


def _calls() -> dict[str, str]:
    """``installed('<name>', '<relative path>')`` calls of main.nf."""
    text = re.sub(r"(?m)//.*$", "", MAIN_NF.read_text(encoding="utf-8"))
    return dict(re.findall(r"installed\(\s*'([\w-]+)'\s*,\s*'([^']+)'\s*\)", text))


@pytest.mark.parametrize("name", ["contaminants", "cog2020_def"])
def test_installed_path_equals_db_sources(name: str) -> None:
    calls = _calls()
    assert name in calls, f"main.nf does not resolve {name}: {calls}"
    assert calls[name] == SOURCES[name].path, (
        f"main.nf resolves {name} at {calls[name]}, db_sources.yaml at {SOURCES[name].path}"
    )
