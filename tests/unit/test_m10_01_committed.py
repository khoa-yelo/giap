"""What M10-01 commits: the manifest (two stages), ``acquisition.json`` and the README section.

Spec: docs/specs/M10-01-pxd078459-download.md, deliverables 3, 5 and 6, acceptance criteria 6, 8-10,
question Q5 (two commits: the reviewed selection with ``sha256 = NA``, then the filled hashes);
R11 (no raw MS data, no file above 5 MB); index question L4 (the record XML is not committed).

These tests read only files tracked in this checkout; nothing reads ``$GIAP_DATA_DIR``.
No statistical claim is made in this module, so it has no calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

from helpers.pxd_fixtures import (
    ACQUISITION_REL,
    MANIFEST_COLUMNS,
    MANIFEST_REL,
    PILOT,
    read_tsv,
)
from helpers.slurm_sandbox import REPO

R11_MAX_BYTES = 5 * 1024 * 1024
RUNS = ["WT_1", "WT_2", "WT_3"]
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _tracked(prefix: str) -> list[str]:
    res = subprocess.run(["git", "-C", str(REPO), "ls-files", "--", prefix],
                         capture_output=True, text=True, check=True)  # fmt: skip
    return [ln for ln in res.stdout.splitlines() if ln]


def _natural(name: str) -> list[object]:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def test_no_committed_file_under_validation_exceeds_5_mb() -> None:
    big = [p for p in _tracked("validation") if (REPO / p).stat().st_size > R11_MAX_BYTES]
    assert not big, big


def test_no_raw_data_or_record_xml_is_committed() -> None:
    raw = re.compile(r"(\.(d|zip|tdf|tdf_bin|mzml|xml|raw|7z|rar|tar|tar\.gz)$)|(\.d/)", re.I)
    bad = [p for p in _tracked("validation/pilot") if raw.search(p)]
    assert not bad, bad


def test_manifest_is_committed_with_three_runs() -> None:
    path = REPO / MANIFEST_REL
    assert str(MANIFEST_REL) in _tracked(str(MANIFEST_REL)), f"{MANIFEST_REL} is not committed"
    header, rows = read_tsv(path)
    assert tuple(header) == MANIFEST_COLUMNS
    assert sorted({r["run"] for r in rows}) == RUNS
    firsts = []
    for r in rows:
        assert r["url"].startswith(("https://", "ftps://")), r["url"]
        assert r["url_record"] not in ("", "NA"), r
        assert r["file_name"] == r["url_record"].rstrip("/").split("/")[-1] or "/" in r["file_name"]
        assert r["selected_because"] not in ("", "NA"), r
        if r["bytes"] != "NA":
            assert r["bytes"].isdigit(), r["bytes"]
        if r["checksum_record"] != "NA":
            assert re.match(r"^(md5|sha1|sha256):[0-9a-f]+$", r["checksum_record"]), r
        if r["run"] not in [f[0] for f in firsts]:
            firsts.append((r["run"], r["file_name"]))
    assert [f[0] for f in firsts] == RUNS, "rows are in selection order"
    names = [f[1] for f in firsts]
    assert names == sorted(names, key=_natural), "selection is by natural sort (decision 1)"


def test_manifest_sha256_is_all_na_or_all_filled() -> None:
    """Q5: stage 1 commits ``sha256 = NA`` everywhere; stage 2 fills every row."""
    _, rows = read_tsv(REPO / MANIFEST_REL)
    values = [r["sha256"] for r in rows]
    assert all(v == "NA" for v in values) or all(HEX64.match(v) for v in values), values


def test_acquisition_json_comes_with_the_filled_manifest() -> None:
    """Stage 2: ``acquisition.json`` is committed only after the hashes are, and agrees with it."""
    path = REPO / ACQUISITION_REL
    if not path.exists():
        return  # stage 1: the download has not run yet
    _, rows = read_tsv(REPO / MANIFEST_REL)
    assert all(HEX64.match(r["sha256"]) for r in rows), "acquisition.json before the hashes"
    doc = json.loads(path.read_text())
    assert [r["run"] for r in doc["runs"]] == RUNS
    assert doc["dataset_mode"] == "diaPASEF", (
        "acceptance 9: otherwise the issue stops on the §17 path (exit 3), nothing is committed"
    )


def test_readme_has_a_raw_data_section() -> None:
    """Acceptance 10 / deliverable 6."""
    text = (PILOT / "README.md").read_text()
    m = re.search(r"^#+\s*Raw data\s*$(.*?)(?=^#{1,2}\s|\Z)", text, re.M | re.S)
    assert m, "no 'Raw data' section in validation/pilot/README.md"
    sec = m.group(1)
    assert "PXD078459" in sec
    assert "proteomecentral.proteomexchange.org" in sec, "the record URL"
    assert re.search(r"20\d\d-\d\d-\d\d", sec), "the retrieval date"
    assert re.search(r"\b[0-9a-f]{64}\b", sec), "the record sha256"
    assert re.search(r"licen[cs]e", sec, re.I), "the data licence the record states"
    assert re.search(r"natural sort", sec, re.I), "the selection rule"
    for run in RUNS:
        assert run in sec
    assert re.search(r"diaPASEF|ddaPASEF", sec), "the acquisition mode"
    assert re.search(r"\b\d{6,}\b", sec), "the job IDs"


def test_manifest_lives_where_the_sandbox_copies_it() -> None:
    """The M3-24 sandbox copies ``validation/``, so the committed manifest gives the defaults."""
    assert Path(REPO / MANIFEST_REL).is_relative_to(REPO / "validation")
