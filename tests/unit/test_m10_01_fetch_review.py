"""M10-01 review round 1 (#290): the run-label guard and a server that ignores byte ranges."""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import ModuleType

import pytest
from helpers.pxd_fixtures import PILOT, load


@pytest.fixture(scope="module")
def pf() -> ModuleType:
    return load(PILOT / "pxd_fetch.py")


def _row(**over: str) -> dict[str, str]:
    body = b"archive bytes\n"
    row = {
        "run": "WT_1",
        "file_name": "WT_1.d.zip",
        "url": "https://example.org/WT_1.d.zip",
        "url_record": "http://example.org/WT_1.d.zip",
        "bytes": str(len(body)),
        "checksum_record": "sha256:" + hashlib.sha256(body).hexdigest(),
        "sha256": "NA",
        "sample": "NA",
        "selected_because": "token:WT",
    }
    return {**row, **over}


@pytest.mark.parametrize("run", ["../../x", "/abs", "a/b", ".hidden", ""])
def test_an_unsafe_run_label_is_a_precondition(pf: ModuleType, run: str) -> None:
    with pytest.raises(pf.Precondition, match="unsafe run label"):
        pf.check_rows([_row(run=run)])


def test_a_refused_range_restarts_the_download(pf: ModuleType, tmp_path: Path) -> None:
    """curl exit 33 (range not supported) keeps the partial file; the retry must start from 0."""
    target = tmp_path / "WT_1.d.zip"
    target.write_bytes(b"archive")  # a partial file from an earlier attempt
    seen: list[bool] = []

    def get(url: str, dest: Path) -> int:
        seen.append(dest.exists())
        if len(seen) == 1:
            return 33
        dest.write_bytes(b"archive bytes\n")
        return 0

    assert pf.fetch(_row(), target, get) is None
    assert seen == [True, False], "the retry must not resume the partial file"
    assert target.read_bytes() == b"archive bytes\n"
