"""PR #236 review (MVP-03c), lead fixes B1-B2 and the decision-20 order.

* B1: an unreadable ``versions.yml`` (stderr in a ``--version`` line, a second document line,
  a non-UTF-8 byte) is skipped like a missing one; the other tasks' tools are still reported.
* B2: a container digest wins over an existing conda lock; ``1.10`` stays the string ``"1.10"``.
* Decision 20: of distinct payloads A, B, A the last in file order (A) is the one callers read.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from helpers.manifest_data import trace_row

from gippa.io import manifest_sources as ms
from gippa.manifest import _event_data

DIGEST = "9" * 64
BAD = {
    "stderr_warning": b'"A:MINIPROT":\n    miniprot: WARNING: x 0.18\n',
    "second_line": b'"A:MINIPROT":\n    miniprot: 0.18\n  stray: [\n',
    "not_utf8": b'"A:MINIPROT":\n    miniprot: 0.18\xff\n',
}


def _task(root: Path, name: str, versions: bytes) -> Path:
    d = root / "work" / name
    d.mkdir(parents=True)
    (d / "versions.yml").write_bytes(versions)
    return d


@pytest.mark.parametrize("kind", sorted(BAD))
def test_unreadable_versions_yml_is_skipped(tmp_path: Path, kind: str) -> None:
    bad = _task(tmp_path, "bad", BAD[kind])
    good = _task(tmp_path, "good", b'"A:PROMOTECH":\n    promotech: 1.0\n')
    rows = [trace_row("A:MINIPROT (g)", workdir=str(bad)),
            trace_row("A:PROMOTECH (c1)", workdir=str(good))]  # fmt: skip
    reps = ms.tool_reports(rows, tmp_path / "work", tmp_path / "locks")
    assert [(r["name"], r["version"]) for r in reps] == [("promotech", "1.0")]


def test_container_digest_wins_over_an_existing_lock(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", b'"ANNOTATE:GENOMAD":\n    genomad: 1.12.0\n')
    locks = tmp_path / "locks"
    locks.mkdir()
    (locks / "genomad.explicit.txt").write_text("@EXPLICIT\nhttps://x/y.conda\n")
    ref = f"docker://quay.io/biocontainers/genomad@sha256:{DIGEST}"
    row = trace_row("ANNOTATE:GENOMAD (genome)", workdir=str(d), container=ref)
    (rep,) = ms.tool_reports([row], tmp_path / "work", locks)
    assert rep["container_digest"] == f"sha256:{DIGEST}"
    assert rep["conda_lock_hash"] is None


def test_version_1_10_stays_a_string(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", b'"NOVEL:TBLASTN":\n    blast: 1.10\n')
    row = trace_row("NOVEL:TBLASTN (x)", workdir=str(d))
    (rep,) = ms.tool_reports([row], tmp_path / "work", tmp_path / "locks")
    assert rep["version"] == "1.10"


def test_last_distinct_payload_comes_last() -> None:
    a, b = {"search_db": "a", "n_queries": 1}, {"search_db": "b", "n_queries": 2}
    events = [{"event": "mapping.search_db_resolved", "data": d} for d in (a, b, a)]
    got = _event_data(SimpleNamespace(events=events), "mapping.search_db_resolved")  # type: ignore[arg-type]
    assert got == [b, a]
