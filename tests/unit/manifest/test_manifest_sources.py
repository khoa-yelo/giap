"""MVP-03a readers and writer in ``giap.io.manifest_sources``
(docs/specs/MVP-03-manifest-finalize.md D-B step 1, D-C ``tools``/``inputs``/``git`` rows, edge
cases "An event line that is not JSON", "versions.yml missing", "Non-git install"; lead
decisions 5 and 9).

Interface assumed (test-author names; the spec names the module and its contents only):

* ``rotate_events(logs_dir: Path) -> None``: ``logs_dir/events.jsonl`` → ``events.previous.jsonl``
  (one generation: an older previous file is replaced); no current file → nothing raised.
* ``read_events(path: Path) -> list[dict]``: one dict per non-blank line, file order; a missing file
  → ``[]``; a line that is not a JSON object raises ``ValueError``.
* ``read_trace(path: Path) -> list[dict[str, str]] | None``: ``trace.tsv`` rows as
  ``csv.DictReader`` dicts; a missing file → ``None``.
* ``tool_reports(trace_rows, work_dir: Path, locks_dir: Path) -> list[dict]``: for each COMPLETED
  or CACHED row, the ``versions.yml`` in the row's task directory (``workdir`` column) gives one
  ``{process, name, version, container_digest, conda_lock_hash}`` per tool line; ``process`` is the
  ``versions.yml`` key (``task.process``). ``container_digest`` = ``sha256:<hex>`` taken from the
  trace ``container`` column (``...@sha256:<hex>`` or the Apptainer cache name
  ``...@sha256-<hex>.img``), else None; when the task had no container (``-``/empty),
  ``conda_lock_hash`` = SHA-256 of ``locks_dir/<tool>.explicit.txt`` when that file exists, else
  None. A missing ``versions.yml``
  is skipped silently; FAILED rows are ignored.
* ``input_row(role: str, path: Path) -> dict``: ``{role, file, path, sha256, size}``.
* ``git_info(root: Path) -> tuple[str | None, bool | None]``: ``(HEAD sha, dirty)``; ``(None,
  None)`` outside a git checkout.
* ``write_manifest(path: Path, manifest) -> None``: writes ``giap.manifest.canonical_json``;
  NaN/±∞ raise ``ValueError``.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest
from helpers.manifest_data import trace_row

DIGEST = "9" * 64


def _ms():  # type: ignore[no-untyped-def]
    from giap.io import manifest_sources

    return manifest_sources


# ----------------------------------------------------------------------------- rotation


def test_rotate_moves_events_to_previous(tmp_path: Path) -> None:
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "events.jsonl").write_text('{"run": 1}\n')
    _ms().rotate_events(logs)
    assert not (logs / "events.jsonl").exists()
    assert (logs / "events.previous.jsonl").read_text() == '{"run": 1}\n'


def test_rotate_keeps_one_generation(tmp_path: Path) -> None:
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "events.previous.jsonl").write_text('{"run": 0}\n')
    (logs / "events.jsonl").write_text('{"run": 1}\n')
    _ms().rotate_events(logs)
    assert (logs / "events.previous.jsonl").read_text() == '{"run": 1}\n'
    assert not list(logs.glob("events*.jsonl.*"))


def test_rotate_without_events_is_a_no_op(tmp_path: Path) -> None:
    _ms().rotate_events(tmp_path / "logs")  # directory absent too
    (tmp_path / "logs").mkdir()
    _ms().rotate_events(tmp_path / "logs")
    assert not (tmp_path / "logs" / "events.jsonl").exists()


# ----------------------------------------------------------------------------- events


def test_read_events_in_file_order(tmp_path: Path) -> None:
    p = tmp_path / "events.jsonl"
    p.write_text('{"event": "b"}\n\n{"event": "a"}\n')
    assert _ms().read_events(p) == [{"event": "b"}, {"event": "a"}]


def test_read_events_missing_file(tmp_path: Path) -> None:
    assert _ms().read_events(tmp_path / "none.jsonl") == []


@pytest.mark.parametrize("bad", ["not json", "[1, 2]", '{"a": 1'])
def test_read_events_bad_line_is_an_error(tmp_path: Path, bad: str) -> None:
    p = tmp_path / "events.jsonl"
    p.write_text(f'{{"event": "a"}}\n{bad}\n')
    with pytest.raises(ValueError):
        _ms().read_events(p)


# ----------------------------------------------------------------------------- trace


def test_read_trace(tmp_path: Path) -> None:
    p = tmp_path / "trace.tsv"
    p.write_text("name\tstatus\trealtime\nPROMOTECH (chunk_1)\tCOMPLETED\t5.5s\n")
    assert _ms().read_trace(p) == [
        {"name": "PROMOTECH (chunk_1)", "status": "COMPLETED", "realtime": "5.5s"}
    ]


def test_read_trace_missing(tmp_path: Path) -> None:
    assert _ms().read_trace(tmp_path / "trace.tsv") is None


# ----------------------------------------------------------------------------- tools


def _task(tmp_path: Path, name: str, versions: str | None) -> Path:
    d = tmp_path / "work" / name
    d.mkdir(parents=True)
    if versions is not None:
        (d / "versions.yml").write_text(versions)
    return d


def test_tool_reports_with_container_digest(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", '"ANNOTATE:GENOMAD":\n    genomad: 1.12.0\n')
    ref = f"docker://quay.io/biocontainers/genomad@sha256:{DIGEST}"
    row = trace_row("ANNOTATE:GENOMAD (genome)", workdir=str(d), container=ref)
    reps = _ms().tool_reports([row], tmp_path / "work", tmp_path / "locks")
    assert reps == [{"process": "ANNOTATE:GENOMAD", "name": "genomad", "version": "1.12.0",
                     "container_digest": f"sha256:{DIGEST}", "conda_lock_hash": None}]  # fmt: skip


def test_tool_reports_apptainer_cache_name(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", '"ANNOTATE:BAKTA":\n    bakta: 1.11.4\n')
    img = f"/g/containers/quay.io-biocontainers-bakta@sha256-{DIGEST}.img"
    row = trace_row("ANNOTATE:BAKTA (genome)", "CACHED", workdir=str(d), container=img)
    (rep,) = _ms().tool_reports([row], tmp_path / "work", tmp_path / "locks")
    assert rep["container_digest"] == f"sha256:{DIGEST}"


def test_tool_reports_conda_lock_hash_without_container(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", '"ANNOTATE:PROMOTECH":\n    promotech: 1.0\n')
    locks = tmp_path / "locks"
    locks.mkdir()
    (locks / "promotech.explicit.txt").write_text("@EXPLICIT\nhttps://x/y.conda\n")
    row = trace_row("ANNOTATE:PROMOTECH (chunk_1)", workdir=str(d), container="-")
    (rep,) = _ms().tool_reports([row], tmp_path / "work", locks)
    want = hashlib.sha256((locks / "promotech.explicit.txt").read_bytes()).hexdigest()
    assert rep["container_digest"] is None
    assert rep["conda_lock_hash"] == want


def test_tool_reports_skip_missing_versions_and_failed_rows(tmp_path: Path) -> None:
    gone = _task(tmp_path, "cleaned", None)
    failed = _task(tmp_path, "failed", '"A:PROMOTECH":\n    promotech: 1.0\n')
    rows = [trace_row("A:BAKTA (g)", workdir=str(gone)),
            trace_row("A:PROMOTECH (c1)", "FAILED", workdir=str(failed))]  # fmt: skip
    assert _ms().tool_reports(rows, tmp_path / "work", tmp_path / "locks") == []


def test_tool_reports_every_tool_line(tmp_path: Path) -> None:
    d = _task(tmp_path, "ab", '"NOVEL:TBLASTN":\n    blast: 2.16.0\n    python: 3.12.1\n')
    row = trace_row("NOVEL:TBLASTN (x)", workdir=str(d))
    reps = _ms().tool_reports([row], tmp_path / "work", tmp_path / "locks")
    assert sorted((r["name"], r["version"]) for r in reps) == [("blast", "2.16.0"),
                                                               ("python", "3.12.1")]  # fmt: skip


# ----------------------------------------------------------------------------- inputs, git


def test_input_row(tmp_path: Path) -> None:
    p = tmp_path / "depth.tsv"
    p.write_bytes(b"contig\tmean_depth\nchr1\t30\n")
    row = _ms().input_row("replicon_depth", p)
    assert row == {"role": "replicon_depth", "file": "depth.tsv", "path": str(p),
                   "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                   "size": p.stat().st_size}  # fmt: skip


def test_git_info_outside_a_checkout(tmp_path: Path) -> None:
    assert _ms().git_info(tmp_path) == (None, None)


@pytest.mark.skipif(shutil.which("git") is None, reason="git not on PATH")
def test_git_info_clean_then_dirty(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.invalid",
           "-c", "commit.gpgsign=false"]  # fmt: skip
    subprocess.run([*git, "init", "-q"], check=True)
    (repo / "a.txt").write_text("a\n")
    subprocess.run([*git, "add", "a.txt"], check=True)
    subprocess.run([*git, "commit", "-q", "-m", "init"], check=True)
    head = subprocess.run([*git, "rev-parse", "HEAD"], check=True, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip
    assert _ms().git_info(repo) == (head, False)
    (repo / "a.txt").write_text("changed\n")
    assert _ms().git_info(repo) == (head, True)


# ----------------------------------------------------------------------------- writer


def test_write_manifest_is_canonical(tmp_path: Path) -> None:
    from giap.manifest import canonical_json

    doc = {"b": 1, "a": ["é", None]}
    p = tmp_path / "run_manifest.json"
    _ms().write_manifest(p, doc)
    assert p.read_text(encoding="utf-8") == canonical_json(doc)
    assert json.loads(p.read_text(encoding="utf-8")) == doc


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_write_manifest_rejects_non_finite(tmp_path: Path, bad: float) -> None:
    with pytest.raises(ValueError):
        _ms().write_manifest(tmp_path / "run_manifest.json", {"x": bad})
