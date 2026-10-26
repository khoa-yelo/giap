"""``update_digests.py --nextflow-config`` and the CLI guards (M11-06, #338 review B1).

Spec docs/specs/M11-06-ghcr-containers.md item 2 (lead decision on ambiguity 5): the pinned set
equals the apptainer profile's conda-off set, so every newly pinned process joins the single
``withName: '…' { conda = '' }`` selector in the same run. Since M11-14p the selector lives in
``workflow/conf/container_common.config``, which the apptainer and docker profiles include; the
script follows ``includeConfig`` from ``workflow/nextflow.config``.
No network: every run edits copies in a temporary directory.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "update_digests.py"
CONFIG = REPO / "workflow" / "conf" / "containers.config"
NEXTFLOW = REPO / "workflow" / "nextflow.config"
COMMON = REPO / "workflow" / "conf" / "container_common.config"

DATE = "2026-10-25"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
SELECTOR = re.compile(r"^\s*withName: '(?P<names>[A-Z0-9_|]+)' \{ conda = '' \}$", re.M)


def _ref(tool: str, hexd: str) -> str:
    return f"ghcr.io/khoa-yelo/giap-{tool}@sha256:{hexd}"


def _run(
    tmp: Path,
    digests: dict[str, str],
    *,
    common_text: str | None = None,
    date: str = DATE,
    commit: str = COMMIT,
    config_text: str | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path, Path]:
    cfg = tmp / "containers.config"
    cfg.write_text(CONFIG.read_text() if config_text is None else config_text)
    nf = tmp / "nextflow.config"
    nf.write_text(NEXTFLOW.read_text())
    conf = tmp / "conf"
    conf.mkdir(exist_ok=True)
    for src in (REPO / "workflow" / "conf").glob("*.config"):
        (conf / src.name).write_text(src.read_text())
    common = conf / "container_common.config"
    common.write_text(COMMON.read_text() if common_text is None else common_text)
    dj = tmp / "digests.json"
    dj.write_text(json.dumps(digests))
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--digests",
            str(dj),
            "--config",
            str(cfg),
            "--nextflow-config",
            str(nf),
            "--date",
            date,
            "--source-commit",
            commit,
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    return proc, cfg, common


def _selector_names(text: str) -> list[str]:
    hits = SELECTOR.findall(text)
    assert len(hits) == 1, hits
    return hits[0].split("|")


def test_real_config_has_one_selector_in_the_shared_file() -> None:
    assert not SELECTOR.findall(NEXTFLOW.read_text())
    _selector_names(COMMON.read_text())


def test_pinned_processes_join_the_conda_off_selector(tmp_path: Path) -> None:
    before = _selector_names(COMMON.read_text())
    proc, _, nf = _run(
        tmp_path, {"OSTIR": _ref("ostir", "a" * 64), "OSTIR_NOVEL": _ref("ostir", "a" * 64)}
    )
    assert proc.returncode == 0, proc.stderr
    after = _selector_names(nf.read_text())
    assert after == [*before, "OSTIR", "OSTIR_NOVEL"]


def test_only_the_selector_line_changes(tmp_path: Path) -> None:
    original = COMMON.read_text()
    proc, _, nf = _run(tmp_path, {"TERMINATORNET": _ref("terminatornet", "b" * 64)})
    assert proc.returncode == 0, proc.stderr
    old, new = original.splitlines(), nf.read_text().splitlines()
    assert len(old) == len(new)
    changed = [i for i, (a, b) in enumerate(zip(old, new, strict=True)) if a != b]
    assert len(changed) == 1 and SELECTOR.match(new[changed[0]])
    assert nf.read_text().endswith("\n") == original.endswith("\n")


def test_second_run_does_not_duplicate_names(tmp_path: Path) -> None:
    digests = {"BACTERMFINDER": _ref("bactermfinder", "c" * 64)}
    proc, _, nf = _run(tmp_path, digests)
    assert proc.returncode == 0, proc.stderr
    once = nf.read_text()
    proc2, _, nf2 = _run(tmp_path, digests, common_text=once)
    assert proc2.returncode == 0, proc2.stderr
    assert nf2.read_text() == once
    assert _selector_names(once).count("BACTERMFINDER") == 1


@pytest.mark.parametrize("copies", [0, 2])
def test_zero_or_two_selectors_exit_1_and_leave_both_files(tmp_path: Path, copies: int) -> None:
    lines = COMMON.read_text().splitlines(keepends=True)
    idx = next(i for i, ln in enumerate(lines) if SELECTOR.match(ln.rstrip("\n")))
    if copies == 0:
        lines[idx] = "// selector removed\n"
    else:
        lines.insert(idx, lines[idx])
    nf_text = "".join(lines)
    cfg_text = CONFIG.read_text()
    proc, cfg, nf = _run(tmp_path, {"OSTIR": _ref("ostir", "d" * 64)}, common_text=nf_text)
    assert proc.returncode == 1
    assert "selector" in proc.stderr
    assert nf.read_text() == nf_text
    assert cfg.read_text() == cfg_text


def test_duplicate_pin_in_config_exits_1(tmp_path: Path) -> None:
    line = f"  withName: 'OSTIR' {{ container = 'docker://{_ref('ostir', 'e' * 64)}' }}"
    text = CONFIG.read_text()
    lines = text.splitlines()
    start = lines.index("process {")
    lines[start + 1 : start + 1] = [line, line]
    bad = "\n".join(lines) + "\n"
    proc, cfg, _ = _run(tmp_path, {"OSTIR": _ref("ostir", "f" * 64)}, config_text=bad)
    assert proc.returncode == 1
    assert "more than once" in proc.stderr
    assert cfg.read_text() == bad


@pytest.mark.parametrize("date", ["2026-13-01", "25-10-2026", "today", ""])
def test_bad_date_exits_1(tmp_path: Path, date: str) -> None:
    proc, cfg, _ = _run(tmp_path, {"OSTIR": _ref("ostir", "a" * 64)}, date=date)
    assert proc.returncode == 1
    assert cfg.read_text() == CONFIG.read_text()


@pytest.mark.parametrize("commit", ["abc", "G" * 40, "0123456789ABCDEF0123", "a" * 41])
def test_bad_source_commit_exits_1(tmp_path: Path, commit: str) -> None:
    proc, cfg, _ = _run(tmp_path, {"OSTIR": _ref("ostir", "a" * 64)}, commit=commit)
    assert proc.returncode == 1
    assert cfg.read_text() == CONFIG.read_text()


def test_tag_lines_follow_tool_order_after_existing_ones(tmp_path: Path) -> None:
    proc, cfg, _ = _run(
        tmp_path,
        {
            "OSTIR": _ref("ostir", "a" * 64),
            "PROMOTECH": _ref("promotech", "b" * 64),
            "PROMOTECH_NOVEL": _ref("promotech", "b" * 64),
        },
    )
    assert proc.returncode == 0, proc.stderr
    lines = cfg.read_text().splitlines()
    header = next(i for i, ln in enumerate(lines) if ln.startswith("// Tags the digests"))
    block = []
    for ln in lines[header + 1 :]:
        if not ln.startswith("//   "):
            break
        block.append(ln)
    ours = [ln for ln in block if "giap-promotech:" in ln or "giap-ostir:" in ln]
    assert [("promotech" in ln) for ln in ours] == [True, False]
    assert block[-2:] == ours
