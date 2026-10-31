"""``scripts/ci/check_roles.py``: the trailer exemption of the bot's digest commit is narrow."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "ci" / "check_roles.py"
PIN_FILE = "workflow/conf/containers.config"
BOT = "github-actions[bot]"
MISSING = "missing Agent-Role trailer"


def _git(path: Path, *args: str, author: str = "t") -> str:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "GITHUB_"))}
    env |= {
        "GIT_AUTHOR_NAME": author,
        "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": author,
        "GIT_COMMITTER_EMAIL": "t@example.invalid",
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    return subprocess.run(
        ["git", "-c", "commit.gpgsign=false", *args],
        cwd=path, env=env, capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip


def _check(tmp_path: Path, files: list[str], author: str) -> subprocess.CompletedProcess[str]:
    """``check_roles.py base`` on one trailer-less commit by ``author`` that changes ``files``."""
    _git(tmp_path, "init", "-q", "-b", "main")
    (tmp_path / "README").write_text("base\n")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-qm", "base")
    _git(tmp_path, "branch", "base")
    for name in files:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("pin\n")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-qm", "ci: pin GHCR digests", author=author)
    return subprocess.run(
        [sys.executable, str(SCRIPT), "base"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )  # fmt: skip


def test_bot_commit_of_only_the_pin_files_needs_no_trailer(tmp_path: Path) -> None:
    files = [PIN_FILE, "workflow/conf/container_common.config", "workflow/conf/images.json"]
    proc = _check(tmp_path, files, BOT)
    assert proc.returncode == 0, proc.stdout


@pytest.mark.parametrize(
    ("files", "author"),
    [
        ([PIN_FILE], "t"),
        ([PIN_FILE, "src/gippa/x.py"], BOT),
        ([PIN_FILE, "tests/unit/test_x.py"], BOT),
        (["workflow/main.nf"], BOT),
    ],
    ids=["another-author", "bot-with-src", "bot-with-tests", "bot-other-workflow-file"],
)
def test_any_other_trailerless_commit_still_fails(
    tmp_path: Path, files: list[str], author: str
) -> None:
    proc = _check(tmp_path, files, author)
    assert proc.returncode == 1
    assert MISSING in proc.stdout
