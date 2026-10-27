"""The ``needs_net`` probe runs only for selected items (tests/conftest.py; PR #115 review).

A child pytest collects a toy test marked ``needs_net`` and ``e2e`` under the repository's
conftest, with ``socket.create_connection`` replaced by a plugin that records each call and
raises ``OSError`` (so nothing ever reaches the network). Deselected by ``-m "not e2e"``: no call.
Selected (positive control): the probe is attempted and the test is skipped as a missing
capability.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

PLUGIN = """
import os, socket

def _blocked(address, *args, **kwargs):
    with open(os.environ["GIPPA_PROBE_LOG"], "a") as fh:
        fh.write(f"{address[0]}\\n")
    raise OSError("network blocked by test")

socket.create_connection = _blocked
"""

TOY = """
import pytest

@pytest.mark.e2e
@pytest.mark.needs_net("probe.invalid")
def test_toy():
    pass
"""


def _run(tmp_path: Path, marker_expr: str) -> tuple[subprocess.CompletedProcess[str], str]:
    plugdir = tmp_path / "plug"
    plugdir.mkdir(exist_ok=True)
    (plugdir / "gippa_block_net.py").write_text(PLUGIN)
    toy_dir = REPO / "tests" / "unit"
    toy = toy_dir / f"test_zz_probe_toy_{os.getpid()}_{tmp_path.name}.py"
    log = tmp_path / f"probe_{marker_expr.replace(' ', '_')}.log"
    env = dict(os.environ, GIPPA_PROBE_LOG=str(log))
    env["PYTHONPATH"] = os.pathsep.join([str(plugdir), env.get("PYTHONPATH", "")])
    toy.write_text(TOY)
    try:
        run = subprocess.run(
            [sys.executable, "-m", "pytest", str(toy), "-m", marker_expr, "-q",
             "-p", "gippa_block_net", "-p", "no:cacheprovider", "-p", "no:xdist", "-rs"],
            cwd=REPO, env=env, capture_output=True, text=True, timeout=300,
        )  # fmt: skip
    finally:
        toy.unlink(missing_ok=True)
    return run, log.read_text() if log.exists() else ""


def test_deselected_needs_net_item_opens_no_socket(tmp_path: Path) -> None:
    run, probes = _run(tmp_path, "not e2e")
    assert "1 deselected" in run.stdout, run.stdout + run.stderr
    assert probes == ""


@pytest.mark.parametrize("expr", ["e2e"])
def test_selected_needs_net_item_is_probed_and_skipped(tmp_path: Path, expr: str) -> None:
    run, probes = _run(tmp_path, expr)
    assert probes.split() == ["probe.invalid"]
    assert "1 skipped" in run.stdout, run.stdout + run.stderr
    assert "needs_net(probe.invalid)" in run.stdout
