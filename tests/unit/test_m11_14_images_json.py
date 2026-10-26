"""``workflow/conf/images.json`` and ``update_digests.py`` writing it (M11-14 item 2, AC 3).

Spec docs/specs/M11-14-container-install.md, "Lock file" and the edge cases "two processes, one
image" and "an image with no process" (fgsea). Offline: every run edits copies in a temporary
directory; the repository's own files are only read.

Interface fixed by these tests (readings chosen where the spec leaves room):

* ``images.json`` is a JSON object ``{tool: entry}``; each entry has exactly the keys
  ``key``, ``tag``, ``ref`` (bare digest ref), ``processes``, ``built_from`` (commit) and
  ``licence``. The file ends with a newline.
* ``update_digests.py`` gains ``--images PATH`` (default: ``images.json`` next to ``--config``)
  and ``--root DIR`` (the tree the keys are computed from; default the repository). For every
  pinned tool it writes the entry: ``key = giap.io.image_key.key(tool, root)``,
  ``tag = tag(TOOLS[tool].version, key)``, ``ref`` from ``digests.json``,
  ``processes`` = the processes whose ``withName`` line in the new ``containers.config`` carries
  that ref (sorted), ``built_from = --source-commit``, ``licence = TOOLS[tool].licence``.
  Entries of other tools are kept. A second run with the same input changes nothing.
* fgsea has no Nextflow process: ``digests.json`` names it ``FGSEA``; it is written only to
  ``images.json`` (``"processes": []``) and adds no ``withName`` line.
* A rejected batch (M11-06 examples 4-5) leaves ``images.json`` untouched (or absent).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "update_digests.py"
KEY_SCRIPT = REPO / "scripts" / "release" / "image_key.py"
IMAGES = REPO / "workflow" / "conf" / "images.json"
CONFIG = REPO / "workflow" / "conf" / "containers.config"

DATE = "2026-10-25"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
ENV = {
    **os.environ,
    "PYTHONPATH": os.pathsep.join([str(REPO / "src"), os.environ.get("PYTHONPATH", "")]),
}
FIELDS = {"key", "tag", "ref", "processes", "built_from", "licence"}
PROCESSES = {
    "PROMOTECH": "promotech",
    "PROMOTECH_NOVEL": "promotech",
    "PROMOTER_CALCULATOR": "promoter_calculator",
    "TERMINATORNET": "terminatornet",
    "BACTERMFINDER": "bactermfinder",
    "OSTIR": "ostir",
    "OSTIR_NOVEL": "ostir",
}


def _ref(tool: str, hexd: str) -> str:
    return f"ghcr.io/khoa-yelo/giap-{tool}@sha256:{hexd}"


def _tree(root: Path) -> Path:
    for rel in ("workflow/containers", "workflow/envs"):
        shutil.copytree(REPO / rel, root / rel)
    (root / "workflow/bin").mkdir(parents=True)
    shutil.copy2(REPO / "workflow/bin/giap_fgsea.R", root / "workflow/bin/giap_fgsea.R")
    (root / "workflow/conf").mkdir(parents=True)
    shutil.copy2(CONFIG, root / "workflow/conf/containers.config")
    return root


def _run(root: Path, digests: dict[str, str]) -> subprocess.CompletedProcess[str]:
    dj = root / "digests.json"
    dj.write_text(json.dumps(digests))
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--digests",
            str(dj),
            "--config",
            str(root / "workflow/conf/containers.config"),
            "--root",
            str(root),
            "--date",
            DATE,
            "--source-commit",
            COMMIT,
        ],
        capture_output=True,
        text=True,
        cwd=root,
        timeout=120,
        env=ENV,
    )


def _images(root: Path) -> dict[str, dict[str, object]]:
    data: dict[str, dict[str, object]] = json.loads(
        (root / "workflow/conf/images.json").read_text()
    )
    return data


def _entry(root: Path, tool: str, ref: str, processes: list[str]) -> dict[str, object]:
    from giap.annotate.tools import TOOLS
    from giap.io.image_key import key, tag

    k = key(tool, root)
    return {
        "key": k,
        "tag": tag(TOOLS[tool].version, k),
        "ref": ref,
        "processes": processes,
        "built_from": COMMIT,
        "licence": TOOLS[tool].licence,
    }


def _pinned_processes(config: str, ref: str) -> list[str]:
    pat = re.compile(r"^\s*withName: '([A-Z0-9_]+)' \{ container = '([^']+)' \}$")
    return sorted(
        m.group(1) for ln in config.splitlines() if (m := pat.match(ln)) and m.group(2) == ref
    )


# --- the bundled file --------------------------------------------------------------------------


def test_repo_images_json_exists_and_is_well_formed() -> None:
    assert IMAGES.is_file(), "workflow/conf/images.json is missing (M11-14 14a)"
    data = json.loads(IMAGES.read_text())
    assert isinstance(data, dict)
    assert IMAGES.read_text().endswith("\n")
    for tool, entry in data.items():
        assert tool != "uniop", "UniOP is never published (R5)"
        assert set(entry) == FIELDS, tool
        assert re.fullmatch(r"[0-9a-f]{12}", entry["key"]), tool
        assert entry["tag"].endswith(f"-{entry['key']}"), tool
        assert re.fullmatch(rf"ghcr\.io/[^/]+/[a-z]+-{tool}@sha256:[0-9a-f]{{64}}", entry["ref"])


def test_repo_images_json_agrees_with_containers_config() -> None:
    data = json.loads(IMAGES.read_text())
    config = CONFIG.read_text()
    for tool, entry in data.items():
        assert entry["processes"] == _pinned_processes(config, entry["ref"]), tool
        if tool == "fgsea":
            assert entry["processes"] == []


# --- update_digests.py writes it ---------------------------------------------------------------


def test_two_processes_one_image_one_entry_listing_both(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    ref = _ref("ostir", "a" * 64)
    proc = _run(root, {"OSTIR": ref, "OSTIR_NOVEL": ref})
    assert proc.returncode == 0, proc.stderr
    assert _images(root) == {"ostir": _entry(root, "ostir", ref, ["OSTIR", "OSTIR_NOVEL"])}
    assert (root / "workflow/conf/images.json").read_text().endswith("\n")


def test_images_json_defaults_to_the_config_directory(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    proc = _run(root, {"PROMOTECH": _ref("promotech", "b" * 64)})
    assert proc.returncode == 0, proc.stderr
    assert (root / "workflow/conf/images.json").is_file()
    entry = _images(root)["promotech"]
    assert entry["processes"] == ["PROMOTECH"]
    assert re.fullmatch(r"git-56251ad9-[0-9a-f]{12}", str(entry["tag"]))
    assert entry["licence"] == "GPL-3.0"


def test_fgsea_entry_has_no_process_and_no_config_line(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    before = (root / "workflow/conf/containers.config").read_text()
    ref = _ref("fgsea", "c" * 64)
    proc = _run(root, {"FGSEA": ref})
    assert proc.returncode == 0, proc.stderr
    assert _images(root) == {"fgsea": _entry(root, "fgsea", ref, [])}
    after = (root / "workflow/conf/containers.config").read_text()
    assert "withName: 'FGSEA'" not in after
    assert [ln for ln in after.splitlines() if "withName" in ln] == [
        ln for ln in before.splitlines() if "withName" in ln
    ]


def test_other_entries_are_kept(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    p_ref = _ref("promotech", "d" * 64)
    assert _run(root, {"PROMOTECH": p_ref, "PROMOTECH_NOVEL": p_ref}).returncode == 0
    first = _images(root)["promotech"]
    proc = _run(root, {"OSTIR": _ref("ostir", "e" * 64)})
    assert proc.returncode == 0, proc.stderr
    data = _images(root)
    assert data["promotech"] == first
    assert set(data) == {"promotech", "ostir"}


def test_second_run_is_byte_identical(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    digests = {"OSTIR": _ref("ostir", "f" * 64), "FGSEA": _ref("fgsea", "1" * 64)}
    assert _run(root, digests).returncode == 0
    once = (root / "workflow/conf/images.json").read_bytes()
    cfg = (root / "workflow/conf/containers.config").read_bytes()
    proc = _run(root, digests)
    assert proc.returncode == 0, proc.stderr
    assert (root / "workflow/conf/images.json").read_bytes() == once
    assert (root / "workflow/conf/containers.config").read_bytes() == cfg


def test_rejected_batch_writes_no_images_json(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    proc = _run(root, {"OSTIR": _ref("ostir", "2" * 64), "UNIOP": _ref("uniop", "3" * 64)})
    assert proc.returncode == 1
    assert not (root / "workflow/conf/images.json").exists()
    proc = _run(root, {"OSTIR": "ghcr.io/khoa-yelo/giap-ostir@sha256:abc"})
    assert proc.returncode == 1
    assert not (root / "workflow/conf/images.json").exists()


def test_rejected_batch_leaves_an_existing_images_json(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    assert _run(root, {"OSTIR": _ref("ostir", "4" * 64)}).returncode == 0
    once = (root / "workflow/conf/images.json").read_bytes()
    proc = _run(root, {"PROMOTECH": _ref("promotech", "5" * 64), "OSTIR": "bad"})
    assert proc.returncode == 1
    assert (root / "workflow/conf/images.json").read_bytes() == once


def test_fgsea_ref_must_be_its_own_ghcr_image(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    proc = _run(root, {"FGSEA": _ref("ostir", "6" * 64)})
    assert proc.returncode == 1
    assert "FGSEA" in proc.stderr
    assert not (root / "workflow/conf/images.json").exists()


def test_check_passes_after_update_digests_and_fails_after_an_input_change(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    digests = {p: _ref(t, format(i + 7, "x") * 64) for i, (p, t) in enumerate(PROCESSES.items())}
    digests["PROMOTECH_NOVEL"] = digests["PROMOTECH"]
    digests["OSTIR_NOVEL"] = digests["OSTIR"]
    digests["FGSEA"] = _ref("fgsea", "e" * 64)
    proc = _run(root, digests)
    assert proc.returncode == 0, proc.stderr
    check = [sys.executable, str(KEY_SCRIPT), "--root", str(root), "--check"]
    ok = subprocess.run(check, capture_output=True, text=True, timeout=120, env=ENV)
    assert ok.returncode == 0, ok.stdout + ok.stderr
    key = _images(root)["terminatornet"]["key"]
    env = root / "workflow/envs/terminatornet.yml"
    env.write_text(env.read_text() + "# touched\n")
    bad = subprocess.run(check, capture_output=True, text=True, timeout=120, env=ENV)
    assert bad.returncode == 1
    assert f"stale image: terminatornet (inputs changed since {key})" in bad.stdout + bad.stderr


@settings(max_examples=8, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    chosen=st.dictionaries(
        st.sampled_from([*sorted(PROCESSES), "FGSEA"]),
        st.text("0123456789abcdef", min_size=64, max_size=64),
        min_size=1,
    )
)
def test_property_images_json_is_consistent_with_containers_config(chosen: dict[str, str]) -> None:
    # same-image processes share one digest (the workflow pins one digest per image)
    by_tool: dict[str, str] = {}
    for p, h in sorted(chosen.items()):
        by_tool.setdefault(PROCESSES.get(p, "fgsea"), h)
    digests = {
        p: _ref(PROCESSES.get(p, "fgsea"), by_tool[PROCESSES.get(p, "fgsea")]) for p in chosen
    }
    with tempfile.TemporaryDirectory() as d:
        root = _tree(Path(d))
        proc = _run(root, digests)
        assert proc.returncode == 0, proc.stderr
        data = _images(root)
        config = (root / "workflow/conf/containers.config").read_text()
        assert set(data) == {PROCESSES.get(p, "fgsea") for p in digests}
        for tool, entry in data.items():
            assert set(entry) == FIELDS
            assert entry["processes"] == _pinned_processes(config, str(entry["ref"]))
            assert entry == _entry(root, tool, str(entry["ref"]), list(entry["processes"]))
