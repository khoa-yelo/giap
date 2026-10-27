"""Docker / apptainer profile parity, read statically (M11-14 item 3b, acceptance 11; PR 14p).

Spec: docs/specs/M11-14-container-install.md, "3b. Profile parity: docker and plain apptainer".
``workflow/nextflow.config`` is parsed as text: each profile body is cut out of ``profiles { }``
and every ``includeConfig`` inside it is inlined (paths relative to the including file, as
Nextflow resolves them), so a block shared through an included file counts for both profiles.

Asserted for ``docker`` and ``apptainer``:

- the DB bind (``-v``/``-B ${params.db_dir}:${params.db_dir}:ro``);
- the ps shim on exactly ``BAKTA|ANTISMASH``
  (``-v``/``-B ${projectDir}/assets/ps/ps:/usr/bin/ps:ro``);
- the per-task ``TMPDIR`` ``beforeScript``;
- the pinned-process conda-off selector, written once (one shared block both profiles include),
  naming exactly the processes pinned in ``containers.config``;
- for docker, ``-u $(id -u):$(id -g)`` and no Apptainer ``-B`` syntax.

``update_digests.py --nextflow-config workflow/nextflow.config`` (as ``containers.yml`` calls it)
keeps editing that one selector, so a newly pinned process loses conda under both profiles.
No network, no Docker; the update runs on copies in a temporary directory.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
NEXTFLOW = WF / "nextflow.config"
CONTAINERS = WF / "conf" / "containers.config"
SCRIPT = REPO / "scripts" / "release" / "update_digests.py"

TMPDIR_SCRIPT = "mkdir -p tmp; export TMPDIR=$PWD/tmp"
DB_BIND = {
    "docker": "-v ${params.db_dir}:${params.db_dir}:ro",
    "apptainer": "-B ${params.db_dir}:${params.db_dir}:ro",
}
PS_SHIM = {
    "docker": "-v ${projectDir}/assets/ps/ps:/usr/bin/ps:ro",
    "apptainer": "-B ${projectDir}/assets/ps/ps:/usr/bin/ps:ro",
}
PROFILES = ("docker", "apptainer")
CONDA_OFF = re.compile(r"withName:\s*'(?P<names>[A-Za-z0-9_|]+)'\s*\{\s*conda\s*=\s*''\s*\}")
CONTAINER_OPTS = re.compile(
    r"withName:\s*'(?P<names>[A-Za-z0-9_|]+)'\s*\{\s*containerOptions\s*=\s*(?P<q>[\"'])(?P<val>.*?)(?P=q)\s*\}"
)
INCLUDE = re.compile(r"""^\s*includeConfig\s+['"](?P<path>[^'"]+)['"]\s*$""", re.M)


# --- a minimal Nextflow-config reader ------------------------------------------------------------


def _strip_comments(text: str) -> str:
    """Drop ``//`` comments (after whitespace or at line start; ``docker://`` stays)."""
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)


def _block(text: str, header: re.Pattern[str]) -> str:
    """The body of the first ``<header> {`` block, braces balanced."""
    m = header.search(text)
    assert m, f"no block matching {header.pattern!r}"
    start = text.index("{", m.start())
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i]
    raise AssertionError(f"unbalanced block {header.pattern!r}")


def _inline_includes(text: str, base: Path, seen: tuple[Path, ...] = ()) -> str:
    def sub(m: re.Match[str]) -> str:
        path = (base / m.group("path")).resolve()
        assert path not in seen, f"includeConfig cycle at {path}"
        assert path.is_file(), f"includeConfig {m.group('path')!r} not found from {base}"
        body = _strip_comments(path.read_text())
        return _inline_includes(body, path.parent, (*seen, path))

    return INCLUDE.sub(sub, text)


def _profile(name: str, nf: Path = NEXTFLOW) -> str:
    """The text of profile ``name`` with its includeConfig files inlined."""
    text = _strip_comments(nf.read_text())
    profiles = _block(text, re.compile(r"^profiles\s*\{", re.M))
    body = _block(profiles, re.compile(rf"(^|\s){re.escape(name)}\s*\{{", re.M))
    return _inline_includes(body, nf.parent)


def _conda_off(text: str) -> set[str]:
    names: set[str] = set()
    for m in CONDA_OFF.finditer(text):
        names |= set(m.group("names").split("|"))
    return names


def _shim(text: str) -> dict[frozenset[str], str]:
    return {
        frozenset(m.group("names").split("|")): m.group("val")
        for m in CONTAINER_OPTS.finditer(text)
    }


def _pinned(containers: Path = CONTAINERS) -> set[str]:
    return set(
        re.findall(r"withName:\s*'([A-Z0-9_]+)'\s*\{\s*container\s*=", containers.read_text())
    )


def _config_files(root: Path) -> list[Path]:
    return [root / "nextflow.config", *sorted((root / "conf").rglob("*.config"))]


# --- the parser itself, on small examples ---------------------------------------------------------


def test_reader_inlines_profile_includes(tmp_path: Path) -> None:
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf" / "shared.config").write_text(
        "process {\n    withName: 'A|B' { conda = '' }\n}\n"
    )
    nf = tmp_path / "nextflow.config"
    nf.write_text(
        "profiles {\n"
        "    one { includeConfig 'conf/shared.config' }\n"
        "    two {\n        x.y = 1   // { not a brace that counts\n"
        "        includeConfig 'conf/shared.config'\n    }\n"
        "    three { z = 2 }\n"
        "}\n"
    )
    assert _conda_off(_profile("one", nf)) == {"A", "B"}
    assert _conda_off(_profile("two", nf)) == {"A", "B"}
    assert _conda_off(_profile("three", nf)) == set()


# --- item 3b table, per profile -------------------------------------------------------------------


@pytest.mark.parametrize("profile", PROFILES)
def test_db_dir_bound_read_only(profile: str) -> None:
    body = _profile(profile)
    assert DB_BIND[profile] in body, f"{profile}: no read-only bind of params.db_dir"
    assert re.search(r"params\.db_dir\s*\?", body), (
        f"{profile}: the DB bind must apply only when params.db_dir is set"
    )


@pytest.mark.parametrize("profile", PROFILES)
def test_ps_shim_on_bakta_and_antismash(profile: str) -> None:
    shims = _shim(_profile(profile))
    key = frozenset({"BAKTA", "ANTISMASH"})
    assert key in shims, f"{profile}: no containerOptions selector for BAKTA|ANTISMASH: {shims}"
    assert shims[key] == PS_SHIM[profile]


@pytest.mark.parametrize("profile", PROFILES)
def test_per_task_tmpdir_before_script(profile: str) -> None:
    body = _profile(profile)
    found = re.findall(r"beforeScript\s*=\s*'([^']*)'", body)
    assert found == [TMPDIR_SCRIPT], f"{profile}: beforeScript {found}"


@pytest.mark.parametrize("profile", PROFILES)
def test_pinned_processes_never_use_conda(profile: str) -> None:
    assert _conda_off(_profile(profile)) == _pinned()


def test_docker_runs_as_the_calling_user() -> None:
    body = _profile("docker")
    assert "-u $(id -u):$(id -g)" in body
    m = re.search(r"runOptions\s*=\s*(.+)", body)
    assert m and "-u $(id -u):$(id -g)" in m.group(1), "the -u must be in docker.runOptions"


def test_docker_profile_uses_docker_syntax_only() -> None:
    body = _profile("docker")
    assert "docker.enabled" in body or re.search(r"docker\s*\{[^}]*enabled\s*=\s*true", body)
    assert not re.search(r"(^|[\s\"'])-B\s", body), "Apptainer -B bind in the docker profile"
    assert "apptainer." not in body


def test_apptainer_profile_keeps_apptainer_syntax() -> None:
    body = _profile("apptainer")
    assert "apptainer.enabled" in body
    assert not re.search(r"(^|[\s\"'])-v\s", body), "Docker -v bind in the apptainer profile"
    assert "-u $(id -u)" not in body


# --- one shared selector ----------------------------------------------------------------------


def test_conda_off_selector_written_once_and_shared() -> None:
    hits = [
        (p.relative_to(REPO), m.group(0))
        for p in _config_files(WF)
        for m in CONDA_OFF.finditer(_strip_comments(p.read_text()))
    ]
    assert len(hits) == 1, f"expected one conda-off selector in workflow/ config files: {hits}"
    for profile in PROFILES:
        assert CONDA_OFF.search(_profile(profile)), f"{profile} does not include the selector"


def test_tmpdir_and_ps_shim_not_duplicated_per_profile() -> None:
    """Shared items stay in sync: the beforeScript line exists once across the config files."""
    texts = [_strip_comments(p.read_text()) for p in _config_files(WF)]
    count = sum(t.count(TMPDIR_SCRIPT) for t in texts)
    assert count == 1, f"the TMPDIR beforeScript is written {count} times; share one block"


def test_conda_profile_alone_keeps_conda_for_pinned_tools() -> None:
    """The shared block is included by the container profiles only, not applied globally."""
    top = _strip_comments(NEXTFLOW.read_text())
    outside = top.replace(_block(top, re.compile(r"^profiles\s*\{", re.M)), "")
    assert not CONDA_OFF.search(_inline_includes(outside, WF)), "conda-off outside the profiles"
    assert _conda_off(_profile("conda")) == set()


def _ref(tool: str, hexd: str) -> str:
    return f"ghcr.io/khoa-yelo/gippa-{tool}@sha256:{hexd}"


def test_update_digests_edits_the_shared_selector(tmp_path: Path) -> None:
    """``--nextflow-config workflow/nextflow.config`` (containers.yml) reaches the one selector."""
    wf = tmp_path / "workflow"
    (wf / "conf").mkdir(parents=True)
    shutil.copy2(NEXTFLOW, wf / "nextflow.config")
    for p in (WF / "conf").rglob("*.config"):
        dest = wf / "conf" / p.relative_to(WF / "conf")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    digests = tmp_path / "digests.json"
    digests.write_text(json.dumps({"OSTIR": _ref("ostir", "a" * 64)}))
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--digests",
            str(digests),
            "--config",
            str(wf / "conf" / "containers.config"),
            "--nextflow-config",
            str(wf / "nextflow.config"),
            "--date",
            "2026-10-25",
            "--source-commit",
            "0123456789abcdef0123456789abcdef01234567",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
    pinned = _pinned(wf / "conf" / "containers.config")
    assert "OSTIR" in pinned
    for profile in PROFILES:
        assert _conda_off(_profile(profile, wf / "nextflow.config")) == pinned, profile
    hits = sum(len(CONDA_OFF.findall(_strip_comments(p.read_text()))) for p in _config_files(wf))
    assert hits == 1
