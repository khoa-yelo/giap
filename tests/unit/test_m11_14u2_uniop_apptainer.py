"""M11-14u2: the apptainer build of UniOP (fakeroot fallback), its record and its preflight row.

docs/specs/M11-14-container-install.md item 1b ("Build per runtime: apptainer (14u2)", "No
fakeroot means both attempts failed", "Permissions", "What is recorded (14u2)"), item 3 (row
"uniop with a local build", apptainer column; the stale-key row), examples 16b, 19d, 20; AC 9
(14u2: "``uniop.def`` has no ``chown`` and no ``curl``/``wget``/URL", "no ``push``"); #396 review S5
(the apptainer path never calls docker; ``local_build`` validates against the schema).

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- Not root: ``apptainer build --fakeroot <sif> <def>``; on a non-zero exit exactly one retry
  ``apptainer build --fakeroot --ignore-fakeroot-command <sif> <def>``. Root (``os.geteuid() ==
  0``): ``apptainer build <sif> <def>`` once. ``<def>`` is the staged ``uniop.def`` (bytes equal to
  the bundled one) and the build runs with the staged context as its working directory, because
  ``%files`` sources are relative to it.
- ``<sif>`` = ``<cache>/localhost-<pkg>-uniop-git-3ed14dfc-<key>.sif``, mode ``0600``.
- ``local_tools.config`` gains ``apptainer { process { withName: 'UNIOP' { container =
  '<absolute sif>'; conda = '' } } }``.
- ``image:uniop``: ``path`` = the ``.sif``; ``local_build`` = ``{"runtime": "apptainer", "ref":
  <sif>, "digest": "sha256:<sif bytes>", "key": key("uniop"), "terms": "flag"}``.
- Both builds failing: ``image:uniop`` ``absent`` with the example 20 message (the last stderr
  line of the second build), exit 2, no apptainer block.
- Preflight (``missing_tool_reasons(..., local_tools=<db>/local_tools.config)``) reads
  ``local_build`` from ``<db>/db_manifest.json`` next to it: present iff the ``.sif`` exists, its
  sha256 equals ``local_build.digest`` and ``local_build.key == key("uniop")``; a stale key's
  reason contains ``local build <old> does not match this release (<new>)``.
"""

from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.uniop_14u1 import (
    ACCEPT,
    SKIP_DBS,
    VERSION,
    FakeNet,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    mode,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
    uniop_key,
)
from helpers.uniop_14u2 import (
    APPTAINER_FAILED,
    GLIBC,
    UNIOP_DEF,
    as_root,
    calls,
    install_fake_runtimes,
    set_local_build,
    sif_path,
    validate_db_manifest,
)

pytestmark = pytest.mark.usefixtures("no_network")

LOCAL_TOOLS = "local_tools.config"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    install_fake_runtimes(tmp_path, monkeypatch)
    old = os.umask(0o022)
    yield
    os.umask(old)


@pytest.fixture
def net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


@pytest.fixture
def archive(net: FakeNet, monkeypatch: pytest.MonkeyPatch) -> bytes:
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    return data


def _uniop(db: Path, cache: Path, *extra: str, profile: str = "apptainer") -> tuple[int, str]:
    return setup_cli(["--db-dir", str(db), "--cache-dir", str(cache), *SKIP_DBS, "--uniop",
                      "--profile", profile, *extra])  # fmt: skip


def _builds(tmp_path: Path) -> list[list[str]]:
    return calls(tmp_path, "apptainer", "build")


def _built(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Example 16b: the first ``--fakeroot`` build fails, the retry succeeds."""
    monkeypatch.setenv("FAKE_APPTAINER_FAIL", "fakeroot")
    db, cache = tmp_path / "db", tmp_path / "cache"
    code, out = _uniop(db, cache, ACCEPT)
    assert code == 0, out
    return db, cache


# ---------------------------------------------------------------- example 16b


def test_example_16b_two_builds_in_order(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, cache = _built(tmp_path, monkeypatch)
    sif = str(sif_path(cache))
    builds = _builds(tmp_path)
    assert len(builds) == 2, builds
    assert builds[0][:-1] == ["build", "--fakeroot", sif], builds
    assert builds[1][:-1] == ["build", "--fakeroot", "--ignore-fakeroot-command", sif], builds
    assert Path(builds[0][-1]).name == "uniop.def" and builds[0][-1] == builds[1][-1]


def test_first_fakeroot_success_builds_once(archive: bytes, tmp_path: Path) -> None:
    code, out = _uniop(tmp_path / "db", tmp_path / "cache", ACCEPT)
    assert code == 0, out
    builds = _builds(tmp_path)
    assert [b[:2] for b in builds] == [["build", "--fakeroot"]], builds
    assert "--ignore-fakeroot-command" not in builds[0]


def test_relative_cache_dir_builds(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--db-dir DB --cache-dir cache`` (relative): the build runs in the context, so the def
    and the sif must be absolute (#401 review B1)."""
    monkeypatch.chdir(tmp_path)
    code, out = _uniop(Path("db"), Path("cache"), ACCEPT)
    assert code == 0, out
    assert sif_path(tmp_path / "cache").is_file(), out


def test_root_builds_without_fakeroot(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    as_root(monkeypatch)
    cache = tmp_path / "cache"
    code, out = _uniop(tmp_path / "db", cache, ACCEPT)
    assert code == 0, out
    (build,) = _builds(tmp_path)
    assert build[:-1] == ["build", str(sif_path(cache))], build


def test_the_def_built_is_the_bundled_one(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _built(tmp_path, monkeypatch)
    seen = tmp_path / "apptainer_state" / "def_seen"
    assert seen.read_bytes() == UNIOP_DEF.read_bytes()


def test_sif_path_and_mode(archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _, cache = _built(tmp_path, monkeypatch)
    sif = sif_path(cache)
    assert sif.is_file()
    assert mode(sif) == 0o600
    assert re.fullmatch(r"localhost-[a-z]+-uniop-git-3ed14dfc-[0-9a-f]{12}\.sif", sif.name)


def test_apptainer_block_names_the_absolute_sif(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.preflight.environment import _local_containers

    db, cache = _built(tmp_path, monkeypatch)
    path = db / LOCAL_TOOLS
    text = path.read_text(encoding="utf-8")
    sif = str(sif_path(cache))
    assert Path(sif).is_absolute()
    assert _local_containers(path, "apptainer") == {"uniop": sif}
    assert _local_containers(path, "docker") == {}
    assert "conda.enabled" not in text
    for line in text.splitlines():
        if re.search(r"\bcontainer\s*=", line):
            assert re.search(r"\bconda\s*=\s*''", line), line


def test_example_16b_record(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, cache = _built(tmp_path, monkeypatch)
    sif = sif_path(cache)
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "ok", item
    assert item["version"] == VERSION
    assert item["checksum"] == f"sha256:{sha256_hex(archive)}"
    assert item["path"] == str(sif)
    assert item["local_build"] == {
        "runtime": "apptainer",
        "ref": str(sif),
        "digest": "sha256:" + hashlib.sha256(sif.read_bytes()).hexdigest(),
        "key": uniop_key(),
        "terms": "flag",
    }
    validate_db_manifest(db)


def test_apptainer_path_never_calls_docker_and_never_pushes(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#396 review S5: with the terms accepted, the apptainer route never touches docker."""
    _built(tmp_path, monkeypatch)
    assert calls(tmp_path, "docker") == []
    assert calls(tmp_path, "apptainer", "push") == []
    assert {c[0] for c in calls(tmp_path, "apptainer")} <= {"build", "exec", "inspect"}


def test_apptainer_build_adds_its_block_and_keeps_docker(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.preflight.environment import _local_containers

    db = tmp_path / "db"
    db.mkdir()
    (db / LOCAL_TOOLS).write_text(
        "profiles {\n    docker { process { withName: 'UNIOP' { container = "
        "'localhost/x-uniop:t'; conda = '' } } }\n}\n"
    )
    code, out = _uniop(db, tmp_path / "cache", ACCEPT)
    assert code == 0, out
    assert _local_containers(db / LOCAL_TOOLS, "docker") == {"uniop": "localhost/x-uniop:t"}
    assert _local_containers(db / LOCAL_TOOLS, "apptainer") == {
        "uniop": str(sif_path(tmp_path / "cache"))
    }


# ---------------------------------------------------------------- example 20: no fakeroot


def test_example_20_both_builds_fail(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    line = "FATAL:   could not use fakeroot: no mapping entry found in /etc/subuid"
    monkeypatch.setenv("FAKE_APPTAINER_FAIL", "all")
    monkeypatch.setenv("FAKE_APPTAINER_FAIL_LINE", line)
    db = tmp_path / "db"
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer",
                           ACCEPT])  # fmt: skip
    assert code == 2, out
    log = calls(tmp_path, "apptainer")
    subs = [c[0] for c in log]
    assert "pull" in subs and subs.count("build") == 2, subs
    assert max(i for i, s in enumerate(subs) if s == "pull") < subs.index("build"), (
        "the other images are installed first"
    )
    items = manifest_items(db)
    assert items["image:bakta"]["status"] == "ok"
    item = items["image:uniop"]
    assert item["status"] == "absent", item
    assert item["message"] == APPTAINER_FAILED.format(line=line), item["message"]
    assert "local_build" not in item
    lt = db / LOCAL_TOOLS
    if lt.exists():
        from gippa.preflight.environment import _local_containers

        assert _local_containers(lt, "apptainer") == {}


def test_glibc_failure_of_both_names_the_last_line(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_APPTAINER_FAIL", "all")
    monkeypatch.setenv("FAKE_APPTAINER_FAIL_LINE", GLIBC)
    db = tmp_path / "db"
    code, out = _uniop(db, tmp_path / "cache", ACCEPT)
    assert code == 2, out
    assert manifest_items(db)["image:uniop"]["message"] == APPTAINER_FAILED.format(line=GLIBC)
    assert not sif_path(tmp_path / "cache").exists()


# ---------------------------------------------------------------- uniop.def


def test_uniop_def_has_no_chown_and_fetches_nothing() -> None:
    text = UNIOP_DEF.read_text(encoding="utf-8")
    body = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert "chown" not in body
    for word in ("curl", "wget", "http://", "https://", "git clone"):
        assert word not in body, word
    assert re.search(r"^Bootstrap:\s*docker\s*$", body, re.MULTILINE)
    assert re.search(r"^From:\s*\S+@sha256:[0-9a-f]{64}\s*$", body, re.MULTILINE)
    assert "uniop-src.tar.gz" in body and "workflow/envs/uniop.yml" in body
    assert "chmod -R a+rX /opt/conda/share/uniop" in body


# ---------------------------------------------------------------- preflight rows (item 3)


def _reasons(db: Path, profile: str) -> dict[str, str]:
    from gippa.annotate.tools import TOOLS
    from gippa.preflight.environment import missing_tool_reasons

    return missing_tool_reasons([TOOLS["uniop"]], profile, db / "containers", None,
                                db / LOCAL_TOOLS)  # fmt: skip


def test_preflight_present_after_the_build(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, _ = _built(tmp_path, monkeypatch)
    assert _reasons(db, "apptainer") == {}


def test_preflight_missing_when_the_sif_bytes_differ(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, cache = _built(tmp_path, monkeypatch)
    sif = sif_path(cache)
    sif.write_bytes(sif.read_bytes() + b"x")
    assert "uniop" in _reasons(db, "apptainer")


def test_preflight_missing_when_the_sif_is_gone(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, cache = _built(tmp_path, monkeypatch)
    sif_path(cache).unlink()
    assert "uniop" in _reasons(db, "apptainer")


def test_preflight_stale_key(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db, _ = _built(tmp_path, monkeypatch)
    set_local_build(db, key="0123456789ab")
    reasons = _reasons(db, "apptainer")
    assert "uniop" in reasons
    want = f"local build 0123456789ab does not match this release ({uniop_key()})"
    assert want in reasons["uniop"], reasons


def test_example_19d_docker_only_build_is_missing_under_apptainer(
    archive: bytes, tmp_path: Path
) -> None:
    """A docker build never satisfies ``-profile apptainer``; nothing tries ``docker://localhost``."""
    db = tmp_path / "db"
    code, out = _uniop(db, tmp_path / "cache", ACCEPT, profile="docker")
    assert code == 0, out
    reasons = _reasons(db, "apptainer")
    assert "uniop" in reasons
    assert "docker://localhost" not in " ".join(reasons.values())
    assert [c for c in calls(tmp_path, "apptainer") if "docker://localhost" in " ".join(c)] == []
