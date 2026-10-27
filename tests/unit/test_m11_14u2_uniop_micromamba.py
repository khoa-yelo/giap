"""M11-14u2: the micromamba build of UniOP (plain ``conda``), the ``farmshare*`` prebuilt check.

docs/specs/M11-14-container-install.md item 1b ("micromamba (14u2), plain conda profile",
"micromamba, farmshare* lists (the prebuilt case)", "local_tools.config", "What is recorded
(14u2)"), item 3 (row "uniop with a local build", ``conda`` column), examples 16c, 16d; AC 9 (14u2:
"a micromamba-only build under ``-profile docker`` is missing (no docker block), with the docker
command"); #396 review S5 (``local_build`` validates against the schema).

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- Plain ``conda``: one ``micromamba create`` with ``-p``/``--prefix``
  ``<cache-dir>/envs/gippa-uniop`` and ``-f``/``--file`` naming a file whose bytes are
  ``workflow/envs/locks/uniop.explicit.txt``;
  then the verified archive is extracted (first path component stripped) into
  ``<prefix>/share/uniop`` with ``GIPPA_COMMIT`` = the pinned commit. Nothing is fetched by the
  builder itself; docker and apptainer are never called.
- ``local_tools.config`` gains ``conda { process { withName: 'UNIOP' { conda = '<prefix>' } } }``
  and never ``conda.enabled``.
- ``image:uniop``: ``path`` = the prefix; ``local_build.runtime`` ``micromamba``, ``ref`` = the
  prefix, ``key`` = key("uniop"), ``terms`` ``flag``.
- ``farmshare*``: no fetch, no build, ``local_tools.config`` untouched; ``local_build`` =
  ``{"runtime": "prebuilt", "ref": "$GIPPA_ROOT/micromamba/envs/gippa-uniop", "key": key("uniop"),
  "terms": ""}``; without ``share/uniop`` exit 2 naming ``sbatch scripts/slurm/build_envs.sbatch
  uniop``.
- Preflight under plain ``conda`` reads the ``conda`` block's prefix (no ``$GIPPA_ROOT`` needed):
  present iff ``<prefix>/share/uniop/GIPPA_COMMIT`` equals the pin and the key is current.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.uniop_14u1 import (
    ACCEPT,
    COMMIT,
    PKG,
    SKIP_DBS,
    FakeNet,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
    uniop_key,
)
from helpers.uniop_14u2 import (
    PREBUILT_HINT,
    UNIOP_LOCK,
    calls,
    install_fake_runtimes,
    prebuilt_root,
    set_local_build,
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


def _uniop(db: Path, cache: Path, *extra: str, profile: str = "conda") -> tuple[int, str]:
    return setup_cli(["--db-dir", str(db), "--cache-dir", str(cache), *SKIP_DBS, "--uniop",
                      "--profile", profile, *extra])  # fmt: skip


def _built(tmp_path: Path) -> tuple[Path, Path]:
    db, cache = tmp_path / "db", tmp_path / "cache"
    code, out = _uniop(db, cache, ACCEPT)
    assert code == 0, out
    return db, cache / "envs" / "gippa-uniop"


def _flag(argv: list[str], *names: str) -> str | None:
    for i, a in enumerate(argv):
        if a in names and i + 1 < len(argv):
            return argv[i + 1]
        for n in names:
            if n.startswith("--") and a.startswith(n + "="):
                return a.split("=", 1)[1]
    return None


# ---------------------------------------------------------------- example 16c


def test_example_16c_one_create_from_the_lock(archive: bytes, tmp_path: Path) -> None:
    _, prefix = _built(tmp_path)
    creates = calls(tmp_path, "micromamba", "create")
    assert len(creates) == 1, creates
    assert Path(_flag(creates[0], "-p", "--prefix") or "") == prefix
    seen = tmp_path / "micromamba_state" / "lock_seen"
    assert seen.read_bytes() == UNIOP_LOCK.read_bytes(), "the env is created from the lock"


def test_example_16c_archive_extracted_with_the_commit(archive: bytes, tmp_path: Path) -> None:
    _, prefix = _built(tmp_path)
    home = prefix / "share" / "uniop"
    assert (home / "GIPPA_COMMIT").read_text(encoding="utf-8").strip() == COMMIT
    assert (home / "src" / "UniOP").is_file(), "first path component stripped (UniOP-<commit>/)"
    assert (home / "README.fixture").read_bytes() == b"not UniOP code\n"


def test_example_16c_conda_block(archive: bytes, tmp_path: Path) -> None:
    db, prefix = _built(tmp_path)
    text = (db / LOCAL_TOOLS).read_text(encoding="utf-8")
    assert "conda.enabled" not in text
    m = re.search(r"\bconda\s*\{\s*process\s*\{\s*withName\s*:\s*'UNIOP'\s*\{\s*conda\s*=\s*"
                  r"'([^']+)'\s*\}\s*\}\s*\}", text)  # fmt: skip
    assert m is not None, text
    assert Path(m.group(1)) == prefix and Path(m.group(1)).is_absolute()
    assert not re.search(r"\b(docker|apptainer)\s*\{", text), text


def test_example_16c_record(archive: bytes, tmp_path: Path) -> None:
    db, prefix = _built(tmp_path)
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "ok", item
    assert item["path"] == str(prefix)
    lb = item["local_build"]
    assert (lb["runtime"], lb["ref"], lb["key"], lb["terms"]) == (
        "micromamba", str(prefix), uniop_key(), "flag",
    )  # fmt: skip
    validate_db_manifest(db)


def test_micromamba_never_calls_a_container_runtime(archive: bytes, tmp_path: Path) -> None:
    _built(tmp_path)
    assert calls(tmp_path, "docker") == [] and calls(tmp_path, "apptainer") == []


def test_micromamba_failure_is_recorded_and_exits_2(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_MAMBA_FAIL", "critical libmamba fixture solve boom")
    db = tmp_path / "db"
    code, out = _uniop(db, tmp_path / "cache", ACCEPT)
    assert code == 2, out
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "absent" and "fixture solve boom" in item["message"], item
    assert not (db / LOCAL_TOOLS).exists() or "conda" not in (db / LOCAL_TOOLS).read_text()


# ---------------------------------------------------------------- preflight rows (item 3)


def _reasons(db: Path, profile: str) -> dict[str, str]:
    from gippa.annotate.tools import TOOLS
    from gippa.preflight.environment import missing_tool_reasons

    return missing_tool_reasons([TOOLS["uniop"]], profile, db / "containers", None,
                                db / LOCAL_TOOLS)  # fmt: skip


def test_preflight_plain_conda_present_without_gippa_root(archive: bytes, tmp_path: Path) -> None:
    db, _ = _built(tmp_path)
    assert _reasons(db, "conda") == {}


def test_preflight_plain_conda_missing_on_a_wrong_commit(archive: bytes, tmp_path: Path) -> None:
    db, prefix = _built(tmp_path)
    (prefix / "share" / "uniop" / "GIPPA_COMMIT").write_text("0" * 40 + "\n")
    assert "uniop" in _reasons(db, "conda")


def test_preflight_plain_conda_stale_key(archive: bytes, tmp_path: Path) -> None:
    db, _ = _built(tmp_path)
    set_local_build(db, key="0123456789ab")
    reasons = _reasons(db, "conda")
    want = f"local build 0123456789ab does not match this release ({uniop_key()})"
    assert want in reasons.get("uniop", ""), reasons


def test_micromamba_only_build_is_missing_under_docker(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 9 (14u2): no docker block, so ``-profile docker`` is missing, naming the docker
    command."""
    from gippa.config import resolve_config
    from gippa.preflight import environment
    from gippa.preflight.run import PreflightInputs, run_preflight

    db, _ = _built(tmp_path)
    assert "uniop" in _reasons(db, "docker")
    monkeypatch.setattr(environment, "probe_uniop", lambda *a, **k: (0, ""))
    g = tmp_path / "g.fa"
    g.write_text(">chr1 [topology=circular]\n" + "ACGT" * 30_000 + "\n")
    p = tmp_path / "p.faa"
    p.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA\n")
    cfg = resolve_config(overrides=["preflight.min_n50_warn=1",
                                    "preflight.min_largest_contig_warn=1",
                                    "annotation.uniop.enabled=true"])  # fmt: skip
    report, _ = run_preflight(PreflightInputs(genome=g, proteins=p), cfg, None, "run",
                              tmp_path / "OUT", db, profile="docker")  # fmt: skip
    (tools,) = [c for c in report.checks if c.rule_id == "environment.tools"]
    want = (f"uniop not installed; run: {PKG} setup --uniop --profile docker "
            f"--db-dir {db.absolute()} --accept-uniop-terms")  # fmt: skip
    assert want in tools.message, tools.message


# ---------------------------------------------------------------- example 16d: prebuilt


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "farmshare", "farmshare_local"])
def test_example_16d_prebuilt_builds_nothing(
    net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, profile: str
) -> None:
    root = prebuilt_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    db = tmp_path / "db"
    db.mkdir()
    before = (
        "profiles {\n    docker { process { withName: 'UNIOP' { container = 'x'; "
        "conda = '' } } }\n}\n"
    )
    (db / LOCAL_TOOLS).write_text(before)
    code, out = _uniop(db, tmp_path / "cache", ACCEPT, profile=profile)
    assert code == 0, out
    assert net.contacted == [], "no fetch for the prebuilt env"
    for prog in ("docker", "apptainer", "micromamba"):
        assert calls(tmp_path, prog) == [], prog
    assert (db / LOCAL_TOOLS).read_text() == before, "local_tools.config unchanged"
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "ok", item
    prefix = root / "micromamba" / "envs" / "gippa-uniop"
    lb = item["local_build"]
    assert (lb["runtime"], lb["ref"], lb["terms"]) == ("prebuilt", str(prefix), ""), lb
    assert lb["key"] == uniop_key()
    validate_db_manifest(db)


def test_example_16d_prebuilt_writes_no_local_tools(
    net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIPPA_ROOT", str(prebuilt_root(tmp_path / "gippa_root")))
    db = tmp_path / "db"
    code, out = _uniop(db, tmp_path / "cache", ACCEPT, profile="farmshare,apptainer")
    assert code == 0, out
    assert not (db / LOCAL_TOOLS).exists()


@pytest.mark.parametrize("share", [False, None], ids=["no-share-uniop", "no-prefix"])
def test_example_16d_without_the_prefix_exits_2(
    net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, share: bool | None
) -> None:
    root = tmp_path / "gippa_root"
    if share is None:
        root.mkdir()
    else:
        prebuilt_root(root, share=share)
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    db = tmp_path / "db"
    code, out = _uniop(db, tmp_path / "cache", ACCEPT, profile="farmshare,apptainer")
    assert code == 2, out
    assert PREBUILT_HINT in out, out
    assert net.contacted == []
    assert calls(tmp_path, "micromamba") == [] and calls(tmp_path, "apptainer") == []
    item = manifest_items(db).get("image:uniop")
    assert item is None or item["status"] != "ok"
