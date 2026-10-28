"""M11-14u2: the ``image:uniop`` record across reruns, the docker preflight row, the run manifest.

docs/specs/M11-14-container-install.md item 1b ("What is recorded (14u2)": the ``local_build``
object and the run manifest's ``uniop`` entry), item 3 (rows "uniop with a local build", docker
column, and "not built, stale or wrong runtime"; "UniOP message"), examples 19, 19b, 19c; edge
case "UniOP: rerun after success"; #396 review S3 (a later ``setup --images`` without the terms
flag, or for a runtime that was not built, keeps a good record) and S5 (the image-ID check after
``docker image inspect``; ``0600`` on the archive copy in the build context; a corrupt archive in
place is deleted before re-fetching; ``local_build`` validates against the schema).

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- Preflight under docker: present iff the local tag exists **and** ``docker image inspect``'s
  ``.Id`` equals ``local_build.digest`` **and** ``local_build.key == key("uniop")``; the stale
  reason contains ``local build <old> does not match this release (<new>)``. The preflight
  message is ``uniop not installed; run: <pkg> setup --uniop --profile <P> --db-dir <db>
  --accept-uniop-terms`` with the run's profile and DB directory.
- A rerun with an ``ok`` item whose key and digest still match builds nothing (idempotent).
- A later UniOP step that does not build (terms not accepted, or no build for this runtime) keeps
  an ``ok`` ``image:uniop`` and its ``local_build`` as they are, while ``local_tools.config``
  names the built image, and exits 0 under ``--images``.
- Run manifest (``build_manifest``): the ``tools`` entry ``uniop`` gains ``local_build: true``,
  ``runtime`` (from ``local_build.runtime``), ``source_sha256`` (the archive checksum of the
  ``image:uniop`` item, ``sha256:<hex>`` or bare hex) and ``container_digest``: docker the image
  ID, apptainer ``sha256:<.sif bytes>`` (``local_build.digest``), micromamba and prebuilt null or
  absent. The schema accepts the entry.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest
from helpers.uniop_14u1 import (
    ACCEPT,
    ARCHIVE_PATH,
    ARCHIVE_URL,
    DB_PATH,
    FAKE_IMAGE_ID,
    PKG,
    SKIP_DBS,
    VERSION,
    FakeNet,
    Route,
    docker_calls,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
    uniop_key,
    uniop_tag,
)
from helpers.uniop_14u2 import (
    RUN_SCHEMA,
    install_fake_runtimes,
    set_local_build,
    validate_db_manifest,
)

pytestmark = pytest.mark.usefixtures("no_network")

LOCAL_TOOLS = "local_tools.config"
OTHER_ID = "sha256:" + "7a" * 32

FAKE_DOCKER_STAT = r"""#!/bin/sh
# records the mode of the staged archive at build time, then runs the 14u1 fake
if [ "$1" = build ]; then
    last=""; for a in "$@"; do last=$a; done
    # GNU stat -c (Linux), else BSD stat -f %Lp (macOS): the octal permission bits
    { stat -c %a "$last/uniop-src.tar.gz" 2>/dev/null ||
      stat -f %Lp "$last/uniop-src.tar.gz"; } >> "$FAKE_ARCHIVE_MODES"
fi
exec "$FAKE_DOCKER_REAL" "$@"
"""


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


def _uniop(db: Path, *extra: str, profile: str = "docker") -> tuple[int, str]:
    return setup_cli(["--db-dir", str(db), *SKIP_DBS, "--uniop", "--profile", profile, *extra])


def _built(tmp_path: Path) -> Path:
    db = tmp_path / "db"
    code, out = _uniop(db, ACCEPT)
    assert code == 0, out
    return db


def _reasons(db: Path, profile: str = "docker") -> dict[str, str]:
    from gippa.annotate.tools import TOOLS
    from gippa.preflight.environment import missing_tool_reasons

    return missing_tool_reasons([TOOLS["uniop"]], profile, db / "containers", None,
                                db / LOCAL_TOOLS)  # fmt: skip


# ---------------------------------------------------------------- S5: docker build gaps


def test_docker_local_build_validates_against_the_schema(archive: bytes, tmp_path: Path) -> None:
    data = validate_db_manifest(_built(tmp_path))
    (item,) = [i for i in data["items"] if i["name"] == "image:uniop"]
    assert set(item["local_build"]) == {"runtime", "ref", "digest", "key", "terms"}


@pytest.mark.parametrize("bad_id", ["", "not-an-id", "sha256:abc", "5e" * 32])
def test_image_id_from_inspect_must_be_a_sha256_digest(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad_id: str
) -> None:
    monkeypatch.setenv("FAKE_IMAGE_ID", bad_id)
    db = tmp_path / "db"
    code, out = _uniop(db, ACCEPT)
    assert code == 2, out
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "absent" and "docker image inspect" in item["message"], item
    assert "local_build" not in item
    assert not (db / LOCAL_TOOLS).exists() or uniop_tag() not in (db / LOCAL_TOOLS).read_text()


def test_build_without_an_inspectable_tag_is_absent(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The build claims success but the tag is not in the store: no record of a phantom image."""
    script = tmp_path / "fakebin" / "docker"
    text = script.read_text().replace('echo "$t" >> "$state/tags"', ": no tag recorded")
    script.write_text(text)
    db = tmp_path / "db"
    code, out = _uniop(db, ACCEPT)
    assert code == 2, out
    assert manifest_items(db)["image:uniop"]["status"] == "absent"


def test_archive_copy_in_the_build_context_is_0600(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bindir = tmp_path / "fakebin"
    real = bindir / "docker.real"
    (bindir / "docker").rename(real)
    (bindir / "docker").write_text(FAKE_DOCKER_STAT)
    (bindir / "docker").chmod(0o755)
    monkeypatch.setenv("FAKE_DOCKER_REAL", str(real))
    modes = tmp_path / "modes.txt"
    monkeypatch.setenv("FAKE_ARCHIVE_MODES", str(modes))
    _built(tmp_path)
    assert modes.read_text().split() == ["600"]


def test_corrupt_archive_in_place_is_deleted_before_refetching(
    net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A corrupt copy is removed even when the re-fetch then fails (offline)."""
    pin_uniop_source(monkeypatch, sha256_hex(fixture_archive()))
    net.routes[("github.com", ARCHIVE_PATH)] = Route(error="[Errno 101] Network is unreachable")
    db = tmp_path / "db"
    (db / "uniop").mkdir(parents=True)
    (db / DB_PATH).write_bytes(b"corrupt bytes")
    code, out = _uniop(db, ACCEPT)
    assert code == 2, out
    assert not (db / DB_PATH).exists(), "the corrupt archive must be deleted"
    assert net.contacted, "a re-fetch was attempted"


def test_corrupt_archive_in_place_is_replaced(archive: bytes, net: FakeNet, tmp_path: Path) -> None:
    db = tmp_path / "db"
    (db / "uniop").mkdir(parents=True)
    (db / DB_PATH).write_bytes(b"corrupt bytes")
    code, out = _uniop(db, ACCEPT)
    assert code == 0, out
    assert (db / DB_PATH).read_bytes() == archive
    assert len(docker_calls(tmp_path, "build")) == 1


# ---------------------------------------------------------------- idempotent rerun


def test_rerun_after_success_builds_nothing(archive: bytes, net: FakeNet, tmp_path: Path) -> None:
    db = _built(tmp_path)
    before = manifest_items(db)["image:uniop"]
    contacted = len(net.contacted)
    code, out = _uniop(db, ACCEPT)
    assert code == 0, out
    assert len(docker_calls(tmp_path, "build")) == 1, "no second build"
    assert len(net.contacted) == contacted, "no second download"
    assert manifest_items(db)["image:uniop"]["local_build"] == before["local_build"]


def test_rerun_rebuilds_when_the_image_id_changed(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = _built(tmp_path)
    monkeypatch.setenv("FAKE_IMAGE_ID", OTHER_ID)  # the tag now names another image
    code, out = _uniop(db, ACCEPT)
    assert code == 0, out
    assert len(docker_calls(tmp_path, "build")) == 2
    assert manifest_items(db)["image:uniop"]["local_build"]["digest"] == OTHER_ID


# ---------------------------------------------------------------- S3: a good record is kept


@pytest.mark.parametrize("profile", ["docker", "apptainer"])
def test_later_images_without_the_flag_keeps_the_good_record(
    archive: bytes, net: FakeNet, tmp_path: Path, profile: str
) -> None:
    db = _built(tmp_path)
    before = manifest_items(db)["image:uniop"]
    lt_before = (db / LOCAL_TOOLS).read_text()
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", profile])
    assert code == 0, out
    after = manifest_items(db)["image:uniop"]
    assert after["status"] == "ok", after
    assert after["local_build"] == before["local_build"]
    assert (db / LOCAL_TOOLS).read_text() == lt_before
    assert len(docker_calls(tmp_path, "build")) == 1


def test_later_uniop_for_an_unbuilt_runtime_without_the_flag_keeps_the_record(
    archive: bytes, tmp_path: Path
) -> None:
    db = _built(tmp_path)
    before = manifest_items(db)["image:uniop"]
    code, out = _uniop(db, profile="apptainer")
    assert code == 2, out  # explicit --uniop without acceptance still exits 2
    after = manifest_items(db)["image:uniop"]
    assert after["status"] == "ok" and after["local_build"] == before["local_build"], after


# ---------------------------------------------------------------- preflight rows (docker)


def test_example_19_present_after_example_16(archive: bytes, tmp_path: Path) -> None:
    assert _reasons(_built(tmp_path)) == {}


def test_example_19c_image_id_differs_from_the_record(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = _built(tmp_path)
    monkeypatch.setenv("FAKE_IMAGE_ID", OTHER_ID)
    assert "uniop" in _reasons(db)


def test_example_19c_local_tag_removed(archive: bytes, tmp_path: Path) -> None:
    db = _built(tmp_path)
    (tmp_path / "docker_state" / "tags").write_text("")
    assert "uniop" in _reasons(db)


def test_example_19c_stale_key(archive: bytes, tmp_path: Path) -> None:
    db = _built(tmp_path)
    set_local_build(db, key="0123456789ab")
    reason = _reasons(db).get("uniop", "")
    assert f"local build 0123456789ab does not match this release ({uniop_key()})" in reason


def test_missing_local_build_counts_as_not_built(archive: bytes, tmp_path: Path) -> None:
    """A docker block without a recorded ``local_build`` (e.g. a failed or foreign record) is
    not a usable build."""
    db = _built(tmp_path)
    path = db / "db_manifest.json"
    data = json.loads(path.read_text())
    data["items"] = [i for i in data["items"] if i["name"] != "image:uniop"]
    path.write_text(json.dumps(data))
    assert "uniop" in _reasons(db)


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


def _preflight_tools(tmp_path: Path, db: Path, profile: str,
                     monkeypatch: pytest.MonkeyPatch) -> Any:  # fmt: skip
    from gippa.config import resolve_config
    from gippa.preflight import environment
    from gippa.preflight.run import PreflightInputs, run_preflight

    monkeypatch.setattr(environment, "probe_uniop", lambda *a, **k: (0, ""))
    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(exist_ok=True)
    g.write_text(">chr1 [topology=circular]\n" + "ACGT" * 30_000 + "\n")
    p = tmp_path / "in" / "p.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    cfg = resolve_config(overrides=["preflight.min_n50_warn=1",
                                    "preflight.min_largest_contig_warn=1",
                                    "annotation.uniop.enabled=true"])  # fmt: skip
    report, _ = run_preflight(PreflightInputs(genome=g, proteins=p), cfg, None, "run",
                              tmp_path / "OUT", db, profile=profile)  # fmt: skip
    (tools,) = [c for c in report.checks if c.rule_id == "environment.tools"]
    return tools


def _missing(check: Any) -> set[str]:
    value = check.value if isinstance(check.value, dict) else {}
    return set(value.get("missing", []))


def test_example_19_preflight_uniop_present(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = _built(tmp_path)
    assert "uniop" not in _missing(_preflight_tools(tmp_path, db, "docker", monkeypatch))


def test_example_19b_nothing_built_names_the_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tools = _preflight_tools(tmp_path, Path("/data/db"), "apptainer", monkeypatch)
    want = (f"uniop not installed; run: {PKG} setup --uniop --profile apptainer --db-dir /data/db "
            "--accept-uniop-terms")  # fmt: skip
    assert "uniop" in _missing(tools) and want in tools.message, tools.message


def test_example_19c_stale_preflight_message(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = _built(tmp_path)
    set_local_build(db, key="0123456789ab")
    tools = _preflight_tools(tmp_path, db, "docker", monkeypatch)
    assert "uniop" in _missing(tools)
    assert f"local build 0123456789ab does not match this release ({uniop_key()})" in tools.message
    want = (f"run: {PKG} setup --uniop --profile docker --db-dir {db.absolute()} "
            "--accept-uniop-terms")  # fmt: skip
    assert want in tools.message, tools.message


# ---------------------------------------------------------------- run manifest uniop entry

ARCHIVE_SHA = "ab" * 32
SIF_DIGEST = "sha256:" + "c4" * 32


def _db_manifest(runtime: str, digest: str) -> dict[str, Any]:
    from helpers.fake_env import setup_manifest

    data = setup_manifest()
    items = [i for i in data["items"] if i["name"] != "image:uniop"]  # type: ignore[attr-defined]
    items.append({"name": "image:uniop", "version": VERSION, "url": ARCHIVE_URL,
                  "checksum": f"sha256:{ARCHIVE_SHA}", "size_bytes": 0,
                  "retrieved": "2026-10-26", "path": "", "status": "ok", "message": "",
                  "local_build": {"runtime": runtime, "ref": "r", "digest": digest,
                                  "key": "0123456789ab",
                                  "terms": "" if runtime == "prebuilt" else "flag"}})  # fmt: skip
    data["items"] = items
    return data


@pytest.mark.parametrize(
    ("runtime", "digest", "want"),
    [
        ("docker", FAKE_IMAGE_ID, FAKE_IMAGE_ID),
        ("apptainer", SIF_DIGEST, SIF_DIGEST),
        ("micromamba", "", None),
        ("prebuilt", "", None),
    ],
)
def test_run_manifest_uniop_entry(runtime: str, digest: str, want: str | None) -> None:
    from helpers.manifest_data import build, validate

    report = {"process": "ANNOTATE:UNIOP", "name": "uniop", "version": "3ed14dfc",
              "container_digest": None, "conda_lock_hash": None}  # fmt: skip
    m = build(tool_reports=[report], db_manifest=_db_manifest(runtime, digest))
    (entry,) = [t for t in m["tools"] if t["name"] == "uniop"]
    assert entry["local_build"] is True, entry
    assert entry["runtime"] == runtime
    assert entry["source_sha256"] in (ARCHIVE_SHA, f"sha256:{ARCHIVE_SHA}"), entry
    assert entry.get("container_digest") == want, entry
    assert "version" in entry
    validate(m)


def test_run_manifest_other_tools_have_no_local_build() -> None:
    from helpers.manifest_data import build

    reports = [
        {"process": "ANNOTATE:UNIOP", "name": "uniop", "version": "3ed14dfc",
         "container_digest": None, "conda_lock_hash": None},
        {"process": "ANNOTATE:BAKTA", "name": "bakta", "version": "1.11.4",
         "container_digest": FAKE_IMAGE_ID, "conda_lock_hash": None},
    ]  # fmt: skip
    m = build(tool_reports=reports, db_manifest=_db_manifest("docker", FAKE_IMAGE_ID))
    bakta = next(t for t in m["tools"] if t["name"] == "bakta")
    assert "local_build" not in bakta and "runtime" not in bakta


def test_run_manifest_schema_accepts_the_local_build_keys() -> None:
    schema = json.loads(RUN_SCHEMA.read_text(encoding="utf-8"))
    props = schema["properties"]["tools"]["items"]["properties"]
    assert {"local_build", "runtime", "source_sha256"} <= set(props), sorted(props)
    assert props["local_build"].get("type") == "boolean" or props["local_build"].get("const")
