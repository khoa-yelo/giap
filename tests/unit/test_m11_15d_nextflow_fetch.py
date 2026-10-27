"""M11-15d: ``gippa setup`` fetches the pinned Nextflow ``-dist`` launcher (spec item 5).

docs/specs/M11-15-portable-clusters.md item 5 (pip path, pin, allowlist, Java in setup),
examples 5.1, 5.2, 5.5, 5.6; AC "setup test with a fake fetcher for 5.1/5.2/5.6 (no network);
per-item redirect test (the UniOP row still refuses the asset host; the ``nextflow`` row refuses
``codeload.github.com``); ``_load_sources`` still rejects a path outside the DB directory".

- The ``nextflow`` row of ``db_sources.yaml``: the official ``nextflow-26.10.0-dist`` asset of the
  ``nextflow-io/nextflow`` GitHub release ``v26.10.0`` (R6), ``path: bin/nextflow``, a sha256 pin,
  and ``redirect_hosts: [<the one release-asset host>]``.
- The redirect allowlist is per item: a row's ``redirect_hosts`` (``DbItem.redirect_hosts``)
  replaces the ``REDIRECT_HOSTS`` default for that row's URL; ``REDIRECT_HOSTS`` is unchanged, so
  the UniOP archive (same host ``github.com``) still redirects only to ``codeload.github.com``.
  ``HttpsFetcher.fetch(url, dest)`` keeps its signature and finds the row by its URL (as
  ``_body_cap`` does).
- ``gippa setup`` installs ``<db_dir>/bin/nextflow`` (mode 0755) when ``nextflow`` is not on
  ``PATH``; ``--no-nextflow`` skips it; ``SetupOptions.nextflow`` defaults to ``True``.
- Java: a warning only in setup.

The network is a fake ``http.client`` (``helpers.uniop_14u1.FakeNet``) behind the real
``HttpsFetcher``, sockets and DNS blocked (``no_network``). The pin is replaced by the sha256 of a
few fake bytes (the real pin is computed by the implementer from the official asset). ``java``
and ``nextflow`` are fake scripts; ``PATH`` drops every host directory holding a ``nextflow``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import os
import re
import stat
import sys
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest
import yaml
from helpers.uniop_14u1 import (
    ARCHIVE_PATH,
    ARCHIVE_URL,
    CODELOAD_PATH,
    FakeNet,
    Route,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    setup_cli,
)

pytestmark = [pytest.mark.usefixtures("no_network"), pytest.mark.real_prereqs]

REPO = Path(__file__).resolve().parents[2]
SOURCES_YAML = REPO / "src" / "gippa" / "data" / "db_sources.yaml"
NF_VERSION = "26.10.0"
NF_URL = (
    f"https://github.com/nextflow-io/nextflow/releases/download/v{NF_VERSION}/"
    f"nextflow-{NF_VERSION}-dist"
)
NF_PATH = f"/nextflow-io/nextflow/releases/download/v{NF_VERSION}/nextflow-{NF_VERSION}-dist"
ASSET_HOSTS = ("release-assets.githubusercontent.com", "objects.githubusercontent.com")
ASSET_PATH = "/github-production-release-asset/43446421/0f5e7c1a-26-10-0-dist"
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models"]
FAKE_DIST = b"#!/bin/sh\n# fake nextflow-26.10.0-dist (test bytes, not Nextflow)\nexit 0\n"
JAVA_21 = 'openjdk version "21.0.1" 2023-10-17\nOpenJDK Runtime Environment (build 21.0.1+12)\n'
JAVA_11 = 'openjdk version "11.0.2" 2019-01-15\nOpenJDK Runtime Environment 18.9\n'
JAVA_HEAD = "Java 17 or newer is needed for Nextflow"


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def raw_rows() -> dict[str, Any]:
    return yaml.safe_load(SOURCES_YAML.read_text(encoding="utf-8"))


def nf_row() -> dict[str, Any]:
    rows = raw_rows()
    assert "nextflow" in rows, "db_sources.yaml has no `nextflow` row (M11-15d item 5)"
    return rows["nextflow"]


def asset_host() -> str:
    hosts = nf_row().get("redirect_hosts")
    assert isinstance(hosts, list) and len(hosts) == 1, (
        f"the nextflow row pins exactly one release-asset host, got {hosts!r}"
    )
    return str(hosts[0])


def other_asset_host() -> str:
    return next(h for h in ASSET_HOSTS if h != asset_host())


def pin(monkeypatch: pytest.MonkeyPatch, **sha_hex: str) -> Any:
    """Replace the checksum of each named row by ``sha256:<hex>`` in every loaded ``gippa``
    module that holds the original ``SOURCES`` mapping."""
    from gippa.setupdb import items

    original = items.SOURCES
    for name in sha_hex:
        assert name in original, f"SOURCES has no {name!r} row"
    patched = MappingProxyType(
        {
            **original,
            **{
                n: dataclasses.replace(original[n], checksum=f"sha256:{h}")
                for n, h in sha_hex.items()
            },
        }
    )
    for mod_name, mod in list(sys.modules.items()):
        if (mod_name == "gippa" or mod_name.startswith("gippa.")) and getattr(
            mod, "SOURCES", None
        ) is original:
            monkeypatch.setattr(mod, "SOURCES", patched)
    return patched


def make_exe(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)
    return path


def fake_java_script(stderr: str, code: int = 0) -> str:
    lines = "".join(f"printf '%s\\n' '{ln}' >&2\n" for ln in stderr.splitlines())
    return f"#!/bin/sh\n{lines}exit {code}\n"


FAKE_NF_ON_PATH = """#!/bin/sh
case "$1" in
-v) echo "nextflow version 25.10.0.10289";;
*) printf '\\n      N E X T F L O W\\n      version 25.10.0 build 10289\\n'
   printf '      created 22-10-2025 15:38 UTC\\n      cite doi:10.1038/nbt.3820\\n';;
esac
exit 0
"""


def host_path_without(*names: str) -> list[str]:
    """The host ``PATH`` minus every directory that holds one of ``names``."""
    out = []
    for d in os.environ.get("PATH", "").split(os.pathsep):
        if d and not any((Path(d) / n).exists() for n in names):
            out.append(d)
    return out


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    for var in ("JAVA_HOME", "JAVA_CMD", "NXF_HOME", "NXF_OFFLINE", "NXF_VER"):
        monkeypatch.delenv(var, raising=False)
    old = os.umask(0o022)
    yield
    os.umask(old)


@pytest.fixture
def fakebin(tmp_path: Path) -> Path:
    d = tmp_path / "fakebin"
    d.mkdir()
    return d


def use_path(monkeypatch: pytest.MonkeyPatch, fakebin: Path, *, host: bool = True) -> None:
    """``fakebin`` first; the host ``PATH`` without any ``nextflow`` (and, with ``host=False``,
    nothing else at all: no real ``java`` either)."""
    rest = host_path_without("nextflow") if host else []
    monkeypatch.setenv("PATH", os.pathsep.join([str(fakebin), *rest]))


@pytest.fixture
def net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


def serve_nextflow(net: FakeNet, body: bytes = FAKE_DIST) -> None:
    """GitHub's 302 from the release download URL to the pinned release-asset host."""
    loc = f"https://{asset_host()}{ASSET_PATH}?sp=r&sv=2018-11-09&sr=b&sig=x"
    net.routes[("github.com", NF_PATH)] = Route(302, {"Location": loc})
    net.routes[(asset_host(), ASSET_PATH)] = Route(
        200, {"Content-Type": "application/octet-stream"}, body
    )


def nf_args(db: Path, *extra: str) -> list[str]:
    return ["--db-dir", str(db), *SKIP_DBS, *extra]


def no_leftovers(db: Path) -> list[str]:
    return [str(p) for p in db.rglob("*") if p.name.endswith((".partial", ".part"))]


# ---------------------------------------------------------------- the db_sources.yaml row


def test_nextflow_row_is_the_official_dist_asset() -> None:
    row = nf_row()
    assert row["url"] == NF_URL, "R6: the nextflow-io/nextflow GitHub release v26.10.0 -dist"
    assert str(row["version"]) == NF_VERSION
    assert row["path"] == "bin/nextflow"
    assert row.get("tool") is None, "fetched over HTTPS by gippa itself"
    assert re.fullmatch(r"sha256:[0-9a-f]{64}", str(row["checksum"])), row["checksum"]
    assert asset_host() in ASSET_HOSTS


def test_nextflow_row_loads_through_load_sources() -> None:
    from gippa.setupdb import SOURCES

    assert "nextflow" in SOURCES, "SOURCES has no nextflow item"
    item = SOURCES["nextflow"]
    assert (item.url, item.host, item.path, item.version) == (
        NF_URL,
        "github.com",
        "bin/nextflow",
        NF_VERSION,
    )
    assert item.checksum == nf_row()["checksum"]
    assert set(item.redirect_hosts) == {asset_host()}


def test_no_other_row_names_the_asset_host() -> None:
    """The per-item allowlist widens only the nextflow row (UniOP keeps codeload only)."""
    rows = raw_rows()
    for name, row in rows.items():
        if name != "nextflow":
            assert asset_host() not in (row.get("redirect_hosts") or []), name
    uniop = rows["uniop_source"].get("redirect_hosts")
    assert uniop in (None, ["codeload.github.com"]), uniop


def test_redirect_hosts_default_is_unchanged() -> None:
    from gippa.setupdb.run import REDIRECT_HOSTS

    assert REDIRECT_HOSTS["github.com"] == frozenset({"codeload.github.com"})
    assert all(asset_host() not in hosts for hosts in REDIRECT_HOSTS.values())


def _load_from(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, rows: dict[str, Any]) -> Any:
    from gippa.setupdb import items

    data = tmp_path / "data"
    data.mkdir(exist_ok=True)
    (data / "db_sources.yaml").write_text(yaml.safe_dump(rows), encoding="utf-8")
    monkeypatch.setattr(items, "files", lambda _pkg: data)
    return items._load_sources()


def _row(path: str, **extra: Any) -> dict[str, Any]:
    return {"version": NF_VERSION, "url": NF_URL, "checksum": "sha256:" + "0" * 64,
            "path": path, "tool": None, "check": "sha256", **extra}  # fmt: skip


def test_load_sources_reads_redirect_hosts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    got = _load_from(tmp_path, monkeypatch, {
        "nextflow": _row("bin/nextflow", redirect_hosts=["release-assets.githubusercontent.com"]),
        "plain": _row("x/y"),
    })  # fmt: skip
    assert set(got["nextflow"].redirect_hosts) == {"release-assets.githubusercontent.com"}
    assert not got["plain"].redirect_hosts, "a row without redirect_hosts uses the default"


@pytest.mark.parametrize(
    "path", ["../bin/nextflow", "/usr/local/bin/nextflow", "bin/../../nextflow", ""]
)
def test_load_sources_still_rejects_a_path_outside_the_db_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, path: str
) -> None:
    safe = _load_from(
        tmp_path, monkeypatch, {"nextflow": _row("bin/nextflow", redirect_hosts=[ASSET_HOSTS[0]])}
    )
    assert set(safe["nextflow"].redirect_hosts) == {ASSET_HOSTS[0]}, "the row itself is valid"
    rows = {"nextflow": _row(path, redirect_hosts=[ASSET_HOSTS[0]])}
    with pytest.raises(ValueError, match="unsafe path"):
        _load_from(tmp_path, monkeypatch, rows)


# ---------------------------------------------------------------- per-item redirects (AC)


def test_nextflow_redirect_to_the_asset_host_is_followed(net: FakeNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher

    serve_nextflow(net)
    HttpsFetcher().fetch(NF_URL, tmp_path / "nf")
    assert (tmp_path / "nf").read_bytes() == FAKE_DIST
    assert net.contacted == [("github.com", NF_PATH), (asset_host(), ASSET_PATH)]


@pytest.mark.parametrize(
    "location",
    [
        "https://codeload.github.com" + CODELOAD_PATH,
        "OTHER_ASSET_HOST",
        "http://ASSET_HOST" + ASSET_PATH,
        "https://ASSET_HOST.evil.example" + ASSET_PATH,
        "https://evil.example/nextflow",
    ],
)
def test_nextflow_row_refuses_every_other_redirect(
    net: FakeNet, tmp_path: Path, location: str
) -> None:
    from gippa.setupdb import HttpsFetcher, SourceError

    if location == "OTHER_ASSET_HOST":
        location = f"https://{other_asset_host()}{ASSET_PATH}"
    location = location.replace("ASSET_HOST", asset_host())
    net.routes[("github.com", NF_PATH)] = Route(302, {"Location": location})
    with pytest.raises(SourceError, match="refused"):
        HttpsFetcher().fetch(NF_URL, tmp_path / "nf")
    assert net.contacted == [("github.com", NF_PATH)], "the new location is never contacted"
    assert not (tmp_path / "nf").exists()


def test_uniop_row_still_refuses_the_asset_host(net: FakeNet, tmp_path: Path) -> None:
    from gippa.setupdb import HttpsFetcher, SourceError

    location = f"https://{asset_host()}{ASSET_PATH}"
    net.routes[("github.com", ARCHIVE_PATH)] = Route(302, {"Location": location})
    with pytest.raises(SourceError, match="refused"):
        HttpsFetcher().fetch(ARCHIVE_URL, tmp_path / "a.tar.gz")
    assert net.contacted == [("github.com", ARCHIVE_PATH)]


def test_uniop_row_still_follows_codeload(net: FakeNet, tmp_path: Path) -> None:
    """Regression guard next to the new row (same host ``github.com``)."""
    from gippa.setupdb import HttpsFetcher

    nf_row()  # fails first while the nextflow row is missing
    data = fixture_archive()
    net.serve_github(data)
    HttpsFetcher().fetch(ARCHIVE_URL, tmp_path / "a.tar.gz")
    assert (tmp_path / "a.tar.gz").read_bytes() == data


# ---------------------------------------------------------------- options


def test_setup_options_fetch_nextflow_by_default(tmp_path: Path) -> None:
    from gippa.setupdb import SetupOptions

    opts = SetupOptions(tmp_path / "db")
    assert opts.nextflow is True


# ---------------------------------------------------------------- example 5.1


@pytest.mark.parametrize("umask", [0o022, 0o077])
def test_example_5_1_downloads_checks_installs_and_records(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, umask: int
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    pin(monkeypatch, nextflow=sha256_hex(FAKE_DIST))
    serve_nextflow(net)
    db = tmp_path / "db"
    os.umask(umask)
    code, out = setup_cli(nf_args(db))
    assert code == 0, out
    assert net.contacted == [("github.com", NF_PATH), (asset_host(), ASSET_PATH)]
    exe = db / "bin" / "nextflow"
    assert exe.read_bytes() == FAKE_DIST
    assert stat.S_IMODE(exe.stat().st_mode) == 0o755
    rec = manifest_items(db)["nextflow"]
    assert (rec["status"], rec["path"], rec["url"], rec["version"]) == (
        "ok",
        "bin/nextflow",
        NF_URL,
        NF_VERSION,
    )
    assert rec["checksum"] == f"sha256:{sha256_hex(FAKE_DIST)}"
    assert rec["size_bytes"] == len(FAKE_DIST)
    assert not no_leftovers(db)
    assert JAVA_HEAD not in out, "Java 21 is new enough: no warning"


def test_rerun_adopts_the_verified_launcher_without_downloading(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    pin(monkeypatch, nextflow=sha256_hex(FAKE_DIST))
    db = tmp_path / "db"
    (db / "bin").mkdir(parents=True)
    make_exe(db / "bin" / "nextflow", FAKE_DIST.decode())
    code, out = setup_cli(nf_args(db))
    assert code == 0, out
    assert net.contacted == []
    assert manifest_items(db)["nextflow"]["status"] == "ok"


# ---------------------------------------------------------------- example 5.2


def test_example_5_2_checksum_mismatch_installs_nothing(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    pinned = sha256_hex(b"the official bytes")
    pin(monkeypatch, nextflow=pinned)
    serve_nextflow(net, FAKE_DIST)
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db))
    assert code == 2, out
    assert not (db / "bin" / "nextflow").exists(), "nothing installed"
    assert not no_leftovers(db)
    rec = manifest_items(db)["nextflow"]
    assert rec["status"] != "ok"
    assert sha256_hex(FAKE_DIST) in rec["message"] and pinned in rec["message"], rec["message"]
    assert "pinned" in rec["message"]
    assert sha256_hex(FAKE_DIST) in out and pinned in out


# ---------------------------------------------------------------- offline / unreachable (R6)


def test_unreachable_github_records_the_host_and_tries_no_mirror(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    pin(monkeypatch, nextflow=sha256_hex(FAKE_DIST))
    net.routes[("github.com", NF_PATH)] = Route(error="[Errno 101] Network is unreachable")
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db))
    assert code == 2, out
    assert net.contacted == [("github.com", NF_PATH)], "R6: no other source is tried"
    assert not (db / "bin" / "nextflow").exists()
    assert not no_leftovers(db)
    rec = manifest_items(db)["nextflow"]
    assert rec["status"] != "ok" and "github.com" in rec["message"], rec


def test_refused_redirect_is_recorded_and_exits_2(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    pin(monkeypatch, nextflow=sha256_hex(FAKE_DIST))
    net.routes[("github.com", NF_PATH)] = Route(
        302, {"Location": "https://codeload.github.com" + CODELOAD_PATH}
    )
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db))
    assert code == 2, out
    assert ("codeload.github.com", CODELOAD_PATH) not in net.contacted
    assert not (db / "bin" / "nextflow").exists()
    assert manifest_items(db)["nextflow"]["status"] != "ok"


# ---------------------------------------------------------------- example 5.5 and --no-nextflow


def test_example_5_5_nextflow_on_path_is_not_downloaded(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    make_exe(fakebin / "nextflow", FAKE_NF_ON_PATH)
    use_path(monkeypatch, fakebin)
    serve_nextflow(net)
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db))
    assert code == 0, out
    assert net.contacted == [], "no download when nextflow is on PATH"
    assert not (db / "bin" / "nextflow").exists()
    # "the record names the found version": setup's output, or a db_manifest item (either way)
    rec = manifest_items(db).get("nextflow")
    assert "25.10.0" in out or (rec is not None and "25.10.0" in rec["version"]), (out, rec)
    assert str(fakebin / "nextflow") in out or rec is not None, "the nextflow found is named"


def test_no_nextflow_skips_the_launcher(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    serve_nextflow(net)
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db, "--no-nextflow"))
    assert code == 0, out
    assert net.contacted == []
    assert not (db / "bin" / "nextflow").exists()


def test_explicit_uniop_never_plans_the_launcher(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Q-K27: ``gippa setup --uniop`` installs UniOP only; with no nextflow on PATH the hosts
    contacted are still exactly GitHub's archive URL and codeload (14u1 example 16)."""
    from helpers.uniop_14u1 import ACCEPT, install_fake_docker, pin_uniop_source, uniop_args

    nf_row()  # 15d: the row exists, so not planning it is a decision, not an absence
    make_exe(fakebin / "java", fake_java_script(JAVA_21))
    use_path(monkeypatch, fakebin)
    install_fake_docker(tmp_path, monkeypatch)
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    serve_nextflow(net)
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 0, out
    assert net.contacted == [("github.com", ARCHIVE_PATH), ("codeload.github.com", CODELOAD_PATH)]
    assert not (db / "bin" / "nextflow").exists()
    assert "nextflow" not in manifest_items(db)


# ---------------------------------------------------------------- example 5.6: Java warns only


COG = ("cog2020_def", "cog2020_fun")


def serve_cog(net: FakeNet) -> dict[str, str]:
    from gippa.setupdb.items import SOURCES

    pins = {}
    for name in COG:
        url = SOURCES[name].url
        body = f"{name}\tfake COG table\n".encode()
        net.routes[("ftp.ncbi.nlm.nih.gov", url.split("ftp.ncbi.nlm.nih.gov", 1)[1])] = Route(
            200, {}, body
        )
        pins[name] = sha256_hex(body)
    return pins


@pytest.mark.parametrize("java", ["none", "stub", "11"])
def test_example_5_6_setup_without_java_warns_and_still_downloads(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, java: str
) -> None:
    if java == "stub":  # the macOS /usr/bin/java without a JDK
        make_exe(fakebin / "java", fake_java_script("Unable to locate a Java Runtime.", 1))
    elif java == "11":
        make_exe(fakebin / "java", fake_java_script(JAVA_11))
    use_path(monkeypatch, fakebin, host=False)  # no host java either
    pin(monkeypatch, nextflow=sha256_hex(FAKE_DIST), **serve_cog(net))
    serve_nextflow(net)
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db, "--cog"))
    assert code == 0, out
    assert JAVA_HEAD in out, "a warning names the Java problem"
    assert "warning" in out.lower()
    assert ("found 11 at" in out) == (java == "11")
    if java != "11":
        assert "not found" in out
    items = manifest_items(db)
    assert all(items[n]["status"] == "ok" for n in (*COG, "nextflow")), items
    assert (db / "bin" / "nextflow").read_bytes() == FAKE_DIST


def test_setup_exit_follows_the_db_results_not_java(
    net: FakeNet, fakebin: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No Java and a failed item: exit 2 for the item, the Java line stays a warning."""
    use_path(monkeypatch, fakebin, host=False)
    pin(monkeypatch, nextflow=sha256_hex(b"other bytes"))
    serve_nextflow(net)
    db = tmp_path / "db"
    code, out = setup_cli(nf_args(db))
    assert code == 2, out
    assert JAVA_HEAD in out
    assert manifest_items(db)["nextflow"]["status"] != "ok"
