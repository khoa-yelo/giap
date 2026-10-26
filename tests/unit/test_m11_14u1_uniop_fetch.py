"""M11-14u1: the UniOP archive fetch, the runtime rule R(P), the image key line, the Dockerfile.

docs/specs/M11-14-container-install.md item 1b ("Command", "Runtime from the profile list",
"Fetch", "Redirect", "Permissions", "Image key", "Dockerfile (14u1)"), examples 16 (fetch part),
16a, 16e, 17, AC 9 (14u1 list; "No network, enforced"); R5, R6, R11.

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- ``db_sources.yaml`` row ``uniop_source`` exactly as item 1b lays it out. The fetch, the
  checksum check and the image key read the pin at call time from the ``SOURCES`` mapping of the
  ``giap.setupdb`` modules (tests patch it: the real archive's sha256 is computed on FarmShare).
- ``SetupOptions`` gains ``uniop: bool = False``; ``plan()`` adds ``uniop_source`` only with it.
- ``giap setup --uniop --profile P [--cache-dir DIR] [--accept-uniop-terms]``. The archive goes to
  ``<db_dir>/uniop/UniOP-<commit>.tar.gz`` (dir ``0700``, file ``0600``) through the real
  ``HttpsFetcher`` (one GET of the pinned URL, GitHub's 302 to codeload followed).
- ``REDIRECT_HOSTS["github.com"] == frozenset({"codeload.github.com"})`` in ``giap.setupdb.run``.
- ``giap.preflight.environment.uniop_runtime(profile)`` is R(P): ``"micromamba"``,
  ``"apptainer"``, ``"docker"``, or ``None`` for no runtime **and** for an ambiguous list; setup
  tells the two apart in its message (``no container or conda profile`` /
  ``ambiguous UniOP runtime for <P>: use one of conda | docker | apptainer``).
- ``key_text("uniop")`` carries one extra line ``uniop-src.tar.gz\\t100644\\t<pinned hex>``,
  sorted with the file lines; ``uncovered_copy_sources("uniop")`` accepts ``uniop-src.tar.gz``.

No test touches the network (``no_network``; ``FakeNet`` replaces ``http.client``) or the real
``$GIAP_DB_DIR``; ``docker`` is a logging fake on ``PATH``.
"""

from __future__ import annotations

import hashlib
import os
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.uniop_14u1 import (
    ACCEPT,
    ARCHIVE_PATH,
    ARCHIVE_URL,
    CODELOAD_PATH,
    COMMIT,
    DB_PATH,
    REPO,
    SKIP_DBS,
    VERSION,
    FakeNet,
    Route,
    docker_calls,
    fixture_archive,
    install_fake_docker,
    install_fake_net,
    isolate,
    manifest_items,
    mode,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
    uniop_source,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.usefixtures("no_network")

DOCKERFILE = REPO / "workflow" / "containers" / "uniop" / "Dockerfile"
UNIOP_DEF = REPO / "workflow" / "containers" / "uniop" / "uniop.def"
TAG_RE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}")


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    install_fake_docker(tmp_path, monkeypatch)
    old = os.umask(0o022)  # the modes below must not come from a strict umask
    yield
    os.umask(old)


@pytest.fixture
def net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


@pytest.fixture
def archive(monkeypatch: pytest.MonkeyPatch) -> bytes:
    """The fixture archive, with the ``uniop_source`` pin set to its sha256."""
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    return data


def uniop_args(db: Path, *extra: str, profile: str = "docker") -> list[str]:
    return ["--db-dir", str(db), *SKIP_DBS, "--uniop", "--profile", profile, *extra]


# ---------------------------------------------------------------- AC 9: no network, enforced


def test_no_network_fixture_blocks_a_real_https_fetch(tmp_path: Path) -> None:
    """Under ``no_network`` the real fetcher (no FakeNet) cannot reach any host."""
    from giap.setupdb import HttpsFetcher, SourceError

    with pytest.raises((SourceError, OSError)):
        HttpsFetcher().fetch(ARCHIVE_URL, tmp_path / "a.tar.gz")
    assert not (tmp_path / "a.tar.gz").exists()


# ---------------------------------------------------------------- the uniop_source row


def test_uniop_source_row_is_laid_out_like_the_others() -> None:
    from giap.annotate.tools import TOOLS

    row = uniop_source()
    assert row.name == "uniop_source"
    assert row.version == VERSION
    assert row.url == ARCHIVE_URL
    assert row.host == "github.com"
    assert re.fullmatch(r"sha256:[0-9a-f]{64}", row.checksum), row.checksum
    assert row.path == DB_PATH
    assert row.tool is None
    assert row.check == "sha256"
    assert row.size_gib == 0  # left at the default
    assert row.version == TOOLS["uniop"].version, "the archive commit is the registry pin"


def test_uniop_source_loads_through_load_sources() -> None:
    from giap.setupdb.items import _load_sources

    assert _load_sources()["uniop_source"] == uniop_source()


@settings(max_examples=60, deadline=None)
@given(
    bakta=st.sampled_from(["light", "full", "skip"]),
    flags=st.fixed_dictionaries(
        {
            k: st.booleans()
            for k in (
                "genomad",
                "antismash",
                "promotech_models",
                "go",
                "kegg",
                "cog",
                "contaminants",
            )
        }
    ),
)
def test_plain_setup_never_plans_the_uniop_archive(bakta: str, flags: dict[str, bool]) -> None:
    from giap.setupdb.items import SetupOptions, plan

    names = [i.name for i in plan(SetupOptions(db_dir=Path("/x"), bakta=bakta, **flags))]  # type: ignore[arg-type]
    assert "uniop_source" not in names


def test_uniop_option_plans_the_archive() -> None:
    from giap.setupdb.items import SetupOptions, plan

    opts = SetupOptions(db_dir=Path("/x"), bakta="skip", genomad=False, promotech_models=False,
                        uniop=True)  # type: ignore[call-arg]  # fmt: skip
    assert [i.name for i in plan(opts)] == ["uniop_source"]


# ---------------------------------------------------------------- R6: the redirect allowlist


def test_redirect_allowlist_adds_only_codeload_for_github() -> None:
    from giap.setupdb.run import REDIRECT_HOSTS

    assert REDIRECT_HOSTS["github.com"] == frozenset({"codeload.github.com"})
    assert REDIRECT_HOSTS["purl.obolibrary.org"] == frozenset(
        {"release.geneontology.org", "current.geneontology.org"}
    ), "the existing entries are unchanged"
    assert set(REDIRECT_HOSTS) == {"github.com", "purl.obolibrary.org"}


def test_github_302_to_codeload_is_followed_without_content_length(
    net: FakeNet, tmp_path: Path
) -> None:
    from giap.setupdb import HttpsFetcher

    data = fixture_archive()
    net.serve_github(data)
    HttpsFetcher().fetch(ARCHIVE_URL, tmp_path / "a.tar.gz")
    assert (tmp_path / "a.tar.gz").read_bytes() == data
    assert net.contacted == [("github.com", ARCHIVE_PATH), ("codeload.github.com", CODELOAD_PATH)]


@pytest.mark.parametrize(
    "location",
    [
        "https://objects.githubusercontent.com/x.tar.gz",
        "https://evil.example/x.tar.gz",
        "http://codeload.github.com" + CODELOAD_PATH,
        "https://codeload.github.com.evil.example" + CODELOAD_PATH,
    ],
)
def test_example_16a_redirect_to_any_other_host_is_refused(
    net: FakeNet, tmp_path: Path, location: str
) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("github.com", ARCHIVE_PATH)] = Route(302, {"Location": location})
    with pytest.raises(SourceError, match="refused"):
        HttpsFetcher().fetch(ARCHIVE_URL, tmp_path / "a.tar.gz")
    assert net.contacted == [("github.com", ARCHIVE_PATH)], "the new location is never contacted"
    assert not (tmp_path / "a.tar.gz").exists()


def test_codeload_is_allowed_only_as_a_github_redirect(net: FakeNet, tmp_path: Path) -> None:
    """The allowlist is per source: a zenodo source may not be sent to codeload."""
    from giap.setupdb import HttpsFetcher, SourceError

    path = "/records/1/files/x.tar.gz"
    net.routes[("zenodo.org", path)] = Route(302, {"Location": "https://codeload.github.com/x"})
    with pytest.raises(SourceError, match="refused"):
        HttpsFetcher().fetch("https://zenodo.org" + path, tmp_path / "x")
    assert ("codeload.github.com", "/x") not in net.contacted


# ---------------------------------------------------------------- example 16: the verified fetch


def test_example_16_one_download_into_the_db_dir(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    net.serve_github(archive)
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 0, out
    assert net.contacted == [("github.com", ARCHIVE_PATH), ("codeload.github.com", CODELOAD_PATH)]
    got = db / DB_PATH
    assert got.read_bytes() == archive
    assert not list((db / "uniop").glob("*.partial")), "nothing partial is left"


def test_archive_and_directory_modes(net: FakeNet, archive: bytes, tmp_path: Path) -> None:
    net.serve_github(archive)
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 0, out
    assert mode(db / "uniop") == 0o700
    assert mode(db / DB_PATH) == 0o600
    assert mode(db / "containers" / "uniop") == 0o700  # <cache>/uniop, the default cache


def test_cache_dir_option_holds_the_staging_context(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    net.serve_github(archive)
    db, cache = tmp_path / "db", tmp_path / "cache"
    code, out = setup_cli(uniop_args(db, ACCEPT, "--cache-dir", str(cache)))
    assert code == 0, out
    assert mode(cache / "uniop") == 0o700
    (build,) = docker_calls(tmp_path, "build")
    assert Path(build[-1]).resolve() == (cache / "uniop" / "context").resolve()
    assert not (db / "containers").exists()


def test_rerun_with_the_verified_archive_does_not_download_again(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    """Reading: an archive already in place that verifies is adopted (as every item is)."""
    net.serve_github(archive)
    db = tmp_path / "db"
    assert setup_cli(uniop_args(db, ACCEPT))[0] == 0
    first = list(net.contacted)
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 0, out
    assert net.contacted == first


# ---------------------------------------------------------------- example 17: checksum mismatch


def test_example_17_checksum_mismatch_deletes_and_builds_nothing(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    bad = bytearray(archive)
    bad[len(bad) // 2] ^= 0x01
    net.serve_github(bytes(bad))
    db = tmp_path / "db"
    db.mkdir()
    before = (
        "profiles {\n    docker { process { withName: 'UNIOP' { container = 'old'; "
        "conda = '' } } }\n}\n"
    )
    (db / "local_tools.config").write_text(before)
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    got, want = sha256_hex(bytes(bad)), sha256_hex(archive)
    flat = " ".join(out.split())
    assert re.search(rf"uniop: archive sha256 (sha256:)?{got} != pinned (sha256:)?{want}", flat), (
        out
    )
    assert not (db / DB_PATH).exists(), "the mismatching file is deleted"
    left = [p for p in (db / "uniop").rglob("*") if p.is_file()] if (db / "uniop").is_dir() else []
    assert left == [], "no copy of the mismatching bytes is left"
    assert docker_calls(tmp_path, "build") == [], "no build command runs"
    assert (db / "local_tools.config").read_text() == before, "local_tools.config unchanged"
    item = manifest_items(db).get("image:uniop")
    assert item is not None and item["status"] != "ok", item


@settings(
    max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(pos=st.integers(min_value=0), bit=st.integers(min_value=0, max_value=7))
def test_any_flipped_bit_fails_the_check_before_any_build(
    tmp_path_factory: pytest.TempPathFactory,
    monkeypatch: pytest.MonkeyPatch,
    pos: int,
    bit: int,
) -> None:
    tmp = tmp_path_factory.mktemp("flip")
    data = fixture_archive()
    bad = bytearray(data)
    bad[pos % len(bad)] ^= 1 << bit
    with monkeypatch.context() as mp:
        pin_uniop_source(mp, sha256_hex(data))
        install_fake_docker(tmp, mp)
        net = install_fake_net(mp)
        net.serve_github(bytes(bad))
        code, out = setup_cli(uniop_args(tmp / "db", ACCEPT))
        assert code == 2, out
        assert "uniop: archive sha256" in " ".join(out.split()), out
        assert docker_calls(tmp, "build") == []
        assert not (tmp / "db" / DB_PATH).exists()


# ---------------------------------------------------------------- offline


def test_offline_records_a_failure_that_preflight_reads(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    net.routes[("github.com", ARCHIVE_PATH)] = Route(error="[Errno 101] Network is unreachable")
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    assert docker_calls(tmp_path, "build") == []
    assert not (db / "local_tools.config").exists()
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "absent", item
    assert "github.com" in item["message"], item
    assert item["version"] == VERSION and item["url"] == ARCHIVE_URL
    from giap.preflight import environment

    record = getattr(environment, "uniop_setup_record", None)  # 14d (#387) reads it
    if record is not None:
        assert record(db) == f"(last setup: {item['message']})"


def test_refused_redirect_is_recorded_and_exits_2(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    net.routes[("github.com", ARCHIVE_PATH)] = Route(302, {"Location": "https://evil.example/x"})
    db = tmp_path / "db"
    code, out = setup_cli(uniop_args(db, ACCEPT))
    assert code == 2, out
    assert ("evil.example", "/x") not in net.contacted
    item = manifest_items(db)["image:uniop"]
    assert item["status"] != "ok" and "refused" in item["message"], item
    assert docker_calls(tmp_path, "build") == []


# ---------------------------------------------------------------- example 16e: R(P)

R_OF_P = [
    ("docker", "docker"),
    ("apptainer", "apptainer"),
    ("conda", "micromamba"),
    ("farmshare", "micromamba"),
    ("farmshare_local,apptainer", "micromamba"),
    ("farmshare,apptainer", "micromamba"),
    ("apptainer,docker", "apptainer"),
    ("conda,docker", None),  # ambiguous: missing
    ("conda,apptainer", None),  # ambiguous: missing
    ("standard", None),  # none: missing
]


def _runtime(profile: str) -> str | None:
    from giap.preflight import environment

    fn = getattr(environment, "uniop_runtime", None)
    if fn is None:
        pytest.fail("giap.preflight.environment.uniop_runtime (R(P)) is not implemented")
    return fn(profile)  # type: ignore[no-any-return]


@pytest.mark.parametrize(("profile", "want"), R_OF_P)
def test_example_16e_runtime_from_the_profile_list(profile: str, want: str | None) -> None:
    assert _runtime(profile) == want


PROFILE_NAMES = ["conda", "farmshare", "farmshare_local", "apptainer", "docker", "standard",
                 "test", "stub"]  # fmt: skip


def _expected_runtime(names: set[str]) -> str | None:
    if names & {"farmshare", "farmshare_local"}:
        return "micromamba"
    if "conda" in names:
        return None if names & {"docker", "apptainer"} else "micromamba"
    if "apptainer" in names:
        return "apptainer"
    return "docker" if "docker" in names else None


@settings(max_examples=200, deadline=None)
@given(st.lists(st.sampled_from(PROFILE_NAMES), min_size=1, max_size=5))
def test_runtime_rule_property(names: list[str]) -> None:
    """The first matching rule wins; order and repeats in the list do not matter."""
    assert _runtime(",".join(names)) == _expected_runtime(set(names))


def test_setup_ambiguous_profile_exits_2_without_fetching(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    code, out = setup_cli(uniop_args(tmp_path / "db", ACCEPT, profile="conda,docker"))
    assert code == 2, out
    flat = " ".join(out.split())
    assert "ambiguous UniOP runtime for conda,docker: use one of conda | docker | apptainer" in flat
    assert net.contacted == []
    assert docker_calls(tmp_path, "build") == []


def test_setup_without_a_runtime_exits_2_without_fetching(
    net: FakeNet, archive: bytes, tmp_path: Path
) -> None:
    code, out = setup_cli(uniop_args(tmp_path / "db", ACCEPT, profile="standard"))
    assert code == 2, out
    assert "no container or conda profile" in " ".join(out.split())
    assert net.contacted == []


# ---------------------------------------------------------------- image key (14a, 14u1)


def _pin_hex() -> str:
    return str(uniop_source().checksum).removeprefix("sha256:")


def test_uniop_key_text_has_the_archive_line() -> None:
    from giap.io.image_key import key_text

    lines = key_text("uniop").splitlines(keepends=True)
    assert f"uniop-src.tar.gz\t100644\t{_pin_hex()}\n" in lines
    paths = [ln.split("\t", 1)[0] for ln in lines]
    assert paths == sorted(paths), "the archive line is sorted with the file lines"
    assert len(paths) == len(set(paths))


def test_uniop_def_is_a_key_input() -> None:
    from giap.io.image_key import key_inputs

    assert UNIOP_DEF.is_file(), "workflow/containers/uniop/uniop.def is new in 14u1"
    assert "workflow/containers/uniop/uniop.def" in key_inputs("uniop")


def test_other_tools_have_no_archive_line() -> None:
    from giap.io.image_key import key_text

    for tool in ("ostir", "promotech", "fgsea"):
        assert "uniop-src.tar.gz" not in key_text(tool), tool


def test_the_key_is_sha256_of_the_text_with_the_archive_line() -> None:
    from giap.io.image_key import key, key_text

    assert key("uniop") == hashlib.sha256(key_text("uniop").encode()).hexdigest()[:12]


@settings(
    max_examples=30, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(a=st.binary(min_size=32, max_size=32), b=st.binary(min_size=32, max_size=32))
def test_a_pin_bump_changes_the_key_and_the_tag(
    monkeypatch: pytest.MonkeyPatch, a: bytes, b: bytes
) -> None:
    from giap.io.image_key import key, key_text, tag

    with monkeypatch.context() as mp:
        pin_uniop_source(mp, a.hex())
        ka, text = key("uniop"), key_text("uniop")
        assert f"uniop-src.tar.gz\t100644\t{a.hex()}\n" in text
    with monkeypatch.context() as mp:
        pin_uniop_source(mp, b.hex())
        kb = key("uniop")
    assert (ka == kb) == (a == b)
    t = tag(VERSION, ka)
    assert t == f"git-{COMMIT[:8]}-{ka}" and TAG_RE.fullmatch(t)


def test_changing_uniop_def_changes_the_key(tmp_path: Path) -> None:
    from giap.io.image_key import key

    root = tmp_path / "root"
    shutil.copytree(REPO / "workflow" / "containers" / "uniop", root / "workflow/containers/uniop")
    (root / "workflow" / "envs" / "locks").mkdir(parents=True)
    shutil.copy(REPO / "workflow/envs/uniop.yml", root / "workflow/envs/uniop.yml")
    before = key("uniop", root)
    d = root / "workflow/containers/uniop/uniop.def"
    assert d.is_file()
    d.write_text(d.read_text() + "\n# changed\n")
    assert key("uniop", root) != before


# ---------------------------------------------------------------- Dockerfile (14u1)


def test_dockerfile_adds_no_url() -> None:
    for line in DOCKERFILE.read_text(encoding="utf-8").splitlines():
        words = line.strip().split()
        if words and words[0].upper() == "ADD":
            assert not re.search(r"\b(https?|ftp)://", line), line
        assert not re.search(r"\b(curl|wget)\b", line) or line.lstrip().startswith("#"), line


def test_dockerfile_copies_the_verified_archive() -> None:
    sources: list[str] = []
    for line in DOCKERFILE.read_text(encoding="utf-8").splitlines():
        words = line.strip().split()
        if words and words[0].upper() == "COPY":
            sources += [w for w in words[1:-1] if not w.startswith("--")]
    assert "uniop-src.tar.gz" in sources
    assert COMMIT in DOCKERFILE.read_text(encoding="utf-8"), "GIAP_COMMIT is still the pin"


def test_dockerfile_header_says_local_build_only() -> None:
    header = " ".join(
        ln.lstrip("# ").strip()
        for ln in DOCKERFILE.read_text(encoding="utf-8").splitlines()
        if ln.startswith("#")
    )
    assert "built only locally by" in header and "setup --uniop" in header, header
    assert "never pushed" in header, header


def test_copy_sources_are_covered_by_the_key() -> None:
    from giap.io.image_key import uncovered_copy_sources

    assert uncovered_copy_sources("uniop") == []
