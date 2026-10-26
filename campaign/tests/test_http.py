"""https-only JSON GET with cache, checked downloads and safe zip extraction (campaign/SPEC.md §2).

Everything is offline: ``get_json`` gets a fake ``opener`` and ``download`` a fake ``run`` that
plays curl by writing the ``-o`` target.
"""

from __future__ import annotations

import email.message
import hashlib
import itertools
import json
import logging
import stat
import subprocess
import tempfile
import urllib.error
import zipfile
from pathlib import Path
from typing import Any

import pytest
from giap_campaign.http import DownloadError, NotFound, download, get_json, safe_extract
from hypothesis import given, settings
from hypothesis import strategies as st

URL = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001"
NCBI_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/562"
FILE_URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD000001/evidence.txt"


# --- fakes --------------------------------------------------------------------------------------


class FakeResponse:
    def __init__(self, body: bytes) -> None:
        self._body = body
        self.status = 200

    def read(self, *args: Any) -> bytes:
        return self._body

    def getcode(self) -> int:
        return 200

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *exc: object) -> None:
        return None


def _http_error(url: str, code: int) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(url, code, f"status {code}", email.message.Message(), None)


class FakeOpener:
    """Replays ``outcomes`` in order (bytes → response, exception → raised); records requests."""

    def __init__(self, *outcomes: bytes | BaseException) -> None:
        self.outcomes = list(outcomes)
        self.requests: list[Any] = []

    @property
    def calls(self) -> int:
        return len(self.requests)

    def __call__(self, request: Any, *args: Any, **kwargs: Any) -> FakeResponse:
        self.requests.append(request)
        if not self.outcomes:
            raise AssertionError("opener called more often than expected")
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, BaseException):
            raise outcome
        return FakeResponse(outcome)


def _url_of(request: Any) -> str:
    return str(getattr(request, "full_url", request))


def _headers_of(request: Any) -> dict[str, str]:
    if isinstance(request, str):
        return {}
    return {k.lower(): v for k, v in request.header_items()}


class FakeCurl:
    """Stands in for ``subprocess.run``: writes ``content`` to the path after ``-o``."""

    def __init__(self, content: bytes | None, returncode: int = 0) -> None:
        self.content = content
        self.returncode = returncode
        self.cmds: list[list[str]] = []

    def __call__(self, cmd: Any, *args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        cmd = [str(c) for c in cmd]
        self.cmds.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1])
        if self.content is not None:
            out.write_bytes(self.content)
        if self.returncode != 0 and kwargs.get("check"):
            raise subprocess.CalledProcessError(self.returncode, cmd)
        return subprocess.CompletedProcess(cmd, self.returncode, "", "")


def _part(dest: Path) -> Path:
    return Path(str(dest) + ".part")


# --- get_json -----------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "url",
    [
        "http://www.ebi.ac.uk/pride/ws/archive/v3/projects",
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/x.json",
        "file:///etc/passwd",
        "www.ebi.ac.uk/pride/ws/archive/v3/projects",
    ],
)
def test_get_json_rejects_non_https(url: str) -> None:
    opener = FakeOpener(b"{}")
    with pytest.raises(ValueError):
        get_json(url, opener=opener)
    assert opener.calls == 0


def test_get_json_parses_body() -> None:
    opener = FakeOpener(b'{"accession": "PXD000001", "n": [1, 2]}')
    assert get_json(URL, opener=opener) == {"accession": "PXD000001", "n": [1, 2]}
    assert opener.calls == 1
    assert _url_of(opener.requests[0]) == URL


def test_get_json_cache_hit_avoids_second_fetch(tmp_path: Path) -> None:
    body = b'[{"accession": "PXD000001"}]'
    opener = FakeOpener(body)
    first = get_json(URL, cache_dir=tmp_path, opener=opener)
    second = get_json(URL, cache_dir=tmp_path, opener=opener)
    assert first == second == [{"accession": "PXD000001"}]
    assert opener.calls == 1
    cached = tmp_path / f"{hashlib.sha256(URL.encode()).hexdigest()}.json"
    assert cached.is_file()
    assert json.loads(cached.read_text(encoding="utf-8")) == first


def test_get_json_cache_is_keyed_by_url(tmp_path: Path) -> None:
    opener = FakeOpener(b'{"a": 1}', b'{"b": 2}')
    assert get_json(URL, cache_dir=tmp_path, opener=opener) == {"a": 1}
    assert get_json(URL + "/files/all", cache_dir=tmp_path, opener=opener) == {"b": 2}
    assert opener.calls == 2


def test_get_json_without_cache_dir_fetches_each_time() -> None:
    opener = FakeOpener(b"1", b"2")
    assert get_json(URL, opener=opener) == 1
    assert get_json(URL, opener=opener) == 2
    assert opener.calls == 2


@pytest.mark.parametrize(
    "transient",
    [
        urllib.error.URLError("connection reset"),
        _http_error(URL, 500),
        _http_error(URL, 503),
        TimeoutError("timed out"),
    ],
    ids=["urlerror", "http500", "http503", "timeout"],
)
def test_get_json_retries_once_then_succeeds(transient: BaseException) -> None:
    opener = FakeOpener(transient, b'{"ok": true}')
    assert get_json(URL, opener=opener) == {"ok": True}
    assert opener.calls == 2


def test_get_json_fails_after_one_retry() -> None:
    opener = FakeOpener(_http_error(URL, 502), _http_error(URL, 502), b"{}")
    with pytest.raises(Exception) as info:
        get_json(URL, opener=opener)
    assert not isinstance(info.value, AssertionError)
    assert opener.calls == 2


def test_get_json_retries_zero_means_single_attempt() -> None:
    opener = FakeOpener(urllib.error.URLError("down"), b"{}")
    with pytest.raises(Exception) as info:
        get_json(URL, retries=0, opener=opener)
    assert not isinstance(info.value, AssertionError)
    assert opener.calls == 1


def test_not_found_is_a_lookup_error() -> None:
    assert issubclass(NotFound, LookupError)


def test_get_json_404_raises_not_found_without_retry() -> None:
    opener = FakeOpener(_http_error(URL, 404), b"{}")
    with pytest.raises(NotFound):
        get_json(URL, opener=opener)
    assert opener.calls == 1


def test_get_json_other_4xx_not_retried() -> None:
    opener = FakeOpener(_http_error(URL, 403), b"{}")
    with pytest.raises(Exception) as info:
        get_json(URL, opener=opener)
    assert not isinstance(info.value, (AssertionError, NotFound))
    assert opener.calls == 1


def test_get_json_errors_are_not_cached(tmp_path: Path) -> None:
    opener = FakeOpener(_http_error(URL, 404), b'{"late": 1}')
    with pytest.raises(NotFound):
        get_json(URL, cache_dir=tmp_path, opener=opener)
    assert get_json(URL, cache_dir=tmp_path, opener=opener) == {"late": 1}
    assert opener.calls == 2


def test_get_json_sends_ncbi_api_key_header_and_never_logs_it(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, tmp_path: Path
) -> None:
    secret = "s3cr3t-key-0123456789"
    monkeypatch.setenv("NCBI_API_KEY", secret)
    opener = FakeOpener(b'{"taxonomy_nodes": []}')
    with caplog.at_level(logging.DEBUG):
        get_json(NCBI_URL, cache_dir=tmp_path, opener=opener)
    assert _headers_of(opener.requests[0]).get("api-key") == secret
    assert secret not in caplog.text
    assert all(secret not in p.name for p in tmp_path.iterdir())


def test_get_json_no_api_key_header_without_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    opener = FakeOpener(b"{}")
    get_json(NCBI_URL, opener=opener)
    assert "api-key" not in _headers_of(opener.requests[0])


# --- download -----------------------------------------------------------------------------------

CONTENT = b"Sequence\tRaw file\tProteins\nPEPTIDEK\trun1\tP1\n"
DIGESTS = {
    "md5": hashlib.md5(CONTENT).hexdigest(),
    "sha1": hashlib.sha1(CONTENT).hexdigest(),
    "sha256": hashlib.sha256(CONTENT).hexdigest(),
}


def test_download_runs_curl_https_only(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    curl = FakeCurl(CONTENT)
    download(FILE_URL, dest, expected_bytes=None, checksum=None, run=curl)
    assert len(curl.cmds) == 1
    cmd = curl.cmds[0]
    assert cmd[0] == "curl"
    assert "-fsSL" in cmd
    pairs = list(itertools.pairwise(cmd))
    assert ("--proto", "=https") in pairs
    assert ("--proto-redir", "=https") in pairs
    assert ("--retry", "1") in pairs
    assert ("-o", str(_part(dest))) in pairs
    assert cmd[-1] == FILE_URL


@pytest.mark.parametrize("algo", ["md5", "sha1", "sha256"])
@pytest.mark.parametrize("upper", [False, True], ids=["lower", "upper"])
def test_download_good_checksum_renames_and_returns_sha256(
    tmp_path: Path, algo: str, upper: bool
) -> None:
    dest = tmp_path / "evidence.txt"
    checksum = DIGESTS[algo].upper() if upper else DIGESTS[algo]
    got = download(
        FILE_URL, dest, expected_bytes=len(CONTENT), checksum=checksum, run=FakeCurl(CONTENT)
    )
    assert got == DIGESTS["sha256"]
    assert dest.read_bytes() == CONTENT
    assert not _part(dest).exists()


@pytest.mark.parametrize("algo", ["md5", "sha1", "sha256"])
def test_download_checksum_mismatch_raises_and_removes_part(tmp_path: Path, algo: str) -> None:
    dest = tmp_path / "evidence.txt"
    wrong = getattr(hashlib, algo)(b"something else").hexdigest()
    with pytest.raises(DownloadError):
        download(FILE_URL, dest, expected_bytes=None, checksum=wrong, run=FakeCurl(CONTENT))
    assert not _part(dest).exists()
    assert not dest.exists()


def test_download_size_mismatch_raises_and_removes_part(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    with pytest.raises(DownloadError):
        download(
            FILE_URL,
            dest,
            expected_bytes=len(CONTENT) + 1,
            checksum=DIGESTS["sha256"],
            run=FakeCurl(CONTENT),
        )
    assert not _part(dest).exists()
    assert not dest.exists()


@pytest.mark.parametrize(
    "checksum",
    ["", "not-a-checksum", "0" * 31, "0" * 33, "g" * 32, "0" * 128],
    ids=["empty", "word", "31hex", "33hex", "32nonhex", "128hex"],
)
def test_download_unrecognised_checksum_is_ignored(tmp_path: Path, checksum: str) -> None:
    dest = tmp_path / "evidence.txt"
    got = download(
        FILE_URL, dest, expected_bytes=len(CONTENT), checksum=checksum, run=FakeCurl(CONTENT)
    )
    assert got == DIGESTS["sha256"]
    assert dest.read_bytes() == CONTENT


def test_download_no_checks_still_returns_sha256(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    got = download(FILE_URL, dest, expected_bytes=None, checksum=None, run=FakeCurl(CONTENT))
    assert got == DIGESTS["sha256"]
    assert dest.read_bytes() == CONTENT


def test_download_existing_valid_dest_is_not_downloaded_again(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    dest.write_bytes(CONTENT)
    curl = FakeCurl(b"different bytes")
    got = download(FILE_URL, dest, expected_bytes=len(CONTENT), checksum=DIGESTS["md5"], run=curl)
    assert curl.cmds == []
    assert got == DIGESTS["sha256"]
    assert dest.read_bytes() == CONTENT


def test_download_existing_invalid_dest_is_replaced(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    dest.write_bytes(b"truncated")
    curl = FakeCurl(CONTENT)
    got = download(FILE_URL, dest, expected_bytes=len(CONTENT), checksum=DIGESTS["sha1"], run=curl)
    assert len(curl.cmds) == 1
    assert got == DIGESTS["sha256"]
    assert dest.read_bytes() == CONTENT


def test_download_curl_failure_raises(tmp_path: Path) -> None:
    dest = tmp_path / "evidence.txt"
    with pytest.raises((DownloadError, subprocess.CalledProcessError)):
        download(
            FILE_URL, dest, expected_bytes=None, checksum=None, run=FakeCurl(None, returncode=22)
        )
    assert not dest.exists()


@settings(max_examples=30, deadline=None)
@given(content=st.binary(min_size=0, max_size=256), algo=st.sampled_from(["md5", "sha1", "sha256"]))
def test_download_property_matching_digest_always_accepted(content: bytes, algo: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        dest = Path(d) / "f.bin"
        checksum = getattr(hashlib, algo)(content).hexdigest()
        got = download(
            FILE_URL, dest, expected_bytes=len(content), checksum=checksum, run=FakeCurl(content)
        )
        assert got == hashlib.sha256(content).hexdigest()
        assert dest.read_bytes() == content


# --- safe_extract -------------------------------------------------------------------------------


def _zip(path: Path, members: dict[str, bytes]) -> Path:
    with zipfile.ZipFile(path, "w") as zf:
        for name, data in members.items():
            zf.writestr(name, data)
    return path


def test_safe_extract_extracts_members(tmp_path: Path) -> None:
    members = {
        "README.md": b"readme\n",
        "ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna": b">c\nACGT\n",
        "ncbi_dataset/data/GCF_000005845.2/genomic.gbff": b"LOCUS x\n//\n",
    }
    archive = _zip(tmp_path / "genome.zip", members)
    out = tmp_path / "out"
    paths = safe_extract(archive, out)
    files = {p.resolve() for p in paths if p.is_file()}
    assert files == {(out / name).resolve() for name in members}
    for name, data in members.items():
        assert (out / name).read_bytes() == data


@pytest.mark.parametrize(
    "bad",
    ["../evil.txt", "ok/../../evil.txt", "/abs_evil_giap_campaign_test.txt"],
    ids=["dotdot", "nested-dotdot", "absolute"],
)
def test_safe_extract_refuses_unsafe_paths(tmp_path: Path, bad: str) -> None:
    archive = _zip(tmp_path / "bad.zip", {"fine.txt": b"x", bad: b"evil"})
    out = tmp_path / "a" / "out"
    out.mkdir(parents=True)
    with pytest.raises(DownloadError):
        safe_extract(archive, out)
    assert not (tmp_path / "a" / "evil.txt").exists()
    assert not (tmp_path / "evil.txt").exists()
    assert not Path("/abs_evil_giap_campaign_test.txt").exists()


def test_safe_extract_refuses_symlinks(tmp_path: Path) -> None:
    archive = tmp_path / "link.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("fine.txt", b"x")
        info = zipfile.ZipInfo("link")
        info.create_system = 3
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        zf.writestr(info, "/etc/passwd")
    out = tmp_path / "out"
    with pytest.raises(DownloadError):
        safe_extract(archive, out)
    assert not (out / "link").is_symlink()
