"""``validation/fetch_regulondb.py`` without the network (M6-11; GIPPA_SPEC.md §3.9, §15.3).

Nothing here opens a socket: the one network call of the TLS path is injected, the GraphQL
transport is stubbed at ``urllib.request.urlopen``, and no RegulonDB bytes are committed (R6/R11)
— the payloads below are tables we wrote.

What the tests hold the fetcher to:

* the SHA-256 pin on the GlobalSign intermediate is refused when it changes, and the message
  names the pin, so a 2029 expiry is obvious rather than mysterious;
* the intermediate is only ever *added* to the default context — verification is never disabled
  (no ``verify=False``, no ``CERT_NONE``, no ``_create_unverified_context``, and the source is
  grepped for those);
* the provenance JSON carries the URL, release, byte count, SHA-256 and every retrieval date;
* writes are atomic and leave no ``.part`` behind.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib
import io
import json
import re
import shutil
import ssl
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
fr: Any = importlib.import_module("validation.fetch_regulondb")

TU_SET = (
    "1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes\n"
    "TU01\tlacZYA\tOP01\tlacZYA\tlacZ;lacY;lacA\n"
)
TERMINATOR_SET = "1)id\t2)leftEndPos\t3)rightEndPos\t4)strand\nTM01\t10\t40\tforward\n"
GENE_SET = (
    "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\n"
    "G1\tlacZ\tb1\t1\t9\tforward\n"
)
PAYLOADS = {"TUSet": TU_SET, "TerminatorSet": TERMINATOR_SET, "GeneProductSet": GENE_SET}

RELEASES = [
    {
        "regulonDBVersion": "14.5.0",
        "ecocycVersion": "29.5",
        "releaseDate": "2026-01-28",
        "genomeVersion": "NC_000913.3",
    },
    {
        "regulonDBVersion": "12.0.0",
        "ecocycVersion": "27.0",
        "releaseDate": "2023-10-23",
        "genomeVersion": "NC_000913.3",
    },
]


class _Response(io.BytesIO):
    """The slice of ``http.client.HTTPResponse`` ``urlopen`` users touch."""

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def _stub_urlopen(captured: list[dict[str, Any]]):  # type: ignore[no-untyped-def]
    def fake(req: Any, *_: object, **kwargs: object) -> _Response:
        body = json.loads(req.data.decode())
        captured.append({"url": req.full_url, "body": body, "context": kwargs.get("context")})
        name = body["variables"].get("fileName")
        if name is None:
            data = {"getDatabaseInfo": RELEASES}
        else:
            data = {
                "getDataOfFile": {
                    "fileName": name,
                    "fileFormat": "rif-version 1",
                    "version": "1.0",
                    "rdbVersion": "14.5",
                    "creationDate": "06-09-2026",
                    "license": "CC BY 4.0",
                    "citation": "RegulonDB",
                    "columnsDetails": "# Columns:",
                    "content": PAYLOADS[name],
                }
            }
        return _Response(json.dumps({"data": data}).encode())

    return fake


# ----------------------------------------------------------------- the pin


def _self_signed_der() -> bytes:
    """A DER blob standing in for the intermediate; only its bytes matter to the pin check."""
    return b"\x30\x82\x00\x04not-a-real-certificate"


def test_the_pin_is_the_sha256_the_backlog_records() -> None:
    assert fr.INTERMEDIATE_SHA256 == (
        "d160e2de4e56cb10b66c1cb0adcb79cf93c78dcd9db30c2018220262c04063f9"
    )
    assert fr.INTERMEDIATE_URL.startswith("http://secure.globalsign.com/cacert/")
    assert fr.ROOT_CN == "GlobalSign Root R46"


def test_verify_intermediate_refuses_a_changed_certificate() -> None:
    with pytest.raises(fr.FetchError) as exc:
        fr.verify_intermediate(_self_signed_der())
    message = str(exc.value)
    assert fr.INTERMEDIATE_SHA256 in message  # the pin is named, so an expiry is obvious
    assert hashlib.sha256(_self_signed_der()).hexdigest() in message


def test_regulondb_context_refuses_pem_from_the_ca_endpoint() -> None:
    with pytest.raises(fr.FetchError, match="PEM, not the DER"):
        fr.regulondb_context(fetch=lambda _url: b"-----BEGIN CERTIFICATE-----\n")


def test_regulondb_context_refuses_a_certificate_that_fails_the_pin() -> None:
    with pytest.raises(fr.FetchError, match="sha256"):
        fr.regulondb_context(fetch=lambda _url: _self_signed_der())


def _pin(monkeypatch: pytest.MonkeyPatch, der: bytes) -> None:
    """Make ``der`` satisfy the pin, so the test can reach the chain check behind it."""
    monkeypatch.setattr(fr, "INTERMEDIATE_SHA256", hashlib.sha256(der).hexdigest())


def _openssl_results(
    monkeypatch: pytest.MonkeyPatch, verify: tuple[int, str], issuer: str
) -> list[list[str]]:
    """Stub ``openssl``: ``verify`` gives ``(code, output)``, ``x509 -issuer`` gives ``issuer``."""
    calls: list[list[str]] = []

    def fake(argv: list[str], **_: object) -> subprocess.CompletedProcess[str]:
        calls.append(argv)
        if argv[1] == "verify":
            return subprocess.CompletedProcess(argv, verify[0], verify[1], "")
        return subprocess.CompletedProcess(argv, 0, issuer, "")

    monkeypatch.setattr(fr.subprocess, "run", fake)
    return calls


def test_verify_intermediate_checks_the_chain_against_pythons_own_trust_store(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Lead decision 6 (supersedes "no -CAfile"): ``openssl verify`` uses Python's own store.

    The handshake trusts ``ssl.get_default_verify_paths()``; the chain check must use the same
    ``-CAfile`` / ``-CApath``, or the two checks could disagree about which roots are trusted.
    """
    der = _self_signed_der()
    _pin(monkeypatch, der)
    calls = _openssl_results(monkeypatch, (1, "error 20 at 0 depth lookup"), "")

    with pytest.raises(fr.FetchError, match="does not chain to a trusted root"):
        fr.verify_intermediate(der)

    argv = calls[0]
    assert argv[1] == "verify"
    paths = ssl.get_default_verify_paths()
    passed = {flag: argv[argv.index(flag) + 1] for flag in ("-CAfile", "-CApath") if flag in argv}
    assert passed, f"openssl verify must name Python's store via -CAfile/-CApath: {argv}"
    if "-CAfile" in passed:
        assert passed["-CAfile"] == paths.cafile
    if "-CApath" in passed:
        assert passed["-CApath"] == paths.capath


def test_verify_intermediate_requires_the_expected_root(monkeypatch: pytest.MonkeyPatch) -> None:
    der = _self_signed_der()
    _pin(monkeypatch, der)
    _openssl_results(monkeypatch, (0, "ok"), "issuer=C = BE, O = Someone Else, CN = Other Root\n")

    with pytest.raises(fr.FetchError, match="not to 'GlobalSign Root R46'"):
        fr.verify_intermediate(der)


def test_verify_intermediate_accepts_a_pinned_certificate_under_the_right_root(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    der = _self_signed_der()
    _pin(monkeypatch, der)
    _openssl_results(
        monkeypatch, (0, "intermediate.pem: OK"), f"issuer=C = BE, CN = {fr.ROOT_CN}\n"
    )

    pem = fr.verify_intermediate(der)

    assert pem.startswith("-----BEGIN CERTIFICATE-----")
    assert ssl.PEM_cert_to_DER_cert(pem) == der


def test_verify_intermediate_without_openssl_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    der = _self_signed_der()
    _pin(monkeypatch, der)
    monkeypatch.setattr(fr.shutil, "which", lambda _name: None)
    with pytest.raises(fr.FetchError, match="openssl is not on PATH"):
        fr.verify_intermediate(der)


def test_verification_is_never_disabled_anywhere_in_the_fetcher() -> None:
    source = (REPO / "validation" / "fetch_regulondb.py").read_text(encoding="utf-8")
    for forbidden in (
        "_create_unverified_context",
        "CERT_NONE",
        "check_hostname = False",
        "verify=False",
        "--insecure",
    ):
        assert forbidden not in source, f"{forbidden} must never appear"
    # the spaced spelling above misses ``check_hostname=False`` (a keyword or a tight assignment)
    assert not re.search(r"check_hostname\s*=\s*False", source)
    assert not re.search(r"verify_mode\s*=\s*(ssl\.)?CERT_(NONE|OPTIONAL)", source)


# -------------------------------------------------------------- the GraphQL


def test_fetch_dataset_returns_the_content_and_the_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: list[dict[str, Any]] = []
    monkeypatch.setattr(fr.urllib.request, "urlopen", _stub_urlopen(captured))
    ctx = ssl.create_default_context()

    ds = fr.fetch_dataset("TUSet", ctx)

    assert ds.name == "TUSet"
    assert ds.content == TU_SET.encode()
    assert ds.sha256 == hashlib.sha256(TU_SET.encode()).hexdigest()
    assert ds.rdb_version == "14.5"
    assert captured[0]["url"] == fr.GRAPHQL_URL
    assert captured[0]["context"] is ctx  # the completed chain is used for every request


def test_graphql_raises_on_an_error_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake(req: Any, *_: object, **__: object) -> _Response:
        return _Response(json.dumps({"errors": [{"message": "boom"}]}).encode())

    monkeypatch.setattr(fr.urllib.request, "urlopen", fake)
    with pytest.raises(fr.FetchError, match="GraphQL error"):
        fr.graphql("{ x }", {}, ssl.create_default_context())


def test_fetch_dataset_raises_on_empty_content(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake(req: Any, *_: object, **__: object) -> _Response:
        return _Response(json.dumps({"data": {"getDataOfFile": {"content": "  "}}}).encode())

    monkeypatch.setattr(fr.urllib.request, "urlopen", fake)
    with pytest.raises(fr.FetchError, match="came back empty"):
        fr.fetch_dataset("TUSet", ssl.create_default_context())


def test_latest_release_picks_the_newest_row() -> None:
    assert fr.latest_release(RELEASES)["regulonDBVersion"] == "14.5.0"
    assert fr.latest_release(RELEASES[::-1])["regulonDBVersion"] == "14.5.0"
    assert fr.latest_release(RELEASES[1])["regulonDBVersion"] == "12.0.0"


def test_latest_release_raises_when_there_is_none() -> None:
    with pytest.raises(fr.FetchError, match="no release"):
        fr.latest_release([])


def test_the_three_datasets_cover_tus_terminators_and_the_gene_coordinates() -> None:
    assert set(fr.DATASETS) == {"TUSet", "TerminatorSet", "GeneProductSet"}


# ------------------------------------------------------------- provenance


def _datasets() -> list[Any]:
    return [
        fr.Dataset(
            name=name,
            content=text.encode(),
            file_format="rif-version 1",
            version="1.0",
            rdb_version="14.5",
            creation_date="06-09-2026",
            license="CC BY 4.0",
            citation="RegulonDB",
            columns_details="# Columns:",
        )
        for name, text in PAYLOADS.items()
    ]


def test_provenance_records_url_release_bytes_and_hashes() -> None:
    now = dt.datetime(2026, 10, 21, 20, 5, 2, tzinfo=dt.UTC)
    meta = fr.provenance(_datasets(), RELEASES[0], now)

    assert meta["url"] == fr.GRAPHQL_URL
    assert meta["regulondb_version"] == "14.5.0"
    assert meta["release_date"] == "2026-01-28"
    assert meta["retrieved"] == "2026-10-21"
    assert meta["retrieved_at"] == "2026-10-21T20:05:02Z"
    assert meta["tls"]["intermediate_sha256"] == fr.INTERMEDIATE_SHA256
    assert meta["tls"]["anchored_to"] == fr.ROOT_CN
    files = {f["dataset"]: f for f in meta["files"]}
    assert set(files) == set(PAYLOADS)
    assert files["TUSet"]["bytes"] == len(TU_SET.encode())
    assert files["TUSet"]["sha256"] == hashlib.sha256(TU_SET.encode()).hexdigest()
    assert files["TUSet"]["file"] == "TUSet.txt"
    json.dumps(meta)  # the record must be JSON-serialisable as written


def test_provenance_keeps_every_retrieval_date(tmp_path: Path) -> None:
    """Spec example 6: fetched twice, identical hashes, both dates recorded."""
    first = fr.provenance(_datasets(), RELEASES[0], dt.datetime(2026, 10, 21, tzinfo=dt.UTC))
    second = fr.provenance(
        _datasets(), RELEASES[0], dt.datetime(2026, 11, 2, tzinfo=dt.UTC), previous=first
    )
    again = fr.provenance(
        _datasets(), RELEASES[0], dt.datetime(2026, 11, 2, tzinfo=dt.UTC), previous=second
    )

    assert first["retrieved_dates"] == ["2026-10-21"]
    assert second["retrieved_dates"] == ["2026-10-21", "2026-11-02"]
    assert again["retrieved_dates"] == ["2026-10-21", "2026-11-02"]  # same day, no duplicate
    assert [f["sha256"] for f in first["files"]] == [f["sha256"] for f in second["files"]]


# ---------------------------------------------------------------- writing


def test_write_atomically_leaves_no_part_file(tmp_path: Path) -> None:
    target = tmp_path / "TUSet.txt"
    fr.write_atomically(target, b"one")
    fr.write_atomically(target, b"two")
    assert target.read_bytes() == b"two"
    assert list(tmp_path.iterdir()) == [target]


def test_release_dir_sanitises_the_release_string(tmp_path: Path) -> None:
    assert fr.release_dir(tmp_path, "14.5.0") == tmp_path / "regulondb" / "14.5.0"
    assert fr.release_dir(tmp_path, "../etc") == tmp_path / "regulondb" / ".._etc"
    assert fr.release_dir(tmp_path, "") == tmp_path / "regulondb" / "unknown"


def test_main_writes_the_release_directory_and_the_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    captured: list[dict[str, Any]] = []
    monkeypatch.setenv("GIPPA_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(fr, "regulondb_context", lambda: ssl.create_default_context())
    monkeypatch.setattr(fr.urllib.request, "urlopen", _stub_urlopen(captured))

    assert fr.main() == 0

    out = tmp_path / "regulondb" / "14.5.0"
    assert (out / "TUSet.txt").read_text() == TU_SET
    assert (out / "TerminatorSet.txt").read_text() == TERMINATOR_SET
    assert (out / "GeneProductSet.txt").read_text() == GENE_SET
    meta = json.loads((out / "provenance.json").read_text())
    assert meta["regulondb_version"] == "14.5.0"
    assert len(meta["retrieved_dates"]) == 1
    assert "RegulonDB 14.5.0" in capsys.readouterr().out
    assert not list(out.glob("*.part"))


def test_main_is_idempotent_and_appends_the_second_date(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIPPA_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(fr, "regulondb_context", lambda: ssl.create_default_context())
    monkeypatch.setattr(fr.urllib.request, "urlopen", _stub_urlopen([]))
    assert fr.main() == 0
    out = tmp_path / "regulondb" / "14.5.0"
    first = json.loads((out / "provenance.json").read_text())

    class _Later(dt.datetime):
        @classmethod
        def now(cls, tz: dt.tzinfo | None = None) -> dt.datetime:
            return dt.datetime(2026, 11, 2, 9, 0, 0, tzinfo=tz)

    monkeypatch.setattr(fr.dt, "datetime", _Later)
    assert fr.main() == 0
    second = json.loads((out / "provenance.json").read_text())

    assert [f["sha256"] for f in first["files"]] == [f["sha256"] for f in second["files"]]
    assert second["retrieved_dates"] == [first["retrieved"], "2026-11-02"]


def test_main_fails_without_gippa_data_dir(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.delenv("GIPPA_DATA_DIR", raising=False)
    assert fr.main() == 1
    assert "GIPPA_DATA_DIR is not set" in capsys.readouterr().err


def test_main_reports_a_broken_chain_instead_of_downgrading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GIPPA_DATA_DIR", str(tmp_path))

    def boom() -> ssl.SSLContext:
        raise fr.FetchError("the GlobalSign intermediate changed")

    monkeypatch.setattr(fr, "regulondb_context", boom)
    assert fr.main() == 1
    err = capsys.readouterr().err
    assert "verified chain" in err and fr.HOST in err
    assert not (tmp_path / "regulondb").exists()


# ------------------------------------------------------------------- R6/R11


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not on PATH")
def test_no_regulondb_bytes_are_committed() -> None:
    """R6/R11 and the §15.3 row: the data lives in ``$GIPPA_DATA_DIR``, never in the repo."""
    names = {"TUSet.txt", "TerminatorSet.txt", "GeneProductSet.txt", "OperonSet.txt"}
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True
    ).stdout.splitlines()
    offenders = [p for p in tracked if Path(p).name in names or Path(p).parts[:1] == ("regulondb",)]
    assert offenders == []
    assert not any(p.startswith("validation/regulondb/") for p in tracked)


# ------------------------------------------- lead decision 6: the returned context

# A throwaway self-signed CA generated for this test only (its key was discarded). It stands in
# for the GlobalSign intermediate so ``load_verify_locations`` sees a real certificate; it is not
# RegulonDB or GlobalSign material.
_TEST_CA_PEM = """-----BEGIN CERTIFICATE-----
MIIBjzCCATWgAwIBAgIUenhfRAXNvOJTtqpcISZ6mh+2Gl4wCgYIKoZIzj0EAwIw
HDEaMBgGA1UEAwwRR0lBUCB0ZXN0LW9ubHkgQ0EwIBcNMjYxMDIxMjIwNDIxWhgP
MjEyNjA5MjcyMjA0MjFaMBwxGjAYBgNVBAMMEUdJQVAgdGVzdC1vbmx5IENBMFkw
EwYHKoZIzj0CAQYIKoZIzj0DAQcDQgAEmgMd0yo6IzYs94NVdu3ZoRKEL9sqnqlF
4GcDPDB0LJy9tdNJSfKRVojv3Op1o0VSuX7AwpQqwweUTSHBR3EncaNTMFEwHQYD
VR0OBBYEFKG+SDLUOtucSYR5N1v+Xo7LliQjMB8GA1UdIwQYMBaAFKG+SDLUOtuc
SYR5N1v+Xo7LliQjMA8GA1UdEwEB/wQFMAMBAf8wCgYIKoZIzj0EAwIDSAAwRQIh
ALM0oNzktFewxfj+XpBCEWBM9WWDkk+m4SpFcrhGbtulAiBExcWdMcTb7np9prfb
esXjFja+BismJNLM2IiKyDinhA==
-----END CERTIFICATE-----
"""
_TEST_CA_DER = ssl.PEM_cert_to_DER_cert(_TEST_CA_PEM)


def _good_openssl(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    """``openssl`` present, the chain verifies and the issuer is the expected root."""
    monkeypatch.setattr(fr.shutil, "which", lambda _name: "/usr/bin/openssl")
    return _openssl_results(monkeypatch, (0, "intermediate.pem: OK"), f"issuer=CN = {fr.ROOT_CN}\n")


def _partial_chain_default_context(monkeypatch: pytest.MonkeyPatch) -> None:
    """Simulate Python >= 3.13, whose default context sets ``VERIFY_X509_PARTIAL_CHAIN``."""
    original = ssl.create_default_context

    def with_partial_chain(*args: Any, **kwargs: Any) -> ssl.SSLContext:
        ctx = original(*args, **kwargs)
        ctx.verify_flags |= ssl.VERIFY_X509_PARTIAL_CHAIN
        return ctx

    monkeypatch.setattr(fr.ssl, "create_default_context", with_partial_chain)


def test_regulondb_context_keeps_verification_on_and_adds_the_pinned_ca(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    baseline = ssl.create_default_context().cert_store_stats()["x509_ca"]
    _pin(monkeypatch, _TEST_CA_DER)
    _good_openssl(monkeypatch)
    _partial_chain_default_context(monkeypatch)
    fetched: list[str] = []

    def fetch(url: str) -> bytes:
        fetched.append(url)
        return _TEST_CA_DER

    ctx = fr.regulondb_context(fetch=fetch)

    assert fetched == [fr.INTERMEDIATE_URL]
    assert ctx.verify_mode == ssl.CERT_REQUIRED
    assert ctx.check_hostname is True
    # the handshake must reach a trusted root, whatever the Python default (lead decision 6)
    assert not ctx.verify_flags & ssl.VERIFY_X509_PARTIAL_CHAIN
    assert ctx.cert_store_stats()["x509_ca"] == baseline + 1


def test_an_openssl_timeout_is_a_fetch_error_not_a_traceback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _pin(monkeypatch, _TEST_CA_DER)
    monkeypatch.setattr(fr.shutil, "which", lambda _name: "/usr/bin/openssl")

    def hang(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        raise subprocess.TimeoutExpired(argv, kwargs.get("timeout") or 60)

    monkeypatch.setattr(fr.subprocess, "run", hang)
    with pytest.raises(fr.FetchError):
        fr.verify_intermediate(_TEST_CA_DER)


def test_main_returns_1_when_openssl_times_out(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GIPPA_DATA_DIR", str(tmp_path))
    _pin(monkeypatch, _TEST_CA_DER)
    monkeypatch.setattr(fr.shutil, "which", lambda _name: "/usr/bin/openssl")

    def hang(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        raise subprocess.TimeoutExpired(argv, 60)

    def intermediate_only(req: Any, *_: object, **__: object) -> _Response:
        assert getattr(req, "full_url", req) == fr.INTERMEDIATE_URL, "no other request expected"
        return _Response(_TEST_CA_DER)

    monkeypatch.setattr(fr.subprocess, "run", hang)
    monkeypatch.setattr(fr.urllib.request, "urlopen", intermediate_only)

    assert fr.main() == 1
    assert fr.HOST in capsys.readouterr().err
    assert not (tmp_path / "regulondb").exists()


# ------------------------------------------------- lead decision 7: retry exhaustion


def test_main_returns_1_and_writes_nothing_when_every_attempt_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GIPPA_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(fr, "regulondb_context", lambda: ssl.create_default_context())
    calls: list[str] = []

    def down(req: Any, *_: object, **__: object) -> _Response:
        calls.append(req.full_url)
        raise OSError("connection refused")

    sleeps: list[float] = []
    monkeypatch.setattr(fr.urllib.request, "urlopen", down)
    monkeypatch.setattr(fr.time, "sleep", sleeps.append)

    assert fr.main() == 1

    assert len(calls) == fr.ATTEMPTS  # one failing request per attempt, then give up
    assert len(sleeps) == fr.ATTEMPTS - 1  # no sleep after the last attempt
    assert not (tmp_path / "regulondb").exists()
    assert list(tmp_path.iterdir()) == []
    err = capsys.readouterr().err
    assert f"attempt {fr.ATTEMPTS}" in err
