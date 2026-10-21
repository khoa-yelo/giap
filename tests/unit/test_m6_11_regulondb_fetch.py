"""``validation/fetch_regulondb.py`` without the network (M6-11; GIAP_SPEC.md §3.9, §15.3).

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
        fr.verify_intermediate(_self_signed_der(), "-----BEGIN CERTIFICATE-----\n")
    message = str(exc.value)
    assert fr.INTERMEDIATE_SHA256 in message  # the pin is named, so an expiry is obvious
    assert hashlib.sha256(_self_signed_der()).hexdigest() in message


def test_regulondb_context_refuses_pem_from_the_ca_endpoint() -> None:
    with pytest.raises(fr.FetchError, match="PEM, not the DER"):
        fr.regulondb_context(fetch=lambda _url: b"-----BEGIN CERTIFICATE-----\n")


def test_regulondb_context_refuses_a_certificate_that_fails_the_pin() -> None:
    with pytest.raises(fr.FetchError, match="sha256"):
        fr.regulondb_context(fetch=lambda _url: _self_signed_der())


@pytest.mark.skipif(shutil.which("openssl") is None, reason="openssl is not on PATH")
def test_verify_intermediate_without_openssl_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    der = _self_signed_der()
    monkeypatch.setattr(fr, "INTERMEDIATE_SHA256", hashlib.sha256(der).hexdigest())
    monkeypatch.setattr(fr.shutil, "which", lambda _name: None)
    with pytest.raises(fr.FetchError, match="openssl is not on PATH"):
        fr.verify_intermediate(der, "-----BEGIN CERTIFICATE-----\n")


def test_trusted_root_pem_raises_when_the_root_is_absent() -> None:
    with pytest.raises(fr.FetchError, match="not in the system trust store"):
        fr.trusted_root_pem("No Such Root CA 1999")


def test_trusted_root_pem_returns_pem_for_a_root_that_is_there() -> None:
    ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    ctx.load_default_certs(ssl.Purpose.SERVER_AUTH)
    names = [
        value
        for cert in ctx.get_ca_certs()
        for rdn in cert.get("subject", ())
        for key, value in rdn
        if key == "commonName"
    ]
    if not names:
        pytest.skip("no system trust store on this host")
    pem = fr.trusted_root_pem(names[0])
    assert pem.startswith("-----BEGIN CERTIFICATE-----")
    assert ssl.PEM_cert_to_DER_cert(pem)


@pytest.mark.skipif(shutil.which("openssl") is None, reason="openssl is not on PATH")
def test_verify_intermediate_runs_openssl_against_the_root(monkeypatch: pytest.MonkeyPatch) -> None:
    """The pinned bytes must still verify against the root; the check is a real subprocess call."""
    der = _self_signed_der()
    monkeypatch.setattr(fr, "INTERMEDIATE_SHA256", hashlib.sha256(der).hexdigest())
    seen: dict[str, Any] = {}

    def fake_run(argv: list[str], **_: object) -> subprocess.CompletedProcess[str]:
        seen["argv"] = argv
        return subprocess.CompletedProcess(argv, 1, "", "error 20 at 0 depth lookup")

    monkeypatch.setattr(fr.subprocess, "run", fake_run)
    with pytest.raises(fr.FetchError, match="does not verify"):
        fr.verify_intermediate(der, "-----BEGIN CERTIFICATE-----\nAAA\n-----END CERTIFICATE-----\n")
    assert seen["argv"][1:3] == ["verify", "-CAfile"]


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
    monkeypatch.setenv("GIAP_DATA_DIR", str(tmp_path))
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
    monkeypatch.setenv("GIAP_DATA_DIR", str(tmp_path))
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


def test_main_fails_without_giap_data_dir(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.delenv("GIAP_DATA_DIR", raising=False)
    assert fr.main() == 1
    assert "GIAP_DATA_DIR is not set" in capsys.readouterr().err


def test_main_reports_a_broken_chain_instead_of_downgrading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setenv("GIAP_DATA_DIR", str(tmp_path))

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
    """R6/R11 and the §15.3 row: the data lives in ``$GIAP_DATA_DIR``, never in the repo."""
    names = {"TUSet.txt", "TerminatorSet.txt", "GeneProductSet.txt", "OperonSet.txt"}
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True
    ).stdout.splitlines()
    offenders = [p for p in tracked if Path(p).name in names or Path(p).parts[:1] == ("regulondb",)]
    assert offenders == []
    assert not any(p.startswith("validation/regulondb/") for p in tracked)
