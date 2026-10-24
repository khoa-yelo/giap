"""``scripts/slurm/fetch_pxd078459.sbatch``: resolve and download against stubs (M10-01).

Spec: docs/specs/M10-01-pxd078459-download.md, "Two jobs, in order", "Download rules", examples
5-8, edge cases, acceptance criterion 3; lead decisions 1 and 4 (https on the same host when HEAD
gives the record's byte size; names and URLs from the official record only, R6).

Every run happens in ``helpers.pxd_fixtures.pxd_sandbox``: a throw-away checkout with a fake
``curl`` that serves registered URLs from files, a fake ``df``, a ``wget`` that refuses, fake
Slurm tools and an unreachable proxy. **Nothing reaches the network or the real data.**

Conventions settled by the test author:

* ``--resolve`` saves the record to ``$GIAP_DATA_DIR/pxd078459/record/PXD078459.xml``, prints (or
  writes under ``record/``) its sha256 and the UTC date, and writes the draft manifest as a
  ``.tsv`` under ``record/``. It never GETs a raw file; HEAD requests are allowed.
* An ``ftp://`` record URI becomes ``url = https://<same host and path>`` in the draft only when a
  HEAD on that URL answers 200 with ``Content-Length`` equal to the record's ``bytes`` (lead
  decision 1); otherwise the resolve exits 2 or the row keeps a non-https ``url`` (which the
  download then refuses).
* The download reads ``validation/pilot/inputs/pxd078459_files.tsv`` of the checkout
  (``PXD078459_FILES``/``PXD078459_BASE_URL`` win when set), and writes the filled manifest
  (``sha256`` column) as a ``.tsv`` with the manifest header directly under the data directory, for
  the agent's second manifest commit (Q5). The repository is never written.
* Unsafe archives (absolute or ``..`` members, symlinks leaving the run folder, an unpacked size
  larger than the free space) are refused with exit 1 (2 accepted for the size case) and the
  message ``unsafe archive member: <member>`` for the path cases.

No statistical claim is made in this module, so it has no calibration test (spec §8.4).
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import pytest
from helpers.pxd_fixtures import (
    BASE,
    GB,
    IPROX,
    MANIFEST_COLUMNS,
    RECORD_URL,
    FileSpec,
    PxdSandbox,
    example1_files,
    make_run,
    out,
    px_xml,
    pxd_sandbox,
    read_tsv,
    sha256_hex,
    tar_run,
    unsafe_tar,
    unsafe_zip,
    zip_bomb,
    zip_run,
)

OFFICIAL_HOSTS = {"proteomecentral.proteomexchange.org", "www.iprox.cn", IPROX}
SMALL_RUN = {"msms": {0: 3, 9: 9}, "dia_groups": 3}


@pytest.fixture
def px(tmp_path: Path) -> PxdSandbox:
    return pxd_sandbox(tmp_path)


def _build(px: PxdSandbox, run: str, *, ext: str = ".d.zip", folder: str | None = None) -> bytes:
    src = px.sb.tmp / "build" / f"{run}.d"
    if not src.exists():
        make_run(src, **SMALL_RUN)
    archive = px.sb.tmp / "build" / f"{run}{ext}"
    folder = folder or f"{run}.d"
    if ext == ".d.zip":
        return zip_run(archive, folder, src)
    return tar_run(archive, folder, src, gz=ext.endswith(".gz"))


def _row(run: str, name: str, body: bytes, **over: object) -> dict[str, object]:
    row: dict[str, object] = {
        "run": run,
        "file_name": name,
        "url": f"{BASE}/{name}",
        "url_record": f"{BASE}/{name}",
        "bytes": len(body),
        "checksum_record": f"md5:{hashlib.md5(body).hexdigest()}",
        "sha256": "NA",
        "sample": "WT",
        "selected_because": "token:WT",
    }
    row.update(over)
    return row


def _three(
    px: PxdSandbox, *, ext: str = ".d.zip", serve: bool = True, **over: object
) -> tuple[list[dict[str, object]], dict[str, bytes]]:
    rows, bodies = [], {}
    for i in (1, 2, 3):
        run = f"WT_{i}"
        name = f"{run}{ext}"
        body = _build(px, run, ext=ext)
        bodies[name] = body
        rows.append(_row(run, name, body, **over))
        if serve:
            px.serve(f"{BASE}/{name}", body)
    px.write_manifest(rows)
    return rows, bodies


def _sums(px: PxdSandbox) -> dict[str, str]:
    path = px.dest / "SHA256SUMS"
    assert path.is_file(), f"{path} was not written"
    listed = {}
    for line in path.read_text().splitlines():
        m = re.match(r"^([0-9a-f]{64})\s+\*?(\S+)$", line.strip())
        assert m, f"unparsable SHA256SUMS line: {line!r}"
        listed[Path(m.group(2)).name] = m.group(1)
    return listed


def _get_counts(px: PxdSandbox) -> Counter[str]:
    return Counter(c["url"] for c in px.gets())


# ------------------------------------------------------------------ download: happy path


def test_downloads_exactly_the_manifest_urls(px: PxdSandbox) -> None:
    rows, _ = _three(px)
    res = px.run()
    assert res.returncode == 0, out(res)
    assert _get_counts(px) == Counter(str(r["url"]) for r in rows), px.gets()
    assert not [c for c in px.calls() if c.get("tool") == "wget"]
    assert not any(c["insecure"] for c in px.calls()), "R6: TLS verification was disabled"
    assert px.repo_is_clean() == "", px.repo_is_clean()


def test_sha256sums_lists_every_archive(px: PxdSandbox) -> None:
    _, bodies = _three(px)
    res = px.run()
    assert res.returncode == 0, out(res)
    assert _sums(px) == {name: sha256_hex(b) for name, b in bodies.items()}


def test_runs_are_unpacked_to_run_d_and_archives_deleted(px: PxdSandbox) -> None:
    _, bodies = _three(px)
    res = px.run()
    assert res.returncode == 0, out(res)
    for i in (1, 2, 3):
        run_d = px.dest / f"WT_{i}.d"
        assert (run_d / "analysis.tdf").is_file(), sorted(p.name for p in px.dest.iterdir())
        assert (run_d / "analysis.tdf_bin").is_file()
    for name in bodies:
        assert not (px.dest / name).exists(), f"{name} kept after a successful unpack (§3.9)"


def test_filled_manifest_is_written_next_to_the_data(px: PxdSandbox) -> None:
    """Q5 / deliverable 3: the second commit copies a manifest with ``sha256`` filled."""
    _, bodies = _three(px)
    res = px.run()
    assert res.returncode == 0, out(res)
    filled = []
    for p in px.dest.glob("*.tsv"):
        header, rows = read_tsv(p)
        if tuple(header) == MANIFEST_COLUMNS:
            filled.append(rows)
    assert len(filled) == 1, sorted(p.name for p in px.dest.iterdir())
    got = {r["file_name"]: r["sha256"] for r in filled[0]}
    assert got == {name: sha256_hex(b) for name, b in bodies.items()}


@pytest.mark.parametrize("ext", [".tar", ".tar.gz"])
def test_tar_archives_are_unpacked(px: PxdSandbox, ext: str) -> None:
    _three(px, ext=ext)
    res = px.run()
    assert res.returncode == 0, out(res)
    for i in (1, 2, 3):
        assert (px.dest / f"WT_{i}.d" / "analysis.tdf").is_file()


def test_nested_folder_is_normalised_to_run_d(px: PxdSandbox) -> None:
    """Edge case: ``WT_1/20260101_EcN_WT_1.d/`` in the archive ends up at ``WT_1.d/``."""
    rows, _ = _three(px)
    body = _build(px, "WT_1", folder="WT_1/20260101_EcN_WT_1.d")
    rows[0] = _row("WT_1", "WT_1.d.zip", body)
    px.serve(f"{BASE}/WT_1.d.zip", body)
    px.write_manifest(rows)
    res = px.run()
    assert res.returncode == 0, out(res)
    tdfs = list((px.dest / "WT_1.d").rglob("analysis.tdf"))
    assert tdfs == [px.dest / "WT_1.d" / "analysis.tdf"], tdfs


@pytest.mark.parametrize("ext", [".d.7z", ".d.rar"])
def test_unsupported_archive_format_exits_2(px: PxdSandbox, ext: str) -> None:
    rows, _ = _three(px)
    rows[1] = _row("WT_2", f"WT_2{ext}", b"7z\xbc\xaf\x27\x1c" + b"\x00" * 64)
    px.serve(str(rows[1]["url"]), b"7z\xbc\xaf\x27\x1c" + b"\x00" * 64)
    px.write_manifest(rows)
    res = px.run()
    assert res.returncode == 2, out(res)


# ------------------------------------------------------------------ scheme (R6, decision 11)


@pytest.mark.parametrize("scheme", ["ftp", "http"])
def test_plain_ftp_or_http_row_exits_2_before_any_download(px: PxdSandbox, scheme: str) -> None:
    rows, _ = _three(px)
    bad = f"{scheme}://{IPROX}/IPX0012345000/IPX0012345001/WT_2.d.zip"
    rows[1]["url"] = bad
    rows[1]["url_record"] = bad
    px.write_manifest(rows)
    res = px.run()
    assert res.returncode == 2, out(res)
    msg = res.stdout + res.stderr
    assert "WT_2.d.zip" in msg and "decision 11" in msg, msg
    assert not px.gets(), f"downloads happened before the scheme check: {px.gets()}"


def test_ftp_record_with_https_url_downloads_https(px: PxdSandbox) -> None:
    """Q1 / lead decision 1: ``url_record`` may be ftp:// when ``url`` is the same-host https."""
    rows, _ = _three(px)
    for r in rows:
        r["url_record"] = str(r["url"]).replace("https://", "ftp://", 1)
    px.write_manifest(rows)
    res = px.run()
    assert res.returncode == 0, out(res)
    assert all(c["url"].startswith("https://") for c in px.gets()), px.gets()


# ------------------------------------------------------------------ manifest precedence


def test_no_manifest_and_no_override_exits_2(px: PxdSandbox) -> None:
    px.remove_manifest()
    res = px.run()
    assert res.returncode == 2, out(res)
    assert "no file list: run --resolve first" in res.stdout + res.stderr
    assert not px.gets()


def test_env_overrides_win_over_the_manifest(px: PxdSandbox) -> None:
    _three(px)
    other = f"https://{IPROX}/override"
    names = ["A_1.d.zip", "A_2.d.zip", "A_3.d.zip"]
    for i, n in enumerate(names, start=1):
        px.serve(f"{other}/{n}", _build(px, f"WT_{i}"))
    px.run(PXD078459_BASE_URL=other, PXD078459_FILES=" ".join(names))
    assert set(_get_counts(px)) == {f"{other}/{n}" for n in names}, px.gets()


# ------------------------------------------------------------------ verification and retry


def test_example6_checksum_mismatch_twice_exits_1(px: PxdSandbox) -> None:
    _, bodies = _three(px)
    px.serve(f"{BASE}/WT_2.d.zip", b"corrupt one", b"corrupt two")
    res = px.run()
    assert res.returncode == 1, out(res)
    assert "checksum mismatch after one retry: WT_2.d.zip" in res.stdout + res.stderr
    assert _get_counts(px)[f"{BASE}/WT_2.d.zip"] == 2, "exactly one retry (R6)"
    sums = _sums(px)
    assert set(sums) == {"WT_1.d.zip", "WT_3.d.zip"}, sums
    assert sums["WT_1.d.zip"] == sha256_hex(bodies["WT_1.d.zip"])
    assert not (px.dest / "WT_2.d.zip").exists(), "a mismatching file is deleted"


def test_checksum_mismatch_once_then_ok(px: PxdSandbox) -> None:
    _, bodies = _three(px)
    px.serve(f"{BASE}/WT_2.d.zip", b"corrupt once", bodies["WT_2.d.zip"])
    res = px.run()
    assert res.returncode == 0, out(res)
    assert _get_counts(px)[f"{BASE}/WT_2.d.zip"] == 2
    assert _sums(px)["WT_2.d.zip"] == sha256_hex(bodies["WT_2.d.zip"])


def test_example7_partial_file_is_resumed(px: PxdSandbox) -> None:
    _, bodies = _three(px)
    px.dest.mkdir(parents=True, exist_ok=True)
    (px.dest / "WT_1.d.zip").write_bytes(bodies["WT_1.d.zip"])
    (px.dest / "WT_2.d.zip").write_bytes(bodies["WT_2.d.zip"])
    half = len(bodies["WT_3.d.zip"]) // 2
    (px.dest / "WT_3.d.zip").write_bytes(bodies["WT_3.d.zip"][:half])
    res = px.run()
    assert res.returncode == 0, out(res)
    assert "already there" in res.stdout + res.stderr
    gets = px.gets()
    assert [c["url"] for c in gets] == [f"{BASE}/WT_3.d.zip"], gets
    assert gets[0]["continue_at"] == "-", f"no --continue-at - in {gets[0]['argv']}"
    assert gets[0]["resumed_from"] == half
    assert _sums(px) == {n: sha256_hex(b) for n, b in bodies.items()}


def test_unpacked_run_with_recorded_sha256_is_skipped(px: PxdSandbox) -> None:
    rows, bodies = _three(px)
    rows[0]["sha256"] = sha256_hex(bodies["WT_1.d.zip"])
    px.write_manifest(rows)
    make_run(px.dest / "WT_1.d", **SMALL_RUN)
    res = px.run()
    assert res.returncode == 0, out(res)
    assert f"{BASE}/WT_1.d.zip" not in _get_counts(px), px.gets()


def test_unpacked_run_without_recorded_sha256_is_downloaded_again(px: PxdSandbox) -> None:
    _three(px)
    make_run(px.dest / "WT_1.d", **SMALL_RUN)
    res = px.run()
    assert res.returncode == 0, out(res)
    assert _get_counts(px)[f"{BASE}/WT_1.d.zip"] == 1


def test_login_page_instead_of_an_archive_exits_1(px: PxdSandbox) -> None:
    """Edge case: HTML where an archive was expected (no record checksum to catch it)."""
    rows, _ = _three(px)
    page = b"<!DOCTYPE html><html><body>Please log in</body></html>\n"
    rows[0].update(bytes="NA", checksum_record="NA")
    px.write_manifest(rows)
    px.serve(f"{BASE}/WT_1.d.zip", page, content_type="text/html")
    res = px.run()
    assert res.returncode == 1, out(res)


def test_http_403_exits_1(px: PxdSandbox) -> None:
    _three(px)
    px.serve(f"{BASE}/WT_1.d.zip", b"forbidden", status=403, content_type="text/html")
    res = px.run()
    assert res.returncode == 1, out(res)
    sums = px.dest / "SHA256SUMS"
    assert not sums.exists() or "WT_1.d.zip" not in _sums(px), sums.read_text()


# ------------------------------------------------------------------ unsafe archives


def _escaped(px: PxdSandbox) -> list[Path]:
    """Files the hostile archives try to plant, wherever they landed (payloads excluded)."""
    found = []
    for p in px.sb.tmp.rglob("*"):
        if px.payloads in p.parents or (px.sb.tmp / "build") in p.parents:
            continue
        if p.name in {"evil.sh", "abs_evil.sh", "escaped.txt"}:
            found.append(p)
    return found


def _links(px: PxdSandbox) -> list[Path]:
    return [p for p in px.dest.rglob("*") if p.is_symlink()] if px.dest.exists() else []


@pytest.mark.parametrize("kind", ["dotdot", "absolute", "zip_symlink", "tar_symlink"])
def test_example8_unsafe_archive_member_exits_1(px: PxdSandbox, kind: str) -> None:
    rows, _ = _three(px)
    archive = px.sb.tmp / "build" / "hostile"
    if kind == "dotdot":
        member, name = "../../evil.sh", "WT_1.d.zip"
        body = unsafe_zip(archive, member)
    elif kind == "absolute":
        member, name = str(px.sb.tmp / "abs_evil.sh"), "WT_1.d.zip"
        body = unsafe_zip(archive, member)
    elif kind == "zip_symlink":
        member, name = "WT_1.d/link", "WT_1.d.zip"
        body = unsafe_zip(archive, member, symlink_to="../../..")
    else:
        member, name = "WT_1.d/link", "WT_1.d.tar"
        body = unsafe_tar(archive, symlink_to="../../..")
    rows[0] = _row("WT_1", name, body)
    px.serve(str(rows[0]["url"]), body)
    px.write_manifest(rows)
    res = px.run()
    assert res.returncode == 1, out(res)
    msg = res.stdout + res.stderr
    assert "unsafe archive member" in msg, msg
    if kind == "dotdot":
        assert "unsafe archive member: ../../evil.sh" in msg, msg
    assert not _escaped(px), _escaped(px)
    assert not _links(px), f"symlinks left in the data directory: {_links(px)}"
    assert not (px.dest / "WT_1.d" / "analysis.tdf").exists(), "nothing is extracted"


def test_archive_larger_than_free_space_is_refused(px: PxdSandbox) -> None:
    """Size limit: the 3x disk margin covers the record's bytes; an archive inflating beyond the
    free space (a zip bomb, or a wrong ``bytes``) is refused before it fills the disk."""
    rows, _ = _three(px)
    bomb = zip_bomb(px.sb.tmp / "build" / "bomb.zip", unpacked_bytes=64 * 1024 * 1024)
    rows[0] = _row("WT_1", "WT_1.d.zip", bomb)
    px.serve(f"{BASE}/WT_1.d.zip", bomb)
    px.write_manifest(rows)
    need = 3 * sum(int(str(r["bytes"])) for r in rows)
    res = px.run(FAKE_DF_AVAIL_BYTES=str(need + 8 * 1024 * 1024))
    assert res.returncode in (1, 2), out(res)
    big = px.dest / "WT_1.d" / "analysis.tdf_bin"
    assert not big.exists() or big.stat().st_size < 64 * 1024 * 1024, "the bomb was inflated"


# ------------------------------------------------------------------ disk


def test_too_little_disk_exits_2_before_the_first_byte(px: PxdSandbox) -> None:
    rows, _ = _three(px)
    for r in rows:
        r["bytes"] = 2 * GB
    px.write_manifest(rows)
    res = px.run(FAKE_DF_AVAIL_BYTES=str(17 * GB))  # needs 3 x 6 GB = 18 GB
    assert res.returncode == 2, out(res)
    assert not px.gets(), px.gets()
    assert px.df_calls(), "free space is read with df (acceptance 3: stubbed df)"


def test_na_bytes_assume_5_gb_per_run(px: PxdSandbox) -> None:
    _three(px, bytes="NA")
    res = px.run(FAKE_DF_AVAIL_BYTES=str(40 * GB))  # needs 3 x 15 GB = 45 GB
    assert res.returncode == 2, out(res)
    assert not px.gets()


def test_na_bytes_with_enough_disk_downloads(px: PxdSandbox) -> None:
    _three(px, bytes="NA")
    res = px.run(FAKE_DF_AVAIL_BYTES=str(60 * 1024**3))
    assert res.returncode == 0, out(res)
    assert len(px.gets()) == 3


# ------------------------------------------------------------------ --resolve


def _resolve_record(px: PxdSandbox, files: list[FileSpec], **kw: Any) -> bytes:
    xml = px_xml(files, **kw)
    px.serve(RECORD_URL, xml, content_type="application/xml")
    return xml


def _draft(px: PxdSandbox) -> list[dict[str, str]]:
    found = []
    record = px.dest / "record"
    for p in record.glob("*.tsv") if record.exists() else []:
        header, rows = read_tsv(p)
        if tuple(header) == MANIFEST_COLUMNS:
            found.append(rows)
    assert len(found) == 1, f"expected one draft manifest under {record}"
    return found[0]


def _no_raw_get(px: PxdSandbox) -> None:
    raw = [c["url"] for c in px.gets() if c["url"] != RECORD_URL]
    assert not raw, f"--resolve downloaded raw data: {raw}"


def test_resolve_saves_the_record_and_writes_a_draft(px: PxdSandbox) -> None:
    px.remove_manifest()
    xml = _resolve_record(px, example1_files())
    res = px.run("--resolve")
    assert res.returncode == 0, out(res)
    assert RECORD_URL in _get_counts(px), px.gets()
    saved = px.dest / "record" / "PXD078459.xml"
    assert saved.read_bytes() == xml
    digest = hashlib.sha256(xml).hexdigest()
    texts = [res.stdout, res.stderr]
    texts += [p.read_text(errors="replace") for p in saved.parent.iterdir() if p != saved]
    assert any(digest in t for t in texts), "the record's sha256 is not recorded"
    assert any(re.search(r"20\d\d-\d\d-\d\d", t) for t in texts), "no UTC retrieval date"
    rows = _draft(px)
    assert [r["run"] for r in rows] == ["WT_1", "WT_2", "WT_3"]
    assert [r["file_name"] for r in rows] == ["WT_1.d.zip", "WT_2.d.zip", "WT_3.d.zip"]
    _no_raw_get(px)
    assert px.repo_is_clean() == "", "the draft is written outside the checkout"


def _ftp_files() -> list[FileSpec]:
    return [
        FileSpec(f.name, f.uri.replace("https://", "ftp://", 1), f.raw, f.size, f.md5)
        for f in example1_files()
    ]


def test_resolve_uses_https_on_the_same_host_when_head_matches(px: PxdSandbox) -> None:
    px.remove_manifest()
    files = _ftp_files()
    _resolve_record(px, files)
    for f in files:
        px.serve(f.uri.replace("ftp://", "https://", 1), None, head={"content_length": f.size})
    res = px.run("--resolve")
    assert res.returncode == 0, out(res)
    rows = {r["file_name"]: r for r in _draft(px)}
    for f in files[:3]:
        assert rows[f.name]["url_record"] == f.uri
        assert rows[f.name]["url"] == f.uri.replace("ftp://", "https://", 1)
    heads = {c["url"] for c in px.heads()}
    assert heads <= {f.uri.replace("ftp://", "https://", 1) for f in files}, heads
    _no_raw_get(px)


@pytest.mark.parametrize("head", [{"content_length": 12345}, {"status": 404}])
def test_resolve_without_a_matching_head_keeps_no_https(
    px: PxdSandbox, head: dict[str, int]
) -> None:
    px.remove_manifest()
    files = _ftp_files()
    _resolve_record(px, files)
    for f in files:
        h = head if f.name == "WT_2.d.zip" else {"content_length": f.size}
        px.serve(f.uri.replace("ftp://", "https://", 1), None, head=h)
    res = px.run("--resolve")
    _no_raw_get(px)
    if res.returncode == 0:
        rows = {r["file_name"]: r for r in _draft(px)}
        assert not rows["WT_2.d.zip"]["url"].startswith("https://"), rows["WT_2.d.zip"]
    else:
        assert res.returncode == 2, out(res)
        assert "WT_2.d.zip" in res.stdout + res.stderr


def test_resolve_fewer_than_three_exits_2_without_a_draft(px: PxdSandbox) -> None:
    px.remove_manifest()
    _resolve_record(px, example1_files()[:2] + example1_files()[3:])
    res = px.run("--resolve")
    assert res.returncode == 2, out(res)
    assert "fewer than three wild-type runs" in res.stdout + res.stderr
    record = px.dest / "record"
    assert not (record.exists() and list(record.glob("*.tsv"))), "example 4: no manifest"
    _no_raw_get(px)


def test_resolve_refuses_a_non_timstof_record(px: PxdSandbox) -> None:
    px.remove_manifest()
    _resolve_record(px, example1_files(), instrument="Q Exactive HF")
    res = px.run("--resolve")
    assert res.returncode == 2, out(res)
    assert "does not describe timsTOF data" in res.stdout + res.stderr


def test_resolve_with_the_record_unreachable_fails(px: PxdSandbox) -> None:
    px.remove_manifest()
    px.serve(RECORD_URL, b"down", b"down", status=503, content_type="text/html")
    res = px.run("--resolve")
    assert res.returncode != 0, out(res)
    record = px.dest / "record"
    assert not (record.exists() and list(record.glob("*.tsv")))


def test_only_official_hosts_are_contacted(px: PxdSandbox) -> None:
    """R6: the record, iProX and the host the record names; never a mirror."""
    px.remove_manifest()
    files = _ftp_files()
    _resolve_record(px, files)
    for f in files:
        px.serve(f.uri.replace("ftp://", "https://", 1), None, head={"content_length": f.size})
    px.run("--resolve")
    _three(px)
    px.run()
    hosts = {urlsplit(c["url"]).hostname for c in px.calls() if c.get("tool") == "curl"}
    assert hosts and hosts <= OFFICIAL_HOSTS, hosts
