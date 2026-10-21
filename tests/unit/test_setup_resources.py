"""`giap setup --go --kegg --cog --contaminants` (issue M2-09; GIAP_SPEC.md §7.9.1, §3.2, §7.2.8).

Examples 1-6 and the edge cases of docs/specs/M2-09-giap-setup-resources.md with fake fetchers.
Every test runs under ``setup_guard`` (M2-08 lead decision 3): no network, no tool, no real
$GIAP_DB_DIR (``GIAP_DB_DIR`` and the tool roots point into ``tmp_path``).

Ambiguities resolved here:
- The GO pin is a release date: ``SOURCES["go_basic"].version == "YYYY-MM-DD"`` and the URL is
  ``.../obo/go/releases/<version>/...``; ``verify_go(path, release)`` takes that date. Tests re-pin
  the release to the spec's ``2026-07-01`` (and the checksum to the fake bytes) so they do not
  depend on the release chosen in ``db_sources.yaml``.
- The live ``goslim_generic.obo`` header is ``data-version: go/releases/<date>/subsets/
  goslim_generic.owl`` (checked against release.geneontology.org, 2025-07-22), not
  ``releases/<date>``; the spec's "same" is read as "same release date", so both forms verify.
- ``purl.obolibrary.org`` answers with a 302 to ``https://release.geneontology.org/...``; §3.2 lists
  ``*.geneontology.org`` as an official GO source, so an https redirect from purl to a
  ``*.geneontology.org`` host is followed. Every other cross-host or non-https redirect is refused.
- KEGG items carry no pinned checksum (``""``); the record's checksum is the sha256 of the bytes
  retrieved and ``version`` is the retrieval date. ``--check`` compares the file with that
  recorded checksum.
- A GO version mismatch is reported on one stderr line ``go_basic: ... version mismatch (expected
  releases/2026-07-01)``; the record is never ``ok`` and the file is not installed.
- antiSMASH (backlog): whether the pinned image already holds the databases is asked with
  ``antismash --check-prereqs`` through the command runner (M2-08 spec "verify rules"); exit 0
  without ``--databases`` means "in the image": record ``ok`` with ``url: "image"``, no download.
- Backlog items already met by the M2-08 code (non-200, redirect cap, tarball deletion, unsafe tar
  members reported by GIAP's own check rather than by ``tarfile``'s ``data`` filter) are pinned
  here so they stay met.
- The live download test of the acceptance list is left out, as in M2-08: tests never download.
"""

from __future__ import annotations

import dataclasses
import hashlib
import http.client
import io
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
from collections.abc import Sequence
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any, NamedTuple

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from setup_guard import BLOCKED, real_https_fetch, setup_guard  # noqa: F401
from test_setupdb import GENOMAD_TAR, NoFetcher, NoRunner, _tar, make_genomad, md5

GIB = 1024**3
REPO = Path(__file__).resolve().parents[2]
RELEASE = "2026-07-01"
OTHER_RELEASE = "2026-09-01"
KO_URL = "https://rest.kegg.jp/link/pathway/ko"
PATHWAY_URL = "https://rest.kegg.jp/list/pathway"
COG_DEF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab"
COG_FUN_URL = "https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab"
GO_NAMES = ("go_basic", "goslim_generic")
KEGG_NAMES = ("kegg_ko_pathway", "kegg_pathway_names")
COG_NAMES = ("cog2020_def", "cog2020_fun")
NEW_NAMES = (*GO_NAMES, *KEGG_NAMES, *COG_NAMES, "contaminants")
SHA256 = re.compile(r"sha256:[0-9a-f]{64}")


class _Usage(NamedTuple):
    total: int
    used: int
    free: int


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Plenty of disk, inside Slurm, every GIAP root (and $GIAP_DB_DIR) in ``tmp_path``."""
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: _Usage(10_000 * GIB, 0, 10_000 * GIB))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    roots = {"GIAP_ROOT": "root", "MAMBA_ROOT_PREFIX": "mamba",
             "NXF_APPTAINER_CACHEDIR": "cache", "GIAP_DB_DIR": "not-the-real-db"}  # fmt: skip
    for var, sub in roots.items():
        monkeypatch.setenv(var, str(tmp_path / "_env" / sub))
    monkeypatch.setenv("GIAP_REPO", str(REPO))


# ---------------------------------------------------------------- payloads


def go_obo(release: str, slim: bool = False) -> bytes:
    """A small OBO file with the header form of the live release files."""
    version = f"go/releases/{release}/subsets/goslim_generic.owl" if slim else f"releases/{release}"
    return (
        f"format-version: 1.2\ndata-version: {version}\n"
        'subsetdef: goslim_generic "Generic GO slim"\nontology: go\n\n'
        "[Term]\nid: GO:0008150\nname: biological_process\nnamespace: biological_process\n"
    ).encode()


KO_BODY = b"ko:K00001\tpath:map00010\n"  # example 3
PATHWAY_BODY = b"map00010\tGlycolysis / Gluconeogenesis\nmap01100\tMetabolic pathways\n"
COG_DEF_BODY = b"COG0001\tH\tGlutamate-1-semialdehyde aminotransferase\tHemL\t\t\t1BS0\n"
COG_FUN_BODY = b"J\tFCCCFC\tTranslation, ribosomal structure and biogenesis\n"


def fasta(n: int) -> bytes:
    return b"".join(
        f">Cont_P{i:05d} contaminant {i}\nMKVLAAGIVG\nLLAA\n".encode() for i in range(n)
    )


# ---------------------------------------------------------------- fakes


def sha256(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def sources() -> Any:
    from giap.setupdb import SOURCES

    return SOURCES


def repin(monkeypatch: pytest.MonkeyPatch, name: str, **fields: Any) -> Any:
    """Replace ``SOURCES[name]`` (every ``giap.setupdb`` module that holds it) by a copy."""
    item = dataclasses.replace(sources()[name], **fields)
    patched = {**sources(), name: item}
    for mod_name, mod in list(sys.modules.items()):
        if mod_name.startswith("giap.setupdb") and hasattr(mod, "SOURCES"):
            monkeypatch.setattr(mod, "SOURCES", patched)
    return item


def go_url(release: str, slim: bool) -> str:
    sub = "subsets/goslim_generic.obo" if slim else "go-basic.obo"
    return f"https://purl.obolibrary.org/obo/go/releases/{release}/{sub}"


def pin_go(monkeypatch: pytest.MonkeyPatch, served: str = RELEASE) -> dict[str, bytes]:
    """Pin GO to ``RELEASE`` with checksums of OBO files that announce ``served``."""
    out = {}
    for name, slim in (("go_basic", False), ("goslim_generic", True)):
        body = go_obo(served, slim)
        item = repin(monkeypatch, name, version=RELEASE, url=go_url(RELEASE, slim),
                     checksum=sha256(body))  # fmt: skip
        out[item.url] = body
    return out


def pin_cog(monkeypatch: pytest.MonkeyPatch) -> dict[str, bytes]:
    out = {}
    for name, body in (("cog2020_def", COG_DEF_BODY), ("cog2020_fun", COG_FUN_BODY)):
        out[repin(monkeypatch, name, checksum=sha256(body)).url] = body
    return out


def pin_contaminants(monkeypatch: pytest.MonkeyPatch, body: bytes) -> dict[str, bytes]:
    return {repin(monkeypatch, "contaminants", checksum=sha256(body)).url: body}


def kegg_routes() -> dict[str, bytes]:
    return {KO_URL: KO_BODY, PATHWAY_URL: PATHWAY_BODY}


def pin_all(monkeypatch: pytest.MonkeyPatch) -> dict[str, bytes]:
    return {**pin_go(monkeypatch), **kegg_routes(), **pin_cog(monkeypatch),
            **pin_contaminants(monkeypatch, fasta(381))}  # fmt: skip


class MapFetcher:
    """Serves bytes by URL; an exception value is raised after half the bytes were written."""

    def __init__(self, routes: dict[str, bytes | Exception]) -> None:
        self.routes = routes
        self.urls: list[str] = []
        self.dests: list[Path] = []

    def fetch(self, url: str, dest: Path) -> None:
        self.urls.append(url)
        self.dests.append(dest)
        value = self.routes.get(url)
        if value is None:
            raise AssertionError(f"unexpected fetch of {url}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(value, Exception):
            dest.write_bytes(b"half a download")
            raise value
        dest.write_bytes(value)


def source_error(host: str, reason: str = "timed out") -> Exception:
    from giap.setupdb import SourceError

    return SourceError(host, reason)


# ---------------------------------------------------------------- helpers


def opts(db: Path, **kw: Any) -> Any:
    from giap.setupdb import SetupOptions

    base: dict[str, Any] = {"db_dir": db, "bakta": "skip", "genomad": False, "antismash": False,
                            "promotech_models": False, "check": False}  # fmt: skip
    return SetupOptions(**(base | kw))


def setup(db: Path, fetcher: Any = None, runner: Any = None, **kw: Any) -> int:
    from giap.setupdb import run_setup

    return int(run_setup(opts(db, **kw), runner or NoRunner(), fetcher or NoFetcher()))


def items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {i["name"]: i for i in data["items"]}


def validate_manifest(db: Path) -> None:
    import jsonschema

    schema = json.loads((REPO / "src/giap/schemas/db_manifest.schema.json").read_text())
    jsonschema.validate(json.loads((db / "db_manifest.json").read_text()), schema)


def files_under(d: Path) -> list[Path]:
    return sorted(p for p in d.rglob("*") if p.is_file()) if d.exists() else []


def leftovers(db: Path) -> list[Path]:
    return [p for p in db.rglob("*") if ".partial" in p.name or p.name.endswith((".tmp", ".gz"))]


def today_dates() -> set[str]:
    return {datetime.now(UTC).date().isoformat(), date.today().isoformat()}


def item_file(db: Path, name: str) -> Path:
    return db / sources()[name].path


# ---------------------------------------------------------------- sources and plan


def test_new_items_are_the_official_sources() -> None:
    from giap.setupdb import DbItem

    s = sources()
    assert set(NEW_NAMES) <= set(s)
    for name in NEW_NAMES:
        assert isinstance(s[name], DbItem) and s[name].tool is None, s[name]
        assert s[name].url.startswith("https://"), s[name]
    go_v = s["go_basic"].version
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", go_v), go_v
    assert s["goslim_generic"].version == go_v, "both GO files come from one release"
    assert s["go_basic"].url == go_url(go_v, slim=False)
    assert s["goslim_generic"].url == go_url(go_v, slim=True)
    assert s["kegg_ko_pathway"].url == KO_URL and s["kegg_pathway_names"].url == PATHWAY_URL
    assert s["cog2020_def"].url == COG_DEF_URL and s["cog2020_fun"].url == COG_FUN_URL
    hosts = {n: s[n].host for n in NEW_NAMES}
    assert {hosts[n] for n in GO_NAMES} == {"purl.obolibrary.org"}
    assert {hosts[n] for n in KEGG_NAMES} == {"rest.kegg.jp"}
    assert {hosts[n] for n in COG_NAMES} == {"ftp.ncbi.nlm.nih.gov"}
    assert hosts["contaminants"] in {"raw.githubusercontent.com", "github.com"}
    assert re.search(r"/[0-9a-f]{40}/", s["contaminants"].url), "pinned by commit SHA"
    paths = {n: s[n].path for n in NEW_NAMES}
    assert paths == {
        "go_basic": "go/go-basic.obo",
        "goslim_generic": "go/goslim_generic.obo",
        "kegg_ko_pathway": "kegg/ko_pathway.tsv",
        "kegg_pathway_names": "kegg/pathway_names.tsv",
        "cog2020_def": "cog/cog-20.def.tab",
        "cog2020_fun": "cog/fun-20.tab",
        "contaminants": "contaminants/contaminants.fasta",
    }


def test_pinned_and_unpinned_checksums() -> None:
    s = sources()
    for name in (*GO_NAMES, *COG_NAMES, "contaminants"):
        assert SHA256.fullmatch(s[name].checksum), (name, s[name].checksum)
    for name in KEGG_NAMES:
        assert s[name].checksum == "", "KEGG REST data are unversioned: recorded, not pinned"


@pytest.mark.parametrize(
    ("flag", "names"),
    [
        ("go", GO_NAMES),
        ("kegg", KEGG_NAMES),
        ("cog", COG_NAMES),
        ("contaminants", ("contaminants",)),
    ],
)
def test_plan_includes_resources(tmp_path: Path, flag: str, names: tuple[str, ...]) -> None:
    from giap.setupdb import plan

    assert {i.name for i in plan(opts(tmp_path, **{flag: True}))} == set(names)


def test_plan_all_resources_and_none(tmp_path: Path) -> None:
    from giap.setupdb import plan

    every = {"go": True, "kegg": True, "cog": True, "contaminants": True}
    assert {i.name for i in plan(opts(tmp_path, **every))} == set(NEW_NAMES)
    assert plan(opts(tmp_path)) == []


def test_resource_flags_are_no_longer_ignored(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    fetcher = MapFetcher(pin_go(monkeypatch))
    setup(tmp_path, fetcher, go=True)
    assert "not handled" not in capsys.readouterr().err
    assert fetcher.urls, "--go fetched nothing"


# ---------------------------------------------------------------- verify_go / verify_kegg_links


def _write(tmp_path: Path, body: bytes, name: str = "f") -> Path:
    p = tmp_path / name
    p.write_bytes(body)
    return p


def test_verify_go_header_matches_release(tmp_path: Path) -> None:
    from giap.setupdb import verify_go

    assert verify_go(_write(tmp_path, go_obo(RELEASE)), RELEASE) is True
    assert verify_go(_write(tmp_path, go_obo(OTHER_RELEASE)), RELEASE) is False


def test_verify_go_accepts_the_goslim_header_form(tmp_path: Path) -> None:
    from giap.setupdb import verify_go

    assert verify_go(_write(tmp_path, go_obo(RELEASE, slim=True)), RELEASE) is True
    assert verify_go(_write(tmp_path, go_obo(OTHER_RELEASE, slim=True)), RELEASE) is False


@pytest.mark.parametrize(
    "body",
    [b"", b"<html><body>502 Bad Gateway</body></html>\n",
     b"format-version: 1.2\nontology: go\n\n[Term]\nid: GO:0008150\n",
     b"format-version: 1.2\ndata-version: releases/2026-07-011\n"],
)  # fmt: skip
def test_verify_go_rejects_files_without_the_release(tmp_path: Path, body: bytes) -> None:
    from giap.setupdb import verify_go

    assert verify_go(_write(tmp_path, body), RELEASE) is False


def test_verify_go_missing_file_is_false(tmp_path: Path) -> None:
    from giap.setupdb import verify_go

    assert verify_go(tmp_path / "nope.obo", RELEASE) is False


_dates = st.dates(min_value=date(2015, 1, 1), max_value=date(2035, 12, 31)).map(date.isoformat)


@settings(max_examples=60, deadline=None)
@given(header=_dates, pinned=_dates, slim=st.booleans())
def test_verify_go_property(header: str, pinned: str, slim: bool) -> None:
    """Property: an OBO file verifies exactly when its header release equals the pin."""
    from giap.setupdb import verify_go

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "go.obo"
        p.write_bytes(go_obo(header, slim))
        assert verify_go(p, pinned) is (header == pinned)


def test_verify_kegg_links_example(tmp_path: Path) -> None:
    from giap.setupdb import verify_kegg_links

    assert verify_kegg_links(_write(tmp_path, KO_BODY)) is True
    live_shape = b"ko:K00001\tpath:map00010\nko:K00001\tpath:ko00010\nko:K00002\tpath:map00010\n"
    assert verify_kegg_links(_write(tmp_path, live_shape)) is True


@pytest.mark.parametrize(
    "body",
    [b"", b"\n", b"<html>Forbidden</html>\n", b"ko:K00001\tpath:ko00010\n",
     b"ko:K0001\tpath:map00010\n", b"ko:K00001 path:map00010\n", b"map00010\tGlycolysis\n"],
)  # fmt: skip
def test_verify_kegg_links_rejects(tmp_path: Path, body: bytes) -> None:
    from giap.setupdb import verify_kegg_links

    assert verify_kegg_links(_write(tmp_path, body)) is False


def test_verify_kegg_links_missing_file_is_false(tmp_path: Path) -> None:
    from giap.setupdb import verify_kegg_links

    assert verify_kegg_links(tmp_path / "nope.tsv") is False


_ko = st.integers(0, 99_999).map(lambda i: f"ko:K{i:05d}")
_map = st.integers(0, 99_999).map(lambda i: f"path:map{i:05d}")


@settings(max_examples=40, deadline=None)
@given(st.lists(st.tuples(_ko, _map), min_size=1, max_size=30))
def test_verify_kegg_links_property(pairs: list[tuple[str, str]]) -> None:
    """Property: any non-empty KO -> map link table verifies."""
    from giap.setupdb import verify_kegg_links

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "ko.tsv"
        p.write_text("".join(f"{k}\t{m}\n" for k, m in pairs))
        assert verify_kegg_links(p) is True


# ---------------------------------------------------------------- count_fasta_records


def test_count_fasta_records_examples(tmp_path: Path) -> None:
    from giap.setupdb import count_fasta_records

    assert count_fasta_records(_write(tmp_path, b"")) == 0
    assert count_fasta_records(_write(tmp_path, fasta(381))) == 381
    assert count_fasta_records(_write(tmp_path, b">a\nMK\n>b\n\n>c desc > x\nM")) == 3


@settings(max_examples=40, deadline=None)
@given(
    st.lists(
        st.tuples(st.text("ABCXYZ_0123", min_size=1, max_size=8),
                  st.lists(st.text("ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=60), max_size=4)),
        max_size=40,
    )
)  # fmt: skip
def test_count_fasta_records_property(records: list[tuple[str, list[str]]]) -> None:
    """Property: the count equals the number of ``>`` header lines, whatever the line wrapping."""
    from giap.setupdb import count_fasta_records

    body = "".join(f">{h}\n" + "".join(f"{s}\n" for s in seq) for h, seq in records)
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "x.fasta"
        p.write_text(body)
        assert count_fasta_records(p) == len(records)


# ---------------------------------------------------------------- example 1 and 2: GO


def test_example1_go_ok(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    routes = pin_go(monkeypatch)
    fetcher = MapFetcher(dict(routes))
    assert setup(tmp_path, fetcher, go=True) == 0
    assert sorted(fetcher.urls) == sorted(routes)
    got = items(tmp_path)
    for name in GO_NAMES:
        rec, item = got[name], sources()[name]
        assert rec["status"] == "ok", rec
        assert RELEASE in rec["version"]
        assert (rec["url"], rec["checksum"], rec["path"]) == (item.url, item.checksum, item.path)
        assert rec["size_bytes"] == len(routes[item.url])
        assert rec["retrieved"] in today_dates()
        assert item_file(tmp_path, name).read_bytes() == routes[item.url]
    assert leftovers(tmp_path) == []
    validate_manifest(tmp_path)


def test_example1_on_the_command_line(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from click.testing import CliRunner

    from giap.cli import cli
    from giap.setupdb import HttpsFetcher

    fetcher = MapFetcher(dict(pin_go(monkeypatch)))
    monkeypatch.setattr(HttpsFetcher, "fetch", lambda _self, url, dest: fetcher.fetch(url, dest))
    db = tmp_path / "D"
    argv = ["setup", "--db-dir", str(db), "--bakta", "skip", "--no-genomad",
            "--no-promotech-models", "--go"]  # fmt: skip
    r = CliRunner().invoke(cli, argv)
    assert r.exit_code == 0, r.output
    assert {items(db)[n]["status"] for n in GO_NAMES} == {"ok"}


def test_example2_go_version_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    fetcher = MapFetcher(dict(pin_go(monkeypatch, served=OTHER_RELEASE)))  # checksums match
    assert setup(tmp_path, fetcher, go=True) == 2
    err = capsys.readouterr().err
    assert re.search(r"go_basic:[^\n]*version mismatch \(expected releases/2026-07-01\)", err), err
    rec = items(tmp_path)["go_basic"]
    assert rec["status"] != "ok", rec
    assert not item_file(tmp_path, "go_basic").exists()
    assert files_under(tmp_path / "go") == []
    assert leftovers(tmp_path) == []


def test_go_checksum_mismatch_is_mismatch_and_names_host(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    routes = pin_go(monkeypatch)
    url = sources()["go_basic"].url
    routes[url] = routes[url] + b"\n[Term]\nid: GO:9999999\n"  # right release, wrong bytes
    assert setup(tmp_path, MapFetcher(dict(routes)), go=True) == 2
    rec = items(tmp_path)["go_basic"]
    assert rec["status"] == "mismatch", rec
    assert "purl.obolibrary.org" in rec["message"] + capsys.readouterr().err
    assert not item_file(tmp_path, "go_basic").exists()


# ---------------------------------------------------------------- example 3 and KEGG


def test_example3_kegg_checksum_recorded_version_is_retrieval_date(tmp_path: Path) -> None:
    assert setup(tmp_path, MapFetcher(dict(kegg_routes())), kegg=True) == 0
    got = items(tmp_path)
    rec = got["kegg_ko_pathway"]
    assert rec["status"] == "ok", rec
    assert rec["checksum"] == sha256(KO_BODY)
    assert rec["version"] in today_dates() and rec["retrieved"] in today_dates()
    assert rec["url"] == KO_URL and rec["size_bytes"] == len(KO_BODY)
    assert item_file(tmp_path, "kegg_ko_pathway").read_bytes() == KO_BODY
    names = got["kegg_pathway_names"]
    assert names["status"] == "ok" and names["checksum"] == sha256(PATHWAY_BODY)
    validate_manifest(tmp_path)


@pytest.mark.parametrize(
    ("name", "body"),
    [("kegg_ko_pathway", b"<html>Service unavailable</html>\n"),
     ("kegg_ko_pathway", b""),
     ("kegg_pathway_names", b"<html>Service unavailable</html>\n"),
     ("kegg_pathway_names", b"path:map00010\tGlycolysis\n"),
     ("kegg_pathway_names", b"")],
)  # fmt: skip
def test_kegg_content_that_does_not_verify_is_not_installed(
    tmp_path: Path, name: str, body: bytes
) -> None:
    routes: dict[str, bytes | Exception] = {**kegg_routes(), sources()[name].url: body}
    assert setup(tmp_path, MapFetcher(routes), kegg=True) == 2
    rec = items(tmp_path)[name]
    assert rec["status"] != "ok", rec
    assert "rest.kegg.jp" in rec["message"]
    assert not item_file(tmp_path, name).exists()
    assert leftovers(tmp_path) == []


def test_kegg_check_uses_the_recorded_checksum(tmp_path: Path) -> None:
    assert setup(tmp_path, MapFetcher(dict(kegg_routes())), kegg=True) == 0
    assert setup(tmp_path, kegg=True, check=True) == 0  # NoFetcher: no server contact
    with item_file(tmp_path, "kegg_ko_pathway").open("ab") as fh:
        fh.write(b"ko:K00002\tpath:map00010\n")  # still well-formed, but not what was retrieved
    assert setup(tmp_path, kegg=True, check=True) == 2
    assert BLOCKED == []


def test_kegg_rerun_adopts_without_fetching(tmp_path: Path) -> None:
    assert setup(tmp_path, MapFetcher(dict(kegg_routes())), kegg=True) == 0
    before = items(tmp_path)
    assert setup(tmp_path, NoFetcher(), kegg=True) == 0
    after = items(tmp_path)
    for name in KEGG_NAMES:
        assert after[name]["status"] == "ok"
        assert after[name]["checksum"] == before[name]["checksum"]


# ---------------------------------------------------------------- example 4: contaminants


def test_example4_contaminants_380_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    fetcher = MapFetcher(dict(pin_contaminants(monkeypatch, fasta(380))))  # checksum matches
    assert setup(tmp_path, fetcher, contaminants=True) == 2
    rec = items(tmp_path)["contaminants"]
    assert rec["status"] == "mismatch", rec
    assert "expected 381 records" in rec["message"] + capsys.readouterr().err
    assert not item_file(tmp_path, "contaminants").exists()
    assert leftovers(tmp_path) == []


def test_contaminants_381_records_ok(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    body = fasta(381)
    assert setup(tmp_path, MapFetcher(dict(pin_contaminants(monkeypatch, body))),
                 contaminants=True) == 0  # fmt: skip
    rec = items(tmp_path)["contaminants"]
    assert rec["status"] == "ok" and rec["checksum"] == sha256(body), rec
    assert item_file(tmp_path, "contaminants").read_bytes() == body


def test_contaminants_checksum_mismatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    url = next(iter(pin_contaminants(monkeypatch, fasta(381))))
    other = fasta(381).replace(b"MKVLAAGIVG", b"MKVLAAGIVA", 1)  # 381 records, other bytes
    assert setup(tmp_path, MapFetcher({url: other}), contaminants=True) == 2
    rec = items(tmp_path)["contaminants"]
    assert rec["status"] == "mismatch", rec
    assert sources()["contaminants"].host in rec["message"]
    assert not item_file(tmp_path, "contaminants").exists()


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(n=st.sampled_from([0, 1, 380, 381, 382, 762]))
def test_contaminant_count_property(monkeypatch: pytest.MonkeyPatch, n: int) -> None:
    """Property: a checksum-valid library installs exactly when it holds 381 records."""
    with tempfile.TemporaryDirectory() as tmp, monkeypatch.context() as mp:
        db = Path(tmp) / "db"
        code = setup(db, MapFetcher(dict(pin_contaminants(mp, fasta(n)))), contaminants=True)
        rec = items(db)["contaminants"]
        assert (code == 0) is (n == 381) and (rec["status"] == "ok") is (n == 381), (n, rec)
        assert item_file(db, "contaminants").exists() is (n == 381)


# ---------------------------------------------------------------- example 5: COG timeout


def test_example5_cog_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    pin_cog(monkeypatch)
    err = source_error("ftp.ncbi.nlm.nih.gov", "timed out")
    fetcher = MapFetcher({COG_DEF_URL: err, COG_FUN_URL: err})
    assert setup(tmp_path, fetcher, cog=True) == 2
    assert "ftp.ncbi.nlm.nih.gov" in capsys.readouterr().err
    got = items(tmp_path)
    for name in COG_NAMES:
        assert got[name]["status"] == "absent", got[name]
        assert "ftp.ncbi.nlm.nih.gov" in got[name]["message"]
    assert files_under(tmp_path / "cog") == [], "nothing left under D/cog/"
    assert leftovers(tmp_path) == []


def test_cog_ok_and_checksum_mismatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    routes: dict[str, bytes | Exception] = {**pin_cog(monkeypatch)}
    routes[COG_FUN_URL] = COG_FUN_BODY + b"K\tX\tTranscription\n"
    assert setup(tmp_path, MapFetcher(routes), cog=True) == 2
    got = items(tmp_path)
    assert got["cog2020_def"]["status"] == "ok"
    assert got["cog2020_def"]["checksum"] == sha256(COG_DEF_BODY)
    assert got["cog2020_fun"]["status"] == "mismatch"
    assert files_under(tmp_path / "cog") == [item_file(tmp_path, "cog2020_def")]


# ---------------------------------------------------------------- shared directories


def test_failing_goslim_leaves_go_basic_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    routes: dict[str, bytes | Exception] = {**pin_go(monkeypatch)}
    routes[sources()["goslim_generic"].url] = source_error("purl.obolibrary.org")
    assert setup(tmp_path, MapFetcher(routes), go=True) == 2
    got = items(tmp_path)
    assert got["go_basic"]["status"] == "ok" and got["goslim_generic"]["status"] == "absent"
    assert files_under(tmp_path / "go") == [item_file(tmp_path, "go_basic")]
    assert leftovers(tmp_path) == []


def test_failing_download_keeps_an_installed_neighbour(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    routes = pin_go(monkeypatch)
    assert setup(tmp_path, MapFetcher(dict(routes)), go=True) == 0
    basic = item_file(tmp_path, "go_basic").read_bytes()
    item_file(tmp_path, "goslim_generic").unlink()
    failing: dict[str, bytes | Exception] = {
        sources()["goslim_generic"].url: source_error("purl.obolibrary.org")
    }
    assert setup(tmp_path, MapFetcher(failing), go=True) == 2
    assert item_file(tmp_path, "go_basic").read_bytes() == basic
    assert items(tmp_path)["go_basic"]["status"] == "ok"
    assert files_under(tmp_path / "go") == [item_file(tmp_path, "go_basic")]


FILE_ITEMS = ("go_basic", "goslim_generic", "cog2020_def", "cog2020_fun")


@settings(
    max_examples=16, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(failing=st.sets(st.sampled_from(FILE_ITEMS)), how=st.sampled_from(["source", "bytes"]))
def test_partial_cleanup_property(
    monkeypatch: pytest.MonkeyPatch, failing: set[str], how: str
) -> None:
    """Property: a failing item leaves no file of its own and never removes a sibling's file."""
    with tempfile.TemporaryDirectory() as tmp, monkeypatch.context() as mp:
        db = Path(tmp) / "db"
        routes: dict[str, bytes | Exception] = {**pin_go(mp), **pin_cog(mp)}
        for name in failing:
            item = sources()[name]
            routes[item.url] = (source_error(item.host) if how == "source"
                                else bytes(routes[item.url]) + b"tampered")  # type: ignore[operator]  # fmt: skip
        code = setup(db, MapFetcher(routes), go=True, cog=True)
        assert code == (2 if failing else 0)
        got = items(db)
        for name in FILE_ITEMS:
            assert (got[name]["status"] == "ok") is (name not in failing), (name, got[name])
            assert item_file(db, name).exists() is (name not in failing), name
        expected = sorted(item_file(db, n) for n in FILE_ITEMS if n not in failing)
        assert sorted(files_under(db / "go") + files_under(db / "cog")) == expected
        assert leftovers(db) == []


# ---------------------------------------------------------------- example 6: --check, offline


def test_example6_check_after_full_install_is_offline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    every = {"go": True, "kegg": True, "cog": True, "contaminants": True}
    fetcher = MapFetcher(dict(pin_all(monkeypatch)))
    assert setup(tmp_path, fetcher, **every) == 0
    assert {items(tmp_path)[n]["status"] for n in NEW_NAMES} == {"ok"}
    validate_manifest(tmp_path)
    manifest_before = (tmp_path / "db_manifest.json").read_bytes()
    assert setup(tmp_path, NoFetcher(), NoRunner(), check=True, **every) == 0
    assert (tmp_path / "db_manifest.json").read_bytes() == manifest_before
    assert BLOCKED == []


def test_example6_on_the_command_line(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from click.testing import CliRunner

    from giap.cli import cli

    every = {"go": True, "kegg": True, "cog": True, "contaminants": True}
    assert setup(tmp_path, MapFetcher(dict(pin_all(monkeypatch))), **every) == 0
    argv = ["setup", "--db-dir", str(tmp_path), "--bakta", "skip", "--no-genomad",
            "--no-promotech-models", "--check", "--go", "--kegg", "--cog",
            "--contaminants"]  # fmt: skip
    r = CliRunner().invoke(cli, argv)
    assert r.exit_code == 0, r.output
    assert BLOCKED == [], "the guard saw a network or tool call"


def test_check_with_nothing_installed_is_absent_offline(tmp_path: Path) -> None:
    every = {"go": True, "kegg": True, "cog": True, "contaminants": True}
    assert setup(tmp_path, NoFetcher(), NoRunner(), check=True, **every) == 2
    assert BLOCKED == []
    assert not (tmp_path / "go").exists() and not (tmp_path / "kegg").exists()


def test_check_detects_a_tampered_pinned_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert setup(tmp_path, MapFetcher(dict(pin_cog(monkeypatch))), cog=True) == 0
    item_file(tmp_path, "cog2020_fun").write_bytes(COG_FUN_BODY + b"X\n")
    assert setup(tmp_path, cog=True, check=True) == 2


# ---------------------------------------------------------------- HttpsFetcher (backlog)


class _Sock:
    def __init__(self, raw: bytes) -> None:
        self.raw = raw

    def makefile(self, *_a: Any, **_k: Any) -> io.BytesIO:
        return io.BytesIO(self.raw)


Route = tuple[int, dict[str, str], bytes]


class FakeNet:
    """Fake ``http.client``: ``routes[(scheme, host, path)] = (status, headers, body)``.

    ``headers`` may override ``Content-Length`` (to fake a truncated body); ``raise_on`` makes a
    request to that key raise the given exception (a timeout)."""

    def __init__(self) -> None:
        self.routes: dict[tuple[str, str, str], Route] = {}
        self.raise_on: dict[tuple[str, str, str], Exception] = {}
        self.contacted: list[tuple[str, str, str]] = []

    def conn(self, scheme: str) -> type:
        net = self

        class Conn(http.client.HTTPConnection):
            def __init__(self, host: str, port: int | None = None, **_kw: Any) -> None:
                super().__init__(host, port)
                self.path = ""

            def connect(self) -> None:
                raise AssertionError("FakeNet never opens a socket")

            def close(self) -> None: ...

            def request(self, method: str, url: str, *_a: Any, **_k: Any) -> None:  # type: ignore[override]
                self.path = url

            def getresponse(self) -> http.client.HTTPResponse:
                key = (scheme, self.host, self.path.split("?")[0])
                net.contacted.append(key)
                if key in net.raise_on:
                    raise net.raise_on[key]
                status, headers, body = net.routes.get(key, (404, {}, b"not found"))
                headers = {"Content-Length": str(len(body))} | headers
                head = "".join(f"{k}: {v}\r\n" for k, v in headers.items())
                raw = f"HTTP/1.1 {status} X\r\n{head}Connection: close\r\n\r\n".encode() + body
                resp = http.client.HTTPResponse(_Sock(raw), method="GET")  # type: ignore[arg-type]
                resp.begin()
                return resp

        return Conn


@pytest.fixture
def net(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    request.getfixturevalue("real_https_fetch")
    fake = FakeNet()
    monkeypatch.setattr(http.client, "HTTPSConnection", fake.conn("https"))
    monkeypatch.setattr(http.client, "HTTPConnection", fake.conn("http"))
    return fake


PURL_PATH = f"/obo/go/releases/{RELEASE}/go-basic.obo"
GO_RELEASE_PATH = f"/{RELEASE}/ontology/go-basic.obo"


def test_purl_redirect_to_geneontology_over_https_is_followed(net: FakeNet, tmp_path: Path) -> None:
    """The live purl answers 302 -> https://release.geneontology.org/... (§3.2 lists both)."""
    from giap.setupdb import HttpsFetcher

    body = go_obo(RELEASE)
    net.routes[("https", "purl.obolibrary.org", PURL_PATH)] = (
        302, {"Location": "https://release.geneontology.org" + GO_RELEASE_PATH}, b"")  # fmt: skip
    net.routes[("https", "release.geneontology.org", GO_RELEASE_PATH)] = (200, {}, body)
    HttpsFetcher().fetch("https://purl.obolibrary.org" + PURL_PATH, tmp_path / "f")
    assert (tmp_path / "f").read_bytes() == body
    assert net.contacted[-1] == ("https", "release.geneontology.org", GO_RELEASE_PATH)


@pytest.mark.parametrize(
    "location",
    ["http://release.geneontology.org" + GO_RELEASE_PATH,
     "https://evil.example" + GO_RELEASE_PATH,
     "https://release.geneontology.org.evil.example" + GO_RELEASE_PATH,
     "https://geneontology.org.evil.example" + GO_RELEASE_PATH,
     "https://evilgeneontology.org" + GO_RELEASE_PATH,
     "ftp://release.geneontology.org" + GO_RELEASE_PATH],
)  # fmt: skip
def test_purl_redirect_elsewhere_is_refused(net: FakeNet, tmp_path: Path, location: str) -> None:
    from urllib.parse import urlparse

    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "purl.obolibrary.org", PURL_PATH)] = (302, {"Location": location}, b"")
    u = urlparse(location)
    net.routes[(u.scheme, u.hostname or "", u.path)] = (200, {}, b"EVIL")
    with pytest.raises(SourceError):
        HttpsFetcher().fetch("https://purl.obolibrary.org" + PURL_PATH, tmp_path / "f")
    assert net.contacted == [("https", "purl.obolibrary.org", PURL_PATH)]
    assert not (tmp_path / "f").exists() or b"EVIL" not in (tmp_path / "f").read_bytes()


def test_kegg_redirect_to_another_host_is_refused(net: FakeNet, tmp_path: Path) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "rest.kegg.jp", "/link/pathway/ko")] = (
        301, {"Location": "https://release.geneontology.org/link/pathway/ko"}, b"")  # fmt: skip
    net.routes[("https", "release.geneontology.org", "/link/pathway/ko")] = (200, {}, KO_BODY)
    with pytest.raises(SourceError):
        HttpsFetcher().fetch(KO_URL, tmp_path / "f")
    assert net.contacted == [("https", "rest.kegg.jp", "/link/pathway/ko")]


@pytest.mark.parametrize("status", [204, 206, 403, 404, 410, 500, 503])
def test_non_200_is_rejected_before_writing(net: FakeNet, tmp_path: Path, status: int) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "rest.kegg.jp", "/link/pathway/ko")] = (status, {}, KO_BODY)
    with pytest.raises(SourceError) as exc:
        HttpsFetcher().fetch(KO_URL, tmp_path / "f")
    assert exc.value.host == "rest.kegg.jp" and str(status) in exc.value.reason
    assert not (tmp_path / "f").exists(), "nothing is written for a non-200 response"


def test_non_200_through_install_leaves_item_absent(net: FakeNet, tmp_path: Path) -> None:
    from giap.setupdb import HttpsFetcher, install

    net.routes[("https", "rest.kegg.jp", "/link/pathway/ko")] = (503, {}, KO_BODY)
    rec = install(sources()["kegg_ko_pathway"], tmp_path, NoRunner(), HttpsFetcher())
    assert rec.status == "absent" and "rest.kegg.jp" in rec.message, rec
    assert files_under(tmp_path) == [] or all(".lock" in p.name for p in files_under(tmp_path))


def test_truncated_body_is_an_error(net: FakeNet, tmp_path: Path) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "rest.kegg.jp", "/link/pathway/ko")] = (
        200, {"Content-Length": "1000"}, KO_BODY)  # fmt: skip
    with pytest.raises(SourceError):
        HttpsFetcher().fetch(KO_URL, tmp_path / "f")


def test_timeout_is_a_source_error_naming_the_host(net: FakeNet, tmp_path: Path) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.raise_on[("https", "ftp.ncbi.nlm.nih.gov", "/pub/COG/COG2020/data/fun-20.tab")] = (
        TimeoutError("timed out"))  # fmt: skip
    with pytest.raises(SourceError) as exc:
        HttpsFetcher().fetch(COG_FUN_URL, tmp_path / "f")
    assert exc.value.host == "ftp.ncbi.nlm.nih.gov"


def test_redirect_loop_is_capped(net: FakeNet, tmp_path: Path) -> None:
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes[("https", "rest.kegg.jp", "/a")] = (302, {"Location": "/b"}, b"")
    net.routes[("https", "rest.kegg.jp", "/b")] = (302, {"Location": "https://rest.kegg.jp/a"}, b"")
    with pytest.raises(SourceError, match="redirect"):
        HttpsFetcher().fetch("https://rest.kegg.jp/a", tmp_path / "f")
    assert 2 <= len(net.contacted) <= 11, net.contacted
    assert not (tmp_path / "f").exists()


@settings(
    max_examples=20, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(hops=st.integers(0, 30))
def test_redirect_chain_property(net: FakeNet, hops: int) -> None:
    """Property: a same-host https chain is followed when short, refused when long; requests are
    bounded either way and never exceed one more than the redirects allowed."""
    from giap.setupdb import HttpsFetcher, SourceError

    net.routes.clear()
    net.contacted.clear()
    for i in range(hops):
        net.routes[("https", "rest.kegg.jp", f"/r{i}")] = (302, {"Location": f"/r{i + 1}"}, b"")
    net.routes[("https", "rest.kegg.jp", f"/r{hops}")] = (200, {}, KO_BODY)
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp) / "f"
        try:
            HttpsFetcher().fetch("https://rest.kegg.jp/r0", dest)
        except SourceError:
            assert hops >= 3, "a few redirects are always followed"
            assert len(net.contacted) <= 11
            assert not dest.exists()
        else:
            assert hops <= 10, "at most 10 redirects"
            assert dest.read_bytes() == KO_BODY
            assert len(net.contacted) == hops + 1


# ---------------------------------------------------------------- geNOMAD archive (backlog)


def _genomad(payload: bytes) -> Any:
    return dataclasses.replace(sources()["genomad"], checksum=md5(payload))


def _archives(root: Path) -> list[Path]:
    return [p for p in root.rglob("*") if p.name.endswith((".tar.gz", ".tgz", ".tar"))]


@pytest.mark.parametrize(
    "payload",
    [GENOMAD_TAR, _tar([("version.txt", b"1.9\n"), ("genomad_db", b"mmseqs")])],
    ids=["nested", "flat"],
)
def test_genomad_tarball_deleted_after_unpacking(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, payload: bytes
) -> None:
    """§3.9: the archive is gone before the unpacked copy is verified and moved into place."""
    import giap.setupdb.run as run_mod
    from giap.setupdb import install

    seen: list[list[Path]] = []
    real = getattr(run_mod, "check_dir", None)
    if real is not None:

        def spy(item: Any, d: Path) -> Any:
            seen.append(_archives(tmp_path))
            return real(item, d)

        monkeypatch.setattr(run_mod, "check_dir", spy)
    rec = install(
        _genomad(payload), tmp_path, NoRunner(), MapFetcher({sources()["genomad"].url: payload})
    )
    assert rec.status == "ok", rec
    assert all(s == [] for s in seen), f"archive still on disk at verification: {seen}"
    assert _archives(tmp_path) == []
    assert (tmp_path / "genomad_db" / "version.txt").read_text().strip() == "1.9"


def _special_tar(kind: str) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        good = tarfile.TarInfo("genomad_db/version.txt")
        good.size = 4
        tar.addfile(good, io.BytesIO(b"1.9\n"))
        info = tarfile.TarInfo("genomad_db/odd")
        if kind == "fifo":
            info.type = tarfile.FIFOTYPE
        elif kind == "chardev":
            info.type, info.devmajor, info.devminor = tarfile.CHRTYPE, 1, 3
        elif kind == "hardlink-out":
            info.type, info.linkname = tarfile.LNKTYPE, "../../escape"
        elif kind == "hardlink-abs":
            info.type, info.linkname = tarfile.LNKTYPE, "/etc/passwd"
        tar.addfile(info)
    return buf.getvalue()


@pytest.mark.parametrize("kind", ["fifo", "chardev", "hardlink-out", "hardlink-abs"])
def test_unsafe_tar_members_are_caught_by_giaps_own_check(tmp_path: Path, kind: str) -> None:
    """The member check itself rejects these (named in the record), not tarfile's data filter."""
    from giap.setupdb import install

    db = tmp_path / "a" / "db"
    db.mkdir(parents=True)
    payload = _special_tar(kind)
    rec = install(
        _genomad(payload), db, NoRunner(), MapFetcher({sources()["genomad"].url: payload})
    )
    assert rec.status != "ok", rec
    assert "unsafe" in rec.message and "odd" in rec.message, rec.message
    assert not (db / "genomad_db").exists()
    assert leftovers(db) == []


# ---------------------------------------------------------------- _provenance (backlog)


def _old_genomad_record(db: Path, status: str) -> None:
    from giap.setupdb import DbRecord, write_manifest

    item = sources()["genomad"]
    write_manifest(db, [DbRecord(item.name, item.version, item.url, item.checksum, 1, "2025-01-02",
                                 item.path, status, "old")])  # type: ignore[arg-type]  # fmt: skip


@pytest.mark.parametrize("status", ["absent", "mismatch"])
def test_previous_record_is_reused_only_when_it_was_ok(tmp_path: Path, status: str) -> None:
    make_genomad(tmp_path)
    _old_genomad_record(tmp_path, status)
    assert setup(tmp_path, genomad=True) == 0
    rec = items(tmp_path)["genomad"]
    assert rec["status"] == "ok", rec
    assert rec["retrieved"] != "2025-01-02" and rec["message"] != "old", rec


def test_previous_ok_record_is_reused(tmp_path: Path) -> None:
    make_genomad(tmp_path)
    _old_genomad_record(tmp_path, "ok")
    assert setup(tmp_path, genomad=True) == 0
    assert items(tmp_path)["genomad"]["retrieved"] == "2025-01-02"


# ---------------------------------------------------------------- antiSMASH in the image (backlog)


class AntismashRunner:
    """``--check-prereqs`` answers from the image (no ``--databases``) or a directory; the
    downloader writes a small database tree."""

    def __init__(self, image_has_dbs: bool) -> None:
        self.image_has_dbs = image_has_dbs
        self.calls: list[list[str]] = []

    def run(self, tool: str, argv: Sequence[str]) -> int:
        argv = list(argv)
        self.calls.append(argv)
        assert tool == "antismash", tool
        if "--check-prereqs" in argv:
            if "--databases" in argv:
                d = Path(argv[argv.index("--databases") + 1])
                return 0 if d.is_dir() and any(d.iterdir()) else 1
            return 0 if self.image_has_dbs else 1
        if any("download" in a for a in argv):
            out = Path(argv[argv.index("--database-dir") + 1])
            (out / "clusterblast").mkdir(parents=True, exist_ok=True)
            (out / "clusterblast" / "proteins.fasta").write_text(">p\nMK\n")
            return 0
        return 1

    def downloads(self) -> list[list[str]]:
        return [c for c in self.calls if any("download" in a for a in c)]


def test_antismash_dbs_in_the_image_recorded_as_image(tmp_path: Path) -> None:
    runner = AntismashRunner(image_has_dbs=True)
    assert setup(tmp_path, NoFetcher(), runner, antismash=True) == 0
    rec = items(tmp_path)["antismash"]
    assert (rec["status"], rec["url"]) == ("ok", "image"), rec
    assert rec["checksum"] == ""
    assert runner.downloads() == [], "no download when the image already has the databases"
    assert any("--check-prereqs" in c for c in runner.calls)
    assert not (tmp_path / sources()["antismash"].path).exists()
    validate_manifest(tmp_path)


def test_antismash_dbs_not_in_the_image_are_downloaded(tmp_path: Path) -> None:
    runner = AntismashRunner(image_has_dbs=False)
    assert setup(tmp_path, NoFetcher(), runner, antismash=True) == 0
    rec = items(tmp_path)["antismash"]
    assert rec["status"] == "ok" and rec["url"] == sources()["antismash"].url, rec
    assert rec["url"].startswith("https://dl.secondarymetabolites.org/")
    assert len(runner.downloads()) == 1
    assert leftovers(tmp_path) == []


def test_antismash_installed_copy_is_adopted_not_image(tmp_path: Path) -> None:
    d = tmp_path / sources()["antismash"].path / "clusterblast"
    d.mkdir(parents=True)
    (d / "proteins.fasta").write_text(">p\nMK\n")
    runner = AntismashRunner(image_has_dbs=True)
    assert setup(tmp_path, NoFetcher(), runner, antismash=True) == 0
    rec = items(tmp_path)["antismash"]
    assert rec["status"] == "ok" and rec["url"] != "image", rec
    assert runner.downloads() == []


# ---------------------------------------------------------------- repository (R5, R11)


def test_notice_names_kegg_academic_use_and_the_contaminant_library() -> None:
    text = (REPO / "NOTICE.md").read_text(encoding="utf-8")
    assert re.search(r"KEGG[^\n]*academic", text)
    assert re.search(r"(?i)contaminant[^\n]*(licen[cs]e|Frankenfield)", text), (
        "NOTICE.md records the contaminant library and its licence as found"
    )


def test_contaminant_library_is_not_vendored() -> None:
    # .nf-test is Nextflow's gitignored work tree: a local `nf-test test workflow/tests` run
    # symlinks the real $GIAP_DB_DIR library into it (mapping_nf_main_review.nf.test), which is
    # staging, not vendoring. Nothing under a gitignored directory can be vendored.
    skip = {".git", "local", ".venv", "node_modules", ".hypothesis", ".mypy_cache", ".ruff_cache",
            ".nf-test"}  # fmt: skip
    hits = []
    for root, dirs, names in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in skip]
        hits += [Path(root, n) for n in names if "contaminant" in n.lower()
                 and n.lower().endswith((".fasta", ".fa", ".faa", ".fasta.gz"))]  # fmt: skip

    # Hand-written test excerpts (a few records) are fine; the real library has 381 records (R11).
    def n_records(path: Path) -> int:
        if path.suffix == ".gz":
            return 10**6
        return path.read_text(encoding="utf-8", errors="replace").count(">")

    vendored = [h for h in hits if n_records(h) > 20]
    assert vendored == [], vendored


def test_giap_setup_sbatch_installs_the_resources() -> None:
    text = (REPO / "scripts/slurm/giap_setup.sbatch").read_text(encoding="utf-8")
    default = next(line for line in text.splitlines() if "set --" in line)
    for flag in ("--go", "--kegg", "--cog", "--contaminants"):
        assert re.search(rf"{flag}(\s|;|$)", default), (flag, default)
