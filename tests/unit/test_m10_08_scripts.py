"""``fetch_a3_assets.sbatch`` and ``mzml_pilot.sbatch`` in a sandbox (M10-08a, deliverables 6-7).

Spec: docs/specs/M10-08-a3-feasibility.md, deliverables 5-7, acceptance criteria 7 and 8, edge
cases "Disk" and "Weights sha256 differs", the job table (steps 2 and 3) and R6/R14.

Every script runs with ``bash`` in ``helpers.a3_kit.a3_sandbox`` (a throw-away checkout, fake
Slurm tools, fake ``curl``/``wget`` that only log and write a stand-in payload, an unreachable
proxy, fake ``micromamba``/``apptainer``/``df``/``sha256sum``/``shasum``). No network, no GPU, no
real raw data, TIMSCONVERT, ProteoWizard or Cascadia. The weights URL is the committed one from
``validation/pilot/a3_assets.tsv``; only its sha256/bytes columns are rewritten in the sandbox
copy to match the fake payload. No statistical claim, so no calibration test (spec §8.4).
"""

from __future__ import annotations

import csv
import gzip
import re
import stat
from pathlib import Path

import pytest
from helpers.a3_kit import (
    A3_ASSETS,
    ASSET_COLUMNS,
    FETCH_SBATCH,
    MZML_SBATCH,
    PWIZ_DIGEST,
    STUB_UNIMOD_XML,
    A3Sandbox,
    a3_sandbox,
    example11_spectra,
    require,
    sha256_hex,
    weights_payload,
    write_mzml,
    wt1_spectra,
)


@pytest.fixture(scope="module")
def good_template(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A WT_1-shaped mzML that passes mzml-check against the committed acquisition.json."""
    return write_mzml(tmp_path_factory.mktemp("tpl") / "good.mzML", wt1_spectra(), indexed=True)


@pytest.fixture
def box(tmp_path: Path, good_template: Path) -> A3Sandbox:
    return a3_sandbox(tmp_path, good_template)


def _out(res) -> str:  # type: ignore[no-untyped-def]
    return res.stdout + res.stderr


def _slurm_seconds(value: str) -> int:
    """``[D-]HH:MM:SS``, ``MM:SS`` or ``MM`` (sbatch --time) in seconds."""
    days, _, hms = value.rpartition("-")
    parts = [int(x) for x in hms.split(":")]
    if len(parts) == 1:
        parts = [0, parts[0], 0]
    elif len(parts) == 2:
        parts = [0, *parts]
    h, m, s = parts
    return (int(days or 0) * 24 + h) * 3600 + m * 60 + s


def _header_value(text: str, key: str) -> str | None:
    m = re.search(rf"^#SBATCH --{key}[= ](\S+)$", text, re.M)
    return m.group(1) if m else None


# ================================================================== fetch_a3_assets.sbatch


def _committed_assets() -> list[dict[str, str]]:
    with require(A3_ASSETS, "deliverable 5").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def _weights_url() -> str:
    return next(r["url"] for r in _committed_assets() if r["asset"] == "cascadia_weights")


def _payload(url: str) -> bytes:
    return weights_payload(url)  # what helpers.a3_kit's fake curl writes (zip magic first)


def _assets_for_fetch(
    box: A3Sandbox, *, weights_sha: str | None = None, unimod_sha: str | None = None
) -> None:
    rows = []
    for r in _committed_assets():
        r = dict(r)
        if r["asset"] == "cascadia_weights":
            data = _payload(r["url"])
            r["sha256"] = sha256_hex(data) if weights_sha is None else weights_sha
            r["bytes"] = str(len(data)) if r["sha256"] else ""
        elif r["asset"] == "unimod_table":
            r["sha256"] = sha256_hex(STUB_UNIMOD_XML) if unimod_sha is None else unimod_sha
            r["bytes"] = str(len(STUB_UNIMOD_XML)) if r["sha256"] else ""
        rows.append(r)
    box.write_assets(rows)


def _weights(box: A3Sandbox) -> Path:
    return box.db / "cascadia" / "cascadia_0.0.7.ckpt"


def _unimod(box: A3Sandbox) -> Path:
    return box.db / "cascadia" / "unimod_tables.xml"


def test_fetch_header_and_no_plain_http() -> None:
    text = require(FETCH_SBATCH, "deliverable 6").read_text()
    assert text.splitlines()[0] == "#!/bin/bash"
    name = _header_value(text, "job-name")
    assert name is not None and name.startswith("giap-")
    assert _header_value(text, "partition") == "normal"
    assert _header_value(text, "cpus-per-task") == "1"
    assert _header_value(text, "mem") == "4G"
    time = _header_value(text, "time")
    assert time is not None and _slurm_seconds(time) <= 2 * 3600, time  # S2: about 2 h at most
    assert _header_value(text, "gres") is None
    assert "set -euo pipefail" in text
    assert "scripts/farmshare_env.sh" in text
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert "http://" not in code
    assert not re.search(r"(^|\s)(-k|--insecure|--no-check-certificate)(\s|$)", code, re.M)


def test_fetch_without_adr_m10_b_exits_2_before_any_request(box: A3Sandbox) -> None:
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", unset=("A3_WEIGHTS_HOST_OK",))
    assert res.returncode == 2, _out(res)
    assert "ADR-M10-B" in _out(res)
    assert box.pilot.fetch_calls() == []
    assert not _weights(box).exists()


def test_fetch_happy_path(box: A3Sandbox) -> None:
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 0, _out(res)
    allowed = {r["url"] for r in _committed_assets()}
    urls = box.urls()
    assert urls, "the weights were not requested"
    assert set(urls) <= allowed, urls
    assert all(u.startswith("https://") for u in urls), urls
    assert _weights_url() in urls
    assert not any("unimod" in u for u in urls), "the Unimod table comes from psims, offline"
    assert not any(c["insecure"] for c in box.pilot.fetch_calls())
    w, u = _weights(box), _unimod(box)
    assert w.read_bytes() == _payload(_weights_url())
    assert u.read_bytes() == STUB_UNIMOD_XML
    for p in (w, u):
        assert stat.S_IMODE(p.stat().st_mode) == 0o444, (p, oct(p.stat().st_mode))
    dirty = [ln for ln in box.pilot.repo_is_clean().splitlines() if "a3_assets.tsv" not in ln]
    assert not dirty, f"the fetch wrote into the checkout: {dirty}"


def test_fetch_weights_sha_mismatch_exits_1_and_leaves_no_file(box: A3Sandbox) -> None:
    _assets_for_fetch(box, weights_sha="0" * 64)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 1, _out(res)
    assert not _weights(box).exists()


def test_fetch_unimod_sha_mismatch_exits_1_and_leaves_no_file(box: A3Sandbox) -> None:
    _assets_for_fetch(box, unimod_sha="f" * 64)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 1, _out(res)
    assert not _unimod(box).exists()


def test_fetch_refuses_unverified_weights(box: A3Sandbox) -> None:
    # R6: never keep an unverified file; a blank sha256 needs --resolve first.
    _assets_for_fetch(box, weights_sha="")
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode in (1, 2), _out(res)
    assert not _weights(box).exists()


def test_fetch_resolve_prints_sha256_and_bytes(box: A3Sandbox) -> None:
    _assets_for_fetch(box, weights_sha="", unimod_sha="")
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", "--resolve", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 0, _out(res)
    data = _payload(_weights_url())
    assert sha256_hex(data) in res.stdout
    assert re.search(rf"\b{len(data)}\b", res.stdout)
    assert sha256_hex(STUB_UNIMOD_XML) in res.stdout
    assert re.search(rf"\b{len(STUB_UNIMOD_XML)}\b", res.stdout)
    assert set(box.urls()) <= {r["url"] for r in _committed_assets()}


def test_fetch_resolve_still_needs_adr_m10_b(box: A3Sandbox) -> None:
    _assets_for_fetch(box, weights_sha="", unimod_sha="")
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", "--resolve", unset=("A3_WEIGHTS_HOST_OK",))
    assert res.returncode == 2, _out(res)
    assert box.pilot.fetch_calls() == []


def test_fetch_pulls_no_image_on_route_b2(box: A3Sandbox) -> None:
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 0, _out(res)
    pulls = [c for c in box.calls("apptainer.jsonl") if c["argv"][:1] == ["pull"]]
    for c in pulls:  # B1 only, and then only by digest
        assert any(PWIZ_DIGEST in a for a in c["argv"]), c


def test_fetch_corrupt_psims_table_fails(box: A3Sandbox) -> None:
    _assets_for_fetch(box)
    env = box.stub_cascadia_env()
    gz = next(env.rglob("unimod_tables.xml.gz"))
    gz.write_bytes(gzip.compress(b"<unimod>tampered</unimod>\n", mtime=0))
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 1, _out(res)
    assert not _unimod(box).exists()


# ================================================================== a3_assets.tsv (deliverable 5)


def test_assets_table_layout() -> None:
    text = require(A3_ASSETS, "deliverable 5").read_text()
    assert tuple(text.splitlines()[0].split("\t")) == ASSET_COLUMNS
    rows = {r["asset"]: r for r in _committed_assets()}
    assert {"cascadia_weights", "unimod_table", "timsconvert", "pytdfsdk"} <= set(rows)
    assert set(rows) <= {"cascadia_weights", "unimod_table", "timsconvert", "pytdfsdk",
                         "pwiz_image"}  # fmt: skip
    for r in rows.values():
        assert r["url"].startswith("https://"), r
        assert r["licence"], r
        assert r["sha256"] in ("", "NA") or re.fullmatch(r"[0-9a-f]{64}", r["sha256"]), r
        assert r["path"] == "" or "GIAP_ROOT" in r["path"] or "GIAP_DB_DIR" in r["path"], r
    w = rows["cascadia_weights"]
    assert w["url"].startswith("https://drive.usercontent.google.com/download?id=")
    assert "export=download" in w["url"] and "confirm=t" in w["url"]
    assert w["path"].endswith("cascadia/cascadia_0.0.7.ckpt")
    assert w["licence"] == "Apache-2.0"
    assert rows["unimod_table"]["path"].endswith("cascadia/unimod_tables.xml")
    assert rows["timsconvert"]["version"] == "3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b"
    assert rows["pytdfsdk"]["version"] == "4808da81d3e37f82da96f392dbb446c3733edacb"
    assert "github.com/gtluu/timsconvert" in rows["timsconvert"]["url"]
    assert "github.com/gtluu/pyTDFSDK" in rows["pytdfsdk"]["url"]
    if "pwiz_image" in rows:
        assert PWIZ_DIGEST in rows["pwiz_image"]["url"] + rows["pwiz_image"]["version"]


# ================================================================== mzml_pilot.sbatch


def _converted(box: A3Sandbox, run: str = "WT_1") -> Path:
    return box.data / "mzml" / f"{run}.mzML"


def _timsconvert_calls(box: A3Sandbox) -> list[list[str]]:
    out = []
    for c in box.calls("micromamba.jsonl"):
        argv = c["argv"]
        if "timsconvert" in argv and "--input" in argv:
            out.append(argv)
    return out


def test_mzml_pilot_header_and_no_network() -> None:
    text = require(MZML_SBATCH, "deliverable 7").read_text()
    assert text.splitlines()[0] == "#!/bin/bash"
    name = _header_value(text, "job-name")
    assert name is not None and name.startswith("giap-")
    assert _header_value(text, "partition") == "normal"
    assert _header_value(text, "cpus-per-task") == "8"
    assert _header_value(text, "mem") == "32000M"  # MaxMemPerCPU 4000M x 8 CPUs
    assert _header_value(text, "time") == "02:00:00"
    assert _header_value(text, "gres") is None
    assert "set -euo pipefail" in text
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    for word in (r"\bcurl\b", r"\bwget\b", r"git\s+clone", r"pip\s+install", r"apptainer\s+pull",
                 r"https?://"):  # fmt: skip
        assert not re.search(word, code), f"mzml_pilot.sbatch must not reach the network: {word}"


def test_mzml_pilot_default_route_timsconvert(box: A3Sandbox) -> None:
    d = box.make_run_dir("WT_1")
    before = box.files_under_root()
    res = box.run("mzml_pilot", "WT_1")
    assert res.returncode == 0, _out(res)
    calls = _timsconvert_calls(box)
    assert len(calls) == 1, box.calls("micromamba.jsonl")
    argv = calls[0]
    assert argv[0] == "run" and "giap-timsconvert" in argv
    t = argv[argv.index("timsconvert") + 1 :]
    assert t.count("--input") == 1
    assert Path(t[t.index("--input") + 1]).resolve() == d.resolve()
    assert t[t.index("--mode") + 1] == "centroid"
    assert t[t.index("--compression") + 1] == "zlib"
    assert "--exclude_mobility" in t
    assert _converted(box).is_file()
    new = box.files_under_root() - before
    mzml_dir = (box.data / "mzml").resolve()
    stray = [p for p in new if mzml_dir not in p.resolve().parents]
    assert not stray, f"writes outside $GIAP_DATA_DIR/pxd078459/mzml/: {stray}"
    assert box.pilot.repo_is_clean() == ""
    assert box.urls() == []


def test_mzml_pilot_runs_the_mzml_check_and_fails_with_it(tmp_path: Path) -> None:
    bad = write_mzml(tmp_path / "bad.mzML", example11_spectra())  # 3 MS1 ≠ WT_1's 734 (M1)
    box = a3_sandbox(tmp_path / "sb", bad)
    box.make_run_dir("WT_1")
    res = box.run("mzml_pilot", "WT_1")
    assert res.returncode == 1, _out(res)
    assert "M1" in _out(res) or any(
        "M1" in p.read_text() for p in (box.data / "mzml").glob("*.json")
    ), "the failed rule is not named anywhere"


def test_mzml_pilot_too_little_disk_exits_2_before_converting(box: A3Sandbox) -> None:
    box.make_run_dir("WT_1", tdf_bin_bytes=1_000_000)
    res = box.run("mzml_pilot", "WT_1", FAKE_DF_AVAIL_BYTES=str(3_000_000))  # < 4 x .d size
    assert res.returncode == 2, _out(res)
    assert _timsconvert_calls(box) == []
    assert not _converted(box).exists()


def test_mzml_pilot_enough_disk_converts(box: A3Sandbox) -> None:
    box.make_run_dir("WT_1", tdf_bin_bytes=1_000_000)
    res = box.run("mzml_pilot", "WT_1", FAKE_DF_AVAIL_BYTES=str(50_000_000))
    assert res.returncode == 0, _out(res)
    assert len(_timsconvert_calls(box)) == 1


def test_mzml_pilot_pwiz_route(box: A3Sandbox) -> None:
    box.make_run_dir("WT_1")
    sif = box.root / "images" / "pwiz.sif"
    sif.parent.mkdir(parents=True)
    sif.write_text("stub sif\n")
    box.run("mzml_pilot", "WT_1", MZML_ROUTE="pwiz", PWIZ_SIF=str(sif))
    calls = [c["argv"] for c in box.calls("apptainer.jsonl") if "msconvert" in c["argv"]]
    real = [a for a in calls if "--help" not in a]
    assert real, box.calls("apptainer.jsonl")
    argv = real[0]
    i = argv.index("msconvert")
    assert argv[i - 1] == "wine"
    assert "--combineIonMobilitySpectra" in argv
    assert any(re.search(r"peakPicking\s+vendor", a) for a in argv), argv
    assert _timsconvert_calls(box) == []


@pytest.mark.parametrize(
    ("args", "env"),
    [((), {}), (("WT_9",), {}), (("WT_1",), {"MZML_ROUTE": "bogus"})],
    ids=["no-run", "unknown-run", "unknown-route"],
)
def test_mzml_pilot_preconditions_exit_2(
    box: A3Sandbox, args: tuple[str, ...], env: dict[str, str]
) -> None:
    box.make_run_dir("WT_1")
    res = box.run("mzml_pilot", *args, **env)
    assert res.returncode == 2, _out(res)
    assert _timsconvert_calls(box) == []


# ================================================================== review round of PR #301


def _code(path: Path) -> str:
    return "\n".join(ln for ln in path.read_text().splitlines() if not ln.lstrip().startswith("#"))


def test_fetch_curl_time_fits_the_job_time_limit() -> None:
    # S2: every curl attempt plus the retry delay must end before Slurm kills the job, so the
    # "retry once later, then Blocked" message can print.
    text = require(FETCH_SBATCH, "deliverable 6").read_text()
    code = _code(FETCH_SBATCH)
    m = re.search(r"(?:--max-time|-m)[ =](\d+)", code)
    assert m, "curl needs --max-time (a bound on each attempt)"
    max_time = int(m.group(1))
    r = re.search(r"--retry[ =](\d+)", code)
    retries = int(r.group(1)) if r else 0
    d = re.search(r"--retry-delay[ =](\d+)", code)
    delay = int(d.group(1)) if d else 0
    time = _header_value(text, "time")
    assert time is not None
    worst = max_time * (retries + 1) + delay * retries
    assert worst <= _slurm_seconds(time), (worst, time)


_HTML_BODIES = {
    "doctype": b"<!DOCTYPE html><html><head><title>Google Drive - Virus scan warning</title>"
    b'</head><body><form action="https://drive.usercontent.google.com/download"></form>'
    b"</body></html>\n",
    "html-tag": b"\n  <html><body>Quota exceeded for this file</body></html>\n",
    "not-zip": b"\x00\x01\x02 neither html nor a zip checkpoint\n",
}


def _serve(box: A3Sandbox, body: bytes) -> dict[str, str]:
    f = box.sb.tmp / "served_body"
    f.write_bytes(body)
    return {"FAKE_CURL_BODY_FILE": str(f), "FAKE_CURL_CONTENT_TYPE": "text/html; charset=utf-8"}


def _dest_files(box: A3Sandbox) -> set[str]:
    d = box.db / "cascadia"
    return {p.name for p in d.iterdir()} if d.is_dir() else set()


@pytest.mark.parametrize("kind", sorted(_HTML_BODIES))
def test_fetch_resolve_rejects_a_page_instead_of_the_weights(box: A3Sandbox, kind: str) -> None:
    # S1: --resolve must never print (and so never pin) the sha256 of a Drive HTML page.
    body = _HTML_BODIES[kind]
    _assets_for_fetch(box, weights_sha="", unimod_sha="")
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", "--resolve", A3_WEIGHTS_HOST_OK="1", **_serve(box, body))
    assert res.returncode == 1, _out(res)
    assert sha256_hex(body) not in _out(res), "the page's sha256 was printed"
    assert _dest_files(box) <= {"unimod_tables.xml"}, _dest_files(box)


@pytest.mark.parametrize("kind", sorted(_HTML_BODIES))
def test_fetch_rejects_a_page_even_when_its_sha256_is_committed(box: A3Sandbox, kind: str) -> None:
    # S1, normal mode: a page pinned by mistake still never becomes the checkpoint.
    body = _HTML_BODIES[kind]
    _assets_for_fetch(box, weights_sha=sha256_hex(body))
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1", **_serve(box, body))
    assert res.returncode == 1, _out(res)
    assert not _weights(box).exists()
    assert _dest_files(box) <= {"unimod_tables.xml"}, _dest_files(box)


@pytest.mark.parametrize("value", ["0", "true", "yes", "11", " 1"])
def test_fetch_host_ok_must_be_exactly_1(box: A3Sandbox, value: str) -> None:
    # N1: only A3_WEIGHTS_HOST_OK=1 opens the weights download.
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK=value)
    assert res.returncode == 2, _out(res)
    assert box.pilot.fetch_calls() == []
    assert not _weights(box).exists()


def test_fetch_unimod_step_runs_without_the_weights_gate(box: A3Sandbox) -> None:
    # N2: the Unimod table is offline; only the Drive weights need ADR-M10-B.
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", unset=("A3_WEIGHTS_HOST_OK",))
    assert res.returncode == 2, _out(res)  # the weights precondition still fails
    assert box.pilot.fetch_calls() == []
    u = _unimod(box)
    assert u.read_bytes() == STUB_UNIMOD_XML
    assert stat.S_IMODE(u.stat().st_mode) == 0o444
    assert not _weights(box).exists()


def test_fetch_resolve_prints_the_unimod_sha_without_the_gate(box: A3Sandbox) -> None:
    _assets_for_fetch(box, weights_sha="", unimod_sha="")
    box.stub_cascadia_env()
    res = box.run("fetch_a3_assets", "--resolve", unset=("A3_WEIGHTS_HOST_OK",))
    assert res.returncode == 2, _out(res)
    assert sha256_hex(STUB_UNIMOD_XML) in res.stdout
    assert box.pilot.fetch_calls() == []


def test_fetch_leaves_another_jobs_part_files_alone(box: A3Sandbox) -> None:
    # N8: .part names are unique per job; a concurrent job's partial files are not touched.
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    dest = box.db / "cascadia"
    dest.mkdir(parents=True, exist_ok=True)
    other = {n: f"another job's partial {n}\n".encode()
             for n in (".cascadia_0.0.7.ckpt.part", ".unimod_tables.xml.part")}  # fmt: skip
    for n, data in other.items():
        (dest / n).write_bytes(data)
    res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1")
    assert res.returncode == 0, _out(res)
    for n, data in other.items():
        assert (dest / n).is_file() and (dest / n).read_bytes() == data, n
    assert _weights(box).read_bytes() == _payload(_weights_url())
    assert _unimod(box).read_bytes() == STUB_UNIMOD_XML


def _curl_output(call: dict[str, object]) -> str:
    argv = [str(a) for a in call["argv"]]  # type: ignore[union-attr]
    for i, a in enumerate(argv):
        if a in ("-o", "--output") and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith("--output="):
            return a.split("=", 1)[1]
    raise AssertionError(f"curl without --output: {argv}")


def test_fetch_part_names_differ_between_jobs(box: A3Sandbox) -> None:
    _assets_for_fetch(box)
    box.stub_cascadia_env()
    for job in ("1001", "1002"):
        res = box.run("fetch_a3_assets", A3_WEIGHTS_HOST_OK="1", SLURM_JOB_ID=job)
        assert res.returncode == 0, _out(res)
    outs = [_curl_output(c) for c in box.pilot.fetch_calls()]
    assert len(outs) == 2 and outs[0] != outs[1], outs
    dest = (box.db / "cascadia").resolve()
    assert all(Path(o).resolve().parent == dest for o in outs), outs  # same filesystem: mv


def test_mzml_pilot_checks_the_run_name_before_converting(box: A3Sandbox) -> None:
    # N3: a run folder that is not in acquisition.json is refused before an hour of conversion.
    box.make_run_dir("WT_9")
    res = box.run("mzml_pilot", "WT_9")
    assert res.returncode == 2, _out(res)
    assert _timsconvert_calls(box) == []
    assert not _converted(box, "WT_9").exists()
