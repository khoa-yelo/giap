"""``scripts/slurm/diann_image.sbatch``: verify Khoa's image and install a read-only copy (M10-02).

Spec: docs/specs/M10-02-diann-searches.md, "Image verification and the read-only copy",
examples 1-3, edge cases "sha mismatch" and "Image missing", acceptance criteria 1 and 11;
index lead decision 2; ADR-0051 "Image handling".

Every run happens in ``helpers.pilot_kit_m10.m10_sandbox``: the "shared image" is a 9-byte
stand-in (mode 0777, like the real one), the sandbox pin file holds its sha256, and ``apptainer``
is a fake that logs its calls. No container runs, the real image is never read or copied.

Conventions settled by the test author:

* ``DIANN_SIF`` is left unset in these tests, so the copy goes to its default
  ``$GIPPA_ROOT/images/diann_2.3.1.sif``.
* The install record ``$GIPPA_ROOT/images/diann_2.3.1.sif.json`` has a ``source`` object (``path``,
  ``size``, ``mode``, ``owner``, ``group``, ``mtime``, ``sha256``) and a ``copy`` object (``path``,
  ``size``, ``mode``, ``sha256``); modes are octal strings as ``stat -c %a`` prints them
  (``"777"``, ``"444"``; a leading ``0`` is accepted).
* "No read of the source" in example 2 is checked by making the source unreadable (mode 000) and
  by the ``cp`` wrapper's log.

No statistical claim is made here, so there is no calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from helpers.pilot_kit_m10 import (
    BANNER,
    DIANN_SIF_SHARED,
    M10_SCRIPTS,
    PIN_NAME,
    STAND_IN,
    STAND_IN_SHA256,
    M10Sandbox,
    m10_sandbox,
    sha256,
)

SCRIPT = M10_SCRIPTS["diann_image"]


@pytest.fixture
def m10(tmp_path: Path) -> M10Sandbox:
    m = m10_sandbox(tmp_path)
    m.sb.env.pop("DIANN_SIF", None)
    return m


def _install(m10: M10Sandbox, **env: str):  # type: ignore[no-untyped-def]
    return m10.run("diann_image", **env)


def _ok(res) -> None:  # type: ignore[no-untyped-def]
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _mode(value: object) -> int:
    return int(str(value), 8)


def _images(m10: M10Sandbox) -> list[str]:
    d = m10.root / "images"
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def _execs(m10: M10Sandbox) -> list[dict[str, object]]:
    return [c for c in m10.calls() if c["argv"] and c["argv"][0] == "exec"]


def _text(path: Path) -> str:
    assert path.is_file(), f"{path} does not exist (M10-02 deliverable 2)"
    return path.read_text()


# ------------------------------------------------------------ example 1: a fresh install


def test_example1_installs_a_read_only_copy_equal_to_the_source(m10: M10Sandbox) -> None:
    res = _install(m10)
    _ok(res)
    assert m10.copy.is_file(), f"no copy at {m10.copy}\n{res.stdout}\n{res.stderr}"
    assert m10.copy.read_bytes() == STAND_IN
    assert m10.copy.stat().st_mode & 0o777 == 0o444, oct(m10.copy.stat().st_mode)
    assert sha256(m10.copy) == STAND_IN_SHA256
    assert m10.shared.read_bytes() == STAND_IN, "the source must never be modified"
    leftovers = [n for n in _images(m10) if ".tmp" in n]
    assert not leftovers, f"temporary file left behind: {leftovers}"


def test_example1_the_images_directory_is_0755(m10: M10Sandbox) -> None:
    _ok(_install(m10))
    assert (m10.root / "images").stat().st_mode & 0o777 == 0o755


def test_example1_writes_the_install_record(m10: M10Sandbox) -> None:
    _ok(_install(m10))
    record = m10.root / "images" / f"{PIN_NAME}.json"
    assert record.is_file(), f"{record} was not written (rule 7)"
    doc = json.loads(record.read_text())
    src, cp = doc["source"], doc["copy"]
    assert src["path"] == str(m10.shared), src
    assert int(src["size"]) == len(STAND_IN)
    assert _mode(src["mode"]) == 0o777, src["mode"]
    for key in ("owner", "group", "mtime"):
        assert src.get(key), f"source.{key} missing: {src}"
    assert src["sha256"] == STAND_IN_SHA256
    assert cp["path"] == str(m10.copy), cp
    assert int(cp["size"]) == len(STAND_IN)
    assert _mode(cp["mode"]) == 0o444, cp["mode"]
    assert cp["sha256"] == STAND_IN_SHA256
    assert STAND_IN_SHA256 in json.dumps(doc.get("pin")), doc.get("pin")
    assert "diann_image.sha256" in json.dumps(doc), "the record must name the pin file"
    for key in ("installed_utc", "job_id", "hostname"):
        assert doc.get(key), f"{key} missing from the install record"
    assert str(doc["job_id"]) == m10.sb.env["SLURM_JOB_ID"]
    assert "khoayelo189/diann:2.3.1" in json.dumps(doc), "the apptainer inspect labels"
    assert "1.4.4" in json.dumps(doc), "the Apptainer version"
    banner = doc["banner"]
    first = banner[0] if isinstance(banner, list) else str(banner).splitlines()[0]
    assert BANNER in first, banner


def test_example1_checks_the_banner_on_the_copy(m10: M10Sandbox) -> None:
    _ok(_install(m10))
    execs = _execs(m10)
    assert execs, "the banner check (apptainer exec … --version) never ran"
    # review S3: the banner is checked on the temporary copy before it is renamed into place,
    # so the image run is the copy or a file beside it, with the pinned bytes
    for c in execs:
        assert Path(str(c["image"])).parent == m10.copy.parent, c["image"]
        assert c["image_sha256"] == STAND_IN_SHA256, "the banner check must run the hashed copy"
    assert all(c["image"] != str(m10.shared) for c in m10.calls()), "the source was executed"


def test_a_banner_without_the_version_of_record_is_exit_2(m10: M10Sandbox) -> None:
    res = _install(m10, FAKE_DIANN_BANNER="DIA-NN 1.9.2 (Data-Independent Acquisition)")
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert BANNER in res.stdout + res.stderr, "the message names the expected banner"


# ------------------------------------------------------------ example 2: idempotent


def test_example2_a_second_run_does_not_read_the_source(m10: M10Sandbox) -> None:
    _ok(_install(m10))
    before = m10.copy.stat()
    cps = len(m10.cp_calls())
    m10.shared.chmod(0o000)  # any read of the source now fails
    try:
        res = _install(m10)
    finally:
        m10.shared.chmod(0o777)
    _ok(res)
    assert "already installed" in res.stdout + res.stderr, res.stdout + res.stderr
    assert str(m10.copy) in res.stdout + res.stderr
    new_cps = m10.cp_calls()[cps:]
    assert not any(str(m10.shared) in a for c in new_cps for a in c), new_cps
    after = m10.copy.stat()
    assert (after.st_ino, after.st_mtime_ns) == (before.st_ino, before.st_mtime_ns)


# ------------------------------------------------------------ example 3: a changed source


def test_example3_a_source_that_no_longer_matches_is_exit_2(m10: M10Sandbox) -> None:
    m10.shared.write_bytes(b"a different image\n")
    other = sha256(m10.shared)
    res = _install(m10)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    msg = res.stdout + res.stderr
    assert other in msg and STAND_IN_SHA256 in msg, f"both hashes must be named:\n{msg}"
    assert "never fetch DIA-NN elsewhere" in msg, msg
    assert not m10.copy.exists(), "a mismatching image was installed"
    assert not [n for n in _images(m10) if n != f"{PIN_NAME}.json"], _images(m10)
    assert not _execs(m10), "apptainer exec ran on an image that failed the pin"


def test_the_copy_not_the_source_is_hashed(m10: M10Sandbox) -> None:
    """Rule 5: the hash that gates the install is the copy's. The cp wrapper shows a copy."""
    _ok(_install(m10))
    cps = m10.cp_calls()
    assert any(str(m10.shared) in c for c in cps), f"the source was not copied with cp: {cps}"
    assert any("--preserve=timestamps" in c for c in cps), cps
    tmp = [a for c in cps for a in c if ".tmp" in a]
    assert tmp and all(Path(t).parent == m10.root / "images" for t in tmp), tmp


# ------------------------------------------------------------ edge cases: preconditions


def test_a_missing_source_is_exit_2_naming_it_and_r6(m10: M10Sandbox) -> None:
    missing = m10.sb.tmp / "nowhere" / PIN_NAME
    res = _install(m10, DIANN_SIF_SHARED=str(missing))
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    msg = res.stdout + res.stderr
    assert str(missing) in msg, msg
    assert "never fetch DIA-NN elsewhere" in msg and "R6" in msg, msg
    assert not m10.copy.exists()
    assert not m10.calls(), "no container call without a verified copy"


def test_a_directory_as_source_is_exit_2(m10: M10Sandbox) -> None:
    d = m10.sb.tmp / "adir.sif"
    d.mkdir()
    res = _install(m10, DIANN_SIF_SHARED=str(d))
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert not m10.copy.exists()


@pytest.mark.parametrize(
    "content",
    [
        "",
        "1946452d  diann_2.3.1.sif\n",
        f"{STAND_IN_SHA256}  other.sif\n",
        f"{STAND_IN_SHA256}  diann_2.3.1.sif\n{STAND_IN_SHA256}  diann_2.3.1.sif\n",
        f"{STAND_IN_SHA256.upper()[:63]}  diann_2.3.1.sif\n",
    ],
    ids=["empty", "short-hash", "wrong-name", "two-lines", "63-hex"],
)
def test_a_malformed_pin_file_is_exit_2(m10: M10Sandbox, content: str) -> None:
    m10.pin_file.write_text(content)
    res = _install(m10)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert not m10.copy.exists()
    assert not _execs(m10)


def test_a_missing_pin_file_is_exit_2(m10: M10Sandbox) -> None:
    m10.pin_file.unlink()
    res = _install(m10)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "diann_image.sha256" in res.stdout + res.stderr
    assert not m10.copy.exists()


# ------------------------------------------------------------ an existing copy is never replaced


def test_an_existing_mismatched_copy_is_exit_2_and_kept(m10: M10Sandbox) -> None:
    m10.copy.parent.mkdir(parents=True)
    m10.copy.write_bytes(b"tampered\n")
    m10.copy.chmod(0o444)
    res = _install(m10)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    msg = res.stdout + res.stderr
    assert str(m10.copy) in msg and "does not match" in msg, msg
    assert m10.copy.read_bytes() == b"tampered\n", "the copy was replaced; it is evidence"
    assert not [c for c in m10.cp_calls() if str(m10.shared) in c], "the source was copied"


def test_an_existing_writable_copy_is_exit_2_and_kept(m10: M10Sandbox) -> None:
    m10.copy.parent.mkdir(parents=True)
    m10.copy.write_bytes(STAND_IN)
    m10.copy.chmod(0o644)
    res = _install(m10)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert m10.copy.stat().st_mode & 0o777 == 0o644, "the script must not fix the copy itself"


@pytest.mark.parametrize("mode", [0o444, 0o440, 0o400])
def test_an_existing_read_only_matching_copy_is_already_installed(
    m10: M10Sandbox, mode: int
) -> None:
    m10.copy.parent.mkdir(parents=True)
    m10.copy.write_bytes(STAND_IN)
    m10.copy.chmod(mode)
    res = _install(m10)
    _ok(res)
    assert "already installed" in res.stdout + res.stderr


def test_the_install_writes_nothing_into_the_repository_and_never_downloads(
    m10: M10Sandbox,
) -> None:
    _ok(_install(m10))
    assert m10.repo_is_clean() == "", m10.repo_is_clean()
    assert not m10.pilot.fetch_calls(), "R6: DIA-NN is never fetched"


# ------------------------------------------------------------ static: house rules (criterion 11)


def _directive(key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", _text(SCRIPT), re.M)
    return (m.group(1) or "") if m else None


def _code() -> str:
    return "\n".join(ln for ln in _text(SCRIPT).splitlines()
                     if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH"))  # fmt: skip


def test_house_rules_of_the_jobs_table() -> None:
    assert _directive("job-name") == "gippa-diann-image"
    assert _directive("partition") == "normal"
    assert _directive("qos") == "normal", "dev allows one job; M10-01/M10-08 may hold it"
    assert _directive("cpus-per-task") == "1"
    assert (_directive("mem") or "").upper().rstrip("B") in {"2G", "2048M"}
    assert _directive("time") in {"00:30:00", "30:00", "30"}
    out = _directive("output") or ""
    assert "/local/logs/" in out or "$GIPPA_ROOT" in out, out
    code = _code()
    assert "scripts/farmshare_env.sh" in code
    assert re.search(r"^set -euo pipefail$", code, re.M)


def test_the_script_reads_the_pin_file_and_defaults() -> None:
    code = _code()
    assert "validation/pilot/diann_image.sha256" in code, "the pin file is the only pin source"
    assert "1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381" not in code, (
        "the pin must not be hard-coded (lead decision 2: Q2)"
    )
    m = re.search(r"DIANN_SIF_SHARED:-([^}]+)\}", code)
    assert m and m.group(1).strip("\"'") == DIANN_SIF_SHARED, "the source default is Khoa's image"
    m = re.search(r"DIANN_SIF:-([^}]+)\}", code)
    assert m and m.group(1).strip("\"'") in {
        "$GIPPA_ROOT/images/diann_2.3.1.sif",
        "${GIPPA_ROOT}/images/diann_2.3.1.sif",
    }, m and m.group(1)
    assert re.search(r"\bmv\s+-n\b", code), "the rename must not overwrite (mv -n, rule 6)"
    assert re.search(r"chmod\s+0?444\b", code), "chmod 0444 (rule 6)"
    for tool in ("wget", "curl", "apptainer pull", "singularity pull", "docker"):
        assert tool not in code, f"R6: {tool} in diann_image.sbatch"
