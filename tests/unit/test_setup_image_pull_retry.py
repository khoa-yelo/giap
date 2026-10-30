"""#504: ``gippa setup`` image pulls retry once, keep the error detail, and say how to retry.

GitHub issue khoa-yelo/gippa#504. On Sherlock (apptainer 1.5.5, unprivileged
build) ``apptainer pull`` died in ``mksquashfs`` (exit 139) for three images. An immediate second
pull succeeded. The message kept only the last stderr line, and nothing told the user that
rerunning setup retries the images. Item 4 of the issue (progress output during long pulls) is out
of scope here.

Interface pinned here (readings chosen where the issue leaves room):

- ``gippa.setupdb.images._pull(runtime, tool, ref, version, cache)`` (the signature is unchanged)
  retries a failed pull once: at most two runs of the runtime per image.
  - The first pull fails and the second succeeds: the record is ``ok``, two runs.
  - Both fail: the record is not ``ok`` (``absent`` as today), two runs, no third.
  - The first succeeds: one run.
  - No partial apptainer image is left before the retry. Real ``apptainer pull`` refuses to
    overwrite an existing file, and the fake does the same.
- The same single retry applies to **docker**. ``_pull`` is shared by both runtimes, and a
  transient registry or daemon error is as likely there; a test pins the docker retry too.
- Not retried, since the outcome is deterministic: a full disk or quota (``is_disk_full``; the
  record stays ``failed`` with the advice of M11-15c1), a docker digest mismatch (``mismatch``),
  and a missing runtime binary (``OSError``).
- Any wait between the attempts goes through ``time.sleep`` (looked up at call time, so these
  tests patch it to a no-op). No wait length is pinned.
- Failure message (item 2): the message of a failed record keeps more than the last line.
  Either (i) the last 5 non-empty stderr lines are all in the message, or (ii) the message names an
  existing file holding the full stderr of the failed pull. The message is bounded: with 200
  stderr lines, at most 20 of them are copied into the message (the manifest stays small).
  Wording is not pinned.
- Summary (item 3): when at least one image is not ``ok`` at the end of ``gippa setup``, one
  ``gippa setup: ...`` line tells the user to rerun ``gippa setup`` to retry the failed images and
  names each failed tool. It comes after the last ``image:<tool>`` line and after every DB item
  line, so it is the last item-related output. There is no such line when every image is ``ok``.
  Setup still exits 2. An image that failed on a full disk or quota is not told to rerun (the same
  disk stays full): a line of its own names it and points to the messages above (review round).

No test runs a real container runtime or reaches a registry. Logging fakes (``apptainer``,
``docker``) sit first on ``PATH`` and play a scripted outcome per call: ``FAKE_OUTCOMES``, for
example ``"fail ok"``, holds the outcomes of calls 1, 2, ... (``ok`` after the list runs out), and
``FAKE_FAIL_FOR`` fails every pull of a matching ref.
"""

from __future__ import annotations

import json
import os
import re
import stat
import tempfile
import time
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.setup_opt_outs import NO_RESOURCE_FLAGS, opt_outs
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REF = "quay.io/biocontainers/genomad@sha256:" + "d5" * 32
BAD = "sha256:" + "ee" * 32
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models", "--no-uniop",
            *opt_outs(*NO_RESOURCE_FLAGS)]  # fmt: skip
TAIL = [
    "INFO:    Converting OCI blobs to SIF format",
    "INFO:    Starting build...",
    "detail-1: Parallel mksquashfs: Using 4 processors",
    "detail-2: Creating 4.0 filesystem on /tmp/build-temp-1/squashfs.img",
    "detail-3: [=====|            ] 120/900  13%",
    "detail-4: Segmentation fault (core dumped)",
    "detail-5: while creating squashfs",
    "FATAL:   While making image from oci registry: error fetching image to cache: while "
    "building SIF from layers: while creating squashfs: /usr/libexec/apptainer/bin/mksquashfs "
    "command failed: exit status 139:",
]

FAKE = r"""#!/bin/sh
prog=$(basename "$0")
printf '%s\n' "$prog $*" >> "$FAKE_RUNTIME_LOG"
n=0
[ -f "$FAKE_COUNT" ] && n=$(cat "$FAKE_COUNT")
n=$((n + 1))
echo "$n" > "$FAKE_COUNT"
outcome=$(echo "$FAKE_OUTCOMES" | awk -v i="$n" '{ print (i <= NF ? $i : "ok") }')
last=""; prev=""
for a in "$@"; do prev=$last; last=$a; done
if [ -n "$FAKE_FAIL_FOR" ]; then
    case "$last" in *"$FAKE_FAIL_FOR"*) outcome=fail;; esac
fi
if [ "$prog $1" = "apptainer pull" ] && [ -e "$prev" ]; then
    echo "FATAL:   Image file already exists: \"$prev\" - will not overwrite" >&2
    exit 255
fi
fail() {
    if [ -n "$FAKE_STDERR_FILE" ]; then cat "$FAKE_STDERR_FILE" >&2; else
        echo "fake $prog: pull failed" >&2; fi
}
case "$prog $1" in
"apptainer pull")
    case "$outcome" in
    ok) printf 'fake sif %s\n' "$last" > "$prev"; exit 0;;
    full) printf 'half a sif\n' > "$prev"
          echo "FATAL:   write $prev: no space left on device" >&2; exit 255;;
    *) printf 'half a sif\n' > "$prev"; fail; exit 255;;
    esac;;
"docker pull")
    case "$outcome" in
    ok) echo "Digest: ${last#*@}"; exit 0;;
    bad) echo "Digest: $FAKE_BAD_DIGEST"; exit 0;;
    full) echo "write /var/lib/docker/x: no space left on device" >&2; exit 1;;
    *) fail; exit 1;;
    esac;;
"docker run")
    exit 1;;
esac
echo "fake $prog: unexpected call: $*" >&2
exit 3
"""


# ---------------------------------------------------------------- helpers and fixtures


def _images() -> Any:
    from gippa.setupdb import images

    return images


def _install(bindir: Path) -> None:
    bindir.mkdir(parents=True, exist_ok=True)
    for exe in ("apptainer", "docker"):
        path = bindir / exe
        path.write_text(FAKE)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _point(work: Path, monkeypatch: pytest.MonkeyPatch, outcomes: str = "") -> None:
    """Fresh call log and counter under ``work``; the scripted outcomes."""
    work.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("FAKE_RUNTIME_LOG", str(work / "runtime.log"))
    monkeypatch.setenv("FAKE_COUNT", str(work / "count"))
    monkeypatch.setenv("FAKE_OUTCOMES", outcomes)


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Never the real $GIPPA_ROOT / $GIPPA_DB_DIR / caches / envs; fake runtimes first on PATH;
    ``time.sleep`` recorded, not slept."""
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_REPO", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "GIPPA_CONTAINER_RUNTIME", "FAKE_FAIL_FOR",
                "FAKE_STDERR_FILE", "FAKE_BAD_DIGEST"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no_mamba"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "cfg"))
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    _install(tmp_path / "fakebin")
    monkeypatch.setenv("PATH", f"{tmp_path / 'fakebin'}{os.pathsep}{os.environ.get('PATH', '')}")
    _point(tmp_path / "w", monkeypatch)
    slept: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda s: slept.append(float(s)))
    return slept


def _calls(prefix: str = "") -> list[list[str]]:
    log = Path(os.environ["FAKE_RUNTIME_LOG"])
    lines = log.read_text().splitlines() if log.exists() else []
    return [ln.split() for ln in lines if ln.startswith(prefix)]


def _stderr(tmp_path: Path, lines: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    f = tmp_path / "stderr.txt"
    f.write_text("".join(f"{ln}\n" for ln in lines))
    monkeypatch.setenv("FAKE_STDERR_FILE", str(f))


def _pull(runtime: str, cache: Path) -> Any:
    return _images()._pull(runtime, "genomad", REF, "1.11.2", cache)


def _full_stderr_kept(message: str, lines: list[str]) -> bool:
    """Reading (ii): the message names an existing file that holds every stderr line."""
    for m in re.finditer(r"(/[^\s'\"(),;]+)", message):
        p = Path(m.group(1).rstrip(".:"))
        if p.is_file():
            text = p.read_text(encoding="utf-8", errors="replace")
            if all(ln in text for ln in lines):
                return True
    return False


# ---------------------------------------------------------------- (a) one retry


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_pull_that_fails_once_then_succeeds_is_ok(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    _point(tmp_path / "w", monkeypatch, "fail ok")
    cache = tmp_path / "cache"
    rec = _pull(runtime, cache)
    assert len(_calls(f"{runtime} pull")) == 2, _calls()
    assert rec.status == "ok", rec
    if runtime == "apptainer":
        img = cache / _images().simple_name(REF)
        assert img.read_text().startswith("fake sif"), "the retry's image, not the partial one"
        assert rec.path == str(img)


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_pull_that_fails_twice_is_not_ok_and_not_tried_a_third_time(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    _point(tmp_path / "w", monkeypatch, "fail fail ok")
    cache = tmp_path / "cache"
    rec = _pull(runtime, cache)
    assert len(_calls(f"{runtime} pull")) == 2, _calls()
    assert rec.status != "ok", rec
    assert rec.status == "absent", rec  # as a single failed pull is today
    assert not (cache / _images().simple_name(REF)).exists(), "no partial image is left"


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_pull_that_succeeds_first_runs_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str, _isolated: list[float]
) -> None:
    _point(tmp_path / "w", monkeypatch, "ok")
    rec = _pull(runtime, tmp_path / "cache")
    assert rec.status == "ok", rec
    assert len(_calls(f"{runtime} pull")) == 1
    assert _isolated == [], "no wait without a failure"


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_full_disk_is_not_retried(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    """A full disk or quota stays full: one run, the M11-15c1 ``failed`` record."""
    _point(tmp_path / "w", monkeypatch, "full ok")
    cache = tmp_path / "cache"
    rec = _pull(runtime, cache)
    assert len(_calls(f"{runtime} pull")) == 1, _calls()
    assert rec.status == "failed", rec
    assert not (cache / _images().simple_name(REF)).exists()


def test_docker_digest_mismatch_is_not_retried(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_BAD_DIGEST", BAD)
    _point(tmp_path / "w", monkeypatch, "bad ok")
    rec = _pull("docker", tmp_path / "cache")
    assert len(_calls("docker pull")) == 1, _calls()
    assert rec.status == "mismatch", rec


def test_missing_runtime_binary_is_not_retried(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, _isolated: list[float]
) -> None:
    (tmp_path / "fakebin" / "apptainer").unlink()
    monkeypatch.setenv("PATH", str(tmp_path / "fakebin"))  # no apptainer anywhere
    rec = _pull("apptainer", tmp_path / "cache")
    assert rec.status == "absent", rec
    assert _isolated == [], "an OSError is not waited on and retried"


@settings(max_examples=30, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(
    runtime=st.sampled_from(["apptainer", "docker"]),
    outcomes=st.lists(st.sampled_from(["ok", "fail"]), min_size=1, max_size=4),
)
def test_retry_property(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str, outcomes: list[str]
) -> None:
    """For any scripted sequence: ok iff one of the first two runs succeeds; runs = 1 if the
    first succeeds, else 2 (never more)."""
    work = Path(tempfile.mkdtemp(dir=tmp_path))
    _point(work, monkeypatch, " ".join(outcomes))
    rec = _pull(runtime, work / "cache")
    first_two = [*outcomes, "ok"][:2]
    assert (rec.status == "ok") == ("ok" in first_two), (outcomes, rec)
    assert len(_calls(f"{runtime} pull")) == (1 if outcomes[0] == "ok" else 2), outcomes


# ---------------------------------------------------------------- (b) the failure message


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_failure_message_keeps_the_stderr_tail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    _stderr(tmp_path, TAIL, monkeypatch)
    _point(tmp_path / "w", monkeypatch, "fail fail")
    rec = _pull(runtime, tmp_path / "cache")
    assert rec.status != "ok"
    assert TAIL[-1].split(":", 1)[0] in rec.message  # the last line, as today
    inline = all(ln.strip() in rec.message for ln in TAIL[-5:])
    assert inline or _full_stderr_kept(rec.message, TAIL), (
        f"neither the last 5 stderr lines nor a log of the full stderr:\n{rec.message}"
    )
    assert REF in rec.message  # still says which pull failed


def test_failure_message_is_bounded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    lines = [f"noise-line-{i:03d}: layer {i} extracted" for i in range(200)] + [TAIL[-1]]
    _stderr(tmp_path, lines, monkeypatch)
    _point(tmp_path / "w", monkeypatch, "fail fail")
    rec = _pull("apptainer", tmp_path / "cache")
    copied = sum(1 for ln in lines if ln in rec.message)
    assert copied <= 20, f"{copied} stderr lines copied into the record message"
    assert copied >= 5 or _full_stderr_kept(rec.message, lines), rec.message


# ---------------------------------------------------------------- (c) the summary at the end


def _setup(args: list[str]) -> Any:
    from gippa.cli import cli

    return CliRunner().invoke(cli, ["setup", *args], catch_exceptions=False)


def _summary(output: str) -> list[int]:
    """Indices of the lines that tell the user to rerun ``gippa setup``."""
    return [i for i, ln in enumerate(output.splitlines())
            if re.search(r"(?i)\brerun\b|\bre-run\b|\brun\b.*\bagain\b", ln)
            and "gippa setup" in ln]  # fmt: skip


def _items(db: Path) -> dict[str, dict[str, Any]]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    return {item["name"]: item for item in data["items"]}


def test_failed_image_ends_setup_with_a_rerun_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_FAIL_FOR", "genomad")
    db = tmp_path / "db"
    res = _setup(["--db-dir", str(db), *SKIP_DBS, "--images", "all", "--profile", "apptainer",
                  "--yes"])  # fmt: skip
    assert res.exit_code == 2, res.output
    assert _items(db)["image:genomad"]["status"] != "ok"
    found = _summary(res.output)
    assert len(found) == 1, f"no single 'rerun gippa setup' summary line:\n{res.output}"
    lines = res.output.splitlines()
    line = lines[found[0]]
    assert line.startswith("gippa setup:"), line
    assert "genomad" in line, f"the summary names the failed image: {line}"
    last_item = max(i for i, ln in enumerate(lines) if "image:" in ln and i != found[0])
    assert found[0] > last_item, res.output


def test_rerun_summary_names_every_failed_image(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_FAIL_FOR", "sha256:")  # every pinned ref fails
    db = tmp_path / "db"
    res = _setup(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer", "--yes"])
    assert res.exit_code == 2, res.output
    failed = sorted(n.removeprefix("image:") for n, i in _items(db).items()
                    if n.startswith("image:") and i["status"] != "ok")  # fmt: skip
    assert failed
    found = _summary(res.output)
    assert len(found) == 1, res.output
    line = res.output.splitlines()[found[0]]
    missing = [t for t in failed if not re.search(rf"\b{re.escape(t)}\b", line)]
    assert not missing, f"summary does not name {missing}: {line}"


def test_rerun_summary_comes_after_the_db_items(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The issue's cascade: the bakta image fails, so ``bakta_light`` cannot be installed; the
    summary is still the last item-related line."""
    monkeypatch.setenv("FAKE_FAIL_FOR", "bakta")
    db = tmp_path / "db"
    args = ["--db-dir", str(db), "--bakta", "light", "--no-genomad", "--no-promotech-models",
            "--no-uniop", *opt_outs(*NO_RESOURCE_FLAGS), "--images", "--profile", "apptainer",
            "--yes"]  # fmt: skip
    res = _setup(args)
    assert res.exit_code == 2, res.output
    lines = res.output.splitlines()
    found = _summary(res.output)
    assert len(found) == 1, res.output
    bakta = [i for i, ln in enumerate(lines) if ln.startswith("gippa setup: bakta_light")]
    assert bakta, f"no bakta_light item line:\n{res.output}"
    assert found[0] > max(bakta), res.output
    assert "bakta" in lines[found[0]]


def test_no_rerun_summary_when_every_image_is_ok(tmp_path: Path) -> None:
    db = tmp_path / "db"
    res = _setup(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer", "--yes"])
    assert res.exit_code == 0, res.output
    assert _summary(res.output) == [], res.output


def test_a_pull_that_recovers_on_retry_needs_no_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End to end: the first pull of setup fails once, its retry works; setup exits 0."""
    _point(tmp_path / "w", monkeypatch, "fail ok")
    db = tmp_path / "db"
    res = _setup(["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", "apptainer", "--yes"])
    assert res.exit_code == 0, res.output
    assert all(i["status"] == "ok" for n, i in _items(db).items() if n.startswith("image:"))
    assert _summary(res.output) == [], res.output


# ---------------------------------------------------------------- review round: PR #506
#
# Interface pinned here, on top of the module docstring:
#
# - A destination image that appears during the wait (another ``gippa setup`` finished it) is
#   ``ok``: no second run, and the file stays as it is.
# - A failed apptainer pull that ran with ``--disable-cache`` (``APPTAINER_CACHEDIR`` unset or
#   empty) has a message naming ``APPTAINER_CACHEDIR``. The argv stays that of example 6.9.
# - In the summary, an image with status ``failed`` (full disk or quota) is on no rerun line; one
#   other ``gippa setup:`` line names it and says ``above``. ``absent`` keeps the rerun line.


def test_image_finished_by_another_setup_during_the_wait_is_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _point(tmp_path / "w", monkeypatch, "fail ok")
    cache = tmp_path / "cache"
    dest = cache / _images().simple_name(REF)
    theirs = "the image another gippa setup finished\n"
    waits: list[float] = []

    def sleep(seconds: float) -> None:
        waits.append(float(seconds))
        dest.write_text(theirs)

    monkeypatch.setattr(time, "sleep", sleep)
    rec = _pull("apptainer", cache)
    assert waits, "the first pull failed, so setup waited"
    assert len(_calls("apptainer pull")) == 1, f"pulled over a finished image: {_calls()}"
    assert rec.status == "ok", rec
    assert rec.path == str(dest), rec
    assert dest.read_text() == theirs, "the finished image is not replaced"


def _argv(cache: Path) -> list[str]:
    return ["apptainer", "pull", "--disable-cache", str(cache / _images().simple_name(REF)),
            f"docker://{REF}"]  # fmt: skip


@pytest.mark.parametrize("value", [None, ""], ids=["unset", "empty"])
def test_failed_pull_without_a_cache_mentions_apptainer_cachedir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, value: str | None
) -> None:
    """#504: every retry rebuilds from nothing under ``--disable-cache``; say what avoids that."""
    if value is not None:
        monkeypatch.setenv("APPTAINER_CACHEDIR", value)
    _point(tmp_path / "w", monkeypatch, "fail fail")
    cache = tmp_path / "cache"
    rec = _pull("apptainer", cache)
    assert rec.status == "absent", rec
    assert _calls("apptainer pull") == [_argv(cache), _argv(cache)], "the argv of example 6.9"
    assert "APPTAINER_CACHEDIR" in rec.message, rec.message


def test_failed_pull_with_a_cache_does_not_mention_apptainer_cachedir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Guard: the user already set it, and the argv has no ``--disable-cache``."""
    monkeypatch.setenv("APPTAINER_CACHEDIR", str(tmp_path / "blobs"))
    _point(tmp_path / "w", monkeypatch, "fail fail")
    rec = _pull("apptainer", tmp_path / "cache")
    assert rec.status == "absent", rec
    assert all("--disable-cache" not in c for c in _calls("apptainer pull")), _calls()
    assert "APPTAINER_CACHEDIR" not in rec.message, rec.message


def test_failed_docker_pull_does_not_mention_apptainer_cachedir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Guard: the hint is about apptainer's layer cache."""
    _point(tmp_path / "w", monkeypatch, "fail fail")
    rec = _pull("docker", tmp_path / "cache")
    assert rec.status == "absent", rec
    assert "APPTAINER_CACHEDIR" not in rec.message, rec.message


def _names(line: str, tool: str) -> bool:
    return re.search(rf"(?<![\w-]){re.escape(tool)}(?![\w-])", line) is not None


def _with_status(db: Path, status: str) -> list[str]:
    return sorted(n.removeprefix("image:") for n, i in _items(db).items()
                  if n.startswith("image:") and i["status"] == status)  # fmt: skip


def _setup_images(tmp_path: Path, runtime: str) -> tuple[Any, Path]:
    db = tmp_path / "db"
    res = _setup(["--db-dir", str(db), *SKIP_DBS, "--images", "all", "--profile", runtime,
                  "--cache-dir", str(tmp_path / "image-cache"), "--yes"])  # fmt: skip
    return res, db


@pytest.mark.parametrize("runtime", ["apptainer", "docker"])
def test_full_disk_summary_does_not_say_rerun(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runtime: str
) -> None:
    """The same disk stays full, so rerunning is the wrong advice; the item line said what to do."""
    _point(tmp_path / "w", monkeypatch, "full")
    res, db = _setup_images(tmp_path, runtime)
    assert res.exit_code == 2, res.output
    full = _with_status(db, "failed")
    assert full, res.output
    assert _summary(res.output) == [], f"rerun advice for a full disk:\n{res.output}"
    lines = res.output.splitlines()
    above = [i for i, ln in enumerate(lines) if ln.startswith("gippa setup:")
             and "above" in ln.lower() and all(_names(ln, t) for t in full)]  # fmt: skip
    assert len(above) == 1, f"no single line naming {full} and the messages above:\n{res.output}"
    last_item = max(i for i, ln in enumerate(lines) if "image:" in ln and i != above[0])
    assert above[0] > last_item, res.output


def test_mixed_failures_keep_the_rerun_hint_only_for_the_plain_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """First image: full disk (one run). Second: both runs fail (``absent``). The rest are ok."""
    _point(tmp_path / "w", monkeypatch, "full fail fail")
    res, db = _setup_images(tmp_path, "apptainer")
    assert res.exit_code == 2, res.output
    full, plain = _with_status(db, "failed"), _with_status(db, "absent")
    assert full and plain, res.output
    found = _summary(res.output)
    assert len(found) == 1, f"no single rerun line:\n{res.output}"
    lines = res.output.splitlines()
    rerun = lines[found[0]]
    assert all(_names(rerun, t) for t in plain), f"the rerun line must name {plain}: {rerun}"
    assert not [t for t in full if _names(rerun, t)], f"{full} hit a full disk: {rerun}"
    above = [ln for i, ln in enumerate(lines) if i != found[0] and ln.startswith("gippa setup:")
             and "above" in ln.lower() and all(_names(ln, t) for t in full)]  # fmt: skip
    assert len(above) == 1, f"no single line naming {full} and the messages above:\n{res.output}"
