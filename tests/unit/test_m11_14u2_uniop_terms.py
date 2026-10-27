"""M11-14u2: the UniOP terms notice, the TTY prompt, ``--accept-uniop-terms`` and the exit codes.

docs/specs/M11-14-container-install.md item 1b ("Terms notice (always printed, one line; 14u2)",
"Default: not accepted unless the user says so", "Exit codes"), examples 18, 18b, 18c; AC 9 (14u2:
"the notice is printed on every path, and the prompt default is no"); Khoa 2026-10-26 (Q-K27): an
explicit ``setup --uniop`` installs UniOP only, ``setup --images`` still includes it.

Interface pinned here (readings chosen where the spec leaves room; listed in the PR):

- The notice is one line, printed (stdout or stderr) whenever the UniOP step runs for a build
  runtime (docker, apptainer, micromamba), accepted or not.
- TTY detection goes through ``giap.setupdb.uniop_build.stdin_is_tty() -> bool`` (a module
  function, looked up at call time, so a test can stand in for a terminal; ``CliRunner`` has
  none). On a TTY without the flag, setup prints ``Build UniOP for your own use? [y/N]`` after the
  notice and reads one answer from stdin: ``y``/``yes`` builds with ``terms: "prompt"``; ``n``,
  Enter or anything else is a refusal. With ``--accept-uniop-terms`` there is no prompt and
  ``terms`` is ``"flag"``. Without a TTY there is never a prompt.
- A refusal or a missing acceptance: ``image:uniop`` ``absent`` with ``terms not accepted; rerun
  with --accept-uniop-terms``; exit 2 under ``--uniop``; under ``--images`` the other images are
  installed, a ``warning`` line names the message, exit 0.
- Q-K27: ``setup --uniop`` calls ``run_setup`` (if at all) with an ``SetupOptions`` whose
  ``plan()`` holds nothing but ``uniop_source`` (no Bakta, geNOMAD or Promotech models), even
  without ``--bakta skip --no-genomad --no-promotech-models``.
"""

from __future__ import annotations

import contextlib
import os
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.uniop_14u1 import (
    ACCEPT,
    SKIP_DBS,
    FakeNet,
    docker_calls,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    pin_uniop_source,
    setup_cli,
    sha256_hex,
)
from helpers.uniop_14u2 import (
    NOTICE,
    PROMPT,
    TERMS_NOT_ACCEPTED,
    calls,
    install_fake_runtimes,
    prebuilt_root,
)

pytestmark = pytest.mark.usefixtures("no_network")


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    install_fake_runtimes(tmp_path, monkeypatch)
    old = os.umask(0o022)
    yield
    os.umask(old)


@pytest.fixture
def net(monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


@pytest.fixture
def archive(net: FakeNet, monkeypatch: pytest.MonkeyPatch) -> bytes:
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    return data


def _tty(monkeypatch: pytest.MonkeyPatch, is_tty: bool = True) -> None:
    from giap.setupdb import uniop_build

    monkeypatch.setattr(uniop_build, "stdin_is_tty", lambda: is_tty, raising=False)


def _invoke(args: list[str], answer: str | None = None) -> tuple[int, str]:
    """``giap setup <args>`` with ``answer`` on stdin (``None``: empty stdin)."""
    from giap.cli import cli

    res = CliRunner().invoke(cli, ["setup", *args], input=answer)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    text = res.output
    with contextlib.suppress(ValueError):
        text += res.stderr
    return res.exit_code, text


def _uniop_args(db: Path, *extra: str, profile: str = "docker") -> list[str]:
    return ["--db-dir", str(db), *SKIP_DBS, "--uniop", "--profile", profile, *extra]


def _images_args(db: Path, *extra: str, profile: str = "docker") -> list[str]:
    return ["--db-dir", str(db), *SKIP_DBS, "--images", "--profile", profile, *extra]


# ---------------------------------------------------------------- the notice


def test_notice_is_one_line_naming_the_commit_and_the_package() -> None:
    assert "\n" not in NOTICE and "@3ed14dfc" in NOTICE


@pytest.mark.parametrize("profile", ["docker", "apptainer", "conda"])
@pytest.mark.parametrize("accept", [True, False], ids=["flag", "no-flag"])
def test_notice_printed_on_every_build_path(
    archive: bytes, tmp_path: Path, profile: str, accept: bool
) -> None:
    db = tmp_path / "db"
    extra = [ACCEPT] if accept else []
    code, out = setup_cli(["--cache-dir", str(tmp_path / "cache"),
                           *_uniop_args(db, *extra, profile=profile)])  # fmt: skip
    assert code == (0 if accept else 2), out
    assert NOTICE in out, out
    assert out.count(NOTICE) == 1


def test_notice_printed_under_images_without_the_flag(archive: bytes, tmp_path: Path) -> None:
    code, out = setup_cli(_images_args(tmp_path / "db"))
    assert code == 0, out
    assert NOTICE in out, out


# ---------------------------------------------------------------- examples 18, 18b (no TTY)


def test_example_18_images_without_tty_or_flag(
    archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(_images_args(db))
    assert code == 0, out
    assert net.contacted == [] and docker_calls(tmp_path, "build") == []
    items = manifest_items(db)
    assert items["image:uniop"]["status"] == "absent"
    assert items["image:uniop"]["message"] == TERMS_NOT_ACCEPTED
    assert items["image:bakta"]["status"] == "ok", "other images pulled"
    assert "warning" in out.lower() and TERMS_NOT_ACCEPTED in out, out
    assert PROMPT not in out, "never a prompt without a TTY"


def test_example_18b_uniop_without_tty_or_flag(
    archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(_uniop_args(db))
    assert code == 2, out
    assert TERMS_NOT_ACCEPTED in out, out
    assert PROMPT not in out
    assert net.contacted == [] and docker_calls(tmp_path, "build") == []
    assert manifest_items(db)["image:uniop"]["message"] == TERMS_NOT_ACCEPTED


def test_no_tty_never_reads_stdin(archive: bytes, tmp_path: Path,
                                  monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    """A 'y' waiting on a non-TTY stdin is not an acceptance (no silent yes in a script)."""
    _tty(monkeypatch, False)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_args(db), answer="y\n")
    assert code == 2, out
    assert docker_calls(tmp_path, "build") == []
    assert manifest_items(db)["image:uniop"]["message"] == TERMS_NOT_ACCEPTED


# ---------------------------------------------------------------- example 18c (TTY)


@pytest.mark.parametrize("answer", ["n\n", "\n", "N\n", "whatever\n", ""],
                         ids=["n", "enter", "N", "other", "eof"])  # fmt: skip
def test_example_18c_tty_refusal_under_uniop_exits_2(
    archive: bytes, net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, answer: str
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_args(db), answer=answer)
    assert code == 2, out
    assert PROMPT in out, out
    assert out.index(NOTICE) < out.index(PROMPT)
    assert net.contacted == [] and docker_calls(tmp_path, "build") == []
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "absent" and item["message"] == TERMS_NOT_ACCEPTED


@pytest.mark.parametrize("answer", ["n\n", "\n"], ids=["n", "enter"])
def test_example_18c_tty_refusal_under_images_exits_0(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, answer: str
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_images_args(db), answer=answer)
    assert code == 0, out
    assert PROMPT in out
    items = manifest_items(db)
    assert items["image:uniop"]["message"] == TERMS_NOT_ACCEPTED
    assert items["image:bakta"]["status"] == "ok"
    assert docker_calls(tmp_path, "build") == []


@pytest.mark.parametrize("answer", ["y\n", "yes\n", "Y\n"])
def test_example_18c_tty_yes_builds_with_terms_prompt(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, answer: str
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_args(db), answer=answer)
    assert code == 0, out
    assert len(docker_calls(tmp_path, "build")) == 1
    item = manifest_items(db)["image:uniop"]
    assert item["status"] == "ok" and item["local_build"]["terms"] == "prompt", item


def test_tty_with_the_flag_does_not_prompt(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_args(db, ACCEPT), answer="n\n")
    assert code == 0, out
    assert PROMPT not in out and NOTICE in out
    assert manifest_items(db)["image:uniop"]["local_build"]["terms"] == "flag"


# ---------------------------------------------------------------- after acceptance: exit 2 in both


def test_accepted_build_failure_exits_2_under_images(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_BUILD_FAIL", "ERROR: fixture boom")
    db = tmp_path / "db"
    code, out = setup_cli(_images_args(db, ACCEPT))
    assert code == 2, out
    assert manifest_items(db)["image:bakta"]["status"] == "ok", "the others installed first"


# ---------------------------------------------------------------- Q-K27


def _record_run_setup(monkeypatch: pytest.MonkeyPatch) -> list[set[str]]:
    import giap.setupdb
    from giap.exitcodes import ExitCode

    planned: list[set[str]] = []

    def fake(opts: Any, runner: Any, fetcher: Any) -> ExitCode:
        planned.append({i.name for i in giap.setupdb.plan(opts)})
        return ExitCode.SUCCESS

    monkeypatch.setattr(giap.setupdb, "run_setup", fake)
    return planned


def test_qk27_uniop_installs_uniop_only(
    archive: bytes, net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Khoa 2026-10-26: ``setup --uniop`` plans no Bakta/geNOMAD/Promotech download, even with
    the default ``--bakta light --genomad --promotech-models``."""
    planned = _record_run_setup(monkeypatch)
    db = tmp_path / "db"
    code, out = setup_cli(["--db-dir", str(db), "--uniop", "--profile", "docker", ACCEPT])
    assert code == 0, out
    assert all(p <= {"uniop_source"} for p in planned), planned
    assert {h for h, _ in net.contacted} <= {"github.com", "codeload.github.com"}
    assert docker_calls(tmp_path, "run") == [], "no ToolRunner download"
    names = set(manifest_items(db))
    assert names <= {"image:uniop", "uniop_source"}, names


def test_qk27_uniop_only_real_run_setup(archive: bytes, net: FakeNet, tmp_path: Path) -> None:
    """The same through the real ``run_setup``: only GitHub is contacted and only UniOP is
    recorded."""
    db = tmp_path / "db"
    code, out = setup_cli(["--db-dir", str(db), "--uniop", "--profile", "docker", ACCEPT])
    assert code == 0, out
    assert {h for h, _ in net.contacted} <= {"github.com", "codeload.github.com"}, net.contacted
    assert set(manifest_items(db)) <= {"image:uniop", "uniop_source"}


def test_qk27_images_still_includes_uniop(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(_images_args(db, ACCEPT))
    assert code == 0, out
    assert manifest_items(db)["image:uniop"]["status"] == "ok"


# ---------------------------------------------------------------- prebuilt: no prompt needed


def test_prebuilt_records_empty_terms(
    net: FakeNet, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIAP_ROOT", str(prebuilt_root(tmp_path / "giap_root")))
    db = tmp_path / "db"
    code, out = setup_cli(_uniop_args(db, ACCEPT, profile="farmshare,apptainer"))
    assert code == 0, out
    assert manifest_items(db)["image:uniop"]["local_build"]["terms"] == ""
    assert calls(tmp_path, "micromamba") == []
