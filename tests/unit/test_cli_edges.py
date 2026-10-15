"""CLI and config edge cases (GIAP_SPEC.md §5.1, §5.2)."""

from __future__ import annotations

from pathlib import Path

import pytest
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools; decision 3)

from giap.exitcodes import ExitCode


def _exit_code(argv: list[str]) -> int:
    from giap.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    return 0 if code is None else int(code)


def _setup_argv(tmp_path: Path, *extra: str) -> list[str]:
    return ["setup", "--db-dir", str(tmp_path / "db"), *extra]


def test_unexpected_exception_becomes_exit_5_with_report_message(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    import giap.cli

    def boom(*_: object, **__: object) -> None:
        raise RuntimeError("kaboom")

    monkeypatch.setattr(giap.cli, "resolve_config", boom)
    assert _exit_code(_setup_argv(tmp_path)) == ExitCode.INTERNAL == 5
    err = capsys.readouterr().err
    assert "report" in err.lower()
    assert "internal error" in err.lower()


# `giap setup` is real since M2-08: tests that only check option parsing skip every resource.
_NO_DOWNLOADS = ["--bakta", "skip", "--no-genomad", "--no-antismash", "--no-promotech-models"]


def test_empty_config_file_is_accepted(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    from giap.config import resolve_config

    cfg = tmp_path / "empty.yaml"
    cfg.write_text("")
    assert resolve_config(cfg) == resolve_config()
    # Accepted by the CLI (not invalid input, 1); every resource skipped, so no download.
    argv = [*_setup_argv(tmp_path, "--config", str(cfg)), *_NO_DOWNLOADS]
    assert _exit_code(argv) == ExitCode.SUCCESS


def test_malformed_yaml_config_exits_1(tmp_path: Path) -> None:
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("threads: [1, 2\nseed: : :\n")
    assert _exit_code(_setup_argv(tmp_path, "--config", str(cfg))) == ExitCode.INVALID_INPUT == 1


def test_non_mapping_config_exits_1(tmp_path: Path) -> None:
    cfg = tmp_path / "list.yaml"
    cfg.write_text("- a\n- b\n")
    assert _exit_code(_setup_argv(tmp_path, "--config", str(cfg))) == ExitCode.INVALID_INPUT


def test_set_into_scalar_key_exits_1(tmp_path: Path) -> None:
    assert _exit_code(_setup_argv(tmp_path, "--set", "seed.x=1")) == ExitCode.INVALID_INPUT


def test_set_type_error_exits_1(tmp_path: Path) -> None:
    assert _exit_code(_setup_argv(tmp_path, "--set", "threads=abc")) == ExitCode.INVALID_INPUT


def test_set_without_equals_exits_1(tmp_path: Path) -> None:
    assert _exit_code(_setup_argv(tmp_path, "--set", "threads")) == ExitCode.INVALID_INPUT


def test_set_nested_key_is_accepted(tmp_path: Path) -> None:
    from giap.config import resolve_config

    cfg = resolve_config(None, ["annotation.promoters.mode=targeted"])
    assert cfg.annotation.promoters.mode == "targeted"
    argv = [*_setup_argv(tmp_path, "--set", "annotation.promoters.mode=targeted"), *_NO_DOWNLOADS]
    assert _exit_code(argv) == ExitCode.SUCCESS


def test_set_value_may_contain_equals() -> None:
    from giap.config import resolve_config

    assert resolve_config(None, ["sample_name=a=b"]).sample_name == "a=b"


def test_set_malformed_yaml_value_exits_1(tmp_path: Path) -> None:
    assert _exit_code(_setup_argv(tmp_path, "--set", "threads=[1")) == ExitCode.INVALID_INPUT
