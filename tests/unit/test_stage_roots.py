"""``gippa.runner.local.stage_roots`` (M4-15 lead decisions 19-21; GIPPA_SPEC.md §6.4).

A root file is copied to ``workdir/<key>``. Only for the key ``preflight``, whose
``workdir/preflight`` is preflight's own work directory, does the file go inside that directory
(decision 19); for any other key an existing directory at ``workdir/<key>`` still fails loudly
(decision 21) instead of silently receiving the file.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _src(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / "src" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p


def test_file_root_copied_to_key(tmp_path: Path) -> None:
    from gippa.runner.local import stage_roots

    work = tmp_path / "work"
    stage_roots(work, {"abundance": _src(tmp_path, "abundance.norm.parquet", "x")})
    assert (work / "abundance").read_text(encoding="utf-8") == "x"


def test_preflight_into_existing_directory(tmp_path: Path) -> None:
    from gippa.runner.local import stage_roots

    work = tmp_path / "work"
    (work / "preflight").mkdir(parents=True)
    (work / "preflight" / "genome.norm.fa").write_text(">c\nA\n", encoding="utf-8")
    stage_roots(work, {"preflight": _src(tmp_path, "preflight.json", '{"a": 1}\n')})
    assert (work / "preflight" / "preflight.json").read_text(encoding="utf-8") == '{"a": 1}\n'
    assert (work / "preflight" / "genome.norm.fa").is_file()  # nothing else touched


@pytest.mark.parametrize("key", ["abundance", "genome", "cog2020_def"])
def test_other_key_onto_directory_fails_loudly(tmp_path: Path, key: str) -> None:
    from gippa.runner.local import stage_roots

    work = tmp_path / "work"
    (work / key).mkdir(parents=True)
    with pytest.raises((IsADirectoryError, ValueError, OSError)):
        stage_roots(work, {key: _src(tmp_path, "file.dat", "x")})
    assert list((work / key).iterdir()) == []
