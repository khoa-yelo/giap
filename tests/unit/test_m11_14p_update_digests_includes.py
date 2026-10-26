"""``update_digests.py`` finds the conda-off selector through ``includeConfig`` (#374 review S1).

Nested includes resolve against the including file's directory, a cycle ends, a missing target
and zero or two selectors across the files are errors.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
SELECTOR = "process {\n    withName: 'BAKTA' { conda = '' }\n}\n"


@pytest.fixture(scope="module")
def ud() -> ModuleType:
    path = REPO / "scripts" / "release" / "update_digests.py"
    spec = importlib.util.spec_from_file_location("update_digests_includes", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _tree(tmp: Path, files: dict[str, str]) -> Path:
    for rel, text in files.items():
        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp / rel).write_text(text)
    return tmp / "nextflow.config"


def test_nested_include_resolves_against_the_including_file(ud: ModuleType, tmp_path: Path) -> None:
    nf = _tree(
        tmp_path,
        {
            "nextflow.config": "profiles {\n docker {\n  includeConfig 'conf/a.config'\n }\n}\n",
            "conf/a.config": "includeConfig 'inner/b.config' // the shared block\n",
            "conf/inner/b.config": SELECTOR,
            "inner/b.config": "// a decoy at the wrong level\n",
        },
    )
    assert ud.selector_file(nf) == (tmp_path / "conf/inner/b.config").resolve()


def test_a_cycle_ends(ud: ModuleType, tmp_path: Path) -> None:
    nf = _tree(
        tmp_path,
        {
            "nextflow.config": "includeConfig 'conf/a.config'\n",
            "conf/a.config": "includeConfig 'b.config'\n" + SELECTOR,
            "conf/b.config": "includeConfig 'a.config'\n",
        },
    )
    assert [p.name for p in ud.config_closure(nf)] == ["nextflow.config", "a.config", "b.config"]
    assert ud.selector_file(nf).name == "a.config"


def test_a_missing_target_is_an_error(ud: ModuleType, tmp_path: Path) -> None:
    nf = _tree(tmp_path, {"nextflow.config": "includeConfig 'conf/nope.config'\n" + SELECTOR})
    with pytest.raises(ud.DigestError, match="not found"):
        ud.selector_file(nf)


@pytest.mark.parametrize("where", [(), ("nextflow.config", "conf/a.config")])
def test_zero_or_one_per_file_in_two_files_is_an_error(
    ud: ModuleType, tmp_path: Path, where: tuple[str, ...]
) -> None:
    files = {"nextflow.config": "includeConfig 'conf/a.config'\n", "conf/a.config": "\n"}
    for rel in where:
        files[rel] += SELECTOR
    nf = _tree(tmp_path, files)
    with pytest.raises(ud.DigestError, match=f"found {len(where)}"):
        ud.selector_file(nf)


def test_a_commented_include_is_not_followed(ud: ModuleType, tmp_path: Path) -> None:
    nf = _tree(tmp_path, {"nextflow.config": "// includeConfig 'conf/nope.config'\n" + SELECTOR})
    assert ud.selector_file(nf) == nf.resolve()
