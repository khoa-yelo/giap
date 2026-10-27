"""``workflow/main.nf``'s ``functionResources()`` resolves the five items at Python's paths (M8-06
item 7, 06b; the MVP-01 pattern of ``test_nf_main_installed_items.py``).

``main.nf`` cannot import ``gippa.setupdb``; ``functionResources()`` builds a ``[key, path]``
channel from ``installed(name, relative)`` calls with literal relative paths. These checks keep
the literals equal to ``gippa/data/db_sources.yaml`` and the names equal to
``gippa.setupdb.installed.FUNCTION_RESOURCES``; ``cogTable()``'s path is pinned with them. The
module switch (``'function' in (params.modules ?: [])``) is pinned textually here and behaviourally
by ``workflow/tests/main_function.nf.test``. No Nextflow needed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from gippa.setupdb.items import SOURCES

ROOT = Path(__file__).resolve().parents[2]
MAIN_NF = ROOT / "workflow" / "main.nf"
FUNCTION_NF = ROOT / "workflow" / "subworkflows" / "local" / "function.nf"
NAMES = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun")


def _text(path: Path) -> str:
    return re.sub(r"(?m)//.*$", "", path.read_text(encoding="utf-8"))


def _body(text: str, name: str) -> str:
    """The brace-balanced body of ``def <name>(...) { ... }``."""
    m = re.search(rf"def\s+{name}\s*\([^)]*\)\s*\{{", text)
    assert m is not None, f"main.nf defines no {name}()"
    depth, i = 1, m.end()
    while depth and i < len(text):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        i += 1
    return text[m.end() : i - 1]


def _installed(body: str) -> dict[str, str]:
    return dict(re.findall(r"installed\(\s*'([\w-]+)'\s*,\s*'([^']+)'\s*\)", body))


def test_function_resources_names_equal_installed_function_resources() -> None:
    from gippa.setupdb.installed import FUNCTION_RESOURCES

    assert tuple(NAMES) == FUNCTION_RESOURCES
    calls = _installed(_body(_text(MAIN_NF), "functionResources"))
    assert sorted(calls) == sorted(NAMES), calls


@pytest.mark.parametrize("name", NAMES)
def test_function_resource_path_equals_db_sources(name: str) -> None:
    calls = _installed(_body(_text(MAIN_NF), "functionResources"))
    assert calls.get(name) == SOURCES[name].path, (
        f"functionResources() resolves {name} at {calls.get(name)}, "
        f"db_sources.yaml at {SOURCES[name].path}"
    )


def test_cog_table_path_equals_db_sources() -> None:
    calls = _installed(_body(_text(MAIN_NF), "cogTable"))
    assert calls == {"cog2020_def": SOURCES["cog2020_def"].path}, calls


def test_main_nf_calls_function_only_under_the_module_switch() -> None:
    """Decision 1 / Q5: the first module switch in workflow/, for ``function`` only."""
    text = _text(MAIN_NF)
    assert re.search(r"include\s*\{\s*FUNCTION\s*\}\s*from\s*'./subworkflows/local/function'",
                     text), "main.nf does not include FUNCTION"  # fmt: skip
    switch = re.search(r"if\s*\(\s*'function'\s+in\s+\(?\s*params\.modules", text)
    assert switch is not None, "main.nf has no `'function' in params.modules` switch"
    call = re.search(r"\bFUNCTION\s*\(", text[switch.end() :])
    assert call is not None, "FUNCTION is not called after the switch"
    between = text[switch.end() : switch.end() + call.start()]
    assert "\n\n" not in between and between.count("\n") <= 2, "FUNCTION( is not the switched call"
    assert "functionResources()" in text[switch.end() : switch.end() + call.end() + 200]
    assert len(re.findall(r"\bFUNCTION\s*\(", text)) == 1, "FUNCTION called more than once"
    for mod in ("mobile", "bgc"):  # Q5: not switched in main.nf yet
        assert f"'{mod}' in" not in text, mod


def test_function_subworkflow_uses_the_function_step_alias() -> None:
    text = _text(FUNCTION_NF)
    assert re.search(r"include\s*\{\s*GIPPA_STEP\s+as\s+FUNCTION_STEP\s*\}", text), text[:400]
    assert re.search(r"workflow\s+FUNCTION\s*\{", text)
    assert re.search(r"stepIn\(\s*'function'", text)
