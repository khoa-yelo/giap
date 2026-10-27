"""The generated references walk the real sources (M11-07 acceptance 3, examples 1, 3, 4).

Every ``gippa`` subcommand and non-hidden option is in ``cli.md``; every ``defaults.yaml`` key in
``configuration.md``; every datapackage resource and every ``schemas/tables/*.schema.json`` in
``tables.md``; the §9.4 vocabulary in ``flags.md``; and the committed pages are current
(``gen_reference.py --check`` exits 0).
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path
from types import ModuleType
from typing import Any

import click
import markdown
import pytest
import yaml

from gippa.cli import cli

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "docs" / "gen_reference.py"
SCHEMAS = REPO / "src" / "gippa" / "schemas"
DEFAULTS = REPO / "src" / "gippa" / "config" / "defaults.yaml"
DOCS = REPO / "docs"


def _load() -> ModuleType:
    assert SCRIPT.is_file(), "scripts/docs/gen_reference.py does not exist (M11-07 item 3)"
    spec = importlib.util.spec_from_file_location("gen_reference_real", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gen_reference_real"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def gr() -> ModuleType:
    return _load()


class _Cells(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.rows: list[list[str]] = []
        self._cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: Any) -> None:
        if tag == "tr":
            self.rows.append([])
        elif tag == "td":
            self._cell = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "td" and self._cell is not None:
            self.rows[-1].append("".join(self._cell).strip())
            self._cell = None

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell.append(data)


def _rows(md: str) -> list[list[str]]:
    p = _Cells()
    p.feed(markdown.markdown(md, extensions=["tables"]))
    return [r for r in p.rows if r]


def _sections(md: str) -> dict[str, str]:
    """``## heading`` → body up to the next ``## `` heading."""
    out: dict[str, str] = {}
    head: str | None = None
    buf: list[str] = []
    for line in md.splitlines():
        if line.startswith("## "):
            if head is not None:
                out[head] = "\n".join(buf) + "\n"
            head, buf = line[3:].strip(), []
        elif head is not None:
            buf.append(line)
    if head is not None:
        out[head] = "\n".join(buf) + "\n"
    return out


# --- cli.md -------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cli_md(gr: ModuleType) -> str:
    text: str = gr.render_cli(cli, "gippa")
    return text


def _visible_commands() -> list[str]:
    ctx = click.Context(cli)
    names = cli.list_commands(ctx)
    return [n for n in names if not cli.commands[n].hidden]


def test_example1_sections(cli_md: str) -> None:
    heads = [h for h in _sections(cli_md) if h.startswith("gippa ")]
    expected = ["compare", "convert", "makedb", "preflight", "report", "run", "setup"]
    assert set(expected) <= {h.removeprefix("gippa ") for h in heads}
    assert heads == [f"gippa {n}" for n in _visible_commands()]  # command order


def test_hidden_commands_absent(cli_md: str) -> None:
    heads = set(_sections(cli_md))
    for name, cmd in cli.commands.items():
        if cmd.hidden:
            assert f"gippa {name}" not in heads


@pytest.mark.parametrize("name", _visible_commands())
def test_every_visible_option_documented(cli_md: str, name: str) -> None:
    body = _sections(cli_md)[f"gippa {name}"]
    for p in cli.commands[name].params:
        if isinstance(p, click.Option):
            if p.hidden:
                continue
            for opt in [*p.opts, *p.secondary_opts]:
                assert f"`{opt}`" in body, (name, opt)
        else:
            assert f"`{p.human_readable_name}`" in body, (name, p.name)


def test_example1_run_rows(cli_md: str) -> None:
    rows = {r[0]: r for r in _rows(_sections(cli_md)["gippa run"])}
    assert rows["--genome"][1:3] == ["FILE", "required"]
    assert rows["--profile"][2] == "farmshare,apptainer"


def test_cli_unicode_written_as_is(cli_md: str) -> None:
    for cmd in cli.commands.values():
        for p in cmd.params:
            text = getattr(p, "help", None) or ""
            if any(ord(c) > 127 for c in text):
                assert text.replace("|", "\\|") in cli_md


# --- configuration.md ---------------------------------------------------------------------------


def _leaves(d: dict[str, Any], prefix: str = "") -> list[tuple[str, Any]]:
    out: list[tuple[str, Any]] = []
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict) and v:
            out.extend(_leaves(v, key + "."))
        else:
            out.append((key, v))
    return out


def _resolve(node: dict[str, Any], defs: dict[str, Any]) -> dict[str, Any]:
    while True:
        if "$ref" in node:
            node = defs[node["$ref"].rsplit("/", 1)[-1]]
        elif "allOf" in node and len(node["allOf"]) == 1:
            node = node["allOf"][0]
        else:
            return node


def _schema_for(schema: dict[str, Any], dotted: str) -> dict[str, Any] | None:
    defs = schema.get("$defs", {})
    node: dict[str, Any] = schema
    for part in dotted.split("."):
        node = _resolve(node, defs)
        props = node.get("properties", {})
        if part not in props:
            return None
        node = props[part]
    return node


@pytest.fixture(scope="module")
def config_md(gr: ModuleType) -> str:
    defaults = yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))
    schema = json.loads((SCHEMAS / "config.schema.json").read_text(encoding="utf-8"))
    text: str = gr.render_configuration(defaults, schema)
    return text


def test_every_defaults_key_documented(config_md: str) -> None:
    defaults = yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))
    schema = json.loads((SCHEMAS / "config.schema.json").read_text(encoding="utf-8"))
    rows = {r[0]: r for r in _rows(config_md)}
    missing = [k for k, _ in _leaves(defaults) if k not in rows]
    assert missing == []
    for key, _ in _leaves(defaults):
        node = _schema_for(schema, key)
        desc = (node or {}).get("description")
        assert rows[key][3] == (" ".join(desc.split()) if desc else "—"), key


def test_example3_seed_row(config_md: str) -> None:
    defaults = yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))
    rows = {r[0]: r for r in _rows(config_md)}
    assert rows["seed"][1] == str(defaults["seed"])
    assert rows["seed"][2] == "integer"


def test_configuration_grouped_by_top_level_section(config_md: str) -> None:
    defaults = yaml.safe_load(DEFAULTS.read_text(encoding="utf-8"))
    sections = _sections(config_md)
    for top, v in defaults.items():
        if isinstance(v, dict) and v:
            body = sections[f"`{top}`"]
            for key, _ in _leaves(v, top + "."):
                assert f"| `{key}` |" in body, key


# --- tables.md and flags.md ---------------------------------------------------------------------


@pytest.fixture(scope="module")
def tables_md(gr: ModuleType) -> str:
    text: str = gr.render_tables(SCHEMAS)
    return text


def _table_schemas() -> list[dict[str, Any]]:
    return [
        json.loads(p.read_text(encoding="utf-8"))
        for p in sorted((SCHEMAS / "tables").glob("*.schema.json"))
    ]


def test_every_resource_and_table_schema_documented(tables_md: str) -> None:
    heads = set(_sections(tables_md))
    package = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    for res in package["resources"]:
        assert f"`{res['name']}`" in heads, res["name"]
    for schema in _table_schemas():
        assert f"`{schema['name']}`" in heads, schema["name"]


@pytest.mark.parametrize("schema", _table_schemas(), ids=lambda s: s["name"])
def test_every_field_has_a_row(tables_md: str, schema: dict[str, Any]) -> None:
    body = _sections(tables_md)[f"`{schema['name']}`"]
    names = [r[0] for r in _rows(body)]
    assert names == [f["name"] for f in schema["fields"]]


def test_example4_genes(tables_md: str) -> None:
    body = _sections(tables_md)["`genes`"]
    assert "- Coordinates: 1-based closed" in body
    assert _rows(body)[0][0] == "locus_tag"


GENE_FLAGS_94 = [
    "partial",
    "crosses_origin",
    "pseudogene",
    "multi_identical",
    "shared_evidence",
    "two_proteins_one_gene",
    "cross_strain_protein",
    "imputed_abundance",
    "tolerant_evidence",
    "strain_difference",
    "nterm_extension_evidence",
    "cai_reference",
    "top3_lt3",
    "detectability_constant",
    "abundance_not_explained_by_sequence",
]


def _pattern_tokens(pattern: str) -> list[str]:
    m = re.match(r"\(\?:(.*?)\)\(\?:;", pattern)
    assert m, pattern
    return [t for t in m.group(1).split("|") if re.fullmatch(r"[a-z0-9_]+", t)]


def _flag_patterns() -> list[tuple[str, str]]:
    out = []
    for schema in _table_schemas():
        for f in schema["fields"]:
            pat = f.get("constraints", {}).get("pattern")
            if f["name"] == "flags" and pat:
                out.append((schema["name"], pat))
    return out


def test_flags_has_the_94_gene_vocabulary(gr: ModuleType) -> None:
    text = gr.render_flags(SCHEMAS)
    body = _sections(text)["`genes`"]
    for flag in GENE_FLAGS_94:
        assert f"`{flag}`" in body, flag


@pytest.mark.parametrize(("table", "pattern"), _flag_patterns(), ids=lambda x: str(x)[:30])
def test_flags_lists_every_schema_token(gr: ModuleType, table: str, pattern: str) -> None:
    body = _sections(gr.render_flags(SCHEMAS))[f"`{table}`"]
    for token in _pattern_tokens(pattern):
        assert f"`{token}`" in body, (table, token)


# --- indexes ------------------------------------------------------------------------------------


def _first_heading(p: Path) -> str:
    for line in p.read_text(encoding="utf-8").splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return p.stem


def test_methods_index_lists_every_note(gr: ModuleType) -> None:
    text = gr.render_index(DOCS / "methods", "methods")
    for p in sorted((DOCS / "methods").glob("*.md")):
        if p.name == "index.md":
            assert "(index.md)" not in text
            continue
        assert f"- [{_first_heading(p)}]({p.name})" in text, p.name


def test_adr_index_lists_every_adr_not_the_template(gr: ModuleType) -> None:
    text = gr.render_index(DOCS / "adr", "adr")
    assert "TEMPLATE" not in text
    adrs = sorted((DOCS / "adr").glob("ADR-*.md"))
    assert len(adrs) > 50
    for p in adrs:
        assert f"- [{_first_heading(p)}]({p.name})" in text, p.name


# --- the committed pages are current ------------------------------------------------------------


@pytest.mark.parametrize(
    "rel",
    [
        "reference/cli.md",
        "reference/configuration.md",
        "reference/tables.md",
        "reference/flags.md",
        "methods/index.md",
        "adr/index.md",
    ],
)
def test_generated_page_committed(rel: str) -> None:
    assert (DOCS / rel).is_file(), f"docs/{rel} is not committed (M11-07: pages live in git)"


def test_check_passes_on_this_checkout() -> None:
    assert SCRIPT.is_file(), "scripts/docs/gen_reference.py does not exist"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
