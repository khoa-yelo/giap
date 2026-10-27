"""Fake plots, small ``OUT/`` trees and an HTML tree reader for the M9-02 tests (the ``plots`` and
``report`` steps, ``gippa.report.build`` and ``gippa report``).

Fake plots register into a swapped ``gippa.viz.base._REGISTRY`` (``fake_registry``), as
tests/unit/viz/test_registry.py does. Every ``gippa.steps`` and ``gippa.viz`` / ``gippa.report``
module is imported first, so a real plot module imported for the first time while the registry is
swapped can never register into the temporary dict (the ``registry`` fixture of
tests/unit/conftest.py does the same for steps).

The HTML reader is a small tree builder over ``html.parser`` (standard library); it keeps SVG
elements, so ids and ``url(#...)`` references inside embedded SVGs can be checked.

Interface pinned by these helpers (M9-02 test author, beyond the spec text):
- every report section is ``<section id="<SECTIONS id>">``; nothing else is a ``<section>``;
- the slot of each plot ID that a section lists carries ``data-plot="Vnn"`` (rendered or not);
- each summary card is an element with ``class="card"`` and ``data-card="<name>"`` holding an
  element of class ``value`` (the number) and one of class ``estimate`` (the error estimate or the
  word ``unestimated``); names: ``genes_detected``, ``peptides_placed``, ``islands``,
  ``operon_flags``, ``novel_orfs``;
- each table shown in the report is ``<table data-table="<PATHS key>">``.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import pkgutil
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

#: the 13 §10.3 section ids, in order (M9-02 behaviour 3)
SECTION_IDS = (
    "summary",
    "inputs-qc",
    "mapping",
    "detection",
    "genome-context",
    "operons",
    "variants",
    "novel-orfs",
    "function",
    "determinants",
    "sensitivity",
    "methods",
    "provenance",
)
#: section id -> the plot IDs it shows, in order (behaviour 3)
SECTION_PLOTS: Mapping[str, tuple[str, ...]] = {
    "summary": (),
    "inputs-qc": ("V17",),
    "mapping": ("V05", "V06"),
    "detection": ("V07", "V18"),
    "genome-context": ("V01", "V02", "V10", "V03"),
    "operons": ("V11",),
    "variants": ("V12",),
    "novel-orfs": ("V13",),
    "function": ("V08", "V14"),
    "determinants": ("V15",),
    "sensitivity": (),
    "methods": (),
    "provenance": (),
}
#: plot ID -> name used by the fake plots (the §10.2 names, so file stems look real)
NAMES: Mapping[str, str] = {
    "V01": "circos",
    "V02": "abundance_landscape",
    "V03": "locus",
    "V05": "mapping_sankey",
    "V06": "mapping_fdr",
    "V07": "obs_vs_exp",
    "V08": "cog_abundance",
    "V10": "correlogram",
    "V11": "operon_pairs",
    "V12": "variants",
    "V13": "novel_orfs",
    "V14": "enrichment",
    "V15": "determinants",
    "V17": "qc",
    "V18": "stat_diagnostics",
}

INPUTS = [
    {"file": "genome.fasta", "role": "genome", "sha256": "a" * 64, "size": 1000},
    {"file": "peptides.tsv", "role": "peptides", "sha256": "b" * 64, "size": 2000},
]


# ---------------------------------------------------------------- registry


def _import_tree(package: str) -> None:
    try:
        pkg = importlib.import_module(package)
    except ModuleNotFoundError:
        return
    paths = getattr(pkg, "__path__", None)
    if paths is None:
        return
    for info in pkgutil.walk_packages(paths, prefix=f"{package}."):
        importlib.import_module(info.name)


def preload() -> None:
    """Import every step, viz and report module, so none registers while a registry is swapped."""
    from gippa.steps import load_step_modules

    load_step_modules()
    for package in ("gippa.viz", "gippa.report"):
        _import_tree(package)


def swap_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Give ``gippa.viz.base`` an empty plot registry for one test (each test module wraps this in
    its own ``fake_registry`` fixture)."""
    preload()
    from gippa.viz import base

    monkeypatch.setattr(base, "_REGISTRY", {})
    return base


def figure(kind: str = "line") -> Any:
    """A small matplotlib figure; two calls give byte-identical SVGs with the same internal ids."""
    from gippa.viz import style

    fig = style.new_figure("single", 1.5)
    ax = fig.add_subplot()
    if kind == "line":
        ax.plot([0, 1, 2], [0.0, 1.0, 0.5])
    else:
        ax.bar([0, 1], [2.0, 3.0])
    ax.set_title("fake plot: 2 of 3 detected")
    return fig


def frame() -> pd.DataFrame:
    return pd.DataFrame({"x": [0, 1, 2], "y": [0.0, 1.0, 0.5]})


def register_fake(
    pid: str,
    *,
    name: str | None = None,
    needs: Sequence[str] = (),
    module: str = "core",
    suffixes: Sequence[str | None] = (None,),
    raises: bool = False,
    calls: list[str] | None = None,
) -> None:
    """Register a fake plot ``pid`` drawing one instance per suffix (identical figures)."""
    from gippa.viz import base

    @base.plot(pid, name or NAMES.get(pid, "fake"), needs=list(needs), module=module)
    def _fake(view: Any, cfg: Any) -> Any:
        if calls is not None:
            calls.append(pid)
        if raises:
            raise RuntimeError(f"fake plot {pid} failed")
        return [
            base.Instance(suffix=s, figure=figure(), data=frame(), title="fake") for s in suffixes
        ]


# ---------------------------------------------------------------- OUT trees


def write_preflight(out: Path, inputs: Sequence[Mapping[str, Any]] = INPUTS) -> Path:
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    body = {
        "schema_version": "1",
        "gippa_version": "0.0.0",
        "inputs": [dict(i) for i in inputs],
        "replicons": [{"id": "chr1", "length": 10000, "topology": "circular"}],
    }
    path.write_text(json.dumps(body, sort_keys=True, indent=1), encoding="utf-8")
    return path


def island_rows(n: int, p_values: Sequence[float] | None = None) -> list[dict[str, object]]:
    from helpers.viz_out_tree import ISLAND_ROW

    rows = []
    for i in range(n):
        row = dict(ISLAND_ROW)
        row["island_id"] = f"ISL{i + 1:04d}"
        row["start"] = 1001 + 10 * i
        row["end"] = 5400 + 10 * i
        if p_values is not None:
            row["p_value"] = p_values[i]
        rows.append(row)
    return rows


def run_out(tmp: Path, *, islands: int = 0, preflight: bool = True) -> Path:
    """``OUT/`` with ``tables/``, a preflight JSON and, optionally, ``islands`` island rows."""
    from helpers.viz_out_tree import make_out, write_schema_table

    out = make_out(tmp)
    if preflight:
        write_preflight(out)
    if islands:
        write_schema_table(out / "tables" / "islands.tsv", "islands", island_rows(islands))
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cfg(*overrides: str) -> Any:
    """Defaults with small, fast plot output (SVG only, 50 dpi) unless overridden."""
    from gippa.config import resolve_config

    return resolve_config(overrides=["report.formats=[svg]", "report.dpi=50", *overrides])


# ---------------------------------------------------------------- HTML


VOID = frozenset(
    {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "source", "wbr"}
)


@dataclass(eq=False)
class Node:
    tag: str
    attrs: dict[str, str]
    parent: Node | None = None
    children: list[Node | str] = field(default_factory=list)

    def text(self) -> str:
        parts: list[str] = []
        for c in self.children:
            parts.append(c if isinstance(c, str) else c.text())
        return "".join(parts)

    def iter(self) -> Iterator[Node]:
        yield self
        for c in self.children:
            if isinstance(c, Node):
                yield from c.iter()

    def find_all(self, pred: Callable[[Node], bool]) -> list[Node]:
        return [n for n in self.iter() if pred(n)]

    def by_tag(self, tag: str) -> list[Node]:
        return self.find_all(lambda n: n.tag == tag)

    def by_attr(self, attr: str, value: str | None = None) -> list[Node]:
        return self.find_all(
            lambda n: attr in n.attrs and (value is None or n.attrs[attr] == value)
        )

    def has_class(self, cls: str) -> bool:
        return cls in self.attrs.get("class", "").split()

    def by_class(self, cls: str) -> list[Node]:
        return self.find_all(lambda n: n.has_class(cls))


class _Builder(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.root = Node("#root", {})
        self.cur = self.root
        self.decls: list[str] = []
        self.pis: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        node = Node(tag, {k: v or "" for k, v in attrs}, self.cur)
        self.cur.children.append(node)
        if tag not in VOID:
            self.cur = node

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.cur.children.append(Node(tag, {k: v or "" for k, v in attrs}, self.cur))

    def handle_endtag(self, tag: str) -> None:
        node: Node | None = self.cur
        while node is not None and node.tag != tag:
            node = node.parent
        if node is not None and node.parent is not None:
            self.cur = node.parent

    def handle_data(self, data: str) -> None:
        self.cur.children.append(data)

    def handle_decl(self, decl: str) -> None:
        self.decls.append(decl)

    def handle_pi(self, data: str) -> None:
        self.pis.append(data)


@dataclass
class Page:
    html: str
    root: Node
    pis: list[str]

    def sections(self) -> list[Node]:
        return self.root.by_tag("section")

    def section(self, sid: str) -> Node:
        found = [s for s in self.sections() if s.attrs.get("id") == sid]
        assert len(found) == 1, f"section {sid!r}: {len(found)} found"
        return found[0]

    def slot(self, pid: str) -> Node:
        found = self.root.by_attr("data-plot", pid)
        assert found, f"no element with data-plot={pid!r}"
        return found[0]

    def ids(self) -> list[str]:
        return [n.attrs["id"] for n in self.root.iter() if "id" in n.attrs]

    def links(self) -> list[str]:
        """Every ``href``, ``src`` and ``xlink:href`` value."""
        out = []
        for n in self.root.iter():
            for attr in ("href", "src", "xlink:href"):
                if attr in n.attrs:
                    out.append(n.attrs[attr])
        return out

    def cards(self) -> dict[str, Node]:
        return {
            n.attrs["data-card"]: n for n in self.root.by_class("card") if "data-card" in n.attrs
        }

    def table(self, key: str) -> Node | None:
        found = self.root.find_all(lambda n: n.tag == "table" and n.attrs.get("data-table") == key)
        return found[0] if found else None


def parse(html: str) -> Page:
    b = _Builder()
    b.feed(html)
    b.close()
    return Page(html, b.root, b.pis)


URL_REF = re.compile(r"url\(\s*['\"]?#([^)'\"\s]+)['\"]?\s*\)")


def local_refs(page: Page) -> list[str]:
    """Every ``#id`` an element refers to: ``href="#x"``, ``xlink:href="#x"`` and ``url(#x)``."""
    refs: list[str] = []
    for n in page.root.iter():
        for attr, value in n.attrs.items():
            if attr in ("href", "xlink:href") and value.startswith("#") and len(value) > 1:
                refs.append(value[1:])
            refs.extend(URL_REF.findall(value))
    for n in page.root.by_tag("style"):
        refs.extend(URL_REF.findall(n.text()))
    return refs


def table_rows(table: Node) -> tuple[list[str], list[list[str]]]:
    """Header cells and body rows (cell texts, stripped) of a ``<table>``."""
    rows = table.by_tag("tr")
    header: list[str] = []
    body: list[list[str]] = []
    for tr in rows:
        ths = [c for c in tr.children if isinstance(c, Node) and c.tag == "th"]
        tds = [c for c in tr.children if isinstance(c, Node) and c.tag == "td"]
        if ths and not tds and not header:
            header = [c.text().strip() for c in ths]
        elif tds:
            body.append([c.text().strip() for c in tds])
    return header, body
