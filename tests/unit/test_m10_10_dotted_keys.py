"""``render`` resolves keys that hold dots (contig accessions such as ``CP171243.1``).

The pilot genome's plasmid subsets are ``plasmid:CP171243.1``; the dotted-path resolver split the
version suffix off and reported the share as missing (step 6 render, 2026-10-25).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def m() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "pilot_metrics_dotted", ROOT / "validation/pilot/metrics.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


DOC = {
    "m3": {"A1": {"subsets": {"plasmid:CP171243.1": {"share": 0.25}, "chromosome": {"share": 0.5}}}}
}


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("m3.A1.subsets.plasmid:CP171243.1.share", 0.25),
        ("m3.A1.subsets.chromosome.share", 0.5),
    ],
)
def test_dotted_key_resolves(m: ModuleType, path: str, value: float) -> None:
    assert m.resolve(DOC, path) == value


def test_longest_key_wins(m: ModuleType) -> None:
    doc = {"a": {"b": 1}, "a.b": 2}
    assert m.resolve(doc, "a.b") == 2


def test_missing_dotted_key_stays_missing(m: ModuleType) -> None:
    assert m.resolve(DOC, "m3.A1.subsets.plasmid:CP171243.2.share") is m._MISSING


def test_list_index_still_works(m: ModuleType) -> None:
    assert m.resolve({"x": [{"y": 3}]}, "x.0.y") == 3
