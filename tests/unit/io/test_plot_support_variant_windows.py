"""The ``variant_windows`` support table's schema and datapackage entry, and the ``plots`` step
writing it (M9-18 "Support"; ADR-0043; the M9-03 support conventions).

A separate file from tests/unit/io/test_plot_support.py and tests/unit/steps/
test_plots_support_step.py (M9-03's tests are not edited here).

Interface assumed beyond the spec text (M9-18 test author):

- ``src/gippa/schemas/plot_support/variant_windows.json`` with fields ``variant_id, gene, aa_pos,
  window_start_aa, ref_window`` (``aa_pos`` / ``window_start_aa`` integers ≥ 1), listed in
  ``schemas/datapackage.json`` as resource ``variant_windows`` with ``role: plot_support``.
- The ``plots`` step writes ``plots/support/variant_windows.tsv`` when it is given both
  ``variants`` (``tables/variants.tsv``) and ``genes`` (C2 ``genes.parquet``); without either,
  no ``variant_windows`` file. ``gippa report`` never writes it (it only reads it back).
"""

from __future__ import annotations

import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers import viz_variants as vv
from helpers.report_fixtures import cfg, run_out, swap_registry

import gippa.schemas
from gippa.steps import registry as reg

SCHEMAS = Path(gippa.schemas.__file__).resolve().parent
NAME = "variant_windows"
COLUMNS = ["variant_id", "gene", "aa_pos", "window_start_aa", "ref_window"]


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _schema() -> dict[str, Any]:
    data: dict[str, Any] = json.loads(
        (SCHEMAS / "plot_support" / f"{NAME}.json").read_text(encoding="utf-8")
    )
    return data


# ---------------------------------------------------------------- schema


def test_schema_fields_follow_the_spec_columns() -> None:
    schema = _schema()
    assert [f["name"] for f in schema["fields"]] == COLUMNS
    assert "NA" in schema.get("missingValues", [])
    fields = {f["name"]: f for f in schema["fields"]}
    for col in ("aa_pos", "window_start_aa"):
        assert fields[col]["type"] == "integer"
        assert fields[col].get("constraints", {}).get("minimum") == 1
    assert fields["ref_window"]["type"] == "string"
    assert not {"p_value", "q_value", "stat", "effect", "null_model"} & set(COLUMNS)


def test_datapackage_lists_variant_windows_with_role_plot_support() -> None:
    pkg = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    res = [r for r in pkg["resources"] if r["name"] == NAME]
    assert len(res) == 1, [r["name"] for r in pkg["resources"]]
    r = res[0]
    assert r["role"] == "plot_support"
    assert r["path"].startswith(f"plots/support/{NAME}.tsv")
    assert r.get("profile") == "tabular-data-resource"
    assert (SCHEMAS / r["schema"]).resolve() == (SCHEMAS / "plot_support" / f"{NAME}.json")


def test_written_table_validates_against_its_schema(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    from gippa.io.plot_support import write_support
    from gippa.io.run_view import RunView
    from gippa.viz.support import variant_windows

    out = vv.build_out(tmp_path, vv.main_variants())
    variants = RunView.from_outdir(out).table("variants")
    path = write_support(NAME, variant_windows(variants, vv.genes_parquet_frame()), tmp_path / "p")
    resource = frictionless.Resource(
        path=path.name,
        basepath=str(path.parent),
        format="tsv",
        schema=frictionless.Schema.from_descriptor(_schema()),
    )
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ---------------------------------------------------------------- the plots step


def _inputs(tmp: Path) -> dict[str, Path]:
    out = vv.build_out(tmp / "src", vv.main_variants())
    genes = tmp / "inputs" / "genes.parquet"
    genes.parent.mkdir(parents=True, exist_ok=True)
    vv.genes_parquet_frame().to_parquet(genes, index=False)
    return {"variants": out / "tables" / "variants.tsv", "genes": genes}


def _run(tmp: Path, extra: dict[str, Path]) -> Path:
    out = run_out(tmp)
    plots = tmp / "work" / "plots"
    inputs = {"preflight": out / "preflight" / "preflight.json", **extra}
    reg.run_step("plots", cfg(), inputs, {"plots": plots})
    return plots


def _tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


def _support(plots: Path) -> list[str]:
    d = plots / "support"
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def test_step_writes_variant_windows_with_variants_and_genes(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    plots = _run(tmp_path, _inputs(tmp_path))
    assert "variant_windows.tsv" in _support(plots)
    rows = _tsv(plots / "support" / "variant_windows.tsv")
    assert rows[0] == COLUMNS
    got = {r[0]: r for r in rows[1:]}
    assert set(got) == {"VAR0001", "VAR0002", "VAR0003", "VAR0004"}
    assert got["VAR0001"] == ["VAR0001", "GSIM_00105", "214", "204", vv.P105[203:224]]
    assert got["VAR0002"] == ["VAR0002", "GSIM_00200", "50", "40", vv.P200[39:60]]


@pytest.mark.parametrize("drop", ["variants", "genes"])
def test_step_without_variants_or_genes_writes_no_variant_windows(
    fake_registry: ModuleType, tmp_path: Path, drop: str
) -> None:
    ins = _inputs(tmp_path)
    del ins[drop]
    plots = _run(tmp_path, ins)
    assert not [n for n in _support(plots) if n.startswith(NAME)]


def test_step_variant_windows_is_byte_identical_across_runs(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    ins = _inputs(tmp_path)
    a = _run(tmp_path / "a", ins)
    b = _run(tmp_path / "b", ins)
    name = "variant_windows.tsv"
    assert (a / "support" / name).read_bytes() == (b / "support" / name).read_bytes()


def test_plots_rendered_by_the_step_see_variant_windows(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.viz import base

    seen: dict[str, Any] = {}

    @base.plot("V90", "support_probe")
    def _probe(view: Any, c: Any) -> Any:
        seen["windows"] = view.support(NAME)
        return base.Skip("probe")

    _run(tmp_path, _inputs(tmp_path))
    assert seen["windows"] is not None
    assert set(seen["windows"]["variant_id"]) == {"VAR0001", "VAR0002", "VAR0003", "VAR0004"}
