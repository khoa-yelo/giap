"""M11-15b ``gippa demo`` end to end on the local executor (spec M11-15 item 2, example 2.1, AC 3).

One ``gippa demo --outdir`` run with no ``$GIPPA_DB_DIR``, an empty ``PATH`` (no Nextflow,
Docker or Apptainer) and sockets that refuse to connect; then a second, independent
``gippa run --executor local`` on the same bundled genome, peptides and layers (gunzipped into
``tmp``, an ``annotation.precomputed`` config written by hand), the §6.4 check pattern of
tests/integration/test_m2_21_local_vs_nextflow.py. Every published table must be byte-identical
(no committed sha256 listing, so table changes and lock bumps do not break this test).

Layout pinned as in tests/unit/test_m11_15b_demo.py: ``gippa/data/demo/layers/bakta.json.gz``
and ``layers/<key>.tsv.gz``.
"""

from __future__ import annotations

import gzip
import json
import re
import socket
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow]

REPO = Path(__file__).resolve().parents[2]
DEMO = REPO / "src" / "gippa" / "data" / "demo"
LAYERS = DEMO / "layers"
BANNER = "Synthetic demo data (SIM100K, gippasim seed 20261014)"
DONE = re.compile(r"^gippa demo: done in (\d+) s; report: (\S.*)$", re.M)
MODULES = ["core", "mobile", "bgc"]


def _layer_keys() -> list[str]:
    from gippa.config.model import PrecomputedConfig

    return [k for k in PrecomputedConfig.model_fields if k != "bakta_json"]


def _plain(path: Path) -> bytes:
    data = path.read_bytes()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


@dataclass(frozen=True)
class Runs:
    demo_out: Path
    demo_output: str
    local_out: Path
    connects: list[Any]


def _independent_config(tmp: Path) -> tuple[Path, Path]:
    """Gunzipped bundled layers and a config naming each; returns (config, bakta json)."""
    lay = tmp / "layers"
    lay.mkdir(parents=True)
    bakta = lay / "bakta.json"
    bakta.write_bytes(_plain(LAYERS / "bakta.json.gz"))
    lines = ["annotation:", "  precomputed:", f"    bakta_json: {json.dumps(str(bakta))}"]
    for key in _layer_keys():
        p = lay / f"{key}.tsv"
        p.write_bytes(_plain(LAYERS / f"{key}.tsv.gz"))
        lines.append(f"    {key}: {json.dumps(str(p))}")
    cfg = tmp / "config.yaml"
    cfg.write_text("\n".join(lines) + "\n")
    return cfg, bakta


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> Runs:
    from gippa.cli import cli

    assert LAYERS.is_dir(), "gippa/data/demo/layers/ missing (M11-15 item 2)"
    root = tmp_path_factory.mktemp("m11_15b")
    empty_bin = root / "bin"
    empty_bin.mkdir()
    connects: list[Any] = []

    def no_network(_self: socket.socket, address: Any) -> None:
        connects.append(address)
        raise OSError("network disabled in this test")

    with pytest.MonkeyPatch.context() as mp:
        for name in ("GIPPA_DB_DIR", "GIPPA_CONTAINER_RUNTIME", "NUMBA_CACHE_DIR"):
            mp.delenv(name, raising=False)
        mp.setenv("PATH", str(empty_bin))
        mp.setattr(socket.socket, "connect", no_network)
        mp.chdir(root)
        demo_out = root / "demo"
        res = CliRunner().invoke(cli, ["demo", "--outdir", str(demo_out)])
        if res.exception is not None and not isinstance(res.exception, SystemExit):
            raise res.exception
        assert res.exit_code == 0, res.output
        cfg, bakta = _independent_config(root / "indep")
        local_out = root / "local"
        args = ["run", "--executor", "local", "--genome", str(DEMO / "sim100k.fna.gz")]
        args += ["--peptides", str(DEMO / "sim100k_peptides.tsv.gz"), "--annotation", str(bakta)]
        args += ["--config", str(cfg), "--modules", ",".join(MODULES), "--outdir", str(local_out)]
        ind = CliRunner().invoke(cli, args)
        if ind.exception is not None and not isinstance(ind.exception, SystemExit):
            raise ind.exception
        assert ind.exit_code == 0, ind.output
    return Runs(demo_out, res.output, local_out, connects)


def _manifest(out: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((out / "run_manifest.json").read_text())
    return data


def _tables(out: Path) -> dict[str, bytes]:
    files = [p for d in ("tables", "annotation") for p in sorted((out / d).rglob("*"))]
    return {str(p.relative_to(out)): p.read_bytes() for p in files if p.is_file()}


def test_example_2_1_report_written_and_path_printed(runs: Runs) -> None:
    m = DONE.search(runs.demo_output)
    assert m is not None, runs.demo_output
    report = Path(m.group(2).strip())
    assert report == (runs.demo_out / "report" / "report.html").resolve()
    assert report.is_file()


def test_example_2_1_within_two_minutes(runs: Runs) -> None:
    m = DONE.search(runs.demo_output)
    assert m is not None, runs.demo_output
    assert int(m.group(1)) <= 120


def test_no_network_used(runs: Runs) -> None:
    assert runs.connects == []


def test_manifest_records_demo_and_completion(runs: Runs) -> None:
    man = _manifest(runs.demo_out)
    assert man["demo"] is True
    assert man["status"] == "completed"
    assert man["executor"] == "local"
    assert sorted(man["config"]["modules"]) == sorted(MODULES)
    assert "function" in man["modules_off"]


def test_independent_run_is_not_marked_demo(runs: Runs) -> None:
    assert _manifest(runs.local_out).get("demo") in (None, False)


def test_report_shows_the_synthetic_data_banner(runs: Runs) -> None:
    html = (runs.demo_out / "report" / "report.html").read_text(encoding="utf-8")
    assert BANNER in html


def test_independent_report_has_no_banner(runs: Runs) -> None:
    html = (runs.local_out / "report" / "report.html").read_text(encoding="utf-8")
    assert BANNER not in html


def test_ac3_tables_byte_identical_to_independent_local_run(runs: Runs) -> None:
    demo, local = _tables(runs.demo_out), _tables(runs.local_out)
    assert "tables/genes.tsv" in demo
    assert sorted(demo) == sorted(local)
    differ = [k for k in demo if demo[k] != local[k]]
    assert differ == []


def test_every_bundled_layer_used_including_rbs_novel(runs: Runs) -> None:
    man = _manifest(runs.demo_out)
    pre = man["config"]["annotation"]["precomputed"]
    for key in ["bakta_json", *_layer_keys()]:
        assert pre.get(key), f"annotation.precomputed.{key} not set in the demo run"
    keys = set(_layer_keys())
    assert not keys & set(man["layers_missing"]), man["layers_missing"]
    reasons = {s["name"]: s["reason"] for s in man["skipped"]}
    assert reasons.get("ostir_novel") == "precomputed:rbs_novel", reasons


def test_no_external_tool_step_ran(runs: Runs) -> None:
    """Every external tool was replaced by its layer (no tool, no container on this host)."""
    man = _manifest(runs.demo_out)
    reasons = {s["name"]: s["reason"] for s in man["skipped"]}
    for tool in ("promotech", "ostir", "transtermhp", "uniop", "promoter_calculator"):
        assert reasons.get(tool, "").startswith("precomputed:"), (tool, reasons.get(tool))


def test_demo_adds_no_warning_of_its_own(runs: Runs) -> None:
    """Spec edge case "the demo under --strict: a warning is a bug". SIM100K itself trips
    genome.contiguity (108 kb < 1 Mb), so the claim checked is: no warning beyond those an
    independent local run of the same inputs raises."""
    demo = set(_manifest(runs.demo_out)["warnings"])
    local = set(_manifest(runs.local_out)["warnings"])
    assert demo <= local, demo - local


def test_demo_preflight_has_no_error(runs: Runs) -> None:
    pf = json.loads((runs.demo_out / "preflight" / "preflight.json").read_text())
    assert [c["rule_id"] for c in pf["checks"] if c["level"] == "error"] == []
