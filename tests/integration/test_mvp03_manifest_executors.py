"""MVP-03a examples 9 and 10: the manifest of a local and a Nextflow run of the same inputs
(docs/specs/MVP-03-manifest-finalize.md "Examples" 9-10, "Goal", "Test ideas" mutation guard;
GIPPA_SPEC.md §6.4 extended to the manifest, §6.7, §7.9.3).

The ``test_m2_21_local_vs_nextflow`` setup: gippasim SIM100K with M3-19 proteomics and every
precomputed layer, peptides only, ``--modules core,mobile,bgc --seed 7``, once with ``--executor
local`` and once with ``--executor nextflow --profile farmshare_local`` (real Nextflow and the
prebuilt GIPPA env under the real ``$GIPPA_ROOT``; no annotation tool runs). The database directory
is the test DB of ``helpers.mapping_steps`` (never the real one).

Mutation guard: ``stable_view`` with one of ``/executor``, ``/outdir``, ``/command_line`` or
``/steps`` dropped from ``VOLATILE`` no longer compares equal, so each of those entries is needed.
(``/inputs/*/path`` is not guarded here: both runs read the same input files.)
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.manifest_data import validate

pytestmark = [pytest.mark.needs_nextflow, pytest.mark.slow]

SEED = 7
GUARDED = ("/executor", "/outdir", "/command_line", "/steps")


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import gippasim
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    return gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim")


def _run(sim: Any, out: Path, *extra: str) -> None:
    from gippa.cli import cli

    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides)]
    args += ["--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml)]
    args += ["--outdir", str(out), "--modules", "core,mobile,bgc", "--seed", str(SEED), *extra]
    res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


@pytest.fixture(scope="module")
def outs(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    mp = pytest.MonkeyPatch()
    try:
        from helpers.mapping_steps import write_test_db

        mp.setenv("GIPPA_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db_mvp03"))))
        mp.setenv("NXF_ANSI_LOG", "false")
        local = tmp_path_factory.mktemp("LOCAL") / "OUT"
        nf = tmp_path_factory.mktemp("NF") / "OUT"
        _run(sim, local, "--executor", "local")
        _run(sim, nf, "--executor", "nextflow", "--profile", "farmshare_local")
    finally:
        mp.undo()
    return {"local": local, "nextflow": nf}


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


# ----------------------------------------------------------------------------- example 9


@pytest.mark.parametrize("executor", ["local", "nextflow"])
def test_example_9_manifest_validates(outs: dict[str, Path], executor: str) -> None:
    m = _manifest(outs[executor])
    validate(m)
    assert m["status"] == "completed" and m["exit_code"] == 0
    assert m["executor"] == executor
    assert m["modules_off"] == ["function"]
    assert m["tools"] == []  # every layer precomputed; gippa itself is left out


def test_example_9_inputs(outs: dict[str, Path]) -> None:
    out = outs["local"]
    m = _manifest(out)
    pf = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    by = {r["role"]: r for r in m["inputs"]}
    assert len(by) == len(m["inputs"]), "a role appears twice"
    genome = next(r for r in pf["inputs"] if r["role"] == "genome")
    assert by["genome"]["sha256"] == genome["sha256"]
    assert "precomputed:bakta_json" in by
    for role, row in by.items():
        if role.startswith("precomputed:"):
            p = Path(row["path"])
            assert row["sha256"] == hashlib.sha256(p.read_bytes()).hexdigest(), role


@pytest.mark.parametrize("executor", ["local", "nextflow"])
def test_example_9_digest_of_the_provenance_core(outs: dict[str, Path], executor: str) -> None:
    from gippa.manifest import canonical_json, provenance_core

    out = outs[executor]
    m = _manifest(out)
    prov = out / "work" / "provenance.json"
    core = json.loads(prov.read_text(encoding="utf-8"))
    assert m["provenance_digest"] == hashlib.sha256(canonical_json(core).encode()).hexdigest()
    assert core == provenance_core(m)


# ----------------------------------------------------------------------------- example 10


def test_example_10_stable_views_equal(outs: dict[str, Path]) -> None:
    from gippa.manifest import stable_view

    a, b = _manifest(outs["local"]), _manifest(outs["nextflow"])
    va, vb = stable_view(a), stable_view(b)
    differ = sorted(k for k in set(va) | set(vb) if va.get(k) != vb.get(k))
    assert not differ, f"local and Nextflow manifests differ on {differ}"


def test_example_10_same_digest(outs: dict[str, Path]) -> None:
    a, b = _manifest(outs["local"]), _manifest(outs["nextflow"])
    assert a["provenance_digest"] == b["provenance_digest"]


def test_example_10_nextflow_steps_have_rss(outs: dict[str, Path]) -> None:
    steps = _manifest(outs["nextflow"])["steps"]
    assert steps
    assert all("peak_rss_bytes" in s for s in steps if s["status"] == "completed"), steps


@pytest.mark.parametrize("pointer", GUARDED)
def test_mutation_guard_each_volatile_entry_is_needed(outs: dict[str, Path], pointer: str) -> None:
    from gippa.manifest import VOLATILE, stable_view

    a, b = _manifest(outs["local"]), _manifest(outs["nextflow"])
    fewer = [p for p in VOLATILE if p != pointer]
    assert stable_view(a, volatile=fewer) != stable_view(b, volatile=fewer), pointer
