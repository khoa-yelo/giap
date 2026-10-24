"""``validation/pilot/te_real.py compute`` on the synthetic pilot (M10-06 acceptance criteria 1
(example 2's exceptions row), 8 (the ``compute`` preconditions) and 10 (``compute`` half)).

Spec: docs/specs/M10-06-te-real-and-te7.md ("Inputs" preconditions, "Outputs", "Interface",
example 2, edge case "Contaminant FASTA changed since the searches"). Inputs: M10-05's synthetic
pilot (``helpers.pilot_metrics_m10.build_pilot``: giapsim SIM100K, schema-valid manifests) plus
an A2 ``genes.tsv`` with r-proteins and the four TE4 genes, a ``variants.tsv`` and a ``Cont_``
contaminant FASTA (``helpers.te_real_m10.build_te_pilot``). Nothing from ``$GIAP_ROOT``.

Interface fixed here: a manifest pins the contaminant library through a ``databases`` entry named
``contaminants`` whose ``checksum`` is ``sha256:<hex>`` (the ``giap setup`` form); exceptions
kinds ``te2_violation`` and ``te3_not_detected`` with ``id`` = the locus tag.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_metrics_m10 import REPO, World, output_files, sha256, sim_world
from helpers.te_real_m10 import (
    build_te_pilot,
    codes,
    edit_manifest,
    run_te_real,
    te_compute_args,
)

PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_06_sim"))


def _compute(p: Any, contam: Path, out: Path) -> Any:
    return run_te_real(*te_compute_args(p, contam, out), cwd=p.root)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _doc(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((out / "te_real.json").read_text(encoding="utf-8"))
    return doc


def _exceptions(out: Path) -> list[dict[str, str]]:
    with (out / "te_real_exceptions.tsv").open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


def test_compute_layout_and_all_pass(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    raw = (out / "te_real.json").read_text(encoding="utf-8")
    doc = json.loads(raw)
    assert raw == json.dumps(doc, sort_keys=True, indent=2) + "\n"
    assert doc["schema"] == "giap-pilot-te-real/1"
    assert re.fullmatch(r"[0-9a-f]{40}", str(doc["spec"]))
    assert set(doc) >= {"schema", "inputs", "spec", "te2", "te3", "te4", "te5", "warnings"}
    for name, rec in doc["inputs"].items():
        assert set(rec) >= {"path", "sha256"}, name
    assert sha256(contam) in {rec["sha256"] for rec in doc["inputs"].values()}
    assert doc["te2"]["status"] == "pass"
    assert doc["te3"]["status"] == "pass"
    assert doc["te4"]["status"] == "pass"
    assert set(doc["te5"]) >= {"A1", "A2", "status"}
    assert doc["te5"]["A1"]["status"] == doc["te5"]["A2"]["status"] == "pass"
    assert doc["te5"]["status"] == "pass"
    assert doc["te5"]["A2"]["n_s"] == 1  # the one Cont_P00761 (trypsin) peptide of the pilot
    # the manifests record no contaminants database: unpinned, warned, not an error
    assert "contaminants_unpinned" in codes(doc["warnings"])
    assert {r["kind"] for r in _exceptions(out)} <= {"te2_violation", "te3_not_detected"}


def test_example_2_through_the_cli(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world, violation=True, ribo_missed=2)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))  # a fail is a finding: exit 0
    doc = _doc(out)
    assert doc["te2"]["status"] == "fail" and doc["te2"]["counts"]["violations"] == 1
    assert doc["te2"]["breakdown"]["flags"] == {"database_derived;mass_equivalent": 1}
    assert doc["te3"]["status"] == "fail"  # 8 / 10 detected
    rows = _exceptions(out)
    assert list(rows[0]) == ["te", "kind", "id", "detail"]
    keys = [(r["te"], r["kind"], r["id"]) for r in rows]
    assert keys == sorted(keys)
    assert sum(r["kind"] == "te2_violation" for r in rows) == 1
    assert sum(r["kind"] == "te3_not_detected" for r in rows) == 2


def test_te_real_json_holds_no_sequence(tmp_path: Path, world: World) -> None:
    from helpers.pilot_metrics_m10 import CONT_PEPTIDE

    p, contam = build_te_pilot(tmp_path / "pilot", world, violation=True)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    raw = (out / "te_real.json").read_text(encoding="utf-8")
    for _, pep in p.peptides:
        assert pep not in raw
    assert CONT_PEPTIDE not in raw
    assert not PEPTIDE_LIKE.search(raw), PEPTIDE_LIKE.search(raw)


def test_missing_variants_is_not_evaluable(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world, variants=False)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    te2 = _doc(out)["te2"]
    assert te2["status"] == "not_evaluable"
    assert "missing:tables/variants.tsv" in te2["reason"]


def test_id_mixups_cross_check_warning(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world, violation=True, id_mixups=0.0)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    doc = _doc(out)
    every = codes(doc["warnings"]) + codes(doc["te2"].get("warnings", []))
    assert "te2_id_mixups_inconsistent" in every


def test_te5_fails_when_one_arm_fails(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world, a1_cont_unmapped=True)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    te5 = _doc(out)["te5"]
    assert te5["A2"]["status"] == "pass"
    assert te5["A1"]["status"] == "fail"
    assert te5["A1"]["counts"]["labelled_unmapped"] == 1
    assert te5["status"] == "fail"


def test_compute_twice_is_byte_identical(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world, violation=True, ribo_missed=1)
    o1, o2 = tmp_path / "m1", tmp_path / "m2"
    _ok(_compute(p, contam, o1))
    _ok(_compute(p, contam, o2))
    f1, f2 = output_files(o1), output_files(o2)
    assert {"te_real.json", "te_real_exceptions.tsv"} <= set(f1)
    assert f1 == f2


# ------------------------------------------------------------------ preconditions (exit 2)


def _pin(checksum: str) -> Any:
    def fn(m: dict[str, Any]) -> None:
        m["databases"] = [
            {
                "name": "contaminants",
                "version": "bddf6da9",
                "checksum": checksum,
                "retrieved": "2026-10-01",
                "size_bytes": 1000,
            }
        ]

    return fn


def test_contaminants_pinned_and_equal_no_warning(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    for out_dir in (p.out_a1, p.out_a2):
        edit_manifest(out_dir, _pin("sha256:" + sha256(contam)))
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    assert "contaminants_unpinned" not in codes(_doc(out)["warnings"])


def test_contaminant_fasta_changed_exit_2(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    for out_dir in (p.out_a1, p.out_a2):
        edit_manifest(out_dir, _pin("sha256:" + "0" * 64))
    out = tmp_path / "metrics"
    res = _compute(p, contam, out)
    assert res.returncode == 2, res.stdout + res.stderr
    assert _nothing_written(out)


def test_unknown_checksum_algorithm_exit_2(tmp_path: Path, world: World) -> None:
    """N7: a contaminants checksum with an algorithm GIAP does not know is a precondition (exit
    2), not an internal inconsistency."""
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    for out_dir in (p.out_a1, p.out_a2):
        edit_manifest(out_dir, _pin("crc99:" + sha256(contam)))
    out = tmp_path / "metrics"
    res = _compute(p, contam, out)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)


@pytest.mark.parametrize("name", ["peptide_hits.tsv", "unmapped.tsv"])
def test_te5_arm_without_its_tables_is_not_evaluable(
    tmp_path: Path, world: World, name: str
) -> None:
    """Test gap 4: an arm whose run has no peptide_hits.tsv or unmapped.tsv cannot be judged;
    the other arm still is, and TE5 is ``not_evaluable`` (no fail anywhere)."""
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    for rel in (f"tables/{name}", f"mapping/{name}"):
        (p.out_a1 / rel).unlink(missing_ok=True)
    out = tmp_path / "metrics"
    _ok(_compute(p, contam, out))
    te5 = _doc(out)["te5"]
    assert te5["A1"]["status"] == "not_evaluable"
    assert any(str(r).startswith("missing") for r in te5["A1"]["reason"]), te5["A1"]
    assert te5["A2"]["status"] == "pass"
    assert te5["status"] == "not_evaluable"


def test_failed_run_exit_2(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    edit_manifest(p.out_a1, lambda m: m.update(exit_code=1, status="failed"))
    out = tmp_path / "metrics"
    res = _compute(p, contam, out)
    assert res.returncode == 2, res.stdout + res.stderr
    assert _nothing_written(out)


def test_peptides_sha_mismatch_exit_2(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    args = te_compute_args(p, contam, tmp_path / "metrics")
    i = args.index("--peptides-a2")
    args[i + 1] = str(p.pep_a1)  # the A1 table is not what the A2 run read
    res = run_te_real(*args, cwd=p.root)
    assert res.returncode == 2, res.stdout + res.stderr
    assert _nothing_written(tmp_path / "metrics")


def test_out_inside_repository_exit_2(tmp_path: Path, world: World) -> None:
    p, contam = build_te_pilot(tmp_path / "pilot", world)
    out = REPO / ".m10_06_compute_out_test"
    try:
        res = _compute(p, contam, out)
        assert res.returncode == 2, res.stdout + res.stderr
        assert _nothing_written(out)
    finally:
        shutil.rmtree(out, ignore_errors=True)
