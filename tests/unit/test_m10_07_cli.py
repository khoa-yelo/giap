"""``validation/pilot/entrapment.py compute`` and the ``metrics.py compute --entrapment`` hook
(M10-07 acceptance criteria 5, 7, 8 and 9; example 10; cross-checks 1-3; the S-amended guard of
criterion 10: no discovered peptide string in ``entrapment.json``).

Spec: docs/specs/M10-07-free-entrapment.md ("Inputs", "Outputs", "Hook into M10-05",
"Cross-checks", "Acceptance criteria") with the binding "Statistician decisions" ("Other required
changes": caveat strings, data rule; S1 scoping; S3 ``ci``). Inputs are the spec's hand-made
example-2 fixture written in ``tmp`` (``tests/helpers/entrapment_m10_07.py``); the M10-05 pilot
for the hook is M10-05's synthetic builder (``tests/helpers/pilot_metrics_m10.py``). This branch
is stacked on ``m10/spec-05`` so ``metrics.py`` exists; only the hook is new.

Interface fixed here (the spec names the inputs, not the flags): ``compute --pep-a1 --k12-fasta
--contaminants --genome --strain-db --runs-manifest --out`` (required) and ``--out-a1``,
``--diann-cmd``, ``--provenance`` (optional), mirroring the spec's input names. ``--out`` is
METRICS_DIR. ``provenance.json`` lists the searched FASTAs as ``inputs.fasta = [{path, sha256}]``
(M10-02 "inputs: each --fasta (path, sha256)"). ``command_step1.txt`` is one argv token per line.
``metrics.py compute`` takes ``--entrapment PATH`` and writes the blocks under ``m8``.
"""

from __future__ import annotations

import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from helpers.entrapment_m10_07 import (
    CAVEATS,
    CONTAMINANTS,
    EXAMPLE2_ROWS,
    EXAMPLE8_ROW,
    K12,
    NISSLE,
    REPO,
    accepted_strings,
    codes,
    coding_text,
    diann_step1_tokens,
    read_json,
    row,
    run_entrapment,
    sha256,
    write_a1,
    write_faa,
    write_genome,
    write_out_a1,
    write_runs_manifest,
)

CFG_A1 = REPO / "validation" / "pilot" / "diann_A1.cfg"


@dataclass
class Inputs:
    root: Path
    pep_a1: Path
    k12: Path
    cont: Path
    genome: Path
    strain_db: Path
    manifest: Path

    def args(self, out: Path, *extra: str, **over: Path) -> list[str]:
        v = {"pep-a1": self.pep_a1, "k12-fasta": self.k12, "contaminants": self.cont,
             "genome": self.genome, "strain-db": self.strain_db,
             "runs-manifest": self.manifest}  # fmt: skip
        v |= {k.replace("_", "-"): p for k, p in over.items()}
        flat = [x for k, p in v.items() for x in (f"--{k}", str(p))]
        return ["compute", *flat, "--out", str(out), *extra]


def build_inputs(root: Path, rows: Any = EXAMPLE2_ROWS, *, tags: bool = True) -> Inputs:
    return Inputs(
        root=root,
        pep_a1=write_a1(root / "inputs" / "a1.tsv.gz", list(rows)),
        k12=write_faa(root / "data" / "UP000000625.fasta", K12),
        cont=write_faa(root / "db" / "contaminants.fasta", CONTAMINANTS, prefix="Cont_P"),
        genome=write_genome(
            root / "nissle.fna", {"nissle": coding_text(*NISSLE)}, ["nissle"] if tags else []
        ),
        strain_db=write_faa(root / "strain_db.faa", ["MSAGLLKVEPTLDEKAAGGYWR"], prefix="NISSLE_"),
        manifest=write_runs_manifest(root / "pxd078459_files.tsv"),
    )


@pytest.fixture
def inp(tmp_path: Path) -> Inputs:
    return build_inputs(tmp_path / "in")


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


def _exit2_nothing(res: Any, out: Path, needle: str | None = None) -> None:
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out), sorted(p.name for p in out.iterdir())
    if needle is not None:
        assert needle in res.stderr + res.stdout


@pytest.fixture
def computed(inp: Inputs, tmp_path: Path) -> tuple[Inputs, Path, dict[str, Any]]:
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out)))
    return inp, out, read_json(out / "entrapment.json")


# ------------------------------------------------------------------ outputs, layout, values


def test_compute_writes_the_two_files(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    _, out, _ = computed
    assert {p.name for p in out.iterdir()} >= {"entrapment.json", "entrapment_discoveries.tsv"}


def test_entrapment_json_layout(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    inp, out, doc = computed
    raw = (out / "entrapment.json").read_text(encoding="utf-8")
    assert raw.endswith("\n")
    assert raw == json.dumps(doc, sort_keys=True, indent=2) + "\n"
    assert doc["schema"] == "giap-pilot-entrapment/1"
    assert re.fullmatch(r"[0-9a-f]{40}", str(doc["spec"])), doc["spec"]
    for key in ("inputs", "settings", "universe", "peptide", "precursor", "per_run",
                "cross_checks", "caveats", "warnings"):  # fmt: skip
        assert key in doc, key
    for name, rec in doc["inputs"].items():
        assert set(rec) >= {"path", "sha256"}, name
    shas = {rec["sha256"] for rec in doc["inputs"].values()}
    assert {sha256(inp.pep_a1), sha256(inp.k12), sha256(inp.cont), sha256(inp.genome)} <= shas
    assert "timestamp" not in raw.lower() and "created" not in doc


def test_settings_block(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    s = computed[2]["settings"]
    assert (s["enzyme"], s["missed"], s["min_len"], s["max_len"], s["var_mods"]) == (
        "trypsin_p", 1, 7, 30, 1)  # fmt: skip
    assert list(s["charges"]) == [1, 4]
    assert list(s["mz"]) == [300.0, 1800.0]
    assert "tolerant" in s


def test_example_2_through_the_cli(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    doc = computed[2]
    u, p, q = doc["universe"], doc["peptide"], doc["precursor"]
    assert (u["shared"], u["k12_only"]) == (4, 2)
    assert u["r"] == pytest.approx(0.5, abs=1e-12)
    assert u["r_prec"] == pytest.approx(1 / 3, abs=1e-12)
    assert (p["N_T"], p["N_E"]) == (3, 1)
    assert p["fdp_combined"] == pytest.approx(0.75, abs=1e-12)
    assert p["assessment"] == "not_assessed"
    assert q["assessment"] == "evidence_of_failure"
    assert q["fdp_combined"] == pytest.approx(0.8, abs=1e-12)
    assert all(b["assessment"] == "not_assessed" for b in doc["per_run"].values())


def test_caveats_are_the_four_fixed_strings(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    assert list(computed[2]["caveats"]) == list(CAVEATS)


def test_discoveries_tsv(computed: tuple[Inputs, Path, dict[str, Any]]) -> None:
    _, out, _ = computed
    lines = (out / "entrapment_discoveries.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[:3] == ["peptide", "bucket", "present_rule"]
    got = {}
    for line in lines[1:]:
        pep, bucket, rule = [*line.split("\t"), ""][:3]
        got[pep.replace("I", "L")] = (bucket, rule)
    assert set(got) == {"MSAGLLK", "VEPTLDEK", "AGNVLTEPWR", "AAGGWWR", "QQWWLLDR"}
    assert got["MSAGLLK"] == ("target", "exact")
    assert got["VEPTLDEK"] == ("target", "exact")
    assert got["AGNVLTEPWR"] == ("target", "mass_equivalent")
    assert got["AAGGWWR"] == ("entrapment", "")
    assert got["QQWWLLDR"][0] == "excluded_contaminant"


def test_no_discovered_peptide_string_in_entrapment_json(
    computed: tuple[Inputs, Path, dict[str, Any]],
) -> None:
    """Data rule (statistician decisions, criterion 10): aggregates only."""
    _, out, _ = computed
    text = (out / "entrapment.json").read_text(encoding="utf-8")
    leaked = sorted(s for s in accepted_strings(EXAMPLE2_ROWS) if s in text)
    assert not leaked, f"row-level DIA-NN content in entrapment.json: {leaked}"
    assert "UniMod" not in text


# ------------------------------------------------------------------ criterion 8: determinism


def test_two_computes_are_byte_identical(inp: Inputs, tmp_path: Path) -> None:
    a, b = tmp_path / "m1", tmp_path / "m2"
    _ok(run_entrapment(*inp.args(a)))
    _ok(run_entrapment(*inp.args(b)))
    for name in ("entrapment.json", "entrapment_discoveries.tsv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


# ------------------------------------------------------------------ criterion 5: settings


def test_diann_cmd_equal_to_the_defaults_is_accepted(inp: Inputs, tmp_path: Path) -> None:
    cmd = tmp_path / "command_step1.txt"
    cmd.write_text(
        "\n".join(diann_step1_tokens(CFG_A1, inp.k12, inp.cont)) + "\n", encoding="utf-8"
    )
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--diann-cmd", str(cmd))))
    s = read_json(out / "entrapment.json")["settings"]
    assert (s["missed"], s["min_len"], s["max_len"]) == (1, 7, 30)


@pytest.mark.parametrize(
    "override",
    [{"--missed-cleavages": "2"}, {"--max-pep-len": "40"}, {"--max-pr-charge": "3"},
     {"--min-pr-mz": "350"}, {"--cut": "K*,R*,!*P"}],
)  # fmt: skip
def test_diann_cmd_differing_from_the_defaults_exits_2(
    inp: Inputs, tmp_path: Path, override: dict[str, str]
) -> None:
    cmd = tmp_path / "command_step1.txt"
    toks = diann_step1_tokens(CFG_A1, inp.k12, inp.cont, **override)
    cmd.write_text("\n".join(toks) + "\n", encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--diann-cmd", str(cmd))), out)


# ------------------------------------------------------------------ criterion 7: preconditions


def test_example_10_non_wild_type_run_exits_2(tmp_path: Path) -> None:
    rows = [*EXAMPLE2_ROWS, row("WT_ROS_1", "MSAGLLK", 2, 0.001)]
    inp = build_inputs(tmp_path / "in", rows)
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out, "WT_ROS_1")


def test_fasta_hash_differing_from_provenance_exits_2(inp: Inputs, tmp_path: Path) -> None:
    prov = tmp_path / "provenance.json"
    good = {"inputs": {"fasta": [{"path": str(inp.k12), "sha256": sha256(inp.k12)},
                                 {"path": str(inp.cont), "sha256": sha256(inp.cont)}]}}  # fmt: skip
    prov.write_text(json.dumps(good), encoding="utf-8")
    _ok(run_entrapment(*inp.args(tmp_path / "ok", "--provenance", str(prov))))
    bad = json.loads(json.dumps(good))
    bad["inputs"]["fasta"][0]["sha256"] = "0" * 64
    prov.write_text(json.dumps(bad), encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--provenance", str(prov))), out)


@pytest.mark.parametrize("name", ["pep_a1", "k12", "cont", "genome", "strain_db", "manifest"])
def test_missing_input_exits_2(inp: Inputs, tmp_path: Path, name: str) -> None:
    getattr(inp, name).unlink()
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out)


def test_genome_without_topology_tags_exits_2(tmp_path: Path) -> None:
    inp = build_inputs(tmp_path / "in", tags=False)
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out)


def test_empty_k12_fasta_exits_2(inp: Inputs, tmp_path: Path) -> None:
    inp.k12.write_text("", encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out)


def test_out_inside_the_checkout_exits_2(inp: Inputs) -> None:
    out = REPO / "validation" / "pilot" / "_m10_07_out_inside_checkout"
    assert not out.exists()
    try:
        _exit2_nothing(run_entrapment(*inp.args(out)), out)
    finally:
        if out.exists():
            shutil.rmtree(out)


# ------------------------------------------------------------------ cross-checks


def test_cross_check_1_exact_route_disagreement_exits_1(inp: Inputs, tmp_path: Path) -> None:
    """A1's peptide_hits.tsv has a route = exact row for AAGGWWR, which is not present here."""
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK", "AAGGWWR"], 0, 3)
    out = tmp_path / "metrics"
    res = run_entrapment(*inp.args(out, "--out-a1", str(a1)))
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def test_cross_checks_pass_quietly(inp: Inputs, tmp_path: Path) -> None:
    """Exact rows only for present peptides, and a genome-inconsistent count (2: AGNVLTEPWR and
    AAGGWWR) >= N_E: exit 0 and no warning beyond the run without --out-a1."""
    base = tmp_path / "base"
    _ok(run_entrapment(*inp.args(base)))
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK"], 2, 4)
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    doc = read_json(out / "entrapment.json")
    assert codes(doc["warnings"]) == codes(read_json(base / "entrapment.json")["warnings"])
    assert doc["cross_checks"]


def test_cross_check_2_n_e_above_the_genome_inconsistent_count_warns(
    inp: Inputs, tmp_path: Path
) -> None:
    """The A1 run lists only MSAGLLK and VEPTIDEK (exact), so its genome-inconsistent count is 0
    while N_E = 1: a warning, exit 0."""
    base = tmp_path / "base"
    _ok(run_entrapment(*inp.args(base)))
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK"], 0, 2)
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    doc = read_json(out / "entrapment.json")
    assert len(doc["warnings"]) > len(read_json(base / "entrapment.json")["warnings"])


def test_cross_check_3_digest_mismatch_warns_and_exits_0(tmp_path: Path) -> None:
    inp = build_inputs(tmp_path / "in", [*EXAMPLE2_ROWS, EXAMPLE8_ROW])
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out)))
    doc = read_json(out / "entrapment.json")
    assert "digest_mismatch" in codes(doc["warnings"])
    assert doc["peptide"]["not_in_universe"] == 1


# ------------------------------------------------------------------ criterion 9: the hook


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from helpers.pilot_metrics_m10 import sim_world

    return sim_world(tmp_path_factory.mktemp("m10_07_sim"))


def _pending_section(report: str) -> str:
    m = re.search(r"^## [^\n]*Pending in pass 2[^\n]*\n(.*?)(?=^## |\Z)", report, re.M | re.S)
    assert m, "the report has no 'Pending in pass 2' section"
    return m.group(1)


def _metrics_and_report(world: Any, root: Path, *extra: str) -> tuple[dict[str, Any], str]:
    from helpers.pilot_metrics_m10 import build_pilot, run_cli

    p = build_pilot(root / "pilot", world)
    out = root / "metrics"
    res = run_cli(*p.compute_args(out, *extra), cwd=p.root)
    _ok(res)
    report = root / "REPORT.md"
    _ok(run_cli("render", "--metrics", str(out / "metrics.json"), "--out", str(report)))
    return read_json(out / "metrics.json"), report.read_text(encoding="utf-8")


def test_hook_copies_the_blocks_under_m8_and_render_resolves_them(
    world: Any, inp: Inputs, tmp_path: Path
) -> None:
    ent_dir = tmp_path / "ent"
    _ok(run_entrapment(*inp.args(ent_dir)))
    ent = read_json(ent_dir / "entrapment.json")
    metrics, report = _metrics_and_report(
        world, tmp_path, "--entrapment", str(ent_dir / "entrapment.json")
    )
    m8 = metrics["m8"]
    for key in ("peptide", "precursor", "universe", "caveats", "cross_checks"):
        assert m8[key] == ent[key], key
    assert "{{" not in report
    for caveat in CAVEATS:
        assert caveat in report, caveat
    assert "metric 8" not in _pending_section(report).lower()
    assert "evidence_of_failure" in report or "evidence of failure" in report.lower()


def test_without_entrapment_metric_8_stays_pending(world: Any, tmp_path: Path) -> None:
    metrics, report = _metrics_and_report(world, tmp_path)
    assert "metric 8" in _pending_section(report).lower()
    m8 = metrics.get("m8")
    assert m8 is None or "peptide" not in m8


def test_template_has_a_metric_8_section_and_an_h2_line_without_verdict() -> None:
    from helpers.pilot_metrics_m10 import TEMPLATE

    text = TEMPLATE.read_text(encoding="utf-8")
    paths = set(re.findall(r"\{\{\s*([^}\s]+)\s*\}\}", text))
    for need in ("m8.peptide.fdp_combined", "m8.peptide.fdp_lower", "m8.precursor.assessment",
                 "m8.precursor.fdp_combined", "m8.universe.r"):  # fmt: skip
        assert need in paths, need
    h2 = [line for line in text.splitlines() if "**H2**" in line]
    assert h2, "no H2 line"
    line = h2[0]
    assert "deferred" in line
    assert "m8.peptide.fdp_combined" in line and "m2.A1.genome_inconsistent" in line
    assert "verdict" not in line.lower()
