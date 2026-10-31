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
--contaminants --genome --strain-db --annotation --runs-manifest --out`` (required; ADR-0068:
``--annotation`` is a Bakta JSON, possibly gzipped) and ``--out-a1``,
``--diann-cmd``, ``--provenance`` (optional), mirroring the spec's input names. ``--out`` is
METRICS_DIR. ``provenance.json`` lists the searched FASTAs as ``inputs.fasta = [{path, sha256}]``
(M10-02 "inputs: each --fasta (path, sha256)"). ``command_step1.txt`` is one argv token per line.
``metrics.py compute`` takes ``--entrapment PATH`` and writes the blocks under ``m8``.

Review round rev-1007 (the spec's binding "Lead decisions after the code review"):

* exit codes: 0 written, 1 cross-check 1 only, 2 a precondition (including an unreadable or
  malformed input: provenance JSON, a missing peptide-table column, an NA charge on an accepted
  row; the message names the file), 4 any unexpected exception (traceback on stderr; listed in
  the module docstring); nothing is written unless 0;
* cross-check 1 looks only at the ``target`` and ``entrapment`` buckets and requires
  ``present_rule == "exact"``; discoveries outside those buckets with an A1 ``route = exact``
  row are counted in ``cross_checks.run.n_exact_route_unclassified``;
* cross-check 2 without ``tables/detection_summary.tsv`` warns ``genome_inconsistent_unavailable``;
  ``cross_checks.run`` keeps ``present_not_exact``, ``exact_rejected``,
  ``genome_inconsistent_count`` and ``n_entrapment``;
* ``--out-a1`` is resolved through ``ATTEMPTS.tsv`` next to it, as M10-05's ``resolve_attempt``;
* provenance: the record whose basename is ``UP000000625.fasta`` vouches for ``--k12-fasta``, the
  one named ``contaminants.fasta`` for ``--contaminants`` (any directory, any order);
* ``--out`` inside the checkout is accepted when ``git check-ignore`` says it is ignored;
* the spec SHA comes from git or the run exits 2 (no placeholder);
* ``metrics.json`` ``m8`` also carries ``warnings``; the report renders the warning codes and the
  ``ci.assumes`` string; ``inputs.entrapment`` is hashed like the other inputs.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from helpers.entrapment_m10_07 import (
    CAVEATS,
    CI_ASSUMES,
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
    write_shared_annotation,
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
    annotation: Path

    def args(self, out: Path, *extra: str, **over: Path) -> list[str]:
        v = {"pep-a1": self.pep_a1, "k12-fasta": self.k12, "contaminants": self.cont,
             "genome": self.genome, "strain-db": self.strain_db, "annotation": self.annotation,
             "runs-manifest": self.manifest}  # fmt: skip
        v |= {k.replace("_", "-"): p for k, p in over.items()}
        flat = [x for k, p in v.items() for x in (f"--{k}", str(p))]
        return ["compute", *flat, "--out", str(out), *extra]


def build_inputs(
    root: Path,
    rows: Any = EXAMPLE2_ROWS,
    *,
    tags: bool = True,
    k12: tuple[str, ...] = K12,
    texts: tuple[str, ...] = NISSLE,
    strain_db: tuple[str, ...] = ("MSAGLLKVEPTLDEKAAGGYWR",),
) -> Inputs:
    return Inputs(
        root=root,
        pep_a1=write_a1(root / "inputs" / "a1.tsv.gz", list(rows)),
        k12=write_faa(root / "data" / "UP000000625.fasta", k12),
        cont=write_faa(root / "db" / "contaminants.fasta", CONTAMINANTS, prefix="Cont_P"),
        genome=write_genome(
            root / "nissle.fna", {"nissle": coding_text(*texts)}, ["nissle"] if tags else []
        ),
        strain_db=write_faa(root / "strain_db.faa", list(strain_db), prefix="NISSLE_"),
        manifest=write_runs_manifest(root / "pxd078459_files.tsv"),
        annotation=write_shared_annotation(
            root / "nissle.bakta.json", texts, ["nissle"] if tags else []
        ),
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
    assert doc["schema"] == "gippa-pilot-entrapment/1"
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


@pytest.mark.parametrize(
    "name", ["pep_a1", "k12", "cont", "genome", "strain_db", "manifest", "annotation"]
)
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
    """The final template's "Not evaluated" section (M10-15 replaced "Pending in pass 2")."""
    m = re.search(r"^## [^\n]*Not evaluated[^\n]*\n(.*?)(?=^## |\Z)", report, re.M | re.S)
    assert m, "the report has no 'Not evaluated' section"
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
    assert "`m8`" not in _pending_section(report)
    assert "evidence_of_failure" in report or "evidence of failure" in report.lower()


def test_without_entrapment_metric_8_stays_pending(world: Any, tmp_path: Path) -> None:
    metrics, report = _metrics_and_report(world, tmp_path)
    assert "`m8`: pending (M10-07)" in _pending_section(report)
    m8 = metrics.get("m8")
    assert m8 is None or "universe" not in m8  # M10-15's pending placeholder has no universe


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


# ------------------------------------------------------------------ rev-1007 review round
# Binding: the spec's "Lead decisions after the code review (rev-1007)" and the statistician's
# S1 (``m8`` carries the warnings and the assumption string) and erratum (the per-run reason).


def _rewrite_a1(path: Path, edit: Any) -> Path:
    """Apply ``edit(header, rows) -> (header, rows)`` to the gzipped TSV at ``path``."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    header, rows = edit(lines[0].split("\t"), [line.split("\t") for line in lines[1:]])
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("\n".join("\t".join(r) for r in [header, *rows]) + "\n")
    return path


def _drop_column(name: str) -> Any:
    def edit(header: list[str], rows: list[list[str]]) -> tuple[list[str], list[list[str]]]:
        i = header.index(name)
        return header[:i] + header[i + 1 :], [r[:i] + r[i + 1 :] for r in rows]

    return edit


def _cross_run(doc: dict[str, Any]) -> dict[str, Any]:
    run = doc["cross_checks"].get("run")
    assert isinstance(run, dict), doc["cross_checks"]
    return run


# B1: cross-check 1 sees only the target and entrapment buckets


def test_b1_exact_route_on_unclassified_discoveries_is_not_a_disagreement(tmp_path: Path) -> None:
    """Example 8 plus A1 ``route = exact`` rows for its Met-excised SAGLLKVEPTIDEK
    (``not_in_universe``) and for the contaminant QQWWLLDR (``excluded_contaminant``): neither
    was classified by ``present``, so neither can disagree. Exit 0, the digest-mismatch warning,
    and the two are counted under ``cross_checks.run.n_exact_route_unclassified``."""
    inp = build_inputs(tmp_path / "in", [*EXAMPLE2_ROWS, EXAMPLE8_ROW])
    exact = ["MSAGLLK", "VEPTIDEK", "SAGLLKVEPTIDEK", "QQWWLLDR"]
    a1 = write_out_a1(tmp_path / "A1", exact, 2, 4)
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    doc = read_json(out / "entrapment.json")
    assert "digest_mismatch" in codes(doc["warnings"])
    assert doc["peptide"]["not_in_universe"] == 1
    assert _cross_run(doc)["n_exact_route_unclassified"] == 2


# S3: cross-check 1 requires the exact rule; strain_db is not agreement


def test_s3_strain_db_present_peptide_with_an_a1_exact_row_exits_1(tmp_path: Path) -> None:
    """Example 9's MTEKLLAGGPR is present here only by the ``strain_db`` clause (the genome reads
    VTEK...). GIPPA's exact route is six-frame only, so an A1 ``route = exact`` row for it means
    the two six-frame matchers disagree: exit 1, nothing written."""
    p4 = "MTEKLLAGGPR"
    rows = [*EXAMPLE2_ROWS, row("WT_1", p4, 2, 0.003)]
    kw: dict[str, Any] = {"k12": (*K12, p4), "texts": (*NISSLE, "VTEKLLAGGPR"),
                          "strain_db": ("MSAGLLKVEPTLDEKAAGGYWR", "MTEKLLAGGPRAAGK")}  # fmt: skip
    inp = build_inputs(tmp_path / "in", rows, **kw)
    ok = write_out_a1(tmp_path / "A1ok", ["MSAGLLK", "VEPTIDEK"], 2, 5)
    _ok(run_entrapment(*inp.args(tmp_path / "base", "--out-a1", str(ok))))
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK", p4], 2, 5)
    out = tmp_path / "metrics"
    res = run_entrapment(*inp.args(out, "--out-a1", str(a1)))
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)


# S4: cross-check 2's codes and breakdown


def test_s4_missing_detection_summary_has_its_own_warning_code(inp: Inputs, tmp_path: Path) -> None:
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK"], 2, 4)
    (a1 / "tables" / "detection_summary.tsv").unlink()
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    got = codes(read_json(out / "entrapment.json")["warnings"])
    assert "genome_inconsistent_unavailable" in got
    assert "n_entrapment_above_genome_inconsistent" not in got


@pytest.mark.parametrize(
    ("exact", "want"),
    [(["MSAGLLK", "VEPTIDEK"], (1, 0)), (["MSAGLLK"], (1, 1))],
)
def test_s4_cross_check_2_breakdown(
    inp: Inputs, tmp_path: Path, exact: list[str], want: tuple[int, int]
) -> None:
    """Example 2's targets: MSAGLLK and VEPTLDEK (exact here), AGNVLTEPWR (mass-equivalent).
    ``present_not_exact`` counts the targets made present by ``strain_db`` or rule 2 (1);
    ``exact_rejected`` the targets exact here without an A1 ``route = exact`` row (0, or 1 when
    A1 lists only MSAGLLK). Mutants M22/M23 (both 0) survived the round-1 tests."""
    a1 = write_out_a1(tmp_path / "A1", exact, 2, 4)
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    run = _cross_run(read_json(out / "entrapment.json"))
    assert (run["present_not_exact"], run["exact_rejected"]) == want
    assert run["genome_inconsistent_count"] == 2
    assert run["n_entrapment"] == 1


def test_n6_n_e_equal_to_the_genome_inconsistent_count_does_not_warn(
    inp: Inputs, tmp_path: Path
) -> None:
    """rev-1007 N6 (M10): the cross-check is N_E <= count; N_E = count = 1 is no violation."""
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK"], 1, 4)
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(a1))))
    doc = read_json(out / "entrapment.json")
    assert _cross_run(doc)["genome_inconsistent_count"] == 1
    assert "n_entrapment_above_genome_inconsistent" not in codes(doc["warnings"])


def test_s4_out_a1_without_peptide_hits_exits_2_naming_both_locations(
    inp: Inputs, tmp_path: Path
) -> None:
    a1 = tmp_path / "A1"
    (a1 / "tables").mkdir(parents=True)
    out = tmp_path / "metrics"
    res = run_entrapment(*inp.args(out, "--out-a1", str(a1)))
    _exit2_nothing(res, out, str(a1))
    assert "tables" in res.stderr and "mapping" in res.stderr, res.stderr


# S5: --out-a1 is resolved the way M10-05 resolves a run directory


def test_s5_out_a1_resolves_through_attempts_tsv(inp: Inputs, tmp_path: Path) -> None:
    """M10-05 (lead decision 1) resolves ``OUT_A1`` through M10-04's ``ATTEMPTS.tsv`` to the
    final attempt. The superseded ``out_a1`` holds an exact row for the entrapment peptide
    AAGGWWR (cross-check 1 would exit 1); the final ``out_a1.r1`` agrees: exit 0."""
    runs = tmp_path / "runs"
    write_out_a1(runs / "out_a1", ["MSAGLLK", "VEPTIDEK", "AAGGWWR"], 2, 4)
    write_out_a1(runs / "out_a1.r1", ["MSAGLLK", "VEPTIDEK"], 2, 4)
    cols = ("arm", "attempt", "outdir", "job", "resume", "git_commit", "state", "exit_code")
    att = [
        cols,
        ("A1", "1", str(runs / "out_a1"), "2001", "fresh", "abc1234", "superseded", "1"),
        ("A1", "2", str(runs / "out_a1.r1"), "2002", "resume", "abc1234", "done", "0"),
    ]
    (runs / "ATTEMPTS.tsv").write_text("\n".join("\t".join(r) for r in att) + "\n", "utf-8")
    out = tmp_path / "metrics"
    _ok(run_entrapment(*inp.args(out, "--out-a1", str(runs / "out_a1"))))


# S1: exit 2 for unreadable inputs (naming the file), exit 4 for anything unexpected


@pytest.mark.parametrize("content", ["{not json", "[]", '{"inputs": []}'])
def test_s1_malformed_provenance_exits_2_naming_it(
    inp: Inputs, tmp_path: Path, content: str
) -> None:
    prov = tmp_path / "provenance.json"
    prov.write_text(content, encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--provenance", str(prov))), out, str(prov))


@pytest.mark.parametrize("column", ["charge", "modified_peptide", "run", "q_value"])
def test_s1_pep_a1_without_a_needed_column_exits_2_naming_it(
    inp: Inputs, tmp_path: Path, column: str
) -> None:
    _rewrite_a1(inp.pep_a1, _drop_column(column))
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out, str(inp.pep_a1))


def test_s1_na_charge_on_an_accepted_row_exits_2_naming_pep_a1(inp: Inputs, tmp_path: Path) -> None:
    def edit(header: list[str], rows: list[list[str]]) -> tuple[list[str], list[list[str]]]:
        i, j = header.index("charge"), header.index("peptide")
        return header, [[*r[:i], "NA", *r[i + 1 :]] if r[j] == "MSAGLLK" else r for r in rows]

    _rewrite_a1(inp.pep_a1, edit)
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out, str(inp.pep_a1))


def test_s1_runs_manifest_without_a_run_column_exits_2_naming_it(
    inp: Inputs, tmp_path: Path
) -> None:
    text = inp.manifest.read_text(encoding="utf-8").replace("run\t", "sample_run\t", 1)
    inp.manifest.write_text(text, encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out)), out, str(inp.manifest))


def test_s1_unexpected_exception_exits_4_with_a_traceback(
    inp: Inputs, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Exit 1 is reserved for cross-check 1 (lead decision S1): a bug exits 4 and prints the
    traceback; nothing is written."""
    from helpers.entrapment_m10_07 import load_entrapment

    mod = load_entrapment()

    def boom(*_: Any, **__: Any) -> Any:
        raise RuntimeError("m10-07-unexpected")

    monkeypatch.setattr(mod, "evaluate", boom)
    out = tmp_path / "metrics"
    try:
        rc = mod.main(inp.args(out))
    except SystemExit as err:
        rc = err.code
    except RuntimeError:
        pytest.fail("the unexpected exception escaped main() (the interpreter exits 1)")
    err_text = capsys.readouterr().err
    assert rc == 4, rc
    assert "Traceback" in err_text and "m10-07-unexpected" in err_text
    assert _nothing_written(out)


def test_s1_docstring_lists_exit_4() -> None:
    from helpers.entrapment_m10_07 import load_entrapment

    doc = load_entrapment().__doc__ or ""
    assert re.search(r"\b4\b[^.;]*unexpected|unexpected[^.;]*\b4\b", doc, re.I | re.S), doc


# S2: provenance matched by basename role


def _prov(path: Path, records: list[tuple[str, str]]) -> Path:
    """M10-02's layout: ``inputs.fasta = [{path, sha256, size, mode}]``; the paths are the
    search's (another directory than the inputs here), so only the basename can match."""
    recs = [{"path": f"/scratch/gippa/db/{name}", "sha256": sha, "size": 1, "mode": "0644"}
            for name, sha in records]  # fmt: skip
    path.write_text(json.dumps({"inputs": {"fasta": recs}}), encoding="utf-8")
    return path


def test_s2_provenance_in_another_directory_and_order_is_accepted(
    inp: Inputs, tmp_path: Path
) -> None:
    prov = _prov(
        tmp_path / "provenance.json",
        [("contaminants.fasta", sha256(inp.cont)), ("UP000000625.fasta", sha256(inp.k12))],
    )
    _ok(run_entrapment(*inp.args(tmp_path / "metrics", "--provenance", str(prov))))


def test_s2_swapped_roles_exit_2(inp: Inputs, tmp_path: Path) -> None:
    """--k12-fasta and --contaminants swapped: both hashes are in the record, but each under the
    other role (rev-1007 S2, reproduced exit 0 in round 1)."""
    prov = _prov(
        tmp_path / "provenance.json",
        [("UP000000625.fasta", sha256(inp.k12)), ("contaminants.fasta", sha256(inp.cont))],
    )
    out = tmp_path / "metrics"
    args = inp.args(out, "--provenance", str(prov), k12_fasta=inp.cont, contaminants=inp.k12)
    _exit2_nothing(run_entrapment(*args), out)


@pytest.mark.parametrize("which", ["k12", "contaminants"])
def test_s2_either_recorded_hash_differing_exits_2(inp: Inputs, tmp_path: Path, which: str) -> None:
    """Mutant M17 (the contaminants hash never checked) survived round 1."""
    k, c = sha256(inp.k12), sha256(inp.cont)
    recs = [("UP000000625.fasta", "0" * 64 if which == "k12" else k),
            ("contaminants.fasta", "0" * 64 if which == "contaminants" else c)]  # fmt: skip
    prov = _prov(tmp_path / "provenance.json", recs)
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--provenance", str(prov))), out, str(prov))


def test_s2_two_copies_of_the_k12_hash_do_not_vouch_for_the_contaminants(
    inp: Inputs, tmp_path: Path
) -> None:
    prov = _prov(
        tmp_path / "provenance.json",
        [("UP000000625.fasta", sha256(inp.k12)), ("contaminants.fasta", sha256(inp.k12))],
    )
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--provenance", str(prov))), out)


def test_s2_a_role_missing_from_provenance_exits_2(inp: Inputs, tmp_path: Path) -> None:
    prov = _prov(tmp_path / "provenance.json", [("UP000000625.fasta", sha256(inp.k12))])
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--provenance", str(prov))), out)


# S7: the --met-excision and --unimod4 settings checks (mutants M18/M19)


def _edit_tokens(tokens: list[str], edit: str) -> list[str]:
    if edit == "add --met-excision":
        return [*tokens, "--met-excision"]
    assert edit == "drop --unimod4"
    assert "--unimod4" in tokens
    return [t for t in tokens if t != "--unimod4"]


@pytest.mark.parametrize("edit", ["add --met-excision", "drop --unimod4"])
def test_s7_diann_cmd_met_excision_or_no_unimod4_exits_2(
    inp: Inputs, tmp_path: Path, edit: str
) -> None:
    cmd = tmp_path / "command_step1.txt"
    toks = _edit_tokens(diann_step1_tokens(CFG_A1, inp.k12, inp.cont), edit)
    cmd.write_text("\n".join(toks) + "\n", encoding="utf-8")
    out = tmp_path / "metrics"
    _exit2_nothing(run_entrapment(*inp.args(out, "--diann-cmd", str(cmd))), out, str(cmd))


# B2: --out inside the checkout is accepted when git-ignored


def test_b2_out_under_an_ignored_directory_of_the_checkout_is_accepted(inp: Inputs) -> None:
    """``farmshare_env.sh`` puts ``GIPPA_ROOT`` at ``$GIPPA_REPO/local`` (git-ignored), so the
    sbatch default ``$GIPPA_ROOT/pilot/metrics`` lies inside the checkout; it must not exit 2
    (precedent: ``diann_convert.sbatch``'s ``git check-ignore``)."""
    base = REPO / "local"
    made_base = not base.exists()
    out = base / f"_m10_07_ignored_{uuid.uuid4().hex[:8]}" / "pilot" / "metrics"
    chk = subprocess.run(["git", "-C", str(REPO), "check-ignore", "-q", str(out)], check=False)
    assert chk.returncode == 0, f"{out} is not git-ignored in this checkout"
    try:
        _ok(run_entrapment(*inp.args(out)))
        assert (out / "entrapment.json").is_file()
    finally:
        shutil.rmtree(out.parents[1], ignore_errors=True)
        if made_base:
            shutil.rmtree(base, ignore_errors=True)


# N4: the spec SHA is never a placeholder


def test_n4_spec_sha_unavailable_exits_2_not_zeros(inp: Inputs, tmp_path: Path) -> None:
    """When git cannot give the spec's commit, M10-05's ``_spec_sha`` raises ``Precondition``;
    round 1 wrote 40 zeros instead."""
    out = tmp_path / "metrics"
    res = run_entrapment(*inp.args(out), env_extra={"GIT_DIR": str(tmp_path / "no-such-git")})
    _exit2_nothing(res, out)


# Statistician S1 (rev-1007): m8 carries the warnings, the report the assumption string


def test_stat_s1_m8_carries_the_warnings_and_the_report_the_assumes_string(
    world: Any, tmp_path: Path
) -> None:
    inp = build_inputs(tmp_path / "in", [*EXAMPLE2_ROWS, EXAMPLE8_ROW])
    ent_dir = tmp_path / "ent"
    _ok(run_entrapment(*inp.args(ent_dir)))
    ent = read_json(ent_dir / "entrapment.json")
    assert "digest_mismatch" in codes(ent["warnings"])
    metrics, report = _metrics_and_report(
        world, tmp_path, "--entrapment", str(ent_dir / "entrapment.json")
    )
    assert metrics["m8"]["warnings"] == ent["warnings"]
    assert "{{" not in report
    assert "digest_mismatch" in report or "digest mismatch" in report.lower()
    assert CI_ASSUMES in report


def test_b3_entrapment_input_is_hashed_like_the_others(
    world: Any, inp: Inputs, tmp_path: Path
) -> None:
    """B3: after the merge of main, ``entrapment`` is one of metrics.py's optional inputs and is
    recorded in ``metrics.json`` ``inputs`` with its sha256, like every other input."""
    ent_dir = tmp_path / "ent"
    _ok(run_entrapment(*inp.args(ent_dir)))
    ent = ent_dir / "entrapment.json"
    metrics, _ = _metrics_and_report(world, tmp_path, "--entrapment", str(ent))
    rec = metrics["inputs"]["entrapment"]
    assert rec["sha256"] == sha256(ent)
    assert Path(rec["path"]) == ent.resolve()
