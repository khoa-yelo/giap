"""M11-14w ``genetic_code_check``: examples 6-10, the ``_decide`` boundary, the ``info`` paths,
AC 1, AC 3 (fixture seeds), AC 4 (determinism, revcomp and record order) and AC 5 (preflight
wiring) (docs/specs/M11-14w-wrong-genetic-code-warning.md; GIAP_SPEC.md §7.0, M1-12).

Interface assumed (spec "Where it runs" items 1, 5 and example 8):
- ``giap.preflight.genetic_code.genetic_code_check(seqs: Sequence[str], table: int) -> Check``;
- ``giap.preflight.genetic_code._decide(bp11, bp4, n4)`` returns a ``Check`` or a level string
  (both accepted here; only the level is read from it);
- ``genetic_code_check`` gets its per-record counts from the module-level ``long_orf_cover``
  looked up at call time (the truncated-ratio wording test patches it; spec example 8 asks for
  ``ratio 0.54`` at bp11 = 549, bp4 = 1000, which no small real sequence gives);
- ``run_preflight`` calls ``genetic_code_check`` as a name it can see at call time in
  ``giap.preflight.run`` or ``giap.preflight.genetic_code`` (AC 5 patches both).

The fixture generator is ``helpers.genetic_code_fixture.synthetic_genome`` (spec test plan).
"""

from __future__ import annotations

import dataclasses
import itertools
import json
import re
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import Any

import jsonschema
import pandas as pd
import pytest
from helpers.genetic_code_fixture import random_dna, revcomp, synthetic_genome
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig, resolve_config
from giap.exitcodes import ExitCode
from giap.preflight.checks import Check

RULE_ID = "genome.genetic_code"
RULE_TEXT = (
    "long ORFs (>= 300 codons, six frames): bp with TGA as stop / bp with TGA as sense < 0.55 "
    "-> warn (needs >= 20 such ORFs)"
)
MSG_TABLE = "not evaluated: annotation.translation_table is not 11"
MSG_FEW = "not evaluated: fewer than 20 open reading frames of at least 300 codons"
WARN_START = "This genome probably uses genetic code 4"
VALUE_KEYS = {"min_codons", "bp_tga_stop", "bp_tga_sense", "n_orfs_tga_sense", "ratio", "threshold"}
ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "src" / "giap" / "schemas" / "preflight.schema.json"
LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"

CODE4 = [(0.28, 0.01), (0.40, 0.005)]
CODE11 = [0.28, 0.50, 0.70]
SEEDS = [1, 2, 3, 4, 5]


def gcc(seqs: Sequence[str], table: int = 11) -> Check:
    from giap.preflight.genetic_code import genetic_code_check

    c = genetic_code_check(seqs, table)
    assert isinstance(c, Check)
    return c


def level_of(result: Any) -> str:
    return result if isinstance(result, str) else str(result.level)


@cache
def code4(seed: int = 7, gc: float = 0.28, rate: float = 0.01) -> str:
    return synthetic_genome(seed, gc, 4, tga_rate=rate)


@cache
def code11(seed: int = 7, gc: float = 0.50) -> str:
    return synthetic_genome(seed, gc, 11)


def assert_shape(c: Check) -> None:
    assert c.rule_id == RULE_ID
    assert c.input == "genome"
    assert c.level in {"ok", "info", "warn"}
    assert c.exit_code == 0
    assert c.rule == RULE_TEXT


def assert_value(c: Check) -> dict[str, Any]:
    v = c.value
    assert isinstance(v, dict)
    assert set(v) == VALUE_KEYS
    assert v["min_codons"] == 300
    assert v["threshold"] == 0.55
    bp11, bp4 = v["bp_tga_stop"], v["bp_tga_sense"]
    assert isinstance(bp11, int) and isinstance(bp4, int) and isinstance(v["n_orfs_tga_sense"], int)
    assert 0 <= bp11 <= bp4
    if bp4 == 0:
        assert v["ratio"] is None
    else:
        assert v["ratio"] == pytest.approx(bp11 / bp4)
    return v


# ---------------------------------------------------------------- examples 6-7


def test_example6_code4_fixture_warns() -> None:
    c = gcc([code4()])
    assert_shape(c)
    v = assert_value(c)
    assert c.level == "warn"
    assert v["ratio"] is not None and v["ratio"] < 0.55
    assert v["n_orfs_tga_sense"] >= 20
    assert c.message.startswith(WARN_START)


def test_warning_wording() -> None:
    seq = code4()
    c = gcc([seq])
    v = assert_value(c)
    bp11, bp4 = v["bp_tga_stop"], v["bp_tga_sense"]
    text = " ".join(c.message.split())
    ratio = f"0.{(100 * bp11) // bp4:02d}"
    assert f"(ratio {ratio} < 0.55)" in text
    m = re.search(r"cover (\d+) % of the genome when TGA is read as a codon but only (\d+) % when "
                  r"it is a stop", text)  # fmt: skip
    assert m, text
    q4, q11 = int(m.group(1)), int(m.group(2))
    assert abs(q4 - 100 * bp4 / len(seq)) <= 1
    assert abs(q11 - 100 * bp11 / len(seq)) <= 1
    for phrase in (
        "(TGA = Trp; Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma)",
        "or another code in which TGA is not a stop",
        "GIAP 1.0 supports genetic code 11 only: genes will be cut at TGA",
        "proteins and peptides that contain W will not map",
        "strain_difference and nterm_extension_evidence flags will be artefacts",
        "This warning does not change any result.",
    ):
        assert phrase in text


def test_example7_code11_fixture_ok() -> None:
    c = gcc([code11()])
    assert_shape(c)
    v = assert_value(c)
    assert c.level == "ok"
    assert c.message == ""
    assert v["ratio"] > 0.9
    assert v["n_orfs_tga_sense"] >= 20


# ---------------------------------------------------------------- example 8: the boundary


def test_example8_decide_boundary_on_integers() -> None:
    from giap.preflight.genetic_code import _decide

    assert level_of(_decide(550, 1000, 20)) == "ok"
    assert level_of(_decide(549, 1000, 20)) == "warn"
    assert level_of(_decide(549, 1000, 19)) == "info"
    assert level_of(_decide(0, 1000, 19)) == "info"
    assert level_of(_decide(1000, 1000, 19)) == "info"
    assert level_of(_decide(0, 18_000, 20)) == "warn"
    assert level_of(_decide(18_000, 18_000, 20)) == "ok"


@settings(max_examples=300, deadline=None)
@given(bp4=st.integers(18_000, 10**9), frac=st.integers(0, 10**6), n4=st.integers(20, 10**5))
def test_decide_is_the_exact_integer_rule(bp4: int, frac: int, n4: int) -> None:
    from giap.preflight.genetic_code import _decide

    bp11 = bp4 * frac // 10**6
    want = "warn" if 20 * bp11 < 11 * bp4 else "ok"
    assert level_of(_decide(bp11, bp4, n4)) == want


def test_example8_ratio_is_truncated_in_the_message(monkeypatch: pytest.MonkeyPatch) -> None:
    import giap.preflight.genetic_code as gc_mod

    def fake(seq: str, stops: Any, min_codons: int) -> tuple[int, int]:
        assert min_codons == 300
        return (549, 20) if "TGA" in set(stops) else (1000, 20)

    monkeypatch.setattr(gc_mod, "long_orf_cover", fake)
    c = gcc(["ACGT" * 500])
    assert c.level == "warn"
    assert "(ratio 0.54 < 0.55)" in " ".join(c.message.split())
    assert "0.55 < 0.55" not in c.message
    v = assert_value(c)
    assert (v["bp_tga_stop"], v["bp_tga_sense"], v["n_orfs_tga_sense"]) == (549, 1000, 20)

    def at_cut(seq: str, stops: Any, min_codons: int) -> tuple[int, int]:
        return (550, 20) if "TGA" in set(stops) else (1000, 20)

    monkeypatch.setattr(gc_mod, "long_orf_cover", at_cut)
    assert gcc(["ACGT" * 500]).level == "ok"


# ---------------------------------------------------------------- examples 9-10, info paths


@pytest.mark.parametrize("table", [4, 1, 25])
@pytest.mark.parametrize("which", ["code4", "code11", "empty"])
def test_example9_other_table_is_info(table: int, which: str) -> None:
    seqs = {"code4": [code4()], "code11": [code11()], "empty": []}[which]
    c = gcc(seqs, table)
    assert_shape(c)
    assert c.level == "info"
    assert c.message == MSG_TABLE


def test_example10_plasmid_with_few_long_orfs_is_info() -> None:
    genes = synthetic_genome(3, 0.5, 11, n_genes=5)
    plasmid = random_dna(11, 20_000) + genes + random_dna(12, 30_000 - len(genes))
    c = gcc([plasmid])
    assert_shape(c)
    v = assert_value(c)
    assert v["n_orfs_tga_sense"] < 20
    assert c.level == "info"
    assert c.message == MSG_FEW


def test_example10_poly_a_is_info_not_ok() -> None:
    c = gcc(["A" * 10_000])
    assert_shape(c)
    v = assert_value(c)
    assert c.level == "info"
    assert c.message == MSG_FEW
    assert v == {"min_codons": 300, "bp_tga_stop": 10_000, "bp_tga_sense": 10_000,
                 "n_orfs_tga_sense": 6, "ratio": 1.0, "threshold": 0.55}  # fmt: skip


def test_no_records_is_info_with_null_ratio() -> None:
    c = gcc([])
    assert_shape(c)
    assert c.level == "info"
    assert c.message == MSG_FEW
    v = assert_value(c)
    assert v["bp_tga_stop"] == v["bp_tga_sense"] == v["n_orfs_tga_sense"] == 0
    assert v["ratio"] is None


@pytest.mark.parametrize("seqs", [["N" * 30_000], ["N" * 5000, "RYKMSW" * 1000], [""]])
def test_only_non_acgt_is_info_with_null_ratio(seqs: list[str]) -> None:
    c = gcc(seqs)
    assert c.level == "info"
    v = assert_value(c)
    assert v["n_orfs_tga_sense"] == 0 and v["ratio"] is None


def test_runs_never_cross_records() -> None:
    whole = gcc(["GCC" * 300])
    split = gcc(["GCC" * 150, "GCC" * 150])
    assert assert_value(whole)["bp_tga_sense"] >= 900
    assert assert_value(split)["bp_tga_sense"] == 0


def test_records_are_summed() -> None:
    from giap.preflight.genetic_code import long_orf_cover

    seqs = [code4(1)[:200_000], code11(2)[:150_000], "A" * 3000]
    v = assert_value(gcc(seqs))
    s11 = frozenset({"TAA", "TAG", "TGA"})
    s4 = frozenset({"TAA", "TAG"})
    per11 = [long_orf_cover(s, s11, 300) for s in seqs]
    per4 = [long_orf_cover(s, s4, 300) for s in seqs]
    assert v["bp_tga_stop"] == sum(b for b, _ in per11)
    assert v["bp_tga_sense"] == sum(b for b, _ in per4)
    assert v["n_orfs_tga_sense"] == sum(n for _, n in per4)


def test_lowercase_input_identical() -> None:
    a, b = gcc([code4()]), gcc([code4().lower()])
    assert a == b


def test_n_gaps_do_not_change_the_decision() -> None:
    """Spec "Thresholds": inserting N gaps does not change the decision (100 N every 10 kb)."""
    for seq, want in ((code4(), "warn"), (code11(), "ok")):
        gapped = ("N" * 100).join(seq[i : i + 10_000] for i in range(0, len(seq), 10_000))
        assert gcc([gapped]).level == want


# ---------------------------------------------------------------- AC 3: fixture seeds


@pytest.mark.parametrize(("gc", "rate"), CODE4)
@pytest.mark.parametrize("seed", SEEDS)
def test_ac3_code4_fixture_warns(seed: int, gc: float, rate: float) -> None:
    c = gcc([code4(seed, gc, rate)])
    assert c.level == "warn", c.value
    assert assert_value(c)["n_orfs_tga_sense"] >= 20


@pytest.mark.parametrize("gc", CODE11)
@pytest.mark.parametrize("seed", SEEDS)
def test_ac3_code11_fixture_does_not_warn(seed: int, gc: float) -> None:
    c = gcc([code11(seed, gc)])
    assert c.level == "ok", c.value


# ---------------------------------------------------------------- AC 4: determinism, symmetry


def _contigs(seq: str, cuts: Sequence[int]) -> list[str]:
    edges = [0, *cuts, len(seq)]
    return [seq[a:b] for a, b in itertools.pairwise(edges)]


def test_ac4_deterministic() -> None:
    seqs = _contigs(code4(2), [90_000, 200_000])
    assert gcc(seqs) == gcc(list(seqs))


@pytest.mark.parametrize("make", [code4, code11])
def test_ac4_revcomp_and_reordering_leave_value_unchanged(make: Any) -> None:
    seqs = _contigs(make(4), [60_000, 61_000, 250_000])
    base = gcc(seqs)
    assert gcc([revcomp(s) for s in seqs]).value == base.value
    assert gcc(list(reversed(seqs))).value == base.value
    assert gcc([seqs[2], seqs[0], seqs[3], seqs[1]]).value == base.value
    assert gcc([revcomp(s) for s in reversed(seqs)]).level == base.level


@settings(max_examples=40, deadline=None)
@given(
    seqs=st.lists(
        st.lists(st.sampled_from(["GCC" * 100, "GCA" * 50, "TGA", "TAA", "N", "A", "ATG"]),
                 max_size=40).map("".join),
        max_size=4,
    ),
    data=st.data(),
)  # fmt: skip
def test_ac4_property_revcomp_and_permutation(seqs: list[str], data: st.DataObject) -> None:
    base = gcc(seqs)
    assert_value(base)
    perm = data.draw(st.permutations(seqs))
    flips = data.draw(st.lists(st.booleans(), min_size=len(seqs), max_size=len(seqs)))
    other = [revcomp(s) if f else s for s, f in zip(perm, flips, strict=True)]
    assert gcc(other).value == base.value
    assert gcc(other).level == base.level


# ---------------------------------------------------------------- AC 1, AC 5: preflight wiring


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("db", "repo"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db"))
    monkeypatch.setenv("GIAP_REPO", str(tmp_path / "repo"))
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "repo"))  # never the real $GIAP_ROOT
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)


def _cfg(*extra: str) -> GiapConfig:
    return resolve_config(overrides=[*LOW, *extra])


def _probe() -> Any:
    from giap.preflight.environment import GIB, EnvProbe, requirements

    found: dict[str, str | None] = {
        r.name: r.min_version or "1.0"
        for r in requirements(["core", "mobile", "bgc", "function"], "genome")
    }
    return EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)


def _inputs(tmp_path: Path, fasta: str, annotation: bool = False) -> Any:
    from giap.preflight.run import PreflightInputs

    g = tmp_path / "g.fa"
    g.write_text(fasta)
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    pep = tmp_path / "peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    ann = None
    if annotation:
        ann = tmp_path / "bakta.json"
        ann.write_text("{}")
    return PreflightInputs(genome=g, proteins=p, peptides=pep, annotation=ann)


def _fasta(seqs: Sequence[str]) -> str:
    return "".join(f">c{i + 1}\n{s}\n" for i, s in enumerate(seqs))


def _run(tmp_path: Path, fasta: str, *extra: str, annotation: bool = False) -> Any:
    from giap.preflight.run import run_preflight

    return run_preflight(_inputs(tmp_path, fasta, annotation), _cfg(*extra), _probe(), "preflight")


def _gc_checks(report: Any) -> list[Check]:
    return [c for c in report.checks if c.rule_id == RULE_ID]


@pytest.mark.parametrize(
    ("seqs", "want"),
    [([code4()], "warn"), ([code11()], "ok"), (["A" * 10_000], "info")],
)
def test_ac1_preflight_has_the_check_once_after_the_genome_checks(
    tmp_path: Path, seqs: list[str], want: str
) -> None:
    report, _ = _run(tmp_path, _fasta(seqs))
    found = _gc_checks(report)
    assert len(found) == 1
    c = found[0]
    assert_shape(c)
    assert c.level == want
    i = report.checks.index(c)
    genome_idx = [j for j, x in enumerate(report.checks) if x.input == "genome" and x is not c]
    assert genome_idx, "no genome checks"
    assert i == max(genome_idx) + 1, [x.rule_id for x in report.checks]
    assert c == gcc(seqs)


def test_ac1_parse_error_gives_info(tmp_path: Path) -> None:
    report, _ = _run(tmp_path, "this is not FASTA\n")
    found = _gc_checks(report)
    assert len(found) == 1
    assert found[0].level == "info"
    assert found[0].exit_code == 0


def test_ac1_translation_table_4_gives_info(tmp_path: Path) -> None:
    report, _ = _run(tmp_path, _fasta([code4()]), "annotation.translation_table=4")
    (c,) = _gc_checks(report)
    assert c.level == "info"
    assert c.message == MSG_TABLE


def test_evaluated_with_precomputed_annotation(tmp_path: Path) -> None:
    report, _ = _run(tmp_path, _fasta([code4()]), annotation=True)
    (c,) = _gc_checks(report)
    assert c.level == "warn"


def test_topology_is_ignored(tmp_path: Path) -> None:
    fasta = _fasta([code4()])
    a, _ = _run(tmp_path, fasta, "annotation.bakta.complete=true")
    b, _ = _run(tmp_path, fasta, "annotation.bakta.complete=false")
    assert _gc_checks(a) == _gc_checks(b)


def test_preflight_json_and_md_carry_the_check(tmp_path: Path) -> None:
    from giap.preflight.run import write_outputs

    report, norm = _run(tmp_path, _fasta([code4()]))
    out = tmp_path / "OUT"
    out.mkdir()
    write_outputs(out, report, norm)
    data = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(data)
    rows = [c for c in data["checks"] if c["rule_id"] == RULE_ID]
    assert len(rows) == 1
    row = rows[0]
    assert row["input"] == "genome" and row["level"] == "warn" and row["exit_code"] == 0
    assert set(row["value"]) == VALUE_KEYS
    assert 0 <= row["value"]["ratio"] < 0.55
    assert RULE_ID in (out / "preflight" / "preflight.md").read_text(encoding="utf-8")


def _ok_check() -> Check:
    value = {"min_codons": 300, "bp_tga_stop": 0, "bp_tga_sense": 0, "n_orfs_tga_sense": 0,
             "ratio": None, "threshold": 0.55}  # fmt: skip
    return Check(RULE_ID, "genome", "ok", RULE_TEXT, value, "")


def _same(a: Any, b: Any) -> bool:
    if isinstance(a, pd.DataFrame) or isinstance(b, pd.DataFrame):
        return isinstance(a, pd.DataFrame) and isinstance(b, pd.DataFrame) and a.equals(b)
    return bool(a == b)


def test_ac5_warning_changes_only_checks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import giap.preflight.genetic_code as gc_mod
    import giap.preflight.run as run_mod
    from giap.preflight.run import exit_code

    fasta = _fasta([code4()])
    warned, norm_w = _run(tmp_path, fasta)
    assert [c.level for c in _gc_checks(warned)] == ["warn"]

    monkeypatch.setattr(gc_mod, "genetic_code_check", lambda seqs, table: _ok_check())
    monkeypatch.setattr(run_mod, "genetic_code_check", lambda seqs, table: _ok_check(),
                        raising=False)  # fmt: skip
    quiet, norm_q = _run(tmp_path, fasta)
    assert [c.level for c in _gc_checks(quiet)] == ["ok"]

    assert exit_code(warned, False) == exit_code(quiet, False) == ExitCode.SUCCESS
    for f in dataclasses.fields(norm_w):
        if f.name == "resources":
            continue
        assert _same(getattr(norm_w, f.name), getattr(norm_q, f.name)), f.name
    assert norm_w.cfg is not None and norm_w.cfg.annotation.translation_table == 11
    rest_w = [c for c in warned.checks if c.rule_id != RULE_ID]
    rest_q = [c for c in quiet.checks if c.rule_id != RULE_ID]
    assert rest_w == rest_q
    for f in dataclasses.fields(warned):
        if f.name != "checks":
            assert _same(getattr(warned, f.name), getattr(quiet, f.name)), f.name


def test_ac5_strict_turns_the_warning_into_exit_3(tmp_path: Path) -> None:
    from giap.preflight.run import exit_code

    report, _ = _run(tmp_path, _fasta([code4()]))
    others = [c for c in report.checks if c.level == "warn" and c.rule_id != RULE_ID]
    assert not others, [c.rule_id for c in others]
    assert exit_code(report, False) == ExitCode.SUCCESS
    assert exit_code(report, True) == ExitCode.WARNINGS_STRICT
    ok, _ = _run(tmp_path, _fasta([code11()]))
    assert exit_code(ok, True) == ExitCode.SUCCESS
