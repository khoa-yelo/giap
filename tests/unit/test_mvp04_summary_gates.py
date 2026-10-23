"""The full ``scripts/dev/mvp04_summary.py``: acceptance gates 1-16 of
docs/specs/MVP-04-ecoli-evidence-run.md ("Acceptance criteria and expected numbers"; lead
decisions 2, 11, 15 and 18). Gate 7 keeps its part-1 tests in tests/unit/test_mvp04_summary.py;
``gate7()`` is unchanged.

Fixtures are hand-built (tests/helpers/mvp04_summary_data.py): ``build`` writes a run in which
every gate passes inside its expected range, and each test edits one file.

Interface (pinned by the test-author; the implementer may add keys, not rename these)
====================================================================================

Constants: ``ALPHA = 0.01`` (``stats.alpha.identification``), ``Q_VALUE_MAX = 0.01``
(``input.q_value_max``), ``FORBIDDEN_WARNINGS = ("decoy_unscored", "decoy_pairing_low",
"q_floor_above_alpha", "no_eligible_genes")``, ``REPORT_SECTIONS = ("detection",
"genome-context", "operons", "variants", "novel-orfs")``, ``RUNTIME_KEY = "mvp04_k12_evidence"``.

``status_of(gate_ok: bool, in_expected: bool) -> str``: ``"pass"`` (gate met, value in its
expected range), ``"outside_expected"`` (gate met, value outside the expected range: reported,
not failed, spec "Gates are what the e2e test asserts"), ``"fail"`` (gate not met).

``gate<N>(outdir: Path, inputs: Path) -> dict`` for N in 1-6, 8, 10-16, ``gate9(inputs: Path,
again: Path | None) -> dict`` and ``gate7_entry(outdir, inputs) -> dict`` (gate 7's entry: the
part-1 ``gate7()`` dict plus the keys below). Each is pure: it reads files, writes nothing and
never raises on a missing output file (the gate then fails). Every entry carries

* ``value``: the gate's main quantity (JSON scalar or list); ``expected``: the spec's expected
  range ``[lo, hi]`` (``None`` = open end) or a categorical value; ``gate``: the gate in words;
* ``gate_ok`` (bool), ``in_expected`` (bool), ``status`` = ``status_of(gate_ok, in_expected)``.

Eligible genes are ``tables/genes.tsv`` rows whose ``detection_status`` is not
``partial_excluded`` (§7.3, D8). A true gene is the ``locus_tag`` of a ``truth/peptides.tsv`` row
with ``is_wrong = false`` (decision 11). Per gate (extra keys in brackets):

1. ``run_manifest.json``: value = ``exit_code``; ok iff exit 0, ``status = completed`` and no
   ``steps[*].status = failed`` [``run_status``, ``failed_steps``]. Missing manifest: fail.
2. value = the one ``detection_tier`` of ``genes.tsv`` (``"mixed"`` if several); expected and gate
   ``decoy_picked`` and no ``FORBIDDEN_WARNINGS`` token found [``warnings_found``, sorted]. A token
   is found when it is a substring of a manifest ``flags``/``warnings`` entry or of a line of
   ``logs/events.jsonl``; ``q_floor_above_alpha`` is also found when the tier ``q_floor`` > ALPHA.
3. value = tier row ``decoy_paired_share`` (the run's tier); gate >= 0.8; expected [0.98, 1].
4. value = T = round(1 / tier ``q_floor``) [``q_floor``]; gate q_floor <= 1/110; expected
   [2300, 3100]. NA q_floor: fail.
5. value = decoy wins: eligible genes with ``gene_q_value`` NA that have a paired decoy, i.e. an
   input ``peptides.tsv`` row with ``is_decoy = true``, ``q_value <= Q_VALUE_MAX`` and
   ``rev_<gene>`` among its ``protein_ids`` [``n_paired_genes``: eligible genes with such a
   decoy]; gate >= 1; expected [20, 300].
6. value = tier ``n_detected`` [``n_eligible``, ``n_truth_detected``]; gate >= 1000 and <= 0.90 x
   eligible; expected [0.85, 1.0] x ``n_truth_detected``.
7. value = ``realised_gene_fdp``; gate <= 0.02; expected [0, 0.01]; NaN (nothing detected): fail.
8. value = eligible true genes / eligible [``n_truth_detected``, ``n_eligible``]; gate and
   expected [0.55, 0.70].
9. value = sorted relative paths that differ between the two giapsim calls (``peptides.tsv`` and
   ``truth/*.tsv``, a file present in one only counts as differing); ok iff empty and both have
   ``peptides.tsv`` and ``truth/{peptides,genes,islands,alt_sites,novel_orfs}.tsv``. ``again =
   None``: fail.
10. over the distinct peptides of ``mapping/peptide_hits.tsv`` and ``mapping/unmapped.tsv``
    (``kind = peptide``), decoy-only input peptides left out: exact = with an ``exact`` hit,
    tolerant = a ``tolerant`` hit and no exact one, unmapped = the rest. value = exact share
    [``shares``: ``exact``/``tolerant``/``unmapped``, summing to 1; ``n_peptides``]; gate
    tolerant share > 0; expected [0.95, 1].
11. ``high`` = ``variants.tsv`` rows with ``confidence = high``; a high row is true when
    (``contig``, ``codon_start``, ``alt_aa``) is a ``truth/alt_sites.tsv`` site (both 1-based
    lowest coordinates). value = ``n_high`` [``precision``, ``recall`` over truth sites with
    ``n_true_peptides >= 2``, ``n_recall_sites``, ``n_rows_missing_flags`` (rows of any confidence
    without both ``strain_difference`` and ``database_derived``), ``vcf_valid``
    (``tables/variants.vcf`` opens and iterates under pysam)]; gate n_high >= 30, precision >=
    0.95, recall >= 0.90, no row missing a flag, VCF valid; expected n_high in [60, 300] and
    precision >= 0.97.
12. a call is a ``novel_orfs.tsv`` row with ``called = true`` and numeric ``class_fdr <= ALPHA``;
    it is true when (contig, strand, stop) matches a truth ORF, where the stop is the lowest
    coordinate of the stop codon (GIAP ``end - 2`` on ``+``, ``end`` on ``-``; truth ``stop``).
    ``classes[c]`` for c in intergenic, cds_outframe, antisense, nterm_extension: ``n_truth``,
    ``n_truth_with_peptide`` (``n_true_peptides >= 1``), ``n_calls``, ``n_true_calls``,
    ``precision`` (None without calls), ``below_1_over_alpha`` (``n_truth_with_peptide < 100``,
    M7-12 decision 26). value = non-vacuous classes (>= 1 call); gate >= 2 and precision >= 0.95
    in each; in_expected iff every class not below 1/alpha has >= 30 calls.
13. planted islands = truth genes whose ``tu_id`` lies in ``first_tu_id..last_tu_id``, direction
    ``high`` for ``effect_sd > 0``; a planted island is recovered by an ``islands.tsv`` row of the
    same ``direction`` whose ``genes`` cover >= 50 % of it; a spurious call shares no gene with any
    planted island. value = ``n_recovered`` [``n_spurious``, ``autocorrelation``:
    ``{"morans_i", "p_value"}`` of the lag-1 ``morans_i`` row under the ``tu_block`` null,
    ``regions_ok``, ``regions_missing_types``]. Regions: every ``region_type`` of
    ``annotation/regions.tsv`` has a ``regions_stats.tsv`` row, and a type with a region of
    ``n_eligible >= 3`` has a row with a non-NA ``effect``. Gate: lag-1 tu_block Moran's I > 0 with
    p <= 0.01, n_recovered >= 4, n_spurious <= 2, regions_ok; expected n_recovered [5, 6].
14. value = OP pairs in ``pairs.tsv`` [``n_pairs``, ``n_tus`` (``operons.tsv`` rows),
    ``tests_run`` (some ``test != none``), ``too_few_op_pairs`` (the flag on any pair),
    ``regulondb_recovery``: share of multi-gene truth TUs whose gene set is an ``operons.tsv``
    ``genes`` set, or ``"not available"`` unless ``provenance.json`` ``tu_source = regulondb``
    (example 9)]; gate tests_run, no too_few_op_pairs, n_tus >= 1000; expected [1500, None].
15. ``report/report.html``: a section of ``REPORT_SECTIONS`` has content when it holds an
    ``<svg`` or a ``<table data-table=``. value = sections with content [``empty_sections``,
    ``q_floor_caveat`` (``q_floor_above_alpha`` appears in the page)]; gate all five and no caveat;
    expected [5, 5]. Missing report: fail.
16. ``logs/trace.tsv`` tasks as ``scripts/measure_runtimes.py`` reads them (CACHED/FAILED out):
    value = wall hours (manifest ``finished_at - started_at``) [``cpu_h`` = sum of realtime x
    %cpu / 100, ``peak_rss_gb`` (GiB), ``over_budget`` (wall > 8 h, CPU > 150 CPU-h or RSS > 64
    GB), ``rows``: the ``measure_runtimes.summarise`` markdown rows under ``RUNTIME_KEY``]; never
    fails (``gate_ok`` true: an overrun goes to the backlog); expected [3, 8].

``summarise(outdir, inputs, inputs_again=None) -> dict``: ``{"gates": {"1": ..., "16": ...},
"all_gates_ok": bool, "outside_expected": [gate numbers as str]}``.
``python scripts/dev/mvp04_summary.py OUTDIR INPUTS [--inputs-again DIR]`` writes it to
``OUTDIR/evidence_summary.json`` as strict JSON (NaN written as null) and exits 0.
"""

from __future__ import annotations

import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers import mvp04_summary_data as S

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
GATE_NUMBERS = [str(n) for n in range(1, 17)]
ENTRY_KEYS = {"value", "expected", "gate", "gate_ok", "in_expected", "status"}


def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary_full", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    m = importlib.util.module_from_spec(spec)
    sys.modules["mvp04_summary_full"] = m
    spec.loader.exec_module(m)
    return m


def gate(n: int, d: S.Dirs) -> dict[str, Any]:
    m = mod()
    if n == 9:
        res = m.gate9(d.inputs, d.again)
    elif n == 7:
        res = m.gate7_entry(d.out, d.inputs)
    else:
        res = getattr(m, f"gate{n}")(d.out, d.inputs)
    assert isinstance(res, dict), res
    missing = ENTRY_KEYS - set(res)
    assert not missing, f"gate {n} entry lacks {sorted(missing)}: {res}"
    assert res["status"] == mod().status_of(res["gate_ok"], res["in_expected"]), res
    return res


@pytest.fixture
def d(tmp_path: Path) -> S.Dirs:
    return S.build(tmp_path)


def _set_q(out: Path, genes: dict[str, str | None]) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["locus_tag"] in genes:
                q = genes[r["locus_tag"]]
                r["gene_q_value"] = "NA" if q is None else q
        return rows

    S.edit_tsv(out / "tables" / "genes.tsv", fn)


# ---------------------------------------------------------------- status and constants


def test_status_vocabulary() -> None:
    m = mod()
    assert m.status_of(True, True) == "pass"
    assert m.status_of(True, False) == "outside_expected"
    assert m.status_of(False, True) == "fail"
    assert m.status_of(False, False) == "fail"


def test_constants_follow_the_spec() -> None:
    m = mod()
    assert m.ALPHA == 0.01 and m.Q_VALUE_MAX == 0.01
    assert set(m.FORBIDDEN_WARNINGS) == {"decoy_unscored", "decoy_pairing_low",
                                         "q_floor_above_alpha", "no_eligible_genes"}  # fmt: skip
    assert tuple(m.REPORT_SECTIONS) == S.REPORT_SECTIONS
    assert m.RUNTIME_KEY == "mvp04_k12_evidence"


def test_gate7_from_part1_is_still_there() -> None:
    m = mod()
    assert callable(m.gate7) and callable(m.gate7_from_rows)


@pytest.mark.parametrize("n", range(1, 17))
def test_the_all_pass_fixture_passes_every_gate_in_range(n: int, d: S.Dirs) -> None:
    res = gate(n, d)
    assert res["status"] == "pass", res


# ---------------------------------------------------------------- gate 1: run completes


def test_gate1_reads_the_manifest(d: S.Dirs) -> None:
    res = gate(1, d)
    assert res["value"] == 0 and res["run_status"] == "completed" and res["failed_steps"] == []


def test_gate1_a_failed_step_fails(d: S.Dirs) -> None:
    def fn(m: dict[str, Any]) -> None:
        m["steps"][3]["status"] = "failed"

    S.edit_json(d.out / "run_manifest.json", fn)
    res = gate(1, d)
    assert res["status"] == "fail" and res["failed_steps"] == ["context_stats"]


def test_gate1_non_zero_exit_fails(d: S.Dirs) -> None:
    def fn(m: dict[str, Any]) -> None:
        m["exit_code"], m["status"] = 1, "failed"

    S.edit_json(d.out / "run_manifest.json", fn)
    res = gate(1, d)
    assert (res["value"], res["status"]) == (1, "fail")


def test_gate1_missing_manifest_fails_without_raising(d: S.Dirs) -> None:
    (d.out / "run_manifest.json").unlink()
    assert gate(1, d)["status"] == "fail"


# ---------------------------------------------------------------- gate 2: tier and warnings


def test_gate2_tier_and_no_warnings(d: S.Dirs) -> None:
    res = gate(2, d)
    assert res["value"] == "decoy_picked" and res["warnings_found"] == []


def test_gate2_example4_mayu_tier_fails(d: S.Dirs) -> None:
    """Example 4: unscored decoys give ``mayu`` with ``decoy_unscored``; gate 2 fails."""

    def tier(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["detection_tier"] = "mayu"
        return rows

    S.edit_tsv(d.out / "tables" / "genes.tsv", tier)
    S.edit_json(d.out / "run_manifest.json",
                lambda m: m["warnings"].append("gene_table:detection.decoy_unscored"))  # fmt: skip
    res = gate(2, d)
    assert res["value"] == "mayu" and res["warnings_found"] == ["decoy_unscored"]
    assert res["status"] == "fail"


@pytest.mark.parametrize("token", ["decoy_unscored", "decoy_pairing_low", "no_eligible_genes"])
def test_gate2_a_forbidden_token_in_the_events_log_fails(d: S.Dirs, token: str) -> None:
    with (d.out / "logs" / "events.jsonl").open("a") as fh:
        fh.write(json.dumps({"event": "gene_table.detection", "level": "WARNING",
                             "step": "gene_table", "data": {"flag": token}}) + "\n")  # fmt: skip
    res = gate(2, d)
    assert res["warnings_found"] == [token] and res["status"] == "fail"


def test_gate2_a_forbidden_manifest_flag_fails(d: S.Dirs) -> None:
    S.edit_json(d.out / "run_manifest.json", lambda m: m["flags"].append("no_eligible_genes"))
    assert gate(2, d)["warnings_found"] == ["no_eligible_genes"]


def test_gate2_q_floor_above_alpha_is_derived_from_the_tier_row(d: S.Dirs) -> None:
    """Example 5's contrast: T = 85 gives q_floor 0.0118 > 0.01 even if no log line says so."""
    S.set_tier_metric(d.out, "q_floor", 1 / 85)
    res = gate(2, d)
    assert res["warnings_found"] == ["q_floor_above_alpha"] and res["status"] == "fail"


def test_gate2_mixed_tiers(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        rows[100]["detection_tier"] = "two_peptide"
        return rows

    S.edit_tsv(d.out / "tables" / "genes.tsv", fn)
    assert gate(2, d)["value"] == "mixed"


# ---------------------------------------------------------------- gates 3 and 4: f and T


def test_gate3_value_and_bands(d: S.Dirs) -> None:
    assert gate(3, d)["value"] == pytest.approx(1.0)
    S.set_tier_metric(d.out, "decoy_paired_share", 0.9)
    res = gate(3, d)
    assert (res["value"], res["status"]) == (pytest.approx(0.9), "outside_expected")
    S.set_tier_metric(d.out, "decoy_paired_share", 0.79)
    assert gate(3, d)["status"] == "fail"
    S.set_tier_metric(d.out, "decoy_paired_share", "NA")
    assert gate(3, d)["status"] == "fail"


def test_gate4_example5(d: S.Dirs) -> None:
    """Example 5: T = 2,600 -> q_floor 0.000385."""
    res = gate(4, d)
    assert res["value"] == 2600 and res["q_floor"] == pytest.approx(1 / 2600)
    S.set_tier_metric(d.out, "q_floor", 1 / 200)
    res = gate(4, d)
    assert (res["value"], res["status"]) == (200, "outside_expected")
    S.set_tier_metric(d.out, "q_floor", 1 / 110)
    assert gate(4, d)["gate_ok"] is True
    S.set_tier_metric(d.out, "q_floor", 1 / 85)
    assert gate(4, d)["status"] == "fail"
    S.set_tier_metric(d.out, "q_floor", "NA")
    assert gate(4, d)["status"] == "fail"


# ---------------------------------------------------------------- gate 5: decoy wins


def test_gate5_counts_decoy_wins_only(d: S.Dirs) -> None:
    """40 wins; the 5 decoys above q_value_max leave the competition (ADR-0047)."""
    res = gate(5, d)
    assert res["value"] == len(S.DECOY_WINS)
    assert res["n_paired_genes"] == len(S.DECOY_WINS)


def test_gate5_a_decoy_on_a_target_won_gene_is_no_win(d: S.Dirs) -> None:
    """A paired decoy on a gene with a finite gene q (the target won) is paired, not a win."""
    S.append_tsv(
        d.inputs / "peptides.tsv",
        [
            {
                "peptide": "QQQQQQQK",
                "protein_ids": "rev_" + S.locus(100),
                "q_value": 0.003,
                "pep": 0.03,
                "score": 4.7,
                "is_decoy": "true",
                "intensity": "NA",
                "run": "run2",
            }
        ],
    )
    res = gate(5, d)
    assert (res["value"], res["n_paired_genes"]) == (40, 41)


def test_gate5_partial_genes_never_count(d: S.Dirs) -> None:
    S.append_tsv(
        d.inputs / "peptides.tsv",
        [
            {
                "peptide": "QQQQQQQK",
                "protein_ids": "rev_" + S.locus(0),
                "q_value": 0.003,
                "pep": 0.03,
                "score": 4.7,
                "is_decoy": "true",
                "intensity": "NA",
                "run": "run2",
            }
        ],
    )
    assert gate(5, d)["value"] == 40


def test_gate5_d4_no_decoy_wins_fails(d: S.Dirs) -> None:
    """D4: D = 0 here means the decoys are not competing (the 1869554 picture)."""
    S.edit_tsv(d.inputs / "peptides.tsv", lambda rows: [r for r in rows if r["is_decoy"] != "true"])
    res = gate(5, d)
    assert (res["value"], res["status"]) == (0, "fail")


def test_gate5_few_wins_are_outside_expected(d: S.Dirs) -> None:
    keep = {"rev_" + S.locus(i) for i in list(S.DECOY_WINS)[:5]}
    S.edit_tsv(
        d.inputs / "peptides.tsv",
        lambda rows: [r for r in rows if r["is_decoy"] != "true" or r["protein_ids"] in keep],
    )
    res = gate(5, d)
    assert (res["value"], res["status"]) == (5, "outside_expected")


# ---------------------------------------------------------------- gates 6-8: detection


def test_gate6_values(d: S.Dirs) -> None:
    res = gate(6, d)
    assert res["value"] == 2400
    assert res["n_eligible"] == S.N_GENES - S.N_PARTIAL
    assert res["n_truth_detected"] == len(S.TRUE_GENES)
    lo, hi = res["expected"]
    assert lo == pytest.approx(0.85 * 2520) and hi == pytest.approx(2520)


def test_gate6_saturated_fails(d: S.Dirs) -> None:
    """Example 8: n_detected >= 0.95 x eligible is saturated (1865149)."""
    S.set_tier_metric(d.out, "n_detected", 3800)
    assert gate(6, d)["status"] == "fail"


def test_gate6_too_few_fails_and_low_is_outside_expected(d: S.Dirs) -> None:
    S.set_tier_metric(d.out, "n_detected", 999)
    assert gate(6, d)["status"] == "fail"
    S.set_tier_metric(d.out, "n_detected", 2000)  # < 0.85 x 2,520 = 2,142
    assert gate(6, d)["status"] == "outside_expected"


def test_gate7_entry_wraps_part1(d: S.Dirs) -> None:
    res = gate(7, d)
    assert res["value"] == pytest.approx(20 / 2400)
    assert res["n_detected"] == 2400 and res["n_false_detected"] == 20
    assert res["false_detected_by_channel"]["protein_null"] == 20


def test_gate7_bands(d: S.Dirs) -> None:
    """FDP 0.015 is inside the gate (0.02) but outside the expected 0.01; 0.025 fails."""

    def detect(extra: int) -> None:
        idx = set(range(2600, 2600 + extra))

        def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
            for i, r in enumerate(rows):
                if i in idx:
                    r["detected"], r["detection_status"] = "true", "detected"
            return rows

        S.edit_tsv(d.out / "tables" / "genes.tsv", fn)

    detect(17)  # 37 / 2,417 = 0.0153
    assert gate(7, d)["status"] == "outside_expected"
    detect(45)  # 65 / 2,445 = 0.0266
    assert gate(7, d)["status"] == "fail"


def test_gate7_nothing_detected_fails(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["detected"] = "false"
        return rows

    S.edit_tsv(d.out / "tables" / "genes.tsv", fn)
    res = gate(7, d)
    assert res["status"] == "fail"
    assert res["value"] is None or math.isnan(res["value"])


def test_gate8_share(d: S.Dirs) -> None:
    res = gate(8, d)
    assert res["value"] == pytest.approx(2520 / 3990)
    assert (res["n_truth_detected"], res["n_eligible"]) == (2520, 3990)


def test_gate8_a_protein_null_never_makes_a_gene_true(d: S.Dirs) -> None:
    """Decision 11: the false-detected genes' truth rows carry their locus tag (is_wrong true)."""
    assert gate(8, d)["n_truth_detected"] == len(S.TRUE_GENES)


def test_gate8_partial_genes_leave_the_share(d: S.Dirs) -> None:
    """D8: a true peptide on a partial gene counts in neither numerator nor denominator."""
    S.append_tsv(
        d.inputs / "truth" / "peptides.tsv",
        [{"peptide": "WWWWWWWK", "is_wrong": "false", "origin": "cds", "locus_tag": S.locus(0)}],
    )
    res = gate(8, d)
    assert (res["n_truth_detected"], res["n_eligible"]) == (2520, 3990)


def test_gate8_example8_saturated_share_fails(d: S.Dirs) -> None:
    rows = [{"peptide": S.peptide(10_000 + i), "is_wrong": "false", "origin": "cds",
             "locus_tag": S.locus(i)} for i in range(2530, 3800)]  # fmt: skip
    S.append_tsv(d.inputs / "truth" / "peptides.tsv", rows)
    res = gate(8, d)
    assert res["value"] > 0.9 and res["status"] == "fail"


# ---------------------------------------------------------------- gate 9: TE7 on the inputs


def test_gate9_identical_calls_pass(d: S.Dirs) -> None:
    assert gate(9, d)["value"] == []


@pytest.mark.parametrize(
    "rel",
    [
        "peptides.tsv",
        "truth/peptides.tsv",
        "truth/alt_sites.tsv",
        "truth/islands.tsv",
        "truth/genes.tsv",
        "truth/novel_orfs.tsv",
    ],
)
def test_gate9_one_changed_byte_fails(d: S.Dirs, rel: str) -> None:
    p = d.again / rel
    b = bytearray(p.read_bytes())
    b[-2] = ord("X") if b[-2] != ord("X") else ord("Y")
    p.write_bytes(bytes(b))
    res = gate(9, d)
    assert res["value"] == [rel] and res["status"] == "fail"


def test_gate9_an_extra_truth_table_differs(d: S.Dirs) -> None:
    (d.again / "truth" / "extra.tsv").write_text("a\n")
    assert gate(9, d)["value"] == ["truth/extra.tsv"]


def test_gate9_a_missing_required_table_fails_even_in_both(d: S.Dirs) -> None:
    (d.inputs / "truth" / "islands.tsv").unlink()
    (d.again / "truth" / "islands.tsv").unlink()
    assert gate(9, d)["status"] == "fail"


def test_gate9_provenance_and_sample_genome_are_not_compared(d: S.Dirs) -> None:
    """Gate 9 names ``peptides.tsv`` and ``truth/*.tsv`` only (provenance carries a git SHA)."""
    (d.again / "provenance.json").write_text("{}\n")
    (d.again / "truth" / "sample_genome.fna").write_text(">x\nACGT\n")
    assert gate(9, d)["status"] == "pass"


def test_gate9_without_a_second_call_fails() -> None:
    res = mod().gate9(Path("/nonexistent/inputs"), None)
    assert res["status"] == "fail"


# ---------------------------------------------------------------- gate 10: mapping


def test_gate10_shares(d: S.Dirs) -> None:
    res = gate(10, d)
    n = 2520 + 50 + 20 + 30
    assert res["n_peptides"] == n
    sh = res["shares"]
    assert sh["exact"] == pytest.approx(2540 / n) and res["value"] == pytest.approx(2540 / n)
    assert sh["tolerant"] == pytest.approx(50 / n)
    assert sh["unmapped"] == pytest.approx(30 / n)


def test_gate10_a_peptide_with_both_routes_is_exact(d: S.Dirs) -> None:
    tol = next(r for r in S.read_tsv(d.out / "mapping" / "peptide_hits.tsv")
               if r["route"] == "tolerant")  # fmt: skip
    S.append_tsv(d.out / "mapping" / "peptide_hits.tsv", [{**tol, "route": "exact", "n_diff": 0}])
    sh = gate(10, d)["shares"]
    assert sh["tolerant"] == pytest.approx(49 / 2620)


def test_gate10_no_tolerant_hit_fails(d: S.Dirs) -> None:
    S.edit_tsv(d.out / "mapping" / "peptide_hits.tsv",
               lambda rows: [r for r in rows if r["route"] != "tolerant"])  # fmt: skip
    assert gate(10, d)["status"] == "fail"


def test_gate10_low_exact_share_is_only_reported(d: S.Dirs) -> None:
    rows = [{"sequence_id": S.peptide(20_000 + k), "sequence": S.peptide(20_000 + k),
             "kind": "peptide", "label": "unmapped", "reason": "no_hit", "contaminant_match": "NA"}
            for k in range(300)]  # fmt: skip
    S.append_tsv(d.out / "mapping" / "unmapped.tsv", rows)
    assert gate(10, d)["status"] == "outside_expected"


# ---------------------------------------------------------------- gate 11: variants


def test_gate11_values(d: S.Dirs) -> None:
    res = gate(11, d)
    assert res["value"] == 60 and res["precision"] == 1.0 and res["recall"] == 1.0
    assert res["n_recall_sites"] == S.N_ALT_MULTI
    assert res["n_rows_missing_flags"] == 0 and res["vcf_valid"] is True


def test_gate11_a_wrong_alt_residue_is_false(d: S.Dirs) -> None:
    """Same codon, other ALT residue: not a match (precision 57/60 = 0.95, outside 0.97)."""

    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows[:3]:
            r["alt_aa"] = "L"
        return rows

    S.edit_tsv(d.out / "tables" / "variants.tsv", fn)
    res = gate(11, d)
    assert res["precision"] == pytest.approx(57 / 60)
    assert res["recall"] == pytest.approx(47 / 50)
    assert res["status"] == "outside_expected"


def test_gate11_recall_counts_only_sites_with_two_true_peptides(d: S.Dirs) -> None:
    """Dropping the calls of the 10 single-peptide sites leaves recall at 1."""
    single = {str(S.alt_codon(k)) for k in range(S.N_ALT_MULTI, S.N_ALT)}
    S.edit_tsv(d.out / "tables" / "variants.tsv",
               lambda rows: [r for r in rows if r["codon_start"] not in single])  # fmt: skip
    res = gate(11, d)
    assert res["recall"] == 1.0 and res["value"] == 50


def test_gate11_low_recall_fails(d: S.Dirs) -> None:
    drop = {str(S.alt_codon(k)) for k in range(6)}  # 44 / 50 = 0.88 < 0.90
    S.edit_tsv(d.out / "tables" / "variants.tsv",
               lambda rows: [r for r in rows if r["codon_start"] not in drop])  # fmt: skip
    res = gate(11, d)
    assert res["recall"] == pytest.approx(44 / 50) and res["status"] == "fail"


def test_gate11_candidates_are_not_high(d: S.Dirs) -> None:
    """The 5 candidate rows would be false; they never enter precision."""
    assert gate(11, d)["precision"] == 1.0


def test_gate11_a_row_without_strain_difference_fails(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        rows[-1]["flags"] = "snv;database_derived"
        return rows

    S.edit_tsv(d.out / "tables" / "variants.tsv", fn)
    res = gate(11, d)
    assert res["n_rows_missing_flags"] == 1 and res["status"] == "fail"


def test_gate11_an_invalid_vcf_fails(d: S.Dirs) -> None:
    (d.out / "tables" / "variants.vcf").write_text("not a vcf\n")
    res = gate(11, d)
    assert res["vcf_valid"] is False and res["status"] == "fail"


def test_gate11_too_few_high_fails(d: S.Dirs) -> None:
    S.edit_tsv(d.out / "tables" / "variants.tsv", lambda rows: rows[:29])
    assert gate(11, d)["status"] == "fail"


# ---------------------------------------------------------------- gate 12: novel ORFs


def test_gate12_classes(d: S.Dirs) -> None:
    res = gate(12, d)
    assert res["value"] == 4
    for c in S.NOVEL_CLASSES:
        k = res["classes"][c]
        assert (k["n_truth"], k["n_truth_with_peptide"]) == (120, 110)
        assert (k["n_calls"], k["n_true_calls"], k["precision"]) == (40, 40, 1.0)
        assert k["below_1_over_alpha"] is False


def test_gate12_stop_convention_on_both_strands(d: S.Dirs) -> None:
    """A call whose ``end`` is the truth ``end`` on the minus strand (the highest coordinate,
    not the stop) does not match: GIAP's minus-strand ``end`` is the lowest coordinate."""
    truth = [r for r in S.read_tsv(d.inputs / "truth" / "novel_orfs.tsv")
             if r["class"] == "antisense"]  # fmt: skip
    minus = [r for r in truth if r["strand"] == "-"][:4]
    plus = [r for r in truth if r["strand"] == "+"][:4]
    wrong = [{**r, "start": int(r["start"]), "end": int(r["end"])} for r in minus]
    shifted = [{**r, "start": int(r["start"]), "end": int(r["end"]) + 3} for r in plus]

    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        ids = {f"GIAPNOVEL_{r['orf_id']}" for r in minus + plus}
        keep = [r for r in rows if r["orf_id"] not in ids]
        bad = [dict(zip(S.NOVEL_HEADER, map(str, S.novel_call(t)), strict=True)) for t in shifted]
        for t in wrong:  # end = highest coordinate on '-': the minus convention broken
            row = dict(zip(S.NOVEL_HEADER, map(str, S.novel_call(t)), strict=True))
            row["end"] = str(t["end"])
            bad.append(row)
        return keep + bad

    S.edit_tsv(d.out / "tables" / "novel_orfs.tsv", fn)
    k = gate(12, d)["classes"]["antisense"]
    assert (k["n_calls"], k["n_true_calls"]) == (40, 32)


def test_gate12_uncalled_or_unestimated_rows_are_no_calls(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["class"] == "intergenic":
                r["called"] = "false"
            if r["class"] == "cds_outframe":
                r["class_fdr"] = "unestimated"
            if r["class"] == "antisense":
                r["class_fdr"] = "0.02"
        return rows

    S.edit_tsv(d.out / "tables" / "novel_orfs.tsv", fn)
    res = gate(12, d)
    assert res["value"] == 1
    assert res["classes"]["intergenic"]["n_calls"] == 0
    assert res["classes"]["intergenic"]["precision"] is None
    assert res["status"] == "fail"


def test_gate12_d3_a_class_below_100_is_reported_not_failed(d: S.Dirs) -> None:
    """D3 / M7-12 decision 26: < 100 truth ORFs with a true peptide is vacuous by construction,
    so its call count does not enter the expected check."""

    def truth(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["class"] == "intergenic":
                r["n_true_peptides"] = "0"
        return rows

    S.edit_tsv(d.inputs / "truth" / "novel_orfs.tsv", truth)
    S.edit_tsv(d.out / "tables" / "novel_orfs.tsv",
               lambda rows: [r for r in rows if r["class"] != "intergenic"])  # fmt: skip
    res = gate(12, d)
    assert res["classes"]["intergenic"]["below_1_over_alpha"] is True
    assert res["value"] == 3 and res["status"] == "pass"


def test_gate12_few_calls_in_a_supplied_class_is_outside_expected(d: S.Dirs) -> None:
    drop = {f"GIAPNOVEL_NOVEL_antisense_{k:04d}" for k in range(20)}
    S.edit_tsv(d.out / "tables" / "novel_orfs.tsv",
               lambda rows: [r for r in rows if r["orf_id"] not in drop])  # fmt: skip
    res = gate(12, d)
    assert res["classes"]["antisense"]["n_calls"] == 20
    assert res["status"] == "outside_expected"


def test_gate12_low_precision_fails(d: S.Dirs) -> None:
    false_calls = [{"orf_id": f"FALSE_{k}", "class": "antisense", "contig": S.CONTIG,
                    "strand": "+", "start_longest": 9_000_001 + 600 * k,
                    "end": 9_000_300 + 600 * k, "class_fdr": "0.005", "called": "true",
                    "flags": "NA"} for k in range(3)]  # fmt: skip
    S.append_tsv(d.out / "tables" / "novel_orfs.tsv", false_calls)
    res = gate(12, d)
    assert res["classes"]["antisense"]["precision"] == pytest.approx(40 / 43)
    assert res["status"] == "fail"


# ---------------------------------------------------------------- gate 13: context


def test_gate13_values(d: S.Dirs) -> None:
    res = gate(13, d)
    assert res["value"] == 6 and res["n_spurious"] == 0
    assert res["autocorrelation"]["morans_i"] == pytest.approx(0.22)
    assert res["autocorrelation"]["p_value"] == pytest.approx(0.001)
    assert res["regions_ok"] is True and res["regions_missing_types"] == []


def test_gate13_example10_half_overlap_recovers_a_down_island(d: S.Dirs) -> None:
    genes = S.island_genes(4)

    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        rows[4]["genes"] = ";".join(genes[:8])  # 8 / 16 = 50 %
        rows[5]["genes"] = ";".join(S.island_genes(5)[:7])  # 7 / 16 < 50 %
        return rows

    S.edit_tsv(d.out / "tables" / "islands.tsv", fn)
    res = gate(13, d)
    assert res["value"] == 5 and res["status"] == "pass"


def test_gate13_wrong_direction_is_no_recovery(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows[3:]:
            r["direction"] = "high"
        return rows

    S.edit_tsv(d.out / "tables" / "islands.tsv", fn)
    res = gate(13, d)
    assert res["value"] == 3 and res["n_spurious"] == 0 and res["status"] == "fail"


def test_gate13_four_recovered_is_outside_expected(d: S.Dirs) -> None:
    S.edit_tsv(d.out / "tables" / "islands.tsv", lambda rows: rows[:4])
    res = gate(13, d)
    assert (res["value"], res["status"]) == (4, "outside_expected")


def test_gate13_spurious_calls(d: S.Dirs) -> None:
    extra = [dict(zip(S.ISLAND_HEADER, map(str, S.island_call(0, [S.locus(3500 + 2 * k)])),
                      strict=True)) for k in range(3)]  # fmt: skip
    S.append_tsv(d.out / "tables" / "islands.tsv", extra)
    res = gate(13, d)
    assert res["n_spurious"] == 3 and res["status"] == "fail"


def test_gate13_uses_the_tu_block_null(d: S.Dirs) -> None:
    """Genes of one TU share their abundance, so the gene-label null is anti-conservative: the gate
    reads the lag-1 ``tu_block`` row."""

    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["null_model"] == "tu_block" and r["lag"] == "1":
                r["p_value"] = "0.2"
        return rows

    S.edit_tsv(d.out / "tables" / "autocorrelation.tsv", fn)
    res = gate(13, d)
    assert res["autocorrelation"]["p_value"] == pytest.approx(0.2) and res["status"] == "fail"


def test_gate13_negative_autocorrelation_fails(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["value"] = "-0.1"
        return rows

    S.edit_tsv(d.out / "tables" / "autocorrelation.tsv", fn)
    assert gate(13, d)["status"] == "fail"


def test_gate13_a_called_region_type_without_rows_fails(d: S.Dirs) -> None:
    S.append_tsv(
        d.out / "annotation" / "regions.tsv",
        [
            {
                "region_id": "bgc_1",
                "region_type": "bgc",
                "contig": S.CONTIG,
                "start": 1,
                "end": 9,
                "source_tool": "antiSMASH",
                "score": "NA",
                "genes": S.locus(3900),
                "notes": "",
            }
        ],
    )
    res = gate(13, d)
    assert res["regions_missing_types"] == ["bgc"] and res["status"] == "fail"


def test_gate13_a_large_region_type_with_only_na_effects_fails(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            if r["region_type"] == "prophage":
                r["effect"] = "NA"
        return rows

    S.edit_tsv(d.out / "tables" / "regions_stats.tsv", fn)
    assert gate(13, d)["regions_ok"] is False


# ---------------------------------------------------------------- gate 14: operons


def test_gate14_values(d: S.Dirs) -> None:
    res = gate(14, d)
    assert (res["value"], res["n_pairs"], res["n_tus"]) == (2000, 3999, 2000)
    assert res["tests_run"] is True and res["too_few_op_pairs"] is False
    assert res["regulondb_recovery"] == pytest.approx(1.0)


def test_gate14_example9_gap_rule_recovery_not_available(d: S.Dirs) -> None:
    S.edit_json(d.inputs / "provenance.json", lambda p: p.update(tu_source="gap_rule"))
    res = gate(14, d)
    assert res["regulondb_recovery"] == "not available" and res["status"] == "pass"


def test_gate14_recovery_is_reported_not_gated(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows[:1000]:
            r["genes"] = r["genes"].split(";")[0]
        return rows

    S.edit_tsv(d.out / "tables" / "operons.tsv", fn)
    res = gate(14, d)
    assert res["regulondb_recovery"] == pytest.approx(0.5) and res["gate_ok"] is True


def test_gate14_all_tests_none_or_too_few_op_pairs_fails(d: S.Dirs) -> None:
    def none(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["test"] = "none"
        return rows

    S.edit_tsv(d.out / "tables" / "pairs.tsv", none)
    assert gate(14, d)["status"] == "fail"


def test_gate14_too_few_op_pairs_flag_fails(d: S.Dirs) -> None:
    def fn(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        rows[0]["flags"] = "too_few_op_pairs;unestimated"
        return rows

    S.edit_tsv(d.out / "tables" / "pairs.tsv", fn)
    res = gate(14, d)
    assert res["too_few_op_pairs"] is True and res["status"] == "fail"


def test_gate14_too_few_tus_fails(d: S.Dirs) -> None:
    S.edit_tsv(d.out / "tables" / "operons.tsv", lambda rows: rows[:999])
    assert gate(14, d)["status"] == "fail"


# ---------------------------------------------------------------- gate 15: report


def test_gate15_values(d: S.Dirs) -> None:
    res = gate(15, d)
    assert res["value"] == 5 and res["empty_sections"] == [] and res["q_floor_caveat"] is False


def test_gate15_a_placeholder_section_fails(d: S.Dirs) -> None:
    content = {s: S.FIGURE for s in S.REPORT_SECTIONS} | {"variants": S.PLACEHOLDER}
    (d.out / "report" / "report.html").write_text(S.report_html(content))
    res = gate(15, d)
    assert res["empty_sections"] == ["variants"] and res["status"] == "fail"


def test_gate15_content_outside_the_sections_does_not_count(d: S.Dirs) -> None:
    """A figure in the summary section does not fill the operons card."""
    content = {s: S.TABLE for s in S.REPORT_SECTIONS if s != "operons"} | {"summary": S.FIGURE}
    (d.out / "report" / "report.html").write_text(S.report_html(content))
    assert gate(15, d)["empty_sections"] == ["operons"]


def test_gate15_q_floor_caveat_fails(d: S.Dirs) -> None:
    page = (d.out / "report" / "report.html").read_text()
    page = page.replace(
        "<h1>GIAP report</h1>", '<h1>GIAP report</h1>\n<p class="note">q_floor_above_alpha: ...</p>'
    )
    (d.out / "report" / "report.html").write_text(page)
    res = gate(15, d)
    assert res["q_floor_caveat"] is True and res["status"] == "fail"


def test_gate15_missing_report_fails(d: S.Dirs) -> None:
    shutil.rmtree(d.out / "report")
    assert gate(15, d)["status"] == "fail"


# ---------------------------------------------------------------- gate 16: runtime


def test_gate16_values(d: S.Dirs) -> None:
    """Cached rows are left out (measure_runtimes); CPU = 1h x 2 + 2 x 0.5h x 4 + 2h x 1 = 8."""
    res = gate(16, d)
    assert res["value"] == pytest.approx(5.0)
    assert res["cpu_h"] == pytest.approx(8.0)
    assert res["peak_rss_gb"] == pytest.approx(12.0)
    assert res["over_budget"] is False
    assert res["rows"] and all("mvp04_k12_evidence" in r for r in res["rows"][2:])


def test_gate16_overrun_is_recorded_never_failed(d: S.Dirs) -> None:
    """D10: an overrun goes to the backlog."""
    S.edit_json(d.out / "run_manifest.json",
                lambda m: m.update(finished_at="2026-10-24T10:30:00+00:00"))  # fmt: skip
    res = gate(16, d)
    assert res["value"] == pytest.approx(10.5)
    assert res["over_budget"] is True and res["gate_ok"] is True
    assert res["status"] == "outside_expected"


def test_gate16_rss_over_64_gb_is_over_budget(d: S.Dirs) -> None:
    S.append_tsv(d.out / "logs" / "trace.tsv", [dict(zip(S.TRACE_HEADER, map(str, S.trace_row(
        6, "ANNOTATE:BAKTA (k12)", "10m", "100.0%", "70 GB")), strict=True))])  # fmt: skip
    res = gate(16, d)
    assert res["over_budget"] is True and res["status"] != "fail"


# ---------------------------------------------------------------- summary and CLI


def test_summarise_has_every_gate(d: S.Dirs) -> None:
    res = mod().summarise(d.out, d.inputs, d.again)
    assert sorted(res["gates"], key=int) == GATE_NUMBERS
    for n, e in res["gates"].items():
        assert set(e) >= ENTRY_KEYS, (n, e)
    assert res["all_gates_ok"] is True and res["outside_expected"] == []
    assert {"n_detected", "n_false_detected", "realised_gene_fdp",
            "false_detected_by_channel"} <= set(res["gates"]["7"])  # fmt: skip


def test_summarise_without_the_second_call_fails_gate9_only(d: S.Dirs) -> None:
    res = mod().summarise(d.out, d.inputs)
    failed = [n for n, e in res["gates"].items() if e["status"] == "fail"]
    assert failed == ["9"] and res["all_gates_ok"] is False


def test_summarise_lists_outside_expected_gates(d: S.Dirs) -> None:
    S.set_tier_metric(d.out, "decoy_paired_share", 0.9)
    res = mod().summarise(d.out, d.inputs, d.again)
    assert res["outside_expected"] == ["3"] and res["all_gates_ok"] is True


def _strict(text: str) -> Any:
    def bad(c: str) -> None:
        raise ValueError(f"non-strict JSON constant {c}")

    return json.loads(text, parse_constant=bad)


def test_cli_writes_strict_json(d: S.Dirs) -> None:
    def none(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        for r in rows:
            r["detected"] = "false"
        return rows

    S.edit_tsv(d.out / "tables" / "genes.tsv", none)  # FDP NaN
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(d.out), str(d.inputs), "--inputs-again", str(d.again)],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    res = _strict((d.out / "evidence_summary.json").read_text())
    assert sorted(res["gates"], key=int) == GATE_NUMBERS
    assert res["gates"]["7"]["value"] is None and res["gates"]["7"]["status"] == "fail"


def test_cli_writes_nothing_else(d: S.Dirs) -> None:
    before = sorted(p for p in d.out.rglob("*"))
    ins = sorted(p for p in d.inputs.rglob("*"))
    subprocess.run([sys.executable, str(SCRIPT), str(d.out), str(d.inputs)], check=True,
                   capture_output=True)  # fmt: skip
    after = sorted(p for p in d.out.rglob("*"))
    assert sorted(set(after) - set(before)) == [d.out / "evidence_summary.json"]
    assert sorted(p for p in d.inputs.rglob("*")) == ins
