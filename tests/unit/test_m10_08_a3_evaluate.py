"""``a3_feasibility.py evaluate`` and ``early-fail``: the L6 verdict on a ``.ssl`` (M10-08b).

Spec: docs/specs/M10-08-a3-feasibility.md, "The failure rule L6", "Computing (a)" steps 1-6,
"Computing (b)" steps 7-14, the verdict JSON fields, exit codes, examples 1-10 and 12, edge cases
(empty input, ties, nonstandard letters, duplicated loci, circular vs linear, minus strand) and
acceptance criteria 1, 2, 3 and 5.

Every ``.ssl``, genome, Bakta JSON and contaminant FASTA is synthetic and written under
``tmp_path`` (``helpers.a3_eval_kit``); nothing reads Nissle, the network, a GPU or
``$GIPPA_ROOT``. The direct-``match_exact`` oracle reuses ``gippa.mapping`` exactly as the spec's
steps 9-11 prescribe. The one statistical claim (the Wilson 95 % interval of ``share``, R4) is
checked here for equality with ``gippa.stats.cutoffs.wilson_ci`` and calibrated through the
evaluator in ``tests/calibration/test_m10_08_share_ci_calibration.py``.
"""

from __future__ import annotations

import ast
import json
import random
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_eval_kit import (
    COUNT_FIELDS,
    JSON_FIELDS,
    REASONS,
    RESIDUES,
    RULE,
    SSL_COLUMNS,
    Case,
    Contig,
    Genome,
    Row,
    call,
    coding,
    fold,
    has_subcommand,
    ladder,
    make_case,
    planted_contig,
    random_dna,
    random_peptides,
    ranked_case,
    revcomp,
    sequence_like_strings,
)
from helpers.a3_kit import A3_FEASIBILITY, is_canonical_json, require, sha256_hex
from helpers.slurm_sandbox import REPO
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.stats.cutoffs import wilson_ci


def _ok(res: Any) -> str:
    return f"rc={res.rc}\nstdout:\n{res.stdout}\nstderr:\n{res.stderr}"


def _j(res: Any) -> dict[str, Any]:
    assert res.j is not None, "evaluate wrote no J\n" + _ok(res)
    j: dict[str, Any] = res.j
    return j


def _close(x: float | None, want: float, tol: float = 5e-5) -> bool:
    return x is not None and abs(x - want) <= tol


def _oracle(genome: Genome, top: list[str]) -> dict[str, Any]:
    """Steps 9-12 by direct calls: target matches, decoy matches, matches per replicon."""
    from gippa.config import resolve_config
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.mapping.sixframe import build_sixframe
    from gippa.rng import for_step

    cfg = resolve_config()
    seqs, reps = genome.as_inputs()
    index = build_sixframe(seqs, reps)
    decoys = build_decoys(index, "shuffle", cfg.mapping.decoys.shuffle_window_aa,
                          for_step("a3_feasibility", cfg.seed))  # fmt: skip
    hits = match_exact(top, index, {"shuffle": decoys}, dict.fromkeys(top, "shuffle"),
                       min_len=7, max_len=30)  # fmt: skip
    tgt = hits[hits["search"] == "target"]
    by_rep = {c: int(g["peptide"].nunique()) for c, g in tgt.groupby("contig")}
    return {
        "target": set(tgt["peptide"]),
        "decoy": set(hits.loc[hits["search"] == "decoy", "peptide"]),
        "by_rep": by_rep,
    }


def _top(case: Case, n_pass: int) -> list[str]:
    return case.peptides[: min(1000, n_pass)]


# ================================================================== examples 1-4: the verdict


def test_example1_pass(tmp_path: Path) -> None:
    # 640 of the top 1,000 planted, plus 100 planted peptides ranked after 1,000 (not counted)
    # and 500 extra rows repeating peptides on other charges (counted once).
    case = ranked_case(tmp_path, 2400, [r for r in range(1000) if r % 25 < 16], seed=1,
                       planted_below_top=100, duplicates=500)  # fmt: skip
    res = case.evaluate()
    assert res.rc == 0, _ok(res)
    j = _j(res)
    assert j["n_rows"] == 2900
    assert j["n_pass"] == 2400 and j["n_top"] == 1000 and j["n_match"] == 640
    assert _close(j["share"], 0.64)
    assert _close(j["share_ci95"][0], 0.6098) and _close(j["share_ci95"][1], 0.6692)
    assert j["verdict"] == "pass" and j["reason"] is None
    assert j["fails"] == {"too_few_peptides": False, "low_genome_match": False}
    oracle = _oracle(_genome_of(case), _top(case, 2400))
    assert j["n_match"] == len(oracle["target"])
    assert _close(j["decoy_share"], len(oracle["decoy"]) / 1000, 1e-9)


def _genome_of(case: Case) -> Genome:
    """Re-read the case's genome files into a :class:`Genome` (for the oracle)."""
    from gippa.io.fasta import read_genome_fasta

    meta = {r["id"]: r for r in json.loads(case.bakta.read_text())["sequences"]}
    return Genome([Contig(r.id, r.sequence, meta[r.id]["topology"], meta[r.id]["type"])
                   for r in read_genome_fasta(case.genome)])  # fmt: skip


def test_example2_too_few_peptides(tmp_path: Path) -> None:
    case = ranked_case(tmp_path, 420, range(300), seed=2)
    res = case.evaluate()
    assert res.rc == 3, _ok(res)
    j = _j(res)
    assert j["n_pass"] == 420 and j["n_top"] == 420 and j["n_match"] == 300
    assert _close(j["share"], 0.7143)
    assert j["fails"] == {"too_few_peptides": True, "low_genome_match": False}
    assert j["verdict"] == "fail" and j["reason"] == "too_few_peptides"


def test_example3_low_genome_match(tmp_path: Path) -> None:
    case = ranked_case(tmp_path, 3000, range(0, 900, 6), seed=3)  # 150 ranks spread over the top
    res = case.evaluate()
    assert res.rc == 3, _ok(res)
    j = _j(res)
    assert j["n_pass"] == 3000 and j["n_top"] == 1000 and j["n_match"] == 150
    assert _close(j["share"], 0.15)
    assert _close(j["share_ci95"][0], 0.1292) and _close(j["share_ci95"][1], 0.1735)
    assert j["fails"] == {"too_few_peptides": False, "low_genome_match": True}
    assert j["verdict"] == "fail" and j["reason"] == "low_genome_match"


@pytest.mark.parametrize(
    ("n_pass", "n_match", "rc", "verdict", "reason", "fails"),
    [
        (500, 100, 0, "pass", None, (False, False)),
        (500, 99, 3, "fail", "low_genome_match", (False, True)),
        (499, 300, 3, "fail", "too_few_peptides", (True, False)),
        (400, 10, 3, "fail", "both", (True, True)),
    ],
    ids=["500-100-pass", "500-99-low-match", "499-too-few", "both"],
)
def test_example4_boundaries(
    tmp_path: Path, n_pass: int, n_match: int, rc: int, verdict: str, reason: str | None,
    fails: tuple[bool, bool],
) -> None:  # fmt: skip
    # Both rules say "fewer than": exactly 500 peptides and exactly 20 % pass.
    case = ranked_case(tmp_path, n_pass, range(n_match), seed=4 + n_pass + n_match)
    res = case.evaluate()
    assert res.rc == rc, _ok(res)
    j = _j(res)
    assert j["n_pass"] == n_pass and j["n_top"] == n_pass and j["n_match"] == n_match
    assert j["verdict"] == verdict and j["reason"] == reason
    assert j["fails"] == {"too_few_peptides": fails[0], "low_genome_match": fails[1]}
    assert j["rule"] == RULE


# ================================================================== example 5: stripping, cut-off


def _empty_genome(rng: random.Random) -> Genome:
    return Genome([Contig("chr1", random_dna(600, rng))])


def test_example5_strip_fold_distinct_strict_cutoff(tmp_path: Path) -> None:
    rows = [Row("PEPTLDEK", 0.95), Row("PEPTIDEK", 0.91), Row("M[+15.99]PEPTLDEK", 0.85),
            Row("[+42.01]-SAMPLER", 0.90), Row("PEPTLDEK", 0.80)]  # fmt: skip
    case = make_case(tmp_path, rows, _empty_genome(random.Random(5)))
    res = case.evaluate()
    j = _j(res)
    assert j["n_rows"] == 5 and j["n_len_7_30"] == 5
    assert j["n_below_cutoff"] == 1, "a score exactly at 0.8 is excluded (strict >)"
    assert j["n_pass"] == 3 and j["n_top"] == 3  # PEPTLDEK, MPEPTLDEK, SAMPLER
    assert j["length_hist"] == {"7": 1, "8": 1, "9": 1}
    assert j["cutoff"] == 0.8 and j["cutoff_source"] == "cli_default"
    assert res.rc == 3, _ok(res)  # 3 < 500


def test_stripping_rule_variants(tmp_path: Path) -> None:
    # (...) groups, unbracketed signed masses and '-' go; I folds to L; all fold to one peptide.
    # A named [...] group (no mass shift) is an error since review S3, so it is not a form here.
    forms = ["PEPTIDEKR", "PEPT(ox)LDEKR", "PEPM+15.995TLDEKR", "[+42.011]-PEPTLDEKR",
             "PEPTLDEK-17.03R", "PEPTL[-17.027]DEKR"]  # fmt: skip
    rows = [Row(f, 0.95 - 0.01 * i) for i, f in enumerate(forms)]
    case = make_case(tmp_path, rows, _empty_genome(random.Random(6)))
    j = _j(case.evaluate())
    assert {fold(f) for f in forms} == {"PEPTLDEKR", "PEPMTLDEKR"}
    assert j["n_pass"] == 2 and j["n_rows"] == 6


def test_length_filter_bounds(tmp_path: Path) -> None:
    rng = random.Random(7)
    seqs = {n: "".join(rng.choice(RESIDUES) for _ in range(n)) for n in (6, 7, 30, 31)}
    rows = [Row(s, 0.95) for s in seqs.values()] + [Row("[+42.01]-" + seqs[7][:6], 0.95)]
    case = make_case(tmp_path, rows, _empty_genome(rng))
    j = _j(case.evaluate())
    assert j["n_rows"] == 5 and j["n_len_7_30"] == 2 and j["n_pass"] == 2
    assert j["length_hist"] == {"7": 1, "30": 1}


def test_cutoff_override(tmp_path: Path) -> None:
    rows = [Row("PEPTLDEK", 0.95), Row("SAMPLERK", 0.90), Row("ANOTHERK", 0.85)]
    case = make_case(tmp_path, rows, _empty_genome(random.Random(8)))
    j = _j(case.evaluate("--cutoff", "0.9"))
    assert j["cutoff"] == 0.9 and j["cutoff_source"] == "override"
    assert j["n_pass"] == 1 and j["n_below_cutoff"] == 2


def test_columns_are_located_by_name(tmp_path: Path) -> None:
    rows = [Row("PEPTLDEK", 0.95), Row("SAMPLERK", 0.85)]
    cols = ("score", "extra", *[c for c in SSL_COLUMNS if c != "score"][::-1])
    a = make_case(tmp_path / "a", rows, _empty_genome(random.Random(9)))
    b = make_case(tmp_path / "b", rows, _empty_genome(random.Random(9)), columns=cols)
    ja, jb = _j(a.evaluate()), _j(b.evaluate())
    for k in COUNT_FIELDS:
        assert ja[k] == jb[k], k


def test_max_score_of_an_il_group_sets_its_rank(tmp_path: Path) -> None:
    # 1,001 distinct peptides; X's own rows rank last (0.8001) but its I-form scores 0.995, so the
    # folded entry keeps the maximum and enters the top 1,000. Only X is in the genome.
    rng = random.Random(10)
    x = "PEPTLDEKLLR"
    others = random_peptides(1000, rng, avoid=[x])
    rows = [Row(p, s) for p, s in zip(others, ladder(1000, top=0.99), strict=True)]
    rows += [Row(x, 0.8001), Row("PEPTIDEKLIR", 0.995)]
    case = make_case(tmp_path, rows, Genome([planted_contig("chr1", [x], rng)]))
    j = _j(case.evaluate())
    assert j["n_pass"] == 1001 and j["n_top"] == 1000
    assert j["n_match"] == 1, "the I/L group kept a lower score than its maximum"


# ================================================================== example 6: ties, determinism


def _example6(d: Path, *, shuffle: int | None = None) -> Case:
    rng = random.Random(11)
    tie = ["AAAKLLLR", "AAALLLLR", "ACDEFGHK", "ADDEFGHK"]
    as_written = {"AAAKLLLR": "AAAKLILR", "AAALLLLR": "AAALLLLR",
                  "ACDEFGHK": "ACDEFGHK", "ADDEFGHK": "ADDEFGHK"}  # fmt: skip
    rest = random_peptides(1196, rng, avoid=tie)
    scores = ladder(1196, top=0.99)
    high, low = rest[:998], rest[998:]
    rows = [Row(p, s) for p, s in zip(high, scores[:998], strict=True)]
    rows += [Row(as_written[t], 0.95) for t in reversed(tie)]  # ranks 999-1,002 (1-based)
    rows += [Row(p, min(s, 0.9) - 0.01) for p, s in zip(low, scores[998:], strict=True)]
    # Of the top 1,000: 500 planted among the first 998 plus the two tie winners.
    plant = [*high[:500], "AAAKLLLR", "AAALLLLR"]
    genome = Genome([planted_contig("chr1", plant, random.Random(12))])
    # rank order check of the construction: every 'high' score is above 0.95
    assert min(scores[:998]) > 0.95 > max(min(s, 0.9) - 0.01 for s in scores[998:])
    if shuffle is not None:
        random.Random(shuffle).shuffle(rows)
    return make_case(d, rows, genome)


def test_example6_ties_broken_by_folded_sequence(tmp_path: Path) -> None:
    j = _j(_example6(tmp_path).evaluate())
    assert j["n_pass"] == 1200 and j["n_top"] == 1000
    assert j["n_match"] == 502, "ranks 999 and 1,000 must be AAAKLLLR and AAALLLLR"


def test_example6_json_is_byte_identical_on_rerun_and_shuffled_rows(tmp_path: Path) -> None:
    a = _example6(tmp_path / "a").evaluate()
    a2 = _example6(tmp_path / "a").evaluate()
    b = _example6(tmp_path / "b", shuffle=13).evaluate()
    assert a.raw is not None and a2.raw is not None and b.raw is not None, _ok(a)
    assert a.raw == a2.raw
    # The genome files are byte-identical too, so genome_sha256 agrees; only the row order differs.
    assert a.raw == b.raw


# ================================================================== JSON layout (criterion 3)


def test_json_has_every_field_sorted_keys_trailing_newline_no_sequence(tmp_path: Path) -> None:
    case = ranked_case(tmp_path, 600, range(0, 600, 3), seed=14)
    res = case.evaluate()
    j = _j(res)
    assert res.raw is not None and is_canonical_json(res.raw), "sorted keys, indent 2, newline"
    assert set(j) >= JSON_FIELDS, sorted(JSON_FIELDS - set(j))
    assert j["schema_version"] == 1 and j["run"] == "WT_1" and j["cascadia_version"] == "0.0.7"
    assert j["genome_sha256"] == sha256_hex(case.genome.read_bytes())
    assert j["evidence"] is None, "evidence is early-fail only"
    assert j["verdict"] in ("pass", "fail")
    assert not sequence_like_strings(j), sequence_like_strings(j)
    assert not any(p in res.raw for p in case.peptides), "a peptide sequence leaked into J"
    lo, hi = j["share_ci95"]
    want = wilson_ci(j["n_match"], j["n_top"])
    assert _close(lo, want[0]) and _close(hi, want[1])
    assert lo <= j["share"] <= hi
    assert all(isinstance(v, int) for v in j["length_hist"].values())
    assert sum(j["length_hist"].values()) == j["n_pass"]
    assert all(7 <= int(k) <= 30 for k in j["length_hist"])


# ================================================================== example 10, errors, exit codes


def test_example10_header_only_is_a_fail_verdict(tmp_path: Path) -> None:
    case = make_case(tmp_path, [], _empty_genome(random.Random(15)))
    res = case.evaluate()
    assert res.rc == 3, _ok(res)
    j = _j(res)
    assert (j["n_rows"], j["n_pass"], j["n_top"]) == (0, 0, 0)
    assert j["share"] is None
    assert j["fails"] == {"too_few_peptides": True, "low_genome_match": None}
    assert j["verdict"] == "fail" and j["reason"] == "too_few_peptides"


def test_missing_or_zero_byte_ssl_is_a_precondition(tmp_path: Path) -> None:
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(16)))
    case.ssl.write_bytes(b"")
    res = case.evaluate()
    assert res.rc == 2 and res.j is None, _ok(res)
    case.ssl.unlink()
    res = case.evaluate()
    assert res.rc == 2 and res.j is None, _ok(res)


@pytest.mark.parametrize("drop", ["sequence", "score"])
def test_missing_column_is_an_error(tmp_path: Path, drop: str) -> None:
    cols = [c for c in SSL_COLUMNS if c != drop]
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(17)),
                     columns=cols)  # fmt: skip
    res = case.evaluate()
    assert res.rc == 1 and res.j is None, _ok(res)


@pytest.mark.parametrize("bad", ["nan", "inf", "-inf", "abc", ""])
def test_unparsable_or_non_finite_score_names_the_line(tmp_path: Path, bad: str) -> None:
    rows = [Row("PEPTLDEK", 0.95), Row("SAMPLERK", 0.9), Row("ANOTHERK", bad)]
    case = make_case(tmp_path, rows, _empty_genome(random.Random(18)))
    res = case.evaluate()
    assert res.rc == 1 and res.j is None, _ok(res)
    err = res.stderr.replace(str(case.ssl), "")
    assert re.search(r"(?<![\d.])4(?![\d.])", err), "name the line (header = line 1): " + err


def test_unknown_option_and_missing_genome_are_preconditions(tmp_path: Path) -> None:
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(19)))
    assert case.evaluate("--no-such-option").rc == 2
    case.genome.unlink()
    res = case.evaluate()
    assert res.rc == 2 and res.j is None, _ok(res)


def test_script_runs_as_a_subprocess(tmp_path: Path) -> None:
    has_subcommand("evaluate")
    case = ranked_case(tmp_path, 520, range(200), seed=20)
    out = tmp_path / "J.json"
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": str(REPO / "src"), "LANG": "C.UTF-8"}
    cmd = [sys.executable, str(A3_FEASIBILITY), *case.argv(out)]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600, check=False)
    assert res.returncode == 0, res.stdout + res.stderr
    assert json.loads(out.read_text())["n_pass"] == 520


# ================================================================== examples 7-9: the genome check


_PEP8 = "WHEKMPNR"


def _wrap_contig(topology: str, rng: random.Random) -> Contig:
    nt = coding(_PEP8)  # 24 nt: the last 12 nt and the first 12 nt of the contig
    return Contig("rep1", nt[12:] + random_dna(301 - 24, rng) + nt[:12], topology)


@pytest.mark.parametrize(("topology", "want"), [("circular", 1), ("linear", 0)])
def test_example7_origin_spanning_peptide(tmp_path: Path, topology: str, want: int) -> None:
    rng = random.Random(21)
    genome = Genome([_wrap_contig(topology, rng)])
    case = make_case(tmp_path, [Row(_PEP8, 0.95)], genome)
    j = _j(case.evaluate())
    assert j["n_match"] == want
    assert j["matches_by_replicon"].get("rep1", 0) == want


def test_example8_minus_strand(tmp_path: Path) -> None:
    rng = random.Random(22)
    pep = "DEFHKMNPQW"  # 30 nt, planted as its reverse complement
    for lead in range(3):  # all three minus-strand phases
        contig = Contig("rep1", random_dna(60 + lead, rng) + revcomp(coding(pep))
                        + random_dna(60, rng))  # fmt: skip
        case = make_case(tmp_path / str(lead), [Row(pep, 0.95)], Genome([contig]))
        assert _j(case.evaluate())["n_match"] == 1, f"phase lead {lead}"


def test_example9_in_frame_stop_never_matches(tmp_path: Path) -> None:
    rng = random.Random(23)
    stop_pep, ok_pep = "MKDEQWHNR", "MKDEQWHNY"
    nt = coding(stop_pep)
    nt_stop = nt[:12] + "TAA" + nt[15:]  # Q (CAG) -> TAA: the peptide "reads" it as Q
    assert stop_pep[4] == "Q"
    contig = Contig("rep1", random_dna(60, rng) + nt_stop + random_dna(9, rng) + coding(ok_pep)
                    + random_dna(60, rng))  # fmt: skip
    case = make_case(tmp_path, [Row(stop_pep, 0.95), Row(ok_pep, 0.94)], Genome([contig]))
    j = _j(case.evaluate())
    assert j["n_match"] == 1


def test_x_codons_and_nonstandard_letters_never_match(tmp_path: Path) -> None:
    rng = random.Random(24)
    pep = "KDEFHMNPQW"
    nt = coding(pep)
    contig = Contig("rep1", random_dna(60, rng) + nt[:6] + "NNN" + nt[9:] + random_dna(60, rng))
    nonstd = ["PEPBTLDEK", "PEPJTLDEK", "PEPOTLDEK", "PEPUTLDEK", "PEPXTLDEK", "PEPZTLDEK"]
    rows = [Row(pep, 0.95)] + [Row(s, 0.9) for s in nonstd]
    case = make_case(tmp_path, rows, Genome([contig]))
    j = _j(case.evaluate())
    assert j["n_pass"] == 7, "nonstandard letters stay in the counts"
    assert j["n_nonstandard"] == 6
    assert j["n_match"] == 0


def test_several_loci_count_once_and_matches_per_replicon(tmp_path: Path) -> None:
    rng = random.Random(25)
    a = random_peptides(4, rng)
    b = random_peptides(2, rng, avoid=a)
    rep_a = planted_contig("repA", [a[0], a[0], a[1], a[2]], rng, minus=[a[0]])
    rep_b = planted_contig("repB", b, rng, topology="circular")
    rep_b.kind = "plasmid"
    rows = [Row(p, 0.95 - 0.001 * i) for i, p in enumerate(a + b)]
    genome = Genome([rep_a, rep_b])
    case = make_case(tmp_path, rows, genome)
    j = _j(case.evaluate())
    assert j["n_match"] == 5  # a[3] is not planted
    assert j["matches_by_replicon"] == {"repA": 3, "repB": 2}
    assert j["matches_by_replicon"] == _oracle(genome, a + b)["by_rep"]


def test_decoy_and_contaminant_shares(tmp_path: Path) -> None:
    rng = random.Random(26)
    poly = "AAAAAAAAA"  # a shuffled poly-Ala window is still poly-Ala: a decoy hit by design
    peps = random_peptides(8, rng)
    cont = {"CONT_KRT": "MSRQ" + peps[6] + "GGSG" + peps[7] + "FSRS"}
    contig = planted_contig("rep1", [*peps[:4]], rng)
    contig.seq += "GCT" * 150 + random_dna(30, rng)
    genome = Genome([contig])
    rows = [Row(p, 0.95 - 0.001 * i) for i, p in enumerate([poly, *peps])]
    case = make_case(tmp_path, rows, genome, contaminants=cont)
    j = _j(case.evaluate())
    oracle = _oracle(genome, [poly, *peps])
    assert poly in oracle["decoy"], "fixture: the poly-Ala window must give a decoy hit"
    assert j["n_top"] == 9 and j["n_match"] == len(oracle["target"]) == 5
    assert _close(j["decoy_share"], len(oracle["decoy"]) / 9, 1e-9)
    assert _close(j["contaminant_share"], 2 / 9, 1e-9)


# ================================================================== module structure (criterion 2)


def _imports(path: Path) -> dict[str, set[str]]:
    tree = ast.parse(path.read_text())
    out: dict[str, set[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            out.setdefault(node.module, set()).update(a.name for a in node.names)
        elif isinstance(node, ast.Import):
            for a in node.names:
                out.setdefault(a.name, set())
    return out


def test_module_reuses_gippa_matching_and_wilson() -> None:
    imp = _imports(require(A3_FEASIBILITY, "deliverable 9"))
    names = {n for mod, ns in imp.items() if mod.startswith("gippa.") for n in ns}
    mods = set(imp)
    for need in ("build_sixframe", "build_decoys", "match_exact", "wilson_ci", "for_step"):
        assert need in names or any(m.endswith(need) for m in mods), f"must import {need}"
    assert "match_contaminant_peptides" in names, "contaminant_share uses the GIPPA matcher"
    assert "parse_bakta_json" in names, "replicons come from the Bakta JSON (step 8)"
    assert not {"ahocorasick", "cascadia", "pyteomics", "re2"} & {m.split(".")[0] for m in mods}
    src = A3_FEASIBILITY.read_text()
    for own in ("def match_exact", "def build_sixframe", "def translate", "def wilson"):
        assert own not in src, f"defines its own matcher or interval: {own}"
    assert "import csv" in src or "from csv" in src, "the .ssl is read with the stdlib csv module"


# ================================================================== early-fail (criterion 5)


@pytest.mark.parametrize("reason", REASONS)
def test_early_fail_writes_the_schema(tmp_path: Path, reason: str) -> None:
    has_subcommand("early-fail")
    out = tmp_path / "J.json"
    text = "Khoa declined the TDF-SDK terms (PROGRESS 2026-10-24)"
    res = call(["early-fail", "--reason", reason, "--evidence", text, "--out", str(out)], out)
    assert res.rc == 3, _ok(res)
    j = _j(res)
    assert res.raw is not None and is_canonical_json(res.raw)
    assert set(j) >= JSON_FIELDS, sorted(JSON_FIELDS - set(j))
    assert j["verdict"] == "fail" and j["reason"] == reason and j["evidence"] == text
    assert all(j[k] is None for k in COUNT_FIELDS), {k: j[k] for k in COUNT_FIELDS}
    assert j["fails"] == {"too_few_peptides": None, "low_genome_match": None}
    assert j["rule"] == RULE and j["schema_version"] == 1


@pytest.mark.parametrize("reason", ["too_few_peptides", "both", "gpu_queue", ""])
def test_early_fail_rejects_other_reasons(tmp_path: Path, reason: str) -> None:
    has_subcommand("early-fail")
    out = tmp_path / "J.json"
    res = call(["early-fail", "--reason", reason, "--evidence", "x", "--out", str(out)], out)
    assert res.rc == 2 and not out.exists(), _ok(res)


# ================================================================== Hypothesis: counts = oracle

_VOCAB_PLANTED = ["WHEKMPNRAL", "DEFHKMNPQW", "TTGGKDDEER", "LLMMPPQQRR"]
_VOCAB_ABSENT = ["PEPTIDEKAA", "SAMPLERKKT", "QWERTYKLMN", "PEPB" + "TLDEKR"]
_VOCAB_EDGE = ["AAAAAA", "ACDEFGH", "W" * 30, "Y" * 31]  # lengths 6, 7, 30, 31
_MODS = ["", "[+15.99]", "(ox)", "+42.011", "-", "[+42.01]-"]
_SCORES = [0.5, 0.8, 0.8000001, 0.81, 0.9, 0.95, 0.99]


@st.composite
def _rows(draw: st.DrawFn) -> list[tuple[str, float]]:
    n = draw(st.integers(0, 25))
    out = []
    for _ in range(n):
        base = draw(st.sampled_from(_VOCAB_PLANTED + _VOCAB_ABSENT + _VOCAB_EDGE))
        if draw(st.booleans()):
            base = base.replace("L", "I")  # an I-form that folds back
        k = draw(st.integers(0, len(base)))
        seq = draw(st.sampled_from(_MODS)) + base[:k] + draw(st.sampled_from(_MODS)) + base[k:]
        out.append((seq, draw(st.sampled_from(_SCORES))))
    return out


_PROP_GENOME = Genome([planted_contig("rep1", _VOCAB_PLANTED, random.Random(27))])


@settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_rows(), perm_seed=st.integers(0, 2**16))
def test_property_counts_match_an_independent_oracle(rows: list[tuple[str, float]],
                                                     perm_seed: int) -> None:  # fmt: skip
    has_subcommand("evaluate")
    # the .ssl holds scores to 6 decimals (write_ssl), so the oracle reads what the file says
    kept = [(fold(s), float(f"{sc:.6f}")) for s, sc in rows]
    in_len = [(f, sc) for f, sc in kept if 7 <= len(f) <= 30]
    passing = {f for f, sc in in_len if sc > 0.8}
    oracle = _oracle(_PROP_GENOME, sorted(passing))
    with tempfile.TemporaryDirectory(prefix="gippa-a3-prop-") as tmp:
        d = Path(tmp)
        case = make_case(d / "a", [Row(s, sc) for s, sc in rows], _PROP_GENOME)
        res = case.evaluate()
        shuffled = list(rows)
        random.Random(perm_seed).shuffle(shuffled)
        case_b = make_case(d / "b", [Row(s, sc) for s, sc in shuffled], _PROP_GENOME)
        res_b = case_b.evaluate()
    j = _j(res)
    assert j["n_rows"] == len(rows)
    assert j["n_len_7_30"] == len(in_len)
    assert j["n_below_cutoff"] == sum(1 for _, sc in in_len if sc <= 0.8)
    assert j["n_pass"] == len(passing) and j["n_top"] == len(passing)
    assert j["n_match"] == len(oracle["target"])
    assert j["n_nonstandard"] == sum(1 for f in passing if set(f) - set(RESIDUES))
    assert {int(k): v for k, v in j["length_hist"].items()} == dict(Counter(map(len, passing)))
    assert {k: v for k, v in j["matches_by_replicon"].items() if v} == oracle["by_rep"]
    if passing:
        assert _close(j["share"], len(oracle["target"]) / len(passing), 1e-9)
    else:
        assert j["share"] is None
    assert res.rc == 3  # always fewer than 500 here
    assert res_b.raw == res.raw, "J depends on the row order"


# ================================================================== the pinned minimal CLI (lead)

NISSLE = REPO / "tests" / "data" / "genomes" / "nissle1917.fna.gz"


def test_minimal_cli_missing_ssl_is_a_precondition(tmp_path: Path) -> None:
    # Lead decision (2026-10-24): `evaluate <RUN.ssl> --run <RUN> --out <J>` as M10-10's
    # cascadia_convert.sbatch calls it; the genome, Bakta JSON and contaminants have defaults.
    has_subcommand("evaluate")
    out = tmp_path / "J.json"
    res = call(["evaluate", str(tmp_path / "WT_1.ssl"), "--run", "WT_1", "--out", str(out)], out)
    assert res.rc == 2 and not out.exists(), _ok(res)


def test_minimal_cli_defaults_to_the_nissle_genome(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    has_subcommand("evaluate")
    db = tmp_path / "db"  # never the real $GIPPA_DB_DIR
    (db / "contaminants").mkdir(parents=True)
    from helpers.a3_eval_kit import write_contaminants, write_ssl

    write_contaminants(db / "contaminants" / "contaminants.fasta")
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    ssl = write_ssl(tmp_path / "WT_1.ssl", [])
    out = tmp_path / "J.json"
    res = call(["evaluate", str(ssl), "--run", "WT_1", "--out", str(out)], out)
    assert res.rc == 3, _ok(res)
    j = _j(res)
    assert j["genome_sha256"] == sha256_hex(NISSLE.read_bytes())
    assert j["n_rows"] == 0 and j["verdict"] == "fail"


# ================================================================== review round 1 (M10-08b)
# Findings of the 08b code review, lead decisions of 2026-10-24: S3, S4a, N2, N6, N7, N8.


def _evaluate_or_fail(case: Case, *extra: str) -> Any:
    """``case.evaluate`` that turns an uncaught exception into a test failure naming it."""
    try:
        return case.evaluate(*extra)
    except Exception as exc:  # a traceback is exit 1, never the spec's exit 2
        pytest.fail(f"evaluate raised {type(exc).__name__}: {exc} (want an exit code)")


@pytest.mark.parametrize("named", ["M[Oxidation]PEPTLDEK", "[Acetyl]-PEPTLDEKR",
                                   "PEPC[UNIMOD:4]TLDEK", "PEPTLDEK[]R"])  # fmt: skip
def test_s3_a_named_modification_without_mass_is_an_input_error(tmp_path: Path, named: str) -> None:
    # Spec edge case "Unimod": the .ssl has modified sequences without mass shifts -> exit 1
    # (Cascadia's convert_mods_to_shifts leaves the name when the local table lacks it).
    rows = [Row("PEPTLDEK", 0.95), Row("SAMPLERK", 0.9), Row(named, 0.9)]
    case = make_case(tmp_path, rows, _empty_genome(random.Random(40)))
    res = _evaluate_or_fail(case)
    assert res.rc == 1 and res.j is None, _ok(res)
    err = res.stderr.replace(str(case.ssl), "")
    assert re.search(r"(?<![\d.])4(?![\d.])", err), "name the line (header = line 1): " + err


def test_s3_input_error_is_the_modules_input_error() -> None:
    has_subcommand("evaluate")
    from helpers.a3_eval_kit import module, write_ssl

    mod = module()
    with tempfile.TemporaryDirectory() as d:
        p = write_ssl(Path(d) / "WT_1.ssl", [Row("M[Oxidation]PEPTLDEK", 0.95)])
        with pytest.raises(mod.InputError):
            mod.tally(mod.read_ssl(p), 0.8)


@pytest.mark.parametrize("shifted", ["M[+15.995]PEPTLDEK", "PEPTLDEK[-17.027]R",
                                     "[+42.011]-PEPTLDEKR", "M[+16]PEPTLDEK"])  # fmt: skip
def test_s3_signed_mass_shifts_are_accepted(tmp_path: Path, shifted: str) -> None:
    case = make_case(tmp_path, [Row(shifted, 0.95)], _empty_genome(random.Random(41)))
    res = case.evaluate()
    assert res.rc == 3, _ok(res)
    assert _j(res)["n_pass"] == 1


def test_s4a_a_contaminant_also_in_the_genome_is_a_non_match(tmp_path: Path) -> None:
    # 08b amendment: a contaminant counts as a non-match even if it also matches the genome.
    rng = random.Random(42)
    peps = random_peptides(6, rng)
    both, only_cont = peps[4], peps[5]
    cont = {"CONT_KRT": "MSRQ" + both + "GGSG" + only_cont + "FSRS"}
    genome = Genome([planted_contig("rep1", peps[:5], rng)])  # peps[4] planted too
    rows = [Row(p, 0.95 - 0.001 * i) for i, p in enumerate(peps)]
    case = make_case(tmp_path, rows, genome, contaminants=cont)
    j = _j(case.evaluate())
    oracle = _oracle(genome, peps)
    assert both in oracle["target"], "fixture: the contaminant peptide must match the genome"
    assert len(oracle["target"]) == 5
    assert j["n_top"] == 6 and j["n_match"] == 4, "the contaminant hit was counted as a match"
    assert _close(j["share"], 4 / 6, 1e-9)
    lo, hi = wilson_ci(4, 6)
    assert _close(j["share_ci95"][0], lo, 1e-9) and _close(j["share_ci95"][1], hi, 1e-9)
    assert j["matches_by_replicon"] == {"rep1": 4}, "matches_by_replicon counts the contaminant"
    assert _close(j["contaminant_share"], 2 / 6, 1e-9), "still counted as a contaminant"


def test_n2_nonstandard_counts_every_distinct_passing_peptide(tmp_path: Path) -> None:
    # 1,000 standard peptides fill the top 1,000; three nonstandard ones rank after it.
    rng = random.Random(43)
    std = random_peptides(1000, rng)
    rows = [Row(p, s) for p, s in zip(std, ladder(1000, top=0.99), strict=True)]
    rows += [Row(s, 0.81) for s in ("PEPBTLDEK", "PEPXTLDEK", "PEPZTLDEK")]
    rows += [Row("PEPUTLDEKR", 0.999)]  # one inside the top
    case = make_case(tmp_path, rows, _empty_genome(rng))
    j = _j(case.evaluate())
    assert j["n_pass"] == 1004 and j["n_top"] == 1000
    assert j["n_nonstandard"] == 4, "n_nonstandard counts the distinct passing set, not the top"


def test_n6_missing_mzml_check_file_is_a_precondition(tmp_path: Path) -> None:
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(44)))
    res = _evaluate_or_fail(case, "--mzml-check", str(tmp_path / "no_such_check.json"))
    assert res.rc == 2 and res.j is None, _ok(res)
    assert "no_such_check.json" in res.stderr, "name the missing file: " + res.stderr


def test_n6_mzml_check_rules_are_copied(tmp_path: Path) -> None:
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(45)))
    check = tmp_path / "WT_1.mzml_check.json"
    rules = {f"M{i}": True for i in range(1, 9)}
    check.write_text(json.dumps({"run": "WT_1", "rules": rules, "pass": True, "failed": []}))
    res = case.evaluate("--mzml-check", str(check))
    assert res.rc == 3, _ok(res)
    assert _j(res)["mzml_check"] == rules


@pytest.mark.parametrize("body", [" ", "\n", "\n\n", "  \t \n \n", "\r\n"])
def test_n6_whitespace_only_ssl_is_a_precondition(tmp_path: Path, body: str) -> None:
    # Spec exit codes: 2 for an empty file; a whitespace-only .ssl holds no header either.
    case = make_case(tmp_path, [Row("PEPTLDEK", 0.95)], _empty_genome(random.Random(46)))
    case.ssl.write_text(body)
    res = _evaluate_or_fail(case)
    assert res.rc == 2 and res.j is None, _ok(res)


def test_n7_without_date_the_json_is_stable_within_a_day(tmp_path: Path) -> None:
    import datetime

    case = ranked_case(tmp_path, 520, range(200), seed=47)
    before = datetime.datetime.now(datetime.UTC).date().isoformat()
    a = case.evaluate(out=tmp_path / "a.json")
    b = case.evaluate(out=tmp_path / "b.json")
    after = datetime.datetime.now(datetime.UTC).date().isoformat()
    assert _j(a)["date"] in (before, after), "the default date is today (UTC)"
    assert _j(b)["date"] == _j(a)["date"] or before != after
    if before == after:
        assert a.raw == b.raw


def test_n7_with_date_the_json_is_byte_identical(tmp_path: Path) -> None:
    case = ranked_case(tmp_path, 520, range(200), seed=48)
    a = case.evaluate("--date", "2026-10-24", out=tmp_path / "a.json")
    b = case.evaluate("--date", "2026-10-24", out=tmp_path / "b.json")
    assert _j(a)["date"] == "2026-10-24"
    assert a.raw == b.raw


def test_n7_readme_evaluate_command_passes_date() -> None:
    text = require(REPO / "validation" / "pilot" / "README.md", "deliverable 12").read_text()
    start = text.find("a3_feasibility.py evaluate ")
    assert start >= 0, "the README has no a3_feasibility.py evaluate command"
    cmd = text[start : text.index("```", start)]
    assert "--date" in cmd, "the README's evaluate command must pass --date:\n" + cmd


def _git_ignored(path: Path) -> bool:
    res = subprocess.run(["git", "-C", str(REPO), "check-ignore", "-q", str(path)],
                         capture_output=True, check=False)  # fmt: skip
    assert res.returncode in (0, 1), res.stderr
    return res.returncode == 0


def test_n8_detail_inside_the_checkout_not_git_ignored_is_refused(tmp_path: Path) -> None:
    import uuid

    detail = REPO / "validation" / "pilot" / f"_a3_detail_test_{uuid.uuid4().hex}.tsv"
    assert not _git_ignored(detail), "fixture: the path must not be git-ignored"
    case = ranked_case(tmp_path, 30, range(10), seed=49)
    try:
        res = _evaluate_or_fail(case, "--detail", str(detail))
        written = detail.exists()
    finally:
        detail.unlink(missing_ok=True)
    assert res.rc == 2 and res.j is None, _ok(res)
    assert not written, "the per-peptide list was written inside the checkout"


def test_n8_detail_outside_the_checkout_is_written(tmp_path: Path) -> None:
    case = ranked_case(tmp_path, 30, range(10), seed=50)
    detail = tmp_path / "gippa_root" / "pilot" / "feasibility_top1000.tsv"
    res = case.evaluate("--detail", str(detail))
    j = _j(res)
    lines = detail.read_text().splitlines()
    assert lines[0].split("\t")[:2] == ["rank", "sequence"]
    assert len(lines) == 1 + j["n_top"]
    match = lines[0].split("\t").index("match")
    assert sum(int(ln.split("\t")[match]) for ln in lines[1:]) == j["n_match"]


def test_n8_detail_in_a_git_ignored_path_of_the_checkout_is_allowed(tmp_path: Path) -> None:
    import shutil
    import uuid

    d = REPO / "local" / f"_a3_detail_test_{uuid.uuid4().hex}"
    detail = d / "feasibility_top1000.tsv"
    assert _git_ignored(detail), "fixture: local/ is git-ignored"
    case = ranked_case(tmp_path, 30, range(10), seed=51)
    try:
        res = case.evaluate("--detail", str(detail))
        written = detail.is_file()
    finally:
        shutil.rmtree(d, ignore_errors=True)
    assert res.rc in (0, 3), _ok(res)
    assert written
