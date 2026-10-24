"""``validation/pilot/reproduce_check.py`` (``parse``, ``static_check``, ``verify``, CLI), the
committed draft ``validation/pilot/reproduce.md`` and the ``pilot_reproduce_check.sbatch`` header.

Spec: docs/specs/M10-14-pilot-reproducibility.md, "The runbook format", "S-check", examples 1-6,
edge cases "Empty input" and "Ties", acceptance criteria 1 (a)-(h), 2 and 4; lead decisions
(2026-10-24): criterion 1 (i) (no 3+-digit numbers outside provenance columns) is dropped.

Fixtures are hand-made runbooks and a synthetic checkout in ``tmp_path``
(``helpers.m10_14_kit``); ``verify`` runs against fake ``$GIAP_ROOT`` files in ``tmp_path``.
Nothing reads the real ``$GIAP_ROOT``, Slurm or the network. The conventions the spec leaves open
are listed in the docstring of ``helpers.m10_14_kit``.

``reproduce_check`` moves files and compares hashes; it makes no statistical claim, so this module
has no calibration test (spec §8.4).
"""

from __future__ import annotations

import fnmatch
import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.m10_14_kit import (
    DIANN_SHA,
    EVIDENCE_S5,
    EVIDENCE_S8_PENDING,
    GENOME_REL,
    GENOME_SHA,
    INPUT_GENOME,
    INPUT_IMAGE,
    MANIFEST_REL,
    OTHER_COMMIT,
    PIN_REL,
    REPRO_SBATCH,
    REPRODUCE_MD,
    STEP_S0,
    STEP_S7,
    F,
    check_house_rules,
    evidence_s8,
    load_reproduce_check,
    run_cli,
    runbook,
    s5,
    sbatch_code,
    sbatch_header,
    sha256_bytes,
    sha256_file,
    synthetic_repo,
)
from helpers.slurm_sandbox import REPO
from hypothesis import given, settings
from hypothesis import strategies as st

HEX64 = re.compile(r"(?<![0-9a-fA-F])[0-9a-f]{64}(?![0-9a-fA-F])")
STATUSES = {"ok", "mismatch", "absent", "not_checkable"}


@pytest.fixture(scope="module")
def rc() -> Any:
    return load_reproduce_check()


def _issues(rc: Any, text: str, repo: Path) -> list[str]:
    return [s for s in rc.static_check(rc.parse(text), repo) if not s.startswith("duplicate_hash:")]


# ---------------------------------------------------------------- parse


def test_parse_reads_header_tables_and_steps(rc: Any) -> None:
    rb = rc.parse(runbook(status="draft", freeze="pending", a3_branch="pending"))
    assert (rb.status, rb.freeze_commit, rb.a3_branch) == ("draft", "pending", "pending")
    assert [r["id"] for r in rb.inputs] == ["genome", "diann_image"]
    assert rb.inputs[0] == {
        "id": "genome",
        "path": GENOME_REL,
        "sha256": GENOME_SHA,
        "recorded_in": MANIFEST_REL,
    }
    assert [s.id for s in rb.steps] == ["S0", "S5", "S7", "S8"]
    s7 = rb.steps[2]
    assert s7.title == "conversion" and s7.issue == "M10-03"
    assert s7.job == {
        "script": "scripts/slurm/diann_convert.sbatch",
        "name": "giap-diann-convert-A1",
        "partition": "normal",
        "qos": "normal",
        "cpus": "2",
        "gpus": "0",
        "mem": "8G",
        "time": "00:30:00",
        "constraint": "-",
        "array": "-",
    }
    assert any("diann_convert.sbatch" in c for c in s7.commands)
    assert s7.reads == ["genome"]
    assert s7.writes == ["$GIAP_ROOT/pilot/inputs/a1.tsv.gz"]
    assert rb.steps[0].job is None, "`Job: none (session)` parses as no job"
    assert rb.steps[3].job["children"] == "nextflow:normalx14cpu"
    assert [e["step"] for e in rb.evidence] == ["S5", "S8"]
    assert rb.evidence[0]["job"] == "1879913"
    assert rb.outputs == []
    assert rb.a3 == {"cascadia": "inactive", "fallback": "inactive"}


def test_parse_reads_outputs_and_final_header(rc: Any) -> None:
    out = f"| a2_table | $GIAP_ROOT/pilot/inputs/a2.tsv.gz | {'d' * 64} | - | yes |"
    rb = rc.parse(
        runbook(status="final", freeze=F, a3_branch="cascadia", cascadia="active", outputs=[out])
    )
    assert (rb.status, rb.freeze_commit, rb.a3_branch) == ("final", F, "cascadia")
    assert rb.outputs == [
        {
            "id": "a2_table",
            "path": "$GIAP_ROOT/pilot/inputs/a2.tsv.gz",
            "sha256": "d" * 64,
            "recorded_in": "-",
            "kept": "yes",
        }
    ]
    assert rb.a3 == {"cascadia": "active", "fallback": "inactive"}


_ID = st.from_regex(r"[a-z][a-z0-9_]{0,11}", fullmatch=True)
_HEX = st.text(alphabet="0123456789abcdef", min_size=64, max_size=64)


@settings(max_examples=40, deadline=None)
@given(rows=st.dictionaries(_ID, _HEX, min_size=1, max_size=6))
def test_property_inputs_round_trip(rows: dict[str, str]) -> None:
    """Any well-formed Inputs table parses back to its ids and hashes, in order."""
    rc = load_reproduce_check()
    lines = [
        f"| {i} | $GIAP_ROOT/x/{i}.bin | {h} | $GIAP_ROOT/x/provenance.json |"
        for i, h in rows.items()
    ]
    rb = rc.parse(runbook(inputs=lines))
    assert [(r["id"], r["sha256"]) for r in rb.inputs] == list(rows.items())


# ---------------------------------------------------------------- static_check: the examples


def test_default_fixture_passes(rc: Any, tmp_path: Path) -> None:
    """Example 1 (S7 verbatim): the script exists, -J is giap-*, every field equals the header."""
    repo = synthetic_repo(tmp_path)
    assert _issues(rc, runbook(), repo) == []


def test_example_2_missing_script(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    step = s5(script="scripts/slurm/diann_search.sbatch").replace(
        "script=scripts/slurm/diann_pilot.sbatch", "script=scripts/slurm/diann_search.sbatch"
    )
    issues = _issues(rc, runbook(steps=[STEP_S0, step, STEP_S7]), repo)
    assert "missing_script:scripts/slurm/diann_search.sbatch (S5)" in issues, issues


def test_example_3_job_class_mismatch(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(mem="64G"), STEP_S7]), repo)
    assert "job_class_mismatch:S5:mem:64G!=64000M" in issues, issues


def test_command_override_is_the_expected_value(rc: Any, tmp_path: Path) -> None:
    """Criterion 1 (c): a field equals the script's #SBATCH value *or the command's override*."""
    repo = synthetic_repo(tmp_path)
    ok = s5(mem="32000M", extra="--mem=32000M ")
    assert _issues(rc, runbook(steps=[STEP_S0, ok, STEP_S7]), repo) == []
    bad = s5(mem="64000M", extra="--mem=32000M ")
    issues = _issues(rc, runbook(steps=[STEP_S0, bad, STEP_S7]), repo)
    assert "job_class_mismatch:S5:mem:64000M!=32000M" in issues, issues


def test_job_name_must_start_with_giap(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    step = s5(jname="diann-A1").replace("name=giap-diann-A1", "name=diann-A1")
    issues = _issues(rc, runbook(steps=[STEP_S0, step, STEP_S7]), repo)
    assert any("diann-A1" in s and "S5" in s for s in issues), issues


def test_example_4_inputs_row(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    assert sha256_file(repo / GENOME_REL) == GENOME_SHA
    assert _issues(rc, runbook(inputs=[INPUT_GENOME]), repo) == []
    bad = INPUT_GENOME.replace(GENOME_SHA, GENOME_SHA[:-1] + "f")
    issues = _issues(rc, runbook(inputs=[bad]), repo)
    assert "hash_mismatch:genome" in issues, issues
    assert f"not_recorded:genome:{MANIFEST_REL}" in issues, issues


def test_committed_input_is_not_its_own_record(rc: Any, tmp_path: Path) -> None:
    """A committed file is not its own record: recorded_in must name a second file."""
    repo = synthetic_repo(tmp_path)
    selfrec = f"| genome | {GENOME_REL} | {GENOME_SHA} | {GENOME_REL} |"
    issues = _issues(rc, runbook(inputs=[selfrec]), repo)
    assert any("genome" in s for s in issues), issues


def test_malformed_inputs_rows(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    upper = INPUT_IMAGE.replace(DIANN_SHA, DIANN_SHA.upper())
    assert any("diann_image" in s for s in _issues(rc, runbook(inputs=[INPUT_GENOME, upper]), repo))
    bad_id = INPUT_IMAGE.replace("| diann_image |", "| DIA-NN image |")
    assert _issues(rc, runbook(inputs=[INPUT_GENOME, bad_id]), repo)
    abspath = INPUT_IMAGE.replace("$GIAP_ROOT/images", "/home/someone/images")
    assert _issues(rc, runbook(inputs=[INPUT_GENOME, abspath]), repo)


def test_example_5_static_reads_the_pin(rc: Any, tmp_path: Path) -> None:
    """The image row passes static (its pin is in diann_image.sha256); a changed pin fails."""
    repo = synthetic_repo(tmp_path)
    assert _issues(rc, runbook(), repo) == []
    (repo / PIN_REL).write_text(f"{'0' * 64}  diann_2.3.1.sif\n")
    issues = _issues(rc, runbook(), repo)
    assert f"not_recorded:diann_image:{PIN_REL}" in issues, issues


def test_example_6_a3_branch_mismatch(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path, verdict="fail")
    text = runbook(a3_branch="cascadia", cascadia="active")
    issues = _issues(rc, text, repo)
    assert "a3_branch_mismatch:verdict=fail,active=cascadia" in issues, issues


def test_a3_branch_agrees(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path, verdict="pass")
    assert _issues(rc, runbook(a3_branch="cascadia", cascadia="active"), repo) == []
    repo2 = synthetic_repo(tmp_path / "b", verdict="fail")
    assert _issues(rc, runbook(a3_branch="fallback", fallback="active"), repo2) == []


def test_a3_header_must_match_the_single_active_subsection(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path, verdict="pass")
    assert _issues(rc, runbook(a3_branch="cascadia", fallback="active"), repo)
    assert _issues(rc, runbook(a3_branch="cascadia", cascadia="active", fallback="active"), repo)


def test_empty_runbook_parts(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    assert "no_steps" in _issues(rc, runbook(steps=[]), repo)
    assert "no_inputs" in _issues(rc, runbook(inputs=[]), repo)


def test_duplicate_step(rc: Any, tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), s5(), STEP_S7]), repo)
    assert any(s.startswith("duplicate_step") for s in issues), issues


def test_duplicate_hash_is_a_warning(tmp_path: Path) -> None:
    """Edge case "Ties": the same hash on two rows with different paths is reported, not failed."""
    repo = synthetic_repo(tmp_path)
    copy = (
        f"| genome_copy | $GIAP_ROOT/pilot/inputs/genome.fna.gz | {GENOME_SHA} | {MANIFEST_REL} |"
    )
    p = tmp_path / "reproduce.md"
    p.write_text(runbook(inputs=[INPUT_GENOME, INPUT_IMAGE, copy]))
    res = run_cli("static", str(p), "--repo", str(repo))
    assert res.returncode == 0, res.stdout + res.stderr
    assert "duplicate_hash:genome,genome_copy" in res.stdout + res.stderr


def test_missing_path_token(rc: Any, tmp_path: Path) -> None:
    """Criterion 1 (a): scripts/ validation/ workflow/ tests/data/ tokens must exist; globs pass."""
    repo = synthetic_repo(tmp_path)
    missing = STEP_S7.replace('" A1', '" A1 validation/pilot/diann_A9.cfg')
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), missing]), repo)
    assert any("validation/pilot/diann_A9.cfg" in s for s in issues), issues
    glob = STEP_S7.replace('" A1', '" A1 validation/pilot/diann_A*.cfg')
    assert _issues(rc, runbook(steps=[STEP_S0, s5(), glob]), repo) == []


def test_stray_64_hex(rc: Any, tmp_path: Path) -> None:
    """Criterion 1 (e): every 64-hex string is an Inputs/Outputs row (or the freeze commit)."""
    repo = synthetic_repo(tmp_path)
    stray = "ab12" * 16
    issues = _issues(rc, runbook(prose=f"An old hash {stray} kept for the record."), repo)
    assert any(stray[:12] in s for s in issues), issues
    out = f"| a1_table | $GIAP_ROOT/pilot/inputs/a1.tsv.gz | {stray} | - | yes |"
    assert _issues(rc, runbook(prose=f"The table hashes to {stray}.", outputs=[out]), repo) == []


@pytest.mark.parametrize("cmd", ["uv run giap --help", "uv sync --frozen", "pip install giap"])
def test_forbidden_commands(rc: Any, tmp_path: Path, cmd: str) -> None:
    repo = synthetic_repo(tmp_path)
    step = STEP_S7.replace("```bash\n", f"```bash\n{cmd}\n")
    issues = _issues(rc, runbook(steps=[STEP_S0, s5(), step]), repo)
    assert any("S7" in s for s in issues), issues


def _final(**kw: Any) -> str:
    base: dict[str, Any] = dict(
        status="final",
        freeze=F,
        a3_branch="fallback",
        fallback="active",
        evidence=[EVIDENCE_S5, evidence_s8()],
    )
    base.update(kw)
    return runbook(**base)


def test_final_runbook(rc: Any, tmp_path: Path) -> None:
    """Criterion 1 (h): final means no pending/TBD/empty Evidence cell, a 40-hex freeze commit,
    and every S8-S12 Evidence commit equal to it; the same pending rows pass as a draft."""
    repo = synthetic_repo(tmp_path, verdict="fail")
    assert _issues(rc, _final(), repo) == []
    pending = [EVIDENCE_S5, EVIDENCE_S8_PENDING]
    assert _issues(rc, runbook(evidence=pending), repo) == []
    assert _issues(rc, _final(evidence=pending), repo)
    assert _issues(
        rc, _final(evidence=[EVIDENCE_S5, evidence_s8().replace("1900001", "TBD")]), repo
    )
    assert _issues(
        rc, _final(evidence=[EVIDENCE_S5, evidence_s8().replace("| 7G |", "|  |")]), repo
    )
    assert _issues(rc, _final(freeze="pending"), repo)
    assert _issues(rc, _final(freeze=F[:7]), repo)
    issues = _issues(rc, _final(evidence=[EVIDENCE_S5, evidence_s8(OTHER_COMMIT)]), repo)
    assert any("S8" in s for s in issues), issues


def test_cli_static_exit_codes(tmp_path: Path) -> None:
    repo = synthetic_repo(tmp_path)
    good, bad = tmp_path / "good.md", tmp_path / "bad.md"
    good.write_text(runbook())
    bad.write_text(runbook(steps=[STEP_S0, s5(mem="64G"), STEP_S7]))
    res = run_cli("static", str(good), "--repo", str(repo))
    assert res.returncode == 0, res.stdout + res.stderr
    res = run_cli("static", str(bad), "--repo", str(repo))
    assert res.returncode == 1, res.stdout + res.stderr
    assert "job_class_mismatch:S5:mem:64G!=64000M" in res.stdout + res.stderr


# ---------------------------------------------------------------- verify (fake $GIAP_ROOT)


def _verify_world(tmp_path: Path) -> tuple[Path, Path, str]:
    repo = synthetic_repo(tmp_path)
    root = tmp_path / "root"
    image = root / "images" / "diann_2.3.1.sif"
    image.parent.mkdir(parents=True)
    image.write_bytes(b"stub sif\n")
    image.chmod(0o444)
    img_sha = sha256_bytes(b"stub sif\n")
    (repo / PIN_REL).write_text(f"{img_sha}  diann_2.3.1.sif\n")
    run = root / "data" / "pxd078459" / "WT_1.d"
    run.mkdir(parents=True)
    (run / "analysis.tdf").write_bytes(b"tdf bytes\n")
    (run / "analysis.tdf_bin").write_bytes(b"tdf_bin bytes\n")
    tdf_sha = sha256_bytes(b"tdf bytes\n")
    prov = root / "pilot" / "A1" / "provenance.json"
    prov.parent.mkdir(parents=True)
    prov.write_text(json.dumps({"raw": [{"run": "WT_1", "analysis_tdf_sha256": tdf_sha}]}))
    rows = [
        INPUT_GENOME,
        f"| diann_image | $GIAP_ROOT/images/diann_2.3.1.sif | {img_sha} | {PIN_REL} |",
        f"| wt1_tdf | $GIAP_ROOT/data/pxd078459/WT_1.d/analysis.tdf | {tdf_sha} | "
        "$GIAP_ROOT/pilot/A1/provenance.json |",
        f"| wt1_tdf_bin | $GIAP_ROOT/data/pxd078459/WT_1.d/analysis.tdf_bin | {'e' * 64} | "
        "$GIAP_ROOT/pilot/A1/provenance.json |",
        f"| a3_table | $GIAP_ROOT/pilot/inputs/a3.tsv.gz | {'f' * 64} | "
        "$GIAP_ROOT/pilot/inputs/a3.convert.json |",
    ]
    return repo, root, runbook(inputs=rows)


def _status(res: dict[str, Any], rid: str) -> str:
    return str(res["inputs"][rid]["status"])


def test_verify_all_ok(rc: Any, tmp_path: Path) -> None:
    repo, root, text = _verify_world(tmp_path)
    res = rc.verify(rc.parse(text), root, repo)
    assert res["failures"] == [], res
    assert {r: _status(res, r) for r in res["inputs"]} == {
        "genome": "ok",
        "diann_image": "ok",
        "wt1_tdf": "ok",
        "wt1_tdf_bin": "not_checkable",
        "a3_table": "absent",
    }
    assert all(v["status"] in STATUSES for v in res["inputs"].values())


def test_verify_altered_image(rc: Any, tmp_path: Path) -> None:
    """Example 5: an altered copy gives ``mismatch`` and a failure (and the CLI exits 1)."""
    repo, root, text = _verify_world(tmp_path)
    image = root / "images" / "diann_2.3.1.sif"
    image.chmod(0o644)
    image.write_bytes(b"altered sif\n")
    image.chmod(0o444)
    res = rc.verify(rc.parse(text), root, repo)
    assert _status(res, "diann_image") == "mismatch"
    assert any("diann_image" in f for f in res["failures"]), res["failures"]
    p, out = tmp_path / "rb.md", tmp_path / "verify.json"
    p.write_text(text)
    cli = run_cli(
        "verify", str(p), "--out", str(out), "--repo", str(repo), env={"GIAP_ROOT": str(root)}
    )
    assert cli.returncode == 1, cli.stdout + cli.stderr
    doc = json.loads(out.read_text())
    assert doc["inputs"]["diann_image"]["status"] == "mismatch"
    assert doc["failures"]


def test_verify_image_must_be_read_only(rc: Any, tmp_path: Path) -> None:
    repo, root, text = _verify_world(tmp_path)
    (root / "images" / "diann_2.3.1.sif").chmod(0o644)
    res = rc.verify(rc.parse(text), root, repo)
    assert any("diann_image" in f for f in res["failures"]), res
    assert _status(res, "diann_image") != "ok"


def test_verify_reads_farmshare_records(rc: Any, tmp_path: Path) -> None:
    """A FarmShare ``recorded_in`` (provenance.json) must contain the runbook's hash."""
    repo, root, text = _verify_world(tmp_path)
    (root / "pilot" / "A1" / "provenance.json").write_text(json.dumps({"raw": []}))
    res = rc.verify(rc.parse(text), root, repo)
    assert any("wt1_tdf" in f for f in res["failures"]), res
    assert _status(res, "wt1_tdf") != "ok"


def test_verify_detects_a_changed_farmshare_file(rc: Any, tmp_path: Path) -> None:
    repo, root, text = _verify_world(tmp_path)
    (root / "data" / "pxd078459" / "WT_1.d" / "analysis.tdf").write_bytes(b"other\n")
    res = rc.verify(rc.parse(text), root, repo)
    assert _status(res, "wt1_tdf") == "mismatch"
    assert any("wt1_tdf" in f for f in res["failures"])


def test_verify_cli_ok(tmp_path: Path) -> None:
    repo, root, text = _verify_world(tmp_path)
    p, out = tmp_path / "rb.md", tmp_path / "out" / "verify.json"
    out.parent.mkdir()
    p.write_text(text)
    cli = run_cli(
        "verify", str(p), "--out", str(out), "--repo", str(repo), env={"GIAP_ROOT": str(root)}
    )
    assert cli.returncode == 0, cli.stdout + cli.stderr
    doc = json.loads(out.read_text())
    assert doc["failures"] == []
    assert doc["inputs"]["genome"]["status"] == "ok"


# ---------------------------------------------------------------- the committed draft runbook

PATH_TOKEN = re.compile(
    r"(?:\$\{?GIAP_REPO\}?/)?((?:scripts|validation|workflow|tests/data)/[^\s\"'`;|)(<>]+)"
)
FENCE = re.compile(r"^```bash\n(.*?)^```", re.M | re.S)
JOB_LINE = re.compile(r"^Job:\s*(.*)$", re.M)


def _committed() -> str:
    assert REPRODUCE_MD.is_file(), (
        f"{REPRODUCE_MD.relative_to(REPO)} does not exist (M10-14 deliverable 1, draft in 14a)"
    )
    return REPRODUCE_MD.read_text()


def _header(text: str) -> dict[str, str]:
    first = text.splitlines()[0]
    m = re.fullmatch(r"<!-- reproduce: (.*) -->", first)
    assert m, f"first line is not the runbook header: {first!r}"
    return dict(kv.split("=", 1) for kv in m.group(1).split())


def _table(text: str, heading: str) -> list[dict[str, str]]:
    m = re.search(rf"^## {heading}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    assert m, f"no `## {heading}` section"
    lines = [ln.strip() for ln in m.group(1).splitlines() if ln.strip().startswith("|")]
    assert len(lines) >= 2, f"`## {heading}` holds no table"
    cols = [c.strip() for c in lines[0].strip("|").split("|")]
    return [
        dict(zip(cols, (c.strip() for c in ln.strip("|").split("|")), strict=False))
        for ln in lines[2:]
    ]


def test_committed_runbook_header() -> None:
    hdr = _header(_committed())
    assert hdr.get("status") in {"draft", "final"}, hdr
    fc = hdr.get("freeze_commit", "")
    assert fc == "pending" or re.fullmatch(r"[0-9a-f]{40}", fc), hdr
    assert hdr.get("a3_branch") in {"cascadia", "fallback", "pending"}, hdr
    if hdr["status"] == "final":
        assert re.fullmatch(r"[0-9a-f]{40}", fc), "final runbook needs the 40-hex freeze commit"


def test_committed_runbook_passes_static(rc: Any) -> None:
    """Criterion 2: ``reproduce_check.py static`` on the committed runbook (as S-check does)."""
    text = _committed()
    assert _issues(rc, text, REPO) == []
    res = run_cli("static", str(REPRODUCE_MD))
    assert res.returncode == 0, res.stdout + res.stderr


def test_committed_runbook_lists_s1_to_s7_in_order() -> None:
    """Test plan: the draft lists S1-S7 (with their pass-1 evidence) and the freeze steps."""
    text = _committed()
    ids = re.findall(r"^### (S\d+|S-[a-z]+) ", text, re.M)
    numeric = [int(i[1:]) for i in ids if i[1:].isdigit()]
    assert numeric == sorted(numeric), f"steps out of order: {ids}"
    assert len(set(ids)) == len(ids), f"duplicate step ids: {ids}"
    for n in range(1, 8):
        assert f"S{n}" in ids, f"S{n} missing from {ids}"
    assert "S-check" in ids and any(i in ids for i in ("S8", "S15")), ids
    ev = {r["step"] for r in _table(text, "Evidence")}
    for n in (1, 2, 4, 5, 6):  # S3 runs in a dev job too; S7 may be re-converted at F
        assert f"S{n}" in ev, f"no pass-1 Evidence row for S{n}"


def test_committed_runbook_paths_exist() -> None:
    """Criterion 1 (a), with the test's own regex (a parser bug must not hide a missing file)."""
    text = _committed()
    chunks = FENCE.findall(text) + JOB_LINE.findall(text)
    missing = set()
    for chunk in chunks:
        for tok in PATH_TOKEN.findall(chunk):
            tok = tok.rstrip(".,:")
            if any(ch in tok for ch in "*?[{$"):
                continue
            if not (REPO / tok).exists():
                missing.add(tok)
    assert not missing, f"paths in reproduce.md that do not exist: {sorted(missing)}"


def test_committed_runbook_sbatch_calls() -> None:
    """Criterion 1 (b): every sbatch call names an existing scripts/slurm/*.sbatch; -J is giap-*."""
    text = _committed()
    calls = [
        ln
        for block in FENCE.findall(text)
        for ln in _logical_lines(block)
        if re.match(r"^\s*(\w+=\S*\s+)*sbatch\b", ln)
    ]
    assert calls, "no sbatch call in any bash block"
    for ln in calls:
        scripts = re.findall(r"(scripts/slurm/[\w.-]+\.sbatch)", ln)
        assert scripts, f"sbatch call names no scripts/slurm/*.sbatch: {ln}"
        assert (REPO / scripts[0]).is_file(), f"missing script {scripts[0]}"
        for name in re.findall(r"(?:-J|--job-name)[= ]\s*\"?([^\s\"]+)", ln):
            assert name.startswith("giap-"), f"job name {name!r} (R14: giap-*): {ln}"


def _logical_lines(block: str) -> list[str]:
    return block.replace("\\\n", " ").splitlines()


def test_committed_runbook_job_lines_match_headers() -> None:
    """Criterion 1 (c), own regex: Job fields equal the script's #SBATCH value unless overridden."""
    text = _committed()
    keys = {
        "partition": ("partition", "p"),
        "qos": ("qos", None),
        "cpus": ("cpus-per-task", "c"),
        "mem": ("mem", None),
        "time": ("time", "t"),
    }
    for m in re.finditer(r"^(### (\S+) .*?)(?=^### |^## |\Z)", text, re.M | re.S):
        body, sid = m.group(1), m.group(2)
        if not re.match(r"S(\d+|-[a-z]+)$", sid):
            continue
        jobs = JOB_LINE.findall(body)
        assert len(jobs) == 1, f"{sid}: exactly one Job: line, found {len(jobs)}"
        if jobs[0].startswith("none"):
            continue
        fields = dict(kv.split("=", 1) for kv in jobs[0].split() if "=" in kv)
        for f in (
            "script",
            "name",
            "partition",
            "qos",
            "cpus",
            "gpus",
            "mem",
            "time",
            "constraint",
            "array",
        ):
            assert f in fields, f"{sid}: Job line lacks {f}= ({jobs[0]})"
        assert fields["name"].startswith("giap-"), f"{sid}: name={fields['name']}"
        if fields["script"] == "-":
            continue  # a dev job without a script (srun): nothing to compare
        script = REPO / fields["script"]
        assert script.is_file(), f"{sid}: {fields['script']}"
        hdr = sbatch_header(script)
        cmds = " ".join(FENCE.findall(body))
        for jf, (long, short) in keys.items():
            pat = rf"--{long}[= ]" + (rf"|(?<!\S)-{short}\s" if short else "")
            if re.search(pat, cmds):
                continue  # overridden on the command line
            assert fields[jf] == hdr.get(long), (
                f"{sid}: Job {jf}={fields[jf]} but {fields['script']} says --{long}={hdr.get(long)}"
            )


def test_committed_runbook_inputs_well_formed() -> None:
    """Criterion 1 (d), own regex: committed files hash to their row; records contain the hash."""
    rows = _table(_committed(), "Inputs")
    assert rows, "empty Inputs table"
    assert len({r["id"] for r in rows}) == len(rows), "duplicate Inputs ids"
    for r in rows:
        assert re.fullmatch(r"[a-z0-9_]+", r["id"]), r
        assert re.fullmatch(r"[0-9a-f]{64}", r["sha256"]), r
        path = r["path"]
        assert path.startswith(("$GIAP_ROOT", "$GIAP_DATA_DIR", "$GIAP_DB_DIR")) or (
            not path.startswith("/") and (REPO / path).is_file()
        ), f"{r['id']}: path {path!r} is neither committed nor on FarmShare"
        recs = [x.strip() for x in r["recorded_in"].split(";") if x.strip()]
        assert recs, f"{r['id']}: empty recorded_in"
        if not path.startswith("$"):
            assert sha256_file(REPO / path) == r["sha256"], f"hash_mismatch:{r['id']}"
            assert any(x != path for x in recs), f"{r['id']}: a file is not its own record"
        for rec in recs:
            if rec.startswith("$"):
                continue
            assert (REPO / rec).is_file(), f"{r['id']}: recorded_in {rec} is not committed"
            assert r["sha256"] in (REPO / rec).read_text(errors="replace"), (
                f"not_recorded:{r['id']}:{rec}"
            )


def test_committed_runbook_has_no_stray_hashes() -> None:
    """Criterion 1 (e): every 64-hex string is an Inputs or Outputs row."""
    text = _committed()
    known = {r["sha256"] for r in _table(text, "Inputs")}
    known |= {r["sha256"] for r in _table(text, "Outputs")}
    stray = set(HEX64.findall(text)) - known
    assert not stray, f"64-hex strings that are no Inputs/Outputs row: {sorted(stray)}"


def test_committed_runbook_no_uv_or_pip() -> None:
    """Criterion 1 (g)."""
    for block in FENCE.findall(_committed()):
        assert not re.search(r"\buv\s+(run|sync)\b|\bpip3?\s|\bpython3? -m pip\b", block), block


def test_committed_runbook_reads_name_known_ids() -> None:
    text = _committed()
    ids = {r["id"] for r in _table(text, "Inputs")} | {r["id"] for r in _table(text, "Outputs")}
    for m in re.finditer(r"^Reads:\s*(.*)$", text, re.M):
        for rid in (x.strip() for x in m.group(1).split(",")):
            assert rid in ids or rid in {"", "-", "none"}, f"Reads names unknown id {rid!r}"


def test_committed_runbook_a3_sections() -> None:
    """Criterion 1 (f): one subsection per branch; the header names the single active one."""
    text = _committed()
    hdr = _header(text)
    found = dict(
        re.findall(r"^### A3-(cascadia|fallback) \(status: (active|inactive)\)", text, re.M)
    )
    assert set(found) == {"cascadia", "fallback"}, found
    active = [k for k, v in found.items() if v == "active"]
    if hdr["a3_branch"] == "pending":
        assert active == [], f"a3_branch=pending but {active} active"
    else:
        assert active == [hdr["a3_branch"]], (hdr, found)
    verdict_p = REPO / "validation" / "pilot" / "metrics" / "a3_feasibility.json"
    if verdict_p.is_file() and hdr["a3_branch"] == "cascadia":
        assert json.loads(verdict_p.read_text())["verdict"] == "pass"


def test_committed_runbook_final_mentions_the_r11_audit() -> None:
    """ "The audit's result (pass, at commit F) is one line in the runbook" (14b)."""
    text = _committed()
    if _header(text)["status"] != "final":
        pytest.skip("draft runbook: the R11 audit line is added with the final runbook (14b)")
    assert "test_m10_14_r11_audit" in text


def test_committed_runbook_names_the_freeze_scripts() -> None:
    """The draft lists S-check and the clean-up with the 14a scripts."""
    text = _committed()
    assert "scripts/slurm/pilot_reproduce_check.sbatch" in text
    assert "scripts/slurm/pilot_cleanup.sbatch" in text
    for name in ("giap-run-freeze-a2", "giap-run-freeze-a1", "giap-run-freeze-a2-rep"):
        assert name in text, name


# ---------------------------------------------------------------- pilot_reproduce_check.sbatch


def test_reproduce_check_sbatch_house_rules() -> None:
    check_house_rules(REPRO_SBATCH, name="giap-pilot-repro-check", cpus=2, mem_gb=8, minutes=60)


def test_reproduce_check_sbatch_runs_static_verify_and_reconversion() -> None:
    """S-check steps 1-3: static, verify into $GIAP_ROOT/pilot/reproduce/verify-<F7>.json, and
    the re-conversion with the M10-03 arguments through the venv python (no uv run)."""
    code = sbatch_code(REPRO_SBATCH)
    assert "reproduce_check.py" in code and re.search(r"\bstatic\b", code)
    assert re.search(r"\bverify\b", code) and "--out" in code
    assert "pilot/reproduce" in code and "verify-" in code
    assert "GIAP_PYTHON" in code and "giap.cli" in code
    assert re.search(r"convert\b.*--from\s+diann|--from diann", code), "giap convert --from diann"
    assert "--search-db" in code and "--abundance-out" in code
    assert "TMPDIR" in code, "re-converted tables go to $TMPDIR, never over the inputs"
    assert "identical" in code and "changed" in code


def test_reproduce_check_sbatch_bash32() -> None:
    import subprocess
    import sys

    assert REPRO_SBATCH.is_file(), "scripts/slurm/pilot_reproduce_check.sbatch does not exist"
    res = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "ci" / "check_bash32.py"), str(REPRO_SBATCH)],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr


def test_glob_tokens_are_excepted_by_the_own_regex() -> None:
    """Self-test of the path regex used on the committed runbook (globs are not files)."""
    toks = PATH_TOKEN.findall('sbatch "$GIAP_REPO/scripts/slurm/x.sbatch" validation/pilot/a*.cfg')
    assert toks == ["scripts/slurm/x.sbatch", "validation/pilot/a*.cfg"]
    assert fnmatch.fnmatch("validation/pilot/a1.cfg", toks[1])
