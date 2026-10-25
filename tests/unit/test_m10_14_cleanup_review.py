"""``pilot_cleanup.sbatch``, review round 1 of M10-14a (docs/specs/M10-14-pilot-reproducibility.md,
"Lead decisions (2026-10-24), review round 1": B1, S1-S4, S6, N1, N2, N5).

Same sandbox as ``test_m10_14_cleanup.py`` (``helpers.m10_14_kit.cleanup_sandbox``); some tests
swap in their own fake ``nextflow`` or ``rm``. No Slurm, no network. No statistical claim, so no
calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.m10_14_kit import (
    CLEANUP,
    F7,
    CleanupSandbox,
    F,
    cleanup_runbook,
    cleanup_sandbox,
    commit_runbook,
    sha256_file,
)
from helpers.slurm_sandbox import assert_honours_giap_repo, write_exe

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="bash script")

VALID_VERIFY = {
    "freeze_commit": F,
    "static_ok": True,
    "failures": [],
    "reconvert": {"a1": "identical", "a2": "identical"},
}


@pytest.fixture
def cs(tmp_path: Path) -> CleanupSandbox:
    assert CLEANUP.is_file(), "scripts/slurm/pilot_cleanup.sbatch does not exist (M10-14)"
    assert_honours_giap_repo(CLEANUP.read_text(), "pilot_cleanup.sbatch")
    return cleanup_sandbox(tmp_path)


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return f"exit {res.returncode}\n--- stdout\n{res.stdout}\n--- stderr\n{res.stderr}"


def _verify_json(cs: CleanupSandbox) -> Path:
    return cs.root / "pilot" / "reproduce" / f"verify-{F7}.json"


def _set_verify(cs: CleanupSandbox, doc: object) -> None:
    _verify_json(cs).write_text(json.dumps(doc))


def _without(key: str) -> dict[str, object]:
    return {k: v for k, v in VALID_VERIFY.items() if k != key}


def _refuses(cs: CleanupSandbox, **env: str) -> subprocess.CompletedProcess[str]:
    """Exit 2, no file under $GIAP_ROOT or outside it removed, no rm, no `nextflow clean -f`."""
    before = cs.snapshot()
    res = cs.run(CLEAN="1", **env)
    assert res.returncode == 2, _out(res)
    gone = set(before) - set(cs.snapshot())
    assert not gone, f"deleted despite exit 2: {sorted(gone)}"
    assert cs.rm_log() == [], cs.rm_log()
    assert not any("-f" in c["argv"] for c in cs.nextflow_calls()), cs.nextflow_calls()
    return res


def _plan_rows(cs: CleanupSandbox) -> list[list[str]]:
    p = cs.cleanup_dir / "plan.tsv"
    assert p.is_file(), f"no plan.tsv in {cs.cleanup_dir}"
    return [ln.split("\t") for ln in p.read_text().splitlines()[1:] if ln.strip()]


# ---------------------------------------------------------------- B1: canonical paths


def _root_variants(cs: CleanupSandbox) -> dict[str, str]:
    root = str(cs.root)
    alias = cs.sb.tmp / "root_alias"
    alias.symlink_to(cs.root)
    return {
        "trailing_slash": root + "/",
        "dot": f"{cs.root.parent}/./{cs.root.name}",
        "double_slash": f"{cs.root.parent}//{cs.root.name}",
        "dotdot": f"{root}/../{cs.root.name}",
        "symlink_alias": str(alias),
    }


VARIANTS = ["trailing_slash", "dot", "double_slash", "dotdot", "symlink_alias"]


@pytest.mark.parametrize("variant", VARIANTS)
def test_b1_protected_input_survives_a_non_canonical_giap_root(
    cs: CleanupSandbox, variant: str
) -> None:
    """Reproduced in the review: GIAP_ROOT=<root>/ made GIAP_DATA_DIR=<root>//data, the plan
    target <root>//data/.../mzml no longer matched the Inputs row, and the mzML was deleted."""
    mz = next(f for f in cs.delete_files if f.suffix == ".mzML")
    row = (
        f"| wt1_mzml | $GIAP_DATA_DIR/{mz.relative_to(cs.root / 'data')} | {sha256_file(mz)} | "
        "$GIAP_ROOT/pilot/A1/provenance.json |"
    )
    commit_runbook(cs, cleanup_runbook(extra_inputs=[row]))
    giap_root = _root_variants(cs)[variant]
    _refuses(cs, GIAP_ROOT=giap_root)
    assert mz.is_file(), f"protected Input deleted with GIAP_ROOT={giap_root!r}"


@pytest.mark.parametrize("variant", VARIANTS)
def test_b1_kept_output_survives_a_non_canonical_giap_root(
    cs: CleanupSandbox, variant: str
) -> None:
    work = cs.work_dirs[2]
    row = f"| pass1_work | $GIAP_ROOT/{work.relative_to(cs.root)} | {'a' * 64} | - | yes |"
    commit_runbook(cs, cleanup_runbook(extra_outputs=[row]))
    giap_root = _root_variants(cs)[variant]
    _refuses(cs, GIAP_ROOT=giap_root)
    assert work.is_dir(), f"kept Output deleted with GIAP_ROOT={giap_root!r}"


@pytest.mark.parametrize("variant", ["trailing_slash", "dot", "double_slash"])
def test_b1_non_canonical_giap_root_still_cleans(cs: CleanupSandbox, variant: str) -> None:
    """Canonicalising must not break the ordinary run: the same targets go, the rest stays."""
    giap_root = _root_variants(cs)[variant]
    before = cs.snapshot()
    res = cs.run(CLEAN="1", GIAP_ROOT=giap_root)
    assert res.returncode == 0, _out(res)
    gone = set(before) - set(cs.snapshot())
    assert gone == {str(p) for p in cs.delete_files}, sorted(gone)


# ---------------------------------------------------------------- S1, M12c: the verify JSON


@pytest.mark.parametrize(
    "doc",
    [
        pytest.param(_without("static_ok"), id="static_ok_missing"),
        pytest.param({**VALID_VERIFY, "static_ok": False}, id="static_ok_false"),
        pytest.param(_without("reconvert"), id="reconvert_missing"),
        pytest.param({**VALID_VERIFY, "reconvert": {}}, id="reconvert_empty"),
        pytest.param(
            {**VALID_VERIFY, "reconvert": {"a1": "failed", "a2": "identical"}},
            id="reconvert_failed",
        ),
    ],
)
def test_s1_refuses_unless_static_ok_and_reconvert_passed(
    cs: CleanupSandbox, doc: dict[str, object]
) -> None:
    _set_verify(cs, doc)
    _refuses(cs)


def test_s1_valid_verify_json_is_accepted(cs: CleanupSandbox) -> None:
    _set_verify(cs, {**VALID_VERIFY, "reconvert": {"a1": "changed", "a2": "absent"}})
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert not cs.work_dirs[0].exists()


# ---------------------------------------------------------------- N1: the freeze commit


def test_n1_freeze_commit_must_be_40_hex(cs: CleanupSandbox) -> None:
    """Mutant M14: a 7-hex F (verify-<F7>.json present and valid) must still refuse."""
    text = cleanup_runbook().replace(f"freeze_commit={F}", f"freeze_commit={F7}")
    assert f"freeze_commit={F7} " in text
    commit_runbook(cs, text)
    _set_verify(cs, {**VALID_VERIFY, "freeze_commit": F7})
    _refuses(cs)


def test_n1_freeze_commit_with_trailing_junk_refuses(cs: CleanupSandbox) -> None:
    text = cleanup_runbook().replace(f"freeze_commit={F}", f"freeze_commit={F}zz")
    commit_runbook(cs, text)
    _refuses(cs)


@pytest.mark.parametrize(
    "doc",
    [
        pytest.param({**VALID_VERIFY, "freeze_commit": F7 + "f" * 33}, id="other_commit"),
        pytest.param({**VALID_VERIFY, "freeze_commit": "pending"}, id="pending"),
        pytest.param(_without("freeze_commit"), id="missing"),
    ],
)
def test_n1_verify_json_freeze_commit_must_equal_f(
    cs: CleanupSandbox, doc: dict[str, object]
) -> None:
    _set_verify(cs, doc)
    _refuses(cs)


# ---------------------------------------------------------------- N2: uncaught errors exit 2


@pytest.mark.parametrize("which", ["verify", "evidence"])
def test_n2_malformed_json_exits_2(cs: CleanupSandbox, which: str) -> None:
    """A JSON that is not an object (``.get`` raises) is a precondition, not "a delete failed"."""
    if which == "verify":
        _verify_json(cs).write_text("[]")
    else:
        ev = cs.root / "runs" / "pilot" / "freeze" / F7 / "evidence" / "A1.json"
        ev.write_text('["not", "an", "object"]')
    _refuses(cs)


# ---------------------------------------------------------------- S2: nextflow clean

_NEXTFLOW_LEDGER = r"""#!{python}
'''Fake nextflow: logs argv, cwd and the DELETED.tsv rows naming cwd at call time; exits
$FAKE_NEXTFLOW_RC (default 0). Deletes nothing.'''
import json, os, sys
cwd = os.path.realpath(os.getcwd())
led = os.path.join(os.environ["GIAP_ROOT"], "pilot", "cleanup", "DELETED.tsv")
rows = []
if os.path.isfile(led):
    for ln in open(led).read().splitlines()[1:]:
        cells = ln.split("\t")
        if cells and os.path.realpath(cells[0]) == cwd:
            rows.append(cells)
with open(os.path.join(os.environ["FAKE_SLURM_DIR"], "nextflow.jsonl"), "a") as fh:
    fh.write(json.dumps({"argv": sys.argv[1:], "cwd": cwd, "ledger": rows}) + "\n")
sys.exit(int(os.environ.get("FAKE_NEXTFLOW_RC", "0")))
"""


def _ledger_nextflow(cs: CleanupSandbox) -> None:
    write_exe(cs.sb.bin / "nextflow", _NEXTFLOW_LEDGER.replace("{python}", sys.executable))


def test_s2_nextflow_clean_is_ledgered_before_it_runs(cs: CleanupSandbox) -> None:
    _ledger_nextflow(cs)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    calls = [c for c in cs.nextflow_calls() if c["argv"][:1] == ["clean"]]
    assert {c["cwd"] for c in calls} >= {os.path.realpath(d) for d in cs.run_dirs}, calls
    for c in calls:
        assert any("nextflow_clean" in "\t".join(r) for r in c["ledger"]), (
            f"`nextflow clean` ran in {c['cwd']} before its DELETED.tsv row: {c['ledger']}"
        )
    for d in cs.run_dirs:
        mine = [r for r in cs.ledger() if os.path.realpath(r["path"]) == os.path.realpath(d)]
        assert mine and all(r["job"] == "4242" for r in mine), mine


def test_s2_failing_nextflow_clean_exits_1(cs: CleanupSandbox) -> None:
    _ledger_nextflow(cs)
    res = cs.run(CLEAN="1", FAKE_NEXTFLOW_RC="1")
    assert res.returncode == 1, _out(res)
    for p in cs.keep_files:
        assert p.is_file(), p


# ---------------------------------------------------------------- S3: plan rows


@pytest.mark.parametrize("sep", ["\n", "\t"], ids=["newline", "tab"])
def test_s3_run_dir_name_with_tab_or_newline_refuses(cs: CleanupSandbox, sep: str) -> None:
    """A glob over a directory name with a newline gave the bash loop a relative ``rel/work``
    row (logged ``absent``); with a tab, a shifted row. Either way: exit 2 before any delete."""
    run = cs.root / "runs" / "pilot" / "pass1" / f"evil{sep}rel"
    (run / ".nextflow").mkdir(parents=True)
    (run / "work" / "ab").mkdir(parents=True)
    (run / "work" / "ab" / "task.out").write_text("task\n")
    _refuses(cs)
    assert (run / "work" / "ab" / "task.out").is_file()


def test_s3_giap_root_itself_is_never_a_target(tmp_path: Path) -> None:
    """Plan step 4's ``$GIAP_DATA_DIR/pxd078459/mzml`` reached through an ancestor link whose
    realpath *is* $ROOT: the old ``case "$real/" in "$ROOT"/*`` accepts it (``*`` matches empty)
    and ``rm -rf`` removes the whole of $GIAP_ROOT."""
    cs = cleanup_sandbox(tmp_path)
    real = cs.sb.tmp / "x" / "mzml"  # the real $GIAP_ROOT; $GIAP_ROOT stays the old path
    real.parent.mkdir()
    cs.root.rename(real)
    cs.root.symlink_to(real)
    pxd = real / "data" / "pxd078459"
    subprocess.run(["/bin/rm", "-rf", str(pxd)], check=True)
    pxd.symlink_to(real.parent)  # $ROOT/data/pxd078459/mzml -> <tmp>/x/mzml == $ROOT
    keep = real / "pilot" / "inputs" / "a1.tsv.gz"
    assert keep.is_file()
    res = cs.run(CLEAN="1")
    assert res.returncode == 2, _out(res)
    assert real.is_dir() and keep.is_file(), "$GIAP_ROOT itself was deleted"
    assert all(w.is_dir() for w in cs.work_dirs), "deleted despite exit 2"


# ---------------------------------------------------------------- S4: gate-digest Outputs rows


def test_s4_gate_digest_row_naming_the_digest_file_does_not_block(cs: CleanupSandbox) -> None:
    run = cs.run_dirs[0]  # $R/out_a2
    digest = run / "run_manifest.json"
    row = (
        f"| a2_gate | $GIAP_ROOT/{digest.relative_to(cs.root)} | {sha256_file(digest)} | - | yes |"
    )
    commit_runbook(cs, cleanup_runbook(extra_outputs=[row]))
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert digest.is_file() and not (run / "work").exists()


def test_s4_kept_directory_holding_a_target_refuses(cs: CleanupSandbox) -> None:
    """Mutant M9b: a ``kept = yes`` row naming the run directory holds ``$R/out_a2/work``."""
    run = cs.run_dirs[0]
    row = f"| a2_run | $GIAP_ROOT/{run.relative_to(cs.root)} | {'b' * 64} | - | yes |"
    commit_runbook(cs, cleanup_runbook(extra_outputs=[row]))
    _refuses(cs)
    assert (run / "work").is_dir()


# ---------------------------------------------------------------- S6: de novo outputs are kept


def test_s6_de_novo_outputs_are_kept_pending_k2(cs: CleanupSandbox) -> None:
    ssl = cs.root / "pilot" / "A3_cascadia" / "WT_1" / "WT_1.ssl"
    mztab = cs.root / "pilot" / "A3_casanovo" / "WT_1" / "WT_1.mztab"
    for p in (ssl, mztab):
        p.parent.mkdir(parents=True)
        p.write_text("de novo\n")
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert ssl.is_file() and mztab.is_file(), "de novo outputs deleted (lead decision S6: keep)"
    rows = {r[0]: r for r in _plan_rows(cs)}
    for p in (ssl, mztab):
        assert str(p) in rows, f"{p} not listed in plan.tsv"
        assert rows[str(p)][-1] == "kept_pending_K2", rows[str(p)]
    assert not any(r["path"] in (str(ssl), str(mztab)) for r in cs.ledger())


# ---------------------------------------------------------------- N5: defence in depth


def test_m4_exec_loop_skips_a_target_that_became_a_symlink(cs: CleanupSandbox) -> None:
    """The work/ is replaced by a link between the plan and its delete (here by `nextflow clean`):
    the exec loop's ``-L`` check leaves the link and ledgers ``skipped_symlink``."""
    outside = cs.outside / "late"
    outside.mkdir()
    (outside / "data.txt").write_text("someone else's\n")
    target = cs.run_dirs[0] / "work"
    fake = rf"""#!/bin/bash
if [ "$(pwd -P)" = "{os.path.realpath(cs.run_dirs[0])}" ]; then
  /bin/rm -rf work && ln -s "{outside}" work
fi
exit 0
"""
    write_exe(cs.sb.bin / "nextflow", fake)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert target.is_symlink(), "the link that replaced work/ was removed"
    assert (outside / "data.txt").is_file()
    rows = [r for r in cs.ledger() if r["path"] == str(target)]
    assert rows and rows[-1]["action"] == "skipped_symlink", rows
    assert not any(r["action"] in ("deleting", "deleted") for r in rows), rows


def test_m10_symlinked_run_dir_is_skipped(cs: CleanupSandbox) -> None:
    """A run directory that is a link (to a run outside $GIAP_ROOT) is not a run of the plan:
    no `nextflow clean` there, its work/ untouched, and the rest of the plan runs (exit 0)."""
    far = cs.outside / "far_run"
    (far / ".nextflow").mkdir(parents=True)
    (far / "work" / "ab").mkdir(parents=True)
    (far / "work" / "ab" / "task.out").write_text("far\n")
    link = cs.root / "runs" / "pilot" / "pass1" / "out_far"
    link.symlink_to(far)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert (far / "work" / "ab" / "task.out").is_file()
    cwds = {os.path.realpath(str(c["cwd"])) for c in cs.nextflow_calls()}
    assert os.path.realpath(far) not in cwds
    assert all(not w.exists() for w in cs.work_dirs)


def test_m15_step4_symlink_is_planned_as_skipped(cs: CleanupSandbox) -> None:
    """An mzml/ that links elsewhere inside $GIAP_ROOT is ``skipped_symlink`` in plan.tsv itself."""
    store = cs.root / "data" / "mzml_store"
    subprocess.run(["/bin/mv", str(cs.mzml_dir), str(store)], check=True)
    cs.mzml_dir.symlink_to(store)
    res = cs.run()
    assert res.returncode == 0, _out(res)
    rows = [r for r in _plan_rows(cs) if r[0] == str(cs.mzml_dir)]
    assert rows and rows[0][-1] == "skipped_symlink", rows
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    assert cs.mzml_dir.is_symlink() and (store / "WT_1.mzML").is_file()


def test_m16_no_sha256_above_10_gb(cs: CleanupSandbox) -> None:
    """A regular-file target over 10 GB (sparse here) is ledgered with ``sha256 = NA``."""
    subprocess.run(["/bin/rm", "-rf", str(cs.mzml_dir)], check=True)
    with open(cs.mzml_dir, "wb") as fh:
        fh.truncate(10 * 1024**3 + 1)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    rows = [r for r in cs.ledger() if r["path"] == str(cs.mzml_dir)]
    assert rows and rows[-1]["action"] == "deleted", rows
    assert rows[-1]["bytes"] == str(10 * 1024**3 + 1), rows
    assert all(r["sha256"] == "NA" for r in rows), rows


def test_m16_small_regular_file_is_hashed(cs: CleanupSandbox) -> None:
    subprocess.run(["/bin/rm", "-rf", str(cs.mzml_dir)], check=True)
    cs.mzml_dir.write_bytes(b"one small mzML\n")
    want = sha256_file(cs.mzml_dir)
    res = cs.run(CLEAN="1")
    assert res.returncode == 0, _out(res)
    rows = [r for r in cs.ledger() if r["path"] == str(cs.mzml_dir)]
    assert rows and rows[-1]["sha256"] == want, rows


_RM_FAIL = r"""#!/bin/bash
for a in "$@"; do
  case "$a" in
    *"$FAKE_RM_FAIL"*) echo "rm: cannot remove '$a': fake failure" >&2; exit 1 ;;
  esac
done
exec /bin/rm "$@"
"""


def test_m17_rm_failure_exits_1(cs: CleanupSandbox) -> None:
    write_exe(cs.sb.bin / "rm", _RM_FAIL)
    victim = cs.work_dirs[1]
    res = cs.run(CLEAN="1", FAKE_RM_FAIL=str(victim))
    assert res.returncode == 1, _out(res)
    rows = [r for r in cs.ledger() if r["path"] == str(victim)]
    assert rows and rows[-1]["action"] == "failed", rows
    assert victim.is_dir()
    others = [w for w in cs.work_dirs if w != victim]
    assert all(not w.exists() for w in others), "one failed rm stops the rest"
