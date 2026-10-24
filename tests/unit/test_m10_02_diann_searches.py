"""``scripts/slurm/diann_pilot.sbatch`` after M10-02: verified copy, acquisition gate, run list,
step-1 marker, fresh step 2, provenance on every exit (docs/specs/M10-02-diann-searches.md).

Spec: "The searches", "provenance.json", examples 4-9, edge cases "Image missing", "Not
diaPASEF", "Timeout", "Partial report", "Runs", acceptance criteria 2-5 and 11.

Everything runs in ``helpers.pilot_kit_m10.m10_sandbox`` with the extended fake ``apptainer``:
nothing executes DIA-NN or a container, and no real data is read. ``DIANN_SIF`` points at a
read-only stand-in whose sha256 is the sandbox pin.

Conventions settled by the test author (where the spec leaves the shape open):

* ``provenance.json`` ``step1``/``step2`` are either an object with ``status`` (or ``state``) and
  ``exit_code``, or a string with ``step<N>_exit_code`` next to it (``helpers.pilot_kit_m10.
  step_state``); ``commands`` is a list of the two argv arrays or an object keyed by step.
* A recorded argv (``commands`` and ``command_step<N>.txt``, one token per line) **ends with** the
  DIA-NN arguments the fake container received; it may start with the ``apptainer exec …``
  prefix and the binary.
* A failure of step 1 (DIA-NN non-zero) is exit 1 like a failure of step 2, and leaves no
  ``step1.done``.
* The pilot repeats the install job's banner check: a banner without ``DIA-NN 2.3.1 Academia`` is
  exit 2 before any search.

No statistical claim is made here, so there is no calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_kit_m10 import (
    BANNER,
    COPY_REL,
    M10_SCRIPTS,
    RUNS,
    STAND_IN,
    STAND_IN_SHA256,
    M10Sandbox,
    acquisition_doc,
    commands,
    diann_args,
    ends_with,
    flag_values,
    m10_sandbox,
    sha256,
    step_state,
    strings,
    write_json,
)

ARMS = ("A1", "A2")
SCRIPT = M10_SCRIPTS["diann_pilot"]


@pytest.fixture
def m10(tmp_path: Path) -> M10Sandbox:
    return m10_sandbox(tmp_path)


def _run(m10: M10Sandbox, arm: str = "A1", **env: str):  # type: ignore[no-untyped-def]
    return m10.run("diann_pilot", arm, **env)


def _ok(res) -> None:  # type: ignore[no-untyped-def]
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _rc(res, want: int) -> str:  # type: ignore[no-untyped-def]
    msg = res.stdout + res.stderr
    assert res.returncode == want, f"exit {res.returncode}, want {want}\n{msg}"
    return msg


def _acq(m10: M10Sandbox, doc: Any, name: str = "acq.json") -> dict[str, str]:
    path = write_json(m10.sb.tmp / name, doc)
    return {"ACQUISITION_JSON": str(path)}


def _speclib(m10: M10Sandbox, arm: str = "A1") -> Path:
    return m10.out(arm) / f"{arm}.predicted.speclib"


# ------------------------------------------------------------ example 4: a full search


@pytest.mark.parametrize("arm", ARMS)
def test_example4_two_calls_on_the_copy_and_the_outputs(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm))
    searches = m10.searches()
    assert len(searches) == 2, f"exactly two DIA-NN calls, got {len(searches)}"
    assert {c["image"] for c in m10.calls() if c["image"]} == {str(m10.stand_in)}
    out = m10.out(arm)
    for name in ("step1.done", "outputs.sha256", "provenance.json", "diann_version.txt",
                 "command_step1.txt", "command_step2.txt", f"{arm}_report.parquet"):  # fmt: skip
        assert (out / name).is_file(), f"{name} missing under {out}"
    prov = m10.provenance(arm)
    assert prov["arm"] == arm
    assert prov["status"] == "ok" and int(prov["exit_code"]) == 0, prov["status"]
    assert step_state(prov, "step1") == ("ran", 0)
    assert step_state(prov, "step2") == ("ran", 0)
    assert not (out / "FAILED").exists()


def test_example4_step1_done_records_the_speclib_hash(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    marker = (m10.out("A1") / "step1.done").read_text().strip()
    m = re.fullmatch(r"([0-9a-f]{64})\s+\*?(\S+)", marker)
    assert m, f"step1.done is not '<sha256>  A1.predicted.speclib': {marker!r}"
    assert Path(m.group(2)).name == "A1.predicted.speclib"
    assert m.group(1) == sha256(_speclib(m10))


def test_example4_outputs_sha256_lists_what_dia_nn_wrote(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    out = m10.out("A1")
    lines = (out / "outputs.sha256").read_text().splitlines()
    listed: dict[str, str] = {}
    for line in lines:
        m = re.fullmatch(r"([0-9a-f]{64})\s+\*?(\S+)", line.strip())
        assert m, f"unparsable outputs.sha256 line {line!r}"
        listed[m.group(2)] = m.group(1)
    names = [Path(n).name for n in listed]
    assert names == sorted(names), "outputs.sha256 is sorted by name"
    assert "A1_report.parquet" in names and "A1_empirical.parquet" in names, names
    for n, h in listed.items():
        p = Path(n) if Path(n).is_absolute() else out / n
        assert p.is_file() and sha256(p) == h, n
    excluded = {"provenance.json", "outputs.sha256", "FAILED"}
    assert not excluded & set(names), names
    assert not [n for n in listed if "fasta/" in n], "fasta/ is excluded"


def test_example4_diann_version_txt_holds_the_banner(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    assert BANNER in (m10.out("A1") / "diann_version.txt").read_text()


def test_the_repository_stays_clean(m10: M10Sandbox) -> None:
    _ok(_run(m10, "A2"))
    assert m10.repo_is_clean() == "", m10.repo_is_clean()


# ------------------------------------------------------------ rule 1: only the verified copy


def test_the_default_image_is_the_copy_under_giap_root(m10: M10Sandbox) -> None:
    copy = m10.root / COPY_REL
    copy.parent.mkdir(parents=True)
    copy.write_bytes(STAND_IN)
    copy.chmod(0o444)
    _ok(m10.run("diann_pilot", "A1", unset=("DIANN_SIF",)))
    images = {c["image"] for c in m10.calls() if c["image"]}
    assert images == {str(copy)}, images
    assert str(m10.shared) not in images, "DIANN_SIF_SHARED is never run"


def test_a_missing_copy_is_exit_2_with_no_container_call(m10: M10Sandbox) -> None:
    res = m10.run("diann_pilot", "A1", unset=("DIANN_SIF",))
    msg = _rc(res, 2)
    assert str(m10.root / COPY_REL) in msg, msg
    assert "DIANN_SIF" in msg and "diann_image.sbatch" in msg, msg
    assert not m10.calls(), "a container call (even --version) before the image checks"
    assert m10.provenance("A1")["status"] == "precondition"


def test_example6_an_altered_copy_is_exit_2_before_any_container_call(m10: M10Sandbox) -> None:
    altered = m10.sb.tmp / "images" / "altered.sif"
    altered.write_bytes(b"altered image\n")
    altered.chmod(0o444)
    msg = _rc(_run(m10, DIANN_SIF=str(altered)), 2)
    assert sha256(altered) in msg and STAND_IN_SHA256 in msg, f"both hashes:\n{msg}"
    assert not m10.calls(), "a container call (even --version) on an unverified image"
    prov = m10.provenance("A1")
    assert prov["status"] == "precondition" and int(prov["exit_code"]) == 2
    assert step_state(prov, "step1")[0] == "not_reached"
    assert step_state(prov, "step2")[0] == "not_reached"


@pytest.mark.parametrize("mode", [0o644, 0o664, 0o446, 0o744])
def test_example7_a_writable_copy_is_exit_2(m10: M10Sandbox, mode: int) -> None:
    m10.stand_in.chmod(mode)
    msg = _rc(_run(m10), 2)
    assert "read-only" in msg and "0444" in msg, msg
    assert not m10.calls()


def test_a_malformed_or_missing_pin_is_exit_2(m10: M10Sandbox) -> None:
    m10.pin_file.write_text("not a pin\n")
    _rc(_run(m10), 2)
    assert not m10.calls()
    m10.pin_file.unlink()
    _rc(_run(m10), 2)
    assert not m10.calls()


def test_the_banner_check_is_repeated(m10: M10Sandbox) -> None:
    _rc(_run(m10, FAKE_DIANN_BANNER="DIA-NN 2.7.0"), 2)
    assert not m10.searches(), "a search ran on an image with the wrong banner"


def test_the_script_never_reads_the_shared_image() -> None:
    code = "\n".join(ln for ln in SCRIPT.read_text().splitlines()
                     if not ln.lstrip().startswith("#"))  # fmt: skip
    assert "DIANN_SIF_SHARED" not in code
    assert "/jordan/" not in code, "the shared path belongs to diann_image.sbatch only"
    assert "validation/pilot/diann_image.sha256" in code
    m = re.search(r"DIANN_SIF:-([^}]+)\}", code)
    assert m and m.group(1).strip("\"'") in {
        "$GIAP_ROOT/images/diann_2.3.1.sif",
        "${GIAP_ROOT}/images/diann_2.3.1.sif",
    }, m and m.group(1)


# ------------------------------------------------------------ rule 2: the acquisition gate


@pytest.mark.parametrize("mode", ["ddaPASEF", "mixed", "ms1_only", "unknown", "inconsistent"])
def test_example8_a_mode_other_than_diapasef_is_exit_3(m10: M10Sandbox, mode: str) -> None:
    msg = _rc(_run(m10, **_acq(m10, acquisition_doc(mode=mode))), 3)
    assert "§17" in msg and "ADR-0051" in msg and mode in msg, msg
    assert not m10.calls(), "no container call when the acquisition mode forbids the settings"
    prov = m10.provenance("A1")
    assert prov["status"] == "acquisition_mode" and int(prov["exit_code"]) == 3


@pytest.mark.parametrize("variant", ["synchro-PASEF", "slice-PASEF", "midia-PASEF"])
def test_a_diapasef_variant_in_notes_is_exit_3(m10: M10Sandbox, variant: str) -> None:
    doc = acquisition_doc(notes={"wt_rep2": f"{variant} window scheme, 12 groups"})
    _rc(_run(m10, **_acq(m10, doc)), 3)
    assert not m10.calls()


def test_plain_notes_do_not_stop_the_search(m10: M10Sandbox) -> None:
    doc = acquisition_doc(notes={"wt_rep1": "24 window groups, m/z 400-1200"})
    _ok(_run(m10, **_acq(m10, doc)))


def test_the_committed_acquisition_file_is_the_default(m10: M10Sandbox) -> None:
    m10.acquisition.write_text(json.dumps(acquisition_doc(mode="ddaPASEF")))
    _rc(_run(m10), 3)


@pytest.mark.parametrize("kind", ["missing", "unparsable", "no-runs-key"])
def test_a_missing_or_unparsable_acquisition_file_is_exit_2(m10: M10Sandbox, kind: str) -> None:
    path = m10.sb.tmp / "acq.json"
    if kind == "unparsable":
        path.write_text("{not json")
    elif kind == "no-runs-key":
        path.write_text(json.dumps({"dataset_mode": "diaPASEF"}))
    _rc(_run(m10, ACQUISITION_JSON=str(path)), 2)
    assert not m10.calls()
    assert m10.provenance("A1")["status"] == "precondition"


# ------------------------------------------------------------ rule 3: the runs


def test_the_runs_are_the_acquisition_runs_in_file_order(m10: M10Sandbox) -> None:
    order = ["wt_rep3", "wt_rep1", "wt_rep2"]
    _ok(_run(m10, **_acq(m10, acquisition_doc(order))))
    step2 = diann_args(m10.searches()[1])
    assert flag_values(step2, "--f") == [str(m10.raw_dir / f"{r}.d") for r in order]


def test_an_unlisted_run_directory_is_ignored_with_a_warning(m10: M10Sandbox) -> None:
    old = m10.raw_dir / "old.d"
    old.mkdir()
    (old / "analysis.tdf").write_text("x")
    (old / "analysis.tdf_bin").write_text("x")
    res = _run(m10)
    _ok(res)
    assert "old.d" in res.stdout + res.stderr, "the ignored run must be named in a warning"
    step2 = diann_args(m10.searches()[1])
    assert str(old) not in flag_values(step2, "--f")
    assert flag_values(step2, "--f") == [str(m10.raw_dir / f"{r}.d") for r in RUNS]


@pytest.mark.parametrize("missing", ["analysis.tdf_bin", "analysis.tdf"])
def test_a_listed_run_without_its_tdf_files_is_exit_2(m10: M10Sandbox, missing: str) -> None:
    (m10.raw_dir / "wt_rep2.d" / missing).unlink()
    msg = _rc(_run(m10), 2)
    assert "wt_rep2" in msg, msg
    assert not m10.searches()


def test_a_listed_run_with_no_directory_is_exit_2(m10: M10Sandbox) -> None:
    doc = acquisition_doc([*RUNS, "wt_rep4"])
    msg = _rc(_run(m10, **_acq(m10, doc)), 2)
    assert "wt_rep4" in msg, msg
    assert not m10.searches()


def test_an_empty_run_list_is_exit_2(m10: M10Sandbox) -> None:
    _rc(_run(m10, **_acq(m10, acquisition_doc([]))), 2)
    assert not m10.searches()


# ------------------------------------------------------------ rule 4: the step-1 marker


def test_example5_a_resubmission_skips_step1_and_starts_step2_fresh(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    m10.reset_calls()
    report = m10.out("A1") / "A1_report.parquet"
    assert report.is_file()
    _ok(_run(m10))
    searches = m10.searches()
    assert len(searches) == 1, f"step 1 must be skipped; got {len(searches)} calls"
    assert flag_values(diann_args(searches[0]), "--lib") == [str(_speclib(m10))]
    assert searches[0]["preexisting"].get(str(report)) is False, (
        "the old report must be deleted before step 2 runs"
    )
    prov = m10.provenance("A1")
    assert step_state(prov, "step1")[0] == "skipped"
    assert step_state(prov, "step2") == ("ran", 0)
    assert prov["status"] == "ok"


def test_a_speclib_without_a_marker_is_rebuilt(m10: M10Sandbox) -> None:
    """Example 9, second half: killed inside step 1, so the library may be truncated."""
    out = m10.out("A1")
    out.mkdir(parents=True)
    _speclib(m10).write_text("half a library")
    (out / "A1.parquet").write_text("half")
    _ok(_run(m10))
    searches = m10.searches()
    assert len(searches) == 2, "step 1 must run again without a marker"
    pre = searches[0]["preexisting"]
    assert pre.get(str(_speclib(m10))) is False, "the leftover speclib is deleted before step 1"
    assert pre.get(str(out / "A1.parquet")) is False, "the leftover step-1 parquet is deleted"
    assert _speclib(m10).read_text() != "half a library"


def test_a_marker_whose_hash_disagrees_rebuilds_the_library(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    _speclib(m10).write_text("truncated by a timeout")
    m10.reset_calls()
    _ok(_run(m10))
    assert len(m10.searches()) == 2, "a stale marker must not let step 1 be skipped"
    marker = (m10.out("A1") / "step1.done").read_text().split()[0]
    assert marker == sha256(_speclib(m10))


def test_a_failed_step1_writes_no_marker(m10: M10Sandbox) -> None:
    _rc(_run(m10, FAKE_DIANN_RC_STEP1="1"), 1)
    out = m10.out("A1")
    assert not (out / "step1.done").exists()
    assert (out / "FAILED").is_file()
    prov = m10.provenance("A1")
    assert prov["status"] == "failed"
    assert step_state(prov, "step1") == ("ran", 1)
    assert step_state(prov, "step2")[0] == "not_reached"


# ------------------------------------------------------------ rule 5: step 2 is always fresh


def test_example9_a_partial_report_and_leftovers_are_removed_before_step2(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    out = m10.out("A1")
    leftovers = [out / "A1_report.parquet", out / "A1_report.stats.tsv",
                 out / "A1_empirical.parquet", out / "FAILED"]  # fmt: skip
    for p in leftovers:
        p.write_text("40 % written")
    (out / "outputs.sha256").write_text("stale\n")
    m10.reset_calls()
    _ok(_run(m10))
    (step2,) = m10.searches()
    pre = step2["preexisting"]
    assert pre.get(str(out / "A1_report.parquet")) is False
    assert pre.get(str(out / "A1_empirical.parquet")) is False
    assert not (out / "FAILED").exists(), "FAILED is removed when a new attempt starts"
    assert (out / "A1_report.stats.tsv").read_text() != "40 % written"
    assert "stale" not in (out / "outputs.sha256").read_text()


# ------------------------------------------------------------ rule 6: failures (criterion 4)


@pytest.mark.parametrize(
    "env", [{"FAKE_DIANN_RC_STEP2": "1"}, {"FAKE_DIANN_NO_REPORT": "1"}], ids=["rc1", "no-report"]
)
def test_a_failed_search_is_exit_1_with_failed_and_no_sums(
    m10: M10Sandbox, env: dict[str, str]
) -> None:
    _rc(_run(m10, **env), 1)
    out = m10.out("A1")
    failed = out / "FAILED"
    assert failed.is_file(), "FAILED must be written"
    text = failed.read_text()
    assert "step2" in text.replace(" ", "").lower() or "step 2" in text.lower(), text
    assert re.search(r"\d{4}-\d{2}-\d{2}T?", text), f"FAILED records the UTC time: {text!r}"
    assert not (out / "outputs.sha256").exists()
    prov = m10.provenance("A1")
    assert prov["status"] == "failed" and int(prov["exit_code"]) == 1
    assert step_state(prov, "step1") == ("ran", 0)
    state, code = step_state(prov, "step2")
    assert state == "ran"
    if "FAKE_DIANN_RC_STEP2" in env:
        assert code == 1
    assert prov["outputs"] is None


def test_an_empty_report_is_a_failure(m10: M10Sandbox) -> None:
    """A zero-byte ``<arm>_report.parquet`` counts as no report (rule 6)."""
    fake = m10.sb.bin / "apptainer"
    fake.write_text(fake.read_text().replace(
        '    touch(out)\n', '    touch(out); open(out, "w").close()\n'))  # fmt: skip
    _rc(_run(m10), 1)
    assert (m10.out("A1") / "FAILED").is_file()


# ------------------------------------------------------------ rule 7: provenance (criterion 5)

JOB_KEYS = ("id", "name", "partition", "cpus", "mem", "node", "cpu_model",
            "start_utc", "end_utc", "elapsed_s")  # fmt: skip
TOP_KEYS = ("arm", "status", "exit_code", "step1", "step2", "job", "image", "apptainer", "cfg",
            "repo", "inputs", "commands", "outputs")  # fmt: skip


def _check_schema(prov: dict[str, Any]) -> None:
    missing = [k for k in TOP_KEYS if k not in prov]
    assert not missing, f"provenance.json lacks {missing}"
    assert prov["status"] in {"ok", "failed", "precondition", "acquisition_mode"}
    assert isinstance(prov["exit_code"], int)
    job = prov["job"]
    assert not [k for k in JOB_KEYS if k not in job], f"job lacks {set(JOB_KEYS) - set(job)}"
    for k in ("path", "sha256", "pin", "mode", "banner"):
        assert k in prov["image"], f"image.{k}"
    for k in ("path", "sha256"):
        assert k in prov["cfg"], f"cfg.{k}"
    assert "version" in prov["apptainer"] or isinstance(prov["apptainer"], str)


@pytest.mark.parametrize("arm", ARMS)
def test_provenance_of_a_successful_search(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm, SLURM_MEM_PER_NODE="65536", SLURM_JOB_NAME=f"giap-diann-{arm}",
             SLURM_JOB_PARTITION="normal"))  # fmt: skip
    prov = m10.provenance(arm)
    _check_schema(prov)
    job = prov["job"]
    assert str(job["id"]) == "4242" and str(job["cpus"]) == "16"
    assert str(job["mem"]) == "65536" and job["name"] == f"giap-diann-{arm}"
    assert job["node"] == "barley-01", "the hostname"
    assert float(job["elapsed_s"]) >= 0
    img = prov["image"]
    assert img["path"] == str(m10.stand_in)
    assert img["sha256"] == STAND_IN_SHA256 and STAND_IN_SHA256 in json.dumps(img["pin"])
    assert int(str(img["mode"]), 8) == 0o444
    assert BANNER in json.dumps(img["banner"])
    cfg = m10.repo / "validation" / "pilot" / f"diann_{arm}.cfg"
    assert prov["cfg"]["sha256"] == sha256(cfg)
    repo_strings = strings(prov["repo"])
    assert m10.sb.head in repo_strings, "git rev-parse HEAD of GIAP_REPO"
    assert any(isinstance(v, bool) for v in prov["repo"].values()), "the dirty flag"
    assert prov["outputs"] and str(prov["outputs"]).endswith("outputs.sha256")


@pytest.mark.parametrize("arm", ARMS)
def test_provenance_commands_equal_the_logged_argv(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm))
    logged = [diann_args(c) for c in m10.searches()]
    recorded = commands(m10.provenance(arm))
    assert len(recorded) == 2, recorded
    for rec, log in zip(recorded, logged, strict=True):
        assert ends_with(rec, log), f"provenance commands differ from the call:\n{rec}\n{log}"
    for n, log in enumerate(logged, start=1):
        tokens = (m10.out(arm) / f"command_step{n}.txt").read_text().splitlines()
        assert ends_with(tokens, log), f"command_step{n}.txt differs from the call"
    assert "K*,R*" in logged[0], "the --cut token reaches DIA-NN unglobbed"


@pytest.mark.parametrize("arm", ARMS)
def test_provenance_input_hashes_equal_the_sandbox_files(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm))
    prov = m10.provenance(arm)
    hashes = set(strings(prov["inputs"]))
    fastas = {f for c in m10.searches() for f in flag_values(diann_args(c), "--fasta")}
    for f in fastas:
        assert f in hashes and sha256(Path(f)) in hashes, f"--fasta {f} path/sha256 not recorded"
    if arm == "A2":
        gz = m10.repo / "validation" / "pilot" / "db_a2" / "strain_db.faa.gz"
        assert sha256(gz) in hashes, "A2 records the committed .gz too"
    acq = prov["inputs"]["acquisition_json"]
    assert acq["sha256"] == sha256(m10.acquisition) and acq["dataset_mode"] == "diaPASEF"
    for r in RUNS:
        tdf = m10.raw_dir / f"{r}.d" / "analysis.tdf"
        assert sha256(tdf) in hashes, f"analysis.tdf of {r}"
        assert r in hashes, f"run label {r}"
    bins = json.dumps(prov["inputs"])
    assert sha256(m10.raw_dir / "wt_rep1.d" / "analysis.tdf_bin") not in bins, (
        "analysis.tdf_bin is sized, not hashed"
    )


def test_provenance_is_written_on_a_failure_too(m10: M10Sandbox) -> None:
    _rc(_run(m10, FAKE_DIANN_RC_STEP2="3"), 1)
    prov = m10.provenance("A1")
    _check_schema(prov)
    assert prov["status"] == "failed"
    assert step_state(prov, "step2") == ("ran", 3)


def test_a_usage_error_is_exit_2_before_anything(m10: M10Sandbox) -> None:
    _rc(m10.run("diann_pilot", "A3"), 2)
    assert not m10.calls()


# ------------------------------------------------------------ house rules (criterion 11)


def _directive(key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", SCRIPT.read_text(), re.M)
    return (m.group(1) or "") if m else None


def test_the_jobs_table_resources() -> None:
    assert (_directive("job-name") or "").startswith("giap-diann")
    assert _directive("partition") == "normal" and _directive("qos") == "normal"
    assert _directive("cpus-per-task") == "16"
    assert (_directive("mem") or "").upper().rstrip("B") in {"64G", "65536M"}
    assert _directive("time") == "24:00:00"


def test_step2_keeps_its_quant_files_under_the_arm(m10: M10Sandbox) -> None:
    """DIA-NN writes ``<run>.d.quant`` next to the raw runs unless ``--temp`` is given, so two
    arms searching the same runs at once read each other's (Slurm 1879914, "library size
    mismatch"): step 2 passes ``--temp <out>/quant``."""
    _ok(_run(m10, "A1"))
    step2 = m10.searches()[-1]["inner"]
    assert "--temp" in step2, step2
    assert step2[step2.index("--temp") + 1] == str(m10.out("A1") / "quant")
