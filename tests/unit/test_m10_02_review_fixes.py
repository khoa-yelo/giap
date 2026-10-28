"""The fixes the PR #293 review asked for (M10-02, the DIA-NN A1/A2 searches).

One section per finding, with the lead's decisions:

* B1: a TERM, INT or HUP during a DIA-NN call is a failure: exit 143/130/129, ``provenance.json``
  ``status: failed`` with that exit code, ``STEP<N>_RC`` never preset to 0, and ``FAILED`` written
  with the code. The job is signalled as Slurm does (``scancel``, TIMEOUT): the whole process
  group, while the fake DIA-NN sleeps.
* S1: ``--mem=64000M`` (16 x FarmShare's MaxMemPerCPU of 4000 MB), so a task gets exactly 16
  CPUs; no M10 job asks for more memory than its CPUs carry, and the searches pass
  ``--threads 16``.
* S2: ``step1.done`` holds, after its first line ``<speclib sha256>  <arm>.predicted.speclib``,
  the job ID, the commit and the sha256s of the cfg, the FASTA files and the image; a skipped
  step 1 copies the job ID and commit into ``provenance.json`` ``step1``.
* S3: ``diann_image.sbatch`` checks the banner before the copy is renamed into place: a failed
  banner installs nothing, and a rerun checks it again.
* S4: a marker whose hashes differ from the current inputs makes step 1 rerun; the unpacked A2
  FASTA is rebuilt when its sidecar ``<fasta>.src.sha256`` (the sha256 of the ``.gz`` it came
  from) is missing or differs, and is unpacked to a temporary file and renamed (tested with a
  ``gzip``/``gunzip`` on PATH that fails partway, so the unpack must go through one of them).
* N1: ``FAILED`` never records exit code 0 (a step that exits 0 but writes no output);
  N2: a broken ``GIPPA_PYTHON`` is named; N5: a step-1 rerun removes the old ``<arm>.log.txt``.
* Test gaps: A1 and A2 get different ``--temp`` directories and nothing appears under
  ``RAW_DIR``; a stale ``quant/`` is emptied before step 2; ``outputs.sha256`` leaves out the
  job's own records (``command_step*.txt`` and the rest).

Everything runs in ``helpers.pilot_kit_m10.m10_sandbox``; nothing executes DIA-NN or a container.
No statistical claim is made here, so there is no calibration test (spec §8.4).
"""

from __future__ import annotations

import gzip
import json
import os
import re
import signal
from pathlib import Path

import pytest
from helpers.pilot_kit_m10 import (
    M10_SCRIPTS,
    STAND_IN_SHA256,
    M10Sandbox,
    commands,
    diann_args,
    flag_values,
    m10_sandbox,
    sha256,
    step_state,
)

ARMS = ("A1", "A2")
MAX_MEM_PER_CPU_MB = 4000  # FarmShare partition normal (sacct/scontrol, review S1)


@pytest.fixture
def m10(tmp_path: Path) -> M10Sandbox:
    return m10_sandbox(tmp_path)


@pytest.fixture
def m10i(tmp_path: Path) -> M10Sandbox:
    """The install tests' sandbox: DIANN_SIF unset, so the copy goes to its default."""
    m = m10_sandbox(tmp_path)
    m.sb.env.pop("DIANN_SIF", None)
    return m


def _run(m10: M10Sandbox, arm: str = "A1", **env: str):  # type: ignore[no-untyped-def]
    return m10.run("diann_pilot", arm, **env)


def _ok(res) -> None:  # type: ignore[no-untyped-def]
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _rc(res, want: int) -> str:  # type: ignore[no-untyped-def]
    msg = res.stdout + res.stderr
    assert res.returncode == want, f"exit {res.returncode}, want {want}\n{msg}"
    return msg


def _failed_code(text: str) -> int:
    m = re.search(r"exit_code[=:\s]+(-?\d+)", text)
    assert m, f"FAILED records no exit_code: {text!r}"
    return int(m.group(1))


def _marker(m10: M10Sandbox, arm: str = "A1") -> str:
    p = m10.out(arm) / "step1.done"
    assert p.is_file(), f"{p} missing"
    return p.read_text()


def _step1_fastas(m10: M10Sandbox) -> list[str]:
    return flag_values(diann_args(m10.searches()[0]), "--fasta")


# ------------------------------------------------------------ B1: a signal is a failure

SIGNALS = [(signal.SIGTERM, 143), (signal.SIGINT, 130), (signal.SIGHUP, 129)]


def _signal_job(m10: M10Sandbox, sig: int, step: str):  # type: ignore[no-untyped-def]
    """Start the A1 search, wait until DIA-NN sleeps in ``step``, signal the job's process group
    (as Slurm does on scancel and TIMEOUT), and return the finished process and its output."""
    proc = m10.popen("diann_pilot", "A1", **{f"FAKE_DIANN_SLEEP_{step.upper()}": "120"})
    try:
        assert m10.wait_for_sleep(proc) == step
        os.killpg(proc.pid, sig)
        out, err = proc.communicate(timeout=120)
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.communicate()
    return proc, out + err


@pytest.mark.parametrize("sig, code", SIGNALS, ids=["TERM", "INT", "HUP"])
def test_b1_a_signal_during_step2_is_a_failure(m10: M10Sandbox, sig: int, code: int) -> None:
    proc, msg = _signal_job(m10, sig, "step2")
    assert proc.returncode in {code, -sig}, f"exit {proc.returncode}, want {code}\n{msg}"
    prov = m10.provenance("A1")
    assert prov["status"] == "failed", f"a {sig!r} kill is recorded as {prov['status']!r}"
    assert int(prov["exit_code"]) == code, prov["exit_code"]
    assert step_state(prov, "step1") == ("ran", 0)
    state, rc = step_state(prov, "step2")
    assert state == "ran" and rc != 0, f"step2 of a killed job records exit code {rc}"
    assert prov["outputs"] is None
    out = m10.out("A1")
    assert not (out / "outputs.sha256").exists()
    failed = out / "FAILED"
    assert failed.is_file(), "a killed step writes FAILED"
    text = failed.read_text()
    assert "step2" in text, text
    assert _failed_code(text) == code, text


def test_b1_a_signal_during_step1_is_a_failure_and_leaves_no_marker(m10: M10Sandbox) -> None:
    proc, msg = _signal_job(m10, signal.SIGTERM, "step1")
    assert proc.returncode in {143, -signal.SIGTERM}, f"exit {proc.returncode}\n{msg}"
    prov = m10.provenance("A1")
    assert prov["status"] == "failed" and int(prov["exit_code"]) == 143, prov["status"]
    state, rc = step_state(prov, "step1")
    assert state == "ran" and rc != 0, f"step1 of a killed job records exit code {rc}"
    assert step_state(prov, "step2")[0] == "not_reached"
    out = m10.out("A1")
    assert not (out / "step1.done").exists()
    failed = out / "FAILED"
    assert failed.is_file(), "a killed step writes FAILED"
    assert "step1" in failed.read_text() and _failed_code(failed.read_text()) == 143


# ------------------------------------------------------------ S1: exactly 16 CPUs


def _directives(path: Path) -> dict[str, str]:
    found = re.findall(r"^#SBATCH\s+--([a-z-]+)[= ](\S+)", path.read_text(), re.M)
    return dict(found)


def _mem_mb(value: str) -> float:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)B?", value.strip(), re.I)
    assert m, f"unparsable --mem {value!r}"
    return float(m.group(1)) * {"": 1 / 1024**2, "K": 1 / 1024, "M": 1, "G": 1024,
                                "T": 1024**2}[m.group(2).upper()]  # fmt: skip


def test_s1_the_search_header_asks_for_16_cpus_and_64000m() -> None:
    d = _directives(M10_SCRIPTS["diann_pilot"])
    assert d.get("cpus-per-task") == "16"
    assert d.get("mem") == "64000M", f"--mem={d.get('mem')}: want 64000M (16 x 4000 MB)"


@pytest.mark.parametrize("name", sorted(M10_SCRIPTS))
def test_s1_no_job_asks_for_more_memory_than_its_cpus_carry(name: str) -> None:
    """With MaxMemPerCPU=4000, Slurm raises the CPU count of a job whose --mem exceeds
    cpus x 4000 MB (review S1: 64G gave 17 CPUs per task and 18 allocated), which breaks the
    R14 bookkeeping."""
    d = _directives(M10_SCRIPTS[name])
    cpus = int(d.get("cpus-per-task", "1"))
    mem = _mem_mb(d["mem"])
    assert mem <= cpus * MAX_MEM_PER_CPU_MB, (
        f"{name}: --mem={d['mem']} ({mem:.0f} MB) > {cpus} x {MAX_MEM_PER_CPU_MB} MB"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_s1_both_calls_pass_threads_16(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm))
    for call in m10.searches():
        assert flag_values(diann_args(call), "--threads") == ["16"]
    for argv in commands(m10.provenance(arm)):
        assert flag_values(argv, "--threads") == ["16"]


# ------------------------------------------------------------ S2: the marker names its job


@pytest.mark.parametrize("arm", ARMS)
def test_s2_step1_done_records_the_job_commit_and_input_hashes(m10: M10Sandbox, arm: str) -> None:
    _ok(_run(m10, arm, SLURM_JOB_ID="1111"))
    text = _marker(m10, arm)
    first = text.splitlines()[0].strip()
    m = re.fullmatch(r"([0-9a-f]{64})\s+\*?(\S+)", first)
    assert m and Path(m.group(2)).name == f"{arm}.predicted.speclib", first
    assert m.group(1) == sha256(m10.out(arm) / f"{arm}.predicted.speclib")
    assert "1111" in text, f"step1.done lacks the job ID:\n{text}"
    assert m10.sb.head in text, f"step1.done lacks the commit:\n{text}"
    cfg = m10.repo / "validation" / "pilot" / f"diann_{arm}.cfg"
    assert sha256(cfg) in text, f"step1.done lacks the cfg sha256:\n{text}"
    assert STAND_IN_SHA256 in text, f"step1.done lacks the image sha256:\n{text}"
    gz = m10.repo / "validation" / "pilot" / "db_a2" / "strain_db.faa.gz"
    for f in _step1_fastas(m10):
        ok = sha256(Path(f)) in text or (f.endswith("strain_db.faa") and sha256(gz) in text)
        assert ok, f"step1.done lacks the sha256 of --fasta {f}:\n{text}"


def test_s2_a_skipped_step1_names_the_job_and_commit_that_built_it(m10: M10Sandbox) -> None:
    _ok(_run(m10, SLURM_JOB_ID="1111"))
    built_at = m10.sb.head
    (m10.repo / "UNRELATED.txt").write_text("a commit that changes no input of step 1\n")
    m10.commit("unrelated change")
    assert m10.sb.head != built_at
    m10.reset_calls()
    _ok(_run(m10, SLURM_JOB_ID="2222"))
    assert len(m10.searches()) == 1, "an unrelated commit must not make step 1 rerun"
    prov = m10.provenance("A1")
    assert step_state(prov, "step1")[0] == "skipped"
    step1 = json.dumps(prov["step1"])
    assert "1111" in step1, f"provenance step1 lacks the job that built the library: {step1}"
    assert built_at in step1, f"provenance step1 lacks the commit that built it: {step1}"
    # compare the job id field, not a substring: a commit hash can contain "2222" (CI flake)
    built_by = prov["step1"]["built_by"]
    assert str(built_by["job_id"]) != "2222" and m10.sb.head not in step1, step1
    assert str(prov["job"]["id"]) == "2222"


# ------------------------------------------------------------ S4: stale inputs rerun step 1


def _change_input(m10: M10Sandbox, what: str, fastas: list[str]) -> None:
    if what == "cfg":
        cfg = m10.repo / "validation" / "pilot" / "diann_A1.cfg"
        cfg.write_text(cfg.read_text() + "# a changed cfg\n")
        m10.commit("change the A1 cfg")
    elif what == "fasta":
        f = Path(fastas[-1])
        f.write_text(f.read_text() + ">extra_protein\nMKKAAAAAAR\n")
    else:
        new = b"another stub sif\n"
        m10.stand_in.chmod(0o644)
        m10.stand_in.write_bytes(new)
        m10.stand_in.chmod(0o444)
        m10.pin_file.write_text(f"{sha256(m10.stand_in)}  diann_2.3.1.sif\n")
        m10.commit("re-pin the image")


@pytest.mark.parametrize("what", ["cfg", "fasta", "image"])
def test_s4_a_marker_from_other_inputs_makes_step1_rerun(m10: M10Sandbox, what: str) -> None:
    _ok(_run(m10))
    fastas = _step1_fastas(m10)
    _change_input(m10, what, fastas)
    m10.reset_calls()
    _ok(_run(m10))
    assert len(m10.searches()) == 2, f"a changed {what} must make step 1 rerun"
    prov = m10.provenance("A1")
    assert step_state(prov, "step1") == ("ran", 0)
    text = _marker(m10)
    if what == "cfg":
        assert sha256(m10.repo / "validation" / "pilot" / "diann_A1.cfg") in text
    elif what == "fasta":
        assert sha256(Path(fastas[-1])) in text
    else:
        assert sha256(m10.stand_in) in text


# ------------------------------------------------------------ S4: the unpacked A2 FASTA


def _gz(m10: M10Sandbox) -> Path:
    return m10.repo / "validation" / "pilot" / "db_a2" / "strain_db.faa.gz"


def _unpacked(m10: M10Sandbox) -> Path:
    return m10.out("A2") / "fasta" / "strain_db.faa"


def _sidecar(m10: M10Sandbox) -> Path:
    f = _unpacked(m10)
    return f.with_name(f.name + ".src.sha256")


def _strain_fastas(m10: M10Sandbox) -> list[str]:
    return [f for c in m10.searches() for f in flag_values(diann_args(c), "--fasta")
            if f.endswith("strain_db.faa")]  # fmt: skip


def test_s4_a_truncated_unpacked_fasta_without_a_sidecar_is_rebuilt(m10: M10Sandbox) -> None:
    f = _unpacked(m10)
    f.parent.mkdir(parents=True)
    f.write_text(">truncated_by_a_timeout\nMK\n")
    _ok(_run(m10, "A2"))
    assert _strain_fastas(m10) == [str(f), str(f)], _strain_fastas(m10)
    assert f.read_bytes() == gzip.decompress(_gz(m10).read_bytes()), "the truncated file is reused"
    assert _sidecar(m10).is_file(), f"{_sidecar(m10)} missing"
    assert sha256(_gz(m10)) in _sidecar(m10).read_text()


def test_s4_a_changed_gz_is_unpacked_again_and_step1_reruns(m10: M10Sandbox) -> None:
    _ok(_run(m10, "A2"))
    new = b">sandbox_protein\nMAAAKAAAAAAAR\n"
    _gz(m10).write_bytes(gzip.compress(new))
    m10.commit("change the committed strain database")
    m10.reset_calls()
    _ok(_run(m10, "A2"))
    assert _unpacked(m10).read_bytes() == new, "the old unpacked FASTA is reused"
    assert sha256(_gz(m10)) in _sidecar(m10).read_text()
    assert len(m10.searches()) == 2, "a changed FASTA must make step 1 rerun"


def test_s4_an_interrupted_unpack_leaves_no_fasta_in_place(m10: M10Sandbox) -> None:
    m10.install_gzip_wrappers()
    res = _run(m10, "A2", FAKE_GZIP_FAIL="1")
    assert (m10.sb.state / "gzip.log").is_file(), "the A2 FASTA is unpacked with gzip/gunzip"
    assert res.returncode != 0, f"a failed unpack must stop the job\n{res.stdout}\n{res.stderr}"
    assert not _unpacked(m10).exists(), "a partial unpack was left at the final path"
    assert not m10.searches()
    _ok(_run(m10, "A2"))
    assert _unpacked(m10).read_bytes() == gzip.decompress(_gz(m10).read_bytes())


# ------------------------------------------------------------ S3: the banner before the install


def _bad_banner() -> dict[str, str]:
    return {"FAKE_DIANN_BANNER": "DIA-NN 1.9.2 (Data-Independent Acquisition)"}


def test_s3_a_failed_banner_installs_nothing(m10i: M10Sandbox) -> None:
    res = m10i.run("diann_image", **_bad_banner())
    _rc(res, 2)
    assert not m10i.copy.exists(), "a copy whose banner failed was installed"
    images = m10i.copy.parent
    left = sorted(p.name for p in images.iterdir() if ".tmp" in p.name) if images.is_dir() else []
    assert not left, f"temporary copy left behind: {left}"


def test_s3_a_rerun_after_a_failed_banner_checks_the_banner_again(m10i: M10Sandbox) -> None:
    _rc(m10i.run("diann_image", **_bad_banner()), 2)
    res = m10i.run("diann_image", **_bad_banner())
    msg = _rc(res, 2)
    assert "already installed" not in msg, msg
    _ok(m10i.run("diann_image"))
    assert m10i.copy.is_file() and m10i.copy.stat().st_mode & 0o777 == 0o444
    assert sha256(m10i.copy) == STAND_IN_SHA256


# ------------------------------------------------------------ N1, N2, N5


@pytest.mark.parametrize(
    "env, step",
    [({"FAKE_DIANN_NO_LIB": "1"}, "step1"), ({"FAKE_DIANN_NO_REPORT": "1"}, "step2")],
    ids=["no-speclib", "no-report"],
)
def test_n1_failed_never_records_exit_code_0(
    m10: M10Sandbox, env: dict[str, str], step: str
) -> None:
    _rc(_run(m10, **env), 1)
    text = (m10.out("A1") / "FAILED").read_text()
    assert step in text, text
    assert _failed_code(text) >= 1, f"FAILED says exit code 0 for a failure: {text!r}"
    assert m10.provenance("A1")["status"] == "failed"


@pytest.mark.parametrize("kind", ["missing", "not-executable", "broken"])
def test_n2_a_broken_gippa_python_is_named(m10: M10Sandbox, kind: str) -> None:
    py = m10.sb.tmp / "badpy" / "python"
    py.parent.mkdir()
    if kind == "not-executable":
        py.write_text("#!/bin/sh\nexit 0\n")
        py.chmod(0o644)
    elif kind == "broken":
        py.write_text("#!/bin/sh\necho 'python: error while loading shared libraries' >&2\n"
                      "exit 127\n")  # fmt: skip
        py.chmod(0o755)
    msg = _rc(_run(m10, GIPPA_PYTHON=str(py)), 2)
    assert "GIPPA_PYTHON" in msg, f"the message must name GIPPA_PYTHON:\n{msg}"
    assert not m10.searches()


def test_n5_a_step1_rerun_removes_the_old_log(m10: M10Sandbox) -> None:
    out = m10.out("A1")
    out.mkdir(parents=True)
    log = out / "A1.log.txt"
    log.write_text("the log of an earlier step 1\n")
    _ok(_run(m10))
    step1 = m10.searches()[0]
    assert step1["preexisting"].get(str(log)) is False, "the old A1.log.txt survives step 1"


# ------------------------------------------------------------ test gaps of the review


def _tree(d: Path) -> list[str]:
    return sorted(str(p.relative_to(d)) for p in d.rglob("*"))


def test_gap_a1_and_a2_have_their_own_temp_and_raw_dir_is_untouched(m10: M10Sandbox) -> None:
    before = _tree(m10.raw_dir)
    temps: dict[str, str] = {}
    for arm in ARMS:
        m10.reset_calls()
        _ok(_run(m10, arm))
        step2 = m10.searches()[-1]
        assert step2["temp"], f"{arm}: step 2 has no --temp"
        temps[arm] = step2["temp"]
        assert Path(temps[arm]).is_relative_to(m10.out(arm)), (arm, temps[arm])
        assert _tree(m10.raw_dir) == before, f"{arm} wrote under RAW_DIR"
    assert temps["A1"] != temps["A2"], temps
    assert list(Path(temps["A1"]).glob("*.quant")), "the fake wrote no .quant under --temp"


def test_gap_a_stale_quant_dir_is_emptied_before_step2(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    stale = m10.out("A1") / "quant" / "old_run.d.quant"
    stale.write_text("from an earlier attempt\n")
    m10.reset_calls()
    _ok(_run(m10))
    (step2,) = m10.searches()
    assert "old_run.d.quant" not in (step2["temp_listing"] or []), step2["temp_listing"]
    assert not stale.exists()


def test_gap_outputs_sha256_leaves_out_the_jobs_own_records(m10: M10Sandbox) -> None:
    _ok(_run(m10))
    out = m10.out("A1")
    names = [re.split(r"\s+\*?", ln.strip(), maxsplit=1)[1]
             for ln in (out / "outputs.sha256").read_text().splitlines()]  # fmt: skip
    base = {Path(n).name for n in names}
    for own in ("command_step1.txt", "command_step2.txt", "step1.done", "diann_version.txt",
                "provenance.json", "FAILED", "outputs.sha256"):  # fmt: skip
        assert own not in base, f"outputs.sha256 lists the job's own {own}"
    assert any("quant/" in n for n in names), "the .quant files are listed"
    assert (out / "command_step1.txt").is_file()
