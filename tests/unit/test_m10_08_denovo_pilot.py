"""``denovo_pilot.sbatch``, the Cascadia branch, in a sandbox (M10-08b, deliverable 8).

Spec: docs/specs/M10-08-a3-feasibility.md, deliverable 8, "Commands (B2)", the smoke-run table
(step 4), edge cases "Weights sha256 differs", "CUDA not available" and acceptance criterion 9.
The M3-24 assertions on this script (``tests/unit/test_m3_24_pilot_kit.py``) must keep passing
unchanged; they are not repeated here.

The script runs with ``bash`` in ``helpers.a3_kit.a3_sandbox`` (fake Slurm tools, fake
``micromamba`` that logs every call with its cwd, an unreachable proxy). The ``giap-cascadia`` env
is a real venv holding the 08a stub ``cascadia``/``pyteomics``/``psims`` packages plus, added
here, a stub ``torch`` (``$STUB_CUDA`` decides ``torch.cuda.is_available()``) and a stub
``cascadia`` executable for ``--help``. A fake ``nvidia-smi`` logs its calls. The real 08a
launcher ``scripts/pilot/cascadia_offline.py`` runs against the stubs. No GPU, no Cascadia, no
network, no real ``$GIAP_ROOT``. No statistical claim, so no calibration test.

Conventions settled by the test author: the run comes from ``$A3_RUN`` or the second argument
(``denovo_pilot.sbatch cascadia WT_2``); the expected sha256 values come from
``$GIAP_REPO/validation/pilot/a3_assets.tsv`` (as in ``fetch_a3_assets.sbatch``); a missing mzML
for a requested run is a precondition (exit 2).
"""

from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_kit import (
    A3_ASSETS,
    STUB_UNIMOD_XML,
    A3Sandbox,
    a3_sandbox,
    example11_spectra,
    require,
    sha256_hex,
    weights_payload,
    write_mzml,
)
from helpers.slurm_sandbox import write_exe

RUNS = ("WT_1", "WT_2", "WT_3")
WEIGHTS = weights_payload("sandbox")

_TORCH = """import os
__version__ = "2.0.1"


class cuda:
    @staticmethod
    def is_available():
        return os.environ.get("STUB_CUDA", "1") == "1"
"""

_CASCADIA_BIN = """#!{python}
import json, os, sys
with open(os.path.join(os.environ["FAKE_SLURM_DIR"], "cascadia_bin.jsonl"), "a") as fh:
    fh.write(json.dumps({"argv": sys.argv[1:]}) + "\\n")
print("usage: cascadia sequence [-h] [-o OUT] [-t SCORE_THRESHOLD] [-b BATCH_SIZE] "
      "spectrum_file model")
"""

# Appended to the 08a stub cascadia/cascadia.py: like Cascadia 0.0.7, a finished run writes
# ``<-o>.ssl`` (header only here). ``$STUB_ECHO`` is printed first (a log line, e.g. a Unimod
# URL); ``$STUB_NO_SSL=1`` finishes with exit 0 but no ``.ssl``. ``$STUB_RC`` still exits first.
_SSL_WRITER = """

_stub_main = main


def main():
    import os, sys
    if os.environ.get("STUB_ECHO"):
        print(os.environ["STUB_ECHO"], flush=True)
    _stub_main()
    if os.environ.get("STUB_NO_SSL") == "1":
        return
    a = sys.argv[1:]
    with open(a[a.index("-o") + 1] + ".ssl", "w") as fh:
        fh.write("file\\tscan\\tcharge\\tsequence\\tscore-type\\tscore\\tretention-time\\n")
"""

_NVIDIA_SMI = """#!/bin/bash
echo "$*" >> "$FAKE_SLURM_DIR/nvidia_smi.log"
echo "NVIDIA L40S, 46068 MiB"
"""


@pytest.fixture(scope="module")
def template(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return write_mzml(tmp_path_factory.mktemp("tpl") / "t.mzML", example11_spectra())


@pytest.fixture
def box(tmp_path: Path, template: Path) -> A3Sandbox:
    return a3_sandbox(tmp_path, template)


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return res.stdout + res.stderr


def _msg(box: A3Sandbox, res: subprocess.CompletedProcess[str]) -> str:
    """The output without sandbox paths (the tmp dir carries the test's name)."""
    return _out(res).replace(str(box.sb.tmp), "<tmp>").lower()


def _ready(
    box: A3Sandbox,
    runs: tuple[str, ...] = RUNS,
    *,
    weights_sha: str | None = None,
    unimod_sha: str | None = None,
) -> None:
    """Env, weights, Unimod table, per-run mzML and a3_assets.tsv rows pinning their sha256."""
    env = box.stub_cascadia_env()
    purelib = subprocess.run(
        [str(env / "bin" / "python"), "-I", "-c",
         "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()  # fmt: skip
    (Path(purelib) / "torch").mkdir(exist_ok=True)
    (Path(purelib) / "torch" / "__init__.py").write_text(_TORCH)
    with (Path(purelib) / "cascadia" / "cascadia.py").open("a") as fh:
        fh.write(_SSL_WRITER)
    import sys

    write_exe(env / "bin" / "cascadia", _CASCADIA_BIN.replace("{python}", sys.executable))
    write_exe(box.sb.bin / "nvidia-smi", _NVIDIA_SMI)
    d = box.db / "cascadia"
    d.mkdir(parents=True, exist_ok=True)
    (d / "cascadia_0.0.7.ckpt").write_bytes(WEIGHTS)
    (d / "unimod_tables.xml").write_bytes(STUB_UNIMOD_XML)
    for run in runs:
        write_mzml(box.data / "mzml" / f"{run}.mzML", example11_spectra())
    with require(A3_ASSETS, "deliverable 5").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        if r["asset"] == "cascadia_weights":
            r["sha256"] = sha256_hex(WEIGHTS) if weights_sha is None else weights_sha
            r["bytes"] = str(len(WEIGHTS))
        elif r["asset"] == "unimod_table":
            r["sha256"] = sha256_hex(STUB_UNIMOD_XML) if unimod_sha is None else unimod_sha
            r["bytes"] = str(len(STUB_UNIMOD_XML))
    box.write_assets(rows)


def _run(box: A3Sandbox, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
    env.setdefault("STUB_LOG", str(box.sb.tmp / "stub_cascadia.json"))
    return box.run("denovo_pilot", "cascadia", *args, **env)  # type: ignore[arg-type]


def _sequence_calls(box: A3Sandbox) -> list[dict[str, Any]]:
    return [c for c in box.calls("micromamba.jsonl") if "sequence" in c["argv"]]


def _same(a: str, b: Path) -> bool:
    return Path(a).resolve() == b.resolve()


# ------------------------------------------------------------------ the call itself


def test_one_run_calls_cascadia_through_the_offline_launcher(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 0, _out(res)
    calls = _sequence_calls(box)
    assert len(calls) == 1, calls
    argv = calls[0]["argv"]
    i = argv.index("sequence")
    assert argv[:3] == ["run", "-n", "giap-cascadia"] or (
        "-n" in argv and argv[argv.index("-n") + 1] == "giap-cascadia"
    ), argv
    pre = argv[:i]
    assert "python" in [Path(a).name for a in pre] and "-I" in pre, argv
    assert pre[-2].endswith("scripts/pilot/cascadia_offline.py"), argv
    assert _same(pre[-1], box.db / "cascadia" / "unimod_tables.xml"), argv
    post = argv[i + 1 :]
    assert _same(post[0], box.data / "mzml" / "WT_1.mzML"), "first positional: the mzML"
    assert _same(post[1], box.db / "cascadia" / "cascadia_0.0.7.ckpt"), "second: the checkpoint"
    assert "-o" in post and _same(post[post.index("-o") + 1],
                                  box.root / "pilot" / "A3_cascadia" / "WT_1" / "WT_1")  # fmt: skip
    assert "-b" in post and post[post.index("-b") + 1] == "32"
    for bad in ("--model", "--num_workers", "--output", "-m", "-t", "--score_threshold"):
        assert bad not in argv, f"{bad} in the Cascadia call: {argv}"
    assert (box.root / "pilot" / "A3_cascadia" / "WT_1").is_dir()


def test_launcher_hands_cascadia_the_sequence_argv(box: A3Sandbox) -> None:
    _ready(box)
    log = box.sb.tmp / "stub_cascadia.json"
    res = _run(box, A3_RUN="WT_1", STUB_LOG=str(log))
    assert res.returncode == 0, _out(res)
    import json

    rec = json.loads(log.read_text())
    assert rec["argv"][0] == "sequence" and len(rec["argv"]) == 7, rec["argv"]
    assert _same(rec["url_in_main"], box.db / "cascadia" / "unimod_tables.xml")


def test_an_inherited_run_variable_is_ignored(box: A3Sandbox) -> None:
    """A login shell may export ``RUN`` for something else (FarmShare's BIOS 270 ~/.bashrc sets it
    to a singularity command); sbatch exports it into the job, so the script reads ``A3_RUN``
    only (lead fix, 2026-10-24)."""
    _ready(box)
    res = _run(box, RUN="singularity run -B /farmshare/x /farmshare/y.sif")
    assert res.returncode == 0, _out(res)
    assert len(_sequence_calls(box)) == 3


def test_working_directory_is_the_runs_tmpdir(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 0, _out(res)
    (call,) = _sequence_calls(box)
    tmpdir = Path(box.sb.env["TMPDIR"])
    cwd = Path(call["cwd"]).resolve()
    assert cwd.name == "WT_1" and tmpdir.resolve() in cwd.parents, call["cwd"]  # see S2 below


def test_default_runs_all_three_in_order(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box)
    assert res.returncode == 0, _out(res)
    calls = _sequence_calls(box)
    assert len(calls) == 3, "one Cascadia call per run"
    for call, run in zip(calls, RUNS, strict=True):
        post = call["argv"][call["argv"].index("sequence") + 1 :]
        assert _same(post[0], box.data / "mzml" / f"{run}.mzML")
        assert _same(post[post.index("-o") + 1], box.root / "pilot" / "A3_cascadia" / run / run)
        cwd = Path(call["cwd"]).resolve()
        assert cwd.name == run and Path(box.sb.env["TMPDIR"]).resolve() in cwd.parents


def test_run_as_the_second_argument(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, "WT_2")
    assert res.returncode == 0, _out(res)
    (call,) = _sequence_calls(box)
    assert _same(call["argv"][call["argv"].index("sequence") + 1], box.data / "mzml" / "WT_2.mzML")


def test_help_gpu_and_threshold_are_logged_first(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 0, _out(res)
    calls = box.calls("micromamba.jsonl")
    first_seq = next(i for i, c in enumerate(calls) if "sequence" in c["argv"])
    before = calls[:first_seq]
    assert any("--help" in c["argv"] for c in before), "cascadia --help is logged first"
    assert any("torch" in " ".join(c["argv"]) and "cuda" in " ".join(c["argv"])
               for c in before), "torch.cuda.is_available() is checked before Cascadia"  # fmt: skip
    assert (box.sb.state / "nvidia_smi.log").is_file(), "nvidia-smi is logged"
    assert "0.8" in _out(res), "the default score threshold 0.8 is recorded in the log"


# ------------------------------------------------------------------ refusals


@pytest.mark.parametrize("which", ["weights", "unimod"])
def test_sha256_mismatch_exits_2_before_cascadia(box: A3Sandbox, which: str) -> None:
    wrong = "0" * 64
    _ready(box, weights_sha=wrong if which == "weights" else None,
           unimod_sha=wrong if which == "unimod" else None)  # fmt: skip
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 2, _out(res)
    assert "sha256" in _msg(box, res), "exit 2 for the sha256 check, not another precondition"
    assert not _sequence_calls(box)


def test_unpinned_weights_exit_2(box: A3Sandbox) -> None:
    # The committed a3_assets.tsv has no sha256 until the fetch is resolved: never run unverified.
    _ready(box, weights_sha="")
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 2, _out(res)
    assert "sha256" in _msg(box, res), "exit 2 for the missing pin, not another precondition"
    assert not _sequence_calls(box)


def test_no_cuda_exits_1_before_cascadia(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_CUDA="0")
    assert res.returncode == 1, _out(res)
    assert not _sequence_calls(box)


def test_missing_mzml_is_a_precondition(box: A3Sandbox) -> None:
    _ready(box, runs=("WT_1",))
    res = _run(box, A3_RUN="WT_3")
    assert res.returncode == 2, _out(res)
    assert "wt_3.mzml" in _msg(box, res), "name the missing input"
    assert not _sequence_calls(box)


def test_cascadia_failure_propagates(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_RC="7")
    assert len(_sequence_calls(box)) == 1, _out(res)
    assert res.returncode != 0, _out(res)


# ------------------------------------------------------------------ review round 1 (M10-08b)
# Findings of the 08b code review, lead decisions of 2026-10-24: S1, S2, S4b, N4, N5.

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "slurm" / "denovo_pilot.sbatch"
MAX_MEM_PER_CPU_MB = 4000  # `scontrol show partition gpu`: DefMemPerCPU = MaxMemPerCPU = 4000
UNIMOD_URL = "http://www.unimod.org/obo/unimod.obo"
_SCRATCH = re.compile(r"giap-a3\.[A-Za-z0-9]{6,}")


def _directive(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}[= ](\S+)", path.read_text(), re.M)
    return m.group(1) if m else None


def _mem_mb(value: str) -> float:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)B?", value.strip(), re.I)
    assert m, f"unparsable --mem {value!r}"
    return float(m.group(1)) * {"": 1, "K": 1 / 1024, "M": 1, "G": 1024,
                                "T": 1024**2}[m.group(2).upper()]  # fmt: skip


def test_s1_header_asks_for_16_cpus_and_64000m() -> None:
    # Slurm raises --cpus-per-task to ceil(mem / MaxMemPerCPU): 64G on 8 CPUs became 17.
    assert _directive(SCRIPT, "cpus-per-task") == "16"
    assert _directive(SCRIPT, "mem") == "64000M"


def _gpu_scripts() -> list[Path]:
    d = SCRIPT.parent
    return sorted(p for p in d.glob("*.sbatch") if _directive(p, "partition") == "gpu")


@pytest.mark.parametrize("path", _gpu_scripts(), ids=lambda p: p.name)
def test_s1_house_rule_mem_fits_max_mem_per_cpu(path: Path) -> None:
    """R14 accounting: the CPUs asked for are the CPUs held (mem <= 4000M x CPUs on `gpu`)."""
    cpus, mem = _directive(path, "cpus-per-task"), _directive(path, "mem")
    assert cpus and mem, f"{path.name}: --cpus-per-task {cpus}, --mem {mem}"
    assert _mem_mb(mem) <= MAX_MEM_PER_CPU_MB * int(cpus), (
        f"{path.name}: --mem={mem} > {MAX_MEM_PER_CPU_MB}M x {cpus} CPUs: Slurm silently "
        "raises the CPU count"
    )


def test_s2_scratch_is_a_fresh_per_job_mktemp_dir(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box)
    assert res.returncode == 0, _out(res)
    tmpdir = Path(box.sb.env["TMPDIR"]).resolve()
    cwds = [Path(c["cwd"]).resolve() for c in _sequence_calls(box)]
    assert [c.name for c in cwds] == list(RUNS)
    roots = {c.parent for c in cwds}
    assert len(roots) == 1, f"one scratch root per job: {roots}"
    (root,) = roots
    assert root.parent == tmpdir and _SCRATCH.fullmatch(root.name), (
        f"the scratch must be mktemp -d $TMPDIR/giap-a3.XXXXXX, not {root}"
    )
    assert not root.exists(), "the scratch root is removed on exit"
    assert not list(tmpdir.glob("giap-a3.*")), "nothing left under $TMPDIR"


def test_s2_a_foreign_tmpdir_run_directory_is_untouched(box: A3Sandbox) -> None:
    _ready(box)
    tmpdir = Path(box.sb.env["TMPDIR"])
    theirs = tmpdir / "WT_1"  # another user's (or another job's) directory
    theirs.mkdir()
    (theirs / "their_scratch.bin").write_bytes(b"not ours\n")
    theirs.chmod(0o555)
    try:
        res = _run(box, A3_RUN="WT_1")
        survived = (theirs / "their_scratch.bin").read_bytes() == b"not ours\n"
        mode = theirs.stat().st_mode & 0o777
    finally:
        theirs.chmod(0o755)
    assert res.returncode == 0, _out(res)
    assert survived, "a pre-existing $TMPDIR/WT_1 was modified"
    assert mode == 0o555
    (call,) = _sequence_calls(box)
    assert Path(call["cwd"]).resolve() != theirs.resolve(), "ran in the foreign directory"


def test_s2_scratch_is_removed_when_cascadia_fails(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_RC="7")
    assert res.returncode == 1, _out(res)
    (call,) = _sequence_calls(box)
    cwd = Path(call["cwd"]).resolve()
    assert _SCRATCH.fullmatch(cwd.parent.name), cwd
    assert not cwd.parent.exists(), "the trap removes the scratch root"


def test_s2_unset_tmpdir_falls_back_to_tmp(box: A3Sandbox) -> None:
    _ready(box)
    res = box.run("denovo_pilot", "cascadia", unset=("TMPDIR",), A3_RUN="WT_1",
                  STUB_LOG=str(box.sb.tmp / "stub_cascadia.json"))  # fmt: skip
    assert res.returncode == 0, _out(res)
    (call,) = _sequence_calls(box)
    cwd = Path(call["cwd"]).resolve()
    assert cwd.parent.parent == Path("/tmp").resolve() and _SCRATCH.fullmatch(cwd.parent.name)
    assert not cwd.parent.exists()


def test_s4b_unimod_url_in_the_log_exits_1(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_ECHO=f"downloading {UNIMOD_URL}")
    assert res.returncode == 1, _out(res)
    assert "unimod" in res.stderr.lower(), "say why on stderr: " + res.stderr
    log = box.root / "pilot" / "A3_cascadia" / "WT_1" / "cascadia.log"
    assert log.is_file() and UNIMOD_URL in log.read_text()


def test_s4b_unimod_url_stops_before_the_next_run(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, STUB_ECHO=UNIMOD_URL)
    assert res.returncode == 1, _out(res)
    assert len(_sequence_calls(box)) == 1, "no further run after the launcher failed"


def test_n4_exit_0_without_an_ssl_is_exit_1(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_NO_SSL="1")
    assert len(_sequence_calls(box)) == 1, _out(res)
    assert res.returncode == 1, _out(res)
    assert "wt_1.ssl" in _msg(box, res), "name the missing .ssl"


def test_n4_a_header_only_ssl_is_still_a_success(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1")
    assert res.returncode == 0, _out(res)
    ssl = box.root / "pilot" / "A3_cascadia" / "WT_1" / "WT_1.ssl"
    assert ssl.is_file() and len(ssl.read_text().splitlines()) == 1


def test_n5_a_crash_with_a_unimod_url_gives_the_unimod_message(box: A3Sandbox) -> None:
    _ready(box)
    res = _run(box, A3_RUN="WT_1", STUB_RC="7", STUB_ECHO=f"lxml: failed to load {UNIMOD_URL}")
    assert res.returncode == 1, _out(res)
    assert "unimod" in res.stderr.lower(), (
        "a Cascadia crash on a Unimod fetch must say the launcher did not take effect: "
        + res.stderr
    )
