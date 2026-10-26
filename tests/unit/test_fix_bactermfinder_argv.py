"""BACTERMFINDER calls genome_scan.py with its four arguments from a staged working directory
(docs/specs/fix-bactermfinder-argv.md).

At the pinned commit (ADR-0008) ``genome_scan.py <fasta> <step> <prefix> <batch>`` reads its models
and ``iLearnPlus`` relative to the working directory and writes its intermediates there; the result
table is ``<prefix><fasta>_mean.csv``. The behavioural tests render the ``script:`` block and run it
in bash against a fake, read-only ``BACTERMFINDER_HOME`` whose ``genome_scan.py`` mimics that I/O.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MAIN_NF = ROOT / "workflow" / "modules" / "local" / "bactermfinder" / "main.nf"
MODELS = [
    f"deep-bactermfinder_3cnn_1d_1cat_reduced_10x_{e}.csv_saved_model.h5"
    for e in ("ENAC", "NCP", "PS2", "binary")
]
ARGS = ["genome.fa", "3", "out", "10000"]  # README example; GIAP_SPEC.md §7.1.1 "per its README"
RESULT = "outgenome.fa_mean.csv"
HEADER = (
    "SampleName,probability_binary,probability_ENAC,probability_PS2,probability_NCP,"
    "probability_mean"
)


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _section(name: str) -> str:
    """The triple-quoted body of the process ``script:`` or ``stub:`` section."""
    text = _strip_comments(MAIN_NF.read_text())
    m = re.search(rf"\b{name}\s*:\s*\"\"\"(.*?)\"\"\"", text, flags=re.S)
    assert m, f'{MAIN_NF}: no {name}: """ ... """ block'
    return m.group(1)


def _script() -> str:
    return _section("script")


def _bash(body: str) -> str:
    """Render a Nextflow script body as bash (``\\$`` -> ``$``, ``${task.process}`` filled)."""
    body = body.replace("${task.process}", "BACTERMFINDER").replace("${meta_id}", "S1")
    assert "${" not in body.replace("\\${", ""), "unexpected Groovy interpolation in script block"
    return textwrap.dedent(body.replace("\\$", "$"))  # Nextflow strips the common indent


def _genome_scan_calls(script: str) -> list[list[str]]:
    calls = []
    for line in script.splitlines():
        if "genome_scan.py" not in line or not re.search(r"\bpython3?\b", line):
            continue
        toks = shlex.split(line.split(">")[0].replace("\\$", "$"))
        i = next(i for i, t in enumerate(toks) if t.endswith("genome_scan.py"))
        calls.append(toks[i + 1 :])
    return calls


# ---------------------------------------------------------------- static checks


def test_genome_scan_gets_four_positional_args() -> None:
    calls = _genome_scan_calls(_script())
    assert len(calls) == 1, f"expected one genome_scan.py call, got {calls}"
    assert calls[0] == ARGS, (
        f"genome_scan.py args {calls[0]} != {ARGS} (<fasta> <step> <prefix> <batch>)"
    )


def test_prefix_has_no_dot_and_genome_is_a_bare_name() -> None:
    fasta, _, prefix, _ = _genome_scan_calls(_script())[0]
    assert "." not in prefix  # the script uses prefix.split(".")[0] for one output name
    assert "/" not in fasta  # embedded in the output names


def test_script_stages_bactermfinder_home_and_never_cds_into_it() -> None:
    s = _script()
    assert re.search(r"\bcd\s+\S*btf_work\b", s), "script must cd into the staged btf_work/ dir"
    assert not re.search(r"\bcd\s+[\"']?\\?\$\{?BACTERMFINDER_HOME", s), (
        "must not cd into the shared (read-only in images) BACTERMFINDER_HOME"
    )
    assert re.search(r"\bln\s+-s[f]*\b[^\n]*BACTERMFINDER_HOME", s), (
        "genome_scan.py, iLearnPlus and the models are symlinked from BACTERMFINDER_HOME"
    )
    for item in ("genome_scan.py", "iLearnPlus"):
        assert re.search(rf"BACTERMFINDER_HOME\S*/{re.escape(item)}", s), f"{item} not staged"
    assert re.search(r"BACTERMFINDER_HOME\S*/\S*\*\S*\.h5|saved_model\.h5", s), "models not staged"


def test_bactermfinder_home_default_and_version_line_kept() -> None:
    s = _script()
    assert "BACTERMFINDER_HOME:-\\$CONDA_PREFIX/share/bactermfinder" in s
    assert "GIAP_COMMIT" in s


def test_result_table_copied_into_output_dir() -> None:
    s = _script()
    assert re.search(rf"\b(cp|mv|ln)\b[^\n]*{re.escape(RESULT)}[^\n]*\bbactermfinder/?", s), (
        f"{RESULT} must land in bactermfinder/"
    )
    assert not re.search(r"\b(cp|mv)\b[^\n]*\*[^\n]*\.csv[^\n]*bactermfinder", s), (
        "no CSV glob into bactermfinder/: the step reader rejects several CSVs"
    )


def test_output_globs_unchanged() -> None:
    text = _strip_comments(MAIN_NF.read_text())
    out = re.search(r"\boutput\s*:(.*?)\bscript\s*:", text, flags=re.S)
    assert out
    assert re.search(r"path\s+'bactermfinder'\s*,\s*emit:\s*dir", out.group(1))
    assert re.search(r"path\s+'versions\.yml'\s*,\s*emit:\s*versions", out.group(1))
    assert "mkdir bactermfinder" in _section("stub")


# ---------------------------------------------------------------- behaviour (fake tool)

FAKE_GENOME_SCAN = r"""
import json, os, sys
from pathlib import Path
here = Path.cwd()
assert Path("iLearnPlus/util/FileProcessing.py").is_file(), "iLearnPlus not in cwd"
for m in MODELS:
    assert Path(m).is_file(), f"model {m} not in cwd"
genome, step, prefix, batch = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
assert len(sys.argv) == 5
Path(os.environ["BTF_RECORD"]).write_text(json.dumps({"argv": sys.argv[1:], "cwd": str(here)}))
assert Path(genome).is_file()
Path("Sample.csv").write_text("x\n")
Path("df_sample.fasta").write_text(">w\nACGT\n")
Path("output_sample").mkdir(exist_ok=True)
Path("output_sample/ENAC-0.csv").write_text("x\n")
Path(prefix.split(".")[0] + f"_{genome}_sliding_windows.csv").write_text("x\n")
for e in ("binary", "ENAC", "PS2", "NCP"):
    Path(f"{prefix}_{e}{genome}.csv").write_text("x\n")
Path(f"{prefix}{genome}_mean.csv").write_text(HEADER + "\n1_chr_0_101_+,0.1,0.1,0.1,0.1,0.1\n")
print("Done")
"""


@pytest.fixture
def ran(tmp_path: Path) -> dict[str, object]:
    home = tmp_path / "share" / "bactermfinder"
    (home / "iLearnPlus" / "util").mkdir(parents=True)
    (home / "iLearnPlus" / "util" / "FileProcessing.py").write_text("")
    for m in MODELS:
        (home / m).write_text("h5")
    (home / "GIAP_COMMIT").write_text("024870447f6288ebc9d29b314e957a6dfdfeb5f6\n")
    (home / "genome_scan.py").write_text(
        f"MODELS = {MODELS!r}\nHEADER = {HEADER!r}\n" + FAKE_GENOME_SCAN
    )
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "python").symlink_to(sys.executable)
    for p in [*sorted(home.rglob("*"), reverse=True), home]:  # read-only, as in an image
        p.chmod(p.stat().st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    task = tmp_path / "task"
    task.mkdir()
    (task / "genome.fa").write_text(">chr\n" + "ACGT" * 60 + "\n")
    record = tmp_path / "record.json"
    env = {
        **os.environ,
        "BACTERMFINDER_HOME": str(home),
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "BTF_RECORD": str(record),
    }
    bash = shutil.which("bash")
    assert bash
    proc = subprocess.run(
        [bash, "-euo", "pipefail", "-c", _bash(_script())],
        cwd=task,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    for p in [home, *home.rglob("*")]:
        p.chmod(p.stat().st_mode | stat.S_IWUSR)
    return {"proc": proc, "task": task, "home": home, "record": record}


def _ok(ran: dict[str, object]) -> subprocess.CompletedProcess[str]:
    proc = ran["proc"]
    assert isinstance(proc, subprocess.CompletedProcess)
    assert proc.returncode == 0, f"script failed:\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    return proc


def test_runs_against_read_only_home_with_readme_args(ran: dict[str, object]) -> None:
    import json

    _ok(ran)
    rec = json.loads(Path(str(ran["record"])).read_text())
    assert rec["argv"] == ARGS
    task = Path(str(ran["task"]))
    assert Path(rec["cwd"]).resolve() == (task / "btf_work").resolve()


def test_output_dir_holds_exactly_the_mean_table(ran: dict[str, object]) -> None:
    _ok(ran)
    out = Path(str(ran["task"])) / "bactermfinder"
    assert out.is_dir()
    csvs = sorted(f.name for f in out.iterdir() if f.is_file() and f.suffix == ".csv")
    assert csvs == [RESULT]
    assert (out / RESULT).read_text().splitlines()[0] == HEADER


def test_step_reader_accepts_the_output_dir(
    ran: dict[str, object], monkeypatch: pytest.MonkeyPatch
) -> None:
    """``_bactermfinder`` raises on several CSVs; with one it hands that file to the parser."""
    from giap.steps import terminators as st

    _ok(ran)
    seen: list[Path] = []
    monkeypatch.setattr(st, "parse_bactermfinder", lambda p, reps: seen.append(p) or [])
    st._bactermfinder(Path(str(ran["task"])) / "bactermfinder", {})
    assert [p.name for p in seen] == [RESULT]


def test_home_untouched_and_versions_written(ran: dict[str, object]) -> None:
    _ok(ran)
    home = Path(str(ran["home"]))
    for junk in ("Sample.csv", "df_sample.fasta", "output_sample", RESULT):
        assert not (home / junk).exists()
    versions = (Path(str(ran["task"])) / "versions.yml").read_text()
    assert "bactermfinder: 024870447f6288ebc9d29b314e957a6dfdfeb5f6" in versions
