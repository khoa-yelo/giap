"""``denovo_pilot.sbatch cascadia`` with ``A3_OUT_DIR`` (M10-10a criterion 4, deliverable 1).

Spec: docs/specs/M10-10-a3-run.md, deliverable 1 and "Submission": ``A3_RUN=WT_k
A3_OUT_DIR=$GIAP_ROOT/pilot/A3 denovo_pilot.sbatch cascadia`` writes under
``$A3_OUT_DIR/<RUN>/<RUN>``; unset, M10-08's ``$GIAP_ROOT/pilot/A3_cascadia/<RUN>/<RUN>`` holds,
and M10-08's WT_1 smoke output is never touched by an M10-10 run.

The Cascadia branch itself (``A3_RUN``, the asset sha256 checks, ``cascadia_offline.py … -o …``) is
M10-08b's deliverable 8; these tests drive it in the M10-08 ``a3_sandbox`` with ``micromamba``,
``nvidia-smi`` and a stub ``giap-cascadia`` env (stub ``cascadia``, ``torch.cuda.is_available()``
true) and read the ``-o`` argument from the fake ``micromamba``'s call log. No GPU, no Cascadia.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest
from helpers.a3_kit import A3Sandbox, a3_sandbox, example11_spectra, write_mzml
from helpers.slurm_sandbox import write_exe

_NVIDIA_SMI = """#!/bin/bash
echo "NVIDIA L40S, 46068 MiB"
"""

_TORCH = "class cuda:\n    @staticmethod\n    def is_available():\n        return True\n"

# Appended to the stub ``cascadia/cascadia.py``: write ``<-o>.ssl`` (header only), as M10-08's
# ``_ready`` does, so a finished run has its output.
_SSL_WRITER = """

_stub_main = main


def main():
    import sys
    _stub_main()
    a = sys.argv[1:]
    with open(a[a.index("-o") + 1] + ".ssl", "w") as fh:
        fh.write("file\\tscan\\tcharge\\tsequence\\tscore-type\\tscore\\tretention-time\\n")
"""

_CASCADIA_BIN = """#!{python}
import sys
if "--help" in sys.argv[1:]:
    print("usage: cascadia sequence SPECTRA MODEL [-o OUT] [-b BATCH] [-t THRESHOLD]")
    sys.exit(0)
from cascadia.cascadia import main
main()
"""


def _sandbox(tmp: Path) -> A3Sandbox:
    template = write_mzml(tmp / "template.mzML", example11_spectra())
    box = a3_sandbox(tmp / "sb", template)
    env = box.stub_cascadia_env()
    site = subprocess.run(
        [str(env / "bin" / "python"), "-I", "-c",
         "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()  # fmt: skip
    (Path(site) / "torch").mkdir()
    (Path(site) / "torch" / "__init__.py").write_text(_TORCH)
    with (Path(site) / "cascadia" / "cascadia.py").open("a") as fh:
        fh.write(_SSL_WRITER)
    write_exe(
        env / "bin" / "cascadia", _CASCADIA_BIN.replace("{python}", str(env / "bin" / "python"))
    )
    write_exe(box.sb.bin / "nvidia-smi", _NVIDIA_SMI)

    mzml = box.data / "mzml" / "WT_1.mzML"
    mzml.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(template, mzml)
    ckpt = box.db / "cascadia" / "cascadia_0.0.7.ckpt"
    unimod = box.db / "cascadia" / "unimod_tables.xml"
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    ckpt.write_bytes(b"PK\x03\x04 stub checkpoint\n")
    unimod.write_bytes(b"<unimod/>\n")

    def row(asset: str, path: Path, rel: str) -> dict[str, str]:
        data = path.read_bytes()
        return {"asset": asset, "version": "stub", "url": "https://example.invalid/x",
                "sha256": hashlib.sha256(data).hexdigest(), "bytes": str(len(data)),
                "licence": "Apache-2.0", "path": rel, "retrieved": "2026-10-24"}  # fmt: skip

    box.write_assets([
        row("cascadia_weights", ckpt, "$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt"),
        row("unimod_table", unimod, "$GIAP_DB_DIR/cascadia/unimod_tables.xml"),
    ])  # fmt: skip
    box.extra_env.update({"STUB_LOG": str(tmp / "stub_cascadia.json"), "A3_RUN": "WT_1"})
    return box


def _o_arguments(box: A3Sandbox) -> list[str]:
    """The ``-o`` value of every Cascadia ``sequence`` call (through the launcher or not)."""
    outs = []
    for call in box.calls("micromamba.jsonl"):
        argv = call["argv"]
        if "sequence" in argv and "-o" in argv:
            outs.append(argv[argv.index("-o") + 1])
    return outs


def _msg(res) -> str:  # type: ignore[no-untyped-def]
    return (
        f"rc={res.returncode}\n--- stdout\n{res.stdout[-3000:]}\n--- stderr\n{res.stderr[-3000:]}"
    )


@pytest.fixture
def box(tmp_path: Path) -> A3Sandbox:
    return _sandbox(tmp_path)


def test_a3_out_dir_set_writes_under_it(box: A3Sandbox) -> None:
    a3 = box.root / "pilot" / "A3"
    res = box.run("denovo_pilot", "cascadia", A3_OUT_DIR=str(a3))
    outs = _o_arguments(box)
    assert len(outs) == 1, f"one Cascadia call for A3_RUN=WT_1, got {outs}\n{_msg(res)}"
    assert outs[0].startswith(f"{a3}/WT_1/WT_1"), outs[0]
    assert res.returncode == 0, _msg(res)


def test_a3_out_dir_unset_keeps_the_m10_08_default(box: A3Sandbox) -> None:
    res = box.run("denovo_pilot", "cascadia", unset=("A3_OUT_DIR",))
    outs = _o_arguments(box)
    assert len(outs) == 1, f"one Cascadia call for A3_RUN=WT_1, got {outs}\n{_msg(res)}"
    assert outs[0].startswith(f"{box.root}/pilot/A3_cascadia/WT_1/WT_1"), outs[0]
    assert res.returncode == 0, _msg(res)


def test_the_m10_08_smoke_output_is_left_alone(box: A3Sandbox) -> None:
    smoke = box.root / "pilot" / "A3_cascadia" / "WT_1" / "WT_1.ssl"
    smoke.parent.mkdir(parents=True)
    smoke.write_text("file\tscan\tcharge\tsequence\tscore-type\tscore\n")
    before = {p: p.read_bytes() for p in smoke.parent.parent.rglob("*") if p.is_file()}
    a3 = box.root / "pilot" / "A3"
    res = box.run("denovo_pilot", "cascadia", A3_OUT_DIR=str(a3))
    assert _o_arguments(box), _msg(res)
    after = {p: p.read_bytes() for p in smoke.parent.parent.rglob("*") if p.is_file()}
    assert after == before, "an M10-10 run changed A3_cascadia/ (M10-08's verdict evidence)"


@pytest.mark.parametrize("run", ["WT_2", "WT_3"])
def test_a3_out_dir_names_each_run(box: A3Sandbox, run: str) -> None:
    mzml = box.data / "mzml"
    shutil.copyfile(mzml / "WT_1.mzML", mzml / f"{run}.mzML")
    a3 = box.root / "pilot" / "A3"
    res = box.run("denovo_pilot", "cascadia", A3_OUT_DIR=str(a3), A3_RUN=run)
    outs = _o_arguments(box)
    assert len(outs) == 1, f"{outs}\n{_msg(res)}"
    assert outs[0].startswith(f"{a3}/{run}/{run}"), outs[0]
