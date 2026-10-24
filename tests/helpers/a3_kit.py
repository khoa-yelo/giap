"""Synthetic mzML, acquisition stubs, stub packages and a script sandbox for the M10-08a tests.

Spec: docs/specs/M10-08-a3-feasibility.md (part 08a: env files, ``PILOT_TOOLS``, the offline
Unimod launcher, ``a3_assets.tsv``, ``fetch_a3_assets.sbatch``, ``mzml_pilot.sbatch`` and the
``mzml-check`` subcommand of ``validation/pilot/a3_feasibility.py``).

Nothing here reaches the network, a GPU, real raw data, real Cascadia, TIMSCONVERT or
ProteoWizard, or the real ``$GIAP_ROOT``:

* :func:`write_mzml` writes a small, hand-made mzML 1.1 file (optionally wrapped in
  ``indexedmzML``). Binary arrays are deliberately **not** valid base64 (``mzml-check`` must never
  decode them; it reads ``defaultArrayLength`` only).
* :func:`write_acquisition` writes an ``acquisition.json`` with the M10-01 layout
  (``{"dataset_mode", "runs": [{run, msms_type_counts, dia_window_rows, window_groups, mz_range,
  gradient_length_s, ...}]}``).
* :func:`write_stub_packages` writes stub ``cascadia`` / ``pyteomics`` / ``psims`` packages.
* :func:`a3_sandbox` builds on ``helpers.pilot_kit.pilot_sandbox`` (throw-away checkout, fake
  Slurm tools, fake ``curl``/``wget``, an unreachable proxy) and adds fake ``micromamba``,
  ``apptainer``/``singularity``, ``df``, ``sha256sum`` and ``shasum``. Every fake is either Python
  (shebang = this interpreter) or bash 3.2-compatible shell.

Conventions settled by the test author (the spec leaves them open; see the report to the lead):

* ``mzml-check`` writes ``J`` as JSON, sorted keys, 2-space indent, trailing newline, with at least
  ``rules`` (``{"M1": bool, …, "M8": bool}``), ``pass`` (all rules true), ``failed`` (sorted names
  of the false rules) and ``empty_ms2_share`` (M8's reported share). ``J`` is written on a pass and
  on a rule failure (exit 1).
* ``validation/pilot/a3_feasibility.py`` has ``main(argv: Sequence[str] | None = None) -> int``
  like the other ``validation/pilot`` scripts (used by the Hypothesis properties in-process).
* ``mzml_pilot.sbatch`` reads the B1 image path from ``PWIZ_SIF``.
* ``fetch_a3_assets.sbatch`` reads the expected sha256 values from
  ``validation/pilot/a3_assets.tsv`` in ``$GIAP_REPO`` and finds psims' vendored
  ``unimod_tables.xml.gz`` inside the built ``giap-cascadia`` env
  (``$MAMBA_ROOT_PREFIX/envs/giap-cascadia``); the sandbox env is a real ``venv`` (no pip)
  holding a stub ``psims`` package, so both ``find``-based and
  ``micromamba run -n giap-cascadia python -c 'import psims …'`` lookups work.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import venv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from helpers.pilot_kit import PilotSandbox, pilot_sandbox
from helpers.pxd_fixtures import _DF
from helpers.slurm_sandbox import REPO, write_exe

PILOT = REPO / "validation" / "pilot"
A3_FEASIBILITY = PILOT / "a3_feasibility.py"
A3_ASSETS = PILOT / "a3_assets.tsv"
ACQUISITION = PILOT / "inputs" / "acquisition.json"
LAUNCHER = REPO / "scripts" / "pilot" / "cascadia_offline.py"
SLURM = REPO / "scripts" / "slurm"
FETCH_SBATCH = SLURM / "fetch_a3_assets.sbatch"
MZML_SBATCH = SLURM / "mzml_pilot.sbatch"
ENVS = REPO / "workflow" / "envs"

RULES = tuple(f"M{i}" for i in range(1, 9))
ASSET_COLUMNS = ("asset", "version", "url", "sha256", "bytes", "licence", "path", "retrieved")

CASCADIA_VERSION = "0.0.7"
CASCADIA_WHEEL_SHA_HEAD, CASCADIA_WHEEL_SHA_TAIL = "50b52043", "f092a3"
TIMSCONVERT_COMMIT = "3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b"
PYTDFSDK_COMMIT = "4808da81d3e37f82da96f392dbb446c3733edacb"
PWIZ_DIGEST = "sha256:34ed787cd4c3182c2d81adcbf7514031e594609ca460abc5d3d114dc2334dae8"
WEIGHTS_URL_PREFIX = "https://drive.usercontent.google.com/download?id="

CENTROID, PROFILE = "MS:1000127", "MS:1000128"

# ------------------------------------------------------------------ files and assertions


def require(path: Path, what: str) -> Path:
    """Fail (not skip) with a clear message while an M10-08a deliverable is missing."""
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist yet (M10-08a: {what})"
    return path


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def is_canonical_json(text: str) -> bool:
    """Sorted keys, 2-space indent, trailing newline (the M10-08 JSON byte convention)."""
    return text == json.dumps(json.loads(text), indent=2, sort_keys=True) + "\n"


# ------------------------------------------------------------------ synthetic mzML


@dataclass(frozen=True)
class Spec:
    """One spectrum of a synthetic mzML file.

    ``windows``: the isolation windows of the MS2 spectrum's precursors, one precursor per entry,
    as (target m/z, lower offset, upper offset); ``None`` in an entry writes a precursor without an
    ``isolationWindow``. ``mode``: ``"centroid"``, ``"profile"`` or ``None`` (no flag).
    """

    ms_level: int
    rt: float
    windows: tuple[tuple[float, float, float] | None, ...] = ()
    mode: str | None = "centroid"
    n_peaks: int = 5


def _cv(acc: str, name: str, value: object = "", unit: tuple[str, str] | None = None) -> str:
    u = ""
    if unit is not None:
        u = f' unitCvRef="UO" unitAccession="{unit[0]}" unitName="{unit[1]}"'
    ref = acc.split(":")[0]
    return f'<cvParam cvRef="{ref}" accession="{acc}" name="{name}" value="{value}"{u}/>'


_MZ_UNIT = ("MS:1000040", "m/z")


def _spectrum_xml(i: int, s: Spec, rt_unit: str) -> str:
    parts = [
        f'<spectrum index="{i}" id="scan={i + 1}" defaultArrayLength="{s.n_peaks}">',
        _cv("MS:1000511", "ms level", s.ms_level),
        _cv("MS:1000579", "MS1 spectrum") if s.ms_level == 1 else _cv("MS:1000580", "MSn spectrum"),
    ]
    if s.mode == "centroid":
        parts.append(_cv(CENTROID, "centroid spectrum"))
    elif s.mode == "profile":
        parts.append(_cv(PROFILE, "profile spectrum"))
    rt = s.rt if rt_unit == "second" else s.rt / 60.0
    unit = ("UO:0000010", "second") if rt_unit == "second" else ("UO:0000031", "minute")
    parts += [
        '<scanList count="1">',
        _cv("MS:1000795", "no combination"),
        "<scan>",
        _cv("MS:1000016", "scan start time", repr(rt), unit),
        "</scan>",
        "</scanList>",
    ]
    if s.ms_level == 2 and s.windows:
        parts.append(f'<precursorList count="{len(s.windows)}">')
        for w in s.windows:
            parts.append("<precursor>")
            if w is not None:
                t, lo, up = w
                parts += [
                    "<isolationWindow>",
                    _cv("MS:1000827", "isolation window target m/z", t, ("MS:1000040", "m/z")),
                    _cv("MS:1000828", "isolation window lower offset", lo, ("MS:1000040", "m/z")),
                    _cv("MS:1000829", "isolation window upper offset", up, ("MS:1000040", "m/z")),
                    "</isolationWindow>",
                ]
            parts += [
                '<selectedIonList count="1"><selectedIon>',
                _cv("MS:1000744", "selected ion m/z", w[0] if w else 500.0, _MZ_UNIT),
                "</selectedIon></selectedIonList>",
                "<activation>",
                _cv("MS:1000133", "collision-induced dissociation"),
                "</activation>",
                "</precursor>",
            ]
        parts.append("</precursorList>")
    # Not valid base64 on purpose: a checker that decodes binary arrays fails on this file.
    junk = "!!not-base64: mzml-check must read defaultArrayLength only!!"
    parts += [
        '<binaryDataArrayList count="2">',
        f'<binaryDataArray encodedLength="{len(junk)}">',
        _cv("MS:1000523", "64-bit float"),
        _cv("MS:1000574", "zlib compression"),
        _cv("MS:1000514", "m/z array", "", ("MS:1000040", "m/z")),
        f"<binary>{junk}</binary>",
        "</binaryDataArray>",
        f'<binaryDataArray encodedLength="{len(junk)}">',
        _cv("MS:1000521", "32-bit float"),
        _cv("MS:1000574", "zlib compression"),
        _cv("MS:1000515", "intensity array", "", ("MS:1000131", "number of detector counts")),
        f"<binary>{junk}</binary>",
        "</binaryDataArray>",
        "</binaryDataArrayList>",
        "</spectrum>",
    ]
    return "".join(parts)


def write_mzml(
    path: Path, spectra: Sequence[Spec], *, indexed: bool = False, rt_unit: str = "second"
) -> Path:
    """Write ``spectra`` as a small mzML 1.1 file (namespace ``http://psi.hupo.org/ms/mzml``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    head = (
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0" id="synthetic">'
        '<cvList count="2">'
        '<cv id="MS" fullName="PSI-MS" URI="https://purl.obolibrary.org/obo/ms.obo"/>'
        '<cv id="UO" fullName="Unit Ontology" URI="https://purl.obolibrary.org/obo/uo.obo"/>'
        "</cvList>"
        "<fileDescription><fileContent>"
        + _cv("MS:1000579", "MS1 spectrum")
        + _cv("MS:1000580", "MSn spectrum")
        + "</fileContent></fileDescription>"
        '<softwareList count="1"><software id="giap_test" version="0"/></softwareList>'
        '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1"/>'
        "</instrumentConfigurationList>"
        '<dataProcessingList count="1"><dataProcessing id="dp"/></dataProcessingList>'
        '<run id="synthetic" defaultInstrumentConfigurationRef="IC1">'
        f'<spectrumList count="{len(spectra)}" defaultDataProcessingRef="dp">'
    )
    with path.open("w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="utf-8"?>\n')
        if indexed:
            fh.write(
                '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            )
        fh.write(head)
        for i, s in enumerate(spectra):
            fh.write(_spectrum_xml(i, s, rt_unit))
            fh.write("\n")
        fh.write("</spectrumList></run></mzML>")
        if indexed:
            fh.write('<indexList count="1"><index name="spectrum"></index></indexList>')
            fh.write("<indexListOffset>0</indexListOffset></indexedmzML>")
        fh.write("\n")
    return path


# Example 11: 4 windows of width 25 over m/z 400-500; 3 cycles of MS1 + 4 MS2; last RT 9.5 s.
EX11_WINDOWS = tuple((400.0 + 25.0 * k + 12.5, 12.5, 12.5) for k in range(4))
_EX11_RT = ((0.5, (1.0, 1.5, 2.0, 2.5)), (3.5, (4.0, 4.5, 5.0, 5.5)), (6.5, (8.0, 8.5, 9.0, 9.5)))


def example11_spectra() -> list[Spec]:
    out: list[Spec] = []
    for ms1_rt, ms2_rts in _EX11_RT:
        out.append(Spec(1, ms1_rt))
        out += [Spec(2, rt, (w,)) for rt, w in zip(ms2_rts, EX11_WINDOWS, strict=True)]
    return out


EX11_ACQ: dict[str, Any] = {
    "msms_type_counts": {"0": 3, "9": 3},
    "window_groups": 2,
    "dia_window_rows": 4,
    "mz_range": [400.0, 500.0],
    "gradient_length_s": 10.0,
}


def ms2_indices(spectra: Sequence[Spec]) -> list[int]:
    return [i for i, s in enumerate(spectra) if s.ms_level == 2]


def with_window(spectra: list[Spec], idx: int, *windows: tuple[float, float, float] | None) -> None:
    spectra[idx] = replace(spectra[idx], windows=tuple(windows))


def write_acquisition(path: Path, runs: dict[str, dict[str, Any]]) -> Path:
    """An ``acquisition.json`` in the M10-01 layout with the given runs."""
    rows = []
    for run, vals in runs.items():
        rows.append({
            "run": run, "source_folder": f"{run}.d", "notes": [], "mode": "diaPASEF",
            "pasef_rows": 0, "metadata": {}, "incomplete_run": False, **vals,
        })  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"dataset_mode": "diaPASEF", "runs": rows}, indent=2) + "\n")
    return path


def wt1_spectra(n_ms2: int = 8804) -> list[Spec]:
    """WT_1-shaped spectra: 734 MS1, ``n_ms2`` MS2 over 24 windows of width 25, last RT ≈ 528 s.

    Checked against the committed ``acquisition.json`` (WT_1: 734 MS1 frames, 8,804 diaPASEF
    frames, 12 window groups, 65 window rows, m/z 343.83-1188.83, 539.86 s).
    """
    windows = [(400.0 + 25.0 * k, 12.5, 12.5) for k in range(24)]
    out: list[Spec] = []
    left = n_ms2
    for c in range(734):
        t0 = 0.72 * c
        out.append(Spec(1, round(t0, 3), n_peaks=50))
        take = min(12, left) if c < 733 else left  # 733 x 12 + 8 = 8,804
        for j in range(take):
            w = windows[(12 * (c % 2)) + (j % 12)]
            out.append(Spec(2, round(t0 + 0.05 * (j + 1), 3), (w,), n_peaks=20))
        left -= take
    assert sum(1 for s in out if s.ms_level == 2) == n_ms2, "wt1_spectra: MS2 count"
    return out


# ------------------------------------------------------------------ stub Python packages

_CASCADIA_INIT = """import os
_log = os.environ.get("STUB_IMPORT_LOG")
if _log:
    with open(_log, "a") as fh:
        fh.write("cascadia imported\\n")
"""

_CASCADIA_MAIN = """import json, os, sys
import pyteomics.mass.unimod as _u
URL_AT_IMPORT = _u._unimod_xml_download_url


def main():
    import pyteomics.mass.unimod as u
    rec = {"argv": list(sys.argv[1:]), "url_in_main": u._unimod_xml_download_url,
           "url_at_import": URL_AT_IMPORT}
    with open(os.environ["STUB_LOG"], "w") as fh:
        json.dump(rec, fh)
    rc = int(os.environ.get("STUB_RC", "0"))
    if rc:
        sys.exit(rc)
"""

DEFAULT_UNIMOD_URL = "http://www.unimod.org/xml/unimod_tables.xml"

STUB_UNIMOD_XML = (
    b'<?xml version="1.0" encoding="UTF-8"?>\n'
    b'<unimod xmlns="http://www.unimod.org/xmlns/schema/unimod_tables_1">'
    b'<modifications><modifications_row record_id="35" code_name="Oxidation"/>'
    b"</modifications></unimod>\n"
)


def write_stub_packages(site: Path, *, unimod_gz: bytes | None = None) -> Path:
    """Stub ``cascadia``, ``pyteomics`` and ``psims`` packages under ``site``."""
    pkgs = {
        "cascadia/__init__.py": _CASCADIA_INIT,
        "cascadia/cascadia.py": _CASCADIA_MAIN,
        "pyteomics/__init__.py": "",
        "pyteomics/mass/__init__.py": "",
        "pyteomics/mass/unimod.py": f'_unimod_xml_download_url = "{DEFAULT_UNIMOD_URL}"\n',
        "psims/__init__.py": '__version__ = "1.4.0"\n',
        "psims/controlled_vocabulary/__init__.py": "",
        "psims/controlled_vocabulary/vendor/__init__.py": "",
    }
    for rel, text in pkgs.items():
        p = site / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    gz = site / "psims" / "controlled_vocabulary" / "vendor" / "unimod_tables.xml.gz"
    gz.write_bytes(unimod_gz if unimod_gz is not None else gzip.compress(STUB_UNIMOD_XML, mtime=0))
    return site


def make_stub_venv(prefix: Path) -> Path:
    """A real venv (no pip) at ``prefix`` with the stub packages in its site-packages."""
    venv.EnvBuilder(with_pip=False, symlinks=True, clear=True).create(prefix)
    out = subprocess.run(
        [str(prefix / "bin" / "python"), "-I", "-c",
         "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    write_stub_packages(Path(out.stdout.strip()))
    return prefix


# ------------------------------------------------------------------ fakes for the sandbox

_MICROMAMBA = r"""#!{python}
'''Fake micromamba: logs every call; ``run -n ENV CMD …`` emulates timsconvert or runs CMD in
the fake env prefix ($MAMBA_ROOT_PREFIX/envs/ENV, bin/ first on PATH). Nothing is solved.'''
import json, os, shutil, sys

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "micromamba.jsonl"), "a") as fh:
    fh.write(json.dumps({"argv": argv, "cwd": os.getcwd()}) + "\n")
if not argv or argv[0] != "run":
    sys.exit(0)
rest, env = argv[1:], None
while rest and rest[0].startswith("-"):
    if rest[0] in ("-n", "--name", "-p", "--prefix") and len(rest) > 1:
        env = rest[1]; rest = rest[2:]; continue
    rest = rest[1:]
root = os.environ.get("MAMBA_ROOT_PREFIX", "")
prefix = env if env and env.startswith("/") else os.path.join(root, "envs", env or "base")
if not rest:
    sys.exit(0)
if os.path.basename(rest[0]) == "timsconvert":
    args = rest[1:]
    if "--help" in args or "--version" in args or not args:
        print("usage: timsconvert --input INPUT --outdir OUTDIR [--mode {raw,centroid,profile}] "
              "[--compression {zlib,none}] [--exclude_mobility]")
        sys.exit(0)
    def val(flag):
        return args[args.index(flag) + 1] if flag in args else None
    src, outdir, outfile = val("--input"), val("--outdir"), val("--outfile")
    if src is None or not os.path.isdir(src):
        print("timsconvert: no such .d: %s" % src, file=sys.stderr); sys.exit(1)
    outdir = outdir or os.path.dirname(os.path.abspath(src))
    os.makedirs(outdir, exist_ok=True)
    name = outfile or os.path.splitext(os.path.basename(src.rstrip("/")))[0] + ".mzML"
    shutil.copyfile(os.environ["FAKE_MZML_TEMPLATE"], os.path.join(outdir, name))
    sys.exit(int(os.environ.get("FAKE_TIMSCONVERT_RC", "0")))
env_vars = dict(os.environ, CONDA_PREFIX=prefix,
                PATH=os.path.join(prefix, "bin") + os.pathsep + os.environ.get("PATH", ""))
cmd = rest
if not os.path.sep in cmd[0]:
    cand = os.path.join(prefix, "bin", cmd[0])
    if os.path.exists(cand):
        cmd = [cand] + cmd[1:]
os.execvpe(cmd[0], cmd, env_vars)
"""

_APPTAINER = r"""#!{python}
'''Fake apptainer/singularity: logs the call; ``exec … IMAGE wine msconvert …`` writes the
template mzML to ``-o DIR`` / ``--outfile NAME``; ``pull`` writes a stub image. Nothing runs.'''
import json, os, shutil, sys

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "apptainer.jsonl"), "a") as fh:
    fh.write(json.dumps({"tool": os.path.basename(sys.argv[0]), "argv": argv}) + "\n")
if not argv or argv[0] in ("--version", "version"):
    print("apptainer version 1.4.4"); sys.exit(0)
if argv[0] == "pull":
    targets = [a for a in argv[1:] if not a.startswith("-") and not "://" in a]
    if targets:
        d = os.path.dirname(targets[0])
        if d:
            os.makedirs(d, exist_ok=True)
        open(targets[0], "w").write("stub sif\n")
    sys.exit(0)
if "msconvert" in argv:
    inner = argv[argv.index("msconvert") + 1:]
    if "--help" in inner or not inner:
        print("msconvert (stub)"); sys.exit(0)
    outdir = inner[inner.index("-o") + 1] if "-o" in inner else os.getcwd()
    name = inner[inner.index("--outfile") + 1] if "--outfile" in inner else "out.mzML"
    os.makedirs(outdir, exist_ok=True)
    shutil.copyfile(os.environ["FAKE_MZML_TEMPLATE"], os.path.join(outdir, name))
sys.exit(0)
"""

_SHA256SUM = r"""#!{python}
'''Fake sha256sum/shasum (portable: macOS has no sha256sum). Computes the real digest.'''
import hashlib, sys

args = [a for a in sys.argv[1:] if a not in ("-b", "--binary", "-t", "--text")]
if args[:2] == ["-a", "256"]:
    args = args[2:]
args = [a for a in args if not (a.startswith("-") and a != "-")] or ["-"]
for a in args:
    data = sys.stdin.buffer.read() if a == "-" else open(a, "rb").read()
    print(hashlib.sha256(data).hexdigest() + "  " + a)
"""


@dataclass
class A3Sandbox:
    pilot: PilotSandbox
    template: Path
    extra_env: dict[str, str] = field(default_factory=dict)

    @property
    def sb(self):  # type: ignore[no-untyped-def]
        return self.pilot.sb

    @property
    def root(self) -> Path:
        return self.sb.root

    @property
    def data(self) -> Path:
        return self.sb.root / "data" / "pxd078459"

    @property
    def db(self) -> Path:
        return self.sb.root / "db"

    @property
    def mamba(self) -> Path:
        return self.sb.root / "micromamba"

    def run(
        self, script: str, *args: str, unset: tuple[str, ...] = (), **env: str
    ) -> subprocess.CompletedProcess[str]:
        """Run ``scripts/slurm/<script>.sbatch`` with bash in the sandbox (never through Slurm)."""
        path = self.sb.repo / "scripts" / "slurm" / f"{script}.sbatch"
        assert path.is_file(), f"scripts/slurm/{script}.sbatch does not exist yet (M10-08a)"
        from helpers.slurm_sandbox import assert_honours_giap_repo

        assert_honours_giap_repo(path.read_text(), path.name)
        return self.sb.run(["bash", str(path), *args], env={**self.extra_env, **env},
                           cwd=self.sb.tmp, unset=unset)  # fmt: skip

    def calls(self, name: str) -> list[dict[str, Any]]:
        return self.sb.jsonl(name)

    def urls(self) -> list[str]:
        return [u for c in self.sb.jsonl("fetch.jsonl") for u in c["urls"]]

    def make_run_dir(self, run: str = "WT_1", tdf_bin_bytes: int = 100_000) -> Path:
        d = self.data / f"{run}.d"
        d.mkdir(parents=True, exist_ok=True)
        (d / "analysis.tdf").write_bytes(b"stub tdf\n")
        (d / "analysis.tdf_bin").write_bytes(b"\0" * tdf_bin_bytes)
        return d

    def files_under_root(self) -> set[Path]:
        skip = {self.sb.root / "logs", self.mamba, self.sb.root / "venv"}
        out = set()
        for p in self.sb.root.rglob("*"):
            if p.is_file() and not any(s in p.parents for s in skip):
                out.add(p)
        return out

    def write_assets(self, rows: Iterable[dict[str, str]]) -> Path:
        path = self.sb.repo / "validation" / "pilot" / "a3_assets.tsv"
        lines = ["\t".join(ASSET_COLUMNS)]
        lines += ["\t".join(r.get(c, "") for c in ASSET_COLUMNS) for r in rows]
        path.write_text("\n".join(lines) + "\n")
        return path

    def stub_cascadia_env(self) -> Path:
        return make_stub_venv(self.mamba / "envs" / "giap-cascadia")


def a3_sandbox(tmp: Path, template: Path) -> A3Sandbox:
    pilot = pilot_sandbox(tmp, sif=False, raws=0)
    sb = pilot.sb
    src = sb.repo / "src"
    if not src.exists():
        src.symlink_to(REPO / "src")  # mzml-check imports giap; untracked, ignored below
        with (sb.repo / ".git" / "info" / "exclude").open("a") as fh:
            fh.write("/src\n")
    py = sys.executable
    write_exe(sb.bin / "micromamba", _MICROMAMBA.replace("{python}", py))
    for name in ("apptainer", "singularity"):
        write_exe(sb.bin / name, _APPTAINER.replace("{python}", py))
    write_exe(sb.bin / "df", _DF.replace("{python}", py))
    for name in ("sha256sum", "shasum"):
        write_exe(sb.bin / name, _SHA256SUM.replace("{python}", py))
    (sb.root / "micromamba" / "envs").mkdir(parents=True, exist_ok=True)
    tmpdir = sb.tmp / "nodetmp"
    tmpdir.mkdir()
    sb.env.update({
        "MAMBA_ROOT_PREFIX": str(sb.root / "micromamba"),
        "FAKE_MZML_TEMPLATE": str(template),
        "TMPDIR": str(tmpdir),
        "PYTHONPATH": str(REPO / "src"),
    })  # fmt: skip
    sb.env.pop("DIANN_SIF", None)
    return A3Sandbox(pilot=pilot, template=template)


def copy_tree(src: Path, dst: Path) -> None:
    shutil.copytree(src, dst, symlinks=True, dirs_exist_ok=True)
