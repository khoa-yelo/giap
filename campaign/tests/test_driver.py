"""Campaign driver: ``plan`` and ``run_job`` (campaign/SPEC.md §10).

Every tool and every download is a fake: ``FakeRun`` records the argv the driver builds with the
real pure builders of ``engines.py`` and writes the output file the real tool would write;
``FakeDownload`` writes the bytes a real download would leave at ``dest``. The genome package is a
real NCBI-style zip and the Bruker runs are real zips holding ``<run>.d/analysis.tdf`` (an sqlite
file with a ``Frames`` table), so ``safe_extract``, ``find_genome_files``, ``gbff_topologies``,
``tag_fasta`` and ``tdf_acquisition`` run for real.

Readings of the spec (stated so the implementer can check them):

* ``run_job(job, *, db_dir, profile, threads, keep_inputs, run, download, env)``: ``job`` is a
  mapping with the §10 ``plan`` columns (as one row of ``jobs.tsv``); ``run(argv, cwd)`` returns
  the exit code; ``download(url, dest, expected_bytes, checksum)`` writes ``dest`` (whose parent
  directory the driver creates) and returns its SHA-256.
* ``env`` keys: ``db_dir``, ``diann_sif``, ``cascadia_model``, ``cascadia_unimod``,
  ``casanovo_model``, ``diann_cfg`` (paths) and ``pins`` (name -> expected SHA-256 hex, ``""`` =
  no pin; the CLI fills it from ``campaign/pins.tsv``).
* Pins: the files a job needs follow from its acquisition (DIA: ``diann_sif``, ``cascadia_model``;
  DDA: ``casanovo_model``). A needed file that is missing, has no pin, an empty pin or a different
  SHA-256 makes ``run_job`` return 2 (environment error, §12) before any search or de novo engine
  runs; files the job does not need are not checked.
* A tool that exits 0 without writing its step's output fails the step.
* ``status.json["acquisition"]`` is a plain string mode (``"DIA"``, ``"DDA"``, ``"mixed"`` or
  ``"unknown"``) and the ``resolve_acquisition`` note is the top-level string
  ``"acquisition_note"`` (``""`` when none) (lead decision, aligned with the aggregate tests).
* Throttle with every ``run_gb`` 0 (sizes unknown): the quotient is unbounded, so ``min(8, …)``
  gives 8.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib
import io
import json
import math
import sqlite3
import zipfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from giap.io.converters.base import ConvertResult
from giap.io.issues import Issue
from giap_campaign import driver
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

ACC = "GCF_000005845.2"
DATASET = "PXD000001"
TAX_ID = 511145
ORGANISM = "Escherichia coli str. K-12 substr. MG1655"
PRIDE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD000001"
PROFILE = "farmshare_local,apptainer"
THREADS = 4

STEPS = (
    "fetch_genome",
    "annotate",
    "makedb",
    "fetch_raw",
    "acquisition",
    "to_mzml",
    "search_db",
    "denovo",
    "convert",
    "giap_db",
    "giap_denovo",
    "compare",
)
PLAN_COLUMNS = [
    "job",
    "dataset",
    "tax_id",
    "organism",
    "genome_accession",
    "raw_format",
    "acquisition_metadata",
    "run_files",
    "run_urls",
    "run_bytes",
    "run_checksums",
    "run_gb",
    "denovo_runs",
    "workdir",
]
SEARCH_KINDS = {"diann1", "diann2", "cascadia", "sage", "casanovo"}

DIANN_CFG = """\
# test settings: one flag per line
[step1]
--fasta-search
--gen-spec-lib
--predictor
--cut K*,R*
--missed-cleavages 1

[step2]
--qvalue 0.01
--matrices
"""

PEPTIDES_TSV = "peptide\tq_value\trun\nPEPTIDEK\t0.001\tA\n"


# ---------------------------------------------------------------- input builders


def _gbff(replicons: Sequence[tuple[str, str]]) -> str:
    out = []
    for acc_ver, topology in replicons:
        name = acc_ver.split(".")[0]
        out.append(
            f"LOCUS       {name:<16} 1000 bp    DNA     {topology:<8} BCT 01-JAN-2020\n"
            "DEFINITION  test replicon.\n"
            f"ACCESSION   {name}\n"
            f"VERSION     {acc_ver}\n"
            "ORIGIN\n"
            "        1 acgtacgtac\n"
            "//\n"
        )
    return "".join(out)


def _fna(replicons: Sequence[tuple[str, str]]) -> str:
    return "".join(f">{acc_ver} test replicon\nACGTACGTACGTACGT\nACGT\n" for acc_ver, _ in replicons)


CIRCULAR = (("NC_000913.3", "circular"),)
WITH_LINEAR = (("NC_000913.3", "circular"), ("NZ_CP009001.1", "linear"))


def genome_zip(replicons: Sequence[tuple[str, str]] = CIRCULAR) -> bytes:
    """An NCBI Datasets package as ``/genome/accession/{acc}/download`` returns it."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("README.md", "NCBI Datasets package\n")
        zf.writestr("ncbi_dataset/data/assembly_data_report.jsonl", json.dumps({"a": ACC}) + "\n")
        zf.writestr("ncbi_dataset/data/dataset_catalog.json", "{}\n")
        zf.writestr(f"ncbi_dataset/data/{ACC}/{ACC}_ASM584v2_genomic.fna", _fna(replicons))
        zf.writestr(f"ncbi_dataset/data/{ACC}/genomic.gbff", _gbff(replicons))
    return buf.getvalue()


def tdf_bytes(scratch: Path, msms_types: Iterable[int]) -> bytes:
    """A Bruker ``analysis.tdf`` (sqlite) whose Frames hold MS1 frames and ``msms_types``."""
    scratch.mkdir(parents=True, exist_ok=True)
    path = scratch / f"analysis_{len(list(scratch.iterdir()))}.tdf"
    with contextlib.closing(sqlite3.connect(path)) as con:
        con.execute(
            "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, Polarity TEXT, "
            "ScanMode INTEGER, MsMsType INTEGER, NumScans INTEGER)"
        )
        rows = []
        for i, t in enumerate([0, *msms_types, 0, *msms_types]):
            rows.append((i + 1, 0.1 * i, "+", 9 if t == 9 else 8, t, 927))
        con.executemany("INSERT INTO Frames VALUES (?, ?, ?, ?, ?, ?)", rows)
        con.commit()
    return path.read_bytes()


def bruker_zip(scratch: Path, run: str, msms_types: Iterable[int]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(f"{run}.d/analysis.tdf", tdf_bytes(scratch, msms_types))
        zf.writestr(f"{run}.d/analysis.tdf_bin", b"\0" * 64)
    return buf.getvalue()


def mzml_text(half_width: float, n_spectra: int = 4) -> str:
    """Minimal mzML: MS2 spectra with isolation offsets ``half_width`` each side."""
    spectra = []
    for i in range(n_spectra):
        spectra.append(
            f'<spectrum index="{i}" id="controllerType=0 controllerNumber=1 scan={i + 1}" '
            'defaultArrayLength="0">\n'
            '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n'
            '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>\n'
            '<precursorList count="1"><precursor><isolationWindow>\n'
            '<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" '
            'value="500.0"/>\n'
            '<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" '
            f'value="{half_width}"/>\n'
            '<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" '
            f'value="{half_width}"/>\n'
            "</isolationWindow></precursor></precursorList>\n"
            "</spectrum>\n"
        )
    return (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">\n'
        f'<run id="r"><spectrumList count="{n_spectra}">\n'
        + "".join(spectra)
        + "</spectrumList></run>\n</mzML>\n"
    )


DDA_HALF_WIDTH = 0.7  # 1.4 Th windows
DIA_HALF_WIDTH = 12.5  # 25 Th windows


# ---------------------------------------------------------------- fakes


def _opt(argv: Sequence[str], flag: str) -> str:
    return argv[argv.index(flag) + 1]


def _opts(argv: Sequence[str], flag: str) -> list[str]:
    return [argv[i + 1] for i, tok in enumerate(argv[:-1]) if tok == flag]


def _at(cwd: Path, value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else Path(cwd) / p


def _same(a: str | Path, b: str | Path) -> bool:
    return Path(a).resolve() == Path(b).resolve()


def classify(argv: Sequence[str]) -> str:
    names = [Path(t).name for t in argv]
    if "bakta" in names:
        return "bakta"
    if any(t.endswith("diann-linux") for t in argv):
        return "diann2" if "--f" in argv else "diann1"
    for i, name in enumerate(names[:-1]):
        if name == "giap" and argv[i + 1] in {"makedb", "convert", "run", "compare"}:
            return f"giap_{argv[i + 1]}"
    if "casanovo" in names or "giap-casanovo" in names:
        return "casanovo"
    if ("cascadia" in names or "giap-cascadia" in names) and "sequence" in argv:
        return "cascadia"
    if "sage" in names or "giap-sage" in names:
        return "sage"
    if any("ThermoRawFileParser" in n for n in names):
        return "thermo"
    if "timsconvert" in names:
        return "timsconvert"
    return "unknown"


@dataclass
class Call:
    kind: str
    argv: list[str]
    cwd: Path
    all_str: bool
    info: dict[str, Any] = field(default_factory=dict)

    def opt(self, flag: str) -> str:
        return _opt(self.argv, flag)

    def path(self, flag: str) -> Path:
        return _at(self.cwd, self.opt(flag))


class FakeRun:
    """``run(argv, cwd) -> int``: records the call and writes what the real tool would."""

    def __init__(
        self,
        events: list[str],
        *,
        fail: Iterable[str] = (),
        no_output: Iterable[str] = (),
        thermo_half_width: float = DDA_HALF_WIDTH,
    ) -> None:
        self.events = events
        self.fail = set(fail)
        self.no_output = set(no_output)
        self.thermo_half_width = thermo_half_width
        self.calls: list[Call] = []

    def __call__(self, argv: Sequence[Any], cwd: str | Path, **_: Any) -> int:
        raw = list(argv)
        args = [str(a) for a in raw]
        call = Call(classify(args), args, Path(cwd), all(isinstance(a, str) for a in raw))
        self.calls.append(call)
        self.events.append(call.kind)
        if call.kind in self.fail:
            return 1
        if call.kind not in self.no_output:
            self._produce(call)
        return 0

    def of(self, kind: str) -> list[Call]:
        return [c for c in self.calls if c.kind == kind]

    def kinds(self) -> list[str]:
        return [c.kind for c in self.calls]

    def _produce(self, c: Call) -> None:  # noqa: C901 - one branch per tool
        a, cwd = c.argv, c.cwd
        if c.kind == "bakta":
            out = c.path("--output")
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{c.opt('--prefix')}.json").write_text('{"sequences": [], "features": []}')
        elif c.kind == "giap_makedb":
            out = c.path("--outdir")
            out.mkdir(parents=True, exist_ok=True)
            (out / "strain_db.faa").write_text(">p1\nMPEPTIDEK\n")
        elif c.kind == "thermo":
            raw, out = c.path("-i"), c.path("-o")
            c.info["input_exists"] = raw.is_file()
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{raw.stem}.mzML").write_text(mzml_text(self.thermo_half_width))
        elif c.kind == "timsconvert":
            d, out = c.path("--input"), c.path("--outdir")
            c.info["input_has_tdf"] = (d / "analysis.tdf").is_file()
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{d.name.removesuffix('.d')}.mzML").write_text(mzml_text(DIA_HALF_WIDTH))
        elif c.kind == "diann1":
            lib = c.path("--out-lib")
            lib.parent.mkdir(parents=True, exist_ok=True)
            lib.write_bytes(b"lib")
            lib.with_suffix(".predicted.speclib").write_bytes(b"predicted")
        elif c.kind == "diann2":
            spectra = [_at(cwd, f) for f in _opts(a, "--f")]
            c.info["f"] = spectra
            c.info["f_ok"] = [
                (p / "analysis.tdf").is_file() if p.name.endswith(".d") else p.is_file()
                for p in spectra
            ]
            out = c.path("--out")
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(b"PAR1 fake report PAR1")
            if "--out-lib" in a:
                c.path("--out-lib").write_bytes(b"empirical")
        elif c.kind == "cascadia":
            i = a.index("sequence")
            c.info["mzml"] = _at(cwd, a[i + 1])
            c.info["mzml_exists"] = c.info["mzml"].is_file()
            c.info["model"] = a[i + 2]
            prefix = c.path("-o")
            prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(f"{prefix}.ssl").write_text("file\tscan\tcharge\tsequence\tscore\n")
        elif c.kind == "casanovo":
            inputs = [_at(cwd, t) for t in a if t.lower().endswith((".mzml", ".mzml.gz"))]
            c.info["inputs"] = inputs
            c.info["inputs_exist"] = all(p.is_file() for p in inputs)
            out = c.path("--output")
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("MTD\tmzTab-version\t1.0.0\n")
        elif c.kind == "sage":
            cfg_path = _at(cwd, a[-1])
            cfg = json.loads(cfg_path.read_text())
            c.info["config"] = cfg
            c.info["mzml_exist"] = all(_at(cwd, p).is_file() for p in cfg["mzml_paths"])
            out = _at(cwd, cfg["output_directory"])
            out.mkdir(parents=True, exist_ok=True)
            (out / "results.sage.tsv").write_text("peptide\tproteins\tfilename\n")
        elif c.kind == "giap_convert":
            out = c.path("--out")
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(PEPTIDES_TSV)
        elif c.kind == "giap_run":
            out = c.path("--outdir")
            (out / "work" / "ab").mkdir(parents=True, exist_ok=True)
            (out / "work" / "ab" / "task.out").write_text("x")
            (out / "run_manifest.json").write_text('{"giap_version": "1.0.0"}')
        elif c.kind == "giap_compare":
            out = c.path("--outdir")
            (out / "tables").mkdir(parents=True, exist_ok=True)
            (out / "tables" / "compare_detection.tsv").write_text("peptide\n")


class FakeDownload:
    """``download(url, dest, expected_bytes, checksum) -> sha256``; writes the known bytes."""

    def __init__(self, events: list[str], genome: bytes, raw: dict[str, bytes]) -> None:
        self.events = events
        self.genome = genome
        self.raw = raw
        self.calls: list[dict[str, Any]] = []

    def __call__(
        self,
        url: str,
        dest: str | Path,
        expected_bytes: int | None = None,
        checksum: str | None = None,
        **_: Any,
    ) -> str:
        url = str(url)
        self.calls.append(
            {"url": url, "dest": Path(dest), "expected_bytes": expected_bytes, "checksum": checksum}
        )
        if "/genome/accession/" in url:
            self.events.append("download_genome")
            data = self.genome
        else:
            self.events.append("download_raw")
            data = self.raw[url]
        Path(dest).write_bytes(data)  # like curl -o: the parent directory must exist
        return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------- scenario


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_env(tmp_path: Path) -> tuple[Path, dict[str, Any]]:
    assets = tmp_path / "assets"
    assets.mkdir()
    files = {
        "diann_sif": assets / "diann_2.3.1.sif",
        "cascadia_model": assets / "cascadia_0.0.7.ckpt",
        "cascadia_unimod": assets / "unimod_tables.xml",
        "casanovo_model": assets / "casanovo_massivekb_5.2.1.ckpt",
    }
    for name, path in files.items():
        path.write_bytes(f"fake {name}".encode())
    cfg = assets / "diann.cfg"
    cfg.write_text(DIANN_CFG)
    db_dir = tmp_path / "db"
    (db_dir / "bakta" / "db-light").mkdir(parents=True)
    env: dict[str, Any] = {name: str(path) for name, path in files.items()}
    env["db_dir"] = str(db_dir)
    env["diann_cfg"] = str(cfg)
    env["pins"] = {name: _sha256(path) for name, path in files.items()}
    return db_dir, env


@dataclass
class Scenario:
    tmp_path: Path
    job: dict[str, Any]
    db_dir: Path
    env: dict[str, Any]
    genome: bytes
    raw: dict[str, bytes]
    thermo_half_width: float = DDA_HALF_WIDTH
    events: list[str] = field(default_factory=list)
    run: FakeRun = field(init=False)
    download: FakeDownload = field(init=False)

    def __post_init__(self) -> None:
        self.reset()

    def reset(self, **run_kwargs: Any) -> None:
        """Fresh fakes (a new process), the same workdir."""
        self.events = []
        run_kwargs.setdefault("thermo_half_width", self.thermo_half_width)
        self.run = FakeRun(self.events, **run_kwargs)
        self.download = FakeDownload(self.events, self.genome, self.raw)

    @property
    def workdir(self) -> Path:
        return Path(self.job["workdir"])

    def go(self, **kwargs: Any) -> int:
        args: dict[str, Any] = {
            "db_dir": self.db_dir,
            "profile": PROFILE,
            "threads": THREADS,
            "keep_inputs": False,
            "run": self.run,
            "download": self.download,
            "env": self.env,
        }
        args.update(kwargs)
        return driver.run_job(self.job, **args)

    def status(self) -> dict[str, Any]:
        return json.loads((self.workdir / "status.json").read_text())

    def wd(self, *parts: str) -> Path:
        return self.workdir.joinpath(*parts)


def _job(
    tmp_path: Path,
    raw_format: str,
    files: dict[str, bytes],
    denovo: Sequence[str],
    acquisition: str,
    checksums: dict[str, str] | None = None,
) -> tuple[dict[str, Any], dict[str, bytes]]:
    names = sorted(files)
    checksums = checksums or {}
    urls = {n: f"{PRIDE}/{n}" for n in names}
    job = {
        "job": 0,
        "dataset": DATASET,
        "tax_id": TAX_ID,
        "organism": ORGANISM,
        "genome_accession": ACC,
        "raw_format": raw_format,
        "acquisition_metadata": acquisition,
        "run_files": ";".join(names),
        "run_urls": ";".join(urls[n] for n in names),
        "run_bytes": ";".join(str(len(files[n])) for n in names),
        "run_checksums": ";".join(checksums.get(n, "") for n in names),
        "run_gb": round(sum(len(files[n]) for n in names) / 1e9, 2),
        "denovo_runs": ";".join(denovo),
        "workdir": str(tmp_path / "campaign" / f"{DATASET}_{TAX_ID}"),
    }
    return job, {urls[n]: files[n] for n in names}


def bruker_scenario(
    tmp_path: Path,
    types: dict[str, Sequence[int]] | None = None,
    denovo: Sequence[str] = ("A", "C"),
    acquisition: str = "DIA",
    replicons: Sequence[tuple[str, str]] = CIRCULAR,
) -> Scenario:
    types = types if types is not None else {"A": (9,), "B": (9,), "C": (9,)}
    scratch = tmp_path / "scratch"
    files = {f"{r}.d.zip": bruker_zip(scratch, r, t) for r, t in types.items()}
    checksums = {}
    if "A.d.zip" in files:
        checksums["A.d.zip"] = hashlib.md5(files["A.d.zip"]).hexdigest()
    if "C.d.zip" in files:
        checksums["C.d.zip"] = hashlib.sha256(files["C.d.zip"]).hexdigest()
    job, raw = _job(
        tmp_path, "bruker_d", files, [f"{r}.d.zip" for r in denovo], acquisition, checksums
    )
    db_dir, env = make_env(tmp_path)
    return Scenario(tmp_path, job, db_dir, env, genome_zip(replicons), raw)


def thermo_scenario(
    tmp_path: Path,
    acquisition: str = "DDA",
    half_width: float = DDA_HALF_WIDTH,
) -> Scenario:
    files = {f"X{i}.raw": f"thermo raw bytes X{i}".encode() for i in (1, 2, 3)}
    job, raw = _job(tmp_path, "thermo_raw", files, ["X1.raw", "X3.raw"], acquisition)
    db_dir, env = make_env(tmp_path)
    return Scenario(tmp_path, job, db_dir, env, genome_zip(), raw, thermo_half_width=half_width)


def _patch(monkeypatch: pytest.MonkeyPatch, module: str, attr: str, value: Any) -> None:
    """Patch ``giap_campaign.<module>.<attr>`` and any copy the driver imported by name."""
    monkeypatch.setattr(importlib.import_module(f"giap_campaign.{module}"), attr, value)
    monkeypatch.setattr(driver, attr, value, raising=False)


def _acq(mode: str) -> Any:
    def fake(*_: Any, **__: Any) -> str:
        return mode

    return fake


def _peptides(source: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": ["PEPTIDEK", "LESSK"],
            "modified_peptide": ["PEPTIDEK", "LESSK"],
            "protein_ids": ["p1", "rev_p2"],
            "q_value": [0.001, 0.5],
            "pep": [0.01, 0.9],
            "score": [3.2, 0.1],
            "is_decoy": [False, True],
            "run": ["X1", "X1"],
            "charge": [2, 2],
            "spectrum_id": ["10", "11"],
            "source": [source, source],
            "search_db": ["same", "same"],
        }
    )


class ConverterSpy:
    def __init__(self, source: str, issues: Sequence[Issue] = ()) -> None:
        self.source = source
        self.issues = list(issues)
        self.paths: list[list[Path]] = []

    def __call__(self, paths: Sequence[str | Path], *_: Any, **__: Any) -> ConvertResult:
        self.paths.append([Path(p) for p in paths])
        return ConvertResult(peptides=_peptides(self.source), abundance=None, issues=self.issues)


@pytest.fixture
def dda_patches(monkeypatch: pytest.MonkeyPatch) -> tuple[ConverterSpy, ConverterSpy]:
    sage_spy, casanovo_spy = ConverterSpy("database"), ConverterSpy("denovo")
    _patch(monkeypatch, "spectra", "mzml_acquisition", _acq("DDA"))
    _patch(monkeypatch, "sage", "convert_sage", sage_spy)
    _patch(monkeypatch, "casanovo", "convert_casanovo", casanovo_spy)
    return sage_spy, casanovo_spy


def _acq_mode(status: dict[str, Any]) -> str:
    """``status.json["acquisition"]`` is a plain string mode (lead decision)."""
    acq = status["acquisition"]
    assert isinstance(acq, str)
    assert isinstance(status["acquisition_note"], str)
    return acq


def _assert_all_done(status: dict[str, Any]) -> None:
    steps = status["steps"]
    for name in STEPS:
        assert steps[name]["state"] == "done", (name, steps.get(name))
    assert not [n for n, s in steps.items() if s["state"] != "done"]


# ================================================================ plan


SHEET_COLUMNS = [
    "dataset",
    "tax_id",
    "organism",
    "below_species",
    "title",
    "submission_type",
    "acquisition",
    "instruments",
    "pubmed_ids",
    "dois",
    "genome_accession",
    "genome_choice",
    "raw_format",
    "n_raw",
    "raw_gb",
    "raw_files",
    "raw_urls",
    "raw_bytes",
    "raw_checksums",
    "formats_seen",
    "flags",
    "status",
]


def sheet_row(
    dataset: str,
    tax_id: int = TAX_ID,
    *,
    sizes: Sequence[int | None] = (10**9, 10**9, 10**9),
    status: str = "ready",
    genome: str = ACC,
    raw_format: str = "thermo_raw",
    acquisition: str = "DDA",
    checksums: Sequence[str] | None = None,
) -> dict[str, Any]:
    names = [f"{dataset}_run{i:02d}.raw" for i in range(len(sizes))]
    checksums = list(checksums) if checksums is not None else [""] * len(sizes)
    total = sum(s or 0 for s in sizes)
    return {
        "dataset": dataset,
        "tax_id": tax_id,
        "organism": f"organism {tax_id}",
        "below_species": "true",
        "title": f"title {dataset}",
        "submission_type": "COMPLETE",
        "acquisition": acquisition,
        "instruments": "Q Exactive",
        "pubmed_ids": "",
        "dois": "",
        "genome_accession": genome,
        "genome_choice": "unique" if genome else "",
        "raw_format": raw_format,
        "n_raw": len(sizes) if raw_format else 0,
        "raw_gb": round(total / 1e9, 2),
        "raw_files": ";".join(names) if raw_format else "",
        "raw_urls": ";".join(f"{PRIDE}/{n}" for n in names) if raw_format else "",
        "raw_bytes": ";".join("" if s is None else str(s) for s in sizes) if raw_format else "",
        "raw_checksums": ";".join(checksums) if raw_format else "",
        "formats_seen": raw_format,
        "flags": "",
        "status": status,
    }


def sheet(rows: Sequence[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=SHEET_COLUMNS)


def even_indices(n: int, k: int) -> list[int]:
    """SPEC §7 select_runs spacing (only used on examples without .5 ties)."""
    k = min(k, n)
    if k == 1:
        return [0]
    return [round(i * (n - 1) / (k - 1)) for i in range(k)]


def test_plan_keeps_ready_rows_with_genome_and_raw_format_in_sheet_order(tmp_path: Path) -> None:
    rows = [
        sheet_row("PXD000009"),
        sheet_row("PXD000003", status="needs_genome", genome=""),
        sheet_row("PXD000004", genome=""),  # ready but no genome: curator cleared it
        sheet_row("PXD000005", raw_format=""),
        sheet_row("PXD000006", status="taxon_error"),
        sheet_row("PXD000002", tax_id=83333),
    ]
    jobs, _ = driver.plan(sheet(rows), tmp_path)
    assert list(jobs["dataset"]) == ["PXD000009", "PXD000002"]
    assert list(jobs["job"]) == [0, 1]
    assert [int(t) for t in jobs["tax_id"]] == [TAX_ID, 83333]


def test_plan_columns_exactly_as_spec(tmp_path: Path) -> None:
    jobs, _ = driver.plan(sheet([sheet_row("PXD000001")]), tmp_path)
    assert list(jobs.columns) == PLAN_COLUMNS


def test_plan_empty_sheet_gives_no_jobs_with_the_columns(tmp_path: Path) -> None:
    jobs, _ = driver.plan(sheet([sheet_row("PXD000001", status="needs_genome")]), tmp_path)
    assert len(jobs) == 0
    assert list(jobs.columns) == PLAN_COLUMNS


def test_plan_identity_columns_and_workdir(tmp_path: Path) -> None:
    row = sheet_row("PXD000007", tax_id=1280, acquisition="DIA", raw_format="thermo_raw")
    jobs, _ = driver.plan(sheet([row]), tmp_path / "root")
    job = jobs.iloc[0]
    assert job["dataset"] == "PXD000007"
    assert job["organism"] == "organism 1280"
    assert job["genome_accession"] == ACC
    assert job["raw_format"] == "thermo_raw"
    assert job["acquisition_metadata"] == "DIA"
    assert Path(job["workdir"]) == tmp_path / "root" / "PXD000007_1280"


@pytest.mark.parametrize(
    ("n", "size_gb", "max_runs", "max_gb", "expected_k"),
    [
        (20, 1, 12, 100.0, 12),  # max_runs binds: indices round(i * 19 / 11)
        (20, 10, 12, 100.0, 10),  # max_gb binds: 12 -> 120 GB, 11 -> 110, 10 -> 100 <= 100
        (5, 1, 12, 100.0, 5),  # fewer files than max_runs: all of them
        (4, 50, 12, 10.0, 1),  # every file too big: still one file, the first
    ],
)
def test_plan_applies_select_runs(
    tmp_path: Path, n: int, size_gb: int, max_runs: int, max_gb: float, expected_k: int
) -> None:
    sizes = [size_gb * 10**9] * n
    checksums = [f"{i:032x}" if i % 3 else "" for i in range(n)]  # empty items kept
    row = sheet_row("PXD000001", sizes=sizes, checksums=checksums)
    jobs, _ = driver.plan(sheet([row]), tmp_path, max_runs=max_runs, max_gb=max_gb)
    job = jobs.iloc[0]
    idx = even_indices(n, expected_k)
    files = row["raw_files"].split(";")
    assert job["run_files"].split(";") == [files[i] for i in idx]
    assert job["run_urls"].split(";") == [row["raw_urls"].split(";")[i] for i in idx]
    assert job["run_bytes"].split(";") == [str(sizes[i]) for i in idx]
    assert job["run_checksums"].split(";") == [checksums[i] for i in idx]
    assert float(job["run_gb"]) == pytest.approx(expected_k * size_gb, abs=0.01)


def test_plan_unknown_sizes_count_zero(tmp_path: Path) -> None:
    row = sheet_row("PXD000001", sizes=[None, None, 10**9])
    jobs, throttle = driver.plan(sheet([row]), tmp_path, max_gb=0.5)
    # sizes None count 0, so taking the first and last (0 + 1 GB) exceeds 0.5; then one file
    assert jobs.iloc[0]["run_files"].split(";") == ["PXD000001_run00.raw"]
    assert jobs.iloc[0]["run_bytes"] == ""
    assert float(jobs.iloc[0]["run_gb"]) == 0.0
    assert throttle == 8


def test_plan_denovo_runs_are_evenly_spaced_selected_runs(tmp_path: Path) -> None:
    row = sheet_row("PXD000001", sizes=[10**9] * 20)
    jobs, _ = driver.plan(sheet([row]), tmp_path, max_runs=12, denovo_runs=3)
    selected = jobs.iloc[0]["run_files"].split(";")
    assert len(selected) == 12
    # k = 3 of 12: indices round(i * 11 / 2) = 0, 6 (5.5 rounds to 6 either way), 11
    assert jobs.iloc[0]["denovo_runs"].split(";") == [selected[0], selected[6], selected[11]]


def test_plan_denovo_runs_one_is_the_first_selected(tmp_path: Path) -> None:
    row = sheet_row("PXD000001", sizes=[10**9] * 7)
    jobs, _ = driver.plan(sheet([row]), tmp_path, denovo_runs=1)
    assert jobs.iloc[0]["denovo_runs"] == jobs.iloc[0]["run_files"].split(";")[0]


def test_plan_denovo_runs_capped_by_selected(tmp_path: Path) -> None:
    row = sheet_row("PXD000001", sizes=[10**9] * 2)
    jobs, _ = driver.plan(sheet([row]), tmp_path, denovo_runs=3)
    assert jobs.iloc[0]["denovo_runs"] == jobs.iloc[0]["run_files"]


@pytest.mark.parametrize(
    ("job_gb", "disk", "expected"),
    [
        ((40, 20), 300.0, 2),  # floor(300 / (3 * 40)) = floor(2.5) = 2
        ((10,), 300.0, 8),  # floor(10) = 10, capped at 8
        ((40,), 10.0, 1),  # floor(0.083) = 0, at least 1
        ((50,), 300.0, 2),  # exactly 2
        ((0,), 300.0, 8),  # sizes unknown: run_gb 0, unbounded quotient -> 8
    ],
)
def test_plan_throttle(tmp_path: Path, job_gb: Sequence[int], disk: float, expected: int) -> None:
    rows = []
    for i, gb in enumerate(job_gb):
        sizes: list[int | None] = [10**9] * gb if gb else [None, None]
        rows.append(sheet_row(f"PXD00000{i + 1}", sizes=sizes))
    _, throttle = driver.plan(sheet(rows), tmp_path, max_runs=100, disk_budget_gb=disk)
    assert throttle == expected
    assert isinstance(throttle, int)


@settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    datasets=st.lists(
        st.lists(st.integers(min_value=0, max_value=30), min_size=1, max_size=25),
        min_size=1,
        max_size=4,
    ),
    max_runs=st.integers(min_value=1, max_value=15),
    max_gb=st.integers(min_value=1, max_value=200),
    denovo=st.integers(min_value=1, max_value=5),
    disk=st.integers(min_value=1, max_value=3000),
)
def test_plan_invariants(
    tmp_path: Path,
    datasets: list[list[int]],
    max_runs: int,
    max_gb: int,
    denovo: int,
    disk: int,
) -> None:
    rows = [
        sheet_row(f"PXD{i:06d}", sizes=[gb * 10**9 for gb in sizes])
        for i, sizes in enumerate(datasets)
    ]
    jobs, throttle = driver.plan(
        sheet(rows),
        tmp_path,
        max_runs=max_runs,
        max_gb=float(max_gb),
        denovo_runs=denovo,
        disk_budget_gb=float(disk),
    )
    assert list(jobs["job"]) == list(range(len(rows)))
    for row, (_, job) in zip(rows, jobs.iterrows(), strict=True):
        files = row["raw_files"].split(";")
        sizes = [int(b) for b in row["raw_bytes"].split(";")]
        selected = job["run_files"].split(";")
        pos = [files.index(f) for f in selected]
        assert pos == sorted(set(pos)), "selected runs keep the name order, no repeats"
        assert 1 <= len(selected) <= min(max_runs, len(files))
        total = sum(sizes[i] for i in pos) / 1e9
        assert len(selected) == 1 or total <= max_gb
        assert float(job["run_gb"]) == pytest.approx(total, abs=0.01)
        dn = job["denovo_runs"].split(";")
        dpos = [selected.index(f) for f in dn]
        assert dpos == sorted(set(dpos))
        assert len(dn) == min(denovo, len(selected))
        assert dn[0] == selected[0]
        if len(dn) >= 2:
            assert dn[-1] == selected[-1]
    max_run_gb = max(float(g) for g in jobs["run_gb"])
    expected = 8 if max_run_gb == 0 else max(1, min(8, math.floor(disk / (3 * max_run_gb))))
    assert throttle == expected


# ================================================================ run_job: DIA, Bruker


DIA_BRUKER_EVENTS = [
    "download_genome",
    "bakta",
    "giap_makedb",
    "download_raw",
    "download_raw",
    "download_raw",
    "timsconvert",
    "timsconvert",
    "diann1",
    "diann2",
    "cascadia",
    "cascadia",
    "giap_convert",
    "giap_convert",
    "giap_run",
    "giap_run",
    "giap_compare",
]


def test_dia_bruker_step_order_and_status(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    assert s.events == DIA_BRUKER_EVENTS
    status = s.status()
    _assert_all_done(status)
    assert _acq_mode(status) == "DIA"
    assert status["acquisition_note"] == ""


def test_every_argv_item_is_a_string(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    assert all(c.all_str for c in s.run.calls), [c.argv for c in s.run.calls if not c.all_str]
    assert "unknown" not in s.run.kinds()


def test_genome_download_url_and_raw_download_arguments(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    genome_call, *raw_calls = s.download.calls
    url = genome_call["url"]
    assert url.startswith(f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{ACC}/download")
    assert "include_annotation_type=GENOME_FASTA" in url
    assert "include_annotation_type=GENOME_GBFF" in url
    assert [c["url"] for c in raw_calls] == s.job["run_urls"].split(";")
    sizes = [int(b) for b in s.job["run_bytes"].split(";")]
    assert [c["expected_bytes"] for c in raw_calls] == sizes
    checks = s.job["run_checksums"].split(";")
    assert raw_calls[0]["checksum"] == checks[0]
    assert raw_calls[1]["checksum"] in (None, "")
    assert raw_calls[2]["checksum"] == checks[2]
    for c in raw_calls:
        assert s.wd("raw") in c["dest"].resolve().parents or c["dest"].parent == s.wd("raw")


def test_genome_headers_tagged_circular_and_bakta_complete(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go(keep_inputs=True) == 0
    fna = s.wd("genome", "genome.fna").read_text()
    headers = [line for line in fna.splitlines() if line.startswith(">")]
    assert headers == [">NC_000913.3 test replicon [topology=circular]"]
    (bakta,) = s.run.of("bakta")
    assert bakta.argv[:5] == ["micromamba", "run", "-n", "giap-bakta", "bakta"]
    assert "--complete" in bakta.argv
    assert "--keep-contig-headers" in bakta.argv
    assert bakta.opt("--prefix") == "genome"
    assert bakta.opt("--threads") == str(THREADS)
    assert _same(bakta.opt("--db"), s.db_dir / "bakta" / "db-light")  # bakta_db_path("auto")
    assert _same(bakta.path("--output"), s.wd("annotation"))
    assert _same(_at(bakta.cwd, bakta.argv[-1]), s.wd("genome", "genome.fna"))
    assert s.wd("annotation", "genome.json").is_file()


def test_linear_replicon_means_no_complete(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path, replicons=WITH_LINEAR)
    assert s.go() == 0
    headers = [
        line
        for line in s.wd("genome", "genome.fna").read_text().splitlines()
        if line.startswith(">")
    ]
    assert headers == [
        ">NC_000913.3 test replicon [topology=circular]",
        ">NZ_CP009001.1 test replicon",
    ]
    (bakta,) = s.run.of("bakta")
    assert "--complete" not in bakta.argv


def test_makedb_argv(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    (mk,) = s.run.of("giap_makedb")
    assert "--six-frame-orfs" in mk.argv
    assert mk.opt("--min-orf-aa") == "20"
    assert _same(mk.path("--annotation"), s.wd("annotation", "genome.json"))
    assert _same(mk.path("--genome"), s.wd("genome", "genome.fna"))
    assert _same(mk.opt("--db-dir"), s.db_dir)
    assert _same(mk.path("--outdir"), s.wd("searchdb"))
    assert s.wd("searchdb", "strain_db.faa").is_file()


def test_dia_bruker_converts_only_the_denovo_runs(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go(keep_inputs=True) == 0
    calls = s.run.of("timsconvert")
    assert [c.path("--input").name for c in calls] == ["A.d", "C.d"]
    assert all(c.info["input_has_tdf"] for c in calls)
    for c in calls:
        assert c.argv[:5] == ["micromamba", "run", "-n", "giap-timsconvert", "timsconvert"]
        assert _same(c.path("--outdir"), s.wd("mzml"))
    assert sorted(p.name for p in s.wd("mzml").iterdir()) == ["A.mzML", "C.mzML"]


def test_diann_two_steps_on_the_d_directories(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    (step1,) = s.run.of("diann1")
    (step2,) = s.run.of("diann2")
    for c in (step1, step2):
        i = c.argv.index("exec")
        assert c.argv[i - 1] == "apptainer"
        assert _same(c.argv[i + 1], s.env["diann_sif"])
        assert c.argv[i + 2] == "/diann-2.3.1/diann-linux"
        assert c.opt("--threads") == str(THREADS)
    # step 1: the [step1] flags of env["diann_cfg"] and the strain database
    assert "--predictor" in step1.argv
    assert "--fasta-search" in step1.argv
    assert step1.opt("--cut") == "K*,R*"
    assert [_at(step1.cwd, f).resolve() for f in _opts(step1.argv, "--fasta")] == [
        s.wd("searchdb", "strain_db.faa").resolve()
    ]
    assert "--f" not in step1.argv
    # step 2: the [step2] flags, every selected run as a .d directory, the report, decoys
    assert step2.opt("--qvalue") == "0.01"
    assert "--matrices" in step2.argv
    assert "--predictor" not in step2.argv
    assert "--report-decoys" in step2.argv
    assert "--lib" in step2.argv
    assert [p.name for p in step2.info["f"]] == ["A.d", "B.d", "C.d"]
    assert all(step2.info["f_ok"])
    assert _same(step2.path("--out"), s.wd("search", "report.parquet"))
    assert s.wd("search", "report.parquet").is_file()


def test_cascadia_per_denovo_run(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    calls = s.run.of("cascadia")
    assert [c.info["mzml"].name for c in calls] == ["A.mzML", "C.mzML"]
    assert all(c.info["mzml_exists"] for c in calls)
    assert all(_same(c.info["model"], s.env["cascadia_model"]) for c in calls)
    assert "giap-cascadia" in calls[0].argv
    assert [c.opt("-b") for c in calls] == ["32", "32"]
    assert sorted(p.name for p in s.wd("denovo").glob("*.ssl")) == ["A.ssl", "C.ssl"]


def test_dia_convert_uses_giap_convert(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    db_conv, dn_conv = s.run.of("giap_convert")
    assert db_conv.opt("--from") == "diann"
    assert db_conv.opt("--search-db") == "same"
    assert any(_same(_at(db_conv.cwd, t), s.wd("search", "report.parquet")) for t in db_conv.argv)
    assert _same(db_conv.path("--out"), s.wd("peptides_db.tsv"))
    assert dn_conv.opt("--from") == "cascadia"
    ssl = sorted(Path(t).name for t in dn_conv.argv if t.endswith(".ssl"))
    assert ssl == ["A.ssl", "C.ssl"]
    assert _same(dn_conv.path("--out"), s.wd("peptides_denovo.tsv"))
    assert s.wd("peptides_db.tsv").is_file()
    assert s.wd("peptides_denovo.tsv").is_file()


def test_giap_run_arms_and_compare(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    db_run, dn_run = s.run.of("giap_run")
    for c, outdir in ((db_run, "giap_db"), (dn_run, "giap_denovo")):
        assert _same(c.path("--genome"), s.wd("genome", "genome.fna"))
        assert _same(c.path("--annotation"), s.wd("annotation", "genome.json"))
        assert _same(c.path("--outdir"), s.wd(outdir))
        assert _same(c.opt("--db-dir"), s.db_dir)
        assert c.opt("--profile") == PROFILE
        assert c.opt("--threads") == str(THREADS)
        assert "--resume" in c.argv
    assert _same(db_run.path("--peptides"), s.wd("peptides_db.tsv"))
    assert "--reference-peptides" not in db_run.argv
    assert _same(dn_run.path("--peptides"), s.wd("peptides_denovo.tsv"))
    assert _same(dn_run.path("--reference-peptides"), s.wd("peptides_db.tsv"))
    (cmp_call,) = s.run.of("giap_compare")
    i = cmp_call.argv.index("compare")
    assert _same(_at(cmp_call.cwd, cmp_call.argv[i + 1]), s.wd("giap_db"))
    assert _same(_at(cmp_call.cwd, cmp_call.argv[i + 2]), s.wd("giap_denovo"))
    assert cmp_call.opt("--labels") == "db,denovo"
    assert _same(cmp_call.path("--outdir"), s.wd("compare"))
    assert s.wd("compare").is_dir()


# ================================================================ run_job: DDA and other inputs


DDA_THERMO_EVENTS = [
    "download_genome",
    "bakta",
    "giap_makedb",
    "download_raw",
    "download_raw",
    "download_raw",
    "thermo",
    "thermo",
    "thermo",
    "sage",
    "casanovo",
    "giap_convert",
    "giap_convert",
    "giap_run",
    "giap_run",
    "giap_compare",
]


def test_dda_thermo_path(tmp_path: Path, dda_patches: tuple[ConverterSpy, ConverterSpy]) -> None:
    sage_spy, casanovo_spy = dda_patches
    s = thermo_scenario(tmp_path)
    assert s.go(keep_inputs=True) == 0
    assert s.events == DDA_THERMO_EVENTS
    status = s.status()
    _assert_all_done(status)
    assert _acq_mode(status) == "DDA"
    # ThermoRawFileParser on every run
    thermo = s.run.of("thermo")
    assert [c.path("-i").name for c in thermo] == ["X1.raw", "X2.raw", "X3.raw"]
    assert all(c.info["input_exists"] for c in thermo)
    for c in thermo:
        i = c.argv.index("-n")
        assert c.argv[i + 1] == "giap-thermorawfileparser"
        assert _same(c.path("-o"), s.wd("mzml"))
        assert c.opt("-f") == "2"
    # Sage on every run's mzML
    (sage,) = s.run.of("sage")
    assert "giap-sage" in sage.argv
    cfg = sage.info["config"]
    assert _same(_at(sage.cwd, cfg["database"]["fasta"]), s.wd("searchdb", "strain_db.faa"))
    assert cfg["database"]["generate_decoys"] is True
    assert cfg["database"]["decoy_tag"] == "rev_"
    assert sorted(Path(p).name for p in cfg["mzml_paths"]) == ["X1.mzML", "X2.mzML", "X3.mzML"]
    assert sage.info["mzml_exist"]
    assert _same(_at(sage.cwd, cfg["output_directory"]), s.wd("search"))
    assert s.wd("search", "results.sage.tsv").is_file()
    # Casanovo on the denovo runs
    (cas,) = s.run.of("casanovo")
    assert "giap-casanovo" in cas.argv
    assert _same(cas.opt("--model"), s.env["casanovo_model"])
    assert _same(cas.path("--output"), s.wd("denovo", "denovo.mztab"))
    assert [p.name for p in cas.info["inputs"]] == ["X1.mzML", "X3.mzML"]
    assert cas.info["inputs_exist"]
    # converters, then giap convert --from generic
    assert [[p.name for p in ps] for ps in sage_spy.paths] == [["results.sage.tsv"]]
    assert [[p.name for p in ps] for ps in casanovo_spy.paths] == [["denovo.mztab"]]
    raw_db = pd.read_csv(s.wd("peptides_db.raw.tsv"), sep="\t", keep_default_na=False)
    assert "peptide" in raw_db.columns
    assert "PEPTIDEK" in set(raw_db["peptide"])
    db_conv, dn_conv = s.run.of("giap_convert")
    assert db_conv.opt("--from") == "generic"
    assert db_conv.opt("--search-db") == "same"
    assert any(_same(_at(db_conv.cwd, t), s.wd("peptides_db.raw.tsv")) for t in db_conv.argv)
    assert _same(db_conv.path("--out"), s.wd("peptides_db.tsv"))
    assert dn_conv.opt("--from") == "generic"
    assert _same(dn_conv.path("--out"), s.wd("peptides_denovo.tsv"))
    dn_inputs = [t for t in dn_conv.argv if t.endswith(".tsv") and t != dn_conv.opt("--out")]
    assert len(dn_inputs) == 1
    assert _at(dn_conv.cwd, dn_inputs[0]).is_file()
    # the arms
    _, dn_run = s.run.of("giap_run")
    assert _same(dn_run.path("--reference-peptides"), s.wd("peptides_db.tsv"))
    (cmp_call,) = s.run.of("giap_compare")
    assert cmp_call.opt("--labels") == "db,denovo"


def test_dda_thermo_acquisition_from_the_real_mzml(tmp_path: Path) -> None:
    """No patch: the mzML ThermoRawFileParser writes has 1.4 Th windows -> DDA (SPEC §7)."""
    s = thermo_scenario(tmp_path, acquisition="unknown")
    s.go(keep_inputs=True)
    assert _acq_mode(s.status()) == "DDA"
    assert s.run.of("sage")
    assert not s.run.of("diann1")


def test_dia_thermo_searches_mzml_with_diann(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _patch(monkeypatch, "spectra", "mzml_acquisition", _acq("DIA"))
    s = thermo_scenario(tmp_path, acquisition="DIA", half_width=DIA_HALF_WIDTH)
    assert s.go() == 0
    assert len(s.run.of("thermo")) == 3
    (step2,) = s.run.of("diann2")
    assert [p.name for p in step2.info["f"]] == ["X1.mzML", "X2.mzML", "X3.mzML"]
    assert all(step2.info["f_ok"])
    assert [c.info["mzml"].name for c in s.run.of("cascadia")] == ["X1.mzML", "X3.mzML"]
    assert not s.run.of("sage")
    assert not s.run.of("casanovo")
    assert _acq_mode(s.status()) == "DIA"


def test_dda_bruker_converts_every_run(
    tmp_path: Path, dda_patches: tuple[ConverterSpy, ConverterSpy]
) -> None:
    s = bruker_scenario(tmp_path, types={"A": (8,), "B": (8,), "C": (8,)}, acquisition="DDA")
    assert s.go() == 0
    assert [c.path("--input").name for c in s.run.of("timsconvert")] == ["A.d", "B.d", "C.d"]
    (sage,) = s.run.of("sage")
    assert sorted(Path(p).name for p in sage.info["config"]["mzml_paths"]) == [
        "A.mzML",
        "B.mzML",
        "C.mzML",
    ]
    (cas,) = s.run.of("casanovo")
    assert [p.name for p in cas.info["inputs"]] == ["A.mzML", "C.mzML"]
    assert not s.run.of("diann1")
    assert _acq_mode(s.status()) == "DDA"


def test_mzml_inputs_are_used_as_they_are(
    tmp_path: Path, dda_patches: tuple[ConverterSpy, ConverterSpy]
) -> None:
    files = {f"S{i}.mzML": mzml_text(DDA_HALF_WIDTH).encode() for i in (1, 2)}
    job, raw = _job(tmp_path, "mzml", files, ["S1.mzML"], "DDA")
    db_dir, env = make_env(tmp_path)
    s = Scenario(tmp_path, job, db_dir, env, genome_zip(), raw)
    assert s.go(keep_inputs=True) == 0
    assert not s.run.of("thermo")
    assert not s.run.of("timsconvert")
    (sage,) = s.run.of("sage")
    assert sorted(Path(p).name for p in sage.info["config"]["mzml_paths"]) == ["S1.mzML", "S2.mzML"]
    assert sage.info["mzml_exist"]
    (cas,) = s.run.of("casanovo")
    assert [p.name for p in cas.info["inputs"]] == ["S1.mzML"]


# ================================================================ acquisition


def test_mixed_acquisition_fails_the_job(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path, types={"A": (9,), "B": (8,), "C": (9,)})
    assert s.go() == 1
    status = s.status()
    assert status["steps"]["acquisition"]["state"] == "failed"
    assert status["steps"]["acquisition"]["message"]
    assert _acq_mode(status) == "mixed"
    assert status["acquisition_note"] == "mixed_acquisition"
    assert not set(s.run.kinds()) & (SEARCH_KINDS | {"timsconvert", "giap_run"})


def test_unknown_acquisition_fails_the_job(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path, types={"A": (), "B": (), "C": ()}, acquisition="unknown")
    assert s.go() == 1
    status = s.status()
    assert status["steps"]["acquisition"]["state"] == "failed"
    assert _acq_mode(status) == "unknown"
    assert status["acquisition_note"] == "undetermined"
    assert not set(s.run.kinds()) & SEARCH_KINDS


def test_undetermined_data_falls_back_to_metadata(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path, types={"A": (), "B": (), "C": ()}, acquisition="DIA")
    assert s.go() == 0
    assert _acq_mode(s.status()) == "DIA"
    assert s.status()["acquisition_note"] == "from_metadata"
    assert s.run.of("diann2")


def test_data_override_metadata(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path, acquisition="DDA")  # the frames say DIA
    assert s.go() == 0
    assert _acq_mode(s.status()) == "DIA"
    assert s.status()["acquisition_note"] == "metadata_disagrees"
    assert s.run.of("diann2")
    assert not s.run.of("sage")


# ================================================================ failures and resuming


def test_second_call_after_success_runs_nothing(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    s.reset()
    assert s.go() == 0
    assert s.run.calls == []
    assert s.download.calls == []


def test_first_failure_stops_and_resume_starts_there(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    s.reset(fail={"diann2"})
    assert s.go() == 1
    status = s.status()
    failed = status["steps"]["search_db"]
    assert failed["state"] == "failed"
    assert isinstance(failed["message"], str)
    assert failed["message"]
    for name in ("fetch_genome", "annotate", "makedb", "fetch_raw", "acquisition", "to_mzml"):
        assert status["steps"][name]["state"] == "done", name
    for name in ("denovo", "convert", "giap_db", "giap_denovo", "compare"):
        assert status["steps"].get(name, {}).get("state") != "done", name
    assert s.run.kinds()[-1] == "diann2"  # nothing after the failure

    s.reset()
    assert s.go() == 0
    kinds = s.run.kinds()
    assert kinds[0] in {"diann1", "diann2"}
    assert not {"bakta", "giap_makedb", "timsconvert"} & set(kinds)
    assert s.download.calls == []
    assert "diann2" in kinds
    assert "giap_compare" in kinds
    _assert_all_done(s.status())


def test_tool_exit_zero_without_output_fails_the_step(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    s.reset(no_output={"bakta"})
    assert s.go() == 1
    assert s.status()["steps"]["annotate"]["state"] == "failed"
    assert "giap_makedb" not in s.run.kinds()


def test_giap_convert_failure_fails_convert(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    s.reset(fail={"giap_convert"})
    assert s.go() == 1
    assert s.status()["steps"]["convert"]["state"] == "failed"
    assert "giap_run" not in s.run.kinds()


def test_converter_error_issue_fails_convert(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dda_patches: tuple[ConverterSpy, ConverterSpy]
) -> None:
    bad = ConverterSpy("database", [Issue("missing_column", "error", "no column 'peptide'")])
    _patch(monkeypatch, "sage", "convert_sage", bad)
    s = thermo_scenario(tmp_path)
    assert s.go() == 1
    assert s.status()["steps"]["convert"]["state"] == "failed"
    assert "giap_run" not in s.run.kinds()


# ================================================================ pins (environment errors)


def test_mismatched_diann_pin_is_an_environment_error(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    s.env["pins"]["diann_sif"] = "0" * 64
    assert s.go() == 2
    assert not set(s.run.kinds()) & SEARCH_KINDS


def test_missing_pin_entry_is_an_environment_error(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    del s.env["pins"]["cascadia_model"]
    assert s.go() == 2
    assert not set(s.run.kinds()) & SEARCH_KINDS


def test_missing_pinned_file_is_an_environment_error(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    Path(s.env["diann_sif"]).unlink()
    assert s.go() == 2
    assert not set(s.run.kinds()) & SEARCH_KINDS


def test_empty_casanovo_pin_fails_a_dda_job(
    tmp_path: Path, dda_patches: tuple[ConverterSpy, ConverterSpy]
) -> None:
    s = bruker_scenario(tmp_path, types={"A": (8,), "B": (8,), "C": (8,)}, acquisition="DDA")
    s.env["pins"]["casanovo_model"] = ""  # as committed in campaign/pins.tsv
    assert s.go() == 2
    assert not set(s.run.kinds()) & SEARCH_KINDS


def test_empty_casanovo_pin_does_not_affect_a_dia_job(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    s.env["pins"]["casanovo_model"] = ""
    assert s.go() == 0


def test_pin_fixed_then_the_job_resumes(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    good = s.env["pins"]["diann_sif"]
    s.env["pins"]["diann_sif"] = good.upper()[::-1]
    assert s.go() == 2
    s.env["pins"]["diann_sif"] = good
    s.reset()
    assert s.go() == 0
    assert "bakta" not in s.run.kinds()
    _assert_all_done(s.status())


# ================================================================ cleanup


def test_cleanup_removes_inputs_and_work(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go() == 0
    assert not s.wd("raw").exists()
    assert not s.wd("mzml").exists()
    assert not s.wd("giap_db", "work").exists()
    assert not s.wd("giap_denovo", "work").exists()
    search_files = sorted(p.name for p in s.wd("search").rglob("*") if p.is_file())
    assert search_files == ["report.parquet"]
    for kept in (
        ("status.json",),
        ("genome", "genome.fna"),
        ("searchdb", "strain_db.faa"),
        ("peptides_db.tsv",),
        ("peptides_denovo.tsv",),
        ("giap_db", "run_manifest.json"),
        ("giap_denovo", "run_manifest.json"),
        ("denovo", "A.ssl"),
    ):
        assert s.wd(*kept).is_file(), kept
    assert s.wd("compare").is_dir()
    assert s.status()["steps"]["cleanup"]["state"] == "done"


def test_cleanup_keeps_the_sage_results(
    tmp_path: Path, dda_patches: tuple[ConverterSpy, ConverterSpy]
) -> None:
    s = thermo_scenario(tmp_path)
    assert s.go() == 0
    assert not s.wd("raw").exists()
    assert not s.wd("mzml").exists()
    search_files = sorted(p.name for p in s.wd("search").rglob("*") if p.is_file())
    assert search_files == ["results.sage.tsv"]


def test_keep_inputs_skips_cleanup(tmp_path: Path) -> None:
    s = bruker_scenario(tmp_path)
    assert s.go(keep_inputs=True) == 0
    assert any(s.wd("raw").rglob("analysis.tdf"))
    assert sorted(p.name for p in s.wd("mzml").iterdir()) == ["A.mzML", "C.mzML"]
    assert s.wd("giap_db", "work").is_dir()
    assert s.wd("giap_denovo", "work").is_dir()
    assert s.status()["steps"].get("cleanup", {}).get("state") != "failed"
