"""Pure command and config builders (SPEC §8)."""

from __future__ import annotations

import inspect
import json
import tempfile
from pathlib import Path

import pytest
from giap.mapping.digest import cut_sites, digest
from giap_campaign.engines import (
    bakta_argv,
    bruker_to_mzml_argv,
    cascadia_argv,
    casanovo_argv,
    diann_argvs,
    giap_compare_argv,
    giap_run_argv,
    makedb_argv,
    sage_argv,
    sage_config,
    thermo_to_mzml_argv,
)
from hypothesis import given, settings
from hypothesis import strategies as st

W = Path("/scratch/giap/campaign/PXD000001_562")
DIANN = "/diann-2.3.1/diann-linux"


def _all_str(argv: list[str]) -> None:
    assert isinstance(argv, list)
    assert all(isinstance(a, str) for a in argv), argv


# ---------------------------------------------------------------- bakta, makedb, converters


def test_bakta_argv_complete() -> None:
    argv = bakta_argv(
        genome=W / "genome/genome.fna",
        outdir=W / "annotation",
        db=Path("/db/bakta/db-light"),
        threads=16,
        complete=True,
    )
    _all_str(argv)
    assert argv == [
        "micromamba",
        "run",
        "-n",
        "giap-bakta",
        "bakta",
        "--db",
        "/db/bakta/db-light",
        "--output",
        f"{W}/annotation",
        "--prefix",
        "genome",
        "--threads",
        "16",
        "--keep-contig-headers",
        "--complete",
        f"{W}/genome/genome.fna",
    ]


def test_bakta_argv_not_complete() -> None:
    argv = bakta_argv(
        genome=W / "genome/genome.fna",
        outdir=W / "annotation",
        db=Path("/db/bakta/db"),
        threads=8,
        complete=False,
    )
    _all_str(argv)
    assert argv == [
        "micromamba",
        "run",
        "-n",
        "giap-bakta",
        "bakta",
        "--db",
        "/db/bakta/db",
        "--output",
        f"{W}/annotation",
        "--prefix",
        "genome",
        "--threads",
        "8",
        "--keep-contig-headers",
        f"{W}/genome/genome.fna",
    ]
    assert "--complete" not in argv


def test_bakta_argv_accepts_str_paths() -> None:
    argv = bakta_argv("g.fna", "out", "/db", 4, False)
    assert argv[-1] == "g.fna"
    assert argv[argv.index("--db") + 1] == "/db"


def test_makedb_argv() -> None:
    argv = makedb_argv(
        annotation_json=W / "annotation/genome.json",
        genome=W / "genome/genome.fna",
        outdir=W / "searchdb",
        db_dir=Path("/db"),
    )
    _all_str(argv)
    assert argv == [
        "giap",
        "makedb",
        "--annotation",
        f"{W}/annotation/genome.json",
        "--genome",
        f"{W}/genome/genome.fna",
        "--six-frame-orfs",
        "--min-orf-aa",
        "20",
        "--db-dir",
        "/db",
        "--outdir",
        f"{W}/searchdb",
    ]


def test_thermo_to_mzml_argv() -> None:
    argv = thermo_to_mzml_argv(raw=W / "raw/Ecoli_WT_1.raw", outdir=W / "mzml")
    _all_str(argv)
    assert argv == [
        "micromamba",
        "run",
        "-n",
        "giap-thermorawfileparser",
        "ThermoRawFileParser",
        "-i",
        f"{W}/raw/Ecoli_WT_1.raw",
        "-o",
        f"{W}/mzml",
        "-f",
        "2",
    ]


def test_bruker_to_mzml_argv() -> None:
    argv = bruker_to_mzml_argv(d_dir=W / "raw/run1.d", outdir=W / "mzml")
    _all_str(argv)
    assert argv == [
        "micromamba",
        "run",
        "-n",
        "giap-timsconvert",
        "timsconvert",
        "--input",
        f"{W}/raw/run1.d",
        "--outdir",
        f"{W}/mzml",
    ]


# ---------------------------------------------------------------- de novo


def test_cascadia_argv() -> None:
    argv = cascadia_argv(
        mzml=W / "mzml/run1.mzML",
        out_prefix=W / "denovo/run1",
        model=Path("/db/cascadia/cascadia_0.0.7.ckpt"),
    )
    _all_str(argv)
    # in the giap-cascadia env, as scripts/slurm/denovo_pilot.sbatch (which may go through the
    # offline launcher): ... sequence <mzml> <model> -o <prefix> -b 32
    assert argv[:4] == ["micromamba", "run", "-n", "giap-cascadia"]
    assert argv[-7:] == [
        "sequence",
        f"{W}/mzml/run1.mzML",
        "/db/cascadia/cascadia_0.0.7.ckpt",
        "-o",
        f"{W}/denovo/run1",
        "-b",
        "32",
    ]
    assert any("cascadia" in a for a in argv[4:-7])


def test_casanovo_argv() -> None:
    mzmls = [W / "mzml/run1.mzML", W / "mzml/run3.mzML", W / "mzml/run2.mzML"]
    argv = casanovo_argv(
        mzmls=mzmls,
        out_mztab=W / "denovo/denovo.mztab",
        model=Path("/db/casanovo/casanovo_massivekb_5.2.1.ckpt"),
    )
    _all_str(argv)
    assert argv == [
        "micromamba",
        "run",
        "-n",
        "giap-casanovo",
        "casanovo",
        "sequence",
        "--model",
        "/db/casanovo/casanovo_massivekb_5.2.1.ckpt",
        "--output",
        f"{W}/denovo/denovo.mztab",
        f"{W}/mzml/run1.mzML",
        f"{W}/mzml/run3.mzML",
        f"{W}/mzml/run2.mzML",
    ]


# ---------------------------------------------------------------- giap run / compare


def test_giap_run_argv_db_arm() -> None:
    argv = giap_run_argv(
        genome=W / "genome/genome.fna",
        annotation=W / "annotation/genome.json",
        peptides=W / "peptides_db.tsv",
        outdir=W / "giap_db",
        db_dir=Path("/db"),
        profile="farmshare_local,apptainer",
        threads=16,
    )
    _all_str(argv)
    assert argv == [
        "giap",
        "run",
        "--genome",
        f"{W}/genome/genome.fna",
        "--annotation",
        f"{W}/annotation/genome.json",
        "--peptides",
        f"{W}/peptides_db.tsv",
        "--outdir",
        f"{W}/giap_db",
        "--db-dir",
        "/db",
        "--profile",
        "farmshare_local,apptainer",
        "--threads",
        "16",
        "--resume",
    ]
    assert "--reference-peptides" not in argv


def test_giap_run_argv_denovo_arm_with_reference() -> None:
    argv = giap_run_argv(
        W / "genome/genome.fna",
        W / "annotation/genome.json",
        W / "peptides_denovo.tsv",
        W / "giap_denovo",
        Path("/db"),
        "farmshare_local,apptainer",
        16,
        reference=W / "peptides_db.tsv",
    )
    _all_str(argv)
    assert argv == [
        "giap",
        "run",
        "--genome",
        f"{W}/genome/genome.fna",
        "--annotation",
        f"{W}/annotation/genome.json",
        "--peptides",
        f"{W}/peptides_denovo.tsv",
        "--reference-peptides",
        f"{W}/peptides_db.tsv",
        "--outdir",
        f"{W}/giap_denovo",
        "--db-dir",
        "/db",
        "--profile",
        "farmshare_local,apptainer",
        "--threads",
        "16",
        "--resume",
    ]


def test_giap_compare_argv() -> None:
    argv = giap_compare_argv(
        db_out=W / "giap_db", denovo_out=W / "giap_denovo", outdir=W / "compare"
    )
    _all_str(argv)
    assert argv == [
        "giap",
        "compare",
        f"{W}/giap_db",
        f"{W}/giap_denovo",
        "--outdir",
        f"{W}/compare",
        "--labels",
        "db,denovo",
    ]


# ---------------------------------------------------------------- DIA-NN

CFG = """\
# DIA-NN settings for the campaign (A2 flags without the strain --fasta line).
# One flag per line; $VAR / ${VAR} are expanded from the environment.

[step1]
--fasta-search
--gen-spec-lib
--predictor
--fasta $GIAP_DB_DIR/contaminants/contaminants.fasta
--cut K*,R*
--missed-cleavages 1
--min-pep-len 7
--max-pep-len 30
--unimod4
--var-mods 1
--var-mod UniMod:35,15.994915,M
--min-pr-charge 1
--max-pr-charge 4
--verbose 1   # quiet

[step2]
# the search itself
--fasta ${GIAP_DB_DIR}/contaminants/contaminants.fasta
--gen-spec-lib
--qvalue 0.01
--matrices
--peptidoforms
--temp ${GIAP_ROOT}/tmp
--verbose 1
"""

STEP1_CFG = [
    "--fasta-search",
    "--gen-spec-lib",
    "--predictor",
    "--fasta",
    "/db/contaminants/contaminants.fasta",
    "--cut",
    "K*,R*",
    "--missed-cleavages",
    "1",
    "--min-pep-len",
    "7",
    "--max-pep-len",
    "30",
    "--unimod4",
    "--var-mods",
    "1",
    "--var-mod",
    "UniMod:35,15.994915,M",
    "--min-pr-charge",
    "1",
    "--max-pr-charge",
    "4",
    "--verbose",
    "1",
]
STEP2_CFG = [
    "--fasta",
    "/db/contaminants/contaminants.fasta",
    "--gen-spec-lib",
    "--qvalue",
    "0.01",
    "--matrices",
    "--peptidoforms",
    "--temp",
    "/root/tmp",
    "--verbose",
    "1",
]


def _is_subsequence(needle: list[str], haystack: list[str]) -> bool:
    it = iter(haystack)
    return all(any(h == n for h in it) for n in needle)


def _value(argv: list[str], flag: str) -> str:
    assert argv.count(flag) == 1, (flag, argv)
    return argv[argv.index(flag) + 1]


def _fasta_values(argv: list[str]) -> list[str]:
    return [argv[i + 1] for i, a in enumerate(argv) if a == "--fasta"]


@pytest.fixture
def diann_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIAP_DB_DIR", "/db")
    monkeypatch.setenv("GIAP_ROOT", "/root")
    monkeypatch.setenv("GIAP_REPO", "/repo")


def _diann(tmp_path: Path, cfg: str = CFG, spectra: list[Path] | None = None) -> list[list[str]]:
    cfg_path = tmp_path / "diann.cfg"
    cfg_path.write_text(cfg, encoding="utf-8")
    if spectra is None:
        spectra = [W / "raw/run1.d", W / "raw/run2.d", W / "raw/run10.d"]
    out = diann_argvs(
        spectra=spectra,
        fasta=W / "searchdb/strain_db.faa",
        outdir=W / "search",
        threads=16,
        cfg_path=cfg_path,
        sif=Path("/images/diann_2.3.1.sif"),
    )
    return [list(a) for a in out]


def _check_prefix(argv: list[str]) -> list[str]:
    """``apptainer exec [-B <bind> ...] <sif> /diann-2.3.1/diann-linux``; returns the DIA-NN args."""
    assert argv[:2] == ["apptainer", "exec"]
    k = argv.index("/images/diann_2.3.1.sif")
    binds = argv[2:k]
    assert len(binds) % 2 == 0
    assert all(b in ("-B", "--bind") for b in binds[::2]), binds
    assert argv[k + 1] == DIANN
    return argv[k + 2 :]


@pytest.mark.usefixtures("diann_env")
def test_diann_argvs_two_steps(tmp_path: Path) -> None:
    out = _diann(tmp_path)
    assert len(out) == 2
    step1, step2 = out
    _all_str(step1)
    _all_str(step2)
    a1 = _check_prefix(step1)
    a2 = _check_prefix(step2)
    # the cfg flags, in file order, split on whitespace, comments dropped, variables expanded
    assert _is_subsequence(STEP1_CFG, a1), a1
    assert _is_subsequence(STEP2_CFG, a2), a2
    assert "#" not in " ".join(step1 + step2)
    assert "quiet" not in step1
    assert not any("$" in a for a in step1 + step2)
    assert not any(a.startswith("[step") for a in step1 + step2)
    # step-2 cfg flags do not leak into step 1 and vice versa
    assert "--qvalue" not in a1 and "--matrices" not in a1
    assert "--predictor" not in a2 and "--cut" not in a2


@pytest.mark.usefixtures("diann_env")
def test_diann_argvs_fasta_threads_and_outputs(tmp_path: Path) -> None:
    step1, step2 = _diann(tmp_path)
    strain = f"{W}/searchdb/strain_db.faa"
    for argv in (step1, step2):
        assert strain in _fasta_values(argv)
        assert "/db/contaminants/contaminants.fasta" in _fasta_values(argv)
        assert _value(argv, "--threads") == "16"
    # step 1 predicts a library (no runs); step 2 searches the runs against it
    assert "--f" not in step1
    lib_out = _value(step1, "--out-lib")
    assert lib_out.startswith(f"{W}/search/")
    lib_in = _value(step2, "--lib")
    assert lib_in.startswith(f"{W}/search/")
    assert _value(step2, "--out") == f"{W}/search/report.parquet"
    assert "--out" not in step1


@pytest.mark.usefixtures("diann_env")
def test_diann_argvs_one_f_per_spectra_file_in_order(tmp_path: Path) -> None:
    spectra = [W / "mzml/b.mzML", W / "mzml/a.mzML", W / "mzml/c.mzML"]
    _, step2 = _diann(tmp_path, spectra=spectra)
    fs = [step2[i + 1] for i, a in enumerate(step2) if a == "--f"]
    assert fs == [str(s) for s in spectra]


@pytest.mark.usefixtures("diann_env")
def test_diann_argvs_report_decoys_always_exactly_once(tmp_path: Path) -> None:
    _, step2 = _diann(tmp_path)  # CFG has no --report-decoys
    assert step2.count("--report-decoys") == 1
    with_decoys = CFG.replace("--qvalue 0.01\n", "--qvalue 0.01\n--report-decoys\n")
    step1, step2 = _diann(tmp_path, cfg=with_decoys)
    assert step2.count("--report-decoys") == 1
    assert step1.count("--report-decoys") <= 1


def test_diann_argvs_expands_both_variable_forms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIAP_DB_DIR", "/other/db")
    monkeypatch.setenv("GIAP_ROOT", "/other/root")
    step1, step2 = _diann(tmp_path)
    assert "/other/db/contaminants/contaminants.fasta" in _fasta_values(step1)  # $VAR
    assert "/other/db/contaminants/contaminants.fasta" in _fasta_values(step2)  # ${VAR}
    assert "/other/root/tmp" in step2


@pytest.mark.usefixtures("diann_env")
def test_diann_argvs_reads_the_repo_pilot_cfg(tmp_path: Path) -> None:
    """The pilot A2 cfg (sections, comments, two --fasta lines) parses the same way."""
    pilot = Path(__file__).resolve().parents[2] / "validation/pilot/diann_A2.cfg"
    if not pilot.is_file():
        pytest.skip("validation/pilot/diann_A2.cfg not in this checkout")
    step1, step2 = _diann(tmp_path, cfg=pilot.read_text(encoding="utf-8"))
    assert "--predictor" in step1 and "--f" not in step1
    assert step2.count("--report-decoys") == 1
    assert "--cut" in step1 and _value(step1, "--cut") == "K*,R*"
    assert "/db/contaminants/contaminants.fasta" in _fasta_values(step2)


_PATHS = st.lists(
    st.text(alphabet="abcdefghij0123456789_", min_size=1, max_size=6).map(
        lambda s: f"/data/{s}.mzML"
    ),
    min_size=1,
    max_size=8,
)


@settings(max_examples=30, deadline=None)
@given(spectra=_PATHS, threads=st.integers(min_value=1, max_value=64))
def test_diann_argvs_property(spectra: list[str], threads: int) -> None:
    literal = (
        CFG.replace("${GIAP_DB_DIR}", "/db")
        .replace("$GIAP_DB_DIR", "/db")
        .replace("${GIAP_ROOT}", "/root")
    )
    with tempfile.TemporaryDirectory() as d:
        cfg_path = Path(d) / "diann.cfg"
        cfg_path.write_text(literal, encoding="utf-8")
        step1, step2 = diann_argvs(
            [Path(s) for s in spectra], Path("/f.faa"), Path("/o"), threads, cfg_path, "/x.sif"
        )
    fs = [step2[i + 1] for i, a in enumerate(step2) if a == "--f"]
    assert fs == spectra
    assert step2.count("--report-decoys") == 1
    assert step1[step1.index("--threads") + 1] == str(threads)
    assert step2[step2.index("--threads") + 1] == str(threads)
    assert step2[step2.index("--out") + 1] == "/o/report.parquet"


# ---------------------------------------------------------------- Sage


def _sage() -> dict:
    return sage_config(
        fasta=W / "searchdb/strain_db.faa",
        mzmls=[W / "mzml/run2.mzML", W / "mzml/run1.mzML"],
        outdir=W / "search",
    )


def test_sage_config_values() -> None:
    cfg = _sage()
    db = cfg["database"]
    assert db["fasta"] == f"{W}/searchdb/strain_db.faa"
    assert db["enzyme"] == {
        "cleave_at": "KR",
        "restrict": None,
        "missed_cleavages": 1,
        "min_len": 7,
        "max_len": 30,
    }
    assert db["static_mods"] == {"C": 57.0215}
    assert db["variable_mods"] == {"M": [15.9949]}
    assert isinstance(db["variable_mods"]["M"], list)
    assert db["decoy_tag"] == "rev_"
    assert db["generate_decoys"] is True
    assert cfg["precursor_tol"] == {"ppm": [-20, 20]}
    assert cfg["fragment_tol"] == {"ppm": [-20, 20]}
    assert cfg["output_directory"] == f"{W}/search"
    assert cfg["mzml_paths"] == [f"{W}/mzml/run2.mzML", f"{W}/mzml/run1.mzML"]


def test_sage_config_enzyme_matches_giap_digest() -> None:
    """Sage digests the strain database as GIAP's own digest does (giap.mapping.digest)."""
    enzyme = _sage()["database"]["enzyme"]
    params = inspect.signature(digest).parameters
    assert enzyme["missed_cleavages"] == params["missed"].default
    assert enzyme["min_len"] == params["min_len"].default
    assert enzyme["max_len"] == params["max_len"].default
    default_enzyme = params["enzyme"].default
    # trypsin_p cuts after K/R even before P: Sage's `restrict: null`
    assert default_enzyme == "trypsin_p"
    assert enzyme["restrict"] is None
    protein = "MAKPAAGRPLLKEDR"
    giap_cuts = cut_sites(protein, default_enzyme)
    sage_cuts = [i + 1 for i, aa in enumerate(protein[:-1]) if aa in enzyme["cleave_at"]]
    assert giap_cuts == [*sage_cuts, len(protein)]


def test_sage_config_is_json_serialisable(tmp_path: Path) -> None:
    cfg = _sage()
    text = json.dumps(cfg)
    assert json.loads(text) == cfg
    loaded = json.loads(text)
    assert loaded["database"]["enzyme"]["restrict"] is None
    assert loaded["database"]["generate_decoys"] is True


def test_sage_argv() -> None:
    argv = sage_argv(W / "search/sage.json")
    _all_str(argv)
    assert argv == ["micromamba", "run", "-n", "giap-sage", "sage", f"{W}/search/sage.json"]
