"""T1 and static checks for the R fgsea engine (M8-17a; ADR-0038; spec §1, §3, A3).

- ``workflow/envs/fgsea.yml`` is clean under ``env_problems`` and pins r-base 4.5.3, fgsea 1.36.2
  (``FGSEA_VERSION``) and data.table 1.18.6.1.
- ``workflow/envs/locks/fgsea.explicit.txt`` holds the same three builds as #192's
  ``golden_r.explicit.txt``, only conda.anaconda.org URLs, and no Bioconductor data package,
  goseq or BiasedUrn (ADR-0038 item 5).
- Lead decision 9: conda-only, no apptainer anywhere in the wrapper.
- R5: ``scripts/ci/check_no_rscript.py src`` passes with the new wrapper in the tree;
  ``workflow/bin/giap_fgsea.R`` is executable with the spec's shebang and installs nothing.
"""

from __future__ import annotations

import ast
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
ENV = REPO / "workflow" / "envs" / "fgsea.yml"
LOCK = REPO / "workflow" / "envs" / "locks" / "fgsea.explicit.txt"
WRAPPER = REPO / "src" / "giap" / "io" / "fgsea_r.py"
SCRIPT = REPO / "workflow" / "bin" / "giap_fgsea.R"
CHECK = REPO / "scripts" / "ci" / "check_no_rscript.py"
RECORD = REPO / "scripts" / "slurm" / "record_fgsea.sbatch"

# The three builds of #192's scripts/golden/golden_r.explicit.txt (spec §1).
GOLDEN_BUILDS = (
    "https://conda.anaconda.org/conda-forge/linux-64/r-base-4.5.3-h502d0c9_4.conda"
    "#efdd36d3bd3c0fa459cfd37c9e4200a8",
    "https://conda.anaconda.org/conda-forge/linux-64/r-data.table-1.18.6.1-r45hbd9402d_0.conda"
    "#34e991670e902d15eaeb84790048974b",
    "https://conda.anaconda.org/bioconda/linux-64/bioconductor-fgsea-1.36.2-r45ha27e39d_0.conda"
    "#28894f5017c6780c72a8165b528fb3e9",
)
FORBIDDEN_EXACT = {
    "bioconductor-goseq",
    "bioconductor-genelendatabase",
    "bioconductor-go.db",
    "bioconductor-data-packages",
    "r-biasedurn",
}
FORBIDDEN_PREFIX = ("bioconductor-org.", "bioconductor-txdb.", "bioconductor-bsgenome")


def _lock_urls() -> list[str]:
    return [
        line.strip()
        for line in LOCK.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    ]


def _package(url: str) -> str:
    """``.../r-data.table-1.18.6.1-r45hbd9402d_0.conda#md5`` → ``r-data.table``."""
    base = url.split("#")[0].rsplit("/", 1)[1]
    base = re.sub(r"\.(conda|tar\.bz2)$", "", base)
    return base.rsplit("-", 2)[0].lower()


# ---------------------------------------------------------------- env file


def test_env_file_passes_env_problems() -> None:
    from giap.annotate.tools import env_problems, read_env_file

    env = read_env_file(ENV)
    assert env.name == "giap-fgsea"
    assert env_problems(env) == []
    assert env.pip == ()


def test_env_pins_match_spec_and_wrapper() -> None:
    from giap.annotate.tools import read_env_file
    from giap.io.fgsea_r import FGSEA_VERSION

    assert FGSEA_VERSION == "1.36.2"
    conda = set(read_env_file(ENV).conda)
    assert conda == {
        "r-base=4.5.3",
        f"bioconductor-fgsea={FGSEA_VERSION}",
        "r-data.table=1.18.6.1",
    }


def test_env_file_has_no_goseq_or_data_package() -> None:
    text = ENV.read_text(encoding="utf-8").lower()
    deps = text.split("dependencies:", 1)[1]
    assert "goseq" not in deps
    assert "biasedurn" not in deps
    assert "pip" not in deps


# ---------------------------------------------------------------- explicit lock


def test_lock_exists_and_is_explicit() -> None:
    text = LOCK.read_text(encoding="utf-8")
    assert "@EXPLICIT" in text.splitlines()
    urls = _lock_urls()
    assert urls
    for url in urls:
        assert url.startswith(
            ("https://conda.anaconda.org/conda-forge/", "https://conda.anaconda.org/bioconda/")
        ), url


@pytest.mark.parametrize("build", GOLDEN_BUILDS)
def test_lock_has_golden_build(build: str) -> None:
    assert build in _lock_urls()


def test_lock_has_one_build_per_package() -> None:
    names = [_package(u) for u in _lock_urls()]
    assert len(names) == len(set(names))
    assert {"r-base", "r-data.table", "bioconductor-fgsea"} <= set(names)


def test_lock_has_no_forbidden_package() -> None:
    names = {_package(u) for u in _lock_urls()}
    bad = sorted(n for n in names if n in FORBIDDEN_EXACT or n.startswith(FORBIDDEN_PREFIX))
    assert bad == []
    assert not any("goseq" in n for n in names)


def test_lock_forbidden_check_catches_data_packages() -> None:
    """The name parser used above really sees the packages it must reject."""
    for url, name in [
        (
            "https://conda.anaconda.org/bioconda/noarch/bioconductor-go.db-3.21.0-r45hdfd78af_0"
            ".tar.bz2#0",
            "bioconductor-go.db",
        ),
        (
            "https://conda.anaconda.org/bioconda/linux-64/bioconductor-goseq-1.62.0-r45hdfd78af_0"
            ".conda#0",
            "bioconductor-goseq",
        ),
        (
            "https://conda.anaconda.org/conda-forge/linux-64/r-biasedurn-2.0.12-r45h_0.conda#0",
            "r-biasedurn",
        ),
    ]:
        assert _package(url) == name


# ---------------------------------------------------------------- conda-only (lead decision 9)


def _string_constants(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    docstrings = {
        id(node.body[0].value)
        for node in ast.walk(tree)
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef))
        and node.body
        and isinstance(node.body[0], ast.Expr)
        and isinstance(node.body[0].value, ast.Constant)
    }
    return [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in docstrings
    ]


def test_wrapper_never_reads_an_image_variable() -> None:
    """Lead decision 9 was reversed by M11-14 14c (the image launcher starts ``apptainer exec``
    or ``docker run`` itself); the image ref still comes from images.json, never from a
    ``GIAP_FGSEA_IMAGE`` override."""
    consts = _string_constants(WRAPPER)
    assert not [c for c in consts if c == "GIAP_FGSEA_IMAGE"]


def test_wrapper_reuses_env_built() -> None:
    """Launcher rule 3 reuses ``preflight.environment._env_built`` (spec §4), it is not copied."""
    text = WRAPPER.read_text(encoding="utf-8")
    assert "_env_built" in text
    assert "def _env_built" not in text


def test_wrapper_does_not_import_rpy2() -> None:
    tree = ast.parse(WRAPPER.read_text(encoding="utf-8"))
    mods = [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names] + [
        n.module or "" for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
    ]
    assert not [m for m in mods if m.split(".")[0] == "rpy2"]


# ---------------------------------------------------------------- R5: check_no_rscript


def test_check_no_rscript_passes_on_src_with_wrapper() -> None:
    assert WRAPPER.is_file()
    proc = subprocess.run(
        [sys.executable, str(CHECK), str(REPO / "src")], capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_check_no_rscript_still_catches_a_mention(tmp_path: Path) -> None:
    (tmp_path / "x.py").write_text('CMD = ["R' + 'script", "a.R"]\n')
    proc = subprocess.run(
        [sys.executable, str(CHECK), str(tmp_path)], capture_output=True, text=True
    )
    assert proc.returncode == 1


# ---------------------------------------------------------------- giap_fgsea.R (static)


def test_script_shebang_and_mode() -> None:
    lines = SCRIPT.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "#!/usr/bin/env -S Rscript --vanilla"
    assert os.access(SCRIPT, os.X_OK)


def test_script_calls_fgsea_multilevel_as_specified() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "fgseaMultilevel(" in text
    assert re.search(r'scoreType\s*=\s*"std"', text)
    assert re.search(r"nproc\s*=\s*1\b", text)
    assert "set.seed(" in text
    assert "withCallingHandlers" in text and "muffleWarning" in text
    assert re.search(r"useBytes\s*=\s*TRUE", text)
    assert "file.rename" in text  # results.tsv.tmp → results.tsv, never a partial file
    for code in (2, 3):
        assert re.search(rf"(quit|q)\(\s*(status\s*=\s*)?{code}\b", text), code


def test_script_installs_and_downloads_nothing() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for bad in (
        "install.packages",
        "BiocManager",
        "goseq",
        "BiasedUrn",
        "http://",
        "https://",
        "download.file",
        "apptainer",
        "system(",
        "system2(",
    ):
        assert bad not in text, bad


def test_script_accepts_exactly_the_spec_arguments() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for arg in (
        "--stats",
        "--sets",
        "--jobs",
        "--sample-size",
        "--eps",
        "--n-perm-simple",
        "--seed",
        "--out",
        "--versions",
    ):
        assert arg in text, arg


# ---------------------------------------------------------------- record_fgsea.sbatch (R14)


def test_record_sbatch_header() -> None:
    text = RECORD.read_text(encoding="utf-8")
    assert text.splitlines()[0] == "#!/bin/bash"
    assert re.search(r"^#SBATCH --job-name=giap-record-fgsea$", text, re.M)
    cpus = re.search(r"^#SBATCH --cpus-per-task=(\d+)$", text, re.M)
    assert cpus and int(cpus.group(1)) <= 2
    assert re.search(r"^#SBATCH --time=00:30:00$", text, re.M)
    assert "scripts/farmshare_env.sh" in text
    for name in ("sim_main", "sim_retest", "ties", "unbalanced", "tiny_p"):
        assert name in text, name
    proc = subprocess.run(["bash", "-n", str(RECORD)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
