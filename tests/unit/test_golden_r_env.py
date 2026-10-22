"""The R golden environment, its batch job, the R script and the manual workflow (M8-12 "What is
added", "Degenerate cases"; GIAP_SPEC.md §3.2, §3.7, §3.8, §11.5, §12, R6, R14).

Files checked: ``scripts/golden/golden_r.yml``, ``scripts/golden/golden_r.explicit.txt`` (skips
until the first Slurm solve commits it), ``scripts/golden/make_golden_r.R``,
``scripts/slurm/golden_r.sbatch`` and ``.github/workflows/golden-r.yml``. The R script's calls are
checked as text (whitespace removed): R is never run by the test suite.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "scripts" / "golden"
ENV_YML = GOLDEN / "golden_r.yml"
LOCK = GOLDEN / "golden_r.explicit.txt"
R_SCRIPT = GOLDEN / "make_golden_r.R"
SBATCH = REPO / "scripts" / "slurm" / "golden_r.sbatch"
WORKFLOW = REPO / ".github" / "workflows" / "golden-r.yml"
PINS = {
    "r-base": "4.5.3",
    "bioconductor-fgsea": "1.36.2",
    "bioconductor-goseq": "1.62.0",
    "r-biasedurn": "2.0.12",
}


def _text(path: Path) -> str:
    assert path.is_file(), f"missing {path.relative_to(REPO)}"
    return path.read_text()


def _env() -> dict[str, object]:
    data = yaml.safe_load(_text(ENV_YML))
    assert isinstance(data, dict)
    return data


def _deps() -> dict[str, str]:
    deps = _env()["dependencies"]
    assert isinstance(deps, list)
    assert all(isinstance(d, str) for d in deps), "no pip: section (conda packages only)"
    out = {}
    for dep in deps:
        m = re.match(r"^([A-Za-z0-9_.\-]+)\s*(?:==?\s*(\S+))?$", dep.strip())
        assert m, dep
        out[m.group(1)] = m.group(2) or ""
    return out


# --- env yml (§3.2: conda-forge and bioconda only) -----------------------------------------


def test_env_name() -> None:
    assert _env()["name"] == "giap-golden-r"


def test_env_channels_are_conda_forge_and_bioconda_only() -> None:
    channels = _env()["channels"]
    assert isinstance(channels, list)
    assert "defaults" not in channels
    assert {"conda-forge", "bioconda"} <= set(channels)
    assert set(channels) <= {"conda-forge", "bioconda", "nodefaults"}
    assert channels.index("conda-forge") < channels.index("bioconda")
    for dep in _env()["dependencies"]:  # type: ignore[union-attr]
        assert "::" not in str(dep) or str(dep).split("::")[0] in {"conda-forge", "bioconda"}


@pytest.mark.parametrize(("package", "version"), sorted(PINS.items()))
def test_env_pins(package: str, version: str) -> None:
    assert _deps().get(package) == version


def test_env_lists_data_table() -> None:
    assert "r-data.table" in _deps()


def test_lock_is_an_explicit_export_from_the_two_channels() -> None:
    if not LOCK.is_file():
        pytest.skip("golden_r.explicit.txt is written by the first giap-golden-r Slurm solve")
    lines = [ln for ln in LOCK.read_text().splitlines() if ln and not ln.startswith("#")]
    assert lines[0] == "@EXPLICIT"
    urls = lines[1:]
    assert urls
    for url in urls:
        assert re.match(r"^https://conda\.anaconda\.org/(conda-forge|bioconda)/", url), url
    for package, version in PINS.items():
        assert any(f"/{package}-{version}-" in u for u in urls), f"{package} {version} not locked"


# --- make_golden_r.R -----------------------------------------------------------------------


def _r_flat() -> str:
    return re.sub(r"\s+", "", _text(R_SCRIPT))


@pytest.mark.parametrize(
    "fragment",
    [
        "fgseaMultilevel(",
        "sampleSize=101",
        "minSize=10",
        "maxSize=500",
        "eps=0",
        "nPermSimple=1000",
        'scoreType="std"',
        "nproc=1",
        "calcGseaStat(",
        "goseq(",
        'method="Wallenius"',
        "use_genes_without_cat=TRUE",
        "bias.data=1",
        "digits=17",
        "sessionInfo()",
        "GOLDEN_SEED",
    ],
)
def test_r_script_calls(fragment: str) -> None:
    assert fragment in _r_flat()


def test_r_script_sets_the_seed_before_each_call() -> None:
    flat = _r_flat()
    assert flat.count("set.seed(") >= 2
    for call in ("fgseaMultilevel(", "goseq("):
        before = flat[: flat.index(call)]
        assert "set.seed(" in before, call


@pytest.mark.parametrize("package", ["fgsea", "goseq", "BiasedUrn"])
def test_r_script_records_package_versions(package: str) -> None:
    assert re.search(rf"packageVersion\(\s*[\"']{package}[\"']\s*\)", _text(R_SCRIPT))


def test_r_script_installs_nothing_from_outside_conda() -> None:
    text = _text(R_SCRIPT)
    for bad in ("install.packages", "BiocManager", "devtools", "remotes::", "download.file"):
        assert bad not in text, bad


def test_r_script_writes_the_versions_file() -> None:
    assert "VERSIONS.txt" in _text(R_SCRIPT)


# --- the Slurm job (§3.8, R14) -------------------------------------------------------------


def test_sbatch_parses_and_has_a_giap_header() -> None:
    text = _text(SBATCH)
    proc = subprocess.run(["bash", "-n", str(SBATCH)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert text.splitlines()[0] == "#!/bin/bash"
    assert re.search(r"^#SBATCH --job-name=giap-golden-r$", text, re.M)
    assert re.search(r"^#SBATCH --output=\S*/local/logs/%x\.%j\.out$", text, re.M)
    cpus = re.search(r"^#SBATCH --cpus-per-task=(\d+)$", text, re.M)
    assert cpus and int(cpus.group(1)) <= 4
    assert re.search(r"^#SBATCH --mem=16G$", text, re.M)
    assert "set -euo pipefail" in text
    assert "scripts/farmshare_env.sh" in text


def test_sbatch_builds_with_a_yml_marker_and_checks_versions() -> None:
    text = _text(SBATCH)
    assert "golden_r.yml" in text
    assert "micromamba" in text
    assert ".giap-built" in text, "the build_envs.sh marker pattern"
    assert "sha256" in text
    assert "packageVersion" in text, "the solved version must equal the pin"
    assert "golden_r.explicit.txt" in text, "the lock is exported"
    assert "make_golden_inputs.py" in text
    assert "make_golden_r.R" in text
    assert not re.search(r"(-c|--channel)\s+defaults\b", text)


# --- the manual workflow (§12) -------------------------------------------------------------


def _workflow() -> dict[object, object]:
    data = yaml.safe_load(_text(WORKFLOW))
    assert isinstance(data, dict)
    return data


def test_workflow_is_dispatch_only() -> None:
    wf = _workflow()
    triggers = wf.get("on", wf.get(True))  # PyYAML reads a bare `on:` key as True
    if isinstance(triggers, str):
        triggers = [triggers]
    assert triggers is not None
    assert set(triggers) == {"workflow_dispatch"}


def test_workflow_uses_the_lock_and_opens_a_pr() -> None:
    text = _text(WORKFLOW)
    assert "setup-micromamba" in text
    assert "golden_r.explicit.txt" in text
    assert "make_golden_inputs.py" in text
    assert "make_golden_r.R" in text
    assert "create-pull-request" in text or "gh pr create" in text


def test_workflow_never_pushes_to_main() -> None:
    text = _text(WORKFLOW)
    for line in text.splitlines():
        if "git push" in line:
            assert not re.search(r"\bmain\b", line), line
    assert not re.search(r"^\s*branch:\s*['\"]?main['\"]?\s*$", text, re.M)
