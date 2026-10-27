"""``scripts/release/bioconda_recipe.py`` (M11-10 definition items 1 and 3, examples 1, 3 and 4,
acceptance 1 and 3).

Interface read from the spec: ``python scripts/release/bioconda_recipe.py --version V --sha256 H
[--doi D] --out DIR`` renders ``packaging/bioconda/meta.yaml.j2`` into ``DIR/gippa/meta.yaml``; a
bad ``--sha256`` exits 1 with ``sha256 must be 64 hex characters``. ``--sha256 auto`` (network) is
never exercised: nothing here contacts PyPI or anaconda.org.

The rendered file must parse with PyYAML (acceptance 3), so the conda-build variable in
``build: script`` is read as the literal string ``{{ PYTHON }} -m pip install ...`` (quoted in
YAML); conda-build substitutes it later.
"""

from __future__ import annotations

import os
import string
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "bioconda_recipe.py"
TEMPLATE = REPO / "packaging" / "bioconda" / "meta.yaml.j2"
SHA = "0123456789abcdef" * 4
PIP_SCRIPT = "{{ PYTHON }} -m pip install . --no-deps --no-build-isolation -vv"


def run_render(out: Path, *args: str, cwd: Path = REPO) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), env.get("PYTHONPATH", "")])
    return subprocess.run(
        [
            sys.executable,
            str(cwd / "scripts" / "release" / "bioconda_recipe.py"),
            *args,
            "--out",
            str(out),
        ],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )


def render_ok(out: Path, *args: str) -> dict[str, Any]:
    proc = run_render(out, *args)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    meta = out / "gippa" / "meta.yaml"
    assert meta.is_file(), f"expected {meta}"
    data = yaml.safe_load(meta.read_text())
    assert isinstance(data, dict)
    return data


def test_script_and_template_exist() -> None:
    assert SCRIPT.is_file()
    assert TEMPLATE.is_file()


def test_template_notes_codon_bias_is_bioconda_only() -> None:
    """Edge case: the template comment records that codon-bias exists only on bioconda."""
    lines = TEMPLATE.read_text().lower().splitlines()
    hits = [i for i, ln in enumerate(lines) if "codon-bias" in ln]
    assert hits, "template should mention codon-bias"
    assert any("bioconda" in ln for i in hits for ln in lines[max(0, i - 2) : i + 3]), (
        "template should note next to codon-bias that it is bioconda-only"
    )


# Example 1 and acceptance 3: the keys of definition item 1.


@pytest.fixture(scope="module")
def recipe(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return render_ok(tmp_path_factory.mktemp("recipe"), "--version", "1.0.0", "--sha256", SHA)


def test_package_section(recipe: dict[str, Any]) -> None:
    assert recipe["package"]["name"] == "gippa"
    assert str(recipe["package"]["version"]) == "1.0.0"


def test_source_is_pypi_sdist_with_sha(recipe: dict[str, Any]) -> None:
    assert recipe["source"]["url"] == "https://pypi.org/packages/source/g/gippa/gippa-1.0.0.tar.gz"
    assert recipe["source"]["sha256"] == SHA


def test_build_section(recipe: dict[str, Any]) -> None:
    build = recipe["build"]
    assert build["noarch"] == "python"
    assert build["number"] == 0
    assert build["script"].strip() == PIP_SCRIPT
    assert build["entry_points"] == ["gippa = gippa.cli:main"]
    # lead decision (#345 review B1): bioconda-utils requires run_exports on new recipes
    assert build["run_exports"] == ['{{ pin_subpackage("gippa", max_pin="x") }}']


def test_host_requirements(recipe: dict[str, Any]) -> None:
    host = [str(h).strip() for h in recipe["requirements"]["host"]]
    assert "python >=3.11" in host
    assert "pip" in host
    assert any(h.split()[0] == "hatchling" for h in host)


def test_run_has_python_nextflow_openjdk(recipe: dict[str, Any]) -> None:
    run = [str(r).strip() for r in recipe["requirements"]["run"]]
    assert "python >=3.11" in run
    assert "nextflow >=25.04" in run  # spec §3 Nextflow floor (#345 review S4)
    assert "openjdk >=17" in run


def test_test_section(recipe: dict[str, Any]) -> None:
    test = recipe["test"]
    assert "gippa --version" in test["commands"]
    assert "gippa --help" in test["commands"]
    assert "pip check" in test["commands"] and "pip" in test["requires"]
    assert set(test["imports"]) >= {"gippa", "gippa.cli"}


def test_about_licence(recipe: dict[str, Any]) -> None:
    about = recipe["about"]
    assert about["license"] == "MIT"
    assert about["license_file"] == "LICENSE"
    assert (REPO / "LICENSE").is_file()
    for key in ("home", "summary", "doc_url", "dev_url"):
        assert isinstance(about.get(key), str) and about[key].strip(), key
    assert "github.com/khoa-yelo/gippa" in about["dev_url"]


def test_extra_maintainers_and_no_doi_by_default(recipe: dict[str, Any]) -> None:
    extra = recipe["extra"]
    assert extra["recipe-maintainers"] == ["khoa-yelo"]
    assert not extra.get("identifiers")


def test_no_defaults_channel(tmp_path: Path) -> None:
    """§3.2: conda-forge and bioconda only; the recipe never names the ``defaults`` channel."""
    render_ok(tmp_path, "--version", "1.0.0", "--sha256", SHA)
    lines = (tmp_path / "gippa" / "meta.yaml").read_text().splitlines()
    code = [ln.split("#", 1)[0] for ln in lines]
    assert not any("defaults" in ln for ln in code)


# Example 4.


def test_doi_identifier(tmp_path: Path) -> None:
    data = render_ok(
        tmp_path, "--version", "1.0.0", "--sha256", SHA, "--doi", "10.5281/zenodo.1234567"
    )
    assert data["extra"]["identifiers"] == ["doi:10.5281/zenodo.1234567"]
    assert data["extra"]["recipe-maintainers"] == ["khoa-yelo"]


# Example 3 and the sha256 placeholder.


@pytest.mark.parametrize("bad", ["xyz", SHA[:-1], SHA + "0", "g" * 64, ""])
def test_bad_sha256_exits_1(tmp_path: Path, bad: str) -> None:
    proc = run_render(tmp_path, "--version", "1.0.0", "--sha256", bad)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "sha256 must be 64 hex characters" in proc.stdout + proc.stderr
    assert not (tmp_path / "gippa" / "meta.yaml").exists()


@pytest.mark.parametrize("bad", ['1.0"x', "../../evil", "", "one"])
def test_invalid_version_exits_1(tmp_path: Path, bad: str) -> None:
    """#345 review S1: a version that is not PEP 440 never reaches the URL or the YAML."""
    proc = run_render(tmp_path, "--version", bad, "--sha256", SHA)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "invalid version" in proc.stdout + proc.stderr
    assert not (tmp_path / "gippa" / "meta.yaml").exists()


def test_version_is_normalised_like_the_sdist(tmp_path: Path) -> None:
    data = render_ok(tmp_path, "--version", "1.0.0-rc1", "--sha256", SHA)
    assert data["package"]["version"] == "1.0.0rc1"
    assert data["source"]["url"].endswith("/gippa-1.0.0rc1.tar.gz")


def test_uppercase_hex_sha_accepted(tmp_path: Path) -> None:
    """64 hex characters in either case is a valid digest (reading: hex is case-insensitive)."""
    data = render_ok(tmp_path, "--version", "1.0.0", "--sha256", SHA.upper())
    assert data["source"]["sha256"] == SHA  # written lowercased


# Properties: version and sha256 placeholders are substituted verbatim.

_versions = st.builds(
    lambda a, b, c, rc: f"{a}.{b}.{c}" + (f"rc{rc}" if rc is not None else ""),
    st.integers(0, 30),
    st.integers(0, 30),
    st.integers(0, 30),
    st.none() | st.integers(1, 9),
)
_shas = st.text(alphabet="0123456789abcdef", min_size=64, max_size=64)


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(version=_versions, sha=_shas)
def test_property_placeholders_substituted(tmp_path: Path, version: str, sha: str) -> None:
    out = tmp_path / f"{version}-{sha[:8]}"
    data = render_ok(out, "--version", version, "--sha256", sha)
    assert str(data["package"]["version"]) == version
    assert data["source"]["sha256"] == sha
    assert (
        data["source"]["url"] == f"https://pypi.org/packages/source/g/gippa/gippa-{version}.tar.gz"
    )
    text = (out / "gippa" / "meta.yaml").read_text()
    assert "{{ version" not in text and "{{ sha256" not in text


_bad_shas = st.one_of(
    st.text(alphabet="0123456789abcdef", max_size=80).filter(lambda s: len(s) != 64),
    st.text(alphabet=string.ascii_letters + string.digits, min_size=64, max_size=64).filter(
        lambda s: any(ch not in string.hexdigits for ch in s)
    ),
).filter(lambda s: s != "auto" and not s.startswith("-"))


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(bad=_bad_shas)
def test_property_bad_sha_rejected(tmp_path: Path, bad: str) -> None:
    out = tmp_path / f"bad-{abs(hash(bad))}"
    proc = run_render(out, "--version", "1.0.0", "--sha256", bad)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "sha256 must be 64 hex characters" in proc.stdout + proc.stderr


def test_all_digit_sha_stays_a_string(tmp_path: Path) -> None:
    """A digest of digits only must survive YAML as the same string (it has to be quoted)."""
    digits = "0123456789" * 6 + "0123"
    data = render_ok(tmp_path, "--version", "1.0.0", "--sha256", digits)
    assert data["source"]["sha256"] == digits


def test_two_part_version_stays_a_string(tmp_path: Path) -> None:
    """``1.10`` must not become the float 1.1 when the recipe is parsed."""
    data = render_ok(tmp_path, "--version", "1.10", "--sha256", SHA)
    assert data["package"]["version"] == "1.10"
