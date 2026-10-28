"""M11-15g: permissive-only step image contents (``licences.allow`` + ``check_image_licences.py``).

Spec docs/specs/M11-15-portable-clusters.md item 8a "Contents rule (MIT)": no GPL/AGPL/
non-commercial/unlicensed Python package; ``scripts/release/check_image_licences.py`` "runs inside
the built image and fails on a Python distribution licence not in
``workflow/containers/step/licences.allow``"; item 8 AC 7 "Licence check script tested on a fake
site-packages (one GPL dist fails and is named)"; §13 (GIPPA stays MIT, Q-K34).

Readings chosen where the spec leaves room (the CLI is the contract; nothing is imported from it):

* ``python scripts/release/check_image_licences.py --allow FILE SITE_PACKAGES...`` reads every
  ``*.dist-info/METADATA`` under the given directories; exit 0 when every distribution's licence
  is allowed, 1 otherwise, naming each offending distribution (and its licence, or that it has
  none) on stdout or stderr, and naming no allowed distribution (a count of the checked ones is
  fine). It uses only the standard library (it runs in the image against
  whatever is installed there).
* A licence comes from ``License-Expression`` (PEP 639), else ``License``, else the
  ``License ::`` classifiers. An SPDX ``AND`` needs every operand allowed. A distribution with no
  licence information is unlicensed and fails.
* ``licences.allow``: lines, ``#`` comments and blanks ignored. A file holding one SPDX identifier
  per line is valid input (the implementer may add richer entries, e.g. per-distribution
  overrides for metadata that carries a full licence text). The committed file names ``MIT`` and
  no GPL-family, AGPL, non-commercial, proprietary, SSPL/BUSL or "unknown" entry.
* The real check: the METADATA of every distribution of ``requirements.lock`` installed in this
  venv, plus the package's own, passes against the committed allowlist (a distribution that a
  marker keeps off this platform is not installed here and not checked).
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
import tomllib
from importlib import metadata
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "check_image_licences.py"
ALLOW = REPO / "workflow" / "containers" / "step" / "licences.allow"
LOCK = REPO / "workflow" / "containers" / "step" / "requirements.lock"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
FORBIDDEN = re.compile(
    r"GPL|non-?commercial|\bNC\b|-NC\b|proprietary|\bSSPL\b|\bBUSL\b|commons.clause|unknown",
    re.IGNORECASE,
)
GPL_CLASSIFIER = "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
MIT_CLASSIFIER = "License :: OSI Approved :: MIT License"


# ---------------------------------------------------------------- helpers


def _dist(site: Path, name: str, version: str = "1.0", **fields: str | list[str]) -> None:
    info = site / f"{name.replace('-', '_')}-{version}.dist-info"
    info.mkdir(parents=True)
    lines = ["Metadata-Version: 2.4", f"Name: {name}", f"Version: {version}"]
    for key, value in fields.items():
        header = key.replace("_", "-")
        for v in value if isinstance(value, list) else [value]:
            lines.append(f"{header}: {v}")
    (info / "METADATA").write_text("\n".join(lines) + "\n\n", encoding="utf-8")


def _check(allow: Path, *sites: Path) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} is missing (M11-15g)"
    return subprocess.run(
        [sys.executable, "-I", str(SCRIPT), "--allow", str(allow), *map(str, sites)],
        capture_output=True,
        text=True,
        timeout=120,
    )


def _out(done: subprocess.CompletedProcess[str]) -> str:
    return done.stdout + done.stderr


def _allow_entries() -> list[str]:
    assert ALLOW.is_file(), f"{ALLOW.relative_to(REPO)} is missing (M11-15g)"
    lines = (ln.split("#", 1)[0].strip() for ln in ALLOW.read_text(encoding="utf-8").splitlines())
    return [ln for ln in lines if ln]


def _locked_names() -> list[str]:
    assert LOCK.is_file(), f"{LOCK.relative_to(REPO)} is missing (M11-15g)"
    return sorted(set(re.findall(r"^([A-Za-z0-9][A-Za-z0-9._\-]*)==", LOCK.read_text(), re.M)))


# ---------------------------------------------------------------- the allowlist


def test_allowlist_exists_and_names_mit() -> None:
    entries = _allow_entries()
    assert entries, "licences.allow is empty"
    assert any(re.search(r"(^|[\s:,])MIT($|[\s,])", e) for e in entries), entries


def test_allowlist_has_no_copyleft_noncommercial_or_unknown_entry() -> None:
    bad = [e for e in _allow_entries() if FORBIDDEN.search(e)]
    assert bad == [], bad


# ---------------------------------------------------------------- the script on a fake site


def test_script_uses_only_the_standard_library() -> None:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} is missing (M11-15g)"
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.add(node.module.split(".")[0])
    assert {n for n in names if n not in sys.stdlib_module_names and n != "__future__"} == set()


def test_permissive_site_passes(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakemit", License_Expression="MIT")
    _dist(site, "fakebsd", License_Expression="BSD-3-Clause")
    _dist(site, "fakeapache", License_Expression="Apache-2.0")
    _dist(site, "fakeclassifier", Classifier=[MIT_CLASSIFIER])
    done = _check(ALLOW, site)
    assert done.returncode == 0, _out(done)


def test_one_gpl_dist_fails_and_is_named(tmp_path: Path) -> None:
    """AC 7."""
    site = tmp_path / "site"
    _dist(site, "fakemit", License_Expression="MIT")
    _dist(site, "fakegpl", License_Expression="GPL-3.0-or-later")
    done = _check(ALLOW, site)
    assert done.returncode == 1, _out(done)
    assert "fakegpl" in _out(done), _out(done)
    assert "GPL-3.0" in _out(done), _out(done)
    assert "fakemit" not in _out(done).replace("fakegpl", ""), _out(done)


@pytest.mark.parametrize(
    "fields",
    [
        {"License_Expression": "AGPL-3.0-only"},
        {"License_Expression": "LGPL-2.1-or-later"},
        {"License_Expression": "MIT AND GPL-2.0-only"},
        {"License_Expression": "CC-BY-NC-4.0"},
        {"License": "GPL-3.0"},
        {"Classifier": [GPL_CLASSIFIER]},
        {},
    ],
    ids=["agpl", "lgpl", "and-gpl", "non-commercial", "license-field", "classifier", "none"],
)
def test_forbidden_or_unlicensed_dist_fails_and_is_named(
    tmp_path: Path, fields: dict[str, str | list[str]]
) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    _dist(site, "fakebad", **fields)
    done = _check(ALLOW, site)
    assert done.returncode == 1, _out(done)
    assert "fakebad" in _out(done), _out(done)


def test_the_allowlist_file_is_the_source(tmp_path: Path) -> None:
    """With an allowlist of MIT only, an Apache-2.0 dist fails; with the committed one it passes."""
    site = tmp_path / "site"
    _dist(site, "fakeapache", License_Expression="Apache-2.0")
    only_mit = tmp_path / "only-mit.allow"
    only_mit.write_text("# test allowlist\nMIT\n", encoding="utf-8")
    done = _check(only_mit, site)
    assert done.returncode == 1 and "fakeapache" in _out(done), _out(done)
    assert _check(ALLOW, site).returncode == 0


def test_every_site_directory_is_checked(tmp_path: Path) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    _dist(a, "fakemit", License_Expression="MIT")
    _dist(b, "fakegpl", License_Expression="GPL-2.0-only")
    done = _check(ALLOW, a, b)
    assert done.returncode == 1 and "fakegpl" in _out(done), _out(done)


LICENCES = {
    "MIT": True,
    "BSD-3-Clause": True,
    "Apache-2.0": True,
    "GPL-3.0-only": False,
    "AGPL-3.0-or-later": False,
    "GPL-2.0-or-later": False,
}


@settings(
    max_examples=15,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(st.lists(st.sampled_from(sorted(LICENCES)), min_size=1, max_size=6))
def test_property_fails_iff_a_forbidden_dist_and_names_exactly_those(
    tmp_path_factory: pytest.TempPathFactory, licences: list[str]
) -> None:
    site = tmp_path_factory.mktemp("site")
    for i, lic in enumerate(licences):
        _dist(site, f"dist{i:03d}x", License_Expression=lic)
    done = _check(ALLOW, site)
    bad = {f"dist{i:03d}x" for i, lic in enumerate(licences) if not LICENCES[lic]}
    assert done.returncode == (1 if bad else 0), _out(done)
    named = set(re.findall(r"dist\d{3}x", _out(done)))
    assert named == bad, (named, bad, _out(done))


# ---------------------------------------------------------------- the real resolved set


def test_the_locked_dependency_set_passes_the_committed_allowlist(tmp_path: Path) -> None:
    site = tmp_path / "site"
    checked = []
    for name in [*_locked_names(), PKG]:
        try:
            dist = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue  # kept off this platform by a marker: not installed here, not in the image
        text = dist.read_text("METADATA")
        assert text, f"{name}: no METADATA"
        info = site / f"{name.replace('-', '_')}-{dist.version}.dist-info"
        info.mkdir(parents=True)
        (info / "METADATA").write_text(text, encoding="utf-8")
        checked.append(name)
    assert len(checked) >= 20, checked
    done = _check(ALLOW, site)
    assert done.returncode == 0, _out(done)
