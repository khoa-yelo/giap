"""M11-15g review round 1, B2 and B3: ``check_image_licences.py`` must not fail open.

Review of PR #436 (B2, B3). The CLI is the contract; nothing is imported from the script.

* B2(a): a site argument that is not a directory is a usage error, exit 2. A run that checks no
  distribution exits 1. ``--require NAME`` (repeatable; the 15n release job passes the package's
  own name) exits 1 and names ``NAME`` when no checked distribution has that name (PEP 503
  normalised). Reading chosen here: the own-package check is a flag, so a fake site without the
  package still passes as in the existing tests.
* B2(b), (c): ``License-Expression`` (PEP 639) is authoritative: when present, nothing else is
  read. Without it, every signal present must pass: a one-line ``License`` that parses as an
  expression, and every ``License ::`` classifier. An unparsable ``License`` is not a signal
  when classifiers exist (``Dual License``, ``2-clause BSD``), and fails alone. ``License`` text
  (one line or a full licence text) matching ``GPL|General Public|Affero|non-?commercial``
  fails unless ``licences.allow`` overrides that distribution by name; renaming a real
  distribution whose metadata needs such an override makes it fail again.
* B3: one case per mutant that survived the first test suite; each checks the exit code and that
  the offender (and no passing distribution) is named.

``GIPPA_LICENCE_SCRIPT`` may point the tests at a mutated copy of the script (mutation runs).
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from importlib import metadata
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = Path(
    os.environ.get("GIPPA_LICENCE_SCRIPT", REPO / "scripts" / "release" / "check_image_licences.py")
)
ALLOW = REPO / "workflow" / "containers" / "step" / "licences.allow"
MIT_CLASSIFIER = "License :: OSI Approved :: MIT License"
BSD_CLASSIFIER = "License :: OSI Approved :: BSD License"
GPL_CLASSIFIER = "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
GPL_TEXT = (
    "Copyright (c) 2020 Fake Authors\n"
    "        This program is free software: you can redistribute it and/or modify\n"
    "        it under the terms of the GNU General Public License as published by\n"
    "        the Free Software Foundation, either version 3 of the License."
)


def _dist(site: Path, name: str, version: str = "1.0", **fields: str | list[str]) -> None:
    info = site / f"{name.replace('-', '_')}-{version}.dist-info"
    info.mkdir(parents=True)
    lines = ["Metadata-Version: 2.4", f"Name: {name}", f"Version: {version}"]
    for key, value in fields.items():
        header = key.replace("_", "-")
        for v in value if isinstance(value, list) else [value]:
            lines.append(f"{header}: {v}")
    (info / "METADATA").write_text("\n".join(lines) + "\n\n", encoding="utf-8")


def _check(*args: str | Path) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), SCRIPT
    return subprocess.run(
        [sys.executable, "-I", str(SCRIPT), *map(str, args)],
        capture_output=True,
        text=True,
        timeout=120,
    )


def _out(done: subprocess.CompletedProcess[str]) -> str:
    return done.stdout + done.stderr


def _site_with(tmp_path: Path, **fields: str | list[str]) -> Path:
    """A site with one permissive control (``fakeok``) and the case under test (``fakecase``)."""
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    _dist(site, "fakecase", **fields)
    return site


def _assert_fails_naming_case(done: subprocess.CompletedProcess[str]) -> None:
    out = _out(done)
    assert done.returncode == 1, out
    assert "fakecase" in out, out
    assert "fakeok" not in out, out


def _assert_passes(done: subprocess.CompletedProcess[str]) -> None:
    out = _out(done)
    assert done.returncode == 0, out
    assert "fakecase" not in out and "fakeok" not in out, out


# ---------------------------------------------------------------- B2(a): fail closed on paths


def test_missing_site_directory_is_a_usage_error(tmp_path: Path) -> None:
    done = _check("--allow", ALLOW, tmp_path / "python3.11" / "site-packages")
    assert done.returncode == 2, _out(done)
    assert "python3.11" in _out(done), _out(done)


def test_site_that_is_a_file_is_a_usage_error(tmp_path: Path) -> None:
    not_dir = tmp_path / "site-packages"
    not_dir.write_text("", encoding="utf-8")
    done = _check("--allow", ALLOW, not_dir)
    assert done.returncode == 2, _out(done)


def test_one_missing_site_among_good_ones_is_a_usage_error(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    done = _check("--allow", ALLOW, site, tmp_path / "typo")
    assert done.returncode == 2, _out(done)


def test_a_run_that_checks_nothing_fails(tmp_path: Path) -> None:
    empty = tmp_path / "site"
    empty.mkdir()
    (empty / "README.txt").write_text("not a dist\n", encoding="utf-8")
    done = _check("--allow", ALLOW, empty)
    assert done.returncode == 1, _out(done)


def test_require_names_a_missing_own_distribution(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    done = _check("--allow", ALLOW, "--require", "gippa", site)
    assert done.returncode == 1, _out(done)
    assert "gippa" in _out(done), _out(done)


def test_require_passes_when_the_distribution_is_checked(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    _dist(site, "Fake_Pkg", License_Expression="MIT")
    done = _check("--allow", ALLOW, "--require", "fake-pkg", "--require", "fakeok", site)
    assert done.returncode == 0, _out(done)


# ---------------------------------------------------------------- B2(b), (c): every signal


@pytest.mark.parametrize(
    "fields",
    [
        {"License": "GNU General Public License v3", "Classifier": [MIT_CLASSIFIER]},
        {"License": "GPLv3", "Classifier": [MIT_CLASSIFIER]},
        {"License": "GNU Affero General Public License", "Classifier": [BSD_CLASSIFIER]},
        {"License": "Free for non-commercial use", "Classifier": [MIT_CLASSIFIER]},
        {"License": "Noncommercial research licence", "Classifier": [MIT_CLASSIFIER]},
        {"License": GPL_TEXT, "Classifier": [BSD_CLASSIFIER]},
        {"License": GPL_TEXT},
    ],
    ids=[
        "gpl-prose",
        "gplv3",
        "affero",
        "non-commercial",
        "noncommercial",
        "gpl-text",
        "text-only",
    ],
)
def test_b2b_gpl_or_noncommercial_licence_prose_fails(
    tmp_path: Path, fields: dict[str, str | list[str]]
) -> None:
    _assert_fails_naming_case(_check("--allow", ALLOW, _site_with(tmp_path, **fields)))


def test_b2c_parsable_licence_field_does_not_hide_a_gpl_classifier(tmp_path: Path) -> None:
    site = _site_with(tmp_path, License="MIT", Classifier=[GPL_CLASSIFIER])
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


def test_b2c_disallowed_licence_field_fails_despite_allowed_classifiers(tmp_path: Path) -> None:
    site = _site_with(tmp_path, License="LicenseRef-Proprietary", Classifier=[MIT_CLASSIFIER])
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


def test_b2c_reproducer_two_distributions_both_named(tmp_path: Path) -> None:
    """The review's reproduction: two such distributions, exit 0 before the fix."""
    site = tmp_path / "site"
    _dist(site, "fakeone", License="MIT", Classifier=[GPL_CLASSIFIER])
    _dist(site, "faketwo", License="BSD-3-Clause", Classifier=[GPL_CLASSIFIER, BSD_CLASSIFIER])
    done = _check("--allow", ALLOW, site)
    assert done.returncode == 1, _out(done)
    assert "fakeone" in _out(done) and "faketwo" in _out(done), _out(done)


def test_licence_expression_is_authoritative(tmp_path: Path) -> None:
    """PEP 639: with License-Expression, the legacy License field and classifiers are not read."""
    site = _site_with(
        tmp_path, License_Expression="MIT", License=GPL_TEXT, Classifier=[GPL_CLASSIFIER]
    )
    _assert_passes(_check("--allow", ALLOW, site))


@pytest.mark.parametrize(
    "fields",
    [
        {"License": "MIT", "Classifier": [MIT_CLASSIFIER]},
        {"License": "Dual License", "Classifier": [BSD_CLASSIFIER]},
        {"License": "2-clause BSD", "Classifier": [BSD_CLASSIFIER]},
        {
            "License": "Copyright (c) 2015 Fake\n        All rights reserved, BSD terms.",
            "Classifier": [BSD_CLASSIFIER],
        },
        {"License": "BSD-3-Clause and Public-Domain"},
    ],
    ids=["agree", "unparsable-with-classifier", "prose-bsd", "text-bsd", "lowercase-and"],
)
def test_permissive_signals_still_pass(tmp_path: Path, fields: dict[str, str | list[str]]) -> None:
    _assert_passes(_check("--allow", ALLOW, _site_with(tmp_path, **fields)))


@pytest.mark.parametrize("real", ["matplotlib", "pandas", "scipy"])
def test_the_gpl_prose_override_is_by_name(tmp_path: Path, real: str) -> None:
    """These locked distributions carry GPL words in their bundled licence texts, so they pass
    only through a by-name override (the locked-set test); renamed, the same METADATA fails."""
    try:
        text = metadata.distribution(real).read_text("METADATA")
    except metadata.PackageNotFoundError:
        pytest.skip(f"{real} is not installed here")
    assert text
    renamed = re.sub(r"(?m)^Name: .*$", "Name: fakecase", text, count=1)
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    info = site / "fakecase-1.0.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(renamed, encoding="utf-8")
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


# ---------------------------------------------------------------- B3: one case per mutant


def test_mutant_all_to_any_dual_classifiers_fail(tmp_path: Path) -> None:
    site = _site_with(tmp_path, Classifier=[MIT_CLASSIFIER, GPL_CLASSIFIER])
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


@pytest.mark.parametrize("expr", ["MIT AND", "MIT OR OR BSD-3-Clause", "(MIT", "MIT, BSD"])
def test_mutant_unparsable_licence_expression_fails(tmp_path: Path, expr: str) -> None:
    site = _site_with(tmp_path, License_Expression=expr, Classifier=[MIT_CLASSIFIER])
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


def test_mutant_unparsable_licence_alone_fails(tmp_path: Path) -> None:
    site = _site_with(tmp_path, License="Custom permissive terms, see LICENSE")
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


def test_mutant_multiline_licence_is_not_read_as_an_expression(tmp_path: Path) -> None:
    """``MIT OR\\n Proprietary`` parses to an allowed OR; it is licence text, not an expression,
    and with no classifier there is no other signal."""
    site = _site_with(tmp_path, License="MIT OR\n        LicenseRef-Proprietary")
    _assert_fails_naming_case(_check("--allow", ALLOW, site))


def test_mutant_or_is_not_and(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="MIT")
    _dist(site, "fakeor", License_Expression="GPL-3.0-only OR MIT")
    _dist(site, "fakecase", License_Expression="GPL-3.0-only OR LicenseRef-Proprietary")
    done = _check("--allow", ALLOW, site)
    _assert_fails_naming_case(done)
    assert "fakeor" not in _out(done), _out(done)


def test_mutant_with_exception_is_matched_whole(tmp_path: Path) -> None:
    site = tmp_path / "site"
    _dist(site, "fakeok", License_Expression="Apache-2.0 WITH LLVM-exception")
    _dist(site, "fakecase", License_Expression="MIT WITH Fake-exception")
    _assert_fails_naming_case(_check("--allow", ALLOW, site))
