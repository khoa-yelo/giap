"""The capability plugin itself (GIPPA_SPEC.md §3.6)."""

from pathlib import Path

import pytest

CONFTEST = (Path(__file__).parents[1] / "conftest.py").read_text()

pytest_plugins = ["pytester"]


def test_missing_capability_skips_with_reason(pytester: pytest.Pytester) -> None:
    pytester.makeconftest(CONFTEST)
    pytester.makepyfile("import pytest\n@pytest.mark.needs_docker\ndef test_x():\n    pass\n")
    pytester.makeini("[pytest]\nmarkers =\n    needs_docker: docker\n")
    result = pytester.runpytest_subprocess("-p", "no:cacheprovider")
    if result.ret == 0 and "1 passed" in result.stdout.str():
        pytest.skip("docker is available here")
    result.assert_outcomes(skipped=1)
    result.stdout.fnmatch_lines(["*missing capability: needs_docker*"])


def test_required_capability_fails_instead_of_skipping(pytester: pytest.Pytester) -> None:
    pytester.makeconftest(CONFTEST)
    pytester.makepyfile("import pytest\n@pytest.mark.needs_docker\ndef test_x():\n    pass\n")
    pytester.makeini("[pytest]\nmarkers =\n    needs_docker: docker\n")
    result = pytester.runpytest_subprocess(
        "-p", "no:cacheprovider", "--require-capabilities=needs_docker"
    )
    if "1 passed" in result.stdout.str():
        pytest.skip("docker is available here")
    result.assert_outcomes(errors=1)
