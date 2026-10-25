"""Render the Bioconda recipe for GIAP from ``pyproject.toml`` (M11-10).

Usage: python scripts/release/bioconda_recipe.py --version V --sha256 H [--doi D] --out DIR

Renders ``packaging/bioconda/meta.yaml.j2`` into ``DIR/giap/meta.yaml``. The runtime dependencies
come from ``[project].dependencies`` of the repository's ``pyproject.toml`` (located relative to
this script), each mapped to its conda name. ``--sha256 auto`` reads the sdist digest from PyPI
(network; only after the release). Exit 1 on a bad sha256, an unmapped dependency or an
environment marker; nothing is written then.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
import urllib.request
from pathlib import Path

import jinja2
from packaging.requirements import InvalidRequirement, Requirement

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = Path("packaging") / "bioconda" / "meta.yaml.j2"
SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")

# PyPI name -> conda name. Every runtime dependency must be listed (a new one fails the render).
CONDA_NAME: dict[str, str] = {
    name: name
    for name in (
        "numpy",
        "pandas",
        "pyarrow",
        "numba",
        "scipy",
        "scikit-learn",
        "statsmodels",
        "lightgbm",
        "pyahocorasick",
        "pydantic",
        "click",
        "pyyaml",
        "jinja2",
        "matplotlib",
        "pycirclize",
        "pygenomeviz",
        "upsetplot",
        "plotly",
        "goatools",
        "codon-bias",  # bioconda only
        "biopython",
    )
}


class RecipeError(ValueError):
    """A dependency or argument that cannot go into the recipe."""


def conda_requirement(dep: str) -> str:
    """``pydantic>=2`` -> ``pydantic >=2`` with the conda name; markers and unknown names fail."""
    try:
        req = Requirement(dep)
    except InvalidRequirement as exc:
        raise RecipeError(f"cannot parse dependency {dep!r}: {exc}") from exc
    if req.marker is not None:
        raise RecipeError(
            f"{req.name} has an environment marker ({req.marker}); "
            "a noarch: python recipe cannot express it"
        )
    name = CONDA_NAME.get(req.name.lower())
    if name is None:
        raise RecipeError(f"no conda name for {req.name} (add it to the map)")
    spec = str(req.specifier)
    return f"{name} {spec}" if spec else name


def check_sha256(sha: str) -> str:
    """The digest lowercased; anything but 64 hex characters fails."""
    if not SHA256.match(sha):
        raise RecipeError("sha256 must be 64 hex characters")
    return sha.lower()


def pypi_sha256(version: str) -> str:
    """The sdist sha256 of ``giap==version`` on PyPI (network)."""
    url = f"https://pypi.org/pypi/giap/{version}/json"
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = json.load(resp)
    for f in data["urls"]:
        if f["packagetype"] == "sdist":
            return str(f["digests"]["sha256"])
    raise RecipeError(f"no sdist for giap {version} on PyPI")


def render(
    pyproject: dict[str, object], template: str, version: str, sha256: str, doi: str | None
) -> str:
    """The recipe text for one version (pure)."""
    project = pyproject["project"]
    assert isinstance(project, dict)
    run = [conda_requirement(d) for d in project["dependencies"]]
    env = jinja2.Environment(undefined=jinja2.StrictUndefined, keep_trailing_newline=True)
    return env.from_string(template).render(
        version=version,
        sha256=sha256,
        doi=doi,
        run=run,
        summary=str(project.get("description", "")).replace('"', "'"),
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--version", required=True)
    ap.add_argument("--sha256", required=True, help="64 hex characters, or 'auto' (PyPI)")
    ap.add_argument("--doi", default=None)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    try:
        sha = pypi_sha256(args.version) if args.sha256 == "auto" else args.sha256
        sha = check_sha256(sha)
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
        text = render(pyproject, (ROOT / TEMPLATE).read_text(), args.version, sha, args.doi)
    except RecipeError as exc:
        print(f"bioconda_recipe: {exc}", file=sys.stderr)
        return 1
    out = args.out / "giap" / "meta.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
