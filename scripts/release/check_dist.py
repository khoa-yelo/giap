"""Check the wheel and sdist in DIST_DIR before a release (M11-01 item 4; R5, R11, R12).

Usage: python scripts/release/check_dist.py DIST_DIR [--version V]

Opens every ``*.whl`` and ``*.tar.gz`` with the standard library (never installs) and checks the
required content, forbidden paths, symlinks, sizes (> 5 MB) and the version (default: the
checkout's ``giap.__version__``). Prints one ``ok <artefact>`` line per clean artefact and one
``FAIL <artefact>: <rule>: <path>`` line per violation. Exit 0 clean, 1 violations, 2 when
DIST_DIR has no wheel or no sdist.
"""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parents[2]
MAX_BYTES = 5 * 1024 * 1024  # R11
FORBIDDEN_TOP = ("local", "tests", "validation", ".github", ".claude", "docs", "benchmarks")
FORBIDDEN_PREFIXES = ("workflow/tests/", "giap/workflow/tests/")
FORBIDDEN_SUFFIXES = (".sif", ".raw", ".d", ".mzml", ".parquet", ".pyc")
UNIOP_ALLOWED = ("workflow/modules/local/uniop/", "workflow/containers/uniop/")
SDIST_FILES = (
    "pyproject.toml",
    "README.md",
    "LICENSE",
    "NOTICE.md",
    "CHANGELOG.md",
    "CITATION.cff",
)
SDIST_TREES = ("src/giap/__init__.py", "workflow/main.nf")
WHEEL_FILES = (
    "giap/__init__.py",
    "giap/cli.py",
    "giap/py.typed",
    "giap/config/defaults.yaml",
    "giap/schemas/datapackage.json",
    "giap/data/calibration.yaml",
    "giap/data/db_sources.yaml",
    "giap/workflow/main.nf",
    "giap/workflow/nextflow.config",
)


@dataclass(frozen=True)
class Member:
    """One archive entry: its archive path, the path below the artefact root, size and kind."""

    path: str
    rel: str
    size: int
    symlink: bool


def checkout_version(repo: Path = REPO) -> str:
    """``__version__`` of ``src/giap/__init__.py`` (read, not imported)."""
    m = re.search(
        r'^__version__\s*=\s*"([^"]+)"', (repo / "src/giap/__init__.py").read_text(), re.M
    )
    if m is None:
        raise SystemExit("no __version__ in src/giap/__init__.py")
    return m.group(1)


def wheel_required(repo: Path = REPO) -> list[str]:
    """Item-4 wheel files plus every ``workflow/envs/*.yml`` and ``workflow/conf/*.config``."""
    envs = sorted(f"giap/workflow/envs/{p.name}" for p in (repo / "workflow/envs").glob("*.yml"))
    confs = sorted(
        f"giap/workflow/conf/{p.name}" for p in (repo / "workflow/conf").glob("*.config")
    )
    return [*WHEEL_FILES, *envs, *confs]


def forbidden(rel: str) -> bool:
    """R5/R11 path rules on a root-relative POSIX path."""
    parts = PurePosixPath(rel).parts
    if not parts:
        return False
    if parts[0] in FORBIDDEN_TOP or rel.startswith(FORBIDDEN_PREFIXES):
        return True
    if "__pycache__" in parts or any(p.lower().endswith(FORBIDDEN_SUFFIXES) for p in parts):
        return True
    if "regulondb" in parts[-1].lower():
        return True
    in_uniop = "uniop" in parts[:-1]
    return in_uniop and not rel.removeprefix("giap/").startswith(UNIOP_ALLOWED)


def member_failures(members: list[Member]) -> Iterator[tuple[str, str]]:
    """(rule, path) for symlinks, forbidden paths and files > 5 MB; one per member at most."""
    for m in members:
        if m.symlink:
            yield "symlink", m.path
        elif forbidden(m.rel):
            yield "forbidden path", m.path
        elif m.size > MAX_BYTES:
            yield "file > 5 MB", m.path


def header(text: str, key: str) -> list[str]:
    """Values of an RFC 822 style ``key:`` header (METADATA, PKG-INFO)."""
    return [ln.split(":", 1)[1].strip() for ln in text.splitlines() if ln.startswith(f"{key}:")]


def version_failures(
    name_v: str, meta_v: list[str], want: str, meta_file: str
) -> list[tuple[str, str]]:
    """One ``version`` failure when the file name or the metadata does not carry ``want``."""
    if name_v == want and meta_v == [want]:
        return []
    return [("version", f"want {want}, file name has {name_v}, {meta_file} has {meta_v}")]


def check_wheel(path: Path, version: str, repo: Path = REPO) -> list[tuple[str, str]]:
    """Every violation of one wheel as (rule, path)."""
    with zipfile.ZipFile(path) as zf:
        infos = [i for i in zf.infolist() if not i.is_dir()]
        members = [
            Member(
                i.filename, i.filename, i.file_size, (i.external_attr >> 16) & 0o170000 == 0o120000
            )
            for i in infos
        ]
        names = {m.path for m in members}
        dist_info = sorted(
            {n.split("/", 1)[0] for n in names if n.split("/", 1)[0].endswith(".dist-info")}
        )
        di = dist_info[0] if dist_info else "giap.dist-info"

        def read(name: str) -> str:
            return zf.read(name).decode(errors="replace") if name in names else ""

        metadata = read(f"{di}/METADATA")
        entry_points = read(f"{di}/entry_points.txt")
    fails = list(member_failures(members))
    fails += [("missing", f) for f in wheel_required(repo) if f not in names]
    if not metadata:
        fails.append(("missing", f"{di}/METADATA"))
    if not re.search(
        r"^\[console_scripts\][^\[]*^giap\s*=\s*giap\.cli:main\s*$", entry_points, re.M | re.S
    ):
        fails.append(("console script", f"{di}/entry_points.txt"))
    licences = header(metadata, "License-Expression") + header(metadata, "License")
    if "MIT" not in licences:
        fails.append(("licence not MIT", f"{di}/METADATA"))
    if not {f"{di}/licenses/LICENSE", f"{di}/LICENSE"} & names:
        fails.append(("missing", f"{di}/licenses/LICENSE"))
    name_v = path.name.split("-")[1] if path.name.count("-") >= 2 else ""
    return fails + version_failures(name_v, header(metadata, "Version"), version, "METADATA")


def check_sdist(path: Path, version: str) -> list[tuple[str, str]]:
    """Every violation of one sdist as (rule, path); paths below the ``giap-V/`` root count."""
    with tarfile.open(path, "r:gz") as tf:
        tinfos = [t for t in tf.getmembers() if not t.isdir()]
        roots = sorted({t.name.split("/", 1)[0] for t in tinfos})
        root = roots[0] if roots else ""
        members = [
            Member(
                t.name,
                t.name.split("/", 1)[1] if "/" in t.name else t.name,
                t.size,
                t.issym() or t.islnk(),
            )
            for t in tinfos
        ]
        pkg_info = f"{root}/PKG-INFO"
        meta = tf.extractfile(pkg_info) if pkg_info in tf.getnames() else None
        metadata = meta.read().decode(errors="replace") if meta else ""
    rels = {m.rel for m in members if "/" in m.path}
    fails = list(member_failures(members))
    fails += [("missing", f"{root}/{f}") for f in (*SDIST_FILES, *SDIST_TREES) if f not in rels]
    name_v = path.name.removeprefix("giap-").removesuffix(".tar.gz")
    return fails + version_failures(name_v, header(metadata, "Version"), version, "PKG-INFO")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dist_dir", type=Path)
    ap.add_argument(
        "--version", help="expected version (default: giap.__version__ of the checkout)"
    )
    args = ap.parse_args(argv)
    version = args.version or checkout_version()
    wheels = sorted(args.dist_dir.glob("*.whl"))
    sdists = sorted(args.dist_dir.glob("*.tar.gz"))
    if not wheels or not sdists:
        missing = " and ".join(k for k, v in (("no wheel", wheels), ("no sdist", sdists)) if not v)
        print(f"{missing} in {args.dist_dir}")
        return 2
    bad = False
    for art in [*wheels, *sdists]:
        fails = check_wheel(art, version) if art.suffix == ".whl" else check_sdist(art, version)
        for rule, where in fails:
            print(f"FAIL {art.name}: {rule}: {where}")
        if not fails:
            print(f"ok {art.name}")
        bad = bad or bool(fails)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
