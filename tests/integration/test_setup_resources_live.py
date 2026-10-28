"""Live download of the M2-09 resources (docs/specs/M2-09, lead decision 2 after code review).

Each group (``--go``, ``--kegg``, ``--cog``, ``--contaminants``) is downloaded from its official
host with the real ``HttpsFetcher`` into ``tmp_path`` (never $GIPPA_DB_DIR), must verify (``ok``,
pinned sha256 where pinned), and ``--check`` must then pass offline. Marked ``needs_net(<host>)``
for every host the group contacts (GO: purl.obolibrary.org redirects to release.geneontology.org)
and ``slow``; it runs only in the weekly Slurm job. Sizes are printed for PROGRESS.md.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from gippa.setupdb import SOURCES, HttpsFetcher, SetupOptions, read_records, run_setup

pytestmark = pytest.mark.slow

GROUPS = {
    "go": ("go_basic", "goslim_generic"),
    "kegg": ("kegg_ko_pathway", "kegg_pathway_names"),
    "cog": ("cog2020_def", "cog2020_fun"),
    "contaminants": ("contaminants",),
}
EXTRA_HOSTS = {"go": ("release.geneontology.org",)}  # official redirect target (decision 1)


class NoRunner:
    def run(self, tool: str, argv: Sequence[str]) -> int:
        raise AssertionError(f"no tool expected: {tool} {list(argv)}")


class NoFetcher:
    def fetch(self, url: str, dest: Path) -> None:
        raise AssertionError(f"--check must not fetch {url}")


def _hosts(flag: str) -> list[str]:
    hosts = {SOURCES[n].host for n in GROUPS[flag]} | set(EXTRA_HOSTS.get(flag, ()))
    return sorted(hosts)


def _opts(db: Path, flag: str, check: bool = False) -> SetupOptions:
    # M11-15c2 turns every item on by default: state "nothing but <flag>" explicitly (no antiSMASH
    # prerequisite check, no Nextflow launcher, only this resource group)
    groups = {g: g == flag for g in GROUPS}
    return SetupOptions(db, bakta="skip", genomad=False, antismash=False, promotech_models=False,
                        nextflow=False, check=check, **groups)  # fmt: skip


@pytest.mark.parametrize(
    "flag",
    [pytest.param(f, marks=[pytest.mark.needs_net(h) for h in _hosts(f)], id=f) for f in GROUPS],
)
def test_resource_downloads_and_verifies(flag: str, tmp_path: Path) -> None:
    db = tmp_path / "db"
    assert int(run_setup(_opts(db, flag), NoRunner(), HttpsFetcher())) == 0
    recs = read_records(db)
    for name in GROUPS[flag]:
        item, rec = SOURCES[name], recs[name]
        assert rec.status == "ok", rec
        f = db / item.path
        assert f.is_file() and rec.size_bytes == f.stat().st_size > 0, rec
        if item.checksum:
            assert rec.checksum == item.checksum, rec
        else:
            assert rec.checksum.startswith("sha256:"), rec
        print(f"SIZE {name} {rec.size_bytes} bytes ({rec.size_bytes / 1e6:.2f} MB)")
    assert [p for p in db.rglob("*") if ".partial" in p.name] == []
    assert int(run_setup(_opts(db, flag, check=True), NoRunner(), NoFetcher())) == 0
