"""``scripts/slurm/build_envs.sbatch`` pull loop on bare refs (M11-09b rule 4, example 10, AC 8).

The script runs with ``bash`` (never ``sbatch``) against a throw-away ``GIAP_REPO``: its own copies
of ``scripts/farmshare_env.sh`` and ``scripts/slurm/build_envs.sbatch``, a stub
``scripts/build_envs.sh`` that exits 0, and a ``containers.config`` written by the test.
``GIAP_ROOT`` is ``tmp_path/root``, so ``farmshare_env.sh`` puts the fake ``apptainer`` in
``$GIAP_ROOT/bin`` first on PATH and ``NXF_APPTAINER_CACHEDIR`` is ``$GIAP_ROOT/containers``. The
fake logs its argv and writes a non-empty target (the loop's presence test is ``[ -s ]``).

A config without pins must pull nothing and still exit 0 (spec edge case "Empty or absent pins":
the job exits with the env-build status, here 0).

The config is the checkout's ``containers.config`` with every ``container = '…'`` value made bare
(the M11-09b form) or prefixed (a legacy config, rule 3), so the test does not depend on which form
the checkout holds.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from giap.annotate.tools import read_containers_config

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "workflow" / "conf" / "containers.config"
PIN_VALUE = re.compile(r"(container\s*=\s*')(?:docker://)?")

FAKE_APPTAINER = r"""#!/bin/bash
echo "apptainer $*" >> "$FAKE_LOG"
if [ "$1" = pull ]; then printf 'fake sif\n' > "$2"; fi
exit 0
"""


def _form(text: str, form: str) -> str:
    return PIN_VALUE.sub(r"\1docker://" if form == "prefixed" else r"\1", text)


def _bare_refs() -> set[str]:
    return {r.removeprefix("docker://") for r in read_containers_config(CONFIG).values()}


def _simple_name(bare: str) -> str:
    """Nextflow 26.10 ``SingularityCache.simpleName`` of ``docker://<bare>``."""
    return bare.replace("/", "-").replace(":", "-") + ".img"


class Box:
    def __init__(self, tmp: Path, form: str) -> None:
        self.repo = tmp / "repo"
        self.root = tmp / "root"
        self.log = tmp / "apptainer.log"
        (self.repo / "scripts" / "slurm").mkdir(parents=True)
        (self.repo / "workflow" / "conf").mkdir(parents=True)
        shutil.copy(REPO / "scripts" / "farmshare_env.sh", self.repo / "scripts")
        shutil.copy(
            REPO / "scripts" / "slurm" / "build_envs.sbatch", self.repo / "scripts" / "slurm"
        )
        (self.repo / "scripts" / "build_envs.sh").write_text("#!/bin/bash\nexit 0\n")
        conf = self.repo / "workflow" / "conf" / "containers.config"
        conf.write_text(_form(CONFIG.read_text(encoding="utf-8"), form), encoding="utf-8")
        bin_dir = self.root / "bin"
        bin_dir.mkdir(parents=True)
        fake = bin_dir / "apptainer"
        fake.write_text(FAKE_APPTAINER)
        fake.chmod(0o755)
        self.log.write_text("")
        self.cache = self.root / "containers"

    def run(self) -> subprocess.CompletedProcess[str]:
        env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG"}}
        env.update(
            GIAP_REPO=str(self.repo),
            GIAP_ROOT=str(self.root),
            FAKE_LOG=str(self.log),
            HOME=str(self.root / "home"),
        )
        proc = subprocess.run(
            ["bash", str(self.repo / "scripts" / "slurm" / "build_envs.sbatch")],
            capture_output=True, text=True, cwd=self.repo, env=env, timeout=120, check=False,
        )  # fmt: skip
        assert proc.returncode == 0, proc.stdout + proc.stderr
        return proc

    def pulls(self) -> list[list[str]]:
        return [
            ln.split()[1:]
            for ln in self.log.read_text().splitlines()
            if ln.split()[1:2] == ["pull"]
        ]


@pytest.mark.parametrize("form", ["bare", "prefixed"])
def test_example10_pulls_each_distinct_image_once_with_a_docker_uri(
    tmp_path: Path, form: str
) -> None:
    box = Box(tmp_path, form)
    box.run()
    refs = _bare_refs()
    assert len(refs) == 6, refs  # BLAST and TBLASTN share one image
    want = sorted(["pull", str(box.cache / _simple_name(r)), f"docker://{r}"] for r in refs)
    assert sorted(box.pulls()) == want
    assert {p.name for p in box.cache.iterdir()} == {_simple_name(r) for r in refs}


@pytest.mark.parametrize("form", ["bare", "prefixed"])
def test_example10_second_run_finds_every_image_present(tmp_path: Path, form: str) -> None:
    box = Box(tmp_path, form)
    box.run()
    assert len(box.pulls()) == 6
    box.log.write_text("")
    proc = box.run()
    assert box.pulls() == []
    present = sorted(ln for ln in proc.stdout.splitlines() if ln.startswith("present: "))
    assert present == sorted(f"present: {box.cache / _simple_name(r)}" for r in _bare_refs())


def test_config_without_pins_pulls_nothing(tmp_path: Path) -> None:
    box = Box(tmp_path, "bare")
    conf = box.repo / "workflow" / "conf" / "containers.config"
    conf.write_text("// no pins\nprocess {\n}\n", encoding="utf-8")
    proc = box.run()
    assert box.pulls() == [] and "present:" not in proc.stdout
