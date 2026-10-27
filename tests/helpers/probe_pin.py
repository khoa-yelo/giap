"""The host the autouse ``_pinned_host_probe`` fixture (tests/conftest.py) pins
``gippa.runner.profiles.probe_host`` to (M11-15a item 1 edge case (d)).

Lead decision on 15a: a FarmShare host, so ``--profile auto`` resolves to ``farmshare,apptainer``
(today's default). Kept here, not in a conftest, so tests import it by a unique module name
(``from conftest import ...`` is ambiguous between tests/conftest.py and tests/unit/conftest.py).
"""

from __future__ import annotations

__all__ = ["PINNED_PROBE"]

#: Keyword arguments of ``gippa.runner.profiles.Probe``.
PINNED_PROBE: dict[str, object] = {
    "on_path": frozenset({"sbatch", "apptainer"}),
    "docker_ok": False,
    "slurm_job_id": None,
    "cluster_name": "farmshare",
    "setup_runtime": None,
}
