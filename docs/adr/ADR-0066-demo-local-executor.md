# ADR-0066: `gippa demo` runs the developer-only local executor on bundled layers
- Status: proposed (lead; Khoa approved M11-15 wave 1 for 1.0, Q-K44)
- Date: 2026-10-27
- Spec sections: §5.1 (`gippa run --executor local` documented under the developer guide), §6.4
  (the local executor is developer-only); docs/specs/M11-15-portable-clusters.md item 2, Q-L26
- Needs Khoa: no (covered by the Q-K44 approval of M11-15 wave 1; recorded because it departs
  from §5.1/§6.4)

## Context
GIPPA_SPEC §5.1 documents `gippa run --executor local` "under Developer guide, not in the user
guide" and §6.4 calls the local executor developer-only: it runs every step in-process, needs
precomputed layers for every external tool and has no containers, retries or resume. M11-15 item 2
adds `gippa demo`, the first command a new user runs, which must work with no Nextflow, no
container runtime, no databases and no network. Only the local executor can do that.

## Decision
1. `gippa demo` (without `--full`) runs the bundled SIM100K genome and peptides through the local
   executor, in-process, with every `annotation.precomputed` layer taken from the bundle
   (`gippa/data/demo/layers/`, gunzipped into `OUT/work/demo_layers/`), modules `core,mobile,bgc`,
   no `$GIPPA_DB_DIR` and `NUMBA_CACHE_DIR` under `OUT/work/` unless the user set it.
2. The user guide documents `gippa demo`, never `gippa run --executor local`; the latter stays a
   developer-only option as §5.1 and §6.4 say. The demo exposes no option that reaches the local
   executor's inputs (no `--config`, `--set`, `--genome`).
3. The demo marks itself: the manifest records `demo: true` (an optional boolean of
   `run_manifest.schema.json`) and the report shows "Synthetic demo data (SIM100K, gippasim seed
   20261014)".
4. `gippa demo --full` is the user-facing install check and goes through `gippa run`'s Nextflow
   path (preflight, images, databases, the real tools).

## Evidence (before → after: metrics, tests, simulations)
tests/integration/test_m11_15b_demo_run.py runs the demo with an empty `PATH` and no network and
checks that every published table is byte-identical to an independent
`gippa run --executor local` on the same layers (the §6.4 check pattern). tests/unit/test_m11_15b_demo.py
covers the CLI contract.

## Consequences
- The local executor becomes part of a user command, so a local-executor regression now breaks the
  demo; the integration test and the existing §6.4 local-vs-Nextflow test guard it.
- The bundled layers must stay in sync with gippasim: `scripts/docs/make_demo_layers.py --check`.
