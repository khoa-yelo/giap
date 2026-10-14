# ADR-0003: `giapsim.simulate` takes an explicit output directory
- Status: accepted
- Date: 2026-10-14
- Spec sections: §11.3
- Needs Khoa: no (test tooling only; no GIAP interface, default or threshold changes)

## Context
§11.3 gives the entry point as `giapsim.simulate(config: SimConfig, seed: int) -> SimResult` and says
it "writes a directory of inputs plus truth", without naming where the directory comes from.

## Decision
The signature is `simulate(config, seed, outdir) -> SimResult`. The caller (usually a pytest
`tmp_path`) owns the directory, and `SimResult` carries the paths of the files written.

## Evidence (before → after: metrics, tests, simulations)
tests/unit/test_sim_basic.py uses this signature.

## Consequences
None for users: `giapsim` is not shipped in the `giap` wheel.
